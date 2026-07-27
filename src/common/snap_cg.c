/*
 * snap_cg.c -- ILU(0)-preconditioned CG for Laplacian snap-back.
 *
 * Solves  (L_graph + diag(w)) v = diag(w) target  per coordinate channel.
 * L_graph = D - A (combinatorial graph Laplacian).
 * ILU(0) factorization on the explicit system CSR, then PCG to convergence.
 * All CG internals in double (CLAUDE.md rule 20).
 */
#include "snap_cg.h"
#include "csr.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/* build_system_csr                                                    */
/* Build explicit CSR for M = L_graph + diag(w) from unweighted adj.  */
/* Off-diagonal: -1.0 for each edge. Diagonal: degree_i + w_i.        */
/* Adjacency columns are sorted, so diagonal inserts at correct pos.   */
/* ------------------------------------------------------------------ */
static void build_system_csr(Arena_T arena, const CSR_T adj,
                              const double *w, size_t nv,
                              int32_t **out_off, int32_t **out_col,
                              double **out_val, size_t *out_nnz)
{
    const int32_t *a_off = CSR_offset(adj);
    const int32_t *a_col = CSR_target(adj);
    int32_t adj_nnz = CSR_nnz(adj);
    size_t total_nnz = (size_t)adj_nnz + nv;

    int32_t *off = (int32_t *)ARENA_ALLOC(arena,
                    (long)(nv + 1) * (long)sizeof(int32_t));
    int32_t *col = (int32_t *)ARENA_ALLOC(arena,
                    (long)total_nnz * (long)sizeof(int32_t));
    double  *val = (double *)ARENA_ALLOC(arena,
                    (long)total_nnz * (long)sizeof(double));

    size_t wp = 0;
    for (size_t i = 0; i < nv; i++) {
        off[i] = (int32_t)wp;
        int32_t degree = a_off[i + 1] - a_off[i];
        int diag_done = 0;

        for (int32_t j = a_off[i]; j < a_off[i + 1]; j++) {
            /* Insert diagonal before first neighbor column > i */
            if (!diag_done && a_col[j] > (int32_t)i) {
                col[wp] = (int32_t)i;
                val[wp] = (double)degree + w[i];
                wp++;
                diag_done = 1;
            }
            col[wp] = a_col[j];
            val[wp] = -1.0;
            wp++;
        }
        /* Diagonal at end if all neighbors < i (or isolated vertex) */
        if (!diag_done) {
            col[wp] = (int32_t)i;
            val[wp] = (double)degree + w[i];
            wp++;
        }
    }
    off[nv] = (int32_t)wp;
    assert(wp == total_nnz);

    *out_off = off;
    *out_col = col;
    *out_val = val;
    *out_nnz = total_nnz;
}

/* ------------------------------------------------------------------ */
/* ilu0_factorize                                                      */
/* In-place ILU(0) on CSR values.  Lower triangle stores L (unit diag */
/* implicit), upper triangle stores U (diagonal stored explicitly).     */
/* merge_buf[nv] must be initialized to -1.                            */
/* ------------------------------------------------------------------ */
static void ilu0_factorize(const int32_t *off, const int32_t *col,
                            double *val, size_t nv,
                            int32_t *diag_pos,
                            int32_t *merge_buf)
{
    /* Locate diagonal position in each row */
    for (size_t i = 0; i < nv; i++) {
        diag_pos[i] = -1;
        for (int32_t j = off[i]; j < off[i + 1]; j++) {
            if (col[j] == (int32_t)i) {
                diag_pos[i] = j;
                break;
            }
        }
        assert(diag_pos[i] >= 0);
    }

    /* IKJ ILU(0) factorization */
    for (size_t i = 0; i < nv; i++) {
        /* Mark column positions of row i */
        for (int32_t j = off[i]; j < off[i + 1]; j++)
            merge_buf[col[j]] = j;

        /* Process lower triangle entries (col < i) */
        for (int32_t j = off[i]; j < diag_pos[i]; j++) {
            int32_t k = col[j];
            double dk = val[diag_pos[k]];
            if (fabs(dk) < 1e-15) continue;
            val[j] /= dk;   /* L[i,k] = A[i,k] / U[k,k] */

            /* Subtract from row i using row k's upper triangle */
            for (int32_t p = diag_pos[k] + 1; p < off[k + 1]; p++) {
                int32_t c = col[p];
                if (merge_buf[c] >= 0)
                    val[merge_buf[c]] -= val[j] * val[p];
            }
        }

        /* Unmark */
        for (int32_t j = off[i]; j < off[i + 1]; j++)
            merge_buf[col[j]] = -1;
    }
}

/* ------------------------------------------------------------------ */
/* ilu0_solve — forward/backward substitution: Lz = r, then Uz = z.   */
/* ------------------------------------------------------------------ */
static void ilu0_solve(const int32_t *off, const int32_t *col,
                        const double *val, size_t nv,
                        const int32_t *diag_pos,
                        const double *r, double *z)
{
    /* Forward: L y = r  (unit diagonal implicit) */
    for (size_t i = 0; i < nv; i++) {
        double s = r[i];
        for (int32_t j = off[i]; j < diag_pos[i]; j++)
            s -= val[j] * z[col[j]];
        z[i] = s;
    }

    /* Backward: U z = y */
    for (size_t i = nv; i-- > 0; ) {
        double s = z[i];
        for (int32_t j = diag_pos[i] + 1; j < off[i + 1]; j++)
            s -= val[j] * z[col[j]];
        double d = val[diag_pos[i]];
        z[i] = (fabs(d) > 1e-15) ? s / d : s;
    }
}

/* ------------------------------------------------------------------ */
/* system_matvec — y = M @ x  using explicit CSR (double precision).   */
/* ------------------------------------------------------------------ */
static void system_matvec(const int32_t *off, const int32_t *col,
                           const double *val, size_t nv,
                           const double *x, double *y)
{
    for (size_t i = 0; i < nv; i++) {
        double s = 0.0;
        for (int32_t j = off[i]; j < off[i + 1]; j++)
            s += val[j] * x[col[j]];
        y[i] = s;
    }
}

/* ------------------------------------------------------------------ */
/* pcg_solve_channel                                                   */
/* Single-channel ILU(0)-preconditioned CG.                            */
/* work: 4*nv doubles (r, z, p, Ap).                                   */
/* ------------------------------------------------------------------ */
static int pcg_solve_channel(const int32_t *off, const int32_t *col,
                              const double *sys_val, size_t nv,
                              const int32_t *diag_pos,
                              const double *ilu_val,
                              const double *b, double *x,
                              int max_iter, double tol,
                              double *work,
                              int *out_iters, double *out_rel_res)
{
    double *r  = work;
    double *z  = work + nv;
    double *p  = work + 2 * nv;
    double *Ap = work + 3 * nv;

    /* r = b - M*x */
    system_matvec(off, col, sys_val, nv, x, Ap);
    double bnorm = 0.0;
    for (size_t i = 0; i < nv; i++) {
        r[i] = b[i] - Ap[i];
        bnorm += b[i] * b[i];
    }
    bnorm = sqrt(bnorm);
    if (bnorm < 1e-30) {
        *out_iters = 0;
        *out_rel_res = 0.0;
        return 0;
    }

    /* z = ILU_solve(r);  p = z;  rz = r . z */
    ilu0_solve(off, col, ilu_val, nv, diag_pos, r, z);
    memcpy(p, z, nv * sizeof(double));
    double rz = 0.0;
    for (size_t i = 0; i < nv; i++) rz += r[i] * z[i];

    int iter = 0;
    for (iter = 0; iter < max_iter; iter++) {
        /* Ap = M * p */
        system_matvec(off, col, sys_val, nv, p, Ap);
        double pAp = 0.0;
        for (size_t i = 0; i < nv; i++) pAp += p[i] * Ap[i];
        if (fabs(pAp) < 1e-30) break;
        double alpha = rz / pAp;

        /* x += alpha*p;  r -= alpha*Ap */
        double rnorm = 0.0;
        for (size_t i = 0; i < nv; i++) {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
            rnorm += r[i] * r[i];
        }
        rnorm = sqrt(rnorm);
        if (rnorm / bnorm < tol) {
            *out_iters = iter + 1;
            *out_rel_res = rnorm / bnorm;
            return 0;
        }

        /* z = ILU_solve(r) */
        ilu0_solve(off, col, ilu_val, nv, diag_pos, r, z);
        double rz_new = 0.0;
        for (size_t i = 0; i < nv; i++) rz_new += r[i] * z[i];
        double beta = (fabs(rz) > 1e-30) ? rz_new / rz : 0.0;

        /* p = z + beta*p */
        for (size_t i = 0; i < nv; i++)
            p[i] = z[i] + beta * p[i];
        rz = rz_new;
    }

    /* Did not converge within max_iter */
    double rnorm = 0.0;
    for (size_t i = 0; i < nv; i++) rnorm += r[i] * r[i];
    *out_iters = iter;
    *out_rel_res = sqrt(rnorm) / bnorm;
    return 0;
}

/* ================================================================== */
/* SnapCG_solve -- public entry point                                  */
/* ================================================================== */
int SnapCG_solve(Arena_T arena, const CSR_T adj,
                 const float *w, const float *target,
                 float *verts, size_t nv,
                 int max_iter, double tol)
{
    assert(arena && adj);
    if (nv == 0) return 0;

    Arena_Mark mark = Arena_save(arena);

    /* Convert weights to double */
    double *dw = (double *)ARENA_ALLOC(arena, (long)nv * (long)sizeof(double));
    for (size_t i = 0; i < nv; i++) dw[i] = (double)w[i];

    /* Build explicit system CSR: M = L_graph + diag(w) */
    int32_t *sys_off = NULL, *sys_col = NULL;
    double *sys_val = NULL;
    size_t sys_nnz = 0;
    build_system_csr(arena, adj, dw, nv,
                     &sys_off, &sys_col, &sys_val, &sys_nnz);

    /* Copy values for ILU (original sys_val preserved for matvec) */
    double *ilu_val = (double *)ARENA_ALLOC(arena,
                       (long)sys_nnz * (long)sizeof(double));
    memcpy(ilu_val, sys_val, sys_nnz * sizeof(double));

    /* ILU(0) factorize */
    int32_t *diag_pos = (int32_t *)ARENA_ALLOC(arena,
                          (long)nv * (long)sizeof(int32_t));
    int32_t *merge_buf = (int32_t *)ARENA_ALLOC(arena,
                           (long)nv * (long)sizeof(int32_t));
    for (size_t i = 0; i < nv; i++) merge_buf[i] = -1;

    ilu0_factorize(sys_off, sys_col, ilu_val, nv, diag_pos, merge_buf);

    /* Work vectors: r, z, p, Ap  +  b, x */
    double *work = (double *)ARENA_ALLOC(arena,
                    (long)(4 * nv) * (long)sizeof(double));
    double *b_vec = (double *)ARENA_ALLOC(arena,
                     (long)nv * (long)sizeof(double));
    double *x_vec = (double *)ARENA_ALLOC(arena,
                     (long)nv * (long)sizeof(double));

    static const char *ch_name[3] = { "z", "y", "x" };

    for (int ch = 0; ch < 3; ch++) {
        /* b[i] = w[i] * target[i*3+ch],  x[i] = verts[i*3+ch] */
        for (size_t i = 0; i < nv; i++) {
            b_vec[i] = dw[i] * (double)target[i * 3 + ch];
            x_vec[i] = (double)verts[i * 3 + ch];
        }

        int iters = 0;
        double rel_res = 0.0;
        pcg_solve_channel(sys_off, sys_col, sys_val, nv, diag_pos,
                          ilu_val, b_vec, x_vec,
                          max_iter, tol, work, &iters, &rel_res);

        /* Copy solution back to verts */
        for (size_t i = 0; i < nv; i++)
            verts[i * 3 + ch] = (float)x_vec[i];

        fprintf(stderr, "      CG ch=%s: %d iters, rel_res=%.2e\n",
                ch_name[ch], iters, rel_res);
    }

    Arena_restore(arena, mark);
    return 0;
}
