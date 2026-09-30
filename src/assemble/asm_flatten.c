/* asm_flatten.c -- intrinsic flattening of one mesh chart (Tutte + ARAP). */
#include "asm_flatten.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../common/csr.h"
#include "../common/eig3.h"
#include "../common/pipeline_constants.h"
#include "../common/uv_guard.h"
#include "../common/ves_platform.h"
#include "../flatten/slim_refine.h"
#include "../flatten/sparse_solve.h"
#include "asm_contacts.h"
#include "asm_metric.h"

#define AF_PI 3.14159265358979323846

/* Wall-clock deadline of this thread's flattening (0 = none): a solve that
 * meets it fails, and a rescue keeps what it certified so far. */
static VES_THREAD_LOCAL double af_deadline_at = 0;
static int af_expired(void) { return af_deadline_at > 0 && ves_clock_sec() >= af_deadline_at; }

/* Set while this thread runs the crumpled rescue: its linear solves factor once with TAUCS
 * (af_solver_*). Everywhere else the historical Jacobi PCG runs, so every chart the clean stage
 * flattens keeps its map bit for bit. */
static VES_THREAD_LOCAL int af_direct = 0;

/* ---- small dense helpers ------------------------------------------------ */

static double af_dot3(const double *a, const double *b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static double af_len3(const double *a) { return sqrt(af_dot3(a, a)); }

static void af_sub3(const double *a, const double *b, double *o)
{
    o[0] = a[0]-b[0]; o[1] = a[1]-b[1]; o[2] = a[2]-b[2];
}

static void af_get3(const float *xyz, int32_t i, double *o)
{
    o[0] = (double)xyz[(size_t)i*3+0];
    o[1] = (double)xyz[(size_t)i*3+1];
    o[2] = (double)xyz[(size_t)i*3+2];
}

/* Angle at vertex a of triangle (a,b,c), radians, clamped to (0, pi). */
static double af_angle(const double *pa, const double *pb, const double *pc)
{
    double u[3], v[3];
    af_sub3(pb, pa, u);
    af_sub3(pc, pa, v);
    double lu = af_len3(u), lv = af_len3(v);
    if (lu <= 0.0 || lv <= 0.0) return AF_PI / 3.0;
    double c = af_dot3(u, v) / (lu * lv);
    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;
    return acos(c);
}

/* ---- sparse symmetric matrix on the CSR_from_faces structure ------------- */

typedef struct AfMat {
    size_t         n;
    const int32_t *off;   /* [n+1] */
    const int32_t *col;   /* [nnz] */
    double        *val;   /* [nnz] off-diagonal values */
    double        *diag;  /* [n] */
} AfMat;

static int32_t af_slot(const AfMat *m, int32_t i, int32_t j)
{
    for (int32_t k = m->off[i]; k < m->off[i+1]; k++)
        if (m->col[k] == j) return k;
    return -1;
}

static void af_matvec(const AfMat *m, const double *x, double *y)
{
    for (size_t i = 0; i < m->n; i++) {
        double s = m->diag[i] * x[i];
        for (int32_t k = m->off[i]; k < m->off[i+1]; k++)
            s += m->val[k] * x[(size_t)m->col[k]];
        y[i] = s;
    }
}

/* Jacobi-preconditioned CG.  x holds the initial guess on entry.  Returns 0
 * on convergence, -1 otherwise. */
static int af_pcg(const AfMat *m, const double *b, double *x,
                  double *r, double *z, double *p, double *q,
                  double tol, int maxit)
{
    size_t n = m->n;
    af_matvec(m, x, q);
    double bnorm = 0.0, rz = 0.0;
    for (size_t i = 0; i < n; i++) {
        r[i] = b[i] - q[i];
        bnorm += b[i]*b[i];
        double d = m->diag[i];
        z[i] = (fabs(d) > 1e-300) ? r[i] / d : r[i];
        p[i] = z[i];
        rz += r[i]*z[i];
    }
    bnorm = sqrt(bnorm);
    if (bnorm <= 0.0) bnorm = 1.0;
    double rnorm = 0.0;
    for (size_t i = 0; i < n; i++) rnorm += r[i]*r[i];
    if (sqrt(rnorm) <= tol * bnorm) return 0;
    for (int it = 0; it < maxit; it++) {
        if (!(it & 63) && af_expired()) return -1;
        af_matvec(m, p, q);
        double pq = 0.0;
        for (size_t i = 0; i < n; i++) pq += p[i]*q[i];
        if (!(pq > 1e-300)) return -1;
        double alpha = rz / pq;
        rnorm = 0.0;
        for (size_t i = 0; i < n; i++) {
            x[i] += alpha * p[i];
            r[i] -= alpha * q[i];
            rnorm += r[i]*r[i];
        }
        if (sqrt(rnorm) <= tol * bnorm) return 0;
        double rz_new = 0.0;
        for (size_t i = 0; i < n; i++) {
            double d = m->diag[i];
            z[i] = (fabs(d) > 1e-300) ? r[i] / d : r[i];
            rz_new += r[i]*z[i];
        }
        double beta = rz_new / rz;
        rz = rz_new;
        for (size_t i = 0; i < n; i++) p[i] = z[i] + beta * p[i];
    }
    return -1;
}

/* Linear solves of one matrix against one or more right-hand sides (2026-09-28). Profiling the
 * PHerc0826 10^3 assembly put 30 of 45 sampled chart-rescue thread stacks in af_pcg/af_matvec: the
 * Tutte, ARAP and LSCM solves of large rescue pieces, ARAP re-solving its fixed Laplacian every
 * iteration. Under af_direct a symmetric matrix is factored once (supernodal Cholesky, the
 * sparse_solve TAUCS path) and each solve is two triangular sweeps; the solution is accepted only
 * when its residual meets the PCG tolerance, else PCG continues from it. A non-symmetric matrix
 * (mean-value Tutte weights), a failed factorization and every solve outside af_direct use the
 * historical Jacobi PCG. */
typedef struct AfSolver { const AfMat *m; SparseFactor_T f; } AfSolver;

static void af_solver_init(Arena_T arena, AfSolver *s, const AfMat *m)
{
    s->m = m; s->f = NULL;
    size_t n = m->n;
    if (!af_direct || n < 2 || n > (size_t)INT32_MAX / 4) return;
    Arena_Mark mark = Arena_save(arena);
    size_t cap = n + (size_t)m->off[n], nt = 0;
    int *rows = ARENA_ALLOC(arena, cap * sizeof *rows), *cols = ARENA_ALLOC(arena, cap * sizeof *cols);
    double *vals = ARENA_ALLOC(arena, cap * sizeof *vals);
    int symmetric = 1;
    for (size_t i = 0; i < n && symmetric; i++) {
        rows[nt] = (int)i; cols[nt] = (int)i; vals[nt] = m->diag[i]; nt++;
        for (int32_t k = m->off[i]; k < m->off[i+1]; k++) {
            int32_t j = m->col[k];
            if (j >= (int32_t)i) continue;                 /* lower triangle, row > col */
            double v = m->val[k];
            int32_t kt = af_slot(m, j, (int32_t)i);
            double t = kt >= 0 ? m->val[kt] : 0.0;
            if (fabs(v - t) > 1e-12 * fmax(fabs(v), fabs(t))) { symmetric = 0; break; }
            if (fabs(v) > 0.0) { rows[nt] = (int)i; cols[nt] = (int)j; vals[nt] = v; nt++; }
        }
    }
    if (symmetric && Sparse_factor_spd((int)n, (int)nt, rows, cols, vals, &s->f) != 0) s->f = NULL;
    Arena_restore(arena, mark);
}

static int af_solver_solve(AfSolver *s, const double *b, double *x,
                           double *r, double *z, double *p, double *q, double tol, int maxit)
{
    if (s->f && Sparse_factor_solve(s->f, b, x) == 0) {
        af_matvec(s->m, x, q);
        double rn = 0.0, bn = 0.0;
        for (size_t i = 0; i < s->m->n; i++) { double d = b[i] - q[i]; rn += d*d; bn += b[i]*b[i]; }
        if (isfinite(rn) && sqrt(rn) <= tol * (bn > 0.0 ? sqrt(bn) : 1.0)) return 0;
    }
    return af_pcg(s->m, b, x, r, z, p, q, tol, maxit);
}

static void af_solver_free(AfSolver *s)
{
    if (s->f) Sparse_factor_free(&s->f);
}

/* ---- boundary structure -------------------------------------------------- */

typedef struct AfTopo {
    size_t   nv, nf;
    CSR_T    adj;
    int32_t *edge_count;   /* [nnz] faces incident to the slot's edge */
    int32_t *bnext;        /* [nv] next boundary vertex, -1 = none */
    uint8_t *on_boundary;  /* [nv] */
    int32_t *loop_id;      /* [nv] boundary loop id or -1 */
    size_t   n_loops;
    int32_t *loop_start;   /* [n_loops] a vertex of each loop */
    size_t  *loop_len;     /* [n_loops] vertex count */
    double  *loop_len3d;   /* [n_loops] 3-D length */
    int32_t  outer;        /* longest loop */
} AfTopo;

/* Returns 0 or an ASM_FLAT_* failure. */
static int af_topology(Arena_T arena, const float *xyz, size_t nv,
                       const int32_t *faces, size_t nf, AfTopo *t)
{
    memset(t, 0, sizeof *t);
    t->nv = nv; t->nf = nf;
    t->adj = CSR_from_faces(arena, faces, nf, nv);
    const int32_t *off = CSR_offset(t->adj);
    const int32_t *col = CSR_target(t->adj);
    size_t nnz = (size_t)CSR_nnz(t->adj);
    t->edge_count = ARENA_CALLOC(arena, nnz, sizeof(int32_t));
    AfMat probe = { nv, off, col, NULL, NULL };
    for (size_t f = 0; f < nf; f++) {
        for (int e = 0; e < 3; e++) {
            int32_t a = faces[f*3 + (size_t)e];
            int32_t b = faces[f*3 + (size_t)((e+1)%3)];
            int32_t ka = af_slot(&probe, a, b);
            int32_t kb = af_slot(&probe, b, a);
            if (ka < 0 || kb < 0) return ASM_FLAT_BAD_INDEX;
            t->edge_count[ka]++;
            t->edge_count[kb]++;
        }
    }
    for (size_t k = 0; k < nnz; k++)
        if (t->edge_count[k] >= 3) return ASM_FLAT_NONMANIFOLD;

    t->bnext = ARENA_ALLOC(arena, nv * sizeof(int32_t));
    t->on_boundary = ARENA_CALLOC(arena, nv, sizeof(uint8_t));
    t->loop_id = ARENA_ALLOC(arena, nv * sizeof(int32_t));
    for (size_t i = 0; i < nv; i++) { t->bnext[i] = -1; t->loop_id[i] = -1; }
    size_t n_bedges = 0;
    for (size_t f = 0; f < nf; f++) {
        for (int e = 0; e < 3; e++) {
            int32_t a = faces[f*3 + (size_t)e];
            int32_t b = faces[f*3 + (size_t)((e+1)%3)];
            int32_t ka = af_slot(&probe, a, b);
            if (t->edge_count[ka] == 1) {
                if (t->bnext[a] != -1) return ASM_FLAT_NONMANIFOLD; /* pinch */
                t->bnext[a] = b;
                t->on_boundary[a] = 1;
                t->on_boundary[b] = 1;
                n_bedges++;
            }
        }
    }
    if (n_bedges == 0) return ASM_FLAT_NO_BOUNDARY;
    /* every boundary vertex must have an outgoing halfedge */
    for (size_t i = 0; i < nv; i++)
        if (t->on_boundary[i] && t->bnext[i] < 0) return ASM_FLAT_NONMANIFOLD;

    /* walk loops */
    size_t cap = 16;
    t->loop_start = ARENA_ALLOC(arena, cap * sizeof(int32_t));
    t->loop_len = ARENA_ALLOC(arena, cap * sizeof(size_t));
    t->loop_len3d = ARENA_ALLOC(arena, cap * sizeof(double));
    for (size_t i = 0; i < nv; i++) {
        if (!t->on_boundary[i] || t->loop_id[i] >= 0) continue;
        if (t->n_loops == cap) {
            size_t ncap = cap * 2;
            int32_t *ls = ARENA_ALLOC(arena, ncap * sizeof(int32_t));
            size_t *ll = ARENA_ALLOC(arena, ncap * sizeof(size_t));
            double *l3 = ARENA_ALLOC(arena, ncap * sizeof(double));
            memcpy(ls, t->loop_start, cap * sizeof(int32_t));
            memcpy(ll, t->loop_len, cap * sizeof(size_t));
            memcpy(l3, t->loop_len3d, cap * sizeof(double));
            t->loop_start = ls; t->loop_len = ll; t->loop_len3d = l3;
            cap = ncap;
        }
        int32_t id = (int32_t)t->n_loops;
        int32_t v = (int32_t)i;
        size_t len = 0;
        double len3 = 0.0;
        while (t->loop_id[v] < 0) {
            t->loop_id[v] = id;
            int32_t w = t->bnext[v];
            double pv[3], pw[3], d[3];
            af_get3(xyz, v, pv);
            af_get3(xyz, w, pw);
            af_sub3(pw, pv, d);
            len3 += af_len3(d);
            len++;
            v = w;
            if (len > nv) return ASM_FLAT_NONMANIFOLD;
        }
        if (v != (int32_t)i) return ASM_FLAT_NONMANIFOLD; /* did not close on its start */
        t->loop_start[t->n_loops] = (int32_t)i;
        t->loop_len[t->n_loops] = len;
        t->loop_len3d[t->n_loops] = len3;
        t->n_loops++;
    }
    t->outer = 0;
    for (size_t l = 1; l < t->n_loops; l++)
        if (t->loop_len3d[l] > t->loop_len3d[(size_t)t->outer]) t->outer = (int32_t)l;
    return ASM_FLAT_OK;
}

/* ---- Tutte embedding ------------------------------------------------------ */

/* Builds the capped Laplacian (real vertices + one virtual vertex per hole),
 * pins the outer loop on a circle, solves u and v.  uv is [nv*2]. */
static int af_tutte(Arena_T arena, const float *xyz, size_t nv,
                    const int32_t *faces, size_t nf, const AfTopo *t,
                    double area3d, double *uv, size_t *out_capped, int cap_mvc)
{
    const int32_t *off = CSR_offset(t->adj);
    const int32_t *col = CSR_target(t->adj);
    size_t nnz = (size_t)CSR_nnz(t->adj);
    size_t n_holes = t->n_loops - 1;
    size_t n = nv + n_holes;

    /* symmetric mean-value weights on the real edges */
    double *w = ARENA_CALLOC(arena, nnz, sizeof(double));
    AfMat probe = { nv, off, col, NULL, NULL };
    double mean_edge = 0.0;
    size_t n_edge = 0;
    for (size_t f = 0; f < nf; f++) {
        int32_t a = faces[f*3], b = faces[f*3+1], c = faces[f*3+2];
        double pa[3], pb[3], pc[3];
        af_get3(xyz, a, pa); af_get3(xyz, b, pb); af_get3(xyz, c, pc);
        double ang_a = af_angle(pa, pb, pc);
        double ang_b = af_angle(pb, pc, pa);
        double ang_c = af_angle(pc, pa, pb);
        double d[3];
        af_sub3(pb, pa, d); double lab = af_len3(d);
        af_sub3(pc, pb, d); double lbc = af_len3(d);
        af_sub3(pa, pc, d); double lca = af_len3(d);
        if (lab <= 0.0 || lbc <= 0.0 || lca <= 0.0) continue;
        mean_edge += lab + lbc + lca; n_edge += 3;
        /* tan(angle/2) at each corner contributes to both edges at that corner */
        double ta = tan(0.5 * ang_a), tb = tan(0.5 * ang_b), tc = tan(0.5 * ang_c);
        if (ta > 50.0) ta = 50.0; if (tb > 50.0) tb = 50.0; if (tc > 50.0) tc = 50.0;
        w[af_slot(&probe, a, b)] += ta / lab;  /* corner a, edge a-b */
        w[af_slot(&probe, a, c)] += ta / lca;  /* corner a, edge a-c */
        w[af_slot(&probe, b, a)] += tb / lab;
        w[af_slot(&probe, b, c)] += tb / lbc;
        w[af_slot(&probe, c, b)] += tc / lbc;
        w[af_slot(&probe, c, a)] += tc / lca;
    }
    if (n_edge == 0) return ASM_FLAT_DEGENERATE;
    /* Mean-value hole caps (cap_mvc): each hole is filled for the embedding only by a fan from a
     * virtual vertex at the loop's 3-D centroid, and the fan triangles are weighted by the same
     * mean-value formula as the real faces, including their share of the loop edges.  The uniform
     * cap weight below pulls every loop vertex to the hub with the strength of a real edge, so a
     * loop of hundreds of vertices collapses toward a point and its faces fold (PHerc0343: 1e-5
     * vox faces, thousands of flipped seeds); geometric weights keep a big hole big. */
    double *capw = NULL;
    if (cap_mvc) {
        capw = ARENA_CALLOC(arena, nv, sizeof(double));
        double *capc = ARENA_CALLOC(arena, nv, sizeof(double));   /* hub-side weight per loop vertex */
        for (size_t l = 0; l < t->n_loops; l++) {
            if ((int32_t)l == t->outer) continue;
            double c[3] = {0, 0, 0};
            int32_t v = t->loop_start[l];
            for (size_t k = 0; k < t->loop_len[l]; k++) {
                double pv[3]; af_get3(xyz, v, pv);
                for (int d = 0; d < 3; d++) c[d] += pv[d] / (double)t->loop_len[l];
                v = t->bnext[v];
            }
            v = t->loop_start[l];
            for (size_t k = 0; k < t->loop_len[l]; k++) {
                int32_t nx = t->bnext[v];
                double pv[3], pn[3], d[3];
                af_get3(xyz, v, pv); af_get3(xyz, nx, pn);
                af_sub3(pn, pv, d); double lvn = af_len3(d);
                af_sub3(c, pv, d); double lvc = af_len3(d);
                af_sub3(c, pn, d); double lnc = af_len3(d);
                if (lvn > 0.0 && lvc > 0.0 && lnc > 0.0) {
                    double tv = tan(0.5 * af_angle(pv, pn, c)), tn = tan(0.5 * af_angle(pn, pv, c));
                    double tc = tan(0.5 * af_angle(c, pv, pn));
                    if (tv > 50.0) tv = 50.0; if (tn > 50.0) tn = 50.0; if (tc > 50.0) tc = 50.0;
                    w[af_slot(&probe, v, nx)] += tv / lvn;
                    w[af_slot(&probe, nx, v)] += tn / lvn;
                    capw[v] += tv / lvc; capw[nx] += tn / lnc;
                    capc[v] += tc / lvc; capc[nx] += tc / lnc;
                }
                v = nx;
            }
        }
        for (size_t i = 0; i < nv; i++) {
            capw[i] = 0.5 * (capw[i] + capc[i]);
            if (t->on_boundary[i] && t->loop_id[i] != t->outer && capw[i] < 1e-9) capw[i] = 1e-9;
        }
    }
    mean_edge /= (double)n_edge;
    /* symmetrize */
    for (int32_t i = 0; i < (int32_t)nv; i++) {
        for (int32_t k = off[i]; k < off[i+1]; k++) {
            int32_t j = col[k];
            if (j < i) continue;
            int32_t kj = af_slot(&probe, j, i);
            double s = 0.5 * (w[k] + w[kj]);
            if (s < 1e-9) s = 1e-9;
            w[k] = s; w[kj] = s;
        }
    }

    /* extended CSR: real rows get their slots plus (maybe) one cap slot;
     * cap rows get the hole loop. */
    int32_t *eoff = ARENA_ALLOC(arena, (n + 1) * sizeof(int32_t));
    size_t enz = nnz;
    for (size_t i = 0; i < nv; i++)
        if (t->on_boundary[i] && t->loop_id[i] != t->outer) enz++;
    for (size_t l = 0; l < t->n_loops; l++)
        if ((int32_t)l != t->outer) enz += t->loop_len[l];
    int32_t *ecol = ARENA_ALLOC(arena, enz * sizeof(int32_t));
    double *eval = ARENA_ALLOC(arena, enz * sizeof(double));
    double *ediag = ARENA_CALLOC(arena, n, sizeof(double));
    double wcap = 1.0 / mean_edge;
    /* cap index for hole loop l: nv + rank among holes */
    int32_t *cap_of_loop = ARENA_ALLOC(arena, t->n_loops * sizeof(int32_t));
    {
        int32_t r = 0;
        for (size_t l = 0; l < t->n_loops; l++)
            cap_of_loop[l] = ((int32_t)l == t->outer) ? -1 : (int32_t)nv + r++;
    }
    size_t pos = 0;
    for (size_t i = 0; i < nv; i++) {
        eoff[i] = (int32_t)pos;
        for (int32_t k = off[i]; k < off[i+1]; k++) {
            ecol[pos] = col[k]; eval[pos] = -w[k]; ediag[i] += w[k]; pos++;
        }
        if (t->on_boundary[i] && t->loop_id[i] != t->outer) {
            double wc = capw ? capw[i] : wcap;
            ecol[pos] = cap_of_loop[t->loop_id[i]]; eval[pos] = -wc; ediag[i] += wc; pos++;
        }
    }
    for (size_t l = 0; l < t->n_loops; l++) {
        if ((int32_t)l == t->outer) continue;
        size_t row = (size_t)cap_of_loop[l];
        eoff[row] = (int32_t)pos;
        int32_t v = t->loop_start[l];
        for (size_t s = 0; s < t->loop_len[l]; s++) {
            double wc = capw ? capw[v] : wcap;
            ecol[pos] = v; eval[pos] = -wc; ediag[row] += wc; pos++;
            v = t->bnext[v];
        }
    }
    eoff[n] = (int32_t)pos;
    AfMat M = { n, eoff, ecol, eval, ediag };

    /* Dirichlet on the outer loop: circle of radius sqrt(area/pi), arc-length
     * parameterized, counter-clockwise in the loop's own order. */
    double *bu = ARENA_CALLOC(arena, n, sizeof(double));
    double *bv = ARENA_CALLOC(arena, n, sizeof(double));
    double *xu = ARENA_CALLOC(arena, n, sizeof(double));
    double *xv = ARENA_CALLOC(arena, n, sizeof(double));
    uint8_t *fixed = ARENA_CALLOC(arena, n, sizeof(uint8_t));
    double R = sqrt(area3d / AF_PI);
    if (!(R > 0.0)) return ASM_FLAT_DEGENERATE;
    {
        int32_t v = t->loop_start[(size_t)t->outer];
        double L = t->loop_len3d[(size_t)t->outer];
        double s = 0.0;
        for (size_t k = 0; k < t->loop_len[(size_t)t->outer]; k++) {
            double ang = 2.0 * AF_PI * s / L;
            xu[v] = R * cos(ang);
            xv[v] = R * sin(ang);
            fixed[v] = 1;
            int32_t wv = t->bnext[v];
            double pv[3], pw[3], d[3];
            af_get3(xyz, v, pv); af_get3(xyz, wv, pw); af_sub3(pw, pv, d);
            s += af_len3(d);
            v = wv;
        }
    }
    /* move fixed columns to the RHS, make fixed rows identity */
    for (size_t i = 0; i < n; i++) {
        if (fixed[i]) {
            ediag[i] = 1.0;
            for (int32_t k = eoff[i]; k < eoff[i+1]; k++) eval[k] = 0.0;
            bu[i] = xu[i]; bv[i] = xv[i];
            continue;
        }
        for (int32_t k = eoff[i]; k < eoff[i+1]; k++) {
            int32_t j = ecol[k];
            if (fixed[j]) {
                bu[i] -= eval[k] * xu[j];
                bv[i] -= eval[k] * xv[j];
                eval[k] = 0.0;
            }
        }
    }
    double *r = ARENA_ALLOC(arena, n * sizeof(double));
    double *z = ARENA_ALLOC(arena, n * sizeof(double));
    double *p = ARENA_ALLOC(arena, n * sizeof(double));
    double *q = ARENA_ALLOC(arena, n * sizeof(double));
    int maxit = (int)(4 * n + 200);
    if (maxit > 20000) maxit = 20000;
    AfSolver solver; af_solver_init(arena, &solver, &M);
    int solved = af_solver_solve(&solver, bu, xu, r, z, p, q, 1e-10, maxit) == 0 &&
                 af_solver_solve(&solver, bv, xv, r, z, p, q, 1e-10, maxit) == 0;
    af_solver_free(&solver);
    if (!solved) return ASM_FLAT_SOLVE;
    for (size_t i = 0; i < nv; i++) { uv[i*2] = xu[i]; uv[i*2+1] = xv[i]; }
    *out_capped = n_holes;
    return ASM_FLAT_OK;
}

/* ---- ARAP local/global -------------------------------------------------- */

typedef struct AfRef {
    double x1, y2, x2;   /* reference frame: p0=(0,0) p1=(x1,0) p2=(x2,y2) */
    double area;
    double cot[3];       /* cot of the angle opposite edge e (e = (v_e, v_{e+1})) */
} AfRef;

static void af_reference(const float *xyz, const int32_t *faces, size_t nf, AfRef *ref)
{
    for (size_t f = 0; f < nf; f++) {
        double p0[3], p1[3], p2[3], e1[3], e2[3];
        af_get3(xyz, faces[f*3], p0);
        af_get3(xyz, faces[f*3+1], p1);
        af_get3(xyz, faces[f*3+2], p2);
        af_sub3(p1, p0, e1);
        af_sub3(p2, p0, e2);
        double l1 = af_len3(e1);
        double d = (l1 > 0.0) ? af_dot3(e2, e1) / l1 : 0.0;
        double cross[3] = {e1[1]*e2[2]-e1[2]*e2[1],e1[2]*e2[0]-e1[0]*e2[2],e1[0]*e2[1]-e1[1]*e2[0]};
        double twice_area = af_len3(cross);
        double h = l1 > 0 ? twice_area/l1 : 0;
        ref[f].x1 = l1; ref[f].x2 = d; ref[f].y2 = h;
        ref[f].area = 0.5 * l1 * h;
        /* All three cotangents have one area denominator. Cap them with
         * one positive triangle weight: clipping each edge independently
         * can turn an obtuse face's positive energy into a negative one.
         * Dot/cross products also retain thin faces whose acos or squared
         * height subtraction rounds away their nonzero area. */
        double e12[3]; af_sub3(p2,p1,e12);
        double c0 = af_dot3(e2,e12), c1 = af_dot3(e1,e2), c2 = -af_dot3(e1,e12);
        double largest = fmax(fabs(c0),fmax(fabs(c1),fabs(c2)));
        double weight = twice_area > 0 && largest > 0 ? fmin(1/twice_area,20/largest) : 0;
        ref[f].cot[0] = weight*c0; ref[f].cot[1] = weight*c1; ref[f].cot[2] = weight*c2;
    }
}

/* 2x2 Jacobian of the map reference->uv for face f. */
static void af_jacobian(const AfRef *rf, const double *uv, const int32_t *fv, double J[4])
{
    double u0 = uv[(size_t)fv[0]*2], v0 = uv[(size_t)fv[0]*2+1];
    double u1 = uv[(size_t)fv[1]*2], v1 = uv[(size_t)fv[1]*2+1];
    double u2 = uv[(size_t)fv[2]*2], v2 = uv[(size_t)fv[2]*2+1];
    /* D = [uv1-uv0, uv2-uv0], X = [x1-x0, x2-x0] = [[x1, x2],[0, y2]]; J = D X^-1 */
    double du1 = u1-u0, dv1 = v1-v0, du2 = u2-u0, dv2 = v2-v0;
    double x1 = rf->x1, x2 = rf->x2, y2 = rf->y2;
    double det = x1 * y2;
    if (!(x1 > 0) || !(y2 > 0) || !isfinite(det)) {
        J[0] = J[1] = J[2] = J[3] = NAN; return;
    }
    /* X^-1 = 1/det * [[y2, -x2],[0, x1]] */
    J[0] = du1/x1;                          /* dU/dx */
    J[1] = (du2-J[0]*x2)/y2;               /* dU/dy */
    J[2] = dv1/x1;                          /* dV/dx */
    J[3] = (dv2-J[2]*x2)/y2;               /* dV/dy */
}

static void af_singular(const double J[4], double *lo, double *hi, double *det)
{
    double a = J[0], b = J[1], c = J[2], d = J[3];
    double s1 = sqrt((a+d)*(a+d) + (b-c)*(b-c));
    double s2 = sqrt((a-d)*(a-d) + (b+c)*(b+c));
    *hi = 0.5 * (s1 + s2);
    *det = a*d - b*c;
    *lo = *hi > 0 ? fabs(*det)/(*hi) : 0;
}

/* closest rotation to J (polar factor), always det +1 */
static void af_polar(const double J[4], double R[4])
{
    double p = J[0] + J[3], q = J[2] - J[1];
    double n = sqrt(p*p + q*q);
    if (n < 1e-12) { R[0] = 1; R[1] = 0; R[2] = 0; R[3] = 1; return; }
    p /= n; q /= n;
    R[0] = p; R[1] = -q; R[2] = q; R[3] = p;
}

static double af_arap_energy(const AfRef *ref, const int32_t *faces, size_t nf,
                              const double *uv)
{
    double energy = 0;
    for (size_t f = 0; f < nf; f++) {
        double J[4], R[4]; af_jacobian(ref+f,uv,faces+3*f,J); af_polar(J,R);
        double px[3] = {0,ref[f].x1,ref[f].x2}, py[3] = {0,0,ref[f].y2};
        for (size_t k = 0; k < 3; k++) {
            size_t j = (k+1)%3, a = (size_t)faces[3*f+k], b = (size_t)faces[3*f+j];
            double x = px[k]-px[j], y = py[k]-py[j];
            double du = uv[2*a]-uv[2*b]-R[0]*x-R[1]*y;
            double dv = uv[2*a+1]-uv[2*b+1]-R[2]*x-R[3]*y;
            energy += ref[f].cot[k]*(du*du+dv*dv);
        }
    }
    return energy;
}

static void af_normalize_frame(double *uv, size_t nv, const int32_t *faces, size_t nf, int mirror);

static int af_encoded_feasible(const double *uv, double *scratch, size_t nv, const int32_t *faces, size_t nf)
{
    memcpy(scratch,uv,2*nv*sizeof(double));
    for (size_t k = 0; k < 2*nv; k++) scratch[k] = (float)scratch[k];
    return UvGuard_mesh(scratch,faces,nf);
}

typedef struct AfBoundaryQuery {
    AsmContacts *query;
    const int32_t *faces;
    uint8_t *blocked;
} AfBoundaryQuery;

/* Query actual triangles incident to every boundary loop. This distinguishes
 * touching edges from material overlapping on the same side of an edge,
 * including collinear boundaries and intersections exactly at a vertex. The
 * numerical tolerance is identical to the complete field audit. */
static AfBoundaryQuery af_boundary_query(Arena_T arena, const AfTopo *topo,
                                         const int32_t *faces, size_t nf, const double *uv)
{
    size_t count = 0;
    for (size_t f = 0; f < nf; f++) for (int k = 0; k < 3; k++)
        if (topo->on_boundary[faces[3*f+k]]) { count++; break; }
    int32_t *boundary = ARENA_ALLOC(arena,3*count*sizeof(*boundary)); count = 0;
    for (size_t f = 0; f < nf; f++) for (int k = 0; k < 3; k++)
        if (topo->on_boundary[faces[3*f+k]]) { memcpy(boundary+3*count++,faces+3*f,3*sizeof(*faces)); break; }
    AfBoundaryQuery result = {AsmContacts_new(arena,topo->nv,boundary,count,uv,1e-7),boundary,NULL};
    return result;
}
static void af_boundary_block(void *context, size_t a, size_t b, double area)
{
    AfBoundaryQuery *q = context; (void)area;
    if (q->blocked) for (int k = 0; k < 3; k++)
        q->blocked[q->faces[3*a+k]] = q->blocked[q->faces[3*b+k]] = 1;
}
static int af_boundary_clear(const AfTopo *topo, const double *uv, AfBoundaryQuery *query, uint8_t *blocked)
{
    AsmContactStats stats; (void)topo; query->blocked = blocked;
    return query->query && !AsmContacts_measure(query->query,uv,NULL,af_boundary_block,query,&stats) && !stats.pairs;
}

/* A single limiting triangle must not freeze the unconstrained solve of an
 * entire chart. Damp its incident vertices, propagate through every affected
 * face, then send the resulting direction through the same objective and
 * continuous-orientation line search. This does not edit or remove faces. */
static int af_local_direction(const double *uv, const double *raw, double *target, double *encoded,
                               double *alpha, uint8_t *blocked, size_t nv, const int32_t *faces, size_t nf, int encoded_ready,
                               const AfTopo *topo, AfBoundaryQuery *boxes)
{
    for (size_t v = 0; v < nv; v++) alpha[v] = 1;
    for (int pass = 0; pass < 48; pass++) {
        for (size_t v = 0; v < nv; v++) for (int d = 0; d < 2; d++) target[2*v+d] = uv[2*v+d]+alpha[v]*(raw[2*v+d]-uv[2*v+d]);
        if (encoded_ready) af_encoded_feasible(target,encoded,nv,faces,nf);
        memset(blocked,0,nv); size_t bad = 0;
        for (size_t f = 0; f < nf; f++) if (!UvGuard_interval(uv,target,faces+3*f,1) ||
            (encoded_ready && !UvGuard_positive(encoded,faces+3*f))) {
            bad++; for (int k = 0; k < 3; k++) blocked[faces[3*f+k]] = 1;
        }
        if (!af_boundary_clear(topo,target,boxes,blocked)) bad++;
        if (encoded_ready && !af_boundary_clear(topo,encoded,boxes,blocked)) bad++;
        if (!bad) return 1;
        for (size_t v = 0; v < nv; v++) if (blocked[v]) alpha[v] *= .5;
    }
    return 0;
}

static int af_arap(Arena_T arena, const float *xyz, size_t nv,
                   const int32_t *faces, size_t nf, const AfTopo *t,
                   double *uv, int max_iters, int *out_iters)
{
    const int32_t *off = CSR_offset(t->adj);
    const int32_t *col = CSR_target(t->adj);
    size_t nnz = (size_t)CSR_nnz(t->adj);
    AfRef *ref = ARENA_ALLOC(arena, nf * sizeof(AfRef));
    af_reference(xyz, faces, nf, ref);
    if (!UvGuard_mesh(uv,faces,nf)) return ASM_FLAT_FOLDED;
    double *target = ARENA_ALLOC(arena,2*nv*sizeof(double));
    double *trial = ARENA_ALLOC(arena,2*nv*sizeof(double));
    double *encoded = ARENA_ALLOC(arena,2*nv*sizeof(double));
    double *raw_target = ARENA_ALLOC(arena,2*nv*sizeof(double));
    double *local_alpha = ARENA_ALLOC(arena,nv*sizeof(double));
    uint8_t *blocked = ARENA_ALLOC(arena,nv);
    AfBoundaryQuery boundary = af_boundary_query(arena,t,faces,nf,uv), *boxes = &boundary;
    if (!af_boundary_clear(t,uv,boxes,NULL)) return ASM_FLAT_FOLDED;
    int encoded_ready = af_encoded_feasible(uv,encoded,nv,faces,nf) && af_boundary_clear(t,encoded,boxes,NULL);
    double energy = af_arap_energy(ref,faces,nf,uv);
    if (!isfinite(energy)) return ASM_FLAT_DEGENERATE;

    /* cotangent Laplacian */
    double *val = ARENA_CALLOC(arena, nnz, sizeof(double));
    double *diag = ARENA_CALLOC(arena, nv, sizeof(double));
    AfMat L = { nv, off, col, val, diag };
    for (size_t f = 0; f < nf; f++) {
        for (int e = 0; e < 3; e++) {
            int32_t a = faces[f*3 + (size_t)e];
            int32_t b = faces[f*3 + (size_t)((e+1)%3)];
            double c = ref[f].cot[e];
            val[af_slot(&L, a, b)] -= c;
            val[af_slot(&L, b, a)] -= c;
            diag[a] += c; diag[b] += c;
        }
    }
    /* pin one vertex (the outer loop start) at its Tutte position */
    int32_t pin = t->loop_start[(size_t)t->outer];
    double pin_u = uv[(size_t)pin*2], pin_v = uv[(size_t)pin*2+1];
    /* soft-free: keep the matrix SPD by a tiny mass term, plus the pin */
    for (size_t i = 0; i < nv; i++) diag[i] += 1e-9;
    {
        diag[pin] = 1.0;
        for (int32_t k = off[pin]; k < off[pin+1]; k++) val[k] = 0.0;
    }
    /* the pin column entries move into the RHS every iteration (they are
     * constant): remember them */
    double *pin_col = ARENA_CALLOC(arena, nv, sizeof(double));
    for (int32_t k = off[pin]; k < off[pin+1]; k++) {
        int32_t j = col[k];
        int32_t kj = af_slot(&L, j, pin);
        pin_col[j] = val[kj];
        val[kj] = 0.0;
    }

    double *bu = ARENA_ALLOC(arena, nv * sizeof(double));
    double *bv = ARENA_ALLOC(arena, nv * sizeof(double));
    double *xu = ARENA_ALLOC(arena, nv * sizeof(double));
    double *xv = ARENA_ALLOC(arena, nv * sizeof(double));
    double *r = ARENA_ALLOC(arena, nv * sizeof(double));
    double *z = ARENA_ALLOC(arena, nv * sizeof(double));
    double *p = ARENA_ALLOC(arena, nv * sizeof(double));
    double *q = ARENA_ALLOC(arena, nv * sizeof(double));
    double *Rt = ARENA_ALLOC(arena, nf * 4 * sizeof(double));
    for (size_t i = 0; i < nv; i++) { xu[i] = uv[i*2]; xv[i] = uv[i*2+1]; }

    double scale = 0.0;
    for (size_t f = 0; f < nf; f++) scale += ref[f].area;
    scale = sqrt(scale > 0.0 ? scale : 1.0);
    int iters = 0;
    AfSolver solver = {NULL, NULL};
    for (int it = 0; it < max_iters; it++) {
        if (af_expired()) break;
        /* local */
        for (size_t f = 0; f < nf; f++) {
            double J[4];
            af_jacobian(&ref[f], uv, &faces[f*3], J);
            af_polar(J, &Rt[f*4]);
        }
        /* global rhs */
        for (size_t i = 0; i < nv; i++) { bu[i] = 0.0; bv[i] = 0.0; }
        for (size_t f = 0; f < nf; f++) {
            double px[3] = { 0.0, ref[f].x1, ref[f].x2 };
            double py[3] = { 0.0, 0.0, ref[f].y2 };
            const double *R = &Rt[f*4];
            for (int e = 0; e < 3; e++) {
                int a = e, b = (e+1)%3;
                int32_t ia = faces[f*3 + (size_t)a], ib = faces[f*3 + (size_t)b];
                double dx = px[a] - px[b], dy = py[a] - py[b];
                double rx = R[0]*dx + R[1]*dy, ry = R[2]*dx + R[3]*dy;
                double c = ref[f].cot[e];
                bu[ia] += c * rx; bv[ia] += c * ry;
                bu[ib] -= c * rx; bv[ib] -= c * ry;
            }
        }
        /* pin */
        for (size_t i = 0; i < nv; i++) {
            if ((int32_t)i == pin) continue;
            bu[i] -= pin_col[i] * pin_u;
            bv[i] -= pin_col[i] * pin_v;
        }
        bu[pin] = pin_u; bv[pin] = pin_v;
        xu[pin] = pin_u; xv[pin] = pin_v;
        int maxit = (int)(2 * nv + 200);
        if (maxit > 20000) maxit = 20000;
        if (it == 0) af_solver_init(arena, &solver, &L);   /* the pinned Laplacian never changes */
        if (af_solver_solve(&solver, bu, xu, r, z, p, q, 1e-9, maxit) != 0 ||
            af_solver_solve(&solver, bv, xv, r, z, p, q, 1e-9, maxit) != 0) {
            af_solver_free(&solver);
            return ASM_FLAT_SOLVE;
        }
        for (size_t i = 0; i < nv; i++) { raw_target[2*i] = xu[i]; raw_target[2*i+1] = xv[i]; }
        if (!af_local_direction(uv,raw_target,target,encoded,local_alpha,blocked,nv,faces,nf,encoded_ready,t,boxes))
            memcpy(target,raw_target,2*nv*sizeof(double));
        double alpha = 1, next_energy = energy;
        int accepted = 0;
        for (int ls = 0; ls < 48; ls++, alpha *= .5) {
            int valid = 1;
            for (size_t f = 0; f < nf; f++)
                if (!UvGuard_interval(uv,target,faces+3*f,alpha)) { valid = 0; break; }
            if (!valid) continue;
            for (size_t k = 0; k < 2*nv; k++) trial[k] = uv[k]+alpha*(target[k]-uv[k]);
            if (!UvGuard_mesh(trial,faces,nf)) continue;
            if (!af_boundary_clear(t,trial,boxes,NULL)) continue;
            next_energy = af_arap_energy(ref,faces,nf,trial);
            if (isfinite(next_energy) && next_energy <= energy+32*DBL_EPSILON*fmax(1,fabs(energy))) {
                int encoded_ok = af_encoded_feasible(trial,encoded,nv,faces,nf) && af_boundary_clear(t,encoded,boxes,NULL);
                if (encoded_ready && !encoded_ok) continue;
                encoded_ready |= encoded_ok;
                accepted = 1; break;
            }
        }
        /* A failed direction cannot replace the last feasible iterate. */
        if (!accepted) break;
        double moved = 0.0;
        for (size_t i = 0; i < nv; i++) {
            double du = trial[i*2] - uv[i*2], dv = trial[i*2+1] - uv[i*2+1];
            double m = du*du + dv*dv;
            if (m > moved) moved = m;
            uv[i*2] = xu[i] = trial[i*2]; uv[i*2+1] = xv[i] = trial[i*2+1];
        }
        energy = next_energy;
        iters = it + 1;
        if (sqrt(moved) < 1e-4 * scale) break;
    }
    af_solver_free(&solver);
    *out_iters = iters;
    return ASM_FLAT_OK;
}

/* ---- frame normalization ---------------------------------------------------- */

static double af_sd_energy(const AfRef *ref, const int32_t *faces, size_t nf, const double *uv, double *bad_area)
{
    double energy = 0; if (bad_area) *bad_area = 0;
    for (size_t f = 0; f < nf; f++) {
        double j[4], lo, hi, det; af_jacobian(ref+f,uv,faces+3*f,j); af_singular(j,&lo,&hi,&det);
        if (!(det > 0) || !(lo > 0) || !isfinite(hi)) return INFINITY;
        double norm = j[0]*j[0]+j[1]*j[1]+j[2]*j[2]+j[3]*j[3];
        energy += ref[f].area*fmax(0,norm*(1+1/(det*det))-4);
        if (bad_area && (lo < .75 || hi > 1.25)) *bad_area += ref[f].area;
    } return energy;
}

/* Guarded ARAP can reach a limiting face before the rest metric is recovered.
 * Use the maintained native SLIM solver for that metric continuation, keeping
 * the same original triangles, boundary, gauge and actual float32 certificate.
 * Its candidate is a proposal; this routine independently checks the committed
 * path and stored metric. No metric substitution, cuts, guides or repacking. */
static void af_metric_continue(Arena_T arena, const float *xyz, size_t nv, const int32_t *faces, size_t nf,
                                double *uv, double area, int max_iters, AsmFlattenStats *st)
{
    if (max_iters <= 0 || af_expired()) return;
    Arena_Mark mark = Arena_save(arena); AfRef *ref = ARENA_ALLOC(arena,nf*sizeof(AfRef)); af_reference(xyz,faces,nf,ref);
    double bad = 0, baseline = af_sd_energy(ref,faces,nf,uv,&bad);
    /* A small damaged patch still fails the individual-face certificate.
     * Its source area must not decide whether it receives any repair. */
    if (!(bad > 0) || !isfinite(baseline)) { Arena_restore(arena,mark); return; }
    double *target = ARENA_ALLOC(arena,2*nv*sizeof(double)), *trial = ARENA_ALLOC(arena,2*nv*sizeof(double)), *encoded = ARENA_ALLOC(arena,2*nv*sizeof(double));
    SlimRefineOpts opts; SlimRefine_defaults(&opts); opts.iterations = max_iters < 60 ? 60 : max_iters;
    opts.pcg_iterations = 240; opts.threads = 1; opts.local_steps = 1; opts.local_step_min_faces = 1;
    opts.guard_boundary = 1; opts.strip_pack = 0; opts.verbose = 0;
    uint8_t *fixed = ARENA_CALLOC(arena,nv,1); fixed[0] = 1; opts.fixed_vertices = fixed;
    SlimRefineStats result; int rc = SlimRefine_run_double(arena,xyz,nv,faces,nf,uv,&opts,target,&result);
    if (!rc && !result.flips_after && !result.boundary_intersections_after) {
        for (int ls = 0; ls < 32; ls++) {
            double alpha = ldexp(1.,-ls); int valid = 1;
            /* The complete SLIM result already followed its guarded curved
             * sequence of steps. A straight shortcut between its endpoints
             * need not be feasible. Only interpolated fallback proposals
             * introduce a new path requiring this additional interval guard. */
            if (ls) for (size_t f = 0; f < nf; f++) if (!UvGuard_interval(uv,target,faces+3*f,alpha)) { valid = 0; break; }
            if (!valid) continue;
            for (size_t k = 0; k < 2*nv; k++) trial[k] = uv[k]+alpha*(target[k]-uv[k]);
            if (!af_encoded_feasible(trial,encoded,nv,faces,nf)) continue;
            double after_bad = 0, energy = af_sd_energy(ref,faces,nf,encoded,&after_bad);
            if (!isfinite(energy) || !(energy < baseline) || after_bad > bad+32*DBL_EPSILON*area) continue;
            Arena_Mark query_mark = Arena_save(arena);
            AsmContacts *query = AsmContacts_new(arena,nv,faces,nf,encoded,1e-7); AsmContactStats contacts;
            int clear = query && !AsmContacts_measure(query,encoded,NULL,NULL,NULL,&contacts) && !contacts.pairs;
            Arena_restore(arena,query_mark);
            if (!clear) continue;
            memcpy(uv,trial,2*nv*sizeof(double)); st->metric_refine_iters = (int)result.iterations_run;
            fprintf(stderr,"[assemble flatten] metric continuation: %zu faces, bad source area %.4g%% -> %.4g%%, alpha %.4g\n",nf,100*bad/area,100*after_bad/area,alpha); break;
        }
    }
    Arena_restore(arena,mark);
}

static void af_normalize_frame(double *uv, size_t nv, const int32_t *faces, size_t nf,
                               int mirror)
{
    (void)faces; (void)nf;
    if (mirror) for (size_t i = 0; i < nv; i++) uv[i*2] = -uv[i*2];
    double cu = 0.0, cv = 0.0;
    for (size_t i = 0; i < nv; i++) { cu += uv[i*2]; cv += uv[i*2+1]; }
    cu /= (double)nv; cv /= (double)nv;
    double sxx = 0.0, sxy = 0.0, syy = 0.0;
    for (size_t i = 0; i < nv; i++) {
        double x = uv[i*2] - cu, y = uv[i*2+1] - cv;
        sxx += x*x; sxy += x*y; syy += y*y;
    }
    double ang = 0.5 * atan2(2.0 * sxy, sxx - syy); /* principal axis angle */
    double c = cos(-ang), s = sin(-ang);
    double m3 = 0.0;
    for (size_t i = 0; i < nv; i++) {
        double x = uv[i*2] - cu, y = uv[i*2+1] - cv;
        double xr = c*x - s*y, yr = s*x + c*y;
        uv[i*2] = xr; uv[i*2+1] = yr;
        m3 += xr*xr*xr;
    }
    if (m3 < 0.0) for (size_t i = 0; i < nv; i++) { uv[i*2] = -uv[i*2]; uv[i*2+1] = -uv[i*2+1]; }
}

/* A stalled circular embedding need not be the only feasible initializer.
 * A projection of the ORIGINAL source geometry is another candidate. Keep
 * its original units, align only its rigid gauge, and certify every encoded
 * triangle before it can replace a map. */
static int af_recovery_candidate(Arena_T arena, const AfRef *ref,
                                  const int32_t *faces, size_t nf, size_t nv,
                                  const double *candidate, double *encoded,
                                  double *current, double *energy, double *bad)
{
    if (!af_encoded_feasible(candidate,encoded,nv,faces,nf)) return 0;
    double after_bad = 0, after = af_sd_energy(ref,faces,nf,encoded,&after_bad);
    if (!isfinite(after) || !(after < *energy) || after_bad > *bad) return 0;
    double total_area = 0, before10 = 0, after10 = 0;
    for (size_t f = 0; f < nf; f++) {
        double before_uv[6], before_j[4], after_j[4], blo,bhi,bdet,alo,ahi,adet;
        int32_t local[3] = {0,1,2};
        for (int k = 0; k < 3; k++) for (int d = 0; d < 2; d++)
            before_uv[2*k+d] = (float)current[2*(size_t)faces[3*f+k]+d];
        af_jacobian(ref+f,before_uv,local,before_j); af_singular(before_j,&blo,&bhi,&bdet);
        af_jacobian(ref+f,encoded,faces+3*f,after_j); af_singular(after_j,&alo,&ahi,&adet);
        total_area += ref[f].area;
        if (blo >= .9 && bhi <= 1.1) before10 += ref[f].area;
        if (alo >= .9 && ahi <= 1.1) after10 += ref[f].area;
        if (blo >= .75 && bhi <= 1.25) {
            if (alo < .75 || ahi > 1.25) return 0;
        } else {
            double old_violation = fmax(.75/blo,bhi/1.25), violation = fmax(.75/alo,ahi/1.25);
            if (violation > old_violation) return 0;
        }
    }
    /* Preserve the independent chart-area quality gate as well as every
     * strict original-face certificate. Lower SD energy alone is not enough. */
    if (after10 < fmin(.95*total_area,before10)) return 0;
    Arena_Mark mark = Arena_save(arena);
    AsmContacts *query = AsmContacts_new(arena,nv,faces,nf,encoded,1e-7); AsmContactStats contacts;
    int clear = query && !AsmContacts_measure(query,encoded,NULL,NULL,NULL,&contacts) && !contacts.pairs;
    Arena_restore(arena,mark);
    if (!clear) return 0;
    memcpy(current,candidate,2*nv*sizeof(double)); *energy = after; *bad = after_bad;
    return 1;
}

static void af_projection_recover(Arena_T arena, const float *xyz, size_t nv,
                                   const int32_t *faces, size_t nf, const AfTopo *topo,
                                   double *uv, int max_iters, AsmFlattenStats *st)
{
    if (max_iters <= 0) return;
    Arena_Mark mark = Arena_save(arena);
    AfRef *ref = ARENA_ALLOC(arena,nf*sizeof(AfRef)); af_reference(xyz,faces,nf,ref);
    double *encoded = ARENA_ALLOC(arena,2*nv*sizeof(double));
    for (size_t i = 0; i < 2*nv; i++) encoded[i] = (float)uv[i];
    double bad = 0, energy = af_sd_energy(ref,faces,nf,encoded,&bad), original_bad = bad;
    /* All faces may pass the broad band while too much original area still
     * fails the normal chart-quality gate. That map also needs recovery. */
    double total_area = 0, within10 = 0;
    for (size_t f = 0; f < nf; f++) {
        double j[4],lo,hi,det; af_jacobian(ref+f,encoded,faces+3*f,j); af_singular(j,&lo,&hi,&det);
        total_area += ref[f].area;
        if (lo >= .9 && hi <= 1.1) within10 += ref[f].area;
    }
    if (!(bad > 0 || within10 < .95*total_area) || !isfinite(energy)) { Arena_restore(arena,mark); return; }
    double center[3] = {0}, covariance[3][3] = {{0}}, eigen[3], basis[3][3];
    for (size_t i = 0; i < nv; i++) for (int d = 0; d < 3; d++) center[d] += xyz[3*i+d]/(double)nv;
    for (size_t i = 0; i < nv; i++) for (int a = 0; a < 3; a++) for (int b = 0; b < 3; b++)
        covariance[a][b] += ((double)xyz[3*i+a]-center[a])*((double)xyz[3*i+b]-center[b]);
    Eig3_sym(covariance,eigen,basis);
    double *candidate = ARENA_CALLOC(arena,2*nv,sizeof(double));
    for (size_t i = 0; i < nv; i++) for (int d = 0; d < 2; d++) for (int k = 0; k < 3; k++)
        candidate[2*i+d] += ((double)xyz[3*i+k]-center[k])*basis[k][d+1];
    if (!UvGuard_positive(candidate,faces)) for (size_t i = 0; i < nv; i++) candidate[2*i] = -candidate[2*i];
    if (!UvGuard_mesh(candidate,faces,nf)) { Arena_restore(arena,mark); return; }
    double cu = 0,cv = 0,tu = 0,tv = 0, dot = 0,cross = 0;
    for (size_t i = 0; i < nv; i++) {
        cu += candidate[2*i]/(double)nv; cv += candidate[2*i+1]/(double)nv;
        tu += uv[2*i]/(double)nv; tv += uv[2*i+1]/(double)nv;
    }
    for (size_t i = 0; i < nv; i++) {
        double u = candidate[2*i]-cu, v = candidate[2*i+1]-cv, a = uv[2*i]-tu, b = uv[2*i+1]-tv;
        dot += u*a+v*b; cross += u*b-v*a;
    }
    double theta = atan2(cross,dot), ct = cos(theta), sn = sin(theta);
    for (size_t i = 0; i < nv; i++) {
        double u = candidate[2*i]-cu, v = candidate[2*i+1]-cv;
        candidate[2*i] = ct*u-sn*v+tu; candidate[2*i+1] = sn*u+ct*v+tv;
    }
    if (!af_encoded_feasible(candidate,encoded,nv,faces,nf)) { Arena_restore(arena,mark); return; }
    Arena_Mark query_mark = Arena_save(arena);
    AsmContacts *query = AsmContacts_new(arena,nv,faces,nf,encoded,1e-7); AsmContactStats contacts;
    int clear = query && !AsmContacts_measure(query,encoded,NULL,NULL,NULL,&contacts) && !contacts.pairs;
    Arena_restore(arena,query_mark);
    if (!clear) { Arena_restore(arena,mark); return; }
    int accepted = af_recovery_candidate(arena,ref,faces,nf,nv,candidate,encoded,uv,&energy,&bad);
    int iters = 0;
    if (af_arap(arena,xyz,nv,faces,nf,topo,candidate,max_iters,&iters) == ASM_FLAT_OK) {
        AsmFlattenStats refinement = {0};
        af_metric_continue(arena,xyz,nv,faces,nf,candidate,st->area3d,max_iters,&refinement);
        if (af_recovery_candidate(arena,ref,faces,nf,nv,candidate,encoded,uv,&energy,&bad)) accepted = 1;
        st->arap_iters += iters; st->metric_refine_iters += refinement.metric_refine_iters;
    }
    if (accepted) fprintf(stderr,"[assemble flatten] source projection recovery: %zu faces, bad source area %.17g -> %.17g; original metric and encoded contact guards pass\n",nf,original_bad,bad);
    Arena_restore(arena,mark);
}

/* Find a local original-metric certificate when aggregate SD minimization
 * leaves a few faces outside the strict band. The caller keeps its old map
 * unless EVERY encoded original face passes, with no chart self-contact. */
typedef struct AfBand {
    const AfRef *ref; const int32_t *faces;
    const size_t *active, *vertices; const int32_t *index;
    size_t na, n;
    const double *seed, *scale;
    double *trial;
} AfBand;

static void af_band_basis(const AfRef *r, double g[6])
{
    g[2] = 1/r->x1; g[3] = -r->x2/(r->x1*r->y2);
    g[4] = 0; g[5] = 1/r->y2;
    g[0] = -g[2]; g[1] = -g[3]-g[5];
}

static double af_band_loss(AfBand *s, const double *x, double *gradient)
{
    double loss = 0;
    memset(gradient,0,2*s->n*sizeof *gradient);
    for (size_t i = 0; i < s->n; i++) for (size_t d = 0; d < 2; d++)
        s->trial[2*s->vertices[i]+d] = s->seed[2*s->vertices[i]+d]+x[2*i+d]/s->scale[i];
    for (size_t k = 0; k < s->na; k++) {
        size_t f = s->active[k]; const int32_t *fv = s->faces+3*f;
        double j[4],lo,hi,det; af_jacobian(s->ref+f,s->trial,fv,j); af_singular(j,&lo,&hi,&det);
        if (!(det > 0) || !(lo > 0) || !isfinite(hi)) return INFINITY;
        double dhi = hi-fmin(1.2498,fmax(.7502,hi));
        double dlo = lo-fmin(1.2498,fmax(.7502,lo));
        loss += dhi*dhi+dlo*dlo;
        if (!dhi && !dlo) continue;
        /* J V diag((sigma-clamp(sigma))/sigma) V^T is J-P(J).
         * The original source frame supplies every gradient column. */
        double a = j[0]*j[0]+j[2]*j[2], b = j[0]*j[1]+j[2]*j[3], c = j[1]*j[1]+j[3]*j[3];
        double theta = .5*atan2(2*b,a-c), ct = cos(theta), sn = sin(theta);
        double h = dhi/hi, l = dlo/lo;
        double q00 = h*ct*ct+l*sn*sn, q01 = (h-l)*ct*sn, q11 = h*sn*sn+l*ct*ct;
        double dj[4] = {2*(j[0]*q00+j[1]*q01),2*(j[0]*q01+j[1]*q11),
                        2*(j[2]*q00+j[3]*q01),2*(j[2]*q01+j[3]*q11)}, g[6];
        af_band_basis(s->ref+f,g);
        for (int v = 0; v < 3; v++) {
            int32_t i = s->index[fv[v]]; if (i < 0) continue;
            gradient[2*(size_t)i] += (dj[0]*g[2*v]+dj[1]*g[2*v+1])/s->scale[i];
            gradient[2*(size_t)i+1] += (dj[2]*g[2*v]+dj[3]*g[2*v+1])/s->scale[i];
        }
    }
    for (size_t i = 0; i < s->n; i++) for (size_t d = 0; d < 2; d++) {
        double move = x[2*i+d]/s->scale[i];
        loss += 1e-8*move*move;
        gradient[2*i+d] += 2e-8*move/s->scale[i];
    }
    return isfinite(loss) ? loss : INFINITY;
}

static double af_band_dot(const double *a, const double *b, size_t n)
{
    double result = 0; for (size_t i = 0; i < n; i++) result += a[i]*b[i]; return result;
}

static int af_band_certificate(Arena_T arena, const AfRef *ref,
                                const int32_t *faces, size_t nf, size_t nv,
                                const double *trial, double *encoded)
{
    if (!af_encoded_feasible(trial,encoded,nv,faces,nf)) return 0;
    for (size_t f = 0; f < nf; f++) {
        double j[4],lo,hi,det; af_jacobian(ref+f,encoded,faces+3*f,j); af_singular(j,&lo,&hi,&det);
        if (!(det > 0) || !(lo >= .75) || !(hi <= 1.25)) return 0;
    }
    Arena_Mark mark = Arena_save(arena);
    AsmContacts *query = AsmContacts_new(arena,nv,faces,nf,encoded,1e-7); AsmContactStats contacts;
    int clear = query && !AsmContacts_measure(query,encoded,NULL,NULL,NULL,&contacts) && !contacts.pairs;
    Arena_restore(arena,mark); return clear;
}

static void af_band_recover(Arena_T arena, const float *xyz, size_t nv,
                             const int32_t *faces, size_t nf, double *uv,
                             int max_iters, AsmFlattenStats *st)
{
    if (max_iters <= 0 || af_expired()) return;
    Arena_Mark mark = Arena_save(arena);
    AfRef *ref = ARENA_ALLOC(arena,nf*sizeof *ref); af_reference(xyz,faces,nf,ref);
    double *seed = ARENA_ALLOC(arena,2*nv*sizeof *seed), *encoded = ARENA_ALLOC(arena,2*nv*sizeof *encoded);
    if (!af_encoded_feasible(uv,seed,nv,faces,nf)) { Arena_restore(arena,mark); return; }
    uint8_t *free_vertex = ARENA_CALLOC(arena,nv,1), *next = ARENA_ALLOC(arena,nv);
    size_t bad = 0;
    for (size_t f = 0; f < nf; f++) {
        double j[4],lo,hi,det; af_jacobian(ref+f,seed,faces+3*f,j); af_singular(j,&lo,&hi,&det);
        if (lo < .75 || hi > 1.25) { bad++; for (int k = 0; k < 3; k++) free_vertex[faces[3*f+k]] = 1; }
    }
    if (!bad) { Arena_restore(arena,mark); return; }
    if (bad > ASM_FLATTEN_BAND_RECOVER_MAX_BAD) {
        /* measured hopeless (pipeline_constants.h): keep the map, the caller's
         * certificate decides, and the next seed or excision tier gets the time */
        fprintf(stderr,"[assemble flatten] local original-metric recovery: %zu faces, %zu strict failures -> skipped (over %d)\n",
                nf,bad,ASM_FLATTEN_BAND_RECOVER_MAX_BAD);
        Arena_restore(arena,mark); return;
    }
    for (int ring = 0; ring < 8; ring++) {
        memcpy(next,free_vertex,nv);
        for (size_t f = 0; f < nf; f++) if (free_vertex[faces[3*f]] || free_vertex[faces[3*f+1]] || free_vertex[faces[3*f+2]])
            for (int k = 0; k < 3; k++) next[faces[3*f+k]] = 1;
        uint8_t *swap = next; next = free_vertex; free_vertex = swap;
    }
    size_t n = 0, na = 0;
    for (size_t i = 0; i < nv; i++) n += free_vertex[i] != 0;
    for (size_t f = 0; f < nf; f++) na += free_vertex[faces[3*f]] || free_vertex[faces[3*f+1]] || free_vertex[faces[3*f+2]];
    /* This local proposal has bounded work. An oversized region retains its
     * original map and all failures remain in the normal complete audit. */
    if (!n || n > 8192 || na > 32768) { Arena_restore(arena,mark); return; }
    size_t *vertices = ARENA_ALLOC(arena,n*sizeof *vertices), *active = ARENA_ALLOC(arena,na*sizeof *active);
    int32_t *index = ARENA_ALLOC(arena,nv*sizeof *index);
    double *scale = ARENA_CALLOC(arena,n,sizeof *scale);
    size_t p = 0;
    for (size_t i = 0; i < nv; i++) { index[i] = -1; if (free_vertex[i]) { vertices[p] = i; index[i] = (int32_t)p++; } }
    p = 0;
    for (size_t f = 0; f < nf; f++) if (free_vertex[faces[3*f]] || free_vertex[faces[3*f+1]] || free_vertex[faces[3*f+2]]) {
        active[p++] = f; double g[6]; af_band_basis(ref+f,g);
        for (int k = 0; k < 3; k++) if (index[faces[3*f+k]] >= 0)
            scale[index[faces[3*f+k]]] += g[2*k]*g[2*k]+g[2*k+1]*g[2*k+1];
    }
    for (size_t i = 0; i < n; i++) scale[i] = sqrt(fmax(scale[i],1e-12));
    size_t count = 2*n; enum { HISTORY = 15 };
    double *x = ARENA_CALLOC(arena,count,sizeof *x), *g = ARENA_ALLOC(arena,count*sizeof *g);
    double *tx = ARENA_ALLOC(arena,count*sizeof *tx), *tg = ARENA_ALLOC(arena,count*sizeof *tg);
    double *direction = ARENA_ALLOC(arena,count*sizeof *direction);
    double *ss = ARENA_ALLOC(arena,HISTORY*count*sizeof *ss), *yy = ARENA_ALLOC(arena,HISTORY*count*sizeof *yy);
    double rho[HISTORY], coefficients[HISTORY]; int history = 0;
    double *trial = ARENA_ALLOC(arena,2*nv*sizeof *trial), *current = ARENA_ALLOC(arena,2*nv*sizeof *current);
    memcpy(trial,seed,2*nv*sizeof *trial); memcpy(current,seed,2*nv*sizeof *current);
    AfBand ctx = {ref,faces,active,vertices,index,na,n,seed,scale,trial};
    double value = af_band_loss(&ctx,x,g); int steps = 0, accepted = 0;
    int limit = max_iters > 100 ? 4000 : 40*max_iters;
    if (limit > ASM_FLATTEN_BAND_RECOVER_MAX_STEPS) limit = ASM_FLATTEN_BAND_RECOVER_MAX_STEPS;
    for (int step = 0; step < limit && isfinite(value); step++) {
        if (af_expired()) break;
        memcpy(direction,g,count*sizeof *direction);
        for (int h = history-1; h >= 0; h--) {
            coefficients[h] = rho[h]*af_band_dot(ss+(size_t)h*count,direction,count);
            for (size_t i = 0; i < count; i++) direction[i] -= coefficients[h]*yy[(size_t)h*count+i];
        }
        if (history) {
            double norm = af_band_dot(yy+(size_t)(history-1)*count,yy+(size_t)(history-1)*count,count);
            double factor = 1/(rho[history-1]*norm);
            for (size_t i = 0; i < count; i++) direction[i] *= factor;
        }
        for (int h = 0; h < history; h++) {
            double beta = rho[h]*af_band_dot(yy+(size_t)h*count,direction,count);
            for (size_t i = 0; i < count; i++) direction[i] += ss[(size_t)h*count+i]*(coefficients[h]-beta);
        }
        for (size_t i = 0; i < count; i++) direction[i] = -direction[i];
        double slope = af_band_dot(direction,g,count);
        if (!(slope < 0) || !isfinite(slope)) {
            history = 0; for (size_t i = 0; i < count; i++) direction[i] = -g[i];
            slope = -af_band_dot(g,g,count);
        }
        if (!(slope < 0) || !isfinite(slope)) break;
        double after = INFINITY; int found = 0;
        for (int line = 0; line < 40; line++) {
            double alpha = ldexp(1.,-line);
            for (size_t i = 0; i < count; i++) tx[i] = x[i]+alpha*direction[i];
            after = af_band_loss(&ctx,tx,tg);
            if (!(after < value) || after > value+1e-4*alpha*slope) continue;
            int path = 1;
            for (size_t k = 0; k < na; k++) if (!UvGuard_interval(current,trial,faces+3*active[k],1)) { path = 0; break; }
            if (path) { found = 1; break; }
        }
        if (!found) break;
        steps++;
        if (af_band_certificate(arena,ref,faces,nf,nv,trial,encoded)) {
            memcpy(uv,encoded,2*nv*sizeof *uv); accepted = 1; break;
        }
        double sy = 0, snorm = 0, ynorm = 0;
        for (size_t i = 0; i < count; i++) {
            double a = tx[i]-x[i], b = tg[i]-g[i]; sy += a*b; snorm += a*a; ynorm += b*b;
        }
        if (isfinite(sy) && sy > 1e-12*sqrt(snorm*ynorm)) {
            if (history == HISTORY) {
                memmove(ss,ss+count,(HISTORY-1)*count*sizeof *ss); memmove(yy,yy+count,(HISTORY-1)*count*sizeof *yy);
                memmove(rho,rho+1,(HISTORY-1)*sizeof *rho); history--;
            }
            for (size_t i = 0; i < count; i++) { ss[(size_t)history*count+i] = tx[i]-x[i]; yy[(size_t)history*count+i] = tg[i]-g[i]; }
            rho[history++] = 1/sy;
        }
        memcpy(x,tx,count*sizeof *x); memcpy(g,tg,count*sizeof *g);
        memcpy(current,trial,2*nv*sizeof *current); value = after;
    }
    st->metric_refine_iters += steps;
    fprintf(stderr,"[assemble flatten] local original-metric recovery: %zu faces, %zu free vertices / %zu active faces, %d steps; %zu strict failures -> %s\n",
            nf,n,na,steps,bad,accepted ? "0 (complete encoded certificate)" : "unchanged (proposal refused)");
    Arena_restore(arena,mark);
}

/* ---- public ---------------------------------------------------------------- */

void AsmFlatten_face_singular(const float *xyz, const int32_t *faces, size_t nf,
                              const float *uv, double *sig_lo, double *sig_hi,
                              int8_t *det_sign)
{
    for (size_t f = 0; f < nf; f++) {
        AfRef rf;
        af_reference(xyz, &faces[f*3], 1, &rf);
        double uvd[6];
        for (int k = 0; k < 3; k++) {
            uvd[k*2] = (double)uv[(size_t)faces[f*3+(size_t)k]*2];
            uvd[k*2+1] = (double)uv[(size_t)faces[f*3+(size_t)k]*2+1];
        }
        int32_t loc[3] = { 0, 1, 2 };
        double J[4], lo, hi, det;
        af_jacobian(&rf, uvd, loc, J);
        af_singular(J, &lo, &hi, &det);
        if (sig_lo) sig_lo[f] = lo;
        if (sig_hi) sig_hi[f] = hi;
        if (det_sign) det_sign[f] = (int8_t)(det > 0.0 ? 1 : (det < 0.0 ? -1 : 0));
    }
}

/* ---- free-boundary conformal seed ----------------------------------------- */

static int32_t af_row_slot(const int32_t *off, const int32_t *col, int32_t row, int32_t c)
{
    for (int32_t k = off[row]; k < off[row+1]; k++) if (col[k] == c) return k;
    return -1;
}

/* Least-squares conformal map (Levy, Petitjean, Ray & Maillot 2002) of the
 * ORIGINAL chart geometry, with no fixed boundary. The Cauchy-Riemann
 * residuals of each face in its own reference frame give a symmetric 2nv
 * system: the u-u and v-v blocks are the cotangent Dirichlet form and the u-v
 * block the (boundary-supported) area term. Two pinned vertices remove the
 * similarity gauge; the map is then scaled to the 3-D area. A developable
 * chart comes out as a similar copy of its isometric flattening whatever its
 * bending, which a circular Tutte seed cannot give a chart that turns back on
 * itself. The same Jacobi-preconditioned CG as the rest of this file solves
 * it (no TAUCS inside the per-cube OpenMP loop). Returns ASM_FLAT_*. */
static int af_lscm(Arena_T arena, const float *xyz, size_t nv, const int32_t *faces, size_t nf,
                   const AfTopo *t, const AfRef *ref, double area3d, double *uv)
{
    const int32_t *off = CSR_offset(t->adj);
    const int32_t *col = CSR_target(t->adj);
    size_t nnz = (size_t)CSR_nnz(t->adj), n = 2*nv;
    if (nv < 3 || nnz == 0) return ASM_FLAT_DEGENERATE;
    /* row i (u_i) and row nv+i (v_i) each hold the vertex's neighbours in the
     * same coordinate, then in the other coordinate; diagonals are separate */
    int32_t *moff = ARENA_ALLOC(arena, (n+1)*sizeof(int32_t));
    int32_t *mcol = ARENA_ALLOC(arena, 4*nnz*sizeof(int32_t));
    double *mval = ARENA_CALLOC(arena, 4*nnz, sizeof(double));
    double *mdiag = ARENA_CALLOC(arena, n, sizeof(double));
    for (size_t i = 0; i < nv; i++) {
        int32_t deg = off[i+1]-off[i];
        moff[i] = 2*off[i];
        moff[nv+i] = (int32_t)(2*nnz) + 2*off[i];
        for (int32_t k = 0; k < deg; k++) {
            int32_t j = col[off[i]+k];
            mcol[moff[i]+k] = j;             mcol[moff[i]+deg+k] = (int32_t)nv+j;
            mcol[moff[nv+i]+k] = (int32_t)nv+j; mcol[moff[nv+i]+deg+k] = j;
        }
    }
    moff[n] = (int32_t)(4*nnz);
    for (size_t f = 0; f < nf; f++) {
        double g[6], area = ref[f].area;
        af_band_basis(ref+f, g);
        for (int a = 0; a < 3; a++) for (int b = 0; b < 3; b++) {
            int32_t i = faces[3*f+(size_t)a], j = faces[3*f+(size_t)b];
            double same = area*(g[2*a]*g[2*b]+g[2*a+1]*g[2*b+1]);
            double cross = area*(g[2*a+1]*g[2*b]-g[2*a]*g[2*b+1]);   /* d2E / du_i dv_j */
            if (i == j) { mdiag[i] += same; mdiag[nv+(size_t)i] += same; continue; }
            int32_t s1 = af_row_slot(moff,mcol,i,j), s2 = af_row_slot(moff,mcol,(int32_t)nv+i,(int32_t)nv+j);
            int32_t s3 = af_row_slot(moff,mcol,i,(int32_t)nv+j), s4 = af_row_slot(moff,mcol,(int32_t)nv+j,i);
            if (s1 < 0 || s2 < 0 || s3 < 0 || s4 < 0) return ASM_FLAT_BAD_INDEX;
            mval[s1] += same; mval[s2] += same; mval[s3] += cross; mval[s4] += cross;
        }
    }
    /* gauge: the vertex farthest from the centroid, and the one farthest from it */
    double c[3] = {0,0,0};
    for (size_t i = 0; i < nv; i++) for (int d = 0; d < 3; d++) c[d] += xyz[3*i+(size_t)d]/(double)nv;
    size_t i0 = 0, i1 = 0; double best = -1;
    for (size_t i = 0; i < nv; i++) {
        double p[3]; af_get3(xyz,(int32_t)i,p); af_sub3(p,c,p);
        double d = af_dot3(p,p); if (d > best) { best = d; i0 = i; }
    }
    double p0[3]; af_get3(xyz,(int32_t)i0,p0); best = -1;
    for (size_t i = 0; i < nv; i++) {
        double p[3]; af_get3(xyz,(int32_t)i,p); af_sub3(p,p0,p);
        double d = af_dot3(p,p); if (d > best) { best = d; i1 = i; }
    }
    if (i0 == i1 || !(best > 0)) return ASM_FLAT_DEGENERATE;
    double *x = ARENA_CALLOC(arena, n, sizeof(double)), *b = ARENA_CALLOC(arena, n, sizeof(double));
    uint8_t *fixed = ARENA_CALLOC(arena, n, 1);
    size_t pins[4] = {i0, nv+i0, i1, nv+i1};
    x[i1] = sqrt(best);
    for (int k = 0; k < 4; k++) fixed[pins[k]] = 1;
    for (size_t i = 0; i < n; i++) {
        if (fixed[i]) {
            mdiag[i] = 1.0; b[i] = x[i];
            for (int32_t k = moff[i]; k < moff[i+1]; k++) mval[k] = 0.0;
            continue;
        }
        for (int32_t k = moff[i]; k < moff[i+1]; k++) if (fixed[mcol[k]]) {
            b[i] -= mval[k]*x[(size_t)mcol[k]]; mval[k] = 0.0;
        }
    }
    AfMat M = { n, moff, mcol, mval, mdiag };
    double *r = ARENA_ALLOC(arena, n*sizeof(double)), *z = ARENA_ALLOC(arena, n*sizeof(double));
    double *p = ARENA_ALLOC(arena, n*sizeof(double)), *q = ARENA_ALLOC(arena, n*sizeof(double));
    int maxit = (int)(4*n + 200);
    if (maxit > 20000) maxit = 20000;
    AfSolver solver; af_solver_init(arena, &solver, &M);
    int solved = af_solver_solve(&solver, b, x, r, z, p, q, 1e-10, maxit) == 0;
    af_solver_free(&solver);
    if (!solved) return ASM_FLAT_SOLVE;
    double area_uv = 0.0;
    for (size_t i = 0; i < nv; i++) { uv[2*i] = x[i]; uv[2*i+1] = x[nv+i]; }
    for (size_t f = 0; f < nf; f++) {
        double J[4], lo, hi, det; af_jacobian(ref+f,uv,faces+3*f,J); af_singular(J,&lo,&hi,&det);
        area_uv += det*ref[f].area;
    }
    if (area_uv < 0) { for (size_t i = 0; i < nv; i++) uv[2*i] = -uv[2*i]; area_uv = -area_uv; }
    if (!(area_uv > 0) || !isfinite(area_uv)) return ASM_FLAT_DEGENERATE;
    double s = sqrt(area3d/area_uv);
    for (size_t i = 0; i < 2*nv; i++) uv[i] *= s;
    return ASM_FLAT_OK;
}

/* Offer the conformal seed to a chart whose map still leaves the stress band.
 * The candidate keeps the chart's own original triangles and units, takes the
 * current rigid gauge, is refined by the same guarded ARAP and metric
 * continuation, and replaces the map only through af_recovery_candidate
 * (lower symmetric-Dirichlet energy, no new strict-band failure, the chart
 * area gate, a clean float32 encoding and no self-contact). A chart inside
 * the band is never touched. */
static void af_conformal_recover(Arena_T arena, const float *xyz, size_t nv,
                                  const int32_t *faces, size_t nf, const AfTopo *topo,
                                  double *uv, double band, int max_iters, AsmFlattenStats *st)
{
    if (max_iters <= 0) return;
    Arena_Mark mark = Arena_save(arena);
    AfRef *ref = ARENA_ALLOC(arena,nf*sizeof(AfRef)); af_reference(xyz,faces,nf,ref);
    double *encoded = ARENA_ALLOC(arena,2*nv*sizeof(double));
    for (size_t i = 0; i < 2*nv; i++) encoded[i] = (float)uv[i];
    size_t stressed = 0;
    for (size_t f = 0; f < nf; f++) {
        double J[4],lo,hi,det; af_jacobian(ref+f,encoded,faces+3*f,J); af_singular(J,&lo,&hi,&det);
        if (lo < 1.0-band || hi > 1.0+band || det <= 0.0) stressed++;
    }
    double bad = 0, energy = af_sd_energy(ref,faces,nf,encoded,&bad), original_bad = bad;
    if (!stressed || !isfinite(energy)) { Arena_restore(arena,mark); return; }
    double *candidate = ARENA_ALLOC(arena,2*nv*sizeof(double));
    if (af_lscm(arena,xyz,nv,faces,nf,topo,ref,st->area3d,candidate) != ASM_FLAT_OK ||
        !UvGuard_mesh(candidate,faces,nf)) { Arena_restore(arena,mark); return; }
    double cu = 0,cv = 0,tu = 0,tv = 0, dot = 0,cross = 0;
    for (size_t i = 0; i < nv; i++) {
        cu += candidate[2*i]/(double)nv; cv += candidate[2*i+1]/(double)nv;
        tu += uv[2*i]/(double)nv; tv += uv[2*i+1]/(double)nv;
    }
    for (size_t i = 0; i < nv; i++) {
        double u = candidate[2*i]-cu, v = candidate[2*i+1]-cv, a = uv[2*i]-tu, w = uv[2*i+1]-tv;
        dot += u*a+v*w; cross += u*w-v*a;
    }
    double theta = atan2(cross,dot), ct = cos(theta), sn = sin(theta);
    for (size_t i = 0; i < nv; i++) {
        double u = candidate[2*i]-cu, v = candidate[2*i+1]-cv;
        candidate[2*i] = ct*u-sn*v+tu; candidate[2*i+1] = sn*u+ct*v+tv;
    }
    if (!af_encoded_feasible(candidate,encoded,nv,faces,nf)) { Arena_restore(arena,mark); return; }
    Arena_Mark query_mark = Arena_save(arena);
    AsmContacts *query = AsmContacts_new(arena,nv,faces,nf,encoded,1e-7); AsmContactStats contacts;
    int clear = query && !AsmContacts_measure(query,encoded,NULL,NULL,NULL,&contacts) && !contacts.pairs;
    Arena_restore(arena,query_mark);
    if (!clear) { Arena_restore(arena,mark); return; }
    int accepted = af_recovery_candidate(arena,ref,faces,nf,nv,candidate,encoded,uv,&energy,&bad);
    int iters = 0;
    if (af_arap(arena,xyz,nv,faces,nf,topo,candidate,max_iters,&iters) == ASM_FLAT_OK) {
        AsmFlattenStats refinement = {0};
        af_metric_continue(arena,xyz,nv,faces,nf,candidate,st->area3d,max_iters,&refinement);
        if (af_recovery_candidate(arena,ref,faces,nf,nv,candidate,encoded,uv,&energy,&bad)) accepted = 1;
        st->arap_iters += iters; st->metric_refine_iters += refinement.metric_refine_iters;
    }
    if (accepted) {
        size_t after = 0;
        for (size_t i = 0; i < 2*nv; i++) encoded[i] = (float)uv[i];
        for (size_t f = 0; f < nf; f++) {
            double J[4],lo,hi,det; af_jacobian(ref+f,encoded,faces+3*f,J); af_singular(J,&lo,&hi,&det);
            if (lo < 1.0-band || hi > 1.0+band || det <= 0.0) after++;
        }
        fprintf(stderr,"[assemble flatten] conformal recovery: %zu faces, stressed faces %zu -> %zu, bad source area %.6g -> %.6g; original metric and encoded contact guards pass\n",
                nf,stressed,after,original_bad,bad);
    }
    Arena_restore(arena,mark);
}

/* Conformal map of the chart with every hole filled, for the embedding only, by a fan from a virtual
 * vertex at the loop's 3-D centroid (wound opposite to the real faces along the loop).  The free
 * outer boundary avoids the circle's distortion of large or elongated charts, and the fans keep
 * holes open (the plain conformal map lets slit-like and tiny holes pinch shut: PHerc0343 seeds had
 * touching hole boundaries on 9 of 10 multi-loop charts).  Writes the real vertices' uv. */
static int af_lscm_capped(Arena_T arena, const float *xyz, size_t nv, const int32_t *faces, size_t nf,
                          const AfTopo *t, double area3d, double *uv)
{
    size_t holes = t->n_loops > 0 ? t->n_loops - 1 : 0, extra_f = 0;
    for (size_t l = 0; l < t->n_loops; l++) if ((int32_t)l != t->outer) extra_f += t->loop_len[l];
    size_t anv = nv + holes, anf = nf + extra_f;
    float *ax = ARENA_ALLOC(arena, anv * 3 * sizeof(float));
    int32_t *af = ARENA_ALLOC(arena, anf * 3 * sizeof(int32_t));
    memcpy(ax, xyz, nv * 3 * sizeof(float));
    memcpy(af, faces, nf * 3 * sizeof(int32_t));
    /* directed boundary edge a->b of a real face: the fan face is (b, a, hub) */
    size_t f = nf, h = nv;
    for (size_t l = 0; l < t->n_loops; l++) {
        if ((int32_t)l == t->outer) continue;
        double c[3] = {0, 0, 0};
        int32_t v = t->loop_start[l];
        for (size_t k = 0; k < t->loop_len[l]; k++) {
            for (int d = 0; d < 3; d++) c[d] += (double)xyz[3*(size_t)v + (size_t)d] / (double)t->loop_len[l];
            v = t->bnext[v];
        }
        for (int d = 0; d < 3; d++) ax[3*h + (size_t)d] = (float)c[d];
        v = t->loop_start[l];
        for (size_t k = 0; k < t->loop_len[l]; k++) {
            int32_t nx = t->bnext[v];
            /* the real face holding boundary edge {v,nx} fixes the direction */
            int dir_vn = 0;
            for (size_t g = 0; g < nf && !dir_vn; g++)
                for (int e = 0; e < 3; e++)
                    if (faces[3*g + (size_t)e] == v && faces[3*g + (size_t)((e + 1) % 3)] == nx) { dir_vn = 1; break; }
            af[3*f] = dir_vn ? nx : v; af[3*f + 1] = dir_vn ? v : nx; af[3*f + 2] = (int32_t)h;
            f++;
            v = nx;
        }
        h++;
    }
    AfTopo at;
    if (af_topology(arena, ax, anv, af, anf, &at) != ASM_FLAT_OK) return ASM_FLAT_NONMANIFOLD;
    AfRef *ref = ARENA_ALLOC(arena, anf * sizeof(AfRef));
    af_reference(ax, af, anf, ref);
    for (size_t g = 0; g < anf; g++) if (!(ref[g].area > 0) || !isfinite(ref[g].area)) return ASM_FLAT_DEGENERATE;
    double *auv = ARENA_ALLOC(arena, anv * 2 * sizeof(double));
    int rc = af_lscm(arena, ax, anv, af, anf, &at, ref, area3d, auv);
    if (rc != ASM_FLAT_OK) return rc;
    memcpy(uv, auv, nv * 2 * sizeof(double));
    return ASM_FLAT_OK;
}

/* ---- pinch points ------------------------------------------------------------- */

typedef struct AfPos { float x, y, z; int32_t v; } AfPos;

static int af_pos_cmp(const void *pa, const void *pb)
{
    const AfPos *a = pa, *b = pb;
    if (a->x != b->x) return a->x < b->x ? -1 : 1;
    if (a->y != b->y) return a->y < b->y ? -1 : 1;
    if (a->z != b->z) return a->z < b->z ? -1 : 1;
    return a->v < b->v ? -1 : (a->v > b->v);
}

/* Boundary vertices that share a 3-D position with another boundary vertex: the mesher splits a
 * bowtie (pinch) vertex into copies at one position, and a slit closed by two zero-area faces
 * leaves the same kind of copies once those faces are dropped.  group[v] = 1 + the group of copies
 * at v's position, 0 elsewhere.  Returns the number of copies. */
static size_t af_pinch_copies(Arena_T arena, const float *xyz, size_t nv, const AfTopo *t, int32_t *group)
{
    for (size_t i = 0; i < nv; i++) group[i] = 0;
    size_t nb = 0;
    for (size_t i = 0; i < nv; i++) nb += t->on_boundary[i] != 0;
    if (nb < 2) return 0;
    Arena_Mark mark = Arena_save(arena);
    AfPos *p = ARENA_ALLOC(arena, nb * sizeof(AfPos));
    nb = 0;
    for (size_t i = 0; i < nv; i++) if (t->on_boundary[i]) {
        p[nb].x = xyz[3*i]; p[nb].y = xyz[3*i+1]; p[nb].z = xyz[3*i+2]; p[nb].v = (int32_t)i; nb++;
    }
    qsort(p, nb, sizeof(AfPos), af_pos_cmp);
    size_t n = 0;
    int32_t groups = 0;
    for (size_t i = 0; i < nb; ) {
        size_t j = i + 1;
        while (j < nb && p[j].x == p[i].x && p[j].y == p[i].y && p[j].z == p[i].z) j++;
        if (j - i > 1) {
            groups++;
            for (size_t k = i; k < j; k++) group[p[k].v] = groups;
            n += j - i;
        }
        i = j;
    }
    Arena_restore(arena, mark);
    return n;
}

/* uv corners of face f, as the contact predicates take them */
static void af_face_uv(const double *uv, const int32_t *faces, int32_t f, double out[6])
{
    for (int c = 0; c < 3; c++) {
        out[2*c] = uv[2*(size_t)faces[3*(size_t)f + (size_t)c]];
        out[2*c+1] = uv[2*(size_t)faces[3*(size_t)f + (size_t)c]+1];
    }
}

/* A near-isometric seed keeps pinch copies together, so their fans touch or overlap and the chart
 * can never pass the contact certificate (PHerc0343: 10 of 13 unflattenable charts; the circular
 * Tutte seed separates copies only by crushing the region between them).  Each group of copies is
 * first gathered at its mean uv (a conformal solve leaves copies a few thousandths of a voxel apart
 * in arbitrary directions, which can put one copy inside another's fan), then every copy moves into
 * its own fan along the mean direction to its incident face centroids.  The step starts at 1e-3 of
 * the group's shortest incident uv edge (inside the stress band, outside the 1e-7 contact tolerance
 * and the float32 lattice) and grows while faces of different copies still overlap, up to a quarter
 * edge.  No move may take an incident face out of strict orientation.  Returns the number of copies
 * moved. */
static size_t af_open_pinches(Arena_T arena, const float *xyz, size_t nv, const int32_t *faces, size_t nf,
                              const AfTopo *t, double *uv)
{
    Arena_Mark mark = Arena_save(arena);
    int32_t *group = ARENA_ALLOC(arena, nv * sizeof(int32_t));
    if (!af_pinch_copies(arena, xyz, nv, t, group)) { Arena_restore(arena, mark); return 0; }
    int32_t *foff = ARENA_CALLOC(arena, nv + 1, sizeof(int32_t));
    for (size_t k = 0; k < nf * 3; k++) foff[faces[k] + 1]++;
    for (size_t i = 0; i < nv; i++) foff[i + 1] += foff[i];
    int32_t *fill = ARENA_ALLOC(arena, nv * sizeof(int32_t));
    memcpy(fill, foff, nv * sizeof(int32_t));
    int32_t *fidx = ARENA_ALLOC(arena, nf * 3 * sizeof(int32_t));
    for (size_t f = 0; f < nf; f++) for (int c = 0; c < 3; c++) fidx[fill[faces[3*f + (size_t)c]]++] = (int32_t)f;
    int32_t ngroups = 0;
    for (size_t i = 0; i < nv; i++) if (group[i] > ngroups) ngroups = group[i];
    /* the copies of each group, contiguous */
    int32_t *goff = ARENA_CALLOC(arena, (size_t)ngroups + 2, sizeof(int32_t));
    for (size_t i = 0; i < nv; i++) if (group[i]) goff[group[i] + 1]++;
    for (int32_t g = 0; g <= ngroups; g++) goff[g + 1] += goff[g];
    int32_t *gfill = ARENA_ALLOC(arena, ((size_t)ngroups + 1) * sizeof(int32_t));
    memcpy(gfill, goff, ((size_t)ngroups + 1) * sizeof(int32_t));
    int32_t *gv = ARENA_ALLOC(arena, ((size_t)goff[ngroups + 1] + 1) * sizeof(int32_t));
    for (size_t i = 0; i < nv; i++) if (group[i]) gv[gfill[group[i]]++] = (int32_t)i;
    double *home = ARENA_ALLOC(arena, 2 * nv * sizeof(double));
    memcpy(home, uv, 2 * nv * sizeof(double));
    size_t moved = 0;
    for (int32_t g = 1; g <= ngroups; g++) {
        int32_t b = goff[g], e = goff[g + 1];
        if (e - b < 2) continue;
        double mu = 0.0, mv = 0.0, shortest = INFINITY;
        for (int32_t k = b; k < e; k++) { mu += home[2*(size_t)gv[k]] / (double)(e - b); mv += home[2*(size_t)gv[k]+1] / (double)(e - b); }
        /* each copy's opening direction, and the group's shortest incident edge */
        double dir[2*8];
        int ncopy = 0;
        for (int32_t k = b; k < e && ncopy < 8; k++, ncopy++) {
            int32_t v = gv[k];
            double du = 0.0, dv = 0.0;
            for (int32_t q = foff[v]; q < foff[v + 1]; q++) {
                const int32_t *fv = faces + 3*(size_t)fidx[q];
                double cu = 0.0, cv = 0.0;
                for (int c = 0; c < 3; c++) {
                    cu += home[2*(size_t)fv[c]] / 3.0; cv += home[2*(size_t)fv[c]+1] / 3.0;
                    if (fv[c] == v) continue;
                    double eu = home[2*(size_t)fv[c]] - mu, ev = home[2*(size_t)fv[c]+1] - mv;
                    double el = sqrt(eu*eu + ev*ev);
                    if (el < shortest) shortest = el;
                }
                cu -= mu; cv -= mv;
                double cl = sqrt(cu*cu + cv*cv);
                if (cl > 0.0) { du += cu / cl; dv += cv / cl; }
            }
            double dl = sqrt(du*du + dv*dv);
            dir[2*ncopy] = dl > 0.0 ? du / dl : 0.0;
            dir[2*ncopy+1] = dl > 0.0 ? dv / dl : 0.0;
        }
        if (!(shortest > 0.0) || !isfinite(shortest)) continue;
        for (double frac = 1e-3; frac <= 0.25 + 1e-12; frac *= 4.0) {
            /* every copy at mean + step * direction, halving a step that would fold a face */
            for (int c = 0; c < ncopy; c++) {
                int32_t v = gv[b + c];
                int placed = 0;
                for (double step = frac * shortest; step > 1e-9 * shortest; step *= 0.5) {
                    uv[2*(size_t)v] = mu + step * dir[2*c]; uv[2*(size_t)v+1] = mv + step * dir[2*c+1];
                    int ok = 1;
                    for (int32_t q = foff[v]; q < foff[v + 1] && ok; q++) ok = UvGuard_positive(uv, faces + 3*(size_t)fidx[q]);
                    if (ok) { placed = 1; break; }
                }
                if (!placed) { uv[2*(size_t)v] = home[2*(size_t)v]; uv[2*(size_t)v+1] = home[2*(size_t)v+1]; }
            }
            /* faces of different copies must not overlap */
            int clear = 1;
            for (int c = 0; c < ncopy && clear; c++)
                for (int d = c + 1; d < ncopy && clear; d++) {
                    int32_t v = gv[b + c], w = gv[b + d];
                    for (int32_t q = foff[v]; q < foff[v + 1] && clear; q++)
                        for (int32_t r = foff[w]; r < foff[w + 1] && clear; r++) {
                            if (fidx[q] == fidx[r]) continue;
                            double fa[6], fb[6];
                            af_face_uv(uv, faces, fidx[q], fa); af_face_uv(uv, faces, fidx[r], fb);
                            if (AsmContacts_pair(fa, fb, 1e-7, NULL, NULL) > 0) clear = 0;
                        }
                }
            if (clear) break;
        }
        for (int c = 0; c < ncopy; c++) {
            int32_t v = gv[b + c];
            if (uv[2*(size_t)v] != home[2*(size_t)v] || uv[2*(size_t)v+1] != home[2*(size_t)v+1]) moved++;
        }
    }
    Arena_restore(arena, mark);
    return moved;
}

/* Seed modes: 0 Tutte with uniform hole caps (the historical seed), 1 capped conformal, 2 Tutte with
 * mean-value hole caps, 3 plain conformal.  The recovery modes (1-3) open pinch points.  Returns 1
 * when the seed is strictly oriented. */
static int af_seed(Arena_T arena, const float *xyz, size_t nv, const int32_t *faces, size_t nf,
                   const AfTopo *topo, double area3d, int mode, double *uv, size_t *capped, int *out_rc)
{
    int rc = ASM_FLAT_OK;
    if (mode == 0 || mode == 2) {
        rc = af_tutte(arena, xyz, nv, faces, nf, topo, area3d, uv, capped, mode == 2);
        if (rc == ASM_FLAT_OK && !UvGuard_positive(uv,faces)) for (size_t i = 0; i < nv; i++) uv[2*i] = -uv[2*i];
    } else if (mode == 1) {
        rc = topo->n_loops > 1 ? af_lscm_capped(arena, xyz, nv, faces, nf, topo, area3d, uv) : ASM_FLAT_OK;
        if (topo->n_loops <= 1) mode = 3;
    }
    if (mode == 3) {
        AfRef *ref = ARENA_ALLOC(arena, nf * sizeof(AfRef));
        af_reference(xyz, faces, nf, ref);
        rc = af_lscm(arena, xyz, nv, faces, nf, topo, ref, area3d, uv);
    }
    *out_rc = rc;
    if (rc != ASM_FLAT_OK || !UvGuard_mesh(uv,faces,nf)) return 0;
    if (mode != 0) af_open_pinches(arena, xyz, nv, faces, nf, topo, uv);
    return 1;
}

static int af_flatten_chart(Arena_T arena, const float *xyz, size_t nv, const int32_t *faces, size_t nf,
                            double band, int max_iters, int conformal, int seed_mode, float *out_uv,
                            AsmFlattenStats *st);

/* The audit's chart-metric certificate on the stored float32 map (asm_audit: strict_bad_faces and
 * bad_metric_charts): every face strictly oriented with singular values in [0.75, 1.25], and at least
 * 95% of the original area within [0.9, 1.1].  A map that fails it cannot be part of a qualified
 * assembly. */
static int af_chart_certified(const float *xyz, const int32_t *faces, size_t nf, const float *uv)
{
    double total = 0.0, within10 = 0.0;
    for (size_t f = 0; f < nf; f++) {
        AfRef r;
        af_reference(xyz, faces + 3*f, 1, &r);
        double corner[6];
        for (int c = 0; c < 3; c++) {
            corner[2*c] = uv[2*(size_t)faces[3*f + (size_t)c]];
            corner[2*c+1] = uv[2*(size_t)faces[3*f + (size_t)c]+1];
        }
        const int32_t local[3] = {0, 1, 2};
        double J[4], lo, hi, det;
        af_jacobian(&r, corner, local, J);
        af_singular(J, &lo, &hi, &det);
        if (!UvGuard_positive(corner, local) || !(det > 0.0) || !AsmMetric_within25(lo, hi)) return 0;
        total += r.area;
        if (AsmMetric_within10(lo, hi)) within10 += r.area;
    }
    return within10 + 32.0 * DBL_EPSILON * total >= 0.95 * total;
}

/* Run the given seed modes, in order, through the full refinement and return the first map that
 * carries the audit's chart-metric certificate.  A valid map outside it (a guarded ARAP stalled on
 * a crushed circular seed leaves faces at 1e-3 of their length) would enter the layout as a bad
 * chart, so it is not a recovery.  On total failure the first mode's statistics are reported. */
/* The seed a requested mode actually builds (af_seed): the capped conformal seed (1) of a chart with a
 * single boundary loop has no hole to cap and IS the conformal seed (3). */
static int af_seed_effective_mode(int mode, size_t n_loops)
{
    return mode == 1 && n_loops <= 1 ? 3 : mode;
}

static int af_flatten_best(Arena_T arena, const float *xyz, size_t nv, const int32_t *faces, size_t nf,
                           double band, int max_iters, const int *modes, int n_modes, float *out_uv,
                           AsmFlattenStats *st)
{
    Arena_Mark mark = Arena_save(arena);
    float *trial = ARENA_ALLOC(arena, nv * 2 * sizeof(float));
    AsmFlattenStats first, ts;
    memset(&first, 0, sizeof first);
    int have = 0, tried[4] = {0, 0, 0, 0};
    /* a mode that builds a seed already refined here is skipped: the refinement is deterministic, so
     * it would only repeat a failure (a disk's capped conformal seed doubled the rescue's work) */
    size_t n_loops = 2;
    {
        Arena_Mark tm = Arena_save(arena);
        AfTopo topo;
        if (af_topology(arena, xyz, nv, faces, nf, &topo) == ASM_FLAT_OK) n_loops = topo.n_loops;
        Arena_restore(arena, tm);
    }
    for (int k = 0; k < n_modes && !have; k++) {
        int eff = af_seed_effective_mode(modes[k], n_loops);
        if (eff >= 0 && eff < 4) {
            if (tried[eff]) continue;
            tried[eff] = 1;
        }
        int rc = af_flatten_chart(arena,xyz,nv,faces,nf,band,max_iters,ASM_FLATTEN_CONFORMAL_RECOVERY,modes[k],trial,&ts);
        if (k == 0) first = ts;
        if (rc != 0 || !af_chart_certified(xyz, faces, nf, trial)) continue;
        memcpy(out_uv, trial, nv * 2 * sizeof(float));
        *st = ts;
        have = 1;
    }
    if (!have) *st = first;
    Arena_restore(arena, mark);
    return have ? 0 : -1;
}

int AsmFlatten_chart(Arena_T arena,
                     const float *xyz, size_t nv,
                     const int32_t *faces, size_t nf,
                     double band, int max_iters,
                     float *out_uv, AsmFlattenStats *st)
{
    /* the historical seed first: every chart it flattens stays bitwise unchanged */
    if (af_flatten_chart(arena,xyz,nv,faces,nf,band,max_iters,ASM_FLATTEN_CONFORMAL_RECOVERY,0,out_uv,st) == 0)
        return 0;
    AsmFlattenStats failed = *st;
    static const int recovery[] = {3, 2, 1};
    if (af_flatten_best(arena,xyz,nv,faces,nf,band,max_iters,recovery,3,out_uv,st) == 0) return 0;
    *st = failed;
    return -1;
}

int AsmFlatten_chart_certified(Arena_T arena, const float *xyz, size_t nv, const int32_t *faces, size_t nf,
                               double band, int max_iters, float *out_uv, AsmFlattenStats *st)
{
    if (af_flatten_chart(arena,xyz,nv,faces,nf,band,max_iters,ASM_FLATTEN_CONFORMAL_RECOVERY,0,out_uv,st) == 0 &&
        af_chart_certified(xyz, faces, nf, out_uv))
        return 0;
    AsmFlattenStats failed = *st;
    static const int recovery[] = {3, 2, 1};
    if (af_flatten_best(arena,xyz,nv,faces,nf,band,max_iters,recovery,3,out_uv,st) == 0) return 0;
    *st = failed;
    if (st->ok) { st->ok = 0; st->fail_reason = ASM_FLAT_FOLDED; st->fail_stage = 6; }
    return -1;
}

static int af_flatten_chart(Arena_T arena, const float *xyz, size_t nv, const int32_t *faces, size_t nf,
                            double band, int max_iters, int conformal, int seed_mode, float *out_uv,
                            AsmFlattenStats *st)
{
    memset(st, 0, sizeof *st);
    st->nv = nv; st->nf = nf;
    if (nv < 3 || nf < 1) { st->fail_reason = ASM_FLAT_EMPTY; return -1; }
    for (size_t k = 0; k < nf*3; k++)
        if (faces[k] < 0 || (size_t)faces[k] >= nv) { st->fail_reason = ASM_FLAT_BAD_INDEX; return -1; }

    Arena_Mark mark = Arena_save(arena);
    int rc = ASM_FLAT_OK;
    AfTopo topo;
    rc = af_topology(arena, xyz, nv, faces, nf, &topo);
    if (rc != ASM_FLAT_OK) { st->fail_reason = rc; Arena_restore(arena, mark); return -1; }
    st->n_boundary_loops = topo.n_loops;

    double area3d = 0.0;
    {
        AfRef *ref = ARENA_ALLOC(arena, nf * sizeof(AfRef));
        af_reference(xyz, faces, nf, ref);
        for (size_t f = 0; f < nf; f++) {
            if (!(ref[f].area > 0) || !isfinite(ref[f].area)) {
                st->fail_reason = ASM_FLAT_DEGENERATE; Arena_restore(arena,mark); return -1;
            }
            area3d += ref[f].area;
        }
    }
    st->area3d = area3d;
    if (!(area3d > 0.0)) { st->fail_reason = ASM_FLAT_DEGENERATE; Arena_restore(arena, mark); return -1; }

    double *uv = ARENA_ALLOC(arena, nv * 2 * sizeof(double));
    size_t capped = 0;
    int seeded = af_seed(arena, xyz, nv, faces, nf, &topo, area3d, seed_mode, uv, &capped, &rc);
    st->n_holes_capped = capped;
    st->seed_mode = seed_mode;
    if (!seeded) {
        st->fail_reason = rc != ASM_FLAT_OK ? rc : ASM_FLAT_FOLDED; st->fail_stage = 1;
        Arena_restore(arena,mark);
        return -1;
    }
    /* Freeze the deterministic seed gauge before optimization. Recomputing a
     * PCA rotation after each trial changes the float32 rounding of untouched
     * limiting faces and defeats incident-vertex damping across the chart. */
    af_normalize_frame(uv,nv,faces,nf,0);

    int iters = 0;
    rc = af_arap(arena, xyz, nv, faces, nf, &topo, uv, max_iters, &iters);
    if (rc != ASM_FLAT_OK) { st->fail_reason = rc; st->fail_stage = 2; Arena_restore(arena, mark); return -1; }
    st->arap_iters = iters;
    af_metric_continue(arena,xyz,nv,faces,nf,uv,area3d,max_iters,st);
    af_projection_recover(arena,xyz,nv,faces,nf,&topo,uv,max_iters,st);
    if (conformal) af_conformal_recover(arena,xyz,nv,faces,nf,&topo,uv,band,max_iters,st);
    af_band_recover(arena,xyz,nv,faces,nf,uv,max_iters,st);

    /* orientation: the majority of faces decides; mirror if negative */
    {
        AfRef *ref = ARENA_ALLOC(arena, nf * sizeof(AfRef));
        af_reference(xyz, faces, nf, ref);
        double pos_area = 0.0, neg_area = 0.0;
        for (size_t f = 0; f < nf; f++) {
            double J[4], lo, hi, det;
            af_jacobian(&ref[f], uv, &faces[f*3], J);
            af_singular(J, &lo, &hi, &det);
            if (det >= 0.0) pos_area += ref[f].area; else neg_area += ref[f].area;
        }
        if (neg_area > pos_area) { st->fail_reason = ASM_FLAT_FOLDED; st->fail_stage = 3; Arena_restore(arena,mark); return -1; }
        /* The cache is float32. Audit exactly those values after the final
         * frame change rather than declaring the unencoded double map safe. */
        for (size_t i = 0; i < nv*2; i++) uv[i] = out_uv[i] = (float)uv[i];
        /* final statistics */
        st->sigma_lo = 1e300; st->sigma_hi = 0.0;
        size_t n_stress = 0;
        double area_uv = 0.0;
        for (size_t f = 0; f < nf; f++) {
            double J[4], lo, hi, det;
            af_jacobian(&ref[f], uv, &faces[f*3], J);
            af_singular(J, &lo, &hi, &det);
            if (det < 0.0) st->n_flipped++;
            if (lo < st->sigma_lo) st->sigma_lo = lo;
            if (hi > st->sigma_hi) st->sigma_hi = hi;
            if (lo < 1.0 - band || hi > 1.0 + band || det <= 0.0) n_stress++;
            area_uv += fabs(det) * ref[f].area;
        }
        st->area_uv = area_uv;
        st->stress_frac = nf ? (double)n_stress / (double)nf : 0.0;
    }
    if (!UvGuard_mesh(uv,faces,nf)) {
        st->fail_reason = ASM_FLAT_FOLDED; st->fail_stage = 4; Arena_restore(arena,mark); return -1;
    }
    AsmContacts *query = AsmContacts_new(arena,nv,faces,nf,uv,1e-7); AsmContactStats contacts; memset(&contacts, 0, sizeof contacts);
    if (!query || AsmContacts_measure(query,uv,NULL,NULL,NULL,&contacts) || contacts.pairs) {
        st->fail_reason = ASM_FLAT_FOLDED; st->fail_stage = 5; st->fail_contacts = query ? contacts.pairs : 0;
        Arena_restore(arena,mark); return -1;
    }
    st->ok = 1;
    Arena_restore(arena, mark);
    return 0;
}

int AsmFlatten_face_degenerate(const float *xyz, const int32_t *face)
{
    AfRef r;
    af_reference(xyz, face, 1, &r);
    return !(r.area > 0) || !isfinite(r.area);
}

size_t AsmFlatten_drop_degenerate(const float *xyz, size_t nv, const int32_t *faces, size_t nf,
                                  int32_t *out_faces, int32_t *out_vmap, size_t *out_nv)
{
    for (size_t i = 0; i < nv; i++) out_vmap[i] = -1;
    size_t kept = 0;
    for (size_t f = 0; f < nf; f++) {
        if (AsmFlatten_face_degenerate(xyz, faces + 3*f)) continue;
        memcpy(out_faces + 3*kept, faces + 3*f, 3 * sizeof(int32_t));
        kept++;
    }
    /* referenced vertices keep their relative order */
    for (size_t k = 0; k < kept * 3; k++) out_vmap[out_faces[k]] = 0;
    size_t n = 0;
    for (size_t i = 0; i < nv; i++) if (out_vmap[i] >= 0) out_vmap[i] = (int32_t)n++;
    for (size_t k = 0; k < kept * 3; k++) out_faces[k] = out_vmap[out_faces[k]];
    *out_nv = n;
    return kept;
}

/* ---- hole cutting ------------------------------------------------------------ */

/* binary min-heap of (key, vertex) for Dijkstra */
typedef struct { double *key; int32_t *v; size_t n; } AfHeap;

static void af_heap_push(AfHeap *h, double k, int32_t v)
{
    size_t i = h->n++;
    while (i > 0) {
        size_t p = (i - 1) / 2;
        if (h->key[p] <= k) break;
        h->key[i] = h->key[p]; h->v[i] = h->v[p]; i = p;
    }
    h->key[i] = k; h->v[i] = v;
}

static int32_t af_heap_pop(AfHeap *h, double *k)
{
    int32_t top = h->v[0];
    *k = h->key[0];
    double lk = h->key[--h->n];
    int32_t lv = h->v[h->n];
    size_t i = 0;
    for (;;) {
        size_t c = 2 * i + 1;
        if (c >= h->n) break;
        if (c + 1 < h->n && h->key[c + 1] < h->key[c]) c++;
        if (h->key[c] >= lk) break;
        h->key[i] = h->key[c]; h->v[i] = h->v[c]; i = c;
    }
    h->key[i] = lk; h->v[i] = lv;
    return top;
}

/* ---- crumpled charts ---------------------------------------------------------- */

/* the root of a with path halving */
static int32_t af_find(int32_t *parent, int32_t a)
{
    while (parent[a] != a) { parent[a] = parent[parent[a]]; a = parent[a]; }
    return a;
}

static void af_union(int32_t *parent, int32_t a, int32_t b)
{
    a = af_find(parent, a); b = af_find(parent, b);
    if (a != b) parent[a < b ? b : a] = a < b ? a : b;
}

/* Multi-source Dijkstra over the chart's edges (3-D length) from the vertices marked in src to the
 * nearest vertex marked in dst.  The faces on the path's vertices are marked in excise and the path's
 * vertices in onpath.  Returns the target vertex, or -1 when none is reachable. */
static int32_t af_join(const AfTopo *t, const float *xyz, size_t nv, const int32_t *faces, size_t nf,
                       const uint8_t *src, const uint8_t *dst, double *dist, int32_t *pred, uint8_t *done,
                       AfHeap *h, uint8_t *onpath, uint8_t *excise)
{
    const int32_t *off = CSR_offset(t->adj), *col = CSR_target(t->adj);
    memset(onpath, 0, nv);
    for (size_t v = 0; v < nv; v++) { dist[v] = INFINITY; pred[v] = -1; done[v] = 0; }
    h->n = 0;
    for (size_t v = 0; v < nv; v++) if (src[v]) { dist[v] = 0.0; af_heap_push(h, 0.0, (int32_t)v); }
    int32_t target = -1;
    while (h->n) {
        double d;
        int32_t a = af_heap_pop(h, &d);
        if (done[a] || d > dist[a]) continue;
        done[a] = 1;
        if (dst[a]) { target = a; break; }
        for (int32_t k = off[a]; k < off[a + 1]; k++) {
            int32_t b = col[k];
            if (done[b]) continue;
            double pa[3], pb[3], e[3];
            af_get3(xyz, a, pa); af_get3(xyz, b, pb); af_sub3(pb, pa, e);
            double nd = d + af_len3(e);
            if (nd < dist[b]) { dist[b] = nd; pred[b] = a; af_heap_push(h, nd, b); }
        }
    }
    if (target < 0) return -1;
    for (int32_t v = target; v >= 0; v = pred[v]) onpath[v] = 1;
    for (size_t f = 0; f < nf; f++)
        if (onpath[faces[3*f]] || onpath[faces[3*f+1]] || onpath[faces[3*f+2]]) excise[f] = 1;
    return target;
}

size_t AsmFlatten_crumple_faces(Arena_T arena, const float *xyz, size_t nv, const int32_t *faces, size_t nf,
                                double cone_rad, int ring, int cut_holes, const uint8_t *seed, uint8_t *excise)
{
    memset(excise, 0, nf);
    Arena_Mark mark = Arena_save(arena);
    AfTopo t;
    if (af_topology(arena, xyz, nv, faces, nf, &t) != ASM_FLAT_OK) { Arena_restore(arena, mark); return 0; }
    /* angle sums, and the cone and pinch vertices */
    double *angle = ARENA_CALLOC(arena, nv, sizeof(double));
    for (size_t f = 0; f < nf; f++) {
        double p[3][3];
        for (int c = 0; c < 3; c++) af_get3(xyz, faces[3*f + (size_t)c], p[c]);
        for (int c = 0; c < 3; c++) angle[faces[3*f + (size_t)c]] += af_angle(p[c], p[(c+1)%3], p[(c+2)%3]);
    }
    int32_t *group = ARENA_ALLOC(arena, nv * sizeof(int32_t));
    af_pinch_copies(arena, xyz, nv, &t, group);
    uint8_t *vmark = ARENA_CALLOC(arena, nv, 1);
    size_t cones = 0;
    for (size_t v = 0; v < nv; v++) {
        int cone = t.on_boundary[v] ? (angle[v] > 2.0*AF_PI + 0.05 || group[v] != 0)
                                    : fabs(2.0*AF_PI - angle[v]) > cone_rad;
        if (cone) { vmark[v] = 1; cones++; }
    }
    size_t seeds = 0;
    for (size_t f = 0; seed && f < nf; f++) seeds += seed[f] != 0;
    if (!cones && !seeds && !(cut_holes && t.n_loops > 1)) { Arena_restore(arena, mark); return 0; }
    for (size_t f = 0; f < nf; f++)
        excise[f] = (uint8_t)(vmark[faces[3*f]] || vmark[faces[3*f+1]] || vmark[faces[3*f+2]] || (seed && seed[f]));
    for (int r = 0; r < ring; r++) {
        for (size_t f = 0; f < nf; f++) if (excise[f]) for (int c = 0; c < 3; c++) vmark[faces[3*f + (size_t)c]] = 1;
        for (size_t f = 0; f < nf; f++)
            if (vmark[faces[3*f]] || vmark[faces[3*f+1]] || vmark[faces[3*f+2]]) excise[f] = 1;
    }
    /* With cut_holes every hole loop is first joined to the outer loop along a spanning tree of
     * shortest paths (a band or holed chart whose hole carries holonomy cannot flatten; the joined
     * region is a disk).  Each excised cluster that does not reach the joined boundary is then joined
     * to it, so no hole is left around net curvature. */
    int32_t *parent = ARENA_ALLOC(arena, nv * sizeof(int32_t));
    for (size_t v = 0; v < nv; v++) parent[v] = (int32_t)v;
    memset(vmark, 0, nv);
    for (size_t f = 0; f < nf; f++) if (excise[f]) {
        af_union(parent, faces[3*f], faces[3*f+1]);
        af_union(parent, faces[3*f], faces[3*f+2]);
        for (int c = 0; c < 3; c++) vmark[faces[3*f + (size_t)c]] = 1;
    }
    double *dist = ARENA_ALLOC(arena, nv * sizeof(double));
    int32_t *pred = ARENA_ALLOC(arena, nv * sizeof(int32_t));
    uint8_t *done = ARENA_ALLOC(arena, nv), *onpath = ARENA_ALLOC(arena, nv), *src = ARENA_ALLOC(arena, nv);
    size_t heap_cap = (size_t)CSR_nnz(t.adj) + nv;
    AfHeap h = { ARENA_ALLOC(arena, heap_cap * sizeof(double)), ARENA_ALLOC(arena, heap_cap * sizeof(int32_t)), 0 };
    uint8_t *joined = ARENA_CALLOC(arena, nv, 1);             /* the boundary a cut may end on */
    for (size_t v = 0; v < nv; v++)
        if (t.on_boundary[v] && (!cut_holes || t.loop_id[v] == t.outer)) joined[v] = 1;
    if (cut_holes) {
        uint8_t *loop_joined = ARENA_CALLOC(arena, t.n_loops, 1);
        loop_joined[t.outer] = 1;
        for (size_t round = 0; round < t.n_loops; round++) {
            memset(src, 0, nv);
            size_t pending = 0;
            for (size_t v = 0; v < nv; v++)
                if (t.on_boundary[v] && !loop_joined[t.loop_id[v]]) { src[v] = 1; pending++; }
            if (!pending) break;
            int32_t target = af_join(&t, xyz, nv, faces, nf, src, joined, dist, pred, done, &h, onpath, excise);
            if (target < 0) break;
            int32_t start = target;
            while (pred[start] >= 0) start = pred[start];
            loop_joined[t.loop_id[start]] = 1;
            for (size_t v = 0; v < nv; v++)
                if (onpath[v] || (t.on_boundary[v] && t.loop_id[v] == t.loop_id[start])) joined[v] = 1;
        }
    }
    uint8_t *touches = ARENA_CALLOC(arena, nv, 1);            /* per cluster root: reaches the joined boundary */
    for (size_t v = 0; v < nv; v++) if (vmark[v] && joined[v]) touches[af_find(parent, (int32_t)v)] = 1;
    for (size_t v0 = 0; v0 < nv; v0++) {
        if (!vmark[v0] || af_find(parent, (int32_t)v0) != (int32_t)v0 || touches[v0]) continue;
        for (size_t v = 0; v < nv; v++) src[v] = (uint8_t)(vmark[v] && af_find(parent, (int32_t)v) == (int32_t)v0);
        af_join(&t, xyz, nv, faces, nf, src, joined, dist, pred, done, &h, onpath, excise);
    }
    /* a vertex left with a split fan (two boundary gaps) is a pinch: clear its fan */
    int32_t *outgoing = ARENA_ALLOC(arena, nv * sizeof(int32_t));
    for (int pass = 0; pass < 64; pass++) {
        Arena_Mark pass_mark = Arena_save(arena);
        size_t kept = 0;
        for (size_t f = 0; f < nf; f++) kept += !excise[f];
        if (!kept) break;
        int32_t *kf = ARENA_ALLOC(arena, kept * 3 * sizeof(int32_t));
        kept = 0;
        for (size_t f = 0; f < nf; f++) if (!excise[f]) memcpy(kf + 3*kept++, faces + 3*f, 3 * sizeof(int32_t));
        CSR_T adj = CSR_from_faces(arena, kf, kept, nv);
        const int32_t *ko = CSR_offset(adj), *kc = CSR_target(adj);
        int32_t *cnt = ARENA_CALLOC(arena, (size_t)CSR_nnz(adj), sizeof(int32_t));
        for (size_t f = 0; f < kept; f++) for (int e = 0; e < 3; e++) {
            int32_t a = kf[3*f + (size_t)e], b = kf[3*f + (size_t)((e+1)%3)];
            for (int32_t k = ko[a]; k < ko[a+1]; k++) if (kc[k] == b) { cnt[k]++; break; }
            for (int32_t k = ko[b]; k < ko[b+1]; k++) if (kc[k] == a) { cnt[k]++; break; }
        }
        for (size_t v = 0; v < nv; v++) outgoing[v] = 0;
        for (size_t f = 0; f < kept; f++) for (int e = 0; e < 3; e++) {
            int32_t a = kf[3*f + (size_t)e], b = kf[3*f + (size_t)((e+1)%3)];
            for (int32_t k = ko[a]; k < ko[a+1]; k++) if (kc[k] == b) { if (cnt[k] == 1) outgoing[a]++; break; }
        }
        int changed = 0;
        for (size_t f = 0; f < nf; f++) if (!excise[f])
            for (int c = 0; c < 3; c++) if (outgoing[faces[3*f + (size_t)c]] > 1) { excise[f] = 1; changed = 1; break; }
        Arena_restore(arena, pass_mark);
        if (!changed) break;
    }
    size_t n = 0;
    for (size_t f = 0; f < nf; f++) n += excise[f];
    Arena_restore(arena, mark);
    return n;
}

typedef struct AfMarkPairs { uint8_t *mark; size_t n; } AfMarkPairs;

static void af_mark_pair(void *context, size_t a, size_t b, double area)
{
    AfMarkPairs *m = context; (void)area;
    if (!m->mark[a]) { m->mark[a] = 1; m->n++; }
    if (!m->mark[b]) { m->mark[b] = 1; m->n++; }
}

/* The faces that stop the chart's conformal seed (pinch copies opened) from being an embedding: the
 * faces it folds, or else the faces it puts in contact.  mark[nf]; returns their number (0 when the
 * seed cannot be solved or is already clean). */
static size_t af_seed_obstruction(Arena_T arena, const float *xyz, size_t nv, const int32_t *faces, size_t nf, uint8_t *mark)
{
    memset(mark, 0, nf);
    Arena_Mark m = Arena_save(arena);
    AfTopo t;
    size_t n = 0;
    if (af_topology(arena, xyz, nv, faces, nf, &t) == ASM_FLAT_OK) {
        AfRef *ref = ARENA_ALLOC(arena, nf * sizeof(AfRef));
        af_reference(xyz, faces, nf, ref);
        double area = 0.0;
        for (size_t f = 0; f < nf; f++) area += ref[f].area;
        double *uv = ARENA_ALLOC(arena, 2 * nv * sizeof(double));
        if (area > 0.0 && af_lscm(arena, xyz, nv, faces, nf, &t, ref, area, uv) == ASM_FLAT_OK) {
            if (!UvGuard_mesh(uv, faces, nf)) {
                for (size_t f = 0; f < nf; f++) if (!UvGuard_positive(uv, faces + 3*f)) { mark[f] = 1; n++; }
            } else {
                af_open_pinches(arena, xyz, nv, faces, nf, &t, uv);
                AsmContacts *q = AsmContacts_new(arena, nv, faces, nf, uv, 1e-7);
                AfMarkPairs mp = { mark, 0 };
                AsmContactStats cs; memset(&cs, 0, sizeof cs);
                if (q) AsmContacts_measure(q, uv, NULL, af_mark_pair, &mp, &cs);
                n = mp.n;
            }
        }
    }
    Arena_restore(arena, m);
    return n;
}

/* One rescue: excise (tiers of rings and hole cuts around the cones and the seed faces), flatten each
 * piece certified; at depth 0 a piece that fails is rescued again from the faces its conformal seed
 * folds or puts in contact.  Writes piece[nf] (-1 = not in a certified piece, else first_piece + k)
 * and the pieces' uv; returns the number of certified pieces. */
static size_t af_rescue(Arena_T arena, const float *xyz, size_t nv, const int32_t *faces, size_t nf,
                        double band, int max_iters, double min_piece_area, int depth, const uint8_t *seed,
                        int32_t first_piece, int32_t *piece, float *uv, AsmFlattenStats *st)
{
    /* rings of excision around the cones, and whether every hole is also cut to the outer loop */
    static const struct { int ring, cut_holes; } tiers[] = { {1, 0}, {1, 1}, {3, 1} };
    Arena_Mark mark = Arena_save(arena);
    double total = 0.0;
    for (size_t f = 0; f < nf; f++) { AfRef r; af_reference(xyz, faces + 3*f, 1, &r); if (r.area > 0) total += r.area; }
    int32_t *trial_piece = ARENA_ALLOC(arena, nf * sizeof(int32_t));
    float *trial_uv = ARENA_ALLOC(arena, nv * 2 * sizeof(float));
    uint8_t *excise = ARENA_ALLOC(arena, nf);
    double best_area = 0.0;
    size_t best_pieces = 0;
    for (size_t f = 0; f < nf; f++) piece[f] = -1;
    for (size_t k = 0; k < sizeof tiers / sizeof tiers[0]; k++) {
        if (af_expired()) break;
        Arena_Mark m2 = Arena_save(arena);
        if (!AsmFlatten_crumple_faces(arena, xyz, nv, faces, nf, 0.3, tiers[k].ring, tiers[k].cut_holes, seed, excise)) {
            Arena_restore(arena, m2);
            continue;
        }
        /* pieces: connected components of the kept faces, largest area first */
        int32_t *parent = ARENA_ALLOC(arena, nv * sizeof(int32_t));
        for (size_t v = 0; v < nv; v++) parent[v] = (int32_t)v;
        for (size_t f = 0; f < nf; f++) if (!excise[f]) {
            af_union(parent, faces[3*f], faces[3*f+1]);
            af_union(parent, faces[3*f], faces[3*f+2]);
        }
        double *carea = ARENA_CALLOC(arena, nv, sizeof(double));
        for (size_t f = 0; f < nf; f++) if (!excise[f]) {
            AfRef r; af_reference(xyz, faces + 3*f, 1, &r);
            carea[af_find(parent, faces[3*f])] += r.area > 0 ? r.area : 0.0;
        }
        for (size_t f = 0; f < nf; f++) trial_piece[f] = -1;
        int32_t *vmap = ARENA_ALLOC(arena, nv * sizeof(int32_t));
        uint8_t *tried = ARENA_CALLOC(arena, nv, 1);
        double certified = 0.0;
        size_t pieces = 0;
        AsmFlattenStats first; memset(&first, 0, sizeof first);
        for (;;) {
            /* the largest untried piece */
            int32_t root = -1;
            for (size_t v = 0; v < nv; v++)
                if (!tried[v] && parent[v] == (int32_t)v && carea[v] >= min_piece_area && (root < 0 || carea[v] > carea[root])) root = (int32_t)v;
            if (root < 0 || af_expired()) break;
            tried[root] = 1;
            Arena_Mark m3 = Arena_save(arena);
            for (size_t v = 0; v < nv; v++) vmap[v] = -1;
            size_t pf = 0, pv = 0;
            for (size_t f = 0; f < nf; f++) if (!excise[f] && af_find(parent, faces[3*f]) == root) {
                pf++;
                for (int c = 0; c < 3; c++) vmap[faces[3*f + (size_t)c]] = 0;
            }
            for (size_t v = 0; v < nv; v++) if (vmap[v] >= 0) vmap[v] = (int32_t)pv++;
            int32_t *lf = ARENA_ALLOC(arena, pf * 3 * sizeof(int32_t)), *lidx = ARENA_ALLOC(arena, pf * sizeof(int32_t));
            float *lx = ARENA_ALLOC(arena, pv * 3 * sizeof(float)), *luv = ARENA_ALLOC(arena, pv * 2 * sizeof(float));
            pf = 0;
            for (size_t f = 0; f < nf; f++) if (!excise[f] && af_find(parent, faces[3*f]) == root) {
                for (int c = 0; c < 3; c++) lf[3*pf + (size_t)c] = vmap[faces[3*f + (size_t)c]];
                lidx[pf++] = (int32_t)f;
            }
            for (size_t v = 0; v < nv; v++) if (vmap[v] >= 0) memcpy(lx + 3*(size_t)vmap[v], xyz + 3*v, 3 * sizeof(float));
            AsmFlattenStats ps;
            if (AsmFlatten_chart_certified(arena, lx, pv, lf, pf, band, max_iters, luv, &ps) == 0) {
                for (size_t q = 0; q < pf; q++) trial_piece[lidx[q]] = first_piece + (int32_t)pieces;
                for (size_t v = 0; v < nv; v++) if (vmap[v] >= 0) memcpy(trial_uv + 2*v, luv + 2*(size_t)vmap[v], 2 * sizeof(float));
                if (!pieces) first = ps;
                pieces++;
                certified += carea[root];
            } else if (depth == 0) {
                /* the piece's own obstruction: what its conformal seed folds or puts in contact */
                uint8_t *obst = ARENA_ALLOC(arena, pf);
                if (af_seed_obstruction(arena, lx, pv, lf, pf, obst)) {
                    int32_t *sub = ARENA_ALLOC(arena, pf * sizeof(int32_t));
                    float *suv = ARENA_ALLOC(arena, pv * 2 * sizeof(float));
                    AsmFlattenStats ss;
                    size_t n = af_rescue(arena, lx, pv, lf, pf, band, max_iters, min_piece_area, 1, obst,
                                         first_piece + (int32_t)pieces, sub, suv, &ss);
                    if (n) {
                        uint8_t *vin = ARENA_CALLOC(arena, pv, 1);
                        for (size_t q = 0; q < pf; q++) if (sub[q] >= 0) {
                            trial_piece[lidx[q]] = sub[q];
                            for (int c = 0; c < 3; c++) vin[lf[3*q + (size_t)c]] = 1;
                            AfRef r; af_reference(lx, lf + 3*q, 1, &r);
                            certified += r.area > 0 ? r.area : 0.0;
                        }
                        for (size_t v = 0; v < nv; v++) if (vmap[v] >= 0 && vin[vmap[v]])
                            memcpy(trial_uv + 2*v, suv + 2*(size_t)vmap[v], 2 * sizeof(float));
                        if (!pieces) first = ss;
                        pieces += n;
                    }
                }
            }
            Arena_restore(arena, m3);
        }
        if (certified > best_area) {
            best_area = certified; best_pieces = pieces; *st = first;
            memcpy(piece, trial_piece, nf * sizeof(int32_t));
            memcpy(uv, trial_uv, nv * 2 * sizeof(float));
        }
        Arena_restore(arena, m2);
        if (best_area >= 0.9 * total) break;
    }
    Arena_restore(arena, mark);
    return best_pieces;
}

int AsmFlatten_rescue_crumpled(Arena_T arena, const float *xyz, size_t nv, const int32_t *faces, size_t nf,
                               double band, int max_iters, double min_piece_area,
                               int32_t *piece, float *uv, size_t *n_pieces, AsmFlattenStats *st)
{
    AsmFlattenStats first; memset(&first, 0, sizeof first);
    int direct = af_direct;
    af_direct = 1;
    size_t n = af_rescue(arena, xyz, nv, faces, nf, band, max_iters, min_piece_area, 0, NULL, 0, piece, uv, &first);
    af_direct = direct;
    *n_pieces = n;
    if (!n) return -1;
    *st = first;
    return 0;
}

/* ---- selftest ---------------------------------------------------------------- */

typedef struct AfFix {
    float   *xyz; int32_t *faces; size_t nv, nf;
} AfFix;

/* grid of (nx+1)*(ny+1) vertices mapped by fn; optional hole removes faces
 * whose centre lies within hole_r of (hole_cx, hole_cy) in grid units */
static void af_fix_grid(Arena_T arena, int nx, int ny, double sx, double sy,
                        void (*fn)(double, double, float *), int hole,
                        double hole_cx, double hole_cy, double hole_r, AfFix *out)
{
    size_t nv = (size_t)(nx+1) * (size_t)(ny+1);
    out->xyz = ARENA_ALLOC(arena, nv * 3 * sizeof(float));
    out->faces = ARENA_ALLOC(arena, (size_t)nx * (size_t)ny * 2 * 3 * sizeof(int32_t));
    for (int j = 0; j <= ny; j++)
        for (int i = 0; i <= nx; i++)
            fn(sx * i / nx, sy * j / ny, &out->xyz[((size_t)j*(size_t)(nx+1)+(size_t)i)*3]);
    size_t nf = 0;
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            if (hole) {
                double cx = i + 0.5, cy = j + 0.5;
                if ((cx-hole_cx)*(cx-hole_cx) + (cy-hole_cy)*(cy-hole_cy) < hole_r*hole_r) continue;
            }
            int32_t a = (int32_t)(j*(nx+1)+i), b = a+1, c = a+(nx+1), d = c+1;
            /* consistent CCW in (x,y) parameter space */
            out->faces[nf*3] = a; out->faces[nf*3+1] = b; out->faces[nf*3+2] = d; nf++;
            out->faces[nf*3] = a; out->faces[nf*3+1] = d; out->faces[nf*3+2] = c; nf++;
        }
    }
    /* compact vertices: drop unreferenced */
    int32_t *remap = ARENA_ALLOC(arena, nv * sizeof(int32_t));
    for (size_t i = 0; i < nv; i++) remap[i] = -1;
    size_t kept = 0;
    for (size_t k = 0; k < nf*3; k++) if (remap[out->faces[k]] < 0) remap[out->faces[k]] = (int32_t)kept++;
    float *xyz2 = ARENA_ALLOC(arena, kept * 3 * sizeof(float));
    for (size_t i = 0; i < nv; i++) if (remap[i] >= 0) memcpy(&xyz2[(size_t)remap[i]*3], &out->xyz[i*3], 3*sizeof(float));
    for (size_t k = 0; k < nf*3; k++) out->faces[k] = remap[out->faces[k]];
    out->xyz = xyz2; out->nv = kept; out->nf = nf;
}

static void af_fn_plane(double x, double y, float *o) { o[0] = 0.0f; o[1] = (float)y; o[2] = (float)x; }
static void af_fn_cyl(double x, double y, float *o)
{   /* x = arc length around a radius-40 cylinder, y = axis */
    double r = 40.0, ang = x / r;
    o[0] = (float)y; o[1] = (float)(r * sin(ang)); o[2] = (float)(r * cos(ang));
}
static void af_fn_cone(double x, double y, float *o)
{   /* developable cone: (s, phi) with s in [20,100], phi = x/100 * 2pi*0.6 */
    double s = 20.0 + y, phi = x / 100.0 * 2.0 * AF_PI * 0.6, k = 0.6;
    o[0] = (float)(s * sqrt(1.0 - k*k)); o[1] = (float)(s * k * sin(phi/k)); o[2] = (float)(s * k * cos(phi/k));
}
static void af_fn_sphere(double x, double y, float *o)
{   /* spherical cap: lon/lat patch on a radius-60 sphere, +/-35 degrees */
    double R = 60.0, lon = (x/100.0 - 0.5) * 1.2, lat = (y/100.0 - 0.5) * 1.2;
    o[0] = (float)(R * sin(lat)); o[1] = (float)(R * cos(lat) * sin(lon)); o[2] = (float)(R * cos(lat) * cos(lon));
}

static void af_fn_uturn(double x, double y, float *o)
{   /* a crushed-scroll layer turning back inside one cube: two 60-vox arms
     * joined by a radius-6 half turn, 40 vox along the turn axis */
    double r = 6.0, arm = 60.0, bend = AF_PI * r, p, q;
    if (x < arm) { p = x - arm; q = r; }
    else if (x < arm + bend) { double a = (x - arm) / r; p = r * sin(a); q = r * cos(a); }
    else { p = -(x - arm - bend); q = -r; }
    o[0] = (float)y; o[1] = (float)q; o[2] = (float)p;
}

static int af_check(const char *name, int cond, int *fails)
{
    if (!cond) { fprintf(stderr, "  asm_flatten selftest FAIL: %s\n", name); (*fails)++; }
    return cond;
}

/* A unit planar grid of nx x ny cells cut along the grid line x = xs for y in [y0, y1]: every
 * vertex strictly inside the cut, and an end on the grid boundary, is split into a left and a right
 * copy at the SAME position (as the mesher leaves a zero-width slit or a split bowtie).  With
 * `close` the cut must be one grid edge pair (y1 = y0 + 2) and is closed again by the two
 * zero-area triangles the mesher uses, so the chart also carries two collinear faces. */
static void af_fix_slit(Arena_T arena, int nx, int ny, int xs, int y0, int y1, int close, AfFix *out)
{
    /* row-major grid (af_fix_grid renumbers vertices by first appearance) */
    size_t row = (size_t)nx + 1, nv = row * (size_t)(ny + 1);
    out->nv = nv; out->nf = (size_t)nx * (size_t)ny * 2;
    out->xyz = ARENA_ALLOC(arena, nv * 3 * sizeof(float));
    out->faces = ARENA_ALLOC(arena, out->nf * 3 * sizeof(int32_t));
    for (int j = 0; j <= ny; j++) for (int i = 0; i <= nx; i++) af_fn_plane((double)i, (double)j, out->xyz + 3*((size_t)j * row + (size_t)i));
    {
        size_t f = 0;
        for (int j = 0; j < ny; j++) for (int i = 0; i < nx; i++) {
            int32_t a = (int32_t)((size_t)j * row + (size_t)i), b = a + 1, c = a + (int32_t)row, d = c + 1;
            int32_t tri[6] = { a, b, d, a, d, c };
            memcpy(out->faces + 3*f, tri, sizeof tri);
            f += 2;
        }
    }
    int32_t *copy = ARENA_ALLOC(arena, nv * sizeof(int32_t));
    for (size_t i = 0; i < nv; i++) copy[i] = -1;
    size_t extra = 0;
    for (int j = y0; j <= y1; j++) {
        int end = (j == y0 || j == y1);
        if (end && j != 0 && j != ny) continue;            /* an interior end is the slit's tip */
        copy[(size_t)j * row + (size_t)xs] = (int32_t)(nv + extra++);
    }
    size_t nf = out->nf + (close ? 2 : 0);
    float *xyz = ARENA_ALLOC(arena, (nv + extra) * 3 * sizeof(float));
    int32_t *faces = ARENA_ALLOC(arena, nf * 3 * sizeof(int32_t));
    memcpy(xyz, out->xyz, nv * 3 * sizeof(float));
    for (size_t i = 0; i < nv; i++) if (copy[i] >= 0) memcpy(xyz + 3*(size_t)copy[i], out->xyz + 3*i, 3 * sizeof(float));
    memcpy(faces, out->faces, out->nf * 3 * sizeof(int32_t));
    for (size_t f = 0; f < out->nf; f++) {
        double cx = 0.0, cy = 0.0;
        for (int c = 0; c < 3; c++) { cx += out->xyz[3*(size_t)faces[3*f + (size_t)c] + 2] / 3.0; cy += out->xyz[3*(size_t)faces[3*f + (size_t)c] + 1] / 3.0; }
        if (cx < (double)xs || cy < (double)y0 || cy > (double)y1) continue;   /* the left side keeps the originals */
        for (int c = 0; c < 3; c++) if (copy[faces[3*f + (size_t)c]] >= 0) faces[3*f + (size_t)c] = copy[faces[3*f + (size_t)c]];
    }
    if (close) {
        /* the slit a - m(left) - b - m'(right): two collinear triangles fill it, wound against the
         * half-edges already present */
        int32_t a = (int32_t)((size_t)y0 * row + (size_t)xs), b = (int32_t)((size_t)y1 * row + (size_t)xs);
        int32_t m = (int32_t)((size_t)(y0 + 1) * row + (size_t)xs), mc = copy[m];
        int am = 0;                                        /* does a real face hold a->m? */
        for (size_t f = 0; f < out->nf && !am; f++)
            for (int e = 0; e < 3; e++)
                if (faces[3*f + (size_t)e] == a && faces[3*f + (size_t)((e + 1) % 3)] == m) am = 1;
        int32_t fill[6] = { a, m, b, a, b, mc };
        if (am) { fill[1] = b; fill[2] = m; fill[4] = mc; fill[5] = b; }
        memcpy(faces + 3*out->nf, fill, sizeof fill);
    }
    out->xyz = xyz; out->faces = faces; out->nv = nv + extra; out->nf = nf;
}

/* A developable cone disk around its apex: angle defect `defect`, `rings` rings of unit spacing,
 * `sectors` vertices per ring, apex at vertex 0. */
static void af_fix_cone_disk(Arena_T arena, int rings, int sectors, double defect, AfFix *out)
{
    size_t ns = (size_t)sectors;
    out->nv = 1 + (size_t)rings * ns;
    out->nf = ns + 2 * (size_t)(rings - 1) * ns;
    out->xyz = ARENA_ALLOC(arena, out->nv * 3 * sizeof(float));
    out->faces = ARENA_ALLOC(arena, out->nf * 3 * sizeof(int32_t));
    double k = (2.0 * AF_PI - defect) / (2.0 * AF_PI);
    out->xyz[0] = out->xyz[1] = out->xyz[2] = 0.0f;
    for (int i = 1; i <= rings; i++) for (int j = 0; j < sectors; j++) {
        double phi = 2.0 * AF_PI * j / sectors, rho = k * i, h = i * sqrt(1.0 - k * k);
        float *p = out->xyz + 3 * (1 + (size_t)(i - 1) * ns + (size_t)j);
        p[0] = (float)h; p[1] = (float)(rho * sin(phi)); p[2] = (float)(rho * cos(phi));
    }
    size_t f = 0;
    for (int j = 0; j < sectors; j++) {
        int32_t a = 1 + j, b = 1 + (j + 1) % sectors;
        out->faces[3*f] = 0; out->faces[3*f+1] = a; out->faces[3*f+2] = b; f++;
    }
    for (int i = 1; i < rings; i++) for (int j = 0; j < sectors; j++) {
        int32_t p = (int32_t)(1 + (size_t)(i - 1) * ns + (size_t)j), q = (int32_t)(1 + (size_t)(i - 1) * ns + (size_t)((j + 1) % sectors));
        int32_t r = p + sectors, s = q + sectors;
        int32_t tri[6] = { p, r, s, p, s, q };
        memcpy(out->faces + 3*f, tri, sizeof tri);
        f += 2;
    }
}

/* Crumpled charts: the excision marks a cone's fan and joins it to the boundary; the rest flattens. */
static int af_crumple_selftest(void)
{
    int fails = 0;
    Arena_T arena = Arena_new();
    AfFix fx;
    af_fix_cone_disk(arena, 12, 48, 1.5, &fx);
    uint8_t *ex = ARENA_ALLOC(arena, fx.nf);
    size_t n = AsmFlatten_crumple_faces(arena, fx.xyz, fx.nv, fx.faces, fx.nf, 0.3, 1, 0, NULL, ex);
    int apex = 1, reaches = 0;
    for (size_t f = 0; f < fx.nf; f++) {
        for (int c = 0; c < 3; c++) {
            if (fx.faces[3*f + (size_t)c] == 0 && !ex[f]) apex = 0;
            if (ex[f] && fx.faces[3*f + (size_t)c] >= (int32_t)(1 + 11 * 48)) reaches = 1;
        }
    }
    af_check("cone fan excised", n > 48 && n < fx.nf / 4 && apex, &fails);
    af_check("cone excision joined to the boundary", reaches, &fails);
    int32_t *piece = ARENA_ALLOC(arena, fx.nf * sizeof(int32_t));
    float *uv = ARENA_ALLOC(arena, fx.nv * 2 * sizeof(float));
    size_t np = 0;
    AsmFlattenStats st;
    int rc = AsmFlatten_rescue_crumpled(arena, fx.xyz, fx.nv, fx.faces, fx.nf, 0.05, 30, 10.0, piece, uv, &np, &st);
    double total = 0.0, kept = 0.0;
    for (size_t f = 0; f < fx.nf; f++) {
        AfRef r; af_reference(fx.xyz, fx.faces + 3*f, 1, &r);
        total += r.area;
        if (piece[f] == 0) kept += r.area;
        if (piece[f] >= 0 && ex[f]) apex = 0;                 /* an excised face never joins a piece */
    }
    fprintf(stderr, "  asm_flatten selftest: cone rescue rc %d pieces %zu, certified %.1f%% of the area, sigma [%.3f, %.3f]\n",
            rc, np, 100.0 * kept / total, st.sigma_lo, st.sigma_hi);
    af_check("cone rescued", rc == 0 && np == 1 && kept >= 0.8 * total && apex, &fails);
    /* the rescued piece's map carries the certificate */
    if (rc == 0) {
        int32_t *vm = ARENA_ALLOC(arena, fx.nv * sizeof(int32_t)), *pf = ARENA_ALLOC(arena, fx.nf * 3 * sizeof(int32_t));
        size_t m = 0, pv = 0;
        for (size_t v = 0; v < fx.nv; v++) vm[v] = -1;
        for (size_t f = 0; f < fx.nf; f++) if (piece[f] == 0) for (int c = 0; c < 3; c++) vm[fx.faces[3*f + (size_t)c]] = 0;
        for (size_t v = 0; v < fx.nv; v++) if (vm[v] >= 0) vm[v] = (int32_t)pv++;
        float *px = ARENA_ALLOC(arena, pv * 3 * sizeof(float)), *pu = ARENA_ALLOC(arena, pv * 2 * sizeof(float));
        for (size_t v = 0; v < fx.nv; v++) if (vm[v] >= 0) {
            memcpy(px + 3*(size_t)vm[v], fx.xyz + 3*v, 3 * sizeof(float));
            memcpy(pu + 2*(size_t)vm[v], uv + 2*v, 2 * sizeof(float));
        }
        for (size_t f = 0; f < fx.nf; f++) if (piece[f] == 0) { for (int c = 0; c < 3; c++) pf[3*m + (size_t)c] = vm[fx.faces[3*f + (size_t)c]]; m++; }
        af_check("cone piece certified", af_chart_certified(px, pf, m, pu), &fails);
    }
    /* a closed band (a full turn of a radius-12 cylinder, 16 high) has no cone and cannot flatten;
     * cutting its hole to the outer loop leaves a certifiable strip */
    {
        const int around = 64, up = 16;
        AfFix band;
        band.nv = (size_t)around * (size_t)(up + 1);
        band.nf = (size_t)around * (size_t)up * 2;
        band.xyz = ARENA_ALLOC(arena, band.nv * 3 * sizeof(float));
        band.faces = ARENA_ALLOC(arena, band.nf * 3 * sizeof(int32_t));
        for (int j = 0; j <= up; j++) for (int i = 0; i < around; i++) {
            double phi = 2.0 * AF_PI * i / around;
            float *p = band.xyz + 3 * ((size_t)j * (size_t)around + (size_t)i);
            p[0] = (float)j; p[1] = (float)(12.0 * sin(phi)); p[2] = (float)(12.0 * cos(phi));
        }
        size_t f = 0;
        for (int j = 0; j < up; j++) for (int i = 0; i < around; i++) {
            int32_t a = j * around + i, b = j * around + (i + 1) % around, c = a + around, d = b + around;
            int32_t tri[6] = { a, b, d, a, d, c };
            memcpy(band.faces + 3*f, tri, sizeof tri);
            f += 2;
        }
        uint8_t *bx = ARENA_ALLOC(arena, band.nf);
        af_check("band has no cone to excise", AsmFlatten_crumple_faces(arena, band.xyz, band.nv, band.faces, band.nf, 0.3, 1, 0, NULL, bx) == 0, &fails);
        size_t cut = AsmFlatten_crumple_faces(arena, band.xyz, band.nv, band.faces, band.nf, 0.3, 1, 1, NULL, bx);
        af_check("band cut from hole to outer loop", cut > 0 && cut < band.nf / 8, &fails);
        int32_t *bp = ARENA_ALLOC(arena, band.nf * sizeof(int32_t));
        float *buv = ARENA_ALLOC(arena, band.nv * 2 * sizeof(float));
        size_t bn = 0;
        AsmFlattenStats bs;
        int brc = AsmFlatten_rescue_crumpled(arena, band.xyz, band.nv, band.faces, band.nf, 0.05, 30, 10.0, bp, buv, &bn, &bs);
        double btotal = 0.0, bkept = 0.0;
        for (size_t g = 0; g < band.nf; g++) {
            AfRef r; af_reference(band.xyz, band.faces + 3*g, 1, &r);
            btotal += r.area;
            if (brc == 0 && bp[g] == 0) bkept += r.area;
        }
        fprintf(stderr, "  asm_flatten selftest: band rescue rc %d pieces %zu, certified %.1f%% of the area\n", brc, bn, 100.0 * bkept / btotal);
        af_check("band rescued by its hole cut", brc == 0 && bn == 1 && bkept >= 0.8 * btotal, &fails);
    }
    /* the conformal seed of a clean sheet has no obstruction to excise */
    af_fix_grid(arena, 8, 6, 8.0, 6.0, af_fn_plane, 0, 0, 0, 0, &fx);
    {
        uint8_t *ob = ARENA_ALLOC(arena, fx.nf);
        af_check("clean sheet seed unobstructed", af_seed_obstruction(arena, fx.xyz, fx.nv, fx.faces, fx.nf, ob) == 0, &fails);
    }
    /* a chart with no cone or pinch is left alone; a pinch copy counts as a cone */
    af_fix_grid(arena, 8, 6, 8.0, 6.0, af_fn_plane, 0, 0, 0, 0, &fx);
    ex = ARENA_ALLOC(arena, fx.nf);
    af_check("flat chart has nothing to excise", AsmFlatten_crumple_faces(arena, fx.xyz, fx.nv, fx.faces, fx.nf, 0.3, 1, 0, NULL, ex) == 0, &fails);
    af_fix_slit(arena, 8, 6, 4, 2, 4, 0, &fx);
    ex = ARENA_ALLOC(arena, fx.nf);
    af_check("pinched chart excises its copies' fans", AsmFlatten_crumple_faces(arena, fx.xyz, fx.nv, fx.faces, fx.nf, 0.3, 1, 0, NULL, ex) > 0, &fails);
    Arena_dispose(&arena);
    return fails;
}

/* Pinch copies and the zero-area slit closure (PHerc0343: 10 of 13 unflattenable charts). */
static int af_pinch_selftest(void)
{
    int fails = 0;
    for (int interior = 0; interior < 2; interior++) {
        Arena_T arena = Arena_new(); AfFix fx; AfTopo topo;
        af_fix_slit(arena, 16, 10, 8, interior ? 2 : 0, interior ? 8 : 5, 0, &fx);
        af_check("slit chart topology", af_topology(arena, fx.xyz, fx.nv, fx.faces, fx.nf, &topo) == ASM_FLAT_OK &&
                 topo.n_loops == (size_t)(interior ? 2 : 1), &fails);
        int32_t *group = ARENA_ALLOC(arena, fx.nv * sizeof(int32_t));
        af_check("slit copies found", af_pinch_copies(arena, fx.xyz, fx.nv, &topo, group) == 10, &fails);
        /* the near-isometric seed keeps the copies together and cannot be certified */
        double *seed = ARENA_ALLOC(arena, fx.nv * 2 * sizeof(double));
        AfRef *ref = ARENA_ALLOC(arena, fx.nf * sizeof(AfRef));
        af_reference(fx.xyz, fx.faces, fx.nf, ref);
        af_check("slit conformal seed", af_lscm(arena, fx.xyz, fx.nv, fx.faces, fx.nf, &topo, ref, 160.0, seed) == ASM_FLAT_OK &&
                 UvGuard_mesh(seed, fx.faces, fx.nf), &fails);
        AsmContacts *q = AsmContacts_new(arena, fx.nv, fx.faces, fx.nf, seed, 1e-7);
        AsmContactStats cs; memset(&cs, 0, sizeof cs);
        AsmContacts_measure(q, seed, NULL, NULL, NULL, &cs);
        size_t before = cs.pairs;
        af_check("slit copies opened", af_open_pinches(arena, fx.xyz, fx.nv, fx.faces, fx.nf, &topo, seed) == 10 &&
                 UvGuard_mesh(seed, fx.faces, fx.nf), &fails);
        AsmContacts *q2 = AsmContacts_new(arena, fx.nv, fx.faces, fx.nf, seed, 1e-7);
        memset(&cs, 0, sizeof cs);
        AsmContacts_measure(q2, seed, NULL, NULL, NULL, &cs);
        fprintf(stderr, "  asm_flatten selftest: %s slit seed contacts %zu -> %zu after opening\n", interior ? "interior" : "boundary", before, cs.pairs);
        af_check("opened slit has no contacts", cs.pairs == 0, &fails);
        float *uv = ARENA_ALLOC(arena, fx.nv * 2 * sizeof(float));
        AsmFlattenStats st;
        static const int recovery[] = {3, 2, 1};
        af_check("slit recovery flattens and certifies",
                 af_flatten_best(arena, fx.xyz, fx.nv, fx.faces, fx.nf, 0.05, 30, recovery, 3, uv, &st) == 0 &&
                 st.seed_mode == 3 && af_chart_certified(fx.xyz, fx.faces, fx.nf, uv) && st.stress_frac < 0.05, &fails);
        af_check("slit chart flattens", AsmFlatten_chart(arena, fx.xyz, fx.nv, fx.faces, fx.nf, 0.05, 30, uv, &st) == 0, &fails);
        /* The conformal solve on PHerc0343 left copies 0.005 vox apart with one inside the other's
         * fan: overlapping material that a fixed step from the copies' own positions did not clear.
         * The exact map with the right-hand copies pushed into the left fan must open cleanly. */
        double *field = ARENA_ALLOC(arena, fx.nv * 2 * sizeof(double));
        for (size_t v = 0; v < fx.nv; v++) { field[2*v] = fx.xyz[3*v+2]; field[2*v+1] = fx.xyz[3*v+1]; }
        size_t first_copy = (size_t)17 * 11;                     /* af_fix_slit appends copies after the grid */
        for (size_t v = first_copy; v < fx.nv; v++) field[2*v] -= 0.003;
        AsmContacts *q3 = AsmContacts_new(arena, fx.nv, fx.faces, fx.nf, field, 1e-7);
        memset(&cs, 0, sizeof cs);
        AsmContacts_measure(q3, field, NULL, NULL, NULL, &cs);
        size_t pushed = cs.pairs;
        af_open_pinches(arena, fx.xyz, fx.nv, fx.faces, fx.nf, &topo, field);
        AsmContacts *q4 = AsmContacts_new(arena, fx.nv, fx.faces, fx.nf, field, 1e-7);
        memset(&cs, 0, sizeof cs);
        AsmContacts_measure(q4, field, NULL, NULL, NULL, &cs);
        af_check("misplaced copies overlap before opening", pushed > 0, &fails);
        af_check("misplaced copies gathered and opened", cs.pairs == 0 && UvGuard_mesh(field, fx.faces, fx.nf), &fails);
        Arena_dispose(&arena);
    }
    {
        /* the mesher's closed slit: DEGENERATE as it is; without its two zero-area faces it is a
         * pinched slit and flattens */
        Arena_T arena = Arena_new(); AfFix fx;
        af_fix_slit(arena, 16, 10, 8, 4, 6, 1, &fx);
        float *uv = ARENA_ALLOC(arena, fx.nv * 2 * sizeof(float));
        AsmFlattenStats st;
        af_check("closed slit is degenerate", AsmFlatten_chart(arena, fx.xyz, fx.nv, fx.faces, fx.nf, 0.05, 30, uv, &st) != 0 &&
                 st.fail_reason == ASM_FLAT_DEGENERATE, &fails);
        int32_t *kf = ARENA_ALLOC(arena, fx.nf * 3 * sizeof(int32_t)), *vm = ARENA_ALLOC(arena, fx.nv * sizeof(int32_t));
        size_t knv = 0, knf = AsmFlatten_drop_degenerate(fx.xyz, fx.nv, fx.faces, fx.nf, kf, vm, &knv);
        af_check("two zero-area faces dropped", knf + 2 == fx.nf && knv == fx.nv &&
                 AsmFlatten_face_degenerate(fx.xyz, fx.faces + 3*(fx.nf - 1)) && !AsmFlatten_face_degenerate(fx.xyz, fx.faces), &fails);
        float *kx = ARENA_ALLOC(arena, knv * 3 * sizeof(float));
        for (size_t i = 0; i < fx.nv; i++) if (vm[i] >= 0) memcpy(kx + 3*(size_t)vm[i], fx.xyz + 3*i, 3 * sizeof(float));
        af_check("chart without zero-area faces flattens", AsmFlatten_chart(arena, kx, knv, kf, knf, 0.05, 30, uv, &st) == 0 &&
                 af_chart_certified(kx, kf, knf, uv), &fails);
        Arena_dispose(&arena);
    }
    {
        /* the certificate refuses a valid map with a face outside the strict band */
        const float xyz[] = {0,0,0, 0,0,1, 0,1,0};
        const int32_t face[] = {0,1,2};
        const float good[] = {0,0, 1,0, 0,1}, squeezed[] = {0,0, 1,0, 0,.7f};
        af_check("certificate accepts isometry", af_chart_certified(xyz, face, 1, good), &fails);
        af_check("certificate refuses a 0.7 face", !af_chart_certified(xyz, face, 1, squeezed), &fails);
    }
    return fails;
}

static int af_projection_recovery_selftest(void)
{
    int fails = 0;
    for (int mild = 0; mild < 2; mild++) for (int rotated = 0; rotated < 2; rotated++) for (int hole = 0; hole < 2; hole++) {
        Arena_T arena = Arena_new(); AfFix patch; AfTopo topo;
        af_fix_grid(arena,12,8,12,8,af_fn_plane,hole,6,4,1.5,&patch);
        double *field = ARENA_ALLOC(arena,2*patch.nv*sizeof(double));
        for (size_t v = 0; v < patch.nv; v++) {
            double x = patch.xyz[3*v+2], y = patch.xyz[3*v+1];
            /* One half already has the right metric; the other is a
             * feasible, compressed map left by a stalled initializer. */
            field[2*v] = x <= 6 ? x : 6+(mild ? .85 : .4)*(x-6); field[2*v+1] = y;
            if (rotated) { patch.xyz[3*v] = (float)(100+x*.6); patch.xyz[3*v+1] = (float)(-200+y); patch.xyz[3*v+2] = (float)(300+x*.8); }
        }
        float *original = ARENA_ALLOC(arena,3*patch.nv*sizeof(float)); memcpy(original,patch.xyz,3*patch.nv*sizeof(float));
        int32_t *original_faces = ARENA_ALLOC(arena,3*patch.nf*sizeof(int32_t)); memcpy(original_faces,patch.faces,3*patch.nf*sizeof(int32_t));
        AfRef *ref = ARENA_ALLOC(arena,patch.nf*sizeof(AfRef)); af_reference(patch.xyz,patch.faces,patch.nf,ref);
        AsmFlattenStats st = {0}; for (size_t f = 0; f < patch.nf; f++) st.area3d += ref[f].area;
        int rc = af_topology(arena,patch.xyz,patch.nv,patch.faces,patch.nf,&topo);
        af_check("projection recovery topology",rc == 0,&fails);
        double before_bad = 0, after_bad = 0; af_sd_energy(ref,patch.faces,patch.nf,field,&before_bad);
        af_check("projection recovery begins outside its applicable metric gate",mild ? before_bad == 0 : before_bad > 0,&fails);
        if (!rc) af_projection_recover(arena,patch.xyz,patch.nv,patch.faces,patch.nf,&topo,field,30,&st);
        for (size_t k = 0; k < 2*patch.nv; k++) field[k] = (float)field[k];
        double energy = af_sd_energy(ref,patch.faces,patch.nf,field,&after_bad);
        af_check("source projection repairs all original faces",after_bad == 0 && isfinite(energy) && energy < 1e-6*st.area3d,&fails);
        af_check("source projection preserves original geometry",!memcmp(original,patch.xyz,3*patch.nv*sizeof(float)) && !memcmp(original_faces,patch.faces,3*patch.nf*sizeof(int32_t)),&fails);
        AsmContacts *query = AsmContacts_new(arena,patch.nv,patch.faces,patch.nf,field,1e-7); AsmContactStats contacts;
        af_check("projection cache has no contacts",query && !AsmContacts_measure(query,field,NULL,NULL,NULL,&contacts) && !contacts.pairs,&fails);
        double *healthy = ARENA_ALLOC(arena,2*patch.nv*sizeof(double)); memcpy(healthy,field,2*patch.nv*sizeof(double));
        if (!rc) af_projection_recover(arena,patch.xyz,patch.nv,patch.faces,patch.nf,&topo,field,30,&st);
        af_check("healthy map is unchanged",!memcmp(healthy,field,2*patch.nv*sizeof(double)),&fails);
        Arena_dispose(&arena);
    }
    {
        /* Projecting a nearly complete cylinder folds its far side. That
         * candidate must not replace a valid, imperfect intrinsic map. */
        Arena_T arena = Arena_new(); AfFix patch; AfTopo topo;
        af_fix_grid(arena,24,8,240,100,af_fn_cyl,0,0,0,0,&patch);
        double *field = ARENA_ALLOC(arena,2*patch.nv*sizeof(double));
        for (int j = 0; j <= 8; j++) for (int i = 0; i <= 24; i++) { size_t v = (size_t)j*25+i; field[2*v] = i*5.; field[2*v+1] = j*12.5; }
        double *saved = ARENA_ALLOC(arena,2*patch.nv*sizeof(double)); memcpy(saved,field,2*patch.nv*sizeof(double));
        AsmFlattenStats st = {0}; int rc = af_topology(arena,patch.xyz,patch.nv,patch.faces,patch.nf,&topo);
        af_check("folded projection control topology",rc == 0,&fails);
        if (!rc) af_projection_recover(arena,patch.xyz,patch.nv,patch.faces,patch.nf,&topo,field,30,&st);
        af_check("folded source projection refused",!memcmp(saved,field,2*patch.nv*sizeof(double)),&fails);
        Arena_dispose(&arena);
    }
    return fails;
}

static int af_band_selftest(void)
{
    int fails = 0;
    {
        Arena_T arena = Arena_new();
        const float xyz[9] = {0,0,0, 1,0,0, 0,1,0}; const int32_t faces[3] = {0,1,2};
        AfRef ref; af_reference(xyz,faces,1,&ref);
        size_t active[1] = {0}, vertices[3] = {0,1,2}; int32_t index[3] = {0,1,2};
        double seed[6] = {0,0, 1.4,.25, .2,.6}, scale[3] = {1,2,3}, trial[6];
        double x[6] = {0}, gradient[6], temporary[6];
        AfBand band = {&ref,faces,active,vertices,index,1,3,seed,scale,trial};
        double value = af_band_loss(&band,x,gradient); int ok = value > 0 && isfinite(value);
        for (int i = 0; i < 6; i++) {
            x[i] = 1e-6; double plus = af_band_loss(&band,x,temporary);
            x[i] = -1e-6; double minus = af_band_loss(&band,x,temporary); x[i] = 0;
            double numerical = (plus-minus)/2e-6;
            if (fabs(numerical-gradient[i]) > 1e-6*(1+fabs(numerical))) ok = 0;
        }
        af_check("local original-metric gradient",ok,&fails);
        Arena_dispose(&arena);
    }
    for (int shifted = 0; shifted < 2; shifted++) for (int hole = 0; hole < 2; hole++) {
        Arena_T arena = Arena_new(); AfFix patch;
        af_fix_grid(arena,12,8,12,8,af_fn_plane,hole,6,4,1.5,&patch);
        double *field = ARENA_ALLOC(arena,2*patch.nv*sizeof *field);
        for (size_t v = 0; v < patch.nv; v++) {
            double x = patch.xyz[3*v+2], y = patch.xyz[3*v+1];
            field[2*v] = (shifted ? 1e6 : 0)+(x <= 6 ? x : 6+.4*(x-6));
            field[2*v+1] = (shifted ? -2e6 : 0)+y;
        }
        float *original = ARENA_ALLOC(arena,3*patch.nv*sizeof *original); memcpy(original,patch.xyz,3*patch.nv*sizeof *original);
        int32_t *original_faces = ARENA_ALLOC(arena,3*patch.nf*sizeof *original_faces); memcpy(original_faces,patch.faces,3*patch.nf*sizeof *original_faces);
        AfRef *ref = ARENA_ALLOC(arena,patch.nf*sizeof *ref); af_reference(patch.xyz,patch.faces,patch.nf,ref);
        double bad = 0; af_sd_energy(ref,patch.faces,patch.nf,field,&bad);
        af_check("local metric control starts damaged",bad > 0,&fails);
        double *prior = ARENA_ALLOC(arena,2*patch.nv*sizeof *prior); memcpy(prior,field,2*patch.nv*sizeof *prior);
        AsmFlattenStats st = {0}; af_band_recover(arena,patch.xyz,patch.nv,patch.faces,patch.nf,field,30,&st);
        double *encoded = ARENA_ALLOC(arena,2*patch.nv*sizeof *encoded);
        int certified = af_band_certificate(arena,ref,patch.faces,patch.nf,patch.nv,field,encoded);
        /* At a million-voxel gauge the float lattice is too coarse for the
         * nominal proposal margin to guarantee recovery. Refusal must keep
         * the full old map; it must not leak a partially repaired iterate. */
        af_check("local recovery certifies the map or preserves coarse-gauge input",certified || (shifted && !memcmp(prior,field,2*patch.nv*sizeof *prior)),&fails);
        af_check("local recovery preserves original geometry",!memcmp(original,patch.xyz,3*patch.nv*sizeof *original) && !memcmp(original_faces,patch.faces,3*patch.nf*sizeof *original_faces),&fails);
        memcpy(encoded,field,2*patch.nv*sizeof *encoded);
        if (certified) af_band_recover(arena,patch.xyz,patch.nv,patch.faces,patch.nf,field,30,&st);
        af_check("healthy local map remains exact",!memcmp(encoded,field,2*patch.nv*sizeof *encoded),&fails);
        Arena_dispose(&arena);
    }
    {
        /* A very sharp cone fan cannot put its apex inside this short rim
         * while keeping every source triangle in the band. Refuse the whole
         * proposal, including every intermediate local optimizer step. */
        Arena_T arena = Arena_new(); AsmFlattenStats st = {0};
        const float xyz[12] = {0,0,10, 1,0,0, -.5f,.8660254f,0, -.5f,-.8660254f,0};
        const int32_t faces[9] = {0,1,2, 0,2,3, 0,3,1};
        double field[8] = {0,0, 1,0, -.5,.8660254, -.5,-.8660254}, before[8]; memcpy(before,field,sizeof before);
        af_band_recover(arena,xyz,4,faces,3,field,1,&st);
        af_check("infeasible local certificate preserves prior map",!memcmp(before,field,sizeof before),&fails);
        Arena_dispose(&arena);
    }
    {
        /* A map outside the band on all 768 faces is past the measured recovery budget
         * (ASM_FLATTEN_BAND_RECOVER_MAX_BAD): no step is taken and the map stays exact. */
        Arena_T arena = Arena_new(); AfFix patch;
        af_fix_grid(arena,24,16,24,16,af_fn_plane,0,0,0,0,&patch);
        double *field = ARENA_ALLOC(arena,2*patch.nv*sizeof *field), *before = ARENA_ALLOC(arena,2*patch.nv*sizeof *before);
        for (size_t v = 0; v < patch.nv; v++) { field[2*v] = .4*patch.xyz[3*v+2]; field[2*v+1] = .4*patch.xyz[3*v+1]; }
        memcpy(before,field,2*patch.nv*sizeof *before);
        AsmFlattenStats st = {0}; af_band_recover(arena,patch.xyz,patch.nv,patch.faces,patch.nf,field,30,&st);
        af_check("hopeless local recovery is not attempted",patch.nf > ASM_FLATTEN_BAND_RECOVER_MAX_BAD &&
                 st.metric_refine_iters == 0 && !memcmp(before,field,2*patch.nv*sizeof *before),&fails);
        Arena_dispose(&arena);
    }
    af_check("capped conformal seed of a disk is the conformal seed",
             af_seed_effective_mode(1,1) == 3 && af_seed_effective_mode(1,2) == 1 && af_seed_effective_mode(2,1) == 2 &&
             af_seed_effective_mode(3,1) == 3 && af_seed_effective_mode(0,1) == 0,&fails);
    return fails;
}

static int af_solver_selftest(void)
{
    /* The rescue's direct solves factor a symmetric matrix once and agree with the historical PCG;
     * a non-symmetric matrix, and every solve outside the rescue, stays on PCG. */
    enum { N = 60 };
    int fails = 0;
    Arena_T arena = Arena_new();
    int32_t *off = ARENA_ALLOC(arena, (N+1)*sizeof *off), *col = ARENA_ALLOC(arena, 2*N*sizeof *col);
    double *val = ARENA_ALLOC(arena, 2*N*sizeof *val), *diag = ARENA_ALLOC(arena, N*sizeof *diag);
    int32_t k = 0;
    for (int i = 0; i < N; i++) {
        off[i] = k; diag[i] = 2.0 + 1e-3*i;
        if (i > 0) { col[k] = i-1; val[k] = -1.0; k++; }
        if (i < N-1) { col[k] = i+1; val[k] = -1.0; k++; }
    }
    off[N] = k;
    AfMat m = { N, off, col, val, diag };
    double b[N], x0[N], x1[N], r[N], z[N], p[N], q[N];
    for (int i = 0; i < N; i++) { b[i] = sin(.3*i) + .1*i; x0[i] = 0.0; x1[i] = 0.0; }
    int pcg = af_pcg(&m, b, x0, r, z, p, q, 1e-13, 10000);
    int prev = af_direct;
    af_direct = 1;
    AfSolver s; af_solver_init(arena, &s, &m);
    int factored = s.f != NULL;
    int direct = af_solver_solve(&s, b, x1, r, z, p, q, 1e-10, 10000);
    af_solver_free(&s);
    double err = 0.0;
    for (int i = 0; i < N; i++) err = fmax(err, fabs(x0[i] - x1[i]));
    af_check("rescue direct solve factors the SPD matrix and matches PCG",
             pcg == 0 && factored && direct == 0 && err < 1e-8, &fails);
    val[0] = -.5;                                  /* (0,1) no longer equals (1,0) */
    af_solver_init(arena, &s, &m);
    af_check("non-symmetric matrix stays on PCG", s.f == NULL, &fails);
    af_solver_free(&s);
    af_direct = 0;
    val[0] = -1.0;
    af_solver_init(arena, &s, &m);
    af_check("outside the rescue nothing is factored", s.f == NULL, &fails);
    af_solver_free(&s);
    af_direct = prev;
    Arena_dispose(&arena);
    return fails;
}

static int af_projection_area_selftest(void)
{
    /* Removing a four-percent severe failure cannot create a new chart
     * quality failure by damaging two percent of previously good area. */
    const float xyz[]={0,0,0,188,0,0,0,1,0, 200,0,0,204,0,0,200,1,0, 220,0,0,228,0,0,220,1,0};
    const int32_t faces[]={0,1,2,3,4,5,6,7,8};
    double uv[]={0,0,188,0,0,1, 200,0,204,0,200,1, 220,0,231.2,0,220,1};
    double candidate[]={0,0,188,0,0,1, 200,0,204.44,0,200,1, 220,0,228.88,0,220,1};
    double before[18],encoded[18];memcpy(before,uv,sizeof uv);
    Arena_T arena=Arena_new();AfRef ref[3];af_reference(xyz,faces,3,ref);
    double bad=0,energy=af_sd_energy(ref,faces,3,uv,&bad);int fails=0;
    af_check("area-quality control begins with four percent strict failures",bad==4,&fails);
    int accepted=af_recovery_candidate(arena,ref,faces,3,9,candidate,encoded,uv,&energy,&bad);
    af_check("projection preserves the narrower chart-area budget",!accepted && !memcmp(before,uv,sizeof uv),&fails);
    Arena_dispose(&arena);return fails;
}
int AsmFlatten_selftest(void)
{
    int fails = af_projection_recovery_selftest()+af_band_selftest()+af_solver_selftest()+af_projection_area_selftest()+af_pinch_selftest()+
                af_crumple_selftest();
    Arena_T arena = Arena_new();
    AfFix fx;
    AsmFlattenStats st;
    float *uv;

    {
        /* Two damaged corner faces occupy less than one percent of this
         * otherwise isometric sheet. The area gate used to skip them. */
        AfFix patch; af_fix_grid(arena,20,20,20,20,af_fn_plane,0,0,0,0,&patch);
        AfRef *ref = ARENA_ALLOC(arena,patch.nf*sizeof(AfRef));
        af_reference(patch.xyz,patch.faces,patch.nf,ref);
        double *field = ARENA_ALLOC(arena,2*patch.nv*sizeof(double));
        double area = 0;
        for (size_t f = 0; f < patch.nf; f++) area += ref[f].area;
        for (size_t v = 0; v < patch.nv; v++) { field[2*v] = patch.xyz[3*v+2]; field[2*v+1] = patch.xyz[3*v+1]; }
        field[2*(patch.nv-1)] += .5;
        double before_bad = 0, after_bad = 0;
        double before = af_sd_energy(ref,patch.faces,patch.nf,field,&before_bad);
        af_check("small damaged patch is below old area gate",before_bad > 0 && before_bad < .01*area,&fails);
        AsmFlattenStats refined = {0};
        af_metric_continue(arena,patch.xyz,patch.nv,patch.faces,patch.nf,field,area,30,&refined);
        for (size_t k = 0; k < 2*patch.nv; k++) field[k] = (float)field[k];
        double after = af_sd_energy(ref,patch.faces,patch.nf,field,&after_bad);
        af_check("small damaged patch repaired",refined.metric_refine_iters > 0 && after < before && after_bad == 0,&fails);
        af_check("small patch stored field remains oriented",UvGuard_mesh(field,patch.faces,patch.nf),&fails);
    }

    /* A strip wrapped through more than a full turn has positive faces but
     * overlaps itself. Boundary feasibility must be an independent guard. */
    {
        AfFix strip; AfTopo topo;
        af_fix_grid(arena,20,1,1,1,af_fn_plane,0,0,0,0,&strip);
        af_check("boundary control topology",af_topology(arena,strip.xyz,strip.nv,strip.faces,strip.nf,&topo) == 0,&fails);
        double *field = ARENA_ALLOC(arena,2*strip.nv*sizeof(double));
        for (size_t v = 0; v < strip.nv; v++) { field[2*v] = strip.xyz[3*v+2]; field[2*v+1] = strip.xyz[3*v+1]; }
        AfBoundaryQuery boundary = af_boundary_query(arena,&topo,strip.faces,strip.nf,field), *boxes = &boundary;
        uint8_t *blocked = ARENA_CALLOC(arena,strip.nv,1);
        for (int wrap = 0; wrap < 2; wrap++) {
            for (size_t v = 0; v < strip.nv; v++) {
                double angle = -(wrap ? 2.5 : .75)*AF_PI*strip.xyz[3*v+2], radius = 1+strip.xyz[3*v+1];
                field[2*v] = radius*cos(angle); field[2*v+1] = radius*sin(angle);
            }
            af_check("wrapped strip faces positive",UvGuard_mesh(field,strip.faces,strip.nf),&fails);
            af_check("wrapped boundary independently checked",af_boundary_clear(&topo,field,boxes,blocked) == !wrap,&fails);
        }
        size_t count = 0; for (size_t v = 0; v < strip.nv; v++) count += blocked[v] != 0;
        af_check("boundary contacts mark incident vertices",count >= 4 && count < strip.nv,&fails);
    }
    {
        /* Two boundary banks can touch without overlapping material. A
         * segment-intersection guard used to freeze such valid proposals. */
        int32_t faces[] = {0,1,2,0,2,3,4,5,6,4,6,7};
        double field[] = {0,0,1,0,1,1,0,1,1,0,2,0,2,1,1,1};
        uint8_t flags[] = {1,1,1,1,1,1,1,1}; AfTopo topo = {0};
        topo.nv = 8; topo.on_boundary = flags;
        AfBoundaryQuery q = af_boundary_query(arena,&topo,faces,4,field);
        af_check("touching boundary banks allowed",af_boundary_clear(&topo,field,&q,NULL),&fails);
        for (size_t v = 4; v < 8; v++) field[2*v] -= .25;
        af_check("overlapping boundary banks refused",!af_boundary_clear(&topo,field,&q,NULL),&fails);
    }

    /* Positive endpoints alone miss an inversion and re-inversion along
     * the proposed step. Also exercise reflection, collapse and translation. */
    {
        int32_t face[] = {0,1,2};
        double a[] = {0,0, 1,0, 0,1}, b[] = {0,0, -1,0, 0,-1};
        af_check("two-root endpoint is oriented",UvGuard_mesh(b,face,1),&fails);
        af_check("two-root step refused",!UvGuard_interval(a,b,face,1),&fails);
        af_check("bounded step before fold accepted",UvGuard_interval(a,b,face,.25),&fails);
        b[5] = 1;
        af_check("reflected endpoint refused",!UvGuard_interval(a,b,face,1),&fails);
        b[2] = 0;
        af_check("collapsed endpoint refused",!UvGuard_interval(a,b,face,1),&fails);
        for (size_t i = 0; i < 6; i++) b[i] = a[i]+1024;
        af_check("rigid translation accepted",UvGuard_interval(a,b,face,1),&fails);
    }

    /* Thin obtuse faces still define a positive metric. Independent cot
     * clipping can make the local energy negative; subtracting squared
     * lengths can also erase the height of a nondegenerate thin face. */
    for (int mirror = 0; mirror < 2; mirror++) for (int thin = 0; thin < 2; thin++) {
        float height = thin ? 1e-8f : .01f;
        float xyz[] = {0,0,0, 1,0,0, -1,height,0}; int32_t face[] = {0,1,2};
        if (mirror) xyz[7] = -xyz[7];
        AfRef reference; af_reference(xyz,face,1,&reference);
        double energy = reference.cot[0]+4*reference.cot[1]+reference.cot[2];
        int ok = energy >= -1e-12 && fabs(reference.area-.5*height) <= 1e-6*height;
        if (!ok) { fprintf(stderr,"  asm_flatten selftest FAIL: obtuse thin face mirror %d thin %d, energy %.9g area %.9g expected %.9g\n",mirror,thin,energy,reference.area,.5*height); fails++; }
    }
    {
        /* CG must reject a negative-curvature direction even when solving
         * the resulting indefinite linear equation is numerically easy. */
        int32_t off[] = {0,0}, col[] = {0}; double diag[] = {-1}, val[] = {0};
        AfMat matrix = {1,off,col,val,diag}; double b[] = {1}, x[] = {0}, r[1], z[1], p[1], q[1];
        af_check("negative curvature refused",af_pcg(&matrix,b,x,r,z,p,q,1e-9,4) != 0,&fails);
    }

    /* 1. planar grid: exact isometry */
    af_fix_grid(arena, 12, 8, 120.0, 80.0, af_fn_plane, 0, 0, 0, 0, &fx);
    uv = ARENA_ALLOC(arena, fx.nv * 2 * sizeof(float));
    af_check("plane flattens", AsmFlatten_chart(arena, fx.xyz, fx.nv, fx.faces, fx.nf, 0.05, 50, uv, &st) == 0, &fails);
    af_check("plane sigma ~1", st.ok && st.sigma_lo > 0.999 && st.sigma_hi < 1.001, &fails);
    af_check("plane no flips", st.n_flipped == 0, &fails);
    af_check("plane area", fabs(st.area_uv - st.area3d) < 1e-3 * st.area3d, &fails);
    af_check("plane one loop", st.n_boundary_loops == 1, &fails);

    /* 2. half cylinder: developable */
    af_fix_grid(arena, 24, 10, 125.0, 100.0, af_fn_cyl, 0, 0, 0, 0, &fx);
    uv = ARENA_ALLOC(arena, fx.nv * 2 * sizeof(float));
    af_check("cylinder flattens", AsmFlatten_chart(arena, fx.xyz, fx.nv, fx.faces, fx.nf, 0.05, 80, uv, &st) == 0, &fails);
    af_check("cylinder sigma", st.ok && st.sigma_lo > 0.99 && st.sigma_hi < 1.01, &fails);
    af_check("cylinder no flips", st.n_flipped == 0, &fails);
    af_check("cylinder stress 0", st.stress_frac == 0.0, &fails);

    /* 3. cone patch: developable */
    af_fix_grid(arena, 24, 12, 100.0, 80.0, af_fn_cone, 0, 0, 0, 0, &fx);
    uv = ARENA_ALLOC(arena, fx.nv * 2 * sizeof(float));
    af_check("cone flattens", AsmFlatten_chart(arena, fx.xyz, fx.nv, fx.faces, fx.nf, 0.05, 80, uv, &st) == 0, &fails);
    af_check("cone sigma", st.ok && st.sigma_lo > 0.98 && st.sigma_hi < 1.02, &fails);
    af_check("cone no flips", st.n_flipped == 0, &fails);

    /* 4. holed rectangle on a cylinder: two loops, one capped */
    af_fix_grid(arena, 24, 12, 125.0, 120.0, af_fn_cyl, 1, 12.0, 6.0, 2.6, &fx);
    uv = ARENA_ALLOC(arena, fx.nv * 2 * sizeof(float));
    af_check("holed flattens", AsmFlatten_chart(arena, fx.xyz, fx.nv, fx.faces, fx.nf, 0.05, 80, uv, &st) == 0, &fails);
    af_check("holed two loops", st.n_boundary_loops == 2 && st.n_holes_capped == 1, &fails);
    af_check("holed sigma", st.ok && st.sigma_lo > 0.98 && st.sigma_hi < 1.02, &fails);
    af_check("holed no flips", st.n_flipped == 0, &fails);

    /* 5. spherical cap: not developable, must be stressed at band 0.02 */
    af_fix_grid(arena, 20, 20, 100.0, 100.0, af_fn_sphere, 0, 0, 0, 0, &fx);
    uv = ARENA_ALLOC(arena, fx.nv * 2 * sizeof(float));
    af_check("sphere flattens", AsmFlatten_chart(arena, fx.xyz, fx.nv, fx.faces, fx.nf, 0.02, 80, uv, &st) == 0, &fails);
    af_check("sphere is stressed", st.ok && st.stress_frac > 0.2, &fails);
    af_check("sphere no flips", st.n_flipped == 0, &fails);
    /* a rescue returns only certified maps */
    af_check("certified rescue is certified or refused",
             AsmFlatten_chart_certified(arena, fx.xyz, fx.nv, fx.faces, fx.nf, 0.02, 80, uv, &st) == 0 ?
             af_chart_certified(fx.xyz, fx.faces, fx.nf, uv) : st.fail_reason != ASM_FLAT_OK, &fails);

    /* 6. closed tetrahedron: no boundary -> refused */
    {
        float tet[12] = { 0,0,0,  0,0,10,  0,10,0,  10,0,0 };
        int32_t tf[12] = { 0,2,1,  0,1,3,  0,3,2,  1,2,3 };
        float tuv[8];
        af_check("tetra refused", AsmFlatten_chart(arena, tet, 4, tf, 4, 0.05, 10, tuv, &st) != 0
                 && st.fail_reason == ASM_FLAT_NO_BOUNDARY, &fails);
    }
    /* 7. determinism: same input twice gives identical uv */
    {
        af_fix_grid(arena, 10, 6, 60.0, 40.0, af_fn_cyl, 0, 0, 0, 0, &fx);
        float *a = ARENA_ALLOC(arena, fx.nv * 2 * sizeof(float));
        float *b = ARENA_ALLOC(arena, fx.nv * 2 * sizeof(float));
        AsmFlattenStats s2;
        AsmFlatten_chart(arena, fx.xyz, fx.nv, fx.faces, fx.nf, 0.05, 50, a, &st);
        AsmFlatten_chart(arena, fx.xyz, fx.nv, fx.faces, fx.nf, 0.05, 50, b, &s2);
        af_check("deterministic", memcmp(a, b, fx.nv * 2 * sizeof(float)) == 0, &fails);
    }
    /* 8. a crushed-scroll turn: the conformal seed recovers the isometric
     * flattening at the production iteration cap; the Tutte-only path is
     * reported for the record */
    {
        af_fix_grid(arena, 72, 16, 120.0 + 6.0 * AF_PI, 40.0, af_fn_uturn, 0, 0, 0, 0, &fx);
        float *a = ARENA_ALLOC(arena, fx.nv * 2 * sizeof(float));
        float *b = ARENA_ALLOC(arena, fx.nv * 2 * sizeof(float));
        AsmFlattenStats base;
        int rb = af_flatten_chart(arena, fx.xyz, fx.nv, fx.faces, fx.nf, 0.05, 30, 0, 0, a, &base);
        int rc = af_flatten_chart(arena, fx.xyz, fx.nv, fx.faces, fx.nf, 0.05, 30, 1, 0, b, &st);
        fprintf(stderr, "  asm_flatten selftest: U-turn Tutte-only rc %d stress %.4f sigma [%.4f, %.4f]; "
                        "with conformal seed rc %d stress %.4f sigma [%.4f, %.4f]\n",
                rb, base.stress_frac, base.sigma_lo, base.sigma_hi, rc, st.stress_frac, st.sigma_lo, st.sigma_hi);
        af_check("U-turn flattens", rc == 0 && st.ok, &fails);
        af_check("U-turn isometric", st.stress_frac == 0.0 && st.sigma_lo > 0.99 && st.sigma_hi < 1.01, &fails);
        af_check("U-turn no flips", st.n_flipped == 0, &fails);
    }
    /* 9. a chart already inside the band is bitwise unchanged by the seed */
    {
        af_fix_grid(arena, 24, 10, 125.0, 100.0, af_fn_cyl, 0, 0, 0, 0, &fx);
        float *a = ARENA_ALLOC(arena, fx.nv * 2 * sizeof(float));
        float *b = ARENA_ALLOC(arena, fx.nv * 2 * sizeof(float));
        AsmFlattenStats s2;
        af_flatten_chart(arena, fx.xyz, fx.nv, fx.faces, fx.nf, 0.05, 80, 0, 0, a, &s2);
        af_flatten_chart(arena, fx.xyz, fx.nv, fx.faces, fx.nf, 0.05, 80, 1, 0, b, &st);
        af_check("conformal seed inert inside the band",
                 s2.stress_frac == 0.0 && memcmp(a, b, fx.nv * 2 * sizeof(float)) == 0, &fails);
    }
    Arena_dispose(&arena);
    if (fails == 0) fprintf(stderr, "  asm_flatten selftest: all passed\n");
    return fails;
}
