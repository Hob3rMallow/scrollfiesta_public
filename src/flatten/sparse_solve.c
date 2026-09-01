#include "sparse_solve.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "taucs.h"

/* taucs.h does not prototype these high-level entry points; declare them. */
extern taucs_ccs_matrix *taucs_ccs_create(int m, int n, int nnz, int flags);
extern void              taucs_ccs_free(taucs_ccs_matrix *m);
extern int  taucs_linsolve(taucs_ccs_matrix *A, void **F, int nrhs,
                           void *X, void *B, char *options[], void *opt_arg[]);

/* qsort comparator for (row,val) pairs packed as parallel arrays is awkward;
 * sort an index permutation of a column segment by row instead. */
typedef struct { int row; double val; } RV;
static int cmp_rv(const void *a, const void *b)
{
    int ra = ((const RV *)a)->row, rb = ((const RV *)b)->row;
    return (ra < rb) ? -1 : (ra > rb ? 1 : 0);
}

/* Build a TAUCS lower-symmetric CCS matrix from COO triplets (row >= col,
 * duplicates summed).  Returns NULL on invalid input or allocation failure. */
static taucs_ccs_matrix *sp_build_ccs(int n, int nt,
                                      const int *rows, const int *cols,
                                      const double *vals)
{
    /* --- Bucket triplets by column. --- */
    int *colcnt = (int *)calloc((size_t)n + 1, sizeof(int));
    if (!colcnt) return NULL;
    for (int t = 0; t < nt; t++) {
        int c = cols[t];
        if (c < 0 || c >= n || rows[t] < c || rows[t] >= n) { free(colcnt); return NULL; }
        colcnt[c]++;
    }
    int *cstart = (int *)malloc(((size_t)n + 1) * sizeof(int));
    if (!cstart) { free(colcnt); return NULL; }
    cstart[0] = 0;
    for (int c = 0; c < n; c++) cstart[c + 1] = cstart[c] + colcnt[c];

    RV *buf = (RV *)malloc((size_t)(nt > 0 ? nt : 1) * sizeof(RV));
    int *cur = (int *)malloc((size_t)n * sizeof(int));
    if (!buf || !cur) { free(colcnt); free(cstart); free(buf); free(cur); return NULL; }
    for (int c = 0; c < n; c++) cur[c] = cstart[c];
    for (int t = 0; t < nt; t++) {
        int c = cols[t];
        buf[cur[c]].row = rows[t];
        buf[cur[c]].val = vals[t];
        cur[c]++;
    }

    /* --- Per column: sort by row, merge duplicates -> CCS. --- */
    taucs_ccs_matrix *A = taucs_ccs_create(n, n, (nt > 0 ? nt : 1),
                              TAUCS_DOUBLE | TAUCS_SYMMETRIC | TAUCS_LOWER);
    if (!A) { free(colcnt); free(cstart); free(buf); free(cur); return NULL; }

    int k = 0;
    for (int c = 0; c < n; c++) {
        A->colptr[c] = k;
        int s = cstart[c], e = cstart[c + 1];
        if (e > s) {
            qsort(buf + s, (size_t)(e - s), sizeof(RV), cmp_rv);
            int i = s;
            while (i < e) {
                int r = buf[i].row;
                double acc = 0.0;
                while (i < e && buf[i].row == r) { acc += buf[i].val; i++; }
                A->rowind[k] = r;
                A->values.d[k] = acc;
                k++;
            }
        }
    }
    A->colptr[n] = k;

    free(colcnt); free(cstart); free(buf); free(cur);
    return A;
}

int Sparse_solve_sym(int n, int nt,
                     const int *rows, const int *cols, const double *vals,
                     const double *b, double *x, SparseMode mode)
{
    if (n <= 0 || nt < 0) return -1;
    if (!rows || !cols || !vals || !b || !x) return -1;

    /* Small active sets legitimately reduce to one column.  The exact scalar
     * solve is simpler and avoids sending a degenerate tree through the
     * multifrontal machinery. */
    if (n == 1) {
        double a = 0.0;
        for (int t = 0; t < nt; t++) {
            if (rows[t] != 0 || cols[t] != 0) return -1;
            a += vals[t];
        }
        if (a == 0.0 || (mode == SPARSE_SPD && a < 0.0)) return -1;
        x[0] = b[0] / a;
        return 0;
    }

    taucs_ccs_matrix *A = sp_build_ccs(n, nt, rows, cols, vals);
    if (!A) return -1;

    /* --- Solve. --- */
    char *opt_spd[]   = { "taucs.factor.LLT=true", "taucs.factor.mf=true", NULL };
    /* TAUCS's legacy f2c GENMMD path is not reliable on very small graphs
     * (the 13-vertex welded-scroll disk reaches it with only a handful of
     * interior unknowns and can fault under /O2).  Ordering cannot materially
     * reduce fill at this scale, so use the exact identity permutation. */
    char *opt_spd_small[] = { "taucs.factor.LLT=true", "taucs.factor.mf=true",
                              "taucs.factor.ordering=identity", NULL };
    char *opt_indef[] = { "taucs.solve.minres=true",
                          "taucs.solve.maxits=10000",
                          "taucs.solve.convergetol=1e-11", NULL };
    void *F = NULL;
    char **solve_opts = (mode == SPARSE_SPD)
                      ? ((n <= 64) ? opt_spd_small : opt_spd)
                      : opt_indef;
    int rc = taucs_linsolve(A, &F, 1, x, (void *)b, solve_opts, NULL);
    if (F) { void *Fp = F; taucs_linsolve(NULL, &Fp, 0, NULL, NULL, NULL, NULL); }
    taucs_ccs_free(A);
    return (rc == TAUCS_SUCCESS) ? 0 : -1;
}

/* ---- factor-once / solve-many ------------------------------------------- */

struct SparseFactor {
    int n;
    taucs_ccs_matrix *A;   /* kept for taucs_linsolve's solve pass */
    void *F;               /* taucs factorization (LLT supernodal) */
    double scalar;         /* n==1 exact path: A==NULL, F==NULL */
};

int Sparse_factor_spd(int n, int nt,
                      const int *rows, const int *cols, const double *vals,
                      SparseFactor_T *out)
{
    if (out == NULL) return -1;
    *out = NULL;
    if (n <= 0 || nt < 0 || !rows || !cols || !vals) return -1;

    struct SparseFactor *f =
        (struct SparseFactor *)calloc(1, sizeof *f);
    if (!f) return -1;
    f->n = n;

    if (n == 1) {
        double a = 0.0;
        for (int t = 0; t < nt; t++) {
            if (rows[t] != 0 || cols[t] != 0) { free(f); return -1; }
            a += vals[t];
        }
        if (a <= 0.0) { free(f); return -1; }
        f->scalar = a;
        *out = f;
        return 0;
    }

    f->A = sp_build_ccs(n, nt, rows, cols, vals);
    if (!f->A) {
        fprintf(stderr,
                "Sparse_factor_spd: CCS build failed (n=%d nt=%d; likely "
                "allocation)\n", n, nt);
        free(f);
        return -1;
    }

    char *opt_factor[]       = { "taucs.factor.LLT=true",
                                 "taucs.factor.mf=true", NULL };
    char *opt_factor_small[] = { "taucs.factor.LLT=true",
                                 "taucs.factor.mf=true",
                                 "taucs.factor.ordering=identity", NULL };
    /* Factor through a one-shot solve against a throwaway rhs, KEEPING the
     * factorization handle: this is the exact call shape the one-shot
     * solver uses (proven at 16-25M unknowns); the factor-only nrhs=0 path
     * failed at large connected systems where this succeeds. */
    double *dummy_b = (double *)calloc((size_t)n, sizeof(double));
    double *dummy_x = (double *)calloc((size_t)n, sizeof(double));
    int rc = -1;
    if (dummy_b == NULL || dummy_x == NULL) {
        fprintf(stderr, "Sparse_factor_spd: scratch allocation failed "
                "(n=%d)\n", n);
        free(dummy_b);
        free(dummy_x);
        taucs_ccs_free(f->A);
        free(f);
        return -1;
    }
    rc = taucs_linsolve(f->A, &f->F, 1, dummy_x, dummy_b,
                        (n <= 64) ? opt_factor_small : opt_factor, NULL);
    free(dummy_b);
    free(dummy_x);
    if (rc != TAUCS_SUCCESS || f->F == NULL) {
        fprintf(stderr,
                "Sparse_factor_spd: taucs_linsolve factor rc=%d F=%s "
                "(n=%d nt=%d)\n", rc, f->F != NULL ? "set" : "null", n, nt);
        if (f->F) { void *Fp = f->F; taucs_linsolve(NULL, &Fp, 0, NULL, NULL, NULL, NULL); }
        taucs_ccs_free(f->A);
        free(f);
        return -1;
    }
    *out = f;
    return 0;
}

int Sparse_factor_solve_multi(SparseFactor_T f, const double *B, double *X,
                              int nrhs)
{
    if (!f || !B || !X || nrhs <= 0) return -1;
    if (f->A == NULL) {           /* scalar path */
        for (int j = 0; j < nrhs; j++) X[j] = B[j] / f->scalar;
        return 0;
    }
    char *opt_solve[] = { "taucs.factor=false", NULL };
    int rc = taucs_linsolve(f->A, &f->F, nrhs, X, (void *)B, opt_solve, NULL);
    return (rc == TAUCS_SUCCESS) ? 0 : -1;
}

int Sparse_factor_solve(SparseFactor_T f, const double *b, double *x)
{
    return Sparse_factor_solve_multi(f, b, x, 1);
}

void Sparse_factor_free(SparseFactor_T *pf)
{
    if (!pf || !*pf) return;
    struct SparseFactor *f = *pf;
    if (f->F) { void *Fp = f->F; taucs_linsolve(NULL, &Fp, 0, NULL, NULL, NULL, NULL); }
    if (f->A) taucs_ccs_free(f->A);
    free(f);
    *pf = NULL;
}

/* ---- selftest ------------------------------------------------------------ */

int Sparse_selftest(void)
{
    /* 6-node SPD: path-graph Laplacian + 0.5 I, two right-hand sides.  The
     * factored solves must agree with the one-shot solver to solver noise. */
    enum { N = 6 };
    int rows[2 * N], cols[2 * N];
    double vals[2 * N];
    int nt = 0;
    for (int i = 0; i < N; i++) {
        double d = 0.5 + (i > 0 ? 1.0 : 0.0) + (i + 1 < N ? 1.0 : 0.0);
        rows[nt] = i; cols[nt] = i; vals[nt] = d; nt++;
        if (i > 0) { rows[nt] = i; cols[nt] = i - 1; vals[nt] = -1.0; nt++; }
    }
    double b1[N], b2[N], x_ref[N], x_fac[N];
    for (int i = 0; i < N; i++) {
        b1[i] = 1.0 + (double)i;
        b2[i] = ((i % 2) == 0) ? 3.0 : -2.0;
    }

    SparseFactor_T f = NULL;
    if (Sparse_factor_spd(N, nt, rows, cols, vals, &f) != 0) {
        fprintf(stderr, "sparse selftest: factor failed\n");
        return -1;
    }
    const double *bs[2] = { b1, b2 };
    for (int k = 0; k < 2; k++) {
        memset(x_ref, 0, sizeof x_ref);
        memset(x_fac, 0, sizeof x_fac);
        if (Sparse_solve_sym(N, nt, rows, cols, vals, bs[k], x_ref,
                             SPARSE_SPD) != 0 ||
            Sparse_factor_solve(f, bs[k], x_fac) != 0) {
            fprintf(stderr, "sparse selftest: solve %d failed\n", k);
            Sparse_factor_free(&f);
            return -1;
        }
        for (int i = 0; i < N; i++) {
            if (fabs(x_ref[i] - x_fac[i]) > 1e-9) {
                fprintf(stderr,
                        "sparse selftest: rhs %d mismatch at %d: %g vs %g\n",
                        k, i, x_ref[i], x_fac[i]);
                Sparse_factor_free(&f);
                return -1;
            }
        }
    }
    Sparse_factor_free(&f);

    /* n==1 scalar path */
    {
        int r0 = 0, c0 = 0;
        double v0 = 4.0, b0 = 10.0, x0 = 0.0;
        SparseFactor_T s = NULL;
        if (Sparse_factor_spd(1, 1, &r0, &c0, &v0, &s) != 0 ||
            Sparse_factor_solve(s, &b0, &x0) != 0 ||
            fabs(x0 - 2.5) > 1e-12) {
            fprintf(stderr, "sparse selftest: scalar path failed\n");
            Sparse_factor_free(&s);
            return -1;
        }
        Sparse_factor_free(&s);
    }
    fprintf(stderr, "sparse selftest: factor reuse OK (2 rhs + scalar)\n");
    return 0;
}
