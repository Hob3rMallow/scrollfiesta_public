#include "sparse_solve.h"

#include <float.h>
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

/* TAUCS' multifrontal path is not reliable for degenerate tiny elimination
 * trees on Win64.  Solve those systems directly; production-sized systems
 * still use TAUCS below.  Partial pivoting also covers the small indefinite
 * systems requested by ABF++. */
static int sparse_solve_small_dense(int n, int nt,
                                    const int *rows, const int *cols,
                                    const double *vals, const double *b,
                                    double *x)
{
    size_t stride = (size_t)n + 1;
    double *aug = (double *)calloc((size_t)n * stride, sizeof(*aug));
    if (aug == NULL) return -1;

    double scale = 0.0;
    for (int t = 0; t < nt; t++) {
        int r = rows[t], c = cols[t];
        double value = vals[t];
        if (r < 0 || r >= n || c < 0 || c > r || !isfinite(value)) {
            free(aug);
            return -1;
        }
        aug[(size_t)r * stride + (size_t)c] += value;
        if (r != c) aug[(size_t)c * stride + (size_t)r] += value;
        if (fabs(value) > scale) scale = fabs(value);
    }
    for (int r = 0; r < n; r++) {
        if (!isfinite(b[r])) {
            free(aug);
            return -1;
        }
        aug[(size_t)r * stride + (size_t)n] = b[r];
    }

    double pivot_floor = DBL_EPSILON * fmax(1.0, scale) * (double)n;
    for (int k = 0; k < n; k++) {
        int pivot = k;
        double pivot_abs = fabs(aug[(size_t)k * stride + (size_t)k]);
        for (int r = k + 1; r < n; r++) {
            double candidate = fabs(aug[(size_t)r * stride + (size_t)k]);
            if (candidate > pivot_abs) {
                pivot = r;
                pivot_abs = candidate;
            }
        }
        if (!(pivot_abs > pivot_floor) || !isfinite(pivot_abs)) {
            free(aug);
            return -1;
        }
        if (pivot != k) {
            for (int c = k; c <= n; c++) {
                double tmp = aug[(size_t)k * stride + (size_t)c];
                aug[(size_t)k * stride + (size_t)c] =
                    aug[(size_t)pivot * stride + (size_t)c];
                aug[(size_t)pivot * stride + (size_t)c] = tmp;
            }
        }
        double diagonal = aug[(size_t)k * stride + (size_t)k];
        for (int r = k + 1; r < n; r++) {
            double factor = aug[(size_t)r * stride + (size_t)k] / diagonal;
            aug[(size_t)r * stride + (size_t)k] = 0.0;
            for (int c = k + 1; c <= n; c++)
                aug[(size_t)r * stride + (size_t)c] -=
                    factor * aug[(size_t)k * stride + (size_t)c];
        }
    }

    for (int r = n - 1; r >= 0; r--) {
        double rhs = aug[(size_t)r * stride + (size_t)n];
        for (int c = r + 1; c < n; c++)
            rhs -= aug[(size_t)r * stride + (size_t)c] * x[c];
        double diagonal = aug[(size_t)r * stride + (size_t)r];
        if (!(fabs(diagonal) > pivot_floor) || !isfinite(rhs)) {
            free(aug);
            return -1;
        }
        x[r] = rhs / diagonal;
    }
    free(aug);
    return 0;
}

int Sparse_solve_sym(int n, int nt,
                     const int *rows, const int *cols, const double *vals,
                     const double *b, double *x, SparseMode mode)
{
    if (n <= 0 || nt < 0) return -1;
    if (!rows || !cols || !vals || !b || !x) return -1;
    if (n <= 32) return sparse_solve_small_dense(n, nt, rows, cols, vals, b, x);



    /* --- Bucket triplets by column. --- */
    int *colcnt = (int *)calloc((size_t)n + 1, sizeof(int));
    if (!colcnt) return -1;
    for (int t = 0; t < nt; t++) {
        int c = cols[t];
        if (c < 0 || c >= n || rows[t] < c || rows[t] >= n) { free(colcnt); return -1; }
        colcnt[c]++;
    }
    int *cstart = (int *)malloc(((size_t)n + 1) * sizeof(int));
    if (!cstart) { free(colcnt); return -1; }
    cstart[0] = 0;
    for (int c = 0; c < n; c++) cstart[c + 1] = cstart[c] + colcnt[c];

    RV *buf = (RV *)malloc((size_t)(nt > 0 ? nt : 1) * sizeof(RV));
    int *cur = (int *)malloc((size_t)n * sizeof(int));
    if (!buf || !cur) { free(colcnt); free(cstart); free(buf); free(cur); return -1; }
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
    if (!A) { free(colcnt); free(cstart); free(buf); free(cur); return -1; }

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

    /* --- Solve. --- */
    char *opt_spd[]   = { "taucs.factor.LLT=true", "taucs.factor.mf=true", NULL };
    char *opt_indef[] = { "taucs.solve.minres=true",
                          "taucs.solve.maxits=10000",
                          "taucs.solve.convergetol=1e-11", NULL };
    void *F = NULL;
    int rc = taucs_linsolve(A, &F, 1, x, (void *)b,
                            (mode == SPARSE_SPD) ? opt_spd : opt_indef, NULL);
    if (F) { void *Fp = F; taucs_linsolve(NULL, &Fp, 0, NULL, NULL, NULL, NULL); }
    taucs_ccs_free(A);
    return (rc == TAUCS_SUCCESS) ? 0 : -1;
}
