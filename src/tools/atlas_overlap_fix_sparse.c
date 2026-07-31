/*
 * atlas_overlap_fix_sparse.c
 *
 * Sparse_solve_sym for the standalone overlap-fix tool: Jacobi-preconditioned
 * conjugate gradients over the COO-symmetric input.  atlas_register's normal
 * equations here are a few-hundred-unknown graph system, far below anything
 * that justifies linking TAUCS (and its 512 MB stack demand) into this tool.
 * All arithmetic is double (CLAUDE.md rule 20).
 */

#include "../flatten/sparse_solve.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

static void coo_matvec(int n, int nt, const int *rows, const int *cols,
                       const double *vals, const double *x, double *y)
{
    memset(y, 0, (size_t)n * sizeof(double));
    for (int t = 0; t < nt; t++) {
        int r = rows[t], c = cols[t];
        y[r] += vals[t] * x[c];
        if (r != c) y[c] += vals[t] * x[r];
    }
}

int Sparse_solve_sym(int n, int nt,
                     const int *rows, const int *cols, const double *vals,
                     const double *b, double *x, SparseMode mode)
{
    (void)mode;   /* CG covers the SPD case; the MINRES caller (ABF++) is
                   * not part of this tool */
    if (n <= 0 || nt < 0) return -1;

    double *scratch = (double *)calloc((size_t)n * 5, sizeof(double));
    if (scratch == NULL) return -1;
    double *r = scratch;
    double *z = scratch + n;
    double *p = scratch + 2 * n;
    double *ap = scratch + 3 * n;
    double *diag = scratch + 4 * n;

    for (int t = 0; t < nt; t++)
        if (rows[t] == cols[t]) diag[rows[t]] += vals[t];

    coo_matvec(n, nt, rows, cols, vals, x, ap);
    double bnorm2 = 0.0;
    for (int i = 0; i < n; i++) {
        r[i] = b[i] - ap[i];
        bnorm2 += b[i] * b[i];
    }
    double tol2 = 1.0e-20 * (bnorm2 > 1.0 ? bnorm2 : 1.0);

    double rz = 0.0;
    for (int i = 0; i < n; i++) {
        z[i] = diag[i] > 1.0e-30 ? r[i] / diag[i] : r[i];
        p[i] = z[i];
        rz += r[i] * z[i];
    }

    int ok = 0;
    int max_iter = 10 * n + 100;
    for (int iter = 0; iter < max_iter; iter++) {
        double rr = 0.0;
        for (int i = 0; i < n; i++) rr += r[i] * r[i];
        if (rr <= tol2) { ok = 1; break; }

        coo_matvec(n, nt, rows, cols, vals, p, ap);
        double pap = 0.0;
        for (int i = 0; i < n; i++) pap += p[i] * ap[i];
        if (!(pap > 0.0)) break;   /* lost positive definiteness */
        double alpha = rz / pap;
        for (int i = 0; i < n; i++) {
            x[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
        }
        double rz_new = 0.0;
        for (int i = 0; i < n; i++) {
            z[i] = diag[i] > 1.0e-30 ? r[i] / diag[i] : r[i];
            rz_new += r[i] * z[i];
        }
        double beta = rz > 0.0 ? rz_new / rz : 0.0;
        for (int i = 0; i < n; i++) p[i] = z[i] + beta * p[i];
        rz = rz_new;
    }

    free(scratch);
    return ok ? 0 : -1;
}
