/*
 * taucs_smoke.c -- Stage-1 smoke test for the in-tree TAUCS + CLAPACK build.
 *
 * Proves the sparse-solver foundation before any ABF++ code:
 *   1. SPD system  : 1D Laplacian, supernodal Cholesky (exercises CLAPACK
 *                    dpotrf_/dgemm_/dsyrk_/dtrsm_ via TAUCS's multifrontal LLT).
 *   2. Indefinite  : symmetric indefinite 3x3, unpreconditioned MINRES
 *                    (the path ABF++'s Newton system will use).
 *
 * Usage: taucs_smoke [--selftest]   (always runs the tests)
 * Exit 0 = all pass, 1 = failure.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "taucs.h"

/* taucs.h does not prototype these high-level entry points; declare them. */
extern taucs_ccs_matrix *taucs_ccs_create(int m, int n, int nnz, int flags);
extern void             taucs_ccs_free(taucs_ccs_matrix *m);
extern int  taucs_linsolve(taucs_ccs_matrix *A, void **F, int nrhs,
                           void *X, void *B, char *options[], void *opt_arg[]);

/* Free a factorization handle (the (A==NULL, nrhs==0) convention). */
static void free_factor(void *F)
{
    if (F) {
        void *Fp = F;
        taucs_linsolve(NULL, &Fp, 0, NULL, NULL, NULL, NULL);
    }
}

/* Test 1: SPD tridiagonal Laplacian (diag 2, off-diag -1), N unknowns.
 * b = A*ones = [1,0,...,0,1], so the solution is the all-ones vector.
 * Solved with the supernodal multifrontal Cholesky (BLAS3 + dpotrf). */
static int test_spd(void)
{
    const int N = 200;
    const int nnz = 2 * N - 1;
    taucs_ccs_matrix *A = taucs_ccs_create(N, N, nnz,
                              TAUCS_DOUBLE | TAUCS_SYMMETRIC | TAUCS_LOWER);
    if (!A) { fprintf(stderr, "  [spd] ccs_create failed\n"); return 1; }

    int k = 0;
    for (int j = 0; j < N; j++) {
        A->colptr[j] = k;
        A->rowind[k] = j;     A->values.d[k] = 2.0;  k++;   /* diagonal */
        if (j < N - 1) { A->rowind[k] = j + 1; A->values.d[k] = -1.0; k++; }
    }
    A->colptr[N] = k;

    double *b = (double *)malloc(sizeof(double) * (size_t)N);
    double *x = (double *)malloc(sizeof(double) * (size_t)N);
    if (!b || !x) { fprintf(stderr, "  [spd] oom\n"); taucs_ccs_free(A); free(b); free(x); return 1; }
    for (int i = 0; i < N; i++) { b[i] = 0.0; x[i] = 0.0; }
    b[0] = 1.0; b[N - 1] = 1.0;

    char *opt[] = { "taucs.factor.LLT=true", "taucs.factor.mf=true", NULL };
    void *F = NULL;
    int rc = taucs_linsolve(A, &F, 1, x, b, opt, NULL);

    double err = 0.0;
    if (rc == TAUCS_SUCCESS) {
        for (int i = 0; i < N; i++) { double e = fabs(x[i] - 1.0); if (e > err) err = e; }
    }
    free_factor(F);
    int fail = (rc != TAUCS_SUCCESS) || (err > 1e-8);
    fprintf(stderr, "  [spd]   N=%d rc=%d max|x-1|=%.3e -> %s\n",
            N, rc, err, fail ? "FAIL" : "ok");
    taucs_ccs_free(A); free(b); free(x);
    return fail ? 1 : 0;
}

/* Test 2: symmetric INDEFINITE 3x3 via unpreconditioned MINRES.
 * A = [[2,1,0],[1,-1,1],[0,1,2]], x=[1,2,3] -> b=[4,2,8]. */
static int test_indef(void)
{
    taucs_ccs_matrix *A = taucs_ccs_create(3, 3, 5,
                              TAUCS_DOUBLE | TAUCS_SYMMETRIC | TAUCS_LOWER);
    if (!A) { fprintf(stderr, "  [indef] ccs_create failed\n"); return 1; }

    int k = 0;
    A->colptr[0] = k; A->rowind[k] = 0; A->values.d[k] = 2.0; k++;
                      A->rowind[k] = 1; A->values.d[k] = 1.0; k++;
    A->colptr[1] = k; A->rowind[k] = 1; A->values.d[k] = -1.0; k++;
                      A->rowind[k] = 2; A->values.d[k] = 1.0; k++;
    A->colptr[2] = k; A->rowind[k] = 2; A->values.d[k] = 2.0; k++;
    A->colptr[3] = k;

    double b[3] = { 4.0, 2.0, 8.0 };
    double x[3] = { 0.0, 0.0, 0.0 };
    const double xexp[3] = { 1.0, 2.0, 3.0 };

    /* No LLT/LU (would fail on an indefinite matrix) -> f->type NONE ->
     * unpreconditioned MINRES on the (ordered) matrix. */
    char *opt[] = { "taucs.solve.minres=true",
                    "taucs.solve.maxits=1000",
                    "taucs.solve.convergetol=1e-12", NULL };
    void *F = NULL;
    int rc = taucs_linsolve(A, &F, 1, x, b, opt, NULL);

    double err = 0.0;
    if (rc == TAUCS_SUCCESS) {
        for (int i = 0; i < 3; i++) { double e = fabs(x[i] - xexp[i]); if (e > err) err = e; }
    }
    free_factor(F);
    int fail = (rc != TAUCS_SUCCESS) || (err > 1e-6);
    fprintf(stderr, "  [indef] rc=%d max|x-xexp|=%.3e -> %s\n",
            rc, err, fail ? "FAIL" : "ok");
    taucs_ccs_free(A);
    return fail ? 1 : 0;
}

int main(int argc, char **argv)
{
    (void)argc; (void)argv;
    fprintf(stderr, "taucs_smoke: TAUCS + CLAPACK solver foundation\n");
    int fails = 0;
    fails += test_spd();
    fails += test_indef();
    fprintf(stderr, "taucs_smoke: %s (%d failure%s)\n",
            fails ? "FAIL" : "PASS", fails, fails == 1 ? "" : "s");
    return fails ? 1 : 0;
}
