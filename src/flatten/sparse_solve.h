#ifndef SPARSE_SOLVE_INCLUDED
#define SPARSE_SOLVE_INCLUDED

/* ============================================================================
 * sparse_solve.h -- thin wrapper over TAUCS for the flatten module.
 *
 * Input is a SYMMETRIC matrix given as lower-triangle COO triplets (row >= col;
 * duplicate (row,col) entries are summed). Builds a TAUCS CCS matrix and solves
 * A x = b:
 *   mode 0 (SPD)        -> supernodal multifrontal Cholesky (CLAPACK BLAS3).
 *   mode 1 (indefinite) -> unpreconditioned MINRES (for ABF++'s Newton system).
 *
 * All arithmetic is double (CLAUDE.md rule 20). Returns 0 on success.
 * ==========================================================================*/

typedef enum {
    SPARSE_SPD = 0,        /* symmetric positive definite -> Cholesky */
    SPARSE_SYM_INDEFINITE  /* symmetric (possibly indefinite) -> MINRES */
} SparseMode;

/* n: matrix dimension. nt: number of COO triplets (lower triangle, row>=col).
 * rows/cols/vals: [nt]. b: [n] (rhs). x: [n] (solution, also initial guess).
 * Returns 0 on success, non-zero on failure. */
int Sparse_solve_sym(int n, int nt,
                     const int *rows, const int *cols, const double *vals,
                     const double *b, double *x, SparseMode mode);

#endif
