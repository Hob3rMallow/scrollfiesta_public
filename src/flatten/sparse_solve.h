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

/* Largest dimension for the one-shot SPD solver's dense fast path (0 disables
 * it). Callers preserving its arithmetic when adding factor reuse can retain
 * one-shot solves at this size; Sparse_factor_spd uses TAUCS for n > 1. */
int Sparse_spd_dense_limit(void);

/* ---- factor-once / solve-many (SPD only) ---------------------------------
 * For loops that solve the SAME matrix against many right-hand sides (ADMM
 * iterations, IRLS rounds, multigrid-free preconditioning): factor once with
 * supernodal multifrontal Cholesky, then each solve is two triangular sweeps.
 * The handle owns the permuted matrix and the TAUCS factorization; free it
 * with Sparse_factor_free.  All functions return 0 on success. */
typedef struct SparseFactor *SparseFactor_T;

int  Sparse_factor_spd(int n, int nt,
                       const int *rows, const int *cols, const double *vals,
                       SparseFactor_T *out);
int  Sparse_factor_solve(SparseFactor_T f, const double *b, double *x);
/* nrhs right-hand sides, column-major contiguous ([n*nrhs]); one call
 * amortizes the solver's per-call permutation across the batch. */
int  Sparse_factor_solve_multi(SparseFactor_T f, const double *B, double *X,
                               int nrhs);
void Sparse_factor_free(SparseFactor_T *f);
/* A caller-owned CCS pattern: the lower triangle (row >= col), columns in
 * order, rows ascending inside each column, every diagonal entry present.
 * The factor keeps its own copy of the pattern, its ordering and its
 * symbolic factor; Sparse_factor_refactor_ccs replaces the values (same
 * pattern, same order) and pays the numeric factorisation only.  For loops
 * that assemble straight into a fixed pattern (SLIM
 * global steps): no triplet sort, no ordering, no pattern-cache traffic. */
int  Sparse_factor_spd_ccs(int n, const int *colptr, const int *rowind,
                           const double *values, SparseFactor_T *out);
int  Sparse_factor_refactor_ccs(SparseFactor_T f, const double *values);
/* (2026-10-27) Repeated SPD solves of ONE triplet pattern (a SLIM's iterations on one chart) through a caller-owned cache
 * (NULL to start; one thread): the first solve makes the triplets' CCS slots, the ordering and the symbolic factor,
 * later solves of the same pattern refactor numerically.  The arithmetic is Sparse_solve_sym(SPARSE_SPD)'s, bit for
 * bit (its dense path at n <= SPARSE_DENSE_MAX is used as is); a different pattern starts the cache over.  0 ok. */
typedef struct SparseRepeat *SparseRepeat_T;
int  Sparse_repeat_solve(SparseRepeat_T *cache, int n, int nt, const int *rows, const int *cols, const double *vals,
                         const double *b, double *x);
void Sparse_repeat_free(SparseRepeat_T *cache);
/* Symbolic-factor pattern cache counters (see sparse_solve.c). */
void Sparse_factor_cache_stats(unsigned long long *hits, unsigned long long *misses);

/* Unit test: factor reuse across multiple RHS agrees with the one-shot
 * solver; exercises the n==1 scalar path.  Returns 0 on pass. */
int Sparse_selftest(void);
/* Concurrent factorizations of side x side grid Laplacians (tool benchmark). */
int Sparse_benchmark(int side, int systems, int threads);

#endif
