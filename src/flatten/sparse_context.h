#ifndef SPARSE_CONTEXT_INCLUDED
#define SPARSE_CONTEXT_INCLUDED

/* One serial caller owns this context for a complete nonlinear solve.
 * Only ordering/symbolic structure persists; every changed matrix receives
 * a new numeric factor. Exact CCS-pattern comparison controls invalidation. */
typedef struct SparseContext *SparseContext_T;
typedef struct SparseContextStats {
    unsigned long long calls, pattern_hits, pattern_builds, invalidations;
    unsigned long long right_hand_sides;
    double assembly_seconds, symbolic_seconds, numeric_seconds, solve_seconds;
} SparseContextStats;

/* A successful factor may serve several right-hand sides. Any attempted
 * refactor invalidates that numeric factor until the new one succeeds. */
int Sparse_context_factor(SparseContext_T *context, int n, int nt,
    const int *rows, const int *cols, const double *values);
int Sparse_context_solve_factored(SparseContext_T context,
    const double *rhs, double *result);
int Sparse_context_solve(SparseContext_T *context, int n, int nt,
    const int *rows, const int *cols, const double *values,
    const double *rhs, double *result);
void Sparse_context_stats(SparseContext_T context, SparseContextStats *stats);
void Sparse_context_free(SparseContext_T *context);
int Sparse_context_selftest(void);

#endif
