#ifndef SPARSE_LOCAL_GLOBAL_INCLUDED
#define SPARSE_LOCAL_GLOBAL_INCLUDED

#include "../common/arena.h"
#include <stddef.h>

/* Convex least squares with local half-space copies. H is SPD and supplied
 * as lower-triangle triplets (duplicates sum). Minimize .5*x'H*x-b'x subject
 * to C*x >= lower. Six slots allow a row on one triangle, including aliases.
 * Geometry and nonlinear acceptance remain the caller's responsibility. */
typedef struct SparseLocalRow {
    int index[6];
    double value[6], lower;
} SparseLocalRow;

typedef struct SparseLocalGlobalStats {
    size_t iterations, factors, solves, penalty_changes, local_refinements;
    int converged, warm;
    double primal_residual, stationarity_residual, complementarity;
    double objective, penalty, linear_residual;
} SparseLocalGlobalStats;

/* Alternating projection / global least squares (scaled ADMM with convex
 * augmented-least-squares refinement). Numeric
 * factors are reused only while H, C and the penalty are unchanged. Row and
 * column equilibration do not alter the problem. Optional starting point and
 * nonnegative multipliers are in ORIGINAL units. Outputs are changed only
 * after original-unit feasibility and optimality certificates pass.
 * Return 0: certified, 1: work limit (unfinished), -1: invalid/numerical error.
 * A work limit, small step or constant loss alone never means convergence. */
int Sparse_local_global(Arena_T arena, size_t n, size_t nt,
    const int *row, const int *col, const double *h, const double *b,
    size_t nr, const SparseLocalRow *constraints,
    const double *initial, const double *initial_multipliers, int limit,
    double *out, double *out_multipliers, SparseLocalGlobalStats *stats);

int Sparse_local_global_selftest(void);

#endif
