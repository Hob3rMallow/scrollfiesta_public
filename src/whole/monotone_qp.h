#ifndef MONOTONE_QP_INCLUDED
#define MONOTONE_QP_INCLUDED

#include <stddef.h>
#include <stdint.h>

#include "../common/arena.h"

/*
 * A convex quadratic written as weighted least-squares rows:
 *
 *   min_x  1/2 sum_r weight_r (a_r^T x - target_r)^2
 *
 * subject to chain-oriented lower bounds
 *
 *   x[hi] - x[lo] >= lower
 *
 * and exact scalar anchors.  The implementation is a primal active-set
 * method specialized to these difference constraints.  Active differences
 * are collapsed into blocks and every numerical solve is reduced SPD
 * Cholesky; no indefinite KKT matrix and no post-solve projection are used.
 */

typedef struct {
    int32_t var;
    double  value;
} MonotoneQpCoeff;

typedef struct {
    size_t  first;
    int32_t count;
    double  target;
    double  weight;
    int32_t kind;
    int32_t owner;
} MonotoneQpRow;

typedef struct {
    int32_t lo;
    int32_t hi;
    double  lower;
    int32_t stroke;
    int32_t edge;
} MonotoneQpBound;

typedef struct {
    int32_t var;
    double  value;
    int32_t component;
} MonotoneQpAnchor;

typedef struct {
    size_t nvar;

    const MonotoneQpRow   *rows;
    size_t                 nrows;
    const MonotoneQpCoeff *coeff;
    size_t                 ncoeff;

    const MonotoneQpBound *bounds;
    size_t                 nbounds;
    const MonotoneQpAnchor *anchors;
    size_t                  nanchors;
} MonotoneQpProblem;

typedef struct {
    int    max_active_iterations;
    double feasibility_tolerance;
    double direction_tolerance;
    double multiplier_tolerance;
    double linear_residual_tolerance;
    int    verbose;
} MonotoneQpOptions;

typedef struct {
    int    iteration;
    int    n_active;
    int    added_bound;
    int    dropped_bound;
    double objective;
    double max_violation;
    double direction_norm;
    double step;
    double reduced_linear_residual;
} MonotoneQpTraceEntry;

typedef struct {
    MonotoneQpTraceEntry *entry;
    size_t capacity;
    size_t count;
} MonotoneQpTrace;

typedef struct {
    int    iterations;
    int    spd_solves;
    int    multiplier_solves;
    int    active_final;
    int    added_total;
    int    dropped_total;
    double objective_initial;
    double objective_final;
    double max_violation;
    double min_slack;
    double stationarity_residual;
    double max_reduced_linear_residual;
} MonotoneQpStats;

void MonotoneQpOptions_default(MonotoneQpOptions *opts);

/*
 * x is both the feasible initial point and the solution.  Anchored entries
 * must already equal their requested values, and all lower bounds must be
 * feasible to opts->feasibility_tolerance.
 *
 * Returns:
 *   0  converged;
 *  -1  malformed problem;
 *  -2  infeasible initial point;
 *  -3  reduced SPD solve failed;
 *  -4  active-set iteration limit.
 */
int MonotoneQp_solve(Arena_T arena,
                     const MonotoneQpProblem *problem,
                     const MonotoneQpOptions *opts,
                     double *x,
                     MonotoneQpTrace *trace,
                     MonotoneQpStats *stats);

/*
 * Exact minimizer of the quadratic objective subject to the ANCHORS ONLY --
 * bounds are ignored (the caller projects afterwards).  One reduced SPD
 * solve; no active set, no iteration limit, no feasible start required.
 * Anchored entries of x are overwritten with their requested values.
 */
int MonotoneQp_solve_ls(Arena_T arena,
                        const MonotoneQpProblem *problem,
                        double *x,
                        double *out_linear_residual);

double MonotoneQp_objective(const MonotoneQpProblem *problem,
                            const double *x);

#endif
