#ifndef VES_ACTIVE_SET_QP_INCLUDED
#define VES_ACTIVE_SET_QP_INCLUDED
#include <stddef.h>

/* ============================================================================
 * active_set_qp.h -- a strictly convex QP with few binding inequalities,
 * solved by a dual active set on ONE sparse Cholesky factor.
 *
 *     minimise  1/2 x'Hx + g'x      subject to  a_i'x <= u_i,  i = 1..m
 *
 * Same calling convention as ConvexQp_solve so the two are interchangeable:
 * H is lower-triangle COO holding TWICE the true Hessian (duplicates sum),
 * `target` is -2g, A is general COO with `upper` bounds, and `result`
 * receives x.  Returns 1 on success, 0 when the problem is infeasible, H is
 * not positive definite after a ridge, or the working set outgrows its cap.
 *
 * Why this and not ADMM: H is factored once (TAUCS supernodal Cholesky) and
 * every iteration is two triangular solves plus a dense system the size of
 * the WORKING SET -- the constraints that actually bind -- so the cost scales
 * with the binding contacts, not with the candidate rows, and the answer is
 * the exact optimum of the linearised model rather than an iterate cut off
 * by an iteration budget.  Goldfarb & Idnani (1983): start at the
 * unconstrained minimiser, add the most violated constraint, take the
 * longer of the primal step to feasibility and the dual step that would
 * zero a multiplier, dropping constraints whose multipliers reach zero;
 * finite termination for positive definite H.
 * ==========================================================================*/

/* Why a solve was refused (ActiveSetQpStats.reason; 0 = solved). */
enum {
    ASQP_REASON_OK = 0,
    ASQP_REASON_INPUT,        /* invalid or non-finite input */
    ASQP_REASON_MEMORY,       /* allocation failed */
    ASQP_REASON_FACTOR,       /* H would not factor, ridged or not */
    ASQP_REASON_SOLVE,        /* a triangular solve failed or went non-finite */
    ASQP_REASON_INFEASIBLE,   /* the rows admit no point */
    ASQP_REASON_STEPS,        /* the step guard was exhausted */
    ASQP_REASON_CAP,          /* the working set outgrew ASQP_MAX_ACTIVE */
    ASQP_REASON_DENSE,        /* the working-set Cholesky lost positive definiteness */
    ASQP_REASON_ILLPOSED,     /* the objective passed the caller's limit (ActiveSetQp_limit) */
    ASQP_REASON_TIME          /* the solve outran ASQP_MAX_SECONDS and was going nowhere */
};

typedef struct ActiveSetQpStats {
    int    iterations;      /* add/drop steps taken */
    int    active;          /* binding constraints at the solution */
    int    dropped;         /* constraints dropped along the way */
    int    ridge_retries;   /* factorisations that needed a diagonal ridge */
    int    reason;          /* ASQP_REASON_* */
    int    warm;            /* rows seeded from the caller's warm start */
    double factor_seconds;  /* time in the Cholesky factorisation */
    double solve_seconds;   /* time in triangular solves */
} ActiveSetQpStats;

const char *ActiveSetQp_reason(int reason);

int ActiveSetQp_solve(int n, int hn, const int *hr, const int *hc, const double *hv,
                      const double *target, int m, int an, const int *ar, const int *ac,
                      const double *av, const double *upper, double *result);

/* Factor once, solve many: the clearance QP re-solves the same H against a
 * growing set of cutting-plane rows.  ActiveSetQp_factor returns 1 and a
 * factor object (stats: factor time, ridge retries, or the reason it
 * refused); ActiveSetQp_solve_factored solves one row set on it.  The
 * one-shot ActiveSetQp_solve is factor + solve + free. */
typedef struct ActiveSetQpFactor *ActiveSetQpFactor_T;
int  ActiveSetQp_factor(int n, int hn, const int *hr, const int *hc, const double *hv,
                        ActiveSetQpFactor_T *out, ActiveSetQpStats *stats);
int  ActiveSetQp_solve_factored(ActiveSetQpFactor_T factor, const double *target, int m, int an,
                                const int *ar, const int *ac, const double *av, const double *upper,
                                double *result, ActiveSetQpStats *stats);
/* As above with a warm start: `rows[0..*count)` are the rows to seed the
 * working set with (those binding at the previous, similar solve); on
 * success the final working set is written back, at most `capacity` rows.
 * Seeding costs one triangular solve per row instead of two per step. */
int  ActiveSetQp_solve_warm(ActiveSetQpFactor_T factor, const double *target, int m, int an,
                            const int *ar, const int *ac, const double *av, const double *upper,
                            double *result, ActiveSetQpStats *stats, int *rows, int *count, int capacity);
void ActiveSetQp_free(ActiveSetQpFactor_T *factor);
/* The QP's value at its UNCONSTRAINED minimiser, -1/2 g'H^-1 g (one solve).
 * Every constrained optimum lies above it, so a value far below the scale
 * of the energy being modelled says the linearisation is ill-posed. */
int  ActiveSetQp_model(ActiveSetQpFactor_T factor, const double *target, double *model);
/* The dual active set raises the objective monotonically toward the
 * constrained optimum, so a running value above `limit` proves the optimum
 * lies above it too: the solve is abandoned there (ASQP_REASON_ILLPOSED)
 * instead of chasing a step the caller will not take.  Default: no limit. */
void ActiveSetQp_limit(ActiveSetQpFactor_T factor, double limit);
/* Absolute ves_clock_sec deadline across factors, warm starts and solves.
 * 0 disables; process-local, so callers must serialize QP use. */
void ActiveSetQp_deadline(double deadline);

/* As above, also reporting how the solve went. `stats` may be NULL. */
int ActiveSetQp_solve_stats(int n, int hn, const int *hr, const int *hc, const double *hv,
                            const double *target, int m, int an, const int *ar, const int *ac,
                            const double *av, const double *upper, double *result,
                            ActiveSetQpStats *stats);

/* KKT checks on hand cases and random problems, plus agreement with the
 * OSQP path where it is linked.  Returns the failure count. */
int ActiveSetQp_selftest(void);

#endif
