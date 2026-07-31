#ifndef ATLAS_STRIP_INCLUDED
#define ATLAS_STRIP_INCLUDED

#include <stddef.h>
#include <stdint.h>

#include "../common/arena.h"
#include "monotone_qp.h"

/*
 * StrokeStrip-style coarse discretization of a scroll atlas.
 *
 * Samples are ordered along constant-axial slice curves ("strokes").
 * Cross-sections collect observations of one latent material ruling.  The
 * unknown vector is [sample u values, one latent intercept q per cross-section].
 */

typedef struct {
    double p[3];
    double tangent[3];
    double s;
    int32_t stroke;
    int32_t ordinal;
    double confidence;
} AtlasStripSample;

typedef struct {
    int32_t first;
    int32_t count;
    int32_t component;
    int32_t resolved;
} AtlasStripStroke;

/*
 * A cross-section observation may lie directly on a sample (value1 < 0) or
 * interpolate between value0/value1 using value_t in [0,1].  Its longitudinal
 * derivative is measured on the oriented edge deriv_lo -> deriv_hi.
 */
typedef struct {
    int32_t value0;
    int32_t value1;
    double  value_t;

    int32_t deriv_lo;
    int32_t deriv_hi;
    double  deriv_length;

    double p[3];
    double tangent[3];
    double dual_width;
    double membership;  /* geometry/context prior or likelihood in [0,1] */
    double base_weight;
    int32_t observation;
} AtlasStripMember;

typedef struct {
    size_t first;
    int32_t count;
    double tangent[3];
    double weight;
    int32_t id;
} AtlasStripCrossSection;

/*
 * A short, geometry-certified continuation between two ordered-stroke
 * endpoints on the same axial slice.  Unlike a cross-section observation,
 * the endpoints are not the same isovalue: the tangent-projected physical
 * gap is retained in both relaxed and final solves,
 *
 *     u[a] - u[b] = target ~= tau dot (p_a - p_b).
 *
 * This is StrokeStrip's fragment-continuation row.  Candidate construction
 * should only emit mutually selected, sub-half-wrap joins; uncertain
 * cross-slice alternatives belong in cross_sections instead.
 */
typedef struct {
    int32_t a;
    int32_t b;
    double target;
    double weight;
    double membership;
    int32_t observation;
} AtlasStripContinuation;

typedef struct {
    const AtlasStripSample *samples;
    size_t nsamples;
    const AtlasStripStroke *strokes;
    size_t nstrokes;

    const AtlasStripMember *members;
    size_t nmembers;
    const AtlasStripCrossSection *cross_sections;
    size_t ncross_sections;

    const AtlasStripContinuation *continuations;
    size_t ncontinuations;

    const MonotoneQpAnchor *anchors;
    size_t nanchors;

    /*
     * Optional per-sample prior targets (u_i ~= prior_u[i], one singleton
     * row per finite entry, weight lambda_prior).  The paper's zero-target
     * registration is only correct for overdrawn strokes; on a coiled sheet
     * an undetermined fragment gauge must resolve toward the geometric
     * spiral, not toward stacking.  NaN entries emit no row (trust-gated
     * upstream).  NULL or lambda_prior == 0 disables.
     */
    const double *prior_u;
} AtlasStripProblem;

typedef enum {
    ATLAS_STRIP_RELAXED = 0,
    ATLAS_STRIP_FINAL = 1
} AtlasStripMode;

enum {
    ATLAS_STRIP_ROW_LENGTH = 1,
    ATLAS_STRIP_ROW_ALIGN = 2,
    ATLAS_STRIP_ROW_LOCAL = 3,
    ATLAS_STRIP_ROW_CONTINUATION = 4,
    ATLAS_STRIP_ROW_PRIOR = 5
};

typedef struct {
    AtlasStripMode mode;
    double lambda_length;
    double lambda_align;
    double lambda_local;
    double lambda_continuation;
    double monotone_fraction;
    double membership_floor;
    double lambda_prior;    /* weight of prior_u singleton rows; 0 = off.
                             * Must dominate the ~1e-5 robust-floor couplings
                             * and be dominated by evidence rows (0.05..24):
                             * it decides which wind an undetermined fragment
                             * sits on, never its shape. */
} AtlasStripOptions;

typedef struct {
    MonotoneQpProblem qp;
    size_t nsamples;
    size_t ncross_sections;
} AtlasStripSystem;

typedef struct {
    double energy_length;
    double energy_align;
    double energy_local;
    double energy_continuation;
    double rms_length_row;
    double rms_align_row;
    double rms_local_speed_error;
    double rms_continuation;
    double max_align_residual;
    double max_local_speed_error;
    double max_continuation_residual;
    double min_monotone_ratio;
    size_t active_members;
    size_t active_continuations;
    double energy_prior;
    double rms_prior_row;
    size_t active_priors;
} AtlasStripMetrics;

/*
 * StrokeStrip-style local/global relaxation.  The first global solve uses a
 * deliberately tiny alignment likelihood, followed by IRLS approximations to
 * the paper's robust L1 initialization.  Subsequent rounds update each
 * member's WTS-adjacency likelihood from its residual to the latent
 * cross-section intercept and solve another monotone QP.
 *
 * likelihood_sigma <= 0 derives sigma from the longest stroke span / 30.
 */
typedef struct {
    int    l1_iterations;
    int    likelihood_iterations;
    double initial_likelihood;
    double l1_epsilon;
    double likelihood_sigma;
    double likelihood_floor;
    double convergence_tolerance;
    int    use_active_set;  /* 1 = the historical MonotoneQp active-set path.
                             * Default 0: anchors-only exact LS solve followed
                             * by per-stroke PAVA projection -- one SPD
                             * factorization per round, no iteration cap, no
                             * rc -2/-3/-4 failure modes. */
} AtlasStripRobustOptions;

enum {
    ATLAS_STRIP_ROBUST_INITIAL = 0,
    ATLAS_STRIP_ROBUST_L1 = 1,
    ATLAS_STRIP_ROBUST_LIKELIHOOD = 2
};

typedef struct {
    int    iteration;
    int    phase;
    int    qp_rc;
    double objective;
    double max_u_change;
    double max_membership_change;
    double residual_rms;
    double max_residual;
    double membership_min;
    double membership_mean;
    double membership_max;
    size_t downweighted_members;
} AtlasStripRobustTraceEntry;

/* Called after each robust coarse solve while x and the memberships used by
 * that solve are still live.  This is intentionally a read-only observation
 * hook: diagnostics must never be able to perturb the alternating solve. */
typedef int (*AtlasStripRobustIterationFn)(
    void *context, const AtlasStripProblem *problem, const double *x,
    const AtlasStripMember *members,
    const AtlasStripRobustTraceEntry *iteration);

typedef struct {
    AtlasStripRobustTraceEntry *entry;
    size_t capacity;
    size_t count;
    AtlasStripRobustIterationFn iteration_fn;
    void *iteration_context;
} AtlasStripRobustTrace;

typedef struct {
    int total_qp_solves;
    int initial_solves;
    int l1_solves;
    int likelihood_solves;
    double likelihood_sigma;
    double max_u_change;
    double max_membership_change;
    size_t downweighted_members;
    MonotoneQpStats final_qp;
} AtlasStripRobustStats;

void AtlasStripOptions_default(AtlasStripOptions *opts);
void AtlasStripRobustOptions_default(AtlasStripRobustOptions *opts);

/*
 * Assemble the quadratic rows and monotonicity bounds.  All returned arrays
 * live in arena.  The caller solves system.qp with MonotoneQp_solve.
 */
int AtlasStrip_build(Arena_T arena,
                     const AtlasStripProblem *problem,
                     const AtlasStripOptions *opts,
                     AtlasStripSystem *system);

/* Initialize q variables in x=[u,q] from the current sample u values. */
void AtlasStrip_initialize_intercepts(const AtlasStripProblem *problem,
                                      const AtlasStripOptions *opts,
                                      double *x);

/*
 * Exact Euclidean projection of the sample entries of x onto the per-stroke
 * monotonicity cone u[i+1]-u[i] >= monotone_fraction*(s[i+1]-s[i]), by
 * weighted pool-adjacent-violators on z_i = u_i - monotone_fraction*s_i.
 * Anchored samples carry a huge weight so a component's gauge pin cannot
 * move.  Intercept entries are untouched.
 */
void AtlasStrip_project_monotone(Arena_T arena,
                                 const AtlasStripProblem *problem,
                                 const AtlasStripOptions *opts,
                                 double *x);

void AtlasStrip_measure(const AtlasStripProblem *problem,
                        const AtlasStripOptions *opts,
                        const double *x,
                        AtlasStripMetrics *metrics);

/* Evaluate a direct or barycentrically interpolated member observation. */
double AtlasStrip_member_value(const AtlasStripMember *member,
                               const double *u);

/*
 * Solve the alternating robust relaxation.  x contains [sample u, q] and its
 * sample entries must be feasible for the monotonicity bounds and exact
 * anchors; q entries are initialized internally.  out_membership receives the
 * final likelihood-weighted member memberships and may be NULL.
 */
int AtlasStrip_solve_robust(Arena_T arena,
                            const AtlasStripProblem *problem,
                            const AtlasStripOptions *strip_opts,
                            const AtlasStripRobustOptions *robust_opts,
                            const MonotoneQpOptions *qp_opts,
                            double *x,
                            double *out_membership,
                            AtlasStripRobustTrace *trace,
                            AtlasStripRobustStats *stats);

#endif
