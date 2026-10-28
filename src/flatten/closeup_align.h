#ifndef CLOSEUP_ALIGN_INCLUDED
#define CLOSEUP_ALIGN_INCLUDED

#include "slim_refine.h"

/* A proposal for ONE complete connected material component. The caller owns
 * extraction/scatter, provenance, coverage qualification and the assembled
 * all-material contact audit. Neither a successful proposal nor target_met
 * qualifies missing evidence, source support, seams or external contacts. */
enum {
    CLOSEUP_ALIGN_PROPOSED = 0,
    CLOSEUP_ALIGN_UNCHANGED = 1,
    CLOSEUP_ALIGN_NO_EVIDENCE = 2,
    CLOSEUP_ALIGN_INVALID = -1,
    CLOSEUP_ALIGN_INFEASIBLE = -2,
    CLOSEUP_ALIGN_SOLVER_FAILED = -3,
    CLOSEUP_ALIGN_TOPOLOGY = -4,
    CLOSEUP_ALIGN_MAX_TRIALS = 26
};

typedef struct CloseupAlignMetric {
    size_t faces, measured_faces, invalid_evidence_faces, inverted_faces;
    size_t protected_bad25_faces;
    double area, measured_area, confidence_mass;
    double fiber_rms_degrees; /* NAN without usable evidence, never zero */
    double sd_mean;          /* physical-area mean SD INCLUDING its minimum 4 */
    double direction_energy; /* sum area*relative_evidence*joint cross-axis
                               squared lengths; no trial direction multiplier */
    double min_det, min_sigma, max_sigma;
    double bad10_area, bad25_area;
    double protected_area, protected_bad10_area, protected_bad25_area;
} CloseupAlignMetric;

typedef struct CloseupAlignOpts {
    /* Required linear_solve/context are reused by both projection and SLIM.
     * iterations/local solve controls are honored. The module uses actual XYZ
     * rest metrics, free boundaries and no packing or hard vertex pins.
     * guide_uv, metric_uv, scaffold_uv and fixed_vertices must be NULL.
     * directions and direction_weight here are ignored: evidence is an explicit
     * immutable argument and the bounded weight policy belongs to this module. */
    SlimRefineOpts slim;
    int field_initialize;       /* optional positive-scale Poisson seed, default 1 */
    double target_rms_degrees;  /* shared close-up target, default 10; no coverage claim */
    int local_global_iterations; /* optional angular refinement; 0 keeps legacy control */
    const double *local_global_bad10_limits; /* optional absolute all/protected area caps; never above 5% */
    /* Optional fixed-matrix lifetime for linear_solve, using linear_context.
     * Both hooks or neither. A successful begin (0) is paired with exactly one
     * end, including on solve failure or exception; a failed/raising begin
     * retains responsibility for its own partial state. The scope covers only
     * the eight projection solves, never nonlinear SLIM. end must not raise.
     * Cached factors must own storage independent of the borrowed arena,
     * whose per-solve scratch is restored after each callback. */
    int (*projection_begin)(void *context);
    void (*projection_end)(void *context);
    /* Optional feasibility filter on the exact restored output coordinates.
     * Return >0 to admit, <=0 to refuse; refusal continues the bounded search.
     * The callback must not modify UV or inputs. Persistent query storage must
     * be independent of the borrowed scratch arena. NULL retains local checks
     * only; the caller still certifies the initial and final assembled surface. */
    void *candidate_context;
    int (*candidate_guard)(void *, size_t, const double *);
} CloseupAlignOpts;

typedef struct CloseupAlignTrial {
    int seed;                  /* 0 old seed, 1 field seed, 2 rigid turn, 3 local/global */
    double direction_weight;
    int solver_status, feasible, scale_nonregression, selected;
    int candidate_guard_pass; /* true without an external hook */
    size_t iterations_run, accepted_steps;
    CloseupAlignMetric metric;
} CloseupAlignTrial;

typedef struct CloseupAlignStats {
    int status, input_feasible, input_manifold_disk;
    int input_topology_supported;
    size_t boundary_loops;
    int field_attempted, field_status, field_feasible;
    int selected_seed, target_met, scale_nonregression;
    size_t field_faces, trials;
    double selected_direction_weight;
    int local_global_attempted,local_global_status,local_global_termination,local_global_loss_converged;
    double local_global_initial_loss,local_global_final_loss;
    CloseupAlignMetric before, field, after;
    CloseupAlignTrial trial[CLOSEUP_ALIGN_MAX_TRIALS];
} CloseupAlignStats;

void CloseupAlign_defaults(CloseupAlignOpts *opts);

/* Read-only physical/Jacobian measurement. orientation is +1 or -1, the
 * expected UV handedness. Both unsigned arms use ONE complementary X/Y
 * assignment, weighted by physical area times the unmodified positive
 * relative evidence. protected_faces[nf] identifies original/DATA faces;
 * NULL measures the aggregate of all faces; it does not imply nonregression
 * on a particular source-only subset. This does not certify injectivity.
 * Return 0, or INVALID without changing *metric. Scratch is always restored;
 * Arena_Failed propagates after restoration. */
int CloseupAlign_measure(Arena_T arena, const float *xyz, size_t nv,
                         const int32_t *faces, size_t nf, const double *uv,
                         const SlimRefineDirection *directions,
                         const uint8_t *protected_faces, int orientation,
                         CloseupAlignMetric *metric);

/* Separate uv_out[nv*2] and stats are caller-owned; no input is written.
 * No orphans, disconnected components or mixed orientation seeds. The native
 * manifold audit, boundary cycles and Euler characteristic require an orientable
 * genus-zero surface with boundary; annuli and multiple holes are supported.
 * Pinches, inconsistent winding and unsupported topology return TOPOLOGY. The
 * normalized UV boundary must have exactly one positive outer cycle and negative
 * hole cycles, in addition to the native simple-boundary and positive-face
 * checks. No cuts are silently introduced to satisfy these conditions.
 * An optional
 * field seed is refused unless it passes the existing native SLIM face and
 * boundary checks. The old-seed SLIM control is always evaluated. Candidate
 * selection requires nonincreasing physical SD and 10/25% bad area, including
 * protected material, before favoring smaller measured joint angular error.
 * Positive field scale updates are bounded; no cuts, welds or topology edits.
 * The input handedness and physical-area UV centroid are preserved. A scalar
 * translation gauge in the linear system fixes only the nullspace, not any
 * material point: every actual boundary vertex can move in the returned map.
 * Return PROPOSED, UNCHANGED or NO_EVIDENCE with complete output; negative
 * status leaves uv_out unchanged. Exceptions also leave both outputs unchanged
 * and restore all borrowed-arena scratch before propagation. */
int CloseupAlign_propose(Arena_T arena, const float *xyz, size_t nv,
                         const int32_t *faces, size_t nf, const double *uv_in,
                         const SlimRefineDirection *directions,
                         const uint8_t *protected_faces,
                         const CloseupAlignOpts *opts, double *uv_out,
                         CloseupAlignStats *stats);

/* A separate rigid alignment pass on an already flattened component. It
 * minimizes the measured joint angular error over a quarter turn, then tests
 * at most CLOSEUP_ALIGN_RIGID_ATTEMPTS exact endpoints (the optimum and halves).
 * Every endpoint must pass the same native output, protected-metric and caller
 * guards. The fixed direction audit objective must also not increase. No
 * nonlinear/linear solve is called; all proposal API input/output contracts,
 * including handedness, centroid, evidence and failure atomicity, still apply.
 * Failure or no admissible improvement never changes the supplied input. */
int CloseupAlign_rotate(Arena_T arena, const float *xyz, size_t nv,
                        const int32_t *faces, size_t nf, const double *uv_in,
                        const SlimRefineDirection *directions,
                        const uint8_t *protected_faces,
                        const CloseupAlignOpts *opts, double *uv_out,
                        CloseupAlignStats *stats);

int CloseupAlign_selftest(void);
int CloseupAlign_rigid_selftest(void);

#endif
