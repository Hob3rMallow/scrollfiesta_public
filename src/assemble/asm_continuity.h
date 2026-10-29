#ifndef ASM_CONTINUITY_INCLUDED
#define ASM_CONTINUITY_INCLUDED
#include "asm_types.h"

enum {
    ASM_CONT_SOURCE = 1u,
    ASM_CONT_VERIFIED = 2u,
    /* Geometry verified before the assembly cut split its two groups.
     * A later pose must preserve this continuation even while SWITCHED. */
    ASM_CONT_CUT = 4u,
    /* An accepted clean2 seam has too few ordered samples for a source
     * certificate. Keep its incumbent placement bundle intact, without
     * promoting the unresolved seam to a verified material join. */
    ASM_CONT_INCUMBENT = 8u,
    /* An explicitly retained original-source repair obligation. This does
     * not assert CT material identity or clear the placement flags. */
    ASM_CONT_REQUIRED = 16u,
    ASM_CONT_RECONSTRUCTED = 32u,
    /* Explicit cut along an original mesh edge within one chart. The field
     * validates both boundary edges, opposite winding, immutable source IDs
     * and identical XYZ before certifying its two endpoint observations. */
    ASM_CONT_INTRACHART_CUT = 64u
};

/* A bounded sparse tangent constraint A*x >= rhs. Each physical seam
 * touches at most two charts with three pose coordinates apiece. */
typedef struct AsmContinuityHalfspace {
    int ids[6], count;
    double values[6], rhs;
} AsmContinuityHalfspace;
int AsmContinuity_project_step(Arena_T arena, int n, int nt,
    const int *rows, const int *cols, const double *values, const double *rhs,
    double *step, const AsmContinuityHalfspace *constraints, size_t count);

typedef struct AsmContinuityStats {
    size_t source_runs, verified_joins, unresolved_hypotheses;
    size_t unsupported_islands, unsupported_charts;
    double unsupported_area;
    size_t unexplained_breaks, metadata_only_joins, wrong_placements;
    int audit_complete, confetti_free;
} AsmContinuityStats;

typedef struct AsmContinuityTarget {
    int32_t chart;
    double u, v, x, y, weight;
    int32_t parent;
} AsmContinuityTarget;

/* Transport original 3-D displacements into boundary-vertex tangent frames.
 * Both seam classification and subsequent continuity use these same frames.
 * Missing frames are reported per observation; no source observations move
 * or disappear. Scratch is released before return. Invalid inputs return -1. */
int AsmContinuity_chart_gaps(Arena_T arena, const AsmChart *chart,
                             const int32_t *vertices, const double *displacements,
                             size_t count, float *gaps, uint8_t *valid);

/* Freeze original physical tangent/face maps before any local UV repair.
 * Scratch is bounded by the largest incident chart. */
int AsmContinuity_prepare(AsmRun *run);
/* Recover omitted per-chart mutual source-boundary hypotheses before any
 * placement deformation. Existing proposals and frozen witnesses are never
 * replaced; new hypotheses remain placement-only until separately accepted. */
int AsmContinuity_reconstruct(AsmRun *run);
/* Provisional clean2 bundles may move together, but their original trim
 * residuals may not worsen during registration, local repair or placement. */
int AsmContinuity_preserves_incumbent(const AsmRun *run);
/* Equivalent transaction check when only the charts in saved can change.
 * Checks every incumbent seam incident to either endpoint in that set. */
int AsmContinuity_preserves_incumbent_trial(const AsmRun *run, const AsmChart *saved, size_t n);
int AsmContinuity_target(const AsmRun *run, const AsmRelation *r,
                         const AsmCorr *c, int reverse, double target[2]);
/* Recomputes residuals from carried source trim edges and current UV/poses,
 * never r->tx/ty or a frozen pre-deformation gap direction. */
int AsmContinuity_measure(const AsmRun *run, const AsmRelation *r,
                         double *rms, size_t *support);

/* Symmetric residual of one ordered witness, in the current chart poses.
 * A relation's aggregate support does not certify every point on its run. */
double AsmContinuity_pair_error(const AsmRun *run, const AsmRelation *r, const AsmCorr *c);
/* Certifies only measured, parity-consistent seams. Hard vetoes persist. */
size_t AsmContinuity_confirm(AsmRun *run);
/* Trial solve only. Caller owns the pose checkpoint and commits after its
 * overlap/rim audit. Every chart outside ids is held exactly fixed. */
int AsmContinuity_close(AsmRun *run, const AsmChart *saved, size_t n, const uint8_t *selected,
                        const AsmContinuityTarget *targets, size_t n_targets);
/* Seam-only recovery when the robust fit sacrifices too many boundary
 * points. Same pose bounds; caller still certifies and audits the result. */
int AsmContinuity_close_boundaries(AsmRun *run, const AsmChart *saved, size_t n, const uint8_t *selected);
/* The saved checkpoint is constant throughout a line search. Its baseline
 * audit can be reused; every trial's current geometry must still be checked.
 * The caller supplies an empty cache for each new checkpoint. */
typedef struct AsmContinuityUvPair {
    int32_t chart_a, vertex_a, chart_b, vertex_b;
} AsmContinuityUvPair;

typedef struct AsmContinuityAuditCache {
    size_t baseline_bad, baseline_queries;
    int baseline_ready;
    uint8_t *baseline_violations;
    Arena_T scratch; /* caller's checkpoint scope; NULL uses the run arena */
    /* Optional caller-owned feedback, indexed by chart id. The geometry
     * callback clears it and marks both charts at a newly violated query. */
    uint8_t *blocked_charts;
    size_t blocked_charts_count;
    /* Optional point feedback, also indexed by chart id. Each non-NULL
     * entry holds nv bytes owned by the caller; NULL denotes fixed material.
     * Both sampled endpoints of an introduced conflict are marked. */
    uint8_t **blocked_vertices;
    size_t blocked_vertices_count;
    /* Optional bounded list of exact newly conflicting source-vertex pairs.
     * The callback resets pairs_count on every audit, including success. */
    AsmContinuityUvPair *pairs;
    size_t pairs_capacity, pairs_count;
    /* Optional unchanged checkpoint of a larger spatial neighbourhood.
     * The callback audits these charts, including fixed neighbours, while
     * the solver's saved array still identifies only its moving charts. */
    const AsmChart *guard_charts;
    size_t guard_count;
} AsmContinuityAuditCache;
typedef int (*AsmContinuityUvAudit)(AsmRun *run, const AsmChart *saved, size_t n, AsmContinuityAuditCache *cache);
/* Registration commits only after the required geometry audit approves
 * the full component, including its fixed chart. Failure restores all
 * original poses and metadata before placement can trust the component. */
size_t AsmContinuity_register_components(AsmRun *run, AsmContinuityUvAudit audit);
/* After local UV repair, retry nearby failing SOURCE seams with rigid
 * poses. Entry holds all chart poses from before ordinary registration,
 * indexed by chart id; its original movement limits still apply. */
size_t AsmContinuity_register_failed(AsmRun *run, const AsmChart *entry, AsmContinuityUvAudit audit);
/* Correct local boundary UV residuals before placement certification. The
 * source positions, faces and measured trim gaps remain fixed. */
size_t AsmContinuity_relax_uv(AsmRun *run, AsmContinuityUvAudit audit, double stress_band, const char *diag_dir);
/* Jointly polish placed source seams in bounded, overlapping chart groups.
 * Each chart stays within eight voxels of its checkpoint at entry to this
 * pass. The required callback audits the supplied spatial neighbourhood. */
size_t AsmContinuity_relax_placed_uv(AsmRun *run, AsmContinuityUvAudit audit, double stress_band, const char *diag_dir);
/* Only saved's charts may change; certified neighbours stay fixed. UV is
 * committed into the caller's buffers, which must belong to its trial. */
int AsmContinuity_relax_boundary_trial(AsmRun *run, const AsmChart *saved, size_t n,
                                       const uint8_t *selected, AsmContinuityUvAudit audit, double stress_band);
/* Diagnostic: stop at the first improving direction that passes every
 * original trial check, before searching additional directions. */
int AsmContinuity_relax_boundary_step(AsmRun *run, const AsmChart *saved, size_t n,
                                      const uint8_t *selected, AsmContinuityUvAudit audit, double stress_band);
size_t AsmContinuity_components(AsmRun *run);
void AsmContinuity_audit(const AsmRun *run, int32_t primary,
                         const int32_t *emitted_chart, const float *xyz, size_t nv,
                         const int32_t *faces, size_t nf,
                         const char *path, AsmContinuityStats *st);
int AsmContinuity_selftest(void);
#endif
