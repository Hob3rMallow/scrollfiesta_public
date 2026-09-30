#ifndef ASM_POSE_INCLUDED
#define ASM_POSE_INCLUDED

#include <stddef.h>
#include <stdint.h>

#include "asm_axis.h"
#include "asm_types.h"

/* ============================================================================
 * asm_pose.h -- rigid SE(2) pose graph over charts.
 *
 * One (x, y, phi) per chart in the layout; each seam relation is a relative
 * pose measurement with information weights. Translations are solved at
 * chart centroids and residuals at the fitted seam centres, so changing a
 * local UV origin cannot change the material layout or robust switches.
 * Output poses retain the original chart frames. Parity bits (charts mirrored
 * relative to each other) are resolved first by a 2-colouring of the
 * relation graph; a parity-inconsistent cycle is a wrong relation and is
 * dropped.  Rotations start from the anchored chordal (complex least
 * squares) solve, translations from the linear solve at those rotations;
 * then Gauss-Newton with IRLS (Cauchy) and Levenberg damping refines all
 * three, and relations whose robust weight collapses are switched off.
 *
 * Linear algebra: TAUCS supernodal Cholesky through sparse_solve.h.  The
 * solver is process-global, so this module is never called from inside a
 * parallel region.
 * ==========================================================================*/

typedef struct AsmPoseOpts {
    int    gn_iters;      /* Gauss-Newton / IRLS iterations */
    double cauchy_c;      /* robust scale in normalized-residual units */
    double switch_off;    /* relations below this robust weight are switched off */
    double damping;       /* Levenberg damping added to the diagonal */
    int    verbose_prior; /* print the prior's chart count (the tool sets it for the pose stage only) */
    double axis_dir[3];   /* mean unit direction of the scroll axis in (z,y,x) over the box (reporting;
                           * the prior itself follows `axis` per chart when one is set) */
    const AsmAxis *axis;  /* the axis polyline: each chart's axial coordinate is the arc length at the
                           * foot of its centroid plus the local tangent's offset (NULL: world z) */
} AsmPoseOpts;

typedef struct AsmPoseStats {
    size_t charts_active, rels_active;
    size_t components;
    size_t parity_conflicts;    /* relations dropped by the parity colouring */
    size_t switched_off;
    size_t z_prior_charts;      /* charts carrying the world-z prior */
    size_t frames_from_anchor, frames_from_mean, frames_none;   /* components whose axial frame direction came from the
                                                                 * anchor chart's gradient / the charts' circular mean
                                                                 * (AP_ZPRIOR_FRAME) / nowhere (no chart under the prior) */
    int    largest_frame;       /* the largest component's frame source: 1 anchor, 2 mean, 0 none */
    size_t vprior_charts;       /* charts under the weak absolute axial prior (v3) */
    double vprior_fit[3];       /* the largest component's v = a + b s + c u (frame units) */
    double vplane_p50_before, vplane_p90_before, vplane_p50_after, vplane_p90_after;   /* |V - plane| over the largest component, vox */
    size_t components_after_switch;   /* components over the relations that survived the robust switching */
    size_t largest_component_charts;
    double largest_component_area;
    size_t degree_hist[6];      /* active charts by active-relation degree 0..4, 5+ */
    double degree_area[6];      /* their 3-D area */
    double resid_p50, resid_p95, resid_max;   /* normalized (sigma units) */
    double xy_p50, xy_p95;                    /* translation residual, vox */
    int    iters;
    double sec;
} AsmPoseStats;

void AsmPose_default_opts(AsmPoseOpts *o);

/* Solves poses for every chart in the layout; sets chart->pose_*, ->placed,
 * ->component and ASM_CHART_MIRROR; sets rel->robust_w / ->residual and
 * flags ASM_REL_SWITCHED / ASM_REL_DROPPED.  Returns 0 on success. */
int AsmPose_solve(AsmRun *run, const AsmPoseOpts *o, AsmPoseStats *st);

/* POST-PLACEMENT REFINE (2026-09-10): a Gauss-Newton / IRLS re-solve over an explicit relation
 * list from the charts' CURRENT poses (warm start), with the `fixed` charts pinned (plus the
 * largest chart of every connected component without a fixed chart, as its gauge).  Touches no
 * component / placed / mirror state; relations whose parity contradicts the charts' mirror bits
 * are skipped; the listed relations get robust_w / residual and ASM_REL_SWITCHED as the solve
 * does; poses are written for the free charts only. */
typedef struct AsmPoseRefineStats {
    size_t charts_active, charts_free, charts_pinned;
    size_t rels, rels_parity_refused, switched;
    int    iters;
    double resid_p50, resid_p95, xy_p50, xy_p95;
    double sec;
} AsmPoseRefineStats;
int AsmPose_refine(AsmRun *run, const AsmPoseOpts *o, const int32_t *rel_idx, size_t n_rel, const uint8_t *fixed, AsmPoseRefineStats *st);

int AsmPose_selftest(void);

#endif
