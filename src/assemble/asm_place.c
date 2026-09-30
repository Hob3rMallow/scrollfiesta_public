/* asm_place.c -- placement of disconnected components in the scroll frame.
 *
 * The sheet is a straight strip: u runs along the winding with the scroll's
 * centre at the LEFT (u increases outward), v along the roll axis (height).
 * That fixes every component's ORIENTATION on its own:
 *   - rotation: the direction of increasing world z in the component's
 *     intrinsic frame becomes +v (exact for closed turns; an arc around a
 *     tilted axis is biased by the tilt x its angular position, so arcs
 *     re-take their rotation from the anchors);
 *   - mirror: the outward turn vector (displacement to the layer neighbour on
 *     the side away from the bending direction) must point to +u; components
 *     without their own layer pairs follow the scroll's chirality, measured
 *     on the components that have both a bending and a turn vote.
 * Placement is then a translation found by VOTING: every layer hit between a
 * placed chart and the component proposes the component's offset for each
 * layer count k = 1..K (k turns along the sheet, with the measured growth of
 * the turn length per layer); the offsets are clustered, and the clusters
 * are tried in order of weight -- fit, gate, tentative placement, audit --
 * so a hit that passed through a hole in the intermediate layer (k = 2 read
 * as k = 1) cannot put the component on top of material that is already
 * there. */
#include "asm_place.h"
#include "asm_metric.h"
#include "asm_continuity.h"
#include "asm_contacts.h"
#include "asm_field.h"
#include "asm_store.h"
#include "asm_audit.h"
#include "asm_report.h"
#include "../common/uv_guard.h"
#include "../flatten/sparse_solve.h"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../common/ves_omp.h"

#include "../common/pipeline_constants.h"
#include "../common/union_find.h"
#include "../common/ves_platform.h"

#define APL_PI 3.14159265358979323846
#define APL_MAX_CLUSTERS 4
/* V FROM THE AXIS.  The papyrus is perpendicular to the local axis, so a sheet's v is arc length
 * along the axis: the primary's own relation v = a + b s + c u (s = the chart's axial coordinate,
 * u for a slightly rotated frame) bounds every placed piece's v within the
 * measured scatter of that model; anchors choose the pose inside this band.
 * A piece whose charts still scatter about that line by more than APL_V_AXIS_TOL
 * (rms, vox) is not that sheet and is refused.  Measured need (4x21x21, 2026-09-09): with 816
 * placed components the layout's v span was 14,469 vox for a 512-vox slab.  0 = off (A/B). */
#define APL_V_FROM_AXIS 1
#define APL_V_AXIS_TOL 30.0
#define APL_V_AXIS_MAX_SHIFT 100.0   /* a piece the anchors put further than this from the sheet's v(s) has bad anchors: refused */
#define APL_V_MODEL_MAX_CORR 0.7     /* the u term is dropped when s and u are this collinear over the primary (the centre read c = -0.22 and moved pieces 1,389 vox) */
#define APL_V_MODEL_B_MIN 0.95       /* the pinning AND the radius turn model are armed only when the primary's v really is arc
                                      * length along the axis: b in [0.95, 1.05].  Measured 2026-09-09: familiar 1.002, 10x 1.008,
                                      * 4x21x21 0.999; the centre's CT-prior table (33.5 deg off +z) reads 0.9135 and, armed at
                                      * [0.85, 1.15], cost 9 points of coverage (72.3 -> 63.1): a table that is not the sheet's
                                      * axis must neither pin nor set hop lengths */
#define APL_V_MODEL_B_MAX 1.05
#define APL_V_MODEL_MAX_C 0.10       /* and when it exceeds a 6-degree frame rotation (the rotation prior holds the frame to ~4 deg): the centre read -0.22 with s-u correlation 0.17 and placed 1 component */
/* TURN LENGTH FROM THE RADIUS.  One turn of the spiral at radius r is 2 pi r; a hop of k wraps
 * from a chart at radius r is 2 pi k (r +- k d / 2) with d the local layer distance.  Armed by an
 * axis point (verdict.umb_y/umb_x or the axis table); the measured/estimated turn model stays
 * the fallback.  The measured growth per wrap read 2.1 vox on the 4x21x21 against 2 pi d = 69,
 * so every hop away from the primary's radius was hundreds of vox wrong.  0 = off (A/B). */
/* SWEEP-LOOP COST (2026-09-10: the 4x21x21 placed 5,680 sweeps in its 900-s budget and left 6.6% of
 * its area never attempted; every sweep re-measured the turns over all layer hits, re-fitted the
 * quadric of every placed chart in apl_stretch_turn, and insertion-sorted the anchors).  Each lever is
 * EXACT (bit-identical placement) and stays behind its constant. */
#define APL_TURN_REMEASURE 1     /* 1 = measure the turns and growth at sweep 0 and only after a sweep that placed something (a sweep that places nothing leaves both unchanged: apl_turn_pair_ok needs both charts placed); 0 = every sweep */
/* 1 = re-measure the base turn t0 whenever the placed material changed, accepting only a
 * MEASURED estimate resting on more turns than the incumbent; 0 = t0 is whatever sweep 0 could
 * see, for the whole run.  The 4x5x5 measures t0 over 568 turns at sweep 0 and never needs
 * this; the 21x5x5 starts with NONE ("t0 curvature over 0 turns", a 3,612-vox guess frozen for
 * 4,148 sweeps) and t0_measured 0 also disables the layer model in apl_collect_anchors.
 * A/B pending (2026-09-17). */
#define APL_T0_REMEASURE 0
/* A hop cluster whose ESTIMATED turns are known no better than this relative error is refused
 * outright instead of fitted.  Every hop gate is a fraction of |T| (the fit gate 5%, the
 * cluster window 15%, the consensus scale 5%), which is right when the turn is known and
 * exactly wrong when it is a guess: the less we know it the wider the gates open.  On the
 * 21x5x5 that puts component 18 (1.18M vox^2) through 42 attempts at rms 39.6 against a gate
 * of 342.4 -- the fit is not what stopped it -- and the deferral ledger then blames the
 * consensus or the ambiguity rather than the turn.  0 = off, the shipped behaviour.
 * The statistic is the WORST relative error among the candidate's estimated turns, so read any
 * bound against what a healthy pile shows: on the 4x5x5, whose 73,971 radius anchors need no
 * estimate at all, the candidates that do fall through to the stretch carry relative errors of
 * 137 to 520 -- the stretch is noise there and the radius model is why it does not matter.
 * A/B pending (2026-09-17); measure the deferred area and the lineage count, not coverage. */
#define APL_HOP_TURN_MAX_REL_ERR 0.0
#define APL_CURV_CACHE 1         /* 1 = the placed side's per-chart curvature in apl_stretch_turn is computed once per placed chart (a placed pose is frozen inside the loop); 0 = every call */
#define APL_CLUSTER_QSORT 1      /* 1 = qsort the anchor offsets (ties by index: the stable insertion sort's order); 0 = the O(n^2) insertion sort */
#define APL_RIM_SQUARED_CULL 1   /* reject distant grid samples before hypot; retain the original test near the reach boundary */
#define APL_PROFILE_EVERY 500    /* sweeps between [place] profile lines (0 = only the final line) */
#define APL_TICK(field) do { double tn_ = ves_clock_sec(); st->field += tn_ - tq; tq = tn_; } while (0)
/* SEAMS FIRST (2026-09-10).  Seam anchors (k = 0: the correspondences of a WEAK / SHORT relation to a
 * placed chart) are exact; hop anchors carry the turn model's error, 100+ vox on the outer wraps.  The
 * offset windows merged them (wu = 0.15 x the turn = 300-750 vox out there) and the weighted fit was
 * the hop fit: on the 10x10x10 dev75 the seam anchors were outvoted 28:1 and 214 seam-anchored
 * components (9.1% of the area) were REVERTED for landing a layer off.  With this on, the seam anchors
 * are clustered among themselves at a tight window, ranked first and fitted with their own gate; the
 * hop clusters follow.  0 = the merged clusters (dev75). */
#define APL_SEAM_FIRST 1
#define APL_SEAM_FRONTIER_FIRST 0 /* Global priority still regresses after recovery fixes (dev155): coverage92.4->85.3, stacked UV1162->4862. Retain inside-out area order. */
#define APL_SWITCHED_SEAM_ANCHORS 1 /* robust pose rejection can separate components without invalidating their measured boundary correspondences */
#define APL_SEAM_CLUSTER_VOX 24.0    /* the seam clusters' window: two seams to different neighbours that disagree stay two hypotheses */
#define APL_MAX_SEAM_CLUSTERS 2
#define APL_MAX_CLUSTERS_TOTAL (APL_MAX_CLUSTERS + 2 * APL_MAX_SEAM_CLUSTERS)
#define APL_SEAM_MIN_ANCHORS 6       /* a seam cluster's minimum members (a SHORT seam carries 6-15 correspondences, under min_anchors 12) */
#define APL_SEAM_GATE_MIN 8.0        /* vox: a seam cluster's rms gate is max(this, APL_SEAM_GATE_RMS_X x its seams' rigid rms), never the turn fraction */
#define APL_SEAM_GATE_RMS_X 4.0
#define APL_SEAM_V_PIN 0             /* 1 = the v-axis pin also SHIFTS seam placements; 0 = it only refuses them (the neighbour they seam to is pinned already) */
/* DEFERRAL (2026-09-10): a small component placed by a HOP cluster is accepted at a gate of 0.05 x the
 * turn (200-600 vox on the outer wraps) whenever the co-location and rim audits find nothing placed
 * within reach -- a piece 100 vox off in empty space passes, and that is the confetti's slits.  With
 * this on such a placement is DEFERRED (status deferred, comp_placed = 2): the retry-on-new-evidence
 * rule re-tries it when a neighbour lands, by then with a rim to judge it or a seam to snap to.  Seam
 * placements are never deferred.  0 = accept as before (dev75). */
#define APL_SMALL_HOP_NEEDS_RIM 1
#define APL_SMALL_AREA_VOX2 4.0e4   /* vox^2: about two cube-sized charts */
#define APL_DEFER_RETRY_RIM 1      /* retry silent placements after the rim grid grows, even when their anchor count stays unchanged */
#define APL_DEFER_ONLY_WITH_SEAMS 0 /* 1 = defer only a piece that holds a seam (accepted-dropped, switched, WEAK, SHORT) to some other chart; a seamless piece keeps the hop placement.  MEASURED 2026-09-10 (dev89 vs dev88): 10x coverage 71.3 -> 76.1 but lineages 133 -> 357 and stacked (u,z) 8.0k -> 28.0k; familiar 86.5 -> 87.6, lineages 22 -> 36; 21x5x5 69.7 -> 74.1, lineages 132 -> 329.  The seamless pieces ARE the floating confetti; OFF by the contiguity-first rule. 0 = defer every small hop-only piece whose audits are silent */
#define APL_DEFER_LAST_CALL 0       /* 1 = the deferred components get one more pass at the end with the deferral off and the STRICT hop gate (the size gate alone, no turn fraction).  MEASURED 2026-09-10 (10x10x10 dev87 vs dev83): coverage 71.3 -> 75.7 but lineages 133 -> 391, stacked (u,z) 8.0k -> 16.3k, within-wrap gap p90 130 -> 472: a hop placement whose anchors agree to 10 vox can still be a turn off, and the strict gate cannot tell.  OFF by the contiguity-first rule; 0 = they stay deferred */
#define APL_OUTLINE_CROP_VOX 8.0 /* vertices this close to the box's top / bottom z trace the crop planes (the outline score) */
#define APL_TURN_FROM_RADIUS 1   /* A/B 2026-09-09: familiar dev66 81% placed / largest 17% (radius) vs dev67 84% / 16% (measured); 10x dev65 61.9% placed (radius) vs dev67 48.4% (measured) */
#define APL_TURN_MIN 150.0    /* one wrap is at least 2 pi x 25 vox: a shorter 'turn' is a fold inside a component */
#define APL_TURN_MAX 12500.0  /* and at most 2 pi x 2000 vox for these scrolls */
/* METRIC RE-SOLVE: every APL_METRIC_EVERY newly placed components (and at the end) the u
 * offset of every placed component is refit to ALL layer hits between placed charts under
 * one turn model T(w) = t0 + g w fitted to the same hits, and the audit grids are rebuilt.
 * Motive (2026-09-09): measured turns of neighbouring placed charts disagree by 17-92% of a
 * turn on the 10x10x10 and the familiar-in-10x agreement is 192 vox p50 -- the greedy chain
 * accumulates each hop's estimate error and no per-anchor estimate can undo that. */
#define APL_METRIC_RESOLVE 0   /* v2 (sign-agreeing pairs only, 0.05 turn per pass, applied only when the robust residual falls)
                                * measured 2026-09-09: familiar neutral (92.2%, 8/8, residual 20 -> 13 vox), 10x10x10 placed
                                * 68.0 -> 59.2%: its pair residual is 250-270 vox p50 (a quarter turn) and rigid per-component
                                * shifts cut it by 5-15 vox a pass while the moved components then fail the audits.  The
                                * inconsistency is inside the placed components' own geometry vs the hits, not between them. */
#define APL_METRIC_EVERY 8
/* The metric re-solve's TURN MODEL alone -- one t0 and one growth fitted by IRLS to every
 * placed pair's own u displacement -- adopted without its component shifts, which stay
 * measured dead above.  The turn at a given radius is a property of the pile, not of the
 * candidate in front of us; today each candidate re-derives it from its own 20-79 layer hits
 * and gets anything from 645 to 9,105 vox (21x5x5, 2026-09-17).  Adopting the model also turns
 * t0_measured on, which is what re-enables the layer model in apl_collect_anchors.  The 4x5x5
 * already measures t0 over 568 real turns, so it has nothing to gain.  A/B pending. */
#define APL_TURN_MODEL_POOL 0
#define APL_TURN_MODEL_MIN_PAIRS 50     /* sign-agreeing placed pairs before the pooled fit may replace t0 */
#define APL_TURN_MODEL_MAX_RESID 0.15   /* and its robust residual must stay under this fraction of a turn */
/* COORDINATE CALIBRATION: a placed pair's hop count and sign are read from its actual u
 * displacement in the layout (round(|du| / T)), and the layer coordinate is re-solved with those
 * edges pinned (weight x4).  Measured motive 2026-09-09: on the 10x10x10 the hits INSIDE the
 * primary are off by a whole turn at the median against every-edge-is-one-hop (the observed
 * displacement is ~0 or ~2 turns for half of them: two-wrap hits across holes), while the
 * primary's own u is right; on the familiar the same residual is 4% of a turn. */
/* the hop count of every layer edge is re-estimated from the coordinate after each robust round (1..3) */
/* the wrap count of a layer pair is the probe's measured crossing order (1 = adjacent) */
#define APL_PAIR_K_FROM_PROBE 1
#define APL_LAYER_KINT 0   /* measured 2026-09-09: familiar neutral (91.6%, 8/8), 10x10x10 60.7 -> 51.9%; the graph's p95
                            * residual stays 2.6 layers either way, so the inconsistency is not the hop count but
                            * the votes themselves (direction and partner) at the core.  Off by the champion rule. */
#define APL_COORD_CALIBRATE 0   /* measured 2026-09-09: familiar 91.7 -> 88.9% coverage (8/8), 10x10x10 60.7 -> 48.6%; the
                                 * pinned placed pairs leave the graph's p95 residual at 2.6 layers: the unplaced
                                 * frontier's two-wrap edges stay unresolved and the calibrated core shifts the
                                 * frontier's coordinates the wrong way.  Off by the champion rule. */
/* TURN FIELD: the hop length for a placed chart without a measured turn is the robust median of
 * the observed one-wrap displacements of ALL placed pairs near it in (u, layer) -- the real turn
 * varies with azimuth by +-20% (oval sections) and every azimuth-blind model (t0 + g w) leaves a
 * quarter-turn residual on the 10x10x10 (measured 250-270 vox p50, 2026-09-09).  Binned in u
 * (APL_TFIELD_DU vox) and layer (APL_TFIELD_DW), rebuilt every APL_METRIC_EVERY placements. */
#define APL_TURN_FIELD 0   /* measured 2026-09-09: familiar +0.1 coverage (8/8), 10x10x10 60.7 -> 59.0%: the local
                            * observed turn does not beat the layer model where the layout itself is a quarter
                            * turn inconsistent.  Off by the champion rule. */
#define APL_TFIELD_DU 250.0
#define APL_TFIELD_DW 0.6
#define APL_TFIELD_MIN_OBS 5
#define APL_METRIC_MAX_SHIFT 0.05   /* per pass, fraction of the mean turn (v1 at 0.30 hit the cap every pass) */
#define APL_TURN_KHOP 0   /* measured 2026-09-09: 10x10x10 +0.9 coverage, 21x5x5 +3.5, familiar -2.4 (loses 8/8), centre -7.7:
                            * the primary's own turns were never two-wrap (t0 1952 -> 1949), the division only
                            * halves genuine later measurements on the small piles.  Off by the champion rule. */
#define APL_LAYER_UNKNOWN (-1e30)
#define APL_BORROW_TURNS 0   /* v4 (layout-window donors, median of up to 8) DEAD too: familiar 92 -> 73%, centre coverage 66 -> 22%.  v1 (any wrap), v2 (same layer coordinate, nearest along u) and v3 (v2 with the side mapped
                                * outward-to-outward) all measured DEAD 2026-09-09: familiar 92 -> 32-37%, centre 74 -> 19-57%;
                                * a borrowed turn wrong by one wrap poisons every anchor from that chart */
#define APL_BORROW_REACH 0.35

void AsmPlace_default_opts(AsmPlaceOpts *o)
{
    o->min_anchors = ASM_PLACE_MIN_ANCHORS;
    o->max_rms = ASM_PLACE_MAX_RMS;
    o->rel_rms = 0.04;
    o->max_layers = 3;
    o->layer_tol = 0.35;
    o->min_turn_hits = 3;
    o->max_contra_ratio = 0.25;
    o->max_rot_refine = 0.6;
    o->max_rim_ratio = 0.3;
    o->crosswrap_anchor_veto = ASM_WRAP_VETO_ANCHOR;
    o->stress_band = ASM_STRESS_BAND;
    o->diag_dir = NULL;
    o->axis = NULL;
}

/* ---- pose helpers -------------------------------------------------------------- */

static void apl_pose(const AsmChart *c, double u, double v, double *gx, double *gy)
{
    if (c->flags & ASM_CHART_MIRROR) u = -u;
    double ct = cos(c->pose_theta), st = sin(c->pose_theta);
    *gx = ct * u - st * v + c->pose_x;
    *gy = st * u + ct * v + c->pose_y;
}

static void apl_vertex_point(const AsmChart *c, int32_t v, double *gx, double *gy)
{
    double p[2]; AsmChart_point(c,(size_t)v,p); *gx=p[0]; *gy=p[1];
}

static void apl_face_point(const AsmChart *c, int32_t f, double l0, double l1, double *gx, double *gy)
{
    const int32_t *fv = &c->faces[(size_t)f*3];
    double l2 = 1.0 - l0 - l1;
    if (c->placed_uv) {
        *gx=l0*c->placed_uv[2*(size_t)fv[0]]+l1*c->placed_uv[2*(size_t)fv[1]]+l2*c->placed_uv[2*(size_t)fv[2]];
        *gy=l0*c->placed_uv[2*(size_t)fv[0]+1]+l1*c->placed_uv[2*(size_t)fv[1]+1]+l2*c->placed_uv[2*(size_t)fv[2]+1];
        return;
    }
    double u = l0 * c->uv[(size_t)fv[0]*2] + l1 * c->uv[(size_t)fv[1]*2] + l2 * c->uv[(size_t)fv[2]*2];
    double v = l0 * c->uv[(size_t)fv[0]*2+1] + l1 * c->uv[(size_t)fv[1]*2+1] + l2 * c->uv[(size_t)fv[2]*2+1];
    apl_pose(c, u, v, gx, gy);
}

/* compose a component transform (theta, t, mirror) onto a chart pose */
static void apl_compose(AsmChart *c, double theta, double tx, double ty, int mirror)
{
    double x = c->pose_x, y = c->pose_y;
    if (mirror) x = -x;
    double ct = cos(theta), st = sin(theta);
    c->pose_x = ct * x - st * y + tx;
    c->pose_y = st * x + ct * y + ty;
    c->pose_theta = theta + (mirror ? -c->pose_theta : c->pose_theta);
    if (mirror) c->flags ^= ASM_CHART_MIRROR;
}

/* Compare the material positions, not the pose parameters. Angle wrapping,
 * mirror parity and the extent of the chart all matter to equivalence. */
static int apl_same_pose(const AsmChart *a, const AsmChart *b, double tolerance)
{
    if (a->id != b->id || a->nv != b->nv || !a->uv || !b->uv || !a->nv) return 0;
    double ct_a = cos(a->pose_theta), sn_a = sin(a->pose_theta);
    double ct_b = cos(b->pose_theta), sn_b = sin(b->pose_theta);
    double ma = a->flags & ASM_CHART_MIRROR ? -1.0 : 1.0;
    double mb = b->flags & ASM_CHART_MIRROR ? -1.0 : 1.0;
    for (size_t v = 0; v < a->nv; v++) {
        double au = ma*a->uv[2*v], av = a->uv[2*v+1];
        double bu = mb*b->uv[2*v], bv = b->uv[2*v+1];
        double dx = ct_a*au-sn_a*av+a->pose_x - (ct_b*bu-sn_b*bv+b->pose_x);
        double dy = sn_a*au+ct_a*av+a->pose_y - (sn_b*bu+ct_b*bv+b->pose_y);
        if (!isfinite(dx) || !isfinite(dy) || dx*dx+dy*dy > tolerance*tolerance) return 0;
    }
    return 1;
}

static int apl_cmp_double(const void *x, const void *y)
{
    double a = *(const double *)x, b = *(const double *)y;
    return (a > b) - (a < b);
}

/* ---- orientation from the scroll structure ------------------------------------- */

static int apl_solve3(double A[9], double b[3], double x[3])
{
    for (int col = 0; col < 3; col++) {
        int piv = col;
        for (int r = col + 1; r < 3; r++) if (fabs(A[r*3+col]) > fabs(A[piv*3+col])) piv = r;
        for (int k = 0; k < 3; k++) { double t = A[col*3+k]; A[col*3+k] = A[piv*3+k]; A[piv*3+k] = t; }
        { double t = b[col]; b[col] = b[piv]; b[piv] = t; }
        if (fabs(A[col*3+col]) < 1e-12) return -1;
        for (int r = 0; r < 3; r++) {
            if (r == col) continue;
            double f = A[r*3+col] / A[col*3+col];
            for (int k = 0; k < 3; k++) A[r*3+k] -= f * A[col*3+k];
            b[r] -= f * b[col];
        }
    }
    for (int k = 0; k < 3; k++) x[k] = b[k] / A[k*3+k];
    return 0;
}

/* Per-chart quadratic fit x(u,v) = c0 + J (u,v) + 1/2 (u,v)^T H (u,v) of the
 * 3-D position over the chart's CURRENT global uv.  Writes the tangents
 * t_u, t_v (3-vectors, world z,y,x order) and the u-curvature vector h_uu. */
static int apl_chart_quadric(const AsmChart *c, double tu[3], double tv[3], double huu[3])
{
    double M[36] = { 0 }, R[3][6] = { { 0 } };
    double cu = 0.0, cv = 0.0;
    size_t n = c->nv;
    for (size_t i = 0; i < n; i++) { double gx, gy; apl_vertex_point(c, (int32_t)i, &gx, &gy); cu += gx; cv += gy; }
    cu /= (double)n; cv /= (double)n;
    for (size_t i = 0; i < n; i++) {
        double gx, gy;
        apl_vertex_point(c, (int32_t)i, &gx, &gy);
        double u = gx - cu, v = gy - cv;
        double phi[6] = { 1.0, u, v, u*u, u*v, v*v };
        for (int r = 0; r < 6; r++) {
            for (int k = 0; k < 6; k++) M[r*6+k] += phi[r] * phi[k];
            for (int d = 0; d < 3; d++) R[d][r] += phi[r] * c->xyz[i*3 + (size_t)d];
        }
    }
    for (int col = 0; col < 6; col++) {
        int piv = col;
        for (int r = col + 1; r < 6; r++) if (fabs(M[r*6+col]) > fabs(M[piv*6+col])) piv = r;
        for (int k = 0; k < 6; k++) { double t = M[col*6+k]; M[col*6+k] = M[piv*6+k]; M[piv*6+k] = t; }
        for (int d = 0; d < 3; d++) { double t = R[d][col]; R[d][col] = R[d][piv]; R[d][piv] = t; }
        if (fabs(M[col*6+col]) < 1e-9) return -1;
        for (int r = 0; r < 6; r++) {
            if (r == col) continue;
            double f = M[r*6+col] / M[col*6+col];
            for (int k = 0; k < 6; k++) M[r*6+k] -= f * M[col*6+k];
            for (int d = 0; d < 3; d++) R[d][r] -= f * R[d][col];
        }
    }
    for (int d = 0; d < 3; d++) {
        tu[d] = R[d][1] / M[1*6+1];
        tv[d] = R[d][2] / M[2*6+2];
        huu[d] = 2.0 * R[d][3] / M[3*6+3];
    }
    return 0;
}

typedef struct AplTurn { double tx, ty; int n; } AplTurn;

typedef struct AplOrient {
    double theta;     /* rotation applied so that +v = increasing world z */
    double grad;      /* |dz/d(u,v)| of the component, 1 = perfectly level */
    double bend;      /* area x curvature weighted (viewer normal . u-curvature): > 0 = viewer side inward */
    double turn_u;    /* area weighted cosine of the outward turn vector against +u: > 0 = u increases outward */
    int    has_turn;  /* turn_u measured from the component's own layer pairs */
    int    mirrored;
} AplOrient;

/* Rotate one component so that +v is the direction of increasing world z. */
static void apl_orient_measure(const AsmRun *run, int32_t comp, AplOrient *o)
{
    size_t nc = run->n_charts;
    double sxx = 0, sxy = 0, syy = 0, sx = 0, sy = 0, sxz = 0, syz = 0, sz = 0, n = 0;
    for (size_t i = 0; i < nc; i++) {
        const AsmChart *c = &run->charts[i];
        if (!AsmChart_in_layout(c) || c->component != comp) continue;
        for (size_t k = 0; k < c->nv; k += 2) {
            double gx, gy;
            apl_vertex_point(c, (int32_t)k, &gx, &gy);
            double z = c->xyz[k*3];
            sxx += gx*gx; sxy += gx*gy; syy += gy*gy; sx += gx; sy += gy; sxz += gx*z; syz += gy*z; sz += z; n += 1;
        }
    }
    o->theta = 0.0; o->grad = 0.0;
    if (n < 6) return;
    double A[9] = { sxx, sxy, sx, sxy, syy, sy, sx, sy, n }, b[3] = { sxz, syz, sz }, x[3];
    if (apl_solve3(A, b, x) != 0) return;
    double ga = x[0], gb = x[1], gl = hypot(ga, gb);
    o->grad = gl;
    if (gl > 0.05) o->theta = 0.5 * APL_PI - atan2(gb, ga);   /* z-gradient -> +v */
}

static void apl_orient_rotate(AsmRun *run, int32_t comp, AplOrient *o)
{
    apl_orient_measure(run,comp,o);
    size_t nc=run->n_charts;
    for (size_t i = 0; i < nc; i++) {
        AsmChart *c = &run->charts[i];
        if (!AsmChart_in_layout(c) || c->component != comp) continue;
        apl_compose(c, o->theta, 0.0, 0.0, 0);
    }
}

/* Two votes in the component's current frame.  bend: the inward direction
 * comes from the configured roll axis, falling back to u-curvature h_uu, and
 * its sign against the viewer normal t_u x t_v (right-handed in x,y,z) says
 * which side the viewer is on.  turn_u: the displacement to the OUTWARD
 * layer neighbour (the side opposite to inward along the chart's own normal)
 * is one turn further along the sheet, so its u component says whether u
 * increases outward -- the scroll's centre must end up at the LEFT.
 * Also records every chart's outward side (+1/-1 along its stored normal,
 * 0 = flat/unknown), which does not depend on the pose. */
/* Local dents can bend away from the roll axis. When an axis is supplied,
 * decide the normal's outward sign from the original material points. A
 * chart that straddles a fold or has mostly tangential normals gives no
 * reliable chart-wide vote and retains the curvature fallback. */
static int8_t apl_axis_inward(const AsmChart *c, const AsmAxis *axis, double inward[3])
{
    if (!axis || !c->xyz || !c->nrm || c->nv < 12) return 0;
    size_t stride=(c->nv+127)/128, samples=0, positive=0, negative=0;
    double sum[3]={0};
    for (size_t v=0; v<c->nv; v+=stride) {
        double p[3],n[3],foot[3],r;
        for (int d=0; d<3; d++) { p[d]=c->xyz[3*v+d]; n[d]=c->nrm[3*v+d]; }
        samples++;
        double nl=sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
        if (!(nl>1e-9) || !isfinite(nl)) continue;
        AsmAxis_project(axis,p,NULL,&r,NULL,foot);
        if (!(r>5) || !isfinite(r)) continue;
        double dot=0;
        for (int d=0; d<3; d++) dot+=n[d]*(p[d]-foot[d])/(nl*r);
        if (fabs(dot)<.25 || !isfinite(dot)) continue;
        if (dot>0) positive++; else negative++;
        for (int d=0; d<3; d++) sum[d]-=(p[d]-foot[d])/(r*r);
    }
    size_t known=positive+negative, majority=positive>negative ? positive : negative;
    if (known<6 || 4*known<3*samples || 10*majority<9*known) return 0;
    for (int d=0; d<3; d++) inward[d]=sum[d]/(double)known;
    return positive>negative ? (int8_t)1 : (int8_t)-1;
}

static void apl_component_out_sides(const AsmRun *run, int32_t comp, const AsmAxis *axis, int8_t *out_side)
{
    for (size_t i = 0; i < run->n_charts; i++) {
        const AsmChart *c = &run->charts[i];
        if (!AsmChart_in_layout(c) || c->component != comp || c->nv < 12) continue;
        double tu[3], tv[3], huu[3];
        int8_t axial=apl_axis_inward(c,axis,huu);
        if (axial) { out_side[c->id]=axial; continue; }
        if (apl_chart_quadric(c, tu, tv, huu) != 0) continue;
        double hl = sqrt(huu[0]*huu[0] + huu[1]*huu[1] + huu[2]*huu[2]);
        if (hl < 1e-6) continue;
        double mn[3] = { 0.0, 0.0, 0.0 };
        for (size_t k = 0; k < c->nv; k++) { mn[0] += c->nrm[k*3]; mn[1] += c->nrm[k*3+1]; mn[2] += c->nrm[k*3+2]; }
        double inward = mn[0]*huu[0] + mn[1]*huu[1] + mn[2]*huu[2];
        if (fabs(inward) < 1e-9) continue;
        out_side[c->id] = inward > 0.0 ? (int8_t)-1 : (int8_t)1;
    }
}

static void apl_orient_votes(const AsmRun *run, int32_t comp, const AsmAxis *axis,
                             const AplTurn *turn, AplOrient *o, const int8_t *out_side)
{
    size_t nc = run->n_charts;
    double bs = 0.0, bw = 0.0, ts = 0.0, tw = 0.0;
    for (size_t i = 0; i < nc; i++) {
        const AsmChart *c = &run->charts[i];
        if (!AsmChart_in_layout(c) || c->component != comp || c->nv < 12) continue;
        double tu[3], tv[3], huu[3];
        if (apl_chart_quadric(c, tu, tv, huu) != 0) continue;
        apl_axis_inward(c,axis,huu);
        double hl = sqrt(huu[0]*huu[0] + huu[1]*huu[1] + huu[2]*huu[2]);
        if (hl < 1e-6) continue;
        /* right-handed cross product in (x, y, z): world arrays are (z, y, x) */
        double ux = tu[2], uy = tu[1], uz = tu[0], vx = tv[2], vy = tv[1], vz = tv[0];
        double Nx = uy*vz - uz*vy, Ny = uz*vx - ux*vz, Nz = ux*vy - uy*vx;
        double Nl = sqrt(Nx*Nx + Ny*Ny + Nz*Nz);
        if (Nl > 1e-9) {
            double dot = (Nx*huu[2] + Ny*huu[1] + Nz*huu[0]) / (Nl * hl);
            bs += c->area3d * hl * dot; bw += c->area3d * hl;
        }
        int8_t os=out_side[c->id];
        if (!os) continue;
        const AplTurn *t = &turn[(size_t)c->id * 2 + (os > 0 ? 1 : 0)];
        if (t->n <= 0) continue;
        double len = hypot(t->tx, t->ty);
        if (len < 1e-6) continue;
        ts += c->area3d * t->tx / len; tw += c->area3d;
    }
    o->bend = bw > 0.0 ? bs / bw : 0.0;
    o->turn_u = tw > 0.0 ? ts / tw : 0.0;
    o->has_turn = tw > 0.0;
}

static void apl_mirror_component(AsmRun *run, int32_t comp)
{
    for (size_t i = 0; i < run->n_charts; i++) {
        AsmChart *c = &run->charts[i];
        if (!AsmChart_in_layout(c) || c->component != comp) continue;
        apl_compose(c, 0.0, 0.0, 0.0, 1);
    }
}

/* ---- turn vectors ---------------------------------------------------------------- */

static int apl_turn_pair_ok(const AsmRun *run, const uint8_t *comp_placed, const AsmLayerPair *lp)
{
    const AsmChart *A = &run->charts[(size_t)lp->a], *B = &run->charts[(size_t)lp->b];
    if (!AsmChart_in_layout(A) || !AsmChart_in_layout(B)) return 0;
    if (comp_placed[A->component] != 1 || comp_placed[B->component] != 1) return 0;
    if (A->component != B->component) return 0;
    return 1;
}

/* Turn vectors: for every placed chart and side, the median 2-D displacement
 * to its layer neighbours inside the same placed component.  A probe from Q
 * that hit P is evidence for P as well (P's displacement to Q is the
 * negative); which side of P that is needs both outward sides: Q found P on
 * Q's side s, so P is outward of Q iff s is Q's outward side, and then Q sits
 * on P's INWARD side. */
/* lidx (may be NULL): the per-chart layer coordinate; a placed pair whose coordinates differ by k
 * layers is a k-wrap adjacency (physical where the spacing is small and the mesh has holes: 53%
 * of the 10x10x10's placed pairs), and its displacement is k turns -- divided out here so a
 * measured turn is always ONE wrap (without this the 10x's t0 read 1,952 and its growth 154
 * per layer, both about double the physical values, and measured turns of neighbouring charts
 * disagreed by up to 92% of a turn). */
static void apl_measure_turns(const AsmRun *run, const uint8_t *comp_placed, const int8_t *out_side, const double *lidx, AplTurn *turn, int min_hits)
{
    Arena_Mark mark = Arena_save(run->arena);
    size_t nc = run->n_charts;
    size_t *off = ARENA_CALLOC(run->arena, nc * 2 + 1, sizeof(size_t));
    double *min_d = ARENA_ALLOC(run->arena, (nc * 2 + 1) * sizeof(double));
    for (size_t i = 0; i < nc * 2; i++) { turn[i].tx = 0.0; turn[i].ty = 0.0; turn[i].n = 0; min_d[i] = 1e300; }
    double *bx = NULL, *by = NULL;
    size_t *fill = NULL;
    /* pass 0: the nearest layer per chart AND side (the two sides of a chart
     * have different gaps; the per-chart distance would reject the far one);
     * pass 1: count hits of pairs within 35% of that; pass 2: fill */
    for (int pass = 0; pass < 3; pass++) {
        for (size_t p = 0; p < run->n_layers; p++) {
            const AsmLayerPair *lp = &run->layers[p];
            if (!apl_turn_pair_ok(run, comp_placed, lp)) continue;
            /* Excluded multi-layer pairs must not set the distance cutoff
             * for adjacent-layer measurements in pass zero. */
            if (APL_PAIR_K_FROM_PROBE && lp->k > 1) continue;
            int rev = 0, side_b = 0;
            if (out_side[lp->a] != 0 && out_side[lp->b] != 0) {
                side_b = (lp->side == out_side[lp->a]) ? -out_side[lp->b] : out_side[lp->b];
                rev = 1;
            }
            size_t sa = (size_t)lp->a * 2 + (lp->side > 0 ? 1 : 0);
            size_t sb = (size_t)lp->b * 2 + (side_b > 0 ? 1 : 0);
            if (pass == 0) {
                if (lp->d_median < min_d[sa]) min_d[sa] = lp->d_median;
                if (rev && lp->d_median < min_d[sb]) min_d[sb] = lp->d_median;
                continue;
            }
            int fwd = lp->d_median <= 1.35 * min_d[sa];
            if (rev) rev = lp->d_median <= 1.35 * min_d[sb];
            if (!fwd && !rev) continue;
            const AsmChart *A = &run->charts[(size_t)lp->a], *B = &run->charts[(size_t)lp->b];
            double kdiv = 1.0;
            if (APL_TURN_KHOP && lidx != NULL && lidx[lp->a] > APL_LAYER_UNKNOWN && lidx[lp->b] > APL_LAYER_UNKNOWN) {
                double dl = fabs(lidx[lp->b] - lidx[lp->a]);
                int kk = (int)floor(dl + 0.5);
                if (kk < 1) kk = 1;
                if (kk > 4) kk = 4;
                kdiv = (double)kk;
            }
            for (int k = 0; k < lp->hit_count; k++) {
                const AsmLayerHit *h = &run->layer_hits[(size_t)lp->hit_first + (size_t)k];
                if (fabs(h->dist - lp->d_median) > 0.25 * lp->d_median + 1.0) continue;
                if (pass == 1) { if (fwd) off[sa + 1]++; if (rev) off[sb + 1]++; continue; }
                double ox, oy, hx, hy;
                apl_vertex_point(A, h->va, &ox, &oy);
                apl_face_point(B, h->fb, h->l0, h->l1, &hx, &hy);
                if (fwd) { size_t q = off[sa] + fill[sa]++; bx[q] = (hx - ox) / kdiv; by[q] = (hy - oy) / kdiv; }
                if (rev) { size_t q = off[sb] + fill[sb]++; bx[q] = (ox - hx) / kdiv; by[q] = (oy - hy) / kdiv; }
            }
        }
        if (pass == 1) {
            for (size_t i = 0; i < nc * 2; i++) off[i + 1] += off[i];
            size_t total = off[nc * 2];
            bx = ARENA_ALLOC(run->arena, (total ? total : 1) * sizeof(double));
            by = ARENA_ALLOC(run->arena, (total ? total : 1) * sizeof(double));
            fill = ARENA_CALLOC(run->arena, nc * 2 + 1, sizeof(size_t));
        }
    }
    for (size_t sl = 0; sl < nc * 2; sl++) {
        size_t n = off[sl + 1] - off[sl];
        if ((int)n < min_hits) continue;
        qsort(bx + off[sl], n, sizeof(double), apl_cmp_double);
        qsort(by + off[sl], n, sizeof(double), apl_cmp_double);
        double tx = bx[off[sl] + n / 2], ty = by[off[sl] + n / 2], len = hypot(tx, ty);
        if (len < APL_TURN_MIN || len > APL_TURN_MAX) continue;   /* a fold or a broken layout, not a turn */
        turn[sl].tx = tx; turn[sl].ty = ty; turn[sl].n = (int)n;
    }
    Arena_restore(run->arena, mark);
}

/* Growth of the turn length per layer, measured on charts with both turns:
 * |T_out| - |T_in|.  The offset-curve rule says 2 pi d for a convex section
 * at uniform layer distance d; measuring it removes the local-d scatter
 * (the layer distance varies 9-18 vox around one turn, which is 25 vox of
 * turn length).  Returns the number of charts measured; *growth falls back
 * to 2 pi x the global layer distance when fewer than 3. */
/* BORROWED TURNS: a placed chart without a measured turn on a side takes the turn of the
 * nearest placed chart (in the layout) that has one on that side, when that chart is within
 * three quarters of its own turn length -- the same wrap neighbourhood, where the turn varies
 * by the growth per layer at most.  Measured motive 2026-09-08: on the familiar 532 of a cut
 * part's 540 anchors used ESTIMATED hop lengths although the primary holds 110 measured turns
 * (the touching charts just lacked their own), and the part landed a wrap off and was reverted.
 * Borrowed entries carry n = -1: hop lengths treat them as measured, the growth and t0
 * statistics (n > 0) do not. */
static int apl_cmp_double_fwd(const void *x, const void *y) { double a = *(const double *)x, b = *(const double *)y; return (a > b) - (a < b); }

/* v4: donors are chosen by the PLACED LAYOUT, not by the layer coordinate -- a placed chart within
 * APL_BORROW_REACH of a turn along u and 300 vox in v of the recipient is on the same wrap by
 * construction of the layout (neighbouring wraps sit one turn apart in u).  Up to 8 such donors
 * contribute; the borrowed turn is their MEDIAN length along their mean direction, so one wrong
 * donor cannot poison the recipient (v1-v3 took the single nearest and were measured dead). */
static size_t apl_borrow_turns(const AsmRun *run, const uint8_t *comp_placed, const int8_t *out_side, const double *lidx, AplTurn *turn)
{
    (void)lidx;
    size_t nc = run->n_charts, borrowed = 0;
    Arena_Mark mark = Arena_save(run->arena);
    int32_t *meas = ARENA_ALLOC(run->arena, (nc ? nc : 1) * sizeof(int32_t));
    double *mx = ARENA_ALLOC(run->arena, (nc ? nc : 1) * sizeof(double)), *my = ARENA_ALLOC(run->arena, (nc ? nc : 1) * sizeof(double));
    for (int want_out = 0; want_out < 2; want_out++) {
        size_t nm = 0;
        for (size_t i = 0; i < nc; i++) {
            if (out_side[i] == 0) continue;
            int sd = want_out ? (out_side[i] > 0 ? 1 : 0) : (out_side[i] > 0 ? 0 : 1);
            if (turn[i * 2 + (size_t)sd].n > 0) { apl_pose(&run->charts[i], 0.0, 0.0, &mx[nm], &my[nm]); meas[nm++] = (int32_t)i; }
        }
        if (nm == 0) continue;
        for (size_t i = 0; i < nc; i++) {
            const AsmChart *c = &run->charts[i];
            if (out_side[i] == 0) continue;
            int side = want_out ? (out_side[i] > 0 ? 1 : 0) : (out_side[i] > 0 ? 0 : 1);
            AplTurn *t = &turn[i * 2 + (size_t)side];
            if (t->n != 0 || !AsmChart_in_layout(c) || c->component < 0 || comp_placed[c->component] != 1) continue;
            double gx, gy;
            apl_pose(c, 0.0, 0.0, &gx, &gy);
            double lens[8], sx = 0.0, sy = 0.0, dist[8];
            size_t nd = 0;
            for (size_t k = 0; k < nm; k++) {
                int32_t d = meas[k];
                int dside = want_out ? (out_side[d] > 0 ? 1 : 0) : (out_side[d] > 0 ? 0 : 1);
                const AplTurn *src = &turn[(size_t)d * 2 + (size_t)dside];
                double len = hypot(src->tx, src->ty);
                double du = fabs(mx[k] - gx), dv = fabs(my[k] - gy);
                if (du > APL_BORROW_REACH * len || dv > 300.0) continue;
                double dd = hypot(du, dv);
                if (nd < 8) { lens[nd] = len; dist[nd] = dd; sx += src->tx / len; sy += src->ty / len; nd++; }
                else {
                    size_t worst = 0;
                    for (size_t q = 1; q < 8; q++) if (dist[q] > dist[worst]) worst = q;
                    if (dd < dist[worst]) { lens[worst] = len; dist[worst] = dd; }   /* direction sum keeps the first 8: fine, they agree */
                }
            }
            if (nd < 2) continue;
            qsort(lens, nd, sizeof(double), apl_cmp_double_fwd);
            double med = lens[nd / 2], dl = hypot(sx, sy);
            if (dl < 1e-9) continue;
            t->tx = med * sx / dl; t->ty = med * sy / dl; t->n = -(int)nd;
            borrowed++;
        }
    }
    Arena_restore(run->arena, mark);
    return borrowed;
}

/* the turn length at the primary's wrap: the median of its charts' measured outward turns, else the
 * turn of its mean curvature (2 pi / kappa), else 2 pi x 150 vox; bounded */
static double apl_turn_t0(const AsmRun *run, const AplTurn *turn, const int8_t *out_side, const int32_t *orig_comp, int32_t primary,
                          int *measured, size_t *ntl_out)
{
    Arena_Mark mt = Arena_save(run->arena);
    size_t nc = run->n_charts;
    double *tl = ARENA_ALLOC(run->arena, (nc ? nc : 1) * sizeof(double));
    size_t ntl = 0; double ksum = 0.0, kw = 0.0, t0;
    for (size_t i = 0; i < nc; i++) {
        const AsmChart *c = &run->charts[i];
        if (!AsmChart_in_layout(c) || orig_comp[i] != primary || out_side[i] == 0) continue;
        const AplTurn *t = &turn[i * 2 + (out_side[i] > 0 ? 1 : 0)];
        if (t->n > 0) tl[ntl++] = hypot(t->tx, t->ty);
        if (c->nv >= 12) {
            double tu[3], tv[3], huu[3];
            if (apl_chart_quadric(c, tu, tv, huu) == 0) {
                double hl = sqrt(huu[0]*huu[0] + huu[1]*huu[1] + huu[2]*huu[2]);
                if (hl > 1e-4 && hl < 0.2) { ksum += c->area3d * hl; kw += c->area3d; }
            }
        }
    }
    *measured = 0;
    if (ntl >= 3) { qsort(tl, ntl, sizeof(double), apl_cmp_double); t0 = tl[ntl / 2]; *measured = 1; }
    else if (kw > 0.0) t0 = 2.0 * APL_PI / (ksum / kw);
    else t0 = 2.0 * APL_PI * 150.0;
    if (t0 < APL_TURN_MIN) t0 = APL_TURN_MIN;
    if (t0 > APL_TURN_MAX) t0 = APL_TURN_MAX;
    *ntl_out = ntl;
    Arena_restore(run->arena, mt);
    return t0;
}

static size_t apl_measure_growth(const AsmRun *run, const AplTurn *turn, const int8_t *out_side, double *growth)
{
    Arena_Mark mark = Arena_save(run->arena);
    size_t nc = run->n_charts, n = 0;
    double *g = ARENA_ALLOC(run->arena, (nc ? nc : 1) * sizeof(double));
    for (size_t c = 0; c < nc; c++) {
        if (out_side[c] == 0) continue;
        const AplTurn *to = &turn[c * 2 + (out_side[c] > 0 ? 1 : 0)];
        const AplTurn *ti = &turn[c * 2 + (out_side[c] > 0 ? 0 : 1)];
        if (to->n <= 0 || ti->n <= 0) continue;
        g[n++] = hypot(to->tx, to->ty) - hypot(ti->tx, ti->ty);
    }
    if (n >= 3) {
        qsort(g, n, sizeof(double), apl_cmp_double);
        *growth = g[n/2];
    } else *growth = 2.0 * APL_PI * run->layer_d_global;
    Arena_restore(run->arena, mark);
    return n;
}

/* Turn vector of chart c on the given side.  Returns 1 measured, 2 estimated
 * from the other side by the growth rule, 0 unknown. */
static int apl_turn_for(const AplTurn *turn, int32_t c, int side, int8_t out_side, double growth, double *tx, double *ty)
{
    size_t s = (size_t)c * 2 + (side > 0 ? 1 : 0), o = (size_t)c * 2 + (side > 0 ? 0 : 1);
    if (turn[s].n != 0) { *tx = turn[s].tx; *ty = turn[s].ty; return 1; }   /* measured (n > 0) or borrowed from the nearest measured chart (n < 0) */
    if (turn[o].n != 0 && out_side != 0) {
        double len = hypot(turn[o].tx, turn[o].ty);
        if (len < 1e-6) return 0;
        double target = side == out_side ? len + growth : len - growth;   /* outward turns are longer */
        if (target < 0.3 * len || target < APL_TURN_MIN || target > APL_TURN_MAX) return 0;
        double f = target / len;
        *tx = -turn[o].tx * f; *ty = -turn[o].ty * f;
        return 2;
    }
    return 0;
}

/* k layers along the given side: k turns, each one growth longer (outward)
 * or shorter (inward) than the previous.  Returns 0 when the turn or the
 * side is unknown (k >= 2 needs the outward side). */
static int apl_turn_k(const AplTurn *turn, int32_t c, int side, int8_t out_side, double growth, int k, double *tx, double *ty)
{
    double t1x, t1y;
    int src = apl_turn_for(turn, c, side, out_side, growth, &t1x, &t1y);
    if (!src) return 0;
    if (k == 1) { *tx = t1x; *ty = t1y; return src; }
    if (src != 1) return 0;   /* an estimated turn carries one hop; two or three turns need a measured one */
    if (out_side == 0) return 0;
    double len = hypot(t1x, t1y);
    if (len < 1e-6) return 0;
    double extra = growth * (double)k * (double)(k - 1) * 0.5;
    double total = side == out_side ? (double)k * len + extra : (double)k * len - extra;
    if (total < 0.3 * len) return 0;
    *tx = t1x * total / len; *ty = t1y * total / len;
    return src;
}

/* ---- translation fit with a bounded rotation refinement ---------------------------- */

/* Where an anchor's turn vector came from.  Separates the families that CARRY
 * x->growth (OPPOSITE, MEASURED_K, LAYER) from those that do not, which is the
 * split the turn-model A/B is judged on.  Until 2026-09-19 this had to be
 * reverse-engineered from Ty == 0.0 plus magnitude matching against 2 pi k r. */
enum {
    APL_TSRC_SEAM = 0,      /* k = 0: an exact seam relation */
    APL_TSRC_MEASURED,      /* apl_turn_k -> 1 at k = 1: a measured turn on that side */
    APL_TSRC_MEASURED_K,    /* apl_turn_k -> 1 at k >= 2: measured + growth * k(k-1)/2 */
    APL_TSRC_OPPOSITE,      /* apl_turn_k -> 2: the other side's turn +/- growth */
    APL_TSRC_RADIUS,        /* apl_turn_radius: 2 pi k r from the axis table */
    APL_TSRC_LAYER,         /* t0 + lidx(P) * growth, summed per hop */
    APL_TSRC_STRETCH,       /* the pooled stretch estimate */
    APL_TSRC_T0,            /* a bare t0, measured or not */
    APL_TSRC_N
};

typedef struct AplAnchor {
    double sx, sy;        /* point of the component, its current frame */
    double tx, ty;        /* where the evidence puts it */
    double w;
    int    k, estimated, reverse;
    int    turn_src;      /* APL_TSRC_*: provenance of the turn vector (Tx, Ty) */
    int32_t chart_p; int8_t side; double Tx, Ty;
    int32_t chart_c;      /* the component's own chart touched by the hit */
    int    cluster;
    double seam_rms;      /* k = 0: the seam relation's rigid rms (its gate); 0 for hops */
    int32_t relation, source_run;
    int measured_order, measured_turn, source_flip;
    double measured_x, measured_y;
    double local_u, local_v;
    size_t observation;  /* same physical ray hit across alternate wrap hypotheses */
    const AsmLayerHit *hit; /* retained vertex/barycentric coordinates for UV trials */
} AplAnchor;

static double apl_fit_translation(const AplAnchor *a, const double *w, size_t n, double max_rot,
                                  double *theta, double *tx, double *ty)
{
    double sw = 0.0, cs0 = 0.0, cs1 = 0.0, cd0 = 0.0, cd1 = 0.0;
    for (size_t i = 0; i < n; i++) { sw += w[i]; cs0 += w[i]*a[i].sx; cs1 += w[i]*a[i].sy; cd0 += w[i]*a[i].tx; cd1 += w[i]*a[i].ty; }
    if (!(sw > 0.0)) return 1e300;
    cs0 /= sw; cs1 /= sw; cd0 /= sw; cd1 /= sw;
    double sxx = 0.0, sxy = 0.0;
    for (size_t i = 0; i < n; i++) {
        double sx = a[i].sx - cs0, sy = a[i].sy - cs1, dx = a[i].tx - cd0, dy = a[i].ty - cd1;
        sxx += w[i] * (sx*dx + sy*dy);
        sxy += w[i] * (sx*dy - sy*dx);
    }
    double th = atan2(sxy, sxx);
    if (th > max_rot) th = max_rot;
    if (th < -max_rot) th = -max_rot;
    double c = cos(th), s = sin(th);
    *theta = th;
    *tx = cd0 - (c*cs0 - s*cs1);
    *ty = cd1 - (s*cs0 + c*cs1);
    double e2 = 0.0;
    for (size_t i = 0; i < n; i++) {
        double px = c*a[i].sx - s*a[i].sy + *tx, py = s*a[i].sx + c*a[i].sy + *ty;
        e2 += w[i] * ((px-a[i].tx)*(px-a[i].tx) + (py-a[i].ty)*(py-a[i].ty));
    }
    return sqrt(e2 / sw);
}

typedef struct AplAxisFit {
    double mx, my, ax, ay, axis_v, slope, band, variance, dot, cross;
} AplAxisFit;

static double apl_axis_excess(double residual, double band)
{
    return copysign(fmax(0.0, fabs(residual)-band), residual);
}

static double apl_axis_fit_cost(const AplAxisFit *f, double theta)
{
    double c = cos(theta), s = sin(theta);
    double rx = c*f->mx-s*f->my, ry = s*f->mx+c*f->my;
    double perpendicular = apl_axis_excess(f->ay-f->axis_v-ry-f->slope*(f->ax-rx), f->band);
    return f->variance-2*(c*f->dot+s*f->cross)+perpendicular*perpendicular/(1+f->slope*f->slope);
}

/* Rigid least squares with the component centroid constrained to a band
 * around V = axis_v + slope*U. A measured height model has uncertainty;
 * fitting it as an exact equality can discard otherwise coherent anchors.
 * Translation in the band has a closed-form optimum for each angle.
 * Optimize the bounded angle jointly; band zero retains the exact plane. */
static double apl_fit_axis_translation(const AplAnchor *a, const double *w, size_t n, double max_rot,
                                       double cu, double cv, double axis_v, double slope, double band,
                                       double *theta, double *tx, double *ty)
{
    AplAxisFit f; memset(&f, 0, sizeof f); f.axis_v = axis_v; f.slope = slope; f.band = fmax(0.0,band);
    double sw = 0;
    for (size_t i = 0; i < n; i++) if (w[i] > 0) {
        sw += w[i]; f.mx += w[i]*(a[i].sx-cu); f.my += w[i]*(a[i].sy-cv);
        f.ax += w[i]*a[i].tx; f.ay += w[i]*a[i].ty;
    }
    if (!(sw > 0)) return 1e300;
    f.mx /= sw; f.my /= sw; f.ax /= sw; f.ay /= sw;
    for (size_t i = 0; i < n; i++) if (w[i] > 0) {
        double x = a[i].sx-cu-f.mx, y = a[i].sy-cv-f.my;
        double u = a[i].tx-f.ax, v = a[i].ty-f.ay;
        f.variance += w[i]*(x*x+y*y+u*u+v*v);
        f.dot += w[i]*(x*u+y*v); f.cross += w[i]*(x*v-y*u);
    }
    f.variance /= sw; f.dot /= sw; f.cross /= sw;
    double limit = fmin(fmax(max_rot, 0.0), APL_PI);
    double best = fmax(-limit, fmin(limit, atan2(f.cross, f.dot)));
    double best_cost = apl_axis_fit_cost(&f, best), step = 2*limit/128;
    for (int j = 0; j <= 128; j++) {
        double angle = -limit+j*step, cost = apl_axis_fit_cost(&f, angle);
        if (cost < best_cost) { best = angle; best_cost = cost; }
    }
    double lo = fmax(-limit, best-step), hi = fmin(limit, best+step);
    const double golden = 0.6180339887498948482;
    double l = hi-golden*(hi-lo), r = lo+golden*(hi-lo);
    double fl = apl_axis_fit_cost(&f, l), fr = apl_axis_fit_cost(&f, r);
    for (int j = 0; j < 40 && hi-lo > 1e-12; j++) {
        if (fl < fr) { hi = r; r = l; fr = fl; l = hi-golden*(hi-lo); fl = apl_axis_fit_cost(&f, l); }
        else { lo = l; l = r; fl = fr; r = lo+golden*(hi-lo); fr = apl_axis_fit_cost(&f, r); }
    }
    double angle = 0.5*(lo+hi);
    if (apl_axis_fit_cost(&f, angle) < best_cost) best = angle;
    double c = cos(best), s = sin(best), rx = c*f.mx-s*f.my, ry = s*f.mx+c*f.my;
    double excess = apl_axis_excess(f.ay-axis_v-ry-slope*(f.ax-rx), f.band)/(1+slope*slope);
    double centroid_u = f.ax-rx+slope*excess, centroid_v = f.ay-ry-excess;
    *theta = best; *tx = centroid_u-c*cu+s*cv;
    *ty = centroid_v-s*cu-c*cv;
    /* Recompute the returned residual directly, avoiding cancellation of
     * the centered quadratic when the final fit is nearly exact. */
    double e2 = 0;
    for (size_t i = 0; i < n; i++) if (w[i] > 0) {
        double x = c*a[i].sx-s*a[i].sy+*tx-a[i].tx;
        double y = s*a[i].sx+c*a[i].sy+*ty-a[i].ty;
        e2 += w[i]*(x*x+y*y);
    }
    return sqrt(e2/sw);
}

/* ---- driver ------------------------------------------------------------------------ */

typedef struct AplCompDiag { int status; int attempts; size_t last_anchors; double last_rms, last_gate, last_rot; int last_k1, last_k2, last_clusters;
                             int placed_by; size_t seam_anchors; unsigned rim_logged; } AplCompDiag;   /* placed_by: 0 none, 1 a seam cluster, 2 a hop cluster */
enum { APL_NEVER = 0, APL_PLACED = 1, APL_REFUSED = 2, APL_REVERTED = 3, APL_PRIMARY = 4, APL_FEW_ANCHORS = 5, APL_DEFERRED = 6 };
static const char *apl_status_names[] = { "never", "placed", "refused", "reverted", "primary", "few_anchors", "deferred" };

typedef struct AplRimGrid AplRimGrid;
typedef struct AplStrict AplStrict;
typedef struct AplConsensus { size_t observations, inliers, model_inliers; double clipped_loss, model_loss; } AplConsensus;
typedef struct AplCtx {
    AsmRun *run; const AsmPlaceOpts *o; const AsmConflictOpts *copts;
    uint8_t *comp_placed; AplTurn *turn; int8_t *out_side; int32_t *orig_comp; double growth;
    AplAnchor *anchors; double *wts; size_t hcap;
    int32_t primary; double *comp_area; AplCompDiag *diag; FILE *fanch; int attempt;
    int stretch_ok; double stretch_T, stretch_s, stretch_d, stretch_rel; size_t stretch_n, stretch_pairs;   /* per-attempt stretch estimate and the relative error of (s-1) */
    double turn_rel;          /* the worst relative error of any turn this candidate's anchors rest on (0 = all measured) */
    double cand_grad;   /* the candidate's own z-gradient magnitude */
    AsmConflictGrid *cgrid; AplRimGrid *rgrid;   /* persistent audits of the placed material */
    const uint8_t *excluded_samples;          /* old samples of a staged relocation only */
    double strip_umin, strip_umax, span_cap;      /* the placed strip's u extent and its physical cap */
    const double *lidx;    /* per chart: layer coordinate relative to the primary's anchor chart (APL_LAYER_UNKNOWN = none) */
    double t0;             /* turn length at the primary's wrap; T(w) = t0 + w x growth */
    double *tfield;        /* per chart: the turn field's hop length at the chart (0 = none) */
    int t0_measured;       /* t0 is the median of measured turns (else a curvature estimate: the layer model is not armed) */
    /* the axis: per chart axial coordinate s and radius r (r only with an axis point) */
    double *s_axis, *r_axis; int have_radius;
    double *cu, *cv;                         /* per chart: the uv centroid in the chart frame (unmirrored) */
    double *tr_top, *tr_bot;                 /* OUTLINE scratch: crop-plane vertex v per bin (CSR by off_*) */
    size_t *off_t, *off_b, *fill_t, *fill_b;
    double vm_a, vm_b, vm_c; int vm_ok;      /* the primary's v = a + b s + c u */
    double vm_corr;                          /* the s-u correlation over the primary (the u term is dropped when collinear) */
    double vm_p90;                           /* the primary's own |v - v(s,u)| p90 about the model, vox */
    int vm_pin;                              /* the pinning is armed (b within the physical band) */
    int diagnostic_axis_policy;              /* case replay only: 1 SOURCE prior, 2 all height priors */
    int diagnostic_first_step;               /* stop local UV search after a certified improvement */
    AplStrict *strict;                       /* opt-in original-face certificate before committing a proposal */
    double uout_x, uout_y; int uout_ok;      /* the layout's outward turn direction (unit) */
    double turn_c; size_t turn_c_n;          /* T(r) = 2 pi r + turn_c: the section's perimeter excess, fitted over turn_c_n charts (0 = none adopted) */
    double turn_slope;                       /* the free slope of |T| on r (2 pi on a circular section); 0 = unfitted */
    size_t anchors_radius;
    AsmPlaceStats *st;                       /* the run's stats: the profile counters */
    double *curv_hl; uint8_t *curv_ok;       /* APL_CURV_CACHE: per placed chart, |h_uu| of its quadric (-1 = the fit failed) */
    size_t na_seam;                          /* seam anchors (k = 0) in the current candidate's evidence */
    int last_silent;                         /* the last audit found insufficient or ambiguous evidence */
    int last_shape_refused;
    const char *last_reason;                 /* which gate refused the last candidate, for stage5_candidates.csv */
    int source_geometry_conflict;            /* this attempt has a supported boundary fit that failed geometry */
    size_t na; int cluster, cluster_kind;
    double evidence_gate;
    uint8_t *seam_selected;
    int trial, source_flip; AsmChart *proposal; size_t nproposal;
    Arena_T proposal_uv;
    int uv_trial;
    size_t n_observations;
    AplConsensus consensus;
    FILE *fevidence, *frim, *fcut;
    double fit_cu, fit_cv, fit_s, fit_area;
    size_t cuts_preserved_source;
} AplCtx;

/* Use the primary's observed height uncertainty, with the same scatter
 * multiplier and shift cap as the final shape gate. A zero-scatter model
 * still specifies an exact plane. */
static double apl_axis_band(const AplCtx *x)
{
    return fmax(0.0,1.5*x->vm_p90);
}

static void apl_uv_centroid(const AsmChart *c, double *u, double *v)
{
    *u = *v = 0;
    for (size_t i = 0; i < c->nv; i++) { *u += c->uv[2*i]; *v += c->uv[2*i+1]; }
    if (c->nv) { *u /= c->nv; *v /= c->nv; }
}

/* Candidate and winning UV buffers have separate lifetimes. Never retain
 * a pointer into a transaction's scratch arena or mutate the original UV
 * while merely comparing two placements. */
static void apl_copy_candidate(AsmChart *out, const AsmChart *in, size_t n,
                                const AsmRun *run, Arena_T *uv_arena)
{
    if (*uv_arena) Arena_free(*uv_arena);
    for (size_t i = 0; i < n; i++) {
        out[i] = in[i];
        if (in[i].uv != run->charts[in[i].id].uv) {
            if (!*uv_arena) *uv_arena = Arena_new();
            out[i].uv = ARENA_ALLOC(*uv_arena,in[i].nv*2*sizeof(float));
            memcpy(out[i].uv,in[i].uv,in[i].nv*2*sizeof(float));
        }
    }
}

static void apl_commit_chart(AplCtx *x, const AsmChart *candidate, float *owned_uv)
{
    AsmChart c = *candidate;
    if (c.uv != owned_uv) {
        memcpy(owned_uv,c.uv,c.nv*2*sizeof(float)); c.uv = owned_uv;
        if (x->cu && x->cv) apl_uv_centroid(&c,&x->cu[c.id],&x->cv[c.id]);
        if (x->curv_ok) x->curv_ok[c.id] = 0;
    }
    x->run->charts[c.id] = c;
}

/* The chart's centroid in the layout: its pose composed with its uv centroid (mirrored as the
 * chart is); with (th, tx, ty) the tentative component transform on top. */
static void apl_layout_centroid(const AplCtx *x, const AsmChart *c, size_t i, double th, double tx, double ty, double *gx, double *gy)
{
    if (c->placed_uv) {
        double u=0,v=0;
        for(size_t k=0;k<c->nv;k++){u+=c->placed_uv[2*k];v+=c->placed_uv[2*k+1];}
        if(c->nv){u/=c->nv;v/=c->nv;}
        *gx=cos(th)*u-sin(th)*v+tx; *gy=sin(th)*u+cos(th)*v+ty;
        return;
    }
    double u = (c->flags & ASM_CHART_MIRROR) ? -x->cu[i] : x->cu[i], v = x->cv[i];
    double ct = cos(c->pose_theta), st = sin(c->pose_theta);
    double px = ct * u - st * v + c->pose_x, py = st * u + ct * v + c->pose_y;
    double c2 = cos(th), s2 = sin(th);
    *gx = c2 * px - s2 * py + tx; *gy = s2 * px + c2 * py + ty;
}

/* Area-weighted v(s,u) from the initial primary, or from committed placed
 * material when recalibrating. A pending candidate never supplies model
 * evidence. The U term needs independent axial and winding extent: a
 * diagonal band alone cannot determine both coefficients. */
static int apl_v_model_fit(AplCtx *x, int placed)
{
    AsmRun *run = x->run;
    double A[9] = { 0 }, B[3] = { 0 };
    double n = 0.0, ws = 0.0, ss = 0.0, su = 0.0, sss = 0.0, suu = 0.0, ssu = 0.0;
    for (size_t i = 0; i < run->n_charts; i++) {
        const AsmChart *c = &run->charts[i];
        if (!AsmChart_in_layout(c) || c->area3d <= 0.0) continue;
        if (placed ? c->component != x->primary || c->placement_state == ASM_PLACE_NONE : x->orig_comp[i] != x->primary) continue;
        double u, v;
        apl_layout_centroid(x, c, i, 0.0, 0.0, 0.0, &u, &v);
        double w = c->area3d, s = x->s_axis[i];
        double f[3] = { 1.0, s, u };
        for (int r = 0; r < 3; r++) { for (int k = 0; k < 3; k++) A[r*3+k] += w * f[r] * f[k]; B[r] += w * f[r] * v; }
        n += 1.0; ws += w; ss += w * s; su += w * u; sss += w * s * s; suu += w * u * u; ssu += w * s * u;
    }
    x->vm_ok = 0; x->vm_corr = 0.0;
    if (n < 5.0) return -1;
    {
        double vs = sss / ws - (ss / ws) * (ss / ws), vu = suu / ws - (su / ws) * (su / ws), cv = ssu / ws - (ss / ws) * (su / ws);
        x->vm_corr = (vs > 1e-9 && vu > 1e-9) ? cv / sqrt(vs * vu) : 1.0;
        int drop_u = fabs(x->vm_corr) > APL_V_MODEL_MAX_CORR || vu < 250.0 * 250.0;
        if (!drop_u) {
            /* the full fit first: a u coefficient beyond a plausible frame rotation drops the term */
            double M[9], R[3];
            memcpy(M, A, sizeof M); memcpy(R, B, sizeof R);
            int ok = 1;
            for (int col = 0; col < 3 && ok; col++) {
                int piv = col;
                for (int r = col + 1; r < 3; r++) if (fabs(M[r*3+col]) > fabs(M[piv*3+col])) piv = r;
                for (int k = 0; k < 3; k++) { double tt = M[col*3+k]; M[col*3+k] = M[piv*3+k]; M[piv*3+k] = tt; }
                { double tt = R[col]; R[col] = R[piv]; R[piv] = tt; }
                if (fabs(M[col*3+col]) < 1e-9) { ok = 0; break; }
                for (int r = 0; r < 3; r++) {
                    if (r == col) continue;
                    double f = M[r*3+col] / M[col*3+col];
                    for (int k = 0; k < 3; k++) M[r*3+k] -= f * M[col*3+k];
                    R[r] -= f * R[col];
                }
            }
            if (!ok || fabs(R[2] / M[8]) > APL_V_MODEL_MAX_C) drop_u = 1;
        }
        if (drop_u) {
            /* no u term: solve the 2x2 for (a, b) */
            double det = A[0] * A[4] - A[1] * A[3];
            if (fabs(det) < 1e-9) return -1;
            x->vm_a = (B[0] * A[4] - A[1] * B[1]) / det; x->vm_b = (A[0] * B[1] - A[3] * B[0]) / det; x->vm_c = 0.0;
            x->vm_ok = 1;
            return 0;
        }
    }
    for (int col = 0; col < 3; col++) {
        int piv = col;
        for (int r = col + 1; r < 3; r++) if (fabs(A[r*3+col]) > fabs(A[piv*3+col])) piv = r;
        for (int k = 0; k < 3; k++) { double t = A[col*3+k]; A[col*3+k] = A[piv*3+k]; A[piv*3+k] = t; }
        { double t = B[col]; B[col] = B[piv]; B[piv] = t; }
        if (fabs(A[col*3+col]) < 1e-9) return -1;
        for (int r = 0; r < 3; r++) {
            if (r == col) continue;
            double f = A[r*3+col] / A[col*3+col];
            for (int k = 0; k < 3; k++) A[r*3+k] -= f * A[col*3+k];
            B[r] -= f * B[col];
        }
    }
    x->vm_a = B[0] / A[0]; x->vm_b = B[1] / A[4]; x->vm_c = B[2] / A[8];
    x->vm_ok = 1;
    return 0;
}

static int apl_v_model(AplCtx *x)
{
    return apl_v_model_fit(x,0);
}

typedef struct AplHeightGroup {
    int32_t component;
    double area, residual, squared;
} AplHeightGroup;

static int apl_height_group_cmp(const void *a, const void *b)
{
    int32_t x = ((const AplHeightGroup *)a)->component, y = ((const AplHeightGroup *)b)->component;
    return (x > y)-(x < y);
}

/* The first component is a seed, not a permanent complete-sheet height
 * calibration. New committed material can improve the extrapolation.
 * Accept a refit only if its scale is physical, its whole placed-area error
 * decreases, and every incumbent component retains the SAME mean-height
 * and scatter limits. The seed uncertainty and all geometry remain fixed. */
static int apl_v_model_update(AplCtx *x)
{
    if (!APL_V_FROM_AXIS || !x->vm_ok || !x->vm_pin) return 0;
    AplCtx fitted = *x;
    if (apl_v_model_fit(&fitted,1) || !fitted.vm_ok || !isfinite(fitted.vm_a) ||
        !isfinite(fitted.vm_b) || !isfinite(fitted.vm_c) ||
        fitted.vm_b < APL_V_MODEL_B_MIN || fitted.vm_b > APL_V_MODEL_B_MAX ||
        fabs(fitted.vm_c) > APL_V_MODEL_MAX_C) return 0;
    Arena_Mark mark = Arena_save(x->run->arena);
    AplHeightGroup *groups = ARENA_ALLOC(x->run->arena,(x->run->n_charts ? x->run->n_charts : 1)*sizeof *groups);
    size_t count = 0;
    double before = 0, after = 0;
    int ok = 1;
    for (size_t i = 0; i < x->run->n_charts; i++) {
        const AsmChart *c = &x->run->charts[i];
        if (!AsmChart_in_layout(c) || c->area3d <= 0 || c->component != x->primary || c->placement_state == ASM_PLACE_NONE) continue;
        double u,v; apl_layout_centroid(x,c,i,0,0,0,&u,&v);
        double was = x->vm_a+x->vm_b*x->s_axis[i]+x->vm_c*u-v;
        double now = fitted.vm_a+fitted.vm_b*x->s_axis[i]+fitted.vm_c*u-v;
        if (!isfinite(was) || !isfinite(now)) { ok = 0; break; }
        groups[count++] = (AplHeightGroup){x->orig_comp[i],c->area3d,c->area3d*now,c->area3d*now*now};
        before += c->area3d*was*was; after += c->area3d*now*now;
    }
    if (!count || !isfinite(before) || !isfinite(after) || !(after < before-fmax(1e-8,1e-8*before))) ok = 0;
    if (ok) qsort(groups,count,sizeof *groups,apl_height_group_cmp);
    for (size_t i = 0; i < count && ok;) {
        size_t j = i+1; double area = groups[i].area, residual = groups[i].residual, squared = groups[i].squared;
        while (j < count && groups[j].component == groups[i].component) {
            area += groups[j].area; residual += groups[j].residual; squared += groups[j].squared; j++;
        }
        double mean = residual/area, scatter = sqrt(fmax(0.0,squared/area-mean*mean));
        if (fabs(mean) > fmax(APL_V_AXIS_MAX_SHIFT,2*x->vm_p90) ||
            scatter > fmax(APL_V_AXIS_TOL,1.5*x->vm_p90)) ok = 0;
        i = j;
    }
    Arena_restore(x->run->arena,mark);
    if (!ok) return 0;
    x->vm_a = fitted.vm_a; x->vm_b = fitted.vm_b; x->vm_c = fitted.vm_c; x->vm_corr = fitted.vm_corr;
    return 1;
}

/* the primary's own scatter about the model (p90), the yardstick of the v-axis gate */
static double apl_v_model_p90(AplCtx *x)
{
    AsmRun *run = x->run;
    Arena_Mark mark = Arena_save(run->arena);
    double *r = ARENA_ALLOC(run->arena, (run->n_charts ? run->n_charts : 1) * sizeof(double));
    size_t n = 0;
    for (size_t i = 0; i < run->n_charts; i++) {
        const AsmChart *c = &run->charts[i];
        if (!AsmChart_in_layout(c) || x->orig_comp[i] != x->primary) continue;
        double gx, gy;
        apl_layout_centroid(x, c, i, 0.0, 0.0, 0.0, &gx, &gy);
        r[n++] = fabs(x->vm_a + x->vm_b * x->s_axis[i] + x->vm_c * gx - gy);
    }
    double p90 = 0.0;
    if (n > 0) { qsort(r, n, sizeof(double), apl_cmp_double); p90 = r[(size_t)((double)(n - 1) * 0.9)]; }
    Arena_restore(run->arena, mark);
    return p90;
}

/* The layout's outward turn direction: the mean unit vector of the primary's measured outward
 * turns; without any, the sign of dr/du over the primary along +u. */
static void apl_uout(AplCtx *x)
{
    AsmRun *run = x->run;
    double sx = 0.0, sy = 0.0;
    size_t n = 0;
    x->uout_ok = 0;
    for (size_t i = 0; i < run->n_charts; i++) {
        const AsmChart *c = &run->charts[i];
        if (!AsmChart_in_layout(c) || x->orig_comp[i] != x->primary || x->out_side[i] == 0) continue;
        const AplTurn *t = &x->turn[i * 2 + (x->out_side[i] > 0 ? 1 : 0)];
        if (t->n <= 0) continue;
        double l = hypot(t->tx, t->ty);
        if (l < 1e-6) continue;
        sx += t->tx / l; sy += t->ty / l; n++;
    }
    if (n >= 3 && hypot(sx, sy) > 0.5 * (double)n) { double l = hypot(sx, sy); x->uout_x = sx / l; x->uout_y = sy / l; x->uout_ok = 1; return; }
    if (!x->have_radius) return;
    double su = 0.0, sr = 0.0, suu = 0.0, sur = 0.0, m = 0.0;
    for (size_t i = 0; i < run->n_charts; i++) {
        const AsmChart *c = &run->charts[i];
        if (!AsmChart_in_layout(c) || x->orig_comp[i] != x->primary) continue;
        double gx, gy;
        apl_layout_centroid(x, c, i, 0.0, 0.0, 0.0, &gx, &gy);
        su += gx; sr += x->r_axis[i]; suu += gx * gx; sur += gx * x->r_axis[i]; m += 1.0;
    }
    if (m < 5.0) return;
    double det = m * suu - su * su;
    if (fabs(det) < 1e-9) return;
    double slope = (m * sur - su * sr) / det;
    x->uout_x = slope >= 0.0 ? 1.0 : -1.0; x->uout_y = 0.0; x->uout_ok = 1;
}

/* The section's perimeter excess.  A wrap is not a circle: its cross-section is a spiral arm
 * with ends, so its length is 2 pi r + c with c > 0, and apl_turn_radius drops the c.  Fit it
 * from the charts that HAVE a measured outward turn: e_i = |T_i| - 2 pi r_i, four rounds of a
 * median with a multiplicative reject, so one mis-gauged wrap (|T| near 2x the truth) cannot
 * drag it.  Leaves turn_c = 0 -- today's behaviour -- unless the fit is well supported.
 * Also stores the free slope of |T| on r (2 pi on a true circular section), which is what
 * ASM_PLACE_GROWTH_PER_WRAP mode 2 carries to the growth. */
static void apl_fit_turn_radius(AplCtx *x)
{
    const AsmRun *run = x->run;
    x->turn_c = 0.0; x->turn_c_n = 0; x->turn_slope = 0.0;
    if (!x->have_radius || x->r_axis == NULL) return;
    Arena_Mark mark = Arena_save(run->arena);
    size_t nc = run->n_charts, n = 0;
    double *e = ARENA_ALLOC(run->arena, (nc ? nc : 1) * sizeof(double));
    double *tl = ARENA_ALLOC(run->arena, (nc ? nc : 1) * sizeof(double));
    double *rr = ARENA_ALLOC(run->arena, (nc ? nc : 1) * sizeof(double));
    for (size_t i = 0; i < nc; i++) {
        const AsmChart *c = &run->charts[i];
        if (!AsmChart_in_layout(c) || x->out_side[i] == 0) continue;
        const AplTurn *t = &x->turn[i * 2 + (x->out_side[i] > 0 ? 1 : 0)];
        if (t->n <= 0) continue;                       /* measured only; a borrowed turn (n < 0) is not evidence */
        double L = hypot(t->tx, t->ty), r = x->r_axis[i];
        if (!(L > APL_TURN_MIN) || !(r > 10.0)) continue;
        tl[n] = L; rr[n] = r; e[n] = L - 2.0 * APL_PI * r; n++;
    }
    if (n < (size_t)ASM_PLACE_TURN_C_MIN_N) { Arena_restore(run->arena, mark); return; }
    double c = 0.0;
    for (int round = 0; round < 4; round++) {
        size_t m = 0;
        for (size_t i = 0; i < n; i++) {
            if (round > 0) {
                double pred = 2.0 * APL_PI * rr[i] + c;
                if (!(pred > 1e-6)) continue;
                double ratio = tl[i] / pred;
                if (ratio < 0.6 || ratio > 1.45) continue;   /* a wrap gauged one turn out, not a noisy one */
            }
            e[m++] = tl[i] - 2.0 * APL_PI * rr[i];
        }
        if (m < (size_t)ASM_PLACE_TURN_C_MIN_N) break;
        qsort(e, m, sizeof(double), apl_cmp_double);
        c = e[m / 2];
        x->turn_c_n = m;
    }
    /* guard: an intercept may not swallow the smallest radius, nor rival a whole turn */
    double rmin = rr[0];
    for (size_t i = 1; i < n; i++) if (rr[i] < rmin) rmin = rr[i];
    if (x->turn_c_n < (size_t)ASM_PLACE_TURN_C_MIN_N ||
        c <= -0.25 * 2.0 * APL_PI * rmin || fabs(c) >= 0.5 * APL_TURN_MAX) { x->turn_c = 0.0; x->turn_c_n = 0; }
    else x->turn_c = c;
    if (n >= 2) {   /* free slope of |T| on r, for the growth's mode 2 */
        double sr = 0.0, st = 0.0, srr = 0.0, srt = 0.0, w = (double)n;
        for (size_t i = 0; i < n; i++) { sr += rr[i]; st += tl[i]; srr += rr[i]*rr[i]; srt += rr[i]*tl[i]; }
        double det = w * srr - sr * sr;
        if (fabs(det) > 1e-9) x->turn_slope = (w * srt - sr * st) / det;
    }
    Arena_restore(run->arena, mark);
}

/* The growth per PHYSICAL WRAP that apl_turn_for and apl_turn_k consume.  Mode 0 is
 * apl_measure_growth untouched (shipped, and provably today); 1 is the offset-curve rule
 * 2 pi d; 2 carries the fitted slope of |T| on r to d, clamped.  See
 * ASM_PLACE_GROWTH_PER_WRAP.  Returns apl_measure_growth's support count either way, so the
 * log line keeps reporting how many charts the shipped estimator rested on. */
static size_t apl_growth_for_hops(const AplCtx *x, const AsmRun *run, const AplTurn *turn,
                                  const int8_t *out_side, double *growth)
{
    size_t n = apl_measure_growth(run, turn, out_side, growth);
    if (!ASM_PLACE_GROWTH_PER_WRAP) return n;
    double d = run->layer_d_global > 0.0 ? run->layer_d_global : 9.5;
    double rule = 2.0 * APL_PI * d;
    if (ASM_PLACE_GROWTH_PER_WRAP == 2 && x->turn_slope > 0.0) {
        double g = x->turn_slope * d, lo = 0.8 * rule, hi = 1.25 * rule;
        *growth = g < lo ? lo : (g > hi ? hi : g);
    } else *growth = rule;
    return n;
}

/* Hop of k wraps from chart P on `outward` (1) or inward (0): 2 pi k (r +- k d / 2) along the
 * layout's outward direction.  Returns 1 with the vector, 0 when the radius model is not armed. */
static int apl_turn_radius(const AplCtx *x, int32_t P, int outward, int k, double *tx, double *ty)
{
    if (!APL_TURN_FROM_RADIUS || !x->have_radius || !x->uout_ok) return 0;
    const AsmRun *run = x->run;
    double d = run->chart_layer_d ? run->chart_layer_d[P] : 0.0;
    if (d <= 0.0) d = run->layer_d_global > 0.0 ? run->layer_d_global : 9.5;
    double r = x->r_axis[P] + (outward ? 1.0 : -1.0) * (double)k * d * 0.5;
    if (r < 10.0) r = 10.0;
    /* each of the k wraps carries one perimeter excess, so the intercept enters k times */
    double L = 2.0 * APL_PI * (double)k * r + (ASM_PLACE_TURN_INTERCEPT ? (double)k * x->turn_c : 0.0);
    if (L < APL_TURN_MIN || L > APL_TURN_MAX * (double)k) return 0;
    double sgn = outward ? 1.0 : -1.0;
    *tx = sgn * x->uout_x * L; *ty = sgn * x->uout_y * L;
    return 1;
}

/* TURN LENGTH FROM THE LAYER STRETCH.  When no placed chart carries a
 * measured turn vector (no component spans two layers: the arbitrary-start /
 * crooked case), the correspondences between a placed arc and the candidate
 * still say how long one turn is: the outer of two concentric arcs is longer
 * by (R + d) / R, so the stretch s of the candidate's u against the placed
 * u over their hits gives R = d / (s - 1) and one turn = 2 pi R (the offset-
 * curve rule, applied at the arc instead of around the turn).  Pooled over
 * every hit between the placed material and the candidate so the fit spans
 * the whole arc.  Returns 1 with *T (signed along +u: positive = the
 * candidate is outward) when s is measurable and the turn is plausible. */
static int apl_stretch_turn(AplCtx *x, int32_t best, double *T, double *s_out, double *dbar_out, size_t *nhits,
                            double *rel_out, size_t *npairs_out)
{
    AsmRun *run = x->run;
    double n = 0.0, su = 0.0, sb = 0.0, suu = 0.0, sub = 0.0, sbb = 0.0, dsum = 0.0, umin = 1e300, umax = -1e300;
    size_t npairs = 0;
    double dirvote = 0.0;   /* > 0: the candidate is OUTWARD of the placed material (bending-side vote) */
    for (size_t p = 0; p < run->n_layers; p++) {
        const AsmLayerPair *lp = &run->layers[p];
        const AsmChart *A = &run->charts[(size_t)lp->a], *B = &run->charts[(size_t)lp->b];
        if (!AsmChart_in_layout(A) || !AsmChart_in_layout(B)) continue;
        int forward = x->comp_placed[A->component] == 1 && B->component == best;
        int reverse = A->component == best && x->comp_placed[B->component] == 1;
        if (!forward && !reverse) continue;
        if (x->out_side[lp->a] != 0) {
            /* forward: B (candidate) sits on A's side lp->side; reverse: B (placed) sits on the
             * candidate A's side lp->side, so the candidate is outward iff that is A's inward side */
            int cand_outward = forward ? (lp->side == x->out_side[lp->a]) : (lp->side != x->out_side[lp->a]);
            dirvote += cand_outward ? (double)lp->hit_count : -(double)lp->hit_count;
        }
        int pair_used = 0;
        for (int h = 0; h < lp->hit_count; h++) {
            const AsmLayerHit *hit = &run->layer_hits[(size_t)lp->hit_first + (size_t)h];
            if (fabs(hit->dist - lp->d_median) > 0.25 * lp->d_median + 1.0) continue;
            double ax, ay, bx, by, up, ub;
            apl_vertex_point(A, hit->va, &ax, &ay);
            apl_face_point(B, hit->fb, hit->l0, hit->l1, &bx, &by);
            if (forward) { up = ax; ub = bx; } else { up = bx; ub = ax; }
            n += 1.0; su += up; sb += ub; suu += up * up; sub += up * ub; sbb += ub * ub; dsum += lp->d_median;
            pair_used = 1;
            if (up < umin) umin = up;
            if (up > umax) umax = up;
        }
        npairs += (size_t)pair_used;   /* the hits of one layer pair are one observation of the stretch, not hit_count of them */
    }
    *nhits = (size_t)n;
    if (npairs_out) *npairs_out = npairs;
    if (rel_out) *rel_out = -1.0;
    if (n < 12.0 || umax - umin < 100.0) return 0;
    double var = suu - su * su / n, cov = sub - su * sb / n;
    if (var <= 1e-9) return 0;
    double sl = cov / var, dbar = dsum / n;
    *s_out = sl; *dbar_out = dbar;
    /* HOW WELL IS THE TURN KNOWN?  T is 2*pi*d/(s-1), a ratio of small differences: on the
     * 21x5x5 s sits within 4% of 1, so a 1% error in s is a 25% error in T, and the
     * per-component turns that come out span 645 to 9,105 vox.  The slope's standard error is
     * free from the same accumulators -- se(s)^2 = sigma^2 / Sxx -- and se(s-1) = se(s),
     * because subtracting u_placed from the response leaves every residual untouched.
     *
     * The degrees of freedom are LAYER PAIRS, not hits.  The hits of one pair share that
     * pair's geometry, its d_median and its two charts' poses, so fifty of them are one look
     * at the stretch, not fifty; counting hits would let a densely probed pair claim a
     * precision it has not got.  se = sqrt(SSE / ((pairs - 2) * Sxx)) is exactly invariant
     * when the same geometry is recorded twice, which is the property the control checks. */
    {
        double syy = sbb - sb * sb / n;
        double sse = syy - cov * cov / var;
        if (sse < 0.0) sse = 0.0;
        double dof = (double)npairs - 2.0;
        double se = dof >= 1.0 ? sqrt(sse / (dof * var)) : -1.0;
        double delta = fabs(sl - 1.0);
        if (rel_out) *rel_out = se < 0.0 ? 1e300 : (delta > 1e-12 ? se / delta : 1e300);
    }
    if (sl < 0.5 || sl > 2.0) return 0;              /* wrong mirror or garbage: no estimate */
    /* the DIRECTION comes from the bending-side vote of the hits (which side of the
     * placed charts the candidate sits on), not from the stretch's sign: s within a
     * few percent of 1 has no reliable sign */
    int outward = dirvote != 0.0 ? dirvote > 0.0 : sl > 1.0;
    /* the curvature of the arcs bounds the turn: the turn through a chart is
     * 2 pi / kappa for a round section and within a factor ~2 of it for an
     * oval; a stretch estimate outside [0.4, 2.5] x that (s within 2% of 1
     * amplifies the hit noise into any R) is replaced by the curvature turn */
    double kw = 0.0, ksum = 0.0;
    for (size_t i = 0; i < run->n_charts; i++) {
        const AsmChart *c = &run->charts[i];
        if (!AsmChart_in_layout(c) || c->nv < 12) continue;
        int on_cand = c->component == best, on_placed = x->comp_placed[c->component] == 1;
        if (outward ? !on_placed : !on_cand) continue;   /* outward: the turn through the PLACED arc */
        double hl;
        if (APL_CURV_CACHE && on_placed && x->curv_ok != NULL && x->curv_ok[i]) { hl = x->curv_hl[i]; if (x->st) x->st->quadric_cached++; }
        else {
            double tu[3], tv[3], huu[3];
            if (x->st) x->st->quadric_calls++;
            hl = apl_chart_quadric(c, tu, tv, huu) != 0 ? -1.0 : sqrt(huu[0]*huu[0] + huu[1]*huu[1] + huu[2]*huu[2]);
            if (APL_CURV_CACHE && on_placed && x->curv_ok != NULL) { x->curv_hl[i] = hl; x->curv_ok[i] = 1; }   /* a placed pose is frozen inside the loop */
        }
        if (hl < 0.0) continue;                          /* the quadric fit failed */
        if (hl < 1e-4 || hl > 0.2) continue;             /* R between 5 and 10,000 vox */
        ksum += c->area3d * hl; kw += c->area3d;
    }
    double t_curv = kw > 0.0 ? 2.0 * APL_PI / (ksum / kw) : 0.0;
    double turn = 0.0;
    if (outward && sl > 1.002) turn = 2.0 * APL_PI * dbar / (sl - 1.0);
    else if (!outward && sl < 0.998) turn = 2.0 * APL_PI * dbar * sl / (1.0 - sl);
    /* An estimate the data does not support falls back to the curvature turn.  At
     * ASM_STRETCH_MAX_REL_ERR 0 that judgement is the shipped fixed bar on |s-1|; positive,
     * it is the measured relative error of (s-1), which is the same statement about T. */
    double rel = rel_out ? *rel_out : -1.0;
    (void)rel;   /* the shipped fixed bar does not read it; the armed branch does */
    int unsupported = ASM_STRETCH_MAX_REL_ERR > 0.0
                    ? (rel < 0.0 || rel > ASM_STRETCH_MAX_REL_ERR)
                    : (fabs(sl - 1.0) < 0.02);
    if (t_curv > 0.0 && (turn <= 0.0 || unsupported || turn < 0.4 * t_curv || turn > 2.5 * t_curv)) turn = t_curv;
    if (turn < APL_TURN_MIN || turn > APL_TURN_MAX) return 0;
    *T = outward ? turn : -turn;
    return 1;
}

/* LAYER COORDINATE: one number per chart, decided GLOBALLY before any placement.
 * Charts joined by an accepted seam sit on the same layer (difference 0); a
 * layer probe from chart a that meets chart b at d vox says b sits
 * round(d / local layer distance) layers to a's outward side (+) or inward
 * side (-).  The weighted least-squares solution over every such edge, with
 * the primary's biggest chart at 0 and three Cauchy reweights to silence the
 * contradicting edges, is a radial layer number: the primary's own turns
 * come out 0, 1, 2, ... and a fragment between two of them lands between.
 * A candidate's hop count and DIRECTION to each placed chart then follow from
 * the coordinate difference, consistently for every anchor, so no wrong local
 * side vote and no chain of estimated turns can run away.  The hop LENGTH
 * stays local (the measured turn at the placed chart, else an estimate). */
typedef struct AplLEdge { int32_t a, b; double d, w; } AplLEdge;   /* x_b - x_a = d */
/* PLACEMENT-GATED CUT: a candidate whose anchors form two clusters that each fit rigidly
 * but disagree on u by more than this (and a quarter turn) is two wraps joined by a seam no
 * co-location audit can see (the two parts do not overlap in the component's own frame);
 * it is cut along the seam relations between the charts the two clusters touch, and both
 * parts go back to the candidate pool.  Measured motive 2026-09-08: 9-10% of the centre's
 * and 21x5x5's area sat in reverted/refused multi-cluster components. */
#define APL_CUT_MIN_SEP 150.0
#define APL_PLACE_BUDGET_SEC 900.0   /* the placement loop's wall-clock budget for a 100-cube pile (the assembly budget is 30 min) */
/* ... scaled with the pile (2026-09-16): the 21x5x5 (525 cubes) and pherc343 (1,000 cubes) both hit the flat
 * 900 s (after 97 and 4,513 sweeps) and left 67% / 41.5% of their area never attempted or deferred at the end. */
#define APL_PLACE_BUDGET_SEC_PER_CUBE 9.0
/* A HOP candidate whose rigid pose has fewer witnesses than this in support is not sent to the joint
 * continuity repair.  2026-09-16 on the 21x5x5: 94% of the sweep loop's time was those repairs
 * (743 of 794 s of continuity solves, up to 467k active vertices), and repairs from zero initial
 * support ended supported 9 times in 85 (11%; pherc343: 74 in 1,753, 4%) while the pile's budget
 * cut the loop at 97 sweeps.  Seam candidates keep their repair: their anchors are exact. */
#define APL_HOP_REPAIR_MIN_SUPPORT 1
/* The boundary UV relax trials (local boundary UV repair after a joint repair, and the supported
 * boundary correction after a seam placement) run only for candidates up to this many vertices.
 * 2026-09-16 on the 21x5x5: each pass over a 447k-vertex candidate cost ~30 s (solve 12 s, the
 * guarded trials' checks 18 s), 3,000 of the placement loop's 4,725 s; a bigger candidate keeps its
 * rigid pose, every audit still runs on it, and the repair stage closes its seams later. */
#define APL_RELAX_TRIAL_MAX_VERTICES 100000
#define APL_HOP_REPAIR_MAX_VERTICES 50000   /* ... and only for candidates larger than this: a small candidate's repair is cheap (the control keeps its recovery) */
static double apl_budget_sec(const AsmRun *run)
{
    double scaled = APL_PLACE_BUDGET_SEC_PER_CUBE * (double)(run ? run->n_cubes : 0);
    return scaled > APL_PLACE_BUDGET_SEC ? scaled : APL_PLACE_BUDGET_SEC;
}
#define APL_SWEEPS_LEGACY_CAP 0      /* A/B: 1 = the flat 512-sweep cap (reviewer bug 1, 2026-09-09) */
/* HOP RETRY.  On the 10x10x10 (dev48) 2,375 of ~3,000 fitted placements were reverted by
 * the audit with a median contradiction of 84% of the candidate's area: the piece sits ON
 * material one layer away, i.e. its hop count is off by one, not its evidence.  When a fit
 * is reverted with a contradiction above APL_HOP_RETRY_MIN_RATIO, the same fit is tried one
 * turn inward and one turn outward (the cluster's mean one-hop vector) under the same audits.
 * 0 = off (A/B, 2026-09-09). */
#define APL_HOP_RETRY 0   /* DEAD (measured 2026-09-09): unconditional retries placed 212 of 4,521 on the 10x10x10 (coverage 53.6 -> 73.9%)
                          * but 189 of them had no rim neighbour at all and the stacked (u,z) cells went 43k -> 118k, cross-wrap 484 -> 1,018;
                          * with the rim confirmation required, 0 of 4,150 pass (1,046 fail the rim audit, ~2,900 land on material
                          * in both directions: those pieces are duplicated material, not misplaced).  Costs 150-280 s of probes. */
#define APL_HOP_RETRY_MIN_RATIO 0.5
#define APL_HOP_RETRY_REQUIRE_RIM 1   /* a retry is accepted only on a positive rim audit (>= 8 rim vertices, none violating) */
/* seam anchors from WEAK (gate-rejected, rigid) relations between placed charts and a candidate */
#define APL_SEAM_ANCHORS 1
/* 1 = every layer hit is one wrap, the layer coordinate gives only its direction */
/* measured turns divided by the pair's layer-coordinate hop count (one wrap per measured turn) */
#define APL_HOP_ONE_WRAP 0   /* measured DEAD 2026-09-09: 10x10x10 coverage 60.7 -> 42.3% (the coordinate's two-hop calls
                              * there are physical: at the core the 6-9 vox spacing lets a probe cross a hole to the wrap
                              * two layers out), familiar/centre/21x5x5 within a point */
/* rms gate floor as a fraction of the anchors' mean hop length (0 = off) */
#define APL_GATE_TURN_FRAC 0.05
/* borrowed turns (a placed chart without a measured turn takes the nearest measured one within
 * 0.75 of its length).  MEASURED DEAD 2026-09-08: familiar 88.4 -> 31.6% placed, centre 73.8 ->
 * 54.8%: the nearest measured turn is often a different wrap's, and a wrong hop length is worse
 * than the curvature estimate.  Kept in-tree, off. */
/* v2 (2026-09-09): donors on the SAME wrap only (layer coordinate within 0.5), nearest along u,
 * within this fraction of the donor's turn length.  Interpolates the measured turns along a
 * wrap into the gaps a candidate fills. */
/* hop length at a placed chart without a measured turn: 1 = t0 + layer(P) x growth (the layer
 * coordinate is per chart), 0 = the candidate's stretch/curvature estimate, else t0 */
#define APL_HOP_LAYER_MODEL 1
#define APL_MAX_CUTS 512
/* a component is a SHEET (usable as the frame of reference, its z prior trusted) when its
 * fitted |grad z| over (u,v) is at least this.  0.7 was set on the familiar box (vertical
 * axis, big components read 0.95-1.02); at the centre crop the axis is 37 deg off vertical and
 * the five biggest components (46 percent of the box) read 0.56-0.68, so the primary fell to a
 * 4 percent fragment and no turn was ever measurable (2026-09-08). */
#define APL_SHEET_GRAD 0.7
#define APL_LAYER_SEAM_W 0.0   /* weight of a seam edge (u difference in turns) in the coordinate solve */
#define APL_LAYER_SEAM_ZERO_W_CONNECTS 0   /* A/B: 1 = zero-weight seam edges still connect the union-find (reviewer bug 2) */

/* the local turn length at chart c: its measured outward turn, else its inward one grown, else the model at its layer */
static double apl_turn_len_at(const AplTurn *turn, int32_t c, int8_t out_side, double t0, double growth, double layer_guess)
{
    if (turn != NULL && out_side != 0) {
        const AplTurn *to = &turn[(size_t)c * 2 + (out_side > 0 ? 1 : 0)], *ti = &turn[(size_t)c * 2 + (out_side > 0 ? 0 : 1)];
        if (to->n != 0) return hypot(to->tx, to->ty);
        if (ti->n != 0) return hypot(ti->tx, ti->ty) + growth;
    }
    double t = t0 + layer_guess * growth;
    if (t < APL_TURN_MIN) t = APL_TURN_MIN;
    if (t > APL_TURN_MAX) t = APL_TURN_MAX;
    return t;
}

static void apl_layer_index(AsmRun *run, const int8_t *out_side, const AplTurn *turn, double t0, double growth,
                            int32_t anchor_chart, double *lidx, const double *kover,
                            size_t *n_seam_edges, size_t *n_layer_edges, size_t *n_known, double *res_p50, double *res_p95)
{
    Arena_T arena = run->arena;
    Arena_Mark mark = Arena_save(arena);
    size_t nc = run->n_charts;
    for (size_t i = 0; i < nc; i++) lidx[i] = APL_LAYER_UNKNOWN;
    *n_seam_edges = 0; *n_layer_edges = 0; *n_known = 0; *res_p50 = 0.0; *res_p95 = 0.0;
    size_t cap = run->n_rels + run->n_layers + 1, ne = 0;
    AplLEdge *E = ARENA_ALLOC(arena, cap * sizeof(AplLEdge));
    size_t seam_first = 0;
    for (size_t r = 0; r < run->n_rels; r++) {
        const AsmRelation *rl = &run->rels[r];
        if (rl->flags & (ASM_REL_DROPPED | ASM_REL_SWITCHED | ASM_REL_CONTACT | ASM_REL_PLACEMENT_ONLY)) continue;
        const AsmChart *A = &run->charts[(size_t)rl->a], *B = &run->charts[(size_t)rl->b];
        if (!AsmChart_in_layout(A) || !AsmChart_in_layout(B)) continue;
        if (rl->robust_w < 0.5) continue;
        /* a seam joins two charts of one component whose relative u is fixed by the component's
         * own layout: their layer difference is that u difference in TURNS.  A seam that said
         * "same layer" would pull every turn of a spiral onto one coordinate (measured 2026-09-08:
         * the primary's ramp compressed to 0.88 of u/T, p90 error 1.08 layers, hops from it wrong).
         * The turn length is local (pass 1: measured or t0; pass 2: the model at the solved layer). */
        double ax, ay, bx, by;
        apl_pose(A, 0.0, 0.0, &ax, &ay);
        apl_pose(B, 0.0, 0.0, &bx, &by);
        (*n_seam_edges)++;
        if (APL_LAYER_SEAM_W <= 0.0 && !APL_LAYER_SEAM_ZERO_W_CONNECTS) continue;   /* a zero-weight edge would still connect the union-find and mark charts 'known' (reviewer, 2026-09-09) */
        E[ne].a = rl->a; E[ne].b = rl->b; E[ne].d = bx - ax; E[ne].w = APL_LAYER_SEAM_W; ne++;   /* d = du for now */
    }
    size_t seam_count = ne;
    double dg = run->layer_d_global > 0.0 ? run->layer_d_global : 10.0;
    for (size_t p = 0; p < run->n_layers; p++) {
        const AsmLayerPair *lp = &run->layers[p];
        if (!AsmChart_in_layout(&run->charts[(size_t)lp->a]) || !AsmChart_in_layout(&run->charts[(size_t)lp->b])) continue;
        if (out_side[lp->a] == 0 || lp->count <= 0) continue;
        /* MEASURED on the familiar's audited layout (2026-09-08): 97% of layer pairs are the
         * ADJACENT wrap whatever the hit distance (8-26 vox is papyrus thickness and gap noise,
         * not wrap count), so a hit is one layer; the distance ratio would call a third of the
         * d >= 20 hits two layers and it is wrong.  13% of the local votes point the wrong way
         * (a mis-oriented chart votes wrong on every hit): the reweighting silences them. */
        (void)dg;
        int k = APL_PAIR_K_FROM_PROBE ? (lp->k >= 1 ? (lp->k > 3 ? 3 : lp->k) : 1) : 1;   /* the probe's crossing order */
        double sgn = lp->side == out_side[lp->a] ? 1.0 : -1.0;
        E[ne].a = lp->a; E[ne].b = lp->b; E[ne].d = sgn * (double)k;
        if (kover != NULL && kover[p] != 0.0) { E[ne].d = kover[p]; E[ne].w = 4.0 * (lp->count < 20 ? (double)lp->count : 20.0) / 4.0; ne++; (*n_layer_edges)++; continue; }
        E[ne].w = (lp->count < 20 ? (double)lp->count : 20.0) / 4.0; ne++;
        (*n_layer_edges)++;
    }
    /* connectivity: the anchor is the primary's chart inside the LARGEST connected set of edges
     * (by the primary's area in it), never simply its biggest chart -- measured 2026-09-08 on the
     * no-winding centre pile: the biggest chart had no layer edge and 1 of 856 charts was reached */
    UnionFind uf = UF_new(arena, (int32_t)nc);
    for (size_t e = 0; e < ne; e++) uf_union(&uf, E[e].a, E[e].b);
    {
        int32_t pcomp = run->charts[(size_t)anchor_chart].component;
        double *set_area = ARENA_CALLOC(arena, nc, sizeof(double));
        for (size_t i = 0; i < nc; i++) {
            const AsmChart *c = &run->charts[i];
            if (!AsmChart_in_layout(c) || c->component != pcomp) continue;
            set_area[uf_find(&uf, (int32_t)i)] += c->area3d;
        }
        int32_t best_root = -1;
        for (size_t i = 0; i < nc; i++) if (best_root < 0 || set_area[i] > set_area[best_root]) best_root = (int32_t)i;
        if (best_root >= 0 && set_area[best_root] > 0.0 && uf_find(&uf, anchor_chart) != best_root) {
            int32_t pick = -1;
            for (size_t i = 0; i < nc; i++) {
                const AsmChart *c = &run->charts[i];
                if (!AsmChart_in_layout(c) || c->component != pcomp || uf_find(&uf, (int32_t)i) != best_root) continue;
                if (pick < 0 || c->area3d > run->charts[(size_t)pick].area3d) pick = (int32_t)i;
            }
            if (pick >= 0) anchor_chart = pick;
        }
    }
    int32_t root = uf_find(&uf, anchor_chart);
    uint8_t *known = ARENA_CALLOC(arena, nc, 1);
    for (size_t i = 0; i < nc; i++) if (AsmChart_in_layout(&run->charts[i]) && uf_find(&uf, (int32_t)i) == root) { known[i] = 1; (*n_known)++; }
    /* weighted least squares by conjugate gradients on the anchored Laplacian, IRLS outside */
    double *x = ARENA_CALLOC(arena, nc, sizeof(double)), *rhs = ARENA_ALLOC(arena, nc * sizeof(double));
    double *rr = ARENA_ALLOC(arena, nc * sizeof(double)), *pp = ARENA_ALLOC(arena, nc * sizeof(double)), *Ap = ARENA_ALLOC(arena, nc * sizeof(double));
    double *w = ARENA_ALLOC(arena, (ne ? ne : 1) * sizeof(double)), *res = ARENA_ALLOC(arena, (ne ? ne : 1) * sizeof(double));
    double *seam_du = ARENA_ALLOC(arena, (seam_count ? seam_count : 1) * sizeof(double));
    for (size_t e = seam_first; e < seam_count; e++) seam_du[e] = E[e].d;
    for (size_t e = 0; e < ne; e++) w[e] = E[e].w;
    for (int round = 0; round < 8; round++) {
        /* seam differences in turns at the local turn length; pass 2 (rounds 4-7) reads the layer solved so far */
        for (size_t e = seam_first; e < seam_count; e++) {
            double lg = round >= 4 && known[E[e].a] ? x[E[e].a] : 0.0;
            double T = apl_turn_len_at(turn, E[e].a, out_side[E[e].a], t0, growth, lg);
            E[e].d = seam_du[e] / T;
        }
        if (round == 4) for (size_t e = 0; e < ne; e++) w[e] = E[e].w;   /* fresh robust weights for pass 2 */
        memset(rhs, 0, nc * sizeof(double));
        for (size_t e = 0; e < ne; e++) {
            if (!known[E[e].a] || !known[E[e].b]) continue;
            rhs[E[e].b] += w[e] * E[e].d; rhs[E[e].a] -= w[e] * E[e].d;
        }
        rhs[anchor_chart] = 0.0;
        /* r = rhs - A x */
        for (size_t i = 0; i < nc; i++) rr[i] = 0.0;
        for (size_t e = 0; e < ne; e++) {
            if (!known[E[e].a] || !known[E[e].b]) continue;
            double t = w[e] * (x[E[e].b] - x[E[e].a]);
            rr[E[e].b] += t; rr[E[e].a] -= t;
        }
        rr[anchor_chart] = x[anchor_chart];
        for (size_t i = 0; i < nc; i++) rr[i] = known[i] ? rhs[i] - rr[i] : 0.0;
        memcpy(pp, rr, nc * sizeof(double));
        double rs = 0.0;
        for (size_t i = 0; i < nc; i++) rs += rr[i] * rr[i];
        double rs0 = rs;
        size_t maxit = nc + 50 < 4000 ? nc + 50 : 4000;
        for (size_t it = 0; it < maxit && rs > 1e-14 * (rs0 > 0.0 ? rs0 : 1.0); it++) {
            for (size_t i = 0; i < nc; i++) Ap[i] = 0.0;
            for (size_t e = 0; e < ne; e++) {
                if (!known[E[e].a] || !known[E[e].b]) continue;
                double t = w[e] * (pp[E[e].b] - pp[E[e].a]);
                Ap[E[e].b] += t; Ap[E[e].a] -= t;
            }
            Ap[anchor_chart] = pp[anchor_chart];
            double pAp = 0.0;
            for (size_t i = 0; i < nc; i++) if (known[i]) pAp += pp[i] * Ap[i];
            if (!(pAp > 0.0)) break;
            double alpha = rs / pAp;
            for (size_t i = 0; i < nc; i++) if (known[i]) { x[i] += alpha * pp[i]; rr[i] -= alpha * Ap[i]; }
            double rs2 = 0.0;
            for (size_t i = 0; i < nc; i++) rs2 += rr[i] * rr[i];
            double beta = rs2 / rs; rs = rs2;
            for (size_t i = 0; i < nc; i++) pp[i] = known[i] ? rr[i] + beta * pp[i] : 0.0;
        }
        /* the hop count of a layer edge is an UNKNOWN the solve decides (measured 2026-09-09: inside
         * the 10x10x10's primary the hits sit at ~0 or ~2 turns of displacement for half the edges
         * while every edge was counted one hop): after each robust round a layer edge takes the
         * integer nearest its current coordinate difference, 1 to 3 wraps, keeping its sign */
        if (APL_LAYER_KINT && round >= 1) {
            for (size_t e = seam_count; e < ne; e++) {
                if (!known[E[e].a] || !known[E[e].b]) continue;
                double dl = fabs(x[E[e].b] - x[E[e].a]);
                int k = (int)floor(dl + 0.5);
                if (k < 1) k = 1;
                if (k > 3) k = 3;
                E[e].d = (E[e].d < 0.0 ? -1.0 : 1.0) * (double)k;
            }
        }
        /* reweight: an edge off by half a layer or more is a contradiction */
        for (size_t e = 0; e < ne; e++) {
            double r = (known[E[e].a] && known[E[e].b]) ? x[E[e].b] - x[E[e].a] - E[e].d : 0.0;
            res[e] = fabs(r);
            w[e] = E[e].w / (1.0 + (r / 0.5) * (r / 0.5));
        }
    }
    for (size_t i = 0; i < nc; i++) lidx[i] = known[i] ? x[i] - x[anchor_chart] : APL_LAYER_UNKNOWN;
    if (ne > 0) {
        qsort(res, ne, sizeof(double), apl_cmp_double);
        *res_p50 = res[ne / 2]; *res_p95 = res[(ne * 95) / 100];
    }
    Arena_restore(arena, mark);
}

/* Collect the evidence for component `best`: every layer hit between a
 * placed chart and the component, forward (placed chart probing) and
 * reverse (the component probing placed material), for k = 1..K. */
static int apl_placement_seam(const AsmRelation *r, int framed, int veto)
{
    /* Cleaning/parity drops and measured contacts remain vetoes. A switched
     * relation is only a placement hypothesis, never readmitted as a join:
     * its boundary fit and the current layout audits must still confirm it. */
    if (r->flags & (ASM_REL_DROPPED | ASM_REL_CONTACT)) return 0;
    /* The unimodality gate says these correspondences do not describe one curve under one
     * rigid motion.  That is negative evidence about where the two charts sit, not an anchor.
     * At ASM_WRAP_VETO_ANCHOR 2 it binds only on an UNFRAMED pile: where the axis frame is
     * armed the layout has a trustworthy geometric model and the audits police a bad pose
     * anyway, and the 4x5x5's 70 such seams are load-bearing evidence (2026-09-18). */
    if ((r->flags & ASM_REL_CROSSWRAP) &&
        (veto == 1 || (veto == 2 && !framed))) return 0;
    if (r->flags & ASM_REL_PLACEMENT_ONLY) return 1;
    return APL_SWITCHED_SEAM_ANCHORS && (r->flags & ASM_REL_SWITCHED) &&
           r->corr_count >= APL_SEAM_MIN_ANCHORS && r->rms <= 4.0;
}

/* Crossing order counts surfaces actually observed by the ray. A missing
 * intermediate surface makes it a lower bound, not an exact wrap count.
 * The measured gap and nearby spacing bound extra hypotheses; they do not
 * certify any of those hypotheses. The competing poses still face the
 * same source, overlap, rim and uniqueness checks. */
static int apl_hop_upper(const AplCtx *x, const AsmLayerPair *lp)
{
    int upper = lp->k > 0 ? lp->k : 1;
    double pitch = x->run->layer_d_global;
    if (x->run->chart_layer_d) for (int side = 0; side < 2; side++) {
        double d = x->run->chart_layer_d[side ? lp->b : lp->a];
        if (d > 0 && (!(pitch > 0) || d < pitch)) pitch = d;
    }
    if (pitch > 0 && lp->d_median > 0) {
        double lower_pitch = pitch*fmax(0.5, 1.0-x->o->layer_tol);
        int possible = (int)fmin((double)x->o->max_layers, floor(lp->d_median/lower_pitch));
        if (possible > upper) upper = possible;
    } else upper = x->o->max_layers;
    return upper < x->o->max_layers ? upper : x->o->max_layers;
}

static size_t apl_collect_anchors(AplCtx *x, int32_t best)
{
    AsmRun *run = x->run;
    const AsmPlaceOpts *o = x->o;
    size_t na = 0;
    x->na_seam = 0;
    x->n_observations = 0;
    x->turn_rel = 0.0;
    x->fit_cu = x->fit_cv = x->fit_s = x->fit_area = 0;
    if (APL_V_FROM_AXIS && x->vm_ok && x->vm_pin && x->cu && x->cv && x->s_axis) {
        for (size_t i = 0; i < run->n_charts; i++) {
            const AsmChart *c = &run->charts[i];
            if (!AsmChart_in_layout(c) || c->component != best) continue;
            double u, v; apl_layout_centroid(x, c, i, 0, 0, 0, &u, &v);
            x->fit_area += c->area3d; x->fit_cu += c->area3d*u; x->fit_cv += c->area3d*v;
            x->fit_s += c->area3d*x->s_axis[i];
        }
        if (x->fit_area > 0) { x->fit_cu /= x->fit_area; x->fit_cv /= x->fit_area; x->fit_s /= x->fit_area; }
    }
    /* SEAM ANCHORS: a gate-rejected but rigidly fitting seam between a placed chart and the
     * candidate gives exact (u,v) targets for the candidate's boundary vertices; the fit, the
     * gates and the audits judge them like any other anchor (weight 1, hop 0). */
    if (APL_SEAM_ANCHORS) for (size_t r = 0; r < run->n_rels; r++) {
        const AsmRelation *R = &run->rels[r];
        if (!apl_placement_seam(R, x->have_radius, x->o->crosswrap_anchor_veto)) continue;
        const AsmChart *A = &run->charts[(size_t)R->a], *B = &run->charts[(size_t)R->b];
        if (!AsmChart_in_layout(A) || !AsmChart_in_layout(B)) continue;
        int fwd = x->comp_placed[A->component] == 1 && B->component == best;
        int rev = A->component == best && x->comp_placed[B->component] == 1;
        if (!fwd && !rev) continue;
        if (na + (size_t)R->corr_count > x->hcap) break;
        for (int32_t k = 0; k < R->corr_count; k++) {
            const AsmCorr *cr = &run->corr[(size_t)R->corr_first + (size_t)k];
            double target[2];
            if (!AsmContinuity_target(run, R, cr, rev, target)) continue;
            double ax, ay, bx, by;
            apl_vertex_point(A, cr->va, &ax, &ay);
            apl_vertex_point(B, cr->vb, &bx, &by);
            AplAnchor *an = &x->anchors[na++];
            if (fwd) { an->sx = bx; an->sy = by; an->tx = ax; an->ty = ay; an->chart_p = R->a; an->chart_c = R->b; }
            else     { an->sx = ax; an->sy = ay; an->tx = bx; an->ty = by; an->chart_p = R->b; an->chart_c = R->a; }
            an->tx = target[0]; an->ty = target[1];
            an->source_flip = (((A->flags ^ B->flags) & ASM_CHART_MIRROR) != ((R->flags & ASM_REL_PARITY) ? ASM_CHART_MIRROR : 0));
            if (an->source_flip) an->sx = -an->sx;
            an->relation = (int32_t)r; an->source_run = cr->run;
            an->hit = NULL;
            an->measured_order = 0; an->measured_turn = 0;
            an->w = 1.0; an->k = 0; an->estimated = 0; an->reverse = rev; an->side = 0; an->Tx = 0.0; an->Ty = 0.0; an->cluster = -1;
            an->turn_src = APL_TSRC_SEAM;
            an->seam_rms = R->rms; x->na_seam++;
        }
    }
    for (size_t p = 0; p < run->n_layers; p++) {
        const AsmLayerPair *lp = &run->layers[p];
        const AsmChart *A = &run->charts[(size_t)lp->a], *B = &run->charts[(size_t)lp->b];
        if (!AsmChart_in_layout(A) || !AsmChart_in_layout(B)) continue;
        int forward = x->comp_placed[A->component] == 1 && B->component == best;
        int reverse = A->component == best && x->comp_placed[B->component] == 1;
        if (!forward && !reverse) continue;
        /* the placed chart's turn toward the component: forward = A's own
         * side; reverse = B's side facing A (B is on A's side lp->side; A is
         * outward of B iff B is inward of A) */
        int32_t P = forward ? lp->a : lp->b;
        int side;
        if (forward) side = lp->side;
        else {
            if (x->out_side[lp->a] == 0 || x->out_side[lp->b] == 0) continue;
            int b_outward_of_a = lp->side == x->out_side[lp->a];
            side = b_outward_of_a ? -x->out_side[lp->b] : x->out_side[lp->b];
        }
        size_t observation = x->n_observations;
        x->n_observations += (size_t)lp->hit_count;
        int first_hop = lp->k > 0 ? lp->k : 1, last_hop = apl_hop_upper(x, lp);
        if (x->out_side[P] != 0) {
          for (int kk = first_hop; kk <= last_hop; kk++) {
            int32_t dk = side == x->out_side[P] ? kk : -kk;
            if (kk > o->max_layers || na + (size_t)lp->hit_count > x->hcap) continue;
            int hop_side = side;
            double tx = 0.0, ty = 0.0;
            int tsrc = apl_turn_k(x->turn, P, hop_side, x->out_side[P], x->growth, kk, &tx, &ty);
            int measured_turn = tsrc == 1;
            int tsource = tsrc == 1 ? (kk > 1 ? APL_TSRC_MEASURED_K : APL_TSRC_MEASURED)
                        : tsrc == 2 ? APL_TSRC_OPPOSITE : APL_TSRC_T0;
            double measured_tx = tx, measured_ty = ty;
            /* Preserve both a measured turn and the one-hop estimate from
             * a measured opposite turn. On a noncircular section a local
             * radius is not the circumference; replacing the opposite-side
             * estimate makes one coherent wrap vote for different poses at
             * different azimuths. Both fallback estimates retain estimated
             * weight and provenance. */
            if (!tsrc && apl_turn_radius(x, P, dk > 0, kk, &tx, &ty)) { tsrc = 2; tsource = APL_TSRC_RADIUS; x->anchors_radius += (size_t)lp->hit_count; }
            if (!tsrc) {
                /* no measured turn at P: a local estimate of one turn, grown per hop.  The
                 * LAYER MODEL reads P's own layer coordinate: T = t0 + layer(P) x growth, the
                 * primary's measured turn at its anchor grown per wrap (the stretch estimate is
                 * the candidate's curvature, noisy on a chart-sized arc). */
                int use_layer = APL_HOP_LAYER_MODEL && x->t0_measured && x->lidx != NULL && x->lidx[P] > APL_LAYER_UNKNOWN;
                double t1 = use_layer ? x->t0 + x->lidx[P] * x->growth
                          : (x->stretch_ok ? fabs(x->stretch_T) : x->t0);
                tsource = use_layer ? APL_TSRC_LAYER : (x->stretch_ok ? APL_TSRC_STRETCH : APL_TSRC_T0);
                if (!x->t0_measured && !x->stretch_ok)
                    t1 += 0.5*x->growth*(dk > 0 ? 1.0 : -1.0); /* integrate between the two measured radii */
                if (APL_TURN_FIELD && x->tfield != NULL && x->tfield[P] > 0.0) t1 = x->tfield[P];   /* the local observed turn */
                /* how well this fallback knows its turn: the layer model rests on a MEASURED t0,
                 * the stretch on its own standard error, and a bare unmeasured t0 on nothing */
                {
                    double arel = use_layer ? 0.0
                                : (x->stretch_ok ? (x->stretch_rel >= 0.0 ? x->stretch_rel : 1e300)
                                : (x->t0_measured ? 0.0 : 1e300));
                    if (arel > x->turn_rel) x->turn_rel = arel;
                }
                if (t1 < APL_TURN_MIN) t1 = APL_TURN_MIN;
                if (t1 > APL_TURN_MAX) t1 = APL_TURN_MAX;
                double sum = 0.0;
                for (int j = 0; j < kk; j++) sum += t1 + (double)j * x->growth * (dk > 0 ? 1.0 : -1.0);
                tx = dk > 0 ? sum : -sum; ty = 0.0; tsrc = 2;
            }
            for (int h = 0; h < lp->hit_count; h++) {
                const AsmLayerHit *hit = &run->layer_hits[(size_t)lp->hit_first + (size_t)h];
                double ax, ay, bx, by;
                apl_vertex_point(A, hit->va, &ax, &ay);
                apl_face_point(B, hit->fb, hit->l0, hit->l1, &bx, &by);
                AplAnchor *an = &x->anchors[na++];
                if (forward) { an->sx = bx; an->sy = by; an->tx = ax + tx; an->ty = ay + ty; }
                else         { an->sx = ax; an->sy = ay; an->tx = bx + tx; an->ty = by + ty; }
                an->w = (tsrc == 1 ? 1.0 : 0.5) / (double)kk;
                an->k = kk; an->estimated = tsrc != 1; an->reverse = reverse; an->turn_src = tsource;
                an->chart_p = P; an->side = (int8_t)side; an->Tx = tx; an->Ty = ty; an->cluster = -1;
                an->chart_c = forward ? lp->b : lp->a; an->seam_rms = 0.0;
                an->relation = -1; an->source_run = -1; an->source_flip = 0;
                an->observation = observation+(size_t)h;
                an->hit = hit;
                an->measured_order = hit->order; an->measured_turn = measured_turn;
                an->measured_x = (forward ? ax : bx) + measured_tx;
                an->measured_y = (forward ? ay : by) + measured_ty;
            }
          }
            continue;
        }
        for (int k = first_hop; k <= last_hop && na + (size_t)lp->hit_count <= x->hcap; k++) {
            double tx, ty;
            int tsrc = apl_turn_k(x->turn, P, side, x->out_side[P], x->growth, k, &tx, &ty);
            int measured_turn = tsrc == 1;
            double measured_tx = tsrc ? tx : 0.0, measured_ty = tsrc ? ty : 0.0;
            int tsource = tsrc == 1 ? (k > 1 ? APL_TSRC_MEASURED_K : APL_TSRC_MEASURED)
                        : tsrc == 2 ? APL_TSRC_OPPOSITE : APL_TSRC_T0;
            if (!tsrc && k == 1 && x->stretch_ok) {
                /* no turn vector anywhere near P: the pooled stretch estimate; the
                 * candidate is on the side the stretch says, whatever the probe side */
                tx = x->stretch_T; ty = 0.0; tsrc = 2; tsource = APL_TSRC_STRETCH;
                double arel = x->stretch_rel >= 0.0 ? x->stretch_rel : 1e300;
                if (arel > x->turn_rel) x->turn_rel = arel;
            }
            if (!tsrc && x->out_side[P] != 0 && apl_turn_radius(x, P, side == x->out_side[P], k, &tx, &ty)) { tsrc = 2; tsource = APL_TSRC_RADIUS; x->anchors_radius += (size_t)lp->hit_count; }
            if (!tsrc) break;
            for (int h = 0; h < lp->hit_count; h++) {
                const AsmLayerHit *hit = &run->layer_hits[(size_t)lp->hit_first + (size_t)h];
                double ax, ay, bx, by;
                apl_vertex_point(A, hit->va, &ax, &ay);
                apl_face_point(B, hit->fb, hit->l0, hit->l1, &bx, &by);
                AplAnchor *an = &x->anchors[na++];
                if (forward) { an->sx = bx; an->sy = by; an->tx = ax + tx; an->ty = ay + ty; }
                else         { an->sx = ax; an->sy = ay; an->tx = bx + tx; an->ty = by + ty; }
                an->w = (tsrc == 1 ? 1.0 : 0.5) / (double)k;
                an->k = k; an->estimated = tsrc == 2; an->reverse = reverse; an->turn_src = tsource;
                an->chart_p = P; an->side = (int8_t)side; an->Tx = tx; an->Ty = ty; an->cluster = -1;
                an->chart_c = forward ? lp->b : lp->a; an->seam_rms = 0.0;
                an->relation = -1; an->source_run = -1; an->source_flip = 0;
                an->observation = observation+(size_t)h;
                an->hit = hit;
                an->measured_order = hit->order; an->measured_turn = measured_turn;
                an->measured_x = (forward ? ax : bx) + measured_tx;
                an->measured_y = (forward ? ay : by) + measured_ty;
            }
        }
    }
    for (size_t i = 0; i < na; i++) {
        AplAnchor *a = &x->anchors[i];
        const AsmChart *c = &run->charts[a->chart_c];
        double sx = a->source_flip ? -a->sx : a->sx;
        double dx = sx-c->pose_x, dy = a->sy-c->pose_y, ct = cos(c->pose_theta), sn = sin(c->pose_theta);
        a->local_u = ct*dx+sn*dy; a->local_v = -sn*dx+ct*dy;
        if (c->flags & ASM_CHART_MIRROR) a->local_u = -a->local_u;
    }
    return na;
}

typedef struct AplDuIdx { double du; size_t idx; } AplDuIdx;
static int apl_cmp_duidx(const void *p, const void *q)
{
    const AplDuIdx *a = p, *b = q;
    if (a->du != b->du) return a->du < b->du ? -1 : 1;
    return (a->idx > b->idx) - (a->idx < b->idx);
}

/* Range-add / maximum tree over possible lower V edges of a window. Each
 * active point contributes to every V window containing it. Sweeping the U
 * edge finds the heaviest rectangle in O(n log n), without averaging two
 * separated V groups into an empty gap. */
static void apl_window_add(double *value, double *lazy, size_t node, size_t lo, size_t hi,
                           size_t first, size_t last, double weight)
{
    if (first >= hi || last <= lo) return;
    if (first <= lo && hi <= last) { value[node] += weight; lazy[node] += weight; return; }
    size_t mid = lo+(hi-lo)/2;
    apl_window_add(value, lazy, node*2, lo, mid, first, last, weight);
    apl_window_add(value, lazy, node*2+1, mid, hi, first, last, weight);
    value[node] = lazy[node]+fmax(value[node*2], value[node*2+1]);
}

static size_t apl_window_edge(const double *value, size_t n)
{
    size_t node = 1, lo = 0, hi = n;
    while (hi-lo > 1) {
        size_t mid = lo+(hi-lo)/2;
        if (value[node*2] >= value[node*2+1]) { node *= 2; hi = mid; }
        else { node = node*2+1; lo = mid; }
    }
    return lo;
}

static size_t apl_window_bound(const double *v, size_t n, double key, int upper)
{
    size_t lo = 0, hi = n;
    while (lo < hi) {
        size_t mid = lo+(hi-lo)/2;
        if (v[mid] < key || (upper && v[mid] == key)) lo = mid+1; else hi = mid;
    }
    return lo;
}

static int apl_cluster_eligible(const AplAnchor *a, int seam_only)
{
    return a->cluster < 0 && a->w > 0 &&
           (seam_only ? a->k == 0 && a->source_flip == seam_only-1 : a->k > 0);
}

/* Greedy disjoint rectangles: U width wu, V width 2*wv (the existing V
 * tolerance). Assign precisely the winning rectangle, never a second,
 * wider window around its centroid. */
static int apl_cluster_pass(AplCtx *x, size_t na, const size_t *order, const double *du, double wu, double wv, int rounds, int seam_only, int kind,
                            int ncl, double *cl_w, double *cl_u, double *cl_v, int *cl_kind)
{
    Arena_Mark mark = Arena_save(x->run->arena);
    double *v = ARENA_ALLOC(x->run->arena, (na ? na : 1)*sizeof(double));
    size_t nv = 0;
    for (size_t i = 0; i < na; i++) if (apl_cluster_eligible(&x->anchors[i], seam_only))
        v[nv++] = x->anchors[i].ty-x->anchors[i].sy;
    if (!nv) { Arena_restore(x->run->arena, mark); return ncl; }
    qsort(v, nv, sizeof(double), apl_cmp_double);
    size_t unique = 0;
    for (size_t i = 0; i < nv; i++) if (!unique || v[i] != v[unique-1]) v[unique++] = v[i];
    nv = unique;
    size_t *first = ARENA_ALLOC(x->run->arena, na*sizeof(size_t));
    size_t *last = ARENA_ALLOC(x->run->arena, na*sizeof(size_t));
    for (size_t i = 0; i < na; i++) if (apl_cluster_eligible(&x->anchors[i], seam_only)) {
        double y = x->anchors[i].ty-x->anchors[i].sy;
        first[i] = apl_window_bound(v, nv, y-2*wv, 0);
        last[i] = apl_window_bound(v, nv, y, 1);
    }
    double *value = ARENA_ALLOC(x->run->arena, 4*nv*sizeof(double));
    double *lazy = ARENA_ALLOC(x->run->arena, 4*nv*sizeof(double));
    for (int round = 0; round < rounds; round++) {
        memset(value, 0, 4*nv*sizeof(double)); memset(lazy, 0, 4*nv*sizeof(double));
        double best_w = 0.0, bottom = 0.0; size_t bi = 0, bj = 0, j = 0;
        for (size_t i = 0; i < na; i++) {
            while (j < na && du[order[j]]-du[order[i]] <= wu) {
                size_t p = order[j++]; const AplAnchor *a = &x->anchors[p];
                if (apl_cluster_eligible(a, seam_only)) apl_window_add(value, lazy, 1, 0, nv, first[p], last[p], a->w);
            }
            if (value[1] > best_w+1e-9) { best_w = value[1]; bi = i; bj = j; bottom = v[apl_window_edge(value, nv)]; }
            size_t p = order[i]; const AplAnchor *a = &x->anchors[p];
            if (apl_cluster_eligible(a, seam_only)) apl_window_add(value, lazy, 1, 0, nv, first[p], last[p], -a->w);
        }
        if (best_w <= 1e-9) break;
        double w2 = 0.0, su = 0.0, sv = 0.0;
        for (size_t q = bi; q < bj; q++) {
            AplAnchor *a = &x->anchors[order[q]];
            if (!apl_cluster_eligible(a, seam_only)) continue;
            double y = a->ty-a->sy;
            if (y >= bottom && y-bottom <= 2*wv) {
                a->cluster = ncl; w2 += a->w; su += a->w*(a->tx-a->sx); sv += a->w*y;
            }
        }
        if (w2 <= 0.0) break;
        cl_w[ncl] = w2; cl_u[ncl] = su / w2; cl_v[ncl] = sv / w2;
        if (cl_kind) cl_kind[ncl] = kind;
        ncl++;
    }
    Arena_restore(x->run->arena, mark);
    return ncl;
}

static double apl_fit_cluster(AplCtx *x, size_t na, int cl, size_t min_members, double max_rot, double *th, double *tx, double *ty,
                              size_t *kept_out, double *res);

/* Cluster the offsets t - s by weight: greedy windows in u (width a fraction
 * of the turn length) and v.  Returns the number of clusters, labels in
 * anchors[].cluster, ranked 0 = heaviest; with APL_SEAM_FIRST the seam
 * clusters (k == 0 anchors, tight window) come first, cl_kind[] = 1. */
static int apl_cluster_anchors(AplCtx *x, size_t na, double *cl_w, double *cl_u, double *cl_v, int *cl_kind)
{
    Arena_Mark mark = Arena_save(x->run->arena);
    double *du = ARENA_ALLOC(x->run->arena, (na ? na : 1) * sizeof(double));
    size_t *order = ARENA_ALLOC(x->run->arena, (na ? na : 1) * sizeof(size_t));
    double tl = 0.0; size_t ntl = 0;
    for (size_t i = 0; i < na; i++) { du[i] = x->anchors[i].tx - x->anchors[i].sx; if (x->anchors[i].k > 0) { tl += hypot(x->anchors[i].Tx, x->anchors[i].Ty) / (double)x->anchors[i].k; ntl++; } }
    tl = ntl ? tl / (double)ntl : 1000.0;
    double wu = fmax(0.15 * tl, 80.0), wv = 100.0;
    for (size_t i = 0; i < na; i++) order[i] = i;
    if (APL_CLUSTER_QSORT) {
        /* index sort by du, ties by index: the order the stable insertion sort gave, in n log n */
        AplDuIdx *key = ARENA_ALLOC(x->run->arena, (na ? na : 1) * sizeof(AplDuIdx));
        for (size_t i = 0; i < na; i++) { key[i].du = du[i]; key[i].idx = i; }
        if (na > 1) qsort(key, na, sizeof(AplDuIdx), apl_cmp_duidx);
        for (size_t i = 0; i < na; i++) order[i] = key[i].idx;
    } else {
        /* the insertion sort (O(n^2): candidates on the big rungs carry tens of thousands of anchors) */
        for (size_t i = 1; i < na; i++) {
            size_t k = order[i]; size_t j = i;
            while (j > 0 && du[order[j-1]] > du[k]) { order[j] = order[j-1]; j--; }
            order[j] = k;
        }
    }
    int ncl = 0;
    if (APL_SEAM_FIRST && x->na_seam > 0) for (int flip = 0; flip < 2; flip++)
        ncl = apl_cluster_pass(x, na, order, du, APL_SEAM_CLUSTER_VOX, APL_SEAM_CLUSTER_VOX, APL_MAX_SEAM_CLUSTERS, flip+1, 1, ncl, cl_w, cl_u, cl_v, cl_kind);
    int first_hop = ncl;
    ncl = apl_cluster_pass(x, na, order, du, wu, wv, APL_MAX_CLUSTERS, 0, 0, ncl, cl_w, cl_u, cl_v, cl_kind);
    if (first_hop < ncl && x->o) {
        /* A rigid rotation spreads t-s along a long component. Merge whole
         * coherent fragments when t-Rs fits the original rectangle and
         * the combined native fit still meets the witness precision.
         * Leave every other group intact: reframing all observations can
         * fabricate a new, weak alternative from an incoherent wrap. */
        AplCtx fit = *x;
        fit.wts = ARENA_ALLOC(x->run->arena,(na ? na : 1)*sizeof(double));
        double *residual = ARENA_ALLOC(x->run->arena,(na ? na : 1)*sizeof(double));
        int *saved = ARENA_ALLOC(x->run->arena,(na ? na : 1)*sizeof(int));
        for (int seed = first_hop; seed < ncl; seed++) {
            double th,tx,ty; size_t kept;
            double rms = apl_fit_cluster(&fit,na,seed,(size_t)x->o->min_anchors,x->o->max_rot_refine,&th,&tx,&ty,&kept,residual);
            if (!(rms <= 16) || !isfinite(th) || fabs(th) < 1e-8) continue;
            double cu = 0, cv = 0, sw = 0;
            for (size_t i = 0; i < na; i++) if (fit.wts[i] > 0) {
                cu += fit.wts[i]*x->anchors[i].sx; cv += fit.wts[i]*x->anchors[i].sy; sw += fit.wts[i];
            }
            if (!(sw > 0)) continue;
            cu /= sw; cv /= sw;
            double ct = cos(th), sn = sin(th);
            double lo[APL_MAX_CLUSTERS_TOTAL][2], hi[APL_MAX_CLUSTERS_TOTAL][2];
            for (int cl = first_hop; cl < ncl; cl++) for (int d = 0; d < 2; d++) { lo[cl][d] = 1e300; hi[cl][d] = -1e300; }
            for (size_t i = 0; i < na; i++) {
                const AplAnchor *a = &x->anchors[i]; int cl = a->cluster;
                if (cl < first_hop || cl >= ncl) continue;
                double u = a->sx-cu, v = a->sy-cv;
                double offset[2] = {a->tx-a->sx+(1-ct)*u+sn*v, a->ty-a->sy-sn*u+(1-ct)*v};
                for (int d = 0; d < 2; d++) { lo[cl][d] = fmin(lo[cl][d],offset[d]); hi[cl][d] = fmax(hi[cl][d],offset[d]); }
            }
            for (int cl = seed+1; cl < ncl; cl++) {
                if (hi[cl][0] < lo[cl][0] ||
                    fmax(hi[seed][0],hi[cl][0])-fmin(lo[seed][0],lo[cl][0]) > wu ||
                    fmax(hi[seed][1],hi[cl][1])-fmin(lo[seed][1],lo[cl][1]) > 2*wv) continue;
                double oth,otx,oty;
                rms = apl_fit_cluster(&fit,na,cl,(size_t)x->o->min_anchors,x->o->max_rot_refine,&oth,&otx,&oty,&kept,residual);
                if (!(rms <= 16)) continue;
                for (size_t i = 0; i < na; i++) {
                    saved[i] = x->anchors[i].cluster;
                    if (saved[i] == cl) x->anchors[i].cluster = seed;
                }
                rms = apl_fit_cluster(&fit,na,seed,(size_t)x->o->min_anchors,x->o->max_rot_refine,&oth,&otx,&oty,&kept,residual);
                if (!(rms <= 16)) {
                    for (size_t i = 0; i < na; i++) x->anchors[i].cluster = saved[i];
                    continue;
                }
                for (int d = 0; d < 2; d++) { lo[seed][d] = fmin(lo[seed][d],lo[cl][d]); hi[seed][d] = fmax(hi[seed][d],hi[cl][d]); }
            }
        }
        int next = first_hop;
        for (int cl = first_hop; cl < ncl; cl++) {
            double sw = 0, su = 0, sv = 0;
            for (size_t i = 0; i < na; i++) if (x->anchors[i].cluster == cl) {
                AplAnchor *a = &x->anchors[i]; a->cluster = next;
                sw += a->w; su += a->w*(a->tx-a->sx); sv += a->w*(a->ty-a->sy);
            }
            if (sw > 0) { cl_w[next] = sw; cl_u[next] = su/sw; cl_v[next] = sv/sw; cl_kind[next] = 0; next++; }
        }
        ncl = next;
    }
    Arena_restore(x->run->arena, mark);
    return ncl;
}

static double apl_fit_cluster_transform(AplCtx *x, size_t na, double rot, double *th, double *tx, double *ty)
{
    int axis = APL_V_FROM_AXIS && x->vm_ok && x->vm_pin && x->fit_area > 0;
    for (size_t i = 0; i < na && axis; i++) if (x->wts[i] > 0 && x->anchors[i].k == 0) axis = 0;
    if (axis) return apl_fit_axis_translation(x->anchors, x->wts, na, rot,
                       x->fit_cu, x->fit_cv, x->vm_a+x->vm_b*x->fit_s, x->vm_c, apl_axis_band(x), th, tx, ty);
    return apl_fit_translation(x->anchors, x->wts, na, rot, th, tx, ty);
}

/* Fit one cluster.  Returns the rms after robust trimming; kept anchors have
 * wts > 0.  res[] receives every anchor's residual (cluster members). */
static double apl_fit_cluster(AplCtx *x, size_t na, int cl, size_t min_members, double max_rot, double *th, double *tx, double *ty,
                              size_t *kept_out, double *res)
{
    size_t members = 0;
    double umin = 1e300, umax = -1e300, vmin = 1e300, vmax = -1e300;
    for (size_t i = 0; i < na; i++) {
        if (x->anchors[i].cluster == cl) {
            const AplAnchor *a = &x->anchors[i]; x->wts[i] = a->w; members++;
            umin = fmin(umin, a->sx); umax = fmax(umax, a->sx);
            vmin = fmin(vmin, a->sy); vmax = fmax(vmax, a->sy);
        }
        else x->wts[i] = 0.0;
        res[i] = 0.0;
    }
    if (members < min_members) { *kept_out = members; return 1e300; }
    /* Either coordinate constrains a rigid rotation. A long vertical run
     * must not get the short-patch rotation limit merely because its U
     * extent is small. The final axis/shape checks still validate the fit. */
    double extent = fmax(umax-umin, vmax-vmin);
    double rot = extent >= 300.0 ? max_rot : extent >= 200.0 ? fmin(max_rot, 0.25) : fmin(max_rot, 0.05);
    if (x->cand_grad < APL_SHEET_GRAD && extent >= 200.0) rot = APL_PI;   /* no usable z prior: the anchors decide */
    double rms = apl_fit_cluster_transform(x, na, rot, th, tx, ty);
    double c = cos(*th), s = sin(*th);
    Arena_Mark mark = Arena_save(x->run->arena);
    double *sorted = ARENA_ALLOC(x->run->arena, members * sizeof(double));
    size_t q = 0;
    for (size_t i = 0; i < na; i++) {
        if (x->wts[i] <= 0.0) continue;
        double px = c*x->anchors[i].sx - s*x->anchors[i].sy + *tx, py = s*x->anchors[i].sx + c*x->anchors[i].sy + *ty;
        res[i] = hypot(px - x->anchors[i].tx, py - x->anchors[i].ty);
        sorted[q++] = res[i];
    }
    qsort(sorted, q, sizeof(double), apl_cmp_double);
    double cut = fmax(3.0 * sorted[q/2], 2.0);
    size_t kept = 0;
    for (size_t i = 0; i < na; i++) if (x->wts[i] > 0.0) { if (res[i] <= cut) kept++; else x->wts[i] = 0.0; }
    if (kept >= min_members) {
        rms = apl_fit_cluster_transform(x, na, rot, th, tx, ty);
        c = cos(*th); s = sin(*th);
        for (size_t i = 0; i < na; i++) {
            if (x->anchors[i].cluster != cl) continue;
            double px = c*x->anchors[i].sx - s*x->anchors[i].sy + *tx, py = s*x->anchors[i].sx + c*x->anchors[i].sy + *ty;
            res[i] = hypot(px - x->anchors[i].tx, py - x->anchors[i].ty);
        }
    } else {
        kept = members;
        for (size_t i = 0; i < na; i++) x->wts[i] = x->anchors[i].cluster == cl ? x->anchors[i].w : 0.0;
    }
    Arena_restore(x->run->arena, mark);
    *kept_out = kept;
    return rms;
}

/* Rim isometry audit.  The layout is isometric, so two material points can
 * never be FARTHER apart in 3-D than in (u,v) (the 3-D distance is bounded by
 * the geodesic).  A chart placed into a hole of the wrong turn sits a layer
 * distance away in 3-D from the hole's rim while touching it in (u,v): the
 * excess d3 - duv over the rim is the signature that the co-location audit
 * cannot see because nothing overlaps.  Returns the fraction of the new
 * component's sampled vertices with placed material within `reach` in (u,v)
 * that violate the bound by more than max(0.5 d_layer, 7) vox; *n_rim gets
 * the number of such rim vertices. */
typedef struct AplRimSample { int64_t key; double u, v; float p[3]; int32_t chart; } AplRimSample;

static int apl_cmp_rim(const void *x, const void *y)
{
    const AplRimSample *a = x, *b = y;
    return (a->key > b->key) - (a->key < b->key);
}

static int64_t apl_rim_key(double u, double v, double cell)
{
    int64_t cu = (int64_t)floor(u / cell) + (1 << 30), cv = (int64_t)floor(v / cell) + (1 << 30);
    return (cu << 32) | (cv & 0xffffffffLL);
}

static double apl_rim_audit(const AsmRun *run, const int32_t *orig_comp, int32_t best, int32_t primary, double reach, size_t *n_rim)
{
    Arena_Mark mark = Arena_save(run->arena);
    size_t nc = run->n_charts, ns = 0, cap = 0;
    const size_t stride_placed = 4, stride_new = 2;
    for (size_t i = 0; i < nc; i++) {
        const AsmChart *c = &run->charts[i];
        if (!AsmChart_in_layout(c) || c->component != primary || orig_comp[i] == best) continue;
        cap += c->nv / stride_placed + 1;
    }
    *n_rim = 0;
    if (cap == 0) { Arena_restore(run->arena, mark); return 0.0; }
    AplRimSample *S = ARENA_ALLOC(run->arena, cap * sizeof(AplRimSample));
    double cell = reach;
    for (size_t i = 0; i < nc; i++) {
        const AsmChart *c = &run->charts[i];
        if (!AsmChart_in_layout(c) || c->component != primary || orig_comp[i] == best) continue;
        for (size_t k = 0; k < c->nv; k += stride_placed) {
            AplRimSample *s = &S[ns++];
            s->chart = c->id;
            apl_vertex_point(c, (int32_t)k, &s->u, &s->v);
            s->p[0] = c->xyz[k*3]; s->p[1] = c->xyz[k*3+1]; s->p[2] = c->xyz[k*3+2];
            s->key = apl_rim_key(s->u, s->v, cell);
        }
    }
    qsort(S, ns, sizeof(AplRimSample), apl_cmp_rim);
    size_t rim = 0, bad = 0;
    for (size_t i = 0; i < nc; i++) {
        const AsmChart *c = &run->charts[i];
        if (!AsmChart_in_layout(c) || orig_comp[i] != best) continue;
        double dl = run->chart_layer_d ? run->chart_layer_d[i] : 0.0;
        if (dl <= 0.0) dl = run->layer_d_global;
        double tol = fmax(0.5 * dl, 7.0);
        for (size_t k = 0; k < c->nv; k += stride_new) {
            double u, v;
            apl_vertex_point(c, (int32_t)k, &u, &v);
            const float *p = &c->xyz[k*3];
            int touched = 0, violated = 0;
            for (int du = -1; du <= 1 && !violated; du++) for (int dv = -1; dv <= 1 && !violated; dv++) {
                int64_t key = apl_rim_key(u + du * cell, v + dv * cell, cell);
                /* lower bound */
                size_t lo = 0, hi = ns;
                while (lo < hi) { size_t mid = (lo + hi) / 2; if (S[mid].key < key) lo = mid + 1; else hi = mid; }
                for (size_t q = lo; q < ns && S[q].key == key; q++) {
                    double duv = hypot(S[q].u - u, S[q].v - v);
                    if (duv > reach) continue;
                    touched = 1;
                    double d3 = sqrt((S[q].p[0]-p[0])*(S[q].p[0]-p[0]) + (S[q].p[1]-p[1])*(S[q].p[1]-p[1]) + (S[q].p[2]-p[2])*(S[q].p[2]-p[2]));
                    if (d3 - duv > tol) { violated = 1; break; }
                }
            }
            if (touched) { rim++; bad += (size_t)violated; }
        }
    }
    Arena_restore(run->arena, mark);
    *n_rim = rim;
    return rim ? (double)bad / (double)rim : 0.0;
}

/* Persistent RIM GRID: samples of the placed material hashed by (u,v) cell so the
 * rim isometry audit costs the candidate's size, not the layout's. */
struct AplRimGrid {
    Arena_T arena; double cell;
    int64_t *keys; int32_t *head; size_t cap, count;
    AplRimSample *s; int32_t *next; size_t n, ncap;
};

static size_t aprg_hash(int64_t key, size_t cap) { return (size_t)(((uint64_t)key * 0x9E3779B97F4A7C15ULL) >> 20) & (cap - 1); }

static void aprg_rehash(AplRimGrid *g, size_t ncap)
{
    int64_t *nk = ARENA_ALLOC(g->arena, ncap * sizeof(int64_t));
    int32_t *nh = ARENA_ALLOC(g->arena, ncap * sizeof(int32_t));
    for (size_t i = 0; i < ncap; i++) nh[i] = -1;
    for (size_t i = 0; i < g->cap; i++) {
        if (g->head[i] < 0) continue;
        size_t h = aprg_hash(g->keys[i], ncap);
        while (nh[h] >= 0) h = (h + 1) & (ncap - 1);
        nk[h] = g->keys[i]; nh[h] = g->head[i];
    }
    g->keys = nk; g->head = nh; g->cap = ncap;
}

static int32_t *aprg_slot(AplRimGrid *g, int64_t key, int create)
{
    if (create && (g->count + 1) * 2 > g->cap) aprg_rehash(g, g->cap * 2);
    size_t h = aprg_hash(key, g->cap), steps = 0;
    for (;;) {
        if (g->head[h] < 0) { if (!create) return NULL; g->keys[h] = key; g->count++; return &g->head[h]; }
        if (g->keys[h] == key) return &g->head[h];
        h = (h + 1) & (g->cap - 1);
        if (++steps > g->cap) { fprintf(stderr, "[rim grid] slot search exhausted: cap %zu count %zu key %lld\n", g->cap, g->count, (long long)key); return NULL; }
    }
}

static AplRimGrid *aprg_new_sized(Arena_T arena, double cell, size_t reserve)
{
    AplRimGrid *g = ARENA_CALLOC(arena, 1, sizeof *g);
    g->arena = arena; g->cell = cell; g->cap = 16;
    while (g->cap < (1u << 16) && g->cap/2 < reserve) g->cap *= 2;
    g->keys = ARENA_ALLOC(arena, g->cap * sizeof(int64_t));
    g->head = ARENA_ALLOC(arena, g->cap * sizeof(int32_t));
    for (size_t i = 0; i < g->cap; i++) g->head[i] = -1;
    g->ncap = 16;
    while (g->ncap < (1u << 16) && g->ncap < reserve) g->ncap *= 2;
    g->s = ARENA_ALLOC(arena, g->ncap * sizeof(AplRimSample));
    g->next = ARENA_ALLOC(arena, g->ncap * sizeof(int32_t));
    return g;
}

static AplRimGrid *aprg_new(Arena_T arena, double cell)
{
    return aprg_new_sized(arena,cell,1u << 16);
}

static void aprg_insert_chart_stride(AplRimGrid *g, const AsmChart *c, size_t stride)
{
    size_t need = c->nv / stride + 1;
    if (g->n + need > g->ncap) {
        size_t ncap = g->ncap; while (ncap < g->n + need) ncap *= 2;
        AplRimSample *ns = ARENA_ALLOC(g->arena, ncap * sizeof(AplRimSample));
        int32_t *nn = ARENA_ALLOC(g->arena, ncap * sizeof(int32_t));
        memcpy(ns, g->s, g->n * sizeof(AplRimSample)); memcpy(nn, g->next, g->n * sizeof(int32_t));
        g->s = ns; g->next = nn; g->ncap = ncap;
    }
    for (size_t k = 0; k < c->nv; k += stride) {
        AplRimSample *s = &g->s[g->n];
        s->chart = c->id;
        apl_vertex_point(c, (int32_t)k, &s->u, &s->v);
        s->p[0] = c->xyz[k*3]; s->p[1] = c->xyz[k*3+1]; s->p[2] = c->xyz[k*3+2];
        s->key = apl_rim_key(s->u, s->v, g->cell);
        int32_t *slot = aprg_slot(g, s->key, 1);
        if (slot == NULL) break;
        g->next[g->n] = *slot; *slot = (int32_t)g->n; g->n++;
    }
}

/* -1 outside the rim, 0 compatible, 1 violating. Most neighbours on a
 * sheet have nearly equal physical and UV distance. D3^2 <= Duv^2+tol^2
 * proves compatibility without evaluating two square roots. Conservative
 * margins leave both numerical boundaries to the original distance test. */
static int aprg_pair_test(double dx, double dy, const float *p, const float *q,
                          double reach, double reach2_safe, double reach2_inner, double tol)
{
    double uv2 = dx*dx+dy*dy;
    if (APL_RIM_SQUARED_CULL && uv2 > reach2_safe) return -1;
    double duv = -1;
    if (!(uv2 < reach2_inner)) {
        duv = hypot(dx,dy);
        if (duv > reach) return -1;
    }
    /* Keep the original float arithmetic for physical distance. Changing
     * it to double here would change decisions at the existing threshold. */
    double d3_squared = (q[0]-p[0])*(q[0]-p[0]) + (q[1]-p[1])*(q[1]-p[1]) + (q[2]-p[2])*(q[2]-p[2]);
    if (d3_squared <= (uv2+tol*tol)*(1.0-64.0*DBL_EPSILON)) return 0;
    if (duv < 0) duv = hypot(dx,dy);
    return sqrt(d3_squared)-duv > tol;
}

static void aprg_insert_chart(AplRimGrid *g, const AsmChart *c)
{
    aprg_insert_chart_stride(g,c,4);
}

/* A guarded UV trial inserts each chart once, in source-vertex order.
 * Map a grid sample back to the exact original reference vertex without
 * enlarging the persistent placement grid's sample representation. */
typedef struct AprgCheckpoint {
    const AsmChart **by_id;
    size_t *first_sample, count, reference_stride;
} AprgCheckpoint;

/* The candidate's charts against the grid: actual violating fraction of
 * rim vertices. With a checkpoint, the optional flags/partners describe
 * only newly violating vertex PAIRS, including at already bad queries. */
static double aprg_probe_record(AplRimGrid *g, const AsmRun *run, const AsmChart *const *charts,
                                size_t n, double reach, size_t *n_rim, const uint8_t *excluded,
                                uint8_t *violations, int32_t *partners, int32_t *partner_samples,
                                const AprgCheckpoint *checkpoint, size_t query_stride)
{
    /* The grid and chart geometry stay fixed throughout a probe. Partition
     * query vertices, including large single charts, and reduce only
     * integer counts. The two supported strides (one and two) both divide
     * the block size. Temporary work descriptors
     * belong to this probe, not the persistent grid's arena. */
    typedef struct AprgRange { const AsmChart *chart; size_t begin, query_base; } AprgRange;
    const size_t block = 512;
    Arena_Mark mark = Arena_save(run->arena);
    size_t nranges = 0;
    for (size_t k = 0; k < n; k++) nranges += (charts[k]->nv+block-1)/block;
    AprgRange *ranges = ARENA_ALLOC(run->arena,(nranges ? nranges : 1)*sizeof *ranges);
    size_t at = 0, query_base = 0;
    for (size_t k = 0; k < n; k++) {
        for (size_t v = 0; v < charts[k]->nv; v += block)
            ranges[at++] = (AprgRange){charts[k],v,query_base};
        query_base += (charts[k]->nv+query_stride-1)/query_stride;
    }
    size_t rim = 0, bad = 0;
    double reach2_safe = reach * reach * (1.0 + 8.0 * DBL_EPSILON);
    double reach2_inner = reach * reach * (1.0 - 64.0 * DBL_EPSILON);
    ptrdiff_t range;
#ifdef _OPENMP
    int threads = omp_get_max_threads();
    if ((size_t)threads > nranges) threads = nranges ? (int)nranges : 1;
#pragma omp parallel for schedule(dynamic,1) reduction(+:rim,bad) num_threads(threads) if(nranges >= 8)
#endif
    for (range = 0; range < (ptrdiff_t)nranges; range++) {
        const AsmChart *c = ranges[range].chart;
        const AsmChart *old = checkpoint && c->id >= 0 && (size_t)c->id < checkpoint->count ? checkpoint->by_id[c->id] : NULL;
        double dl = run->chart_layer_d ? run->chart_layer_d[c->id] : 0.0;
        if (dl <= 0.0) dl = run->layer_d_global;
        double tol = fmax(0.5 * dl, 7.0);
        size_t end = ranges[range].begin+block; if (end > c->nv) end = c->nv;
        for (size_t v = ranges[range].begin; v < end; v += query_stride) {
            double u, w;
            apl_vertex_point(c, (int32_t)v, &u, &w);
            const float *p = &c->xyz[v*3];
            double old_u = 0, old_w = 0;
            int old_valid = old && v < old->nv && old->uv && old->xyz;
            if (old_valid) apl_vertex_point(old,(int32_t)v,&old_u,&old_w);
            int touched = 0, violated = 0, actual_bad = 0; int32_t partner = -1, partner_sample = -1;
            for (int du = -1; du <= 1 && !violated; du++) for (int dv = -1; dv <= 1 && !violated; dv++) {
                int32_t *slot = aprg_slot(g, apl_rim_key(u + du * g->cell, w + dv * g->cell, g->cell), 0);
                if (slot == NULL) continue;
                size_t walk = 0;
                for (int32_t q = *slot; q >= 0; q = g->next[q]) {
                    const AplRimSample *s = &g->s[q];
                    if (++walk > g->n + 1) { fprintf(stderr, "[rim grid] list cycle at node %d (n %zu)\n", q, g->n); break; }
                    if (excluded && excluded[s->chart]) continue;
                    double dx = s->u - u, dy = s->v - w;
                    int result = aprg_pair_test(dx,dy,p,s->p,reach,reach2_safe,reach2_inner,tol);
                    if (result < 0) continue;
                    touched = 1;
                    if (result) {
                        actual_bad = 1;
                        if (checkpoint && old_valid && s->chart >= 0 && (size_t)s->chart < checkpoint->count) {
                            const AsmChart *ref = checkpoint->by_id[s->chart];
                            size_t first = checkpoint->first_sample[s->chart];
                            if (ref && ref->nv && ref->uv && ref->xyz && (size_t)q >= first &&
                                (size_t)q-first <= (ref->nv-1)/checkpoint->reference_stride) {
                                size_t vertex = checkpoint->reference_stride*((size_t)q-first); double ru, rv;
                                apl_vertex_point(ref,(int32_t)vertex,&ru,&rv);
                                if (aprg_pair_test(ru-old_u,rv-old_w,old->xyz+3*v,ref->xyz+3*vertex,
                                                   reach,reach2_safe,reach2_inner,tol) > 0) continue;
                            }
                        }
                        /* Missing checkpoint identity cannot certify a new
                         * conflict as pre-existing. Reject it conservatively. */
                        violated = 1; partner = s->chart; partner_sample = q; break;
                    }
                }
            }
            if (violations) violations[ranges[range].query_base+v/query_stride] = (uint8_t)violated;
            if (partners) partners[ranges[range].query_base+v/query_stride] = partner;
            if (partner_samples) partner_samples[ranges[range].query_base+v/query_stride] = partner_sample;
            if (touched) { rim++; bad += (size_t)actual_bad; }
        }
    }
    Arena_restore(run->arena,mark);
    *n_rim = rim;
    return rim ? (double)bad / (double)rim : 0.0;
}

static double aprg_probe_excluding(AplRimGrid *g, const AsmRun *run, const AsmChart *const *charts,
                                    size_t n, double reach, size_t *n_rim, const uint8_t *excluded)
{
    return aprg_probe_record(g,run,charts,n,reach,n_rim,excluded,NULL,NULL,NULL,NULL,2);
}

static double aprg_probe(AplRimGrid *g, const AsmRun *run, const AsmChart *const *charts, size_t n, double reach, size_t *n_rim)
{
    return aprg_probe_excluding(g,run,charts,n,reach,n_rim,NULL);
}

/* Bounded examples from an already rejected SOURCE proposal. Record actual
 * trial UV and physical points, so a failure can be checked independently
 * without reconstructing modified UV from the original chart cache. These
 * are examples, not a census; they never participate in placement decisions. */
static void apl_write_rim_failure(AplCtx *x, int32_t component, const AsmChart *const *charts, size_t n)
{
    unsigned bit = x->uv_trial ? 2u : 1u;
    if (!x->frim || !x->diag || x->cluster_kind != 1 || (x->diag[component].rim_logged & bit)) return;
    x->diag[component].rim_logged |= bit;
    AplRimGrid *g = x->rgrid;
    size_t rows = 0, queries = 0, examined = 0;
    const double reach = 24, outer = 24*24*(1.0+8.0*DBL_EPSILON), inner = 24*24*(1.0-64.0*DBL_EPSILON);
    for (size_t k = 0; k < n && rows < 32 && queries < 2048 && examined < 32768; k++) {
        const AsmChart *c = charts[k];
        double dl = x->run->chart_layer_d ? x->run->chart_layer_d[c->id] : 0;
        if (dl <= 0) dl = x->run->layer_d_global;
        double tol = fmax(.5*dl,7.0);
        for (size_t v = 0; v < c->nv && rows < 32 && queries < 2048 && examined < 32768; v += 2) {
            queries++;
            double u, w; apl_vertex_point(c,(int32_t)v,&u,&w);
            const float *p = c->xyz+3*v;
            int recorded = 0;
            for (int du = -1; du <= 1 && !recorded && examined < 32768; du++)
            for (int dv = -1; dv <= 1 && !recorded && examined < 32768; dv++) {
                int32_t *slot = aprg_slot(g,apl_rim_key(u+du*g->cell,w+dv*g->cell,g->cell),0);
                if (!slot) continue;
                for (int32_t q = *slot; q >= 0 && examined < 32768; q = g->next[q]) {
                    examined++;
                    const AplRimSample *s = &g->s[q];
                    if (x->excluded_samples && x->excluded_samples[s->chart]) continue;
                    double dx = s->u-u, dy = s->v-w;
                    if (aprg_pair_test(dx,dy,p,s->p,reach,outer,inner,tol) != 1) continue;
                    double physical = sqrt((s->p[0]-p[0])*(s->p[0]-p[0])+(s->p[1]-p[1])*(s->p[1]-p[1])+(s->p[2]-p[2])*(s->p[2]-p[2]));
                    const AsmChart *partner = &x->run->charts[s->chart];
                    fprintf(x->frim,"%d,%d,%d,%d,%d,%zu,%d,%d,%d,%.17g,%.17g,%.17g,%.17g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.17g,%.17g,%.17g\n",
                            x->attempt,component,x->cluster,x->uv_trial,c->id,v,s->chart,q,partner->placement_state,
                            u,w,s->u,s->v,p[0],p[1],p[2],s->p[0],s->p[1],s->p[2],hypot(dx,dy),physical,tol);
                    rows++; recorded = 1; break;
                }
            }
        }
    }
}

/* A UV polish keeps every chart pose fixed and every vertex within eight
 * voxels of saved[].uv. If a rim pair still violates after increasing its
 * UV distance by that entire movement allowance, no permitted polish can
 * fix this query. Count against ALL candidate queries, an upper bound on
 * the final touched denominator. This can only prove failure; inconclusive
 * candidates still take the original solve and full audits. */
static int aprg_polish_blocked(AplRimGrid *g, const AsmRun *run, const AsmChart *saved, size_t n,
                               const uint8_t *excluded, double max_ratio, size_t *guaranteed, size_t *queries)
{
    *guaranteed = *queries = 0;
    if (!g || !g->n || !isfinite(max_ratio) || max_ratio < 0 || max_ratio >= 1) return 0;
    for (size_t k = 0; k < n; k++) {
        if (saved[k].id < 0 || (size_t)saved[k].id >= run->n_charts ||
            run->charts[saved[k].id].nv != saved[k].nv || !saved[k].uv ||
            !run->charts[saved[k].id].uv || !run->charts[saved[k].id].xyz) return 0;
        *queries += (saved[k].nv+1)/2;
    }
    for (size_t k = 0; k < n; k++) {
        const AsmChart *c = &run->charts[saved[k].id];
        double dl = run->chart_layer_d ? run->chart_layer_d[c->id] : 0;
        if (dl <= 0) dl = run->layer_d_global;
        double tol = fmax(.5*dl,7.0);
        for (size_t v = 0; v < c->nv; v += 2) {
            /* Earlier joint repair may already have consumed part of the
             * original disk. A later polish can move back across that disk,
             * so eight voxels from the CURRENT point is not a valid bound. */
            double consumed = hypot((double)c->uv[2*v]-saved[k].uv[2*v],
                                    (double)c->uv[2*v+1]-saved[k].uv[2*v+1]);
            double allowance = 8.0+consumed+1e-4, reach = 24.0-allowance;
            if (!(reach > 0)) continue;
            double u, w; apl_vertex_point(c,(int32_t)v,&u,&w);
            const float *p = c->xyz+3*v;
            int blocked = 0;
            for (int du = -1; du <= 1 && !blocked; du++) for (int dv = -1; dv <= 1 && !blocked; dv++) {
                int32_t *slot = aprg_slot(g,apl_rim_key(u+du*g->cell,w+dv*g->cell,g->cell),0);
                if (!slot) continue;
                for (int32_t q = *slot; q >= 0; q = g->next[q]) {
                    const AplRimSample *s = &g->s[q];
                    if (excluded && excluded[s->chart]) continue;
                    if (aprg_pair_test(s->u-u,s->v-w,p,s->p,reach,
                                      reach*reach*(1+8*DBL_EPSILON),reach*reach*(1-64*DBL_EPSILON),tol+allowance) > 0) {
                        blocked = 1; break;
                    }
                }
            }
            *guaranteed += (size_t)blocked;
            if (*guaranteed >= 8 && (double)*guaranteed > max_ratio*(double)*queries) return 1;
        }
    }
    return 0;
}

/* Joint repair changes relative chart poses. The persistent grid contains
 * only previously placed material, so it cannot see a repair folding the
 * candidate onto itself. Compare that internal geometry before and after
 * repair with the same spatial isometry test. A common rigid transform
 * leaves this count unchanged. */
static size_t apl_internal_record(AsmRun *run, const AsmChart *saved, size_t n, int current,
                                   uint8_t *violations, int32_t *partners, int32_t *partner_vertices, int preserve_pairs)
{
    /* A single chart can contract or fold during UV repair. Relative rigid
     * registration needs two charts, but this audit also certifies UV edits. */
    if (!n) return 0;
    Arena_Mark mark = Arena_save(run->arena);
    /* A local UV edit can change an odd-indexed vertex independently.
     * Guarded post-placement repair must therefore audit every vertex on
     * both sides; a sparse query/reference set cannot certify that edit. */
    const size_t query_stride = preserve_pairs ? 1 : 2, reference_stride = preserve_pairs ? 1 : 4;
    /* Most trial components are small. Avoid initializing a 65,536-slot
     * persistent-layout table for each of their before/after audits.
     * Sampling, cell size, neighbour search, and growth remain identical. */
    size_t reserve = 0;
    for (size_t k = 0; k < n; k++) {
        size_t need = saved[k].nv/reference_stride+1;
        if (need >= (1u << 16)-reserve) { reserve = 1u << 16; break; }
        reserve += need;
    }
    AplRimGrid *grid = aprg_new_sized(run->arena,24.0,reserve);
    const AsmChart **charts = ARENA_ALLOC(run->arena, n*sizeof(*charts));
    AprgCheckpoint checkpoint = {0};
    if (current && (preserve_pairs || partner_vertices)) {
        checkpoint.count = run->n_charts;
        checkpoint.reference_stride = reference_stride;
        checkpoint.by_id = ARENA_CALLOC(run->arena,run->n_charts,sizeof *checkpoint.by_id);
        checkpoint.first_sample = ARENA_CALLOC(run->arena,run->n_charts,sizeof *checkpoint.first_sample);
    }
    for (size_t k = 0; k < n; k++) {
        charts[k] = current ? &run->charts[saved[k].id] : &saved[k];
        if (checkpoint.by_id) {
            checkpoint.by_id[saved[k].id] = &saved[k];
            checkpoint.first_sample[saved[k].id] = grid->n;
        }
        aprg_insert_chart_stride(grid,charts[k],reference_stride);
    }
    size_t touched = 0;
    double bad = aprg_probe_record(grid,run,charts,n,24.0,&touched,NULL,violations,partners,partner_vertices,
                                  current && preserve_pairs ? &checkpoint : NULL,query_stride);
    if (partner_vertices) for (size_t k = 0, at = 0; k < n; k++) for (size_t v = 0; v < saved[k].nv; v += query_stride, at++) {
        int32_t sample = partner_vertices[at]; partner_vertices[at] = -1;
        if (sample < 0 || (size_t)sample >= grid->n || !checkpoint.by_id) continue;
        int32_t id = grid->s[sample].chart;
        if (id < 0 || (size_t)id >= checkpoint.count || !checkpoint.by_id[id]) continue;
        size_t first = checkpoint.first_sample[id], count = checkpoint.by_id[id]->nv;
        if (count && (size_t)sample >= first && (size_t)sample-first <= (count-1)/reference_stride)
            partner_vertices[at] = (int32_t)(reference_stride*((size_t)sample-first));
    }
    size_t count = (size_t)floor(bad*(double)touched+0.5);
    Arena_restore(run->arena, mark);
    return count;
}

static size_t apl_internal_bad(AsmRun *run, const AsmChart *saved, size_t n, int current)
{
    return apl_internal_record(run,saved,n,current,NULL,NULL,NULL,0);
}

/* Preserve the identity of every previously sound query vertex. Counts
 * alone let a repair trade one overlap for another while claiming that no
 * new internal contradiction was introduced. The cache belongs to the
 * caller's original checkpoint and stores that audit's query vertices.
 * A placed UV neighbourhood also preserves reference-vertex identity:
 * existing overlaps must not license additional conflicting pairs. */
static size_t apl_internal_new(AsmRun *run, const AsmChart *saved, size_t n,
                               AsmContinuityAuditCache *cache, size_t *after_bad)
{
    size_t queries = 0; *after_bad = 0;
    cache->pairs_count = 0;
    int preserve_pairs = cache->guard_charts != NULL;
    const size_t query_stride = preserve_pairs ? 1 : 2;
    if (cache->blocked_charts) {
        if (cache->blocked_charts_count != run->n_charts) return SIZE_MAX;
        memset(cache->blocked_charts,0,run->n_charts);
    }
    if (cache->blocked_vertices) {
        if (cache->blocked_vertices_count != run->n_charts) return SIZE_MAX;
        for (size_t i = 0; i < run->n_charts; i++) if (cache->blocked_vertices[i])
            memset(cache->blocked_vertices[i],0,run->charts[i].nv);
    }
    for (size_t i = 0; i < n; i++) {
        if (saved[i].id < 0 || (size_t)saved[i].id >= run->n_charts || run->charts[saved[i].id].nv != saved[i].nv) return SIZE_MAX;
        queries += (saved[i].nv+query_stride-1)/query_stride;
    }
    if (!cache->baseline_ready) {
        Arena_T scratch = cache->scratch ? cache->scratch : run->arena;
        cache->baseline_violations = ARENA_CALLOC(scratch,queries ? queries : 1,1);
        cache->baseline_queries = queries;
        cache->baseline_bad = apl_internal_record(run,saved,n,0,cache->baseline_violations,NULL,NULL,preserve_pairs);
        cache->baseline_ready = 1;
    }
    if (cache->baseline_queries != queries || !cache->baseline_violations) return SIZE_MAX;
    Arena_Mark mark = Arena_save(run->arena);
    uint8_t *current = ARENA_CALLOC(run->arena,queries ? queries : 1,1);
    int report_pairs = cache->pairs && cache->pairs_capacity;
    int32_t *partners = cache->blocked_charts || cache->blocked_vertices || report_pairs ? ARENA_ALLOC(run->arena,(queries ? queries : 1)*sizeof *partners) : NULL;
    int32_t *vertices = cache->blocked_vertices || report_pairs ? ARENA_ALLOC(run->arena,(queries ? queries : 1)*sizeof *vertices) : NULL;
    *after_bad = apl_internal_record(run,saved,n,1,current,partners,vertices,preserve_pairs);
    size_t introduced = 0;
    for (size_t k = 0, q = 0; k < n; k++) for (size_t v = 0; v < saved[k].nv; v += query_stride, q++) {
        if (!current[q] || (!preserve_pairs && cache->baseline_violations[q])) continue;
        introduced++;
        if (cache->blocked_charts) {
            cache->blocked_charts[saved[k].id] = 1;
            if (partners[q] >= 0 && (size_t)partners[q] < run->n_charts) cache->blocked_charts[partners[q]] = 1;
        }
        if (cache->blocked_vertices) {
            if (cache->blocked_vertices[saved[k].id]) cache->blocked_vertices[saved[k].id][v] = 1;
            int32_t id = partners[q], vertex = vertices[q];
            if (id >= 0 && (size_t)id < run->n_charts && vertex >= 0 && (size_t)vertex < run->charts[id].nv && cache->blocked_vertices[id])
                cache->blocked_vertices[id][vertex] = 1;
        }
        if (report_pairs && cache->pairs_count < cache->pairs_capacity) {
            int32_t id = partners[q], vertex = vertices[q];
            if (id >= 0 && (size_t)id < run->n_charts && vertex >= 0 && (size_t)vertex < run->charts[id].nv)
                cache->pairs[cache->pairs_count++] = (AsmContinuityUvPair){saved[k].id,(int32_t)v,id,vertex};
        }
    }
    Arena_restore(run->arena,mark); return introduced;
}

typedef struct AplSelfContact {
    const AsmChart *saved;
    AsmContinuityAuditCache *cache;
    const AsmContacts *query;
    const double *before;
    size_t introduced;
} AplSelfContact;

static void apl_self_contact(void *context, size_t a, size_t b, double area)
{
    AplSelfContact *x = context;
    const AsmChart *c = x->saved;
    double before=AsmContacts_pair_faces(x->query,x->before,a,b,NULL,NULL);
    if (isfinite(before) && before > 0 && area <= before+fmax(1e-10,1e-10*before)) return;
    x->introduced++;
    if (x->cache->blocked_charts) x->cache->blocked_charts[c->id] = 1;
    if (x->cache->blocked_vertices && x->cache->blocked_vertices[c->id]) {
        for (int k = 0; k < 3; k++) {
            x->cache->blocked_vertices[c->id][c->faces[3*a+k]] = 1;
            x->cache->blocked_vertices[c->id][c->faces[3*b+k]] = 1;
        }
    }
}

/* A sub-voxel self-overlap can preserve every local face metric and every
 * vertex-distance inequality. Check the actual triangles of every changed
 * chart as well. Compare each contact with the same original face pair:
 * reducing one existing contact cannot license a new pair or enlarge one.
 * Local UVs avoid roundoff from an arbitrary rigid placement or reflection. */
static int apl_self_contacts_preserved(AsmRun *run, const AsmChart *saved, size_t n, AsmContinuityAuditCache *cache)
{
    for (size_t i = 0; i < n; i++) {
        const AsmChart *old = &saved[i], *c = &run->charts[old->id];
        if (c->nf < 2) continue;
        if (!old->uv || !c->uv || c->nv != old->nv || c->nf != old->nf || c->faces != old->faces ||
            c->nv > SIZE_MAX/(2*sizeof(double))) return 0;
        if (c->uv == old->uv || !memcmp(c->uv,old->uv,2*c->nv*sizeof(float))) continue;
        Arena_Mark mark = Arena_save(run->arena);
        double *uv = ARENA_ALLOC(run->arena,2*c->nv*sizeof *uv);
        double *before=ARENA_ALLOC(run->arena,2*c->nv*sizeof *before);
        for (size_t v = 0; v < 2*c->nv; v++) uv[v] = c->uv[v];
        for(size_t v=0;v<2*c->nv;v++)before[v]=old->uv[v];
        AsmContacts *query = AsmContacts_new(run->arena,c->nv,c->faces,c->nf,uv,AsmContacts_material_tolerance());
        AsmContactStats stats;
        AplSelfContact context={old,cache,query,before,0};
        int ok = query && !AsmContacts_measure(query,uv,NULL,apl_self_contact,&context,&stats) && !context.introduced;
        Arena_restore(run->arena,mark);
        if (!ok) return 0;
    }
    return 1;
}

static int apl_uv_audit(AsmRun *run, const AsmChart *saved, size_t n, AsmContinuityAuditCache *cache)
{
    size_t after;
    cache->pairs_count = 0;
    if (!AsmContinuity_preserves_incumbent_trial(run,saved,n)) return 0;
    if (!!cache->guard_charts != !!cache->guard_count) return 0;
    if (cache->guard_charts) { saved = cache->guard_charts; n = cache->guard_count; }
    return apl_internal_new(run,saved,n,cache,&after) == 0 && apl_self_contacts_preserved(run,saved,n,cache);
}

/* Registration changes relative poses just like a later joint repair.
 * Its fixed chart is part of the audit: a two-chart component must not
 * become overlapping merely because its moving subset has only one chart. */
static int apl_registration_audit(AsmRun *run, const AsmChart *saved, size_t n, AsmContinuityAuditCache *cache)
{
    if (!AsmContinuity_preserves_incumbent_trial(run,saved,n)) return 0;
    size_t after, introduced = apl_internal_new(run,saved,n,cache,&after);
    if (introduced) fprintf(stderr,"  [registration geometry] component %d, %zu charts: bad rim vertices %zu -> %zu, %zu %s; trial rejected\n",
                            saved[0].component,n,cache->baseline_bad,after,introduced,
                            cache->guard_charts ? "queries with new conflicting pairs" : "previously sound vertices violated");
    return introduced == 0;
}

/* Tentatively place `best` with (th, tx, ty); audit; revert on contradiction.
 * Returns 1 placed, 0 reverted. */
/* ---- placement-gated cut ----------------------------------------------------------- */

/* Two clusters that each fit but disagree: cut component `best` along the seam relations
 * between the charts the two clusters touch.  Returns the new component id of the smaller
 * side, or -1 when no cut applies.  Charts are labelled by a multi-source BFS over the active
 * relations from the two clusters' chart sets (nearest seed set wins, ties to the heavier
 * cluster). Conflicting layer votes do not authorize tearing an already
 * measured source continuation. Only a partition that crosses no sound
 * source seam or retained incumbent may commit; otherwise the component
 * remains one placement candidate. Failed seams remain recoverable. */
/* One row per CUT PROPOSAL, whichever precondition refused it.  A cut is the only way a
 * component whose anchors form two clusters a wrap apart can ever be placed, and on the
 * 21x5x5 all 370 proposals die; without this the summary line says only how many. */
typedef struct AplCutRow {
    int32_t component; int c1, c2;
    double w1, w2, rms1, rms2, u1, u2, sep, min_sep, t0;
    size_t n1, n2, m1, m2;
    double area1, area2;
    size_t crossing, crossing_sound, sound_max_support;
    double sound_min_rms;
    long first_sound_rel;
} AplCutRow;

static void apl_write_cut(const AplCtx *x, const AplCutRow *r, const char *outcome)
{
    if (!x->fcut) return;
    fprintf(x->fcut, "%d,%d,%s,%d,%d,%.6g,%.6g,%.6g,%.6g,%.2f,%.2f,%.2f,%.2f,%.2f,"
                     "%zu,%zu,%zu,%zu,%.6e,%.6e,%zu,%zu,%.6g,%zu,%ld\n",
            x->attempt, r->component, outcome, r->c1, r->c2, r->w1, r->w2, r->rms1, r->rms2,
            r->u1, r->u2, r->sep, r->min_sep, r->t0,
            r->n1, r->n2, r->m1, r->m2, r->area1, r->area2,
            r->crossing, r->crossing_sound, r->sound_min_rms, r->sound_max_support, r->first_sound_rel);
}

static int32_t apl_cut_component(AplCtx *x, int32_t best, size_t na, int ncl, const double *cl_w, const double *cl_u, const double *cl_rms,
                                 double gate, const int32_t *rel_off, const int32_t *rel_adj, int32_t *ncomp_io, size_t ncomp_cap,
                                 double *comp_area, int32_t *orig_comp, AplOrient *ori, AplCompDiag *diag, int32_t *comp_ring, const double *lidx,
                                 size_t *rels_deferred)
{
    AsmRun *run = x->run;
    size_t nc = run->n_charts;
    /* the proposal's ledger row: every refusal reports what it saw (APL_CUT_LEDGER) */
    AplCutRow row; memset(&row, 0, sizeof row);
    row.component = best; row.c1 = row.c2 = -1; row.first_sound_rel = -1;
    row.sound_min_rms = -1.0;
    if ((size_t)*ncomp_io + 1 >= ncomp_cap) { apl_write_cut(x, &row, "component_cap"); return -1; }
    /* the two heaviest clusters that fit, far enough apart in u */
    int c1 = -1, c2 = -1;
    for (int cl = 0; cl < ncl; cl++) {
        if (cl_rms[cl] > gate) continue;
        if (c1 < 0 || cl_w[cl] > cl_w[c1]) { c2 = c1; c1 = cl; }
        else if (c2 < 0 || cl_w[cl] > cl_w[c2]) c2 = cl;
    }
    row.c1 = c1; row.c2 = c2; row.t0 = x->t0;
    if (c1 >= 0) { row.w1 = cl_w[c1]; row.rms1 = cl_rms[c1]; row.u1 = cl_u[c1]; }
    if (c2 >= 0) { row.w2 = cl_w[c2]; row.rms2 = cl_rms[c2]; row.u2 = cl_u[c2]; }
    if (c1 < 0 || c2 < 0) { apl_write_cut(x, &row, "no_two_fitting_clusters"); return -1; }
    double sep = fabs(cl_u[c1] - cl_u[c2]), min_sep = fmax(APL_CUT_MIN_SEP, 0.25 * x->t0);
    row.sep = sep; row.min_sep = min_sep;
    if (sep < min_sep) { apl_write_cut(x, &row, "clusters_too_close"); return -1; }
    Arena_Mark mark = Arena_save(run->arena);
    int8_t *label = ARENA_CALLOC(run->arena, nc, 1);     /* 0 none, 1 side of c1, 2 side of c2, 3 contested */
    for (size_t i = 0; i < na; i++) {
        const AplAnchor *a = &x->anchors[i];
        if (a->cluster != c1 && a->cluster != c2) continue;
        int8_t l = a->cluster == c1 ? 1 : 2;
        int32_t ch = a->chart_c;
        if (ch < 0 || (size_t)ch >= nc || run->charts[(size_t)ch].component != best) continue;
        if (label[ch] == 0) label[ch] = l;
        else if (label[ch] != l) label[ch] = 3;
    }
    size_t n1 = 0, n2 = 0;
    for (size_t i = 0; i < nc; i++) { if (label[i] == 3) label[i] = 0; else if (label[i] == 1) n1++; else if (label[i] == 2) n2++; }
    row.n1 = n1; row.n2 = n2;
    if (n1 == 0 || n2 == 0) { apl_write_cut(x, &row, "no_distinct_charts"); Arena_restore(run->arena, mark); return -1; }
    /* multi-source BFS over active relations within the component: nearest seed set */
    int32_t *queue = ARENA_ALLOC(run->arena, nc * sizeof(int32_t));
    size_t qh = 0, qt = 0;
    for (size_t i = 0; i < nc; i++) if (label[i] != 0 && run->charts[i].component == best) queue[qt++] = (int32_t)i;
    while (qh < qt) {
        int32_t u = queue[qh++];
        for (int32_t e = rel_off[u]; e < rel_off[u + 1]; e++) {
            const AsmRelation *R = &run->rels[(size_t)rel_adj[e]];
            if (R->flags & (ASM_REL_DROPPED | ASM_REL_CONTACT | ASM_REL_SWITCHED | ASM_REL_PLACEMENT_ONLY)) continue;
            int32_t v = R->a == u ? R->b : R->a;
            if (!AsmChart_in_layout(&run->charts[(size_t)v]) || run->charts[(size_t)v].component != best) continue;
            if (label[v] != 0) continue;
            label[v] = label[u];
            queue[qt++] = v;
        }
    }
    size_t m1 = 0, m2 = 0; double a1 = 0.0, a2 = 0.0;
    for (size_t i = 0; i < nc; i++) {
        if (!AsmChart_in_layout(&run->charts[i]) || run->charts[i].component != best) continue;
        if (label[i] == 0) label[i] = 1;   /* unreachable (should not happen in a connected component): the heavier side */
        if (label[i] == 1) { m1++; a1 += run->charts[i].area3d; } else { m2++; a2 += run->charts[i].area3d; }
    }
    row.m1 = m1; row.m2 = m2; row.area1 = a1; row.area2 = a2;
    if (m1 == 0 || m2 == 0) { apl_write_cut(x, &row, "partition_empty"); Arena_restore(run->arena, mark); return -1; }
    /* the smaller side becomes the new component */
    int8_t moved = a2 <= a1 ? 2 : 1;
    /* The old cut committed before either resulting placement was known.
     * A fragment could then land on another wrap while its seam neighbour
     * remained unplaced, making the promised later seam closure impossible.
     * Preserve source continuations at the partition boundary itself. */
    /* Survey EVERY crossing relation before deciding.  The verdict is the same -- any sound
     * source continuation across the partition refuses the cut -- but stopping at the first
     * one told us nothing about the other 369 proposals the 21x5x5 makes. */
    size_t crossing = 0, crossing_sound = 0, sound_support = 0;
    double sound_min_rms = -1.0;
    long first_sound = -1;
    for (size_t r = 0; r < run->n_rels; r++) {
        const AsmRelation *R = &run->rels[r];
        if (run->charts[R->a].component != best || run->charts[R->b].component != best ||
            label[R->a] == label[R->b]) continue;
        crossing++;
        double boundary_rms = -1.0; size_t boundary_support = 0;
        if (!((R->continuity & ASM_CONT_INCUMBENT) ||
              AsmContinuity_measure(run,R,&boundary_rms,&boundary_support))) continue;
        crossing_sound++;
        if (first_sound < 0) first_sound = (long)r;
        if (boundary_rms >= 0.0 && (sound_min_rms < 0.0 || boundary_rms < sound_min_rms)) sound_min_rms = boundary_rms;
        if (boundary_support > sound_support) sound_support = boundary_support;
    }
    row.crossing = crossing; row.crossing_sound = crossing_sound;
    row.sound_min_rms = sound_min_rms; row.sound_max_support = sound_support; row.first_sound_rel = first_sound;
    if (crossing_sound) {
        x->cuts_preserved_source++;
        apl_write_cut(x, &row, "source_preserved");
        Arena_restore(run->arena,mark); return -1;
    }
    apl_write_cut(x, &row, "cut");
    int32_t newc = (*ncomp_io)++;
    size_t deferred = 0;
    for (size_t r = 0; r < run->n_rels; r++) {
        AsmRelation *R = &run->rels[r];
        if (R->flags & (ASM_REL_DROPPED | ASM_REL_CONTACT | ASM_REL_PLACEMENT_ONLY)) continue;
        if (run->charts[(size_t)R->a].component != best || run->charts[(size_t)R->b].component != best) continue;
        if (label[R->a] != label[R->b]) {
            double boundary_rms; size_t boundary_support;
            if ((R->continuity & ASM_CONT_VERIFIED) && AsmContinuity_measure(run, R, &boundary_rms, &boundary_support))
                R->continuity |= ASM_CONT_CUT;
            R->flags |= ASM_REL_SWITCHED;
            R->continuity &= ~(uint32_t)ASM_CONT_VERIFIED;
            deferred++;
        }
    }
    comp_area[newc] = 0.0; comp_area[best] = 0.0;
    comp_ring[newc] = 1000000; comp_ring[best] = 1000000;
    for (size_t i = 0; i < nc; i++) {
        AsmChart *c = &run->charts[i];
        if (!AsmChart_in_layout(c) || c->component != best) continue;
        if (label[i] == moved) { c->component = newc; orig_comp[i] = newc; }
        int32_t cc = c->component;
        comp_area[cc] += c->area3d;
        if (lidx[i] > APL_LAYER_UNKNOWN) { int32_t ring = (int32_t)floor(fabs(lidx[i]) + 0.5); if (ring < comp_ring[cc]) comp_ring[cc] = ring; }
    }
    ori[newc] = ori[best];
    memset(&diag[newc], 0, sizeof diag[newc]);
    *rels_deferred += deferred;
    fprintf(stderr, "  [place] CUT component %d: clusters %d (u %.0f, weight %.1f, rms %.1f) and %d (u %.0f, weight %.1f, rms %.1f) disagree by %.0f vox; "
                    "seeds %zu + %zu charts, parts %zu charts / %.3e vox^2 and %zu charts / %.3e vox^2, %zu seam relations deferred -> new component %d\n",
            best, c1 + 1, cl_u[c1], cl_w[c1], cl_rms[c1], c2 + 1, cl_u[c2], cl_w[c2], cl_rms[c2], sep, n1, n2,
            moved == 2 ? m1 : m2, moved == 2 ? a1 : a2, moved == 2 ? m2 : m1, moved == 2 ? a2 : a1, deferred, newc);
    Arena_restore(run->arena, mark);
    return newc;
}

/* ---- metric re-solve --------------------------------------------------------------- */

typedef struct AplObs { int32_t ca, cb; double du, dv, w, k, s, wt; } AplObs;   /* one placed pair: b is s*k wraps from a */

static size_t apl_metric_resolve(AplCtx *x, const uint8_t *comp_placed, const int32_t *orig_comp, int32_t ncomp,
                                 const double *lidx, int32_t primary, Arena_T *garena, const AsmConflictOpts *copts, size_t *n_obs_out,
                                 int model_only)
{
    AsmRun *run = x->run;
    Arena_T arena = run->arena;
    size_t nc = run->n_charts;
    *n_obs_out = 0;
    if (lidx == NULL) return 0;
    Arena_Mark mark = Arena_save(arena);
    /* 1. one observation per placed pair: the median hit displacement in the layout */
    AplObs *obs = ARENA_ALLOC(arena, (run->n_layers ? run->n_layers : 1) * sizeof(AplObs));
    size_t no = 0;
    double *tmpu = ARENA_ALLOC(arena, 64 * sizeof(double)), *tmpv = ARENA_ALLOC(arena, 64 * sizeof(double));
    for (size_t p = 0; p < run->n_layers; p++) {
        const AsmLayerPair *lp = &run->layers[p];
        const AsmChart *A = &run->charts[(size_t)lp->a], *B = &run->charts[(size_t)lp->b];
        if (!AsmChart_in_layout(A) || !AsmChart_in_layout(B)) continue;
        int32_t ca = orig_comp[lp->a], cb = orig_comp[lp->b];
        if (ca < 0 || cb < 0 || comp_placed[ca] != 1 || comp_placed[cb] != 1) continue;
        if (lidx[lp->a] <= APL_LAYER_UNKNOWN || lidx[lp->b] <= APL_LAYER_UNKNOWN) continue;
        double dl = lidx[lp->b] - lidx[lp->a];
        if (fabs(dl) < 0.5) continue;
        int k = (int)floor(fabs(dl) + 0.5);
        if (k > 4) continue;
        size_t nh = 0;
        for (int h = 0; h < lp->hit_count && nh < 64; h++) {
            const AsmLayerHit *hit = &run->layer_hits[(size_t)lp->hit_first + (size_t)h];
            double ax, ay, bx, by;
            apl_vertex_point(A, hit->va, &ax, &ay);
            apl_face_point(B, hit->fb, hit->l0, hit->l1, &bx, &by);
            tmpu[nh] = bx - ax; tmpv[nh] = by - ay; nh++;
        }
        if (nh < 2) continue;
        qsort(tmpu, nh, sizeof(double), apl_cmp_double);
        qsort(tmpv, nh, sizeof(double), apl_cmp_double);
        AplObs *o = &obs[no++];
        o->ca = ca; o->cb = cb; o->du = tmpu[nh / 2]; o->dv = tmpv[nh / 2];
        o->w = 0.5 * (lidx[lp->a] + lidx[lp->b]); o->k = (double)k; o->s = dl < 0.0 ? -1.0 : 1.0;
        o->wt = (double)(lp->hit_count < 16 ? lp->hit_count : 16);
    }
    *n_obs_out = no;
    if (no < 20) { Arena_restore(arena, mark); return 0; }
    /* 2. the turn model on the observations whose sign agrees with the coordinate: T = |du| / k ~ t0 + g w */
    double t0 = x->t0, g = x->growth;
    {
        double *wr = ARENA_ALLOC(arena, no * sizeof(double));
        for (size_t i = 0; i < no; i++) wr[i] = 1.0;
        for (int it = 0; it < 6; it++) {
            double sw = 0, sx = 0, sy = 0, sxx = 0, sxy = 0;
            for (size_t i = 0; i < no; i++) {
                const AplObs *o = &obs[i];
                if (o->du * o->s <= 0.0) continue;
                double tt = fabs(o->du) / o->k, ww = wr[i] * o->wt;
                sw += ww; sx += ww * o->w; sy += ww * tt; sxx += ww * o->w * o->w; sxy += ww * o->w * tt;
            }
            if (sw <= 0.0) break;
            double det = sw * sxx - sx * sx;
            double gg = fabs(det) > 1e-9 ? (sw * sxy - sx * sy) / det : 0.0;
            if (gg < 0.0) gg = 0.0;
            if (gg > 400.0) gg = 400.0;
            double tt0 = (sy - gg * sx) / sw;
            if (tt0 < APL_TURN_MIN) tt0 = APL_TURN_MIN;
            if (tt0 > APL_TURN_MAX) tt0 = APL_TURN_MAX;
            t0 = tt0; g = gg;
            double scale = 0.15 * t0;
            for (size_t i = 0; i < no; i++) {
                const AplObs *o = &obs[i];
                double r = (o->du * o->s > 0.0) ? fabs(o->du) / o->k - (t0 + g * o->w) : 4.0 * scale;
                wr[i] = 1.0 / (1.0 + (r / scale) * (r / scale));
            }
        }
    }
    if (model_only) {
        /* THE TURN MODEL, POOLED, WITHOUT THE COMPONENT SHIFTS.  Part 3 below refits every
         * placed component's u offset and is MEASURED DEAD (10x10x10 placed 68.0 -> 59.2%:
         * the inconsistency sits inside the placed components' own geometry, so rigid shifts
         * move sound material into failing audits).  The model fitted above is a different
         * object and carries none of that risk: ONE t0 and ONE growth over every placed pair,
         * where the per-candidate stretch reads 20-79 hits of a single component and returns
         * turns from 645 to 9,105 vox.  Adopt the model alone. */
        double *res = ARENA_ALLOC(arena, no * sizeof(double));
        size_t nr = 0;
        for (size_t i = 0; i < no; i++) {
            const AplObs *o = &obs[i];
            if (o->du * o->s <= 0.0) continue;
            res[nr++] = fabs(fabs(o->du) / o->k - (t0 + g * o->w));
        }
        double med = 0.0;
        if (nr) { qsort(res, nr, sizeof(double), apl_cmp_double); med = res[nr/2]; }
        int sound = nr >= (size_t)APL_TURN_MODEL_MIN_PAIRS && t0 > APL_TURN_MIN && t0 < APL_TURN_MAX &&
                    med <= APL_TURN_MODEL_MAX_RESID * t0;
        if (sound) { x->t0 = t0; x->growth = g; x->t0_measured = 1; }
        fprintf(stderr, "  [place] pooled turn model over %zu placed pairs (%zu sign-agreeing): T(w) = %.0f + w x %.1f vox, residual p50 %.0f vox (%.0f%% of a turn): %s\n",
                no, nr, t0, g, med, t0 > 0.0 ? 100.0 * med / t0 : 0.0, sound ? "adopted" : "refused");
        Arena_restore(arena, mark);
        return 0;
    }
    /* 3. component offsets in u: r = (du + d_cb - d_ca) - s k T(w), IRLS + conjugate gradients, primary fixed */
    double *d = ARENA_CALLOC(arena, (size_t)ncomp, sizeof(double));
    double *rhs = ARENA_ALLOC(arena, (size_t)ncomp * sizeof(double)), *rr = ARENA_ALLOC(arena, (size_t)ncomp * sizeof(double));
    double *pp = ARENA_ALLOC(arena, (size_t)ncomp * sizeof(double)), *Ap = ARENA_ALLOC(arena, (size_t)ncomp * sizeof(double));
    double *wgt = ARENA_ALLOC(arena, no * sizeof(double));
    uint8_t *active = ARENA_CALLOC(arena, (size_t)ncomp, 1);
    uint8_t *use = ARENA_CALLOC(arena, no, 1);   /* v2: a pair enters the offset solve only when its sign agrees with the coordinate */
    for (size_t i = 0; i < no; i++) {
        wgt[i] = obs[i].wt;
        use[i] = (obs[i].ca != obs[i].cb && obs[i].du * obs[i].s > 0.0) ? 1 : 0;
        if (use[i]) { active[obs[i].ca] = 1; active[obs[i].cb] = 1; }
    }
    active[primary] = 0;   /* fixed */
    double tbar = t0 + g * 0.0;
    for (int round = 0; round < 4; round++) {
        memset(rhs, 0, (size_t)ncomp * sizeof(double));
        for (size_t i = 0; i < no; i++) {
            const AplObs *o = &obs[i];
            if (!use[i]) continue;
            double target = o->s * o->k * (t0 + g * o->w) - o->du;   /* d_cb - d_ca should equal this */
            if (active[o->cb]) rhs[o->cb] += wgt[i] * target;
            if (active[o->ca]) rhs[o->ca] -= wgt[i] * target;
        }
        /* r = rhs - A d */
        for (int32_t c = 0; c < ncomp; c++) rr[c] = 0.0;
        for (size_t i = 0; i < no; i++) {
            const AplObs *o = &obs[i];
            if (!use[i]) continue;
            double t = wgt[i] * (d[o->cb] - d[o->ca]);
            if (active[o->cb]) rr[o->cb] += t;
            if (active[o->ca]) rr[o->ca] -= t;
        }
        for (int32_t c = 0; c < ncomp; c++) rr[c] = active[c] ? rhs[c] - rr[c] : 0.0;
        memcpy(pp, rr, (size_t)ncomp * sizeof(double));
        double rs = 0.0;
        for (int32_t c = 0; c < ncomp; c++) rs += rr[c] * rr[c];
        double rs0 = rs;
        for (int it = 0; it < 500 && rs > 1e-12 * (rs0 > 0.0 ? rs0 : 1.0); it++) {
            for (int32_t c = 0; c < ncomp; c++) Ap[c] = 0.0;
            for (size_t i = 0; i < no; i++) {
                const AplObs *o = &obs[i];
                if (!use[i]) continue;
                double t = wgt[i] * (pp[o->cb] - pp[o->ca]);
                if (active[o->cb]) Ap[o->cb] += t;
                if (active[o->ca]) Ap[o->ca] -= t;
            }
            /* a tiny ridge keeps components without an active neighbour in place */
            for (int32_t c = 0; c < ncomp; c++) if (active[c]) Ap[c] += 1e-6 * pp[c];
            double pAp = 0.0;
            for (int32_t c = 0; c < ncomp; c++) if (active[c]) pAp += pp[c] * Ap[c];
            if (!(pAp > 0.0)) break;
            double alpha = rs / pAp;
            for (int32_t c = 0; c < ncomp; c++) if (active[c]) { d[c] += alpha * pp[c]; rr[c] -= alpha * Ap[c]; }
            double rs2 = 0.0;
            for (int32_t c = 0; c < ncomp; c++) rs2 += rr[c] * rr[c];
            double beta = rs2 / rs; rs = rs2;
            for (int32_t c = 0; c < ncomp; c++) pp[c] = active[c] ? rr[c] + beta * pp[c] : 0.0;
        }
        double scale = 0.08 * t0;
        for (size_t i = 0; i < no; i++) {
            const AplObs *o = &obs[i];
            double r = (o->du + d[o->cb] - d[o->ca]) - o->s * o->k * (t0 + g * o->w);
            wgt[i] = o->wt / (1.0 + (r / scale) * (r / scale));
        }
    }
    /* 4. damp and apply */
    double cap = APL_METRIC_MAX_SHIFT * tbar, dmax = 0.0, dsum = 0.0; size_t moved = 0, nact = 0;
    for (int32_t c = 0; c < ncomp; c++) {
        if (!active[c]) { d[c] = 0.0; continue; }
        if (d[c] > cap) d[c] = cap;
        if (d[c] < -cap) d[c] = -cap;
        nact++; dsum += fabs(d[c]);
        if (fabs(d[c]) > dmax) dmax = fabs(d[c]);
        if (fabs(d[c]) > 1.0) moved++;
    }
    /* v2: the pass is applied only when the median |residual| over the pairs it used falls */
    double r0 = 0.0, r1 = 0.0;
    {
        size_t nu = 0;
        double *ra = ARENA_ALLOC(arena, (no ? no : 1) * sizeof(double)), *rb = ARENA_ALLOC(arena, (no ? no : 1) * sizeof(double));
        for (size_t i = 0; i < no; i++) {
            const AplObs *o = &obs[i];
            if (!use[i]) continue;
            double model = o->s * o->k * (t0 + g * o->w);
            ra[nu] = fabs(o->du - model);
            rb[nu] = fabs(o->du + d[o->cb] - d[o->ca] - model);
            nu++;
        }
        if (nu > 0) { qsort(ra, nu, sizeof(double), apl_cmp_double); qsort(rb, nu, sizeof(double), apl_cmp_double); r0 = ra[nu / 2]; r1 = rb[nu / 2]; }
    }
    int applied = r1 + 1e-9 < r0;
    if (applied) {
        for (size_t i = 0; i < nc; i++) {
            AsmChart *c = &run->charts[i];
            if (!AsmChart_in_layout(c)) continue;
            int32_t oc = orig_comp[i];
            if (oc < 0 || oc >= ncomp || !active[oc]) continue;
            c->pose_x += d[oc];
        }
        /* the model the later candidates use */
        if (no >= 50) { x->t0 = t0; x->growth = g; x->t0_measured = 1; }
    } else moved = 0;
    fprintf(stderr, "  [place] metric re-solve: %zu placed pairs, turn model T(w) = %.0f + w x %.1f vox, %zu components refit (shift mean %.0f max %.0f vox, %zu moved > 1 vox), residual p50 %.1f -> %.1f: %s\n",
            no, t0, g, nact, nact ? dsum / (double)nact : 0.0, dmax, moved, r0, r1, applied ? "applied" : "not applied");
    Arena_restore(arena, mark);
    /* 5. the audit grids hold the old positions: rebuild them from every placed chart */
    if (copts && garena != NULL && *garena != NULL && moved > 0) {
        if (x->cgrid != NULL) AsmConflictGrid_dispose(x->cgrid);
        Arena_dispose(garena);
        *garena = Arena_new();
        x->cgrid = AsmConflictGrid_new(*garena, copts, run);
        x->rgrid = aprg_new(*garena, 24.0);
        x->strip_umin = 1e300; x->strip_umax = -1e300;
        for (size_t i = 0; i < nc; i++) {
            const AsmChart *c = &run->charts[i];
            if (!AsmChart_in_layout(c) || orig_comp[i] < 0 || comp_placed[orig_comp[i]] != 1) continue;
            AsmConflictGrid_insert(x->cgrid, c); aprg_insert_chart(x->rgrid, c);
            if (c->pose_x < x->strip_umin) x->strip_umin = c->pose_x;
            if (c->pose_x > x->strip_umax) x->strip_umax = c->pose_x;
        }
    }
    return moved;
}

/* ---- coordinate calibration ------------------------------------------------------- */

/* per layer pair: the signed hop count read from the layout when both charts are placed
 * (0 = not placed / not readable) */
static size_t apl_kover_build(AplCtx *x, const uint8_t *comp_placed, const int32_t *orig_comp, const double *lidx, double *kover)
{
    AsmRun *run = x->run;
    Arena_T arena = run->arena;
    size_t n = 0;
    Arena_Mark mark = Arena_save(arena);
    double *tmp = ARENA_ALLOC(arena, 64 * sizeof(double));
    for (size_t p = 0; p < run->n_layers; p++) {
        const AsmLayerPair *lp = &run->layers[p];
        kover[p] = 0.0;
        const AsmChart *A = &run->charts[(size_t)lp->a], *B = &run->charts[(size_t)lp->b];
        if (!AsmChart_in_layout(A) || !AsmChart_in_layout(B)) continue;
        int32_t ca = orig_comp[lp->a], cb = orig_comp[lp->b];
        if (ca < 0 || cb < 0 || comp_placed[ca] != 1 || comp_placed[cb] != 1) continue;
        size_t nh = 0;
        for (int h = 0; h < lp->hit_count && nh < 64; h++) {
            const AsmLayerHit *hit = &run->layer_hits[(size_t)lp->hit_first + (size_t)h];
            double ax, ay, bx, by;
            apl_vertex_point(A, hit->va, &ax, &ay);
            apl_face_point(B, hit->fb, hit->l0, hit->l1, &bx, &by);
            tmp[nh++] = bx - ax;
        }
        if (nh < 2) continue;
        qsort(tmp, nh, sizeof(double), apl_cmp_double);
        double du = tmp[nh / 2];
        /* the turn at the inner of the two: the previous coordinate says which layer that is */
        double w = 0.0;
        if (lidx != NULL && lidx[lp->a] > APL_LAYER_UNKNOWN && lidx[lp->b] > APL_LAYER_UNKNOWN) w = fmin(lidx[lp->a], lidx[lp->b]);
        double T = x->t0 + w * x->growth;
        if (T < APL_TURN_MIN) T = APL_TURN_MIN;
        double q = fabs(du) / T;
        if (q < 0.5 || q > 4.5) continue;
        int k = (int)floor(q + 0.5);
        kover[p] = du < 0.0 ? -(double)k : (double)k;
        n++;
    }
    Arena_restore(arena, mark);
    return n;
}

/* per component: the nearest layer to the primary, from the coordinate */
static void apl_comp_ring_update(const AsmRun *run, const double *lidx, const int32_t *orig_comp, int32_t ncomp, int32_t *comp_ring)
{
    for (int32_t c = 0; c < ncomp; c++) comp_ring[c] = 1000000;
    for (size_t i = 0; i < run->n_charts; i++) {
        if (lidx[i] <= APL_LAYER_UNKNOWN || !AsmChart_in_layout(&run->charts[i])) continue;
        int32_t ring = (int32_t)floor(fabs(lidx[i]) + 0.5), oc = orig_comp[i];
        if (oc >= 0 && oc < ncomp && ring < comp_ring[oc]) comp_ring[oc] = ring;
    }
}

/* ---- turn field ------------------------------------------------------------------ */

typedef struct AplTObs { double u, w, t; } AplTObs;
static int apl_cmp_tobs_cell(const void *x, const void *y)
{
    const AplTObs *a = x, *b = y;
    long ca = (long)floor(a->w / APL_TFIELD_DW), cb = (long)floor(b->w / APL_TFIELD_DW);
    if (ca != cb) return ca < cb ? -1 : 1;
    long ua = (long)floor(a->u / APL_TFIELD_DU), ub = (long)floor(b->u / APL_TFIELD_DU);
    if (ua != ub) return ua < ub ? -1 : 1;
    return (a->t > b->t) - (a->t < b->t);
}

/* rebuild x->tfield: for every placed chart, the median observed one-wrap turn of the placed pairs
 * in the 3x3 (u, layer) cells around it, else 0 */
static size_t apl_turn_field_build(AplCtx *x, const uint8_t *comp_placed, const int32_t *orig_comp, const double *lidx)
{
    AsmRun *run = x->run;
    Arena_T arena = run->arena;
    size_t nc = run->n_charts, filled = 0;
    for (size_t i = 0; i < nc; i++) x->tfield[i] = 0.0;
    if (lidx == NULL) return 0;
    Arena_Mark mark = Arena_save(arena);
    AplTObs *obs = ARENA_ALLOC(arena, (run->n_layers ? run->n_layers : 1) * sizeof(AplTObs));
    size_t no = 0;
    double *tmp = ARENA_ALLOC(arena, 64 * sizeof(double));
    for (size_t p = 0; p < run->n_layers; p++) {
        const AsmLayerPair *lp = &run->layers[p];
        const AsmChart *A = &run->charts[(size_t)lp->a], *B = &run->charts[(size_t)lp->b];
        if (!AsmChart_in_layout(A) || !AsmChart_in_layout(B)) continue;
        int32_t ca = orig_comp[lp->a], cb = orig_comp[lp->b];
        if (ca < 0 || cb < 0 || comp_placed[ca] != 1 || comp_placed[cb] != 1) continue;
        if (lidx[lp->a] <= APL_LAYER_UNKNOWN || lidx[lp->b] <= APL_LAYER_UNKNOWN) continue;
        double dl = lidx[lp->b] - lidx[lp->a];
        if (fabs(dl) < 0.5) continue;
        int k = (int)floor(fabs(dl) + 0.5);
        if (k > 4) continue;
        double sgn = dl < 0.0 ? -1.0 : 1.0;
        size_t nh = 0; double usum = 0.0;
        for (int h = 0; h < lp->hit_count && nh < 64; h++) {
            const AsmLayerHit *hit = &run->layer_hits[(size_t)lp->hit_first + (size_t)h];
            double ax, ay, bx, by;
            apl_vertex_point(A, hit->va, &ax, &ay);
            apl_face_point(B, hit->fb, hit->l0, hit->l1, &bx, &by);
            tmp[nh++] = bx - ax; usum += ax;
        }
        if (nh < 2) continue;
        qsort(tmp, nh, sizeof(double), apl_cmp_double);
        double du = tmp[nh / 2];
        if (du * sgn <= 0.0) continue;   /* the sign must agree with the coordinate */
        double t = fabs(du) / (double)k;
        if (t < APL_TURN_MIN || t > APL_TURN_MAX) continue;
        obs[no].u = usum / (double)nh; obs[no].w = 0.5 * (lidx[lp->a] + lidx[lp->b]) - 0.5 * sgn * (double)k; /* the INNER layer's coordinate */
        obs[no].w = lidx[lp->a] < lidx[lp->b] ? lidx[lp->a] : lidx[lp->b];
        obs[no].t = t; no++;
    }
    if (no < (size_t)APL_TFIELD_MIN_OBS) { Arena_restore(arena, mark); return 0; }
    qsort(obs, no, sizeof(AplTObs), apl_cmp_tobs_cell);
    /* cell index: first observation of each (w cell, u cell) run */
    size_t ncell = 0;
    for (size_t i = 0; i < no; i++) if (i == 0 || apl_cmp_tobs_cell(&obs[i - 1], &obs[i]) != 0 || (long)floor(obs[i-1].w / APL_TFIELD_DW) != (long)floor(obs[i].w / APL_TFIELD_DW) || (long)floor(obs[i-1].u / APL_TFIELD_DU) != (long)floor(obs[i].u / APL_TFIELD_DU)) ncell++;
    size_t *cstart = ARENA_ALLOC(arena, (ncell + 1) * sizeof(size_t));
    long *cw = ARENA_ALLOC(arena, (ncell ? ncell : 1) * sizeof(long)), *cu = ARENA_ALLOC(arena, (ncell ? ncell : 1) * sizeof(long));
    size_t q = 0;
    for (size_t i = 0; i < no; i++) {
        long wi = (long)floor(obs[i].w / APL_TFIELD_DW), ui = (long)floor(obs[i].u / APL_TFIELD_DU);
        if (i == 0 || wi != cw[q - 1] || ui != cu[q - 1]) { cstart[q] = i; cw[q] = wi; cu[q] = ui; q++; }
    }
    cstart[q] = no; ncell = q;
    double *vals = ARENA_ALLOC(arena, (no ? no : 1) * sizeof(double));
    for (size_t i = 0; i < nc; i++) {
        const AsmChart *c = &run->charts[i];
        if (!AsmChart_in_layout(c) || orig_comp[i] < 0 || comp_placed[orig_comp[i]] != 1 || lidx[i] <= APL_LAYER_UNKNOWN) continue;
        double gx, gy;
        apl_pose(c, 0.0, 0.0, &gx, &gy);
        long wi = (long)floor(lidx[i] / APL_TFIELD_DW), ui = (long)floor(gx / APL_TFIELD_DU);
        size_t nv = 0;
        for (size_t k = 0; k < ncell; k++) {
            if (cw[k] < wi - 1 || cw[k] > wi + 1 || cu[k] < ui - 1 || cu[k] > ui + 1) continue;
            for (size_t j = cstart[k]; j < cstart[k + 1]; j++) vals[nv++] = obs[j].t;
        }
        if (nv < (size_t)APL_TFIELD_MIN_OBS) continue;
        qsort(vals, nv, sizeof(double), apl_cmp_double);
        x->tfield[i] = vals[nv / 2];
        filled++;
    }
    Arena_restore(arena, mark);
    return filled;
}

/* A frozen-neighbour closure uses the actual boundary targets. The three
 * unknowns are a rigid correction of the candidate; certified charts never
 * move. Centering keeps the normal system well conditioned at scroll scale. */
static int apl_close_seam(AplCtx *x, double *th, double *tx, double *ty)
{
    if (x->cluster_kind != 1) return 1;
    double cx = 0, cy = 0, sw = 0;
    for (size_t i = 0; i < x->na; i++) if (x->wts[i] > 0 && x->anchors[i].k == 0) {
        cx += x->wts[i]*x->anchors[i].sx; cy += x->wts[i]*x->anchors[i].sy; sw += x->wts[i];
    }
    if (sw < APL_SEAM_MIN_ANCHORS) return 0;
    cx /= sw; cy /= sw;
    double gx = cos(*th)*cx-sin(*th)*cy+*tx, gy = sin(*th)*cx+cos(*th)*cy+*ty;
    const int rows[6] = {0,1,1,2,2,2}, cols[6] = {0,0,1,0,1,2};
    for (int it = 0; it < 5; it++) {
        double H[6] = {0}, b[3] = {0}, delta[3] = {0};
        double c = cos(*th), s = sin(*th);
        for (size_t i = 0; i < x->na; i++) {
            const AplAnchor *a = &x->anchors[i];
            if (x->wts[i] <= 0 || a->k != 0) continue;
            double u = a->sx-cx, v = a->sy-cy, px = c*u-s*v, py = s*u+c*v;
            double ex = a->tx-px-gx, ey = a->ty-py-gy, w = x->wts[i];
            H[0] += w; H[2] += w; H[3] -= w*py; H[4] += w*px; H[5] += w*(px*px+py*py);
            b[0] += w*ex; b[1] += w*ey; b[2] += w*(-py*ex+px*ey);
        }
        if (Sparse_solve_sym(3, 6, rows, cols, H, b, delta, SPARSE_SPD) != 0) return 0;
        if (!isfinite(delta[2]) || fabs(*th+delta[2]) > x->o->max_rot_refine) return 0;
        gx += delta[0]; gy += delta[1]; *th += delta[2];
        if (hypot(delta[0],delta[1]) < 1e-6 && fabs(delta[2]) < 1e-8) break;
    }
    *tx = gx-cos(*th)*cx+sin(*th)*cy; *ty = gy-sin(*th)*cx-cos(*th)*cy;
    return 1;
}

/* Measure the actual trial geometry. A pose which already meets these
 * constraints is a valid incumbent; optional refinement must not veto it
 * merely because an unconstrained Gauss-Newton step is too large. */
static int apl_trial_support(AplCtx *x, int32_t best, const AsmContinuityTarget *targets,
                             size_t n_targets, size_t *support, int32_t *parent, double *support_rms, int report)
{
    AsmRun *run = x->run;
    *support = 0; *parent = -1; *support_rms = 0;
    double e2 = 0;
    if (x->cluster_kind == 0) for (size_t t = 0; t < n_targets; t++) {
        const AsmContinuityTarget *a = &targets[t];
        if (run->charts[a->parent].placement_state == ASM_PLACE_NONE) continue;
        double gx, gy; apl_pose(&run->charts[a->chart], a->u, a->v, &gx, &gy);
        double e = hypot(gx-a->x, gy-a->y);
        if (e <= 16.0) { (*support)++; *parent = a->parent; e2 += e*e; }
    }
    for (size_t r = 0; r < run->n_rels; r++) {
        const AsmRelation *R = &run->rels[r];
        if (!(R->continuity & ASM_CONT_SOURCE) || (R->flags & (ASM_REL_CONTACT | ASM_REL_DROPPED))) continue;
        if (!AsmRel_is_join(R) && !(R->continuity & ASM_CONT_CUT) && (x->cluster_kind != 1 || !x->seam_selected[r])) continue;
        int ca = x->orig_comp[R->a] == best, cb = x->orig_comp[R->b] == best;
        if (!ca && !cb) continue;
        const AsmChart *fixed = &run->charts[ca ? R->b : R->a];
        if (ca != cb && (fixed->placement_state == ASM_PLACE_NONE || fixed->component != x->primary)) continue;
        double rms; size_t n;
        if (!AsmContinuity_measure(run, R, &rms, &n)) {
            x->last_reason = "source_continuation_disagrees";
            if (report) fprintf(stderr, "  [place] source continuation %d-%d disagrees (%.2f vox, %zu pairs); transaction deferred\n", R->a, R->b, rms, n);
            return 0;
        }
        if (ca != cb && x->cluster_kind == 1) { *support += n; *parent = fixed->id; e2 += n*rms*rms; }
    }
    if (*support) *support_rms = sqrt(e2 / *support);
    if (*support < (x->cluster_kind == 1 ? (size_t)APL_SEAM_MIN_ANCHORS : (size_t)x->o->min_anchors))
        { x->last_reason = "trial_support_short"; return 0; }
    return 1;
}

/* Axis statistics of the actual candidate geometry, after composing its
 * pose and any boundary correction. The original poses identify its charts. */
static double apl_candidate_axis(const AplCtx *x, const AsmChart *saved, size_t n, double *rms)
{
    double sw = 0, sr = 0, sr2 = 0;
    *rms = 0;
    if (!APL_V_FROM_AXIS || !x->vm_ok || !x->vm_pin) return 0;
    for (size_t k = 0; k < n; k++) {
        const AsmChart *c = &x->run->charts[saved[k].id];
        double u, v;
        if (c->uv != saved[k].uv) {
            double cu, cv; apl_uv_centroid(c,&cu,&cv); apl_pose(c,cu,cv,&u,&v);
        } else apl_layout_centroid(x, c, c->id, 0, 0, 0, &u, &v);
        double e = x->vm_a + x->vm_b*x->s_axis[c->id] + x->vm_c*u - v;
        sw += c->area3d; sr += c->area3d*e; sr2 += c->area3d*e*e;
    }
    if (!(sw > 0)) return 0;
    double shift = sr/sw;
    *rms = sqrt(fmax(0.0, sr2/sw-shift*shift));
    return shift;
}

/* Validate the final transaction, not the rotation-limited proposal that
 * preceded seam closure. A repair may fix a bad proposal or move a good
 * one outside the permitted shape; both must be judged here. */
static int apl_candidate_shape(AplCtx *x, const AsmChart *saved, size_t n, double axis_pin)
{
    if (!n) return 0;
    if (APL_V_FROM_AXIS && x->vm_ok && x->vm_pin) {
        double vrms, dv = apl_candidate_axis(x, saved, n, &vrms);
        double anchor_shift = 0, anchor_weight = 0;
        if (x->cluster_kind == 0) for (size_t i = 0; i < x->na; i++) if (x->wts[i] > 0 && x->anchors[i].k > 0) {
            const AplAnchor *a = &x->anchors[i]; double u, v;
            apl_pose(&x->run->charts[a->chart_c], a->local_u, a->local_v, &u, &v);
            anchor_shift += x->wts[i]*(a->ty-v); anchor_weight += x->wts[i];
        }
        if (anchor_weight > 0) anchor_shift /= anchor_weight;
        if (!isfinite(vrms) || !isfinite(dv) || !isfinite(anchor_shift)) return 0;
        if (vrms > fmax(APL_V_AXIS_TOL, 1.5*x->vm_p90) ||
            fmax(fmax(fabs(dv), fabs(axis_pin)), fabs(anchor_shift)) > fmax(APL_V_AXIS_MAX_SHIFT, 2*x->vm_p90)) {
            fprintf(stderr, "  [place]   final v-axis gate refuses: scatter %.1f vox rms, residual shift %.1f, applied pin %.1f, anchor shift %.1f vox\n", vrms, dv, axis_pin, anchor_shift);
            if (x->diagnostic_axis_policy == 2 ||
                (x->diagnostic_axis_policy == 1 && x->cluster_kind == 1)) {
                fprintf(stderr, "  [place axis diagnostic] checking the remaining gates, kind %d; no production policy change\n", x->cluster_kind);
            } else {
                if (x->st) x->st->v_refused++;
                return 0;
            }
        }
    }
    double ulo = 1e300, uhi = -1e300, vlo = 1e300, vhi = -1e300, plo = 1e300, phi = -1e300;
    for (size_t k = 0; k < n; k++) {
        const AsmChart *c = &x->run->charts[saved[k].id];
        if (!isfinite(c->pose_x) || !isfinite(c->pose_y) || !isfinite(c->pose_theta)) return 0;
        ulo = fmin(ulo, c->pose_x); uhi = fmax(uhi, c->pose_x);
        vlo = fmin(vlo, c->pose_y); vhi = fmax(vhi, c->pose_y);
        plo = fmin(plo, saved[k].pose_y); phi = fmax(phi, saved[k].pose_y);
    }
    double span = fmax(uhi, x->strip_umax)-fmin(ulo, x->strip_umin);
    if (x->span_cap > 0 && span > x->span_cap) {
        fprintf(stderr, "  [place]   final span gate refuses: %.0f vox (cap %.0f)\n", span, x->span_cap);
        return 0;
    }
    double tv = 0, tw = 0;
    for (size_t i = 0; i < x->na; i++) if (x->wts[i] > 0) { tv += x->wts[i]*x->anchors[i].ty; tw += x->wts[i]; }
    if (tw > 0) {
        double extent = phi-plo+300, tvm = tv/tw;
        if (vlo > tvm+extent || vhi < tvm-extent || vhi-vlo > 2*extent+600) {
            fprintf(stderr, "  [place]   final swing gate refuses: body v %.0f..%.0f, anchors v %.0f, own extent %.0f\n", vlo, vhi, tvm, phi-plo);
            return 0;
        }
    }
    return 1;
}

/* Audit of the entire incident layer evidence. A ray is counted once,
 * after trying its possible wrap counts. Cluster membership and fit trimming
 * cannot hide contrary observations or multiply an observation's support. */
static AplConsensus apl_candidate_consensus(const AplCtx *x)
{
    AplConsensus result = {0};
    if (!x->n_observations) return result;
    Arena_Mark mark = Arena_save(x->run->arena);
    double *error = ARENA_ALLOC(x->run->arena, x->n_observations*sizeof(double));
    double *model = ARENA_ALLOC(x->run->arena, x->n_observations*sizeof(double));
    size_t *ref = ARENA_ALLOC(x->run->arena, x->n_observations*sizeof(size_t));
    for (size_t i = 0; i < x->n_observations; i++) { error[i] = model[i] = 1e300; ref[i] = SIZE_MAX; }
    for (size_t i = 0; i < x->na; i++) {
        const AplAnchor *a = &x->anchors[i];
        if (a->k <= 0 || a->measured_order <= 0 || a->measured_order > a->k || a->observation >= x->n_observations) continue;
        if (x->run->charts[a->chart_p].placement_state == ASM_PLACE_NONE) continue;
        const AsmChart *c = &x->run->charts[a->chart_c]; double u, v;
        if (x->uv_trial && a->hit) {
            if (a->reverse) apl_vertex_point(c,a->hit->va,&u,&v);
            else apl_face_point(c,a->hit->fb,a->hit->l0,a->hit->l1,&u,&v);
        } else apl_pose(c,a->local_u,a->local_v,&u,&v);
        double e = hypot(u-a->tx, v-a->ty);
        double scale = fmax(16.0, 0.05*hypot(a->Tx,a->Ty)/(double)a->k);
        error[a->observation] = fmin(error[a->observation], e);
        model[a->observation] = fmin(model[a->observation], e/scale);
        ref[a->observation] = i;
    }
    for (size_t i = 0; i < x->n_observations; i++) if (error[i] < 1e299) {
        double e = error[i]/16.0, z = model[i];
        result.observations++; result.inliers += error[i] <= 16.0;
        result.model_inliers += z <= 1.0;
        result.clipped_loss += fmin(1.0, e*e);
        result.model_loss += z <= 1.0 ? 0.5*z*z : z-0.5;
    }
    if (x->fevidence) {
        typedef struct AplEvidence { size_t n, inliers; double loss; } AplEvidence;
        size_t nc = x->run->n_charts;
        AplEvidence *groups = ARENA_CALLOC(x->run->arena, 2*nc, sizeof(AplEvidence));
        int32_t component = -1;
        for (size_t i = 0; i < x->n_observations; i++) if (ref[i] != SIZE_MAX) {
            const AplAnchor *a = &x->anchors[ref[i]];
            component = x->orig_comp[a->chart_c];
            for (int side = 0; side < 2; side++) {
                AplEvidence *g = &groups[(side ? nc : 0)+(size_t)(side ? a->chart_c : a->chart_p)];
                g->n++; g->inliers += error[i] <= 16.0;
                double z = model[i]; g->loss += z <= 1.0 ? 0.5*z*z : z-0.5;
            }
        }
        for (size_t i = 0; i < 2*nc; i++) if (groups[i].n)
            fprintf(x->fevidence, "%d,%d,%d,%s,%zu,%zu,%zu,%.9g\n", x->attempt, component, x->cluster,
                    i < nc ? "parent" : "candidate", i%nc, groups[i].n, groups[i].inliers, groups[i].loss);
    }
    Arena_restore(x->run->arena, mark);
    return result;
}

static int apl_relax_boundary_trial(AplCtx *x, const AsmChart *saved, size_t n)
{
    return x->diagnostic_first_step
        ? AsmContinuity_relax_boundary_step(x->run,saved,n,x->seam_selected,apl_uv_audit,x->o->stress_band)
        : AsmContinuity_relax_boundary_trial(x->run,saved,n,x->seam_selected,apl_uv_audit,x->o->stress_band);
}

static int apl_place_and_audit_once(AplCtx *x, int32_t best, double th, double tx, double ty,
                                    double *out_ratio, int require_rim, int polish, int *did_polish)
{
    *did_polish = 0;
    if (out_ratio) *out_ratio = 0.0;
    x->last_silent = 0; x->last_shape_refused = 0; x->last_reason = "accepted";
    x->uv_trial = 0;
    if (x->proposal_uv) Arena_free(x->proposal_uv);
    int32_t parent = -1;
    size_t support = 0;
    double support_rms = 0;
    size_t min_support = x->cluster_kind == 1 ? APL_SEAM_MIN_ANCHORS : (size_t)x->o->min_anchors;
    AsmRun *run = x->run;
    Arena_T arena = run->arena;
    size_t nc = run->n_charts;
    Arena_Mark m2 = Arena_save(arena);
    size_t nsave = 0;
    for (size_t i = 0; i < nc; i++) if (AsmChart_in_layout(&run->charts[i]) && run->charts[i].component == best) nsave++;
    AsmChart *saved = ARENA_ALLOC(arena, (nsave ? nsave : 1) * sizeof(AsmChart));
    AsmContinuityTarget *targets = ARENA_ALLOC(arena, (x->na ? x->na : 1)*sizeof(AsmContinuityTarget));
    size_t n_targets = 0;
    if (x->cluster_kind == 0) for (size_t i = 0; i < x->na; i++) {
        const AplAnchor *a = &x->anchors[i];
        if (a->cluster != x->cluster || x->wts[i] <= 0 || a->k <= 0 || a->measured_order <= 0 || a->measured_order > a->k) continue;
        targets[n_targets++] = (AsmContinuityTarget){ a->chart_c, a->local_u, a->local_v, a->tx, a->ty, 0.05*x->wts[i], a->chart_p };
    }
    size_t q = 0;
    for (size_t i = 0; i < nc; i++) {
        AsmChart *c = &run->charts[i];
        if (!AsmChart_in_layout(c) || c->component != best) continue;
        saved[q++] = *c;
        apl_compose(c, th, tx, ty, x->source_flip);
        c->component = x->primary;
    }
    double axis_pin = 0;
    if (APL_V_FROM_AXIS && x->vm_ok && x->vm_pin) {
        if (x->cluster_kind == 1 && !APL_SEAM_V_PIN) { if (x->st) x->st->v_seam_kept++; }
        else {
            double vrms;
            axis_pin = apl_axis_excess(apl_candidate_axis(x, saved, nsave, &vrms), apl_axis_band(x));
            for (size_t k = 0; k < nsave; k++) run->charts[saved[k].id].pose_y += axis_pin;
            if (x->st) {
                x->st->v_pinned++;
                x->st->v_shift_max = fmax(x->st->v_shift_max, fabs(axis_pin));
            }
        }
    }
    int ok = apl_trial_support(x, best, targets, n_targets, &support, &parent, &support_rms, 0);
    size_t cand_vertices = 0;
    for (size_t k = 0; k < nsave; k++) cand_vertices += saved[k].nv;
    if (!ok && x->cluster_kind == 0 && support < (size_t)APL_HOP_REPAIR_MIN_SUPPORT && cand_vertices > (size_t)APL_HOP_REPAIR_MAX_VERTICES) {
        fprintf(stderr, "  [place] hop candidate: component %d, initial support %zu/%zu below %d; no joint repair\n", best, support, min_support, APL_HOP_REPAIR_MIN_SUPPORT);
    } else if (!ok) {
        fprintf(stderr, "  [place] joint repair: component %d, initial support %zu/%zu, %zu charts, %zu targets\n", best, support, min_support, nsave, n_targets);
        double *repair_start = ARENA_ALLOC(arena,3*(nsave ? nsave : 1)*sizeof *repair_start);
        for (size_t k = 0; k < nsave; k++) {
            const AsmChart *c = &run->charts[saved[k].id];
            repair_start[3*k] = c->pose_x; repair_start[3*k+1] = c->pose_y; repair_start[3*k+2] = c->pose_theta;
        }
        AsmChart *fit_start = NULL;
        if (x->cluster_kind == 1) {
            fit_start = ARENA_ALLOC(arena, (nsave ? nsave : 1)*sizeof(AsmChart));
            for (size_t k = 0; k < nsave; k++) fit_start[k] = run->charts[saved[k].id];
        }
        int solved = AsmContinuity_close(run, saved, nsave, x->seam_selected, targets, n_targets);
        /* A SOURCE repair may stop after earlier useful moves. Certify its
         * actual pose independently, including the original cumulative
         * bounds: a failed solve can also leave a move beyond those bounds.
         * Model-based hop proposals still require numerical completion;
         * accepting their partial solves regresses final overlap (dev223).
         * An accepted partial SOURCE repair receives every geometry audit. */
        int repair_bounded = 1;
        for (size_t k = 0; k < nsave && repair_bounded; k++) {
            const AsmChart *c = &run->charts[saved[k].id];
            if (!isfinite(c->pose_x) || !isfinite(c->pose_y) || !isfinite(c->pose_theta) ||
                hypot(c->pose_x-repair_start[3*k],c->pose_y-repair_start[3*k+1]) > 256 ||
                fabs(c->pose_theta-repair_start[3*k+2]) > .25) repair_bounded = 0;
        }
        if (!solved && x->cluster_kind != 1) x->last_reason = "repair_unsolved";
        else if (!repair_bounded) x->last_reason = "repair_unbounded";
        ok = (solved || x->cluster_kind == 1) && repair_bounded &&
             apl_trial_support(x, best, targets, n_targets, &support, &parent, &support_rms, solved);
        if (!solved && ok) fprintf(stderr,"  [place] joint solver stopped for component %d; bounded current SOURCE pose passes support, checking geometry\n",best);
        if (solved) {
            if (!ok && fit_start) {
                /* Restart at the same proposal so the fallback has the same
                 * total movement bounds and cannot compound two repairs. */
                for (size_t k = 0; k < nsave; k++) run->charts[saved[k].id] = fit_start[k];
                fprintf(stderr, "  [place] boundary recovery: component %d, %zu charts\n", best, nsave);
                int boundary_solved = AsmContinuity_close_boundaries(run, saved, nsave, x->seam_selected);
                if (boundary_solved)
                    ok = apl_trial_support(x, best, targets, n_targets, &support, &parent, &support_rms, 1);
                if (!ok && boundary_solved) {
                    size_t bytes = 0;
                    for (size_t k = 0; k < nsave; k++) bytes += saved[k].nv*2*sizeof(float);
                    if (bytes <= 268435456 && cand_vertices > (size_t)APL_RELAX_TRIAL_MAX_VERTICES)
                        fprintf(stderr,"  [place] local boundary UV repair skipped: component %d, %zu vertices above %d\n",best,cand_vertices,APL_RELAX_TRIAL_MAX_VERTICES);
                    else if (bytes <= 268435456) {
                        for (size_t k = 0; k < nsave; k++) {
                            AsmChart *c = &run->charts[saved[k].id];
                            c->uv = ARENA_ALLOC(arena,c->nv*2*sizeof(float));
                            memcpy(c->uv,saved[k].uv,c->nv*2*sizeof(float));
                        }
                        AsmChart *previous = ARENA_ALLOC(arena,nsave*sizeof *previous);
                        float *previous_uv = ARENA_ALLOC(arena,bytes);
                        for (int uv_pass = 0; uv_pass < 3; uv_pass++) {
                            size_t offset = 0;
                            for (size_t k = 0; k < nsave; k++) {
                                previous[k] = run->charts[saved[k].id];
                                previous[k].uv = previous_uv+offset;
                                memcpy(previous[k].uv,run->charts[saved[k].id].uv,saved[k].nv*2*sizeof(float));
                                offset += saved[k].nv*2;
                            }
                            if (!apl_relax_boundary_trial(x,saved,nsave)) break;
                            int next_ok = apl_trial_support(x,best,targets,n_targets,&support,&parent,&support_rms,1);
                            int within_bound = 1;
                            if (!next_ok) {
                                within_bound = AsmContinuity_close_boundaries(run,saved,nsave,x->seam_selected);
                                if (within_bound) next_ok = apl_trial_support(x,best,targets,n_targets,&support,&parent,&support_rms,1);
                                for (size_t k = 0; k < nsave && within_bound; k++) {
                                    const AsmChart *c = &run->charts[saved[k].id], *start = &fit_start[k];
                                    if (hypot(c->pose_x-start->pose_x,c->pose_y-start->pose_y) > 256 ||
                                        fabs(c->pose_theta-start->pose_theta) > 0.25) {
                                        fprintf(stderr,"  [place] boundary repair exceeds original movement bound: chart %d\n",c->id);
                                        within_bound = 0;
                                    }
                                }
                            }
                            if (!within_bound || (ok && !next_ok)) {
                                for (size_t k = 0; k < nsave; k++) {
                                    AsmChart *c = &run->charts[saved[k].id]; float *owned_uv = c->uv;
                                    memcpy(owned_uv,previous[k].uv,c->nv*2*sizeof(float));
                                    *c = previous[k]; c->uv = owned_uv;
                                }
                                if (ok) apl_trial_support(x,best,targets,n_targets,&support,&parent,&support_rms,0);
                                break;
                            }
                            x->uv_trial = 1; ok = next_ok;
                            fprintf(stderr,"  [place] local boundary UV repair: component %d, pass %d, source constraints %s\n",best,uv_pass+1,ok ? "pass" : "fail");
                        }
                    }
                }
            }
        }
        if (ok) {
            /* Common pose composition and axis pinning are rigid. The
             * original checkpoint has the same internal geometry, so its
             * audit is needed only after a successful repair. */
            AsmContinuityAuditCache internal = {0}; internal.scratch = arena;
            size_t internal_after, introduced = apl_internal_new(run,saved,nsave,&internal,&internal_after);
            if (introduced) {
                fprintf(stderr, "  [place] joint repair creates internal contradictions: %zu -> %zu vertices, %zu previously sound vertices violated; transaction deferred\n", internal.baseline_bad, internal_after,introduced);
                ok = 0; x->last_reason = "internal_contradiction";
            }
        }
        if (!ok) x->last_silent = 1;
    }
    if (ok && polish && x->cluster_kind == 1 && nsave && x->copts && x->cgrid) {
        size_t guaranteed, queries;
        if (aprg_polish_blocked(x->rgrid,run,saved,nsave,x->excluded_samples,x->o->max_rim_ratio,&guaranteed,&queries)) {
            /* Keep the ordinary audits of the original supported candidate.
             * In particular they must still compute out_ratio for callers. */
            polish = 0;
            fprintf(stderr,"  [place] UV polish cannot clear rim: component %d, at least %zu of %zu queries violate throughout the original movement disk; auditing original candidate\n",
                    best,guaranteed,queries);
        }
    }
    if (ok && polish && x->cluster_kind == 1 && nsave) {
        /* Passing the broad placement gate does not minimize a boundary's
         * error. Correct already supported seams as well as failed ones;
         * otherwise their residual becomes a slit or stretched raster join.
         * Stage independent UV buffers so any later audit can roll back.
         * The outer wrapper retries the original incumbent if this optional
         * improvement fails its full geometry audit. */
        size_t bytes = 0;
        for (size_t k = 0; k < nsave; k++) bytes += saved[k].nv*2*sizeof(float);
        if (bytes <= 268435456 && cand_vertices > (size_t)APL_RELAX_TRIAL_MAX_VERTICES)
            fprintf(stderr,"  [place] supported boundary correction skipped: component %d, %zu vertices above %d\n",best,cand_vertices,APL_RELAX_TRIAL_MAX_VERTICES);
        else if (bytes <= 268435456) {
            AsmChart *incumbent = ARENA_ALLOC(arena,nsave*sizeof *incumbent);
            for (size_t k = 0; k < nsave; k++) {
                AsmChart *c = &run->charts[saved[k].id]; incumbent[k] = *c;
                c->uv = ARENA_ALLOC(arena,c->nv*2*sizeof(float));
                memcpy(c->uv,incumbent[k].uv,c->nv*2*sizeof(float));
            }
            for (int pass = 0; pass < 3; pass++) {
                if (!apl_relax_boundary_trial(x,saved,nsave)) break;
                *did_polish = 1; x->uv_trial = 1;
            }
            if (*did_polish) {
                ok = apl_trial_support(x,best,targets,n_targets,&support,&parent,&support_rms,1);
                fprintf(stderr,"  [place] supported boundary correction: component %d, source constraints %s, RMS %.3f\n",best,ok ? "pass" : "fail",support_rms);
            } else for (size_t k = 0; k < nsave; k++) run->charts[saved[k].id] = incumbent[k];
        }
    }
    if (ok && !AsmContinuity_preserves_incumbent_trial(run,saved,nsave)) { ok = 0; x->last_reason = "incumbent_changed"; }
    if (ok && !apl_candidate_shape(x, saved, nsave, axis_pin)) { ok = 0; x->last_shape_refused = 1; x->last_reason = "shape_refused"; }
    if (ok) {
        x->consensus = apl_candidate_consensus(x);
        /* A small selected fit cannot establish a wrap while contradicting
         * most of the incident evidence. Test unique physical observations
         * over all permitted wrap counts, with the existing model uncertainty.
         * Independent SOURCE support takes precedence over this layer model.
         * Two poses consistent with most observations remain ambiguous; this
         * gate does not select the pose with the largest supporting count.
         * Run before any persistent grid insertion, including direct commits. */
        if (x->cluster_kind == 0 && x->n_observations &&
            x->consensus.model_inliers <= x->consensus.observations/2) {
            ok = 0; x->last_silent = 1; x->last_reason = "consensus_minority";
            fprintf(stderr,"  [place] component %d: hop model agrees with only %zu/%zu unique incident observations; deferred\n",
                    best,x->consensus.model_inliers,x->consensus.observations);
        }
    }
    int source_supported = ok && x->cluster_kind == 1;
    if (ok && x->copts && x->cgrid != NULL) {
        /* the candidate's charts in their tentative poses against the persistent grids */
        const AsmChart **cand = ARENA_ALLOC(arena, (nsave ? nsave : 1) * sizeof(AsmChart *));
        for (size_t k = 0; k < nsave; k++) cand[k] = &run->charts[saved[k].id];
        int32_t p1 = -1; double m1 = 0.0;
        double tp0 = ves_clock_sec();
        double m = AsmConflictGrid_probe_excluding(x->cgrid, cand, nsave, x->excluded_samples, &p1, &m1);
        double tp1 = ves_clock_sec();
        double ratio = x->comp_area[best] > 0.0 ? m / x->comp_area[best] : 0.0;
        if (out_ratio) *out_ratio = ratio;
        if (ratio > x->o->max_contra_ratio) {
            ok = 0; x->last_reason = "contradiction_mass";
            fprintf(stderr, "  [place]   audit rejects: contradiction mass %.3e = %.2f of its area; partner chart %d (component %d, %.2e)\n",
                    m, ratio, p1, p1 >= 0 ? x->orig_comp[p1] : -1, m1);
        }
        if (ok) {
            size_t n_rim = 0;
            double bad = aprg_probe_excluding(x->rgrid, run, cand, nsave, 24.0, &n_rim, x->excluded_samples);
            if (n_rim >= 8 && bad > x->o->max_rim_ratio) {
                ok = 0; x->last_reason = "rim_isometry";
                apl_write_rim_failure(x,best,cand,nsave);
                fprintf(stderr, "  [place]   rim audit rejects: %.0f%% of %zu rim vertices are farther in 3-D than in (u,v) by more than half a layer\n", 100.0 * bad, n_rim);
            } else if (require_rim == 1 && (n_rim < 8 || bad > 0.1)) {
                ok = 0; x->last_reason = "rim_unconfirmed";
                fprintf(stderr, "  [place]   rim audit does not confirm the retry: %zu rim vertices, %.0f%% violating\n", n_rim, 100.0 * bad);
            } else if (support < min_support) {
                ok = 0; x->last_silent = 1; x->last_reason = "insufficient_support";
                fprintf(stderr, "  [place]   insufficient positive support: %zu measured pairs, %zu rim vertices; deferred\n", support, n_rim);
            } else if (n_rim >= 8 && bad > 0.0)
                fprintf(stderr, "  [place]   rim audit: %.0f%% of %zu rim vertices violate\n", 100.0 * bad, n_rim);
        }
        double tp2 = ves_clock_sec();
        if (ok) {
            for (size_t k = 0; k < nsave; k++) {
                AsmChart *c = &run->charts[saved[k].id];
                c->placement_state = x->cluster_kind == 1 ? ASM_PLACE_SEAM : ASM_PLACE_LAYER;
                c->placement_parent = parent; c->placement_support = support;
                c->placement_rms = support_rms;
            }
            if (!x->trial) {
                AsmContinuity_confirm(run);
                for (size_t k = 0; k < nsave; k++) { AsmConflictGrid_insert(x->cgrid, cand[k]); aprg_insert_chart(x->rgrid, cand[k]); }
            }
        }
        double tp3 = ves_clock_sec();
        if (x->st) { x->st->t_probe += tp1 - tp0; x->st->t_rim += tp2 - tp1; x->st->t_insert += tp3 - tp2; }
        fprintf(stderr, "  [place]   audit timing: conflict probe %.2f s, rim probe %.2f s, insert %.2f s (grid cells %zu)\n",
                tp1 - tp0, tp2 - tp1, tp3 - tp2, AsmConflictGrid_cells(x->cgrid));
    }
    if (!ok && source_supported) x->source_geometry_conflict = 1;
    if (ok && !x->cgrid) {
        ok = support >= min_support;
        if (!ok) { x->last_silent = 1; x->last_reason = "insufficient_support"; }
        else for (size_t k = 0; k < nsave; k++) {
            AsmChart *c = &run->charts[saved[k].id];
            c->placement_state = x->cluster_kind == 1 ? ASM_PLACE_SEAM : ASM_PLACE_LAYER;
            c->placement_parent = parent; c->placement_support = support;
            c->placement_rms = support_rms;
        }
    }
    if (ok && x->trial) {
        x->nproposal = nsave;
        for (size_t k = 0; k < nsave; k++) {
            x->proposal[k] = run->charts[saved[k].id];
            if (x->proposal[k].uv != saved[k].uv) {
                if (!x->proposal_uv) x->proposal_uv = Arena_new();
                x->proposal[k].uv = ARENA_ALLOC(x->proposal_uv,saved[k].nv*2*sizeof(float));
                memcpy(x->proposal[k].uv,run->charts[saved[k].id].uv,saved[k].nv*2*sizeof(float));
            }
        }
    }
    if (ok && !x->trial) for (size_t k = 0; k < nsave; k++) apl_commit_chart(x,&run->charts[saved[k].id],saved[k].uv);
    if (!ok || x->trial) for (size_t k = 0; k < nsave; k++) run->charts[saved[k].id] = saved[k];
    Arena_restore(arena, m2);
    return ok;
}

static int apl_audit_proposal(AplCtx *x, int32_t best, double th, double tx, double ty, double *out_ratio, int require_rim)
{
    /* Failing geometry at a supported physical boundary does not supply
     * evidence for a different wrap. A hop previously ignored the selected
     * weak SOURCE proposals and could occupy a later continuation's space.
     * Keep this conflict pending; another SOURCE pose can still resolve it. */
    if (x->cluster_kind == 0 && x->source_geometry_conflict) {
        if (out_ratio) *out_ratio = 0;
        x->last_silent = 1; x->last_shape_refused = 0; x->last_reason = "source_geometry_conflict";
        fprintf(stderr,"  [place] component %d: supported boundary fit has an unresolved geometry conflict; alternate hop deferred\n",best);
        return 0;
    }
    int refined = 0;
    int ok = apl_place_and_audit_once(x,best,th,tx,ty,out_ratio,require_rim,1,&refined);
    if (!ok && refined) {
        fprintf(stderr,"  [place] boundary correction failed geometry audit; retrying original supported candidate\n");
        ok = apl_place_and_audit_once(x,best,th,tx,ty,out_ratio,require_rim,0,&refined);
    }
    return ok;
}

static int apl_place_and_audit(AplCtx *x, int32_t best, double th, double tx, double ty, double *out_ratio, int require_rim)
{
    /* Preserve the input pose independently of optional closure, which can
     * fail after partial steps or succeed at a pose that fails geometry.
     * Each proposal receives the same SOURCE, shape and conflict audits,
     * with and without UV polish when needed. Failed trials roll back. */
    double closed_th = th, closed_tx = tx, closed_ty = ty;
    if (apl_close_seam(x,&closed_th,&closed_tx,&closed_ty)) {
        if (apl_audit_proposal(x,best,closed_th,closed_tx,closed_ty,out_ratio,require_rim)) return 1;
        if (closed_th == th && closed_tx == tx && closed_ty == ty) return 0;
        fprintf(stderr,"  [place] seam closure failed audit for component %d; checking original proposal\n",best);
    } else fprintf(stderr,"  [place] seam closure unavailable for component %d; checking original proposal\n",best);
    return apl_audit_proposal(x,best,th,tx,ty,out_ratio,require_rim);
}

/* A hop estimate is revisable when subsequently placed material supplies
 * an actual boundary continuation. Reopen the original hop group together
 * with its already attached source pieces. Holding those descendants fixed
 * would pin the very hop being corrected. Preserve every passing source
 * connection and never reopen a bundle connected to the primary root.
 * Neither failure nor ambiguity changes the incumbent. */
static int apl_reconsider_hop(AplCtx *x, int32_t best, double deadline)
{
    AsmRun *run = x->run;
    if (best == x->primary || x->comp_placed[best] != 1) return 0;
    size_t n = 0;
    for (size_t i = 0; i < run->n_charts; i++) {
        const AsmChart *c = &run->charts[i];
        if (!AsmChart_in_layout(c) || x->orig_comp[i] != best) continue;
        if (c->component != x->primary || c->placement_state != ASM_PLACE_LAYER) return 0;
        n++;
    }
    if (!n) return 0;
    Arena_Mark mark = Arena_save(run->arena);
    size_t ng = 0;
    for (size_t i = 0; i < run->n_charts; i++) if (x->orig_comp[i] >= 0 && (size_t)x->orig_comp[i] >= ng) ng = (size_t)x->orig_comp[i]+1;
    int32_t *first = ARENA_ALLOC(run->arena,ng*sizeof *first);
    for (size_t i = 0; i < ng; i++) first[i] = -1;
    UnionFind bundle = UF_new(run->arena,(int32_t)run->n_charts);
    for (size_t i = 0; i < run->n_charts; i++) {
        int32_t g = x->orig_comp[i];
        if (g < 0 || !AsmChart_in_layout(&run->charts[i])) continue;
        if (first[g] < 0) first[g] = (int32_t)i; else uf_union(&bundle,first[g],(int32_t)i);
    }
    for (size_t i = 0; i < run->n_rels; i++) {
        const AsmRelation *r = &run->rels[i];
        const AsmChart *a = &run->charts[r->a], *b = &run->charts[r->b];
        if (!AsmChart_in_layout(a) || !AsmChart_in_layout(b) ||
            a->component != x->primary || b->component != x->primary ||
            a->placement_state == ASM_PLACE_NONE || b->placement_state == ASM_PLACE_NONE) continue;
        double rms; size_t support;
        if (AsmContinuity_measure(run,r,&rms,&support)) uf_union(&bundle,r->a,r->b);
    }
    uint8_t *excluded = ARENA_CALLOC(run->arena,run->n_charts,1);
    uint8_t *was_hop = ARENA_CALLOC(run->arena,ng,1);
    int32_t root = uf_find(&bundle,first[best]);
    double bundle_area = 0; n = 0;
    for (size_t i = 0; i < run->n_charts; i++) {
        const AsmChart *c = &run->charts[i];
        if (!AsmChart_in_layout(c) || uf_find(&bundle,(int32_t)i) != root) continue;
        if (c->component != x->primary || c->placement_state == ASM_PLACE_NONE || c->placement_state == ASM_PLACE_ROOT) {
            Arena_restore(run->arena,mark); return 0;
        }
        excluded[i] = 1; n++; bundle_area += c->area3d;
        if (c->placement_state == ASM_PLACE_LAYER) was_hop[x->orig_comp[i]] = 1;
    }
    memset(x->seam_selected,0,run->n_rels);
    size_t broken = 0;
    for (size_t i = 0; i < run->n_rels; i++) {
        const AsmRelation *r = &run->rels[i];
        if (!(r->continuity & ASM_CONT_SOURCE) || (r->flags & (ASM_REL_DROPPED | ASM_REL_CONTACT))) continue;
        int a = excluded[r->a], b = excluded[r->b];
        if (!a && !b) continue;
        const AsmChart *A = &run->charts[r->a], *B = &run->charts[r->b];
        if (!AsmChart_in_layout(A) || !AsmChart_in_layout(B) ||
            A->placement_state == ASM_PLACE_NONE || B->placement_state == ASM_PLACE_NONE ||
            A->component != x->primary || B->component != x->primary) continue;
        x->seam_selected[i] = 1;
        double rms; size_t support;
        if (a != b && !AsmContinuity_measure(run,r,&rms,&support)) broken++;
    }
    if (!broken) { Arena_restore(run->arena,mark); return 0; }
    AsmChart *saved = ARENA_ALLOC(run->arena,n*sizeof *saved);
    AsmChart *winner = ARENA_ALLOC(run->arena,n*sizeof *winner);
    size_t q = 0;
    for (size_t i = 0; i < run->n_charts; i++) {
        AsmChart *c = &run->charts[i];
        if (excluded[i]) saved[q++] = *c;
    }
    /* A cut label can contain several disconnected pieces. Every moved
     * piece needs a source path to fixed material; one seam cannot certify
     * the rest of a group merely because they share a placement label.
     * All selected paths must subsequently pass the actual geometry gate. */
    UnionFind uf = UF_new(run->arena,(int32_t)run->n_charts);
    for (size_t i = 0; i < run->n_rels; i++) if (x->seam_selected[i]) uf_union(&uf,run->rels[i].a,run->rels[i].b);
    uint8_t *anchored = ARENA_CALLOC(run->arena,run->n_charts,1);
    for (size_t i = 0; i < run->n_charts; i++) if (!excluded[i] && AsmChart_in_layout(&run->charts[i]) &&
        run->charts[i].component == x->primary && run->charts[i].placement_state != ASM_PLACE_NONE)
        anchored[uf_find(&uf,(int32_t)i)] = 1;
    for (size_t i = 0; i < n; i++) if (!anchored[uf_find(&uf,saved[i].id)]) {
        Arena_restore(run->arena,mark); return 0;
    }
    int32_t *original_groups = x->orig_comp;
    int32_t *trial_groups = ARENA_ALLOC(run->arena,run->n_charts*sizeof *trial_groups);
    memcpy(trial_groups,original_groups,run->n_charts*sizeof *trial_groups);
    double original_area = x->comp_area[best];
    x->orig_comp = trial_groups; x->comp_area[best] = bundle_area;
    for (size_t i = 0; i < n; i++) {
        AsmChart *c = &run->charts[saved[i].id]; c->component = best; c->placement_state = ASM_PLACE_NONE;
        trial_groups[c->id] = best;
    }
    x->comp_placed[best] = 0; x->excluded_samples = excluded;
    x->trial = 1; x->stretch_ok = 0; x->attempt++;
    x->na = apl_collect_anchors(x,best);
    double weight[APL_MAX_CLUSTERS_TOTAL], u[APL_MAX_CLUSTERS_TOTAL], v[APL_MAX_CLUSTERS_TOTAL];
    int kind[APL_MAX_CLUSTERS_TOTAL];
    int ncl = apl_cluster_anchors(x,x->na,weight,u,v,kind);
    double *res = ARENA_ALLOC(run->arena,(x->na ? x->na : 1)*sizeof *res);
    Arena_T winner_uv = NULL;
    size_t winner_n = 0; int ambiguous = 0, exhausted = 0;
    fprintf(stderr,"  [place] reconsider hop: component %d, %zu charts, %zu broken source seams\n",best,n,broken);
    for (int cl = 0; cl < ncl; cl++) {
        if (deadline > 0 && ves_clock_sec() > deadline) { exhausted = 1; break; }
        if (kind[cl] != 1) continue;
        x->source_flip = 0;
        for (size_t i = 0; i < x->na; i++) if (x->anchors[i].cluster == cl && x->anchors[i].k == 0) {
            x->source_flip = x->anchors[i].source_flip; break;
        }
        double th, tx, ty; size_t kept;
        double rms = apl_fit_cluster(x,x->na,cl,APL_SEAM_MIN_ANCHORS,x->o->max_rot_refine,&th,&tx,&ty,&kept,res);
        if (rms >= 1e299) continue;
        x->cluster = cl; x->cluster_kind = 1;
        /* seam_selected deliberately includes every incident source seam,
         * not just the cluster which proposes this starting pose. */
        if (!apl_place_and_audit(x,best,th,tx,ty,NULL,2)) continue;
        if (!winner_n) {
            winner_n = x->nproposal;
            apl_copy_candidate(winner,x->proposal,winner_n,run,&winner_uv);
        } else {
            if (winner_n != x->nproposal) ambiguous = 1;
            else for (size_t k = 0; k < winner_n; k++)
                if (!apl_same_pose(&winner[k],&x->proposal[k],8)) { ambiguous = 1; break; }
        }
    }
    if (deadline > 0 && ves_clock_sec() > deadline) exhausted = 1;
    for (size_t i = 0; i < n; i++) run->charts[saved[i].id] = saved[i];
    x->orig_comp = original_groups; x->comp_area[best] = original_area;
    x->comp_placed[best] = 1; x->excluded_samples = NULL; x->source_flip = 0;
    int accepted = winner_n && !ambiguous && !exhausted;
    if (accepted) {
        for (size_t i = 0; i < winner_n; i++) {
            winner[i].placement_group = run->charts[winner[i].id].placement_group;
            apl_commit_chart(x,&winner[i],run->charts[winner[i].id].uv);
            /* The cached winding-direction curvature was measured in the
             * old pose, which may have rotated even with unchanged UV. */
            if (x->curv_ok) x->curv_ok[winner[i].id] = 0;
        }
        AsmContinuity_confirm(run);
        for (size_t g = 0; g < ng; g++) if (was_hop[g]) {
            if (x->diag) x->diag[g].placed_by = 1;
            if (x->st) { x->st->seam_placed++; if (x->st->hop_placed) x->st->hop_placed--; }
        }
        fprintf(stderr,"  [place] source recovery: component %d closes %zu previously broken seams, %zu charts, %zu witnesses at %.3f vox RMS\n",
                best,broken,winner_n,winner[0].placement_support,winner[0].placement_rms);
    } else fprintf(stderr,"  [place] reconsidered hop kept unchanged: component %d%s\n",best,ambiguous ? " (ambiguous)" : "");
    if (exhausted && x->st) x->st->budget_exhausted = 1;
    if (winner_uv) Arena_dispose(&winner_uv);
    Arena_restore(run->arena,mark);
    return accepted;
}

static void apl_rebuild_grids(AplCtx *x, Arena_T *arena)
{
    if (!x->copts || !arena || !*arena) return;
    AsmConflictGrid_dispose(x->cgrid); Arena_dispose(arena); *arena = Arena_new();
    x->cgrid = AsmConflictGrid_new(*arena,x->copts,x->run);
    x->rgrid = aprg_new(*arena,24);
    x->strip_umin = 1e300; x->strip_umax = -1e300;
    for (size_t i = 0; i < x->run->n_charts; i++) {
        const AsmChart *c = &x->run->charts[i];
        if (!AsmChart_in_layout(c) || c->component != x->primary || c->placement_state == ASM_PLACE_NONE) continue;
        AsmConflictGrid_insert(x->cgrid,c); aprg_insert_chart(x->rgrid,c);
        x->strip_umin = fmin(x->strip_umin,c->pose_x); x->strip_umax = fmax(x->strip_umax,c->pose_x);
    }
}

static void apl_profile_line(const AsmPlaceStats *st, size_t sweeps, double total, const char *prefix)
{
    double gates = st->t_clusterloop - st->t_fit - st->t_audit;
    double known = st->t_measure + st->t_growth + st->t_setup + st->t_count + st->t_mirror + st->t_stretch + st->t_collect + st->t_cluster + st->t_clusterloop + st->t_cut;
    fprintf(stderr, "%s: %zu sweeps in %.1f s (%.1f/s) | measure turns %.1f s (%zu re-measurements) growth %.1f setup %.1f candidates %.1f mirror-fit %.1f stretch %.1f (quadric %zu fitted, %zu cached) collect %.1f cluster %.1f fit %.1f gates %.1f audit %.1f (probe %.1f rim %.1f insert %.1f) cut %.1f other %.1f\n",
            prefix, sweeps, total, total > 0.0 ? (double)sweeps / total : 0.0, st->t_measure, st->remeasures, st->t_growth, st->t_setup, st->t_count, st->t_mirror, st->t_stretch,
            st->quadric_calls, st->quadric_cached, st->t_collect, st->t_cluster, st->t_fit, gates < 0.0 ? 0.0 : gates, st->t_audit, st->t_probe, st->t_rim, st->t_insert, st->t_cut, total - known);
}

/* A direct boundary continuation must get its attempt before a hop-only arc
 * can occupy its destination. Ring order and area still break ties within
 * each evidence class; failed candidates wait for new evidence as before. */
static int32_t apl_next_component(int32_t ncomp, const uint8_t *placed,
                                  const int32_t *anchors, const int32_t *seams,
                                  const int32_t *tried, const int32_t *ring,
                                  const double *area, int min_anchors)
{
    int32_t best = -1;
    for (int32_t c = 0; c < ncomp; c++) {
        int has_seam = APL_SEAM_FIRST && seams[c] >= APL_SEAM_MIN_ANCHORS;
        if (placed[c] == 1 || (anchors[c] < min_anchors && !has_seam)) continue;
        if (placed[c] == 2 && anchors[c] <= tried[c]) continue;
        if (best < 0) { best = c; continue; }
        int best_seam = APL_SEAM_FIRST && seams[best] >= APL_SEAM_MIN_ANCHORS;
        if (APL_SEAM_FRONTIER_FIRST && (APL_SEAM_FRONTIER_FIRST == 1 || ring[c] == ring[best]) && has_seam != best_seam) {
            if (has_seam) best = c;
            continue;
        }
        if (ring[c] < ring[best] || (ring[c] == ring[best] && area[c] > area[best])) best = c;
    }
    return best;
}

static void apl_record_poses(const AsmRun *run, const AsmPlaceOpts *o, const char *stage, int append)
{
    if (!o->diag_dir) return;
    char path[2048]; snprintf(path, sizeof path, "%s/stage5_poses.csv", o->diag_dir);
    FILE *fp = fopen(path, append ? "ab" : "wb"); if (!fp) return;
    if (!append) fputs("stage,chart,component,u,v,theta,mirror\n", fp);
    for (size_t i = 0; i < run->n_charts; i++) {
        const AsmChart *c = &run->charts[i]; if (!AsmChart_in_layout(c)) continue;
        fprintf(fp, "%s,%zu,%d,%.9g,%.9g,%.12g,%d\n", stage, i, c->component,
                c->pose_x, c->pose_y, c->pose_theta, (c->flags & ASM_CHART_MIRROR) != 0);
    }
    fclose(fp);
    if (!run->continuity_ready) return;
    snprintf(path, sizeof path, "%s/stage5_source_residuals.csv", o->diag_dir);
    int first = strcmp(stage, "registered_rigid") == 0;
    fp = fopen(path, first ? "wb" : "ab"); if (!fp) return;
    if (first) fputs("stage,chart_a,chart_b,component_a,component_b,flags,source,passes,pairs,rms\n", fp);
    for (size_t i = 0; i < run->n_rels; i++) {
        const AsmRelation *r = &run->rels[i];
        if (!AsmChart_in_layout(&run->charts[r->a]) || !AsmChart_in_layout(&run->charts[r->b])) continue;
        double rms; size_t support; int passes = AsmContinuity_measure(run, r, &rms, &support);
        fprintf(fp, "%s,%d,%d,%d,%d,%u,%d,%d,%zu,%.9g\n", stage, r->a, r->b,
                run->charts[r->a].component, run->charts[r->b].component, r->flags,
                (r->continuity & ASM_CONT_SOURCE) != 0, passes, support, rms);
    }
    fclose(fp);
}

#include "asm_place_strict.inc"

/* A pose which fails admission cannot make a sound pose ambiguous. Apply the
 * complete original-face gates before comparing resumed placement candidates,
 * including every incident source hypothesis and every fixed sheet face. */
static int apl_retry_candidate_qualified(const AplCtx *x,int32_t component)
{
    Arena_T arena=Arena_new();AsmRun trial=*x->run;trial.arena=arena;
    trial.charts=ARENA_ALLOC(arena,trial.n_charts*sizeof(AsmChart));
    memcpy(trial.charts,x->run->charts,trial.n_charts*sizeof(AsmChart));
    for(size_t k=0;k<x->nproposal;k++)trial.charts[x->proposal[k].id]=x->proposal[k];
    AsmField *f=AsmField_new(arena,&trial);int okay=f!=NULL;
    size_t bad_charts=0,failed_seams=0;AsmContactStats contacts={0};
    if(f){
        uint8_t *active=ARENA_CALLOC(arena,f->nf,1);
        for(size_t c=0;c<trial.n_charts;c++)if(x->orig_comp[c]==component){
            if(f->face_offset[c]==SIZE_MAX){bad_charts++;continue;}
            double area=0,good=0;int valid=1;
            for(size_t k=0;k<trial.charts[c].nf;k++){
                size_t t=f->face_offset[c]+k;active[t]=1;double j[4],lo,hi,det;
                AsmField_face(f,t,f->uv,j,&lo,&hi,&det);
                if(!UvGuard_positive(f->uv,f->faces+3*t) || !AsmMetric_within25(lo,hi))valid=0;
                area+=f->metric[t].area;if(AsmMetric_within10(lo,hi))good+=f->metric[t].area;
            }
            bad_charts+=!valid || !(area>0) || good<.95*area-32*DBL_EPSILON*area;
        }
        for(size_t i=0;i<f->n_seams;i++){
            const AsmFieldSeam *s=f->seams+i;const AsmRelation *r=trial.rels+s->relation;
            if((r->flags&ASM_REL_CONTACT) || !(r->continuity&(ASM_CONT_SOURCE|ASM_CONT_CUT)) ||
                (x->orig_comp[s->a]!=component && x->orig_comp[s->b]!=component))continue;
            AsmSeamMeasure m;AsmField_seam(f,i,f->uv,&m);failed_seams+=!m.pass;
        }
        AsmContacts_active_faces(f->contacts,active);
        okay=!AsmContacts_measure(f->contacts,f->uv,NULL,NULL,NULL,&contacts) &&
             !bad_charts && !failed_seams && !contacts.pairs;
    }
    fprintf(stderr,"  [place retry] original field: %zu bad charts, %zu failed source seams, %zu affected overlap pairs; %s\n",
            bad_charts,failed_seams,contacts.pairs,okay?"qualified":"refused");
    Arena_dispose(&arena);return okay;
}

/* One stage5_candidates.csv row per candidate, accepted or refused, with the gate that
 * judged it and the reason it failed.  The consensus columns are written only when the
 * audit reached them; a candidate refused on its fit carries zeros there. */
static void apl_write_candidate(FILE *fp, const AplCtx *x, int32_t best, double area, int cl, int kind,
                                int accepted, const char *reason, size_t kept, double rms, double gate,
                                size_t k1, size_t k2, size_t est, double mean_turn, double turn_rel, int measured)
{
    if (!fp) return;
    int chart = -1, mirror = 0;
    double u = 0.0, v = 0.0, srms = 0.0;
    size_t nprop = 0, support = 0;
    if (accepted && x->nproposal) {
        const AsmChart *c = &x->proposal[0];
        double cu, cv;
        apl_uv_centroid(c, &cu, &cv); apl_pose(c, cu, cv, &u, &v);
        chart = c->id; mirror = (c->flags & ASM_CHART_MIRROR) != 0;
        nprop = x->nproposal; support = c->placement_support; srms = c->placement_rms;
    }
    fprintf(fp, "%d,%d,%d,%d,%d,%s,%zu,%zu,%.6g,%zu,%zu,%zu,%.9g,%.9g,%d,%.6f,%.6f,%d,%zu,%.6g,%.6g,%zu,%zu,%zu,%.6g,%.6g,%.6e\n",
            x->attempt, best, cl, kind, accepted, reason, nprop, support, srms,
            measured ? x->consensus.observations : (size_t)0,
            measured ? x->consensus.inliers : (size_t)0,
            measured ? x->consensus.model_inliers : (size_t)0,
            measured ? x->consensus.clipped_loss : 0.0,
            measured ? x->consensus.model_loss : 0.0,
            chart, u, v, mirror, kept, rms, gate, k1, k2, est, mean_turn, turn_rel, area);
}

/* The part of placement that runs before any placement search: source
 * continuity, rigid registration and UV relaxation of each clean2 component,
 * and the continuity components that mark measured seams VERIFIED. */
static int apl_preamble(AsmRun *run, const AsmPlaceOpts *o)
{
    apl_record_poses(run, o, "clean2", 0);
    if (AsmContinuity_prepare(run) != 0) return -1;
    Arena_T registration_arena = Arena_new();
    AsmChart *registration_entry = ARENA_ALLOC(registration_arena,(run->n_charts ? run->n_charts : 1)*sizeof *registration_entry);
    memcpy(registration_entry,run->charts,run->n_charts*sizeof *registration_entry);
    AsmContinuity_register_components(run,apl_registration_audit);
    apl_record_poses(run, o, "registered_rigid", 1);
    AsmContinuity_relax_uv(run, apl_uv_audit, o->stress_band, o->diag_dir);
    AsmContinuity_register_failed(run,registration_entry,apl_registration_audit);
    Arena_dispose(&registration_arena);
    apl_record_poses(run, o, "registered", 1);
    AsmContinuity_components(run);
    return 0;
}

static int apl_run(AsmRun *run, const AsmPlaceOpts *o, const AsmConflictOpts *copts, AsmPlaceStats *st,
                   const uint8_t *retained, int diagnostic_policy, double diagnostic_seconds)
{
    memset(st, 0, sizeof *st);
    if (!retained && apl_preamble(run, o) != 0) return -1;
    double t0 = ves_clock_sec();
    Arena_T arena = run->arena;
    Arena_Mark mark = Arena_save(arena);
    size_t nc = run->n_charts;
    for (size_t i = 0; i < nc; i++) {
        if (retained && retained[i]) continue;
        run->charts[i].placement_group = run->charts[i].component;
        run->charts[i].placement_state = ASM_PLACE_NONE;
        run->charts[i].placement_parent = -1;
        run->charts[i].placement_support = 0;
        run->charts[i].placement_rms = 0;
    }
    int32_t ncomp = 0;
    for (size_t i = 0; i < nc; i++) if (AsmChart_in_layout(&run->charts[i]) && run->charts[i].component + 1 > ncomp) ncomp = run->charts[i].component + 1;
    if (ncomp == 0) { Arena_restore(arena, mark); return 0; }
    size_t ncomp_cap = (size_t)ncomp + APL_MAX_CUTS + 8;   /* cuts append components */
    double *comp_area = ARENA_CALLOC(arena, ncomp_cap, sizeof(double));
    for (size_t i = 0; i < nc; i++) {
        const AsmChart *c = &run->charts[i];
        if (!AsmChart_in_layout(c) || c->component < 0) continue;
        comp_area[c->component] += c->area3d;
        st->area_in += c->area3d;
    }
    int32_t primary = 0;
    if (!retained) for (int32_t c = 1; c < ncomp; c++) if (comp_area[c] > comp_area[primary]) primary = c;
    st->components_in = (size_t)ncomp;

    uint8_t *comp_placed = ARENA_CALLOC(arena, ncomp_cap, sizeof(uint8_t));
    AplTurn *turn = ARENA_ALLOC(arena, nc * 2 * sizeof(AplTurn));
    int8_t *out_side = ARENA_CALLOC(arena, nc ? nc : 1, sizeof(int8_t));
    AplCompDiag *diag = ARENA_CALLOC(arena, ncomp_cap, sizeof(AplCompDiag));
    FILE *fanch = NULL, *fcandidate = NULL, *fevidence = NULL, *frim = NULL, *fcut = NULL;
    if (o->diag_dir) {
        char path[2048];
        snprintf(path, sizeof path, "%s/stage5_anchors.csv", o->diag_dir);
        fanch = fopen(path, "wb");
        /* turn_source is the APL_TSRC_* enum (0 seam, 1 measured, 2 measured+growth k>=2,
         * 3 opposite+/-growth, 4 radius, 5 layer model, 6 stretch, 7 bare t0).  Before
         * 2026-09-19 this column was a 2-way `estimated ? 2 : 1`, so older CSVs do NOT
         * carry these codes -- the growth-carrying families (2, 3, 5) were indistinguishable. */
        if (fanch) fputs("attempt,component,area,cluster,chart_p,chart_c,side,k,turn_source,reverse,Tx,Ty,sx,sy,tx,ty,residual,kept\n", fanch);
        snprintf(path, sizeof path, "%s/stage5_candidates.csv", o->diag_dir);
        fcandidate = fopen(path, "wb");
        /* Every candidate, accepted or not, with the gate that judged it and the reason it
         * failed.  Until 2026-09-17 only the survivors were written, so the deferrals that
         * dominate the big rungs (the consensus minority, the ambiguity) existed in stderr
         * alone and could not be mined. */
        if (fcandidate) fputs("attempt,component,cluster,kind,accepted,reason,charts,selected_support,selected_rms,"
                              "observations,all_inliers,model_inliers,clipped_loss,model_loss,chart,u_mean,v_mean,mirror,"
                              "anchors,fit_rms,gate,k1,k2,estimated_turns,mean_turn,turn_rel,area\n", fcandidate);
        snprintf(path, sizeof path, "%s/stage5_evidence.csv", o->diag_dir);
        fevidence = fopen(path, "wb");
        if (fevidence) fputs("attempt,component,cluster,side,chart,observations,inliers,model_loss\n", fevidence);
        /* Every CUT PROPOSAL and which of its preconditions refused it.  The 21x5x5 proposes
         * 370 and commits none; the summary line says only how many the source preservation
         * vetoed, never what the crossing seams looked like. */
        snprintf(path, sizeof path, "%s/stage5_cut_proposals.csv", o->diag_dir);
        fcut = fopen(path, "wb");
        if (fcut) fputs("attempt,component,outcome,c1,c2,w1,w2,rms1,rms2,u1,u2,sep,min_sep,t0,"
                        "n1,n2,m1,m2,area1,area2,crossing,crossing_sound,sound_min_rms,sound_max_support,first_sound_rel\n", fcut);
        snprintf(path, sizeof path, "%s/stage5_rim_failures.csv", o->diag_dir);
        frim = fopen(path, "wb");
        if (frim) fputs("attempt,component,cluster,uv_trial,chart,vertex,partner,sample,partner_state,u,v,partner_u,partner_v,z,y,x,partner_z,partner_y,partner_x,uv_distance,physical_distance,tolerance\n",frim);
    }

    /* 1. orient every component in the scroll frame: +v = up (world z), and
     * u increasing OUTWARD so the scroll's centre is at the left.  Components
     * with internal layer pairs measure that directly from their outward turn
     * vector; the others follow the scroll's chirality = the sign relation
     * between the bending vote and the turn vote, measured on the components
     * that have both. */
    AplOrient *ori = ARENA_CALLOC(arena, ncomp_cap, sizeof(AplOrient));
    double chirality = 0.0;
    for (int32_t c = 0; c < ncomp; c++) {
        if (comp_area[c] <= 0.0) continue;
        if(retained && c==primary) apl_orient_measure(run,c,&ori[c]);
        else apl_orient_rotate(run, c, &ori[c]);
        /* the frame of reference must be a SHEET: the largest component whose
         * intrinsic v follows world z (the crumpled core wall of a centre crop is
         * the biggest piece and reads a z-gradient of 0.6) */
        if (!retained && ori[c].grad >= APL_SHEET_GRAD && (ori[primary].grad < APL_SHEET_GRAD || comp_area[c] > comp_area[primary])) primary = c;
        apl_component_out_sides(run, c, o->axis, out_side);
        comp_placed[c] = 1;
        apl_measure_turns(run, comp_placed, out_side, NULL, turn, o->min_turn_hits);
        comp_placed[c] = 0;
        apl_orient_votes(run, c, o->axis, turn, &ori[c], out_side);
        if (ori[c].has_turn) chirality += comp_area[c] * ori[c].bend * ori[c].turn_u;
    }
    for (int32_t c = 0; c < ncomp; c++) {
        if (comp_area[c] <= 0.0) continue;
        /* a weak turn vote (a fold, a piece facing itself around a corner)
         * is no evidence: such components take the chirality rule and are
         * re-checked against their anchors at placement */
        if (ori[c].has_turn && fabs(ori[c].turn_u) < 0.6) ori[c].has_turn = 0;
        int mirror = ori[c].has_turn ? (ori[c].turn_u < 0.0) : (chirality != 0.0 && chirality * ori[c].bend < 0.0);
        if(retained && c==primary) mirror=0;
        if (mirror) { apl_mirror_component(run, c); ori[c].mirrored = 1; st->oriented_mirrored++; }
    }
    st->area_primary = comp_area[primary];
    fprintf(stderr, "  [place] primary component %d oriented: rotation %.3f rad, z-gradient %.3f, bend %+.3f, outward turn along u %+.3f%s%s; scroll chirality %+.0f\n",
            primary, ori[primary].theta, ori[primary].grad, ori[primary].bend, ori[primary].turn_u,
            ori[primary].has_turn ? "" : " (unmeasured)", ori[primary].mirrored ? ", mirrored" : "",
            chirality > 0.0 ? 1.0 : chirality < 0.0 ? -1.0 : 0.0);
    {
        size_t n_os = 0;
        for (size_t c = 0; c < nc; c++) n_os += out_side[c] != 0;
        fprintf(stderr, "  [place] outward side known for %zu of %zu charts\n", n_os, nc);
    }

    int32_t *orig_comp = ARENA_ALLOC(arena, (nc ? nc : 1) * sizeof(int32_t));
    for (size_t i = 0; i < nc; i++) orig_comp[i] = run->charts[i].component;
    if (o->diag_dir) {
        char path[2048];
        snprintf(path, sizeof path, "%s/stage5_layers.csv", o->diag_dir);
        FILE *fl = fopen(path, "wb");
        if (fl) {
            fputs("chart_a,chart_b,comp_a,comp_b,side,out_side_a,count,d_median,hit_count,area_a,area_b,k,nsign\n", fl);
            for (size_t p = 0; p < run->n_layers; p++) {
                const AsmLayerPair *lp = &run->layers[p];
                const AsmChart *A = &run->charts[(size_t)lp->a], *B = &run->charts[(size_t)lp->b];
                fprintf(fl, "%d,%d,%d,%d,%d,%d,%d,%.2f,%d,%.4e,%.4e,%d,%d\n", lp->a, lp->b, A->component, B->component, lp->side, out_side[lp->a],
                        lp->count, lp->d_median, lp->hit_count, A->area3d, B->area3d, (int)lp->k, (int)lp->nsign);
            }
            fclose(fl);
        }
    }

    comp_placed[primary] = 1;
    for (size_t i = 0; i < nc; i++) if ((!retained || !retained[i]) && AsmChart_in_layout(&run->charts[i]) && run->charts[i].component == primary)
        run->charts[i].placement_state = ASM_PLACE_ROOT;
    diag[primary].status = APL_PRIMARY;
    st->placed = 1;
    st->area_placed = comp_area[primary];
    double *rms_list = ARENA_ALLOC(arena, ncomp_cap * sizeof(double));
    size_t n_rms = 0;
    /* Scratch is bounded by incident evidence of one candidate, including
     * every alternate hop. A subsequent cut can only reduce that bound. */
    size_t *incident = ARENA_CALLOC(arena, ncomp_cap, sizeof(size_t));
    size_t hop_count = (size_t)(o->max_layers > 0 ? o->max_layers : 1);
    for (size_t r = 0; r < run->n_rels; r++) {
        const AsmRelation *R = &run->rels[r];
        int32_t a = run->charts[R->a].component, b = run->charts[R->b].component;
        if (R->corr_count <= 0) continue;
        if (a >= 0) incident[a] += (size_t)R->corr_count;
        if (b >= 0 && b != a) incident[b] += (size_t)R->corr_count;
    }
    for (size_t p = 0; p < run->n_layers; p++) {
        const AsmLayerPair *L = &run->layers[p];
        int32_t a = run->charts[L->a].component, b = run->charts[L->b].component;
        if (L->hit_count <= 0) continue;
        size_t n = (size_t)L->hit_count * hop_count;
        if (a >= 0) incident[a] += n;
        if (b >= 0 && b != a) incident[b] += n;
    }
    size_t hcap = 1;
    for (int32_t c = 0; c < ncomp; c++) if (c != primary && incident[c] > hcap) hcap = incident[c];
    fprintf(stderr, "  [place] candidate scratch: %zu anchors (%.1f MiB)\n", hcap, (double)hcap*(sizeof(AplAnchor)+2*sizeof(double))/(1024.0*1024.0));
    AplCtx x;
    memset(&x, 0, sizeof x);
    x.diagnostic_axis_policy = diagnostic_policy > 0 && diagnostic_policy <= 2 ? diagnostic_policy : 0;
    x.diagnostic_first_step = diagnostic_policy == 3 || diagnostic_policy == 4;
    if (diagnostic_policy == 4) x.strict = aps_new(arena, run, o->diag_dir);
    x.st = st;
    x.seam_selected = ARENA_CALLOC(arena, run->n_rels ? run->n_rels : 1, 1);
    x.trial = 1; x.proposal = ARENA_ALLOC(arena, (nc ? nc : 1)*sizeof(AsmChart));
    AsmChart *winner = ARENA_ALLOC(arena, (nc ? nc : 1)*sizeof(AsmChart));
    uint8_t *winner_seams = x.strict ? ARENA_CALLOC(arena, run->n_rels ? run->n_rels : 1, 1) : NULL;
    Arena_T winner_uv = NULL;
    x.curv_hl = ARENA_CALLOC(arena, nc ? nc : 1, sizeof(double));
    x.curv_ok = ARENA_CALLOC(arena, nc ? nc : 1, sizeof(uint8_t));
    double *lidx = ARENA_ALLOC(arena, (nc ? nc : 1) * sizeof(double));
    x.tfield = ARENA_CALLOC(arena, nc ? nc : 1, sizeof(double));
    double *kover = ARENA_CALLOC(arena, run->n_layers ? run->n_layers : 1, sizeof(double));
    int32_t coord_anchor = -1;
    int32_t *comp_ring = ARENA_ALLOC(arena, ncomp_cap * sizeof(int32_t));   /* per component: nearest layer to the primary, 1000000 = unknown */
    for (size_t i = 0; i < nc; i++) lidx[i] = APL_LAYER_UNKNOWN;
    for (int32_t c = 0; c < ncomp; c++) comp_ring[c] = 1000000;
    x.lidx = lidx;
    Arena_T garena = NULL;   /* the persistent grids live in their OWN arena: the placement's
                              * scoped save/restore on the run arena would free their growth */
    if (copts) {
        double tg = ves_clock_sec();
        garena = Arena_new();
        x.cgrid = AsmConflictGrid_new(garena, copts, run);
        x.rgrid = aprg_new(garena, 24.0);
        for (size_t i = 0; i < nc; i++) {
            const AsmChart *c = &run->charts[i];
            if (!AsmChart_in_layout(c) || c->component != primary) continue;
            AsmConflictGrid_insert(x.cgrid, c); aprg_insert_chart(x.rgrid, c);
        }
        fprintf(stderr, "  [place] persistent audit grids built from the primary: %zu cells, %.1f s\n", AsmConflictGrid_cells(x.cgrid), ves_clock_sec() - tg);
    }
    {
        /* the strip's u extent (chart origins) and the physical cap on its span */
        float lo[3] = { 1e30f, 1e30f, 1e30f }, hi[3] = { -1e30f, -1e30f, -1e30f };
        x.strip_umin = 1e300; x.strip_umax = -1e300;
        for (size_t i = 0; i < nc; i++) {
            const AsmChart *ch = &run->charts[i];
            if (!AsmChart_in_layout(ch)) continue;
            for (int d = 0; d < 3; d++) { if (ch->bbox_lo[d] < lo[d]) lo[d] = ch->bbox_lo[d]; if (ch->bbox_hi[d] > hi[d]) hi[d] = ch->bbox_hi[d]; }
            if (ch->component != primary) continue;
            if (ch->pose_x < x.strip_umin) x.strip_umin = ch->pose_x;
            if (ch->pose_x > x.strip_umax) x.strip_umax = ch->pose_x;
        }
        double box_diag = hypot((double)hi[1] - lo[1], (double)hi[2] - lo[2]);
        double dl = run->layer_d_global > 0.0 ? run->layer_d_global : 10.0;
        x.span_cap = (box_diag / dl + 2.0) * APL_TURN_MAX;
        fprintf(stderr, "  [place] strip span cap %.0f vox (box diagonal %.0f / layer %.1f = %.0f wraps x %.0f)\n", x.span_cap, box_diag, dl, box_diag / dl, APL_TURN_MAX);
    }
    x.run = run; x.o = o; x.copts = copts; x.comp_placed = comp_placed; x.turn = turn; x.out_side = out_side; x.orig_comp = orig_comp;
    {
        /* the axis, per chart: s = the table's arc length at the centroid's foot, r = the distance
         * from the polyline (world z and no radius without a table) */
        x.s_axis = ARENA_ALLOC(arena, (nc ? nc : 1) * sizeof(double));
        x.r_axis = ARENA_ALLOC(arena, (nc ? nc : 1) * sizeof(double));
        x.have_radius = o->axis != NULL;
        x.cu = ARENA_CALLOC(arena, nc ? nc : 1, sizeof(double));
        x.cv = ARENA_CALLOC(arena, nc ? nc : 1, sizeof(double));
        for (size_t i = 0; i < nc; i++) {
            const AsmChart *c = &run->charts[i];
            if (c->uv != NULL && c->nv > 0) {
                double su = 0.0, sv = 0.0;
                for (size_t v = 0; v < c->nv; v++) { su += c->uv[v*2]; sv += c->uv[v*2+1]; }
                x.cu[i] = su / (double)c->nv; x.cv[i] = sv / (double)c->nv;
            }
            if (o->axis) {
                double cen[3] = { c->centroid[0], c->centroid[1], c->centroid[2] };
                AsmAxis_project(o->axis, cen, &x.s_axis[i], &x.r_axis[i], NULL, NULL);
            } else { x.s_axis[i] = c->centroid[0]; x.r_axis[i] = 0.0; }
        }
        x.primary = primary;
        x.anchors_radius = 0; x.vm_ok = 0; x.vm_pin = 0; x.uout_ok = 0; x.uout_x = 1.0; x.uout_y = 0.0;
        if (APL_V_FROM_AXIS) {
            apl_v_model(&x);
            x.vm_pin = x.vm_ok && x.vm_b >= APL_V_MODEL_B_MIN && x.vm_b <= APL_V_MODEL_B_MAX;
            x.vm_p90 = x.vm_ok ? apl_v_model_p90(&x) : 0.0;
            if (x.have_radius && x.vm_ok && !x.vm_pin) {
                x.have_radius = 0;   /* the same check disarms the radius turn model: its r comes from the same table */
                fprintf(stderr, "  [place] turn length from the radius: DISARMED with the pinning (b %.4f outside [%.2f, %.2f]: the table is not this sheet's axis); the measured turn model places\n", x.vm_b, APL_V_MODEL_B_MIN, APL_V_MODEL_B_MAX);
            }
            char pin_state[128];
            if (x.vm_pin) snprintf(pin_state, sizeof pin_state, "ARMED");
            else snprintf(pin_state, sizeof pin_state, "OFF (b outside [%.2f, %.2f]: the table is not this sheet's axis)", APL_V_MODEL_B_MIN, APL_V_MODEL_B_MAX);
            if (x.vm_ok) fprintf(stderr, "  [place] v from the axis: the primary's v = %.1f + %.4f s + %.4f u (s = axial coordinate; s-u correlation %.2f; the primary's own scatter p90 %.1f vox); tolerance %.0f vox rms, max correction %.0f vox; pinning %s\n", x.vm_a, x.vm_b, x.vm_c, x.vm_corr, x.vm_p90, fmax(APL_V_AXIS_TOL, 1.5 * x.vm_p90), fmax(APL_V_AXIS_MAX_SHIFT, 2.0 * x.vm_p90), pin_state);
            else fprintf(stderr, "  [place] v from the axis: no model (primary too small), placement keeps the anchors' v\n");
        }
        apl_uout(&x);
        if (x.have_radius) {
            double rmin = 1e300, rmax = -1e300;
            for (size_t i = 0; i < nc; i++) if (AsmChart_in_layout(&run->charts[i])) { if (x.r_axis[i] < rmin) rmin = x.r_axis[i]; if (x.r_axis[i] > rmax) rmax = x.r_axis[i]; }
            fprintf(stderr, "  [place] turn length from the radius: axis %s (%zu rows), radii %.0f..%.0f vox (turns %.0f..%.0f), outward direction %s(%.2f, %.2f)\n",
                    o->axis->is_line ? "line" : "table", o->axis->n_rows, rmin, rmax, 2.0 * APL_PI * rmin, 2.0 * APL_PI * rmax, x.uout_ok ? "" : "UNKNOWN ", x.uout_x, x.uout_y);
        }
    }
    x.anchors = ARENA_ALLOC(arena, hcap * sizeof(AplAnchor));
    x.wts = ARENA_ALLOC(arena, hcap * sizeof(double));
    x.hcap = hcap; x.primary = primary; x.comp_area = comp_area; x.diag = diag; x.fanch = fanch; x.fevidence = fevidence; x.frim = frim; x.fcut = fcut;
    double *res = ARENA_ALLOC(arena, hcap * sizeof(double));
    int32_t *anchor_count = ARENA_ALLOC(arena, ncomp_cap * sizeof(int32_t));
    int32_t *tried_at = ARENA_CALLOC(arena, ncomp_cap, sizeof(int32_t));
    size_t *tried_placed = ARENA_CALLOC(arena, ncomp_cap, sizeof(size_t));
    size_t *reconsidered_at = ARENA_CALLOC(arena,ncomp_cap,sizeof(size_t));
    uint8_t *rim_deferred = ARENA_CALLOC(arena, ncomp_cap, sizeof(uint8_t));
    int32_t *seam_count = ARENA_CALLOC(arena, ncomp_cap, sizeof(int32_t));   /* seam-anchor correspondences per unplaced component */
    uint8_t *has_seams = ARENA_CALLOC(arena, ncomp_cap, sizeof(uint8_t));    /* APL_DEFER_ONLY_WITH_SEAMS: the component holds a seam to another component */
    for (size_t r = 0; r < run->n_rels; r++) {
        const AsmRelation *R = &run->rels[r];
        if (R->flags & ASM_REL_CONTACT) continue;
        const AsmChart *A = &run->charts[(size_t)R->a], *B = &run->charts[(size_t)R->b];
        if (!AsmChart_in_layout(A) || !AsmChart_in_layout(B) || A->component < 0 || B->component < 0 || A->component == B->component) continue;
        if (A->component < (int32_t)ncomp_cap) has_seams[A->component] = 1;
        if (B->component < (int32_t)ncomp_cap) has_seams[B->component] = 1;
    }
    /* chart adjacency over the seam relations, for the placement-gated cut (dropped relations are skipped at walk time) */
    int32_t *rel_off = ARENA_CALLOC(arena, nc + 1, sizeof(int32_t));
    int32_t *rel_adj = ARENA_ALLOC(arena, (run->n_rels ? run->n_rels : 1) * 2 * sizeof(int32_t));
    {
        for (size_t r = 0; r < run->n_rels; r++) { rel_off[run->rels[r].a + 1]++; rel_off[run->rels[r].b + 1]++; }
        for (size_t i = 0; i < nc; i++) rel_off[i + 1] += rel_off[i];
        int32_t *fill = ARENA_CALLOC(arena, nc ? nc : 1, sizeof(int32_t));
        for (size_t r = 0; r < run->n_rels; r++) { rel_adj[rel_off[run->rels[r].a] + fill[run->rels[r].a]++] = (int32_t)r; rel_adj[rel_off[run->rels[r].b] + fill[run->rels[r].b]++] = (int32_t)r; }
    }
    size_t n_cuts = 0;
    size_t n_growth = 0;
    x.growth = 2.0 * APL_PI * run->layer_d_global;

    /* one candidate per sweep (retries and cuts included): the cap follows the component count
     * and a wall-clock budget, not a constant -- at 512 the 21x5x5 / 10x10x10 / 4x21x21 left
     * 645 / 367 / 1,234 components with enough anchors unattempted (12 / 5 / 13% of the area;
     * reviewer, 2026-09-09), which the ledger then called 'never' */
    int max_sweeps = APL_SWEEPS_LEGACY_CAP ? 512 : 512 + 16 * ncomp;
    int placed_dirty = 1;   /* APL_TURN_REMEASURE: the turns' inputs changed (sweep 0, or a placement) */
    size_t t0_turns = 0;    /* how many measured turns the incumbent t0 rests on (APL_T0_REMEASURE) */
    int height_dirty = 0;   /* a changed global model is evidence even without new direct anchors */
    double t_sweep0 = ves_clock_sec();
    const double budget_sec = diagnostic_policy >= 0 && diagnostic_seconds > 0
        ? fmin(apl_budget_sec(run), diagnostic_seconds) : apl_budget_sec(run);
    if (x.strict) x.strict->deadline = t_sweep0 + budget_sec;
    size_t source_recoveries = 0;
    int last_call = 0;      /* APL_DEFER_LAST_CALL: the second pass over the deferred components */
    for (int phase = 0; phase < (APL_DEFER_LAST_CALL ? 2 : 1); phase++) {
    if (phase == 1) {
        size_t nd = 0;
        for (int32_t c = 0; c < ncomp; c++) if (comp_placed[c] != 1 && diag[c].status == APL_DEFERRED) { comp_placed[c] = 0; tried_at[c] = 0; nd++; }
        if (nd == 0) break;
        last_call = 1;
        fprintf(stderr, "  [place] last call: %zu deferred components retried with the strict hop gate (no turn fraction) and the deferral off\n", nd);
    }
    for (int sweep = 0; sweep < max_sweeps; sweep++) {
        if (ves_clock_sec() - t_sweep0 > budget_sec) { st->budget_exhausted = 1; fprintf(stderr, "  [place] wall-clock budget %.0f s reached after %d sweeps\n", budget_sec, sweep); break; }
        st->sweeps++;
        double tq = ves_clock_sec();
        if (APL_PROFILE_EVERY > 0 && sweep > 0 && (sweep % APL_PROFILE_EVERY) == 0) apl_profile_line(st, (size_t)sweep, tq - t_sweep0, "  [place] profile");
        int remeasured_now = 0;
        if (!APL_TURN_REMEASURE || placed_dirty) {
            apl_measure_turns(run, comp_placed, out_side, lidx, turn, o->min_turn_hits);
            APL_TICK(t_measure);
            apl_fit_turn_radius(&x);   /* the section's perimeter excess and the slope of |T| on r */
            n_growth = apl_growth_for_hops(&x, run, turn, out_side, &x.growth);
            APL_TICK(t_growth);
            if (placed_dirty && apl_v_model_update(&x)) {
                height_dirty = 1;
                fprintf(stderr,"  [place] recalibrated height from placed material: b %.6f c %.6f, unchanged uncertainty band %.1f vox\n",
                        x.vm_b,x.vm_c,apl_axis_band(&x));
            }
            placed_dirty = 0; st->remeasures++; remeasured_now = 1;
        }
        /* t0 was measured ONCE, at sweep 0, from whatever turns existed then.  On a pile that
         * starts with none -- the 21x5x5 logs "t0 curvature over 0 turns" -- that guess is
         * frozen for all 4,148 sweeps, and t0_measured 0 also disables the layer model in
         * apl_collect_anchors.  Follow the evidence: re-measure whenever the placed material
         * changed, and accept only a MEASURED estimate resting on more turns than the
         * incumbent, so the model can sharpen but never soften.  The LAYER COORDINATE stays
         * solved once; re-solving it would move the ring order under the sweep loop. */
        if (APL_T0_REMEASURE && sweep > 0 && remeasured_now) {
            size_t ntl2 = 0; int meas2 = 0;
            double t0b = apl_turn_t0(run, turn, out_side, orig_comp, primary, &meas2, &ntl2);
            if (meas2 && ntl2 > t0_turns) {
                fprintf(stderr, "  [place] turn model re-measured at sweep %d: T(w) = %.0f + w x %.1f vox over %zu turns (was %.0f, %s over %zu)\n",
                        sweep, t0b, x.growth, ntl2, x.t0, x.t0_measured ? "measured" : "curvature", t0_turns);
                x.t0 = t0b; x.t0_measured = 1; t0_turns = ntl2;
            }
        }
        size_t n_borrowed = APL_BORROW_TURNS ? apl_borrow_turns(run, comp_placed, out_side, lidx, turn) : 0;
        if (sweep == 0) fprintf(stderr, "  [place] borrowed turns: %zu chart sides took the nearest measured turn\n", n_borrowed);
        if (sweep == 0) {
            fprintf(stderr, "  [place] turn growth per layer %.1f vox measured on %zu charts (offset-curve rule 2 pi d = %.1f) [growth mode %d]\n",
                    x.growth, n_growth, 2.0 * APL_PI * run->layer_d_global, ASM_PLACE_GROWTH_PER_WRAP);
            fprintf(stderr, "  [place] section perimeter excess: T(r) = 2 pi r + %.0f vox over %zu charts, free slope %.3f (2 pi = %.3f, ratio %.3f) [intercept mode %d, %s]\n",
                    x.turn_c, x.turn_c_n, x.turn_slope, 2.0 * APL_PI,
                    x.turn_slope > 0.0 ? x.turn_slope / (2.0 * APL_PI) : 0.0,
                    ASM_PLACE_TURN_INTERCEPT,
                    (ASM_PLACE_TURN_INTERCEPT && x.turn_c_n) ? "APPLIED to the radius turn" : "not applied");
            /* the turn length at the primary's wrap: measured outward turns of its charts, else the curvature */
            size_t ntl = 0;
            x.t0 = apl_turn_t0(run, turn, out_side, orig_comp, primary, &x.t0_measured, &ntl);
            t0_turns = x.t0_measured ? ntl : 0;
            fprintf(stderr, "  [place] turn model: T(w) = %.0f + w x %.1f vox (t0 %s over %zu turns)\n", x.t0, x.growth, x.t0_measured ? "measured" : "curvature", ntl);
            /* the layer coordinate needs the turn model: solved here, once */
        {
            int32_t anchor = -1;
            for (size_t i = 0; i < nc; i++) {
                const AsmChart *c = &run->charts[i];
                if (!AsmChart_in_layout(c) || c->component != primary) continue;
                if (anchor < 0 || c->area3d > run->charts[(size_t)anchor].area3d) anchor = (int32_t)i;
            }
            size_t n_seam = 0, n_layer = 0, n_known = 0;
            double p50 = 0.0, p95 = 0.0, lmin = 0.0, lmax = 0.0;
            coord_anchor = anchor;
            if (anchor >= 0) apl_layer_index(run, out_side, turn, x.t0, x.growth, anchor, lidx, NULL, &n_seam, &n_layer, &n_known, &p50, &p95);
            else for (size_t i = 0; i < nc; i++) lidx[i] = APL_LAYER_UNKNOWN;
            for (size_t i = 0; i < nc; i++) {
                if (lidx[i] <= APL_LAYER_UNKNOWN || !AsmChart_in_layout(&run->charts[i])) continue;
                if (lidx[i] < lmin) lmin = lidx[i];
                if (lidx[i] > lmax) lmax = lidx[i];
                int32_t ring = (int32_t)floor(fabs(lidx[i]) + 0.5), oc = orig_comp[i];
                if (ring < comp_ring[oc]) comp_ring[oc] = ring;
            }
            fprintf(stderr, "  [place] layer coordinate: %zu seam edges + %zu layer edges, %zu of %zu charts reached, layers %.1f..%.1f, edge residual p50 %.2f p95 %.2f layers\n",
                    n_seam, n_layer, n_known, nc, lmin, lmax, p50, p95);
            if (APL_COORD_CALIBRATE && coord_anchor >= 0) {
                size_t ncal = apl_kover_build(&x, comp_placed, orig_comp, lidx, kover);
                size_t ns2 = 0, nl2 = 0, nk2 = 0; double q50 = 0.0, q95 = 0.0;
                apl_layer_index(run, out_side, turn, x.t0, x.growth, coord_anchor, lidx, kover, &ns2, &nl2, &nk2, &q50, &q95);
                apl_comp_ring_update(run, lidx, orig_comp, ncomp, comp_ring);
                fprintf(stderr, "  [place] coordinate calibrated on %zu placed pairs (hop counts from the layout): %zu charts reached, edge residual p50 %.2f p95 %.2f layers\n", ncal, nk2, q50, q95);
            }
        }
        if (APL_TURN_KHOP) {
            /* the coordinate now tells each placed pair's hop count: measure again per wrap */
            apl_measure_turns(run, comp_placed, out_side, lidx, turn, o->min_turn_hits);
            apl_fit_turn_radius(&x);
            n_growth = apl_growth_for_hops(&x, run, turn, out_side, &x.growth);
            size_t ntl2 = 0;
            x.t0 = apl_turn_t0(run, turn, out_side, orig_comp, primary, &x.t0_measured, &ntl2);
            fprintf(stderr, "  [place] turn model per wrap: T(w) = %.0f + w x %.1f vox (t0 %s over %zu turns, growth on %zu charts)\n",
                    x.t0, x.growth, x.t0_measured ? "measured" : "curvature", ntl2, n_growth);
        }
        }
        APL_TICK(t_setup);
        /* evidence per unplaced component: hits from or to placed charts with a known turn */
        memset(anchor_count, 0, (size_t)ncomp * sizeof(int32_t));
        memset(seam_count, 0, (size_t)ncomp * sizeof(int32_t));
        if (APL_SEAM_ANCHORS) for (size_t r = 0; r < run->n_rels; r++) {
            const AsmRelation *R = &run->rels[r];
            if (!apl_placement_seam(R, x.have_radius, o->crosswrap_anchor_veto)) continue;
            const AsmChart *A = &run->charts[(size_t)R->a], *B = &run->charts[(size_t)R->b];
            if (!AsmChart_in_layout(A) || !AsmChart_in_layout(B)) continue;
            if (comp_placed[A->component] == 1 && comp_placed[B->component] != 1) { anchor_count[B->component] += R->corr_count; seam_count[B->component] += R->corr_count; }
            else if (comp_placed[A->component] != 1 && comp_placed[B->component] == 1) { anchor_count[A->component] += R->corr_count; seam_count[A->component] += R->corr_count; }
        }
        for (size_t p = 0; p < run->n_layers; p++) {
            const AsmLayerPair *lp = &run->layers[p];
            const AsmChart *A = &run->charts[(size_t)lp->a], *B = &run->charts[(size_t)lp->b];
            if (!AsmChart_in_layout(A) || !AsmChart_in_layout(B)) continue;
            double tx = 0.0, ty = 0.0;
            if (comp_placed[A->component] == 1 && comp_placed[B->component] != 1) {
                if (lidx[lp->a] > APL_LAYER_UNKNOWN && lidx[lp->b] > APL_LAYER_UNKNOWN && fabs(lidx[lp->b] - lidx[lp->a]) >= 0.5) anchor_count[B->component] += lp->hit_count;
                else if (apl_turn_for(turn, lp->a, lp->side, out_side[lp->a], x.growth, &tx, &ty)) anchor_count[B->component] += lp->hit_count;
                else if (turn[(size_t)lp->a * 2].n == 0 && turn[(size_t)lp->a * 2 + 1].n == 0) anchor_count[B->component] += lp->hit_count;   /* stretch-estimable */
            } else if (comp_placed[A->component] != 1 && comp_placed[B->component] == 1) {
                if (lidx[lp->a] > APL_LAYER_UNKNOWN && lidx[lp->b] > APL_LAYER_UNKNOWN && fabs(lidx[lp->b] - lidx[lp->a]) >= 0.5) { anchor_count[A->component] += lp->hit_count; continue; }
                if (turn[(size_t)lp->b * 2].n == 0 && turn[(size_t)lp->b * 2 + 1].n == 0) { anchor_count[A->component] += lp->hit_count; continue; }
                if (out_side[lp->a] == 0 || out_side[lp->b] == 0) continue;
                int side = lp->side == out_side[lp->a] ? -out_side[lp->b] : out_side[lp->b];
                if (apl_turn_for(turn, lp->b, side, out_side[lp->b], x.growth, &tx, &ty)) anchor_count[A->component] += lp->hit_count;
            }
        }
        APL_TICK(t_count);
        int32_t best = apl_next_component(ncomp, comp_placed, anchor_count, seam_count,
                                          tried_at, comp_ring, comp_area, o->min_anchors);
        if (best < 0 && height_dirty) {
            /* Finish the normal frontier first, then retry each failed
             * component under the new model. Anchor counts alone cannot
             * detect this evidence change. No unchanged-model retry loop. */
            size_t retry = 0;
            for (int32_t c = 0; c < ncomp; c++) if (comp_placed[c] == 2) {
                comp_placed[c] = 0; tried_at[c] = 0; retry++;
            }
            height_dirty = 0;
            if (retry) {
                fprintf(stderr,"  [place] height-model retry: %zu failed components after committed calibration changed (all gates retained)\n",retry);
                best = apl_next_component(ncomp, comp_placed, anchor_count, seam_count,
                                           tried_at, comp_ring, comp_area, o->min_anchors);
            }
        }
        if (best < 0 && APL_DEFER_RETRY_RIM && !last_call) {
            /* Anchor counts only track direct graph neighbours. A new chart
             * can add the missing rim support without adding any anchors.
             * Retry those silent attempts once per changed placement grid;
             * all original fit, contradiction and rim gates remain armed. */
            size_t retry = 0;
            for (int32_t c = 0; c < ncomp; c++) {
                if (comp_placed[c] != 1 && rim_deferred[c] && tried_placed[c] < st->placed) {
                    comp_placed[c] = 0; tried_at[c] = 0; retry++;
                }
            }
            if (retry > 0) {
                fprintf(stderr, "  [place] rim retry: %zu deferred components after new material was placed (all gates retained)\n", retry);
                best = apl_next_component(ncomp, comp_placed, anchor_count, seam_count,
                                           tried_at, comp_ring, comp_area, o->min_anchors);
            }
        }
        if (best < 0) {
            /* A later seam may disprove an earlier hop. Revisit those
             * incumbents before declaring the remaining material stranded. */
            if(retained) break;
            int recovered = 0;
            size_t epoch = st->placed+source_recoveries+1;
            for (int32_t c = 0; c < ncomp; c++) {
                if (ves_clock_sec()-t_sweep0 > budget_sec) { st->budget_exhausted = 1; break; }
                if (reconsidered_at[c] == epoch) continue;
                reconsidered_at[c] = epoch; x.cand_grad = ori[c].grad;
                if (apl_reconsider_hop(&x,c,t_sweep0+budget_sec)) {
                    source_recoveries++; recovered = 1;
                    apl_rebuild_grids(&x,&garena);
                    for (int32_t j = 0; j < ncomp; j++) if (comp_placed[j] != 1) {
                        comp_placed[j] = 0; tried_at[j] = 0;
                    }
                    placed_dirty = 1; break;
                }
            }
            if (recovered) continue;
            break;
        }
        tried_at[best] = anchor_count[best];
        tried_placed[best] = st->placed;
        rim_deferred[best] = 0;
        x.cand_grad = ori[best].grad;
        x.attempt++;
        diag[best].attempts++;
        /* a component whose mirror came from the chirality prior (no turn
         * vector of its own) may have a wrong bending vote (a fold, a concave
         * stretch): try both mirror states and keep the one the anchors fit
         * decisively better -- a wrong mirror fits at half the extent, never
         * close */
        if (!ori[best].has_turn) {
            double best_rms[2] = { 1e300, 1e300 };
            for (int state = 0; state < 2; state++) {
                if (state == 1) apl_mirror_component(run, best);
                x.stretch_ok = apl_stretch_turn(&x, best, &x.stretch_T, &x.stretch_s, &x.stretch_d, &x.stretch_n, &x.stretch_rel, &x.stretch_pairs);
                size_t na2 = apl_collect_anchors(&x, best);
                if (na2 >= (size_t)o->min_anchors || (APL_SEAM_FIRST && x.na_seam >= (size_t)APL_SEAM_MIN_ANCHORS)) {
                    double w2[APL_MAX_CLUSTERS_TOTAL], u2[APL_MAX_CLUSTERS_TOTAL], v2[APL_MAX_CLUSTERS_TOTAL]; int kd2[APL_MAX_CLUSTERS_TOTAL];
                    int n2 = apl_cluster_anchors(&x, na2, w2, u2, v2, kd2);
                    for (int cl = 0; cl < n2; cl++) {
                        double th2, tx2, ty2; size_t kept2 = 0;
                        double r2 = apl_fit_cluster(&x, na2, cl, kd2[cl] == 1 ? (size_t)APL_SEAM_MIN_ANCHORS : (size_t)o->min_anchors, o->max_rot_refine, &th2, &tx2, &ty2, &kept2, res);
                        if (r2 < best_rms[state]) best_rms[state] = r2;
                    }
                }
            }
            /* now mirrored; keep it only if decisively better */
            if (!(best_rms[1] < 0.5 * best_rms[0])) apl_mirror_component(run, best);
            else {
                ori[best].mirrored = !ori[best].mirrored;
                st->mirrored++;
                fprintf(stderr, "  [place] component %d: anchors prefer the mirrored state (rms %.1f vs %.1f)\n", best, best_rms[1], best_rms[0]);
            }
        }
        APL_TICK(t_mirror);
        x.stretch_ok = apl_stretch_turn(&x, best, &x.stretch_T, &x.stretch_s, &x.stretch_d, &x.stretch_n, &x.stretch_rel, &x.stretch_pairs);
        if (x.stretch_ok)
            fprintf(stderr, "  [place] component %d: turn length from the layer stretch: s %.4f +- %.4f over %zu hits on %zu layer pairs, gap %.1f vox -> one turn %.0f vox, (s-1) known to %.0f%% (%s)\n",
                    best, x.stretch_s, x.stretch_rel >= 0.0 ? x.stretch_rel * fabs(x.stretch_s - 1.0) : -1.0,
                    x.stretch_n, x.stretch_pairs, x.stretch_d, fabs(x.stretch_T),
                    x.stretch_rel >= 0.0 ? 100.0 * x.stretch_rel : -1.0,
                    x.stretch_T > 0.0 ? "candidate outward (bending-side vote)" : "candidate inward (bending-side vote)");
        APL_TICK(t_stretch);
        size_t na = apl_collect_anchors(&x, best);
        APL_TICK(t_collect);
        diag[best].seam_anchors = x.na_seam;
        if (na < (size_t)o->min_anchors && !(APL_SEAM_FIRST && x.na_seam >= (size_t)APL_SEAM_MIN_ANCHORS)) { comp_placed[best] = 2; diag[best].status = APL_FEW_ANCHORS; diag[best].last_anchors = na; continue; }
        double cl_w[APL_MAX_CLUSTERS_TOTAL], cl_u[APL_MAX_CLUSTERS_TOTAL], cl_v[APL_MAX_CLUSTERS_TOTAL], cl_rms[APL_MAX_CLUSTERS_TOTAL];
        int cl_kind[APL_MAX_CLUSTERS_TOTAL];
        int ncl = apl_cluster_anchors(&x, na, cl_w, cl_u, cl_v, cl_kind);
        APL_TICK(t_cluster);
        for (int cl = 0; cl < APL_MAX_CLUSTERS_TOTAL; cl++) { cl_rms[cl] = 1e300; if (cl >= ncl) cl_kind[cl] = 0; }
        diag[best].last_clusters = ncl;
        double gate = o->max_rms + o->rel_rms * sqrt(comp_area[best]);
        double mean_turn = 0.0;
        {
            /* the anchor scatter of ESTIMATED hops scales with the turn length, not with the
             * component's size: on the 10x10x10 (turns 2,000-5,000 vox) 43 sheets holding 10.7%
             * of the area were refused at rms 100+ against gates of 50-60 while the audits never
             * got to judge them.  The gate is at least this fraction of the anchors' mean hop. */
            double tsum = 0.0; size_t tn = 0;
            for (size_t i = 0; i < na; i++) { tsum += hypot(x.anchors[i].Tx, x.anchors[i].Ty); tn++; }
            if (tn > 0) mean_turn = tsum / (double)tn;
            if (tn > 0 && APL_GATE_TURN_FRAC > 0.0) gate = fmax(gate, APL_GATE_TURN_FRAC * tsum / (double)tn);
        }
        int placed_now = 0, status = APL_REFUSED;
        double tcl0 = ves_clock_sec();
        size_t winner_n = 0, winner_kept = 0;
        int winner_kind = -1, ambiguous = 0;
        double winner_rms = 0;
        AplCompDiag winner_diag; memset(&winner_diag, 0, sizeof winner_diag);
        x.source_geometry_conflict = 0;
        for (int cl = 0; cl < ncl; cl++) {
            x.source_flip = 0;
            for (size_t i = 0; i < na; i++) if (x.anchors[i].cluster == cl && x.anchors[i].k == 0) { x.source_flip = x.anchors[i].source_flip; break; }
            double th, tx, ty; size_t kept = 0;
            double tf0 = ves_clock_sec();
            double rms = apl_fit_cluster(&x, na, cl, cl_kind[cl] == 1 ? (size_t)APL_SEAM_MIN_ANCHORS : (size_t)o->min_anchors, o->max_rot_refine, &th, &tx, &ty, &kept, res);
            st->t_fit += ves_clock_sec() - tf0;
            cl_rms[cl] = rms;
            /* a seam cluster's gate is its seams' own precision, never the turn fraction; the last call's hop gate is the size gate alone */
            double gate_cl = last_call && cl_kind[cl] == 0 ? o->max_rms + o->rel_rms * sqrt(comp_area[best]) : gate;
            if (cl_kind[cl] == 1) {
                double smax = 0.0;
                for (size_t i = 0; i < na; i++) if (x.anchors[i].cluster == cl && x.anchors[i].seam_rms > smax) smax = x.anchors[i].seam_rms;
                gate_cl = fmax(APL_SEAM_GATE_MIN, APL_SEAM_GATE_RMS_X * smax);
            }
            size_t k1 = 0, k2 = 0, est = 0, rev = 0;
            for (size_t i = 0; i < na; i++) { if (x.wts[i] <= 0.0) continue; if (x.anchors[i].k == 1) k1++; else k2++; est += x.anchors[i].estimated; rev += x.anchors[i].reverse; }
            if (rms >= 1e299) {
                fprintf(stderr, "  [place] component %d (area %.3e) cluster %d/%d (weight %.1f, offset %.0f,%.0f): %zu anchors, too few\n",
                        best, comp_area[best], cl + 1, ncl, cl_w[cl], cl_u[cl], cl_v[cl], kept);
                apl_write_candidate(fcandidate, &x, best, comp_area[best], cl, cl_kind[cl], 0, "too_few_anchors",
                                    kept, -1.0, gate_cl, k1, k2, est, mean_turn, x.turn_rel, 0);
                continue;
            }
            fprintf(stderr, "  [place] component %d (area %.3e) cluster %d/%d (%s, weight %.1f, offset %.0f,%.0f): %zu anchors (k1 %zu k2+ %zu, estimated-T %zu, reverse %zu) rms %.2f gate %.2f rot %.3f -> %s\n",
                    best, comp_area[best], cl + 1, ncl, cl_kind[cl] == 1 ? "seam" : "hop", cl_w[cl], cl_u[cl], cl_v[cl], kept, k1, k2, est, rev, rms, gate_cl, th,
                    rms > gate_cl ? (cl_kind[cl] == 1 ? "repair" : "refused") : "fit");
            if (fanch) {
                for (size_t i = 0; i < na; i++) {
                    const AplAnchor *a = &x.anchors[i];
                    if (a->cluster != cl) continue;
                    fprintf(fanch, "%d,%d,%.6e,%d,%d,%d,%d,%d,%d,%d,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%d\n", x.attempt, best, comp_area[best], cl,
                            a->chart_p, a->chart_c, a->side, a->k, a->turn_src, a->reverse, a->Tx, a->Ty, a->sx, a->sy, a->tx, a->ty, res[i], x.wts[i] > 0.0);
                }
            }
            diag[best].last_anchors = kept; diag[best].last_rms = rms; diag[best].last_gate = gate_cl; diag[best].last_rot = th;
            diag[best].last_k1 = (int)k1; diag[best].last_k2 = (int)k2;
            if (rms > gate_cl && cl_kind[cl] != 1) {
                status = APL_REFUSED;
                apl_write_candidate(fcandidate, &x, best, comp_area[best], cl, cl_kind[cl], 0, "rms_over_gate",
                                    kept, rms, gate_cl, k1, k2, est, mean_turn, x.turn_rel, 0);
                continue;
            }
            /* A hop pose is a claim about WHERE the next wrap sits, and the turn is that
             * distance.  When the turn itself is not known to better than this, the pose
             * cannot be judged: fitting it against a gate that is a fraction of a turn nobody
             * measured is how a 342-vox gate comes to pass an rms of 39.6 (21x5x5, 2026-09-17).
             * Refuse it and say so, rather than deferring it later as ambiguous. */
            if (APL_HOP_TURN_MAX_REL_ERR > 0.0 && cl_kind[cl] != 1 && est > 0 &&
                x.turn_rel > APL_HOP_TURN_MAX_REL_ERR) {
                status = APL_REFUSED;
                fprintf(stderr, "  [place] component %d cluster %d/%d: the turn its %zu estimated anchors rest on is known to %.0f%%; refused\n",
                        best, cl + 1, ncl, est, 100.0 * x.turn_rel);
                apl_write_candidate(fcandidate, &x, best, comp_area[best], cl, cl_kind[cl], 0, "turn_unknown",
                                    kept, rms, gate_cl, k1, k2, est, mean_turn, x.turn_rel, 0);
                continue;
            }
            double ratio = 0.0;
            double ta0 = ves_clock_sec();
            x.na = na; x.cluster = cl; x.cluster_kind = cl_kind[cl]; x.evidence_gate = gate_cl;
            memset(x.seam_selected, 0, run->n_rels); x.source_flip = 0;
            for (size_t i = 0; i < na; i++) if (x.anchors[i].cluster == cl && x.anchors[i].k == 0 && x.anchors[i].relation >= 0) {
                x.seam_selected[x.anchors[i].relation] = 1;
                x.source_flip = x.anchors[i].source_flip;
            }
            int audit_ok = apl_place_and_audit(&x, best, th, tx, ty, &ratio, 2);
            if(audit_ok && retained && diagnostic_policy < 0 && !apl_retry_candidate_qualified(&x,best))audit_ok=0;
            if (audit_ok && x.strict && !aps_candidate(&x, best)) audit_ok = 0;
            st->t_audit += ves_clock_sec() - ta0;
            if (audit_ok) {
                fprintf(stderr, "  [place] candidate supported: component %d cluster %d, %zu witnesses at %.2f vox rms, %zu fitted anchors\n",
                        best, cl+1, x.proposal[0].placement_support, x.proposal[0].placement_rms, kept);
                apl_write_candidate(fcandidate, &x, best, comp_area[best], cl, cl_kind[cl], 1, "accepted",
                                    kept, rms, gate_cl, k1, k2, est, mean_turn, x.turn_rel, 1);
                int replace = winner_n == 0 || cl_kind[cl] > winner_kind;
                if (!replace && cl_kind[cl] == winner_kind) {
                    if (winner_n != x.nproposal) ambiguous = 1;
                    else for (size_t k = 0; k < winner_n; k++)
                        if (!apl_same_pose(&winner[k], &x.proposal[k], 8.0)) { ambiguous = 1; break; }
                }
                if (replace) {
                    winner_n = x.nproposal; winner_kind = cl_kind[cl]; ambiguous = 0;
                    apl_copy_candidate(winner,x.proposal,winner_n,run,&winner_uv);
                    if (winner_seams) memcpy(winner_seams, x.seam_selected, run->n_rels);
                    winner_kept = kept; winner_rms = rms; winner_diag = diag[best];
                }
            } else if (x.last_shape_refused) {
                apl_write_candidate(fcandidate, &x, best, comp_area[best], cl, cl_kind[cl], 0,
                                    x.last_reason ? x.last_reason : "shape_refused", kept, rms, gate_cl, k1, k2, est, mean_turn, x.turn_rel, 1);
                status = APL_REFUSED;
            } else if (x.last_silent) {
                apl_write_candidate(fcandidate, &x, best, comp_area[best], cl, cl_kind[cl], 0,
                                    x.last_reason ? x.last_reason : "deferred", kept, rms, gate_cl, k1, k2, est, mean_turn, x.turn_rel, 1);
                status = APL_DEFERRED; st->deferred++; rim_deferred[best] = 1;
            } else {
                apl_write_candidate(fcandidate, &x, best, comp_area[best], cl, cl_kind[cl], 0,
                                    x.last_reason ? x.last_reason : "reverted", kept, rms, gate_cl, k1, k2, est, mean_turn, x.turn_rel, 1);
                status = APL_REVERTED; st->reverted_by_audit++;
                if (APL_HOP_RETRY && ratio > APL_HOP_RETRY_MIN_RATIO) {
                    /* the cluster's mean ONE-hop vector (anchors carry the k-hop vector) */
                    double hx = 0.0, hy = 0.0, hw = 0.0;
                    for (size_t i = 0; i < na; i++) {
                        const AplAnchor *a = &x.anchors[i];
                        if (x.wts[i] <= 0.0 || a->cluster != cl || a->k < 1) continue;
                        hx += x.wts[i] * a->Tx / (double)a->k; hy += x.wts[i] * a->Ty / (double)a->k; hw += x.wts[i];
                    }
                    if (hw > 0.0 && hypot(hx, hy) > 1.0) {
                        hx /= hw; hy /= hw;
                        for (int dir = -1; dir <= 1 && !placed_now; dir += 2) {
                            double r2 = 0.0;
                            st->hop_retries++;
                            if (apl_place_and_audit(&x, best, th, tx + (double)dir * hx, ty + (double)dir * hy, &r2, APL_HOP_RETRY_REQUIRE_RIM)) {
                                placed_now = 1; st->hop_retry_placed++;
                                st->anchors_used += kept; rms_list[n_rms++] = rms;
                                fprintf(stderr, "  [place]   hop retry %+d ACCEPTED (one turn %.0f,%.0f vox; the fit's landing contradicted %.0f%% of its area)\n", dir, hx, hy, 100.0 * ratio);
                            } else
                                fprintf(stderr, "  [place]   hop retry %+d rejected (contradiction %.0f%%)\n", dir, 100.0 * r2);
                        }
                    }
                }
            }
        }
        if (winner_n && !ambiguous) {
            for (size_t k = 0; k < winner_n; k++) apl_commit_chart(&x,&winner[k],run->charts[winner[k].id].uv);
            if (x.strict) {
                aps_require_selected(run, winner_seams);
                for (size_t k = 0; k < winner_n; k++) if (!aps_box(run->charts + winner[k].id, x.strict->boxes + winner[k].id)) x.strict->valid = 0;
            }
            AsmContinuity_confirm(run);
            if (x.cgrid) for (size_t k = 0; k < winner_n; k++) {
                const AsmChart *c = &run->charts[winner[k].id];
                AsmConflictGrid_insert(x.cgrid, c); aprg_insert_chart(x.rgrid, c);
            }
            placed_now = 1; st->anchors_used += winner_kept; rms_list[n_rms++] = winner_rms;
            diag[best] = winner_diag; diag[best].placed_by = winner_kind == 1 ? 1 : 2;
            if (winner_kind == 1) st->seam_placed++; else st->hop_placed++;
        } else if (ambiguous) {
            fprintf(stderr, "  [place] multiple supported poses; candidate remains deferred\n");
            status = APL_DEFERRED; rim_deferred[best] = 1; st->deferred++;
        }
        x.source_flip = 0;
        st->t_clusterloop += ves_clock_sec() - tcl0;
        if (!retained && !placed_now && ncl >= 2 && n_cuts < APL_MAX_CUTS) {
            double tc0 = ves_clock_sec();
            double cl_rms_hop[APL_MAX_CLUSTERS_TOTAL];   /* a seam cluster against a 150-vox-wrong hop cluster is not two wraps */
            for (int cl = 0; cl < APL_MAX_CLUSTERS_TOTAL; cl++) cl_rms_hop[cl] = cl_kind[cl] == 1 ? 1e300 : cl_rms[cl];
            int32_t newc = apl_cut_component(&x, best, na, ncl, cl_w, cl_u, cl_rms_hop, gate, rel_off, rel_adj, &ncomp, ncomp_cap,
                                             comp_area, orig_comp, ori, diag, comp_ring, lidx, &st->cut_rels_deferred);
            st->t_cut += ves_clock_sec() - tc0;
            if (newc >= 0) {
                n_cuts++; st->cuts++;
                comp_placed[best] = 0; comp_placed[newc] = 0;
                tried_at[best] = 0; tried_at[newc] = 0;
                diag[best].status = APL_NEVER; diag[best].attempts = 0;
                continue;   /* both parts are candidates again */
            }
        }
        if (!placed_now) { comp_placed[best] = 2; diag[best].status = status; continue; }
        comp_placed[best] = 1;
        placed_dirty = 1;
        diag[best].status = APL_PLACED;
        if (!retained && APL_METRIC_RESOLVE && st->placed > 0 && (st->placed % APL_METRIC_EVERY) == 0) {
            size_t nobs = 0;
            apl_metric_resolve(&x, comp_placed, orig_comp, ncomp, lidx, primary, &garena, copts, &nobs, 0);
        }
        /* the pooled TURN MODEL alone, without the measured-dead component shifts */
        if (!retained && APL_TURN_MODEL_POOL && st->placed > 0 && (st->placed % APL_METRIC_EVERY) == 0) {
            size_t nobs = 0;
            apl_metric_resolve(&x, comp_placed, orig_comp, ncomp, lidx, primary, &garena, copts, &nobs, 1);
        }
        if (APL_COORD_CALIBRATE && coord_anchor >= 0 && st->placed > 0 && (st->placed % APL_METRIC_EVERY) == 0) {
            size_t ncal = apl_kover_build(&x, comp_placed, orig_comp, lidx, kover);
            size_t ns2 = 0, nl2 = 0, nk2 = 0; double q50 = 0.0, q95 = 0.0;
            apl_layer_index(run, out_side, turn, x.t0, x.growth, coord_anchor, lidx, kover, &ns2, &nl2, &nk2, &q50, &q95);
            apl_comp_ring_update(run, lidx, orig_comp, ncomp, comp_ring);
            if ((st->placed % (APL_METRIC_EVERY * 4)) == 0)
                fprintf(stderr, "  [place] coordinate re-calibrated on %zu placed pairs: %zu charts reached, edge residual p50 %.2f p95 %.2f layers\n", ncal, nk2, q50, q95);
        }
        if (APL_TURN_FIELD && (st->placed % APL_METRIC_EVERY) == 0) {
            size_t nf = apl_turn_field_build(&x, comp_placed, orig_comp, lidx);
            if (st->placed == 0 || (st->placed % (APL_METRIC_EVERY * 4)) == 0) fprintf(stderr, "  [place] turn field: %zu placed charts carry a local observed turn\n", nf);
        }
        for (size_t i = 0; i < nc; i++) {
            const AsmChart *ch = &run->charts[i];
            if (!AsmChart_in_layout(ch) || orig_comp[i] != best) continue;
            if (ch->pose_x < x.strip_umin) x.strip_umin = ch->pose_x;
            if (ch->pose_x > x.strip_umax) x.strip_umax = ch->pose_x;
        }
        st->placed++;
        st->area_placed += comp_area[best];
        if (last_call) st->last_call_placed++;
    }
    }
    st->t_loop = ves_clock_sec() - t_sweep0;
    apl_profile_line(st, st->sweeps, st->t_loop, "[assemble place] sweep profile");
    for (int32_t c = 0; c < ncomp; c++) if (comp_placed[c] != 1 && comp_area[c] > 0.0) st->unplaced++;
    /* CONFETTI LEDGER (2026-09-10): the cube-sized pieces of the outer wraps, by status and by the seam
     * evidence the unplaced ones still hold */
    {
        int32_t *nch = ARENA_CALLOC(arena, (size_t)(ncomp ? ncomp : 1), sizeof(int32_t));
        uint8_t *evid = ARENA_CALLOC(arena, (size_t)(ncomp ? ncomp : 1), sizeof(uint8_t));   /* bit 1 accepted-but-dropped / switched, 2 weak, 4 short */
        for (size_t i = 0; i < nc; i++) if (AsmChart_in_layout(&run->charts[i]) && orig_comp[i] >= 0 && orig_comp[i] < ncomp) nch[orig_comp[i]]++;
        for (size_t r = 0; r < run->n_rels; r++) {
            const AsmRelation *R = &run->rels[r];
            if (R->flags & ASM_REL_CONTACT) continue;
            const AsmChart *A = &run->charts[(size_t)R->a], *B = &run->charts[(size_t)R->b];
            if (!AsmChart_in_layout(A) || !AsmChart_in_layout(B)) continue;
            int32_t ca = orig_comp[R->a], cb = orig_comp[R->b];
            if (ca < 0 || cb < 0 || ca >= ncomp || cb >= ncomp || ca == cb) continue;
            int pa = comp_placed[ca] == 1, pb = comp_placed[cb] == 1;
            if (pa == pb) continue;
            uint8_t bit = (R->flags & ASM_REL_SHORT) ? 4 : (R->flags & ASM_REL_WEAK) ? 2 : (R->flags & (ASM_REL_DROPPED | ASM_REL_SWITCHED)) ? 1 : 0;
            if (bit) evid[pa ? cb : ca] |= bit;
        }
        for (int32_t c = 0; c < ncomp; c++) {
            if (comp_area[c] <= 0.0) continue;
            int s = diag[c].status; if (s < 0 || s > 6) s = 0;
            if (nch[c] == 1) { st->singles++; st->singles_area += comp_area[c]; st->singles_by_status[s]++; st->singles_area_by_status[s] += comp_area[c]; }
            if (s == APL_PLACED && diag[c].last_rms > 30.0) { st->placed_rms30++; st->placed_rms30_area += comp_area[c]; }
            if (s == APL_REVERTED) st->reverted_area += comp_area[c];
            if (s == APL_DEFERRED) { st->deferred_comps++; st->deferred_area += comp_area[c]; }
            if (comp_placed[c] != 1 && s != APL_PRIMARY && diag[c].attempts == 0) {
                st->never_attempted++; st->never_attempted_area += comp_area[c];
            }
            if (comp_placed[c] != 1 && s != APL_PRIMARY) {
                if (evid[c]) { st->evid_comps++; st->evid_area += comp_area[c]; if (evid[c] & 1) st->evid_dropped++; if (evid[c] & 2) st->evid_weak++; if (evid[c] & 4) st->evid_short++; }
                else { st->noevid_comps++; st->noevid_area += comp_area[c]; }
            }
        }
    }
    /* Placement-time UV repair holds established neighbours fixed. Once
     * all placements are finished, both sides can share a seam correction;
     * nearby fixed material participates in every geometry audit. The old
     * grids are no longer queried after this point. */
    apl_record_poses(run,o,"placed_before_uv",1);
    if(!retained){
        AsmContinuity_relax_placed_uv(run,apl_uv_audit,o->stress_band,o->diag_dir);
        AsmContinuity_confirm(run);
    }
    /* Re-level the whole placed strip: the primary's own z-gradient frame is
     * biased by its partial outer turns (an arc around a tilted axis reads a
     * spurious dz/du), and every arc inherited that frame through its
     * anchors, so the strip drifts in v along u.  One rotation over all
     * placed material puts +v back on world z. */
    if(!retained) {
        AplOrient lev;
        apl_orient_rotate(run, primary, &lev);
        st->relevel_theta = lev.theta;
        fprintf(stderr, "  [place] re-levelled the placed strip by %.4f rad (z-gradient %.3f)\n", lev.theta, lev.grad);
    }
    if (n_rms) {
        qsort(rms_list, n_rms, sizeof(double), apl_cmp_double);
        st->rms_p50 = rms_list[n_rms/2]; st->rms_max = rms_list[n_rms-1];
    }
    st->growth_measured = x.growth;
    apl_record_poses(run, o, "placed", 1);
    st->growth_charts = n_growth;
    if (fanch) fclose(fanch);
    if (fcandidate) fclose(fcandidate);
    if (fevidence) fclose(fevidence);
    if (frim) fclose(frim);
    if (fcut) fclose(fcut);
    if (o->diag_dir) {
        char path[2048];
        snprintf(path, sizeof path, "%s/stage5_components.csv", o->diag_dir);
        FILE *fc = fopen(path, "wb");
        if (fc) {
            fputs("component,area,status,attempts,clusters,last_anchors,last_k1,last_k2,last_rms,last_gate,last_rot,mirrored,turn_u,bend,z_gradient,"
                  "cz,cy,cx,z0,z1,y0,y1,x0,x1,u_mean,v_mean,placed_by,seam_anchors\n", fc);
            for (int32_t c = 0; c < ncomp; c++) {
                if (comp_area[c] <= 0.0) continue;
                /* area-weighted 3-D centroid, 3-D bbox and mean layout position of the component's charts */
                double cz = 0, cy = 0, cx = 0, um = 0, vm = 0, aw = 0;
                float lo[3] = { 1e30f, 1e30f, 1e30f }, hi[3] = { -1e30f, -1e30f, -1e30f };
                for (size_t i = 0; i < nc; i++) {
                    const AsmChart *ch = &run->charts[i];
                    if (!AsmChart_in_layout(ch) || orig_comp[i] != c) continue;
                    double gx, gy;
                    apl_pose(ch, 0.0, 0.0, &gx, &gy);
                    cz += ch->area3d * ch->centroid[0]; cy += ch->area3d * ch->centroid[1]; cx += ch->area3d * ch->centroid[2];
                    um += ch->area3d * gx; vm += ch->area3d * gy; aw += ch->area3d;
                    for (int d = 0; d < 3; d++) { if (ch->bbox_lo[d] < lo[d]) lo[d] = ch->bbox_lo[d]; if (ch->bbox_hi[d] > hi[d]) hi[d] = ch->bbox_hi[d]; }
                }
                if (aw > 0.0) { cz /= aw; cy /= aw; cx /= aw; um /= aw; vm /= aw; }
                fprintf(fc, "%d,%.6e,%s,%d,%d,%zu,%d,%d,%.2f,%.2f,%.4f,%d,%.3f,%.3f,%.3f,%.1f,%.1f,%.1f,%.0f,%.0f,%.0f,%.0f,%.0f,%.0f,%.1f,%.1f,%d,%zu\n",
                        c, comp_area[c], apl_status_names[diag[c].status], diag[c].attempts,
                        diag[c].last_clusters, diag[c].last_anchors, diag[c].last_k1, diag[c].last_k2, diag[c].last_rms, diag[c].last_gate, diag[c].last_rot,
                        ori[c].mirrored, ori[c].turn_u, ori[c].bend, ori[c].grad, cz, cy, cx, lo[0], hi[0], lo[1], hi[1], lo[2], hi[2], um, vm,
                        diag[c].placed_by, diag[c].seam_anchors);
            }
            fclose(fc);
        }
        snprintf(path, sizeof path, "%s/stage5_charts.csv", o->diag_dir);
        FILE *fch = fopen(path, "wb");
        if (fch) {
            fputs("chart,cube,orig_component,status,area,nv,cz,cy,cx,u_mean,v_mean,mirror,out_side,layer,s,r,phi\n", fch);
            for (size_t i = 0; i < nc; i++) {
                const AsmChart *ch = &run->charts[i];
                if (!AsmChart_in_layout(ch)) continue;
                double gx, gy;
                apl_pose(ch, 0.0, 0.0, &gx, &gy);
                int32_t oc = orig_comp[i];
                /* the axis frame of the centroid (s along the axis, r from it, phi about it); r = -1 without an axis */
                double s_ax = 0.0, r_ax = -1.0, phi_ax = 0.0;
                if (o->axis) {
                    double cen[3] = { ch->centroid[0], ch->centroid[1], ch->centroid[2] }, fr[3];
                    AsmAxis_frame(o->axis, cen, fr);
                    s_ax = fr[0]; r_ax = hypot(fr[1], fr[2]); phi_ax = atan2(fr[2], fr[1]);
                }
                fprintf(fch, "%d,%d,%d,%s,%.4e,%zu,%.1f,%.1f,%.1f,%.1f,%.1f,%d,%d,%.3f,%.1f,%.1f,%.4f\n", ch->id, ch->cube, oc,
                        oc >= 0 && oc < ncomp ? apl_status_names[diag[oc].status] : "none", ch->area3d, ch->nv,
                        ch->centroid[0], ch->centroid[1], ch->centroid[2], gx, gy, (ch->flags & ASM_CHART_MIRROR) != 0, out_side[i],
                        lidx[i] > APL_LAYER_UNKNOWN ? lidx[i] : -999.0, s_ax, r_ax, phi_ax);
            }
            fclose(fch);
        }
        fprintf(stderr, "  [place] cuts: %zu components split at disagreeing anchor clusters (%zu seam relations deferred), %d components now\n", st->cuts, st->cut_rels_deferred, ncomp);
        fprintf(stderr, "  [place] source continuations preserved against %zu conflicting layer-cut proposals\n",x.cuts_preserved_source);
        /* the biggest unplaced components, the coverage that is missing */
        int32_t top[10]; size_t ntop = 0;
        for (int32_t c = 0; c < ncomp; c++) {
            if (comp_placed[c] == 1 || comp_area[c] <= 0.0) continue;
            size_t pos = ntop;
            while (pos > 0 && comp_area[top[pos-1]] < comp_area[c]) pos--;
            if (pos >= 10) continue;
            for (size_t k = (ntop < 10 ? ntop : 9); k > pos; k--) top[k] = top[k-1];
            top[pos] = c;
            if (ntop < 10) ntop++;
        }
        for (size_t k = 0; k < ntop; k++) {
            int32_t c = top[k];
            size_t from_placed = 0, with_turn = 0, from_unplaced = 0;
            int32_t src[3] = { -1, -1, -1 }; size_t srcn[3] = { 0, 0, 0 };
            for (size_t p = 0; p < run->n_layers; p++) {
                const AsmLayerPair *lp = &run->layers[p];
                const AsmChart *A = &run->charts[(size_t)lp->a], *B = &run->charts[(size_t)lp->b];
                if (!AsmChart_in_layout(A) || !AsmChart_in_layout(B) || B->component != c) continue;
                if (comp_placed[A->component] != 1) { from_unplaced += (size_t)lp->hit_count; continue; }
                from_placed += (size_t)lp->hit_count;
                double tx = 0.0, ty = 0.0;
                if (apl_turn_for(turn, lp->a, lp->side, out_side[lp->a], x.growth, &tx, &ty)) with_turn += (size_t)lp->hit_count;
                else if (turn[(size_t)lp->a * 2].n == 0 && turn[(size_t)lp->a * 2 + 1].n == 0) with_turn += (size_t)lp->hit_count;
                int32_t oc = orig_comp[lp->a];
                int slot = -1;
                for (int s2 = 0; s2 < 3; s2++) if (src[s2] == oc) slot = s2;
                if (slot < 0) for (int s2 = 0; s2 < 3; s2++) if (src[s2] < 0) { src[s2] = oc; slot = s2; break; }
                if (slot >= 0) srcn[slot] += (size_t)lp->hit_count;
            }
            fprintf(stderr, "  [place] unplaced #%zu: component %d area %.3e %s (attempts %d, clusters %d, last anchors %zu, rms %.1f / gate %.1f) | layer hits from placed %zu (with a turn vector %zu) from unplaced %zu | sources: comp %d x%zu, comp %d x%zu, comp %d x%zu\n",
                    k + 1, c, comp_area[c], apl_status_names[diag[c].status], diag[c].attempts, diag[c].last_clusters, diag[c].last_anchors, diag[c].last_rms, diag[c].last_gate,
                    from_placed, with_turn, from_unplaced, src[0], srcn[0], src[1], srcn[1], src[2], srcn[2]);
        }
    }
    st->anchors_radius = x.anchors_radius;
    st->axis_armed = !(x.vm_ok && !x.vm_pin);
    st->axis_b = x.vm_ok ? x.vm_b : 0.0;
    /* WOBBLE: the layout against the sheet's own v(s): the primary's residual, the placed area
     * beyond 30 / 100 vox of the line, the spine's bin-to-bin steps, the v span */
    /* the strip was re-levelled after placement: refit the primary's plane in the final frame
     * before measuring anything against it (a 0.006 rad relevel moves v by 350 vox at u 60,000) */
    if (x.vm_ok) apl_v_model(&x);
    if (x.vm_ok) {
        Arena_Mark wm = Arena_save(arena);
        double *res = ARENA_ALLOC(arena, (nc ? nc : 1) * sizeof(double));
        double *wres = ARENA_ALLOC(arena, (nc ? nc : 1) * sizeof(double));
        size_t np = 0, npl = 0;
        double w30 = 0.0, w100 = 0.0, wpl = 0.0, vmin = 1e300, vmax = -1e300, umin = 1e300, umax = -1e300;
        for (size_t i = 0; i < nc; i++) {
            const AsmChart *c = &run->charts[i];
            if (!AsmChart_in_layout(c) || c->component != primary) continue;
            double gx, gy;
            apl_layout_centroid(&x, c, i, 0.0, 0.0, 0.0, &gx, &gy);
            double r = x.vm_a + x.vm_b * x.s_axis[i] + x.vm_c * gx - gy;
            if (gy < vmin) vmin = gy;
            if (gy > vmax) vmax = gy;
            if (gx < umin) umin = gx;
            if (gx > umax) umax = gx;
            if (orig_comp[i] == primary) res[np++] = fabs(r);
            else { wres[npl] = fabs(r); npl++; wpl += c->area3d; if (fabs(r) > 30.0) w30 += c->area3d; if (fabs(r) > 100.0) w100 += c->area3d; }
        }
        if (np > 0) { qsort(res, np, sizeof(double), apl_cmp_double); st->wobble_primary_p50 = res[np / 2]; st->wobble_primary_p90 = res[(size_t)((double)(np - 1) * 0.9)]; }
        st->placed_area_beyond30 = wpl > 0.0 ? w30 / wpl : 0.0;
        st->placed_area_beyond100 = wpl > 0.0 ? w100 / wpl : 0.0;
        st->v_span = vmax > vmin ? vmax - vmin : 0.0;
        /* the spine: per 2,000-vox u bin, the area-weighted median residual of on-line material */
        if (umax > umin) {
            size_t nb = (size_t)((umax - umin) / 2000.0) + 1;
            double *steps = ARENA_ALLOC(arena, (nb ? nb : 1) * sizeof(double));
            double prev = 0.0; int have_prev = 0; size_t ns = 0;
            for (size_t b = 0; b < nb; b++) {
                double lo = umin + 2000.0 * (double)b, hi = lo + 2000.0;
                size_t m = 0; double warea = 0.0;
                for (size_t i = 0; i < nc; i++) {
                    const AsmChart *c = &run->charts[i];
                    if (!AsmChart_in_layout(c) || c->component != primary) continue;
                    double gx, gy;
                    apl_layout_centroid(&x, c, i, 0.0, 0.0, 0.0, &gx, &gy);
                    if (gx < lo || gx >= hi) continue;
                    double r = x.vm_a + x.vm_b * x.s_axis[i] + x.vm_c * gx - gy;
                    if (fabs(r) > 100.0) continue;
                    res[m++] = r; warea += c->area3d;
                }
                if (m < 3 || warea < 5.0e4) continue;
                qsort(res, m, sizeof(double), apl_cmp_double);
                double med = res[m / 2];
                if (have_prev) steps[ns++] = fabs(med - prev);
                prev = med; have_prev = 1;
            }
            if (ns > 0) { qsort(steps, ns, sizeof(double), apl_cmp_double); st->spine_step_p50 = steps[ns / 2]; st->spine_step_p90 = steps[(size_t)((double)(ns - 1) * 0.9)]; st->spine_bins = ns + 1; }
        }
        /* OUTLINE: the trace of the box's top and bottom crop planes in the sheet.  In the axis
         * frame the planes z = const are tilted against the sheet, so their trace waves by
         * 2 r tan(tilt) per turn while the sheet itself is straight (the spine above); the wave is
         * reported beside that prediction so the eye can tell a crop edge from a wobble.  The
         * trace, not the material's edge: a chart missing at the top of a bin moves the edge by a
         * cube (the familiar's edge read 384-537 vox p2p against a 102-vox prediction). */
        if (umax > umin) {
            size_t nb = (size_t)((umax - umin) / 2000.0) + 1;
            double rsum = 0.0, rw = 0.0, zlo = 1e300, zhi = -1e300;
            for (size_t i = 0; i < nc; i++) {
                const AsmChart *c = &run->charts[i];
                if (!AsmChart_in_layout(c) || c->component != primary || c->uv == NULL || c->xyz == NULL) continue;
                for (size_t k = 0; k < c->nv; k++) {
                    double z = (double)c->xyz[k*3];
                    if (z < zlo) zlo = z;
                    if (z > zhi) zhi = z;
                }
                if (x.have_radius) { rsum += c->area3d * x.r_axis[i]; rw += c->area3d; }
            }
            st->outline_r_mean = rw > 0.0 ? rsum / rw : 0.0;
            if (zhi > zlo + 2.0 * APL_OUTLINE_CROP_VOX) {
                /* per bin, the layout v of the vertices on each crop plane: counts then medians */
                size_t *cnt_t = ARENA_CALLOC(arena, nb, sizeof(size_t)), *cnt_b = ARENA_CALLOC(arena, nb, sizeof(size_t));
                size_t tot_t = 0, tot_b = 0;
                for (int pass = 0; pass < 2; pass++) {
                    for (size_t i = 0; i < nc; i++) {
                        const AsmChart *c = &run->charts[i];
                        if (!AsmChart_in_layout(c) || c->component != primary || c->uv == NULL || c->xyz == NULL) continue;
                        double ct = cos(c->pose_theta), stt = sin(c->pose_theta);
                        double mir = (c->flags & ASM_CHART_MIRROR) ? -1.0 : 1.0;
                        for (size_t k = 0; k < c->nv; k++) {
                            double z = (double)c->xyz[k*3];
                            int on_top = z >= zhi - APL_OUTLINE_CROP_VOX, on_bot = z <= zlo + APL_OUTLINE_CROP_VOX;
                            if (!on_top && !on_bot) continue;
                            double u = mir * c->uv[k*2], v = c->uv[k*2+1];
                            double gx = ct * u - stt * v + c->pose_x, gy = stt * u + ct * v + c->pose_y;
                            if (gx < umin || gx >= umax) continue;
                            size_t b = (size_t)((gx - umin) / 2000.0);
                            if (b >= nb) b = nb - 1;
                            if (pass == 0) { if (on_top) cnt_t[b]++; else cnt_b[b]++; }
                            else if (on_top) x.tr_top[x.off_t[b] + x.fill_t[b]++] = gy;
                            else x.tr_bot[x.off_b[b] + x.fill_b[b]++] = gy;
                        }
                    }
                    if (pass == 0) {
                        x.off_t = ARENA_ALLOC(arena, (nb + 1) * sizeof(size_t)); x.off_b = ARENA_ALLOC(arena, (nb + 1) * sizeof(size_t));
                        x.fill_t = ARENA_CALLOC(arena, nb, sizeof(size_t)); x.fill_b = ARENA_CALLOC(arena, nb, sizeof(size_t));
                        x.off_t[0] = 0; x.off_b[0] = 0;
                        for (size_t b = 0; b < nb; b++) { x.off_t[b+1] = x.off_t[b] + cnt_t[b]; x.off_b[b+1] = x.off_b[b] + cnt_b[b]; }
                        tot_t = x.off_t[nb]; tot_b = x.off_b[nb];
                        x.tr_top = ARENA_ALLOC(arena, (tot_t ? tot_t : 1) * sizeof(double));
                        x.tr_bot = ARENA_ALLOC(arena, (tot_b ? tot_b : 1) * sizeof(double));
                    }
                }
                double tmin = 1e300, tmax = -1e300, bmin = 1e300, bmax = -1e300;
                size_t nbt = 0, nbb = 0;
                for (size_t b = 0; b < nb; b++) {
                    if (cnt_t[b] >= 8) {
                        double *seg = x.tr_top + x.off_t[b];
                        qsort(seg, cnt_t[b], sizeof(double), apl_cmp_double);
                        double med = seg[cnt_t[b] / 2];
                        if (med < tmin) tmin = med;
                        if (med > tmax) tmax = med;
                        nbt++;
                    }
                    if (cnt_b[b] >= 8) {
                        double *seg = x.tr_bot + x.off_b[b];
                        qsort(seg, cnt_b[b], sizeof(double), apl_cmp_double);
                        double med = seg[cnt_b[b] / 2];
                        if (med < bmin) bmin = med;
                        if (med > bmax) bmax = med;
                        nbb++;
                    }
                }
                st->outline_top_p2p = nbt > 1 ? tmax - tmin : 0.0;
                st->outline_bot_p2p = nbb > 1 ? bmax - bmin : 0.0;
                st->outline_bins_top = nbt; st->outline_bins_bot = nbb;
            }
            if (o->axis != NULL && zhi >= zlo) {
                double dir[3] = { 1.0, 0.0, 0.0 };
                AsmAxis_mean_dir(o->axis, zlo, zhi, dir);
                double cz = fabs(dir[0]);
                if (cz > 1.0) cz = 1.0;
                st->axis_tilt_deg = acos(cz) * 180.0 / APL_PI;
                st->outline_pred_p2p = 2.0 * st->outline_r_mean * tan(acos(cz));
            }
        }
        Arena_restore(arena, wm);
    }
    st->v_model[0] = x.vm_ok ? x.vm_a : 0.0; st->v_model[1] = x.vm_ok ? x.vm_b : 0.0; st->v_model[2] = x.vm_ok ? x.vm_c : 0.0;
    st->sec = ves_clock_sec() - t0;
    if (x.cgrid != NULL) AsmConflictGrid_dispose(x.cgrid);
    if (garena != NULL) Arena_dispose(&garena);
    if (x.proposal_uv) Arena_dispose(&x.proposal_uv);
    if (winner_uv) Arena_dispose(&winner_uv);
    int strict_rc = 0;
    if (x.strict) {
        fprintf(stderr, "[place strict] %zu/%zu candidate certificates passed in %.3f s\n", x.strict->accepted, x.strict->checked, x.strict->seconds);
        if (x.strict->log && fclose(x.strict->log)) x.strict->io = 0;
        strict_rc = !x.strict->io || !x.strict->valid;
    }
    Arena_restore(arena, mark);
    return strict_rc ? -1 : 0;
}

int AsmPlace_run(AsmRun *run,const AsmPlaceOpts *o,const AsmConflictOpts *copts,AsmPlaceStats *st)
{return apl_run(run,o,copts,st,NULL,-1,0);}

/* Retry one complete missing component against the actual repaired sheet.
 * The ordinary placement machinery owns all pose/evidence decisions. This
 * wrapper only supplies an immutable global context and restricts the scope;
 * the admission caller still has to qualify every original face and seam. */
int AsmPlace_retry_component(AsmRun *run,const uint8_t *core,const AsmPlaceOpts *o,
                             const AsmConflictOpts *copts,AsmPlaceStats *st)
{
    if(!run || !core || !o || !copts || !st || !run->continuity_ready)return -1;
    Arena_T arena=Arena_new();size_t nc=run->n_charts,count=0,fixed=0;
    AsmRun trial=*run;trial.arena=arena;
    trial.charts=ARENA_ALLOC(arena,nc*sizeof(AsmChart));memcpy(trial.charts,run->charts,nc*sizeof(AsmChart));
    trial.rels=ARENA_ALLOC(arena,run->n_rels*sizeof(AsmRelation));memcpy(trial.rels,run->rels,run->n_rels*sizeof(AsmRelation));
    uint8_t *retained=ARENA_CALLOC(arena,nc,1);int result=0;
    for(size_t c=0;c<nc;c++){
        AsmChart *ch=trial.charts+c;
        if(core[c]){
            if(AsmChart_registered(ch) || !AsmChart_in_layout(ch) || ch->placed_uv){result=-1;goto done;}
            ch->component=1;ch->placed=1;ch->placement_state=ASM_PLACE_NONE;count++;
            ch->uv=ARENA_ALLOC(arena,2*ch->nv*sizeof(float));memcpy(ch->uv,run->charts[c].uv,2*ch->nv*sizeof(float));
        }else if(AsmChart_registered(ch)){
            retained[c]=1;fixed++;ch->component=0;
        }else{ch->flags|=ASM_CHART_EXTRAS;ch->component=-1;}
    }
    if(!count || !fixed)goto done;
    /* A SOURCE continuation to an omitted chart must use anchored admission,
     * and a partial missing component must never masquerade as a layer arc. */
    for(size_t i=0;i<run->n_rels;i++){
        const AsmRelation *r=run->rels+i;
        if(!(r->flags&ASM_REL_CONTACT) && (r->continuity&(ASM_CONT_SOURCE|ASM_CONT_CUT)) &&
            core[r->a]!=core[r->b])goto done;
    }
    if(apl_run(&trial,o,copts,st,retained,-1,0)){result=-1;goto done;}
    for(size_t c=0;c<nc;c++)if(core[c] && !AsmChart_registered(trial.charts+c))goto done;
    for(size_t c=0;c<nc;c++)if(retained[c]){
        AsmChart after=trial.charts[c];after.component=run->charts[c].component;
        if(memcmp(&after,run->charts+c,sizeof after)){result=-1;goto done;}
    }
    for(size_t c=0;c<nc;c++)if(core[c]){
        AsmChart after=trial.charts[c];
        after.uv=ARENA_ALLOC(run->arena,2*after.nv*sizeof(float));
        memcpy(after.uv,trial.charts[c].uv,2*after.nv*sizeof(float));
        run->charts[c]=after;
    }
    result=1;
done:
    Arena_dispose(&arena);return result;
}

/* ---- selftest --------------------------------------------------------------------- */

static int apl_registration_selftest(void)
{
    int fails = 0;
    for (int mirror = 0; mirror < 2; mirror++) for (int collision = 0; collision < 2; collision++) {
        AsmRun run = {0}; run.arena = Arena_new(); run.layer_d_global = 12;
        float xyz[2][243], uv[2][162], normals[243]; uint8_t boundary[81]; int32_t faces[384];
        for (int j = 0; j < 9; j++) for (int u = 0; u < 9; u++) {
            int v = 9*j+u; boundary[v] = j == 0 || j == 8 || u == 0 || u == 8;
            normals[3*v] = normals[3*v+1] = 0; normals[3*v+2] = 1;
            for (int side = 0; side < 2; side++) {
                xyz[side][3*v] = (float)(16*j);
                xyz[side][3*v+1] = (float)(5*u+(side && !collision ? 44 : 0));
                xyz[side][3*v+2] = (float)(side && collision ? 12 : 0);
                uv[side][2*v] = (float)(5*u); uv[side][2*v+1] = (float)(16*j);
            }
            if (j < 8 && u < 8) { int32_t f[] = {v,v+1,v+10,v,v+10,v+9}; memcpy(faces+6*(8*j+u),f,sizeof f); }
        }
        for (int side = 0; side < 2; side++) {
            AsmChart c = {0}; c.id = side; c.nv = 81; c.nf = 128; c.area3d = side ? 5120 : 5121;
            c.xyz = xyz[side]; c.uv = uv[side]; c.nrm = normals; c.faces = faces; c.boundary = boundary;
            c.flags = mirror ? ASM_CHART_MIRROR : 0; c.pose_x = side*(mirror ? -200 : 200);
            AsmRun_push_chart(&run,&c);
        }
        AsmRelation rel = {0}; rel.a = 0; rel.b = 1; rel.corr_count = rel.n_corr = 9; rel.rms = 0.1;
        for (int j = 0; j < 9; j++) { AsmCorr c = {9*j+(collision ? 0 : 8),9*j}; AsmRun_push_corr(&run,&c); }
        AsmRun_push_rel(&run,&rel); AsmContinuity_prepare(&run);
        AsmChart saved[2] = {run.charts[0],run.charts[1]};
        AsmRelation saved_rel = run.rels[0]; AsmCorr saved_corr[9]; memcpy(saved_corr,run.corr,sizeof saved_corr);
        size_t before = apl_internal_bad(&run,saved,2,0);
        size_t kept = AsmContinuity_register_components(&run,apl_registration_audit);
        size_t after = apl_internal_bad(&run,saved,2,1);
        int ok = kept == (size_t)!collision && after <= before && !memcmp(&run.charts[0],&saved[0],sizeof saved[0]) &&
            !memcmp(&run.rels[0],&saved_rel,sizeof saved_rel) && !memcmp(run.corr,saved_corr,sizeof saved_corr);
        if (collision) ok = ok && !memcmp(run.charts,saved,sizeof saved);
        else {
            double rms; size_t support;
            ok = ok && fabs(run.charts[1].pose_x-(mirror ? -44 : 44)) < 1e-5 &&
                 AsmContinuity_measure(&run,&run.rels[0],&rms,&support) && rms < 1e-5 && support == 9;
        }
        if (!ok) { fprintf(stderr,"  asm_place selftest FAIL: registration collision %d mirror %d, kept %zu, bad %zu -> %zu\n",collision,mirror,kept,before,after); fails++; }
        Arena_dispose(&run.arena);
    }
    return fails;
}

static int apl_registration_reject_test(AsmRun *run, const AsmChart *saved, size_t n, AsmContinuityAuditCache *cache)
{
    (void)run; (void)saved; (void)n; (void)cache; return 0;
}

static int apl_registration_planes_selftest(void)
{
    /* The unconstrained fit closes a weak seam by compressing the gap to
     * a fourth chart. A feasible fit shares that motion with the fourth
     * chart and preserves all four SOURCE seams and every original pair. */
    int fails = 0;
    for (int mirror = 0; mirror < 2; mirror++) for (int stalled = 0; stalled < 2; stalled++) {
        AsmRun run = {0}; run.arena = Arena_new();
        float uv[4][12] = {{0}}, xyz[4][18] = {{0}}, original_uv[4][12], original_xyz[4][18];
        int32_t faces[] = {0,1,3,0,3,2,2,3,5,2,5,4};
        for (int i = 0; i < 4; i++) {
            for (int v = 0; v < 6; v++) {
                uv[i][2*v] = (float)(2*(v%2)-1); uv[i][2*v+1] = (float)(2*(v/2)-2);
                xyz[i][3*v] = uv[i][2*v]+(i == 3 && !stalled ? 20 : 0); xyz[i][3*v+1] = uv[i][2*v+1];
            }
            AsmChart c = {0}; c.id = i; c.nv = 6; c.nf = 4; c.uv = uv[i]; c.xyz = xyz[i]; c.faces = faces;
            c.area3d = i ? 100 : 101; c.pose_theta = .4; c.flags = mirror ? ASM_CHART_MIRROR : 0;
            if (i == 3 && !stalled) { c.pose_x = (mirror ? -18 : 18)*cos(.4); c.pose_y = (mirror ? -18 : 18)*sin(.4); }
            AsmRun_push_chart(&run,&c);
        }
        for (int i = 0; i < 4; i++) {
            AsmRelation r = {0}; r.a = i == 1 ? 1 : 0; r.b = i == 0 ? 1 : i == 3 ? 3 : 2;
            r.corr_first = (int32_t)run.n_corr; r.corr_count = r.n_corr = 6; r.continuity = ASM_CONT_SOURCE;
            r.flags = i == 2 ? ASM_REL_WEAK : 0;
            for (int v = 0; v < 6; v++) {
                AsmCorr c = {0}; c.va = c.vb = v; c.valid = 3;
                c.gap_a[0] = c.gap_b[0] = i == 2 ? 10 : i == 3 ? 18 : 0;
                if (i == 3 && stalled) {
                    /* A sub-voxel, already passing seam asks the final
                     * chart for a rotation beyond the .25-radian bound.
                     * A shorter bounded step may move it while closing
                     * the weak seam; the complete audit must still pass. */
                    double u = uv[0][2*v], w = uv[0][2*v+1];
                    c.gap_a[0] = c.gap_b[0] = (float)((cos(.3)-1)*u-sin(.3)*w);
                    c.gap_a[1] = c.gap_b[1] = (float)(sin(.3)*u+(cos(.3)-1)*w);
                }
                AsmRun_push_corr(&run,&c);
            }
            AsmRun_push_rel(&run,&r);
        }
        AsmChart saved[4]; AsmRelation rels[4]; AsmCorr corr[24];
        memcpy(saved,run.charts,sizeof saved); memcpy(rels,run.rels,sizeof rels); memcpy(corr,run.corr,sizeof corr);
        memcpy(original_uv,uv,sizeof uv); memcpy(original_xyz,xyz,sizeof xyz);
        AsmContinuityAuditCache cache = {0}; cache.guard_charts = saved; cache.guard_count = 4;
        int ok = apl_uv_audit(&run,saved,4,&cache) && !cache.baseline_bad;
        size_t kept = AsmContinuity_register_failed(&run,saved,apl_registration_audit);
        ok = ok && kept == 1 && apl_uv_audit(&run,saved,4,&cache) &&
             !memcmp(run.charts,saved,sizeof saved[0]) && !memcmp(run.rels,rels,sizeof rels) &&
             !memcmp(run.corr,corr,sizeof corr) && !memcmp(uv,original_uv,sizeof uv) && !memcmp(xyz,original_xyz,sizeof xyz);
        for (int i = 0; i < 4; i++) {
            double rms; size_t support;
            if (!AsmContinuity_measure(&run,&run.rels[i],&rms,&support)) ok = 0;
        }
        for (int i=1; i<4; i++) ok = ok &&
            hypot(run.charts[i].pose_x-saved[i].pose_x,run.charts[i].pose_y-saved[i].pose_y)<=256 &&
            fabs(run.charts[i].pose_theta-saved[i].pose_theta)<=.25;
        if (stalled) {
            memcpy(run.charts,saved,sizeof saved);
            ok = ok && !AsmContinuity_register_failed(&run,saved,apl_registration_reject_test) && !memcmp(run.charts,saved,sizeof saved);
        }
        if (!ok) { fprintf(stderr,"  asm_place selftest FAIL: constrained registration mirror %d stalled %d, kept %zu\n",mirror,stalled,kept); fails++; }
        Arena_dispose(&run.arena);
    }
    return fails;
}

static int apl_seam_incumbent_selftest(void)
{
    int fails = 0;
    for (int first_step = 0; first_step < 2; first_step++)
    for (int mirror = 0; mirror < 2; mirror++) for (int test = 0; test < 2; test++) {
        AsmRun run = {0}; run.arena = Arena_new(); run.layer_d_global = 12; run.continuity_ready = 1;
        float uv[3][36], xyz[3][54], normals[54]; int32_t faces[48]; uint8_t boundary[18];
        AsmChart charts[3] = {{0}}, original[3], proposal[1]; AsmCorr corr[18] = {{0}};
        double sign = mirror ? -1 : 1, ct = cos(test ? .08 : .5), sn = sin(test ? .08 : .5);
        AplAnchor anchors[54] = {{0}}; double weights[54] = {0};
        for (int v = 0; v < 18; v++) {
            double u = (test ? 0 : -5)+10.0*v/17, w = (v%2 ? 1 : -1)*(test ? .01 : .2);
            uv[1][2*v] = (float)u; uv[1][2*v+1] = (float)w;
            uv[0][2*v] = (float)(ct*u-sn*w); uv[0][2*v+1] = (float)(sn*u+ct*w);
            uv[2][2*v] = (float)u; uv[2][2*v+1] = (float)(24.05+w);
            for (int side = 0; side < 3; side++) {
                xyz[side][3*v] = (float)(1000+sign*uv[0][2*v]); xyz[side][3*v+1] = uv[0][2*v+1]; xyz[side][3*v+2] = side == 2 ? 100 : 0;
            }
            normals[3*v] = normals[3*v+1] = 0; normals[3*v+2] = 1; boundary[v] = 1;
            corr[v].va = corr[v].vb = v; corr[v].valid = 3; corr[v].run = 0;
            AplAnchor *a = &anchors[v]; a->chart_p = 0; a->chart_c = 1; a->relation = 0;
            a->sx = sign*uv[1][2*v]; a->sy = uv[1][2*v+1];
            a->local_u = uv[1][2*v]; a->local_v = uv[1][2*v+1];
            a->tx = 1000+sign*uv[0][2*v]; a->ty = uv[0][2*v+1]; a->w = weights[v] = 1;
            if (v < 16) { faces[3*v] = v; faces[3*v+1] = v+(v%2 ? 2 : 1); faces[3*v+2] = v+(v%2 ? 1 : 2); }
        }
        for (int i = 0; i < 3; i++) {
            charts[i].id = charts[i].component = i; charts[i].nv = 18; charts[i].nf = 16; charts[i].area3d = 4;
            charts[i].uv = uv[i]; charts[i].xyz = xyz[i]; charts[i].nrm = normals; charts[i].faces = faces; charts[i].boundary = boundary;
            charts[i].flags = mirror ? ASM_CHART_MIRROR : 0;
        }
        charts[0].pose_x = 1000; charts[0].placement_state = ASM_PLACE_ROOT;
        charts[2].pose_x = 1000; charts[2].component = 0; charts[2].placement_state = ASM_PLACE_ROOT;
        memcpy(original,charts,sizeof original);
        AsmRelation rel = {0}; rel.a = 0; rel.b = 1; rel.flags = ASM_REL_WEAK; rel.continuity = ASM_CONT_SOURCE;
        rel.corr_count = rel.n_corr = 18; rel.rms = 1;
        run.charts = charts; run.n_charts = test ? 3 : 2; run.rels = &rel; run.n_rels = 1; run.corr = corr; run.n_corr = 18;
        AsmPlaceOpts opts; AsmPlace_default_opts(&opts); opts.max_rot_refine = .1;
        AsmConflictOpts copts; AsmConflict_default_opts(&copts);
        int32_t groups[] = {0,1,2}; uint8_t selected[] = {1}; double area[] = {4,4,4};
        AplCtx x = {0}; x.run = &run; x.o = &opts; x.copts = &copts; x.orig_comp = groups; x.primary = 0;
        x.trial = 1; x.proposal = proposal; x.cluster_kind = 1; x.anchors = anchors; x.wts = weights; x.na = 18;
        x.seam_selected = selected; x.comp_area = area; x.diagnostic_first_step = first_step;
        Arena_T grid_arena = Arena_new();
        x.cgrid = AsmConflictGrid_new(grid_arena,&copts,&run); x.rgrid = aprg_new(grid_arena,24);
        AsmConflictGrid_insert(x.cgrid,&charts[0]); aprg_insert_chart(x.rgrid,&charts[0]);
        if (test) { AsmConflictGrid_insert(x.cgrid,&charts[2]); aprg_insert_chart(x.rgrid,&charts[2]); }
        charts[1].pose_x = 1000; charts[1].component = 0;
        double rms; size_t support, touched; const AsmChart *query = &charts[1];
        int ok = AsmContinuity_measure(&run,&rel,&rms,&support) && support == 18 && rms < 2 &&
                 aprg_probe(x.rgrid,&run,&query,1,24,&touched) == 0 && touched >= 8;
        charts[1] = original[1];
        double th = 0, tx = 1000, ty = 0;
        /* Both original proposals pass SOURCE and the real geometry audit.
         * The first closure exceeds its .1-radian bound. The second closes
         * successfully, but pulls a thin nondegenerate strip within the
         * rim reach of physically distant held material. All its original
         * vertices, including unsampled ones, are outside that reach. */
        int closed = apl_close_seam(&x,&th,&tx,&ty);
        ok = ok && closed == test;
        if (test) {
            apl_compose(&charts[1],th,tx,ty,0); charts[1].component = 0;
            ok = ok && aprg_probe(x.rgrid,&run,&query,1,24,&touched) > opts.max_rim_ratio && touched >= 8;
            charts[1] = original[1];
        }
        int accepted = apl_place_and_audit(&x,1,0,1000,0,NULL,2);
        ok = ok && accepted && x.nproposal == 1 && proposal[0].pose_x == 1000 && proposal[0].pose_theta == 0 &&
             proposal[0].placement_state == ASM_PLACE_SEAM && proposal[0].placement_support == 18 &&
             !memcmp(charts,original,sizeof original);
        /* Physical SOURCE support is independent of a contrary layer model. */
        for (int v = 0; v < 18; v++) {
            anchors[18+v] = anchors[v]; anchors[18+v].k = anchors[18+v].measured_order = 1;
            anchors[18+v].cluster = 1; anchors[18+v].tx += 2000;
            anchors[18+v].Tx = 1000; anchors[18+v].observation = (size_t)v;
        }
        x.na = 36; x.n_observations = 18;
        int source_precedes = apl_place_and_audit(&x,1,0,1000,0,NULL,2);
        ok = ok && source_precedes && x.consensus.observations == 18 && !x.consensus.model_inliers &&
             proposal[0].placement_state == ASM_PLACE_SEAM && !memcmp(charts,original,sizeof original);
        x.na = 18; x.n_observations = 0;
        /* A genuinely unsupported or geometrically conflicting proposal
         * must still fail, leaving both original chart structures exact. */
        if (!test) ok = ok && !apl_place_and_audit(&x,1,0,1100,0,NULL,2) && !memcmp(charts,original,sizeof original);
        for (int v = 0; v < 18; v++) xyz[1][3*v+2] = 1000;
        int rejected = !apl_place_and_audit(&x,1,0,1000,0,NULL,2);
        ok = ok && rejected && !memcmp(charts,original,sizeof original);
        if (test) {
            /* A real boundary has matching source vertices, but a third
             * placed chart makes its proposed layout geometrically
             * inconsistent. A separate, clear hop fits eighteen model
             * anchors exactly. It cannot erase that boundary conflict. */
            for (int v = 0; v < 18; v++) {
                xyz[1][3*v+2] = 0;
                uv[2][2*v+1] -= 14;
            }
            apl_rebuild_grids(&x,&grid_arena);
            x.source_geometry_conflict = 0;
            int seam_refused = !apl_place_and_audit(&x,1,0,1000,0,NULL,2);
            ok = ok && seam_refused && x.source_geometry_conflict && !memcmp(charts,original,sizeof original);
            x.cluster_kind = 0; selected[0] = 0;
            for (int v = 0; v < 18; v++) {
                anchors[v].k = anchors[v].measured_order = 1;
                anchors[v].tx = 2000+sign*uv[1][2*v]; anchors[v].ty = uv[1][2*v+1];
                anchors[v].Tx = 1000; anchors[v].observation = (size_t)v;
            }
            x.n_observations = 18;
            int hop_refused = !apl_place_and_audit(&x,1,0,2000,0,NULL,2);
            /* the ledger's reason must name the gate that refused, not a generic "deferred" */
            ok = ok && hop_refused && x.last_silent && !memcmp(charts,original,sizeof original) &&
                 x.last_reason && !strcmp(x.last_reason, "source_geometry_conflict");
            /* Control: absent conflicting SOURCE evidence, the same exact
             * model fit passes all the ordinary geometric audits. */
            x.source_geometry_conflict = 0;
            int clear_hop = apl_place_and_audit(&x,1,0,2000,0,NULL,2);
            ok = ok && clear_hop && proposal[0].placement_state == ASM_PLACE_LAYER &&
                 proposal[0].placement_support == 18 && x.consensus.observations == 18 &&
                 x.consensus.model_inliers == 18 && !memcmp(charts,original,sizeof original);
            /* Eighteen exact selected anchors still fit with eighteen
             * contrary, unselected observations. Alternate wrap hypotheses
             * for the first eighteen rays must not inflate their support. */
            for (int v = 0; v < 18; v++) {
                anchors[18+v] = anchors[v]; anchors[18+v].cluster = 1;
                anchors[18+v].k = 2; anchors[18+v].Tx = 2000; anchors[18+v].tx += 1000;
                anchors[36+v] = anchors[v]; anchors[36+v].cluster = 2;
                anchors[36+v].observation = 18+(size_t)v; anchors[36+v].tx += 1000;
            }
            x.na = 54; x.n_observations = 36;
            int contrary_refused = !apl_place_and_audit(&x,1,0,2000,0,NULL,2);
            ok = ok && contrary_refused && x.last_silent && x.consensus.observations == 36 &&
                 x.consensus.model_inliers == 18 && !memcmp(charts,original,sizeof original) &&
                 x.last_reason && !strcmp(x.last_reason, "consensus_minority");
            size_t cells = AsmConflictGrid_cells(x.cgrid);
            x.trial = 0;
            int direct_refused = !apl_place_and_audit(&x,1,0,2000,0,NULL,2);
            ok = ok && direct_refused && !memcmp(charts,original,sizeof original) &&
                 AsmConflictGrid_cells(x.cgrid) == cells;
            /* Check rollback even for grid cells that already existed. */
            charts[1].pose_x = 2000;
            ok = ok && aprg_probe(x.rgrid,&run,&query,1,24,&touched) == 0 && !touched;
            charts[1] = original[1]; x.trial = 1;
            x.na = 53; x.n_observations = 35;
            int majority_passes = apl_place_and_audit(&x,1,0,2000,0,NULL,2);
            ok = ok && majority_passes && x.consensus.observations == 35 && x.consensus.model_inliers == 18 &&
                 !memcmp(charts,original,sizeof original);
            AsmChart first_pose = proposal[0];
            /* If every ray permits either wrap, both poses remain valid and
             * distinct for the existing ambiguity check; support is not a
             * winner score. An exactly tied model above stays deferred. */
            x.na = 36; x.n_observations = 18; x.cluster = 1;
            for (int v = 0; v < 18; v++) { weights[v] = 0; weights[18+v] = 1; }
            int alternate_passes = apl_place_and_audit(&x,1,0,3000,0,NULL,2);
            ok = ok && alternate_passes && x.consensus.observations == 18 && x.consensus.model_inliers == 18 &&
                 !apl_same_pose(&first_pose,&proposal[0],8) && !memcmp(charts,original,sizeof original);
            if (!contrary_refused || !direct_refused || !majority_passes || !alternate_passes)
                fprintf(stderr,"  asm_place selftest detail: contrary %d direct %d majority %d ambiguous alternative %d\n",
                        contrary_refused,direct_refused,majority_passes,alternate_passes);
            if (!seam_refused || !hop_refused || !clear_hop)
                fprintf(stderr,"  asm_place selftest detail: source conflict %d hop deferred %d clear control %d\n",seam_refused,hop_refused,clear_hop);
        }
        if (!ok) { fprintf(stderr,"  asm_place selftest FAIL: supported seam incumbent mirror %d test %d, closed %d, accepted %d, unsafe rejected %d\n",mirror,test,closed,accepted,rejected); fails++; }
        AsmConflictGrid_dispose(x.cgrid); Arena_dispose(&grid_arena);
        if (x.proposal_uv) Arena_dispose(&x.proposal_uv);
        Arena_dispose(&run.arena);
    }
    return fails;
}

static int apl_stopped_repair_selftest(void)
{
    int fails = 0;
    for (int mirror = 0; mirror < 2; mirror++) for (int total_bound = 0; total_bound < 2; total_bound++) {
        AsmRun run = {0}; run.arena = Arena_new(); run.layer_d_global = 12; run.continuity_ready = 1;
        AsmChart charts[4] = {{0}}, original[4], proposal[3]; AsmRelation rels[4] = {{0}}; AsmCorr corr[24] = {{0}};
        float uv[4][12], xyz[4][18], normals[18]; uint8_t boundary[6];
        int32_t faces[] = {0,1,3,0,3,2,2,3,5,2,5,4};
        for (int i = 0; i < 4; i++) {
            for (int v = 0; v < 6; v++) {
                uv[i][2*v] = (float)(2*(v%2)-1); uv[i][2*v+1] = (float)(2*(v/2)-2);
                xyz[i][3*v] = uv[i][2*v]; xyz[i][3*v+1] = uv[i][2*v+1]; xyz[i][3*v+2] = 0;
                normals[3*v] = normals[3*v+1] = 0; normals[3*v+2] = 1; boundary[v] = 1;
            }
            AsmChart *c = &charts[i]; c->id = i; c->component = i ? 1 : 0; c->nv = 6; c->nf = 4; c->area3d = 8;
            c->uv = uv[i]; c->xyz = xyz[i]; c->nrm = normals; c->faces = faces; c->boundary = boundary;
            c->pose_theta = .4; c->flags = mirror ? ASM_CHART_MIRROR : 0;
        }
        charts[0].placement_state = ASM_PLACE_ROOT;
        for (int i = 0; i < 4; i++) {
            AsmRelation *r = &rels[i]; r->a = i == 1 ? 1 : 0; r->b = i == 0 ? 1 : i == 3 ? 3 : 2;
            r->corr_first = 6*i; r->corr_count = r->n_corr = 6; r->continuity = ASM_CONT_SOURCE;
            r->flags = i == 2 ? ASM_REL_WEAK : 0;
            for (int v = 0; v < 6; v++) {
                AsmCorr *c = &corr[6*i+v]; c->va = c->vb = v; c->valid = 3;
                c->gap_a[0] = c->gap_b[0] = i == 2 ? 10 : 0;
                if (i == 3) {
                    double angle = total_bound ? .26 : .3, u = uv[0][2*v], w = uv[0][2*v+1];
                    c->gap_a[0] = c->gap_b[0] = (float)((cos(angle)-1)*u-sin(angle)*w);
                    c->gap_a[1] = c->gap_b[1] = (float)(sin(angle)*u+(cos(angle)-1)*w);
                }
            }
        }
        run.charts = charts; run.n_charts = 4; run.rels = rels; run.n_rels = 4; run.corr = corr; run.n_corr = 24;
        memcpy(original,charts,sizeof original);
        int32_t groups[] = {0,1,1,1}; uint8_t selected[] = {1,1,1,1}; double area[] = {8,24};
        AsmPlaceOpts opts; AsmPlace_default_opts(&opts); AsmConflictOpts copts; AsmConflict_default_opts(&copts);
        AplCtx x = {0}; x.run = &run; x.o = &opts; x.copts = &copts; x.orig_comp = groups; x.primary = 0;
        x.trial = 1; x.proposal = proposal; x.cluster_kind = 1; x.seam_selected = selected; x.comp_area = area;
        Arena_T grid_arena = Arena_new(); x.cgrid = AsmConflictGrid_new(grid_arena,&copts,&run); x.rgrid = aprg_new(grid_arena,24);
        AsmConflictGrid_insert(x.cgrid,&charts[0]); aprg_insert_chart(x.rgrid,&charts[0]);
        /* A full step for the last chart exceeds either the single-step or
         * cumulative rotation limit. Backtracking can retain a useful fit;
         * every chart must stay inside the original bounds and every SOURCE
         * seam and original geometry pair must survive the final audit. */
        int solved = AsmContinuity_close(&run,original+1,3,selected,NULL,0);
        int within = 1, seams = 1;
        for (int i = 1; i < 4; i++) if (hypot(charts[i].pose_x-original[i].pose_x,charts[i].pose_y-original[i].pose_y) > 256 ||
                                       fabs(charts[i].pose_theta-original[i].pose_theta) > .25) within = 0;
        for (int i = 0; i < 4; i++) { double rms; size_t support; if (!AsmContinuity_measure(&run,&rels[i],&rms,&support)) seams = 0; }
        int ok = seams && within;
        memcpy(charts,original,sizeof original);
        int polished = 0;
        int accepted = apl_place_and_audit_once(&x,1,0,0,0,NULL,2,0,&polished);
        ok = ok && accepted && !polished && !memcmp(charts,original,sizeof original);
        if (accepted) ok = ok && x.nproposal == 3 && proposal[0].placement_state == ASM_PLACE_SEAM &&
                           proposal[0].placement_support >= 6;
        if (accepted) for (int i=0; i<3; i++) ok = ok &&
            hypot(proposal[i].pose_x-original[i+1].pose_x,proposal[i].pose_y-original[i+1].pose_y)<=256 &&
            fabs(proposal[i].pose_theta-original[i+1].pose_theta)<=.25;
        for (int i = 1; i < 4; i++) for (int v = 0; v < 6; v++) xyz[i][3*v+2] = 1000;
        int rejected = !apl_place_and_audit_once(&x,1,0,0,0,NULL,2,0,&polished);
        ok = ok && rejected && !memcmp(charts,original,sizeof original);
        if (!ok) { fprintf(stderr,"  asm_place selftest FAIL: bounded repair mirror %d total bound %d, solved %d, seams %d, within %d, accepted %d, unsafe rejected %d\n",
                          mirror,total_bound,solved,seams,within,accepted,rejected); fails++; }
        AsmConflictGrid_dispose(x.cgrid); Arena_dispose(&grid_arena);
        if (x.proposal_uv) Arena_dispose(&x.proposal_uv);
        Arena_dispose(&run.arena);
    }
    return fails;
}

static int apl_reconsider_selftest(void)
{
    int fails = 0;
    for (int mirror = 0; mirror < 2; mirror++) for (int test = 0; test < 8; test++) {
        AsmRun run = {0}; run.arena = Arena_new(); run.layer_d_global = 12;
        float xyz[3][243], uv[3][162], original_uv[3][162], normals[243];
        int32_t faces[384]; uint8_t boundary[81];
        for (int j = 0; j < 9; j++) for (int u = 0; u < 9; u++) {
            int v = j*9+u; boundary[v] = j == 0 || j == 8 || u == 0 || u == 8;
            normals[3*v] = normals[3*v+1] = 0; normals[3*v+2] = 1;
            for (int side = 0; side < 3; side++) {
                xyz[side][3*v] = (float)(j*16);
                xyz[side][3*v+1] = (float)(u*5+((side == 1 || (side == 2 && test == 1)) ? 44 : 0));
                if (side == 2 && test >= 6) xyz[side][3*v+1] = (float)(u*5+88);
                xyz[side][3*v+2] = side == 2 && test == 1 ? 12.0f : 0;
                uv[side][2*v] = (float)(u*5); uv[side][2*v+1] = (float)(j*16);
            }
            if (j < 8 && u < 8) {
                int32_t f[6] = {v,v+1,v+10,v,v+10,v+9}; memcpy(faces+(j*8+u)*6,f,sizeof f);
            }
        }
        double sign = mirror ? -1.0 : 1.0;
        int nch = test == 1 || test == 2 || test >= 5 ? 3 : 2;
        for (int i = 0; i < nch; i++) {
            AsmChart c = {0}; c.id = i; c.component = 0; c.placed = 1;
            c.nv = 81; c.nf = 128; c.xyz = xyz[i]; c.uv = uv[i]; c.faces = faces; c.nrm = normals; c.boundary = boundary;
            c.area3d = 5120; c.flags = mirror ? ASM_CHART_MIRROR : 0;
            c.pose_x = sign*(i == 1 ? 1244 : (i == 2 ? (test == 1 ? 44 : 200) : 0));
            c.placement_state = i == 1 ? (test == 4 ? ASM_PLACE_SEAM : ASM_PLACE_LAYER) : ASM_PLACE_ROOT;
            if (test == 5 && i == 2) c.placement_state = ASM_PLACE_LAYER;
            if (test >= 6 && i == 2) { c.pose_x = sign*1288; c.placement_state = test == 6 ? ASM_PLACE_SEAM : ASM_PLACE_ROOT; }
            c.placement_group = i;
            AsmRun_push_chart(&run,&c);
        }
        for (int i = 0; i < (test == 2 || test >= 6 ? 2 : 1); i++) {
            AsmRelation r = {0}; r.a = i ? 2 : 0; r.b = 1; r.flags = ASM_REL_WEAK; r.rms = 0.1;
            if (i && test >= 6) { r.a = 1; r.b = 2; r.flags = 0; }
            r.corr_first = (int32_t)run.n_corr; r.corr_count = r.n_corr = 9;
            for (int j = 0; j < 9; j++) { AsmCorr c = {j*9+8,j*9}; AsmRun_push_corr(&run,&c); }
            AsmRun_push_rel(&run,&r);
        }
        AsmContinuity_prepare(&run);
        AsmChart saved[3]; memcpy(saved,run.charts,run.n_charts*sizeof *saved);
        AsmRelation rels[2]; memcpy(rels,run.rels,run.n_rels*sizeof *rels);
        memcpy(original_uv,uv,sizeof uv);
        AsmPlaceOpts opts; AsmPlace_default_opts(&opts);
        AsmConflictOpts copts; AsmConflict_default_opts(&copts);
        AsmPlaceStats stats = {0}; AsmChart proposal[3]; AplCompDiag diag[3] = {{0}};
        AplCtx x = {0}; x.run = &run; x.o = &opts; x.copts = &copts; x.st = &stats; x.primary = 0; x.trial = 1; x.proposal = proposal;
        int32_t groups[3] = {0,1,2}; uint8_t placed[3] = {1,1,1}, selected[2] = {0}, curv[3] = {1,1,1};
        double area[3] = {5120,5120,5120}, cu[3], cv[3], weights[18]; AplAnchor anchors[18];
        if (test == 5) { groups[2] = 1; area[1] += area[2]; }
        x.orig_comp = groups; x.comp_placed = placed; x.comp_area = area; x.seam_selected = selected; x.diag = diag;
        x.cu = cu; x.cv = cv; x.curv_ok = curv; x.anchors = anchors; x.wts = weights; x.hcap = 18; x.cand_grad = 1;
        for (int i = 0; i < nch; i++) apl_uv_centroid(&run.charts[i],&cu[i],&cv[i]);
        Arena_T grid_arena = Arena_new(); apl_rebuild_grids(&x,&grid_arena);
        int accepted = apl_reconsider_hop(&x,1,test == 3 ? ves_clock_sec()-1 : 0);
        int ok = accepted == (test == 0 || test == 6) && !x.excluded_samples && placed[1] == 1 && !memcmp(uv,original_uv,sizeof uv);
        ok = ok && x.orig_comp == groups && groups[0] == 0 && groups[1] == 1 && groups[2] == (test == 5 ? 1 : 2) &&
             area[1] == (test == 5 ? 10240 : 5120) && area[2] == 5120;
        if (accepted) {
            double rms; size_t support;
            ok = ok && fabs(run.charts[1].pose_x-sign*44) < 1e-5 &&
                 run.charts[1].placement_state == ASM_PLACE_SEAM && curv[1] == 0 &&
                 AsmContinuity_measure(&run,&run.rels[0],&rms,&support) && support == 9 && rms < 1e-5 &&
                 (run.rels[0].continuity & ASM_CONT_VERIFIED) && !memcmp(&run.charts[0],&saved[0],sizeof *saved);
            if (test == 6) ok = ok && fabs(run.charts[2].pose_x-sign*88) < 1e-5 && curv[2] == 0 &&
                run.charts[2].placement_group == 2 && run.charts[2].placement_state == ASM_PLACE_SEAM &&
                AsmContinuity_measure(&run,&run.rels[1],&rms,&support) && rms < 1e-5 && support == 9;
            apl_rebuild_grids(&x,&grid_arena);
            const AsmChart *old = &saved[1]; size_t touched;
            aprg_probe(x.rgrid,&run,&old,1,24,&touched);
            ok = ok && touched == 0; /* no stale samples remain at the old hop */
        } else {
            ok = ok && !memcmp(run.charts,saved,run.n_charts*sizeof *saved) &&
                 !memcmp(run.rels,rels,run.n_rels*sizeof *rels) && curv[1] == 1;
            if (test == 3) ok = ok && stats.budget_exhausted;
        }
        if (!ok) { fprintf(stderr,"  asm_place selftest FAIL: reconsider hop test %d mirror %d, accepted %d\n",test,mirror,accepted); fails++; }
        AsmConflictGrid_dispose(x.cgrid); Arena_dispose(&grid_arena);
        if (x.proposal_uv) Arena_dispose(&x.proposal_uv);
        Arena_dispose(&run.arena);
    }
    return fails;
}

static int apl_uv_transaction_selftest(double amplitude)
{
    AsmRun run = {0}; run.arena = Arena_new();
    float xyz[2][243], uv[2][162], original[2][162]; int32_t faces[384]; uint8_t boundary[81];
    for (int j = 0; j < 9; j++) for (int u = 0; u < 9; u++) {
        int v = j*9+u; boundary[v] = j == 0 || j == 8 || u == 0 || u == 8;
        for (int side = 0; side < 2; side++) {
            xyz[side][v*3] = (float)(j*16); xyz[side][v*3+1] = (float)(u*5+side*44); xyz[side][v*3+2] = 0;
            uv[side][v*2] = (float)(u*5); uv[side][v*2+1] = (float)(j*16);
        }
        if (j < 8 && u < 8) {
            int32_t f[6] = {v,v+1,v+10,v,v+10,v+9}; memcpy(faces+(j*8+u)*6,f,sizeof f);
        }
    }
    for (int side = 0; side < 2; side++) {
        AsmChart c = {0}; c.id = c.component = side; c.nv = 81; c.nf = 128;
        c.xyz = xyz[side]; c.uv = uv[side]; c.faces = faces; c.boundary = boundary; c.area3d = 5120;
        c.pose_x = side*44; c.placement_state = side ? ASM_PLACE_NONE : ASM_PLACE_ROOT;
        AsmRun_push_chart(&run,&c);
    }
    AsmRelation r = {0}; r.a = 0; r.b = 1; r.corr_count = r.n_corr = 9; r.rms = 0.1; r.flags = ASM_REL_SWITCHED;
    for (int j = 0; j < 9; j++) { AsmCorr c = {j*9+8,j*9}; AsmRun_push_corr(&run,&c); }
    AsmRun_push_rel(&run,&r); AsmContinuity_prepare(&run);
    for (int j = 0; j < 9; j++) for (int u = 0; u < 9; u++)
        uv[1][(j*9+u)*2] += (float)(amplitude*cos(j*APL_PI/4)*(1-u/8.0));
    double initial_rms; size_t initial_support;
    int initially_passes = AsmContinuity_measure(&run,&run.rels[0],&initial_rms,&initial_support);
    memcpy(original,uv,sizeof uv);
    AsmChart saved[2] = {run.charts[0],run.charts[1]}, proposal[1], winner[1];
    AsmPlaceOpts opts; AsmPlace_default_opts(&opts);
    AplCtx x = {0}; x.run = &run; x.o = &opts; x.primary = 0; x.trial = 1; x.proposal = proposal; x.cluster_kind = 1;
    int32_t groups[2] = {0,1}; uint8_t placed[2] = {1,0}, selected[1] = {1}, curv[2] = {1,1};
    AplAnchor anchors[9]; double weights[9], cu[2], cv[2], old_cu[2], old_cv[2];
    x.orig_comp = groups; x.comp_placed = placed; x.seam_selected = selected;
    x.anchors = anchors; x.wts = weights; x.hcap = 9; x.cu = cu; x.cv = cv; x.curv_ok = curv;
    for (int side = 0; side < 2; side++) apl_uv_centroid(&run.charts[side],&cu[side],&cv[side]);
    memcpy(old_cu,cu,sizeof cu); memcpy(old_cv,cv,sizeof cv);
    int ok = initial_support == 9 && initially_passes == (amplitude <= 8) && apl_collect_anchors(&x,1) == 9; x.na = 9;
    for (int i = 0; i < 9; i++) weights[i] = 1;
    ok = ok && apl_place_and_audit(&x,1,0,0,0,NULL,2) && x.uv_trial && x.nproposal == 1;
    ok = ok && !memcmp(run.charts,saved,sizeof saved) && !memcmp(uv,original,sizeof uv) &&
         !memcmp(cu,old_cu,sizeof cu) && !memcmp(cv,old_cv,sizeof cv) && curv[1] == 1;
    Arena_T winner_uv = NULL;
    if (ok) {
        apl_copy_candidate(winner,proposal,1,&run,&winner_uv);
        float winning_uv[162]; memcpy(winning_uv,winner[0].uv,sizeof winning_uv);
        if (initially_passes) {
            /* Put the axis limit between the supported incumbent and its
             * current UV winner. That winner must be forbidden while the
             * incumbent remains available. A different feasible correction
             * may win; a new solver direction need not return identical UVs. */
            int refined = 0;
            ok = apl_place_and_audit_once(&x,1,0,0,0,NULL,2,0,&refined) && !refined;
            AsmChart base = proposal[0];
            double u, v, bu, bv, pu, pv, axial[2] = {0};
            apl_uv_centroid(&base,&u,&v); apl_pose(&base,u,v,&bu,&bv);
            apl_uv_centroid(&winner[0],&u,&v); apl_pose(&winner[0],u,v,&pu,&pv);
            double change = (pv-0.1*pu)-(bv-0.1*bu);
            ok = ok && fabs(change) > 1e-5;
            x.vm_ok = x.vm_pin = 1; x.vm_c = 0.1; x.s_axis = axial;
            x.vm_a = bv-0.1*bu-copysign(APL_V_AXIS_MAX_SHIFT-0.25*fabs(change),change);
            ok = ok && fabs(x.vm_a+0.1*pu-pv) > APL_V_AXIS_MAX_SHIFT &&
                 apl_place_and_audit_once(&x,1,0,0,0,NULL,2,0,&refined) && !refined &&
                 apl_same_pose(&base,&proposal[0],1e-6) && !memcmp(proposal[0].uv,original[1],sizeof original[1]);
            ok = ok && apl_place_and_audit(&x,1,0,0,0,NULL,2) &&
                 !memcmp(run.charts,saved,sizeof saved) && !memcmp(uv,original,sizeof uv) &&
                 !memcmp(winner[0].uv,winning_uv,sizeof winning_uv);
            if (ok) {
                double au, av; apl_uv_centroid(&proposal[0],&u,&v); apl_pose(&proposal[0],u,v,&au,&av);
                AsmMetricFace metric[128]; AsmMetricStats before, after;
                ok = fabs(x.vm_a+0.1*au-av) <= APL_V_AXIS_MAX_SHIFT && !AsmMetric_prepare(&proposal[0],metric);
                if (ok) {
                    AsmMetric_measure(&proposal[0],metric,original[1],&before);
                    AsmMetric_measure(&proposal[0],metric,proposal[0].uv,&after);
                    ok = AsmMetric_preserved(&before,&after);
                }
            }
            x.vm_ok = x.vm_pin = 0; x.s_axis = NULL;
        }
        /* The next candidate fails after correction. It must neither leak
         * UV changes nor destroy the previously saved winning geometry. */
        x.span_cap = 1;
        ok = ok && !apl_place_and_audit(&x,1,0,0,0,NULL,2) &&
             !memcmp(run.charts,saved,sizeof saved) && !memcmp(uv,original,sizeof uv) &&
             !memcmp(winner[0].uv,winning_uv,sizeof winning_uv);
        apl_commit_chart(&x,&winner[0],uv[1]);
        double rms; size_t support;
        ok = ok && AsmContinuity_measure(&run,&run.rels[0],&rms,&support) && support == 9 && rms < initial_rms &&
             run.charts[1].uv == uv[1] && memcmp(uv[1],original[1],sizeof uv[1]) && curv[1] == 0 &&
             !memcmp(&run.charts[0],&saved[0],sizeof saved[0]) && !memcmp(uv[0],original[0],sizeof uv[0]);
        double mean_u, mean_v; apl_uv_centroid(&run.charts[1],&mean_u,&mean_v);
        ok = ok && cu[1] == mean_u && cv[1] == mean_v;
        for (size_t v = 0; v < 81; v++) {
            /* Multiple accepted corrections share one total UV bound. The
             * opposite edge, outside the boundary band, stays exact. */
            ok = ok && hypot((double)uv[1][2*v]-original[1][2*v],
                             (double)uv[1][2*v+1]-original[1][2*v+1]) <= 8;
            if (v%9 == 8) ok = ok && !memcmp(uv[1]+2*v,original[1]+2*v,2*sizeof(float));
        }
    }
    if (!ok) fprintf(stderr,"  asm_place selftest FAIL: local UV trial rollback, winner lifetime or commit (amplitude %.1f)\n",amplitude);
    if (x.proposal_uv) Arena_dispose(&x.proposal_uv);
    if (winner_uv) Arena_dispose(&winner_uv);
    Arena_dispose(&run.arena); return !ok;
}

static int apl_rim_pair_selftest(void)
{
    const double reach = 24;
    const double radii[] = {0,1,7,16,nextafter(24,0),24,nextafter(24,INFINITY),32};
    size_t counts[3] = {0};
    for (int i = 0; i < 65536; i++) {
        double angle = i*0.6180339887498948, r = radii[i%8], tol = 7+(i%5);
        double dx = r*cos(angle), dy = r*sin(angle);
        if (i%7 == 0) { dx = r; dy = 0; }
        float p[3], q[3];
        double physical = i%3 == 0 ? r+tol : (i%3 == 1 ? r : 48.0*(i%101)/100.0);
        for (int k = 0; k < 3; k++) p[k] = i%4 == 0 ? 0 : (float)((k+1)*4096+i%17);
        q[0] = p[0]+(float)(physical*cos(angle));
        q[1] = p[1]+(float)(physical*sin(angle)); q[2] = p[2];
        if (i%11 == 0) { q[0] = p[0]+(float)physical; q[1] = p[1]; }
        if (i%13 == 0) q[0] = nextafterf(q[0],i%2 ? INFINITY : -INFINITY);
        double duv = hypot(dx,dy);
        double d3 = sqrt((q[0]-p[0])*(q[0]-p[0])+(q[1]-p[1])*(q[1]-p[1])+(q[2]-p[2])*(q[2]-p[2]));
        int reference = duv > reach ? -1 : d3-duv > tol;
        int actual = aprg_pair_test(dx,dy,p,q,reach,reach*reach*(1+8*DBL_EPSILON),reach*reach*(1-64*DBL_EPSILON),tol);
        if (reference != actual) {
            fprintf(stderr,"  asm_place selftest FAIL: rim pair %d, UV %.17g physical %.17g tol %.17g: %d vs %d\n",i,duv,d3,tol,actual,reference);
            return 1;
        }
        counts[reference+1]++;
    }
    return !counts[0] || !counts[1] || !counts[2];
}

static int apl_rotated_cluster_selftest(void)
{
    int fails = 0;
    for (int sign = -1; sign <= 1; sign += 2) for (int shifted = 0; shifted < 2; shifted++)
    for (int distinct_rotation = 0; distinct_rotation < 2; distinct_rotation++) {
        AsmRun run = {0}; run.arena = Arena_new();
        AsmPlaceOpts opts; AsmPlace_default_opts(&opts);
        AplAnchor anchors[128] = {{0}}, original[128];
        for (int branch = 0; branch < 2; branch++) for (int v = 0; v < 64; v++) {
            AplAnchor *a = &anchors[64*branch+v];
            double u = 1600.0*(v/2)/31, w = v%2 ? 80 : -80;
            double th = sign*.25*(branch && distinct_rotation ? -1 : 1);
            a->sx = u+(shifted ? 1e6 : 0); a->sy = w-(shifted ? 2e6 : 0);
            a->tx = cos(th)*u-sin(th)*w+1000*(branch+1);
            a->ty = sin(th)*u+cos(th)*w;
            a->w = branch ? .5 : 1; a->k = branch+1; a->measured_order = 1;
            a->Tx = 1000*(branch+1); a->cluster = -1;
        }
        memcpy(original,anchors,sizeof original);
        AplCtx x = {0}; x.run = &run; x.o = &opts; x.anchors = anchors; x.cand_grad = 1;
        double weight[APL_MAX_CLUSTERS_TOTAL], u[APL_MAX_CLUSTERS_TOTAL], v[APL_MAX_CLUSTERS_TOTAL];
        int kind[APL_MAX_CLUSTERS_TOTAL];
        int count = apl_cluster_anchors(&x,128,weight,u,v,kind);
        int ok = distinct_rotation ? count >= 2 : count == 2;
        size_t second = 0;
        for (int i = 0; i < 128; i++) {
            if (i < 64) ok = ok && anchors[i].cluster == 0;
            else {
                ok = ok && anchors[i].cluster != 0;
                if (anchors[i].cluster >= 0) second++;
                if (!distinct_rotation) ok = ok && anchors[i].cluster == 1;
            }
            AplAnchor saved = anchors[i]; saved.cluster = original[i].cluster;
            ok = ok && !memcmp(&saved,&original[i],sizeof saved);
        }
        /* A genuinely different rotation at another wrap remains a
         * separate supported alternative. No pose or witness is edited. */
        ok = ok && second >= 12;
        if (!ok) {
            fprintf(stderr,"  asm_place selftest FAIL: rotated clustering sign %d shifted %d distinct %d: %d clusters, second support %zu\n",
                    sign,shifted,distinct_rotation,count,second);
            fails++;
        }
        Arena_dispose(&run.arena);
    }
    for (int sign = -1; sign <= 1; sign += 2) {
        /* An incoherent alternative must retain its original membership.
         * Rotating its offsets by the good component's fitted angle would
         * split the observations into artificial, smaller hypotheses. */
        AsmRun run = {0}; run.arena = Arena_new();
        AsmPlaceOpts opts; AsmPlace_default_opts(&opts);
        AplAnchor anchors[96] = {{0}}, original[96];
        for (int i = 0; i < 96; i++) {
            AplAnchor *a = anchors+i;
            double u = 1600.0*((i < 64 ? i : 2*(i-64))/2)/31, v = i%2 ? 80 : -80;
            a->sx = u; a->sy = v; a->w = i < 64 ? 1 : .5;
            a->k = i < 64 ? 1 : 2; a->Tx = 1000*a->k; a->measured_order = 1; a->cluster = -1;
            a->tx = i < 64 ? cos(sign*.25)*u-sin(sign*.25)*v+1000 : u+2500+(i%2 ? 20 : -20);
            a->ty = i < 64 ? sin(sign*.25)*u+cos(sign*.25)*v : v+(i%2 ? -60 : 60);
        }
        memcpy(original,anchors,sizeof original);
        AplCtx x = {0}; x.run = &run; x.o = &opts; x.anchors = anchors; x.cand_grad = 1;
        double w[APL_MAX_CLUSTERS_TOTAL], u[APL_MAX_CLUSTERS_TOTAL], v[APL_MAX_CLUSTERS_TOTAL]; int kind[APL_MAX_CLUSTERS_TOTAL];
        int n = apl_cluster_anchors(&x,96,w,u,v,kind), ok = n == 2;
        for (int i = 0; i < 96; i++) {
            ok = ok && anchors[i].cluster == (i < 64 ? 0 : 1);
            AplAnchor saved = anchors[i]; saved.cluster = original[i].cluster;
            ok = ok && !memcmp(&saved,&original[i],sizeof saved);
        }
        if (!ok) { fprintf(stderr,"  asm_place selftest FAIL: incoherent alternative regrouped, sign %d (%d clusters)\n",sign,n); fails++; }
        Arena_dispose(&run.arena);
    }
    return fails;
}

static int apl_v_model_update_selftest(void)
{
    AsmRun run = {0}; run.arena = Arena_new();
    double cu[13] = {0}, cv[13] = {0}, axis[13]; int32_t original[13];
    float uv[2] = {0,0};
    for (int i = 0; i < 13; i++) {
        AsmChart c = {0}; c.nv = 1; c.uv = uv; c.area3d = 1;
        c.component = 0; original[i] = i < 8 ? 0 : 1;
        axis[i] = ((i/2)%2) ? 100 : -100;
        c.pose_x = i < 8 ? ((i%2) ? 500 : -500) : ((i%2) ? 2500 : 1500);
        c.pose_y = axis[i]+(i < 8 ? 0 : ((i%2) ? 45 : 35));
        c.placement_state = i < 8 ? ASM_PLACE_ROOT : ASM_PLACE_LAYER;
        if (i == 12) {
            /* Even an enormous pending candidate in the primary's trial
             * frame cannot contribute to a committed height calibration. */
            c.area3d = 1e12; c.pose_x = 10000; c.pose_y = -10000;
            c.placement_state = ASM_PLACE_NONE; axis[i] = 0;
        }
        AsmRun_push_chart(&run,&c);
    }
    AplCtx x = {0}; x.run = &run; x.primary = 0; x.orig_comp = original;
    x.cu = cu; x.cv = cv; x.s_axis = axis;
    int ok = !apl_v_model(&x) && x.vm_ok && fabs(x.vm_a) < 1e-9 && fabs(x.vm_b-1) < 1e-9 && fabs(x.vm_c) < 1e-9;
    x.vm_pin = 1; x.vm_p90 = 20;
    AsmChart before[13]; memcpy(before,run.charts,sizeof before);
    ok = ok && apl_v_model_update(&x) && fabs(x.vm_a-100.0/41) < 1e-8 &&
         fabs(x.vm_b-1) < 1e-9 && fabs(x.vm_c-67.0/4100) < 1e-9 && x.vm_p90 == 20 &&
         !memcmp(before,run.charts,sizeof before) && !apl_v_model_update(&x);
    /* An unphysical axial scale cannot replace the incumbent calibration. */
    for (int i = 8; i < 12; i++) run.charts[i].pose_y += axis[i];
    double a = x.vm_a, b = x.vm_b, c = x.vm_c;
    ok = ok && !apl_v_model_update(&x) && x.vm_a == a && x.vm_b == b && x.vm_c == c;
    for (int i = 8; i < 12; i++) run.charts[i].pose_y -= axis[i];
    /* Reducing area-weighted error cannot invalidate an already accepted
     * small component at the old mean-height limit. */
    x.vm_a = x.vm_c = 0; x.vm_b = 1;
    run.charts[12].area3d = .001; run.charts[12].pose_y = -99;
    run.charts[12].placement_state = ASM_PLACE_SEAM; original[12] = 2;
    ok = ok && !apl_v_model_update(&x) && x.vm_a == 0 && x.vm_b == 1 && x.vm_c == 0 && x.vm_p90 == 20;
    if (!ok) fprintf(stderr,"  asm_place selftest FAIL: committed height recalibration and incumbent bounds\n");
    Arena_dispose(&run.arena);
    return !ok;
}

static int apl_self_contact_audit_selftest(void)
{
    int fails = 0;
    for (int mirror = 0; mirror < 2; mirror++) for(int deep=0;deep<2;deep++) {
        float scale=deep?128.f:1.f;
        AsmRun run = {0}; run.arena = Arena_new(); run.layer_d_global = 10;
        float uv[12] = {0,0,2,0,0,2, 1.05f,1.05f,3.05f,1.05f,1.05f,3.05f};
        for(int k=0;k<12;k++)uv[k]*=scale;
        float trial[12], smaller[12], xyz[18];
        int32_t faces[6] = {0,1,2,3,4,5};
        for (int v = 0; v < 6; v++) { xyz[3*v] = uv[2*v]; xyz[3*v+1] = uv[2*v+1]; xyz[3*v+2] = 0; }
        AsmChart chart = {0}; chart.nv = 6; chart.nf = 2; chart.uv = uv; chart.xyz = xyz; chart.faces = faces;
        chart.flags = mirror ? ASM_CHART_MIRROR : 0;
        AsmRun_push_chart(&run,&chart);
        AsmChart saved = run.charts[0];
        run.charts[0].pose_theta = .41; run.charts[0].pose_x = 3000; run.charts[0].pose_y = -4000;
        uint8_t blocked[1] = {0}, vertices[6] = {0}; uint8_t *vertex_rows[1] = {vertices};
        AsmContinuityAuditCache cache = {0};
        cache.blocked_charts = blocked; cache.blocked_charts_count = 1;
        cache.blocked_vertices = vertex_rows; cache.blocked_vertices_count = 1;
        int ok = apl_uv_audit(&run,&saved,1,&cache);
        memcpy(trial,uv,sizeof uv);
        for (int j = 6; j < 12; j++) trial[j] -= .1f*scale;
        run.charts[0].uv = trial;
        /* The overlapping right triangle has leg .1*scale and inscribed
         * diameter (2-sqrt(2))*.1*scale: below 4, or 7.50 vox when deep. */
        ok = ok && apl_internal_bad(&run,&saved,1,1) == 0 && apl_uv_audit(&run,&saved,1,&cache)==!deep && blocked[0]==deep;
        for (int v = 0; v < 6; v++) ok = ok && vertices[v]==deep;
        run.charts[0].uv = uv;
        ok = ok && apl_uv_audit(&run,&saved,1,&cache) && !blocked[0];
        /* An existing pair may shrink, but it may not grow. */
        saved.uv = trial;
        memcpy(smaller,uv,sizeof uv);
        for (int j = 6; j < 12; j++) smaller[j] -= .075f*scale;
        run.charts[0].uv = smaller;
        AsmContinuityAuditCache existing = {0};
        ok = ok && apl_uv_audit(&run,&saved,1,&existing);
        for (int j = 6; j < 12; j++) smaller[j] = uv[j]-.15f*scale;
        ok = ok && apl_uv_audit(&run,&saved,1,&existing)==!deep;
        if (!ok) { fprintf(stderr,"  asm_place selftest FAIL: triangle self-contact preservation mirror %d\n",mirror); fails++; }
        Arena_dispose(&run.arena);
    }
    return fails;
}

static int apl_single_chart_uv_audit_selftest(void)
{
    int fails = 0;
    for (int mirror = 0; mirror < 2; mirror++) {
        AsmRun run = {0}; run.arena = Arena_new(); run.layer_d_global = 10;
        float xyz[243], uv[162], contracted[162];
        for (int j = 0; j < 9; j++) for (int u = 0; u < 9; u++) {
            int v = 9*j+u;
            xyz[3*v] = uv[2*v] = (float)(4*u);
            xyz[3*v+1] = uv[2*v+1] = (float)(4*j); xyz[3*v+2] = 0;
            /* Positive orientation, at most 2x local distortion, and every
             * vertex stays within the actual eight-voxel repair bound. */
            contracted[2*v] = .5f*uv[2*v]+8; contracted[2*v+1] = uv[2*v+1];
        }
        AsmChart c = {0}; c.nv = 81; c.xyz = xyz; c.uv = uv;
        c.flags = mirror ? ASM_CHART_MIRROR : 0; AsmRun_push_chart(&run,&c);
        AsmChart saved = run.charts[0];
        run.charts[0].pose_theta = .4; run.charts[0].pose_x = 100; run.charts[0].pose_y = -20;
        AsmContinuityAuditCache cache = {0};
        int ok = apl_uv_audit(&run,&saved,1,&cache);
        for (int v = 0; v < 81; v++)
            ok = ok && hypot(contracted[2*v]-uv[2*v],contracted[2*v+1]-uv[2*v+1]) <= 8;
        run.charts[0].uv = contracted;
        size_t bad = apl_internal_bad(&run,&saved,1,1);
        ok = ok && bad > 0 && !apl_uv_audit(&run,&saved,1,&cache);
        run.charts[0].uv = uv;
        ok = ok && apl_uv_audit(&run,&saved,1,&cache);
        if (!ok) { fprintf(stderr,"  asm_place selftest FAIL: single-chart UV contraction mirror %d, bad %zu\n",mirror,bad); fails++; }
        Arena_dispose(&run.arena);
    }
    return fails;
}

static int apl_guarded_uv_audit_selftest(void)
{
    int fails = 0;
    for (int mirror = 0; mirror < 2; mirror++) {
        AsmRun run = {0}; run.arena = Arena_new(); run.layer_d_global = 10;
        float xyz[2][3] = {{0,0,0},{40,0,0}}, uv[2][2] = {{0,0},{40,0}}, trial[2] = {16,0};
        for (int i = 0; i < 2; i++) {
            AsmChart c = {0}; c.id = i; c.nv = 1; c.xyz = xyz[i]; c.uv = uv[i];
            c.flags = mirror ? ASM_CHART_MIRROR : 0; c.pose_theta = .4;
            AsmRun_push_chart(&run,&c);
        }
        AsmChart saved[2] = {run.charts[0],run.charts[1]};
        uint8_t blocked[2] = {0};
        AsmContinuityAuditCache guard = {0}; guard.guard_charts = saved; guard.guard_count = 2;
        guard.blocked_charts = blocked; guard.blocked_charts_count = 2;
        int ok = apl_uv_audit(&run,saved,1,&guard) && guard.baseline_queries == 2 && !guard.baseline_bad;
        /* This is a move across the diameter of a disk centred at8. The
         * fixed neighbour starts40 vox away and enters the24-voxel rim. */
        run.charts[0].uv = trial;
        AsmContinuityAuditCache alone = {0};
        ok = ok && apl_uv_audit(&run,saved,1,&alone) && !apl_uv_audit(&run,saved,1,&guard) && blocked[0] && blocked[1];
        run.charts[0].uv = uv[0];
        ok = ok && apl_uv_audit(&run,saved,1,&guard) && !blocked[0] && !blocked[1];
        guard.guard_count = 0;
        ok = ok && !apl_uv_audit(&run,saved,1,&guard);
        guard.guard_charts = NULL; guard.guard_count = 2;
        ok = ok && !apl_uv_audit(&run,saved,1,&guard);
        if (!ok) { fprintf(stderr,"  asm_place selftest FAIL: fixed neighbour of a UV repair, mirror %d\n",mirror); fails++; }
        Arena_dispose(&run.arena);
    }
    return fails;
}

static int apl_guarded_unsampled_vertex_selftest(void)
{
    int fails = 0;
    for (int mirror = 0; mirror < 2; mirror++) for (int query = 0; query < 2; query++) for (int reference = 0; reference < 2; reference++) {
        AsmRun run = {0}; run.arena = Arena_new(); run.layer_d_global = 10;
        float uv[2][6], xyz[2][9], trial[6]; int32_t face[3] = {0,1,2};
        for (int i = 0; i < 2; i++) {
            int tip = i ? reference : query;
            for (int v = 0; v < 3; v++) {
                float u = v == tip ? (i ? 25.0f : 0.0f) : (i ? -100.0f : 100.0f);
                float w = v == tip ? 0.0f : (v == 2 ? -100.0f : 100.0f);
                uv[i][2*v] = xyz[i][3*v] = u; uv[i][2*v+1] = xyz[i][3*v+1] = w;
                xyz[i][3*v+2] = (float)(40*i);
            }
            AsmChart c = {0}; c.id = i; c.nv = 3; c.nf = 1; c.xyz = xyz[i]; c.uv = uv[i]; c.faces = face;
            c.flags = mirror ? ASM_CHART_MIRROR : 0; c.pose_theta = .4; c.pose_x = 100; c.pose_y = -50;
            AsmRun_push_chart(&run,&c);
        }
        AsmChart saved[2] = {run.charts[0],run.charts[1]};
        uint8_t blocked[2] = {0}, points[2][5] = {{0xa5,0,0,0,0x5a},{0xa5,0,0,0,0x5a}};
        uint8_t *rows[2] = {points[0]+1,points[1]+1};
        AsmContinuityAuditCache guard = {0}; guard.guard_charts = saved; guard.guard_count = 2;
        guard.blocked_charts = blocked; guard.blocked_charts_count = 2;
        guard.blocked_vertices = rows; guard.blocked_vertices_count = 2;
        AsmContinuityUvPair pairs[3], sentinel;
        memset(&sentinel,0xa5,sizeof sentinel); for (int i = 0; i < 3; i++) pairs[i] = sentinel;
        guard.pairs = pairs+1; guard.pairs_capacity = 1;
        int ok = apl_uv_audit(&run,saved,2,&guard) && !guard.baseline_bad;
        ok = ok && !guard.pairs_count;
        /* Only these two tips meet the rim. Moving the query two voxels
         * introduces a collision even when either endpoint's source index
         * is absent from the old even-query/every-fourth-reference sample. */
        memcpy(trial,uv[0],sizeof trial); trial[2*query] += 2;
        run.charts[0].uv = trial;
        int accepted = apl_uv_audit(&run,saved,2,&guard);
        ok = ok && !accepted && blocked[0] && blocked[1] && rows[0][query] && rows[1][reference];
        ok = ok && guard.pairs_count == 1 && pairs[1].chart_a == 0 && pairs[1].vertex_a == query &&
             pairs[1].chart_b == 1 && pairs[1].vertex_b == reference &&
             !memcmp(pairs,&sentinel,sizeof sentinel) && !memcmp(pairs+2,&sentinel,sizeof sentinel);
        for (int i = 0; i < 2; i++) {
            ok = ok && points[i][0] == 0xa5 && points[i][4] == 0x5a;
            for (int v = 0; v < 3; v++) ok = ok && rows[i][v] == (v == (i ? reference : query));
        }
        /* Pair feedback also works without either optional bitmap. A zero
         * capacity must not write the supplied buffer or retain old counts. */
        guard.blocked_charts = NULL; guard.blocked_vertices = NULL;
        ok = ok && !apl_uv_audit(&run,saved,2,&guard) && guard.pairs_count == 1;
        guard.pairs_capacity = 0; pairs[1] = sentinel;
        ok = ok && !apl_uv_audit(&run,saved,2,&guard) && !guard.pairs_count && !memcmp(pairs+1,&sentinel,sizeof sentinel);
        guard.pairs_capacity = 1; guard.blocked_charts = blocked; guard.blocked_vertices = rows;
        run.charts[0].uv = uv[0];
        ok = ok && apl_uv_audit(&run,saved,2,&guard) && !blocked[0] && !blocked[1];
        ok = ok && !guard.pairs_count;
        if (!ok) {
            fprintf(stderr,"  asm_place selftest FAIL: unsampled UV collision, mirror %d query %d reference %d, accepted %d\n",
                    mirror,query,reference,accepted); fails++;
        }
        Arena_dispose(&run.arena);
    }
    return fails;
}

static int apl_guarded_pair_identity_selftest(void)
{
    int fails = 0;
    for (int mirror = 0; mirror < 2; mirror++) {
        AsmRun run = {0}; run.arena = Arena_new(); run.layer_d_global = 10;
        float uv[2][10] = {{0}}, xyz[2][15] = {{0}}, trial[10];
        for (int i = 0; i < 2; i++) {
            for (int v = 0; v < 5; v++) {
                uv[i][2*v] = xyz[i][3*v] = v == 4 ? 25.0f : 0.0f;
                xyz[i][3*v+2] = (float)(40*i);
            }
            AsmChart c = {0}; c.id = i; c.nv = 5; c.xyz = xyz[i]; c.uv = uv[i];
            c.flags = mirror ? ASM_CHART_MIRROR : 0; c.pose_theta = .4;
            c.pose_x = 100; c.pose_y = -50;
            AsmRun_push_chart(&run,&c);
        }
        AsmChart saved[2] = {run.charts[0],run.charts[1]};
        AsmContinuityAuditCache legacy = {0}, guard = {0};
        uint8_t blocked[4] = {0xa5,0xcc,0xcc,0x5a};
        uint8_t points[2][7]; memset(points,0xcc,sizeof points);
        for (int i = 0; i < 2; i++) { points[i][0] = 0xa5; points[i][6] = 0x5a; }
        uint8_t *point_rows[2] = {points[0]+1,points[1]+1};
        guard.guard_charts = saved; guard.guard_count = 2;
        guard.blocked_charts = blocked+1; guard.blocked_charts_count = 2;
        guard.blocked_vertices = point_rows; guard.blocked_vertices_count = 2;
        size_t before = apl_internal_bad(&run,saved,2,0);
        int ok = before == 6 && apl_uv_audit(&run,saved,2,&guard);
        /* All six sampled queries already overlap the other chart. Moving
         * its second reference sample by two voxels adds bad vertex pairs
         * inside that SAME chart pair, without changing any query flag. */
        memcpy(trial,uv[1],sizeof trial); trial[8] = 23;
        run.charts[1].uv = trial;
        size_t after = apl_internal_bad(&run,saved,2,1);
        int old_accept = apl_uv_audit(&run,saved,2,&legacy);
        int guarded_accept = apl_uv_audit(&run,saved,2,&guard);
        ok = ok && after == before && old_accept && !guarded_accept &&
             blocked[0] == 0xa5 && blocked[1] && blocked[2] && blocked[3] == 0x5a;
        /* Both original query indices and the second reference sample's
         * source index must survive the compact grid representation. */
        for (int i = 0; i < 2; i++) for (int v = 0; v < 5; v++)
            ok = ok && point_rows[i][v] == (i ? v == 4 : v < 4);
        run.charts[1].uv = uv[1];
        ok = ok && apl_uv_audit(&run,saved,2,&guard) && !blocked[1] && !blocked[2];
        for (int i = 0; i < 2; i++) {
            ok = ok && points[i][0] == 0xa5 && points[i][6] == 0x5a;
            for (int v = 0; v < 5; v++) ok = ok && !point_rows[i][v];
        }
        if (!ok) {
            fprintf(stderr,"  asm_place selftest FAIL: additional bad reference pairs at already bad queries, mirror %d, bad %zu -> %zu, accepted legacy %d guarded %d\n",
                    mirror,before,after,old_accept,guarded_accept); fails++;
        }
        Arena_dispose(&run.arena);
    }
    return fails;
}

static int apl_polish_bound_selftest(void)
{
    int fails = 0;
    for (int mirror = 0; mirror < 2; mirror++) {
        AsmRun run = {0}; run.arena = Arena_new(); run.layer_d_global = 14;
        float ref_uv[2] = {0,0}, ref_xyz[3] = {0,0,0};
        float uv[34] = {0}, original[34] = {0}, xyz[51] = {0};
        float far_uv[62], far_xyz[93] = {0};
        for (int i = 0; i < 17; i++) xyz[3*i] = 80;
        for (int i = 0; i < 31; i++) { far_uv[2*i] = 1000; far_uv[2*i+1] = 0; }
        for (int i = 0; i < 3; i++) {
            AsmChart c = {0}; c.id = i; c.flags = mirror ? ASM_CHART_MIRROR : 0;
            c.nv = i == 0 ? 1 : i == 1 ? 17 : 31;
            c.uv = i == 0 ? ref_uv : i == 1 ? uv : far_uv;
            c.xyz = i == 0 ? ref_xyz : i == 1 ? xyz : far_xyz;
            c.pose_theta = .43; c.pose_x = 100; c.pose_y = -50;
            AsmRun_push_chart(&run,&c);
        }
        AplRimGrid *grid = aprg_new(run.arena,24);
        aprg_insert_chart(grid,&run.charts[0]);
        AsmChart saved[2] = {run.charts[1],run.charts[2]}; saved[0].uv = original;
        const AsmChart *candidate = &run.charts[1];
        size_t guaranteed, queries, touched;
        int ok = aprg_polish_blocked(grid,&run,saved,1,NULL,.4,&guaranteed,&queries) &&
                 guaranteed == 8 && queries == 9;
        /* Independently probe an entire circle of allowed eight-voxel
         * corrections. Every query remains bad as the bound predicted. */
        for (int j = 0; j < 64 && ok; j++) {
            double angle = 2*APL_PI*j/64;
            for (int i = 0; i < 17; i++) { uv[2*i] = (float)(8*cos(angle)); uv[2*i+1] = (float)(8*sin(angle)); }
            ok = aprg_probe(grid,&run,&candidate,1,24,&touched) == 1 && touched == 9;
        }
        memset(uv,0,sizeof uv);
        /* Untouched queries can later enter the rim and enlarge its
         * denominator. An initial touched fraction alone is insufficient. */
        ok = ok && !aprg_polish_blocked(grid,&run,saved,2,NULL,.4,&guaranteed,&queries) && queries == 25;
        uint8_t excluded[3] = {1,0,0};
        ok = ok && !aprg_polish_blocked(grid,&run,saved,1,excluded,.4,&guaranteed,&queries);
        run.charts[1].nv = saved[0].nv = 13;
        ok = ok && !aprg_polish_blocked(grid,&run,saved,1,NULL,.4,&guaranteed,&queries);
        run.charts[1].nv = saved[0].nv = 17;
        /* A prior eight-voxel move toward the rim can be reversed. The
         * legal endpoint at25 is outside the24-voxel rim; do not prune it. */
        for (int i = 0; i < 17; i++) { original[2*i] = 17; uv[2*i] = 9; }
        ok = ok && !aprg_polish_blocked(grid,&run,saved,1,NULL,.4,&guaranteed,&queries);
        for (int i = 0; i < 17; i++) uv[2*i] = 25;
        ok = ok && aprg_probe(grid,&run,&candidate,1,24,&touched) == 0 && !touched;
        /* The physical-distance boundary also has a feasible endpoint:
         * fifteen minus eight equals the unchanged seven-voxel tolerance. */
        memset(uv,0,sizeof uv); memset(original,0,sizeof original);
        for (int i = 0; i < 17; i++) xyz[3*i] = 15;
        ok = ok && !aprg_polish_blocked(grid,&run,saved,1,NULL,.4,&guaranteed,&queries);
        /* Roundoff in the common rotation may straddle exact equality;
         * translating the original disk by one voxel gives a strict
         * geometric pass at its legal eight-voxel endpoint. */
        for (int i = 0; i < 17; i++) { original[2*i] = 1; uv[2*i] = 9; }
        ok = ok && aprg_probe(grid,&run,&candidate,1,24,&touched) == 0 && touched == 9;
        if (!ok) { fprintf(stderr,"  asm_place selftest FAIL: conservative UV polish rim bound, mirror %d\n",mirror); fails++; }
        Arena_dispose(&run.arena);
    }
    return fails;
}

static int apl_internal_overlap_exchange_selftest(void)
{
    int fails = 0;
    for (int mirror = 0; mirror < 2; mirror++) {
        AsmRun run = {0}; run.arena = Arena_new(); run.layer_d_global = 10;
        float uv[50], xyz[4][75];
        for (int v = 0; v < 25; v++) {
            uv[2*v] = (float)(2*(v%5)); uv[2*v+1] = (float)(2*(v/5));
            for (int i = 0; i < 4; i++) {
                xyz[i][3*v] = (float)(12*i);
                xyz[i][3*v+1] = uv[2*v]; xyz[i][3*v+2] = uv[2*v+1];
            }
        }
        for (int i = 0; i < 4; i++) {
            AsmChart c = {0}; c.id = i; c.nv = 25; c.xyz = xyz[i]; c.uv = uv;
            c.pose_x = i < 2 ? 0 : i == 2 ? 100 : 200;
            c.flags = mirror ? ASM_CHART_MIRROR : 0; AsmRun_push_chart(&run,&c);
        }
        AsmChart saved[4]; memcpy(saved,run.charts,sizeof saved);
        size_t before = apl_internal_bad(&run,saved,4,0);
        /* Remove A/B's overlap but create an equal-sized overlap on C/D.
         * Aggregate bad counts do not establish that previously sound
         * material stayed sound. Both acceptance paths must reject this. */
        run.charts[1].pose_x = 50; run.charts[3].pose_x = 100;
        size_t after = apl_internal_bad(&run,saved,4,1);
        AsmContinuityAuditCache uv_cache = {0}, registration_cache = {0};
        uint8_t blocked[] = {0xa5,0xcc,0xcc,0xcc,0xcc,0x5a};
        uv_cache.blocked_charts = blocked+1; uv_cache.blocked_charts_count = 4;
        int uv_accepted = apl_uv_audit(&run,saved,4,&uv_cache);
        int registration_accepted = apl_registration_audit(&run,saved,4,&registration_cache);
        int ok = before > 0 && after == before && !uv_accepted && !registration_accepted &&
                 blocked[0] == 0xa5 && blocked[1] == 0 && blocked[2] == 0 &&
                 blocked[3] == 1 && blocked[4] == 1 && blocked[5] == 0x5a;
        /* Removing the original overlap without creating a new one is a
         * valid improvement, using the same original checkpoint/cache. */
        run.charts[3].pose_x = 200;
        ok = ok && !apl_internal_bad(&run,saved,4,1) && apl_uv_audit(&run,saved,4,&uv_cache) &&
             apl_registration_audit(&run,saved,4,&registration_cache);
        for (int i = 1; i <= 4; i++) ok = ok && !blocked[i];
        if (!ok) {
            fprintf(stderr,"  asm_place selftest FAIL: overlap exchanged between regions, mirror %d, bad %zu -> %zu, accepted UV %d registration %d\n",
                    mirror,before,after,uv_accepted,registration_accepted); fails++;
        }
        Arena_dispose(&run.arena);
    }
    return fails;
}

static int apl_parallel_rim_selftest(void)
{
    enum { N = 65, NV = N*N };
    AsmRun run = {0}; run.arena = Arena_new(); run.layer_d_global = 12;
    float *uv = ARENA_ALLOC(run.arena,2*NV*sizeof *uv);
    float *xyz[2];
    for (int i = 0; i < 2; i++) xyz[i] = ARENA_ALLOC(run.arena,3*NV*sizeof(float));
    for (int v = 0; v < NV; v++) {
        uv[2*v] = (float)(2*(v%N)); uv[2*v+1] = (float)(2*(v/N));
        for (int i = 0; i < 2; i++) {
            xyz[i][3*v] = uv[2*v]; xyz[i][3*v+1] = uv[2*v+1];
            xyz[i][3*v+2] = (float)(i && v/N >= N/2 ? 12 : 0);
        }
    }
    for (int i = 0; i < 2; i++) {
        AsmChart c = {0}; c.id = c.component = i; c.nv = NV; c.xyz = xyz[i]; c.uv = uv;
        AsmRun_push_chart(&run,&c);
    }
    AplRimGrid *grid = aprg_new(run.arena,24);
    aprg_insert_chart(grid,&run.charts[0]);
    const AsmChart *candidate = &run.charts[1]; int32_t orig[2] = {0,1};
    const double offsets[] = {0,23.999999999999,24,64,300};
    int fails = 0; size_t touched = 0;
    for (size_t i = 0; i < sizeof offsets/sizeof offsets[0]; i++) {
        run.charts[1].pose_x = offsets[i]; size_t ref_n, actual_n;
        double reference = apl_rim_audit(&run,orig,1,0,24,&ref_n);
        double actual = aprg_probe(grid,&run,&candidate,1,24,&actual_n);
        touched += actual_n;
        if (ref_n != actual_n || reference != actual) {
            fprintf(stderr,"  asm_place selftest FAIL: parallel rim at %.17g: %zu/%.17g vs %zu/%.17g\n",
                    offsets[i],actual_n,actual,ref_n,reference); fails++;
        }
    }
    uint8_t excluded[2] = {1,0}; size_t count = 0;
    run.charts[1].pose_x = 0;
    /* Record two odd-sized charts in one parallel probe. The field has an
     * independent analytic reference: only the raised half of chart1 can
     * violate this flat chart0 grid. Sentinels check both allocation ends. */
    size_t queries = (NV+1)/2, bad_count = 0, recorded_n = 0;
    uint8_t *record = ARENA_ALLOC(run.arena,2*queries+2);
    memset(record,0xcc,2*queries+2); record[0] = 0xa5; record[2*queries+1] = 0x5a;
    int32_t *partners = ARENA_ALLOC(run.arena,(2*queries+2)*sizeof *partners);
    for (size_t i = 0; i < 2*queries+2; i++) partners[i] = -2;
    partners[0] = 12345; partners[2*queries+1] = 54321;
    const AsmChart *both[] = {&run.charts[0],&run.charts[1]};
    double recorded = aprg_probe_record(grid,&run,both,2,24,&recorded_n,NULL,record+1,partners+1,NULL,NULL,2);
    int field_ok = recorded_n == 2*queries && record[0] == 0xa5 && record[2*queries+1] == 0x5a &&
                   partners[0] == 12345 && partners[2*queries+1] == 54321;
    for (int i = 0; i < 2; i++) for (int v = 0; v < NV; v += 2) {
        int expected = i == 1 && v/N >= N/2;
        field_ok = field_ok && record[1+i*queries+(size_t)v/2] == expected;
        field_ok = field_ok && partners[1+i*queries+(size_t)v/2] == (expected ? 0 : -1);
        bad_count += (size_t)expected;
    }
    field_ok = field_ok && recorded == (double)bad_count/(2*queries);
    if (!field_ok) { fprintf(stderr,"  asm_place selftest FAIL: parallel rim violation field or bounds\n"); fails++; }
    if (!touched || aprg_probe_excluding(grid,&run,&candidate,1,24,&count,excluded) != 0 || count != 0) {
        fprintf(stderr,"  asm_place selftest FAIL: parallel rim exclusion or empty support\n"); fails++;
    }
    Arena_dispose(&run.arena); return fails;
}

static int apl_axis_side_selftest(void)
{
    /* A dent curves away from the axis while its stored normal still points
     * outward. Curvature alone reverses every layer order in this patch.
     * World tilt, offset, UV parity and stored-normal parity must not change
     * the physical choice, or let the orientation pass overwrite it. */
    int fails=0;
    for (int world=0; world<2; world++) for (int mirror=0; mirror<2; mirror++) for (int parity=0; parity<2; parity++) {
        AsmRun run={0}; run.arena=Arena_new();
        enum { N=9, NV=N*N };
        float xyz[3*NV],normal[3*NV],uv[2*NV];
        double angle=world ? .4 : 0, ca=cos(angle),sa=sin(angle);
        double origin[3]={world ? 4800 : 0,world ? 3400 : 0,world ? 2900 : 0};
        double dir[3]={ca,0,sa}; int sign=parity ? -1 : 1;
        AsmChart c={0}; c.id=0; c.nv=NV; c.area3d=1600;
        c.xyz=xyz; c.nrm=normal; c.uv=uv; c.flags=mirror ? ASM_CHART_MIRROR : 0;
        for (int j=0; j<N; j++) for (int i=0; i<N; i++) {
            int v=j*N+i; double u=5*(i-4),z=5*(j-4),x=200+.01*u*u;
            double ny=-.02*u,nx=1,nl=hypot(ny,nx);
            double p[3]={ca*z-sa*x,u,sa*z+ca*x};
            double n[3]={-sa*nx,ny,ca*nx};
            for (int d=0; d<3; d++) {
                xyz[3*v+d]=(float)(origin[d]+p[d]); normal[3*v+d]=(float)(sign*n[d]/nl);
                c.centroid[d]+=xyz[3*v+d]/(double)NV;
            }
            uv[2*v]=(float)u; uv[2*v+1]=(float)z;
        }
        AsmAxis *axis=AsmAxis_line(run.arena,origin,dir,origin[0]-500,origin[0]+500);
        run.charts=&c; run.n_charts=1;
        int8_t side[1]={0};
        apl_component_out_sides(&run,0,NULL,side);
        int ok=side[0]==-sign;
        apl_component_out_sides(&run,0,axis,side);
        ok=ok && side[0]==sign;
        AplTurn turn[2]={{200,0,16},{200,0,16}}; AplOrient orientation={0};
        apl_orient_votes(&run,0,axis,turn,&orientation,side);
        ok=ok && side[0]==sign && orientation.has_turn && orientation.turn_u>.99 &&
            (mirror ? orientation.bend>.9 : orientation.bend<-.9);
        /* Opposing normals cannot furnish a chart-wide radial sign. */
        for (int v=0; v<NV/2; v++) for (int d=0; d<3; d++) normal[3*v+d]=-normal[3*v+d];
        double inward[3]={0};
        ok=ok && !apl_axis_inward(&c,axis,inward);
        if (!ok) { fprintf(stderr,"  asm_place selftest FAIL: dented sheet outward side world %d mirror %d normal %d, side %d bend %.6g\n",world,mirror,parity,side[0],orientation.bend); fails++; }
        Arena_dispose(&run.arena);
    }
    return fails;
}

/* THE TURN LENGTH AND HOW WELL IT IS KNOWN.  Two concentric cylindrical arcs, the placed one
 * at radius 240 and the candidate 12 vox outside it, each a 3 x 20 grid whose uv IS its arc
 * length.  The stretch is then s = 252/240 = 1.05 and one turn is 2*pi*12/0.05 = 2*pi*240 =
 * 1508 vox, which the chart curvature independently agrees with.  Three cases pin the
 * estimator and its new standard error:
 *   clean      the turn is right and (s-1) is known to ~0%
 *   noisy      a 5-vox sawtooth on the candidate's u leaves the same slope but makes the
 *              estimate worthless -- se(s-1) is comparable to (s-1) itself
 *   duplicated the SAME geometry with every hit recorded twice must not look better known:
 *              one layer pair's hits are one look at the geometry, which is the whole reason
 *              for the sqrt(hits/pairs) design effect */
static int apl_stretch_error_selftest(void)
{
    enum { ROWS = 3, COLS = 20, NV = ROWS * COLS, NF = 2 * (ROWS - 1) * (COLS - 1), NH = COLS - 1, PAIRS = 4 };
    int fails = 0;
    double rel_clean = -1.0, rel_noisy = -1.0, rel_dup = -1.0, t_clean = 0.0;
    const double ra = 240.0, d = 12.0, du = 20.0, dv = 20.0;
    const double s_true = (ra + d) / ra;
    for (int cs = 0; cs < 3; cs++) {
        AsmRun run; memset(&run, 0, sizeof run); run.arena = Arena_new();
        AsmChart charts[2]; memset(charts, 0, sizeof charts);
        float uv[2][2*NV], xyz[2][3*NV];
        int32_t faces[3*NF];
        for (int k = 0; k < 2; k++) {
            double r = k ? ra + d : ra;
            for (int row = 0; row < ROWS; row++) for (int col = 0; col < COLS; col++) {
                int v = row * COLS + col;
                double arc = du * (double)col * (k ? s_true : 1.0);   /* uv IS arc length on each arc */
                double noise = (k && cs != 0) ? ((col & 1) ? 5.0 : -5.0) : 0.0;
                uv[k][2*v] = (float)(arc + noise); uv[k][2*v+1] = (float)(dv * (double)row);
                double th = arc / r;                                   /* the 3-D arc the uv flattens */
                xyz[k][3*v]   = (float)(dv * (double)row);
                xyz[k][3*v+1] = (float)(r * sin(th));
                xyz[k][3*v+2] = (float)(r * cos(th));
            }
        }
        for (int row = 0, f = 0; row < ROWS - 1; row++) for (int col = 0; col < COLS - 1; col++) {
            int a = row * COLS + col, b = a + 1, c = a + COLS, e = c + 1;
            faces[3*f] = a; faces[3*f+1] = b; faces[3*f+2] = c; f++;
            faces[3*f] = b; faces[3*f+1] = e; faces[3*f+2] = c; f++;
        }
        for (int k = 0; k < 2; k++) {
            charts[k].id = k; charts[k].component = k; charts[k].nv = NV; charts[k].nf = NF;
            charts[k].uv = uv[k]; charts[k].xyz = xyz[k]; charts[k].faces = faces; charts[k].area3d = 1000.0;
        }
        int dup = cs == 2;
        /* PAIRS layer pairs over the same chart pair, each covering its own span of the arc:
         * the standard error's degrees of freedom are these, not the hits inside them */
        AsmLayerPair pairs[PAIRS]; memset(pairs, 0, sizeof pairs);
        AsmLayerHit hits[2*NH]; memset(hits, 0, sizeof hits);
        int nh_written = 0;
        for (int q = 0; q < PAIRS; q++) {
            int lo = q * NH / PAIRS, hi = (q + 1) * NH / PAIRS;
            pairs[q].a = 0; pairs[q].b = 1; pairs[q].side = 1; pairs[q].k = 1; pairs[q].d_median = (float)d;
            pairs[q].hit_first = nh_written;
            for (int i = lo; i < hi; i++) for (int rep = 0; rep < (dup ? 2 : 1); rep++) {
                AsmLayerHit *h = &hits[nh_written++];        /* face 2*i of row 0 has vertex i first */
                h->a = 0; h->b = 1; h->side = 1; h->order = 1;
                h->va = (int32_t)i; h->fb = (int32_t)(2*i); h->l0 = 1.0f; h->l1 = 0.0f;
                h->dist = (float)d;
            }
            pairs[q].hit_count = nh_written - pairs[q].hit_first;
        }
        run.charts = charts; run.n_charts = 2; run.layers = pairs; run.n_layers = PAIRS;
        run.layer_hits = hits; run.n_layer_hits = (size_t)nh_written; run.layer_d_global = d;
        uint8_t placed[2] = { 1, 0 };
        int8_t sides[2] = { 1, 1 };
        AplCtx x; memset(&x, 0, sizeof x);
        x.run = &run; x.comp_placed = placed; x.out_side = sides;
        double T = 0.0, sl = 0.0, dbar = 0.0, rel = -1.0; size_t nh = 0, np = 0;
        int got = apl_stretch_turn(&x, 1, &T, &sl, &dbar, &nh, &rel, &np);
        if (!got || np != (size_t)PAIRS || nh != (size_t)nh_written) {
            fprintf(stderr, "  asm_place stretch FAIL case %d: got %d, s %.4f, %zu hits on %zu pairs\n", cs, got, sl, nh, np);
            fails++; Arena_dispose(&run.arena); continue;
        }
        if (cs == 0) { rel_clean = rel; t_clean = fabs(T); }
        else if (cs == 1) rel_noisy = rel;
        else rel_dup = rel;
        Arena_dispose(&run.arena);
    }
    if (fabs(t_clean - 2.0 * APL_PI * ra) > 2.0) {
        fprintf(stderr, "  asm_place stretch FAIL: clean turn %.1f, want %.1f\n", t_clean, 2.0 * APL_PI * ra);
        fails++;
    }
    if (!(rel_clean >= 0.0 && rel_clean < 0.01)) {
        fprintf(stderr, "  asm_place stretch FAIL: a noiseless fit must know (s-1) exactly, relative error %.4f\n", rel_clean);
        fails++;
    }
    if (!(rel_noisy > 0.5)) {
        fprintf(stderr, "  asm_place stretch FAIL: a 5-vox sawtooth must leave (s-1) unknown, relative error %.4f\n", rel_noisy);
        fails++;
    }
    if (!(rel_dup >= 0.0 && fabs(rel_dup - rel_noisy) <= 1e-9 * fmax(1.0, rel_noisy))) {
        fprintf(stderr, "  asm_place stretch FAIL: recording every hit twice changed the relative error %.6f -> %.6f\n", rel_noisy, rel_dup);
        fails++;
    }
    if (!fails) fprintf(stderr, "  asm_place stretch error: turn %.0f vox, (s-1) known to %.2f%% clean, %.0f%% noisy, %.0f%% duplicated (ok)\n",
                        t_clean, 100.0*rel_clean, 100.0*rel_noisy, 100.0*rel_dup);
    return fails;
}

/* A seam the unimodality gate rejected is negative evidence about where two
 * charts sit.  Whether it may still anchor a placement is ASM_WRAP_VETO_ANCHOR's
 * decision, and this control pins the predicate to that switch either way -- a
 * WEAK seam without the flag always anchors, a DROPPED one never does. */
/* THE SECTION'S PERIMETER EXCESS.  A wrap of radius r whose cross-section is a circle plus a
 * straight flap of length C has perimeter 2 pi r + C, not 2 pi r.  Build three such wraps at
 * r, r+d, r+2d with measured outward turns, and check that:
 *   1. apl_fit_turn_radius recovers C and a free slope of 2 pi;
 *   2. apl_turn_radius agrees with the COMPILED ASM_PLACE_TURN_INTERCEPT -- so this is a
 *      regression lock on mode 0 and a correctness test on mode 1, and it passes in every
 *      cell of the A/B;
 *   3. apl_growth_for_hops agrees with the compiled ASM_PLACE_GROWTH_PER_WRAP;
 *   4. one wrap mis-gauged by a whole turn does not drag the intercept;
 *   5. degenerate input leaves the intercept at 0, with no NaN and no division by zero. */
static int apl_fit_turn_radius_selftest(void)
{
    enum { NC = 40 };
    const double R0 = 300.0, d = 12.0, C = 600.0;
    int fails = 0;
    for (int cs = 0; cs < 3; cs++) {
        AsmRun run; memset(&run, 0, sizeof run); run.arena = Arena_new();
        AsmChart charts[NC]; memset(charts, 0, sizeof charts);
        AplTurn turn[2 * NC]; memset(turn, 0, sizeof turn);
        int8_t sides[NC]; double r_axis[NC];
        float uv[6] = { 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 1.0f };   /* AsmChart_in_layout wants a uv */
        size_t nc = cs == 2 ? 4 : NC;               /* case 2: too few charts to adopt anything */
        for (size_t i = 0; i < nc; i++) {
            double r = R0 + (double)(i % 3) * d;    /* three concentric wraps */
            charts[i].id = (int32_t)i; charts[i].component = 0; charts[i].nv = 3; charts[i].nf = 1;
            charts[i].area3d = 100.0; charts[i].uv = uv;
            sides[i] = 1; r_axis[i] = r;
            double L = 2.0 * APL_PI * r + C;
            /* case 1: one chart in six is gauged a whole turn out -- the reject must shed it */
            if (cs == 1 && (i % 6) == 0) L *= 2.0;
            turn[i * 2 + 1].tx = L; turn[i * 2 + 1].ty = 0.0; turn[i * 2 + 1].n = 4;
        }
        run.charts = charts; run.n_charts = nc; run.layer_d_global = d;
        AplCtx x; memset(&x, 0, sizeof x);
        x.run = &run; x.turn = turn; x.out_side = sides; x.r_axis = r_axis;
        x.have_radius = 1; x.uout_ok = 1; x.uout_x = 1.0; x.uout_y = 0.0;
        apl_fit_turn_radius(&x);
        if (cs == 2) {
            if (x.turn_c != 0.0 || x.turn_c_n != 0) {
                fprintf(stderr, "  asm_place turn-c FAIL case 2: adopted c %.1f over %zu charts on too little support\n", x.turn_c, x.turn_c_n);
                fails++;
            }
            Arena_dispose(&run.arena); continue;
        }
        if (!(fabs(x.turn_c - C) <= 0.05 * C)) {
            fprintf(stderr, "  asm_place turn-c FAIL case %d: c %.1f, want %.1f +/- 5%%\n", cs, x.turn_c, C);
            fails++;
        }
        if (!(x.turn_c_n >= (size_t)ASM_PLACE_TURN_C_MIN_N)) {
            fprintf(stderr, "  asm_place turn-c FAIL case %d: support %zu under the %d bar\n", cs, x.turn_c_n, ASM_PLACE_TURN_C_MIN_N);
            fails++;
        }
        if (cs == 0 && !(fabs(x.turn_slope - 2.0 * APL_PI) <= 0.03 * 2.0 * APL_PI)) {
            fprintf(stderr, "  asm_place turn-c FAIL: free slope %.4f, want %.4f +/- 3%%\n", x.turn_slope, 2.0 * APL_PI);
            fails++;
        }
        if (cs == 0) {
            /* the k = 2 outward hop from chart 0, against the compiled mode */
            double tx = 0.0, ty = 0.0;
            int got = apl_turn_radius(&x, 0, 1, 2, &tx, &ty);
            double rbar = r_axis[0] + 2.0 * d * 0.5;
            double want = 2.0 * (2.0 * APL_PI * rbar) + (ASM_PLACE_TURN_INTERCEPT ? 2.0 * C : 0.0);
            if (!got || !(fabs(fabs(tx) - want) <= 0.01 * want) || ty != 0.0) {
                fprintf(stderr, "  asm_place turn-c FAIL: apl_turn_radius k=2 gave %.1f (got %d, ty %.3f), want %.1f at intercept mode %d\n",
                        fabs(tx), got, ty, want, ASM_PLACE_TURN_INTERCEPT);
                fails++;
            }
            /* the growth actually handed to the hops, against the compiled mode */
            double growth = 0.0;
            apl_growth_for_hops(&x, &run, turn, sides, &growth);
            double rule = 2.0 * APL_PI * d;
            int growth_ok = ASM_PLACE_GROWTH_PER_WRAP == 0 ? growth > 0.0
                          : ASM_PLACE_GROWTH_PER_WRAP == 1 ? fabs(growth - rule) < 1e-9
                          : growth >= 0.8 * rule - 1e-9 && growth <= 1.25 * rule + 1e-9;
            if (!growth_ok || growth != growth) {
                fprintf(stderr, "  asm_place turn-c FAIL: growth %.3f at mode %d (2 pi d = %.3f)\n",
                        growth, ASM_PLACE_GROWTH_PER_WRAP, rule);
                fails++;
            }
        }
        Arena_dispose(&run.arena);
    }
    {   /* no charts at all, and a chart at r = 0: no adoption, no NaN */
        AsmRun run; memset(&run, 0, sizeof run); run.arena = Arena_new();
        AsmChart one; memset(&one, 0, sizeof one);
        AplTurn t2[2]; memset(t2, 0, sizeof t2);
        int8_t side = 1; double r0 = 0.0;
        float uv1[6] = { 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 1.0f };
        one.id = 0; one.nv = 3; one.nf = 1; one.area3d = 1.0; one.uv = uv1;
        t2[1].tx = 1000.0; t2[1].n = 2;
        run.charts = &one; run.n_charts = 1; run.layer_d_global = d;
        AplCtx x; memset(&x, 0, sizeof x);
        x.run = &run; x.turn = t2; x.out_side = &side; x.r_axis = &r0; x.have_radius = 1; x.uout_ok = 1; x.uout_x = 1.0;
        apl_fit_turn_radius(&x);
        if (x.turn_c != 0.0 || x.turn_c_n != 0 || x.turn_slope != x.turn_slope) {
            fprintf(stderr, "  asm_place turn-c FAIL degenerate: c %.1f n %zu slope %.3f\n", x.turn_c, x.turn_c_n, x.turn_slope);
            fails++;
        }
        x.have_radius = 0;
        apl_fit_turn_radius(&x);
        if (x.turn_c != 0.0) { fprintf(stderr, "  asm_place turn-c FAIL: fitted without a radius model\n"); fails++; }
        Arena_dispose(&run.arena);
    }
    if (!fails) fprintf(stderr, "  asm_place section perimeter excess (intercept mode %d, growth mode %d): ok\n",
                        ASM_PLACE_TURN_INTERCEPT, ASM_PLACE_GROWTH_PER_WRAP);
    return fails;
}

static int apl_crosswrap_evidence_selftest(void)
{
    int fails = 0;
    AsmRelation r; memset(&r, 0, sizeof r); r.corr_count = 32; r.rms = 0.5;
    r.flags = ASM_REL_WEAK;
    if (apl_placement_seam(&r, 0, ASM_WRAP_VETO_ANCHOR) != 1) { fprintf(stderr, "  asm_place crosswrap FAIL: a WEAK seam must anchor\n"); fails++; }
    r.flags = ASM_REL_WEAK | ASM_REL_CROSSWRAP;
    /* unframed: the veto binds at 1 and at 2 */
    if (apl_placement_seam(&r, 0, ASM_WRAP_VETO_ANCHOR) != (ASM_WRAP_VETO_ANCHOR ? 0 : 1)) {
        fprintf(stderr, "  asm_place crosswrap FAIL: on an unframed pile a CROSSWRAP seam anchors against ASM_WRAP_VETO_ANCHOR %d\n", ASM_WRAP_VETO_ANCHOR);
        fails++;
    }
    /* framed: only the unconditional setting may take the anchor away */
    if (apl_placement_seam(&r, 1, ASM_WRAP_VETO_ANCHOR) != (ASM_WRAP_VETO_ANCHOR == 1 ? 0 : 1)) {
        fprintf(stderr, "  asm_place crosswrap FAIL: on a framed pile a CROSSWRAP seam does not follow ASM_WRAP_VETO_ANCHOR %d\n", ASM_WRAP_VETO_ANCHOR);
        fails++;
    }
    /* the whole 3 x 2 table, since the policy is now per pile and the constant is only its
     * default: the veto binds at 1 always, at 2 only unframed, and at 0 never */
    {
        int want[3][2] = { { 1, 1 }, { 0, 0 }, { 0, 1 } };   /* [veto][framed] -> anchors? */
        int v = 0, fr = 0;
        for (v = 0; v < 3; v++) for (fr = 0; fr < 2; fr++)
            if (apl_placement_seam(&r, fr, v) != want[v][fr]) {
                fprintf(stderr, "  asm_place crosswrap FAIL: veto %d on a %s pile gave %d, want %d\n",
                        v, fr ? "framed" : "unframed", apl_placement_seam(&r, fr, v), want[v][fr]);
                fails++;
            }
    }
    r.flags = ASM_REL_DROPPED | ASM_REL_CROSSWRAP;
    if (apl_placement_seam(&r, 0, ASM_WRAP_VETO_ANCHOR) != 0) { fprintf(stderr, "  asm_place crosswrap FAIL: a DROPPED seam must never anchor\n"); fails++; }
    return fails;
}

#include "asm_place_strict_test.inc"

int AsmPlace_selftest(void)
{
    /* Both an initially failing seam and an already passing but distorted
     * seam must improve without weakening the transaction's rollback. */
    int fails = apl_fit_turn_radius_selftest() + apl_crosswrap_evidence_selftest() + apl_stretch_error_selftest() + apl_rotated_cluster_selftest() + apl_v_model_update_selftest() + apl_self_contact_audit_selftest() + apl_axis_side_selftest() + apl_stopped_repair_selftest() + apl_seam_incumbent_selftest() + apl_registration_planes_selftest() + apl_guarded_unsampled_vertex_selftest() + apl_guarded_pair_identity_selftest() + apl_guarded_uv_audit_selftest() + apl_polish_bound_selftest() + apl_internal_overlap_exchange_selftest() + apl_parallel_rim_selftest() + apl_single_chart_uv_audit_selftest() + apl_rim_pair_selftest() + apl_registration_selftest() + apl_reconsider_selftest() + apl_uv_transaction_selftest(10) + apl_uv_transaction_selftest(5);
    {
        /* Incompatible model targets must not override a measured source
         * boundary. An unsuccessful repair restores the complete chart. */
        AsmRun run; memset(&run, 0, sizeof run); run.arena = Arena_new();
        AsmChart charts[2]; memset(charts, 0, sizeof charts); float uv[12], xyz[18];
        for (int i = 0; i < 6; i++) {
            uv[2*i] = xyz[3*i] = (float)(4*i); uv[2*i+1] = xyz[3*i+1] = xyz[3*i+2] = 0;
        }
        for (int i = 0; i < 2; i++) { charts[i].id = charts[i].component = i; charts[i].uv = uv; charts[i].xyz = xyz; charts[i].nv = 6; charts[i].area3d = 1000; }
        charts[0].placement_state = ASM_PLACE_ROOT; charts[1].pose_y = 12;
        AsmCorr corr[6]; memset(corr, 0, sizeof corr);
        for (int i = 0; i < 6; i++) { corr[i].va = corr[i].vb = i; corr[i].valid = 3; corr[i].run = 0; }
        AsmRelation rel; memset(&rel, 0, sizeof rel); rel.b = 1; rel.corr_count = 6; rel.continuity = ASM_CONT_SOURCE;
        run.charts = charts; run.n_charts = 2; run.rels = &rel; run.n_rels = 1; run.corr = corr; run.n_corr = 6;
        AplAnchor anchors[92]; memset(anchors, 0, sizeof anchors); double weights[92];
        for (int i = 0; i < 92; i++) {
            anchors[i].chart_c = 1; anchors[i].k = anchors[i].measured_order = 1; anchors[i].w = weights[i] = 1;
            anchors[i].local_u = anchors[i].tx = i%12; anchors[i].ty = i < 12 ? 10 : 200;
        }
        AsmPlaceOpts opts; AsmPlace_default_opts(&opts);
        uint8_t selected[1] = {0}; int32_t groups[2] = {0,1}; AsmChart proposal[1];
        AplCtx x; memset(&x, 0, sizeof x); x.run = &run; x.o = &opts; x.anchors = anchors; x.wts = weights; x.na = 92;
        x.orig_comp = groups; x.seam_selected = selected; x.trial = 1; x.proposal = proposal;
        AsmChart before = charts[1];
        int ok = !apl_place_and_audit(&x, 1, 0, 0, 0, NULL, 2) && memcmp(&charts[1], &before, sizeof before) == 0;
        if (!ok) { fprintf(stderr, "  asm_place selftest FAIL: incompatible model repair cannot override source\n"); fails++; }
        /* Agreement is an acceptance condition, not a prerequisite for
         * attempting a repair. These twelve observations begin outside
         * the support radius but agree after measured source closure. */
        x.na = 12; charts[1].pose_y = 28; before = charts[1];
        for (int i = 0; i < 12; i++) anchors[i].local_u = anchors[i].tx = 4*(i%6);
        ok = apl_place_and_audit(&x, 1, 0, 0, 0, NULL, 2) && proposal[0].placement_support == 12 &&
             fabs(proposal[0].pose_y) < 2 && memcmp(&charts[1], &before, sizeof before) == 0;
        if (!ok) { fprintf(stderr, "  asm_place selftest FAIL: unsupported proposal can recover measured agreement\n"); fails++; }
        Arena_dispose(&run.arena);
    }
    {
        AplAnchor anchors[16]; memset(anchors, 0, sizeof anchors); double weights[16];
        double cu = 350, cv = -80, slope = 0.02, truth = 0.2;
        double axis_v = sin(truth)*cu+cos(truth)*cv-20-slope*(cos(truth)*cu-sin(truth)*cv+1500);
        for (int i = 0; i < 16; i++) {
            weights[i] = 1; anchors[i].sx = cu+100; anchors[i].sy = cv+4*i-30;
            anchors[i].tx = cos(truth)*anchors[i].sx-sin(truth)*anchors[i].sy+1500;
            anchors[i].ty = sin(truth)*anchors[i].sx+cos(truth)*anchors[i].sy-20;
        }
        double theta, tx, ty;
        double rms = apl_fit_axis_translation(anchors, weights, 16, 0.6, cu, cv, axis_v, slope, 0, &theta, &tx, &ty);
        int ok = rms < 1e-4 && fabs(theta-truth) < 1e-5;
        /* A biased vertical model target tempts a free fit to translate
         * away from the axis. Fitting on the plane can use a rigid rotation
         * to improve the evidence, without deforming any chart. */
        for (int i = 0; i < 16; i++) anchors[i].ty += 20;
        double free_th, free_tx, free_ty;
        apl_fit_translation(anchors, weights, 16, 0.6, &free_th, &free_tx, &free_ty);
        free_ty = axis_v+slope*(cos(free_th)*cu-sin(free_th)*cv+free_tx)-sin(free_th)*cu-cos(free_th)*cv;
        double old_e2 = 0;
        for (int i = 0; i < 16; i++) {
            double x = cos(free_th)*anchors[i].sx-sin(free_th)*anchors[i].sy+free_tx-anchors[i].tx;
            double y = sin(free_th)*anchors[i].sx+cos(free_th)*anchors[i].sy+free_ty-anchors[i].ty;
            old_e2 += x*x+y*y;
        }
        rms = apl_fit_axis_translation(anchors, weights, 16, 0.6, cu, cv, axis_v, slope, 0, &theta, &tx, &ty);
        double mean_u = cos(theta)*cu-sin(theta)*cv+tx, mean_v = sin(theta)*cu+cos(theta)*cv+ty;
        ok = ok && rms < 0.5*sqrt(old_e2/16) && fabs(mean_v-axis_v-slope*mean_u) < 1e-8 && fabs(theta) <= 0.6;
        if (!ok) { fprintf(stderr, "  asm_place selftest FAIL: rigid fit with the height constraint\n"); fails++; }
        /* When the measured model admits a band, coherent observations
         * inside it must retain their exact rigid fit. Translating the
         * chart's arbitrary UV origin cannot change that fit. */
        ok = 1;
        for (int shifted = 0; shifted < 2; shifted++) {
            double du = shifted ? 10000 : 0, dv = shifted ? -20000 : 0;
            AplAnchor moved[16]; memcpy(moved,anchors,sizeof moved);
            for (int i = 0; i < 16; i++) { moved[i].sx += du; moved[i].sy += dv; }
            rms = apl_fit_axis_translation(moved,weights,16,0.6,cu+du,cv+dv,axis_v,slope,30,&theta,&tx,&ty);
            double eu = 1500-cos(truth)*du+sin(truth)*dv, ev = -sin(truth)*du-cos(truth)*dv;
            ok = ok && rms < 1e-5 && fabs(theta-truth) < 1e-6 && fabs(tx-eu) < 1e-4 && fabs(ty-ev) < 1e-4;
        }
        /* Evidence beyond the band is still constrained, jointly in angle
         * and translation. The optimum lies on its nearest boundary. */
        for (int i = 0; i < 16; i++) anchors[i].ty += 35;
        rms = apl_fit_axis_translation(anchors,weights,16,0.6,cu,cv,axis_v,slope,30,&theta,&tx,&ty);
        mean_u = cos(theta)*cu-sin(theta)*cv+tx; mean_v = sin(theta)*cu+cos(theta)*cv+ty;
        double plane_rms = apl_fit_axis_translation(anchors,weights,16,0.6,cu,cv,axis_v,slope,0,&free_th,&free_tx,&free_ty);
        ok = ok && rms > 0 && rms < plane_rms && fabs(mean_v-axis_v-slope*mean_u-30) < 1e-8 && fabs(theta) <= 0.6;
        if (!ok) { fprintf(stderr,"  asm_place selftest FAIL: measured height uncertainty and bounded fit\n"); fails++; }
    }
    {
        /* A vertical column constrains rotation as well as a horizontal
         * one. The old U-only extent clamps the correct 0.2-radian fit. */
        AsmRun run; memset(&run, 0, sizeof run); run.arena = Arena_new();
        AplAnchor anchors[16]; memset(anchors, 0, sizeof anchors);
        double weights[16], residuals[16], th, tx, ty; size_t kept;
        AplCtx x; memset(&x, 0, sizeof x); x.run = &run; x.anchors = anchors; x.wts = weights; x.cand_grad = 1;
        for (int i = 0; i < 16; i++) {
            anchors[i].sx = 100; anchors[i].sy = 40*i; anchors[i].w = 1; anchors[i].k = 1;
            anchors[i].tx = cos(0.2)*anchors[i].sx-sin(0.2)*anchors[i].sy+1500;
            anchors[i].ty = sin(0.2)*anchors[i].sx+cos(0.2)*anchors[i].sy-20;
        }
        double rms = apl_fit_cluster(&x, 16, 0, 12, 0.6, &th, &tx, &ty, &kept, residuals);
        if (rms > 1e-8 || kept != 16 || fabs(th-0.2) > 1e-8 || fabs(tx-1500) > 1e-8 || fabs(ty+20) > 1e-8) {
            fprintf(stderr, "  asm_place selftest FAIL: two-dimensional rotation evidence\n"); fails++;
        }
        Arena_dispose(&run.arena);
    }
    {
        AsmRun run; memset(&run, 0, sizeof run); run.arena = Arena_new();
        AsmChart charts[2]; memset(charts, 0, sizeof charts);
        float uv[8] = {0,0,20,0,0,20,20,20}, xyz[2][12];
        for (int k = 0; k < 2; k++) {
            charts[k].id = k; charts[k].uv = uv; charts[k].nv = 4; charts[k].xyz = xyz[k]; charts[k].pose_x = 200*k;
            for (int v = 0; v < 4; v++) { xyz[k][3*v] = uv[2*v]+200*k; xyz[k][3*v+1] = uv[2*v+1]; xyz[k][3*v+2] = 0; }
        }
        run.charts = charts; run.n_charts = 2; run.layer_d_global = 10;
        int ok = apl_internal_bad(&run, charts, 2, 0) == 0;
        charts[1].pose_x = 0;
        size_t overlap = apl_internal_bad(&run, charts, 2, 0);
        ok = ok && overlap > 0;
        for (int k = 0; k < 2; k++) apl_compose(&charts[k], 0.7, 300, -100, 1);
        ok = ok && apl_internal_bad(&run, charts, 2, 0) == overlap;
        if (!ok) { fprintf(stderr, "  asm_place selftest FAIL: internal repair collision / rigid invariance\n"); fails++; }
        Arena_dispose(&run.arena);
    }
    {
        /* A closer multi-layer probe must not suppress an eligible adjacent
         * probe at another location on the same chart. Check both directions
         * and ensure the excluded pair never supplies a measured turn. */
        AsmRun run = {0}; run.arena = Arena_new();
        AsmChart charts[3] = {0}; float uv[6] = {0,0,10,0,0,10}; int32_t face[3] = {0,1,2};
        for (int i = 0; i < 3; i++) {
            charts[i].id = i; charts[i].uv = uv; charts[i].nv = 3;
            charts[i].faces = face; charts[i].nf = 1; charts[i].pose_x = 1000*i;
        }
        AsmLayerPair pairs[2] = {0}; AsmLayerHit hits[12] = {0};
        for (int i = 0; i < 2; i++) {
            pairs[i].k = (int8_t)(i+1); pairs[i].d_median = i ? 5 : 20;
            pairs[i].hit_first = 6*i; pairs[i].hit_count = 6;
            for (int j = 0; j < 6; j++) {
                AsmLayerHit *h = &hits[6*i+j]; h->va = j%3;
                h->l0 = j%3 == 0; h->l1 = j%3 == 1; h->dist = (float)pairs[i].d_median;
                h->order = (int8_t)(i+1);
            }
        }
        run.charts = charts; run.n_charts = 3; run.layers = pairs; run.n_layers = 2;
        run.layer_hits = hits; run.n_layer_hits = 12;
        uint8_t placed[1] = {1}; int8_t sides[3] = {1,1,1}; AplTurn turns[6]; int ok = 1;
        for (int reverse = 0; reverse < 2; reverse++) {
            for (int i = 0; i < 2; i++) {
                pairs[i].a = reverse ? i+1 : 0; pairs[i].b = reverse ? 0 : i+1;
                pairs[i].side = reverse ? -1 : 1;
            }
            apl_measure_turns(&run,placed,sides,NULL,turns,3);
            ok = ok && turns[1].n == 6 && fabs(turns[1].tx-1000) < 1e-9 && fabs(turns[1].ty) < 1e-9 &&
                 turns[2].n == 6 && fabs(turns[2].tx+1000) < 1e-9 && turns[4].n == 0;
        }
        pairs[0].k = 2;
        apl_measure_turns(&run,placed,sides,NULL,turns,3);
        for (int i = 0; i < 6; i++) ok = ok && turns[i].n == 0;
        if (!ok) { fprintf(stderr,"  asm_place selftest FAIL: eligible adjacent turn cutoff\n"); fails++; }
        Arena_dispose(&run.arena);
    }
    {
        /* Measured turn vectors survive an armed, disagreeing radius model,
         * in either probing direction. Missing measurements use the model
         * with estimated provenance. */
        AsmRun run; memset(&run, 0, sizeof run);
        AsmChart charts[2]; memset(charts, 0, sizeof charts);
        float uv[6] = {0,0,10,0,0,10}; int32_t face[3] = {0,1,2};
        for (int i = 0; i < 2; i++) {
            charts[i].id = charts[i].component = i; charts[i].uv = uv; charts[i].nv = 3;
            charts[i].faces = face; charts[i].nf = 1; charts[i].area3d = 50;
        }
        AsmLayerPair lp; memset(&lp, 0, sizeof lp); lp.b = 1; lp.side = lp.k = 1; lp.hit_count = 1; lp.d_median = 10;
        AsmLayerHit hit; memset(&hit, 0, sizeof hit); hit.l0 = 1; hit.order = 1;
        run.charts = charts; run.n_charts = 2; run.layers = &lp; run.n_layers = 1; run.layer_hits = &hit; run.layer_d_global = 10;
        AsmPlaceOpts opts; AsmPlace_default_opts(&opts); opts.max_layers = 1;
        uint8_t placed[2] = {1,0}; int8_t side[2] = {1,1}; double radius[2] = {100,110};
        AplTurn turns[4]; memset(turns, 0, sizeof turns); turns[1] = (AplTurn){1000,20,12};
        AplAnchor anchors[4]; memset(anchors, 0, sizeof anchors);
        AplCtx x; memset(&x, 0, sizeof x); x.run = &run; x.o = &opts; x.comp_placed = placed;
        x.out_side = side; x.turn = turns; x.anchors = anchors; x.hcap = 4;
        x.have_radius = x.uout_ok = 1; x.uout_x = 1; x.r_axis = radius;
        int ok = 1;
        for (int rev = 0; rev < 2; rev++) {
            lp.a = rev; lp.b = 1-rev; lp.side = rev ? -1 : 1;
            ok = ok && apl_collect_anchors(&x, 1) == 1 && anchors[0].Tx == 1000 && anchors[0].Ty == 20 &&
                 anchors[0].measured_turn && !anchors[0].estimated && anchors[0].reverse == rev;
        }
        memset(turns, 0, sizeof turns);
        ok = ok && apl_collect_anchors(&x, 1) == 1 && fabs(anchors[0].Tx-2*APL_PI*105) < 1e-8 &&
             !anchors[0].measured_turn && anchors[0].estimated;
        if (!ok) { fprintf(stderr, "  asm_place selftest FAIL: measured turn precedence and model provenance\n"); fails++; }
        /* Parallel offsets of a convex ellipse differ in perimeter by
         * 2*pi*d, despite different local radii at its two axes. The turn
         * between neighboring offsets uses their mean perimeter. A measured
         * opposite turn therefore predicts the next turn consistently at
         * both azimuths; 2*pi*local_radius does not. */
        double perimeter = 0;
        for (int i = 0; i < 4096; i++) {
            double a = 2*APL_PI*((double)i+0.5)/4096;
            perimeter += hypot(240*sin(a),120*cos(a))*2*APL_PI/4096;
        }
        x.growth = 2*APL_PI*run.layer_d_global;
        ok = 1;
        for (int axis = 0; axis < 2; axis++) for (int direction = -1; direction <= 1; direction += 2)
        for (int rotated = 0; rotated < 2; rotated++) for (int rev = 0; rev < 2; rev++) {
            double theta = rotated ? 0.37 : 0;
            double ct = cos(theta), sn = sin(theta);
            double opposite = -direction*(perimeter-direction*x.growth/2);
            double expected = direction*(perimeter+direction*x.growth/2);
            radius[0] = axis ? 240 : 120;
            x.uout_x = ct; x.uout_y = sn;
            memset(turns,0,sizeof turns);
            turns[direction > 0 ? 0 : 1] = (AplTurn){opposite*ct,opposite*sn,12};
            lp.a = rev; lp.b = 1-rev; lp.side = rev ? -direction : direction;
            ok = ok && apl_collect_anchors(&x,1) == 1 &&
                fabs(anchors[0].Tx-expected*ct) < 1e-8 && fabs(anchors[0].Ty-expected*sn) < 1e-8 &&
                anchors[0].estimated && !anchors[0].measured_turn && anchors[0].reverse == rev;
        }
        if (!ok) { fprintf(stderr,"  asm_place selftest FAIL: opposite turn on a noncircular section\n"); fails++; }
    }
    {
        AsmRun run; memset(&run, 0, sizeof run); run.arena = Arena_new();
        AplAnchor a[24]; double du[24], w[8], u[8], v[8]; size_t order[24]; int kind[8];
        AplCtx x; memset(&x, 0, sizeof x); x.run = &run; x.anchors = a;
        memset(a, 0, sizeof a);
        const double px[4] = {0,7,10,14}, pw[4] = {1,3,3,2};
        for (int i = 0; i < 4; i++) {
            a[i].tx = du[i] = px[i]; a[i].w = pw[i]; a[i].k = 1; a[i].cluster = -1; order[i] = (size_t)i;
        }
        int n = apl_cluster_pass(&x, 4, order, du, 10, 10, 1, 0, 0, 0, w, u, v, kind);
        int ok = n == 1 && w[0] == 8 && a[0].cluster == -1;
        memset(a, 0, sizeof a);
        for (int i = 0; i < 6; i++) {
            a[i].tx = du[i] = i; a[i].ty = i < 3 ? -1000 : 1000;
            a[i].w = 1; a[i].k = 1; a[i].cluster = -1; order[i] = (size_t)i;
        }
        n = apl_cluster_pass(&x, 6, order, du, 10, 10, 2, 0, 0, 0, w, u, v, kind);
        ok = ok && n == 2 && w[0] == 3 && w[1] == 3 && a[0].cluster != a[5].cluster;
        /* Compare against exhaustive rectangles on small, weighted point
         * sets. This checks selection, excluded seams and tied coordinates. */
        uint32_t seed = 71293;
        for (int trial = 0; trial < 40; trial++) {
            memset(a, 0, sizeof a);
            AplDuIdx key[24];
            for (int i = 0; i < 24; i++) {
                seed = seed*1664525u+1013904223u; a[i].tx = du[i] = (int)(seed%41)-20;
                seed = seed*1664525u+1013904223u; a[i].ty = (int)(seed%101)-50;
                seed = seed*1664525u+1013904223u; a[i].w = 1+seed%7;
                a[i].k = i%5 != 0; a[i].cluster = -1; key[i] = (AplDuIdx){du[i], (size_t)i};
            }
            qsort(key, 24, sizeof key[0], apl_cmp_duidx);
            for (int i = 0; i < 24; i++) order[i] = key[i].idx;
            double expected = 0;
            for (int i = 0; i < 24; i++) for (int j = 0; j < 24; j++) {
                double score = 0;
                for (int k = 0; k < 24; k++) if (a[k].k > 0 && a[k].tx >= a[i].tx &&
                    a[k].tx-a[i].tx <= 10 && a[k].ty >= a[j].ty && a[k].ty-a[j].ty <= 14) score += a[k].w;
                expected = fmax(expected, score);
            }
            n = apl_cluster_pass(&x, 24, order, du, 10, 7, 1, 0, 0, 0, w, u, v, kind);
            ok = ok && n == 1 && w[0] == expected;
        }
        if (!ok) { fprintf(stderr, "  asm_place selftest FAIL: bounded two-dimensional anchor clusters\n"); fails++; }
        Arena_dispose(&run.arena);
    }
    {
        AsmRun run; memset(&run, 0, sizeof run); run.arena = Arena_new();
        AsmChart charts[2]; memset(charts, 0, sizeof charts);
        charts[0].placement_state = ASM_PLACE_ROOT; charts[1].pose_x = 1000;
        run.charts = charts; run.n_charts = 2;
        AplAnchor anchors[3]; memset(anchors, 0, sizeof anchors);
        for (int i = 0; i < 3; i++) {
            anchors[i].k = 1; anchors[i].measured_order = 1; anchors[i].chart_c = 1;
            anchors[i].Tx = 1000; anchors[i].tx = 1000; anchors[i].cluster = i;
        }
        anchors[1].k = 2; anchors[1].tx = anchors[1].Tx = 2000; /* alternate for the same ray */
        anchors[2].observation = 1; anchors[2].local_u = 10; anchors[2].tx = 1010;
        AplCtx x; memset(&x, 0, sizeof x); x.run = &run; x.anchors = anchors; x.na = 3; x.n_observations = 2;
        AplConsensus good = apl_candidate_consensus(&x);
        charts[1].pose_x += 40;
        AplConsensus uncertain = apl_candidate_consensus(&x);
        charts[1].pose_x += 460;
        AplConsensus bad = apl_candidate_consensus(&x);
        if (good.observations != 2 || good.inliers != 2 || good.model_inliers != 2 || good.model_loss > 1e-10 ||
            uncertain.observations != 2 || uncertain.inliers != 0 || uncertain.model_inliers != 2 ||
            bad.observations != 2 || bad.inliers != 0 || bad.model_inliers != 0 || bad.model_loss <= good.model_loss) {
            fprintf(stderr, "  asm_place selftest FAIL: independent full-evidence consensus\n"); fails++;
        }
        Arena_dispose(&run.arena);
    }
    {
        /* A short boundary limits the proposal rotation, while its long
         * connected body fails the axial prior in that provisional pose.
         * The measured seam recovers the exact pose before the final gate. */
        AsmRun run; memset(&run, 0, sizeof run); run.arena = Arena_new();
        float small[12], long_uv[24];
        for (int i = 0; i < 6; i++) {
            small[2*i] = long_uv[2*i] = (float)(4*i);
            small[2*i+1] = long_uv[2*i+1] = 0;
            long_uv[12+2*i] = (float)(1000+4*i); long_uv[13+2*i] = 0;
        }
        AsmChart charts[3]; memset(charts, 0, sizeof charts);
        for (int i = 0; i < 3; i++) {
            charts[i].id = i; charts[i].component = i ? 1 : 0;
            charts[i].uv = i == 1 ? long_uv : small; charts[i].nv = i == 1 ? 12 : 6;
            charts[i].area3d = 1000;
        }
        charts[0].placement_state = ASM_PLACE_ROOT;
        charts[1].pose_theta = charts[2].pose_theta = -0.2;
        charts[2].pose_x = 1000*cos(-0.2)-4*sin(-0.2);
        charts[2].pose_y = 1000*sin(-0.2)+4*cos(-0.2);
        AsmChart saved[2] = {charts[1], charts[2]}, proposal[2];
        AsmRelation rels[2]; memset(rels, 0, sizeof rels);
        AsmCorr corr[12]; memset(corr, 0, sizeof corr);
        for (int j = 0; j < 2; j++) {
            rels[j].a = j; rels[j].b = j+1; rels[j].corr_first = 6*j; rels[j].corr_count = 6;
            rels[j].continuity = ASM_CONT_SOURCE;
            for (int i = 0; i < 6; i++) {
                AsmCorr *c = &corr[6*j+i]; c->va = 6*j+i; c->vb = i;
                c->valid = 3; c->run = 0; c->gap_a[1] = c->gap_b[1] = 4;
            }
        }
        rels[0].flags = ASM_REL_SWITCHED;
        run.charts = charts; run.n_charts = 3; run.rels = rels; run.n_rels = 2; run.corr = corr; run.n_corr = 12;
        AsmPlaceOpts o; AsmPlace_default_opts(&o);
        AsmPlaceStats stats; memset(&stats, 0, sizeof stats);
        AplCtx x; memset(&x, 0, sizeof x); x.run = &run; x.o = &o; x.st = &stats;
        uint8_t placed[2] = {1,0}, selected[2] = {1,0}; int32_t groups[3] = {0,1,1};
        double cu[3] = {10,510,10}, cv[3] = {0,0,0}, axis[3] = {0,4,8}, weights[6];
        AplAnchor anchors[6]; memset(anchors, 0, sizeof anchors);
        x.comp_placed = placed; x.orig_comp = groups; x.primary = 0; x.trial = 1; x.proposal = proposal;
        x.seam_selected = selected; x.anchors = anchors; x.wts = weights; x.hcap = 6; x.cluster_kind = 1;
        x.cu = cu; x.cv = cv; x.s_axis = axis; x.vm_ok = x.vm_pin = 1; x.vm_b = 1;
        int ok = apl_collect_anchors(&x, 1) == 6;
        x.na = 6;
        for (int i = 0; i < 6; i++) weights[i] = 1;
        ok = ok && !apl_candidate_shape(&x, saved, 2, 0);
        ok = ok && apl_place_and_audit(&x, 1, 0, 0, 0, NULL, 2) && x.nproposal == 2 &&
             fabs(proposal[0].pose_theta) < 1e-6 && fabs(proposal[0].pose_y-4) < 1e-6 &&
             fabs(proposal[1].pose_x-1000) < 1e-5 && fabs(proposal[1].pose_y-8) < 1e-5 &&
             memcmp(&charts[1], &saved[0], sizeof(AsmChart)) == 0 &&
             memcmp(&charts[2], &saved[1], sizeof(AsmChart)) == 0;
        if (ok) {
            charts[1] = proposal[0]; charts[2] = proposal[1];
            charts[2].pose_y += 200; /* a later solve may not evade the final gate */
            ok = !apl_candidate_shape(&x, saved, 2, 0);
            /* A diagnostic may expose the next gate without changing the
             * ordinary policy. A displaced physical seam still fails its
             * own certificate, and non-finite geometry is never a prior. */
            x.diagnostic_axis_policy = 1;
            ok = ok && apl_candidate_shape(&x, saved, 2, 0);
            size_t support = 0; int32_t parent = -1; double error = 0;
            ok = ok && !apl_trial_support(&x, 1, NULL, 0, &support, &parent, &error, 0);
            x.cluster_kind = 0;
            ok = ok && !apl_candidate_shape(&x, saved, 2, 0);
            x.diagnostic_axis_policy = 2;
            ok = ok && apl_candidate_shape(&x, saved, 2, 0);
            double finite_y = charts[2].pose_y;
            charts[2].pose_y = NAN;
            ok = ok && !apl_candidate_shape(&x, saved, 2, 0);
            charts[2].pose_y = finite_y;
            x.diagnostic_axis_policy = 0;
            ok = ok && !apl_candidate_shape(&x, saved, 2, 0);
        }
        if (!ok) { fprintf(stderr, "  asm_place selftest FAIL: final geometry gates after seam repair\n"); fails++; }
        Arena_dispose(&run.arena);
    }
    {
        /* Alternate wrap votes cannot split a sound source-connected piece.
         * A physically failed seam may separate, and historical cuts remain
         * recoverable. Upstream hard vetoes remain hard. */
        AsmRun run; memset(&run, 0, sizeof run); run.arena = Arena_new();
        AsmChart charts[3]; memset(charts, 0, sizeof charts);
        float uv[12];
        for (int i = 0; i < 6; i++) { uv[2*i] = (float)(4*i); uv[2*i+1] = 0; }
        for (int i = 0; i < 3; i++) {
            charts[i].id = i; charts[i].component = i ? 1 : 0;
            charts[i].uv = uv; charts[i].nv = 6; charts[i].area3d = 1000;
        }
        charts[2].pose_y = 4;
        AsmCorr corr[6]; memset(corr, 0, sizeof corr);
        for (int i = 0; i < 6; i++) {
            corr[i].va = corr[i].vb = i; corr[i].valid = 3; corr[i].run = 0;
            corr[i].gap_a[1] = corr[i].gap_b[1] = 4;
        }
        AsmRelation rel; memset(&rel, 0, sizeof rel);
        rel.a = 1; rel.b = 2; rel.corr_count = 6;
        rel.continuity = ASM_CONT_SOURCE | ASM_CONT_VERIFIED;
        run.charts = charts; run.n_charts = 3; run.rels = &rel; run.n_rels = 1;
        run.corr = corr; run.n_corr = 6; run.continuity_ready = 1;
        AsmPlaceOpts o; AsmPlace_default_opts(&o);
        AplAnchor anchors[6]; memset(anchors, 0, sizeof anchors);
        anchors[0].chart_c = 1; anchors[0].cluster = 0;
        anchors[1].chart_c = 2; anchors[1].cluster = 1;
        AplCtx x; memset(&x, 0, sizeof x); x.run = &run; x.o = &o; x.anchors = anchors; x.t0 = 1000; x.hcap = 6;
        double weight[2] = {10,9}, offset[2] = {0,1000}, rms[2] = {1,1};
        int32_t rel_off[4] = {0,0,1,2}, rel_adj[2] = {0,0}, groups[3] = {0,1,1};
        int32_t rings[4] = {0,1,0,0}, ncomp = 2;
        double areas[4] = {1000,2000,0,0}, layer[3] = {0,1,1};
        AplOrient ori[4]; memset(ori, 0, sizeof ori);
        AplCompDiag diag[4]; memset(diag, 0, sizeof diag);
        size_t deferred = 0;
        int newc = apl_cut_component(&x, 1, 2, 2, weight, offset, rms, 8,
                                    rel_off, rel_adj, &ncomp, 4, areas, groups,
                                    ori, diag, rings, layer, &deferred);
        int ok = newc == -1 && ncomp == 2 && deferred == 0 && x.cuts_preserved_source == 1 &&
                 charts[1].component == 1 && charts[2].component == 1 &&
                 groups[1] == 1 && groups[2] == 1 && areas[1] == 2000 && areas[2] == 0 &&
                 rel.flags == 0 && rel.continuity == (ASM_CONT_SOURCE | ASM_CONT_VERIFIED);
        charts[2].pose_y += 20; /* Actual disagreement, not just a layer vote. */
        newc = apl_cut_component(&x, 1, 2, 2, weight, offset, rms, 8,
                                rel_off, rel_adj, &ncomp, 4, areas, groups,
                                ori, diag, rings, layer, &deferred);
        ok = ok && newc == 2 && ncomp == 3 && deferred == 1 &&
                 !AsmRel_is_join(&rel) && (rel.flags & ASM_REL_SWITCHED) &&
                 !(rel.flags & ASM_REL_DROPPED) && rel.continuity == ASM_CONT_SOURCE;
        charts[2].pose_y -= 20;
        /* Keep the recovery regression for checkpoints made by older cuts. */
        rel.continuity |= ASM_CONT_CUT;
        uint8_t placed[4] = {1,0,0,0}; x.comp_placed = placed;
        charts[1].component = 0; charts[1].placement_state = ASM_PLACE_ROOT;
        ok = ok && apl_collect_anchors(&x, newc) == 6 && x.na_seam == 6;
        /* A selected hop fits perfectly, but tears the retained boundary
         * apart by one turn. It cannot bypass the cut's source constraint. */
        x.orig_comp = groups; charts[2].pose_x = 1000;
        AsmContinuityTarget targets[12];
        for (int i = 0; i < 12; i++) targets[i] = (AsmContinuityTarget){2, (double)i, 0, 1000+i, 4, 1, 1};
        size_t support; int32_t parent; double error;
        ok = ok && !apl_trial_support(&x, newc, targets, 12, &support, &parent, &error, 0);
        charts[2].pose_x = 0;
        charts[2].component = 0; charts[2].placement_state = ASM_PLACE_SEAM;
        ok = ok && AsmContinuity_confirm(&run) == 1 && AsmRel_is_join(&rel) &&
             (rel.continuity & ASM_CONT_VERIFIED) && (rel.flags & ASM_REL_READMIT);
        for (int veto = 0; veto < 2; veto++) {
            uint32_t hard = veto ? ASM_REL_CONTACT : ASM_REL_DROPPED;
            rel.flags = ASM_REL_SWITCHED | hard; rel.continuity = ASM_CONT_SOURCE;
            ok = ok && AsmContinuity_confirm(&run) == 0 && (rel.flags & hard) && !AsmRel_is_join(&rel);
        }
        if (!ok) { fprintf(stderr, "  asm_place selftest FAIL: source-preserving cut / boundary recovery / hard veto\n"); fails++; }
        Arena_dispose(&run.arena);
    }
    {
        float uv[] = {0,0, 1000,0, 0,10};
        AsmChart a; memset(&a, 0, sizeof a); a.uv = uv; a.nv = 3;
        AsmChart b = a; b.pose_theta = 2*APL_PI;
        int ok = apl_same_pose(&a, &b, 8.0);
        b = a; b.flags = ASM_CHART_MIRROR;
        ok = ok && !apl_same_pose(&a, &b, 8.0);
        b = a; b.pose_theta = 0.01; /* 10 vox at the far end, despite a small angle */
        ok = ok && !apl_same_pose(&a, &b, 8.0);
        uv[2] = 10; b.pose_theta = 0.1; /* about one voxel, despite the old 128*angle proxy */
        ok = ok && apl_same_pose(&a, &b, 8.0);
        b = a; b.pose_x = 8.0;
        ok = ok && apl_same_pose(&a, &b, 8.0);
        b.pose_x = 8.000001;
        ok = ok && !apl_same_pose(&a, &b, 8.0);
        if (!ok) { fprintf(stderr, "  asm_place selftest FAIL: material pose equivalence\n"); fails++; }
    }
    {
        AsmRun run; memset(&run, 0, sizeof run); run.layer_d_global = 10;
        double spacing[2] = {10,10}; run.chart_layer_d = spacing;
        AsmPlaceOpts o; AsmPlace_default_opts(&o);
        AplCtx x; memset(&x, 0, sizeof x); x.run = &run; x.o = &o;
        AsmLayerPair lp; memset(&lp, 0, sizeof lp); lp.a = 0; lp.b = 1; lp.k = 1; lp.d_median = 10;
        int ok = apl_hop_upper(&x, &lp) == 1;
        lp.d_median = 20; /* first observed surface, with an intermediate hole */
        ok = ok && apl_hop_upper(&x, &lp) >= 2;
        lp.k = 2; lp.d_median = 10; /* spacing must not erase an observed crossing */
        ok = ok && apl_hop_upper(&x, &lp) >= 2;
        lp.d_median = 100;
        ok = ok && apl_hop_upper(&x, &lp) == o.max_layers;
        if (!ok) { fprintf(stderr, "  asm_place selftest FAIL: missing-layer hop hypotheses\n"); fails++; }
    }
    {
        /* A short observed arc supports the incumbent within measurement
         * tolerance, while its unconstrained correction exceeds the trust
         * region. Keep the incumbent, and prove that silence/wrong crossing
         * order still fail for both small and large pieces. */
        AsmRun run; memset(&run, 0, sizeof run); run.arena = Arena_new();
        AsmChart charts[2]; memset(charts, 0, sizeof charts);
        float uv[24];
        for (int i = 0; i < 12; i++) { uv[2*i] = (float)(4*i-22); uv[2*i+1] = 0; }
        for (int i = 0; i < 2; i++) { charts[i].id = i; charts[i].component = i; charts[i].uv = uv; charts[i].nv = 12; }
        charts[0].placement_state = ASM_PLACE_ROOT;
        run.charts = charts; run.n_charts = 2;
        AsmChart original = charts[1], proposal[1];
        AplAnchor anchors[12]; memset(anchors, 0, sizeof anchors);
        double weights[12]; int32_t groups[2] = {0,1};
        AsmContinuityTarget targets[12];
        for (int i = 0; i < 12; i++) {
            AplAnchor *a = &anchors[i];
            a->sx = a->local_u = uv[2*i]; a->tx = 1000+cos(0.5)*a->sx; a->ty = sin(0.5)*a->sx;
            a->chart_p = 0; a->chart_c = 1; a->k = a->measured_order = 1; a->w = weights[i] = 1;
            targets[i] = (AsmContinuityTarget){1,a->local_u,0,a->tx,a->ty,0.05,0};
        }
        charts[1].pose_x = 1000;
        int ok = !AsmContinuity_close(&run, &original, 1, NULL, targets, 12);
        charts[1] = original;
        AsmPlaceOpts o; AsmPlace_default_opts(&o);
        AplCtx x; memset(&x, 0, sizeof x);
        x.run = &run; x.o = &o; x.anchors = anchors; x.wts = weights; x.na = 12;
        x.orig_comp = groups; x.primary = 0; x.trial = 1; x.proposal = proposal;
        for (int size = 0; size < 2; size++) {
            charts[1].area3d = original.area3d = size ? 1e6 : 1e3;
            ok = ok && apl_place_and_audit(&x, 1, 0, 1000, 0, NULL, 2);
            ok = ok && x.nproposal == 1 && proposal[0].pose_x == 1000 && proposal[0].pose_theta == 0 &&
                 proposal[0].placement_state == ASM_PLACE_LAYER && proposal[0].placement_support == 12 &&
                 memcmp(&charts[1], &original, sizeof original) == 0;
            for (int i = 0; i < 12; i++) anchors[i].measured_order = 2;
            ok = ok && !apl_place_and_audit(&x, 1, 0, 1000, 0, NULL, 2) && memcmp(&charts[1], &original, sizeof original) == 0;
            for (int i = 0; i < 12; i++) { anchors[i].measured_order = 1; weights[i] = 0; }
            ok = ok && !apl_place_and_audit(&x, 1, 0, 1000, 0, NULL, 2) && memcmp(&charts[1], &original, sizeof original) == 0;
            for (int i = 0; i < 12; i++) weights[i] = 1;
        }
        if (!ok) { fprintf(stderr, "  asm_place selftest FAIL: supported incumbent / evidence / rollback\n"); fails++; }
        Arena_dispose(&run.arena);
    }
    {
        AsmRelation r; memset(&r, 0, sizeof r);
        r.corr_count = 12; r.rms = 1.0; r.flags = ASM_REL_SWITCHED;
        int ok = apl_placement_seam(&r, 0, 0) == APL_SWITCHED_SEAM_ANCHORS;
        r.flags |= ASM_REL_DROPPED;
        ok = ok && !apl_placement_seam(&r, 0, 0);
        r.flags = ASM_REL_SWITCHED | ASM_REL_CONTACT;
        ok = ok && !apl_placement_seam(&r, 0, 0);
        r.flags = ASM_REL_SWITCHED; r.corr_count = 5;
        ok = ok && !apl_placement_seam(&r, 0, 0);
        r.corr_count = 12; r.rms = 5.0;
        ok = ok && !apl_placement_seam(&r, 0, 0);
        r.flags = 0; r.rms = 1.0;
        ok = ok && !apl_placement_seam(&r, 0, 0);
        if (!ok) { fprintf(stderr, "  asm_place selftest FAIL: switched-seam evidence vetoes\n"); fails++; }
    }
    {
        /* A small seam continuation competes with a larger, nearer hop arc.
         * It gets first access, but a failed unchanged seam cannot starve
         * the arc, and a sub-minimum seam has no scheduling privilege. */
        uint8_t placed[3] = {1, 0, 0};
        int prior_failures = fails;
        int32_t anchors[3] = {0, 80, 6}, seams[3] = {0, 0, 6};
        int32_t tried[3] = {0, 0, 0}, ring[3] = {0, 1, 2};
        double area[3] = {1e6, 1e5, 1e3};
        int expected = APL_SEAM_FIRST && APL_SEAM_FRONTIER_FIRST == 1 ? 2 : 1;
        if (apl_next_component(3, placed, anchors, seams, tried, ring, area, 12) != expected) fails++;
        ring[2] = 1;
        expected = APL_SEAM_FIRST && APL_SEAM_FRONTIER_FIRST ? 2 : 1;
        if (apl_next_component(3, placed, anchors, seams, tried, ring, area, 12) != expected) fails++;
        placed[2] = 2; tried[2] = 6;
        if (apl_next_component(3, placed, anchors, seams, tried, ring, area, 12) != 1) fails++;
        anchors[2] = 7; seams[2] = 7;
        if (apl_next_component(3, placed, anchors, seams, tried, ring, area, 12) != expected) fails++;
        placed[2] = 0; anchors[2] = 5; seams[2] = 5;
        if (apl_next_component(3, placed, anchors, seams, tried, ring, area, 12) != 1) fails++;
        placed[1] = 1;
        if (apl_next_component(3, placed, anchors, seams, tried, ring, area, 12) != -1) fails++;
        if (fails != prior_failures) fprintf(stderr, "  asm_place selftest FAIL: seam frontier scheduling\n");
    }
    /* SEAMS FIRST (2026-09-10): 20 seam anchors (k = 0, offset 0,0) and 40 hop anchors (k = 1, offset
     * 100,0, one-turn 1,000 vox): with APL_SEAM_FIRST the first cluster is the seam cluster (kind 1,
     * weight 20) and the hops follow; without it one merged cluster (the 150-vox window holds both) */
    {
        AsmRun run; memset(&run, 0, sizeof run); run.arena = Arena_new();
        AplCtx x; memset(&x, 0, sizeof x); x.run = &run;
        AplAnchor an[60]; memset(an, 0, sizeof an);
        for (int i = 0; i < 60; i++) {
            int seam = i < 20;
            an[i].sx = 10.0 * (i % 10); an[i].sy = 5.0 * (i / 10);
            an[i].tx = an[i].sx + (seam ? 0.0 : 100.0); an[i].ty = an[i].sy;
            an[i].w = 1.0; an[i].k = seam ? 0 : 1; an[i].Tx = seam ? 0.0 : 1000.0; an[i].cluster = -1; an[i].seam_rms = seam ? 1.0 : 0.0;
        }
        x.anchors = an; x.na_seam = 20;
        double w[APL_MAX_CLUSTERS_TOTAL], u[APL_MAX_CLUSTERS_TOTAL], v[APL_MAX_CLUSTERS_TOTAL]; int kind[APL_MAX_CLUSTERS_TOTAL];
        int ncl = apl_cluster_anchors(&x, 60, w, u, v, kind);
        int ok = APL_SEAM_FIRST ? (ncl == 2 && kind[0] == 1 && fabs(w[0] - 20.0) < 1e-9 && fabs(u[0]) < 1e-9 && kind[1] == 0 && fabs(w[1] - 40.0) < 1e-9 && fabs(u[1] - 100.0) < 1e-9)
                                : (ncl == 1 && fabs(w[0] - 60.0) < 1e-9);
        if (!ok) { fprintf(stderr, "  asm_place selftest FAIL: seam-first clustering: %d clusters, kind %d w %.1f u %.1f\n", ncl, ncl ? kind[0] : -1, ncl ? w[0] : 0.0, ncl ? u[0] : 0.0); fails++; }
        Arena_dispose(&run.arena);
    }
    /* three layers of a curved sheet on cylinders of radius 100, 110, 120 (10 vox
     * apart); u runs around the cylinder, v along its axis (world z).  Component 0
     * = layers 0 and 1 (gives the measured turn vector), component 1 = layer 2,
     * left in an arbitrary (mirrored, rotated) frame to be oriented and placed. */
    AsmRun run; memset(&run, 0, sizeof run); run.arena = Arena_new();
    enum { N = 5 };
    float uv[N*N*2];
    int32_t faces[(N-1)*(N-1)*2*3];
    float xyz[3][N*N*3], nrm[3][N*N*3];
    int nf = 0;
    for (int j = 0; j < N; j++) for (int i = 0; i < N; i++) { uv[(j*N+i)*2] = 10.0f*i; uv[(j*N+i)*2+1] = 10.0f*j; }
    for (int j = 0; j < N-1; j++) for (int i = 0; i < N-1; i++) {
        int a = j*N+i, b = a+1, c = a+N, d = c+1;
        faces[nf*3] = a; faces[nf*3+1] = b; faces[nf*3+2] = d; nf++;
        faces[nf*3] = a; faces[nf*3+1] = d; faces[nf*3+2] = c; nf++;
    }
    for (int L = 0; L < 3; L++) {
        double R = 100.0 + 10.0 * L;
        for (int j = 0; j < N; j++) for (int i = 0; i < N; i++) {
            double ang = (10.0 * i) / R;
            float *p = &xyz[L][(j*N+i)*3], *n = &nrm[L][(j*N+i)*3];
            p[0] = (float)(10.0 * j);
            p[1] = (float)(R * sin(ang));
            p[2] = (float)(R * cos(ang));
            n[0] = 0.0f; n[1] = (float)sin(ang); n[2] = (float)cos(ang);
        }
    }
    for (int L = 0; L < 3; L++) {
        AsmChart c; memset(&c, 0, sizeof c);
        c.id = L; c.nv = N*N; c.nf = (size_t)nf; c.uv = uv; c.faces = faces; c.nrm = nrm[L]; c.xyz = xyz[L];
        c.area3d = 1600.0; c.placed = 1;
        c.component = (L < 2) ? 0 : 1;
        if (L < 2) { c.pose_x = 2.0 * APL_PI * (100.0 + 10.0 * L) * L; c.pose_y = 0.0; c.pose_theta = 0.0; }
        else { c.pose_x = 1234.0; c.pose_y = -77.0; c.pose_theta = 0.7; c.flags |= ASM_CHART_MIRROR; }   /* to be undone */
        AsmRun_push_chart(&run, &c);
    }
    double T0 = run.charts[1].pose_x - run.charts[0].pose_x;
    run.chart_layer_d = ARENA_CALLOC(run.arena, 3, sizeof(double));
    for (int i = 0; i < 3; i++) run.chart_layer_d[i] = 10.0;
    run.layer_d_global = 10.0;
    int spec[3][3] = { {0, 1, 1}, {1, 0, -1}, {1, 2, 1} };
    for (int pi = 0; pi < 3; pi++) {
        AsmLayerPair lp; memset(&lp, 0, sizeof lp);
        lp.a = spec[pi][0]; lp.b = spec[pi][1]; lp.side = (int8_t)spec[pi][2]; lp.count = 16; lp.d_median = 10.0;
        lp.k = 1;
        lp.hit_first = (int32_t)run.n_layer_hits; lp.hit_count = 0;
        for (int j = 0; j < N-1; j++) for (int i = 0; i < N-1; i++) {
            AsmLayerHit h; memset(&h, 0, sizeof h); h.a = lp.a; h.va = j*N+i; h.b = lp.b; h.side = lp.side; h.dist = 10.0f; h.order = 1;
            h.fb = (j*(N-1)+i)*2; h.l0 = 1.0f; h.l1 = 0.0f;
            AsmRun_push_layer_hit(&run, &h);
            lp.hit_count++;
        }
        AsmRun_push_layer(&run, &lp);
    }
    AsmPlaceOpts o; AsmPlace_default_opts(&o);
    o.min_anchors = 8;
    AsmPlaceStats st;
    AsmPlace_run(&run, &o, NULL, &st);
    const AsmChart *c2 = &run.charts[2];
    double gx, gy;
    apl_pose(c2, 0.0, 0.0, &gx, &gy);
    /* expected: chart 2 at chart 1 + T_out with T_out = T * (1 + 2 pi d / |T|) along +u
     * (only one chart has both turns, so the growth is the offset-curve rule) */
    /* Component labels are not certificates: the two reference arcs have no
     * source seam in this fixture, so all three arcs are independently
     * placed. Integrate the two turns at their mean radii, 105 and 115.
     * The five-column curvature fit has about 0.4% discretization error. */
    double ex = 2.0*APL_PI*(105.0+115.0);
    if (st.placed != 3 || c2->component != 0 || fabs(gx - ex) > 8.0 || fabs(gy) > 2.0 || fabs(c2->pose_theta) > 0.02 ||
        (c2->flags & ASM_CHART_MIRROR) || st.oriented_mirrored + st.mirrored != 1 ||
        c2->placement_state != ASM_PLACE_LAYER || c2->placement_support < 8 || c2->placement_parent < 0) {
        fprintf(stderr, "  asm_place selftest FAIL: placed %zu comp %d at (%.2f, %.2f) expected (%.2f, 0) theta %.3f mirror %d (mirrored %zu)\n",
                st.placed, c2->component, gx, gy, ex, c2->pose_theta, (c2->flags & ASM_CHART_MIRROR) != 0, st.oriented_mirrored); fails++;
    }
    {
        /* u must increase outward: the outer layer (chart 1) sits at larger u than chart 0 */
        double u0, v0, u1, v1;
        apl_pose(&run.charts[0], 0.0, 0.0, &u0, &v0);
        apl_pose(&run.charts[1], 0.0, 0.0, &u1, &v1);
        if (!(u1 > u0 + 100.0)) { fprintf(stderr, "  asm_place selftest FAIL: u does not increase outward (%.1f -> %.1f)\n", u0, u1); fails++; }
    }
    {
        double ax, ay, bx, by;
        apl_vertex_point(&run.charts[0], 0, &ax, &ay);
        apl_vertex_point(&run.charts[0], N*(N-1), &bx, &by);
        if (!(by > ay + 30.0) || fabs(bx - ax) > 2.0) { fprintf(stderr, "  asm_place selftest FAIL: v not along z (%.1f,%.1f)->(%.1f,%.1f)\n", ax, ay, bx, by); fails++; }
    }
    {
        /* Compare the persistent grid with the independent sorted-sample
         * audit while sliding the outer cylinder past the placed rim. Include
         * exact and near-exact reach boundaries, a wrong-wrap overlap, and
         * an empty destination. Both the touched count and verdict must agree. */
        const double offsets[] = {-48.0, -24.000000000001, -24.0, -23.999999999999,
                                  0.0, 23.999999999999, 24.0, 24.000000000001, 48.0, 5000.0};
        int32_t orig[3] = {0, 0, 1};
        AsmChart *cand = &run.charts[2];
        const AsmChart *cand_list[1] = {cand};
        AplRimGrid *grid = aprg_new(run.arena, 24.0);
        aprg_insert_chart(grid, &run.charts[0]);
        aprg_insert_chart(grid, &run.charts[1]);
        size_t touched = 0;
        for (size_t j = 0; j < sizeof offsets / sizeof offsets[0]; j++) {
            cand->pose_x = run.charts[1].pose_x + offsets[j];
            size_t ref_n = 0, grid_n = 0;
            double ref_bad = apl_rim_audit(&run, orig, 1, 0, 24.0, &ref_n);
            double grid_bad = aprg_probe(grid, &run, cand_list, 1, 24.0, &grid_n);
            touched += grid_n;
            if (grid_n != ref_n || fabs(grid_bad - ref_bad) > 1e-12) {
                fprintf(stderr, "  asm_place selftest FAIL: rim grid at offset %.12g: %zu/%.6f vs %zu/%.6f\n",
                        offsets[j], grid_n, grid_bad, ref_n, ref_bad); fails++;
            }
        }
        if (touched == 0) { fprintf(stderr, "  asm_place selftest FAIL: rim comparison had no support\n"); fails++; }
    }
    {
        /* A resumed sheet has authoritative double coordinates. Its old
         * intrinsic poses deliberately point elsewhere. Both the layer
         * observations and persistent audit grids must use the repaired map,
         * and neither a successful retry nor a refusal may move its context. */
        double current[2][2*N*N],saved_current[2][2*N*N];
        for(int c=0;c<2;c++){
            for(int v=0;v<N*N;v++){
                double p[2];AsmChart_point(run.charts+c,(size_t)v,p);
                current[c][2*v]=p[0]+257.125;current[c][2*v+1]=p[1]-91.75;
            }
            run.charts[c].placed_uv=current[c];run.charts[c].pose_x=-5000-100*c;
            run.charts[c].pose_y=8000;run.charts[c].pose_theta=.31;
        }
        memcpy(saved_current,current,sizeof current);
        run.charts[2].placement_state=ASM_PLACE_NONE;run.charts[2].component=1;
        AsmChart saved[3];memcpy(saved,run.charts,sizeof saved);
        uint8_t core[3]={0,0,1};AsmConflictOpts co;AsmConflict_default_opts(&co);
        size_t layers=run.n_layers;run.n_layers=0;
        int refused=AsmPlace_retry_component(&run,core,&o,&co,&st);
        int okay=!refused && !memcmp(saved,run.charts,sizeof saved) && !memcmp(saved_current,current,sizeof current);
        run.n_layers=layers;
        int accepted=AsmPlace_retry_component(&run,core,&o,&co,&st);
        double p[2];AsmChart_point(run.charts+2,0,p);
        okay=okay && accepted==1 && fabs(p[0]-(ex+257.125))<8 && fabs(p[1]+91.75)<2 &&
             !memcmp(saved,run.charts,2*sizeof(AsmChart)) && !memcmp(saved_current,current,sizeof current);
        if(!okay){fprintf(stderr,"  asm_place selftest FAIL: repaired-field retry accepted %d at %.9g,%.9g\n",accepted,p[0],p[1]);fails++;}
    }
    Arena_dispose(&run.arena);
    fails += aps_selftest();
    if (fails == 0) fprintf(stderr, "  asm_place selftest: all passed\n");
    return fails;
}
