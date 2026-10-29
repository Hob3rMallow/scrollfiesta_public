#ifndef ASM_PLACE_INCLUDED
#define ASM_PLACE_INCLUDED

#include <stddef.h>
#include <stdint.h>

#include "asm_axis.h"
#include "asm_types.h"
#include "asm_conflict.h"

/* ============================================================================
 * asm_place.h -- placement of disconnected components by measured turn
 * vectors.
 *
 * In a box crop every outer wrap is an arc cut by the box faces, so the
 * seam graph has one component per arc plus the closed inner spiral.  The
 * spiral's own layer pairs measure the TURN VECTOR: the 2-D displacement
 * from a chart to its layer neighbour (one circumference along u, plus the
 * tilt in v).  An unplaced component sitting k layers from placed material
 * gets anchor targets  origin + k * T  from every layer hit that originates
 * on a placed chart, and is fitted rigidly (reflection allowed: components
 * carry independent handedness).  Where only one side of a chart is
 * measured, the other side follows the offset-curve rule: the perimeter of
 * a convex section grows by 2 pi d per layer of distance d, which needs no
 * axis and no radius.  Placement proceeds outward, ring by ring.
 * ==========================================================================*/

typedef struct AsmPlaceOpts {
    int    min_anchors;
    double max_rms;        /* vox, after robust reweighting */
    double rel_rms;        /* plus this times sqrt(component area): flattening drift over the extent */
    int    max_layers;     /* k allowed */
    double layer_tol;      /* |dist / d_layer - k| tolerance */
    int    min_turn_hits;  /* hits needed to measure a turn vector */
    double max_contra_ratio; /* a placement whose contradiction mass exceeds this fraction of its area is reverted */
    double max_rot_refine;   /* rotation allowed on top of the scroll-frame orientation, rad */
    double max_rim_ratio;    /* placement rejected when this fraction of its rim vertices violates isometry against placed material */
    int crosswrap_anchor_veto; /* ASM_WRAP_VETO_ANCHOR for this pile: 0 = a cross-wrap seam still anchors, 1 = never,
                                * 2 = never where the axis frame is disarmed.  PER PILE because the measurement is:
                                * 21x5x5 +1 gate (coherence 49.7%% -> 82.9%%), pherc343 -13 points of coverage,
                                * 4x5x5 untouched at 2 (its frame is armed).  `place.crosswrap_anchor_veto`. */
    double stress_band;      /* same reporting band as cleaning, after local UV repair */
    const char *diag_dir;    /* writes stage5_anchors.csv + stage5_components.csv here (NULL = none) */
    /* the scroll axis polyline: every chart's axial coordinate s (its v follows the primary's
     * v(s, u)) and radius r (one turn = 2 pi r) come from projecting its centroid; NULL = world z,
     * no radius model */
    const AsmAxis *axis;
} AsmPlaceOpts;

typedef struct AsmPlaceStats {
    int budget_exhausted;
    size_t components_in;
    size_t placed, unplaced;
    size_t mirrored;
    size_t oriented_mirrored;   /* components mirrored by the bending rule */
    size_t sweeps;
    size_t anchors_used;
    double area_in, area_placed, area_primary;
    size_t reverted_by_audit;
    size_t hop_retries, hop_retry_placed;   /* fits retried one turn either way after landing on another wrap, and those accepted */
    size_t v_pinned, v_refused;   /* placements whose v was set from the axial coordinate, and those the residual gate refused */
    double v_shift_max;           /* the largest |v correction| applied, vox */
    size_t anchors_radius;        /* anchors whose hop length came from the radius model */
    /* Components the sweep loop NEVER ATTEMPTED: unplaced with zero attempts, because
     * apl_next_component never saw min_anchors of them (anchors only come from layer pairs
     * whose other end is already placed).  This is frontier STARVATION, and it is a wholly
     * different failure from a refusal or a deferral -- measured 2026-09-19 on the 10x10x10 at
     * 1,037 components / 3.226e7 vox^2 / 38.4% of the pile, which is 69% of all unplaced area
     * and was invisible in every summary the lane printed. */
    size_t never_attempted;
    double never_attempted_area;
    int    axis_armed;            /* 0 when the arc-length gate disarmed the axis (b outside the band): the table is not this sheet's axis */
    double axis_b;                /* the primary's fitted b (0 without a model) */
    double v_model[3];            /* the primary's v = a + b s + c u */
    /* WOBBLE: the layout against v(s) */
    double wobble_primary_p50, wobble_primary_p90;   /* |v - v(s,u)| of the primary's own charts, vox */
    double placed_area_beyond30, placed_area_beyond100;   /* share of placed area further than 30 / 100 vox from the line */
    double spine_step_p50, spine_step_p90;           /* bin-to-bin step of the spine (median residual per 2,000-vox u bin) */
    size_t spine_bins;
    double v_span;                                   /* v extent of the placed layout */
    double outline_top_p2p, outline_bot_p2p;         /* peak-to-peak of the box's top / bottom crop plane's trace in v over 2,000-vox u bins, vox */
    size_t outline_bins_top, outline_bins_bot;       /* bins holding crop-plane vertices */
    double outline_pred_p2p;                         /* 2 r tan(tilt): the wave the box's crop planes make in the axis frame */
    double outline_r_mean;                           /* area-weighted mean radius of the primary's charts, vox */
    double axis_tilt_deg;                            /* the axis's mean direction off +z over the primary's z range */
    size_t cuts, cut_rels_deferred;  /* placement-gated cuts; their measured seams remain recoverable */
    double rms_p50, rms_max;
    double relevel_theta;       /* final rotation putting +v on world z over the whole placed strip */
    double growth_measured;     /* turn length growth per layer, vox */
    size_t growth_charts;       /* charts it was measured on (0 = offset-curve rule) */
    double sec;
    /* SWEEP PROFILE (2026-09-10): seconds per phase of the sweep loop and the counts behind the exact
     * speed levers (the 4x21x21 placed 5,680 sweeps in its 900-s budget, 6.6% of the area never attempted) */
    double t_loop, t_measure, t_growth, t_setup, t_count, t_mirror, t_stretch, t_collect, t_cluster, t_clusterloop, t_fit, t_audit, t_probe, t_rim, t_insert, t_cut;
    size_t remeasures, quadric_calls, quadric_cached;
    /* SEAMS FIRST and the CONFETTI LEDGER (2026-09-10) */
    size_t seam_placed, hop_placed;          /* components placed by a seam cluster / a hop cluster */
    size_t v_seam_kept;                      /* seam placements the v-axis pin left at the seam's v */
    size_t singles; double singles_area;     /* single-chart components at the end (post-cut ids) */
    size_t singles_by_status[7]; double singles_area_by_status[7];   /* never, placed, refused, reverted, primary, few_anchors, deferred */
    size_t deferred;                         /* small hop-only placements refused because the audits were silent (attempts) */
    size_t deferred_comps; double deferred_area;   /* components still deferred at the end */
    size_t last_call_placed;                       /* deferred components the last call placed under the strict gate */
    size_t placed_rms30; double placed_rms30_area;   /* placed with a fit rms over 30 vox */
    double reverted_area;
    size_t evid_comps, evid_dropped, evid_weak, evid_short; double evid_area;   /* unplaced components holding a seam to placed material */
    size_t noevid_comps; double noevid_area;
} AsmPlaceStats;

void AsmPlace_default_opts(AsmPlaceOpts *o);

/* Merges every placeable component into the primary component's frame.
 * Charts of placed components get composed poses and the primary's
 * component id; unplaced components keep their own ids.  Returns 0. */
int AsmPlace_run(AsmRun *run, const AsmPlaceOpts *o, const AsmConflictOpts *copts, AsmPlaceStats *st);

/* One unanchored complete SOURCE component, with every registered chart fixed.
 * 1 = supported proposal, 0 = refused without mutation, -1 = invalid/error.
 * The caller must still enforce native material-admission qualification. */
int AsmPlace_retry_component(AsmRun *run, const uint8_t *core, const AsmPlaceOpts *o,
                            const AsmConflictOpts *copts, AsmPlaceStats *st);

int AsmPlace_selftest(void);

#endif
