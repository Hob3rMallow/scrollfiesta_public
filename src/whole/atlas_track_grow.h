#ifndef ATLAS_TRACK_GROW_INCLUDED
#define ATLAS_TRACK_GROW_INCLUDED

#include <stddef.h>
#include <stdint.h>

#include "../common/arena.h"
#include "atlas_ribbon_fit.h"

/*
 * Evidence-first atlas construction.
 *
 * The input is the raw constant-v section evidence emitted before ribbon
 * registration.  Locally compatible samples form two independent edge
 * families:
 *
 *   horizontal -- increasing section arclength along one evidence spiral;
 *   vertical   -- the same physical evidence continued to the next v row.
 *
 * Closed four-edge cells form the trusted core.  A bounded second pass can
 * complete a missing edge when three sides already grow from that core, the
 * discarded one-way evidence supports the edge, and a local affine UV->XYZ
 * fit predicts all four observed corners within tolerance.  A separate bank
 * bridge can span a longer ribbon gap only when one intact core cell on each
 * side supplies eight jointly consistent UV->XYZ observations. Point contacts
 * never create faces.  Parameterization retains the checkpoint's one GLOBAL
 * packed U gauge; graph and cross-component constraints preserve differences
 * in that gauge instead of independently flattening components or collapsing
 * them onto winding phase.  Snapped winding remains semantic metadata for the
 * organizer and machine-atlas export.
 */

typedef struct {
    double horizontal_max_gap;      /* candidate search along source U (24) */
    double horizontal_stretch;      /* maximum chord / source-U gap (1.6) */
    double horizontal_slack;        /* absolute chord allowance (2.5) */
    double horizontal_min_fraction; /* minimum chord / source-U gap (.35) */
    double vertical_search_u;       /* adjacent-row search window in U (12) */
    double vertical_radius;         /* max cross-plane displacement (8) */
    double tangent_dot_min;         /* oriented tangent agreement (.35) */
    int physical_neighbours;        /* discover H/V edges in local XYZ rather
                                     * than the checkpoint's source-U window */
    const int32_t *direct_relation_a; /* optional sorted direct chart pairs;
                                       * cross-chart neighbours must match */
    const int32_t *direct_relation_b;
    size_t direct_relation_count;
    int preserve_layout_parameterization; /* topology recovery must retain the
                                            * chart layout's complete U field */
    int vertical_winding_delta_max; /* physical V link winding tolerance */
    double vertical_chain_residual_max; /* max within-chain-pair U mismatch;
                                         * 0 disables the topology cut */
    double quad_normal_dot_min;     /* minimum two-triangle normal cosine (0) */
    int min_component_quads;        /* discard weaker cycle islands (64) */
    int gap_fill_rounds;            /* trusted-front completion rounds (2) */
    double gap_fit_rms_max;         /* local affine UV->XYZ RMS gate (1.5) */
    double gap_fit_max_max;         /* local affine UV->XYZ max gate (4) */
    double uv_bridge_max_u;         /* maximum bank-to-bank U gap; 0 disables */
    double uv_bridge_fit_rms_max;   /* eight-observation affine RMS gate (.75) */
    double uv_bridge_fit_max_max;   /* eight-observation affine max gate (2) */
    double parameter_stay_weight;   /* packed-U anchor weight during solve (.01) */
    double parameter_bridge_weight; /* unsupported evidence-edge weight (.25) */
    int parameter_iterations;       /* evidence-coordinate sweeps (160) */
} AtlasTrackGrowOptions;

typedef struct {
    int32_t id;
    size_t quads;
    size_t core_quads;
    size_t inferred_quads;
    int selected;
    size_t first_vertex;
    size_t vertices;
    size_t first_face;
    size_t faces;
    double local_u_span;
    double local_v_span;
    double global_u0;
    double global_v0;
    double atlas_u0;
    double atlas_v0;
    int atlas_band;
    int32_t winding_min;
    int32_t winding_max;
    uint64_t rank_min;
    uint64_t rank_max;
} AtlasTrackGrowComponent;

typedef struct {
    double *xyz;             /* [nv*3], source-space (z,y,x) */
    double *local_uv;        /* [nv*2], component-local view of coupled solve */
    double *global_uv;       /* [nv*2], globally coupled checkpoint/depth coords */
    double *uv;              /* [nv*2], global atlas shifted to a zero origin */
    int32_t *faces;          /* [nf*3], zero-based */
    int32_t *source_target;  /* [nv], AtlasRibbonTarget index */
    int32_t *vertex_component; /* [nv], sheet component */
    int32_t *face_component;   /* [nf], sheet component */
    uint8_t *face_inferred;    /* [nf], UV->XYZ gap-completion provenance */
    int32_t *vertex_winding; /* [nv], snapped physical winding identity */
    uint64_t *vertex_rank;   /* [nv], wind-ordered source chart rank */
    size_t nv;
    size_t nf;

    AtlasTrackGrowComponent *component;
    size_t ncomponents;

    size_t targets;
    size_t accepted_targets;
    size_t conflict_targets;
    size_t horizontal_candidates;
    size_t horizontal_relation_rejects;
    size_t horizontal_edges;
    size_t vertical_candidates;
    size_t vertical_relation_rejects;
    size_t vertical_winding_rejects;
    size_t vertical_edges;
    size_t vertical_chain_checks;
    size_t vertical_chain_pairs;
    size_t vertical_chain_cuts;
    double vertical_chain_residual_rms;
    double vertical_chain_residual_max;
    size_t vertical_global_checks;
    size_t vertical_global_cuts;
    int vertical_global_rounds;
    double vertical_global_residual_rms;
    double vertical_global_residual_max;
    size_t quad_candidates;
    size_t quad_twist_rejects;
    size_t cycle_quads;
    size_t gap_fill_candidates;
    size_t gap_fill_evidence_rejects;
    size_t gap_fill_fit_rejects;
    size_t gap_fill_topology_rejects;
    size_t gap_fill_quads;
    size_t gap_fill_pair_quads;
    size_t gap_fill_horizontal_edges;
    size_t gap_fill_vertical_edges;
    size_t uv_bridge_candidates;
    size_t uv_bridge_evidence_rejects;
    size_t uv_bridge_fit_rejects;
    size_t uv_bridge_topology_rejects;
    size_t uv_bridge_quads;
    size_t uv_bridge_horizontal_edges;
    double uv_bridge_fit_rms;
    double uv_bridge_fit_max;
    double uv_bridge_u_span_max;
    size_t selected_core_quads;
    size_t selected_inferred_quads;
    double core_fit_rms;
    double core_fit_max;
    double gap_fill_fit_rms;
    double gap_fill_fit_max;
    size_t cycle_nodes;
    size_t selected_quads;
    size_t selected_components;
    size_t selected_cycle_nodes;
    size_t selected_rows;
    size_t selected_charts;
    int32_t selected_row_min;
    int32_t selected_row_max;
    int32_t selected_column_min;
    int32_t selected_column_max;

    double source_u_span;
    double local_u_span_max;
    double atlas_u_span;
    double atlas_v_span;
    size_t atlas_placed_overlaps;
    size_t atlas_winding_overlaps;
    int atlas_bands;
    int winding_direction;
    int32_t winding_min;
    int32_t winding_max;
    double surface_area;
    double horizontal_residual_rms;
    double horizontal_residual_max;
    double vertical_residual_rms;
    double vertical_residual_max;
    size_t vertical_residual_hist[10]; /* <=.25,.5,1,2,4,8,16,32,64,>64 */
    size_t parameter_bridge_edges;
    size_t parameter_horizontal_bridges;
    size_t parameter_vertical_bridges;
    size_t parameter_point_bridges;
    size_t parameter_islands;
    double parameter_bridge_residual_rms;
    double parameter_bridge_residual_max;
    double parameter_anchor_rms;
    double parameter_anchor_max;
    double parameter_max_change;
    int parameter_iterations;
} AtlasTrackGrowResult;

void AtlasTrackGrowOptions_default(AtlasTrackGrowOptions *opts);

int AtlasTrackGrow_build(
    Arena_T arena,
    const AtlasRibbonObservationSet *evidence,
    const ScaffoldCalib *calibration,
    const AtlasTrackGrowOptions *opts,
    AtlasTrackGrowResult *out);

int AtlasTrackGrow_selftest(void);

#endif
