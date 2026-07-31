#ifndef ATLAS_BOXCUT_LAYOUT_INCLUDED
#define ATLAS_BOXCUT_LAYOUT_INCLUDED

#include <stddef.h>
#include <stdint.h>

#include "../common/arena.h"
#include "atlas_overlap_audit.h"

typedef struct {
    const int32_t *faces;
    size_t nfaces;
    size_t nvertices;
    const double *base_u;
    const double *v;
    const float *registered_u;
    const float *phi;
    const double *face_radius;
    const int32_t *face_label;
    size_t ncharts;
    const int32_t *adjacency_face0;
    const int32_t *adjacency_face1;
    const double *adjacency_weight;
    size_t intrinsic_adjacency_edges;
    size_t adjacency_edges;
    const AtlasOverlapAudit *overlap;
    const uint64_t *allowed_overlap_key;
    size_t allowed_overlap_keys;
    double radius_gate;
    double average_uv_edge_length;
    double coherence_weight_scale;
} AtlasBoxcutLayoutProblem;

typedef struct {
    int32_t chart0;
    int32_t chart1;
    size_t pairs;
    size_t positive_pairs;
    size_t negative_pairs;
    double target_median;
    double target_mad;
    double sign_agreement;
    double separation_lower;
    double qp_weight;
    int eligible;
    int selected;
    double solved_delta;
    double residual;
} AtlasBoxcutLayoutRelation;

typedef struct {
    int32_t chart0;
    int32_t chart1;
    size_t intrinsic_edges;
    size_t seam_edges;
    double source_weight;
    double qp_weight;
    double solved_delta;
    double residual;
} AtlasBoxcutLayoutCoherence;

typedef struct {
    int32_t chart0;
    int32_t chart1;
    size_t overlap_pairs;
    double chart1_left_upper;
    double chart1_right_lower;
} AtlasBoxcutLayoutDisjunction;

typedef struct {
    size_t observations;
    size_t relations;
    size_t eligible_relations;
    size_t selected_relations;
    size_t cycle_relations_ignored;
    size_t mixed_direction_relations;
    size_t collision_bounds;
    size_t collision_disjunctions;
    size_t collision_relations_rejected;
    size_t collision_outer_iterations;
    size_t collision_bounds_added;
    size_t slab_records;
    size_t slab_count;
    size_t slab_bounds;
    double slab_height;
    size_t initial_exact_overlap_pairs;
    size_t initial_cross_overlap_pairs;
    size_t final_exact_overlap_pairs;
    size_t final_cross_overlap_pairs;
    size_t radial_order_pairs_satisfied;
    size_t target_pairs_satisfied;
    size_t moved_charts;
    size_t moved_faces;
    size_t coherence_adjacency_edges;
    size_t coherence_intrinsic_edges;
    size_t coherence_seam_edges;
    size_t coherence_direct_overlap_edges_excluded;
    size_t coherence_relations;
    size_t coherence_relations_preserved;
    size_t coherence_edges_preserved;
    double radial_order_pair_fraction;
    double target_pair_fraction;
    double relation_residual_rms;
    double relation_residual_max;
    double coherence_residual_rms;
    double coherence_weighted_residual_rms;
    double coherence_residual_median;
    double coherence_residual_p95;
    double coherence_residual_max;
    double maximum_abs_shift;
    int qp_iterations;
    int qp_active_bounds;
    double qp_objective;
    double qp_stationarity;
    double qp_max_violation;
} AtlasBoxcutLayoutStats;

/* Build connected components of the chart-pair forbidden-relative-shift
 * interval union induced by one complete exact-overlap audit. Each component
 * is an independent left/right disjunction. Thresholds are expressed in the
 * original base-U chart shifts, so rows from successive audits can be merged. */
int AtlasBoxcutLayout_collect_disjunctions(
    Arena_T arena, const AtlasBoxcutLayoutProblem *problem,
    const AtlasOverlapAudit *overlap,
    AtlasBoxcutLayoutDisjunction **out_disjunction,
    size_t *out_ndisjunction);

/*
 * Solve one rigid U translation per BoxCutter chart.  Exact overlap pairs
 * supply local inner->outer order observations.  The solve keeps all coherent
 * relations as weighted rows in one convex QP, adds directed separation
 * bounds, and gives every chart an area-weighted weak zero-shift prior.
 * Original-graph boundaries cut by the multicut add zero-relative-shift
 * coherence rows unless the adjacent faces themselves overlap.  Thus a
 * locally good chart moves rigidly and trustworthy neighbours prefer to move
 * as a block.  Face-chart boundary duplication is deliberately left to the
 * caller/exporter because it is a representation operation, not a solve.
 */
int AtlasBoxcutLayout_solve(
    Arena_T arena, const AtlasBoxcutLayoutProblem *problem,
    double **out_chart_shift,
    AtlasBoxcutLayoutRelation **out_relation, size_t *out_nrelation,
    AtlasBoxcutLayoutCoherence **out_coherence, size_t *out_ncoherence,
    AtlasBoxcutLayoutDisjunction **out_disjunction,
    size_t *out_ndisjunction,
    AtlasBoxcutLayoutStats *stats);

int AtlasBoxcutLayout_selftest(void);

#endif
