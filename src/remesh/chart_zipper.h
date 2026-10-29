#ifndef CHART_ZIPPER_H
#define CHART_ZIPPER_H

/*
 * chart_zipper.h -- atomic source-chart zipper for cube seams.
 *
 * The input is a disjoint union of certified source charts.  Each accepted
 * transaction consumes one contiguous boundary interval from two charts on
 * opposite sides of an interior cube plane and appends one monotone triangle
 * strip.  Source faces and vertices are never deleted, split, or moved.
 *
 * ChartZipper deliberately does not fill unmatched paths or repair the result.
 * Hole filling and the exact topology certificate remain separate downstream
 * stages, where they operate on real boundary loops rather than on BPA patch
 * fragments.
 */

#include <stddef.h>
#include <stdint.h>

#include "../common/arena.h"
#include "ball_pivot.h"              /* BpaBridgeGate */

typedef struct {
    float cube_size;
    float band;
    float max_cross_edge;
    float normal_dot_min;
    float min_coverage;
    float min_support_length;
    float ambiguity_ratio;
    float min_arc_ratio;
    float skinny_weight;
    float min_triangle_altitude;
    float conflict_gap;
    float conflict_parallel_angle_deg;
    size_t min_hits;
    size_t ball_results;
    int trace;
} ChartZipperParams;

typedef struct {
    size_t planes;
    size_t boundary_edges;
    size_t seam_boundary_edges;
    size_t chains;
    size_t open_chains;
    size_t closed_chains_opened;
    size_t closed_chains_rejected;
    size_t branched_chains_rejected;

    size_t candidates;
    size_t ambiguous_candidates;
    size_t duplicate_pair_candidates;
    size_t cycle_candidates;
    size_t overlap_candidates;
    size_t geometry_candidates;
    size_t conflict_transactions;
    size_t conflict_pairs;
    size_t postcompact_conflict_transactions;
    size_t postcompact_conflict_pairs;
    size_t balanced_candidates;
    size_t trimmed_candidates;
    size_t accepted_transactions;
    size_t bridge_faces;
    int embedded_certificate;
} ChartZipperStats;

/* One independently selectable zipper transaction.  face_first/face_count
 * index the array returned by ChartZipper_process_with_transactions.  The
 * ranking fields are the exact deterministic ordering keys used internally;
 * streaming weld planners can therefore merge transactions from many seam
 * shards without loading any source triangles. */
typedef struct {
    int32_t chart_a;
    int32_t chart_b;
    size_t face_first;
    size_t face_count;
    /* Required winding relationship between the stored bridge strip and each
     * unflipped source chart. +1 means the stored strip already opposes the
     * chart boundary; -1 means that strip or chart must be flipped. Streaming
     * workers populate these after the zipper returns. */
    int8_t orient_a;
    int8_t orient_b;
    uint16_t flags;
    double full_coverage;
    double span_coverage;
    double mean_gap;
    double support;
    double score;
} ChartZipperTransaction;

/* A lightweight same-sheet hypothesis produced before any zipper strip is
 * triangulated.  Samples are ordered boundary-vertex correspondences and use
 * the caller's input vertex indices.  They are sufficient to measure carried
 * material-coordinate offsets as a function of v without creating weld
 * geometry. */
typedef struct {
    int32_t vertex_a;
    int32_t vertex_b;
} ChartZipperMatchSample;

typedef struct {
    int32_t chart_a;
    int32_t chart_b;
    int32_t port_a;
    int32_t port_b;
    size_t interval_a_first;
    size_t interval_a_last;
    size_t interval_b_first;
    size_t interval_b_last;
    size_t sample_first;
    size_t sample_count;
    double full_coverage;
    double span_coverage;
    double mean_gap;
    double support;
    double score;
    double order_coverage;
} ChartZipperMatch;

enum {
    /* Member of the established pair-local topology/geometry-certified set.
     * Streaming planners retain this subset for physical materialization while
     * using every transaction as a possible relation-graph edge. */
    CHART_ZIPPER_TRANSACTION_PHYSICAL = 1u
};

void ChartZipper_default_params(ChartZipperParams *params);

int ChartZipper_process(Arena_T arena,
                        const float *verts, size_t nv,
                        const int32_t *faces, size_t nf,
                        const int32_t *chart_component,
                        size_t n_chart_components,
                        const BpaBridgeGate *gate,
                        const ChartZipperParams *params,
                        int32_t **out_faces, size_t *out_nf,
                        size_t *out_n_bridge,
                        ChartZipperStats *stats);

/* Transaction-reporting variant used by the on-disk streaming weld.  The
 * report is arena allocated and contains only strips that passed the same
 * topology and exact embedded-geometry gates as ChartZipper_process. */
int ChartZipper_process_with_transactions(
                        Arena_T arena,
                        const float *verts, size_t nv,
                        const int32_t *faces, size_t nf,
                        const int32_t *chart_component,
                        size_t n_chart_components,
                        const BpaBridgeGate *gate,
                        const ChartZipperParams *params,
                        int32_t **out_faces, size_t *out_nf,
                        size_t *out_n_bridge,
                        ChartZipperStats *stats,
                        ChartZipperTransaction **out_transactions,
                        size_t *out_n_transactions);

/* Candidate-graph variant used by streaming seam workers.  Unlike the
 * materializing operation above, this reports every geometrically constructible
 * boundary transaction before chart-pair, graph-cycle, shared-boundary-port,
 * or inter-transaction conflict selection.  Source faces remain the prefix of
 * out_faces and every transaction owns one contiguous appended face range.
 * The caller must classify conflicts and choose a compatible subset before
 * materializing these alternatives together. */
int ChartZipper_enumerate_transactions(
                        Arena_T arena,
                        const float *verts, size_t nv,
                        const int32_t *faces, size_t nf,
                        const int32_t *chart_component,
                        size_t n_chart_components,
                        const BpaBridgeGate *gate,
                        const ChartZipperParams *params,
                        int32_t **out_faces, size_t *out_nf,
                        size_t *out_n_bridge,
                        ChartZipperStats *stats,
                        ChartZipperTransaction **out_transactions,
                        size_t *out_n_transactions);

/* Match-only boundary operation.  This stops after ordered boundary-chain
 * discovery and reciprocal geometric scoring.  It does not run zipper DP,
 * append faces, audit triangle intersections, or modify topology. */
int ChartZipper_match(
                        Arena_T arena,
                        const float *verts, size_t nv,
                        const int32_t *faces, size_t nf,
                        const int32_t *chart_component,
                        size_t n_chart_components,
                        const BpaBridgeGate *gate,
                        const ChartZipperParams *params,
                        ChartZipperMatch **out_matches,
                        size_t *out_n_matches,
                        ChartZipperMatchSample **out_samples,
                        size_t *out_n_samples,
                        ChartZipperStats *stats);

/* Native fixtures include a two-chart join, a distant non-match, and the
 * important one-large-chart/two-micro-chart absorption case. */
int ChartZipper_selftest(void);

#endif /* CHART_ZIPPER_H */
