#ifndef BOUNDARY_ARC_WELD_INCLUDED
#define BOUNDARY_ARC_WELD_INCLUDED

#include <stddef.h>
#include <stdint.h>

#include "../common/arena.h"

/*
 * Merge bit-coincident vertices only when they are evidence for one physical
 * boundary zipper.  A merge is admitted by either of two topology tests:
 *
 *   (1) local link closure -- the two vertices are opposite-oriented boundary
 *       neighbours of one existing fan (a one-vertex slit), or
 *   (2) reciprocal arc correspondence -- unambiguous, opposite-oriented,
 *       normal-compatible anchors persist along both boundary chains for a
 *       minimum arclength.
 *
 * `merge_group` is a necessary provenance partition: vertices in different
 * pre-weld connected components never merge.  It is deliberately not used as
 * proof that two vertices belong to the same layer; the arc tests provide that
 * missing local evidence.
 */
typedef struct {
    double normal_dot_min;       /* compatible sheet normals, default 0.50 */
    double tangent_dot_max;      /* directed boundary tangents, default -0.50 */
    double ambiguity_margin;     /* best-vs-runner-up score gap, default 0.10 */
    double max_anchor_gap;       /* max arclength between anchors, default 32 */
    double min_support_length;   /* supported chain length, default 12 */
    size_t min_anchors;          /* anchors in a supported chain, default 2 */
} BoundaryArcWeldParams;

typedef struct {
    size_t boundary_vertices;
    size_t spatial_candidates;
    size_t reciprocal_pairs;
    size_t ambiguous_vertices;
    size_t arc_links;
    size_t arc_clusters;
    size_t accepted_clusters;
    size_t local_link_pairs;
    size_t accepted_pairs;
    size_t filter_tests;
    size_t filter_accepts;
    size_t chart_clusters_considered;
    size_t chart_clusters_accepted;
    size_t chart_clusters_rejected_self;
    size_t chart_clusters_rejected_nondisk;
    size_t chart_clusters_rejected_nonpath;
    size_t chart_clusters_rejected_sparse;
    size_t chart_clusters_rejected_cycle;
} BoundaryArcWeldStats;

/* Optional chart-level transaction policy.  chart_component[v] addresses
 * chart_is_disk[] and describes the CURRENT (post-split) chart.  merge_group
 * still describes the persistent PRE-split lineage.  A policy-enabled weld
 * admits only complete edge-for-edge boundary zippers between two disk charts,
 * and chooses a maximum-support forest of chart joins.  Consequently a good
 * zipper is not rolled back because an unrelated candidate elsewhere is bad,
 * while a second attachment that would close a chart cycle is rejected. */
typedef struct {
    const int32_t *chart_component;   /* [nv], -1 for unused */
    const uint8_t *chart_is_disk;     /* [n_chart_components] */
    size_t n_chart_components;
    int cross_chart_forest;
    int require_edge_zipper;
} BoundaryArcWeldPolicy;

void BoundaryArcWeld_default_params(BoundaryArcWeldParams *p);

int BoundaryArcWeld_process(Arena_T arena,
                            const float *verts, size_t nv,
                            int32_t *faces, size_t nf,
                            float eps,
                            const int32_t *merge_group,
                            const BoundaryArcWeldParams *params,
                            float **out_verts, size_t *out_nv,
                            size_t *out_nf,
                            BoundaryArcWeldStats *stats);

int BoundaryArcWeld_process_with_policy(
                            Arena_T arena,
                            const float *verts, size_t nv,
                            int32_t *faces, size_t nf,
                            float eps,
                            const int32_t *merge_group,
                            const BoundaryArcWeldParams *params,
                            const BoundaryArcWeldPolicy *policy,
                            float **out_verts, size_t *out_nv,
                            size_t *out_nf,
                            BoundaryArcWeldStats *stats);

/* Synthetic positive/negative regressions for the boundary correspondence. */
int BoundaryArcWeld_selftest(void);

#endif
