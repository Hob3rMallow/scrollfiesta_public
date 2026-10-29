#ifndef CHART_BRIDGE_FOREST_INCLUDED
#define CHART_BRIDGE_FOREST_INCLUDED

#include <stddef.h>
#include <stdint.h>

#include "../common/arena.h"

typedef struct {
    size_t input_charts;
    size_t bridge_faces_in;
    size_t bridge_faces_kept;
    size_t patches;
    size_t accepted_patches;
    size_t accepted_join_patches;
    size_t accepted_growth_patches;
    size_t accepted_hole_patches;
    size_t promotion_rounds;
    size_t promoted_patches;
    size_t promoted_faces;
    size_t accepted_source_shoulder_faces;
    size_t rejected_self_or_mixed;
    size_t rejected_patch_topology;
    size_t rejected_attachment;
    size_t rejected_cycle;
    size_t rejected_edge_conflict;
    size_t rejected_orientation;
} ChartBridgeForestStats;

typedef struct {
    size_t candidate_patches;
    size_t constraints_used;
    size_t constraints_conflicted;
    size_t charts_flipped;
    size_t faces_flipped;
} ChartBridgeOrientationStats;

/* Reconcile the global winding of disconnected source charts from complete raw
 * bridge patches before transaction selection.  A bridge patch supplies a
 * parity constraint between its two chart-boundary paths.  A maximum-area
 * consistent constraint forest is solved, then the lower-face-count side of
 * each parity tree is flipped in the source prefix.  Connectivity and vertex
 * identity never change; ChartBridgeForest_filter still performs every live
 * manifold and disk-topology admission check. */
int ChartBridgeForest_orient_source(Arena_T arena,
                                    const float *verts, size_t nv,
                                    int32_t *faces, size_t nf,
                                    size_t input_nf,
                                    const int32_t *vertex_chart,
                                    size_t n_charts,
                                    ChartBridgeOrientationStats *stats);

/* Filter the appended [input_nf,nf) BPA bridge faces as transactions between
 * original input charts.  A bridge patch must itself be a disk, attach to each
 * source chart along one boundary interval, and join two different components
 * of a maximum-area chart forest.  Adjacent patches then grow that forest only
 * across a single connected boundary interval, retaining maximal safe seam
 * coverage without adding a second attachment/handle.
 *
 * Selection is staged.  After safe chart joins, rejected raw faces are
 * reclassified against the coarser chart forest.  Thus a three-chart corner
 * face whose first two charts have already been safely joined becomes an
 * ordinary two-chart face on the next round; it is never admitted as a broad
 * three-way complex.  Every promoted patch goes through the same live
 * edge/link/orientation/disk checks.  Input faces are always retained. */
int ChartBridgeForest_filter(Arena_T arena,
                             const float *verts, size_t nv,
                             const int32_t *faces, size_t nf,
                             size_t input_nf,
                             const int32_t *vertex_chart,
                             size_t n_charts,
                             int32_t **out_faces, size_t *out_nf,
                             ChartBridgeForestStats *stats);

/* As above, with an explicit bridge-relative range of source-to-bridge
 * shoulder faces.  These faces were created while repairing the combined raw
 * source/bridge fan geometry.  A one-chart face in this range may inherit the
 * unique adjacent two-chart pair and is then judged as part of that complete
 * patch transaction.  Untagged one-chart BPA faces remain rejected. */
int ChartBridgeForest_filter_with_support(
                             Arena_T arena,
                             const float *verts, size_t nv,
                             const int32_t *faces, size_t nf,
                             size_t input_nf,
                             const int32_t *vertex_chart,
                             size_t n_charts,
                             size_t support_bridge_begin,
                             size_t support_bridge_end,
                             int32_t **out_faces, size_t *out_nf,
                             ChartBridgeForestStats *stats);

int ChartBridgeForest_selftest(void);

#endif
