#ifndef ATLAS_BOXCUT_INCLUDED
#define ATLAS_BOXCUT_INCLUDED

#include <stddef.h>
#include <stdint.h>

#include "../common/arena.h"
#include "atlas_overlap_audit.h"
#include "atlas_seam_audit.h"

typedef struct {
    const float *vertices;
    size_t nvertices;
    const int32_t *faces;
    size_t nfaces;
    const double *u;
    const double *v;
    const float *axis_point;
    const float *axis_dir;
    const AtlasOverlapAudit *overlap;
    /* Optional exact-overlap-pair whitelist, aligned with overlap->pairs.
     * Only affirmative external delamination classification belongs here. */
    const uint8_t *allowed_overlap_pair;
    const AtlasSeamAudit *seam;
    const uint8_t *eligible_seam_bundle;
    double radius_gate;
} AtlasBoxcutProblem;

typedef struct {
    size_t faces;
    size_t intrinsic_adjacency_edges;
    size_t seam_adjacency_edges;
    size_t original_graph_edges;
    size_t exact_overlap_pairs;
    size_t radial_repulsive_edges;
    size_t near_radius_pairs_ignored;
    size_t allowed_delamination_pairs;
    size_t repulsive_edges_cut;
    size_t repulsive_edges_joined;
    size_t charts;
    double average_uv_edge_length;
    double radius_gate;
    /* Arena-owned original graph, retained so downstream scroll placement can
     * grow coherent local support patches without rebuilding topology. */
    int32_t *adjacency_face0;
    int32_t *adjacency_face1;
    double *adjacency_weight;
} AtlasBoxcutStats;

/*
 * BoxCutter-style chart partition on an existing UV atlas. Faces are graph
 * nodes; shared mesh edges and verified cross-cube seam edges are attractive
 * original-graph edges, while exact UV intersections whose 3-D radii differ
 * by more than radius_gate are hard repulsive lifted edges.  A zero
 * radius_gate is strict mode and makes every supplied positive-area overlap
 * repulsive, including pairs with equal face-mean radius.  The returned labels
 * cover every input face and never mutate or delete mesh data.
 */
int AtlasBoxcut_partition(Arena_T arena,
                          const AtlasBoxcutProblem *problem,
                          int32_t **out_face_label,
                          AtlasBoxcutStats *stats);

#endif
