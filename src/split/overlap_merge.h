#ifndef OVERLAP_MERGE_INCLUDED
#define OVERLAP_MERGE_INCLUDED

#include <stddef.h>
#include <stdint.h>

#include "../common/arena.h"
#include "../common/mesh_types.h"

/*
 * Conservative local alternative to an overlap multicut.
 *
 * The multicut remains the detector and default result.  This module only
 * proposes a single-sheet replacement when the final cut has exactly two
 * labels and all repulsive edges form one small, connected event.  Faces on
 * both labels are grown by equal topological rings, projected to one local
 * best-fit plane, coalesced away from the fixed frontier, and re-triangulated
 * as the union of the two projected patches.
 *
 * This module constructs topology only.  The initial 3-D lift may still cross
 * itself while it is being fitted; overlap_sep.c must run the production
 * overlap audit on the RAW-fitted candidate before accepting the merge.
 */

typedef struct {
    int     min_rings;
    int     max_rings;
    int     local_search_rings;
    size_t  max_overlap_pairs;
    size_t  max_seed_faces;
    size_t  max_patch_faces;
    double  max_patch_fraction;
    double  max_seed_diameter;
    double  max_patch_diameter;
    double  max_plane_rms;
    double  contact_distance;
    double  merge_distance;
} OverlapMergeConfig;

typedef struct {
    const char *reason;
    int         rings;
    size_t      overlap_pairs;
    size_t      seed_faces;
    size_t      patch_faces_removed;
    size_t      patch_faces_added;
    size_t      patch_face_start;
    size_t      patch_vertices_coalesced;
    size_t      interface_edges;
    size_t      boundary_loops_before;
    size_t      boundary_loops_after;
    double      seed_diameter;
    double      patch_diameter;
    double      plane_rms;
    float       plane_normal[3];
    float       plane_center[3];
    const int32_t *exterior_source_faces;
    const uint8_t *vertex_role; /* [out.nv], OVERLAP_MERGE_VERTEX_* */
    const float   *chart_uv;    /* [out.nv*2], common repair chart */
} OverlapMergeStats;

enum {
    OVERLAP_MERGE_VERTEX_EXTERIOR = 0,
    OVERLAP_MERGE_VERTEX_ANCHOR   = 1,
    OVERLAP_MERGE_VERTEX_REPAIR   = 2
};

void OverlapMerge_default_config(OverlapMergeConfig *config);

/*
 * Returns 0 and fills `out` when a candidate clears every topological gate.
 * The caller must fit and geometrically audit it before replacing the ordinary
 * multicut split.  Returns -1 when no safe topology can be constructed.
 */
int OverlapMerge_try(Arena_T arena,
                     const ComponentMesh *mesh,
                     const int32_t *face_labels,
                     int32_t num_labels,
                     const int32_t *adj_fa,
                     const int32_t *adj_fb,
                     size_t n_adj,
                     const int32_t *ovl_fa,
                     const int32_t *ovl_fb,
                     size_t n_ovl,
                     const OverlapMergeConfig *config,
                     ComponentMesh *out,
                     OverlapMergeStats *stats);

#endif
