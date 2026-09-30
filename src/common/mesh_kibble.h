#ifndef MESH_KIBBLE_INCLUDED
#define MESH_KIBBLE_INCLUDED

#include "mesh_bin.h"

typedef struct {
    size_t components, removed_components, orphan_vertices;
    size_t input_vertices, input_faces, kept_vertices, kept_faces;
    double input_area, removed_area;
} MeshKibbleStats;

/* Drop disconnected components with physical triangle area < min_area.
 * Keep EVERY component meeting the cutoff, without welding or ranking sheets.
 * Drop unreferenced vertices as well. Preserve surviving vertex/face order and
 * UV exactly. Output is arena-owned; input is never modified. min_area == 0
 * explicitly disables cleanup and returns a shallow copy of input.
 *
 * For a cube pile, apply independently to each original cube, before any axis
 * warp or subregion assembly. Absolute area makes shared cubes identical in
 * overlapping canonical blocks. Small clipped boundary patches can be removed;
 * callers must record the policy and removed area rather than infer semantics.
 */
int MeshKibble_filter(Arena_T arena, const MeshBinData *in, double min_area,
                       MeshBinData *out, MeshKibbleStats *stats);
int MeshKibble_selftest(void);

#endif
