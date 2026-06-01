#ifndef VERT_WELD_INCLUDED
#define VERT_WELD_INCLUDED

#include <stddef.h>
#include <stdint.h>
#include "arena.h"

/*
 * Weld_verts — merge vertices that lie within `eps` of each other
 * (transitively), remap face indices, and drop now-degenerate triangles
 * (any face whose remapped indices collide).
 *
 * Determinism: union-find with "smaller index becomes root" rule, plus
 * a sorted spatial hash, gives bit-identical output for bit-identical
 * input. Two cubes that share the same MC + MLS verts in their overlap
 * region produce identical welded positions there, provided no weld
 * group chains into the non-overlap area — which can't happen as long
 * as `eps` is much smaller than the halo width (the per-vertex weld
 * neighborhood is at most ~3·eps across).
 *
 * Inputs / outputs:
 *   verts[nv*3]      : input positions (read-only).
 *   in_normals[nv*3] : optional input normals; if non-NULL, averaged
 *                      within each weld group and renormalized.
 *   faces[nf*3]      : in/out face list. Remapped in place, then
 *                      compacted to drop degenerate triangles.
 *                      Capacity must be at least `nf*3` int32_t.
 *   eps              : merge tolerance (Euclidean voxels).
 *
 *   *out_verts       : newly arena-allocated [out_nv*3] floats
 *                      (centroid of each weld group, in input order
 *                      of the lowest-indexed group member).
 *   *out_nv          : merged vertex count.
 *   *out_normals     : optional; if non-NULL, newly arena-allocated
 *                      [out_nv*3] floats with averaged unit normals.
 *   *out_nf          : face count after dropping degenerates.
 */
void Weld_verts(Arena_T arena,
                const float *verts, size_t nv,
                const float *in_normals,
                int32_t *faces, size_t nf, size_t *out_nf,
                float eps,
                float **out_verts, size_t *out_nv,
                float **out_normals);

#endif
