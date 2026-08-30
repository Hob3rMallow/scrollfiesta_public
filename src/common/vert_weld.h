#ifndef VERT_WELD_INCLUDED
#define VERT_WELD_INCLUDED

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>
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
 *   guard_orient     : when true AND nf>0, refuse to merge a coincident
 *                      pair whose area-weighted winding normals oppose
 *                      (dot < 0) — i.e. fusing opposite-facing (recto vs
 *                      verso) surfaces, which would create a non-orientable
 *                      `same_dir` edge that a winding-repair BFS cannot fix.
 *                      Assumes the input is consistently wound (so an
 *                      opposing dot is exactly a would-be same_dir fusion);
 *                      a degenerate/near-zero normal fails OPEN (merges).
 *                      Ignored on the pure point-dedup path (nf==0).
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
                bool guard_orient,
                float **out_verts, size_t *out_nv,
                float **out_normals);

/* Component-scoped variant.  A candidate pair is mergeable only when
 * merge_group[i] == merge_group[j].  This prevents two already-disconnected,
 * same-facing sheets from being joined.  It is not by itself sufficient for a
 * connected folded/spiral sheet; use Weld_verts_filtered to add a local-pair
 * guard in that case.  Pass NULL to recover Weld_verts semantics. */
void Weld_verts_grouped(Arena_T arena,
                        const float *verts, size_t nv,
                        const float *in_normals,
                        int32_t *faces, size_t nf, size_t *out_nf,
                        float eps,
                        bool guard_orient,
                        const int32_t *merge_group,
                        float **out_verts, size_t *out_nv,
                        float **out_normals);

/* Optional pair-level refinement of Weld_verts_grouped.  The callback is
 * evaluated after the group test and must return true for a spatial candidate
 * to be merged.  It is called once per unordered pair (i < j).  This lets a
 * caller impose a local topological condition which cannot be represented by
 * a single component label. */
typedef bool (*WeldPairFilter)(size_t i, size_t j, void *context);

void Weld_verts_filtered(Arena_T arena,
                         const float *verts, size_t nv,
                         const float *in_normals,
                         int32_t *faces, size_t nf, size_t *out_nf,
                         float eps,
                         bool guard_orient,
                         const int32_t *merge_group,
                         WeldPairFilter pair_filter,
                         void *pair_filter_context,
                         float **out_verts, size_t *out_nv,
                         float **out_normals);

#endif
