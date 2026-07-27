#ifndef MESH_TRIM_INCLUDED
#define MESH_TRIM_INCLUDED

#include <stdint.h>
#include <stddef.h>
#include "arena.h"

/*
 * Mesh_trim_to_owned_box -- filter a mesh to keep only faces whose THREE
 * vertices ALL lie in the half-open box [lo, hi) per axis (Z, Y, X order).
 * Used in halo-overlap grid stitching: each cube's pipeline runs on a padded
 * volume and produces a mesh extending into neighboring cubes' territory;
 * this trim keeps only the faces owned by the central cube. Called with an
 * INSET box [inset, cube_D-inset) so adjacent cubes' charts do NOT touch
 * (a clean 2*inset gap the seam bridge spans); the strict all-vertex test
 * makes the inset exact (centroid-only let a vertex poke ~edge/2 past).
 *
 * Unreferenced vertices are compacted out and face indices remapped.
 *
 * If `pin_mask_in` is non-NULL it is taken to have length nv_in. When
 * `out_pin_mask` is also non-NULL the function emits a compacted pin_mask
 * (length *out_nv) indexed by the new vertex indices. Pass NULL for
 * either to opt out; both must be set to propagate pins. The trim never
 * moves vertices, so propagation is a pure index remap.
 *
 * All output arrays are arena-allocated. Returns 0 on success.
 */
int Mesh_trim_to_owned_box(Arena_T arena,
                           const float *verts_in, size_t nv_in,
                           const int32_t *faces_in, size_t nf_in,
                           const uint8_t *pin_mask_in,
                           float owned_lo, float owned_hi,
                           float **out_verts, size_t *out_nv,
                           int32_t **out_faces, size_t *out_nf,
                           uint8_t **out_pin_mask);

#endif
