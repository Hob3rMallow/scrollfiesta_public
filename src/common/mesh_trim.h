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

/*
 * Mesh_trim_cut_to_owned_box -- cut-at-plane trim: keep the part of the mesh
 * inside [lo, hi]^3 (all axes), CLIPPING triangles that cross the box planes
 * (Sutherland-Hodgman against the 6 axis-aligned half-spaces) instead of
 * dropping them. The drop-only trim above loses up to one full edge length of
 * surface per side, which at coarse (CVT, ~10-15 vox edges) gouges edge-deep
 * notches into the seam rim and starves the seam bridge; the cut trim leaves
 * every cube's boundary EXACTLY on its inset planes, so the inter-cube gap is
 * a uniform 2*inset regardless of mesh density.
 *
 * Numerical contract:
 *   - Vertices within snap_eps of a plane are clamped ONTO it first (kills
 *     the double->float jitter that put CVT boundary sites 1e-5 outside and
 *     cost a whole coarse ring under the drop-only rule).
 *   - New cut vertices are placed exactly on the cutting plane (interpolate,
 *     then overwrite the plane-axis coordinate with the plane value).
 *   - A crossing within snap_eps of an endpoint reuses that endpoint, so no
 *     sub-eps sliver strips are minted.
 *   - Cut points on a shared input edge are SHARED between both incident
 *     faces (hashed per undirected edge + plane), so the cut boundary is
 *     watertight -- no cracks, no duplicated boundary edges.
 *
 * Unlike the drop-only trim this MOVES verts (by <= snap_eps) and CREATES
 * verts; there is no pin-mask variant. Clipped polygons are fan-triangulated;
 * degenerate fans are dropped. Unreferenced verts compacted, indices remapped.
 * out_faces_cut (nullable) receives the number of input faces that were
 * actually clipped (crossed a plane). All outputs arena-allocated. Returns 0.
 */
int Mesh_trim_cut_to_owned_box(Arena_T arena,
                               const float *verts_in, size_t nv_in,
                               const int32_t *faces_in, size_t nf_in,
                               float owned_lo, float owned_hi,
                               float snap_eps,
                               float **out_verts, size_t *out_nv,
                               int32_t **out_faces, size_t *out_nf,
                               size_t *out_faces_cut);

#endif
