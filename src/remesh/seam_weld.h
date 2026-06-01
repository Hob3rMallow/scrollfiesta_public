#ifndef SEAM_WELD_INCLUDED
#define SEAM_WELD_INCLUDED

#include <stddef.h>
#include <stdint.h>
#include "../common/arena.h"

/*
 * Seam-weld a concatenated multi-cube mesh by re-running BPA across each
 * cube-boundary seam ("treat the halo as a new BPA weld").
 *
 * The BPA mesh path does NOT pin seam vertices, so adjacent cubes never
 * produce bit-identical seam points -- a hash-join weld finds ~0 shared
 * verts and the seam stays open (two overlapping, independently
 * triangulated sheets). This routine fixes that directly:
 *
 *   1. Per-vertex normals are recomputed from the faces (area-weighted),
 *      since the final QEM mesh carries none.
 *   2. Axis-aligned cube-boundary planes (integer multiples of cube_size
 *      that have mesh within `band` on BOTH sides) are detected.
 *   3. The open boundary half-edges lying IN a seam plane (both endpoints
 *      within `band`, nearly parallel to the plane) are fed as the initial
 *      front to BallPivot_bridge; only those edge verts form the BPA
 *      sub-cloud, so the rolling ball cannot wander onto the existing
 *      surface and 3-fan an interior edge.
 *   4. BPA rolls across the gap at an escalating, capped radius; the two
 *      sides' fronts meet, and the directed-front glue keeps the bridge
 *      faces wound consistently with the existing surface. (No peelback:
 *      the boundaries are bridged directly -- see seam_weld.c.)
 *
 * No new vertices are created: bridge faces index into the SAME verts[]
 * array. The output is the original faces plus the new bridge faces, in
 * one arena-allocated array.
 *
 * verts[nv*3]  float (z,y,x), world coords (as emitted by the welder)
 * faces[nf*3]  int32, 0-based
 * cube_size    voxels per cube edge (e.g. 128)
 * rho          BPA pivot radius floor (voxels; e.g. 1.5)
 * band         half-width (vox) for selecting seam boundary edges
 *
 * Returns 0 on success. On success *out_faces / *out_nf hold the combined
 * (original + bridge) face list. *out_n_bridge (may be NULL) receives the
 * number of bridge faces appended.
 */
int SeamWeld_bridge(Arena_T arena,
                    const float *verts, size_t nv,
                    const int32_t *faces, size_t nf,
                    float cube_size, float rho, float band,
                    int32_t **out_faces, size_t *out_nf,
                    size_t *out_n_bridge);

#endif
