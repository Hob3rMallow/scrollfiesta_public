#ifndef NORMAL_ORIENT_INCLUDED
#define NORMAL_ORIENT_INCLUDED

#include <stddef.h>
#include "../common/arena.h"

/*
 * Consistent tangent-plane orientation — Hoppe, DeRose, Duchamp, McDonald,
 * Stuetzle, "Surface Reconstruction from Unorganized Points", SIGGRAPH 1992,
 * Section 3.3.
 *
 * Reorients per-point normals IN PLACE so that neighbouring normals agree in
 * sign. The PCA normal AXIS produced by MLS (mls_project.c) is already correct;
 * only its SIGN is wrong — mls_project picks the sign per-point by a fixed
 * (1,1,1)-hemisphere rule, which flips wherever a curved sheet's normal crosses
 * the plane perpendicular to (1,1,1). BPA's per-vertex normal gate depends on
 * that sign, so the flips strand pivots (the pinholes / large strips).
 *
 * Method (Hoppe §3.3):
 *   1. Riemannian graph: connect each point to neighbours within `radius`.
 *   2. Edge weight 1 - |n_i . n_j|  (small where normals are near-parallel).
 *   3. Minimum spanning tree of that graph (favours propagating across low
 *      curvature, deferring ambiguous high-curvature edges).
 *   4. Seed: the extreme point (max z), normal forced toward +z; then DFS/BFS
 *      the MST, flipping any child whose normal opposes its parent.
 * Each graph component is oriented independently (seeded from its own extreme).
 *
 * PRECONDITION: a single-layer surface. On a doubled (recto/verso) envelope the
 * |.| weight treats the two anti-parallel faces as aligned and the MST bridges
 * them — not handled here. The LOP output is single-layer (verified: local
 * thickness ~0.004 vox), so this holds for our pipeline.
 *
 *   verts[nv*3]    float (z, y, x), read only
 *   normals[nv*3]  float, reoriented in place
 *   radius         neighbour ball radius (voxels); ~3x point spacing connects
 *                  the graph
 *   out_flips      (optional) number of normals whose sign was flipped
 *
 * Uses `arena` only for transient scratch (KD-tree, union-find), restored
 * before return. Returns 0 on success.
 */
int NormalOrient_consistent(Arena_T arena,
                            const float *verts,
                            float *normals,
                            size_t nv,
                            float radius,
                            size_t *out_flips);

#endif
