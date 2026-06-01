#ifndef MLS_PROJECT_INCLUDED
#define MLS_PROJECT_INCLUDED

#include <stddef.h>
#include <stdint.h>
#include "arena.h"

/*
 * MLS_project_verts — Locally Optimal Projection (LOP, μ=0) variant of
 * MLS-midpoint surface fitting (Levin 1998/2003; Lipman, Cohen-Or, Levin,
 * Tal-Ezer SIGGRAPH 2007).
 *
 * For each vertex v_i:
 *     c_i = Σ_j θ(|v_i − v_j|) · v_j  /  Σ_j θ(|v_i − v_j|)
 *     n_i = smallest-eigenvector of the same weighted covariance
 *     v_i' = c_i
 *
 * θ is the Wendland C² compact-support kernel
 *     θ(r) = (1 − r/R)^4 · (4r/R + 1)   if r < R, else 0
 *
 * Halo determinism contract:
 *   - Spatial-hash cells are computed as floor((coord + cell_origin) / R).
 *     Two cubes that pass the SAME `cell_origin` in some shared frame
 *     (e.g. world coords) get cell grids that ALIGN across cubes.
 *     If determinism doesn't matter, pass {0,0,0}.
 *   - Within each cell, entries are sorted by position (z,y,x) — not by
 *     vert index — so the accumulation order at a query is identical
 *     across cubes as long as both cubes see bit-identical positions
 *     in the query's R-neighbourhood.
 *
 * Inputs / outputs may alias (in-place projection is supported).
 *
 * `out_verts`   : [nv*3] float, projected positions (z,y,x).
 * `out_normals` : [nv*3] float, unit normals from the weighted covariance.
 *                 May be NULL if normals are not wanted.
 *
 * For vertices with fewer than MLS_MIN_NEIGHBOURS neighbours within
 * `radius_vox`, the original position and a zero normal are written —
 * the welding pass that follows will treat such verts normally.
 *
 * Allocates scratch via Arena_save/restore on `arena`; out arrays are
 * caller-provided (allocate with ARENA_ALLOC of size nv*3*sizeof(float)
 * before calling).
 */
void MLS_project_verts(Arena_T arena,
                       const float *verts, size_t nv,
                       float radius_vox,
                       const float cell_origin[3],
                       float *out_verts,
                       float *out_normals);

#endif
