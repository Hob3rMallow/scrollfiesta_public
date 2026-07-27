#ifndef ORIENT_MESH_INCLUDED
#define ORIENT_MESH_INCLUDED

#include <stddef.h>
#include <stdint.h>
#include "../common/arena.h"

/*
 * Deterministic per-component winding orientation for a triangle mesh.
 *
 * Ball-Pivoting assigns each triangle's winding from its seed vertex's
 * MLS normal sign. Those normals carry the fixed (1,1,1)-hemisphere sign
 * bias of mls_project.c, which flips wherever the true surface normal
 * crosses perpendicular to (1,1,1) — i.e. across a wrapping scroll sheet.
 * The result is an edge-manifold mesh whose interior edges are ~6% wound
 * the SAME direction by their two faces (inconsistent winding).
 *
 * This pass fixes that WITHOUT relying on the biased normals for winding:
 * a BFS over the dual graph (faces adjacent across exactly-2-face edges)
 * recomputes each face's shared-edge direction from the LIVE (already-
 * flipped) `faces` array at visit time, and flips any neighbour that winds
 * a shared edge the same way as its parent. Reading live windings is what
 * makes it correct — a one-shot precomputed-edge BFS goes stale the moment
 * a parent face is itself flipped (the bug in grid_weld's repair_winding).
 *
 * After a component is internally consistent there remain two global
 * orientations; if `normals` is non-NULL, each component's global sign is
 * anchored by the majority of dot(face_geom_normal, avg_vertex_normal), so
 * the kept orientation aligns with the MLS normal field. The anchor is a
 * world-coordinate property, hence identical across neighbouring cubes that
 * share the same geometry (needed for cross-cube seam agreement).
 *
 * `faces` is modified IN PLACE (winding flips swap entries [1] and [2]).
 * Non-manifold edges (>2 incident faces) are not crossed by the BFS, so
 * odd-cycle / bowtie defects survive as a small residual same-direction
 * count — those are cleared upstream by pinhole splitting, not here.
 *
 * Inputs:
 *   verts[nv*3]    float (z, y, x) — read for the normal anchor only
 *   nv             vertex count
 *   normals[nv*3]  per-vertex unit normals (z, y, x); NULL = skip anchor
 *   faces[nf*3]    int32 0-based indices, MODIFIED in place
 *   nf             face count
 *
 * Outputs (optional, may be NULL):
 *   *out_flips             total face winding flips applied
 *   *out_components        connected components over manifold edges
 *   *out_residual_same_dir same-direction interior edges remaining (==0 for
 *                          a fully consistent orientable manifold)
 *
 * Returns 0 on success.
 */
int OrientMesh_consistent(Arena_T arena,
                          const float *verts, size_t nv,
                          const float *normals,
                          int32_t *faces, size_t nf,
                          size_t *out_flips,
                          size_t *out_components,
                          size_t *out_residual_same_dir);

#endif
