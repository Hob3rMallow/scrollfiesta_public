#ifndef MESH_MANIFOLD_INCLUDED
#define MESH_MANIFOLD_INCLUDED

#include <stddef.h>
#include <stdint.h>
#include "arena.h"

/*
 * Combinatorial 2-manifold audit for a triangle mesh. Pure topology -- only the
 * face index array and the vertex count matter; coordinates are irrelevant.
 * Shared by the manifold_guard repair pass and the VESUVIUS_DEBUG pipeline
 * assertion. (The standalone src/tools/manifold_check.c tool keeps its own
 * dependency-free copy; both are locked to the same five cases by selftest.)
 *
 *   EDGE audit: bucket each undirected edge by incident-face count --
 *       1 = boundary, 2 = manifold pair (opposite winding) or same_dir
 *       (same winding: a doubled/flipped face), >=3 = non-manifold edge.
 *   VERTEX audit: a vertex whose incident triangles form >=2 edge-disconnected
 *       fans is a pinch (bowtie). Detected by per-vertex union-find over
 *       incident faces joined through shared vertex-incident edges.
 *
 * MANIFOLD iff nm_edges == 0 && nm_verts == 0 (boundary is allowed).
 */

typedef struct {
    size_t nv, nf;
    size_t boundary_edges;   /* edges with exactly 1 incident face */
    size_t manifold_edges;   /* 2 faces, opposite winding (a good pair)   */
    size_t same_dir_edges;   /* 2 faces, SAME winding (doubled/flipped)   */
    size_t nm_edges;         /* >=3 faces on one edge (non-manifold edge) */
    size_t nm_verts;         /* pinch / bowtie vertices                   */
} MeshManifoldStats;

/* Audit faces[0..3*nf) (indices in [0,nv)). Uses `arena` for scratch only and
 * restores it before returning (leaves no residue). */
MeshManifoldStats MeshManifold_audit(Arena_T arena, size_t nv,
                                     const int32_t *faces, size_t nf);

static inline int MeshManifold_ok(const MeshManifoldStats *s) {
    return s->nm_edges == 0 && s->nm_verts == 0;
}

/* Self-test (5 cases: closed tetra, open quad, 3-fan NM edge, bowtie vertex,
 * doubled face). Returns 0 on success, non-zero = failure count. */
int MeshManifold_selftest(Arena_T arena);

#endif
