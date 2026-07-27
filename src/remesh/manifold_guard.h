#ifndef MANIFOLD_GUARD_INCLUDED
#define MANIFOLD_GUARD_INCLUDED

#include <stddef.h>
#include "../common/arena.h"
#include "../common/mesh_types.h"

/*
 * Final per-cube manifold guard: make every component a 2-manifold (with
 * boundary), as a safety net after trim. Three steps, each monotone toward
 * manifold (none can introduce a new non-manifold edge or pinch):
 *
 *   1. Resolve non-manifold EDGES (>2 incident faces): keep the
 *      opposite-winding pair with the largest area (a clean manifold pair, or
 *      one face if all share a winding) and DELETE the rest. Deletion only
 *      removes incidences, so it can never create a >2-face edge; it can leave
 *      a boundary nick (acceptable) and can orphan a vertex into a pinch (step 2
 *      catches it). Must run first: split_pinch_verts assumes edge-manifold input.
 *   2. Split pinch (bowtie) VERTICES via PinholeFill_split_pinches.
 *   3. Re-assert consistent winding (OrientMesh_consistent), clearing same_dir.
 *
 * The BPA Case-2 vertex-manifold guard (ball_pivot.c) prevents almost all of
 * this at the source; the guard exists to mop up trim-orphaned pinches and any
 * residue, and to make "manifold per-cube output" an invariant rather than a hope.
 *
 * Operates in place on the arena-allocated meshes. `reorient` controls step 3:
 * pass 1 to re-assert consistent winding (per-cube use, where nothing oriented
 * the trimmed mesh), or 0 to skip it (weld use, where the caller already ran
 * its own orientation pass and resolve+split preserve winding). Returns 0 on
 * success.
 */

typedef struct {
    size_t nm_edges_resolved;   /* undirected edges found with >2 faces */
    size_t faces_deleted;       /* faces removed to resolve them */
    size_t pinch_splits;        /* pinch vertices duplicated */
    size_t orient_flips;        /* faces re-wound (0 if reorient==0) */
} ManifoldGuardStats;

int ManifoldGuard_process(Arena_T arena, ComponentMesh *meshes, size_t n_meshes,
                          int reorient, ManifoldGuardStats *out);

/* Self-test: a 3-fan NM edge, a bowtie vertex, and a clean disk -> the first two
 * become manifold, the disk is untouched. Returns 0 on success (failure count). */
int ManifoldGuard_selftest(Arena_T arena);

#endif
