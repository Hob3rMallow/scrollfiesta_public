#ifndef FOLD_CLEANUP_INCLUDED
#define FOLD_CLEANUP_INCLUDED

#include <stddef.h>
#include "../common/arena.h"
#include "../common/mesh_types.h"

/*
 * Remove "fold flaps" from Ball-Pivoting output.
 *
 * BPA occasionally lays a small triangle back over the sheet it just built:
 * two faces share an interior edge, traverse it the SAME direction (a winding
 * inconsistency) and have opposing normals (dot < 0). The surface creases flat
 * back on itself. A BFS winding-repair cannot fix this -- flipping one face
 * trades the winding error for a geometric overlap, and the conflict just
 * moves to an adjacent edge (these survive as residual same_dir edges).
 *
 * For each such fold this deletes the flap: of the two faces, the one whose
 * normal least agrees with its OTHER edge-neighbours (the low-coherence
 * outlier) is removed. Deletion is iterated to a fixpoint because removing one
 * flap can expose or heal an adjacent fold. Faces are compacted in place
 * (arena-reallocated); no vertices are added or moved.
 *
 * Run AFTER orientation and BEFORE PinholeFill_process: deleting a flap leaves
 * a small boundary notch that the pinhole pass recloses with consistent
 * winding. Counts are optional.
 *
 * Returns 0 on success.
 */
int FoldCleanup_process(Arena_T arena,
                        ComponentMesh *meshes, size_t n_meshes,
                        int max_passes,
                        size_t *out_faces_removed);

#endif
