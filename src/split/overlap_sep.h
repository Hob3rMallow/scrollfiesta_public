#ifndef OVERLAP_SEP_INCLUDED
#define OVERLAP_SEP_INCLUDED

#include "../common/arena.h"
#include "../common/mesh_types.h"

/*
 * OverlapSep_process -- Step 3 entry point (Overlap-Detection Multicut).
 *
 * Separate overlapping sheets by projecting to 2D via PCA normal,
 * detecting triangle overlaps, and solving a lifted multicut to find
 * minimal-length cuts separating overlapping regions (Limper 2018 S4.1).
 *
 * Parameters:
 *   arena       - arena for all allocations
 *   mesh        - input component mesh
 *   sheet_count - expected number of sheets from oracle (0 = auto-detect)
 *   n_threads   - thread count (reserved for future parallelization)
 *   timeout_sec - max wall time for this component
 *   out_meshes  - output: array of resulting sub-meshes
 *   out_count   - output: number of resulting sub-meshes
 *
 * Returns 0 on success (>= 2 sub-meshes produced).
 * Returns -1 when separation is not needed or fails
 * (timeout, degenerate, no overlaps, single-sheet). In that case,
 * *out_meshes points to a copy of the original mesh,
 * *out_count = 1.
 */
int OverlapSep_process(Arena_T              arena,
                       const ComponentMesh  *mesh,
                       int                   sheet_count,
                       int                   n_threads,
                       double                timeout_sec,
                       ComponentMesh       **out_meshes,
                       size_t               *out_count);

/* Set a directory for debug OBJ dumps at each internal phase.
 * Pass NULL to disable (default). When set, overlap_sep writes:
 *   {dir}/overlap_input.obj          - input mesh
 *   {dir}/overlap_proj2d.obj         - 2D projection (z=0)
 *   {dir}/overlap_face_labels.obj    - face label coloring
 *   {dir}/overlap_vert_labels.obj    - vertex label coloring
 *   {dir}/overlap_piece_{i}.obj      - each output piece
 */
void OverlapSep_set_debug_dir(const char *dir);

#endif
