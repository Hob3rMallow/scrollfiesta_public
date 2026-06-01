#ifndef BRIDGE_CUT_INCLUDED
#define BRIDGE_CUT_INCLUDED

#include "../common/arena.h"
#include "../common/mesh_types.h"

/*
 * BridgeCut_process — Step 2 entry point.
 *
 * Recursively detect and sever thin bridges in a multi-sheet mesh
 * using Edmonds-Karp max-flow with vertex splitting.
 *
 * arena       — arena for all allocations (save/restore internally)
 * mesh        — input component mesh
 * grid_size   — oracle UV grid size (320)
 * timeout_sec — max wall time for this component (0 = no timeout)
 * out_meshes  — output: arena-allocated array of sub-meshes
 * out_count   — output: number of sub-meshes
 *
 * Returns 0 on success. Output always contains >= 1 mesh.
 */
int BridgeCut_process(Arena_T               arena,
                      const ComponentMesh   *mesh,
                      int                    grid_size,
                      double                 timeout_sec,
                      ComponentMesh        **out_meshes,
                      size_t                *out_count);

#endif
