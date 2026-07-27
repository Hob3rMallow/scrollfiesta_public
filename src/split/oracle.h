#ifndef ORACLE_INCLUDED
#define ORACLE_INCLUDED

#include "../common/arena.h"
#include "../common/mesh_types.h"

/*
 * Oracle_count_sheets — Step 1 entry point.
 *
 * Rasterizes the component mesh into a temporary volume, casts rays along
 * the PCA normal on a UV grid, and counts sheets by gap detection.
 *
 * Returns sheet count >= 1. Returns 1 for degenerate inputs.
 * All scratch memory arena-allocated; arena state unchanged on return.
 */
int Oracle_count_sheets(Arena_T              arena,
                        const ComponentMesh *mesh,
                        int                  grid_size);

#endif
