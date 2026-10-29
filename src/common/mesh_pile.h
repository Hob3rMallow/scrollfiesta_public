#ifndef MESH_PILE_INCLUDED
#define MESH_PILE_INCLUDED

#include <stddef.h>

/* ============================================================================
 * mesh_pile.h -- selection of a per-cube mesh PILE from a dump tree.
 *
 * The pipeline's per-cube dumps live at
 *   <dump>/<cube_id>/<cube_id>_all_obj/<cube_id>_step12_final_all.vmesh
 * with cube_id = z#####_y#####_x##### (the cube's source-voxel origin, which
 * doubles as its world offset).  Every consumer of a pile (the quadribbon
 * lane, the axis tracker, the block driver) must select the SAME files in
 * the SAME order, or the shared-cube registration between tiles and the
 * axis curve the tiles are unwrapped about stop being exact.  This module
 * is that one selector: recursive scan, cube-id parse, optional origin
 * filter, dedup by cube id, deterministic (z, y, x) order.
 * ==========================================================================*/

#define MESH_PILE_MAX_PATH 1024

typedef struct {
    char path[MESH_PILE_MAX_PATH];
    char cube_id[24];      /* "" for a plain world-frame vmesh without an id */
    long oz, oy, ox;       /* cube origin (world offset), 0 without an id */
    int  has_id;
} MeshPileEntry;

/* Origin filter: nonzero keeps the cube.  ctx is the caller's. */
typedef int (*MeshPile_Filter)(void *ctx, long oz, long oy, long ox);

/* Parse "z#####_y#####_x#####" at the head of name.  Returns 0 and fills
 * id (NUL-terminated, 21 bytes used) and the origin; -1 if name is not a
 * cube id. */
int MeshPile_parse_cube_id(const char *name, char id[24],
                           long *oz, long *oy, long *ox);

/* Scan dir recursively (depth <= 6).  At depth 0 every *.vmesh is a
 * candidate (plain world-frame meshes allowed); below depth 0 only
 * *_final_all.vmesh per-cube dumps are.  Cube-id files failing the filter
 * or already present (dedup by id) are skipped.  Entries land in pile[]
 * (capacity cap) sorted by (oz, oy, ox, path).  Returns 0, or -1 when the
 * directory cannot be read or the pile overflows cap. */
int MeshPile_scan(const char *dir, MeshPileEntry *pile, size_t cap,
                  size_t *n, MeshPile_Filter filter, void *ctx);

/* Sort pile[0..n) by (oz, oy, ox, path): the streaming order. */
void MeshPile_sort(MeshPileEntry *pile, size_t n);

int MeshPile_selftest(void);

#endif
