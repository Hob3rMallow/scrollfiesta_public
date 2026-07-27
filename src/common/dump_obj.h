#ifndef DUMP_OBJ_INCLUDED
#define DUMP_OBJ_INCLUDED

#include "arena.h"
#include "mesh_types.h"
#include <stddef.h>

/* Extract cube ID from input path: basename without extension.
 * e.g. "../nnunet-preds/102536988.tif" -> "102536988" */
void DumpObj_extract_cube_id(const char *input_path, char *out, size_t out_sz);

/* Create a directory (like mkdir -p for one level).
 * Returns 0 on success or if dir already exists, -1 on error. */
int DumpObj_ensure_dir(const char *path);

/* Write individual component OBJs + combined all.obj.
 * Arena used for scratch allocation (freed via save/restore). */
void DumpObj_write_meshes(Arena_T arena, const char *dir,
                          const char *cube_id, const char *step_tag,
                          const ComponentMesh *meshes, size_t count);

/* Same but for an array of ComponentMesh pointers (the accumulator). */
void DumpObj_write_mesh_ptrs(Arena_T arena, const char *dir,
                             const char *cube_id, const char *step_tag,
                             const ComponentMesh **meshes, size_t count);

/* Single-mesh write with cube-origin shift. For diagnostic dumps inside
 * step modules that don't go through DumpObj_write_meshes. Verts are
 * shifted by the origin parsed from cube_id (z{vz}_y{vy}_x{vx}) so the
 * dumped OBJ sits at its final welded world position. If cube_id is NULL
 * or doesn't match the format, the shift is zero and verts are written
 * as-is. Arena used for the shifted-verts scratch allocation. */
void DumpObj_write_one_world(Arena_T arena, const char *path,
                             const char *cube_id,
                             const float *verts, size_t nv,
                             const int32_t *faces, size_t nf,
                             const float color[3]);

/* Write an oriented point cloud (verts + per-vertex normals, no faces).
 * Emits "v z y x" then "vn nz ny nx" lines. Verts are world-shifted via
 * the cube_id origin (same convention as DumpObj_write_one_world).
 * `normals` may be NULL — then only `v` lines are written. */
void DumpObj_write_points_world(Arena_T arena, const char *path,
                                const char *cube_id,
                                const float *verts, const float *normals,
                                size_t nv);

#endif
