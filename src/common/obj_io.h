#ifndef OBJ_IO_INCLUDED
#define OBJ_IO_INCLUDED

#include <stdint.h>
#include <stddef.h>
#include "arena.h"

/* Write mesh as text OBJ.
 * Vertices as "v z y x\n".
 * Faces as "f a+1 b+1 c+1\n" (OBJ is 1-indexed).
 * Returns 0 on success, -1 on failure. */
int ObjIO_write(const char *path, const float *verts, size_t nv,
                const int32_t *faces, size_t nf);

/* Write mesh as text OBJ with per-vertex RGB color.
 * Vertices as "v z y x r g b\n" (OBJ extended format).
 * color[3] is a single RGB triplet applied to all vertices.
 * MeshLab and Blender both support this format.
 * Returns 0 on success, -1 on failure. */
int ObjIO_write_colored(const char *path, const float *verts, size_t nv,
                        const int32_t *faces, size_t nf,
                        const float color[3]);

/* Write mesh with two vertex colors split at nv_split.
 * Vertices [0, nv_split) get color1, [nv_split, nv) get color2.
 * Returns 0 on success, -1 on failure. */
int ObjIO_write_twotone(const char *path, const float *verts, size_t nv,
                        const int32_t *faces, size_t nf,
                        size_t nv_split,
                        const float color1[3], const float color2[3]);

/* Write mesh with per-vertex RGB colors.
 * colors[i*3+{0,1,2}] = R,G,B for vertex i.
 * Returns 0 on success, -1 on failure. */
int ObjIO_write_per_vertex_color(const char *path,
                                 const float *verts, size_t nv,
                                 const int32_t *faces, size_t nf,
                                 const float *colors);

/* Read mesh from text OBJ. Arena-allocates output arrays.
 * Returns 0 on success, -1 on failure. */
int ObjIO_read(Arena_T arena, const char *path,
               float **out_verts, size_t *out_nv,
               int32_t **out_faces, size_t *out_nf);

/* Write mesh as text OBJ carrying per-vertex texture (UV) coordinates.
 * Vertices as "v z y x\n", UVs as "vt u v\n" (uv[i*2+0]=u, uv[i*2+1]=v),
 * faces as "f a/a b/b c/c\n" (1-indexed, shared vertex/uv index).
 * This is the interchange form consumed by ScrollPrize's vc_obj2tifxyz.
 * Returns 0 on success, -1 on failure. */
int ObjIO_write_uv(const char *path, const float *verts, size_t nv,
                   const int32_t *faces, size_t nf, const float *uv);

#endif
