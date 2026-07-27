#ifndef PLY_IO_INCLUDED
#define PLY_IO_INCLUDED

#include <stdint.h>
#include <stddef.h>
#include "arena.h"

/* Write oriented point cloud as binary little-endian PLY.
 * verts[n*3] = positions, normals[n*3] = unit normals.
 * Returns 0 on success, -1 on failure. */
int PlyIO_write_points(const char *path,
                       const float *verts, const float *normals, size_t n);

/* Read triangle mesh from binary PLY (PoissonRecon output).
 * Arena-allocates verts, faces, and optionally density.
 * Returns 0 on success, -1 on failure. */
int PlyIO_read_mesh(Arena_T arena, const char *path,
                    float **out_verts, size_t *out_nv,
                    int32_t **out_faces, size_t *out_nf,
                    float **out_density);

/* Defect class enum carried in the winding-mesh PLY's uchar 'defect'
 * property and as the colored-OBJ palette key. */
typedef enum {
    PLY_DEFECT_OK                = 0,
    PLY_DEFECT_NEAR_AXIS         = 1,  /* projection magnitude < R_min */
    PLY_DEFECT_FRAGMENT          = 2,  /* tiny disconnected component */
    PLY_DEFECT_ISOLATED_FROM_AXIS = 3, /* unreachable via valid BFS edges */
    PLY_DEFECT_DISAGREE          = 4   /* (w_max - w_min) > 0.5 */
} PlyDefectClass;

/* Write a triangle mesh + per-vertex {winding, confidence, defect}
 * to a binary little-endian PLY (the winding diagnostic's output
 * format).
 *
 * verts: [nv * 3] in project (z, y, x) convention. The writer emits
 *        them in (x, y, z) on disk so MeshLab loads without an axis
 *        flip.
 * faces: [nf * 3], int32.
 * winding, confidence: [nv] each (float).
 * defect: [nv] (uchar), values from PlyDefectClass.
 *
 * Returns 0 on success, -1 on I/O failure.
 */
int PlyIO_write_winding_mesh(const char *path,
                              const float *verts, size_t nv,
                              const int32_t *faces, size_t nf,
                              const float *winding,
                              const float *confidence,
                              const uint8_t *defect);

#endif
