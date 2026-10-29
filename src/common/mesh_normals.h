#ifndef MESH_NORMALS_INCLUDED
#define MESH_NORMALS_INCLUDED

#include <stddef.h>
#include <stdint.h>

/* Area-weighted, normalized per-vertex winding normals.  This is a purely
 * geometric operation: it has no RAW-volume, texture, or TIFF dependency.
 * Undefined vertices receive the zero vector.  The returned [nv*3] buffer is
 * malloc-owned by the caller (including the nv == 0 case). */
float *MeshNormals_compute(const float *verts, size_t nv,
                           const int32_t *faces, size_t nf);

#endif
