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

/* The same normals, bit for bit, for a mesh stored as nrun vertex-disjoint runs (run r: vertices
 * [vofs[r], vofs[r+1]), faces [fofs[r], fofs[r+1]), vofs[nrun] == nv), computed run-parallel. A face
 * that leaves its run falls back to MeshNormals_compute over fofs[nrun] faces. */
float *MeshNormals_compute_runs(const float *verts, size_t nv, const int32_t *faces,
                                const size_t *vofs, const size_t *fofs, size_t nrun);

/* One run of MeshNormals_compute_runs: normal[v0 .. v1) from faces [f0, f1), which name only those vertices (the
 * caller's guarantee), with the same arithmetic. */
void MeshNormals_run(const float *verts, const int32_t *faces, size_t f0, size_t f1, size_t v0, size_t v1,
                     float *normal);

/* MeshNormals_compute_runs and MeshNormals_run against MeshNormals_compute; 0 on a pass. */
int MeshNormals_selftest(void);

#endif
