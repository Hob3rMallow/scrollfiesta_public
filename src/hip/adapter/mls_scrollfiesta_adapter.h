/* Adapter boundary informed by ScrollFiesta MIT source, Copyright (c) 2026 Nicholas Vining.
 * No upstream source is copied into the validated clean HIP kernel.
 * See ABI_MAPPING.md for file:line observations and inferences.
 */
#ifndef MLS_SCROLLFIESTA_ADAPTER_H
#define MLS_SCROLLFIESTA_ADAPTER_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ScrollFiesta-compatible ABI. The upstream Arena_T argument is scratch-only
 * for the CPU implementation; this adapter constructs the clean HIP arena from
 * the supplied support/query cloud for each call. Buffers use upstream z,y,x
 * ordering. */
void MLS_project_verts(void* upstream_arena,
                       const float* verts_zyx, size_t nv,
                       float radius_vox,
                       const float cell_origin_zyx[3],
                       float* out_verts_zyx,
                       float* out_normals_zyx);

/* Complete caller-level LOP operator. Each pass rebuilds support from the
 * previous pass output, matching mesh_extract (5) and mesh_resplit (20). */
int mls_scrollfiesta_lop_project(const float* verts_zyx, size_t nv,
                                 float radius_vox,
                                 const float cell_origin_zyx[3],
                                 uint32_t passes,
                                 float* out_verts_zyx,
                                 float* out_normals_zyx);

#ifdef __cplusplus
}
#endif

#endif /* MLS_SCROLLFIESTA_ADAPTER_H */
