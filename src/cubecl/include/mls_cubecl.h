#ifndef HERCULANEUM_MLS_CUBECL_H
#define HERCULANEUM_MLS_CUBECL_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

enum mls_cubecl_backend {
    MLS_BACKEND_AUTO = 0,
    MLS_BACKEND_RUST_CPU = 1,
    MLS_BACKEND_CUBECL_CPU = 2,
    MLS_BACKEND_CUBECL_HIP = 3,
};

/* Exact ScrollFiesta symbol: one projection per call. `verts` is both query
 * and support and contains nv z-y-x triples. The grid is rebuilt from it on
 * every call, so ScrollFiesta retains its existing 5/20-pass caller loops.
 * `arena` is scratch state and is intentionally ignored. `out_normals` may be
 * NULL; the input may alias either output, but the two output buffers must not
 * overlap. Select the backend with the MLS_BACKEND environment variable. */
void MLS_project_verts(
    void *arena,
    const float *verts,
    size_t nv,
    float radius_vox,
    const float *cell_origin,
    float *out_verts,
    float *out_normals);

/* Explicit static-support API for standalone xyz benchmarks/imported grids.
 * `arena` must come from an MLS_cubecl_arena_create_* constructor. */
void MLS_cubecl_project_with_arena(
    void *arena,
    const float *verts_xyz,
    size_t nv,
    float radius_vox,
    float *out_verts_xyz,
    float *out_normals_xyz);

/* Build a Rust-owned arena and counting-sort grid from xyz point triples. */
void *MLS_cubecl_arena_create_from_points(
    const float *points_xyz,
    size_t npoints,
    float cell_size,
    const float *cell_origin_or_null);

/* Import an already-built counting-sort cell grid. Data are copied, so the
 * caller may release its input buffers after this function returns. */
void *MLS_cubecl_arena_create_from_grid(
    const float *points_xyz,
    size_t npoints,
    const uint32_t *cell_offsets,
    size_t cell_offsets_len,
    const uint32_t *point_indices,
    size_t point_indices_len,
    const uint32_t dims_xyz[3],
    const float cell_origin[3],
    float cell_size);

void MLS_cubecl_arena_destroy(void *arena);
int32_t MLS_cubecl_arena_set_backend(void *arena, uint32_t backend);
uint32_t MLS_cubecl_available_backends(void);
const char *MLS_cubecl_last_error(void);
const char *MLS_cubecl_backend_name(uint32_t backend);
int32_t MLS_cubecl_backend_from_name(const char *name);

#ifdef __cplusplus
}
#endif

#endif
