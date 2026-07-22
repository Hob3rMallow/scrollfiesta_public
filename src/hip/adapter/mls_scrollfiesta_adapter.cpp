/* Adapter boundary informed by ScrollFiesta MIT source, Copyright (c) 2026 Nicholas Vining.
 * No upstream source is copied into the validated clean HIP kernel.
 * See ABI_MAPPING.md for file:line observations and inferences.
 */
#include "mls_scrollfiesta_adapter.h"
#include "../impl/mls_project_clean.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <new>
#include <vector>

extern "C" void mls_clean_project_verts(void* arena,
                                         const float* verts_xyz, size_t nv,
                                         float radius_vox,
                                         const float* cell_origin_xyz,
                                         float* out_verts_xyz,
                                         float* out_normals_xyz);

namespace {

bool finite3(const float* p) {
    return p && std::isfinite(p[0]) && std::isfinite(p[1]) && std::isfinite(p[2]);
}

bool triple_count_ok(size_t n) {
    return n <= std::numeric_limits<size_t>::max() / (3 * sizeof(float));
}

void zyx_to_xyz(const float* zyx, float* xyz, size_t n);
void xyz_to_zyx(const float* xyz, float* zyx, size_t n);

bool project_one_pass(const float* verts_zyx, size_t nv, float radius_vox,
                      const float cell_origin_zyx[3], float* out_verts_zyx,
                      float* out_normals_zyx) {
    if (nv == 0) return true;
    if (!verts_zyx || !out_verts_zyx || !triple_count_ok(nv) ||
        !std::isfinite(radius_vox) || !(radius_vox > 0.0f) ||
        !finite3(cell_origin_zyx)) {
        return false;
    }

    std::vector<float> support_xyz;
    std::vector<float> out_xyz;
    std::vector<float> normals_xyz;
    try {
        support_xyz.resize(nv * 3);
        out_xyz.resize(nv * 3);
        normals_xyz.resize(nv * 3);
    } catch (const std::bad_alloc&) {
        return false;
    }

    zyx_to_xyz(verts_zyx, support_xyz.data(), nv);
    const float clean_origin_xyz[3] = {
        -cell_origin_zyx[2],
        -cell_origin_zyx[1],
        -cell_origin_zyx[0],
    };

    mls_clean_arena* clean_arena = nullptr;
    const float cell_size_vox = radius_vox;
    mls_clean_status status = mls_clean_arena_create(support_xyz.data(), nv,
                                                     cell_size_vox,
                                                     clean_origin_xyz,
                                                     &clean_arena);
    if (status != MLS_CLEAN_STATUS_OK || !clean_arena) return false;

    mls_clean_project_verts(clean_arena, support_xyz.data(), nv, radius_vox,
                            clean_origin_xyz, out_xyz.data(), normals_xyz.data());
    status = mls_clean_last_status();
    mls_clean_arena_destroy(clean_arena);
    if (status != MLS_CLEAN_STATUS_OK &&
        status != MLS_CLEAN_STATUS_DEGENERATE_NEIGHBORHOOD) {
        return false;
    }

    xyz_to_zyx(out_xyz.data(), out_verts_zyx, nv);
    if (out_normals_zyx) xyz_to_zyx(normals_xyz.data(), out_normals_zyx, nv);
    return true;
}

void zyx_to_xyz(const float* zyx, float* xyz, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        xyz[3 * i + 0] = zyx[3 * i + 2];
        xyz[3 * i + 1] = zyx[3 * i + 1];
        xyz[3 * i + 2] = zyx[3 * i + 0];
    }
}

void xyz_to_zyx(const float* xyz, float* zyx, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        zyx[3 * i + 0] = xyz[3 * i + 2];
        zyx[3 * i + 1] = xyz[3 * i + 1];
        zyx[3 * i + 2] = xyz[3 * i + 0];
    }
}

} // namespace

extern "C" void MLS_project_verts(void* upstream_arena,
                                  const float* verts_zyx, size_t nv,
                                  float radius_vox,
                                  const float cell_origin_zyx[3],
                                  float* out_verts_zyx,
                                  float* out_normals_zyx) {
    (void)upstream_arena;
    (void)project_one_pass(verts_zyx, nv, radius_vox, cell_origin_zyx,
                           out_verts_zyx, out_normals_zyx);
}

extern "C" int mls_scrollfiesta_lop_project(
    const float* verts_zyx, size_t nv, float radius_vox,
    const float cell_origin_zyx[3], uint32_t passes,
    float* out_verts_zyx, float* out_normals_zyx) {
    if (passes == 0 || (nv != 0 && (!verts_zyx || !out_verts_zyx)) ||
        !triple_count_ok(nv)) {
        return -1;
    }
    if (nv == 0) return 0;

    std::vector<float> current;
    std::vector<float> next;
    std::vector<float> normals;
    try {
        current.assign(verts_zyx, verts_zyx + nv * 3);
        next.resize(nv * 3);
        normals.resize(nv * 3);
    } catch (const std::bad_alloc&) {
        return -1;
    }

    for (uint32_t pass = 0; pass < passes; ++pass) {
        if (!project_one_pass(current.data(), nv, radius_vox, cell_origin_zyx,
                              next.data(), normals.data())) {
            return -1;
        }
        current.swap(next);
    }

    std::copy(current.begin(), current.end(), out_verts_zyx);
    if (out_normals_zyx) {
        std::copy(normals.begin(), normals.end(), out_normals_zyx);
    }
    return 0;
}
