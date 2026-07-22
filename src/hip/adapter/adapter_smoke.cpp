#include "mls_scrollfiesta_adapter.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

int main() {
    const int side = 24;
    std::vector<float> verts;
    verts.reserve(side * side * 3);
    for (int y = 0; y < side; ++y) {
        for (int x = 0; x < side; ++x) {
            const float z = 0.03f * std::sin(0.2f * float(x)) + 0.02f * std::cos(0.3f * float(y));
            verts.push_back(z);
            verts.push_back(float(y));
            verts.push_back(float(x));
        }
    }
    std::vector<float> out(verts.size(), -999.0f);
    std::vector<float> normals(verts.size(), -999.0f);
    const float origin[3] = {0.0f, 0.0f, 0.0f};
    MLS_project_verts(nullptr, verts.data(), verts.size() / 3, 4.0f, origin,
                      out.data(), normals.data());

    std::vector<float> lop_out(verts.size());
    std::vector<float> lop_normals(verts.size());
    if (mls_scrollfiesta_lop_project(verts.data(), verts.size() / 3, 4.0f,
                                     origin, 5, lop_out.data(),
                                     lop_normals.data()) != 0) {
        return 3;
    }
    std::vector<float> manual_current = verts;
    std::vector<float> manual_next(verts.size());
    std::vector<float> manual_normals(verts.size());
    for (int pass = 0; pass < 5; ++pass) {
        MLS_project_verts(nullptr, manual_current.data(), verts.size() / 3,
                          4.0f, origin, manual_next.data(), manual_normals.data());
        manual_current.swap(manual_next);
    }
    float lop_max_delta = 0.0f;
    for (size_t i = 0; i < lop_out.size(); ++i) {
        lop_max_delta = std::max(lop_max_delta,
                                 std::fabs(lop_out[i] - manual_current[i]));
    }

    size_t finite = 0;
    double normal_norm_sum = 0.0;
    for (size_t i = 0; i < out.size() / 3; ++i) {
        const float z = out[3 * i + 0], y = out[3 * i + 1], x = out[3 * i + 2];
        const float nz = normals[3 * i + 0], ny = normals[3 * i + 1], nx = normals[3 * i + 2];
        if (std::isfinite(z) && std::isfinite(y) && std::isfinite(x) &&
            std::isfinite(nz) && std::isfinite(ny) && std::isfinite(nx)) {
            ++finite;
            normal_norm_sum += std::sqrt(double(nz) * nz + double(ny) * ny + double(nx) * nx);
        }
    }
    std::printf("finite=%zu total=%zu normal_norm_mean=%.6f lop5_max_delta=%.9g first_out_zyx=%.6f,%.6f,%.6f\n",
                finite, out.size() / 3, normal_norm_sum / double(out.size() / 3),
                lop_max_delta, out[0], out[1], out[2]);
    return finite == out.size() / 3 && lop_max_delta == 0.0f ? 0 : 2;
}
