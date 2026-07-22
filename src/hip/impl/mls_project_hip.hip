#include "mls_project_clean.h"

#include <hip/hip_runtime.h>

#include <algorithm>
#include <atomic>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <limits>
#include <mutex>
#include <new>
#include <vector>

constexpr uint64_t kArenaMagic = UINT64_C(0x434c4e4d4c534152); /* CLNMLSAR */
constexpr int kBlockSize = 256;
constexpr int kJacobiIterations = 100;
constexpr int kMaxCellSpan = 64;

struct Point4 {
    float x, y, z, unused;
};

struct CellSlot {
    int32_t x, y, z;
    uint32_t start;
    uint32_t count; /* zero marks an empty hash slot */
};

struct HostEntry {
    int32_t x, y, z;
    uint32_t ordinal;
    Point4 point;
};

static_assert(sizeof(Point4) == 16, "Point4 must be one vector load");
static_assert(sizeof(CellSlot) == 20, "Unexpected CellSlot padding");

namespace {

thread_local mls_clean_status g_last_status = MLS_CLEAN_STATUS_OK;

__host__ __device__ inline uint32_t cell_hash(int32_t x, int32_t y, int32_t z) {
    uint32_t h = static_cast<uint32_t>(x) * UINT32_C(0x8da6b343);
    h ^= static_cast<uint32_t>(y) * UINT32_C(0xd8163841);
    h ^= static_cast<uint32_t>(z) * UINT32_C(0xcb1ab31f);
    h ^= h >> 16;
    h *= UINT32_C(0x7feb352d);
    h ^= h >> 15;
    h *= UINT32_C(0x846ca68b);
    return h ^ (h >> 16);
}

inline bool finite_float(float x) { return std::isfinite(x); }

inline bool host_cell_coord(float x, float origin, float inv_cell, int32_t* out) {
    const float c = std::floor((x - origin) * inv_cell);
    if (!std::isfinite(c) || c < -2147483648.0f || c >= 2147483648.0f) return false;
    *out = static_cast<int32_t>(c);
    return true;
}

inline bool mul_overflow_size(size_t a, size_t b, size_t* out) {
    if (a != 0 && b > std::numeric_limits<size_t>::max() / a) return true;
    *out = a * b;
    return false;
}

inline bool pointer_range_ok(const void* p, size_t bytes) {
    const uintptr_t begin = reinterpret_cast<uintptr_t>(p);
    return begin <= std::numeric_limits<uintptr_t>::max() - bytes;
}

inline bool ranges_overlap(const void* a, const void* b, size_t bytes) {
    const uintptr_t ab = reinterpret_cast<uintptr_t>(a);
    const uintptr_t bb = reinterpret_cast<uintptr_t>(b);
    return ab < bb + bytes && bb < ab + bytes;
}

class ScopedDevice {
public:
    explicit ScopedDevice(int wanted) : previous_(-1), ok_(false), changed_(false) {
        if (hipGetDevice(&previous_) != hipSuccess) return;
        if (previous_ != wanted) {
            if (hipSetDevice(wanted) != hipSuccess) return;
            changed_ = true;
        }
        ok_ = true;
    }
    ~ScopedDevice() {
        if (changed_) (void)hipSetDevice(previous_);
    }
    bool ok() const { return ok_; }
private:
    int previous_;
    bool ok_;
    bool changed_;
};

} /* namespace */

struct mls_clean_arena {
    uint64_t magic = 0;
    size_t point_count = 0;
    size_t query_capacity = 0;
    float cell_size_vox = 0.0f;
    float inv_cell_size = 0.0f;
    float index_origin[3] = {0.0f, 0.0f, 0.0f};
    int device_id = -1;
    uint32_t cell_table_size = 0;
    Point4* d_points = nullptr;
    CellSlot* d_cells = nullptr;
    float* d_queries = nullptr;
    float* d_out_verts = nullptr;
    float* d_out_normals = nullptr;
    int* d_kernel_status = nullptr;
    uint32_t* d_problem_count = nullptr;
    hipStream_t stream = nullptr;
    std::atomic<int> last_status{MLS_CLEAN_STATUS_OK};
    std::atomic<size_t> last_problem_count{0};
    mutable std::mutex mutex;
};

namespace {

inline bool valid_arena(const mls_clean_arena* arena) {
    return arena != nullptr && arena->magic == kArenaMagic;
}

void record_status(mls_clean_arena* arena, mls_clean_status status, size_t problems = 0) {
    g_last_status = status;
    if (valid_arena(arena)) {
        arena->last_status.store(static_cast<int>(status), std::memory_order_relaxed);
        arena->last_problem_count.store(problems, std::memory_order_relaxed);
    }
}

void release_device_storage(mls_clean_arena* arena) {
    if (!arena) return;
    if (arena->d_queries) (void)hipFree(arena->d_queries);
    if (arena->d_out_verts) (void)hipFree(arena->d_out_verts);
    if (arena->d_out_normals) (void)hipFree(arena->d_out_normals);
    if (arena->d_kernel_status) (void)hipFree(arena->d_kernel_status);
    if (arena->d_problem_count) (void)hipFree(arena->d_problem_count);
    if (arena->d_cells) (void)hipFree(arena->d_cells);
    if (arena->d_points) (void)hipFree(arena->d_points);
    if (arena->stream) (void)hipStreamDestroy(arena->stream);
    arena->d_queries = nullptr;
    arena->d_out_verts = nullptr;
    arena->d_out_normals = nullptr;
    arena->d_kernel_status = nullptr;
    arena->d_problem_count = nullptr;
    arena->d_cells = nullptr;
    arena->d_points = nullptr;
    arena->stream = nullptr;
}

mls_clean_status ensure_query_capacity(mls_clean_arena* arena, size_t needed) {
    if (needed <= arena->query_capacity) return MLS_CLEAN_STATUS_OK;
    size_t capacity = arena->query_capacity ? arena->query_capacity : 1024;
    while (capacity < needed) {
        if (capacity > std::numeric_limits<size_t>::max() / 2) {
            capacity = needed;
            break;
        }
        capacity *= 2;
    }
    size_t values = 0, bytes = 0;
    if (mul_overflow_size(capacity, size_t{3}, &values) ||
        mul_overflow_size(values, sizeof(float), &bytes)) {
        return MLS_CLEAN_STATUS_SIZE_OVERFLOW;
    }

    float* q = nullptr;
    float* ov = nullptr;
    float* on = nullptr;
    if (hipMalloc(reinterpret_cast<void**>(&q), bytes) != hipSuccess ||
        hipMalloc(reinterpret_cast<void**>(&ov), bytes) != hipSuccess ||
        hipMalloc(reinterpret_cast<void**>(&on), bytes) != hipSuccess) {
        if (q) (void)hipFree(q);
        if (ov) (void)hipFree(ov);
        if (on) (void)hipFree(on);
        return MLS_CLEAN_STATUS_OUT_OF_MEMORY;
    }
    if (arena->d_queries) (void)hipFree(arena->d_queries);
    if (arena->d_out_verts) (void)hipFree(arena->d_out_verts);
    if (arena->d_out_normals) (void)hipFree(arena->d_out_normals);
    arena->d_queries = q;
    arena->d_out_verts = ov;
    arena->d_out_normals = on;
    arena->query_capacity = capacity;
    return MLS_CLEAN_STATUS_OK;
}

__device__ inline bool device_finite(float x) {
    return __builtin_isfinite(x);
}

__device__ inline bool device_finite_double(double x) {
    return __builtin_isfinite(x);
}

__device__ inline bool device_finite3(float x, float y, float z) {
    return device_finite(x) && device_finite(y) && device_finite(z);
}

__device__ inline bool device_cell_coord(float x, float origin, float inv_cell,
                                         int32_t* out) {
    const float c = floorf((x - origin) * inv_cell);
    if (!device_finite(c) || c < -2147483648.0f || c >= 2147483648.0f) return false;
    *out = static_cast<int32_t>(c);
    return true;
}

__device__ inline const CellSlot* lookup_cell(const CellSlot* table,
                                              uint32_t table_size,
                                              int32_t x, int32_t y, int32_t z) {
    const uint32_t mask = table_size - 1;
    uint32_t slot = cell_hash(x, y, z) & mask;
    for (uint32_t probe = 0; probe < table_size; ++probe) {
        const CellSlot* cell = table + slot;
        if (cell->count == 0) return nullptr;
        if (cell->x == x && cell->y == y && cell->z == z) return cell;
        slot = (slot + 1) & mask;
    }
    return nullptr;
}

__device__ inline float wendland_weight(float dx, float dy, float dz,
                                        float radius2, float radius) {
    const float d2 = dx * dx + dy * dy + dz * dz;
    if (!(d2 < radius2)) return 0.0f;
    const float q = sqrtf(d2) / radius;
    if (!(q < 1.0f)) return 0.0f;
    const float t = 1.0f - q;
    const float t2 = t * t;
    return (t2 * t2) * (4.0f * q + 1.0f);
}

__device__ inline void compensated_add(float value, float* sum, float* correction) {
    const float adjusted = value - *correction;
    const float next = *sum + adjusted;
    *correction = (next - *sum) - adjusted;
    *sum = next;
}

__device__ bool centroid_pass(const Point4* points, const CellSlot* table,
                              uint32_t table_size, float inv_cell,
                              const float origin[3], int span,
                              float radius2, float radius,
                              float x, float y, float z,
                              float* weight, float* cx, float* cy, float* cz,
                              uint32_t* positive_samples) {
    int32_t bx, by, bz;
    if (!device_cell_coord(x, origin[0], inv_cell, &bx) ||
        !device_cell_coord(y, origin[1], inv_cell, &by) ||
        !device_cell_coord(z, origin[2], inv_cell, &bz)) return false;
    float sw = 0.0f, sx = 0.0f, sy = 0.0f, sz = 0.0f;
    float cw = 0.0f, cx_sum = 0.0f, cy_sum = 0.0f, cz_sum = 0.0f;
    uint32_t used = 0;
    for (int oz = -span; oz <= span; ++oz) {
        const int64_t gz64 = static_cast<int64_t>(bz) + oz;
        if (gz64 < INT32_MIN || gz64 > INT32_MAX) continue;
        for (int oy = -span; oy <= span; ++oy) {
            const int64_t gy64 = static_cast<int64_t>(by) + oy;
            if (gy64 < INT32_MIN || gy64 > INT32_MAX) continue;
            for (int ox = -span; ox <= span; ++ox) {
                const int64_t gx64 = static_cast<int64_t>(bx) + ox;
                if (gx64 < INT32_MIN || gx64 > INT32_MAX) continue;
                const CellSlot* cell = lookup_cell(table, table_size,
                    static_cast<int32_t>(gx64), static_cast<int32_t>(gy64),
                    static_cast<int32_t>(gz64));
                if (!cell) continue;
                const uint32_t end = cell->start + cell->count;
                for (uint32_t j = cell->start; j < end; ++j) {
                    const Point4 p = points[j];
                    const float w = wendland_weight(p.x - x, p.y - y, p.z - z,
                                                    radius2, radius);
                    if (w > 0.0f) ++used;
                    compensated_add(w, &sw, &cw);
                    compensated_add(w * (p.x - x), &sx, &cx_sum);
                    compensated_add(w * (p.y - y), &sy, &cy_sum);
                    compensated_add(w * (p.z - z), &sz, &cz_sum);
                }
            }
        }
    }
    if (!(sw > 0.0f) || !device_finite3(sx, sy, sz) || !device_finite(sw)) return false;
    const float iw = 1.0f / sw;
    *weight = sw;
    *cx = x + sx * iw;
    *cy = y + sy * iw;
    *cz = z + sz * iw;
    *positive_samples = used;
    return device_finite3(*cx, *cy, *cz);
}

__device__ bool covariance_pass(const Point4* points, const CellSlot* table,
                                uint32_t table_size, float inv_cell,
                                const float origin[3], int span,
                                float radius2, float radius,
                                float x, float y, float z, float weight,
                                float cx, float cy, float cz, float cov[6]) {
    int32_t bx, by, bz;
    if (!device_cell_coord(x, origin[0], inv_cell, &bx) ||
        !device_cell_coord(y, origin[1], inv_cell, &by) ||
        !device_cell_coord(z, origin[2], inv_cell, &bz)) return false;
    float xx = 0.0f, xy = 0.0f, xz = 0.0f;
    float yy = 0.0f, yz = 0.0f, zz = 0.0f;
    float cxx = 0.0f, cxy = 0.0f, cxz = 0.0f;
    float cyy = 0.0f, cyz = 0.0f, czz = 0.0f;
    for (int oz = -span; oz <= span; ++oz) {
        const int64_t gz64 = static_cast<int64_t>(bz) + oz;
        if (gz64 < INT32_MIN || gz64 > INT32_MAX) continue;
        for (int oy = -span; oy <= span; ++oy) {
            const int64_t gy64 = static_cast<int64_t>(by) + oy;
            if (gy64 < INT32_MIN || gy64 > INT32_MAX) continue;
            for (int ox = -span; ox <= span; ++ox) {
                const int64_t gx64 = static_cast<int64_t>(bx) + ox;
                if (gx64 < INT32_MIN || gx64 > INT32_MAX) continue;
                const CellSlot* cell = lookup_cell(table, table_size,
                    static_cast<int32_t>(gx64), static_cast<int32_t>(gy64),
                    static_cast<int32_t>(gz64));
                if (!cell) continue;
                const uint32_t end = cell->start + cell->count;
                for (uint32_t j = cell->start; j < end; ++j) {
                    const Point4 p = points[j];
                    const float w = wendland_weight(p.x - x, p.y - y, p.z - z,
                                                    radius2, radius);
                    const float dx = p.x - cx, dy = p.y - cy, dz = p.z - cz;
                    compensated_add(w * dx * dx, &xx, &cxx);
                    compensated_add(w * dx * dy, &xy, &cxy);
                    compensated_add(w * dx * dz, &xz, &cxz);
                    compensated_add(w * dy * dy, &yy, &cyy);
                    compensated_add(w * dy * dz, &yz, &cyz);
                    compensated_add(w * dz * dz, &zz, &czz);
                }
            }
        }
    }
    const float iw = 1.0f / weight;
    cov[0] = xx * iw; cov[1] = xy * iw; cov[2] = xz * iw;
    cov[3] = yy * iw; cov[4] = yz * iw; cov[5] = zz * iw;
    return device_finite3(cov[0], cov[1], cov[2]) &&
           device_finite3(cov[3], cov[4], cov[5]);
}

__device__ bool smallest_eigenvector(const float cov[6], float* nx, float* ny,
                                     float* nz) {
    /* Solve in ScrollFiesta's native z,y,x basis. Jacobi rotation order is
     * basis-sensitive near repeated eigenvalues, including raw vector sign. */
    double a[3][3] = {{cov[5], cov[4], cov[2]},
                      {cov[4], cov[3], cov[1]},
                      {cov[2], cov[1], cov[0]}};
    double v[3][3] = {{1.0, 0.0, 0.0},
                      {0.0, 1.0, 0.0},
                      {0.0, 0.0, 1.0}};
    for (int iteration = 0; iteration < kJacobiIterations; ++iteration) {
        double max_off = 0.0;
        int p = 0, q = 1;
        for (int row = 0; row < 3; ++row) {
            for (int column = row + 1; column < 3; ++column) {
                const double candidate = fabs(a[row][column]);
                if (candidate > max_off) {
                    max_off = candidate;
                    p = row;
                    q = column;
                }
            }
        }
        if (max_off < 1.0e-15) break;
        const double app = a[p][p], aqq = a[q][q], apq = a[p][q];
        const double difference = aqq - app;
        const double angle = fabs(difference) < 1.0e-30
            ? 0.78539816339744830961
            : 0.5 * atan2(2.0 * apq, difference);
        const double c = cos(angle), s = sin(angle);
        double next[3][3];
        for (int row = 0; row < 3; ++row)
            for (int column = 0; column < 3; ++column)
                next[row][column] = a[row][column];
        next[p][p] = c * c * app - 2.0f * s * c * apq + s * s * aqq;
        next[q][q] = s * s * app + 2.0f * s * c * apq + c * c * aqq;
        next[p][q] = next[q][p] = 0.0f;
        for (int row = 0; row < 3; ++row) {
            if (row == p || row == q) continue;
            next[row][p] = next[p][row] = c * a[row][p] - s * a[row][q];
            next[row][q] = next[q][row] = s * a[row][p] + c * a[row][q];
        }
        for (int row = 0; row < 3; ++row)
            for (int column = 0; column < 3; ++column)
                a[row][column] = next[row][column];
        for (int row = 0; row < 3; ++row) {
            const double vp = v[row][p], vq = v[row][q];
            v[row][p] = c * vp - s * vq;
            v[row][q] = s * vp + c * vq;
        }
    }
    if (!device_finite_double(a[0][0]) || !device_finite_double(a[1][1]) ||
        !device_finite_double(a[2][2])) return false;
    int least = 0;
    if (a[1][1] < a[least][least]) least = 1;
    if (a[2][2] < a[least][least]) least = 2;

    double z = v[0][least], y = v[1][least], x = v[2][least];
    const double norm2 = x * x + y * y + z * z;
    if (!(norm2 > 0.0) || !device_finite_double(norm2)) return false;
    const double inv_norm = 1.0 / sqrt(norm2);
    x *= inv_norm; y *= inv_norm; z *= inv_norm;
    const double sign_dot = (x + y) + z;
    if (sign_dot < -1.0e-9) {
        x = -x; y = -y; z = -z;
    } else if (sign_dot < 1.0e-9) {
        int largest = 0;
        if (fabs(y) > fabs(x)) largest = 1;
        if (fabs(z) > fabs(largest == 0 ? x : y)) largest = 2;
        const double sign = largest == 0 ? x : (largest == 1 ? y : z);
        if (sign < 0.0) { x = -x; y = -y; z = -z; }
    }
    if (!device_finite_double(x) || !device_finite_double(y) ||
        !device_finite_double(z)) return false;
    *nx = static_cast<float>(x);
    *ny = static_cast<float>(y);
    *nz = static_cast<float>(z);
    return true;
}

__device__ inline void fail_query(size_t i, float ox, float oy, float oz,
                                  float* out_verts, float* out_normals,
                                  int code, int* status, uint32_t* problems) {
    out_verts[3 * i] = ox;
    out_verts[3 * i + 1] = oy;
    out_verts[3 * i + 2] = oz;
    out_normals[3 * i] = 0.0f;
    out_normals[3 * i + 1] = 0.0f;
    out_normals[3 * i + 2] = 0.0f;
    atomicMax(status, code);
    atomicAdd(problems, 1u);
}

__global__ void mls_project_kernel(const Point4* points, const CellSlot* cells,
                                   uint32_t cell_table_size, float inv_cell,
                                   float ox, float oy, float oz,
                                   const float* queries, size_t nv,
                                   float radius, int cell_span,
                                   uint32_t iterations,
                                   float* out_verts, float* out_normals,
                                   int* status, uint32_t* problems) {
    const size_t i = static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= nv) return;
    const float original_x = queries[3 * i];
    const float original_y = queries[3 * i + 1];
    const float original_z = queries[3 * i + 2];
    if (!device_finite3(original_x, original_y, original_z)) {
        fail_query(i, 0.0f, 0.0f, 0.0f, out_verts, out_normals,
                   MLS_CLEAN_STATUS_NONFINITE_INPUT, status, problems);
        return;
    }
    const float origin[3] = {ox, oy, oz};
    const float radius2 = radius * radius;
    float x = original_x, y = original_y, z = original_z;
    float nx = 0.0f, ny = 0.0f, nz = 0.0f;
    for (uint32_t iteration = 0; iteration < iterations; ++iteration) {
        float weight, cx, cy, cz;
        uint32_t positive_samples = 0;
        if (!centroid_pass(points, cells, cell_table_size, inv_cell, origin,
                           cell_span, radius2, radius, x, y, z,
                           &weight, &cx, &cy, &cz, &positive_samples) ||
            positive_samples < 4) {
            fail_query(i, original_x, original_y, original_z, out_verts, out_normals,
                       MLS_CLEAN_STATUS_DEGENERATE_NEIGHBORHOOD, status, problems);
            return;
        }
        float cov[6];
        if (!covariance_pass(points, cells, cell_table_size, inv_cell, origin,
                             cell_span, radius2, radius, x, y, z, weight,
                             cx, cy, cz, cov) ||
            !smallest_eigenvector(cov, &nx, &ny, &nz)) {
            fail_query(i, original_x, original_y, original_z, out_verts, out_normals,
                       MLS_CLEAN_STATUS_DEGENERATE_NEIGHBORHOOD, status, problems);
            return;
        }
        const float signed_distance = nx * (x - cx) + ny * (y - cy) + nz * (z - cz);
        x -= signed_distance * nx;
        y -= signed_distance * ny;
        z -= signed_distance * nz;
        if (!device_finite3(x, y, z)) {
            fail_query(i, original_x, original_y, original_z, out_verts, out_normals,
                       MLS_CLEAN_STATUS_NUMERICAL_FAILURE, status, problems);
            return;
        }
    }
    out_verts[3 * i] = x;
    out_verts[3 * i + 1] = y;
    out_verts[3 * i + 2] = z;
    out_normals[3 * i] = nx;
    out_normals[3 * i + 1] = ny;
    out_normals[3 * i + 2] = nz;
}

} /* namespace */

extern "C" mls_clean_status mls_clean_arena_create(
    const float* support_points_xyz, size_t point_count, float cell_size_vox,
    const float index_origin[3], mls_clean_arena** out_arena) {
    g_last_status = MLS_CLEAN_STATUS_INVALID_ARGUMENT;
    if (!out_arena) return g_last_status;
    *out_arena = nullptr;
    if (!support_points_xyz || point_count == 0 || !index_origin ||
        !finite_float(cell_size_vox) || !(cell_size_vox > 0.0f) ||
        !finite_float(index_origin[0]) || !finite_float(index_origin[1]) ||
        !finite_float(index_origin[2])) return g_last_status;
    if (point_count > UINT32_MAX) {
        g_last_status = MLS_CLEAN_STATUS_SIZE_OVERFLOW;
        return g_last_status;
    }
    size_t value_count = 0;
    if (mul_overflow_size(point_count, size_t{3}, &value_count)) {
        g_last_status = MLS_CLEAN_STATUS_SIZE_OVERFLOW;
        return g_last_status;
    }

    std::vector<HostEntry> entries;
    std::vector<Point4> sorted_points;
    std::vector<CellSlot> unique_cells;
    std::vector<CellSlot> table;
    try {
        entries.resize(point_count);
        sorted_points.resize(point_count);
    } catch (const std::bad_alloc&) {
        g_last_status = MLS_CLEAN_STATUS_OUT_OF_MEMORY;
        return g_last_status;
    }
    const float inv_cell = 1.0f / cell_size_vox;
    if (!finite_float(inv_cell)) {
        g_last_status = MLS_CLEAN_STATUS_INVALID_ARGUMENT;
        return g_last_status;
    }
    for (size_t i = 0; i < point_count; ++i) {
        const float x = support_points_xyz[3 * i];
        const float y = support_points_xyz[3 * i + 1];
        const float z = support_points_xyz[3 * i + 2];
        if (!finite_float(x) || !finite_float(y) || !finite_float(z)) {
            g_last_status = MLS_CLEAN_STATUS_NONFINITE_INPUT;
            return g_last_status;
        }
        HostEntry& e = entries[i];
        if (!host_cell_coord(x, index_origin[0], inv_cell, &e.x) ||
            !host_cell_coord(y, index_origin[1], inv_cell, &e.y) ||
            !host_cell_coord(z, index_origin[2], inv_cell, &e.z)) {
            g_last_status = MLS_CLEAN_STATUS_SIZE_OVERFLOW;
            return g_last_status;
        }
        e.ordinal = static_cast<uint32_t>(i);
        e.point = Point4{x, y, z, 0.0f};
    }
    std::sort(entries.begin(), entries.end(), [](const HostEntry& a, const HostEntry& b) {
        if (a.z != b.z) return a.z < b.z;
        if (a.y != b.y) return a.y < b.y;
        if (a.x != b.x) return a.x < b.x;
        /* Match ScrollFiesta's stable within-cell accumulation order. The
         * clean representation is x,y,z, so compare point z,y,x here. */
        if (a.point.z != b.point.z) return a.point.z < b.point.z;
        if (a.point.y != b.point.y) return a.point.y < b.point.y;
        if (a.point.x != b.point.x) return a.point.x < b.point.x;
        return a.ordinal < b.ordinal;
    });
    try {
        unique_cells.reserve(point_count);
    } catch (const std::bad_alloc&) {
        g_last_status = MLS_CLEAN_STATUS_OUT_OF_MEMORY;
        return g_last_status;
    }
    for (size_t i = 0; i < point_count;) {
        const size_t begin = i;
        const int32_t x = entries[i].x, y = entries[i].y, z = entries[i].z;
        do {
            sorted_points[i] = entries[i].point;
            ++i;
        } while (i < point_count && entries[i].x == x && entries[i].y == y &&
                 entries[i].z == z);
        unique_cells.push_back(CellSlot{x, y, z, static_cast<uint32_t>(begin),
                                        static_cast<uint32_t>(i - begin)});
    }
    if (unique_cells.size() > (size_t{1} << 30)) {
        g_last_status = MLS_CLEAN_STATUS_SIZE_OVERFLOW;
        return g_last_status;
    }
    size_t table_size = 1;
    const size_t target = unique_cells.size() * 2;
    while (table_size < target) table_size <<= 1;
    try {
        table.resize(table_size);
    } catch (const std::bad_alloc&) {
        g_last_status = MLS_CLEAN_STATUS_OUT_OF_MEMORY;
        return g_last_status;
    }
    const uint32_t mask = static_cast<uint32_t>(table_size - 1);
    for (const CellSlot& cell : unique_cells) {
        uint32_t slot = cell_hash(cell.x, cell.y, cell.z) & mask;
        while (table[slot].count != 0) slot = (slot + 1) & mask;
        table[slot] = cell;
    }

    mls_clean_arena* arena = new (std::nothrow) mls_clean_arena;
    if (!arena) {
        g_last_status = MLS_CLEAN_STATUS_OUT_OF_MEMORY;
        return g_last_status;
    }
    if (hipGetDevice(&arena->device_id) != hipSuccess) {
        delete arena;
        g_last_status = MLS_CLEAN_STATUS_HIP_ERROR;
        return g_last_status;
    }
    arena->point_count = point_count;
    arena->cell_size_vox = cell_size_vox;
    arena->inv_cell_size = inv_cell;
    arena->index_origin[0] = index_origin[0];
    arena->index_origin[1] = index_origin[1];
    arena->index_origin[2] = index_origin[2];
    arena->cell_table_size = static_cast<uint32_t>(table_size);
    size_t point_bytes = 0, table_bytes = 0;
    (void)mul_overflow_size(point_count, sizeof(Point4), &point_bytes);
    (void)mul_overflow_size(table_size, sizeof(CellSlot), &table_bytes);
    hipError_t err = hipStreamCreateWithFlags(&arena->stream, hipStreamNonBlocking);
    if (err == hipSuccess)
        err = hipMalloc(reinterpret_cast<void**>(&arena->d_points), point_bytes);
    if (err == hipSuccess)
        err = hipMalloc(reinterpret_cast<void**>(&arena->d_cells), table_bytes);
    if (err == hipSuccess)
        err = hipMalloc(reinterpret_cast<void**>(&arena->d_kernel_status), sizeof(int));
    if (err == hipSuccess)
        err = hipMalloc(reinterpret_cast<void**>(&arena->d_problem_count), sizeof(uint32_t));
    if (err == hipSuccess)
        err = hipMemcpy(arena->d_points, sorted_points.data(), point_bytes,
                        hipMemcpyHostToDevice);
    if (err == hipSuccess)
        err = hipMemcpy(arena->d_cells, table.data(), table_bytes,
                        hipMemcpyHostToDevice);
    if (err != hipSuccess) {
        release_device_storage(arena);
        delete arena;
        g_last_status = (err == hipErrorOutOfMemory) ? MLS_CLEAN_STATUS_OUT_OF_MEMORY
                                                     : MLS_CLEAN_STATUS_HIP_ERROR;
        return g_last_status;
    }
    arena->magic = kArenaMagic;
    *out_arena = arena;
    g_last_status = MLS_CLEAN_STATUS_OK;
    return g_last_status;
}

extern "C" void mls_clean_arena_destroy(mls_clean_arena* arena) {
    if (!valid_arena(arena)) {
        g_last_status = MLS_CLEAN_STATUS_WRONG_ARENA;
        return;
    }
    bool device_ok = false;
    {
        std::lock_guard<std::mutex> lock(arena->mutex);
        arena->magic = 0;
        ScopedDevice device(arena->device_id);
        device_ok = device.ok();
        if (device_ok) release_device_storage(arena);
    }
    g_last_status = device_ok ? MLS_CLEAN_STATUS_OK : MLS_CLEAN_STATUS_HIP_ERROR;
    delete arena;
}

extern "C" void mls_clean_project_verts_iterations(
                                  void* opaque_arena,
                                  const float* verts, size_t nv,
                                  float radius_vox,
                                  const float* cell_origin,
                                  uint32_t iterations,
                                  float* out_verts,
                                  float* out_normals) {
    mls_clean_arena* arena = static_cast<mls_clean_arena*>(opaque_arena);
    if (!valid_arena(arena)) {
        g_last_status = MLS_CLEAN_STATUS_WRONG_ARENA;
        return;
    }
    std::lock_guard<std::mutex> lock(arena->mutex);
    if (!cell_origin || iterations == 0 || !finite_float(radius_vox) ||
        !(radius_vox > 0.0f) ||
        !finite_float(cell_origin[0]) || !finite_float(cell_origin[1]) ||
        !finite_float(cell_origin[2])) {
        record_status(arena, MLS_CLEAN_STATUS_INVALID_ARGUMENT);
        return;
    }
    if (cell_origin[0] != arena->index_origin[0] ||
        cell_origin[1] != arena->index_origin[1] ||
        cell_origin[2] != arena->index_origin[2]) {
        record_status(arena, MLS_CLEAN_STATUS_ORIGIN_MISMATCH);
        return;
    }
    const float span_f = std::ceil(radius_vox * arena->inv_cell_size);
    if (!finite_float(span_f) || span_f < 1.0f || span_f > kMaxCellSpan) {
        record_status(arena, MLS_CLEAN_STATUS_INVALID_ARGUMENT);
        return;
    }
    if (nv == 0) {
        record_status(arena, MLS_CLEAN_STATUS_OK);
        return;
    }
    if (!verts || !out_verts || !out_normals ||
        (reinterpret_cast<uintptr_t>(verts) % alignof(float)) != 0 ||
        (reinterpret_cast<uintptr_t>(out_verts) % alignof(float)) != 0 ||
        (reinterpret_cast<uintptr_t>(out_normals) % alignof(float)) != 0) {
        record_status(arena, MLS_CLEAN_STATUS_INVALID_ARGUMENT);
        return;
    }
    size_t values = 0, bytes = 0;
    if (mul_overflow_size(nv, size_t{3}, &values) ||
        mul_overflow_size(values, sizeof(float), &bytes) ||
        !pointer_range_ok(verts, bytes) || !pointer_range_ok(out_verts, bytes) ||
        !pointer_range_ok(out_normals, bytes)) {
        record_status(arena, MLS_CLEAN_STATUS_SIZE_OVERFLOW);
        return;
    }
    if (ranges_overlap(verts, out_verts, bytes) ||
        ranges_overlap(verts, out_normals, bytes) ||
        ranges_overlap(out_verts, out_normals, bytes)) {
        record_status(arena, MLS_CLEAN_STATUS_INVALID_ARGUMENT);
        return;
    }
    if (nv > static_cast<size_t>(UINT32_MAX) * kBlockSize) {
        record_status(arena, MLS_CLEAN_STATUS_SIZE_OVERFLOW);
        return;
    }
    ScopedDevice device(arena->device_id);
    if (!device.ok()) {
        record_status(arena, MLS_CLEAN_STATUS_HIP_ERROR);
        return;
    }
    mls_clean_status capacity_status = ensure_query_capacity(arena, nv);
    if (capacity_status != MLS_CLEAN_STATUS_OK) {
        record_status(arena, capacity_status);
        return;
    }
    hipError_t err = hipMemcpyAsync(arena->d_queries, verts, bytes,
                                    hipMemcpyHostToDevice, arena->stream);
    if (err == hipSuccess)
        err = hipMemsetAsync(arena->d_kernel_status, 0, sizeof(int), arena->stream);
    if (err == hipSuccess)
        err = hipMemsetAsync(arena->d_problem_count, 0, sizeof(uint32_t), arena->stream);
    if (err != hipSuccess) {
        record_status(arena, MLS_CLEAN_STATUS_HIP_ERROR);
        return;
    }
    const dim3 block(kBlockSize);
    const dim3 grid(static_cast<uint32_t>((nv + kBlockSize - 1) / kBlockSize));
    hipLaunchKernelGGL(mls_project_kernel, grid, block, 0, arena->stream,
                      arena->d_points, arena->d_cells, arena->cell_table_size,
                      arena->inv_cell_size, arena->index_origin[0],
                      arena->index_origin[1], arena->index_origin[2],
                      arena->d_queries, nv, radius_vox, static_cast<int>(span_f),
                      iterations,
                      arena->d_out_verts, arena->d_out_normals,
                      arena->d_kernel_status, arena->d_problem_count);
    err = hipGetLastError();
    int kernel_status = MLS_CLEAN_STATUS_OK;
    uint32_t problem_count = 0;
    if (err == hipSuccess)
        err = hipMemcpyAsync(out_verts, arena->d_out_verts, bytes,
                             hipMemcpyDeviceToHost, arena->stream);
    if (err == hipSuccess)
        err = hipMemcpyAsync(out_normals, arena->d_out_normals, bytes,
                             hipMemcpyDeviceToHost, arena->stream);
    if (err == hipSuccess)
        err = hipMemcpyAsync(&kernel_status, arena->d_kernel_status, sizeof(int),
                             hipMemcpyDeviceToHost, arena->stream);
    if (err == hipSuccess)
        err = hipMemcpyAsync(&problem_count, arena->d_problem_count, sizeof(uint32_t),
                             hipMemcpyDeviceToHost, arena->stream);
    if (err == hipSuccess) err = hipStreamSynchronize(arena->stream);
    if (err != hipSuccess) {
        record_status(arena, MLS_CLEAN_STATUS_HIP_ERROR);
        return;
    }
    record_status(arena, static_cast<mls_clean_status>(kernel_status), problem_count);
}

extern "C" void MLS_project_verts(void* arena,
                                  const float* verts, size_t nv,
                                  float radius_vox,
                                  const float* cell_origin,
                                  float* out_verts,
                                  float* out_normals) {
    mls_clean_project_verts_iterations(arena, verts, nv, radius_vox, cell_origin,
                                       1u, out_verts, out_normals);
}

extern "C" mls_clean_status mls_clean_last_status(void) {
    return g_last_status;
}

extern "C" mls_clean_status mls_clean_arena_get_info(
    const mls_clean_arena* arena, mls_clean_arena_info* out_info) {
    if (!valid_arena(arena) || !out_info) {
        g_last_status = valid_arena(arena) ? MLS_CLEAN_STATUS_INVALID_ARGUMENT
                                          : MLS_CLEAN_STATUS_WRONG_ARENA;
        return g_last_status;
    }
    std::lock_guard<std::mutex> lock(arena->mutex);
    out_info->struct_size = sizeof(*out_info);
    out_info->point_count = arena->point_count;
    out_info->query_capacity = arena->query_capacity;
    out_info->last_problem_query_count =
        arena->last_problem_count.load(std::memory_order_relaxed);
    out_info->cell_size_vox = arena->cell_size_vox;
    out_info->index_origin[0] = arena->index_origin[0];
    out_info->index_origin[1] = arena->index_origin[1];
    out_info->index_origin[2] = arena->index_origin[2];
    out_info->device_id = arena->device_id;
    out_info->last_status = static_cast<mls_clean_status>(
        arena->last_status.load(std::memory_order_relaxed));
    g_last_status = MLS_CLEAN_STATUS_OK;
    return g_last_status;
}

extern "C" const char* mls_clean_status_string(mls_clean_status status) {
    switch (status) {
        case MLS_CLEAN_STATUS_OK: return "ok";
        case MLS_CLEAN_STATUS_INVALID_ARGUMENT: return "invalid argument";
        case MLS_CLEAN_STATUS_NONFINITE_INPUT: return "non-finite input";
        case MLS_CLEAN_STATUS_ORIGIN_MISMATCH: return "cell origin does not match arena";
        case MLS_CLEAN_STATUS_SIZE_OVERFLOW: return "size or grid-coordinate overflow";
        case MLS_CLEAN_STATUS_OUT_OF_MEMORY: return "out of memory";
        case MLS_CLEAN_STATUS_HIP_ERROR: return "HIP runtime error";
        case MLS_CLEAN_STATUS_DEGENERATE_NEIGHBORHOOD: return "degenerate neighborhood";
        case MLS_CLEAN_STATUS_NUMERICAL_FAILURE: return "numerical failure";
        case MLS_CLEAN_STATUS_WRONG_ARENA: return "not a clean-side arena";
        default: return "unknown clean-side status";
    }
}
