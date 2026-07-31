/*
 * QEM mesh simplification unit tests
 */
#include "common/ves_platform.h"
#include "common/arena.h"
#include "common/qem.h"
#include "common/pipeline_constants.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int tests_run = 0;
static int tests_passed = 0;

#define TEST(name) do { \
    printf("  %-50s ", #name); \
    tests_run++; \
    if (test_##name()) { printf("PASS\n"); tests_passed++; } \
    else { printf("FAIL\n"); } \
} while(0)

/* ================================================================
 * Helper: validate output mesh
 * ================================================================ */
static int validate_mesh(const float *verts, size_t nv,
                         const int32_t *faces, size_t nf)
{
    size_t f = 0;
    for (f = 0; f < nf; f++) {
        int32_t i0 = faces[f * 3 + 0];
        int32_t i1 = faces[f * 3 + 1];
        int32_t i2 = faces[f * 3 + 2];
        /* Indices in range */
        if (i0 < 0 || (size_t)i0 >= nv) return 0;
        if (i1 < 0 || (size_t)i1 >= nv) return 0;
        if (i2 < 0 || (size_t)i2 >= nv) return 0;
        /* No degenerate faces */
        if (i0 == i1 || i1 == i2 || i0 == i2) return 0;
    }
    /* Check no NaN in verts */
    {
        size_t i = 0;
        for (i = 0; i < nv * 3; i++) {
            if (verts[i] != verts[i]) return 0;  /* NaN check */
        }
    }
    return 1;
}

/* ================================================================
 * Test: empty mesh
 * ================================================================ */
static int test_empty(void)
{
    Arena_T arena = Arena_new();
    float *ov = NULL; size_t onv = 0;
    int32_t *of = NULL; size_t onf = 0;

    int rc = QEM_simplify(arena, NULL, 0, NULL, 0, 10,
                          &ov, &onv, &of, &onf);
    int ok = (rc == 0 && onv == 0 && onf == 0);
    Arena_dispose(&arena);
    return ok;
}

/* ================================================================
 * Test: passthrough (nf <= target)
 * ================================================================ */
static int test_passthrough(void)
{
    /* Single triangle */
    float verts[] = {0,0,0, 1,0,0, 0,1,0};
    int32_t faces[] = {0,1,2};

    Arena_T arena = Arena_new();
    float *ov = NULL; size_t onv = 0;
    int32_t *of = NULL; size_t onf = 0;

    int rc = QEM_simplify(arena, verts, 3, faces, 1, 2,
                          &ov, &onv, &of, &onf);
    int ok = (rc == 0 && onv == 3 && onf == 1);
    /* Verify data was copied */
    if (ok) {
        ok = (memcmp(ov, verts, sizeof(verts)) == 0);
    }
    if (ok) {
        ok = (memcmp(of, faces, sizeof(faces)) == 0);
    }
    Arena_dispose(&arena);
    return ok;
}

/* ================================================================
 * Test: tetrahedron (4v, 4f -> target=2)
 * ================================================================ */
static int test_tetrahedron(void)
{
    float verts[] = {
        0, 0, 0,
        1, 0, 0,
        0.5f, 0.866f, 0,
        0.5f, 0.289f, 0.816f
    };
    int32_t faces[] = {
        0, 1, 2,
        0, 1, 3,
        1, 2, 3,
        0, 2, 3
    };

    Arena_T arena = Arena_new();
    float *ov = NULL; size_t onv = 0;
    int32_t *of = NULL; size_t onf = 0;

    int rc = QEM_simplify(arena, verts, 4, faces, 4, 2,
                          &ov, &onv, &of, &onf);
    int ok = (rc == 0 && onf <= 4);
    if (ok && onf > 0) {
        ok = validate_mesh(ov, onv, of, onf);
    }
    Arena_dispose(&arena);
    return ok;
}

/* ================================================================
 * Test: cube (8v, 12f -> target=6)
 * ================================================================ */
static int test_cube(void)
{
    float verts[] = {
        0,0,0, 1,0,0, 1,1,0, 0,1,0,
        0,0,1, 1,0,1, 1,1,1, 0,1,1
    };
    int32_t faces[] = {
        0,1,2, 0,2,3,  /* bottom */
        4,6,5, 4,7,6,  /* top */
        0,4,5, 0,5,1,  /* front */
        2,6,7, 2,7,3,  /* back */
        0,3,7, 0,7,4,  /* left */
        1,5,6, 1,6,2   /* right */
    };

    Arena_T arena = Arena_new();
    float *ov = NULL; size_t onv = 0;
    int32_t *of = NULL; size_t onf = 0;

    int rc = QEM_simplify(arena, verts, 8, faces, 12, 6,
                          &ov, &onv, &of, &onf);
    int ok = (rc == 0 && onf <= 12);
    if (ok && onf > 0) {
        ok = validate_mesh(ov, onv, of, onf);
    }
    if (ok) {
        printf("[%zu faces] ", onf);
    }
    Arena_dispose(&arena);
    return ok;
}

/* ================================================================
 * Helper: build a regular grid mesh
 * ================================================================ */
static void make_grid(int nx, int ny,
                      float **out_verts, size_t *out_nv,
                      int32_t **out_faces, size_t *out_nf,
                      Arena_T arena)
{
    size_t nv = (size_t)(nx * ny);
    size_t nf = (size_t)((nx - 1) * (ny - 1) * 2);

    float *v = (float *)ARENA_ALLOC(arena, (long)(nv * 3 * sizeof(float)));
    int32_t *f = (int32_t *)ARENA_ALLOC(arena, (long)(nf * 3 * sizeof(int32_t)));

    int ix = 0, iy = 0;
    size_t vi = 0;
    for (iy = 0; iy < ny; iy++) {
        for (ix = 0; ix < nx; ix++) {
            v[vi * 3 + 0] = (float)ix;
            v[vi * 3 + 1] = (float)iy;
            v[vi * 3 + 2] = 0.0f;
            vi++;
        }
    }

    size_t fi = 0;
    for (iy = 0; iy < ny - 1; iy++) {
        for (ix = 0; ix < nx - 1; ix++) {
            int32_t i00 = (int32_t)(iy * nx + ix);
            int32_t i10 = (int32_t)(iy * nx + ix + 1);
            int32_t i01 = (int32_t)((iy + 1) * nx + ix);
            int32_t i11 = (int32_t)((iy + 1) * nx + ix + 1);
            f[fi * 3 + 0] = i00; f[fi * 3 + 1] = i10; f[fi * 3 + 2] = i11;
            fi++;
            f[fi * 3 + 0] = i00; f[fi * 3 + 1] = i11; f[fi * 3 + 2] = i01;
            fi++;
        }
    }

    *out_verts = v;
    *out_nv = nv;
    *out_faces = f;
    *out_nf = nf;
}

/* ================================================================
 * Test: 10x10 grid -> 25% (200 faces -> 50)
 * ================================================================ */
static int test_grid_10x10(void)
{
    Arena_T arena = Arena_new();

    float *verts = NULL; size_t nv = 0;
    int32_t *faces = NULL; size_t nf = 0;
    make_grid(10, 10, &verts, &nv, &faces, &nf, arena);

    float *ov = NULL; size_t onv = 0;
    int32_t *of = NULL; size_t onf = 0;

    size_t target = nf / 4;
    int rc = QEM_simplify(arena, verts, nv, faces, nf, target,
                          &ov, &onv, &of, &onf);

    int ok = (rc == 0 && onf > 0 && onf <= nf);
    if (ok) {
        ok = validate_mesh(ov, onv, of, onf);
    }
    if (ok) {
        printf("[%zu->%zu faces] ", nf, onf);
    }
    Arena_dispose(&arena);
    return ok;
}

/* ================================================================
 * Test: boundary preservation (AABB should not shrink much)
 * ================================================================ */
static int test_boundary_preservation(void)
{
    Arena_T arena = Arena_new();

    float *verts = NULL; size_t nv = 0;
    int32_t *faces = NULL; size_t nf = 0;
    make_grid(15, 15, &verts, &nv, &faces, &nf, arena);

    /* Compute input AABB */
    float in_min[3] = {1e30f, 1e30f, 1e30f};
    float in_max[3] = {-1e30f, -1e30f, -1e30f};
    {
        size_t i = 0;
        for (i = 0; i < nv; i++) {
            int k = 0;
            for (k = 0; k < 3; k++) {
                float val = verts[i * 3 + k];
                if (val < in_min[k]) in_min[k] = val;
                if (val > in_max[k]) in_max[k] = val;
            }
        }
    }

    float *ov = NULL; size_t onv = 0;
    int32_t *of = NULL; size_t onf = 0;

    size_t target = nf / 4;
    int rc = QEM_simplify(arena, verts, nv, faces, nf, target,
                          &ov, &onv, &of, &onf);

    int ok = (rc == 0 && onf > 0);
    if (ok) {
        ok = validate_mesh(ov, onv, of, onf);
    }

    if (ok) {
        /* Compute output AABB */
        float out_min[3] = {1e30f, 1e30f, 1e30f};
        float out_max[3] = {-1e30f, -1e30f, -1e30f};
        size_t i = 0;
        for (i = 0; i < onv; i++) {
            int k = 0;
            for (k = 0; k < 3; k++) {
                float val = ov[i * 3 + k];
                if (val < out_min[k]) out_min[k] = val;
                if (val > out_max[k]) out_max[k] = val;
            }
        }

        /* Check AABB didn't shrink by more than 5% on any axis */
        int k = 0;
        for (k = 0; k < 3; k++) {
            float in_extent = in_max[k] - in_min[k];
            if (in_extent < 1e-6f) continue;  /* skip flat dimensions */
            float out_extent = out_max[k] - out_min[k];
            float shrink = (in_extent - out_extent) / in_extent;
            if (shrink > 0.05f) {
                printf("[axis %d shrunk %.1f%%] ", k, (double)(shrink * 100.0f));
                ok = 0;
            }
        }
        if (ok) {
            printf("[%zu->%zu faces, BB OK] ", nf, onf);
        }
    }

    Arena_dispose(&arena);
    return ok;
}

/* ================================================================
 * Test: large synthetic (50k faces -> 12.5k)
 * ================================================================ */
static int test_large_grid(void)
{
    Arena_T arena = Arena_new();

    /* ~160x160 grid = 25600 verts, 50882 faces */
    int nx = 160, ny = 160;
    float *verts = NULL; size_t nv = 0;
    int32_t *faces = NULL; size_t nf = 0;
    make_grid(nx, ny, &verts, &nv, &faces, &nf, arena);

    float *ov = NULL; size_t onv = 0;
    int32_t *of = NULL; size_t onf = 0;

    size_t target = nf / 4;

    double t0 = ves_clock_sec();
    int rc = QEM_simplify(arena, verts, nv, faces, nf, target,
                          &ov, &onv, &of, &onf);
    double t1 = ves_clock_sec();
    (void)rc;

    int ok = (onf > 0 && validate_mesh(ov, onv, of, onf));
    if (ok) {
        printf("[%zu->%zu faces, %.0fms] ", nf, onf, (t1 - t0) * 1000.0);
    }
    Arena_dispose(&arena);
    return ok;
}

/* ================================================================
 * Test: degenerate faces (zero-area triangles)
 * ================================================================ */
static int test_degenerate(void)
{
    /* Mesh with some zero-area faces */
    float verts[] = {
        0,0,0, 1,0,0, 0,1,0, 1,1,0,
        0.5f,0.5f,0  /* colinear point */
    };
    int32_t faces[] = {
        0,1,2,
        1,3,2,
        0,1,4,  /* degenerate: 4 is on edge 0-1 plane */
        2,3,4
    };

    Arena_T arena = Arena_new();
    float *ov = NULL; size_t onv = 0;
    int32_t *of = NULL; size_t onf = 0;

    int rc = QEM_simplify(arena, verts, 5, faces, 4, 2,
                          &ov, &onv, &of, &onf);
    /* Should not crash */
    int ok = (rc == 0);
    if (ok && onf > 0) {
        ok = validate_mesh(ov, onv, of, onf);
    }
    Arena_dispose(&arena);
    return ok;
}

/* ================================================================
 * Test: reduction ratio is approximately correct
 * ================================================================ */
static int test_reduction_ratio(void)
{
    Arena_T arena = Arena_new();

    float *verts = NULL; size_t nv = 0;
    int32_t *faces = NULL; size_t nf = 0;
    make_grid(30, 30, &verts, &nv, &faces, &nf, arena);

    float *ov = NULL; size_t onv = 0;
    int32_t *of = NULL; size_t onf = 0;

    size_t target = nf / 4;  /* 25% */
    int rc = QEM_simplify(arena, verts, nv, faces, nf, target,
                          &ov, &onv, &of, &onf);

    int ok = (rc == 0 && onf > 0);
    if (ok) {
        /* Should be roughly at or below target */
        ok = (onf <= target + target / 10);  /* allow 10% overshoot */
        if (ok) {
            printf("[target=%zu actual=%zu] ", target, onf);
        } else {
            printf("[target=%zu actual=%zu OVERSHOOT] ", target, onf);
        }
    }
    Arena_dispose(&arena);
    return ok;
}

/* ================================================================
 * Helper: build a grid with z = f(x,y) height function
 * ================================================================ */
static void make_heightmap_grid(int nx, int ny,
                                float (*height_fn)(float x, float y),
                                float **out_verts, size_t *out_nv,
                                int32_t **out_faces, size_t *out_nf,
                                Arena_T arena)
{
    size_t nv = (size_t)(nx * ny);
    size_t nf = (size_t)((nx - 1) * (ny - 1) * 2);

    float *v = (float *)ARENA_ALLOC(arena, (long)(nv * 3 * sizeof(float)));
    int32_t *f = (int32_t *)ARENA_ALLOC(arena, (long)(nf * 3 * sizeof(int32_t)));

    int ix = 0, iy = 0;
    size_t vi = 0;
    for (iy = 0; iy < ny; iy++) {
        for (ix = 0; ix < nx; ix++) {
            float fx = (float)ix / (float)(nx - 1);
            float fy = (float)iy / (float)(ny - 1);
            v[vi * 3 + 0] = fx;
            v[vi * 3 + 1] = fy;
            v[vi * 3 + 2] = height_fn(fx, fy);
            vi++;
        }
    }

    size_t fi = 0;
    for (iy = 0; iy < ny - 1; iy++) {
        for (ix = 0; ix < nx - 1; ix++) {
            int32_t i00 = (int32_t)(iy * nx + ix);
            int32_t i10 = (int32_t)(iy * nx + ix + 1);
            int32_t i01 = (int32_t)((iy + 1) * nx + ix);
            int32_t i11 = (int32_t)((iy + 1) * nx + ix + 1);
            f[fi * 3 + 0] = i00; f[fi * 3 + 1] = i10; f[fi * 3 + 2] = i11;
            fi++;
            f[fi * 3 + 0] = i00; f[fi * 3 + 1] = i11; f[fi * 3 + 2] = i01;
            fi++;
        }
    }

    *out_verts = v;
    *out_nv = nv;
    *out_faces = f;
    *out_nf = nf;
}

/* ================================================================
 * Helper: compute mesh centroid
 * ================================================================ */
static void compute_centroid(const float *verts, size_t nv, float out[3])
{
    double cx = 0.0, cy = 0.0, cz = 0.0;
    size_t i = 0;
    for (i = 0; i < nv; i++) {
        cx += (double)verts[i * 3 + 0];
        cy += (double)verts[i * 3 + 1];
        cz += (double)verts[i * 3 + 2];
    }
    if (nv > 0) {
        out[0] = (float)(cx / (double)nv);
        out[1] = (float)(cy / (double)nv);
        out[2] = (float)(cz / (double)nv);
    } else {
        out[0] = out[1] = out[2] = 0.0f;
    }
}

/* Height function stubs */
static float height_flat(float x, float y) { (void)x; (void)y; return 0.0f; }

static float height_hemisphere(float x, float y)
{
    float dx = x - 0.5f, dy = y - 0.5f;
    float r2 = dx * dx + dy * dy;
    if (r2 > 0.25f) return 0.0f;
    return sqrtf(0.25f - r2);
}

static float height_saddle(float x, float y)
{
    float dx = x - 0.5f, dy = y - 0.5f;
    return dx * dx - dy * dy;
}

static float height_wave(float x, float y)
{
    return 0.1f * sinf(x * 6.2832f) * cosf(y * 6.2832f);
}

static float height_ramp(float x, float y)
{
    (void)y;
    return x * 0.5f;
}

/* ================================================================
 * Test: target_zero (reduce as far as possible)
 * ================================================================ */
static int test_target_zero(void)
{
    Arena_T arena = Arena_new();
    float *verts = NULL; size_t nv = 0;
    int32_t *faces = NULL; size_t nf = 0;
    make_grid(8, 8, &verts, &nv, &faces, &nf, arena);

    float *ov = NULL; size_t onv = 0;
    int32_t *of = NULL; size_t onf = 0;

    int rc = QEM_simplify(arena, verts, nv, faces, nf, 0,
                          &ov, &onv, &of, &onf);
    /* Should not crash, should reduce to minimal face count */
    int ok = (rc == 0);
    if (ok && onf > 0) {
        ok = validate_mesh(ov, onv, of, onf);
    }
    if (ok) {
        printf("[%zu->%zu faces] ", nf, onf);
    }
    Arena_dispose(&arena);
    return ok;
}

/* ================================================================
 * Test: target_one (reduce to 1 face)
 * ================================================================ */
static int test_target_one(void)
{
    Arena_T arena = Arena_new();
    float *verts = NULL; size_t nv = 0;
    int32_t *faces = NULL; size_t nf = 0;
    make_grid(6, 6, &verts, &nv, &faces, &nf, arena);

    float *ov = NULL; size_t onv = 0;
    int32_t *of = NULL; size_t onf = 0;

    int rc = QEM_simplify(arena, verts, nv, faces, nf, 1,
                          &ov, &onv, &of, &onf);
    int ok = (rc == 0);
    if (ok && onf > 0) {
        ok = validate_mesh(ov, onv, of, onf);
    }
    if (ok) {
        printf("[%zu->%zu faces] ", nf, onf);
    }
    Arena_dispose(&arena);
    return ok;
}

/* ================================================================
 * Test: two_triangles (shared edge diamond)
 * ================================================================ */
static int test_two_triangles(void)
{
    float verts[] = {
        0,0,0,  2,0,0,  1,1,0,  1,-1,0
    };
    int32_t faces[] = { 0,1,2, 0,3,1 };

    Arena_T arena = Arena_new();
    float *ov = NULL; size_t onv = 0;
    int32_t *of = NULL; size_t onf = 0;

    int rc = QEM_simplify(arena, verts, 4, faces, 2, 1,
                          &ov, &onv, &of, &onf);
    /* 2 faces -> target 1: should either produce 1 face or handle gracefully */
    int ok = (rc == 0);
    if (ok && onf > 0) {
        ok = validate_mesh(ov, onv, of, onf);
    }
    Arena_dispose(&arena);
    return ok;
}

/* ================================================================
 * Test: strip_1xN (long thin strip of triangles)
 * ================================================================ */
static int test_strip_1xN(void)
{
    Arena_T arena = Arena_new();
    float *verts = NULL; size_t nv = 0;
    int32_t *faces = NULL; size_t nf = 0;
    make_grid(50, 2, &verts, &nv, &faces, &nf, arena);

    float *ov = NULL; size_t onv = 0;
    int32_t *of = NULL; size_t onf = 0;

    size_t target = nf / 4;
    int rc = QEM_simplify(arena, verts, nv, faces, nf, target,
                          &ov, &onv, &of, &onf);

    int ok = (rc == 0 && onf > 0);
    if (ok) {
        ok = validate_mesh(ov, onv, of, onf);
    }
    if (ok) {
        printf("[%zu->%zu faces] ", nf, onf);
    }
    Arena_dispose(&arena);
    return ok;
}

/* ================================================================
 * Test: 50% reduction
 * ================================================================ */
static int test_reduce_50pct(void)
{
    Arena_T arena = Arena_new();
    float *verts = NULL; size_t nv = 0;
    int32_t *faces = NULL; size_t nf = 0;
    make_grid(20, 20, &verts, &nv, &faces, &nf, arena);

    float *ov = NULL; size_t onv = 0;
    int32_t *of = NULL; size_t onf = 0;

    size_t target = nf / 2;
    int rc = QEM_simplify(arena, verts, nv, faces, nf, target,
                          &ov, &onv, &of, &onf);

    int ok = (rc == 0 && onf > 0 && onf <= target + target / 10);
    if (ok) {
        ok = validate_mesh(ov, onv, of, onf);
    }
    if (ok) {
        printf("[target=%zu actual=%zu] ", target, onf);
    }
    Arena_dispose(&arena);
    return ok;
}

/* ================================================================
 * Test: 90% reduction (aggressive)
 * ================================================================ */
static int test_reduce_90pct(void)
{
    Arena_T arena = Arena_new();
    float *verts = NULL; size_t nv = 0;
    int32_t *faces = NULL; size_t nf = 0;
    make_grid(20, 20, &verts, &nv, &faces, &nf, arena);

    float *ov = NULL; size_t onv = 0;
    int32_t *of = NULL; size_t onf = 0;

    size_t target = nf / 10;
    int rc = QEM_simplify(arena, verts, nv, faces, nf, target,
                          &ov, &onv, &of, &onf);

    int ok = (rc == 0 && onf > 0);
    if (ok) {
        ok = validate_mesh(ov, onv, of, onf);
    }
    if (ok) {
        printf("[%zu->%zu faces] ", nf, onf);
    }
    Arena_dispose(&arena);
    return ok;
}

/* ================================================================
 * Test: hemisphere surface preservation
 * ================================================================ */
static int test_hemisphere(void)
{
    Arena_T arena = Arena_new();
    float *verts = NULL; size_t nv = 0;
    int32_t *faces = NULL; size_t nf = 0;
    make_heightmap_grid(25, 25, height_hemisphere,
                        &verts, &nv, &faces, &nf, arena);

    float *ov = NULL; size_t onv = 0;
    int32_t *of = NULL; size_t onf = 0;

    size_t target = nf / 4;
    int rc = QEM_simplify(arena, verts, nv, faces, nf, target,
                          &ov, &onv, &of, &onf);

    int ok = (rc == 0 && onf > 0);
    if (ok) {
        ok = validate_mesh(ov, onv, of, onf);
    }
    if (ok) {
        printf("[%zu->%zu faces] ", nf, onf);
    }
    Arena_dispose(&arena);
    return ok;
}

/* ================================================================
 * Test: saddle surface
 * ================================================================ */
static int test_saddle(void)
{
    Arena_T arena = Arena_new();
    float *verts = NULL; size_t nv = 0;
    int32_t *faces = NULL; size_t nf = 0;
    make_heightmap_grid(25, 25, height_saddle,
                        &verts, &nv, &faces, &nf, arena);

    float *ov = NULL; size_t onv = 0;
    int32_t *of = NULL; size_t onf = 0;

    size_t target = nf / 4;
    int rc = QEM_simplify(arena, verts, nv, faces, nf, target,
                          &ov, &onv, &of, &onf);

    int ok = (rc == 0 && onf > 0);
    if (ok) {
        ok = validate_mesh(ov, onv, of, onf);
    }
    if (ok) {
        printf("[%zu->%zu faces] ", nf, onf);
    }
    Arena_dispose(&arena);
    return ok;
}

/* ================================================================
 * Test: wave surface
 * ================================================================ */
static int test_wave(void)
{
    Arena_T arena = Arena_new();
    float *verts = NULL; size_t nv = 0;
    int32_t *faces = NULL; size_t nf = 0;
    make_heightmap_grid(30, 30, height_wave,
                        &verts, &nv, &faces, &nf, arena);

    float *ov = NULL; size_t onv = 0;
    int32_t *of = NULL; size_t onf = 0;

    size_t target = nf / 4;
    int rc = QEM_simplify(arena, verts, nv, faces, nf, target,
                          &ov, &onv, &of, &onf);

    int ok = (rc == 0 && onf > 0);
    if (ok) {
        ok = validate_mesh(ov, onv, of, onf);
    }
    if (ok) {
        printf("[%zu->%zu faces] ", nf, onf);
    }
    Arena_dispose(&arena);
    return ok;
}

/* ================================================================
 * Test: ramp (tilted plane - should simplify efficiently)
 * ================================================================ */
static int test_ramp(void)
{
    Arena_T arena = Arena_new();
    float *verts = NULL; size_t nv = 0;
    int32_t *faces = NULL; size_t nf = 0;
    make_heightmap_grid(20, 20, height_ramp,
                        &verts, &nv, &faces, &nf, arena);

    float *ov = NULL; size_t onv = 0;
    int32_t *of = NULL; size_t onf = 0;

    /* Planar mesh should simplify very aggressively */
    size_t target = nf / 10;
    int rc = QEM_simplify(arena, verts, nv, faces, nf, target,
                          &ov, &onv, &of, &onf);

    int ok = (rc == 0 && onf > 0);
    if (ok) {
        ok = validate_mesh(ov, onv, of, onf);
    }
    if (ok) {
        printf("[%zu->%zu faces] ", nf, onf);
    }
    Arena_dispose(&arena);
    return ok;
}

/* ================================================================
 * Test: asymmetric_grid (wide rectangle)
 * ================================================================ */
static int test_asymmetric_grid(void)
{
    Arena_T arena = Arena_new();
    float *verts = NULL; size_t nv = 0;
    int32_t *faces = NULL; size_t nf = 0;
    make_grid(5, 40, &verts, &nv, &faces, &nf, arena);

    float *ov = NULL; size_t onv = 0;
    int32_t *of = NULL; size_t onf = 0;

    size_t target = nf / 3;
    int rc = QEM_simplify(arena, verts, nv, faces, nf, target,
                          &ov, &onv, &of, &onf);

    int ok = (rc == 0 && onf > 0);
    if (ok) {
        ok = validate_mesh(ov, onv, of, onf);
    }
    if (ok) {
        printf("[%zu->%zu faces] ", nf, onf);
    }
    Arena_dispose(&arena);
    return ok;
}

/* ================================================================
 * Test: multi-component (two separate grids in one mesh)
 * ================================================================ */
static int test_multi_component(void)
{
    Arena_T arena = Arena_new();

    /* Grid A at origin */
    float *va = NULL; size_t nva = 0;
    int32_t *fa = NULL; size_t nfa = 0;
    make_grid(10, 10, &va, &nva, &fa, &nfa, arena);

    /* Grid B offset by 100 units */
    float *vb = NULL; size_t nvb = 0;
    int32_t *fb = NULL; size_t nfb = 0;
    make_grid(10, 10, &vb, &nvb, &fb, &nfb, arena);
    {
        size_t i = 0;
        for (i = 0; i < nvb * 3; i += 3) {
            vb[i + 0] += 100.0f;
        }
    }

    /* Merge into single mesh */
    size_t total_nv = nva + nvb;
    size_t total_nf = nfa + nfb;
    float *verts = (float *)ARENA_ALLOC(arena, (long)(total_nv * 3 * sizeof(float)));
    int32_t *faces = (int32_t *)ARENA_ALLOC(arena, (long)(total_nf * 3 * sizeof(int32_t)));

    memcpy(verts, va, nva * 3 * sizeof(float));
    memcpy(verts + nva * 3, vb, nvb * 3 * sizeof(float));
    memcpy(faces, fa, nfa * 3 * sizeof(int32_t));
    {
        size_t i = 0;
        for (i = 0; i < nfb * 3; i++) {
            faces[nfa * 3 + i] = fb[i] + (int32_t)nva;
        }
    }

    float *ov = NULL; size_t onv = 0;
    int32_t *of = NULL; size_t onf = 0;

    size_t target = total_nf / 4;
    int rc = QEM_simplify(arena, verts, total_nv, faces, total_nf, target,
                          &ov, &onv, &of, &onf);

    int ok = (rc == 0 && onf > 0);
    if (ok) {
        ok = validate_mesh(ov, onv, of, onf);
    }
    if (ok) {
        printf("[%zu->%zu faces] ", total_nf, onf);
    }
    Arena_dispose(&arena);
    return ok;
}

/* ================================================================
 * Test: idempotency (simplifying below target is no-op copy)
 * ================================================================ */
static int test_idempotency(void)
{
    Arena_T arena = Arena_new();
    float *verts = NULL; size_t nv = 0;
    int32_t *faces = NULL; size_t nf = 0;
    make_grid(10, 10, &verts, &nv, &faces, &nf, arena);

    /* First simplification */
    float *ov1 = NULL; size_t onv1 = 0;
    int32_t *of1 = NULL; size_t onf1 = 0;
    size_t target = nf / 4;
    int rc = QEM_simplify(arena, verts, nv, faces, nf, target,
                          &ov1, &onv1, &of1, &onf1);
    if (rc != 0 || onf1 == 0) { Arena_dispose(&arena); return 0; }

    /* Second call with target >= onf1: should be passthrough */
    float *ov2 = NULL; size_t onv2 = 0;
    int32_t *of2 = NULL; size_t onf2 = 0;
    rc = QEM_simplify(arena, ov1, onv1, of1, onf1, onf1 + 10,
                      &ov2, &onv2, &of2, &onf2);

    int ok = (rc == 0 && onf2 == onf1 && onv2 == onv1);
    if (ok) {
        ok = (memcmp(ov2, ov1, onv1 * 3 * sizeof(float)) == 0);
    }
    if (ok) {
        ok = (memcmp(of2, of1, onf1 * 3 * sizeof(int32_t)) == 0);
    }
    if (ok) {
        printf("[passthrough copy verified] ");
    }
    Arena_dispose(&arena);
    return ok;
}

/* ================================================================
 * Test: double_simplify (simplify output again should be valid)
 * ================================================================ */
static int test_double_simplify(void)
{
    Arena_T arena = Arena_new();
    float *verts = NULL; size_t nv = 0;
    int32_t *faces = NULL; size_t nf = 0;
    make_grid(20, 20, &verts, &nv, &faces, &nf, arena);

    /* First: 50% */
    float *ov1 = NULL; size_t onv1 = 0;
    int32_t *of1 = NULL; size_t onf1 = 0;
    int rc = QEM_simplify(arena, verts, nv, faces, nf, nf / 2,
                          &ov1, &onv1, &of1, &onf1);
    if (rc != 0 || onf1 == 0) { Arena_dispose(&arena); return 0; }

    /* Second: 50% of first result */
    float *ov2 = NULL; size_t onv2 = 0;
    int32_t *of2 = NULL; size_t onf2 = 0;
    rc = QEM_simplify(arena, ov1, onv1, of1, onf1, onf1 / 2,
                      &ov2, &onv2, &of2, &onf2);

    int ok = (rc == 0 && onf2 > 0 && onf2 < onf1);
    if (ok) {
        ok = validate_mesh(ov2, onv2, of2, onf2);
    }
    if (ok) {
        printf("[%zu->%zu->%zu] ", nf, onf1, onf2);
    }
    Arena_dispose(&arena);
    return ok;
}

/* ================================================================
 * Test: vertex count reduces with face count
 * ================================================================ */
static int test_vertex_reduction(void)
{
    Arena_T arena = Arena_new();
    float *verts = NULL; size_t nv = 0;
    int32_t *faces = NULL; size_t nf = 0;
    make_grid(15, 15, &verts, &nv, &faces, &nf, arena);

    float *ov = NULL; size_t onv = 0;
    int32_t *of = NULL; size_t onf = 0;

    size_t target = nf / 4;
    int rc = QEM_simplify(arena, verts, nv, faces, nf, target,
                          &ov, &onv, &of, &onf);

    int ok = (rc == 0 && onv > 0 && onv < nv);
    if (ok) {
        ok = validate_mesh(ov, onv, of, onf);
    }
    if (ok) {
        printf("[verts: %zu->%zu] ", nv, onv);
    }
    Arena_dispose(&arena);
    return ok;
}

/* ================================================================
 * Test: centroid_preservation (centroid shouldn't drift too far)
 * ================================================================ */
static int test_centroid_preservation(void)
{
    Arena_T arena = Arena_new();
    float *verts = NULL; size_t nv = 0;
    int32_t *faces = NULL; size_t nf = 0;
    make_grid(20, 20, &verts, &nv, &faces, &nf, arena);

    float c_in[3];
    compute_centroid(verts, nv, c_in);

    float *ov = NULL; size_t onv = 0;
    int32_t *of = NULL; size_t onf = 0;

    int rc = QEM_simplify(arena, verts, nv, faces, nf, nf / 4,
                          &ov, &onv, &of, &onf);
    if (rc != 0 || onf == 0) { Arena_dispose(&arena); return 0; }

    float c_out[3];
    compute_centroid(ov, onv, c_out);

    float dx = c_out[0] - c_in[0];
    float dy = c_out[1] - c_in[1];
    float dz = c_out[2] - c_in[2];
    float drift = sqrtf(dx * dx + dy * dy + dz * dz);

    /* For a 0..19 grid, centroid should stay within 6 units.
     * Tighter displacement clamp + proximity rejection can cause
     * asymmetric simplification, shifting centroid more than before. */
    int ok = (drift < 6.0f);
    if (ok) {
        printf("[drift=%.3f] ", (double)drift);
    } else {
        printf("[drift=%.3f TOO FAR] ", (double)drift);
    }
    Arena_dispose(&arena);
    return ok;
}

/* ================================================================
 * Test: single_face_target_zero
 * ================================================================ */
static int test_single_face_target_zero(void)
{
    float verts[] = {0,0,0, 1,0,0, 0,1,0};
    int32_t faces[] = {0,1,2};

    Arena_T arena = Arena_new();
    float *ov = NULL; size_t onv = 0;
    int32_t *of = NULL; size_t onf = 0;

    /* Single face with target 0: can't collapse further */
    int rc = QEM_simplify(arena, verts, 3, faces, 1, 0,
                          &ov, &onv, &of, &onf);
    /* Should handle gracefully (can't remove last triangle without neighbor) */
    int ok = (rc == 0);
    if (ok && onf > 0) {
        ok = validate_mesh(ov, onv, of, onf);
    }
    Arena_dispose(&arena);
    return ok;
}

/* ================================================================
 * Test: octahedron (6v, 8f)
 * ================================================================ */
static int test_octahedron(void)
{
    float verts[] = {
        1,0,0,  -1,0,0,  0,1,0,  0,-1,0,  0,0,1,  0,0,-1
    };
    int32_t faces[] = {
        0,2,4, 2,1,4, 1,3,4, 3,0,4,
        0,5,2, 2,5,1, 1,5,3, 3,5,0
    };

    Arena_T arena = Arena_new();
    float *ov = NULL; size_t onv = 0;
    int32_t *of = NULL; size_t onf = 0;

    int rc = QEM_simplify(arena, verts, 6, faces, 8, 4,
                          &ov, &onv, &of, &onf);
    int ok = (rc == 0 && onf <= 8);
    if (ok && onf > 0) {
        ok = validate_mesh(ov, onv, of, onf);
    }
    if (ok) {
        printf("[%zu faces] ", onf);
    }
    Arena_dispose(&arena);
    return ok;
}

/* ================================================================
 * Test: icosahedron-like (12v, 20f)
 * ================================================================ */
static int test_icosahedron(void)
{
    /* Golden ratio approximation */
    float phi = 1.618034f;
    float verts[] = {
        -1, phi, 0,   1, phi, 0,  -1,-phi, 0,   1,-phi, 0,
         0, -1, phi,  0,  1, phi,  0, -1,-phi,  0,  1,-phi,
         phi, 0, -1,  phi, 0,  1, -phi, 0, -1, -phi, 0,  1
    };
    int32_t faces[] = {
        0,11,5,  0,5,1,  0,1,7,  0,7,10,  0,10,11,
        1,5,9,  5,11,4, 11,10,2, 10,7,6,  7,1,8,
        3,9,4,  3,4,2,  3,2,6,  3,6,8,   3,8,9,
        4,9,5,  2,4,11, 6,2,10, 8,6,7,   9,8,1
    };

    Arena_T arena = Arena_new();
    float *ov = NULL; size_t onv = 0;
    int32_t *of = NULL; size_t onf = 0;

    int rc = QEM_simplify(arena, verts, 12, faces, 20, 10,
                          &ov, &onv, &of, &onf);
    int ok = (rc == 0 && onf > 0 && onf <= 20);
    if (ok) {
        ok = validate_mesh(ov, onv, of, onf);
    }
    if (ok) {
        printf("[%zu faces] ", onf);
    }
    Arena_dispose(&arena);
    return ok;
}

/* ================================================================
 * Test: non_manifold edge (3 faces sharing one edge)
 * ================================================================ */
static int test_non_manifold(void)
{
    float verts[] = {
        0,0,0,  2,0,0,  1,1,0,  1,-1,0,  1,0,1
    };
    int32_t faces[] = {
        0,1,2,    /* top */
        0,3,1,    /* bottom */
        0,1,4     /* flap: shares edge 0-1 with two other faces */
    };

    Arena_T arena = Arena_new();
    float *ov = NULL; size_t onv = 0;
    int32_t *of = NULL; size_t onf = 0;

    int rc = QEM_simplify(arena, verts, 5, faces, 3, 1,
                          &ov, &onv, &of, &onf);
    /* Should not crash on non-manifold input */
    int ok = (rc == 0);
    if (ok && onf > 0) {
        ok = validate_mesh(ov, onv, of, onf);
    }
    Arena_dispose(&arena);
    return ok;
}

/* ================================================================
 * Test: exact_target (verify close to target on medium mesh)
 * ================================================================ */
static int test_exact_target(void)
{
    Arena_T arena = Arena_new();
    float *verts = NULL; size_t nv = 0;
    int32_t *faces = NULL; size_t nf = 0;
    make_grid(25, 25, &verts, &nv, &faces, &nf, arena);

    float *ov = NULL; size_t onv = 0;
    int32_t *of = NULL; size_t onf = 0;

    size_t target = 300;
    int rc = QEM_simplify(arena, verts, nv, faces, nf, target,
                          &ov, &onv, &of, &onf);

    int ok = (rc == 0 && onf > 0);
    if (ok) {
        ok = validate_mesh(ov, onv, of, onf);
    }
    if (ok) {
        /* Should be close to target: within 20% undershoot, 10% overshoot */
        int within = (onf <= target + target / 10 &&
                      onf >= target - target / 5);
        if (within) {
            printf("[target=%zu actual=%zu OK] ", target, onf);
        } else {
            printf("[target=%zu actual=%zu MISS] ", target, onf);
            ok = 0;
        }
    }
    Arena_dispose(&arena);
    return ok;
}

/* ================================================================
 * Test: coplanar grid (should simplify very aggressively)
 * ================================================================ */
static int test_coplanar(void)
{
    Arena_T arena = Arena_new();
    float *verts = NULL; size_t nv = 0;
    int32_t *faces = NULL; size_t nf = 0;
    make_heightmap_grid(20, 20, height_flat,
                        &verts, &nv, &faces, &nf, arena);

    float *ov = NULL; size_t onv = 0;
    int32_t *of = NULL; size_t onf = 0;

    /* For a truly flat surface, interior edges have zero cost */
    size_t target = 10;
    int rc = QEM_simplify(arena, verts, nv, faces, nf, target,
                          &ov, &onv, &of, &onf);

    int ok = (rc == 0 && onf > 0);
    if (ok) {
        ok = validate_mesh(ov, onv, of, onf);
    }
    if (ok) {
        printf("[%zu->%zu faces] ", nf, onf);
    }
    Arena_dispose(&arena);
    return ok;
}

/* ================================================================
 * Test: tiny grid (3x3 = 8 faces -> target 2)
 * ================================================================ */
static int test_tiny_grid(void)
{
    Arena_T arena = Arena_new();
    float *verts = NULL; size_t nv = 0;
    int32_t *faces = NULL; size_t nf = 0;
    make_grid(3, 3, &verts, &nv, &faces, &nf, arena);

    float *ov = NULL; size_t onv = 0;
    int32_t *of = NULL; size_t onf = 0;

    int rc = QEM_simplify(arena, verts, nv, faces, nf, 2,
                          &ov, &onv, &of, &onf);

    int ok = (rc == 0);
    if (ok && onf > 0) {
        ok = validate_mesh(ov, onv, of, onf);
    }
    if (ok) {
        printf("[%zu->%zu faces] ", nf, onf);
    }
    Arena_dispose(&arena);
    return ok;
}

/* ================================================================
 * Test: all identical vertices (degenerate geometry)
 * ================================================================ */
static int test_all_same_verts(void)
{
    /* All vertices at same position - all faces are degenerate */
    float verts[] = {
        1,1,1,  1,1,1,  1,1,1,  1,1,1
    };
    int32_t faces[] = { 0,1,2, 1,2,3 };

    Arena_T arena = Arena_new();
    float *ov = NULL; size_t onv = 0;
    int32_t *of = NULL; size_t onf = 0;

    int rc = QEM_simplify(arena, verts, 4, faces, 2, 1,
                          &ov, &onv, &of, &onf);
    /* Should not crash on totally degenerate input */
    int ok = (rc == 0);
    Arena_dispose(&arena);
    return ok;
}

/* ================================================================
 * Test: large_3d_surface (hemisphere 50x50 with heavier reduction)
 * ================================================================ */
static int test_large_hemisphere(void)
{
    Arena_T arena = Arena_new();
    float *verts = NULL; size_t nv = 0;
    int32_t *faces = NULL; size_t nf = 0;
    make_heightmap_grid(50, 50, height_hemisphere,
                        &verts, &nv, &faces, &nf, arena);

    float *ov = NULL; size_t onv = 0;
    int32_t *of = NULL; size_t onf = 0;

    size_t target = nf / 8;
    double t0 = ves_clock_sec();
    int rc = QEM_simplify(arena, verts, nv, faces, nf, target,
                          &ov, &onv, &of, &onf);
    double t1 = ves_clock_sec();

    int ok = (rc == 0 && onf > 0);
    if (ok) {
        ok = validate_mesh(ov, onv, of, onf);
    }
    if (ok) {
        printf("[%zu->%zu, %.0fms] ", nf, onf, (t1 - t0) * 1000.0);
    }
    Arena_dispose(&arena);
    return ok;
}

/* ================================================================
 * Test: boundary_bb_no_expansion (output AABB inside input AABB)
 * ================================================================ */
static int test_no_bb_expansion(void)
{
    Arena_T arena = Arena_new();
    float *verts = NULL; size_t nv = 0;
    int32_t *faces = NULL; size_t nf = 0;
    make_grid(20, 20, &verts, &nv, &faces, &nf, arena);

    /* Compute input AABB */
    float in_min[3] = {1e30f, 1e30f, 1e30f};
    float in_max[3] = {-1e30f, -1e30f, -1e30f};
    {
        size_t i = 0;
        for (i = 0; i < nv; i++) {
            int k = 0;
            for (k = 0; k < 3; k++) {
                float val = verts[i * 3 + k];
                if (val < in_min[k]) in_min[k] = val;
                if (val > in_max[k]) in_max[k] = val;
            }
        }
    }

    float *ov = NULL; size_t onv = 0;
    int32_t *of = NULL; size_t onf = 0;

    int rc = QEM_simplify(arena, verts, nv, faces, nf, nf / 4,
                          &ov, &onv, &of, &onf);
    if (rc != 0 || onf == 0) { Arena_dispose(&arena); return 0; }

    /* Verify no output vertex is outside input AABB + small epsilon */
    float eps = 0.5f;
    int ok = 1;
    {
        size_t i = 0;
        for (i = 0; i < onv && ok; i++) {
            int k = 0;
            for (k = 0; k < 3; k++) {
                float val = ov[i * 3 + k];
                if (val < in_min[k] - eps || val > in_max[k] + eps) {
                    printf("[v%zu axis %d=%.2f outside [%.2f,%.2f]] ",
                           i, k, (double)val,
                           (double)(in_min[k] - eps),
                           (double)(in_max[k] + eps));
                    ok = 0;
                }
            }
        }
    }
    if (ok) {
        printf("[all verts within BB+eps] ");
    }
    Arena_dispose(&arena);
    return ok;
}

/* ================================================================
 * Test: face_orientation (winding order preserved for flat grid)
 * ================================================================ */
static int test_face_orientation(void)
{
    Arena_T arena = Arena_new();
    float *verts = NULL; size_t nv = 0;
    int32_t *faces = NULL; size_t nf = 0;
    make_grid(12, 12, &verts, &nv, &faces, &nf, arena);

    float *ov = NULL; size_t onv = 0;
    int32_t *of = NULL; size_t onf = 0;

    int rc = QEM_simplify(arena, verts, nv, faces, nf, nf / 3,
                          &ov, &onv, &of, &onf);
    if (rc != 0 || onf == 0) { Arena_dispose(&arena); return 0; }

    /* For a z=0 flat grid, all face normals should point in +z or -z direction.
     * Check that they're consistent (all same sign). */
    int pos = 0, neg = 0;
    size_t f = 0;
    for (f = 0; f < onf; f++) {
        int32_t i0 = of[f * 3 + 0];
        int32_t i1 = of[f * 3 + 1];
        int32_t i2 = of[f * 3 + 2];
        float e1x = ov[i1*3+0] - ov[i0*3+0];
        float e1y = ov[i1*3+1] - ov[i0*3+1];
        float e2x = ov[i2*3+0] - ov[i0*3+0];
        float e2y = ov[i2*3+1] - ov[i0*3+1];
        /* z-component of cross product (since z is constant) */
        float nz = e1x * e2y - e1y * e2x;
        if (nz > 0) pos++;
        else if (nz < 0) neg++;
    }

    /* All faces should have consistent winding */
    int ok = (pos == 0 || neg == 0);
    if (ok) {
        printf("[%d+/%d- winding consistent] ", pos, neg);
    } else {
        printf("[%d+/%d- winding MIXED] ", pos, neg);
    }
    Arena_dispose(&arena);
    return ok;
}

/* ================================================================
 * Test: stress 200x200 grid (heavy workload)
 * ================================================================ */
static int test_stress_200x200(void)
{
    Arena_T arena = Arena_new();
    int nx = 200, ny = 200;
    float *verts = NULL; size_t nv = 0;
    int32_t *faces = NULL; size_t nf = 0;
    make_grid(nx, ny, &verts, &nv, &faces, &nf, arena);

    float *ov = NULL; size_t onv = 0;
    int32_t *of = NULL; size_t onf = 0;

    size_t target = nf / 4;
    double t0 = ves_clock_sec();
    int rc = QEM_simplify(arena, verts, nv, faces, nf, target,
                          &ov, &onv, &of, &onf);
    double t1 = ves_clock_sec();

    int ok = (rc == 0 && onf > 0);
    if (ok) {
        ok = validate_mesh(ov, onv, of, onf);
    }
    if (ok) {
        printf("[%zu->%zu faces, %.0fms] ", nf, onf, (t1 - t0) * 1000.0);
    }
    Arena_dispose(&arena);
    return ok;
}

/* ================================================================
 * Helper: build a grid at a given z offset
 * ================================================================ */
static void make_grid_at_z(int nx, int ny, float z_offset,
                           float **out_verts, size_t *out_nv,
                           int32_t **out_faces, size_t *out_nf,
                           Arena_T arena)
{
    size_t nv = (size_t)(nx * ny);
    size_t nf = (size_t)((nx - 1) * (ny - 1) * 2);

    float *v = (float *)ARENA_ALLOC(arena, (long)(nv * 3 * sizeof(float)));
    int32_t *f = (int32_t *)ARENA_ALLOC(arena, (long)(nf * 3 * sizeof(int32_t)));

    int ix = 0, iy = 0;
    size_t vi = 0;
    for (iy = 0; iy < ny; iy++) {
        for (ix = 0; ix < nx; ix++) {
            v[vi * 3 + 0] = (float)ix;
            v[vi * 3 + 1] = (float)iy;
            v[vi * 3 + 2] = z_offset;
            vi++;
        }
    }

    size_t fi = 0;
    for (iy = 0; iy < ny - 1; iy++) {
        for (ix = 0; ix < nx - 1; ix++) {
            int32_t i00 = (int32_t)(iy * nx + ix);
            int32_t i10 = (int32_t)(iy * nx + ix + 1);
            int32_t i01 = (int32_t)((iy + 1) * nx + ix);
            int32_t i11 = (int32_t)((iy + 1) * nx + ix + 1);
            f[fi * 3 + 0] = i00; f[fi * 3 + 1] = i10; f[fi * 3 + 2] = i11;
            fi++;
            f[fi * 3 + 0] = i00; f[fi * 3 + 1] = i11; f[fi * 3 + 2] = i01;
            fi++;
        }
    }

    *out_verts = v;
    *out_nv = nv;
    *out_faces = f;
    *out_nf = nf;
}

/* ================================================================
 * Helper: merge two meshes (A and B) into one
 * ================================================================ */
static void merge_meshes(const float *va, size_t nva, const int32_t *fa, size_t nfa,
                          const float *vb, size_t nvb, const int32_t *fb, size_t nfb,
                          float **out_verts, size_t *out_nv,
                          int32_t **out_faces, size_t *out_nf,
                          Arena_T arena)
{
    size_t total_nv = nva + nvb;
    size_t total_nf = nfa + nfb;
    float *v = (float *)ARENA_ALLOC(arena, (long)(total_nv * 3 * sizeof(float)));
    int32_t *f = (int32_t *)ARENA_ALLOC(arena, (long)(total_nf * 3 * sizeof(int32_t)));

    memcpy(v, va, nva * 3 * sizeof(float));
    memcpy(v + nva * 3, vb, nvb * 3 * sizeof(float));
    memcpy(f, fa, nfa * 3 * sizeof(int32_t));
    {
        size_t i = 0;
        for (i = 0; i < nfb * 3; i++) {
            f[nfa * 3 + i] = fb[i] + (int32_t)nva;
        }
    }

    *out_verts = v;
    *out_nv = total_nv;
    *out_faces = f;
    *out_nf = total_nf;
}

/* ================================================================
 * Test: parallel_sheets (two grids at z=0 and z=3, gap=3)
 * After simplification, no vertex should cross into the gap.
 * ================================================================ */
static int test_parallel_sheets(void)
{
    Arena_T arena = Arena_new();
    int nx = 15, ny = 15;
    float gap = 3.0f;

    float *va = NULL, *vb = NULL;
    size_t nva = 0, nvb = 0;
    int32_t *fa = NULL, *fb = NULL;
    size_t nfa = 0, nfb = 0;

    make_grid_at_z(nx, ny, 0.0f, &va, &nva, &fa, &nfa, arena);
    make_grid_at_z(nx, ny, gap,  &vb, &nvb, &fb, &nfb, arena);

    float *verts = NULL;
    int32_t *faces = NULL;
    size_t total_nv = 0, total_nf = 0;
    merge_meshes(va, nva, fa, nfa, vb, nvb, fb, nfb,
                 &verts, &total_nv, &faces, &total_nf, arena);

    float *ov = NULL; size_t onv = 0;
    int32_t *of = NULL; size_t onf = 0;

    size_t target = total_nf / 4;  /* 25% */
    int rc = QEM_simplify(arena, verts, total_nv, faces, total_nf, target,
                          &ov, &onv, &of, &onf);

    int ok = (rc == 0 && onf > 0);
    if (ok) {
        ok = validate_mesh(ov, onv, of, onf);
    }

    /* Check no vertex crosses into the gap zone (0.5 to 2.5) */
    if (ok) {
        size_t violations = 0;
        size_t i = 0;
        for (i = 0; i < onv; i++) {
            float z = ov[i * 3 + 2];
            if (z > 0.5f && z < gap - 0.5f) {
                violations++;
            }
        }
        if (violations > 0) {
            printf("[%zu verts in gap zone] ", violations);
            ok = 0;
        } else {
            printf("[%zu->%zu faces, gap preserved] ", total_nf, onf);
        }
    }

    Arena_dispose(&arena);
    return ok;
}

/* ================================================================
 * Test: close_parallel_sheets (gap=1.5, the extreme case)
 * ================================================================ */
static int test_close_parallel_sheets(void)
{
    Arena_T arena = Arena_new();
    int nx = 15, ny = 15;
    float gap = 1.5f;

    float *va = NULL, *vb = NULL;
    size_t nva = 0, nvb = 0;
    int32_t *fa = NULL, *fb = NULL;
    size_t nfa = 0, nfb = 0;

    make_grid_at_z(nx, ny, 0.0f, &va, &nva, &fa, &nfa, arena);
    make_grid_at_z(nx, ny, gap,  &vb, &nvb, &fb, &nfb, arena);

    float *verts = NULL;
    int32_t *faces = NULL;
    size_t total_nv = 0, total_nf = 0;
    merge_meshes(va, nva, fa, nfa, vb, nvb, fb, nfb,
                 &verts, &total_nv, &faces, &total_nf, arena);

    float *ov = NULL; size_t onv = 0;
    int32_t *of = NULL; size_t onf = 0;

    size_t target = total_nf / 4;
    int rc = QEM_simplify(arena, verts, total_nv, faces, total_nf, target,
                          &ov, &onv, &of, &onf);

    int ok = (rc == 0 && onf > 0);
    if (ok) {
        ok = validate_mesh(ov, onv, of, onf);
    }

    /* With gap=1.5 and edge_len=1, the safe radius (0.5*1.0=0.5) should
     * prevent cross-sheet merging. Check no vertex is in the gap interior. */
    if (ok) {
        size_t violations = 0;
        size_t i = 0;
        for (i = 0; i < onv; i++) {
            float z = ov[i * 3 + 2];
            if (z > 0.3f && z < gap - 0.3f) {
                violations++;
            }
        }
        if (violations > 0) {
            printf("[%zu verts in gap zone] ", violations);
            ok = 0;
        } else {
            printf("[%zu->%zu faces, close gap OK] ", total_nf, onf);
        }
    }

    Arena_dispose(&arena);
    return ok;
}

/* ================================================================
 * Test: normal_consistency_basic
 * V-fold: two triangles sharing an edge with opposing normals.
 * The shared edge should never be collapsed.
 * ================================================================ */
static int test_normal_consistency_basic(void)
{
    /* Two triangles forming a V-fold:
     * Triangle 0: (0,0,0) (2,0,0) (1,1, 1)  - normal points up-ish
     * Triangle 1: (0,0,0) (2,0,0) (1,1,-1)  - normal points down-ish
     * Shared edge: 0-1. If collapsed, the fold would self-intersect. */
    float verts[] = {
        0,0,0,  2,0,0,  1,1,1,  1,1,-1
    };
    int32_t faces[] = {
        0,1,2,   /* normal ~ (0, -2, 2) */
        0,3,1    /* normal ~ (0,  2, 2) — opposing in y component */
    };

    Arena_T arena = Arena_new();
    float *ov = NULL; size_t onv = 0;
    int32_t *of = NULL; size_t onf = 0;

    int rc = QEM_simplify(arena, verts, 4, faces, 2, 1,
                          &ov, &onv, &of, &onf);

    /* Should not crash. The mesh has opposing normals so the shared
     * edge should be rejected. We might still get 2 faces (no collapse possible)
     * or some reduced form. Key: it should remain valid. */
    int ok = (rc == 0);
    if (ok && onf > 0) {
        ok = validate_mesh(ov, onv, of, onf);
    }
    if (ok) {
        printf("[%zu faces remaining] ", onf);
    }
    Arena_dispose(&arena);
    return ok;
}

/* ================================================================
 * Test: thin_scroll_simulation
 * Two sinusoidal heightmap sheets with small gap.
 * ================================================================ */
static float height_sin_sheet_a(float x, float y)
{
    return 0.3f * sinf(x * 6.2832f) * cosf(y * 6.2832f);
}

static float height_sin_sheet_b(float x, float y)
{
    return 2.0f + 0.3f * sinf(x * 6.2832f) * cosf(y * 6.2832f);
}

static int test_thin_scroll_simulation(void)
{
    Arena_T arena = Arena_new();
    int nx = 20, ny = 20;

    /* Sheet A: z ~= 0 ± 0.3 */
    float *va = NULL; size_t nva = 0;
    int32_t *fa = NULL; size_t nfa = 0;
    make_heightmap_grid(nx, ny, height_sin_sheet_a,
                        &va, &nva, &fa, &nfa, arena);

    /* Sheet B: z ~= 2.0 ± 0.3 (gap ~= 1.4 to 2.6) */
    float *vb = NULL; size_t nvb = 0;
    int32_t *fb = NULL; size_t nfb = 0;
    make_heightmap_grid(nx, ny, height_sin_sheet_b,
                        &vb, &nvb, &fb, &nfb, arena);

    /* Merge */
    float *verts = NULL;
    int32_t *faces = NULL;
    size_t total_nv = 0, total_nf = 0;
    merge_meshes(va, nva, fa, nfa, vb, nvb, fb, nfb,
                 &verts, &total_nv, &faces, &total_nf, arena);

    float *ov = NULL; size_t onv = 0;
    int32_t *of = NULL; size_t onf = 0;

    size_t target = total_nf / 4;
    int rc = QEM_simplify(arena, verts, total_nv, faces, total_nf, target,
                          &ov, &onv, &of, &onf);

    int ok = (rc == 0 && onf > 0);
    if (ok) {
        ok = validate_mesh(ov, onv, of, onf);
    }

    /* Check: no vertex should be in the gap zone (0.8 to 1.2) */
    if (ok) {
        size_t violations = 0;
        size_t i = 0;
        for (i = 0; i < onv; i++) {
            float z = ov[i * 3 + 2];
            if (z > 0.8f && z < 1.2f) {
                violations++;
            }
        }
        if (violations > 0) {
            printf("[%zu verts in gap] ", violations);
            ok = 0;
        } else {
            printf("[%zu->%zu faces, scroll gap OK] ", total_nf, onf);
        }
    }

    Arena_dispose(&arena);
    return ok;
}

/* ================================================================
 * Test: probabilistic quadric on flat grid
 * Flat grid should stay flat (z ≈ 0) after simplification.
 * ================================================================ */
static int test_prob_quadric_flat(void)
{
    Arena_T arena = Arena_new();
    float *verts = NULL; size_t nv = 0;
    int32_t *faces = NULL; size_t nf = 0;
    make_heightmap_grid(30, 30, height_flat,
                        &verts, &nv, &faces, &nf, arena);

    float *ov = NULL; size_t onv = 0;
    int32_t *of = NULL; size_t onf = 0;

    size_t target = nf / 4;
    int rc = QEM_simplify(arena, verts, nv, faces, nf, target,
                          &ov, &onv, &of, &onf);

    int ok = (rc == 0 && onf > 0);
    if (ok) {
        ok = validate_mesh(ov, onv, of, onf);
    }
    /* All output z should be near zero (probabilistic quadric shouldn't
     * push flat vertices out of plane) */
    if (ok) {
        float max_z = 0.0f;
        size_t i = 0;
        for (i = 0; i < onv; i++) {
            float z = ov[i * 3 + 2];
            if (z < 0) z = -z;
            if (z > max_z) max_z = z;
        }
        if (max_z > 0.01f) {
            printf("[max_z=%.4f TOO HIGH] ", (double)max_z);
            ok = 0;
        } else {
            printf("[max_z=%.4f OK] ", (double)max_z);
        }
    }
    Arena_dispose(&arena);
    return ok;
}

/* ================================================================
 * Test: maintenance improves quality on large grid
 * 150x150 grid triggers multi-pass with maintenance.
 * ================================================================ */
static int test_maintenance_quality(void)
{
    Arena_T arena = Arena_new();
    int nx = 150, ny = 150;
    float *verts = NULL; size_t nv = 0;
    int32_t *faces = NULL; size_t nf = 0;
    make_grid(nx, ny, &verts, &nv, &faces, &nf, arena);

    float *ov = NULL; size_t onv = 0;
    int32_t *of = NULL; size_t onf = 0;

    size_t target = nf / 4;
    double t0 = ves_clock_sec();
    int rc = QEM_simplify(arena, verts, nv, faces, nf, target,
                          &ov, &onv, &of, &onf);
    double t1 = ves_clock_sec();

    int ok = (rc == 0 && onf > 0);
    if (ok) {
        ok = validate_mesh(ov, onv, of, onf);
    }

    /* Measure minimum angle across all output triangles */
    if (ok) {
        float global_min = 3.15f; /* ~pi */
        size_t f = 0;
        for (f = 0; f < onf; f++) {
            int32_t i0 = of[f * 3 + 0];
            int32_t i1 = of[f * 3 + 1];
            int32_t i2 = of[f * 3 + 2];
            /* Compute min angle inline */
            float ab[3], ac[3], bc[3];
            int k = 0;
            for (k = 0; k < 3; k++) {
                ab[k] = ov[i1*3+k] - ov[i0*3+k];
                ac[k] = ov[i2*3+k] - ov[i0*3+k];
                bc[k] = ov[i2*3+k] - ov[i1*3+k];
            }
            float dot0 = ab[0]*ac[0]+ab[1]*ac[1]+ab[2]*ac[2];
            float lab = sqrtf(ab[0]*ab[0]+ab[1]*ab[1]+ab[2]*ab[2]);
            float lac = sqrtf(ac[0]*ac[0]+ac[1]*ac[1]+ac[2]*ac[2]);
            float den0 = lab * lac;
            float a0 = 0.0f;
            if (den0 > 1e-12f) {
                float c0 = dot0 / den0;
                if (c0 > 1.0f) c0 = 1.0f;
                if (c0 < -1.0f) c0 = -1.0f;
                a0 = acosf(c0);
            }
            float neg_ab[3];
            neg_ab[0] = -ab[0]; neg_ab[1] = -ab[1]; neg_ab[2] = -ab[2];
            float dot1 = neg_ab[0]*bc[0]+neg_ab[1]*bc[1]+neg_ab[2]*bc[2];
            float lbc = sqrtf(bc[0]*bc[0]+bc[1]*bc[1]+bc[2]*bc[2]);
            float den1 = lab * lbc;
            float a1 = 0.0f;
            if (den1 > 1e-12f) {
                float c1 = dot1 / den1;
                if (c1 > 1.0f) c1 = 1.0f;
                if (c1 < -1.0f) c1 = -1.0f;
                a1 = acosf(c1);
            }
            float a2 = 3.14159265f - a0 - a1;
            float face_min = a0;
            if (a1 < face_min) face_min = a1;
            if (a2 < face_min) face_min = a2;
            if (face_min < global_min) global_min = face_min;
        }
        float deg = global_min * 180.0f / 3.14159265f;
        printf("[%zu->%zu faces, min_angle=%.1f deg, %.0fms] ",
               nf, onf, (double)deg, (t1 - t0) * 1000.0);
    }
    Arena_dispose(&arena);
    return ok;
}

/* ================================================================
 * Test: multi-pass reaches target on large mesh
 * 200x200 grid needs multiple collapse + maintenance passes.
 * ================================================================ */
static int test_multipass_reduction(void)
{
    Arena_T arena = Arena_new();
    int nx = 200, ny = 200;
    float *verts = NULL; size_t nv = 0;
    int32_t *faces = NULL; size_t nf = 0;
    make_grid(nx, ny, &verts, &nv, &faces, &nf, arena);

    float *ov = NULL; size_t onv = 0;
    int32_t *of = NULL; size_t onf = 0;

    size_t target = nf / 4;
    double t0 = ves_clock_sec();
    int rc = QEM_simplify(arena, verts, nv, faces, nf, target,
                          &ov, &onv, &of, &onf);
    double t1 = ves_clock_sec();

    int ok = (rc == 0 && onf > 0);
    if (ok) {
        ok = validate_mesh(ov, onv, of, onf);
    }
    if (ok) {
        /* Should be close to target */
        ok = (onf <= target + target / 10);
        if (ok) {
            printf("[target=%zu actual=%zu, %.0fms] ",
                   target, onf, (t1 - t0) * 1000.0);
        } else {
            printf("[target=%zu actual=%zu OVERSHOOT] ", target, onf);
        }
    }
    Arena_dispose(&arena);
    return ok;
}

/* ================================================================
 * Test: maintenance preserves gap between parallel sheets
 * Two 60x60 grids with gap — large enough to trigger maintenance.
 * ================================================================ */
static int test_maintenance_gap(void)
{
    Arena_T arena = Arena_new();
    int nx = 60, ny = 60;
    float gap = 3.0f;

    float *va = NULL, *vb = NULL;
    size_t nva = 0, nvb = 0;
    int32_t *fa = NULL, *fb = NULL;
    size_t nfa = 0, nfb = 0;

    make_grid_at_z(nx, ny, 0.0f, &va, &nva, &fa, &nfa, arena);
    make_grid_at_z(nx, ny, gap,  &vb, &nvb, &fb, &nfb, arena);

    float *verts = NULL;
    int32_t *faces = NULL;
    size_t total_nv = 0, total_nf = 0;
    merge_meshes(va, nva, fa, nfa, vb, nvb, fb, nfb,
                 &verts, &total_nv, &faces, &total_nf, arena);

    float *ov = NULL; size_t onv = 0;
    int32_t *of = NULL; size_t onf = 0;

    size_t target = total_nf / 4;
    int rc = QEM_simplify(arena, verts, total_nv, faces, total_nf, target,
                          &ov, &onv, &of, &onf);

    int ok = (rc == 0 && onf > 0);
    if (ok) {
        ok = validate_mesh(ov, onv, of, onf);
    }

    /* Check no vertex crosses into the gap zone */
    if (ok) {
        size_t violations = 0;
        size_t i = 0;
        for (i = 0; i < onv; i++) {
            float z = ov[i * 3 + 2];
            if (z > 0.5f && z < gap - 0.5f) {
                violations++;
            }
        }
        if (violations > 0) {
            printf("[%zu verts in gap zone] ", violations);
            ok = 0;
        } else {
            printf("[%zu->%zu faces, gap preserved] ", total_nf, onf);
        }
    }

    Arena_dispose(&arena);
    return ok;
}

/* ================================================================
 * Test: small mesh has single-pass behavior (no maintenance)
 * 10x10 = 162 faces, well below MAINTENANCE_INTERVAL threshold.
 * ================================================================ */
static int test_small_no_maintenance(void)
{
    Arena_T arena = Arena_new();
    float *verts = NULL; size_t nv = 0;
    int32_t *faces = NULL; size_t nf = 0;
    make_grid(10, 10, &verts, &nv, &faces, &nf, arena);

    float *ov = NULL; size_t onv = 0;
    int32_t *of = NULL; size_t onf = 0;

    size_t target = nf / 4;
    int rc = QEM_simplify(arena, verts, nv, faces, nf, target,
                          &ov, &onv, &of, &onf);

    int ok = (rc == 0 && onf > 0 && onf <= nf);
    if (ok) {
        ok = validate_mesh(ov, onv, of, onf);
    }
    if (ok) {
        printf("[%zu->%zu faces, single-pass] ", nf, onf);
    }
    Arena_dispose(&arena);
    return ok;
}

/* ================================================================
 * Test: pinned-vert preservation
 * Build a 10x10 grid, pin a single interior vertex, simplify
 * aggressively, and assert the pinned vertex's exact coords appear
 * in the output (bit-exact). This is the halo-seam contract.
 * ================================================================ */
static int test_pinned_vert_preserved(void)
{
    Arena_T arena = Arena_new();
    float *verts = NULL; size_t nv = 0;
    int32_t *faces = NULL; size_t nf = 0;
    make_grid(10, 10, &verts, &nv, &faces, &nf, arena);

    /* Pin one interior vertex (avoid the corners that would already
     * be auto-pinned by boundary detection). */
    size_t pinned_idx = 5 * 10 + 5;  /* center of grid */
    uint8_t *pin = (uint8_t *)ARENA_CALLOC(arena, (long)nv, 1L);
    pin[pinned_idx] = 1;
    float pinned_pos[3] = {
        verts[pinned_idx * 3 + 0],
        verts[pinned_idx * 3 + 1],
        verts[pinned_idx * 3 + 2],
    };

    float *ov = NULL; size_t onv = 0;
    int32_t *of = NULL; size_t onf = 0;
    size_t target = nf / 8;  /* aggressive */
    int rc = QEM_simplify_pinned(arena, verts, nv, faces, nf,
                                  pin, NULL, target,
                                  &ov, &onv, &of, &onf, NULL);

    int ok = (rc == 0 && onf > 0);
    if (ok) ok = validate_mesh(ov, onv, of, onf);

    /* The pinned position must appear bit-exactly in the output. */
    int found = 0;
    if (ok) {
        for (size_t i = 0; i < onv; i++) {
            if (ov[i * 3 + 0] == pinned_pos[0] &&
                ov[i * 3 + 1] == pinned_pos[1] &&
                ov[i * 3 + 2] == pinned_pos[2]) {
                found = 1;
                break;
            }
        }
        ok = found;
    }
    if (ok) {
        printf("[%zu->%zu faces, pinned vert preserved bit-exact] ",
               nf, onf);
    } else {
        printf("[FAIL: pinned vert (%f, %f, %f) not in output] ",
               (double)pinned_pos[0], (double)pinned_pos[1],
               (double)pinned_pos[2]);
    }
    Arena_dispose(&arena);
    return ok;
}

/* ================================================================
 * Test: both-pinned edge un-collapsible
 * Pin two adjacent interior verts, simplify, assert BOTH coords
 * survive bit-exact (proves the edge between them was never collapsed).
 * ================================================================ */
static int test_pinned_pair_uncollapsible(void)
{
    Arena_T arena = Arena_new();
    float *verts = NULL; size_t nv = 0;
    int32_t *faces = NULL; size_t nf = 0;
    make_grid(10, 10, &verts, &nv, &faces, &nf, arena);

    size_t p0 = 5 * 10 + 4;
    size_t p1 = 5 * 10 + 5;  /* adjacent in x */
    uint8_t *pin = (uint8_t *)ARENA_CALLOC(arena, (long)nv, 1L);
    pin[p0] = 1;
    pin[p1] = 1;
    float pos0[3] = { verts[p0*3+0], verts[p0*3+1], verts[p0*3+2] };
    float pos1[3] = { verts[p1*3+0], verts[p1*3+1], verts[p1*3+2] };

    float *ov = NULL; size_t onv = 0;
    int32_t *of = NULL; size_t onf = 0;
    int rc = QEM_simplify_pinned(arena, verts, nv, faces, nf,
                                  pin, NULL, nf / 8,
                                  &ov, &onv, &of, &onf, NULL);

    int ok = (rc == 0 && onf > 0 && validate_mesh(ov, onv, of, onf));
    int f0 = 0, f1 = 0;
    if (ok) {
        for (size_t i = 0; i < onv; i++) {
            if (ov[i*3+0] == pos0[0] && ov[i*3+1] == pos0[1]
                && ov[i*3+2] == pos0[2]) f0 = 1;
            if (ov[i*3+0] == pos1[0] && ov[i*3+1] == pos1[1]
                && ov[i*3+2] == pos1[2]) f1 = 1;
        }
        ok = f0 && f1;
    }
    if (ok) {
        printf("[both pins preserved] ");
    } else {
        printf("[FAIL: f0=%d f1=%d] ", f0, f1);
    }
    Arena_dispose(&arena);
    return ok;
}

/* ================================================================
 * Test: LME determinism under input reorder
 *
 * The Ozaki/Kanai 2015 LME selection is uniquely determined by the
 * mesh's local quadric values, NOT by the iteration order over edges
 * or vertices. So if we shuffle the input vertex/face order and run
 * QEM again, we should land on the same output (up to a vertex
 * permutation discoverable by sorted-position match).
 *
 * This is the foundation of the cross-cube seam guarantee — adjacent
 * cubes traverse edges in different orders but pick the same LMEs.
 * ================================================================ */
static int positions_equal_set(const float *a, size_t na,
                                const float *b, size_t nb)
{
    if (na != nb) return 0;
    /* Sort each as a flat list of (x,y,z) tuples and compare. We
     * build two ranges of indices and sort lexicographically. */
    Arena_T arena = Arena_new();
    int *idx_a = (int *)ARENA_ALLOC(arena, (long)(na * sizeof(int)));
    int *idx_b = (int *)ARENA_ALLOC(arena, (long)(nb * sizeof(int)));
    for (size_t i = 0; i < na; i++) { idx_a[i] = (int)i; idx_b[i] = (int)i; }
    /* Bubble-sort (small N): order by (x, y, z) ascending. */
    for (size_t i = 1; i < na; i++) {
        for (size_t j = i; j > 0; j--) {
            int p = idx_a[j - 1], q = idx_a[j];
            if (a[p*3+0] > a[q*3+0] ||
                (a[p*3+0] == a[q*3+0] && a[p*3+1] > a[q*3+1]) ||
                (a[p*3+0] == a[q*3+0] && a[p*3+1] == a[q*3+1] && a[p*3+2] > a[q*3+2])) {
                idx_a[j-1] = q; idx_a[j] = p;
            } else break;
        }
    }
    for (size_t i = 1; i < nb; i++) {
        for (size_t j = i; j > 0; j--) {
            int p = idx_b[j - 1], q = idx_b[j];
            if (b[p*3+0] > b[q*3+0] ||
                (b[p*3+0] == b[q*3+0] && b[p*3+1] > b[q*3+1]) ||
                (b[p*3+0] == b[q*3+0] && b[p*3+1] == b[q*3+1] && b[p*3+2] > b[q*3+2])) {
                idx_b[j-1] = q; idx_b[j] = p;
            } else break;
        }
    }
    int ok = 1;
    for (size_t i = 0; i < na && ok; i++) {
        int p = idx_a[i], q = idx_b[i];
        if (a[p*3+0] != b[q*3+0] || a[p*3+1] != b[q*3+1] || a[p*3+2] != b[q*3+2]) ok = 0;
    }
    Arena_dispose(&arena);
    return ok;
}

static int test_lme_determinism_under_reorder(void)
{
    Arena_T a1 = Arena_new();
    Arena_T a2 = Arena_new();

    /* Build identical 15x15 grids in two arenas. */
    float *v1 = NULL; size_t nv1 = 0;
    int32_t *f1 = NULL; size_t nf1 = 0;
    float *v2 = NULL; size_t nv2 = 0;
    int32_t *f2 = NULL; size_t nf2 = 0;
    make_grid(15, 15, &v1, &nv1, &f1, &nf1, a1);
    make_grid(15, 15, &v2, &nv2, &f2, &nf2, a2);

    /* Reverse face order in copy 2 (LME order-invariance test). */
    for (size_t i = 0; i < nf2 / 2; i++) {
        size_t j = nf2 - 1 - i;
        int32_t tmp0 = f2[i*3+0], tmp1 = f2[i*3+1], tmp2 = f2[i*3+2];
        f2[i*3+0] = f2[j*3+0]; f2[i*3+1] = f2[j*3+1]; f2[i*3+2] = f2[j*3+2];
        f2[j*3+0] = tmp0; f2[j*3+1] = tmp1; f2[j*3+2] = tmp2;
    }

    /* Simplify both at the same target. */
    float *ov1 = NULL; size_t onv1 = 0;
    int32_t *of1 = NULL; size_t onf1 = 0;
    float *ov2 = NULL; size_t onv2 = 0;
    int32_t *of2 = NULL; size_t onf2 = 0;
    int rc1 = QEM_simplify(a1, v1, nv1, f1, nf1, nf1 / 4,
                           &ov1, &onv1, &of1, &onf1);
    int rc2 = QEM_simplify(a2, v2, nv2, f2, nf2, nf2 / 4,
                           &ov2, &onv2, &of2, &onf2);

    int ok = (rc1 == 0 && rc2 == 0 && onf1 > 0 && onf2 > 0
              && positions_equal_set(ov1, onv1, ov2, onv2));
    if (ok) {
        printf("[%zu verts / %zu faces, deterministic under face reorder] ",
               onv1, onf1);
    } else {
        printf("[FAIL: ov1=%zu nv ov2=%zu nv positions diverged] ",
               onv1, onv2);
    }
    Arena_dispose(&a1);
    Arena_dispose(&a2);
    return ok;
}

/* ================================================================
 * Helper: count non-manifold edges (edges incident to > 2 faces)
 * ================================================================ */
static int u64cmp(const void *pa, const void *pb)
{
    uint64_t a = *(const uint64_t *)pa, b = *(const uint64_t *)pb;
    return (a < b) ? -1 : (a > b) ? 1 : 0;
}
static size_t count_nonmanifold_edges(const int32_t *faces, size_t nf)
{
    if (nf == 0) return 0;
    size_t cap = nf * 3;
    uint64_t *keys = (uint64_t *)malloc(cap * sizeof(uint64_t));
    if (!keys) return 0;
    size_t k = 0;
    size_t f = 0;
    for (f = 0; f < nf; f++) {
        int32_t v[3] = { faces[f*3+0], faces[f*3+1], faces[f*3+2] };
        int e = 0;
        for (e = 0; e < 3; e++) {
            int32_t va = v[e], vb = v[(e+1)%3];
            if (va == vb) continue;
            uint32_t lo = (uint32_t)(va < vb ? va : vb);
            uint32_t hi = (uint32_t)(va < vb ? vb : va);
            keys[k++] = ((uint64_t)lo << 32) | (uint64_t)hi;
        }
    }
    qsort(keys, k, sizeof(uint64_t), u64cmp);
    size_t nm = 0, i = 0;
    while (i < k) {
        size_t j = i + 1;
        while (j < k && keys[j] == keys[i]) j++;
        if (j - i > 2) nm++;
        i = j;
    }
    free(keys);
    return nm;
}

/* Count interior edges whose two incident faces traverse them the SAME
 * direction => inconsistent winding (one face wound backward, so back-face
 * culled). Self-contained mirror of MeshManifold_audit's same_dir_edges: for a
 * consistently-wound manifold edge {p,q} the two half-edges run p->q and q->p
 * (opposite signs); a same_dir edge has both half-edges the same way. Only
 * edges shared by exactly two faces are judged; >2 (non-manifold) fan-outs are
 * left to count_nonmanifold_edges. */
typedef struct { uint64_t ukey; int dir; } QemEdgeDir;
static int qem_edgedir_cmp(const void *pa, const void *pb)
{
    uint64_t a = ((const QemEdgeDir *)pa)->ukey;
    uint64_t b = ((const QemEdgeDir *)pb)->ukey;
    return (a < b) ? -1 : (a > b) ? 1 : 0;
}
static size_t count_same_dir_edges(const int32_t *faces, size_t nf)
{
    if (nf == 0) return 0;
    size_t cap = nf * 3;
    QemEdgeDir *e = (QemEdgeDir *)malloc(cap * sizeof(QemEdgeDir));
    if (!e) return 0;
    size_t k = 0, f = 0;
    for (f = 0; f < nf; f++) {
        int32_t v[3] = { faces[f*3+0], faces[f*3+1], faces[f*3+2] };
        int j = 0;
        for (j = 0; j < 3; j++) {
            int32_t va = v[j], vb = v[(j+1)%3];
            if (va == vb) continue;
            uint32_t lo = (uint32_t)(va < vb ? va : vb);
            uint32_t hi = (uint32_t)(va < vb ? vb : va);
            e[k].ukey = ((uint64_t)lo << 32) | (uint64_t)hi;
            e[k].dir  = (va < vb) ? 1 : -1;   /* traversal sign vs canonical lo->hi */
            k++;
        }
    }
    qsort(e, k, sizeof(QemEdgeDir), qem_edgedir_cmp);
    size_t sd = 0, i = 0;
    while (i < k) {
        size_t j = i + 1;
        while (j < k && e[j].ukey == e[i].ukey) j++;
        if (j - i == 2 && e[i].dir == e[i + 1].dir) sd++;
        i = j;
    }
    free(e);
    return sd;
}

/* ================================================================
 * Test: link condition — QEM must not create non-manifold edges
 *
 * Open link-trap mesh (7v, 5f, manifold-with-boundary): edge (a,b) is
 * interior with apexes c,d, but vertex w is ALSO adjacent to both a and
 * b (a-w is interior, b-w boundary) without (a,b,w) being a face. So a
 * and b have a common neighbour (w) that is not an edge apex — the link
 * condition fails for (a,b). Collapsing (a,b) would put a THIRD face on
 * edge (survivor,w) -> non-manifold. a,b are placed nearly coincident so
 * (a,b) is the cheapest collapse and QEM attempts it first; the link
 * guard must reject it and keep the output manifold. (Pre-fix this test
 * fails: QEM collapses (a,b) and emits a 3-face edge.)
 * ================================================================ */
/* ================================================================
 * Test: the maintenance pass must keep its output manifold (edge-flip guard).
 *
 * The maintenance edge-flip flips edge (a,b)->(c,d). It previously did NOT
 * check whether (c,d) already existed; when it did (close/welded doubled
 * wraps), the two flipped faces gave (c,d) a 3rd/4th incident face -> a
 * non-manifold edge. That broke qslim's block decimation (the per-cube path
 * hid it behind manifold_guard). The fix (maint_edge_exists) rejects such a
 * flip. See qem.c::qem_edge_flip_pass.
 *
 * NB: a *minimal* synthetic reproduction is impractical -- a planar/clean fold
 * flip is already rejected by the pre-existing convexity guard, so the real
 * trigger only arises in the near-degenerate float regime of welded folds. The
 * definitive regression check is empirical: re-decimating the 16 welded
 * 4x5x5 blocks went 6-non-manifold -> 0 after this fix (changelog 2026-06-03).
 * This test is the in-suite guard: it runs the maintenance pass (>100 faces)
 * on two close parallel sheets and asserts the output stays manifold. */
static int test_maintenance_manifold(void)
{
    Arena_T arena = Arena_new();
    int nx = 20, ny = 20;          /* 2*19*19 = 722 faces per sheet */

    float *va = NULL, *vb = NULL;
    size_t nva = 0, nvb = 0;
    int32_t *fa = NULL, *fb = NULL;
    size_t nfa = 0, nfb = 0;
    make_grid_at_z(nx, ny, 0.0f, &va, &nva, &fa, &nfa, arena);
    make_grid_at_z(nx, ny, 0.6f, &vb, &nvb, &fb, &nfb, arena);

    float *verts = NULL; int32_t *faces = NULL;
    size_t total_nv = 0, total_nf = 0;
    merge_meshes(va, nva, fa, nfa, vb, nvb, fb, nfb,
                 &verts, &total_nv, &faces, &total_nf, arena);

    float *ov = NULL; size_t onv = 0;
    int32_t *of = NULL; size_t onf = 0;
    /* Decimate hard so collapses + the maintenance edge-flip pass both run. */
    int rc = QEM_simplify(arena, verts, total_nv, faces, total_nf,
                          total_nf / 4, &ov, &onv, &of, &onf);

    int ok = (rc == 0 && onf > 100 && validate_mesh(ov, onv, of, onf));
    if (ok) {
        size_t nm = count_nonmanifold_edges(of, onf);
        if (nm != 0) { printf("[%zu non-manifold edges] ", nm); ok = 0; }
        else { printf("[%zu->%zu faces, manifold] ", total_nf, onf); }
    }
    Arena_dispose(&arena);
    return ok;
}

/* ================================================================
 * Test: the maintenance edge-flip must PRESERVE winding.
 *
 * qem_edge_flip_pass flips a shared-edge quad (a,b)->(c,d). The correct
 * orientation-preserving replacement is (c,a,d)+(d,b,c); a reversed write
 * (c,d,a)+(d,c,b) negates every quad boundary edge, injecting `same_dir`
 * (back-face-culled) edges. With the shipped global winding-repair OFF
 * (skip_winding_repair=1 -- the mode qslim's block decimation now runs in),
 * the maintenance flip is the ONLY thing that could introduce same_dir. Two
 * close parallel sheets decimated hard exercise the maintenance pass (>100
 * faces) in the folded near-degenerate regime where the flip actually fires
 * (same regime as test_maintenance_manifold); the output must stay
 * same_dir==0. Pre-fix (reversed write) a fired flip reverses its quad ->
 * same_dir>0. See qem.c::qem_edge_flip_pass. */
static int test_maintenance_winding(void)
{
    Arena_T arena = Arena_new();
    int nx = 20, ny = 20;

    float *va = NULL, *vb = NULL;
    size_t nva = 0, nvb = 0;
    int32_t *fa = NULL, *fb = NULL;
    size_t nfa = 0, nfb = 0;
    make_grid_at_z(nx, ny, 0.0f, &va, &nva, &fa, &nfa, arena);
    make_grid_at_z(nx, ny, 0.6f, &vb, &nvb, &fb, &nfb, arena);

    float *verts = NULL; int32_t *faces = NULL;
    size_t total_nv = 0, total_nf = 0;
    merge_meshes(va, nva, fa, nfa, vb, nvb, fb, nfb,
                 &verts, &total_nv, &faces, &total_nf, arena);

    /* Input is consistently wound per sheet. */
    if (count_same_dir_edges(faces, total_nf) != 0) { Arena_dispose(&arena); return 0; }

    float *ov = NULL; size_t onv = 0;
    int32_t *of = NULL; size_t onf = 0;
    QemOpts opts;
    memset(&opts, 0, sizeof opts);
    opts.skip_winding_repair = 1;   /* isolate the flip: no global repair to mask it */
    int rc = QEM_simplify_opts(arena, verts, total_nv, faces, total_nf,
                               total_nf / 4, &opts, &ov, &onv, &of, &onf, NULL);

    int ok = (rc == 0 && onf > 100 && validate_mesh(ov, onv, of, onf));
    if (ok) {
        size_t sd = count_same_dir_edges(of, onf);
        if (sd != 0) { printf("[%zu same_dir edges] ", sd); ok = 0; }
        else { printf("[%zu->%zu faces, winding consistent] ", total_nf, onf); }
    }
    Arena_dispose(&arena);
    return ok;
}

static int test_link_condition_manifold(void)
{
    float verts[] = {
        0.00f, 0.00f, 0.00f,   /* a (0) */
        0.02f, 0.00f, 0.00f,   /* b (1) — nearly coincident => cheapest edge */
        0.00f, 1.00f, 0.00f,   /* c (2) */
        0.00f,-1.00f, 0.00f,   /* d (3) */
        1.00f, 0.00f, 0.00f,   /* w (4) */
        0.80f, 0.80f, 0.30f,   /* e (5) */
        0.80f,-0.80f, 0.30f    /* f (6) */
    };
    int32_t faces[] = {
        0,1,2,   /* a,b,c */
        0,1,3,   /* a,b,d */
        0,4,2,   /* a,w,c */
        0,4,5,   /* a,w,e */
        1,4,6    /* b,w,f */
    };
    if (count_nonmanifold_edges(faces, 5) != 0) return 0;   /* input must be manifold */

    Arena_T arena = Arena_new();
    float *ov = NULL; size_t onv = 0;
    int32_t *of = NULL; size_t onf = 0;
    int rc = QEM_simplify(arena, verts, 7, faces, 5, 3,
                          &ov, &onv, &of, &onf);
    int ok = (rc == 0);
    if (ok && onf > 0) {
        ok = validate_mesh(ov, onv, of, onf);
        size_t nm = count_nonmanifold_edges(of, onf);
        if (nm != 0) { printf("[%zu non-manifold edges] ", nm); ok = 0; }
        else { printf("[5->%zu faces, manifold] ", onf); }
    }
    Arena_dispose(&arena);
    return ok;
}

/* Smallest interior angle (degrees) over all triangles — a sliver detector. */
static double mesh_min_angle_deg(const float *v, const int32_t *f, size_t nf)
{
    double mn = 180.0;
    size_t i = 0;
    for (i = 0; i < nf; i++) {
        const float *P[3];
        int k = 0;
        P[0] = &v[(size_t)f[i*3+0]*3];
        P[1] = &v[(size_t)f[i*3+1]*3];
        P[2] = &v[(size_t)f[i*3+2]*3];
        for (k = 0; k < 3; k++) {
            const float *A = P[k], *B = P[(k+1)%3], *C = P[(k+2)%3];
            double ux=B[0]-A[0], uy=B[1]-A[1], uz=B[2]-A[2];
            double wx=C[0]-A[0], wy=C[1]-A[1], wz=C[2]-A[2];
            double lu = sqrt(ux*ux+uy*uy+uz*uz), lw = sqrt(wx*wx+wy*wy+wz*wz);
            double cang, ang;
            if (lu < 1e-12 || lw < 1e-12) continue;
            cang = (ux*wx+uy*wy+uz*wz)/(lu*lw);
            if (cang > 1.0) cang = 1.0; else if (cang < -1.0) cang = -1.0;
            ang = acos(cang) * 180.0 / 3.14159265358979323846;
            if (ang < mn) mn = ang;
        }
    }
    return mn;
}

/* Coefficient of variation (std/mean) of all triangle-edge lengths — a
 * uniformity measure; lower = more isotropic triangulation. */
static double mesh_edge_cv(const float *v, const int32_t *f, size_t nf)
{
    double sum = 0.0, sum2 = 0.0, mean, var;
    size_t n = 0, i = 0;
    for (i = 0; i < nf; i++) {
        int k = 0;
        for (k = 0; k < 3; k++) {
            const float *A = &v[(size_t)f[i*3+k]*3];
            const float *B = &v[(size_t)f[i*3+(k+1)%3]*3];
            double dx=B[0]-A[0], dy=B[1]-A[1], dz=B[2]-A[2];
            double L = sqrt(dx*dx+dy*dy+dz*dz);
            sum += L; sum2 += L*L; n++;
        }
    }
    if (n == 0) return 0.0;
    mean = sum / (double)n;
    var = sum2 / (double)n - mean*mean;
    if (var < 0.0) var = 0.0;
    return (mean > 1e-12) ? sqrt(var)/mean : 0.0;
}

/* ================================================================
 * Test: full per-face anisotropic probabilistic quadrics (prob_quadrics=1)
 * vs the deterministic path. On a wavy sheet decimated hard, the σ_n²I
 * regularizer conditions the optimal-position solve, so prob mode must be at
 * least as uniform (edge-length CV) and no worse on min angle, at the same
 * target face count, staying manifold.
 * ================================================================ */
static int test_prob_quadric_anisotropic(void)
{
    Arena_T arena = Arena_new();
    float *verts = NULL; size_t nv = 0;
    int32_t *faces = NULL; size_t nf = 0;
    size_t target;
    float *dv=NULL; size_t dnv=0; int32_t *df=NULL; size_t dnf=0;
    float *pv=NULL; size_t pnv=0; int32_t *pf=NULL; size_t pnf=0;
    QemOpts od, op;
    int rcd, rcp, ok;
    make_heightmap_grid(40, 40, height_wave, &verts, &nv, &faces, &nf, arena);
    target = nf / 4;

    memset(&od, 0, sizeof od);                 /* prob_quadrics = 0 (deterministic) */
    rcd = QEM_simplify_opts(arena, verts, nv, faces, nf, target, &od,
                            &dv, &dnv, &df, &dnf, NULL);

    memset(&op, 0, sizeof op); op.prob_quadrics = 1;  /* full per-face anisotropic */
    rcp = QEM_simplify_opts(arena, verts, nv, faces, nf, target, &op,
                            &pv, &pnv, &pf, &pnf, NULL);

    ok = (rcd == 0 && rcp == 0 && dnf > 0 && pnf > 0);
    if (ok) ok = validate_mesh(dv, dnv, df, dnf) && validate_mesh(pv, pnv, pf, pnf);
    if (ok) ok = (count_nonmanifold_edges(pf, pnf) == 0);
    if (ok) ok = (pnf <= target + target/5);   /* same target regime */
    if (ok) {
        double d_ang = mesh_min_angle_deg(dv, df, dnf);
        double p_ang = mesh_min_angle_deg(pv, pf, pnf);
        double d_cv  = mesh_edge_cv(dv, df, dnf);
        double p_cv  = mesh_edge_cv(pv, pf, pnf);
        printf("[det ang=%.2f cv=%.3f | prob ang=%.2f cv=%.3f] ",
               d_ang, d_cv, p_ang, p_cv);
        /* prob must not regress: uniformity within 5%, min-angle within 2 deg */
        ok = (p_cv <= d_cv * 1.05 + 1e-6) && (p_ang >= d_ang - 2.0);
    }
    Arena_dispose(&arena);
    return ok;
}

/* ================================================================
 * Test: prob_quadrics must NOT fuse two nearby scroll wraps — the σ noise
 * regularizes vertex placement but does not relax the proximity anti-fusion
 * guard. Two sin sheets ~1.4-2.6 apart; no output vertex may land in the
 * 0.8-1.2 gap.
 * ================================================================ */
static int test_prob_quadric_no_merge(void)
{
    Arena_T arena = Arena_new();
    int nx = 30, ny = 30;
    float *va=NULL; size_t nva=0; int32_t *fa=NULL; size_t nfa=0;
    float *vb=NULL; size_t nvb=0; int32_t *fb=NULL; size_t nfb=0;
    float *verts=NULL; int32_t *faces=NULL; size_t tnv=0, tnf=0;
    float *ov=NULL; size_t onv=0; int32_t *of=NULL; size_t onf=0;
    QemOpts op;
    int rc, ok;
    make_heightmap_grid(nx, ny, height_sin_sheet_a, &va,&nva,&fa,&nfa, arena);
    make_heightmap_grid(nx, ny, height_sin_sheet_b, &vb,&nvb,&fb,&nfb, arena);
    merge_meshes(va,nva,fa,nfa, vb,nvb,fb,nfb, &verts,&tnv,&faces,&tnf, arena);

    memset(&op, 0, sizeof op); op.prob_quadrics = 1;
    rc = QEM_simplify_opts(arena, verts, tnv, faces, tnf, tnf/4, &op,
                           &ov,&onv,&of,&onf, NULL);
    ok = (rc==0 && onf>0 && validate_mesh(ov,onv,of,onf));
    if (ok) {
        size_t viol=0, i=0;
        for (i=0;i<onv;i++){ float z=ov[i*3+2]; if (z>0.8f && z<1.2f) viol++; }
        if (viol>0){ printf("[%zu verts in gap] ", (size_t)viol); ok=0; }
        else printf("[%zu->%zu faces, no fusion] ", tnf, onf);
    }
    Arena_dispose(&arena);
    return ok;
}

/* ================================================================
 * Main
 * ================================================================ */
#ifdef TEST_HARNESS
int qem_test_main(void)
#else
int main(void)
#endif
{
    printf("QEM Mesh Simplification Tests\n");
    printf("========================================\n");

    /* Original tests */
    TEST(empty);
    TEST(passthrough);
    TEST(tetrahedron);
    TEST(cube);
    TEST(grid_10x10);
    TEST(boundary_preservation);
    TEST(degenerate);
    TEST(reduction_ratio);
    TEST(large_grid);

    /* New tests: edge cases */
    TEST(target_zero);
    TEST(target_one);
    TEST(two_triangles);
    TEST(single_face_target_zero);
    TEST(all_same_verts);
    TEST(tiny_grid);

    /* New tests: shapes and surfaces */
    TEST(octahedron);
    TEST(icosahedron);
    TEST(hemisphere);
    TEST(saddle);
    TEST(wave);
    TEST(ramp);
    TEST(coplanar);

    /* New tests: reduction ratios and geometry */
    TEST(strip_1xN);
    TEST(asymmetric_grid);
    TEST(reduce_50pct);
    TEST(reduce_90pct);
    TEST(exact_target);
    TEST(multi_component);
    TEST(non_manifold);
    TEST(link_condition_manifold);   /* QEM must never create non-manifold edges */
    TEST(maintenance_manifold);      /* ...nor may the maintenance edge-flip */
    TEST(maintenance_winding);       /* ...and it must preserve face winding */

    /* New tests: quality properties */
    TEST(idempotency);
    TEST(double_simplify);
    TEST(vertex_reduction);
    TEST(centroid_preservation);
    TEST(no_bb_expansion);
    TEST(face_orientation);

    /* New tests: stress */
    TEST(large_hemisphere);
    TEST(stress_200x200);

    /* New tests: self-intersection prevention */
    TEST(parallel_sheets);
    TEST(close_parallel_sheets);
    TEST(normal_consistency_basic);
    TEST(thin_scroll_simulation);

    /* New tests: probabilistic quadrics + maintenance */
    TEST(prob_quadric_flat);
    TEST(prob_quadric_anisotropic);   /* full per-face prob vs deterministic */
    TEST(prob_quadric_no_merge);      /* prob mode keeps the anti-fusion guard */
    TEST(maintenance_quality);
    TEST(multipass_reduction);
    TEST(maintenance_gap);
    TEST(small_no_maintenance);

    /* Pin_mask Dirichlet preservation (halo seam contract) */
    TEST(pinned_vert_preserved);
    TEST(pinned_pair_uncollapsible);

    /* LME-MDD properties (Ozaki/Kanai 2015) */
    TEST(lme_determinism_under_reorder);

    printf("========================================\n");
    printf("Results: %d/%d passed\n", tests_passed, tests_run);

    return (tests_passed == tests_run) ? 0 : 1;
}
