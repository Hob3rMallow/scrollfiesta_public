/* test_common.c - Exercise all common modules. */
#include "common/ves_platform.h"
#include "common/arena.h"
#include "common/except.h"
#include "common/mesh_types.h"
#include "common/pipeline_constants.h"
#include "common/union_find.h"
#include "common/pca.h"
#include "common/csr.h"
#include "common/kdtree.h"
#include "common/bfs.h"
#include "common/tiff_io.h"
#include "common/halo_loader.h"
#include "common/ply_io.h"
#include "common/obj_io.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

static int tests_run = 0;
static int tests_passed = 0;

#define TEST(name) \
    do { \
        tests_run++; \
        printf("  %-40s ", #name); \
        fflush(stdout); \
    } while (0)

#define PASS() \
    do { \
        tests_passed++; \
        printf("[PASS]\n"); \
    } while (0)

#define FAIL(msg) \
    do { \
        printf("[FAIL] %s\n", msg); \
    } while (0)

#define ASSERT_EQ(a, b, msg) \
    do { \
        if ((a) != (b)) { \
            FAIL(msg); \
            return; \
        } \
    } while (0)

#define ASSERT_NEAR(a, b, eps, msg) \
    do { \
        if (fabsf((float)(a) - (float)(b)) > (eps)) { \
            printf("[FAIL] %s (got %f, expected %f)\n", \
                   (msg), (double)(a), (double)(b)); \
            return; \
        } \
    } while (0)

/* ===== Arena Tests ===== */

static void test_arena_basic(void)
{
    TEST(arena_basic);

    Arena_T arena = Arena_new();
    assert(arena != NULL);

    /* Allocate some memory */
    float *data = (float *)ARENA_ALLOC(arena, 100L * (long)sizeof(float));
    assert(data != NULL);
    for (int i = 0; i < 100; i++) {
        data[i] = (float)i;
    }
    assert(data[50] > 49.5f && data[50] < 50.5f);

    /* Calloc: should be zero-initialized */
    int *zeros = (int *)ARENA_CALLOC(arena, 50L, (long)sizeof(int));
    assert(zeros != NULL);
    for (int i = 0; i < 50; i++) {
        assert(zeros[i] == 0);
    }

    Arena_dispose(&arena);
    assert(arena == NULL);
    PASS();
}

static void test_arena_save_restore(void)
{
    TEST(arena_save_restore);

    Arena_T arena = Arena_new();

    /* Allocate some persistent data */
    float *persistent = (float *)ARENA_ALLOC(arena, 100L * (long)sizeof(float));
    persistent[0] = 42.0f;

    /* Save mark */
    Arena_Mark mark = Arena_save(arena);

    /* Allocate scratch data */
    float *scratch = (float *)ARENA_ALLOC(arena, 1000L * (long)sizeof(float));
    scratch[0] = 99.0f;
    (void)scratch;

    /* Restore — scratch is now invalid */
    Arena_restore(arena, mark);

    /* Persistent data should still be valid */
    assert(persistent[0] > 41.5f && persistent[0] < 42.5f);

    /* Can allocate again in the reclaimed space */
    float *reuse = (float *)ARENA_ALLOC(arena, 500L * (long)sizeof(float));
    assert(reuse != NULL);
    reuse[0] = 123.0f;

    Arena_dispose(&arena);
    PASS();
}

static void test_arena_large_alloc(void)
{
    TEST(arena_large_alloc);

    Arena_T arena = Arena_new();

    /* Allocate more than default chunk size (64MB) */
    size_t big_size = 80 * 1024 * 1024; /* 80 MB */
    uint8_t *big = (uint8_t *)ARENA_ALLOC(arena, (long)big_size);
    assert(big != NULL);
    big[0] = 1;
    big[big_size - 1] = 255;

    Arena_dispose(&arena);
    PASS();
}

static void test_arena_free_reuse(void)
{
    TEST(arena_free_reuse);

    Arena_T arena = Arena_new();

    float *a = (float *)ARENA_ALLOC(arena, 100L * (long)sizeof(float));
    a[0] = 1.0f;

    /* Free all but keep arena alive */
    Arena_free(arena);

    /* Allocate again — reuse freed chunks */
    float *b = (float *)ARENA_ALLOC(arena, 100L * (long)sizeof(float));
    assert(b != NULL);
    b[0] = 2.0f;

    Arena_dispose(&arena);
    PASS();
}

/* ===== Exception Tests ===== */

static void test_except_try_catch(void)
{
    TEST(except_try_catch);

    int caught = 0;
    TRY
        RAISE(IO_Failed);
    EXCEPT(IO_Failed)
        caught = 1;
    END_TRY;

    ASSERT_EQ(caught, 1, "IO_Failed should be caught");
    PASS();
}

static void test_except_finally(void)
{
    TEST(except_finally);

    int finalized = 0;
    int caught = 0;
    TRY
        RAISE(Arena_Failed);
    EXCEPT(Arena_Failed)
        caught = 1;
    FINALLY
        finalized = 1;
    END_TRY;

    ASSERT_EQ(caught, 1, "Arena_Failed should be caught");
    ASSERT_EQ(finalized, 1, "FINALLY should run");
    PASS();
}

/* ===== Mesh Types Tests ===== */

static void test_mesh_types(void)
{
    TEST(mesh_types_valid);

    ComponentMesh cm;
    float verts[9] = {0,0,0, 1,0,0, 0,1,0};
    int32_t faces[3] = {0, 1, 2};
    cm.verts = verts;
    cm.faces = faces;
    cm.nv = 3;
    cm.nf = 1;
    cm.comp_id = 1;
    cm.self = &cm;

    assert(ComponentMesh_valid(&cm));

    /* Invalid: wrong self */
    ComponentMesh bad = cm;
    bad.self = NULL;
    assert(!ComponentMesh_valid(&bad));

    PASS();
}

/* ===== Union-Find Tests ===== */

static void test_union_find(void)
{
    TEST(union_find);

    Arena_T arena = Arena_new();

    UnionFind uf = UF_new(arena, 10);
    ASSERT_EQ(uf.count, 10, "initial count should be 10");

    uf_union(&uf, 0, 1);
    uf_union(&uf, 2, 3);
    uf_union(&uf, 0, 2);

    /* 0,1,2,3 should be in same component */
    ASSERT_EQ(uf_find(&uf, 0), uf_find(&uf, 1), "0 and 1 same comp");
    ASSERT_EQ(uf_find(&uf, 0), uf_find(&uf, 2), "0 and 2 same comp");
    ASSERT_EQ(uf_find(&uf, 0), uf_find(&uf, 3), "0 and 3 same comp");

    /* 4 should be separate */
    assert(uf_find(&uf, 4) != uf_find(&uf, 0));

    ASSERT_EQ(uf.count, 7, "count should be 7");

    Arena_dispose(&arena);
    PASS();
}

/* ===== PCA Tests ===== */

static void test_pca_normal(void)
{
    TEST(pca_normal);

    /* Points in the XY plane (z=0) — normal should be along Z */
    float points[] = {
        0,0,0,  1,0,0,  0,1,0,  1,1,0,
        0.5f,0.5f,0,  0.2f,0.8f,0,  0.8f,0.2f,0,
        0.3f,0.7f,0,  0.7f,0.3f,0,  0.4f,0.6f,0
    };
    float normal[3] = {0};
    float centroid[3] = {0};

    int ret = PCA_normal(points, 10, normal, centroid);
    ASSERT_EQ(ret, 0, "PCA should succeed");

    /* Points have coords (z,y,x). All x values (coord 2) are 0,
     * so normal should be along x-axis: approximately (0,0,1). */
    float abs2 = fabsf(normal[2]);

    /* Smallest eigenvalue direction: along x (dim 2) since all x=0 */
    assert(abs2 > 0.9f);

    /* Sign convention: largest-magnitude component positive */
    assert(normal[2] > 0.0f);

    /* Check centroid */
    assert(fabsf(centroid[1]) < 1.0f);
    assert(fabsf(centroid[2]) < 1.0f);

    PASS();
}

static void test_pca_project(void)
{
    TEST(pca_project);

    float points[] = {1,0,0, 2,0,0, 3,0,0};
    float dir[3] = {1, 0, 0};
    float proj[3] = {0};
    float mn = 0, mx = 0;

    PCA_project(points, 3, dir, proj, &mn, &mx);

    ASSERT_NEAR(proj[0], 1.0f, 1e-5f, "proj[0]");
    ASSERT_NEAR(proj[1], 2.0f, 1e-5f, "proj[1]");
    ASSERT_NEAR(proj[2], 3.0f, 1e-5f, "proj[2]");
    ASSERT_NEAR(mn, 1.0f, 1e-5f, "min");
    ASSERT_NEAR(mx, 3.0f, 1e-5f, "max");

    PASS();
}

static void test_pca_orthonormal_basis(void)
{
    TEST(pca_orthonormal_basis);

    float normal[3] = {0, 0, 1};
    float u[3] = {0}, v[3] = {0};

    PCA_orthonormal_basis(normal, u, v);

    /* u and v should be perpendicular to normal */
    float dot_un = u[0]*normal[0] + u[1]*normal[1] + u[2]*normal[2];
    float dot_vn = v[0]*normal[0] + v[1]*normal[1] + v[2]*normal[2];
    float dot_uv = u[0]*v[0] + u[1]*v[1] + u[2]*v[2];

    ASSERT_NEAR(dot_un, 0.0f, 1e-5f, "u.n should be 0");
    ASSERT_NEAR(dot_vn, 0.0f, 1e-5f, "v.n should be 0");
    ASSERT_NEAR(dot_uv, 0.0f, 1e-5f, "u.v should be 0");

    /* u and v should be unit vectors */
    float len_u = sqrtf(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
    float len_v = sqrtf(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    ASSERT_NEAR(len_u, 1.0f, 1e-5f, "|u| should be 1");
    ASSERT_NEAR(len_v, 1.0f, 1e-5f, "|v| should be 1");

    PASS();
}

/* ===== CSR Tests ===== */

static void test_csr_from_faces(void)
{
    TEST(csr_from_faces);

    Arena_T arena = Arena_new();

    /* Triangle: vertices 0,1,2 */
    int32_t faces[] = {0, 1, 2};
    CSR_T csr = CSR_from_faces(arena, faces, 1, 3);

    ASSERT_EQ(CSR_nrows(csr), 3, "nrows should be 3");

    /* Each vertex should have 2 neighbors */
    const int32_t *off = CSR_offset(csr);
    ASSERT_EQ(off[1] - off[0], 2, "vertex 0 should have 2 neighbors");
    ASSERT_EQ(off[2] - off[1], 2, "vertex 1 should have 2 neighbors");
    ASSERT_EQ(off[3] - off[2], 2, "vertex 2 should have 2 neighbors");

    CSR_validate(csr);

    Arena_dispose(&arena);
    PASS();
}

static void test_csr_matvec(void)
{
    TEST(csr_matvec);

    Arena_T arena = Arena_new();

    /* Two triangles sharing an edge: (0,1,2) and (0,2,3) */
    int32_t faces[] = {0,1,2, 0,2,3};
    CSR_T csr = CSR_from_faces(arena, faces, 2, 4);

    /* Test matvec with unit weights (adjacency sum) */
    /* Vertex 0 has neighbors: 1, 2, 3 */
    /* Vertex 1 has neighbors: 0, 2 */
    /* Vertex 2 has neighbors: 0, 1, 3 */
    /* Vertex 3 has neighbors: 0, 2 */
    float x[4] = {1.0f, 2.0f, 3.0f, 4.0f};
    float y[4] = {0};
    CSR_matvec(csr, x, y);

    /* y[0] = x[1] + x[2] + x[3] = 2+3+4 = 9 */
    ASSERT_NEAR(y[0], 9.0f, 1e-5f, "y[0]");
    /* y[1] = x[0] + x[2] = 1+3 = 4 */
    ASSERT_NEAR(y[1], 4.0f, 1e-5f, "y[1]");

    CSR_validate(csr);

    Arena_dispose(&arena);
    PASS();
}

static void test_csr_empty(void)
{
    TEST(csr_empty);

    Arena_T arena = Arena_new();

    CSR_T csr = CSR_from_faces(arena, NULL, 0, 0);
    ASSERT_EQ(CSR_nrows(csr), 0, "nrows should be 0");
    ASSERT_EQ(CSR_nnz(csr), 0, "nnz should be 0");

    CSR_validate(csr);

    Arena_dispose(&arena);
    PASS();
}

/* ===== KD-Tree Tests ===== */

static void test_kdtree_nearest(void)
{
    TEST(kdtree_nearest);

    Arena_T arena = Arena_new();

    float points[] = {
        0,0,0,   1,0,0,   0,1,0,   0,0,1,
        1,1,0,   1,0,1,   0,1,1,   1,1,1
    };
    KDTree_T tree = KDTree_new(arena, points, 8);

    /* Query near (1,1,1) — should find index 7 */
    float query[3] = {0.9f, 0.9f, 0.9f};
    float dist_sq = 0;
    size_t idx = KDTree_nearest(tree, query, &dist_sq);
    ASSERT_EQ((int)idx, 7, "nearest to (0.9,0.9,0.9) should be vertex 7");

    /* Query exactly at point 0 */
    float query2[3] = {0, 0, 0};
    idx = KDTree_nearest(tree, query2, &dist_sq);
    ASSERT_EQ((int)idx, 0, "nearest to origin should be vertex 0");
    ASSERT_NEAR(dist_sq, 0.0f, 1e-10f, "dist should be 0");

    Arena_dispose(&arena);
    PASS();
}

static void test_kdtree_ball_query(void)
{
    TEST(kdtree_ball_query);

    Arena_T arena = Arena_new();

    /* Points along a line: 0,1,2,...,9 */
    float points[30];
    for (int i = 0; i < 10; i++) {
        points[i * 3 + 0] = (float)i;
        points[i * 3 + 1] = 0.0f;
        points[i * 3 + 2] = 0.0f;
    }

    KDTree_T tree = KDTree_new(arena, points, 10);

    /* Ball query centered at (5,0,0) with radius 1.5 — should find 4,5,6 */
    float center[3] = {5.0f, 0.0f, 0.0f};
    float radius_sq = 1.5f * 1.5f; /* 2.25 */
    int32_t results[10];
    size_t count = KDTree_ball_query(tree, center, radius_sq, results, 10);

    ASSERT_EQ((int)count, 3, "should find 3 points");

    /* Verify the found indices are 4, 5, 6 (in any order) */
    int found[10] = {0};
    for (size_t i = 0; i < count; i++) {
        found[results[i]] = 1;
    }
    assert(found[4] && found[5] && found[6]);

    Arena_dispose(&arena);
    PASS();
}

/* ===== BFS Tests ===== */

static void test_bfs_kring(void)
{
    TEST(bfs_kring);

    Arena_T arena = Arena_new();

    /* Build a path graph: 0-1-2-3-4 */
    /* Faces: (0,1,2), (2,3,4) — creates edges 0-1, 0-2, 1-2, 2-3, 2-4, 3-4 */
    int32_t faces[] = {0,1,2, 2,3,4};
    CSR_T csr = CSR_from_faces(arena, faces, 2, 5);

    int32_t out[5];
    uint32_t dist[5];

    /* 1-ring from vertex 0 */
    size_t count = BFS_kring(csr, 0, 1, out, 5, dist);

    /* Should find vertex 0 and its direct neighbors (1, 2) */
    assert(count >= 2);
    assert(dist[0] == 0);
    assert(dist[1] <= 1);
    assert(dist[2] <= 1);

    Arena_dispose(&arena);
    PASS();
}

/* ===== TIFF I/O Tests ===== */

static void test_tiff_round_trip(void)
{
    TEST(tiff_round_trip);

    Arena_T arena = Arena_new();

    /* Create a small test volume */
    int D = 4, H = 8, W = 8;
    size_t vol_size = (size_t)D * (size_t)H * (size_t)W;
    uint8_t *vol = (uint8_t *)ARENA_CALLOC(arena, (long)vol_size, 1L);

    /* Write a pattern */
    for (int z = 0; z < D; z++) {
        for (int y = 0; y < H; y++) {
            for (int x = 0; x < W; x++) {
                vol[(size_t)z * (size_t)H * (size_t)W +
                    (size_t)y * (size_t)W + (size_t)x] =
                    (uint8_t)((z + y + x) & 0xFF);
            }
        }
    }

    /* Save */
    const char *path = "/tmp/test_vesuvius_tiff.tif";
    int ret = TiffIO_save(path, vol, D, H, W);
    ASSERT_EQ(ret, 0, "save should succeed");

    /* Load back */
    uint8_t *loaded = NULL;
    int ld = 0, lh = 0, lw = 0;
    ret = TiffIO_load(arena, path, &loaded, &ld, &lh, &lw);
    ASSERT_EQ(ret, 0, "load should succeed");
    ASSERT_EQ(ld, D, "D mismatch");
    ASSERT_EQ(lh, H, "H mismatch");
    ASSERT_EQ(lw, W, "W mismatch");

    /* Compare */
    for (size_t i = 0; i < vol_size; i++) {
        if (vol[i] != loaded[i]) {
            FAIL("volume data mismatch");
            Arena_dispose(&arena);
            return;
        }
    }

    /* Cleanup temp file */
    remove(path);

    Arena_dispose(&arena);
    PASS();
}

/* ===== HaloLoader Tests ===== */

/*
 * Build a small synthetic 2-cube grid in a temp dir: self cube fills with
 * 0x11, +Y neighbor fills with 0x22. Load with halo=1 and verify:
 *   - central 2^3 region of padded buffer is all 0x11
 *   - the +Y shell (py == 3, the row beyond self's Y=1 plane) is all 0x22
 *   - everything else is zero (no neighbors loaded for other 25 offsets)
 *   - origin reflects -halo offset
 */
static void test_halo_loader_basic(void)
{
    TEST(halo_loader_basic);

    Arena_T arena = Arena_new();
    int cube_size = 2;
    int halo = 1;
    int p = cube_size + 2 * halo;  /* 4 */

    /* Pick a temp dir; both Linux and Windows define TMP/TEMP. */
    const char *tmpdir = getenv("TMP");
    if (!tmpdir) tmpdir = getenv("TEMP");
    if (!tmpdir) tmpdir = "/tmp";

    char dir[512];
    snprintf(dir, sizeof(dir), "%s/vesuvius_halo_test", tmpdir);
    /* mkdir; ignore failure if it already exists */
    (void)ves_mkdir(dir);

    /* Self cube at origin (0, 0, 0) filled with 0x11. */
    uint8_t self_vol[2 * 2 * 2];
    memset(self_vol, 0x11, sizeof(self_vol));
    char self_path[1024];
    snprintf(self_path, sizeof(self_path), "%s/z00000_y00000_x00000.tif", dir);
    int rc = TiffIO_save(self_path, self_vol, cube_size, cube_size, cube_size);
    ASSERT_EQ(rc, 0, "save self TIFF");

    /* +Y neighbor at (0, 2, 0) filled with 0x22. */
    uint8_t nbr_vol[2 * 2 * 2];
    memset(nbr_vol, 0x22, sizeof(nbr_vol));
    char nbr_path[1024];
    snprintf(nbr_path, sizeof(nbr_path), "%s/z00000_y00002_x00000.tif", dir);
    rc = TiffIO_save(nbr_path, nbr_vol, cube_size, cube_size, cube_size);
    ASSERT_EQ(rc, 0, "save +Y neighbor TIFF");

    /* Load with halo=1. */
    uint8_t *padded = NULL;
    int p_out = 0;
    int64_t origin[3] = {0, 0, 0};
    rc = HaloLoader_load(arena, dir, "z00000_y00000_x00000",
                         cube_size, halo, &padded, &p_out, origin);
    ASSERT_EQ(rc, 0, "HaloLoader_load");
    ASSERT_EQ(p_out, p, "p_size");
    ASSERT_EQ((int)origin[0], -halo, "origin z");
    ASSERT_EQ((int)origin[1], -halo, "origin y");
    ASSERT_EQ((int)origin[2], -halo, "origin x");

    /* Center 2^3 = padded[halo..halo+2)^3 should be all 0x11. */
    for (int z = halo; z < halo + cube_size; z++) {
        for (int y = halo; y < halo + cube_size; y++) {
            for (int x = halo; x < halo + cube_size; x++) {
                uint8_t v = padded[(size_t)z * (size_t)p * (size_t)p + (size_t)y * (size_t)p + (size_t)x];
                if (v != 0x11) {
                    FAIL("center region should be 0x11");
                    Arena_dispose(&arena);
                    return;
                }
            }
        }
    }

    /* +Y shell at py = halo + cube_size = 3 (only existing row in [3, 4)),
     * for pz in [halo, halo+2) = [1, 3), px in [halo, halo+2) = [1, 3),
     * should be 0x22 (from the +Y neighbor's first Y-row). */
    for (int z = halo; z < halo + cube_size; z++) {
        for (int x = halo; x < halo + cube_size; x++) {
            int py = halo + cube_size;  /* = 3 */
            uint8_t v = padded[(size_t)z * (size_t)p * (size_t)p + (size_t)py * (size_t)p + (size_t)x];
            if (v != 0x22) {
                FAIL("+Y shell should be 0x22");
                Arena_dispose(&arena);
                return;
            }
        }
    }

    /* -Y shell (py=0) should be all zeros (no -Y neighbor file present). */
    for (int z = 0; z < p; z++) {
        for (int x = 0; x < p; x++) {
            uint8_t v = padded[(size_t)z * (size_t)p * (size_t)p + (size_t)0 * (size_t)p + (size_t)x];
            if (v != 0) {
                FAIL("-Y shell should be 0 (no neighbor)");
                Arena_dispose(&arena);
                return;
            }
        }
    }

    remove(self_path);
    remove(nbr_path);
    Arena_dispose(&arena);
    PASS();
}

/*
 * Halo=0 degenerates to a TiffIO_load of the self cube only; padded volume
 * == cube_size^3 and equals the input.
 */
static void test_halo_loader_zero_halo(void)
{
    TEST(halo_loader_zero_halo);

    Arena_T arena = Arena_new();
    int cube_size = 3;

    const char *tmpdir = getenv("TMP");
    if (!tmpdir) tmpdir = getenv("TEMP");
    if (!tmpdir) tmpdir = "/tmp";

    char dir[512];
    snprintf(dir, sizeof(dir), "%s/vesuvius_halo_test", tmpdir);
    (void)ves_mkdir(dir);

    uint8_t vol[3 * 3 * 3];
    for (int i = 0; i < 27; i++) vol[i] = (uint8_t)(i + 1);
    char path[1024];
    snprintf(path, sizeof(path), "%s/z00128_y00256_x00384.tif", dir);
    int rc = TiffIO_save(path, vol, cube_size, cube_size, cube_size);
    ASSERT_EQ(rc, 0, "save TIFF");

    uint8_t *padded = NULL;
    int p_out = 0;
    int64_t origin[3] = {0, 0, 0};
    rc = HaloLoader_load(arena, dir, "z00128_y00256_x00384",
                         cube_size, 0, &padded, &p_out, origin);
    ASSERT_EQ(rc, 0, "HaloLoader_load");
    ASSERT_EQ(p_out, cube_size, "p_size = cube_size");
    ASSERT_EQ((int)origin[0], 128, "origin z");
    ASSERT_EQ((int)origin[1], 256, "origin y");
    ASSERT_EQ((int)origin[2], 384, "origin x");
    for (int i = 0; i < 27; i++) {
        if (padded[i] != (uint8_t)(i + 1)) {
            FAIL("padded value mismatch");
            Arena_dispose(&arena);
            return;
        }
    }

    remove(path);
    Arena_dispose(&arena);
    PASS();
}

/* ===== OBJ I/O Tests ===== */

static void test_obj_round_trip(void)
{
    TEST(obj_round_trip);

    Arena_T arena = Arena_new();

    float verts[] = {0,0,0, 1,0,0, 0,1,0, 1,1,0};
    int32_t faces[] = {0,1,2, 1,3,2};
    const char *path = "/tmp/test_vesuvius.obj";

    int ret = ObjIO_write(path, verts, 4, faces, 2);
    ASSERT_EQ(ret, 0, "write should succeed");

    float *out_verts = NULL;
    int32_t *out_faces = NULL;
    size_t nv = 0, nf = 0;

    ret = ObjIO_read(arena, path, &out_verts, &nv, &out_faces, &nf);
    ASSERT_EQ(ret, 0, "read should succeed");
    ASSERT_EQ((int)nv, 4, "vertex count");
    ASSERT_EQ((int)nf, 2, "face count");

    ASSERT_NEAR(out_verts[0], 0.0f, 1e-4f, "v0.z");
    ASSERT_NEAR(out_verts[3], 1.0f, 1e-4f, "v1.z");

    ASSERT_EQ(out_faces[0], 0, "f0[0]");
    ASSERT_EQ(out_faces[1], 1, "f0[1]");
    ASSERT_EQ(out_faces[2], 2, "f0[2]");

    remove(path);
    Arena_dispose(&arena);
    PASS();
}

/* ===== Pipeline Constants Test ===== */

static void test_pipeline_constants(void)
{
    TEST(pipeline_constants);

    /* Verify key constants are defined and have expected values */
    ASSERT_EQ(MIN_CC_SIZE, 500, "MIN_CC_SIZE");
    ASSERT_EQ(MAX_COMPONENTS, 20, "MAX_COMPONENTS");
    ASSERT_EQ(SEED_RING, 20, "SEED_RING");
    ASSERT_EQ(MAX_FLOW_LIMIT, 100, "MAX_FLOW_LIMIT");
    ASSERT_EQ(BRIDGE_MAX_DEPTH, 10, "BRIDGE_MAX_DEPTH");
    ASSERT_EQ(PARALLEL_VERTEX_CUTOFF, 1000, "PARALLEL_VERTEX_CUTOFF");

    PASS();
}

/* ===== Main ===== */

#ifdef TEST_HARNESS
int test_common_main(void)
#else
int main(void)
#endif
{
    printf("=== test_common ===\n\n");

    printf("[Arena]\n");
    test_arena_basic();
    test_arena_save_restore();
    test_arena_large_alloc();
    test_arena_free_reuse();

    printf("\n[Except]\n");
    test_except_try_catch();
    test_except_finally();

    printf("\n[MeshTypes]\n");
    test_mesh_types();

    printf("\n[UnionFind]\n");
    test_union_find();

    printf("\n[PCA]\n");
    test_pca_normal();
    test_pca_project();
    test_pca_orthonormal_basis();

    printf("\n[CSR]\n");
    test_csr_from_faces();
    test_csr_matvec();
    test_csr_empty();

    printf("\n[KDTree]\n");
    test_kdtree_nearest();
    test_kdtree_ball_query();

    printf("\n[BFS]\n");
    test_bfs_kring();

    printf("\n[TIFF I/O]\n");
    test_tiff_round_trip();

    printf("\n[HaloLoader]\n");
    test_halo_loader_basic();
    test_halo_loader_zero_halo();

    printf("\n[OBJ I/O]\n");
    test_obj_round_trip();

    printf("\n[Constants]\n");
    test_pipeline_constants();

    printf("\n=== Results: %d/%d passed ===\n", tests_passed, tests_run);

    return (tests_passed == tests_run) ? 0 : 1;
}
