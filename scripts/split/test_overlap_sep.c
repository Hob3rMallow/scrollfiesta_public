/*
 * test_overlap_sep.c -- Unit + integration tests for the overlap separator.
 *
 * Tests:
 *   1. 2D triangle overlap detection (known overlapping/non-overlapping pairs)
 *   2. Full OverlapSep_process on real OBJ mesh with debug OBJ dumps
 *   3. Connectivity check: each output piece should be a single CC
 */
#include "split/overlap_sep.h"
#include "split/oracle.h"
#include "common/arena.h"
#include "common/obj_io.h"
#include "common/pca.h"
#include "common/mesh_types.h"
#include "common/union_find.h"
#include "common/csr.h"
#include "common/pipeline_constants.h"

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef _WIN32
#include <direct.h>
#define MKDIR(d) _mkdir(d)
#else
#include <sys/stat.h>
#define MKDIR(d) mkdir(d, 0755)
#endif

static int g_pass = 0, g_fail = 0;

#define CHECK(cond, msg) do { \
    if (cond) { g_pass++; printf("    PASS: %s\n", msg); } \
    else { g_fail++; printf("    FAIL: %s\n", msg); } \
} while(0)

/* ------------------------------------------------------------------ */
/* Half-edge connected components via Union-Find on faces              */
/* ------------------------------------------------------------------ */

typedef struct {
    int32_t v0, v1;
    int32_t face_idx;
} HEdge;

static int cmp_hedge(const void *a, const void *b)
{
    const HEdge *ea = (const HEdge *)a;
    const HEdge *eb = (const HEdge *)b;
    if (ea->v0 != eb->v0) return (ea->v0 < eb->v0) ? -1 : 1;
    if (ea->v1 != eb->v1) return (ea->v1 < eb->v1) ? -1 : 1;
    return 0;
}

static int count_face_components(Arena_T arena,
                                 const int32_t *faces, size_t nf)
{
    if (nf == 0) return 0;

    Arena_Mark mark = Arena_save(arena);

    size_t n_he = nf * 3;
    HEdge *edges = (HEdge *)ARENA_ALLOC(arena,
                     (long)(n_he * sizeof(HEdge)));

    for (size_t f = 0; f < nf; f++) {
        for (int e = 0; e < 3; e++) {
            int32_t a = faces[f * 3 + e];
            int32_t b = faces[f * 3 + ((e + 1) % 3)];
            size_t idx = f * 3 + (size_t)e;
            edges[idx].v0 = (a < b) ? a : b;
            edges[idx].v1 = (a < b) ? b : a;
            edges[idx].face_idx = (int32_t)f;
        }
    }

    qsort(edges, n_he, sizeof(HEdge), cmp_hedge);

    UnionFind uf = UF_new(arena, (int32_t)nf);

    size_t i = 0;
    while (i < n_he) {
        size_t j = i + 1;
        while (j < n_he &&
               edges[j].v0 == edges[i].v0 &&
               edges[j].v1 == edges[i].v1) {
            uf_union(&uf, edges[i].face_idx, edges[j].face_idx);
            j++;
        }
        i = j;
    }

    int result = uf.count;
    Arena_restore(arena, mark);
    return result;
}

static int32_t *label_face_components(Arena_T arena,
                                      const int32_t *faces, size_t nf,
                                      size_t **out_sizes,
                                      size_t *out_count)
{
    if (nf == 0 || nf > (size_t)INT32_MAX) return NULL;
    size_t n_he = nf * 3;
    HEdge *edges = (HEdge *)ARENA_ALLOC(
        arena, (long)(n_he * sizeof(HEdge)));
    for (size_t f = 0; f < nf; f++) {
        for (int e = 0; e < 3; e++) {
            int32_t a = faces[f * 3 + (size_t)e];
            int32_t b = faces[f * 3 + (size_t)((e + 1) % 3)];
            HEdge *h = &edges[f * 3 + (size_t)e];
            h->v0 = a < b ? a : b;
            h->v1 = a < b ? b : a;
            h->face_idx = (int32_t)f;
        }
    }
    qsort(edges, n_he, sizeof(HEdge), cmp_hedge);
    UnionFind uf = UF_new(arena, (int32_t)nf);
    for (size_t i = 0; i < n_he;) {
        size_t j = i + 1;
        while (j < n_he && edges[j].v0 == edges[i].v0 &&
               edges[j].v1 == edges[i].v1) {
            uf_union(&uf, edges[i].face_idx, edges[j].face_idx);
            j++;
        }
        i = j;
    }
    int32_t *root_label = (int32_t *)ARENA_ALLOC(
        arena, (long)(nf * sizeof(int32_t)));
    int32_t *label = (int32_t *)ARENA_ALLOC(
        arena, (long)(nf * sizeof(int32_t)));
    size_t *size = (size_t *)ARENA_CALLOC(
        arena, (long)nf, (long)sizeof(size_t));
    for (size_t i = 0; i < nf; i++) root_label[i] = -1;
    size_t ncomponent = 0;
    for (size_t f = 0; f < nf; f++) {
        int32_t root = uf_find(&uf, (int32_t)f);
        if (root_label[root] < 0)
            root_label[root] = (int32_t)ncomponent++;
        label[f] = root_label[root];
        size[label[f]]++;
    }
    *out_sizes = size;
    *out_count = ncomponent;
    return label;
}

typedef struct {
    int32_t component0, component1;
    size_t pairs;
} AuditComponentPair;

static int compare_u64_key(const void *pa, const void *pb)
{
    uint64_t a = *(const uint64_t *)pa;
    uint64_t b = *(const uint64_t *)pb;
    return a < b ? -1 : (a > b ? 1 : 0);
}

static int compare_audit_pair_count(const void *pa, const void *pb)
{
    const AuditComponentPair *a = (const AuditComponentPair *)pa;
    const AuditComponentPair *b = (const AuditComponentPair *)pb;
    if (a->pairs != b->pairs) return a->pairs > b->pairs ? -1 : 1;
    if (a->component0 != b->component0)
        return a->component0 < b->component0 ? -1 : 1;
    return a->component1 < b->component1 ? -1 :
           (a->component1 > b->component1 ? 1 : 0);
}

static int run_overlap_audit(const char *obj_path, const char *csv_path)
{
    Arena_T arena = Arena_new();
    float *verts = NULL;
    int32_t *faces = NULL;
    size_t nv = 0, nf = 0;
    if (ObjIO_read(arena, obj_path, &verts, &nv, &faces, &nf) != 0) {
        fprintf(stderr, "overlap audit: failed to load %s\n", obj_path);
        Arena_dispose(&arena);
        return 1;
    }
    ComponentMesh mesh;
    memset(&mesh, 0, sizeof mesh);
    mesh.verts = verts;
    mesh.faces = faces;
    mesh.nv = nv;
    mesh.nf = nf;
    PCA_normal(verts, nv, mesh.pca_normal, mesh.centroid);
    mesh.self = &mesh;

    OverlapSepPairSet overlap;
    if (OverlapSep_detect_pairs(arena, &mesh, &overlap) != 0) {
        fprintf(stderr, "overlap audit: detector failed\n");
        Arena_dispose(&arena);
        return 1;
    }
    size_t *component_size = NULL, ncomponent = 0;
    int32_t *component = label_face_components(
        arena, faces, nf, &component_size, &ncomponent);
    if (component == NULL || ncomponent > UINT32_MAX) {
        Arena_dispose(&arena);
        return 1;
    }

    uint64_t *key = (uint64_t *)ARENA_ALLOC(
        arena, (long)((overlap.count ? overlap.count : 1) * sizeof(uint64_t)));
    size_t same_pairs = 0;
    for (size_t i = 0; i < overlap.count; i++) {
        int32_t a = component[overlap.face0[i]];
        int32_t b = component[overlap.face1[i]];
        if (a > b) { int32_t t = a; a = b; b = t; }
        if (a == b) same_pairs++;
        key[i] = ((uint64_t)(uint32_t)a << 32) | (uint32_t)b;
    }
    qsort(key, overlap.count, sizeof(uint64_t), compare_u64_key);
    AuditComponentPair *bucket = (AuditComponentPair *)ARENA_ALLOC(
        arena, (long)((overlap.count ? overlap.count : 1) *
                      sizeof(AuditComponentPair)));
    size_t nbucket = 0, self_components = 0;
    for (size_t first = 0; first < overlap.count;) {
        size_t last = first + 1;
        while (last < overlap.count && key[last] == key[first]) last++;
        bucket[nbucket].component0 = (int32_t)(key[first] >> 32);
        bucket[nbucket].component1 = (int32_t)(key[first] & UINT32_MAX);
        bucket[nbucket].pairs = last - first;
        if (bucket[nbucket].component0 == bucket[nbucket].component1)
            self_components++;
        nbucket++;
        first = last;
    }
    qsort(bucket, nbucket, sizeof(AuditComponentPair),
          compare_audit_pair_count);

    printf("Overlap pair audit: %s\n", obj_path);
    printf("  mesh: %zu vertices, %zu faces, %zu face-components\n",
           nv, nf, ncomponent);
    printf("  exact overlap pairs: %zu\n", overlap.count);
    printf("  same-component: %zu (%.2f%%) across %zu components\n",
           same_pairs, overlap.count ? 100.0 * (double)same_pairs /
                                       (double)overlap.count : 0.0,
           self_components);
    printf("  cross-component: %zu (%.2f%%) across %zu component pairs\n",
           overlap.count - same_pairs,
           overlap.count ? 100.0 * (double)(overlap.count - same_pairs) /
                           (double)overlap.count : 0.0,
           nbucket - self_components);
    printf("  top component pairs (c0,c1,pairs,faces0,faces1):\n");
    size_t top = nbucket < 24 ? nbucket : 24;
    for (size_t i = 0; i < top; i++) {
        const AuditComponentPair *b = &bucket[i];
        printf("    %d,%d,%zu,%zu,%zu%s\n", b->component0,
               b->component1, b->pairs, component_size[b->component0],
               component_size[b->component1],
               b->component0 == b->component1 ? ",SELF" : "");
    }

    if (csv_path != NULL) {
        FILE *fp = fopen(csv_path, "wb");
        if (fp == NULL) {
            fprintf(stderr, "overlap audit: cannot write %s\n", csv_path);
            Arena_dispose(&arena);
            return 1;
        }
        fprintf(fp, "component0,component1,overlap_pairs,component0_faces,"
                    "component1_faces,same_component\n");
        for (size_t i = 0; i < nbucket; i++) {
            const AuditComponentPair *b = &bucket[i];
            fprintf(fp, "%d,%d,%zu,%zu,%zu,%d\n", b->component0,
                    b->component1, b->pairs,
                    component_size[b->component0],
                    component_size[b->component1],
                    b->component0 == b->component1);
        }
        fclose(fp);
        printf("  wrote %s\n", csv_path);
    }
    Arena_dispose(&arena);
    return 0;
}

/* ------------------------------------------------------------------ */
/* Test 1: Synthetic mesh overlap detection                            */
/* ------------------------------------------------------------------ */

static void test_geometry_basics(const char *out_dir)
{
    printf("  Test 1: Synthetic mesh overlap detection\n");

    Arena_T arena = Arena_new();

    /* 2-face mesh: too small, should skip */
    float verts[] = {
        0.0f, 0.0f, 0.0f,
        2.0f, 0.0f, 0.0f,
        1.0f, 2.0f, 0.0f,
        0.5f, 0.5f, 0.0f,
        2.5f, 0.5f, 0.0f,
        1.5f, 2.5f, 0.0f,
    };
    int32_t faces[] = { 0, 1, 2, 3, 4, 5 };

    ComponentMesh cm;
    memset(&cm, 0, sizeof(cm));
    cm.verts = verts;
    cm.faces = faces;
    cm.nv = 6;
    cm.nf = 2;
    cm.pca_normal[0] = 0.0f;
    cm.pca_normal[1] = 0.0f;
    cm.pca_normal[2] = 1.0f;
    cm.centroid[0] = 1.0f;
    cm.centroid[1] = 1.0f;
    cm.centroid[2] = 0.0f;
    cm.self = &cm;

    ComponentMesh *out = NULL;
    size_t out_count = 0;
    int rc = OverlapSep_process(arena, &cm, 2, 1, 10.0, &out, &out_count);

    CHECK(rc == -1, "2-face mesh returns -1 (too small)");
    CHECK(out_count == 1, "2-face mesh returns 1 unchanged copy");

    /* Test 1b: Two overlapping 10x10 grids (324 faces total) */
    int grid_n = 10;
    size_t grid_nv = (size_t)(grid_n * grid_n);
    size_t grid_nf = (size_t)((grid_n - 1) * (grid_n - 1) * 2);
    size_t total_nv = grid_nv * 2;
    size_t total_nf = grid_nf * 2;

    float *big_verts = (float *)ARENA_ALLOC(arena,
        (long)(total_nv * 3 * sizeof(float)));
    int32_t *big_faces = (int32_t *)ARENA_ALLOC(arena,
        (long)(total_nf * 3 * sizeof(int32_t)));

    for (int iy = 0; iy < grid_n; iy++) {
        for (int ix = 0; ix < grid_n; ix++) {
            size_t vi = (size_t)(iy * grid_n + ix);
            big_verts[vi * 3 + 0] = (float)ix;
            big_verts[vi * 3 + 1] = (float)iy;
            big_verts[vi * 3 + 2] = 0.0f;
        }
    }
    for (int iy = 0; iy < grid_n; iy++) {
        for (int ix = 0; ix < grid_n; ix++) {
            size_t vi = grid_nv + (size_t)(iy * grid_n + ix);
            big_verts[vi * 3 + 0] = (float)ix + 0.5f;
            big_verts[vi * 3 + 1] = (float)iy + 0.5f;
            big_verts[vi * 3 + 2] = 0.0f;
        }
    }

    size_t fi = 0;
    for (int base = 0; base < 2; base++) {
        int32_t off = (int32_t)(base * (int32_t)grid_nv);
        for (int iy = 0; iy < grid_n - 1; iy++) {
            for (int ix = 0; ix < grid_n - 1; ix++) {
                int32_t v00 = off + (int32_t)(iy * grid_n + ix);
                int32_t v10 = v00 + 1;
                int32_t v01 = v00 + (int32_t)grid_n;
                int32_t v11 = v01 + 1;
                big_faces[fi * 3 + 0] = v00;
                big_faces[fi * 3 + 1] = v10;
                big_faces[fi * 3 + 2] = v01;
                fi++;
                big_faces[fi * 3 + 0] = v10;
                big_faces[fi * 3 + 1] = v11;
                big_faces[fi * 3 + 2] = v01;
                fi++;
            }
        }
    }

    ComponentMesh cm2;
    memset(&cm2, 0, sizeof(cm2));
    cm2.verts = big_verts;
    cm2.faces = big_faces;
    cm2.nv = total_nv;
    cm2.nf = total_nf;
    cm2.pca_normal[0] = 0.0f;
    cm2.pca_normal[1] = 0.0f;
    cm2.pca_normal[2] = 1.0f;
    cm2.centroid[0] = 5.0f;
    cm2.centroid[1] = 5.0f;
    cm2.centroid[2] = 0.0f;
    cm2.self = &cm2;

    printf("    synthetic grid: nv=%zu nf=%zu\n", total_nv, total_nf);

    /* Dump input synthetic mesh */
    if (out_dir) {
        char path[1024];
        snprintf(path, sizeof(path), "%s/synth_input.obj", out_dir);
        ObjIO_write(path, big_verts, total_nv, big_faces, total_nf);
        printf("    wrote %s\n", path);
    }

    /* Enable debug OBJ dumps for synthetic test */
    if (out_dir) {
        char subdir[1024];
        snprintf(subdir, sizeof(subdir), "%s/synth", out_dir);
        MKDIR(subdir);
        OverlapSep_set_debug_dir(subdir);
    }

    ComponentMesh *out2 = NULL;
    size_t out_count2 = 0;
    rc = OverlapSep_process(arena, &cm2, 2, 1, 30.0, &out2, &out_count2);

    OverlapSep_set_debug_dir(NULL);

    printf("    result: rc=%d, pieces=%zu\n", rc, out_count2);
    CHECK(rc == 0, "overlapping grids are separated");
    CHECK(out_count2 >= 2, "at least 2 pieces produced");

    if (out_count2 >= 2) {
        for (size_t p = 0; p < out_count2; p++) {
            printf("    piece %zu: nv=%zu nf=%zu\n",
                   p, out2[p].nv, out2[p].nf);
        }
    }

    Arena_dispose(&arena);
}

/* ------------------------------------------------------------------ */
/* Test 2: Full pipeline on real OBJ mesh                              */
/* ------------------------------------------------------------------ */

static void test_real_mesh(const char *obj_path, const char *out_dir,
                           int expected_sheets, const char *label,
                           double timeout_sec)
{
    printf("  Test: %s (%s, expect >= %d pieces)\n",
           label, obj_path, expected_sheets);

    Arena_T arena = Arena_new();

    float *verts = NULL;
    size_t nv = 0;
    int32_t *faces = NULL;
    size_t nf = 0;

    int rc = ObjIO_read(arena, obj_path, &verts, &nv, &faces, &nf);
    if (rc != 0) {
        printf("    SKIP: failed to load %s\n", obj_path);
        Arena_dispose(&arena);
        return;
    }
    printf("    loaded: %zu verts, %zu faces\n", nv, nf);

    ComponentMesh cm;
    memset(&cm, 0, sizeof(cm));
    cm.verts = verts;
    cm.faces = faces;
    cm.nv = nv;
    cm.nf = nf;
    PCA_normal(verts, nv, cm.pca_normal, cm.centroid);
    cm.self = &cm;

    printf("    PCA normal: [%.3f, %.3f, %.3f]\n",
           (double)cm.pca_normal[0], (double)cm.pca_normal[1],
           (double)cm.pca_normal[2]);

    /* Detector 1: Oracle UV-grid raycast sheet count. */
    int oracle_sheets = Oracle_count_sheets(arena, &cm, ORACLE_UV_GRID_SIZE);
    printf("    DETECT oracle_sheets = %d  (nf=%zu; overlap-pair count is the "
           "'phase 3 (overlaps)' line below, density = pairs/nf)\n",
           oracle_sheets, nf);

    /* Enable debug OBJ dumps */
    if (out_dir) {
        char subdir[1024];
        snprintf(subdir, sizeof(subdir), "%s/%s", out_dir, label);
        MKDIR(subdir);
        OverlapSep_set_debug_dir(subdir);
    }

    ComponentMesh *out_meshes = NULL;
    size_t out_count = 0;
    rc = OverlapSep_process(arena, &cm, expected_sheets, 1, timeout_sec,
                            &out_meshes, &out_count);

    OverlapSep_set_debug_dir(NULL);

    printf("    OverlapSep_process: rc=%d, pieces=%zu\n", rc, out_count);

    {
        char msg[256];
        snprintf(msg, sizeof(msg), "%s: produced >= %d pieces (got %zu)",
                 label, expected_sheets, out_count);
        CHECK(out_count >= (size_t)expected_sheets, msg);
    }

    /* Check connectivity of each piece */
    int conn_ok = 1;
    for (size_t p = 0; p < out_count; p++) {
        int cc = count_face_components(arena, out_meshes[p].faces,
                                       out_meshes[p].nf);
        printf("    piece %zu: nv=%zu nf=%zu CC=%d\n",
               p, out_meshes[p].nv, out_meshes[p].nf, cc);
        if (cc != 1) conn_ok = 0;
    }

    if (out_count >= 2) {
        char msg[256];
        snprintf(msg, sizeof(msg), "%s: all pieces are single CCs", label);
        CHECK(conn_ok, msg);
    }

    Arena_dispose(&arena);
}

/* ------------------------------------------------------------------ */
/* Main                                                                */
/* ------------------------------------------------------------------ */

int main(int argc, char *argv[])
{
    const char *out_dir = "output/overlap_test";

    if (argc > 2 && strcmp(argv[1], "--audit") == 0)
        return run_overlap_audit(argv[2], argc > 3 ? argv[3] : NULL);

    /* CLI: run the current OverlapSep on one OBJ to reproduce/diagnose a split
     * failure -- test_overlap_sep --obj <mesh.obj> [expected_sheets] [timeout_s] */
    if (argc > 2 && strcmp(argv[1], "--obj") == 0) {
        int exp = (argc > 3) ? atoi(argv[3]) : 2;
        double tmo = (argc > 4) ? atof(argv[4]) : 30.0;
        MKDIR("output"); MKDIR(out_dir);
        printf("CLI overlap test: obj=%s expected=%d timeout=%.1fs\n", argv[2], exp, tmo);
        test_real_mesh(argv[2], out_dir, exp, "cli", tmo);
        printf("Results: %d passed, %d failed\n", g_pass, g_fail);
        return g_fail > 0 ? 1 : 0;
    }

    if (argc > 1) out_dir = argv[1];

    /* Create output directories */
    MKDIR("output");
    MKDIR(out_dir);

    printf("Overlap separator tests:\n");
    printf("========================\n");
    printf("  output dir: %s\n\n", out_dir);

    test_geometry_basics(out_dir);
    printf("\n");

    /* Test all overlap test sheets */
    test_real_mesh("test_sheets/overlap_split_1.obj", out_dir,
                   2, "overlap_split_1", 120.0);
    printf("\n");
    test_real_mesh("test_sheets/overlap_split_2.obj", out_dir,
                   3, "overlap_split_2", 120.0);
    printf("\n");

    /* overlap_split_3: leftover from overlap_split_1.
     * Test if it splits on its own into 2 parts.
     * If single-round doesn't split, try a second round on each piece. */
    test_real_mesh("test_sheets/overlap_split_3.obj", out_dir,
                   2, "overlap_split_3", 120.0);

    printf("\n========================\n");
    printf("Results: %d passed, %d failed\n", g_pass, g_fail);
    return g_fail > 0 ? 1 : 0;
}
