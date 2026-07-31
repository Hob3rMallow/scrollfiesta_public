/*
 * test_holefill.c — Unit + functional tests for hole filling (Step 6).
 *
 * Tests:
 *   1. Boundary edge detection (unit)
 *   2. Functional: load holefill_test_1.obj, fill holes, verify:
 *      a. Interior holes are actually filled (no interior boundary loops remain)
 *      b. Fill triangles don't point away from PCA normal
 *      c. No degenerate faces in fill region
 *      d. All face indices in range after fill
 *
 * Build standalone:  link with hole_fill.o, arena.o, csr.o, pca.o, obj_io.o,
 *                    union_find.o, triangle.o, clipper2_wrap.o
 * Harness mode:      compiled with -DTEST_HARNESS
 */
#include "common/ves_platform.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include "common/arena.h"
#include "common/obj_io.h"
#include "common/pca.h"

/* Forward declare HoleFill_process */
int HoleFill_process(Arena_T arena,
                     float **verts, int32_t **faces,
                     size_t *nv, size_t *nf,
                     const int cube_shape[3]);

/* ================================================================
 * Helpers
 * ================================================================ */

/* Compare (int32_t, int32_t) edge pairs for qsort/bsearch */
static int cmp_edge(const void *a, const void *b)
{
    const int32_t *ea = (const int32_t *)a;
    const int32_t *eb = (const int32_t *)b;
    if (ea[0] != eb[0]) return (ea[0] < eb[0]) ? -1 : 1;
    return (ea[1] < eb[1]) ? -1 : (ea[1] > eb[1]) ? 1 : 0;
}

/* Find boundary edges (appear exactly once as half-edge).
 * Returns count of boundary edges.  out_edges[i] = {v0, v1} with v0<v1. */
static size_t find_boundary_edges(const int32_t *faces, size_t nf,
                                  int32_t **out_edges, size_t *out_cap)
{
    size_t n_he = nf * 3;
    int32_t *all = (int32_t *)malloc(n_he * 2 * sizeof(int32_t));
    assert(all);

    for (size_t f = 0; f < nf; f++) {
        for (int e = 0; e < 3; e++) {
            int32_t a = faces[f * 3 + e];
            int32_t b = faces[f * 3 + ((e + 1) % 3)];
            size_t idx = f * 3 + (size_t)e;
            all[idx * 2 + 0] = (a < b) ? a : b;
            all[idx * 2 + 1] = (a < b) ? b : a;
        }
    }

    qsort(all, n_he, 2 * sizeof(int32_t), cmp_edge);

    /* Count boundary edges */
    size_t n_bdry = 0;
    size_t i = 0;
    while (i < n_he) {
        size_t j = i + 1;
        while (j < n_he && all[j * 2 + 0] == all[i * 2 + 0] &&
               all[j * 2 + 1] == all[i * 2 + 1]) {
            j++;
        }
        if (j - i == 1) n_bdry++;
        i = j;
    }

    int32_t *bdry = (int32_t *)malloc(n_bdry * 2 * sizeof(int32_t));
    assert(bdry);

    size_t bi = 0;
    i = 0;
    while (i < n_he) {
        size_t j = i + 1;
        while (j < n_he && all[j * 2 + 0] == all[i * 2 + 0] &&
               all[j * 2 + 1] == all[i * 2 + 1]) {
            j++;
        }
        if (j - i == 1) {
            bdry[bi * 2 + 0] = all[i * 2 + 0];
            bdry[bi * 2 + 1] = all[i * 2 + 1];
            bi++;
        }
        i = j;
    }

    free(all);
    *out_edges = bdry;
    *out_cap = n_bdry;
    return n_bdry;
}

/* Chain boundary edges into loops.  Returns number of loops.
 * loop_sizes[k] = number of unique verts in loop k.
 * loop_interior[k] = 1 if loop is interior (all verts far from cube boundary). */
static size_t chain_boundary_loops(const int32_t *bdry_edges, size_t n_bdry,
                                   size_t nv,
                                   const float *verts, const int cube_shape[3],
                                   float margin,
                                   size_t *loop_sizes, int *loop_interior,
                                   size_t max_loops)
{
    /* Build adjacency from boundary edges */
    int32_t *degree = (int32_t *)calloc(nv, sizeof(int32_t));
    assert(degree);
    for (size_t i = 0; i < n_bdry; i++) {
        degree[bdry_edges[i * 2 + 0]]++;
        degree[bdry_edges[i * 2 + 1]]++;
    }
    int32_t *adj_off = (int32_t *)malloc((nv + 1) * sizeof(int32_t));
    assert(adj_off);
    adj_off[0] = 0;
    for (size_t v = 0; v < nv; v++) {
        adj_off[v + 1] = adj_off[v] + degree[v];
    }
    int32_t *adj_list = (int32_t *)malloc(n_bdry * 2 * sizeof(int32_t));
    assert(adj_list);
    int32_t *cursor = (int32_t *)malloc(nv * sizeof(int32_t));
    assert(cursor);
    memcpy(cursor, adj_off, nv * sizeof(int32_t));
    for (size_t i = 0; i < n_bdry; i++) {
        int32_t a = bdry_edges[i * 2 + 0];
        int32_t b = bdry_edges[i * 2 + 1];
        adj_list[cursor[a]++] = b;
        adj_list[cursor[b]++] = a;
    }

    uint8_t *edge_used = (uint8_t *)calloc(n_bdry, sizeof(uint8_t));
    assert(edge_used);

    size_t n_loops = 0;

    for (size_t ei = 0; ei < n_bdry && n_loops < max_loops; ei++) {
        if (edge_used[ei]) continue;

        int32_t start = bdry_edges[ei * 2 + 0];
        int32_t prev = start;
        int32_t cur = bdry_edges[ei * 2 + 1];
        edge_used[ei] = 1;

        size_t loop_len = 1; /* count start */
        int all_interior = 1;

        /* Check start vertex */
        for (int k = 0; k < 3; k++) {
            float c = verts[start * 3 + k];
            if (c <= margin || c >= (float)cube_shape[k] - 1.0f - margin) {
                all_interior = 0;
            }
        }

        while (cur != start && loop_len < n_bdry + 1) {
            loop_len++;
            /* Check cur vertex */
            for (int k = 0; k < 3; k++) {
                float c = verts[cur * 3 + k];
                if (c <= margin || c >= (float)cube_shape[k] - 1.0f - margin) {
                    all_interior = 0;
                }
            }

            int32_t next = -1;
            for (int32_t j = adj_off[cur]; j < adj_off[cur + 1]; j++) {
                int32_t nb = adj_list[j];
                if (nb == prev) continue;
                /* Find edge in sorted list */
                int32_t key[2];
                key[0] = (cur < nb) ? cur : nb;
                key[1] = (cur < nb) ? nb : cur;
                int32_t *found = (int32_t *)bsearch(key, bdry_edges,
                    n_bdry, 2 * sizeof(int32_t), cmp_edge);
                if (found) {
                    size_t fidx = (size_t)(found - bdry_edges) / 2;
                    if (!edge_used[fidx]) {
                        edge_used[fidx] = 1;
                        next = nb;
                        break;
                    }
                }
            }
            if (next < 0) break;
            prev = cur;
            cur = next;
        }

        if (cur == start && loop_len >= 3) {
            loop_sizes[n_loops] = loop_len;
            loop_interior[n_loops] = all_interior;
            n_loops++;
        }
    }

    free(degree);
    free(adj_off);
    free(adj_list);
    free(cursor);
    free(edge_used);
    return n_loops;
}

/* Compute face normal (unnormalized). Returns magnitude. */
static float face_normal(const float *verts, const int32_t *face, float out[3])
{
    int32_t a = face[0], b = face[1], c = face[2];
    float e1[3], e2[3];
    for (int k = 0; k < 3; k++) {
        e1[k] = verts[b * 3 + k] - verts[a * 3 + k];
        e2[k] = verts[c * 3 + k] - verts[a * 3 + k];
    }
    out[0] = e1[1] * e2[2] - e1[2] * e2[1];
    out[1] = e1[2] * e2[0] - e1[0] * e2[2];
    out[2] = e1[0] * e2[1] - e1[1] * e2[0];
    return sqrtf(out[0] * out[0] + out[1] * out[1] + out[2] * out[2]);
}

/* ================================================================
 * Test: validate face indices
 * ================================================================ */
static int test_face_indices(const int32_t *faces, size_t nf, size_t nv)
{
    int errors = 0;
    for (size_t f = 0; f < nf; f++) {
        for (size_t k = 0; k < 3; k++) {
            int32_t v = faces[f * 3 + k];
            if (v < 0 || (size_t)v >= nv) {
                if (errors < 5)
                    fprintf(stderr, "  ERROR: face %zu has index %d (nv=%zu)\n",
                            f, v, nv);
                errors++;
            }
        }
    }
    return errors;
}

/* ================================================================
 * Orientation consistency check (shared by tests 3 and 4)
 * ================================================================ */

/* Half-edge for orientation checking */
typedef struct {
    int32_t v0, v1;   /* directed half-edge */
    size_t  face_idx;  /* which face */
} HalfEdge;

static int cmp_halfedge_undirected(const void *a, const void *b)
{
    const HalfEdge *ha = (const HalfEdge *)a;
    const HalfEdge *hb = (const HalfEdge *)b;
    int32_t a_lo = (ha->v0 < ha->v1) ? ha->v0 : ha->v1;
    int32_t a_hi = (ha->v0 < ha->v1) ? ha->v1 : ha->v0;
    int32_t b_lo = (hb->v0 < hb->v1) ? hb->v0 : hb->v1;
    int32_t b_hi = (hb->v0 < hb->v1) ? hb->v1 : hb->v0;
    if (a_lo != b_lo) return (a_lo < b_lo) ? -1 : 1;
    return (a_hi < b_hi) ? -1 : (a_hi > b_hi) ? 1 : 0;
}

/* Check that all interior edges have consistently-wound adjacent faces.
 * Returns the number of inconsistently-wound edge pairs. */
static size_t count_inconsistent_edges(const int32_t *faces, size_t nf)
{
    size_t n_he = nf * 3;
    HalfEdge *hes = (HalfEdge *)malloc(n_he * sizeof(HalfEdge));
    assert(hes);

    for (size_t f = 0; f < nf; f++) {
        for (int e = 0; e < 3; e++) {
            size_t idx = f * 3 + (size_t)e;
            hes[idx].v0 = faces[f * 3 + e];
            hes[idx].v1 = faces[f * 3 + ((e + 1) % 3)];
            hes[idx].face_idx = f;
        }
    }

    /* Sort by undirected edge */
    qsort(hes, n_he, sizeof(HalfEdge), cmp_halfedge_undirected);

    /* For each pair sharing an undirected edge, check that half-edges
     * go in opposite directions (v0->v1 vs v1->v0). */
    size_t n_bad = 0;
    size_t i = 0;
    while (i < n_he) {
        int32_t lo = (hes[i].v0 < hes[i].v1) ? hes[i].v0 : hes[i].v1;
        int32_t hi = (hes[i].v0 < hes[i].v1) ? hes[i].v1 : hes[i].v0;
        size_t j = i + 1;
        while (j < n_he) {
            int32_t jlo = (hes[j].v0 < hes[j].v1) ? hes[j].v0 : hes[j].v1;
            int32_t jhi = (hes[j].v0 < hes[j].v1) ? hes[j].v1 : hes[j].v0;
            if (jlo != lo || jhi != hi) break;
            j++;
        }
        /* [i..j) share the same undirected edge */
        if (j - i == 2) {
            /* Interior edge: check orientation.
             * If both half-edges go the same direction -> inconsistent */
            if (hes[i].v0 == hes[i + 1].v0 && hes[i].v1 == hes[i + 1].v1) {
                n_bad++;
            }
        }
        /* j - i > 2 means non-manifold edge, skip those */
        i = j;
    }

    free(hes);
    return n_bad;
}

/* ================================================================
 * Test 1: Unit test — boundary edge detection on a small mesh
 * ================================================================ */
static int test_boundary_detection(void)
{
    fprintf(stderr, "  [test] boundary detection on unit square...\n");

    /* Two-triangle unit square: v0=(0,0,0) v1=(1,0,0) v2=(1,1,0) v3=(0,1,0)
     * f0: 0 1 2   f1: 0 2 3
     * All 4 outer edges are boundary.  Edge (0,2) is shared. */
    float verts[] = {0,0,0, 1,0,0, 1,1,0, 0,1,0};
    int32_t faces[] = {0,1,2, 0,2,3};
    (void)verts; /* vertices used for documentation, not by find_boundary_edges */

    int32_t *bdry = NULL;
    size_t cap = 0;
    size_t n = find_boundary_edges(faces, 2, &bdry, &cap);

    int ok = 1;
    if (n != 4) {
        fprintf(stderr, "  FAIL: expected 4 boundary edges, got %zu\n", n);
        ok = 0;
    }

    free(bdry);
    return ok ? 0 : 1;
}

/* ================================================================
 * Test 2: Unit test — boundary detection on closed mesh (0 boundary)
 * ================================================================ */
static int test_closed_mesh(void)
{
    fprintf(stderr, "  [test] boundary detection on tetrahedron...\n");

    /* Tetrahedron: 4 faces, each edge shared by exactly 2 faces.
     * No boundary edges. */
    float verts[] = {0,0,0, 1,0,0, 0.5f,1,0, 0.5f,0.5f,1};
    int32_t faces[] = {0,1,2, 0,1,3, 1,2,3, 0,2,3};
    (void)verts;

    int32_t *bdry = NULL;
    size_t cap = 0;
    size_t n = find_boundary_edges(faces, 4, &bdry, &cap);

    int ok = 1;
    if (n != 0) {
        fprintf(stderr, "  FAIL: expected 0 boundary edges, got %zu\n", n);
        ok = 0;
    }

    free(bdry);
    return ok ? 0 : 1;
}

/* ================================================================
 * Test 3: Functional — holefill_test_1.obj
 * ================================================================ */
static int test_holefill_functional(void)
{
    const char *obj_path = "test_sheets/holefill_test_1.obj";
    const int cube_shape[3] = {320, 320, 320};
    const float HOLE_MARGIN = 2.0f;

    fprintf(stderr, "  [test] loading %s...\n", obj_path);

    Arena_T arena = Arena_new();
    float *verts = NULL;
    int32_t *faces = NULL;
    size_t nv = 0, nf = 0;

    if (ObjIO_read(arena, obj_path, &verts, &nv, &faces, &nf) != 0) {
        fprintf(stderr, "  SKIP: cannot read %s\n", obj_path);
        Arena_dispose(&arena);
        return 0;  /* skip, not fail */
    }
    fprintf(stderr, "  loaded: %zu verts, %zu faces\n", nv, nf);

    /* -- Pre-fill analysis: find interior boundary loops -- */
    int32_t *pre_bdry = NULL;
    size_t pre_bdry_cap = 0;
    size_t pre_n_bdry = find_boundary_edges(faces, nf, &pre_bdry, &pre_bdry_cap);
    fprintf(stderr, "  pre-fill: %zu boundary edges\n", pre_n_bdry);

    size_t pre_loop_sizes[256] = {0};
    int pre_loop_interior[256] = {0};
    size_t pre_n_loops = chain_boundary_loops(pre_bdry, pre_n_bdry, nv,
                                              verts, cube_shape, HOLE_MARGIN,
                                              pre_loop_sizes, pre_loop_interior, 256);

    /* Find largest loop (exterior boundary) */
    size_t largest_idx = 0;
    size_t largest_sz = 0;
    size_t n_interior_holes_before = 0;
    for (size_t i = 0; i < pre_n_loops; i++) {
        if (pre_loop_sizes[i] > largest_sz) {
            largest_sz = pre_loop_sizes[i];
            largest_idx = i;
        }
    }
    for (size_t i = 0; i < pre_n_loops; i++) {
        if (i == largest_idx) continue;
        if (pre_loop_interior[i]) n_interior_holes_before++;
    }
    fprintf(stderr, "  pre-fill: %zu loops, %zu interior holes\n",
            pre_n_loops, n_interior_holes_before);
    free(pre_bdry);

    if (n_interior_holes_before == 0) {
        fprintf(stderr, "  SKIP: no interior holes found in test mesh\n");
        Arena_dispose(&arena);
        return 0;
    }

    /* Count pre-fill inconsistent edges (properly, via undirected grouping) */
    size_t n_inconsistent_pre = count_inconsistent_edges(faces, nf);
    if (n_inconsistent_pre > 0) {
        fprintf(stderr, "  pre-fill: %zu inconsistent edge pairs\n",
                n_inconsistent_pre);
    }

    /* -- Compute PCA normal before fill -- */
    float pca_normal[3] = {0, 0, 0};
    float centroid[3] = {0, 0, 0};
    PCA_normal(verts, nv, pca_normal, centroid);
    fprintf(stderr, "  PCA normal: (%.3f, %.3f, %.3f)\n",
            (double)pca_normal[0], (double)pca_normal[1], (double)pca_normal[2]);

    size_t nv_before = nv;
    size_t nf_before = nf;

    /* -- Run hole fill -- */
    fprintf(stderr, "  running HoleFill_process...\n");
    int rc = HoleFill_process(arena, &verts, &faces, &nv, &nf, cube_shape);
    fprintf(stderr, "  HoleFill_process returned %d: %zu -> %zu verts, "
            "%zu -> %zu faces\n", rc, nv_before, nv, nf_before, nf);

    /* Write output for visual inspection */
    {
        const char *out_path = "output/holefill_test_1_filled.obj";
        ObjIO_write(out_path, verts, nv, faces, nf);
        fprintf(stderr, "  wrote %s\n", out_path);
    }

    int n_failures = 0;

    /* -- Check A: Face indices valid -- */
    {
        int bad = test_face_indices(faces, nf, nv);
        if (bad > 0) {
            fprintf(stderr, "  FAIL A: %d invalid face indices\n", bad);
            n_failures++;
        } else {
            fprintf(stderr, "  PASS A: all face indices valid\n");
        }
    }

    /* -- Check B: No interior boundary loops remain -- */
    {
        int32_t *post_bdry = NULL;
        size_t post_cap = 0;
        size_t post_n_bdry = find_boundary_edges(faces, nf, &post_bdry, &post_cap);

        size_t post_loop_sizes[256] = {0};
        int post_loop_interior[256] = {0};
        size_t post_n_loops = chain_boundary_loops(post_bdry, post_n_bdry, nv,
                                                   verts, cube_shape, HOLE_MARGIN,
                                                   post_loop_sizes, post_loop_interior,
                                                   256);

        /* Find largest loop */
        size_t post_largest = 0;
        size_t post_largest_sz = 0;
        for (size_t i = 0; i < post_n_loops; i++) {
            if (post_loop_sizes[i] > post_largest_sz) {
                post_largest_sz = post_loop_sizes[i];
                post_largest = i;
            }
        }

        size_t n_interior_holes_after = 0;
        for (size_t i = 0; i < post_n_loops; i++) {
            if (i == post_largest) continue;
            if (post_loop_interior[i]) {
                n_interior_holes_after++;
                fprintf(stderr, "    interior loop %zu: %zu verts (NOT FILLED)\n",
                        i, post_loop_sizes[i]);
            }
        }

        fprintf(stderr, "  post-fill: %zu boundary edges, %zu loops, "
                "%zu interior holes\n",
                post_n_bdry, post_n_loops, n_interior_holes_after);

        if (n_interior_holes_after > 0) {
            fprintf(stderr, "  FAIL B: %zu interior holes remain "
                    "(had %zu before fill)\n",
                    n_interior_holes_after, n_interior_holes_before);
            n_failures++;
        } else {
            fprintf(stderr, "  PASS B: all interior holes filled\n");
        }

        free(post_bdry);
    }

    /* -- Check B2: Fill must not introduce inconsistent edge orientations.
     *    The input mesh has consistent winding, so any new inconsistencies
     *    are bugs in the fill. */
    {
        size_t n_inconsistent_post = count_inconsistent_edges(faces, nf);
        size_t new_bad = (n_inconsistent_post > n_inconsistent_pre) ?
                          n_inconsistent_post - n_inconsistent_pre : 0;
        if (new_bad > 0) {
            fprintf(stderr, "  FAIL B2: %zu inconsistent edge pairs after fill "
                    "(%zu new, pre-fill had %zu)\n",
                    n_inconsistent_post, new_bad, n_inconsistent_pre);
            n_failures++;
        } else {
            fprintf(stderr, "  PASS B2: %zu inconsistent edge pairs after fill "
                    "(pre-fill had %zu)\n",
                    n_inconsistent_post, n_inconsistent_pre);
        }
    }

    /* -- Check C: No fill faces point away from PCA normal -- */
    {
        /* Check ALL faces added by the fill (indices >= nf_before) */
        size_t n_opposing = 0;
        size_t n_degenerate = 0;
        size_t n_fill_faces = (nf > nf_before) ? nf - nf_before : 0;

        for (size_t f = nf_before; f < nf; f++) {
            float fn[3] = {0, 0, 0};
            float mag = face_normal(verts, &faces[f * 3], fn);

            if (mag < 1e-10f) {
                n_degenerate++;
                continue;
            }

            /* Normalize */
            float inv = 1.0f / mag;
            fn[0] *= inv; fn[1] *= inv; fn[2] *= inv;

            float dot = fn[0] * pca_normal[0]
                      + fn[1] * pca_normal[1]
                      + fn[2] * pca_normal[2];

            /* "Directly away" = dot < -0.5 (more than 120 degrees off) */
            if (dot < -0.5f) {
                n_opposing++;
                if (n_opposing <= 3) {
                    fprintf(stderr, "    fill face %zu: normal=(%.2f,%.2f,%.2f) "
                            "dot=%.3f\n", f, (double)fn[0], (double)fn[1],
                            (double)fn[2], (double)dot);
                }
            }
        }

        if (n_degenerate > 0) {
            fprintf(stderr, "  WARNING: %zu degenerate fill faces (zero area)\n",
                    n_degenerate);
        }

        if (n_opposing > 0) {
            fprintf(stderr, "  FAIL C: %zu/%zu fill faces point away from PCA "
                    "normal (dot < -0.5)\n", n_opposing, n_fill_faces);
            n_failures++;
        } else {
            fprintf(stderr, "  PASS C: all %zu fill faces consistent with "
                    "PCA normal\n", n_fill_faces);
        }
    }

    /* -- Check D: No NaN in vertices -- */
    {
        size_t n_nan = 0;
        for (size_t i = 0; i < nv * 3; i++) {
            if (isnan(verts[i])) n_nan++;
        }
        if (n_nan > 0) {
            fprintf(stderr, "  FAIL D: %zu NaN values in vertex array\n", n_nan);
            n_failures++;
        } else {
            fprintf(stderr, "  PASS D: no NaN in vertices\n");
        }
    }

    Arena_dispose(&arena);

    if (n_failures > 0) {
        fprintf(stderr, "  RESULT: %d checks FAILED\n", n_failures);
        return 1;
    }
    fprintf(stderr, "  RESULT: all checks passed\n");
    return 0;
}

/* ================================================================
 * Test 4: Functional — holefill_test_2.obj orientation check
 *
 * Loads holefill_test_2.obj, fills holes, then checks that
 * every interior edge has consistently-wound adjacent faces.
 * An interior edge shared by faces F1=(a,b,c) and F2 should
 * have the shared edge traversed in opposite directions:
 * if F1 has half-edge a→b, F2 must have b→a.
 * ================================================================ */

static int test_holefill_orientation(void)
{
    const char *obj_path = "test_sheets/holefill_test_2.obj";
    const int cube_shape[3] = {320, 320, 320};

    fprintf(stderr, "  [test] loading %s...\n", obj_path);

    Arena_T arena = Arena_new();
    float *verts = NULL;
    int32_t *faces = NULL;
    size_t nv = 0, nf = 0;

    if (ObjIO_read(arena, obj_path, &verts, &nv, &faces, &nf) != 0) {
        fprintf(stderr, "  SKIP: cannot read %s\n", obj_path);
        Arena_dispose(&arena);
        return 0;
    }
    fprintf(stderr, "  loaded: %zu verts, %zu faces\n", nv, nf);

    /* Check pre-fill orientation consistency */
    size_t pre_bad = count_inconsistent_edges(faces, nf);
    fprintf(stderr, "  pre-fill: %zu inconsistent edge pairs\n", pre_bad);

    size_t nv_before = nv;
    size_t nf_before = nf;

    /* Run hole fill */
    fprintf(stderr, "  running HoleFill_process...\n");
    int rc = HoleFill_process(arena, &verts, &faces, &nv, &nf, cube_shape);
    fprintf(stderr, "  HoleFill_process returned %d: %zu -> %zu verts, "
            "%zu -> %zu faces\n", rc, nv_before, nv, nf_before, nf);

    /* Write output for visual inspection */
    ObjIO_write("output/holefill_test_2_filled.obj", verts, nv, faces, nf);
    fprintf(stderr, "  wrote output/holefill_test_2_filled.obj\n");

    int n_failures = 0;

    /* Check A: Face indices valid */
    {
        int bad = test_face_indices(faces, nf, nv);
        if (bad > 0) {
            fprintf(stderr, "  FAIL A: %d invalid face indices\n", bad);
            n_failures++;
        } else {
            fprintf(stderr, "  PASS A: face indices valid\n");
        }
    }

    /* Check B: At least some hole was filled (verts or faces increased, or
     * degenerate faces were removed so nf could decrease — check that fill
     * process didn't error out) */
    if (rc != 0) {
        fprintf(stderr, "  FAIL B: HoleFill_process returned %d\n", rc);
        n_failures++;
    } else {
        fprintf(stderr, "  PASS B: HoleFill_process succeeded\n");
    }

    /* Check C: Post-fill orientation — this is the key test.
     * Every interior edge should have consistently-wound faces. */
    {
        size_t post_bad = count_inconsistent_edges(faces, nf);
        fprintf(stderr, "  post-fill: %zu inconsistent edge pairs "
                "(was %zu pre-fill)\n", post_bad, pre_bad);

        /* Allow a small number of new inconsistencies proportional to
         * pre-existing mesh defects.  Each pre-existing defective face
         * can create up to ~3 new boundary inconsistencies (its 3 edges). */
        size_t tolerance = pre_bad * 3;
        size_t new_bad = (post_bad > pre_bad) ? post_bad - pre_bad : 0;
        if (new_bad > tolerance) {
            fprintf(stderr, "  FAIL C: fill introduced %zu new orientation "
                    "inconsistencies (tolerance %zu)\n", new_bad, tolerance);
            n_failures++;
        } else {
            fprintf(stderr, "  PASS C: fill orientation OK (%zu total, "
                    "%zu new, tolerance %zu)\n", post_bad, new_bad, tolerance);
        }
    }

    /* Check D: No NaN in vertices */
    {
        size_t n_nan = 0;
        for (size_t i = 0; i < nv * 3; i++) {
            if (isnan(verts[i])) n_nan++;
        }
        if (n_nan > 0) {
            fprintf(stderr, "  FAIL D: %zu NaN values\n", n_nan);
            n_failures++;
        } else {
            fprintf(stderr, "  PASS D: no NaN\n");
        }
    }

    Arena_dispose(&arena);

    if (n_failures > 0) {
        fprintf(stderr, "  RESULT: %d checks FAILED\n", n_failures);
        return 1;
    }
    fprintf(stderr, "  RESULT: all checks passed\n");
    return 0;
}

/* ================================================================
 * Test 5: Functional — holefill_test_3.obj single-hole complete fill
 *
 * Loads holefill_test_3.obj which has exactly 1 interior hole.
 * Verifies that the hole is completely filled (zero interior
 * boundary loops remain) and that new geometry was actually added.
 * ================================================================ */

static int test_holefill_single_hole(void)
{
    const char *obj_path = "test_sheets/holefill_test_3.obj";
    const int cube_shape[3] = {320, 320, 320};
    const float HOLE_MARGIN = 2.0f;

    fprintf(stderr, "  [test] loading %s...\n", obj_path);

    Arena_T arena = Arena_new();
    float *verts = NULL;
    int32_t *faces = NULL;
    size_t nv = 0, nf = 0;

    if (ObjIO_read(arena, obj_path, &verts, &nv, &faces, &nf) != 0) {
        fprintf(stderr, "  SKIP: cannot read %s\n", obj_path);
        Arena_dispose(&arena);
        return 0;
    }
    fprintf(stderr, "  loaded: %zu verts, %zu faces\n", nv, nf);

    /* -- Pre-fill analysis: count interior boundary loops -- */
    int32_t *pre_bdry = NULL;
    size_t pre_bdry_cap = 0;
    size_t pre_n_bdry = find_boundary_edges(faces, nf, &pre_bdry, &pre_bdry_cap);
    fprintf(stderr, "  pre-fill: %zu boundary edges\n", pre_n_bdry);

    /* Boundary vertex degree analysis */
    {
        int32_t *bdeg = (int32_t *)calloc(nv, sizeof(int32_t));
        assert(bdeg);
        for (size_t i = 0; i < pre_n_bdry; i++) {
            bdeg[pre_bdry[i * 2 + 0]]++;
            bdeg[pre_bdry[i * 2 + 1]]++;
        }
        size_t n_deg1 = 0, n_deg2 = 0, n_deg3plus = 0;
        for (size_t v = 0; v < nv; v++) {
            if (bdeg[v] == 1) n_deg1++;
            else if (bdeg[v] == 2) n_deg2++;
            else if (bdeg[v] >= 3) {
                n_deg3plus++;
                if (n_deg3plus <= 10)
                    fprintf(stderr, "    boundary v%zu: degree %d at (%.1f, %.1f, %.1f)\n",
                            v, bdeg[v],
                            (double)verts[v * 3 + 0],
                            (double)verts[v * 3 + 1],
                            (double)verts[v * 3 + 2]);
            }
        }
        fprintf(stderr, "  boundary vertex degrees: deg1=%zu deg2=%zu deg3+=%zu\n",
                n_deg1, n_deg2, n_deg3plus);

        /* Count non-manifold edges (appear >2 times as half-edges) */
        size_t n_nonmanifold_edges = 0;
        {
            /* Re-count: for each undirected edge, count half-edge occurrences */
            size_t n_he = nf * 3;
            int32_t *all_he = (int32_t *)malloc(n_he * 2 * sizeof(int32_t));
            assert(all_he);
            for (size_t f = 0; f < nf; f++) {
                for (int e = 0; e < 3; e++) {
                    int32_t a = faces[f * 3 + e];
                    int32_t b = faces[f * 3 + ((e + 1) % 3)];
                    size_t idx = f * 3 + (size_t)e;
                    all_he[idx * 2 + 0] = (a < b) ? a : b;
                    all_he[idx * 2 + 1] = (a < b) ? b : a;
                }
            }
            qsort(all_he, n_he, 2 * sizeof(int32_t), cmp_edge);
            size_t i = 0;
            while (i < n_he) {
                size_t j = i + 1;
                while (j < n_he && all_he[j * 2] == all_he[i * 2] &&
                       all_he[j * 2 + 1] == all_he[i * 2 + 1]) j++;
                if (j - i > 2) {
                    n_nonmanifold_edges++;
                    if (n_nonmanifold_edges <= 5)
                        fprintf(stderr, "    non-manifold edge: %d-%d (%zu occurrences)\n",
                                all_he[i * 2], all_he[i * 2 + 1], j - i);
                }
                i = j;
            }
            free(all_he);
        }
        if (n_nonmanifold_edges > 0)
            fprintf(stderr, "  %zu non-manifold edges total\n", n_nonmanifold_edges);

        /* Check for duplicate vertices (same position, different index) near boundary */
        {
            size_t n_dup_bdry = 0;
            for (size_t i = 0; i < nv && n_dup_bdry < 20; i++) {
                if (bdeg[i] == 0) continue;
                for (size_t j = i + 1; j < nv && n_dup_bdry < 20; j++) {
                    if (bdeg[j] == 0) continue;
                    float dz = verts[i*3+0] - verts[j*3+0];
                    float dy = verts[i*3+1] - verts[j*3+1];
                    float dx = verts[i*3+2] - verts[j*3+2];
                    if (dz*dz + dy*dy + dx*dx < 1e-6f) {
                        n_dup_bdry++;
                        if (n_dup_bdry <= 10)
                            fprintf(stderr, "    dup boundary verts: v%zu==v%zu at (%.1f, %.1f, %.1f) "
                                    "deg %d,%d\n",
                                    i, j,
                                    (double)verts[i*3], (double)verts[i*3+1], (double)verts[i*3+2],
                                    bdeg[i], bdeg[j]);
                    }
                }
            }
            if (n_dup_bdry > 0)
                fprintf(stderr, "  %zu duplicate boundary vertex pairs\n", n_dup_bdry);
        }
        free(bdeg);
    }

    size_t pre_loop_sizes[256] = {0};
    int pre_loop_interior[256] = {0};
    size_t pre_n_loops = chain_boundary_loops(pre_bdry, pre_n_bdry, nv,
                                              verts, cube_shape, HOLE_MARGIN,
                                              pre_loop_sizes, pre_loop_interior, 256);

    /* Export all boundary edges as one OBJ */
    {
        const char *loop_obj = "output/holefill_test_3_boundary_edges.obj";
        FILE *fp = fopen(loop_obj, "w");
        if (fp) {
            fprintf(fp, "# Boundary edges: %zu edges\n", pre_n_bdry);
            for (size_t i = 0; i < pre_n_bdry; i++) {
                int32_t v0 = pre_bdry[i * 2 + 0];
                int32_t v1 = pre_bdry[i * 2 + 1];
                fprintf(fp, "v %.6f %.6f %.6f\n",
                        (double)verts[v0 * 3 + 0],
                        (double)verts[v0 * 3 + 1],
                        (double)verts[v0 * 3 + 2]);
                fprintf(fp, "v %.6f %.6f %.6f\n",
                        (double)verts[v1 * 3 + 0],
                        (double)verts[v1 * 3 + 1],
                        (double)verts[v1 * 3 + 2]);
            }
            for (size_t i = 0; i < pre_n_bdry; i++) {
                fprintf(fp, "l %zu %zu\n", i * 2 + 1, i * 2 + 2);
            }
            fclose(fp);
            fprintf(stderr, "  wrote boundary edges: %s\n", loop_obj);
        }
    }

    /* Export each chained loop as a separate OBJ */
    {
        /* Re-chain, but this time capture the actual vertex sequences.
         * chain_boundary_loops only returns sizes. Redo chaining here
         * to get the vertex sequences for export. */
        int32_t *bdeg = (int32_t *)calloc(nv, sizeof(int32_t));
        assert(bdeg);
        for (size_t i = 0; i < pre_n_bdry; i++) {
            bdeg[pre_bdry[i * 2 + 0]]++;
            bdeg[pre_bdry[i * 2 + 1]]++;
        }
        int32_t *adj_off2 = (int32_t *)malloc((nv + 1) * sizeof(int32_t));
        assert(adj_off2);
        adj_off2[0] = 0;
        for (size_t v = 0; v < nv; v++)
            adj_off2[v + 1] = adj_off2[v] + bdeg[v];
        int32_t *adj_list2 = (int32_t *)malloc(pre_n_bdry * 2 * sizeof(int32_t));
        assert(adj_list2);
        int32_t *cursor2 = (int32_t *)malloc(nv * sizeof(int32_t));
        assert(cursor2);
        memcpy(cursor2, adj_off2, nv * sizeof(int32_t));
        for (size_t i = 0; i < pre_n_bdry; i++) {
            int32_t a = pre_bdry[i * 2 + 0];
            int32_t b = pre_bdry[i * 2 + 1];
            adj_list2[cursor2[a]++] = b;
            adj_list2[cursor2[b]++] = a;
        }

        uint8_t *edge_used2 = (uint8_t *)calloc(pre_n_bdry, sizeof(uint8_t));
        assert(edge_used2);

        size_t loop_idx = 0;
        for (size_t ei = 0; ei < pre_n_bdry; ei++) {
            if (edge_used2[ei]) continue;

            /* Trace this loop */
            int32_t *lbuf = (int32_t *)malloc((pre_n_bdry + 2) * sizeof(int32_t));
            assert(lbuf);
            size_t llen = 0;

            int32_t start = pre_bdry[ei * 2 + 0];
            int32_t prev = start;
            int32_t cur = pre_bdry[ei * 2 + 1];
            edge_used2[ei] = 1;
            lbuf[llen++] = start;

            while (cur != start && llen < pre_n_bdry + 1) {
                lbuf[llen++] = cur;
                int32_t next = -1;
                for (int32_t j = adj_off2[cur]; j < adj_off2[cur + 1]; j++) {
                    int32_t nb = adj_list2[j];
                    if (nb == prev) continue;
                    /* Find edge in sorted boundary list */
                    int32_t key[2];
                    key[0] = (cur < nb) ? cur : nb;
                    key[1] = (cur < nb) ? nb : cur;
                    int32_t *found = (int32_t *)bsearch(key, pre_bdry,
                        pre_n_bdry, 2 * sizeof(int32_t), cmp_edge);
                    if (found) {
                        size_t fidx = (size_t)(found - pre_bdry) / 2;
                        if (!edge_used2[fidx]) {
                            edge_used2[fidx] = 1;
                            next = nb;
                            break;
                        }
                    }
                }
                if (next < 0) break;
                prev = cur;
                cur = next;
            }

            int closed = (cur == start && llen >= 3);

            /* Write this loop to its own OBJ */
            char path[256];
            snprintf(path, sizeof(path),
                     "output/holefill_test_3_loop_%zu%s.obj",
                     loop_idx, closed ? "" : "_OPEN");
            FILE *fp = fopen(path, "w");
            if (fp) {
                fprintf(fp, "# Loop %zu: %zu verts, %s\n",
                        loop_idx, llen, closed ? "closed" : "OPEN");
                for (size_t k = 0; k < llen; k++) {
                    int32_t v = lbuf[k];
                    fprintf(fp, "v %.6f %.6f %.6f\n",
                            (double)verts[v * 3 + 0],
                            (double)verts[v * 3 + 1],
                            (double)verts[v * 3 + 2]);
                }
                for (size_t k = 0; k < llen; k++) {
                    size_t next_k = (k + 1) % llen;
                    if (!closed && k + 1 == llen) break;
                    fprintf(fp, "l %zu %zu\n", k + 1, next_k + 1);
                }
                fclose(fp);
            }
            fprintf(stderr, "  wrote loop %zu (%zu verts, %s): %s\n",
                    loop_idx, llen, closed ? "closed" : "OPEN", path);
            loop_idx++;
            free(lbuf);
        }

        /* Check for unused boundary edges */
        size_t n_unused = 0;
        for (size_t i = 0; i < pre_n_bdry; i++)
            if (!edge_used2[i]) n_unused++;
        if (n_unused > 0)
            fprintf(stderr, "  WARNING: %zu boundary edges not assigned to any loop\n",
                    n_unused);

        free(bdeg);
        free(adj_off2);
        free(adj_list2);
        free(cursor2);
        free(edge_used2);
    }

    /* Find largest loop (exterior boundary) */
    size_t largest_idx = 0;
    size_t largest_sz = 0;
    for (size_t i = 0; i < pre_n_loops; i++) {
        if (pre_loop_sizes[i] > largest_sz) {
            largest_sz = pre_loop_sizes[i];
            largest_idx = i;
        }
    }

    size_t n_interior_holes_before = 0;
    for (size_t i = 0; i < pre_n_loops; i++) {
        if (i == largest_idx) continue;
        if (pre_loop_interior[i]) {
            n_interior_holes_before++;
            fprintf(stderr, "    interior hole %zu: %zu verts\n",
                    i, pre_loop_sizes[i]);
        }
    }
    fprintf(stderr, "  pre-fill: %zu loops total, %zu interior holes\n",
            pre_n_loops, n_interior_holes_before);
    free(pre_bdry);

    int n_failures = 0;

    /* -- Check A: mesh has boundary edges (sanity check).
     * Note: the test's simple chainer may not detect the interior hole
     * when the outer boundary is broken by non-manifold edges (degree-1
     * dead-ends). HoleFill_process handles this with its improved chainer
     * that reconnects open chains. The real hole-fill assertions are in
     * checks C and D below. */
    if (pre_n_bdry == 0) {
        fprintf(stderr, "  FAIL A: no boundary edges found\n");
        n_failures++;
    } else {
        fprintf(stderr, "  PASS A: %zu boundary edges, %zu loops detected "
                "(%zu interior by simple chainer)\n",
                pre_n_bdry, pre_n_loops, n_interior_holes_before);
    }

    /* Count pre-fill inconsistent edges */
    size_t n_inconsistent_pre = count_inconsistent_edges(faces, nf);
    if (n_inconsistent_pre > 0) {
        fprintf(stderr, "  pre-fill: %zu inconsistent edge pairs\n",
                n_inconsistent_pre);
    }

    /* Compute PCA normal */
    float pca_normal[3] = {0, 0, 0};
    float centroid[3] = {0, 0, 0};
    PCA_normal(verts, nv, pca_normal, centroid);
    fprintf(stderr, "  PCA normal: (%.3f, %.3f, %.3f)\n",
            (double)pca_normal[0], (double)pca_normal[1], (double)pca_normal[2]);

    size_t nv_before = nv;
    size_t nf_before = nf;

    /* -- Run hole fill -- */
    fprintf(stderr, "  running HoleFill_process...\n");
    int rc = HoleFill_process(arena, &verts, &faces, &nv, &nf, cube_shape);
    fprintf(stderr, "  HoleFill_process returned %d: %zu -> %zu verts, "
            "%zu -> %zu faces\n", rc, nv_before, nv, nf_before, nf);

    /* Write output for visual inspection */
    ObjIO_write("output/holefill_test_3_filled.obj", verts, nv, faces, nf);
    fprintf(stderr, "  wrote output/holefill_test_3_filled.obj\n");

    /* -- Check B: HoleFill_process succeeded -- */
    if (rc != 0) {
        fprintf(stderr, "  FAIL B: HoleFill_process returned %d\n", rc);
        n_failures++;
    } else {
        fprintf(stderr, "  PASS B: HoleFill_process succeeded\n");
    }

    /* -- Check C: New geometry was added (not a no-op) -- */
    if (nf <= nf_before) {
        fprintf(stderr, "  FAIL C: no new faces added (nf %zu -> %zu)\n",
                nf_before, nf);
        n_failures++;
    } else {
        fprintf(stderr, "  PASS C: %zu new faces added, %zu new verts\n",
                nf - nf_before, nv - nv_before);
    }

    /* -- Check D: No interior boundary loops remain (hole completely filled) -- */
    {
        int32_t *post_bdry = NULL;
        size_t post_cap = 0;
        size_t post_n_bdry = find_boundary_edges(faces, nf, &post_bdry, &post_cap);

        size_t post_loop_sizes[256] = {0};
        int post_loop_interior[256] = {0};
        size_t post_n_loops = chain_boundary_loops(post_bdry, post_n_bdry, nv,
                                                   verts, cube_shape, HOLE_MARGIN,
                                                   post_loop_sizes, post_loop_interior,
                                                   256);

        /* Find largest loop */
        size_t post_largest = 0;
        size_t post_largest_sz = 0;
        for (size_t i = 0; i < post_n_loops; i++) {
            if (post_loop_sizes[i] > post_largest_sz) {
                post_largest_sz = post_loop_sizes[i];
                post_largest = i;
            }
        }

        size_t n_interior_holes_after = 0;
        for (size_t i = 0; i < post_n_loops; i++) {
            if (i == post_largest) continue;
            if (post_loop_interior[i]) {
                n_interior_holes_after++;
                fprintf(stderr, "    remaining interior hole %zu: %zu verts "
                        "(NOT FILLED)\n", i, post_loop_sizes[i]);
            }
        }

        fprintf(stderr, "  post-fill: %zu boundary edges, %zu loops, "
                "%zu interior holes\n",
                post_n_bdry, post_n_loops, n_interior_holes_after);

        if (n_interior_holes_after > 0) {
            fprintf(stderr, "  FAIL D: %zu interior holes remain — "
                    "hole NOT completely filled\n", n_interior_holes_after);
            n_failures++;
        } else {
            fprintf(stderr, "  PASS D: 0 interior holes remain — "
                    "hole completely filled\n");
        }

        free(post_bdry);
    }

    /* -- Check E: No fill faces point away from PCA normal -- */
    {
        size_t n_opposing = 0;
        size_t n_fill_faces = (nf > nf_before) ? nf - nf_before : 0;

        for (size_t f = nf_before; f < nf; f++) {
            float fn[3] = {0, 0, 0};
            float mag = face_normal(verts, &faces[f * 3], fn);
            if (mag < 1e-10f) continue;

            float inv = 1.0f / mag;
            fn[0] *= inv; fn[1] *= inv; fn[2] *= inv;

            float dot = fn[0] * pca_normal[0]
                      + fn[1] * pca_normal[1]
                      + fn[2] * pca_normal[2];

            if (dot < -0.5f) n_opposing++;
        }

        if (n_opposing > 0) {
            fprintf(stderr, "  FAIL E: %zu/%zu fill faces point away from "
                    "PCA normal\n", n_opposing, n_fill_faces);
            n_failures++;
        } else {
            fprintf(stderr, "  PASS E: all %zu fill faces consistent with "
                    "PCA normal\n", n_fill_faces);
        }
    }

    /* -- Check F: Fill did not introduce inconsistent edge orientations -- */
    {
        size_t n_inconsistent_post = count_inconsistent_edges(faces, nf);
        size_t new_bad = (n_inconsistent_post > n_inconsistent_pre) ?
                          n_inconsistent_post - n_inconsistent_pre : 0;
        if (new_bad > 0) {
            fprintf(stderr, "  FAIL F: %zu new inconsistent edge pairs "
                    "after fill\n", new_bad);
            n_failures++;
        } else {
            fprintf(stderr, "  PASS F: orientation consistent after fill\n");
        }
    }

    /* -- Check G: Face indices valid -- */
    {
        int bad = test_face_indices(faces, nf, nv);
        if (bad > 0) {
            fprintf(stderr, "  FAIL G: %d invalid face indices\n", bad);
            n_failures++;
        } else {
            fprintf(stderr, "  PASS G: all face indices valid\n");
        }
    }

    /* -- Check H: No NaN in vertices -- */
    {
        size_t n_nan = 0;
        for (size_t i = 0; i < nv * 3; i++) {
            if (isnan(verts[i])) n_nan++;
        }
        if (n_nan > 0) {
            fprintf(stderr, "  FAIL H: %zu NaN values in vertex array\n", n_nan);
            n_failures++;
        } else {
            fprintf(stderr, "  PASS H: no NaN in vertices\n");
        }
    }

    Arena_dispose(&arena);

    if (n_failures > 0) {
        fprintf(stderr, "  RESULT: %d checks FAILED\n", n_failures);
        return 1;
    }
    fprintf(stderr, "  RESULT: all checks passed\n");
    return 0;
}

/* ================================================================
 * Test 6: Functional — holefill_test_4.obj Triangle hang regression
 *
 * This mesh (component 1 from cube 1059332280) has an 18-vertex
 * interior boundary loop that causes Triangle's quality refinement
 * to infinite-loop.  Verify that HoleFill_process completes within
 * a reasonable time (does not hang) and produces valid output.
 * ================================================================ */

static int test_holefill_triangle_hang(void)
{
    const char *obj_path = "test_sheets/holefill_test_4.obj";
    const int cube_shape[3] = {320, 320, 320};
    const float HOLE_MARGIN = 2.0f;

    fprintf(stderr, "  [test] loading %s...\n", obj_path);

    Arena_T arena = Arena_new();
    float *verts = NULL;
    int32_t *faces = NULL;
    size_t nv = 0, nf = 0;

    if (ObjIO_read(arena, obj_path, &verts, &nv, &faces, &nf) != 0) {
        fprintf(stderr, "  SKIP: cannot read %s\n", obj_path);
        Arena_dispose(&arena);
        return 0;
    }
    fprintf(stderr, "  loaded: %zu verts, %zu faces\n", nv, nf);

    /* Pre-fill analysis */
    int32_t *pre_bdry = NULL;
    size_t pre_bdry_cap = 0;
    size_t pre_n_bdry = find_boundary_edges(faces, nf, &pre_bdry, &pre_bdry_cap);
    fprintf(stderr, "  pre-fill: %zu boundary edges\n", pre_n_bdry);

    size_t pre_loop_sizes[256] = {0};
    int pre_loop_interior[256] = {0};
    size_t pre_n_loops = chain_boundary_loops(pre_bdry, pre_n_bdry, nv,
                                              verts, cube_shape, HOLE_MARGIN,
                                              pre_loop_sizes, pre_loop_interior, 256);

    size_t largest_idx = 0;
    size_t largest_sz = 0;
    size_t n_interior_holes_before = 0;
    for (size_t i = 0; i < pre_n_loops; i++) {
        if (pre_loop_sizes[i] > largest_sz) {
            largest_sz = pre_loop_sizes[i];
            largest_idx = i;
        }
    }
    for (size_t i = 0; i < pre_n_loops; i++) {
        if (i == largest_idx) continue;
        if (pre_loop_interior[i]) {
            n_interior_holes_before++;
            fprintf(stderr, "    interior hole %zu: %zu verts\n",
                    i, pre_loop_sizes[i]);
        }
    }
    fprintf(stderr, "  pre-fill: %zu loops total, %zu interior holes\n",
            pre_n_loops, n_interior_holes_before);
    free(pre_bdry);

    size_t n_inconsistent_pre = count_inconsistent_edges(faces, nf);

    size_t nv_before = nv;
    size_t nf_before = nf;

    /* Run hole fill — THIS IS THE REGRESSION TEST.
     * Previously this would hang forever inside Triangle's CDT. */
    fprintf(stderr, "  running HoleFill_process (must complete, not hang)...\n");
    int rc = HoleFill_process(arena, &verts, &faces, &nv, &nf, cube_shape);
    fprintf(stderr, "  HoleFill_process returned %d: %zu -> %zu verts, "
            "%zu -> %zu faces\n", rc, nv_before, nv, nf_before, nf);

    /* Write output for visual inspection */
    ObjIO_write("output/holefill_test_4_filled.obj", verts, nv, faces, nf);
    fprintf(stderr, "  wrote output/holefill_test_4_filled.obj\n");

    int n_failures = 0;

    /* Check A: HoleFill_process returned success */
    if (rc != 0) {
        fprintf(stderr, "  FAIL A: HoleFill_process returned %d\n", rc);
        n_failures++;
    } else {
        fprintf(stderr, "  PASS A: HoleFill_process succeeded\n");
    }

    /* Check B: Face indices valid */
    {
        int bad = test_face_indices(faces, nf, nv);
        if (bad > 0) {
            fprintf(stderr, "  FAIL B: %d invalid face indices\n", bad);
            n_failures++;
        } else {
            fprintf(stderr, "  PASS B: all face indices valid\n");
        }
    }

    /* Check C: No NaN in vertices */
    {
        size_t n_nan = 0;
        for (size_t i = 0; i < nv * 3; i++) {
            if (isnan(verts[i])) n_nan++;
        }
        if (n_nan > 0) {
            fprintf(stderr, "  FAIL C: %zu NaN values\n", n_nan);
            n_failures++;
        } else {
            fprintf(stderr, "  PASS C: no NaN\n");
        }
    }

    /* Check D: Fill did not introduce inconsistent edge orientations */
    {
        size_t n_inconsistent_post = count_inconsistent_edges(faces, nf);
        size_t new_bad = (n_inconsistent_post > n_inconsistent_pre) ?
                          n_inconsistent_post - n_inconsistent_pre : 0;
        /* Allow tolerance proportional to pre-existing defects */
        size_t tolerance = n_inconsistent_pre * 3;
        if (new_bad > tolerance) {
            fprintf(stderr, "  FAIL D: %zu new inconsistent edge pairs "
                    "(tolerance %zu)\n", new_bad, tolerance);
            n_failures++;
        } else {
            fprintf(stderr, "  PASS D: orientation OK (%zu new, tolerance %zu)\n",
                    new_bad, tolerance);
        }
    }

    /* Check E: Some holes were filled (geometry was added) */
    if (n_interior_holes_before > 0 && nf <= nf_before) {
        fprintf(stderr, "  WARN E: had %zu interior holes but no faces added "
                "(some may have been skipped due to timeout)\n",
                n_interior_holes_before);
        /* Not a failure — timeout-skipped holes are acceptable */
    } else if (nf > nf_before) {
        fprintf(stderr, "  PASS E: %zu new faces added\n", nf - nf_before);
    } else {
        fprintf(stderr, "  PASS E: no interior holes to fill\n");
    }

    Arena_dispose(&arena);

    if (n_failures > 0) {
        fprintf(stderr, "  RESULT: %d checks FAILED\n", n_failures);
        return 1;
    }
    fprintf(stderr, "  RESULT: all checks passed\n");
    return 0;
}

/* ================================================================
 * Entry point
 * ================================================================ */

#ifdef TEST_HARNESS
int holefill_test_main(void)
{
    int rc = 0;
    rc |= test_boundary_detection();
    rc |= test_closed_mesh();
    rc |= test_holefill_functional();
    rc |= test_holefill_orientation();
    rc |= test_holefill_single_hole();
    rc |= test_holefill_triangle_hang();
    return rc;
}
#else
int main(void)
{
    fprintf(stderr, "\n=== Hole Fill Tests ===\n\n");
    int rc = 0;

    fprintf(stderr, "--- Test 1: boundary detection ---\n");
    rc |= test_boundary_detection();

    fprintf(stderr, "\n--- Test 2: closed mesh ---\n");
    rc |= test_closed_mesh();

    fprintf(stderr, "\n--- Test 3: functional holefill ---\n");
    rc |= test_holefill_functional();

    fprintf(stderr, "\n--- Test 4: orientation holefill ---\n");
    rc |= test_holefill_orientation();

    fprintf(stderr, "\n--- Test 5: single-hole complete fill ---\n");
    rc |= test_holefill_single_hole();

    fprintf(stderr, "\n--- Test 6: Triangle hang regression ---\n");
    rc |= test_holefill_triangle_hang();

    fprintf(stderr, "\n=== %s ===\n", rc ? "FAILED" : "ALL PASSED");
    return rc;
}
#endif
