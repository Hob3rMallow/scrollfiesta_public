/*
 * weld_recoarsen_test -- unit tests for WeldCleanup_recoarsen_seam.
 *
 * Scene A: a fine (1.5-vox) open rectangular grid patch crossing a synthetic
 * seam plane. Asserts: the open boundary-edge set is BIT-IDENTICAL after
 * recoarsening (the hierarchical-composability guarantee), in-band faces
 * coarsen substantially, faces with no in-band vertex are untouched, the mesh
 * stays one connected component and edge-manifold, and stats are monotone.
 * Scene B: no planes -> exact no-op. Scene C: the max_collapse_len cap binds
 * even when collapse_below exceeds it.
 *
 * Verts are (z,y,x); the seam plane is x = 128.
 */
#include "../../src/remesh/weld_cleanup.h"
#include "../../src/common/pipeline_constants.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int g_fail = 0;
#define CHECK(cond, ...) do { \
    if (!(cond)) { g_fail++; fprintf(stderr, "  FAIL: " __VA_ARGS__); fprintf(stderr, "\n"); } \
} while (0)

static int i64_cmp(const void *pa, const void *pb)
{
    int64_t a = *(const int64_t *)pa, b = *(const int64_t *)pb;
    return a < b ? -1 : (a > b ? 1 : 0);
}

/* Sorted array of run-1 (boundary) undirected edges. Caller frees. */
static int64_t *boundary_edges(const int32_t *f, size_t nf, size_t *out_n)
{
    size_t ne = nf * 3, n = 0, i;
    int64_t *keys = (int64_t *)malloc((ne ? ne : 1) * sizeof(int64_t));
    int64_t *bnd = (int64_t *)malloc((ne ? ne : 1) * sizeof(int64_t));
    for (size_t t = 0; t < nf; t++) {
        for (int e = 0; e < 3; e++) {
            int32_t a = f[t*3+e], b = f[t*3+(e+1)%3];
            int32_t u = a < b ? a : b, v = a < b ? b : a;
            keys[t*3+(size_t)e] = ((int64_t)u << 32) | (int64_t)(uint32_t)v;
        }
    }
    qsort(keys, ne, sizeof(int64_t), i64_cmp);
    for (i = 0; i < ne; ) {
        size_t j = i + 1;
        while (j < ne && keys[j] == keys[i]) j++;
        if (j - i == 1) bnd[n++] = keys[i];
        i = j;
    }
    free(keys);
    *out_n = n;
    return bnd;
}

static int max_edge_run(const int32_t *f, size_t nf)
{
    size_t ne = nf * 3, i;
    int maxr = 0;
    int64_t *keys = (int64_t *)malloc((ne ? ne : 1) * sizeof(int64_t));
    for (size_t t = 0; t < nf; t++) {
        for (int e = 0; e < 3; e++) {
            int32_t a = f[t*3+e], b = f[t*3+(e+1)%3];
            int32_t u = a < b ? a : b, v = a < b ? b : a;
            keys[t*3+(size_t)e] = ((int64_t)u << 32) | (int64_t)(uint32_t)v;
        }
    }
    qsort(keys, ne, sizeof(int64_t), i64_cmp);
    for (i = 0; i < ne; ) {
        size_t j = i + 1;
        while (j < ne && keys[j] == keys[i]) j++;
        if ((int)(j - i) > maxr) maxr = (int)(j - i);
        i = j;
    }
    free(keys);
    return maxr;
}

/* Connected components over faces (shared-vertex adjacency via union-find). */
static size_t component_count(const int32_t *f, size_t nf, size_t nv)
{
    int32_t *par = (int32_t *)malloc((nv ? nv : 1) * sizeof(int32_t));
    uint8_t *used = (uint8_t *)calloc(nv ? nv : 1, 1);
    size_t i, n = 0;
    for (i = 0; i < nv; i++) par[i] = (int32_t)i;
    for (size_t t = 0; t < nf; t++) {
        for (int e = 0; e < 3; e++) {
            int32_t a = f[t*3+e], b = f[t*3+(e+1)%3];
            int32_t ra = a, rb = b;
            used[a] = used[b] = 1;
            while (par[ra] != ra) ra = par[ra];
            while (par[rb] != rb) rb = par[rb];
            if (ra != rb) par[ra] = rb;
        }
    }
    for (i = 0; i < nv; i++) {
        if (!used[i]) continue;
        int32_t r = (int32_t)i;
        while (par[r] != r) r = par[r];
        if (r == (int32_t)i) n++;
    }
    free(par); free(used);
    return n;
}

int main(void)
{
    Arena_T a = Arena_new();
    const SeamPlane plane = { 2, 128.0 };   /* axis x, coord 128 */
    int fail_at = 0;

    /* ---- Scene A: fine 1.5-vox grid patch, x in [120,135], y in [0,36] ---- */
    {
        const int NY = 25, NX = 11;         /* verts per side */
        const float SP = 1.5f;
        size_t nv = (size_t)NY * NX, nf = (size_t)(NY-1) * (NX-1) * 2;
        float *V = (float *)ARENA_ALLOC(a, (long)(nv*3*sizeof(float)));
        int32_t *F = (int32_t *)ARENA_ALLOC(a, (long)(nf*3*sizeof(int32_t)));
        for (int i = 0; i < NY; i++) for (int j = 0; j < NX; j++) {
            size_t id = (size_t)i*NX + (size_t)j;
            V[id*3+0] = 50.0f;
            V[id*3+1] = (float)i * SP;
            V[id*3+2] = 120.0f + (float)j * SP;
        }
        size_t fi = 0;
        for (int i = 0; i+1 < NY; i++) for (int j = 0; j+1 < NX; j++) {
            int32_t v00 = (int32_t)((size_t)i*NX+(size_t)j);
            int32_t v01 = v00 + 1;
            int32_t v10 = (int32_t)((size_t)(i+1)*NX+(size_t)j);
            int32_t v11 = v10 + 1;
            F[fi*3+0]=v00; F[fi*3+1]=v10; F[fi*3+2]=v11; fi++;
            F[fi*3+0]=v00; F[fi*3+1]=v11; F[fi*3+2]=v01; fi++;
        }

        size_t nb_before, nb_after;
        int64_t *bnd_before = boundary_edges(F, nf, &nb_before);

        /* faces with NO vert in band (band 6 of x=128 => x in [122,134]) */
        size_t out_faces_before = 0;
        int64_t out_sig_before = 0;
        for (size_t t = 0; t < nf; t++) {
            int allout = 1;
            for (int e = 0; e < 3; e++) {
                float x = V[(size_t)F[t*3+e]*3+2];
                if (x >= 122.0f && x <= 134.0f) { allout = 0; break; }
            }
            if (allout) { out_faces_before++;
                out_sig_before += F[t*3+0]*3 + F[t*3+1]*5 + F[t*3+2]*7; }
        }

        ComponentMesh cm; memset(&cm, 0, sizeof cm);
        cm.verts = V; cm.faces = F; cm.nv = nv; cm.nf = nf;
        cm.comp_id = 1; cm.self = &cm;

        /* Explicit params: this scene asserts the band-restriction MECHANISM
         * (in-band collapses, out-of-band untouched), so pin band=6 regardless
         * of the tuned default (which widened to cover the graded rim+ramp). */
        WeldRecoarsenParams ap; WeldCleanup_default_recoarsen_params(&ap);
        ap.band = 6.0; ap.collapse_below = 3.0;
        WeldRecoarsenStats st;
        int rc = WeldCleanup_recoarsen_seam(a, &cm, &plane, 1, &ap, &st);
        CHECK(rc == 0, "A: rc=%d", rc);
        CHECK(st.faces_in == nf && st.faces_out == cm.nf && cm.nf < nf,
              "A: faces %zu -> %zu (monotone decrease expected)", nf, cm.nf);
        CHECK(st.n_collapses > 0, "A: no collapses happened");

        /* 1. boundary-edge set bit-identical */
        int64_t *bnd_after = boundary_edges(cm.faces, cm.nf, &nb_after);
        CHECK(nb_before == nb_after, "A: boundary edge count %zu -> %zu",
              nb_before, nb_after);
        if (nb_before == nb_after)
            CHECK(memcmp(bnd_before, bnd_after, nb_before*sizeof(int64_t)) == 0,
                  "A: boundary edge SET changed");
        free(bnd_before); free(bnd_after);

        /* 2. in-band coarsening: interior in-band face count dropped >= 40% */
        {
            size_t inband_after = 0, inband_before = 0;
            for (size_t t = 0; t < nf; t++) {
                int has = 0;
                for (int e = 0; e < 3; e++) {
                    float x = V[(size_t)F[t*3+e]*3+2];
                    if (x >= 122.0f && x <= 134.0f) { has = 1; break; }
                }
                if (has) inband_before++;
            }
            for (size_t t = 0; t < cm.nf; t++) {
                int has = 0;
                for (int e = 0; e < 3; e++) {
                    float x = cm.verts[(size_t)cm.faces[t*3+e]*3+2];
                    if (x >= 122.0f && x <= 134.0f) { has = 1; break; }
                }
                if (has) inband_after++;
            }
            CHECK(inband_after * 10 <= inband_before * 6,
                  "A: in-band faces %zu -> %zu (< 40%% reduction)",
                  inband_before, inband_after);
            fprintf(stderr, "  A: in-band faces %zu -> %zu, total %zu -> %zu, "
                    "%zu collapses %zu flips\n",
                    inband_before, inband_after, nf, cm.nf,
                    st.n_collapses, st.n_flips);
        }

        /* 3. faces with no in-band vert untouched (index-identical set) */
        {
            size_t out_faces_after = 0;
            int64_t out_sig_after = 0;
            for (size_t t = 0; t < cm.nf; t++) {
                int allout = 1;
                for (int e = 0; e < 3; e++) {
                    float x = cm.verts[(size_t)cm.faces[t*3+e]*3+2];
                    if (x >= 122.0f && x <= 134.0f) { allout = 0; break; }
                }
                if (allout) { out_faces_after++;
                    out_sig_after += cm.faces[t*3+0]*3 + cm.faces[t*3+1]*5 + cm.faces[t*3+2]*7; }
            }
            CHECK(out_faces_after == out_faces_before &&
                  out_sig_after == out_sig_before,
                  "A: out-of-band faces changed (%zu/%lld -> %zu/%lld)",
                  out_faces_before, (long long)out_sig_before,
                  out_faces_after, (long long)out_sig_after);
        }

        /* 4. one component, edge-manifold, no degenerate faces */
        CHECK(component_count(cm.faces, cm.nf, cm.nv) == 1, "A: component split");
        CHECK(max_edge_run(cm.faces, cm.nf) <= 2, "A: non-manifold edge");
        for (size_t t = 0; t < cm.nf; t++) {
            int32_t x = cm.faces[t*3+0], y = cm.faces[t*3+1], z = cm.faces[t*3+2];
            CHECK(x != y && y != z && x != z, "A: degenerate face %zu", t);
        }
        fprintf(stderr, "Scene A (band collapse + boundary preservation): %s\n",
                g_fail > fail_at ? "FAIL" : "ok");
        fail_at = g_fail;
    }

    /* ---- Scene B: no planes -> exact no-op ---- */
    {
        float V[12] = { 50,0,126,  50,2,126,  50,0,129,  50,2,129 };
        int32_t F[6] = { 0,1,2,  2,1,3 };
        int32_t F0[6]; memcpy(F0, F, sizeof F0);
        ComponentMesh cm; memset(&cm, 0, sizeof cm);
        cm.verts = V; cm.faces = F; cm.nv = 4; cm.nf = 2;
        cm.comp_id = 1; cm.self = &cm;
        WeldRecoarsenStats st;
        int rc = WeldCleanup_recoarsen_seam(a, &cm, NULL, 0, NULL, &st);
        CHECK(rc == 0 && cm.nf == 2 && st.n_collapses == 0,
              "B: no-op violated (nf=%zu coll=%zu)", cm.nf, st.n_collapses);
        CHECK(memcmp(cm.faces, F0, sizeof F0) == 0, "B: faces mutated");
        fprintf(stderr, "Scene B (no planes -> no-op): %s\n",
                g_fail > fail_at ? "FAIL" : "ok");
        fail_at = g_fail;
    }

    /* ---- Scene C: max_collapse_len cap binds even with a big collapse_below.
     * Two triangles in-band with all edges ~5.5 vox: candidate by length
     * (below=6.0) but over the 5.0 cap -> zero collapses. ---- */
    {
        float V[12] = { 50,0,125,  50,5.5f,125,  50,2.75f,130,  50,8.25f,130 };
        int32_t F[6] = { 0,1,2,  2,1,3 };
        ComponentMesh cm; memset(&cm, 0, sizeof cm);
        cm.verts = V; cm.faces = F; cm.nv = 4; cm.nf = 2;
        cm.comp_id = 1; cm.self = &cm;
        WeldRecoarsenParams p; WeldCleanup_default_recoarsen_params(&p);
        p.collapse_below = 6.0;   /* > cap 5.0 */
        WeldRecoarsenStats st;
        int rc = WeldCleanup_recoarsen_seam(a, &cm, &plane, 1, &p, &st);
        CHECK(rc == 0 && st.n_collapses == 0 && cm.nf == 2,
              "C: cap violated (coll=%zu nf=%zu)", st.n_collapses, cm.nf);
        fprintf(stderr, "Scene C (collapse cap binds): %s\n",
                g_fail > fail_at ? "FAIL" : "ok");
    }

    Arena_dispose(&a);
    if (g_fail) { fprintf(stderr, "WELD_RECOARSEN TESTS: %d FAILURE(S)\n", g_fail); return 1; }
    fprintf(stderr, "ALL WELD_RECOARSEN TESTS PASSED\n");
    return 0;
}
