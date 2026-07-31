/*
 * seam_band_cvt_test -- unit tests for SeamBandCvt_process (+ the pinned CVT
 * underneath it).
 *
 * Scene A: one continuous sheet crossing the seam plane x=128, coarse
 * (13-vox) away from the plane and fine + JITTERED (sliver-ridden, the
 * collapse-scar simulacrum) within ~8 vox of it. Asserts: >= 1 patch
 * accepted; every original vertex bit-identical; the GLOBAL boundary-edge
 * set unchanged (no cracks, no new holes, perimeter untouched -- the
 * conformity contract); global edge-manifold; Euler characteristic
 * unchanged; band mean min-angle strictly improved and decent; new-vert
 * src provenance valid.
 * Scene B: no planes -> alias no-op.
 *
 * Verts are (z,y,x).
 */
#include "../../src/remesh/seam_band_cvt.h"
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

static int64_t ekey(int32_t a, int32_t b)
{
    int32_t u = a < b ? a : b, v = a < b ? b : a;
    return ((int64_t)u << 32) | (int64_t)(uint32_t)v;
}

/* sorted run-1 edge keys; caller frees. NOTE: index-based, so it is only
 * comparable across meshes when shared verts keep their indices -- which is
 * exactly the conformity contract under test. */
static size_t bnd_edges(const int32_t *F, size_t nf, int64_t **out)
{
    size_t ne = nf*3, nb = 0, i;
    int64_t *k = (int64_t *)malloc((ne?ne:1)*sizeof(int64_t));
    int64_t *b = (int64_t *)malloc((ne?ne:1)*sizeof(int64_t));
    for (size_t t = 0; t < nf; t++)
        for (int e = 0; e < 3; e++)
            k[t*3+(size_t)e] = ekey(F[t*3+e], F[t*3+(e+1)%3]);
    qsort(k, ne, sizeof(int64_t), i64_cmp);
    for (i = 0; i < ne; ) {
        size_t j = i+1;
        while (j < ne && k[j] == k[i]) j++;
        if (j-i == 1) b[nb++] = k[i];
        i = j;
    }
    free(k);
    *out = b;
    return nb;
}

static int max_run(const int32_t *F, size_t nf)
{
    size_t ne = nf*3, i; int m = 0;
    int64_t *k = (int64_t *)malloc((ne?ne:1)*sizeof(int64_t));
    for (size_t t = 0; t < nf; t++)
        for (int e = 0; e < 3; e++)
            k[t*3+(size_t)e] = ekey(F[t*3+e], F[t*3+(e+1)%3]);
    qsort(k, ne, sizeof(int64_t), i64_cmp);
    for (i = 0; i < ne; ) {
        size_t j = i+1;
        while (j < ne && k[j] == k[i]) j++;
        if ((int)(j-i) > m) m = (int)(j-i);
        i = j;
    }
    free(k);
    return m;
}

static long euler_chi(const int32_t *F, size_t nf, size_t nv)
{
    size_t ne = nf*3, i, uniq = 0, used_v = 0;
    int64_t *k = (int64_t *)malloc((ne?ne:1)*sizeof(int64_t));
    uint8_t *vu = (uint8_t *)calloc(nv?nv:1, 1);
    for (size_t t = 0; t < nf; t++)
        for (int e = 0; e < 3; e++) {
            k[t*3+(size_t)e] = ekey(F[t*3+e], F[t*3+(e+1)%3]);
            vu[F[t*3+e]] = 1;
        }
    qsort(k, ne, sizeof(int64_t), i64_cmp);
    for (i = 0; i < ne; ) {
        size_t j = i+1;
        while (j < ne && k[j] == k[i]) j++;
        uniq++;
        i = j;
    }
    for (i = 0; i < nv; i++) if (vu[i]) used_v++;
    free(k); free(vu);
    return (long)used_v - (long)uniq + (long)nf;
}

static double tri_min_angle_deg(const float *V, int32_t a, int32_t b, int32_t c)
{
    const float *P[3] = { &V[(size_t)a*3], &V[(size_t)b*3], &V[(size_t)c*3] };
    double best = 1e9;
    for (int i = 0; i < 3; i++) {
        const float *p = P[i], *q = P[(i+1)%3], *r = P[(i+2)%3];
        double u[3] = { q[0]-p[0], q[1]-p[1], q[2]-p[2] };
        double w[3] = { r[0]-p[0], r[1]-p[1], r[2]-p[2] };
        double lu = sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);
        double lw = sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]);
        if (lu < 1e-12 || lw < 1e-12) return 0.0;
        double cs = (u[0]*w[0]+u[1]*w[1]+u[2]*w[2])/(lu*lw);
        if (cs > 1.0) cs = 1.0; if (cs < -1.0) cs = -1.0;
        double ang = acos(cs) * 180.0 / 3.14159265358979323846;
        if (ang < best) best = ang;
    }
    return best;
}

static double band_mean_min_angle(const float *V, const int32_t *F, size_t nf,
                                  double plane_x, double band)
{
    double s = 0.0; size_t n = 0;
    for (size_t t = 0; t < nf; t++) {
        int in = 0;
        for (int e = 0; e < 3; e++)
            if (fabs((double)V[(size_t)F[t*3+e]*3+2] - plane_x) <= band) { in = 1; break; }
        if (!in) continue;
        s += tri_min_angle_deg(V, F[t*3+0], F[t*3+1], F[t*3+2]);
        n++;
    }
    return n ? s / (double)n : 0.0;
}

int main(void)
{
    Arena_T a = Arena_new();
    const SeamPlane plane = { 2, 128.0 };

    /* ---- Scene A ---- */
    {
        /* column x positions: coarse away from 128, fine near it */
        const double colx[] = { 76, 89, 102, 115, 121, 124.5, 126.5, 128.5,
                                130.5, 133, 136, 141, 154, 167, 180 };
        const int NC = (int)(sizeof colx / sizeof colx[0]);
        const int NR = 9;                       /* y rows, 13-vox pitch */
        size_t nv = (size_t)NR * NC, nf = (size_t)(NR-1) * (NC-1) * 2;
        float *V = (float *)ARENA_ALLOC(a, (long)(nv*3*sizeof(float)));
        int32_t *F = (int32_t *)ARENA_ALLOC(a, (long)(nf*3*sizeof(int32_t)));
        uint32_t rng = 777u;
        for (int i = 0; i < NR; i++) for (int j = 0; j < NC; j++) {
            size_t id = (size_t)i*NC + (size_t)j;
            double jit = 0.0;
            /* jitter y in the fine zone (interior rows only) to make slivers.
             * Amplitude +-1.0 vox against 2-3.5-vox columns: ugly anisotropy
             * (mean min-angle ~20 deg) but no near-coincident verts -- a tile
             * border through here must still yield usable pinned sites (two
             * pins closer than the RVD can resolve is a LEGITIMATE rejection,
             * covered by the fail-closed gates, not what this scene tests). */
            if (fabs(colx[j] - 128.0) < 9.0 && i > 0 && i < NR-1) {
                rng ^= rng << 13; rng ^= rng >> 17; rng ^= rng << 5;
                jit = ((double)(rng >> 8) * (1.0/16777216.0) - 0.5) * 2.0;
            }
            V[id*3+0] = 50.0f;
            V[id*3+1] = (float)(13.0*i + jit);
            V[id*3+2] = (float)colx[j];
        }
        size_t fi = 0;
        for (int i = 0; i+1 < NR; i++) for (int j = 0; j+1 < NC; j++) {
            int32_t v00 = (int32_t)((size_t)i*NC+(size_t)j), v01 = v00+1;
            int32_t v10 = (int32_t)((size_t)(i+1)*NC+(size_t)j), v11 = v10+1;
            F[fi*3+0]=v00; F[fi*3+1]=v10; F[fi*3+2]=v11; fi++;
            F[fi*3+0]=v00; F[fi*3+1]=v11; F[fi*3+2]=v01; fi++;
        }

        int64_t *bnd_before = NULL;
        size_t nb_before = bnd_edges(F, nf, &bnd_before);
        long chi_before = euler_chi(F, nf, nv);
        double q_before = band_mean_min_angle(V, F, nf, 128.0, 14.0);

        SeamBandCvtParams p; SeamBandCvt_default_params(&p);
        p.min_faces = 32; p.target_h = 10.0;
        SeamBandCvtStats st;
        float *OV; int32_t *OF; int32_t *SRC;
        size_t onv, onf, nnew;
        int rc = SeamBandCvt_process(a, V, nv, F, nf, &plane, 1, &p,
                                     &OV, &onv, &OF, &onf, &SRC, &nnew, &st);
        CHECK(rc == 0, "A: rc=%d", rc);
        fprintf(stderr, "  A: patches=%zu acc=%zu rej[rc=%zu pin=%zu bnd=%zu "
                "chi=%zu nm=%zu comp=%zu] small=%zu clean=%zu\n",
                st.patches, st.accepted, st.rej_rc, st.rej_pin, st.rej_bnd,
                st.rej_chi, st.rej_manifold, st.rej_comp,
                st.skipped_small, st.skipped_clean);
        CHECK(st.accepted >= 1, "A: no patch accepted (patches=%zu rejected=%zu)",
              st.patches, st.rejected);

        /* originals never move */
        for (size_t v = 0; v < nv; v++)
            CHECK(OV[v*3+0]==V[v*3+0] && OV[v*3+1]==V[v*3+1] && OV[v*3+2]==V[v*3+2],
                  "A: original vert %zu moved", v);
        for (size_t i = 0; i < nnew; i++)
            CHECK(SRC[i] >= 0 && (size_t)SRC[i] < nv, "A: bad src[%zu]=%d",
                  i, SRC[i]);

        /* global boundary set unchanged; manifold; chi unchanged */
        {
            int64_t *bnd_after = NULL;
            size_t nb_after = bnd_edges(OF, onf, &bnd_after);
            CHECK(nb_after == nb_before, "A: boundary edges %zu -> %zu",
                  nb_before, nb_after);
            if (nb_after == nb_before)
                CHECK(memcmp(bnd_before, bnd_after, nb_before*sizeof(int64_t)) == 0,
                      "A: boundary edge SET changed (crack or new hole)");
            free(bnd_after);
        }
        CHECK(max_run(OF, onf) <= 2, "A: non-manifold edge");
        CHECK(euler_chi(OF, onf, onv) == chi_before, "A: Euler chi %ld -> %ld",
              chi_before, euler_chi(OF, onf, onv));

        /* quality: band mean min-angle strictly improved and respectable */
        {
            double q_after = band_mean_min_angle(OV, OF, onf, 128.0, 14.0);
            CHECK(q_after > q_before, "A: band quality not improved (%.1f -> %.1f)",
                  q_before, q_after);
            CHECK(q_after >= 30.0, "A: band quality still poor (%.1f deg)", q_after);
            fprintf(stderr, "  A: band mean min-angle %.1f -> %.1f deg, "
                    "faces %zu -> %zu, +%zu verts, %zu/%zu patches\n",
                    q_before, q_after, st.faces_band_in, st.faces_band_out,
                    st.verts_added, st.accepted, st.patches);
        }
        fprintf(stderr, "Scene A (jittered band -> CVT): %s\n", g_fail ? "FAIL" : "ok");
    }

    /* ---- Scene B: no planes -> alias no-op ---- */
    {
        float V[9] = { 50,0,120, 50,13,120, 50,0,133 };
        int32_t F[3] = { 0,1,2 };
        float *OV; int32_t *OF; int32_t *SRC; size_t onv, onf, nnew;
        int rc = SeamBandCvt_process(a, V, 3, F, 1, NULL, 0, NULL,
                                     &OV, &onv, &OF, &onf, &SRC, &nnew, NULL);
        CHECK(rc == 0 && OV == V && OF == F && onv == 3 && onf == 1 && nnew == 0,
              "B: not an alias no-op");
        fprintf(stderr, "Scene B (no planes -> no-op): ok\n");
    }

    Arena_dispose(&a);
    if (g_fail) { fprintf(stderr, "SEAM_BAND_CVT TESTS: %d FAILURE(S)\n", g_fail); return 1; }
    fprintf(stderr, "ALL SEAM_BAND_CVT TESTS PASSED\n");
    return 0;
}
