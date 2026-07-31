/*
 * mesh_trim_cut_test -- unit tests for Mesh_trim_cut_to_owned_box
 * (cut-at-plane trim). Synthetic triangles against the box [1, 127]^3,
 * eps = TRIM_CUT_SNAP_EPS_VOX. Verts are (z, y, x).
 *
 * Covers: all-inside passthrough, fully-outside drop, 1-out and 2-out clips
 * with exact area conservation, cut verts bitwise ON the plane, shared-edge
 * cut vertex dedup (watertight, run-2 interior edge), pass-0 jitter snap
 * (the CVT double->float ring-loss case), on-plane endpoint reuse (no
 * duplicate corners), 2-plane and 3-plane corner clips, winding preservation,
 * empty input.
 */
#include "../../src/common/mesh_trim.h"
#include "../../src/common/pipeline_constants.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int g_fail = 0;
#define CHECK(cond, ...) do { \
    if (!(cond)) { g_fail++; fprintf(stderr, "  FAIL: " __VA_ARGS__); fprintf(stderr, "\n"); } \
} while (0)

static const float LO = 1.0f, HI = 127.0f;
static const float EPS = TRIM_CUT_SNAP_EPS_VOX;

static double tri_area(const float *v, const int32_t *f, size_t t)
{
    const float *a = &v[(size_t)f[t*3+0]*3];
    const float *b = &v[(size_t)f[t*3+1]*3];
    const float *c = &v[(size_t)f[t*3+2]*3];
    double e1[3] = { b[0]-a[0], b[1]-a[1], b[2]-a[2] };
    double e2[3] = { c[0]-a[0], c[1]-a[1], c[2]-a[2] };
    double cr0 = e1[1]*e2[2] - e1[2]*e2[1];
    double cr1 = e1[2]*e2[0] - e1[0]*e2[2];
    double cr2 = e1[0]*e2[1] - e1[1]*e2[0];
    return 0.5 * sqrt(cr0*cr0 + cr1*cr1 + cr2*cr2);
}

static double total_area(const float *v, const int32_t *f, size_t nf)
{
    double s = 0.0;
    for (size_t t = 0; t < nf; t++) s += tri_area(v, f, t);
    return s;
}

/* z-component sign of the face normal (planar z=const tests): >0 = CCW in (y,x)
 * with our (z,y,x) layout and cross(e1,e2) as in tri_area. */
static double face_nz(const float *v, const int32_t *f, size_t t)
{
    const float *a = &v[(size_t)f[t*3+0]*3];
    const float *b = &v[(size_t)f[t*3+1]*3];
    const float *c = &v[(size_t)f[t*3+2]*3];
    double e1[3] = { b[0]-a[0], b[1]-a[1], b[2]-a[2] };
    double e2[3] = { c[0]-a[0], c[1]-a[1], c[2]-a[2] };
    return e1[1]*e2[2] - e1[2]*e2[1];
}

static int i64_cmp(const void *pa, const void *pb)
{
    int64_t a = *(const int64_t *)pa, b = *(const int64_t *)pb;
    return a < b ? -1 : (a > b ? 1 : 0);
}

/* max undirected-edge multiplicity + run count of a specific edge */
static void edge_audit(const int32_t *f, size_t nf, size_t nv,
                       int *out_max_run, int32_t qu, int32_t qv, int *out_q_run)
{
    (void)nv;
    size_t ne = nf * 3;
    int64_t *keys = (int64_t *)malloc((ne ? ne : 1) * sizeof(int64_t));
    for (size_t t = 0; t < nf; t++) {
        for (int e = 0; e < 3; e++) {
            int32_t a = f[t*3+e], b = f[t*3+(e+1)%3];
            int32_t u = a < b ? a : b, v = a < b ? b : a;
            keys[t*3+(size_t)e] = ((int64_t)u << 32) | (int64_t)(uint32_t)v;
        }
    }
    qsort(keys, ne, sizeof(int64_t), i64_cmp);
    int max_run = 0, q_run = 0;
    int64_t qkey = qu >= 0
        ? ((int64_t)(qu < qv ? qu : qv) << 32) | (int64_t)(uint32_t)(qu < qv ? qv : qu)
        : -1;
    for (size_t i = 0; i < ne; ) {
        size_t j = i + 1;
        while (j < ne && keys[j] == keys[i]) j++;
        int run = (int)(j - i);
        if (run > max_run) max_run = run;
        if (qkey >= 0 && keys[i] == qkey) q_run = run;
        i = j;
    }
    free(keys);
    *out_max_run = max_run;
    if (out_q_run) *out_q_run = q_run;
}

static size_t count_verts_at(const float *v, size_t nv, float z, float y, float x)
{
    size_t n = 0;
    for (size_t i = 0; i < nv; i++) {
        if (fabsf(v[i*3+0]-z) < 1e-6f && fabsf(v[i*3+1]-y) < 1e-6f &&
            fabsf(v[i*3+2]-x) < 1e-6f) n++;
    }
    return n;
}

static void expect_inside_box(const float *v, size_t nv, const char *tag)
{
    for (size_t i = 0; i < nv; i++) {
        for (int a = 0; a < 3; a++) {
            CHECK(v[i*3+(size_t)a] >= LO - 1e-6f && v[i*3+(size_t)a] <= HI + 1e-6f,
                  "%s: vert %zu axis %d = %g outside [%g,%g]",
                  tag, i, a, (double)v[i*3+(size_t)a], (double)LO, (double)HI);
        }
    }
}

int main(void)
{
    Arena_T a = Arena_new();
    float *ov; size_t onv; int32_t *of; size_t onf; size_t ncut;
    int rc;

    /* ---- T1: all inside -> passthrough, 0 cut ---- */
    {
        float V[9] = { 5,10,10,  5,20,10,  5,10,20 };
        int32_t F[3] = { 0,1,2 };
        rc = Mesh_trim_cut_to_owned_box(a, V, 3, F, 1, LO, HI, EPS,
                                        &ov, &onv, &of, &onf, &ncut);
        CHECK(rc == 0 && onf == 1 && onv == 3 && ncut == 0,
              "T1 passthrough: rc=%d nf=%zu nv=%zu cut=%zu", rc, onf, onv, ncut);
        CHECK(fabs(total_area(ov, of, onf) - total_area(V, F, 1)) < 1e-9,
              "T1 area changed");
        fprintf(stderr, "T1 all-inside passthrough: %s\n", g_fail ? "FAIL" : "ok");
    }

    /* ---- T2: fully outside one plane -> dropped ---- */
    {
        float V[9] = { 5,10,-5,  5,20,-4,  5,10,-3 };   /* x all < LO */
        int32_t F[3] = { 0,1,2 };
        rc = Mesh_trim_cut_to_owned_box(a, V, 3, F, 1, LO, HI, EPS,
                                        &ov, &onv, &of, &onf, &ncut);
        CHECK(rc == 0 && onf == 0 && onv == 0, "T2 drop: rc=%d nf=%zu", rc, onf);
        fprintf(stderr, "T2 fully-outside drop: ok\n");
    }

    /* ---- T3: one vert out (x < LO) -> quad (2 tris); exact area; cut verts
     *          bitwise on plane; winding preserved ---- */
    {
        /* A=(5,10,0) out, B=(5,14,4) in, C=(5,6,4) in.
         * Crossings at x=1: AB -> (5,11,1), CA -> (5,9,1).
         * area(orig)=16, clipped corner=1 -> inside area 15. */
        float V[9] = { 5,10,0,  5,14,4,  5,6,4 };
        int32_t F[3] = { 0,1,2 };
        double nz_in = face_nz(V, F, 0);
        rc = Mesh_trim_cut_to_owned_box(a, V, 3, F, 1, LO, HI, EPS,
                                        &ov, &onv, &of, &onf, &ncut);
        CHECK(rc == 0 && onf == 2 && onv == 4 && ncut == 1,
              "T3 quad: rc=%d nf=%zu nv=%zu cut=%zu", rc, onf, onv, ncut);
        CHECK(fabs(total_area(ov, of, onf) - 15.0) < 1e-6,
              "T3 area %.9f != 15", total_area(ov, of, onf));
        CHECK(count_verts_at(ov, onv, 5, 11, 1) == 1 &&
              count_verts_at(ov, onv, 5, 9, 1) == 1, "T3 cut verts wrong");
        for (size_t i = 0; i < onv; i++)
            if (ov[i*3+2] < 2.0f)   /* the two cut verts */
                CHECK(ov[i*3+2] == LO, "T3 cut vert not bitwise on plane: %g",
                      (double)ov[i*3+2]);
        for (size_t t = 0; t < onf; t++)
            CHECK(face_nz(ov, of, t) * nz_in > 0.0, "T3 winding flipped");
        expect_inside_box(ov, onv, "T3");
        fprintf(stderr, "T3 one-out clip: ok (area 15.0 exact, on-plane cuts)\n");
    }

    /* ---- T4: two verts out -> single smaller triangle ---- */
    {
        /* A=(5,10,4) in, B=(5,14,0) out, C=(5,6,0) out.
         * Crossings x=1: AB t=3/4 -> (5,13,1), CA from C t=1/4 -> (5,7,1).
         * area(orig)=16; inside tri (A, cutAB, cutCA): base 6 at x=1, apex 3
         * away -> area 9. */
        float V[9] = { 5,10,4,  5,14,0,  5,6,0 };
        int32_t F[3] = { 0,1,2 };
        rc = Mesh_trim_cut_to_owned_box(a, V, 3, F, 1, LO, HI, EPS,
                                        &ov, &onv, &of, &onf, &ncut);
        CHECK(rc == 0 && onf == 1 && onv == 3 && ncut == 1,
              "T4: rc=%d nf=%zu nv=%zu", rc, onf, onv);
        CHECK(fabs(total_area(ov, of, onf) - 9.0) < 1e-6,
              "T4 area %.9f != 9", total_area(ov, of, onf));
        fprintf(stderr, "T4 two-out clip: ok (area 9.0 exact)\n");
    }

    /* ---- T5: two faces sharing a crossing edge -> ONE shared cut vert,
     *          interior cut edge has run 2 (watertight) ---- */
    {
        /* Shared edge B-C crosses x=1 at (5,10,1).
         * A=(5,9,4) in, B=(5,8,-1) out, C=(5,12,3) in, E=(5,13,2) in.
         * t1=(A,B,C), t2=(C,B,E): consistent winding (B->C vs C->B). */
        float V[12] = { 5,9,4,  5,8,-1,  5,12,3,  5,13,2 };
        int32_t F[6] = { 0,1,2,  2,1,3 };
        rc = Mesh_trim_cut_to_owned_box(a, V, 4, F, 2, LO, HI, EPS,
                                        &ov, &onv, &of, &onf, &ncut);
        CHECK(rc == 0 && ncut == 2, "T5: rc=%d cut=%zu", rc, ncut);
        CHECK(count_verts_at(ov, onv, 5, 10, 1) == 1,
              "T5 shared cut vert duplicated (%zu copies)",
              count_verts_at(ov, onv, 5, 10, 1));
        /* find ids of the shared cut vert and C, check edge run == 2 */
        {
            int32_t cut_id = -1, c_id = -1;
            for (size_t i = 0; i < onv; i++) {
                if (count_verts_at(&ov[i*3], 1, 5, 10, 1)) cut_id = (int32_t)i;
                if (count_verts_at(&ov[i*3], 1, 5, 12, 3)) c_id = (int32_t)i;
            }
            CHECK(cut_id >= 0 && c_id >= 0, "T5 verts not found");
            int max_run = 0, q_run = 0;
            edge_audit(of, onf, onv, &max_run, cut_id, c_id, &q_run);
            CHECK(max_run <= 2, "T5 non-manifold edge (run %d)", max_run);
            CHECK(q_run == 2, "T5 shared cut edge run %d != 2 (crack!)", q_run);
        }
        fprintf(stderr, "T5 shared-edge cut: ok (1 shared vert, run-2 cut edge)\n");
    }

    /* ---- T6: pass-0 jitter snap -- vert at LO - 1e-5 (the CVT float-rounding
     *          case) is clamped ONTO the plane; whole face kept, 0 cut ---- */
    {
        float V[9] = { 5,10,LO-1e-5f,  5,20,10,  5,10,20 };
        int32_t F[3] = { 0,1,2 };
        rc = Mesh_trim_cut_to_owned_box(a, V, 3, F, 1, LO, HI, EPS,
                                        &ov, &onv, &of, &onf, &ncut);
        CHECK(rc == 0 && onf == 1 && onv == 3 && ncut == 0,
              "T6 snap: rc=%d nf=%zu cut=%zu", rc, onf, ncut);
        CHECK(count_verts_at(ov, onv, 5, 10, LO) == 1, "T6 vert not snapped to plane");
        fprintf(stderr, "T6 jitter snap (ring-loss fix): ok\n");
    }

    /* ---- T7: on-plane endpoint reuse -- A exactly on plane, B out, C in:
     *          no duplicate corner, single clean triangle ---- */
    {
        float V[9] = { 5,14,LO,  5,12,-1,  5,8,3 };   /* A on-plane, B out, C in */
        int32_t F[3] = { 0,1,2 };
        rc = Mesh_trim_cut_to_owned_box(a, V, 3, F, 1, LO, HI, EPS,
                                        &ov, &onv, &of, &onf, &ncut);
        CHECK(rc == 0 && onf == 1 && onv == 3,
              "T7: rc=%d nf=%zu nv=%zu (dup corner?)", rc, onf, onv);
        CHECK(count_verts_at(ov, onv, 5, 10, 1) == 1,   /* B-C cut at t=0.5 */
              "T7 cut vert missing");
        fprintf(stderr, "T7 on-plane endpoint reuse: ok\n");
    }

    /* ---- T8: 2-plane corner clip (y-lo and x-lo) ---- */
    {
        float V[9] = { 5,-2,-2,  5,4,3,  5,3,4 };
        int32_t F[3] = { 0,1,2 };
        double nz_in = face_nz(V, F, 0);
        rc = Mesh_trim_cut_to_owned_box(a, V, 3, F, 1, LO, HI, EPS,
                                        &ov, &onv, &of, &onf, &ncut);
        CHECK(rc == 0 && onf >= 2 && ncut == 1, "T8: rc=%d nf=%zu", rc, onf);
        expect_inside_box(ov, onv, "T8");
        CHECK(total_area(ov, of, onf) > 0.0 &&
              total_area(ov, of, onf) < total_area(V, F, 1), "T8 area bounds");
        for (size_t t = 0; t < onf; t++)
            CHECK(face_nz(ov, of, t) * nz_in > 0.0, "T8 winding flipped");
        fprintf(stderr, "T8 2-plane corner clip: ok (%zu tris)\n", onf);
    }

    /* ---- T9: 3-plane corner clip at the hi corner ---- */
    {
        float V[9] = { 126,126,126,  129,126,125,  126,129,129 };
        int32_t F[3] = { 0,1,2 };
        rc = Mesh_trim_cut_to_owned_box(a, V, 3, F, 1, LO, HI, EPS,
                                        &ov, &onv, &of, &onf, &ncut);
        CHECK(rc == 0 && onf >= 1 && ncut == 1, "T9: rc=%d nf=%zu", rc, onf);
        expect_inside_box(ov, onv, "T9");
        {
            int max_run = 0;
            edge_audit(of, onf, onv, &max_run, -1, -1, NULL);
            CHECK(max_run <= 2, "T9 non-manifold (run %d)", max_run);
        }
        fprintf(stderr, "T9 3-plane corner clip: ok (%zu tris)\n", onf);
    }

    /* ---- T10: empty input ---- */
    {
        rc = Mesh_trim_cut_to_owned_box(a, NULL, 0, NULL, 0, LO, HI, EPS,
                                        &ov, &onv, &of, &onf, &ncut);
        CHECK(rc == 0 && onf == 0 && onv == 0, "T10 empty: rc=%d", rc);
        fprintf(stderr, "T10 empty input: ok\n");
    }

    Arena_dispose(&a);
    if (g_fail) { fprintf(stderr, "MESH_TRIM_CUT TESTS: %d FAILURE(S)\n", g_fail); return 1; }
    fprintf(stderr, "ALL MESH_TRIM_CUT TESTS PASSED\n");
    return 0;
}
