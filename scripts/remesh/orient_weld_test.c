/* ============================================================================
 * orient_weld_test -- unit tests for OrientWeld_components(_axis), focused on
 * the RADIAL fallback for detached components (2026-07-10 orientation audit:
 * detached shells kept the per-cube (1,1,1) anchor's azimuth-dependent sign
 * because the spatial vote has no ≤radius pairs to see).
 *
 * Geometry: cylinder-arc sheets around axis z through (y,x)=(0,0), built with
 * outward winding by construction; flipping reverses the stored winding only.
 * ==========================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#include "../../src/common/arena.h"
#include "../../src/remesh/orient_weld.h"

#define CHECK(cond, msg) \
    do { if (!(cond)) { fprintf(stderr, "  FAIL: %s\n", (msg)); nfail++; } \
         else { fprintf(stderr, "  ok: %s\n", (msg)); } } while (0)

static void *xgrow(void *p, size_t need, size_t *cap, size_t elem)
{
    size_t ncap = *cap;
    if (need <= ncap) return p;
    ncap = (ncap == 0) ? 1024 : ncap;
    while (ncap < need) ncap *= 2;
    p = realloc(p, ncap * elem);
    if (p == NULL) { fprintf(stderr, "OOM\n"); exit(1); }
    *cap = ncap;
    return p;
}

/* Arc sheet: verts (z,y,x) = (zofs + j, rad cos th, rad sin th),
 * j in [0,nz], th in [th0,th1] over nth+1 columns. Outward winding unless
 * flip. Appends; returns first new vertex index via *out_base. */
static void build_arc(float rad, double th0, double th1, int nth, int nz,
                      float zofs, int flip,
                      float **verts, int32_t **faces, size_t *nv, size_t *nf,
                      size_t *cap_v, size_t *cap_f, size_t *out_base)
{
    int i = 0, j = 0;
    size_t base = *nv;
    *out_base = base;
    for (j = 0; j <= nz; j++) {
        for (i = 0; i <= nth; i++) {
            double th = th0 + (th1 - th0) * (double)i / (double)nth;
            *verts = (float *)xgrow(*verts, (*nv + 1) * 3, cap_v, sizeof(float));
            (*verts)[*nv * 3 + 0] = zofs + (float)j;
            (*verts)[*nv * 3 + 1] = rad * (float)cos(th);
            (*verts)[*nv * 3 + 2] = rad * (float)sin(th);
            (*nv)++;
        }
    }
    for (j = 0; j < nz; j++) {
        for (i = 0; i < nth; i++) {
            int32_t a = (int32_t)(base + (size_t)j * (size_t)(nth + 1) + (size_t)i);
            int32_t b = a + 1;
            int32_t c = a + nth + 1;
            int32_t d = c + 1;
            int32_t t1[3], t2[3];
            t1[0] = a; t1[1] = b; t1[2] = c;
            t2[0] = b; t2[1] = d; t2[2] = c;
            if (flip) {
                int32_t t = t1[1]; t1[1] = t1[2]; t1[2] = t;
                t = t2[1]; t2[1] = t2[2]; t2[2] = t;
            }
            *faces = (int32_t *)xgrow(*faces, (*nf + 2) * 3, cap_f, sizeof(int32_t));
            memcpy(&(*faces)[*nf * 3], t1, sizeof t1); (*nf)++;
            memcpy(&(*faces)[*nf * 3], t2, sizeof t2); (*nf)++;
        }
    }
}

/* Area-weighted sum of face-normal . r_hat over faces whose first vertex is
 * in [v0, v1) -- positive = outward-wound (axis z through y=x=0). */
static double radial_sign_of(const float *verts, const int32_t *faces,
                             size_t nf, size_t v0, size_t v1)
{
    double s = 0.0;
    size_t f = 0;
    for (f = 0; f < nf; f++) {
        size_t a = (size_t)faces[f * 3 + 0];
        size_t b = (size_t)faces[f * 3 + 1];
        size_t c = (size_t)faces[f * 3 + 2];
        double e1[3], e2[3], gn[3], cen_y = 0.0, cen_x = 0.0, rl = 0.0;
        int k = 0;
        if (a < v0 || a >= v1) continue;
        for (k = 0; k < 3; k++) {
            e1[k] = (double)verts[b * 3 + k] - (double)verts[a * 3 + k];
            e2[k] = (double)verts[c * 3 + k] - (double)verts[a * 3 + k];
        }
        gn[0] = e1[1] * e2[2] - e1[2] * e2[1];
        gn[1] = e1[2] * e2[0] - e1[0] * e2[2];
        gn[2] = e1[0] * e2[1] - e1[1] * e2[0];
        cen_y = ((double)verts[a * 3 + 1] + (double)verts[b * 3 + 1]
                 + (double)verts[c * 3 + 1]) / 3.0;
        cen_x = ((double)verts[a * 3 + 2] + (double)verts[b * 3 + 2]
                 + (double)verts[c * 3 + 2]) / 3.0;
        rl = sqrt(cen_y * cen_y + cen_x * cen_x);
        if (rl < 1.0) continue;
        s += (gn[1] * cen_y + gn[2] * cen_x) / rl;   /* z comp of r_hat is 0 */
    }
    return s;
}

int main(void)
{
    int nfail = 0;
    const float AXP[3] = { 0.0f, 0.0f, 0.0f };
    const float AXD[3] = { 1.0f, 0.0f, 0.0f };
    Arena_T arena = Arena_new();

    fprintf(stderr, "[selftest] orient_weld\n");

    /* t1: detached flipped arc + axis -> radial fallback flips it */
    {
        float *v = NULL; int32_t *fc = NULL;
        size_t nv = 0, nf = 0, cv = 0, cf = 0, baseA = 0, baseB = 0;
        size_t flipped = 0, radial = 0;
        build_arc(50.0f, 0.2, 2.2, 24, 8, 0.0f, 0, &v, &fc, &nv, &nf,
                  &cv, &cf, &baseA);                      /* main, outward */
        build_arc(50.0f, 3.4, 4.4, 10, 4, 200.0f, 1, &v, &fc, &nv, &nf,
                  &cv, &cf, &baseB);                      /* far, FLIPPED */
        CHECK(radial_sign_of(v, fc, nf, baseB, nv) < 0.0,
              "t1 setup: far comp wound inward");
        CHECK(OrientWeld_components_axis(arena, v, nv, fc, nf, 3.0f,
                                         AXP, AXD, &flipped, &radial) == 0,
              "t1 pass runs");
        CHECK(flipped == 1 && radial == 1,
              "t1 one flip, decided by radial fallback");
        CHECK(radial_sign_of(v, fc, nf, baseB, nv) > 0.0,
              "t1 far comp now outward (matches anchor)");
        CHECK(radial_sign_of(v, fc, nf, baseA, baseB) > 0.0,
              "t1 anchor untouched");
        free(v); free(fc);
    }

    /* t2: same setup, NO axis (legacy API) -> untouched (documents the bug
     * the fallback fixes) */
    {
        float *v = NULL; int32_t *fc = NULL;
        size_t nv = 0, nf = 0, cv = 0, cf = 0, baseA = 0, baseB = 0;
        size_t flipped = 0;
        build_arc(50.0f, 0.2, 2.2, 24, 8, 0.0f, 0, &v, &fc, &nv, &nf,
                  &cv, &cf, &baseA);
        build_arc(50.0f, 3.4, 4.4, 10, 4, 200.0f, 1, &v, &fc, &nv, &nf,
                  &cv, &cf, &baseB);
        OrientWeld_components(arena, v, nv, fc, nf, 3.0f, &flipped);
        CHECK(flipped == 0 && radial_sign_of(v, fc, nf, baseB, nv) < 0.0,
              "t2 no axis -> detached comp left backward (legacy)");
        free(v); free(fc);
    }

    /* t3: CLOSE flipped arc (2 vox radial gap) -> spatial vote flips it,
     * radial fallback NOT consulted even with axis */
    {
        float *v = NULL; int32_t *fc = NULL;
        size_t nv = 0, nf = 0, cv = 0, cf = 0, baseA = 0, baseB = 0;
        size_t flipped = 0, radial = 0;
        build_arc(50.0f, 0.2, 2.2, 24, 8, 0.0f, 0, &v, &fc, &nv, &nf,
                  &cv, &cf, &baseA);
        build_arc(52.0f, 0.3, 2.1, 20, 8, 0.0f, 1, &v, &fc, &nv, &nf,
                  &cv, &cf, &baseB);                      /* 2 vox out, FLIPPED */
        OrientWeld_components_axis(arena, v, nv, fc, nf, 3.0f,
                                   AXP, AXD, &flipped, &radial);
        CHECK(flipped == 1 && radial == 0,
              "t3 spatial vote flips near comp (fallback not used)");
        CHECK(radial_sign_of(v, fc, nf, baseB, nv) > 0.0,
              "t3 near comp now outward");
        free(v); free(fc);
    }

    /* t4: detached arc already consistent + axis -> counted radial, no flip */
    {
        float *v = NULL; int32_t *fc = NULL;
        size_t nv = 0, nf = 0, cv = 0, cf = 0, baseA = 0, baseB = 0;
        size_t flipped = 0, radial = 0;
        build_arc(50.0f, 0.2, 2.2, 24, 8, 0.0f, 0, &v, &fc, &nv, &nf,
                  &cv, &cf, &baseA);
        build_arc(50.0f, 3.4, 4.4, 10, 4, 200.0f, 0, &v, &fc, &nv, &nf,
                  &cv, &cf, &baseB);                      /* far, consistent */
        OrientWeld_components_axis(arena, v, nv, fc, nf, 3.0f,
                                   AXP, AXD, &flipped, &radial);
        CHECK(flipped == 0 && radial == 1,
              "t4 consistent detached comp: radial-decided, no flip");
        free(v); free(fc);
    }

    /* t6: WEAK spatial evidence must not pre-empt the radial vote. A large
     * flipped far component gets a tiny "contact sliver": one triangle whose
     * two near verts sit 2 vox off the anchor with normals AGREEING with it
     * (a glancing touch voting "keep"). pairs ~20 << 100 -> radial override
     * flips the component anyway (the rank-8 regression from the 4x5x5). */
    {
        float *v = NULL; int32_t *fc = NULL;
        size_t nv = 0, nf = 0, cv = 0, cf = 0, baseA = 0, baseB = 0;
        size_t flipped = 0, radial = 0, t1i = 0, t2i = 0;
        build_arc(50.0f, 0.2, 2.2, 200, 50, 0.0f, 0, &v, &fc, &nv, &nf,
                  &cv, &cf, &baseA);              /* anchor (20k f), outward */
        build_arc(50.0f, 3.4, 5.4, 150, 30, 200.0f, 1, &v, &fc, &nv, &nf,
                  &cv, &cf, &baseB);              /* far (9k f), FLIPPED */
        /* contact sliver: two verts near the anchor + one far B vert */
        t1i = nv; t2i = nv + 1;
        v = (float *)xgrow(v, (nv + 2) * 3, &cv, sizeof(float));
        v[t1i*3+0] = 4.0f; v[t1i*3+1] = 52.0f*(float)cos(1.0);
        v[t1i*3+2] = 52.0f*(float)sin(1.0);
        v[t2i*3+0] = 5.0f; v[t2i*3+1] = v[t1i*3+1]; v[t2i*3+2] = v[t1i*3+2];
        nv += 2;
        fc = (int32_t *)xgrow(fc, (nf + 1) * 3, &cf, sizeof(int32_t));
        fc[nf*3+0] = (int32_t)baseB;      /* B's first vert (far away) */
        fc[nf*3+1] = (int32_t)t2i;        /* winding chosen so the sliver */
        fc[nf*3+2] = (int32_t)t1i;        /* normal AGREES with the anchor */
        nf++;
        CHECK(radial_sign_of(v, fc, nf, baseB, nv) < 0.0,
              "t6 setup: big far comp wound inward");
        OrientWeld_components_axis(arena, v, nv, fc, nf, 3.0f,
                                   AXP, AXD, &flipped, &radial);
        CHECK(flipped == 1 && radial == 1,
              "t6 weak glancing contact -> radial override, one flip");
        CHECK(radial_sign_of(v, fc, nf, baseB, nv) > 0.0,
              "t6 big comp now outward despite agreeing sliver");
        free(v); free(fc);
    }

    /* t5: degenerate inputs */
    {
        size_t flipped = 7, radial = 7;
        CHECK(OrientWeld_components_axis(arena, NULL, 0, NULL, 0, 3.0f,
                                         AXP, AXD, &flipped, &radial) == 0
              && flipped == 0 && radial == 0,
              "t5 empty input clean");
    }

    Arena_dispose(&arena);
    fprintf(stderr, "[selftest] %s (%d failures)\n",
            nfail == 0 ? "ALL PASS" : "FAILURES", nfail);
    return nfail == 0 ? 0 : 1;
}
