/* normal_orient_test.c — unit tests for NormalOrient_consistent (Hoppe 1992 §3.3,
 * src/remesh/normal_orient.c).
 *
 * Builds synthetic flat sheets (coords z,y,x; normals along +/-z) with normals
 * deliberately sign-flipped, and checks:
 *   case 1  a half-flipped sheet becomes internally CONSISTENT, flipping only the
 *           minority, and the majority orientation is preserved;
 *   case 2  an ALREADY-consistent field is left UNCHANGED (0 flips) — the
 *           regression guard for the old +z-seed forcing that globally negated
 *           an already-correct -z field and perturbed BPA;
 *   case 3  two spatially-separated sheets are each oriented independently.
 *
 * Exit 0 = all pass, 1 = a failure.
 */
#include "remesh/normal_orient.h"
#include "common/arena.h"

#include <stdio.h>
#include <stdlib.h>

static int g_fails = 0;
#define CHECK(cond, msg) do { \
    if (!(cond)) { printf("  FAIL: %s\n", (msg)); g_fails++; } \
    else         { printf("  ok:   %s\n", (msg)); } \
} while (0)

/* A flat z=const sheet is consistent iff all its z-normals share one sign. */
static int uniform_zsign(const float *N, size_t lo, size_t hi)
{
    int pos = 0, neg = 0;
    for (size_t i = lo; i < hi; i++) {
        if (N[i*3+0] > 0.0f) pos++;
        else if (N[i*3+0] < 0.0f) neg++;
    }
    return (pos == 0 || neg == 0);
}

int main(void)
{
    printf("=== normal_orient_test ===\n");

    /* ---- case 1: flat sheet, columns x>=12 flipped (160 of 400 = minority) ---- */
    {
        const int GY = 20, GX = 20;
        size_t nv = (size_t)GY * GX;
        float *V = (float *)malloc(nv * 3 * sizeof(float));
        float *N = (float *)malloc(nv * 3 * sizeof(float));
        for (int iy = 0; iy < GY; iy++) for (int ix = 0; ix < GX; ix++) {
            size_t i = (size_t)iy * GX + ix;
            V[i*3+0] = 0.0f; V[i*3+1] = (float)iy; V[i*3+2] = (float)ix;
            float s = (ix >= 12) ? -1.0f : 1.0f;          /* flip the x>=12 minority */
            N[i*3+0] = s; N[i*3+1] = 0.0f; N[i*3+2] = 0.0f;
        }
        Arena_T a = Arena_new();
        size_t flips = 0;
        NormalOrient_consistent(a, V, N, nv, 2.0f, &flips);
        Arena_dispose(&a);
        printf("[case1] half-flipped flat sheet: flips=%zu (expect 160)\n", flips);
        CHECK(uniform_zsign(N, 0, nv), "sheet is internally consistent after orient");
        CHECK(flips == 160, "flipped exactly the 160 minority normals");
        CHECK(N[0] > 0.0f, "majority orientation (+z) preserved");
        free(V); free(N);
    }

    /* ---- case 2: already-consistent all -z -> UNCHANGED (regression guard) ---- */
    {
        const int GY = 16, GX = 16;
        size_t nv = (size_t)GY * GX;
        float *V = (float *)malloc(nv * 3 * sizeof(float));
        float *N = (float *)malloc(nv * 3 * sizeof(float));
        for (int iy = 0; iy < GY; iy++) for (int ix = 0; ix < GX; ix++) {
            size_t i = (size_t)iy * GX + ix;
            V[i*3+0] = 0.0f; V[i*3+1] = (float)iy; V[i*3+2] = (float)ix;
            N[i*3+0] = -1.0f; N[i*3+1] = 0.0f; N[i*3+2] = 0.0f;   /* all -z, consistent */
        }
        Arena_T a = Arena_new();
        size_t flips = 0;
        NormalOrient_consistent(a, V, N, nv, 2.0f, &flips);
        Arena_dispose(&a);
        printf("[case2] already-consistent -z field: flips=%zu (expect 0)\n", flips);
        CHECK(flips == 0, "already-consistent field left UNCHANGED (no global negation)");
        CHECK(N[0] < 0.0f, "kept original -z orientation");
        free(V); free(N);
    }

    /* ---- case 3: two sheets 10 vox apart, each half-flipped -> each consistent ---- */
    {
        const int GY = 10, GX = 10;
        size_t per = (size_t)GY * GX, nv = 2 * per;
        float *V = (float *)malloc(nv * 3 * sizeof(float));
        float *N = (float *)malloc(nv * 3 * sizeof(float));
        for (int s = 0; s < 2; s++)
            for (int iy = 0; iy < GY; iy++) for (int ix = 0; ix < GX; ix++) {
                size_t i = (size_t)s * per + (size_t)iy * GX + ix;
                V[i*3+0] = (s == 0) ? 0.0f : 10.0f;       /* z=0 and z=10, > radius apart */
                V[i*3+1] = (float)iy; V[i*3+2] = (float)ix;
                float sign = (ix >= 5) ? -1.0f : 1.0f;
                N[i*3+0] = sign; N[i*3+1] = 0.0f; N[i*3+2] = 0.0f;
            }
        Arena_T a = Arena_new();
        size_t flips = 0;
        NormalOrient_consistent(a, V, N, nv, 2.0f, &flips);
        Arena_dispose(&a);
        printf("[case3] two separated half-flipped sheets: flips=%zu\n", flips);
        CHECK(uniform_zsign(N, 0, per),   "sheet A consistent after orient");
        CHECK(uniform_zsign(N, per, nv),  "sheet B consistent after orient");
        free(V); free(N);
    }

    printf("=== normal_orient_test %s (%d failure%s) ===\n",
           g_fails ? "FAILED" : "PASSED", g_fails, g_fails == 1 ? "" : "s");
    return g_fails ? 1 : 0;
}
