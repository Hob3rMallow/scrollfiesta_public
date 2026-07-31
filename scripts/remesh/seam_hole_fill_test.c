/*
 * seam_hole_fill_test.c -- unit tests for src/remesh/seam_hole_fill.c.
 *
 * Gates (synthetic meshes with known verdicts):
 *   G1  single-triangle seam puncture (triangular annulus, inner hole straddling
 *       X=128) -> filled with exactly 1 triangle, mesh stays manifold.
 *   G2  same hole moved fully to one side of the seam -> NOT a straddle -> not filled.
 *   G3  hole extent above the cap -> rejected (skip_extent), not filled.
 *   G4  phase gate: verts spanning > tol turns about the umbilicus (two wraps) ->
 *       rejected (skip_phase), not filled.
 *   G5  quad seam puncture (square annulus) -> filled with 2 triangles, manifold.
 *   G6  trivial input (empty mesh, lone triangle) -> no crash, nothing filled.
 *   G7  idempotence: a second pass over the G1 result fills nothing.
 *
 * Build: seam_hole_fill_test.vcxproj (links seam_hole_fill.c + arena.c + except.c).
 * Run:   seam_hole_fill_test.exe   (exit 0 = all pass, 1 = a gate failed)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "common/arena.h"
#include "common/mesh_types.h"
#include "remesh/seam_hole_fill.h"

/* ---- manifold check on a flat face list (<=2 faces/edge, no degenerate) ---- */
typedef struct { int32_t v0, v1; } UE;
static int ue_cmp(const void *a, const void *b)
{
    const UE *x = a, *y = b;
    if (x->v0 != y->v0) return x->v0 < y->v0 ? -1 : 1;
    if (x->v1 != y->v1) return x->v1 < y->v1 ? -1 : 1;
    return 0;
}
static int manifold_ok(const int32_t *f, size_t nf)
{
    size_t i, n = nf * 3;
    for (i = 0; i < nf; i++) {
        int32_t a = f[i*3+0], b = f[i*3+1], c = f[i*3+2];
        if (a == b || b == c || a == c) return 0;
    }
    if (nf == 0) return 1;
    UE *e = (UE *)malloc(n * sizeof(UE));
    for (i = 0; i < nf; i++) {
        int32_t t[3] = { f[i*3+0], f[i*3+1], f[i*3+2] };
        for (int k = 0; k < 3; k++) {
            int32_t a = t[k], b = t[(k+1)%3];
            if (a > b) { int32_t s = a; a = b; b = s; }
            e[i*3+(size_t)k].v0 = a; e[i*3+(size_t)k].v1 = b;
        }
    }
    qsort(e, n, sizeof(UE), ue_cmp);
    int ok = 1; i = 0;
    while (i < n) {
        size_t j = i + 1;
        while (j < n && e[j].v0 == e[i].v0 && e[j].v1 == e[i].v1) j++;
        if (j - i > 2) { ok = 0; break; }
        i = j;
    }
    free(e);
    return ok;
}

/* Build a triangular annulus: inner triangle i0,i1,i2 (the hole) surrounded by
 * outer triangle o0,o1,o2. 6 faces, two boundary loops (inner + outer). Verts:
 * [0..2]=inner, [3..5]=outer. */
static void build_tri_annulus(Arena_T ar, ComponentMesh *cm,
                              const float inner[9], const float outer[9])
{
    memset(cm, 0, sizeof *cm);
    cm->nv = 6; cm->nf = 6;
    cm->verts = (float *)ARENA_ALLOC(ar, (long)(cm->nv * 3 * sizeof(float)));
    cm->faces = (int32_t *)ARENA_ALLOC(ar, (long)(cm->nf * 3 * sizeof(int32_t)));
    for (int k = 0; k < 9; k++) { cm->verts[k] = inner[k]; cm->verts[9+k] = outer[k]; }
    /* i=0..2, o=3..5. Annulus faces wound consistently; inner (0,1,2) left open. */
    int32_t F[18] = {
        3,4,0,  4,1,0,     /* edge 0-1 */
        4,5,1,  5,2,1,     /* edge 1-2 */
        5,3,2,  3,0,2      /* edge 2-0 */
    };
    memcpy(cm->faces, F, sizeof F);
    cm->comp_id = 1; cm->self = cm;
}

/* Square annulus: inner quad 0..3 (hole) + outer quad 4..7. 8 faces. */
static void build_quad_annulus(Arena_T ar, ComponentMesh *cm,
                               const float inner[12], const float outer[12])
{
    memset(cm, 0, sizeof *cm);
    cm->nv = 8; cm->nf = 8;
    cm->verts = (float *)ARENA_ALLOC(ar, (long)(cm->nv * 3 * sizeof(float)));
    cm->faces = (int32_t *)ARENA_ALLOC(ar, (long)(cm->nf * 3 * sizeof(int32_t)));
    for (int k = 0; k < 12; k++) { cm->verts[k] = inner[k]; cm->verts[12+k] = outer[k]; }
    int32_t F[24];
    int fi = 0;
    for (int k = 0; k < 4; k++) {
        int in0 = k, in1 = (k+1)%4, ou0 = 4+k, ou1 = 4+(k+1)%4;
        F[fi*3+0]=ou0; F[fi*3+1]=ou1; F[fi*3+2]=in0; fi++;
        F[fi*3+0]=ou1; F[fi*3+1]=in1; F[fi*3+2]=in0; fi++;
    }
    memcpy(cm->faces, F, sizeof F);
    cm->comp_id = 1; cm->self = cm;
}

int main(void)
{
    Arena_T ar = Arena_new();
    int fails = 0;
    SeamHoleFillParams P; SeamHoleFill_default_params(&P);
    SeamHoleFillStats st;

    /* G1: single-triangle seam puncture straddling X=128 (phase gate off). */
    {
        float inner[9] = { 0,0,127.5f,  2,0,128.5f,  0,2,128.0f };
        float outer[9] = { -1,-1,127.0f,  3,-1,129.0f, -1,3,128.0f };
        ComponentMesh cm; build_tri_annulus(ar, &cm, inner, outer);
        size_t nf0 = cm.nf;
        SeamHoleFill_process(ar, &cm, &P, &st);
        int ok = (st.filled == 1 && st.tris_added == 1 && cm.nf == nf0 + 1 &&
                  manifold_ok(cm.faces, cm.nf));
        printf("[G1] tri seam puncture: filled=%zu tris+=%zu nf %zu->%zu manifold=%d -> %s\n",
               st.filled, st.tris_added, nf0, cm.nf, manifold_ok(cm.faces, cm.nf),
               ok ? "ok" : "FAIL");
        if (!ok) fails++;
    }

    /* G2: same hole shifted fully to x>128, and z/y parked near 64 (mid-cube) so
     * no boundary loop straddles any 128-multiple plane -> no candidate, no fill. */
    {
        float inner[9] = { 64,64,130.5f,  66,64,131.5f,  64,66,131.0f };
        float outer[9] = { 63,63,130.0f,  67,63,132.0f,  63,67,131.0f };
        ComponentMesh cm; build_tri_annulus(ar, &cm, inner, outer);
        SeamHoleFill_process(ar, &cm, &P, &st);
        int ok = (st.filled == 0 && st.seam_candidates == 0);
        printf("[G2] non-straddle: filled=%zu seam_cand=%zu -> %s\n",
               st.filled, st.seam_candidates, ok ? "ok" : "FAIL");
        if (!ok) fails++;
    }

    /* G3: hole extent above the (phase-off) cap -> skip_extent, not filled. */
    {
        float inner[9] = { 0,0,126.5f,  4,0,129.5f,  0,4,128.0f };   /* diag ~6.4 */
        float outer[9] = { -1,-1,126.0f,  5,-1,130.0f, -1,5,128.0f };
        ComponentMesh cm; build_tri_annulus(ar, &cm, inner, outer);
        SeamHoleFill_process(ar, &cm, &P, &st);
        int ok = (st.filled == 0 && st.skip_extent >= 1);
        printf("[G3] over-extent: filled=%zu skip_extent=%zu -> %s\n",
               st.filled, st.skip_extent, ok ? "ok" : "FAIL");
        if (!ok) fails++;
    }

    /* G4: phase gate. umbilicus at origin, pitch 1.0; the hole verts span a big
     * radial range (~4 turns) -> different wraps -> skip_phase. */
    {
        SeamHoleFillParams Q = P;
        /* umbilicus off-origin so the phase gate ARMS (it requires umb != 0). */
        Q.umb_y = 0.0; Q.umb_x = 100.0; Q.pitch = 1.0; Q.wind_tol_turns = 0.25;
        float inner[9] = { 0,0,127.0f,  0,0,131.0f,  0,4,129.0f };
        float outer[9] = { -1,-1,126.0f,  1,-1,132.0f, -1,5,129.0f };
        ComponentMesh cm; build_tri_annulus(ar, &cm, inner, outer);
        SeamHoleFill_process(ar, &cm, &Q, &st);
        int ok = (st.filled == 0 && st.skip_phase >= 1);
        printf("[G4] phase gate: filled=%zu skip_phase=%zu -> %s\n",
               st.filled, st.skip_phase, ok ? "ok" : "FAIL");
        if (!ok) fails++;
    }

    /* G5: quad seam puncture -> 2 fill triangles, manifold. Inner quad lies in the
     * tilted plane x = 127.5 + 0.5*z (planar, convex) so the fan is clean. */
    {
        float inner[12] = { 0,0,127.5f,  2,0,128.5f,  2,2,128.5f,  0,2,127.5f };
        float outer[12] = { -1,-1,127.0f,  3,-1,129.0f,  3,3,129.0f, -1,3,127.0f };
        ComponentMesh cm; build_quad_annulus(ar, &cm, inner, outer);
        size_t nf0 = cm.nf;
        SeamHoleFill_process(ar, &cm, &P, &st);
        int ok = (st.filled == 1 && st.tris_added == 2 && cm.nf == nf0 + 2 &&
                  manifold_ok(cm.faces, cm.nf));
        printf("[G5] quad seam puncture: filled=%zu tris+=%zu nf %zu->%zu manifold=%d -> %s\n",
               st.filled, st.tris_added, nf0, cm.nf, manifold_ok(cm.faces, cm.nf),
               ok ? "ok" : "FAIL");
        if (!ok) fails++;
    }

    /* G6: trivial input -- empty mesh + lone triangle -> no crash, nothing filled. */
    {
        ComponentMesh cm; memset(&cm, 0, sizeof cm); cm.self = &cm;
        cm.verts = NULL; cm.faces = NULL; cm.nv = 0; cm.nf = 0;
        SeamHoleFill_process(ar, &cm, &P, &st);
        int ok1 = (st.filled == 0);
        float V[9] = { 0,0,127.5f, 2,0,128.5f, 0,2,128.0f };
        int32_t Fl[3] = { 0,1,2 };
        ComponentMesh cm2; memset(&cm2, 0, sizeof cm2); cm2.self = &cm2;
        cm2.verts = (float *)ARENA_ALLOC(ar, sizeof V); memcpy(cm2.verts, V, sizeof V);
        cm2.faces = (int32_t *)ARENA_ALLOC(ar, sizeof Fl); memcpy(cm2.faces, Fl, sizeof Fl);
        cm2.nv = 3; cm2.nf = 1;
        size_t lone_nf0 = cm2.nf;
        SeamHoleFill_process(ar, &cm2, &P, &st);
        /* a lone triangle is a single simple 3-loop that is the SOLE boundary of
         * its component; capping it would double the face into a degenerate
         * bubble. The bubble guard must refuse it. */
        int ok2 = (st.filled == 0 && st.skip_bubble >= 1 && cm2.nf == lone_nf0);
        printf("[G6] trivial: empty filled=%zu(exp0) lone filled=%zu skip_bubble=%zu -> %s\n",
               (size_t)(ok1?0:1), st.filled, st.skip_bubble, (ok1 && ok2) ? "ok" : "FAIL");
        if (!(ok1 && ok2)) fails++;
    }

    /* G7: idempotence -- refilling the G1 result closes nothing more. */
    {
        float inner[9] = { 0,0,127.5f,  2,0,128.5f,  0,2,128.0f };
        float outer[9] = { -1,-1,127.0f,  3,-1,129.0f, -1,3,128.0f };
        ComponentMesh cm; build_tri_annulus(ar, &cm, inner, outer);
        SeamHoleFill_process(ar, &cm, &P, &st);
        size_t nf1 = cm.nf;
        SeamHoleFill_process(ar, &cm, &P, &st);
        int ok = (st.filled == 0 && cm.nf == nf1);
        printf("[G7] idempotent: 2nd pass filled=%zu nf=%zu(==%zu) -> %s\n",
               st.filled, cm.nf, nf1, ok ? "ok" : "FAIL");
        if (!ok) fails++;
    }

    Arena_dispose(&ar);
    printf("=== seam_hole_fill_test %s (%d failure%s) ===\n",
           fails ? "FAILED" : "PASSED", fails, fails == 1 ? "" : "s");
    return fails ? 1 : 0;
}
