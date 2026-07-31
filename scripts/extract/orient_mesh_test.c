/* Unit test for OrientMesh_consistent (src/remesh/orient_mesh.c).
 *
 * Builds a flat NxN triangulated sheet with a known-consistent winding,
 * then deliberately mis-winds subsets of faces and asserts the pass:
 *   - drives same-direction interior edges to 0 (orientable),
 *   - with per-vertex normals, flips each component to align the normals,
 *   - with NULL normals, keeps the MAJORITY orientation even when face 0
 *     (the BFS seed) is the mis-wound one.
 *
 * Plain C99. Links arena.c + except.c. Built by orient_mesh_test.vcxproj
 * (Debug x64) and run after every compilation per CLAUDE.md. */
#include "../../src/remesh/orient_mesh.h"
#include "../../src/common/arena.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define N 12                       /* grid side (verts) */
#define NV (N * N)
#define NCELL ((N - 1) * (N - 1))
#define NF (2 * NCELL)

static int failures = 0;
#define CHECK(cond, msg) do { \
    if (!(cond)) { printf("  FAIL: %s\n", (msg)); failures++; } \
    else        { printf("  ok:   %s\n", (msg)); } \
} while (0)

/* Build a sheet in the y-x plane at z=0; verts (z,y,x). */
static void build_sheet(float *verts, int32_t *faces) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            int v = i * N + j;
            verts[v * 3 + 0] = 0.0f;        /* z */
            verts[v * 3 + 1] = (float)i;    /* y */
            verts[v * 3 + 2] = (float)j;    /* x */
        }
    }
    int f = 0;
    for (int i = 0; i < N - 1; i++) {
        for (int j = 0; j < N - 1; j++) {
            int v00 = i * N + j, v01 = i * N + j + 1;
            int v10 = (i + 1) * N + j, v11 = (i + 1) * N + j + 1;
            faces[f * 3 + 0] = v00; faces[f * 3 + 1] = v01; faces[f * 3 + 2] = v11; f++;
            faces[f * 3 + 0] = v00; faces[f * 3 + 1] = v11; faces[f * 3 + 2] = v10; f++;
        }
    }
}

static void flip_face(int32_t *faces, int f) {
    int32_t t = faces[f * 3 + 1];
    faces[f * 3 + 1] = faces[f * 3 + 2];
    faces[f * 3 + 2] = t;
}

/* Count interior edges whose two faces traverse the shared edge the SAME
 * direction (winding inconsistency). O(NF^2) — fine for this small sheet. */
static int count_same_dir(const int32_t *faces, int nf) {
    int same = 0;
    for (int a = 0; a < nf; a++) {
        for (int ea = 0; ea < 3; ea++) {
            int a0 = faces[a * 3 + ea], a1 = faces[a * 3 + (ea + 1) % 3];
            for (int b = a + 1; b < nf; b++) {
                for (int eb = 0; eb < 3; eb++) {
                    int b0 = faces[b * 3 + eb], b1 = faces[b * 3 + (eb + 1) % 3];
                    if (a0 == b0 && a1 == b1) same++;        /* same dir */
                }
            }
        }
    }
    return same;
}

/* Geometric face normal z-component sign (sheet lies in y-x plane). */
static float face_nz(const float *verts, const int32_t *faces, int f) {
    int a = faces[f * 3 + 0], b = faces[f * 3 + 1], c = faces[f * 3 + 2];
    float e1y = verts[b * 3 + 1] - verts[a * 3 + 1];
    float e1x = verts[b * 3 + 2] - verts[a * 3 + 2];
    float e2y = verts[c * 3 + 1] - verts[a * 3 + 1];
    float e2x = verts[c * 3 + 2] - verts[a * 3 + 2];
    return e1y * e2x - e1x * e2y; /* the (z,y,x) n0 slot for a planar sheet */
}

static int all_same_nz_sign(const float *verts, const int32_t *faces, int nf) {
    int pos = 0, neg = 0;
    for (int f = 0; f < nf; f++) {
        float nz = face_nz(verts, faces, f);
        if (nz > 0.0f) pos++; else if (nz < 0.0f) neg++;
    }
    return (pos == 0 || neg == 0);
}

int main(void) {
    Arena_T arena = Arena_new();
    float *verts = (float *)ARENA_ALLOC(arena, (long)(NV * 3 * sizeof(float)));
    int32_t *faces = (int32_t *)ARENA_ALLOC(arena, (long)(NF * 3 * sizeof(int32_t)));
    int32_t *gold = (int32_t *)ARENA_ALLOC(arena, (long)(NF * 3 * sizeof(int32_t)));

    build_sheet(verts, gold);
    float gold_nz = face_nz(verts, gold, 0); /* reference winding sign */
    printf("sheet: %d verts, %d faces; reference face0 nz=%.1f\n",
           NV, NF, (double)gold_nz);

    /* --- Test 0: trivial inputs --- */
    {
        size_t fl = 99, cp = 99, rs = 99;
        int rc = OrientMesh_consistent(arena, verts, NV, NULL, faces, 0,
                                       &fl, &cp, &rs);
        CHECK(rc == 0 && fl == 0 && cp == 0 && rs == 0, "nf=0 returns clean");
    }

    /* --- Test 1: clean sheet stays clean (0 flips, 0 residual) --- */
    {
        memcpy(faces, gold, sizeof(int32_t) * NF * 3);
        CHECK(count_same_dir(faces, NF) == 0, "gold sheet has 0 same-dir");
        size_t fl = 0, cp = 0, rs = 0;
        OrientMesh_consistent(arena, verts, NV, NULL, faces, NF, &fl, &cp, &rs);
        CHECK(cp == 1, "single connected component");
        CHECK(rs == 0, "clean sheet: 0 residual same-dir");
        CHECK(memcmp(faces, gold, sizeof(int32_t) * NF * 3) == 0,
              "clean sheet unchanged (majority anchor)");
    }

    /* --- Test 2: mis-wind a minority subset, NULL normals --- */
    {
        memcpy(faces, gold, sizeof(int32_t) * NF * 3);
        int flipped = 0;
        for (int f = 0; f < NF; f += 3) { flip_face(faces, f); flipped++; }
        CHECK(count_same_dir(faces, NF) > 0, "minority mis-wind makes same-dir");
        CHECK(flipped * 2 < NF, "flipped set is a minority");
        size_t fl = 0, cp = 0, rs = 0;
        OrientMesh_consistent(arena, verts, NV, NULL, faces, NF, &fl, &cp, &rs);
        CHECK(rs == 0, "after orient: 0 residual same-dir");
        CHECK(all_same_nz_sign(verts, faces, NF), "all faces same normal sign");
        CHECK(memcmp(faces, gold, sizeof(int32_t) * NF * 3) == 0,
              "majority anchor restores gold winding");
    }

    /* --- Test 3: mis-wind face 0 only (the BFS seed), NULL normals --- */
    {
        memcpy(faces, gold, sizeof(int32_t) * NF * 3);
        flip_face(faces, 0);
        size_t fl = 0, cp = 0, rs = 0;
        OrientMesh_consistent(arena, verts, NV, NULL, faces, NF, &fl, &cp, &rs);
        CHECK(rs == 0, "seed-flipped: 0 residual after orient");
        CHECK(memcmp(faces, gold, sizeof(int32_t) * NF * 3) == 0,
              "majority anchor ignores mis-wound seed, restores gold");
    }

    /* --- Test 4: flip ALL faces, anchor with normals to original sign --- */
    {
        memcpy(faces, gold, sizeof(int32_t) * NF * 3);
        for (int f = 0; f < NF; f++) flip_face(faces, f);
        /* per-vertex normals all point in the gold-winding normal direction */
        float *nrm = (float *)ARENA_ALLOC(arena, (long)(NV * 3 * sizeof(float)));
        for (int v = 0; v < NV; v++) {
            nrm[v * 3 + 0] = (gold_nz > 0.0f) ? 1.0f : -1.0f;
            nrm[v * 3 + 1] = 0.0f;
            nrm[v * 3 + 2] = 0.0f;
        }
        size_t fl = 0, cp = 0, rs = 0;
        OrientMesh_consistent(arena, verts, NV, nrm, faces, NF, &fl, &cp, &rs);
        CHECK(rs == 0, "all-flipped: 0 residual after orient");
        CHECK(memcmp(faces, gold, sizeof(int32_t) * NF * 3) == 0,
              "normal anchor flips component back to gold winding");
    }

    Arena_dispose(&arena);
    if (failures == 0) { printf("ALL ORIENT_MESH TESTS PASSED\n"); return 0; }
    printf("%d ORIENT_MESH TEST(S) FAILED\n", failures);
    return 1;
}
