/* Unit test for FoldCleanup_process (src/remesh/fold_cleanup.c).
 *
 * A "fold" is an interior edge {u,v} whose two faces both traverse it u->v
 * (a winding inconsistency a BFS cannot flip away). FoldCleanup deletes the
 * low-coherence outlier face of each such edge. Four checks:
 *
 *   A. Clean consistently-wound fan -> NO face removed (no false positive).
 *   B. dot<0 geometric fold flap (apex folded to the far side, opposing
 *      normal) -> the flap is removed; result is winding-consistent + manifold.
 *   C. dot>=0 winding defect whose worse face is a clear low-coherence outlier
 *      (coher 0 < FOLD_OUTLIER_COHER) -> outlier removed; result consistent.
 *   D. dot>=0 BALANCED knot (both faces high-coherence) -> LEFT untouched;
 *      the same_dir edge survives (we never cut a balanced non-orientable
 *      patch). This mirrors the policy gate in fold_cleanup.c.
 *
 * Plain C99. Links fold_cleanup.c + arena.c + except.c. Built + run after
 * every compilation per CLAUDE.md.
 *
 * cl /nologo /W4 /Zi /Fe:output\tools\fold_cleanup_test.exe /Fo:output\tools\ \
 *    scripts\step0-mesh-extract\fold_cleanup_test.c src\remesh\fold_cleanup.c \
 *    src\common\arena.c src\common\except.c
 */
#include "../../src/remesh/fold_cleanup.h"
#include "../../src/common/mesh_types.h"
#include "../../src/common/arena.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int failures = 0;
#define CHECK(cond, msg) do { \
    if (!(cond)) { printf("  FAIL: %s\n", (msg)); failures++; } \
    else        { printf("  ok:   %s\n", (msg)); } \
} while (0)

/* ---- tiny O(nf^2) diagnostics (meshes here are a handful of faces) ---- */

/* interior edges whose two faces traverse the shared edge the SAME direction */
static int count_same_dir(const int32_t *faces, int nf) {
    int same = 0;
    for (int a = 0; a < nf; a++) {
        for (int ea = 0; ea < 3; ea++) {
            int a0 = faces[a * 3 + ea], a1 = faces[a * 3 + (ea + 1) % 3];
            for (int b = a + 1; b < nf; b++) {
                for (int eb = 0; eb < 3; eb++) {
                    int b0 = faces[b * 3 + eb], b1 = faces[b * 3 + (eb + 1) % 3];
                    if (a0 == b0 && a1 == b1) { same++; }
                }
            }
        }
    }
    return same;
}

/* any undirected edge shared by >2 faces => non-manifold */
static int has_nonmanifold_edge(const int32_t *faces, int nf) {
    for (int f = 0; f < nf; f++) {
        for (int e = 0; e < 3; e++) {
            int a = faces[f * 3 + e], b = faces[f * 3 + (e + 1) % 3];
            int lo = a < b ? a : b, hi = a < b ? b : a;
            int occ = 0;
            for (int g = 0; g < nf; g++) {
                for (int eg = 0; eg < 3; eg++) {
                    int c = faces[g * 3 + eg], d = faces[g * 3 + (eg + 1) % 3];
                    int clo = c < d ? c : d, chi = c < d ? d : c;
                    if (clo == lo && chi == hi) { occ++; }
                }
            }
            if (occ > 2) { return 1; }
        }
    }
    return 0;
}

static void set_mesh(ComponentMesh *cm, float *verts, int32_t *faces,
                     size_t nv, size_t nf) {
    memset(cm, 0, sizeof(*cm));
    cm->verts = verts; cm->faces = faces; cm->nv = nv; cm->nf = nf;
    cm->comp_id = 1; cm->self = cm;
}

/* common 5-vertex / 3-triangle consistently-wound planar fan around hub 0.
 * verts are (z,y,x); all z=0 so face normals lie on the z axis. */
static void build_fan(float *verts, int32_t *faces) {
    float pos[5][3] = {
        {0,0,0},   /* 0 hub  c  */
        {0,0,1},   /* 1      p1 */
        {0,1,1},   /* 2      p2 */
        {0,2,1},   /* 3      p3 */
        {0,3,1}    /* 4      p4 */
    };
    for (int i = 0; i < 5; i++) {
        verts[i*3+0] = pos[i][0]; verts[i*3+1] = pos[i][1]; verts[i*3+2] = pos[i][2];
    }
    int32_t f[3 * 3] = { 0,1,2,  0,2,3,  0,3,4 };
    memcpy(faces, f, sizeof(f));
}

int main(void) {
    Arena_T arena = Arena_new();

    /* =========================================================
     * A — clean fan: nothing is a fold, nothing removed
     * ========================================================= */
    {
        float *verts = (float *)ARENA_ALLOC(arena, (long)(5 * 3 * sizeof(float)));
        int32_t *faces = (int32_t *)ARENA_ALLOC(arena, (long)(3 * 3 * sizeof(int32_t)));
        build_fan(verts, faces);

        ComponentMesh cm; set_mesh(&cm, verts, faces, 5, 3);
        CHECK(count_same_dir(cm.faces, (int)cm.nf) == 0, "A: clean fan has 0 same_dir");

        size_t removed = 9;
        FoldCleanup_process(arena, &cm, 1, 8, &removed);
        CHECK(removed == 0, "A: clean fan -> 0 faces removed");
        CHECK(cm.nf == 3, "A: face count unchanged");
        CHECK(count_same_dir(cm.faces, (int)cm.nf) == 0, "A: still 0 same_dir");
        CHECK(!has_nonmanifold_edge(cm.faces, (int)cm.nf), "A: still manifold");
    }

    /* =========================================================
     * B — dot<0 geometric fold flap is deleted
     * Flap (1,2,5) reuses directed edge 1->2 (same as face 0) with apex 5 on
     * the far side & below the plane: opposing normal (dot<0). Its other two
     * edges are boundary (coher 0); face 0 has a neighbour (coher 1) -> the
     * flap is the outlier and is cut.
     * ========================================================= */
    {
        float *verts = (float *)ARENA_ALLOC(arena, (long)(6 * 3 * sizeof(float)));
        int32_t *faces = (int32_t *)ARENA_ALLOC(arena, (long)(4 * 3 * sizeof(int32_t)));
        build_fan(verts, faces);                 /* verts 0..4, faces 0..2 */
        verts[5*3+0] = -1.0f; verts[5*3+1] = 0.5f; verts[5*3+2] = 2.0f; /* apex q */
        faces[3*3+0] = 1; faces[3*3+1] = 2; faces[3*3+2] = 5;           /* flap */

        ComponentMesh cm; set_mesh(&cm, verts, faces, 6, 4);
        CHECK(count_same_dir(cm.faces, (int)cm.nf) == 1, "B: fold present (1 same_dir)");

        size_t removed = 9;
        FoldCleanup_process(arena, &cm, 1, 8, &removed);
        CHECK(removed == 1, "B: exactly the flap removed");
        CHECK(cm.nf == 3, "B: face count back to fan");
        CHECK(count_same_dir(cm.faces, (int)cm.nf) == 0, "B: fold gone (0 same_dir)");
        CHECK(!has_nonmanifold_edge(cm.faces, (int)cm.nf), "B: manifold after");
    }

    /* =========================================================
     * C — dot>=0 winding defect, worse face is low-coherence outlier
     * Flap (1,2,5) reuses 1->2, apex 5 coplanar and on the SAME side as the
     * hub: aligned normal (dot>=0) but it is an isolated flap (coher 0 <
     * FOLD_OUTLIER_COHER) -> removed.
     * ========================================================= */
    {
        float *verts = (float *)ARENA_ALLOC(arena, (long)(6 * 3 * sizeof(float)));
        int32_t *faces = (int32_t *)ARENA_ALLOC(arena, (long)(4 * 3 * sizeof(int32_t)));
        build_fan(verts, faces);
        verts[5*3+0] = 0.0f; verts[5*3+1] = 0.5f; verts[5*3+2] = 0.5f; /* coplanar, hub side */
        faces[3*3+0] = 1; faces[3*3+1] = 2; faces[3*3+2] = 5;

        ComponentMesh cm; set_mesh(&cm, verts, faces, 6, 4);
        CHECK(count_same_dir(cm.faces, (int)cm.nf) == 1, "C: winding defect present");

        size_t removed = 9;
        FoldCleanup_process(arena, &cm, 1, 8, &removed);
        CHECK(removed == 1, "C: low-coherence outlier removed");
        CHECK(cm.nf == 3, "C: face count back to fan");
        CHECK(count_same_dir(cm.faces, (int)cm.nf) == 0, "C: defect gone (0 same_dir)");
        CHECK(!has_nonmanifold_edge(cm.faces, (int)cm.nf), "C: manifold after");
    }

    /* =========================================================
     * D — dot>=0 BALANCED knot is left alone
     * Two faces A,B share edge {0,1} both 0->1 (same_dir), aligned normals,
     * and EACH has an aligned neighbour across its other edge so both have
     * coher 1.0 >= FOLD_OUTLIER_COHER. Policy: leave balanced patches.
     * ========================================================= */
    {
        float pos[6][3] = {
            {0,1,0},   /* 0 U */
            {0,1,3},   /* 1 V */
            {0,0,1},   /* 2 A_apex      (hub side, y=0) */
            {0,-1,1},  /* 3 A_nb_apex   */
            {0,0,2},   /* 4 B_apex      (hub side, y=0) */
            {0,-1,2}   /* 5 B_nb_apex   */
        };
        float *verts = (float *)ARENA_ALLOC(arena, (long)(6 * 3 * sizeof(float)));
        for (int i = 0; i < 6; i++) {
            verts[i*3+0]=pos[i][0]; verts[i*3+1]=pos[i][1]; verts[i*3+2]=pos[i][2];
        }
        /* A=(0,1,2) B=(0,1,4) share {0,1} both 0->1; A_nb=(0,2,3) shares {0,2}
         * with A, B_nb=(0,4,5) shares {0,4} with B (both opposite-dir => not
         * folds, raise coherence). */
        int32_t f[4 * 3] = { 0,1,2,  0,1,4,  0,2,3,  0,4,5 };
        int32_t *faces = (int32_t *)ARENA_ALLOC(arena, (long)(4 * 3 * sizeof(int32_t)));
        memcpy(faces, f, sizeof(f));

        ComponentMesh cm; set_mesh(&cm, verts, faces, 6, 4);
        CHECK(count_same_dir(cm.faces, (int)cm.nf) == 1, "D: one balanced same_dir edge");

        size_t removed = 9;
        FoldCleanup_process(arena, &cm, 1, 8, &removed);
        CHECK(removed == 0, "D: balanced knot left untouched (0 removed)");
        CHECK(cm.nf == 4, "D: face count unchanged");
        CHECK(count_same_dir(cm.faces, (int)cm.nf) == 1, "D: same_dir edge survives");
    }

    Arena_dispose(&arena);
    if (failures == 0) { printf("ALL FOLD_CLEANUP TESTS PASSED\n"); return 0; }
    printf("%d FOLD_CLEANUP TEST(S) FAILED\n", failures);
    return 1;
}
