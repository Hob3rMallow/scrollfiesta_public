/* Unit test for PinholeFill_process (src/remesh/pinhole_fill.c).
 *
 * Three independent checks:
 *   A. Single missing triangle in a flat sheet is closed by exactly one fill
 *      triangle, the mesh stays winding-consistent (0 same-dir edges), and the
 *      filled hole uses precisely the three hole vertices.
 *   B. A bowtie / pinch vertex (two triangle fans meeting at one vertex with no
 *      shared edge) is split into two vertices, leaving zero pinch vertices.
 *      Run with everything pinned so phase 2 cannot also fill the fan rims —
 *      this isolates the phase-1 split.
 *   C. A 3-edge hole whose vertices are pinned is SKIPPED when respect_pins=1,
 *      and the same hole fills when respect_pins=0.
 *
 * Plain C99. Links pinhole_fill.c + arena.c + except.c. Built by
 * pinhole_fill_test.vcxproj (Debug x64) and run after every compilation
 * per CLAUDE.md. */
#include "../../src/remesh/pinhole_fill.h"
#include "../../src/common/mesh_types.h"
#include "../../src/common/arena.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#define N 5                        /* sheet side (verts) */
#define NV (N * N)
#define NCELL ((N - 1) * (N - 1))
#define NF (2 * NCELL)

static int failures = 0;
#define CHECK(cond, msg) do { \
    if (!(cond)) { printf("  FAIL: %s\n", (msg)); failures++; } \
    else        { printf("  ok:   %s\n", (msg)); } \
} while (0)

/* ---- small mesh diagnostics (O(nf^2): meshes here are tiny) ---- */

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

/* unique undirected edges incident to exactly one face */
static int count_boundary_edges(const int32_t *faces, int nf) {
    int nb = 0;
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
            if (occ == 1) { nb++; } /* boundary edge occurs once total */
        }
    }
    return nb;
}

/* any interior edge shared by >2 faces => non-manifold */
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

static int uf_find_local(int *uf, int a) {
    while (uf[a] != a) { uf[a] = uf[uf[a]]; a = uf[a]; }
    return a;
}

/* vertices whose incident faces fall into >=2 edge-connected fans (pinch) */
static int count_pinch_verts(const int32_t *faces, int nf, int nv) {
    int pinch = 0;
    for (int v = 0; v < nv; v++) {
        int inc[256], ni = 0;
        for (int f = 0; f < nf && ni < 256; f++) {
            if (faces[f * 3] == v || faces[f * 3 + 1] == v ||
                faces[f * 3 + 2] == v) { inc[ni++] = f; }
        }
        if (ni < 2) { continue; }
        int uf[256];
        for (int i = 0; i < ni; i++) { uf[i] = i; }
        for (int i = 0; i < ni; i++) {
            for (int j = i + 1; j < ni; j++) {
                int shared = 0;
                for (int ki = 0; ki < 3; ki++) {
                    int wi = faces[inc[i] * 3 + ki];
                    if (wi == v) { continue; }
                    for (int kj = 0; kj < 3; kj++) {
                        int wj = faces[inc[j] * 3 + kj];
                        if (wj != v && wi == wj) { shared = 1; }
                    }
                }
                if (shared) {
                    int ra = uf_find_local(uf, i), rb = uf_find_local(uf, j);
                    if (ra != rb) { uf[ra] = rb; }
                }
            }
        }
        int fans = 0;
        for (int i = 0; i < ni; i++) { if (uf_find_local(uf, i) == i) { fans++; } }
        if (fans >= 2) { pinch++; }
    }
    return pinch;
}

/* ---- builders ---- */

static void build_sheet(float *verts, int32_t *faces) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            int v = i * N + j;
            verts[v * 3 + 0] = 0.0f;
            verts[v * 3 + 1] = (float)i;
            verts[v * 3 + 2] = (float)j;
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

static void set_mesh(ComponentMesh *cm, float *verts, int32_t *faces,
                     size_t nv, size_t nf, uint8_t *pin) {
    memset(cm, 0, sizeof(*cm));
    cm->verts = verts; cm->faces = faces; cm->nv = nv; cm->nf = nf;
    cm->pin_mask = pin; cm->vert_normals = NULL;
    cm->comp_id = 1; cm->self = cm;
}

int main(void) {
    Arena_T arena = Arena_new();

    /* =========================================================
     * Test A — single missing triangle fills with right winding
     * ========================================================= */
    {
        float *verts = (float *)ARENA_ALLOC(arena, (long)(NV * 3 * sizeof(float)));
        int32_t *gold = (int32_t *)ARENA_ALLOC(arena, (long)(NF * 3 * sizeof(int32_t)));
        build_sheet(verts, gold);

        /* remove first triangle of interior cell (i=2,j=2) -> a 3-hole */
        int cell = 2 * (N - 1) + 2;
        int tri = cell * 2;
        int h0 = gold[tri * 3 + 0], h1 = gold[tri * 3 + 1], h2 = gold[tri * 3 + 2];

        int32_t *faces = (int32_t *)ARENA_ALLOC(arena,
                              (long)((NF - 1) * 3 * sizeof(int32_t)));
        int w = 0;
        for (int f = 0; f < NF; f++) {
            if (f == tri) { continue; }
            faces[w * 3 + 0] = gold[f * 3 + 0];
            faces[w * 3 + 1] = gold[f * 3 + 1];
            faces[w * 3 + 2] = gold[f * 3 + 2];
            w++;
        }

        ComponentMesh cm;
        set_mesh(&cm, verts, faces, NV, NF - 1, NULL);

        int perim = 4 * (N - 1);
        CHECK(count_same_dir(cm.faces, (int)cm.nf) == 0, "A: holed sheet winding consistent");
        CHECK(count_boundary_edges(cm.faces, (int)cm.nf) == perim + 3,
              "A: hole adds 3 boundary edges");
        CHECK(count_pinch_verts(cm.faces, (int)cm.nf, (int)cm.nv) == 0,
              "A: no pinch verts before");

        size_t sp = 9, fl = 9, ad = 9, sk = 9;
        PinholeFill_process(arena, &cm, 1, 0, &sp, &fl, &ad, &sk);

        CHECK(sp == 0, "A: no pinch splits");
        CHECK(fl == 1 && ad == 1 && sk == 0, "A: one loop filled, one tri added");
        CHECK(cm.nf == (size_t)NF, "A: face count back to original");
        CHECK(count_boundary_edges(cm.faces, (int)cm.nf) == perim,
              "A: hole closed (only outer boundary remains)");
        CHECK(count_same_dir(cm.faces, (int)cm.nf) == 0,
              "A: fill triangle wound consistently (0 same-dir)");
        CHECK(!has_nonmanifold_edge(cm.faces, (int)cm.nf),
              "A: no non-manifold edge introduced");

        /* the added triangle must use exactly the three hole verts */
        int32_t *t = &cm.faces[(NF - 1) * 3];
        int uses_all = ((t[0] == h0 || t[0] == h1 || t[0] == h2) &&
                        (t[1] == h0 || t[1] == h1 || t[1] == h2) &&
                        (t[2] == h0 || t[2] == h1 || t[2] == h2) &&
                        t[0] != t[1] && t[1] != t[2] && t[0] != t[2]);
        CHECK(uses_all, "A: fill triangle uses the 3 hole vertices");
    }

    /* =========================================================
     * Test B — bowtie pinch vertex is split (phase 1 isolated)
     * ========================================================= */
    {
        /* C=0 shared by two disjoint 2-triangle fans (A: 1,2,3  B: 4,5,6) */
        float *verts = (float *)ARENA_ALLOC(arena, (long)(7 * 3 * sizeof(float)));
        float pos[7][3] = {
            {0,0,0}, {0,1,0}, {0,1,1}, {0,1,2},   /* C, A1, A2, A3 */
            {0,-1,0}, {0,-1,1}, {0,-1,2}          /* B1, B2, B3 */
        };
        for (int i = 0; i < 7; i++) {
            verts[i * 3 + 0] = pos[i][0];
            verts[i * 3 + 1] = pos[i][1];
            verts[i * 3 + 2] = pos[i][2];
        }
        int32_t f0[4 * 3] = {
            0, 1, 2,   0, 2, 3,    /* fan A, shares edge C-A2 */
            0, 4, 5,   0, 5, 6     /* fan B, shares edge C-B2 */
        };
        int32_t *faces = (int32_t *)ARENA_ALLOC(arena, (long)(4 * 3 * sizeof(int32_t)));
        memcpy(faces, f0, sizeof(f0));
        uint8_t *pin = (uint8_t *)ARENA_ALLOC(arena, (long)(7 * sizeof(uint8_t)));
        memset(pin, 1, 7);                         /* pin all => phase 2 no-op */

        ComponentMesh cm;
        set_mesh(&cm, verts, faces, 7, 4, pin);

        CHECK(count_pinch_verts(cm.faces, (int)cm.nf, (int)cm.nv) == 1,
              "B: one pinch vertex before");

        size_t sp = 9, fl = 9, ad = 9, sk = 9;
        PinholeFill_process(arena, &cm, 1, 1, &sp, &fl, &ad, &sk); /* respect_pins=1 */

        CHECK(sp == 1, "B: one pinch split");
        CHECK(cm.nv == 8, "B: vertex count grew by one");
        CHECK(fl == 0 && ad == 0, "B: phase 2 skipped (all pinned)");
        CHECK(count_pinch_verts(cm.faces, (int)cm.nf, (int)cm.nv) == 0,
              "B: zero pinch verts after split");
        CHECK(cm.nf == 4, "B: face count unchanged");
    }

    /* =========================================================
     * Test C — pinned 3-hole skipped, then fills when unpinned
     * ========================================================= */
    {
        float *verts = (float *)ARENA_ALLOC(arena, (long)(NV * 3 * sizeof(float)));
        int32_t *gold = (int32_t *)ARENA_ALLOC(arena, (long)(NF * 3 * sizeof(int32_t)));
        build_sheet(verts, gold);
        int cell = 2 * (N - 1) + 2;
        int tri = cell * 2;
        int h0 = gold[tri * 3 + 0], h1 = gold[tri * 3 + 1], h2 = gold[tri * 3 + 2];

        int32_t *faces = (int32_t *)ARENA_ALLOC(arena,
                              (long)((NF - 1) * 3 * sizeof(int32_t)));
        int w = 0;
        for (int f = 0; f < NF; f++) {
            if (f == tri) { continue; }
            faces[w * 3 + 0] = gold[f * 3 + 0];
            faces[w * 3 + 1] = gold[f * 3 + 1];
            faces[w * 3 + 2] = gold[f * 3 + 2];
            w++;
        }
        uint8_t *pin = (uint8_t *)ARENA_CALLOC(arena, (long)NV, (long)sizeof(uint8_t));
        pin[h0] = 1; pin[h1] = 1; pin[h2] = 1;

        ComponentMesh cm;
        set_mesh(&cm, verts, faces, NV, NF - 1, pin);

        int perim = 4 * (N - 1);

        /* respect pins -> hole left open */
        size_t sp = 9, fl = 9, ad = 9, sk = 9;
        PinholeFill_process(arena, &cm, 1, 1, &sp, &fl, &ad, &sk);
        CHECK(fl == 0 && ad == 0 && sk == 1, "C: pinned 3-hole skipped");
        CHECK(cm.nf == (size_t)(NF - 1), "C: face count unchanged when skipped");
        CHECK(count_boundary_edges(cm.faces, (int)cm.nf) == perim + 3,
              "C: hole still open when pinned");

        /* ignore pins -> hole fills */
        sp = 9; fl = 9; ad = 9; sk = 9;
        PinholeFill_process(arena, &cm, 1, 0, &sp, &fl, &ad, &sk);
        CHECK(fl == 1 && ad == 1 && sk == 0, "C: same hole fills when unpinned");
        CHECK(cm.nf == (size_t)NF, "C: face count restored after fill");
        CHECK(count_boundary_edges(cm.faces, (int)cm.nf) == perim,
              "C: hole closed when unpinned");
    }

    /* =========================================================
     * Test D — a 3-loop with a coincident vertex pair is NOT
     * filled (would be a zero-area sliver), it is skipped.
     * ========================================================= */
    {
        float *verts = (float *)ARENA_ALLOC(arena, (long)(NV * 3 * sizeof(float)));
        int32_t *gold = (int32_t *)ARENA_ALLOC(arena, (long)(NF * 3 * sizeof(int32_t)));
        build_sheet(verts, gold);
        int cell = 2 * (N - 1) + 2;
        int tri = cell * 2;
        int h0 = gold[tri * 3 + 0], h1 = gold[tri * 3 + 1];

        int32_t *faces = (int32_t *)ARENA_ALLOC(arena,
                              (long)((NF - 1) * 3 * sizeof(int32_t)));
        int w = 0;
        for (int f = 0; f < NF; f++) {
            if (f == tri) { continue; }
            faces[w * 3 + 0] = gold[f * 3 + 0];
            faces[w * 3 + 1] = gold[f * 3 + 1];
            faces[w * 3 + 2] = gold[f * 3 + 2];
            w++;
        }
        /* collapse hole vert h1 onto h0 -> coincident pair on the 3-loop */
        verts[h1 * 3 + 0] = verts[h0 * 3 + 0];
        verts[h1 * 3 + 1] = verts[h0 * 3 + 1];
        verts[h1 * 3 + 2] = verts[h0 * 3 + 2];

        ComponentMesh cm;
        set_mesh(&cm, verts, faces, NV, NF - 1, NULL);

        size_t sp = 9, fl = 9, ad = 9, sk = 9;
        PinholeFill_process(arena, &cm, 1, 0, &sp, &fl, &ad, &sk);
        CHECK(fl == 0 && ad == 0, "D: coincident-vertex loop not filled");
        CHECK(sk == 1, "D: coincident-vertex loop counted as skipped");
        CHECK(cm.nf == (size_t)(NF - 1), "D: no sliver triangle added");
    }

    /* =========================================================
     * Test E — the boundary loop of an existing 2-triangle island
     * must NOT be fanned: doing so re-covers both triangles (a
     * doubled-triangle pocket) and pushes their shared edge to a
     * 4-fan. This is the real grid_weld seam failure; the manifold
     * guard must SKIP the loop, leaving the island untouched.
     * ========================================================= */
    {
        /* flat quad split into 2 tris sharing edge {0,2}; 4 boundary edges */
        float pos[4][3] = { {0,0,0}, {0,0,1}, {0,1,1}, {0,1,0} };
        float *verts = (float *)ARENA_ALLOC(arena, (long)(4 * 3 * sizeof(float)));
        for (int i = 0; i < 4; i++) {
            verts[i * 3 + 0] = pos[i][0];
            verts[i * 3 + 1] = pos[i][1];
            verts[i * 3 + 2] = pos[i][2];
        }
        int32_t f0[2 * 3] = { 0, 1, 2,   0, 2, 3 };
        int32_t *faces = (int32_t *)ARENA_ALLOC(arena, (long)(2 * 3 * sizeof(int32_t)));
        memcpy(faces, f0, sizeof(f0));

        ComponentMesh cm;
        set_mesh(&cm, verts, faces, 4, 2, NULL);

        CHECK(count_same_dir(cm.faces, (int)cm.nf) == 0, "E: island winding consistent before");
        CHECK(!has_nonmanifold_edge(cm.faces, (int)cm.nf), "E: island manifold before");

        size_t sp = 9, fl = 9, ad = 9, sk = 9;
        PinholeFill_process(arena, &cm, 1, 0, &sp, &fl, &ad, &sk);

        CHECK(fl == 0 && ad == 0, "E: island back-loop not filled");
        CHECK(sk >= 1, "E: island back-loop counted as skipped");
        CHECK(cm.nf == 2, "E: face count unchanged (no doubled triangles)");
        CHECK(!has_nonmanifold_edge(cm.faces, (int)cm.nf),
              "E: no 4-fan non-manifold edge introduced");
    }

    /* =========================================================
     * Test F — no-merger diameter gate: the SAME single missing
     * triangle as test A, but scaled so the hole is wider than
     * PINHOLE_MAX_DIAM_VOX, is REFUSED (a wide opening might bridge
     * two wraps). It stays open; nothing is filled.
     * ========================================================= */
    {
        float *verts = (float *)ARENA_ALLOC(arena, (long)(NV * 3 * sizeof(float)));
        int32_t *gold = (int32_t *)ARENA_ALLOC(arena, (long)(NF * 3 * sizeof(int32_t)));
        build_sheet(verts, gold);
        /* scale the grid x5: unit cell -> 5 vox, hole diameter 5*sqrt(2)~7.07 */
        for (int v = 0; v < NV; v++) {
            verts[v * 3 + 1] *= 5.0f; verts[v * 3 + 2] *= 5.0f;
        }
        int cell = 2 * (N - 1) + 2;
        int tri = cell * 2;

        int32_t *faces = (int32_t *)ARENA_ALLOC(arena,
                              (long)((NF - 1) * 3 * sizeof(int32_t)));
        int w = 0;
        for (int f = 0; f < NF; f++) {
            if (f == tri) { continue; }
            faces[w * 3 + 0] = gold[f * 3 + 0];
            faces[w * 3 + 1] = gold[f * 3 + 1];
            faces[w * 3 + 2] = gold[f * 3 + 2];
            w++;
        }
        ComponentMesh cm;
        set_mesh(&cm, verts, faces, NV, NF - 1, NULL);
        int perim = 4 * (N - 1);

        size_t sp = 9, fl = 9, ad = 9, sk = 9;
        PinholeFill_process(arena, &cm, 1, 0, &sp, &fl, &ad, &sk);
        CHECK(fl == 0 && ad == 0, "F: wide hole refused (diameter gate)");
        CHECK(sk == 1, "F: wide hole counted as skipped");
        CHECK(cm.nf == (size_t)(NF - 1), "F: no fill triangle added");
        CHECK(count_boundary_edges(cm.faces, (int)cm.nf) == perim + 3,
              "F: wide hole still open");
    }

    /* =========================================================
     * Test G — PinholeFill closes ONLY exact single-triangle (3-loop)
     * holes; a larger loop is left to HoleFill's CDT path. A 10-vertex
     * annulus (inner rim = a decagon hole of diameter 4.0; outer rim =
     * a wide decagon of diameter 7.0). NEITHER rim is a 3-loop, so the
     * inner 10-loop is left OPEN (not fan-filled) and the wide outer
     * rim is refused by the diameter gate. Face count unchanged.
     * (Historical note: an older pinhole fan-filled loops up to
     * PINHOLE_MAX_LOOP; that was removed — 4+ loops now go to HoleFill.)
     * ========================================================= */
    {
        enum { M = 10 };
        const float inR = 2.0f, outR = 3.5f;   /* diams 4.0 and 7.0 */
        float *verts = (float *)ARENA_ALLOC(arena, (long)(2 * M * 3 * sizeof(float)));
        for (int i = 0; i < M; i++) {
            float a = (float)(2.0 * 3.14159265358979323846 * i / M);
            verts[i * 3 + 0] = 0.0f;
            verts[i * 3 + 1] = inR * sinf(a);
            verts[i * 3 + 2] = inR * cosf(a);
            verts[(M + i) * 3 + 0] = 0.0f;
            verts[(M + i) * 3 + 1] = outR * sinf(a);
            verts[(M + i) * 3 + 2] = outR * cosf(a);
        }
        int32_t *faces = (int32_t *)ARENA_ALLOC(arena,
                              (long)(2 * M * 3 * sizeof(int32_t)));
        int nf = 0;
        for (int i = 0; i < M; i++) {
            int j = (i + 1) % M, Mi = M + i, Mj = M + j;
            faces[nf * 3 + 0] = i;  faces[nf * 3 + 1] = j;  faces[nf * 3 + 2] = Mj; nf++;
            faces[nf * 3 + 0] = i;  faces[nf * 3 + 1] = Mj; faces[nf * 3 + 2] = Mi; nf++;
        }
        ComponentMesh cm;
        set_mesh(&cm, verts, faces, 2 * M, (size_t)nf, NULL);
        CHECK(count_same_dir(cm.faces, (int)cm.nf) == 0, "G: annulus winding consistent before");
        CHECK(count_boundary_edges(cm.faces, (int)cm.nf) == 2 * M,
              "G: annulus has inner+outer rim (2M boundary edges)");

        size_t sp = 9, fl = 9, ad = 9, sk = 9;
        PinholeFill_process(arena, &cm, 1, 0, &sp, &fl, &ad, &sk);
        CHECK(fl == 0 && ad == 0, "G: neither rim is a 3-loop -> nothing filled");
        CHECK(sk == 1, "G: wide outer rim refused (diameter gate)");
        CHECK(cm.nf == (size_t)(2 * M), "G: face count unchanged");
        CHECK(count_boundary_edges(cm.faces, (int)cm.nf) == 2 * M,
              "G: both rims remain open (4+ loops left for HoleFill)");
    }

    /* =========================================================
     * Test H — dead-zone single triangle: the same lone missing
     * triangle as A, scaled x4 so its diameter ~5.66 vox falls in the
     * OLD dead zone (4.5 < d < 7-vox inter-wrap clearance). A 3-loop's
     * three edges already exist, so capping it adds no new edge and
     * cannot merge two wraps at any diameter -> it must now FILL
     * (PINHOLE_TRI_DIAM_VOX). Regression guard for the seam-hole fix.
     * ========================================================= */
    {
        float *verts = (float *)ARENA_ALLOC(arena, (long)(NV * 3 * sizeof(float)));
        int32_t *gold = (int32_t *)ARENA_ALLOC(arena, (long)(NF * 3 * sizeof(int32_t)));
        build_sheet(verts, gold);
        for (int v = 0; v < NV; v++) {           /* diameter = 4*sqrt(2) ~ 5.66 */
            verts[v * 3 + 1] *= 4.0f; verts[v * 3 + 2] *= 4.0f;
        }
        int cell = 2 * (N - 1) + 2;
        int tri = cell * 2;
        int32_t *faces = (int32_t *)ARENA_ALLOC(arena,
                              (long)((NF - 1) * 3 * sizeof(int32_t)));
        int w = 0;
        for (int f = 0; f < NF; f++) {
            if (f == tri) { continue; }
            faces[w * 3 + 0] = gold[f * 3 + 0];
            faces[w * 3 + 1] = gold[f * 3 + 1];
            faces[w * 3 + 2] = gold[f * 3 + 2];
            w++;
        }
        ComponentMesh cm;
        set_mesh(&cm, verts, faces, NV, NF - 1, NULL);
        int perim = 4 * (N - 1);

        size_t sp = 9, fl = 9, ad = 9, sk = 9;
        PinholeFill_process(arena, &cm, 1, 0, &sp, &fl, &ad, &sk);
        CHECK(fl == 1 && ad == 1 && sk == 0,
              "H: dead-zone (~5.66 vox) single triangle now fills");
        CHECK(cm.nf == (size_t)NF, "H: face count back to original");
        CHECK(count_boundary_edges(cm.faces, (int)cm.nf) == perim,
              "H: hole closed (only outer boundary remains)");
        CHECK(!has_nonmanifold_edge(cm.faces, (int)cm.nf),
              "H: no non-manifold edge introduced");
    }

    Arena_dispose(&arena);
    if (failures == 0) { printf("ALL PINHOLE_FILL TESTS PASSED\n"); return 0; }
    printf("%d PINHOLE_FILL TEST(S) FAILED\n", failures);
    return 1;
}
