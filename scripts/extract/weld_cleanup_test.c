/*
 * weld_cleanup_test.c -- unit tests for src/remesh/weld_cleanup.c.
 *
 * Gates (all on synthetic meshes with known verdicts):
 *   G1  end-to-end (flip + collapse): an interior needle is removed, mesh stays
 *       manifold + degenerate-free, boundary vertices are not moved.
 *   G2  collapse isolated (flip_max_rounds=0): the needle is removed by collapse
 *       alone -> faces 8 -> 6, manifold.
 *   G3  length guard: with max_collapse_len below the needle's short edge AND
 *       flips off, the collapse is refused -> nothing removed (this is the
 *       wrap-safety guard: an over-long edge is never collapsed).
 *   G4  clean mesh: a regular fan has no targets -> zero flips, zero collapses,
 *       faces unchanged.
 *   G5  trivial input: empty mesh and a lone triangle -> no crash, unchanged.
 *   G6  closed mesh interior collapse: an octahedron with one short equator edge
 *       collapses that edge and stays a closed manifold (link condition ACCEPT).
 *   G7  link-condition REJECT: a triangular bipyramid's shared edge has a
 *       non-apex common neighbour, so collapsing it would be non-manifold; with
 *       flips off the collapse must be refused and the mesh left manifold.
 *
 * Build: weld_cleanup_test.vcxproj (links weld_cleanup.c + arena.c + except.c).
 * Run:   weld_cleanup_test.exe   (exit 0 = all pass, 1 = a gate failed)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "common/arena.h"
#include "common/mesh_types.h"
#include "remesh/weld_cleanup.h"

/* ---- manifold / degeneracy checks on a flat face list ---- */
typedef struct { int32_t v0, v1; } UE;
static int ue_cmp(const void *a, const void *b)
{
    const UE *x=(const UE*)a, *y=(const UE*)b;
    if (x->v0!=y->v0) return x->v0<y->v0?-1:1;
    if (x->v1!=y->v1) return x->v1<y->v1?-1:1;
    return 0;
}
/* 1 if every undirected edge is in <=2 faces AND no face is degenerate. */
static int manifold_ok(const int32_t *f, size_t nf)
{
    size_t i, n=nf*3;
    UE *e; int ok=1;
    for (i=0;i<nf;i++){
        int32_t a=f[i*3+0], b=f[i*3+1], c=f[i*3+2];
        if (a==b||b==c||a==c) return 0;   /* degenerate face */
    }
    if (nf==0) return 1;
    e=(UE*)malloc(n*sizeof(UE));
    for (i=0;i<nf;i++){
        int32_t t[3]={f[i*3+0],f[i*3+1],f[i*3+2]}; int k;
        for (k=0;k<3;k++){
            int32_t a=t[k], b=t[(k+1)%3];
            if (a>b){ int32_t s=a; a=b; b=s; }
            e[i*3+(size_t)k].v0=a; e[i*3+(size_t)k].v1=b;
        }
    }
    qsort(e, n, sizeof(UE), ue_cmp);
    i=0;
    while (i<n){
        size_t j=i+1;
        while (j<n && e[j].v0==e[i].v0 && e[j].v1==e[i].v1) j++;
        if (j-i > 2){ ok=0; break; }
        i=j;
    }
    free(e);
    return ok;
}
static int uses_vert(const int32_t *f, size_t nf, int32_t v)
{
    size_t i; for (i=0;i<nf*3;i++) if (f[i]==v) return 1; return 0;
}

/* ---- mesh builders (arena-allocated verts/faces) ---- */
static void set_v(float *V, int i, float z, float y, float x)
{ V[i*3+0]=z; V[i*3+1]=y; V[i*3+2]=x; }

/* Hexagon fan split into two near-coincident centres C0,C1 (0.05 apart) bridged
 * by two needle triangles. r0..r5 ring (boundary), C0=6 C1=7 (interior). 8 faces,
 * 2 of them slivers (the bridges). Coords in the z=0 plane (here "z y x" = x,y,0
 * conceptually; only relative geometry matters). */
static void build_fan_needle(Arena_T arena, ComponentMesh *cm)
{
    float *V=(float*)ARENA_ALLOC(arena,(long)(8*3*sizeof(float)));
    int32_t *F=(int32_t*)ARENA_ALLOC(arena,(long)(8*3*sizeof(int32_t)));
    int32_t faces[8*3] = {
        6,3,4,  6,4,5,  6,5,0,     /* C0 fan (lower arc r3..r0) */
        7,0,1,  7,1,2,  7,2,3,     /* C1 fan (upper arc r0..r3) */
        6,0,7,  7,3,6              /* bridges (needles) at r0 and r3 */
    };
    set_v(V,0, 5.0f,        0.0f,      0.0f);   /* r0   0deg  */
    set_v(V,1, 2.5f,        4.330127f, 0.0f);   /* r1  60deg  */
    set_v(V,2,-2.5f,        4.330127f, 0.0f);   /* r2 120deg  */
    set_v(V,3,-5.0f,        0.0f,      0.0f);   /* r3 180deg  */
    set_v(V,4,-2.5f,       -4.330127f, 0.0f);   /* r4 240deg  */
    set_v(V,5, 2.5f,       -4.330127f, 0.0f);   /* r5 300deg  */
    set_v(V,6, 0.0f,        0.0f,      0.0f);   /* C0 centre  */
    set_v(V,7, 0.0f,        0.05f,     0.0f);   /* C1 centre (0.05 from C0) */
    memcpy(F, faces, sizeof faces);
    memset(cm,0,sizeof *cm);
    cm->verts=V; cm->faces=F; cm->nv=8; cm->nf=8; cm->comp_id=1; cm->self=cm;
}

/* Clean regular hexagon fan: one centre (6), 6 ring (boundary), 6 faces, no
 * slivers. */
static void build_fan_clean(Arena_T arena, ComponentMesh *cm)
{
    float *V=(float*)ARENA_ALLOC(arena,(long)(7*3*sizeof(float)));
    int32_t *F=(int32_t*)ARENA_ALLOC(arena,(long)(6*3*sizeof(int32_t)));
    int32_t faces[6*3] = { 6,0,1, 6,1,2, 6,2,3, 6,3,4, 6,4,5, 6,5,0 };
    set_v(V,0, 5.0f,        0.0f,      0.0f);
    set_v(V,1, 2.5f,        4.330127f, 0.0f);
    set_v(V,2,-2.5f,        4.330127f, 0.0f);
    set_v(V,3,-5.0f,        0.0f,      0.0f);
    set_v(V,4,-2.5f,       -4.330127f, 0.0f);
    set_v(V,5, 2.5f,       -4.330127f, 0.0f);
    set_v(V,6, 0.0f,        0.0f,      0.0f);
    memcpy(F, faces, sizeof faces);
    memset(cm,0,sizeof *cm);
    cm->verts=V; cm->faces=F; cm->nv=7; cm->nf=6; cm->comp_id=1; cm->self=cm;
}

/* Octahedron with equator verts e0,e1 pulled close (short interior edge e0-e1).
 * top=0 bot=1 e0=2 e1=3 e2=4 e3=5. 8 faces, closed, all verts interior. The two
 * faces on edge e0-e1 are slivers; collapsing e0-e1 is link-valid (common nbrs
 * = {top,bot} = the apexes). */
static void build_octa_sliver(Arena_T arena, ComponentMesh *cm)
{
    float *V=(float*)ARENA_ALLOC(arena,(long)(6*3*sizeof(float)));
    int32_t *F=(int32_t*)ARENA_ALLOC(arena,(long)(8*3*sizeof(int32_t)));
    int32_t faces[8*3] = {
        0,2,3, 0,3,4, 0,4,5, 0,5,2,    /* top fan over equator */
        1,3,2, 1,4,3, 1,5,4, 1,2,5     /* bottom fan (reversed) */
    };
    set_v(V,0, 0.0f, 0.0f, 5.0f);    /* top */
    set_v(V,1, 0.0f, 0.0f,-5.0f);    /* bot */
    set_v(V,2, 5.0f, 0.0f, 0.0f);    /* e0 */
    set_v(V,3, 5.0f, 0.05f,0.0f);    /* e1 (needle-close to e0: edge 0.05) */
    set_v(V,4,-5.0f, 0.0f, 0.0f);    /* e2 */
    set_v(V,5, 0.0f,-5.0f, 0.0f);    /* e3 */
    memcpy(F, faces, sizeof faces);
    memset(cm,0,sizeof *cm);
    cm->verts=V; cm->faces=F; cm->nv=6; cm->nf=8; cm->comp_id=1; cm->self=cm;
}

/* Triangular bipyramid: a=0 b=1 c=2 top=3 bot=4. Edge a-b is shared by
 * (top,a,b),(bot,b,a); c is a common neighbour of a,b but NOT an apex of that
 * edge -> collapsing a-b is non-manifold and the link condition must refuse it.
 * a,b set close so the two a-b faces register as slivers. */
static void build_bipyramid(Arena_T arena, ComponentMesh *cm)
{
    float *V=(float*)ARENA_ALLOC(arena,(long)(5*3*sizeof(float)));
    int32_t *F=(int32_t*)ARENA_ALLOC(arena,(long)(6*3*sizeof(int32_t)));
    int32_t faces[6*3] = {
        3,0,1, 3,1,2, 3,2,0,    /* top fan */
        4,1,0, 4,2,1, 4,0,2     /* bottom fan */
    };
    set_v(V,0, 0.0f, 0.0f, 0.0f);    /* a */
    set_v(V,1, 0.05f,0.0f, 0.0f);    /* b (a-b = 0.05) */
    set_v(V,2, 2.5f, 4.0f, 0.0f);    /* c */
    set_v(V,3, 1.5f, 1.5f, 3.0f);    /* top */
    set_v(V,4, 1.5f, 1.5f,-3.0f);    /* bot */
    memcpy(F, faces, sizeof faces);
    memset(cm,0,sizeof *cm);
    cm->verts=V; cm->faces=F; cm->nv=5; cm->nf=6; cm->comp_id=1; cm->self=cm;
}

#define CHECK(cond, msg) do { if (!(cond)) { printf("    FAIL: %s\n", msg); fail++; } } while (0)

int main(void)
{
    int fail = 0;
    Arena_T arena = Arena_new();

    /* G1: end-to-end. */
    {
        ComponentMesh cm; WeldCleanupParams p; WeldCleanupStats s;
        Arena_Mark m = Arena_save(arena);
        build_fan_needle(arena, &cm);
        WeldCleanup_default_params(&p);
        WeldCleanup_process(arena, &cm, &p, &s);
        printf("  G1 end-to-end: in=%zu out=%zu flips=%zu coll=%zu targets %zu->%zu\n",
               s.faces_in, s.faces_out, s.n_flips, s.n_collapses, s.targets_in, s.targets_out);
        CHECK(s.targets_in == 2, "G1 expected 2 sliver targets in");
        CHECK(s.targets_out == 0, "G1 slivers not all removed");
        CHECK(manifold_ok(cm.faces, cm.nf), "G1 result not manifold/degenerate-free");
        CHECK(s.n_collapses >= 1, "G1 expected at least one collapse (C0-C1)");
        /* boundary ring verts must still be present (never moved/removed) */
        CHECK(uses_vert(cm.faces, cm.nf, 0) && uses_vert(cm.faces, cm.nf, 3),
              "G1 boundary ring vertex dropped");
        Arena_restore(arena, m);
    }

    /* G2: collapse isolated (no flips). */
    {
        ComponentMesh cm; WeldCleanupParams p; WeldCleanupStats s;
        Arena_Mark m = Arena_save(arena);
        build_fan_needle(arena, &cm);
        WeldCleanup_default_params(&p);
        p.flip_max_rounds = 0;
        WeldCleanup_process(arena, &cm, &p, &s);
        printf("  G2 collapse-only: out=%zu flips=%zu coll=%zu targets->%zu\n",
               s.faces_out, s.n_flips, s.n_collapses, s.targets_out);
        CHECK(s.n_flips == 0, "G2 flips should be disabled");
        CHECK(s.n_collapses >= 1, "G2 collapse did not fire");
        CHECK(s.faces_out == 6, "G2 expected 8->6 faces");
        CHECK(s.targets_out == 0, "G2 slivers remain");
        CHECK(manifold_ok(cm.faces, cm.nf), "G2 result not manifold");
        Arena_restore(arena, m);
    }

    /* G3: length guard blocks the collapse (wrap-safety). */
    {
        ComponentMesh cm; WeldCleanupParams p; WeldCleanupStats s;
        Arena_Mark m = Arena_save(arena);
        build_fan_needle(arena, &cm);
        WeldCleanup_default_params(&p);
        p.flip_max_rounds = 0;
        p.max_collapse_len = 0.01;   /* below the 0.05 needle edge -> refuse */
        WeldCleanup_process(arena, &cm, &p, &s);
        printf("  G3 length-guard: out=%zu coll=%zu targets->%zu\n",
               s.faces_out, s.n_collapses, s.targets_out);
        CHECK(s.n_collapses == 0, "G3 collapse should be blocked by length cap");
        CHECK(s.faces_out == 8, "G3 nothing should be removed");
        CHECK(s.targets_out == 2, "G3 slivers should remain");
        Arena_restore(arena, m);
    }

    /* G4: clean mesh untouched. */
    {
        ComponentMesh cm; WeldCleanupParams p; WeldCleanupStats s;
        Arena_Mark m = Arena_save(arena);
        build_fan_clean(arena, &cm);
        WeldCleanup_default_params(&p);
        WeldCleanup_process(arena, &cm, &p, &s);
        printf("  G4 clean: out=%zu flips=%zu coll=%zu targets %zu->%zu\n",
               s.faces_out, s.n_flips, s.n_collapses, s.targets_in, s.targets_out);
        CHECK(s.targets_in == 0, "G4 clean mesh should have no targets");
        CHECK(s.n_collapses == 0, "G4 should not collapse a clean mesh");
        CHECK(s.faces_out == 6, "G4 face count must be unchanged");
        CHECK(manifold_ok(cm.faces, cm.nf), "G4 result not manifold");
        Arena_restore(arena, m);
    }

    /* G5: trivial inputs. */
    {
        ComponentMesh cm; WeldCleanupParams p;
        Arena_Mark m = Arena_save(arena);
        int32_t one[3] = {0,1,2};
        float v3[9] = {0,0,0, 1,0,0, 0,1,0};
        WeldCleanup_default_params(&p);
        memset(&cm,0,sizeof cm); cm.self=&cm;       /* empty */
        CHECK(WeldCleanup_process(arena, &cm, &p, NULL) == 0, "G5 empty crashed");
        memset(&cm,0,sizeof cm); cm.verts=v3; cm.faces=one; cm.nv=3; cm.nf=1; cm.self=&cm;
        CHECK(WeldCleanup_process(arena, &cm, &p, NULL) == 0, "G5 single-tri crashed");
        CHECK(cm.nf == 1, "G5 single triangle should survive");
        printf("  G5 trivial: empty + single-triangle handled, nf=%zu\n", cm.nf);
        Arena_restore(arena, m);
    }

    /* G6: closed-mesh interior collapse (link ACCEPT). */
    {
        ComponentMesh cm; WeldCleanupParams p; WeldCleanupStats s;
        Arena_Mark m = Arena_save(arena);
        build_octa_sliver(arena, &cm);
        WeldCleanup_default_params(&p);
        p.flip_max_rounds = 0;       /* isolate collapse on the closed mesh */
        WeldCleanup_process(arena, &cm, &p, &s);
        printf("  G6 octa: out=%zu coll=%zu targets %zu->%zu manifold=%d\n",
               s.faces_out, s.n_collapses, s.targets_in, s.targets_out,
               manifold_ok(cm.faces, cm.nf));
        CHECK(s.targets_in >= 1, "G6 expected sliver faces on the short edge");
        CHECK(s.n_collapses >= 1, "G6 valid interior collapse should fire");
        CHECK(s.faces_out == 6, "G6 octahedron 8->6 after collapsing e0-e1");
        CHECK(manifold_ok(cm.faces, cm.nf), "G6 closed mesh broke manifoldness");
        Arena_restore(arena, m);
    }

    /* G7: link-condition REJECT (would-be non-manifold). */
    {
        ComponentMesh cm; WeldCleanupParams p; WeldCleanupStats s;
        Arena_Mark m = Arena_save(arena);
        build_bipyramid(arena, &cm);
        WeldCleanup_default_params(&p);
        p.flip_max_rounds = 0;       /* force the collapse path */
        WeldCleanup_process(arena, &cm, &p, &s);
        printf("  G7 bipyramid: out=%zu coll=%zu targets %zu->%zu manifold=%d\n",
               s.faces_out, s.n_collapses, s.targets_in, s.targets_out,
               manifold_ok(cm.faces, cm.nf));
        CHECK(s.targets_in >= 1, "G7 expected a sliver on a-b");
        CHECK(s.n_collapses == 0, "G7 link condition must REFUSE the a-b collapse");
        CHECK(s.faces_out == 6, "G7 mesh must be left intact");
        CHECK(manifold_ok(cm.faces, cm.nf), "G7 mesh must stay manifold");
        Arena_restore(arena, m);
    }

    Arena_dispose(&arena);
    printf("WELD CLEANUP SELFTEST %s\n", fail ? "FAILED" : "PASSED");
    return fail ? 1 : 0;
}
