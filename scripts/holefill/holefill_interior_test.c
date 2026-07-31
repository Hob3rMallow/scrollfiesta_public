/*
 * holefill_interior_test.c -- unit tests for HoleFill_process_ex(interior_only).
 *
 * Validates the GEOMETRIC interior/exterior classifier (signed area vs surface
 * normal): interior holes get filled, outer perimeters and open bays do not --
 * even when the perimeter is small enough to pass the diameter / vertex caps, so
 * the signed-area test is the only thing preventing the sheet from being capped.
 *
 *   G1  flat sheet (7x7) with one interior square hole: the hole is filled, the
 *       outer perimeter is left open (boundary edges 28 -> 24). Manifold.
 *   G2  flat sheet, NO hole: nothing classified interior, faces unchanged, the
 *       perimeter is never filled.
 *   G3  TWO disjoint sheets (5x5), each with an interior hole: BOTH holes filled,
 *       BOTH perimeters left open. This is the multi-component case the old
 *       "skip the single largest loop" heuristic gets wrong; the geometric test
 *       recognises each component's own perimeter.
 *
 * Build: holefill_interior_test.vcxproj. Run: exit 0 = pass, 1 = a gate failed.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/arena.h"
#include "common/obj_io.h"
#include "holefill/hole_fill.h"

/* Exposed (non-static) from hole_fill.c for direct unit testing -- no public
 * prototype, so declare it here. Decomposes a self-pinching "<><>" boundary
 * loop into simple sub-cycles. */
extern size_t split_pinched_loop(Arena_T arena,
                                 int32_t *loop_verts, size_t n,
                                 const float *verts, float eps,
                                 int32_t **out_sub, size_t *out_len,
                                 size_t max_subs);

/* Exposed (non-static) from hole_fill.c for direct unit testing. Collapses
 * near-coincident interior (Steiner) verts of a freshly-built fill and drops the
 * resulting needle triangles; commits only if the result stays manifold. */
extern int prune_degenerate_fill(float *verts, int32_t *faces,
                                 size_t *p_nv, size_t *p_nf,
                                 size_t n_boundary,
                                 float eps_vox, float area_eps);

/* ---- boundary-edge count + manifold check on a flat face list ---- */
typedef struct { int32_t v0, v1; } UE;
static int ue_cmp(const void *a, const void *b)
{
    const UE *x=(const UE*)a, *y=(const UE*)b;
    if (x->v0!=y->v0) return x->v0<y->v0?-1:1;
    if (x->v1!=y->v1) return x->v1<y->v1?-1:1;
    return 0;
}
/* number of undirected edges that appear in exactly one face */
static size_t count_boundary_edges(const int32_t *f, size_t nf)
{
    size_t n=nf*3, i, nb=0;
    UE *e;
    if (nf==0) return 0;
    e=(UE*)malloc(n*sizeof(UE));
    for (i=0;i<nf;i++){
        int32_t t[3]={f[i*3+0],f[i*3+1],f[i*3+2]}; int k;
        for (k=0;k<3;k++){
            int32_t a=t[k], b=t[(k+1)%3];
            if (a>b){ int32_t s=a; a=b; b=s; }
            e[i*3+(size_t)k].v0=a; e[i*3+(size_t)k].v1=b;
        }
    }
    qsort(e,n,sizeof(UE),ue_cmp);
    i=0;
    while (i<n){ size_t j=i+1; while (j<n && e[j].v0==e[i].v0 && e[j].v1==e[i].v1) j++;
                 if (j-i==1) nb++; i=j; }
    free(e);
    return nb;
}
static int manifold_ok(const int32_t *f, size_t nf)
{
    size_t n=nf*3, i; UE *e; int ok=1;
    for (i=0;i<nf;i++){ int32_t a=f[i*3+0],b=f[i*3+1],c=f[i*3+2];
                        if (a==b||b==c||a==c) return 0; }
    if (nf==0) return 1;
    e=(UE*)malloc(n*sizeof(UE));
    for (i=0;i<nf;i++){ int32_t t[3]={f[i*3+0],f[i*3+1],f[i*3+2]}; int k;
        for (k=0;k<3;k++){ int32_t a=t[k],b=t[(k+1)%3]; if(a>b){int32_t s=a;a=b;b=s;}
                           e[i*3+(size_t)k].v0=a; e[i*3+(size_t)k].v1=b; } }
    qsort(e,n,sizeof(UE),ue_cmp);
    i=0;
    while (i<n){ size_t j=i+1; while (j<n && e[j].v0==e[i].v0 && e[j].v1==e[i].v1) j++;
                 if (j-i>2){ ok=0; break; } i=j; }
    free(e);
    return ok;
}

/* Build an NxN vertex grid in the z=0 plane (verts "z y x" = 0, j, i+xoff),
 * triangulated CCW; optionally drop the center cell to make a 4-vert interior
 * hole. Appends into V/F at *pnv/*pnf, returns via updated counts. */
static void add_grid(float *V, int32_t *F, size_t *pnv, size_t *pnf,
                     int N, float xoff, int hole)
{
    size_t base = *pnv;
    int i, j;
    for (i=0;i<N;i++) for (j=0;j<N;j++){
        size_t idx = base + (size_t)(i*N+j);
        V[idx*3+0]=0.0f; V[idx*3+1]=(float)j; V[idx*3+2]=(float)i + xoff;
    }
    int hc = N/2;   /* center cell (hc,hc) */
    for (i=0;i<N-1;i++) for (j=0;j<N-1;j++){
        if (hole && i==hc && j==hc) continue;   /* leave a hole */
        int32_t a=(int32_t)(base+(size_t)(i*N+j));
        int32_t b=(int32_t)(base+(size_t)((i+1)*N+j));
        int32_t c=(int32_t)(base+(size_t)((i+1)*N+(j+1)));
        int32_t d=(int32_t)(base+(size_t)(i*N+(j+1)));
        F[*pnf*3+0]=a; F[*pnf*3+1]=b; F[*pnf*3+2]=c; (*pnf)++;
        F[*pnf*3+0]=a; F[*pnf*3+1]=c; F[*pnf*3+2]=d; (*pnf)++;
    }
    *pnv = base + (size_t)(N*N);
}

/* Build an NxN vertex grid in the z=0 plane (as add_grid) but drop a whole
 * [hlo,hhi) x [hlo,hhi) BLOCK of cells, leaving one large square interior hole.
 * The hole's boundary loop has 4*(hhi-hlo) verts and spans (hhi-hlo) voxels --
 * sized by the caller to exceed the old MAX_LOOP_VERTS / diameter caps so the
 * removal of those caps (interior_only) is what lets it fill. */
static void add_grid_block_hole(float *V, int32_t *F, size_t *pnv, size_t *pnf,
                                int N, float xoff, int hlo, int hhi)
{
    size_t base = *pnv;
    int i, j;
    for (i=0;i<N;i++) for (j=0;j<N;j++){
        size_t idx = base + (size_t)(i*N+j);
        V[idx*3+0]=0.0f; V[idx*3+1]=(float)j; V[idx*3+2]=(float)i + xoff;
    }
    for (i=0;i<N-1;i++) for (j=0;j<N-1;j++){
        if (i>=hlo && i<hhi && j>=hlo && j<hhi) continue;   /* drop the block */
        int32_t a=(int32_t)(base+(size_t)(i*N+j));
        int32_t b=(int32_t)(base+(size_t)((i+1)*N+j));
        int32_t c=(int32_t)(base+(size_t)((i+1)*N+(j+1)));
        int32_t d=(int32_t)(base+(size_t)(i*N+(j+1)));
        F[*pnf*3+0]=a; F[*pnf*3+1]=b; F[*pnf*3+2]=c; (*pnf)++;
        F[*pnf*3+0]=a; F[*pnf*3+1]=c; F[*pnf*3+2]=d; (*pnf)++;
    }
    *pnv = base + (size_t)(N*N);
}

#define CHECK(cond,msg) do { if(!(cond)){ printf("    FAIL: %s\n", msg); fail++; } } while(0)

int main(int argc, char **argv)
{
    int fail = 0;
    Arena_T arena = Arena_new();

    /* G1: single sheet, one interior hole. */
    {
        Arena_Mark m = Arena_save(arena);
        float *V = (float*)ARENA_ALLOC(arena, (long)(64*3*sizeof(float)));
        int32_t *F = (int32_t*)ARENA_ALLOC(arena, (long)(128*3*sizeof(int32_t)));
        size_t nv=0, nf=0;
        add_grid(V, F, &nv, &nf, 7, 0.0f, 1);
        size_t bnd_before = count_boundary_edges(F, nf);
        size_t nl=0, ni=0, nfil=0;
        HoleFill_process_ex(arena, &V, &F, &nv, &nf, NULL, 1, &nl, &ni, &nfil);
        size_t bnd_after = count_boundary_edges(F, nf);
        printf("  G1 sheet+hole: bnd %zu->%zu  loops=%zu interior=%zu filled=%zu\n",
               bnd_before, bnd_after, nl, ni, nfil);
        CHECK(bnd_before == 28, "G1 expected 24 perimeter + 4 hole boundary edges");
        CHECK(ni == 1, "G1 should classify exactly one interior hole");
        CHECK(nfil == 1, "G1 interior hole should be filled");
        CHECK(bnd_after == 24, "G1 hole closed but perimeter must stay open (24)");
        CHECK(manifold_ok(F, nf), "G1 result not manifold");
        Arena_restore(arena, m);
    }

    /* G2: single sheet, NO hole -> nothing interior, perimeter untouched. */
    {
        Arena_Mark m = Arena_save(arena);
        float *V = (float*)ARENA_ALLOC(arena, (long)(64*3*sizeof(float)));
        int32_t *F = (int32_t*)ARENA_ALLOC(arena, (long)(128*3*sizeof(int32_t)));
        size_t nv=0, nf=0;
        add_grid(V, F, &nv, &nf, 7, 0.0f, 0);
        size_t nf_before = nf, bnd_before = count_boundary_edges(F, nf);
        size_t nl=0, ni=0, nfil=0;
        HoleFill_process_ex(arena, &V, &F, &nv, &nf, NULL, 1, &nl, &ni, &nfil);
        printf("  G2 sheet no-hole: bnd=%zu loops=%zu interior=%zu filled=%zu faces %zu->%zu\n",
               bnd_before, nl, ni, nfil, nf_before, nf);
        CHECK(ni == 0, "G2 no interior holes expected");
        CHECK(nfil == 0, "G2 nothing should be filled (perimeter is not interior)");
        CHECK(nf == nf_before, "G2 face count must be unchanged");
        CHECK(count_boundary_edges(F, nf) == bnd_before, "G2 perimeter must stay open");
        Arena_restore(arena, m);
    }

    /* G3: TWO disjoint sheets, each with a hole (multi-component). */
    {
        Arena_Mark m = Arena_save(arena);
        float *V = (float*)ARENA_ALLOC(arena, (long)(64*3*sizeof(float)));
        int32_t *F = (int32_t*)ARENA_ALLOC(arena, (long)(128*3*sizeof(int32_t)));
        size_t nv=0, nf=0;
        add_grid(V, F, &nv, &nf, 5, 0.0f,   1);   /* sheet A */
        add_grid(V, F, &nv, &nf, 5, 100.0f, 1);   /* sheet B, far away */
        size_t bnd_before = count_boundary_edges(F, nf);
        size_t nl=0, ni=0, nfil=0;
        HoleFill_process_ex(arena, &V, &F, &nv, &nf, NULL, 1, &nl, &ni, &nfil);
        size_t bnd_after = count_boundary_edges(F, nf);
        printf("  G3 two sheets+holes: bnd %zu->%zu  loops=%zu interior=%zu filled=%zu\n",
               bnd_before, bnd_after, nl, ni, nfil);
        CHECK(bnd_before == 40, "G3 expected 2*(16 perimeter + 4 hole)");
        CHECK(ni == 2, "G3 should classify BOTH holes interior (multi-component)");
        CHECK(nfil == 2, "G3 both holes should be filled");
        CHECK(bnd_after == 32, "G3 both holes closed, BOTH perimeters open (2*16)");
        CHECK(manifold_ok(F, nf), "G3 result not manifold");
        Arena_restore(arena, m);
    }

    /* G4: pinch-split. A figure-8 "<><>" loop -- two 4-cycles joined at a
     * coincident-but-DISTINCT vertex (verts 0 and 4 share position (1,1), as a
     * weld produces from two cubes' seam verts). split_pinched_loop must peel it
     * into two simple 4-cycles (this is the exact config that hung Triangle). */
    {
        Arena_Mark m = Arena_save(arena);
        /* verts as "z y x"; z=0. slot4 coincides with slot0 at (1,1). */
        float V[8*3] = {
            0,1,1,   0,0,1,   0,0,0,   0,1,0,    /* square A around (1,1) */
            0,1,1,   0,1,2,   0,2,2,   0,2,1     /* square B, vert4==(1,1) */
        };
        int32_t loop[8] = {0,1,2,3,4,5,6,7};
        int32_t *subs[8]; size_t sublen[8];
        size_t ns = split_pinched_loop(arena, loop, 8, V, 0.05f, subs, sublen, 8);
        printf("  G4 figure-8 pinch: split into %zu sub-loops (len %zu, %zu)\n",
               ns, ns>0?sublen[0]:0, ns>1?sublen[1]:0);
        CHECK(ns == 2, "G4 figure-8 must split into exactly 2 sub-loops");
        CHECK(ns == 2 && sublen[0] == 4 && sublen[1] == 4,
              "G4 each sub-loop must be a 4-cycle");
        /* A non-pinched simple loop must pass through unchanged (alias, ns==1). */
        int32_t simple[4] = {0,1,2,3};
        int32_t *ss[4]; size_t sl[4];
        size_t ns2 = split_pinched_loop(arena, simple, 4, V, 0.05f, ss, sl, 4);
        CHECK(ns2 == 1 && sl[0] == 4 && ss[0] == simple,
              "G4 simple loop must pass through unchanged (1 sub-loop, aliased)");
        Arena_restore(arena, m);
    }

    /* G5: ONE LARGE interior hole that exceeds BOTH old caps -- a 140x140 sheet
     * with a 132x132 cell block dropped from the middle (border 4 cells wide).
     * The hole boundary is 4*132 = 528 verts (> MAX_LOOP_VERTS 500) and spans
     * ~132 voxels (> HOLEFILL_MAX_DIAM_VOX 48). The OLD code skipped it as "too
     * large"/"diameter" and left the hole gaping (exactly the comp-007 bug); the
     * interior_only path now size-gates nothing, so the winding test fills it and
     * the single-largest skip leaves the outer perimeter open. */
    {
        Arena_Mark m = Arena_save(arena);
        const int N = 140, HLO = 4, HHI = 136;   /* hole side = 132 cells */
        float   *V = (float*)  ARENA_ALLOC(arena, (long)((size_t)N*N*3*sizeof(float)));
        int32_t *F = (int32_t*)ARENA_ALLOC(arena, (long)((size_t)N*N*2*3*sizeof(int32_t)));
        size_t nv=0, nf=0;
        add_grid_block_hole(V, F, &nv, &nf, N, 0.0f, HLO, HHI);
        size_t bnd_before = count_boundary_edges(F, nf);
        size_t outer = (size_t)(4*(N-1));          /* 556 */
        size_t inner = (size_t)(4*(HHI-HLO));      /* 528 */
        size_t nl=0, ni=0, nfil=0;
        HoleFill_process_ex(arena, &V, &F, &nv, &nf, NULL, 1, &nl, &ni, &nfil);
        size_t bnd_after = count_boundary_edges(F, nf);
        printf("  G5 big interior hole (%zu-vert boundary): bnd %zu->%zu  "
               "loops=%zu interior=%zu filled=%zu\n",
               inner, bnd_before, bnd_after, nl, ni, nfil);
        CHECK(bnd_before == outer + inner, "G5 expected outer perimeter + big hole");
        CHECK(inner > 500, "G5 hole must exceed the old vertex cap to be meaningful");
        CHECK(ni == 1, "G5 should classify exactly one interior hole (the big one)");
        CHECK(nfil == 1, "G5 big interior hole MUST be filled (size cap removed)");
        CHECK(bnd_after == outer, "G5 big hole closed, outer perimeter left open");
        CHECK(manifold_ok(F, nf), "G5 result not manifold");
        Arena_restore(arena, m);
    }

    /* G6: prune_degenerate_fill collapses near-coincident interior (Steiner)
     * verts and drops the resulting needle triangles -- the thin-slit-fill
     * degeneracy seen on the Z=4480 weld seam (8 sliver tris, a 0.06-vox edge,
     * a duplicate Steiner pair). A manifold hexagon disk fanned around TWO
     * interior points 0.05 vox apart (sliver edge c0-c1) must return as a clean
     * 6-tri fan around ONE merged interior point: manifold, boundary untouched.
     * The negative control (interior points 1.0 vox apart) must not change. */
    {
        float V[8*3] = {
            0.f, 0.000f,  2.000f,   /* b0 */
            0.f, 1.732f,  1.000f,   /* b1 */
            0.f, 1.732f, -1.000f,   /* b2 */
            0.f, 0.000f, -2.000f,   /* b3 */
            0.f,-1.732f, -1.000f,   /* b4 */
            0.f,-1.732f,  1.000f,   /* b5 */
            0.f, 0.000f,  0.000f,   /* c0 = idx 6 (kept rep)  */
            0.f, 0.030f,  0.040f,   /* c1 = idx 7, 0.05 vox from c0 */
        };
        int32_t F[8*3] = {
            0,1,6, 1,2,6, 2,3,6,    /* c0 fan over b0..b3 */
            3,4,7, 4,5,7, 5,0,7,    /* c1 fan over b3..b0 */
            3,6,7, 0,7,6,           /* bridge across the c0-c1 seam (the slivers) */
        };
        size_t nv = 8, nf = 8;
        int chg = prune_degenerate_fill(V, F, &nv, &nf, 6, 0.25f, 1e-4f);
        printf("  G6 prune sliver fill: changed=%d  nv 8->%zu  nf 8->%zu  manifold=%s\n",
               chg, nv, nf, manifold_ok(F, nf) ? "yes" : "NO");
        CHECK(chg == 1, "G6 prune should fire on the 0.05-vox Steiner pair");
        CHECK(nv == 7, "G6 coincident Steiner verts should merge to one (8->7)");
        CHECK(nf == 6, "G6 two sliver triangles should drop (8->6)");
        CHECK(manifold_ok(F, nf), "G6 pruned fill must stay manifold");
        {
            const float B0[6*3] = {
                0.f,0.f,2.f, 0.f,1.732f,1.f, 0.f,1.732f,-1.f,
                0.f,0.f,-2.f, 0.f,-1.732f,-1.f, 0.f,-1.732f,1.f };
            int bnd_ok = 1;
            for (int i = 0; i < 6*3; i++) if (V[i] != B0[i]) bnd_ok = 0;
            CHECK(bnd_ok, "G6 boundary verts must not move");
        }
        CHECK(V[6*3+0]==0.f && V[6*3+1]==0.f && V[6*3+2]==0.f,
              "G6 merged Steiner should sit at the kept rep c0");

        /* negative control: interior points 1.0 vox apart -> nothing to collapse */
        {
            float V2[8*3] = {
                0.f, 0.000f,  2.000f,  0.f, 1.732f,  1.000f,  0.f, 1.732f, -1.000f,
                0.f, 0.000f, -2.000f,  0.f,-1.732f, -1.000f,  0.f,-1.732f,  1.000f,
                0.f,-0.350f,  0.000f,  0.f, 0.350f,  0.000f, /* c0,c1 ~0.7 vox apart */
            };
            int32_t F2[8*3] = {
                0,1,6, 1,2,6, 2,3,6,  3,4,7, 4,5,7, 5,0,7,  3,6,7, 0,7,6 };
            size_t nv2 = 8, nf2 = 8;
            int chg2 = prune_degenerate_fill(V2, F2, &nv2, &nf2, 6, 0.25f, 1e-4f);
            printf("  G6 negative control: changed=%d nv=%zu nf=%zu\n", chg2, nv2, nf2);
            CHECK(chg2 == 0, "G6 control: well-separated Steiner must NOT collapse");
            CHECK(nv2 == 8 && nf2 == 8, "G6 control: counts must be unchanged");
        }
    }

    /* Optional real-data check: `holefill_interior_test <in.obj> [out.obj]`
     * loads a single welded component (e.g. the dumped comp-007 with its giant
     * unfilled hole), runs the interior_only fill, and reports how many boundary
     * edges closed. Writes the filled mesh to out.obj if a second arg is given.
     * This is a diagnostic (not a hard gate) -- it never fails the selftest. */
    if (argc >= 2) {
        Arena_Mark m = Arena_save(arena);
        float *V = NULL; int32_t *F = NULL; size_t nv=0, nf=0;
        printf("\n  [real-data] loading %s ...\n", argv[1]);
        if (ObjIO_read(arena, argv[1], &V, &nv, &F, &nf) != 0) {
            printf("  [real-data] FAILED to read %s\n", argv[1]);
        } else {
            size_t bnd_before = count_boundary_edges(F, nf);
            size_t nv0 = nv, nf0 = nf, nl=0, ni=0, nfil=0;
            printf("  [real-data] in: %zu verts %zu faces, %zu boundary edges\n",
                   nv0, nf0, bnd_before);
            HoleFill_process_ex(arena, &V, &F, &nv, &nf, NULL, 1, &nl, &ni, &nfil);
            size_t bnd_after = count_boundary_edges(F, nf);
            printf("  [real-data] loops=%zu interior=%zu filled=%zu\n", nl, ni, nfil);
            printf("  [real-data] out: %zu verts %zu faces, %zu boundary edges "
                   "(closed %zu)\n", nv, nf, bnd_after,
                   bnd_before > bnd_after ? bnd_before - bnd_after : 0);
            printf("  [real-data] manifold=%s\n", manifold_ok(F, nf) ? "yes" : "NO");
            if (argc >= 3) {
                if (ObjIO_write(argv[2], V, nv, F, nf) == 0)
                    printf("  [real-data] wrote %s\n", argv[2]);
                else
                    printf("  [real-data] FAILED to write %s\n", argv[2]);
            }
        }
        Arena_restore(arena, m);
    }

    Arena_dispose(&arena);
    printf("HOLEFILL INTERIOR SELFTEST %s\n", fail ? "FAILED" : "PASSED");
    return fail ? 1 : 0;
}
