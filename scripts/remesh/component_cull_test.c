/* component_cull_test.c -- unit tests for ComponentCull_by_area (kibble removal).
 *
 * Builds a mesh holding a big sheet + a tiny disconnected sheet (2 connectivity-
 * components) and checks the area filter drops the tiny one; that a small but
 * independently supplied input sheet is NOT compared with the whole-cube total;
 * and that nothing is culled when all
 * components clear the threshold. Exit 0 = pass.
 */
#include "remesh/component_cull.h"
#include "common/arena.h"
#include "common/mesh_types.h"
#include "common/pipeline_constants.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int g_fails = 0;
#define CHECK(cond, msg) do { \
    if (!(cond)) { printf("  FAIL: %s\n", (msg)); g_fails++; } \
    else         { printf("  ok:   %s\n", (msg)); } \
} while (0)

/* Fill an R x C unit-spaced flat grid at z=zc into V/F at offsets. Area = (R-1)*(C-1). */
static void build_grid(float *V, int32_t *F, int voff, int foff,
                       int R, int C, float zc, size_t *pnv, size_t *pnf)
{
    int nv = 0;
    for (int r = 0; r < R; r++) for (int c = 0; c < C; c++) {
        int i = voff + nv;
        V[i*3+0]=zc; V[i*3+1]=(float)r; V[i*3+2]=(float)c; nv++;
    }
    int nf = 0;
    for (int r = 0; r < R-1; r++) for (int c = 0; c < C-1; c++) {
        int a=voff+r*C+c, b=voff+r*C+(c+1), d=voff+(r+1)*C+c, e=voff+(r+1)*C+(c+1);
        F[(foff+nf)*3+0]=a; F[(foff+nf)*3+1]=b; F[(foff+nf)*3+2]=d; nf++;
        F[(foff+nf)*3+0]=b; F[(foff+nf)*3+1]=e; F[(foff+nf)*3+2]=d; nf++;
    }
    *pnv=(size_t)nv; *pnf=(size_t)nf;
}

int main(void)
{
    printf("=== component_cull_test ===\n");
    Arena_T arena = Arena_new();

    /* big 16x16 (area 225) + tiny 3x3 (area 4): 4/229 = 1.7% < 2% -> tiny is kibble */
    const int BR=16, BC=16, TR=3, TC=3;
    size_t b_nv=(size_t)BR*BC, b_nf=(size_t)2*(BR-1)*(BC-1);
    size_t t_nv=(size_t)TR*TC, t_nf=(size_t)2*(TR-1)*(TC-1);
    size_t nv=b_nv+t_nv, nf=b_nf+t_nf;
    float   *V=(float*)malloc(nv*3*sizeof(float));
    int32_t *F=(int32_t*)malloc(nf*3*sizeof(int32_t));
    size_t x0,x1,x2,x3;
    build_grid(V,F, 0,            0,            BR,BC, 0.0f,  &x0,&x1);
    build_grid(V,F, (int)b_nv,    (int)b_nf,    TR,TC, 50.0f, &x2,&x3);

    /* ---- test 1: one mesh, two CCs -> tiny culled ---- */
    printf("[one mesh, big + tiny CC]\n");
    {
        ComponentMesh in; memset(&in,0,sizeof in);
        in.verts=V; in.faces=F; in.nv=nv; in.nf=nf; in.comp_id=1; in.self=&in;
        ComponentMesh *out=NULL; size_t n_out=0, n_cull=0;
        int rc = ComponentCull_by_area(arena, &in, 1, KIBBLE_AREA_FRAC, &out, &n_out, &n_cull);
        printf("    n_out=%zu n_culled=%zu\n", n_out, n_cull);
        CHECK(rc==0, "ComponentCull_by_area returns 0");
        CHECK(n_out==1 && n_cull==1, "tiny CC (<2% area) is culled, big kept");
        CHECK(n_out==1 && out[0].nv==b_nv, "the surviving component is the big sheet");
    }

    /* ---- test 2: two semantic input sheets -> each parent's main CC kept ---- */
    printf("[two parent sheets, big + small]\n");
    {
        ComponentMesh m[2]; memset(m,0,sizeof m);
        m[0].verts=V;            m[0].faces=F;          m[0].nv=b_nv; m[0].nf=b_nf; m[0].self=&m[0];
        int32_t *Ft=(int32_t*)malloc(t_nf*3*sizeof(int32_t));
        for(size_t i=0;i<t_nf*3;i++) Ft[i]=F[b_nf*3+i]-(int32_t)b_nv;
        m[1].verts=V+b_nv*3;     m[1].faces=Ft;         m[1].nv=t_nv; m[1].nf=t_nf; m[1].self=&m[1];
        ComponentMesh *out=NULL; size_t n_out=0, n_cull=0;
        int rc = ComponentCull_by_area(arena, m, 2, KIBBLE_AREA_FRAC, &out, &n_out, &n_cull);
        CHECK(rc==0 && n_out==2 && n_cull==0,
              "small semantic input sheet survives parent-relative cull");
        free(Ft);
    }

    /* ---- test 3: big alone -> nothing culled ---- */
    printf("[big alone]\n");
    {
        ComponentMesh in; memset(&in,0,sizeof in);
        in.verts=V; in.faces=F; in.nv=b_nv; in.nf=b_nf; in.comp_id=1; in.self=&in;
        ComponentMesh *out=NULL; size_t n_out=0, n_cull=0;
        int rc = ComponentCull_by_area(arena, &in, 1, KIBBLE_AREA_FRAC, &out, &n_out, &n_cull);
        CHECK(rc==0 && n_out==1 && n_cull==0, "single sheet: nothing culled");
    }

    Arena_dispose(&arena);
    free(V); free(F);
    printf("=== %s (%d failure%s) ===\n", g_fails==0?"PASS":"FAIL", g_fails, g_fails==1?"":"s");
    return g_fails==0 ? 0 : 1;
}
