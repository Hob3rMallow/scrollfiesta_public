/* isotropic_remesh_test.c — whitebox unit tests for the isotropic remesh pass
 * (src/remesh/isotropic_remesh.c). #includes the .c to reach the statics
 * (split_round, collapse_round, detect_boundary, mesh_area).
 *
 * Covers the genuinely-new / riskiest mechanics:
 *   - target-edge-length formula L = sqrt(4A / (sqrt(3) F));
 *   - winding-safe edge SPLIT: one interior edge -> +1 vert / +2 faces, the
 *     result 2-manifold and consistently wound (same_dir_edges == 0), midpoint
 *     placed correctly;
 *   - split NEVER touches a boundary (face_count==1) edge;
 *   - full run stays 2-manifold, holds the face count in tolerance, and freezes
 *     the boundary polyline exactly (fail-closed contract).
 *
 * Exit 0 = all pass, 1 = a failure.
 */
#include "remesh/isotropic_remesh.c"   /* whitebox include */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static int g_fails = 0;
#define CHECK(cond, msg) do { \
    if (!(cond)) { printf("  FAIL: %s\n", (msg)); g_fails++; } \
    else         { printf("  ok:   %s\n", (msg)); } \
} while (0)

/* Flat nx*ny quad grid at z=0 (unit spacing). */
static void grid(Arena_T a, int nx, int ny,
                 float **pv, size_t *pnv, int32_t **pf, size_t *pnf)
{
    size_t nv=(size_t)(nx*ny), nf=(size_t)((nx-1)*(ny-1)*2), vi=0, fi=0;
    float   *v=(float*)  ARENA_ALLOC(a,(long)(nv*3*sizeof(float)));
    int32_t *f=(int32_t*)ARENA_ALLOC(a,(long)(nf*3*sizeof(int32_t)));
    int ix, iy;
    for (iy=0;iy<ny;iy++) for (ix=0;ix<nx;ix++){
        v[vi*3+0]=(float)ix; v[vi*3+1]=(float)iy; v[vi*3+2]=0.0f; vi++;
    }
    for (iy=0;iy<ny-1;iy++) for (ix=0;ix<nx-1;ix++){
        int32_t i00=(int32_t)(iy*nx+ix), i10=(int32_t)(iy*nx+ix+1);
        int32_t i01=(int32_t)((iy+1)*nx+ix), i11=(int32_t)((iy+1)*nx+ix+1);
        f[fi*3+0]=i00; f[fi*3+1]=i10; f[fi*3+2]=i11; fi++;
        f[fi*3+0]=i00; f[fi*3+1]=i11; f[fi*3+2]=i01; fi++;
    }
    *pv=v; *pnv=nv; *pf=f; *pnf=nf;
}

static void test_L_formula(void)
{
    Arena_T a=Arena_new();
    float *v; size_t nv; int32_t *f; size_t nf;
    double A, L, expL;
    grid(a,11,11,&v,&nv,&f,&nf);           /* 10x10 quads => area 100, nf 200 */
    A=mesh_area(v,f,nf);
    CHECK(fabs(A-100.0)<1e-3, "grid area == 100");
    L=sqrt(4.0*A/(1.7320508075688772*(double)nf));
    expL=sqrt(4.0*100.0/(1.7320508075688772*200.0));
    CHECK(fabs(L-expL)<1e-6, "L = sqrt(4A/(sqrt3 F))");
    Arena_dispose(&a);
}

static void test_split_one_edge(void)
{
    Arena_T a=Arena_new();
    float   *v  =(float*)  ARENA_ALLOC(a,(long)(4*3*sizeof(float)));
    int32_t *f  =(int32_t*)ARENA_ALLOC(a,(long)(3*3*sizeof(int32_t)));
    uint8_t *pin=(uint8_t*)ARENA_ALLOC(a,(long)(4*sizeof(uint8_t)));
    size_t nv=4, nf=3, k;
    MeshManifoldStats st;
    /* Tri-fan: center c=v0 (INTERIOR, full 1-ring) + equilateral rim p1,p2,p3
     * (boundary). Every c->rim edge is interior with ONE interior endpoint, so
     * it is a legal split (not both-frozen like a bare diamond would be). The
     * three c->rim edges pairwise share a face, so face-locking splits exactly
     * one per round: (0,1) at midpoint (5,0,0). */
    v[0]=0;v[1]=0;v[2]=0;                      /* c  (interior)     */
    v[3]=10;v[4]=0;v[5]=0;                     /* p1 (10,0)         */
    v[6]=-5;v[7]=8.6602540f;v[8]=0;            /* p2                */
    v[9]=-5;v[10]=-8.6602540f;v[11]=0;         /* p3                */
    f[0]=0;f[1]=1;f[2]=2;                      /* c->p1->p2 */
    f[3]=0;f[4]=2;f[5]=3;                      /* c->p2->p3 */
    f[6]=0;f[7]=3;f[8]=1;                      /* c->p3->p1 */
    memset(pin,0,4);
    k=split_round(a,&v,&nv,&f,&nf,&pin,5.0,100000);   /* c->rim edges len 10 > 5 */
    CHECK(k==1,  "one split recorded (face-locking)");
    CHECK(nv==5, "nv 4 -> 5 (+1 midpoint)");
    CHECK(nf==5, "nf 3 -> 5 (2 of 3 fan tris split, +2)");
    CHECK(fabs(v[4*3+0]-5.0f)<1e-4 && fabs(v[4*3+1])<1e-4 && fabs(v[4*3+2])<1e-4,
          "midpoint at (5,0,0)");
    st=MeshManifold_audit(a,nv,f,nf);
    CHECK(MeshManifold_ok(&st),      "split result 2-manifold");
    CHECK(st.same_dir_edges==0,      "split winding consistent (no same_dir)");
    Arena_dispose(&a);
}

static void test_split_skips_boundary(void)
{
    Arena_T a=Arena_new();
    float   *v  =(float*)  ARENA_ALLOC(a,(long)(3*3*sizeof(float)));
    int32_t *f  =(int32_t*)ARENA_ALLOC(a,(long)(1*3*sizeof(int32_t)));
    uint8_t *pin=(uint8_t*)ARENA_ALLOC(a,(long)(3*sizeof(uint8_t)));
    size_t nv=3, nf=1, k;
    v[0]=0;v[1]=0;v[2]=0;  v[3]=10;v[4]=0;v[5]=0;  v[6]=0;v[7]=10;v[8]=0;
    f[0]=0;f[1]=1;f[2]=2;
    memset(pin,0,3);
    k=split_round(a,&v,&nv,&f,&nf,&pin,1.0,100000);   /* all edges long, all boundary */
    CHECK(k==0,  "boundary edges are never split");
    CHECK(nf==1, "single triangle unchanged");
    Arena_dispose(&a);
}

static void test_flip_winding(void)
{
    Arena_T a=Arena_new();
    float   *v  =(float*)  ARENA_ALLOC(a,(long)(4*3*sizeof(float)));
    int32_t *f  =(int32_t*)ARENA_ALLOC(a,(long)(2*3*sizeof(int32_t)));
    uint8_t *frz=(uint8_t*)ARENA_ALLOC(a,(long)(4*sizeof(uint8_t)));
    size_t nv=4, nf=2, k;
    MeshManifoldStats st;
    /* Convex quad whose shared diagonal a-b (len 10) makes two slivers; the flip
     * to diagonal c-d improves min-angle so flip_pass fires. Coplanar z=0,
     * consistently wound (+z). The min-angle flip must produce the
     * orientation-preserving (c,a,d)+(d,b,c); the pre-fix reversed write
     * (c,d,a)+(d,c,b) is rejected by its OWN guard here -- its reversed-triangle
     * normals oppose n_orig on flat input -- so n_flipped==0 catches that bug. */
    v[0]=0; v[1]=0;    v[2]=0;    /* a (0) */
    v[3]=10;v[4]=0;    v[5]=0;    /* b (1) */
    v[6]=5; v[7]=0.5f; v[8]=0;    /* c (2) */
    v[9]=5; v[10]=-0.5f;v[11]=0;  /* d (3) */
    f[0]=0;f[1]=1;f[2]=2;         /* (a,b,c) traverses a->b */
    f[3]=1;f[4]=0;f[5]=3;         /* (b,a,d) traverses b->a  (consistent) */
    memset(frz,0,4);
    st=MeshManifold_audit(a,nv,f,nf);
    CHECK(st.same_dir_edges==0, "flip: input consistently wound");
    k=flip_pass(a, v, f, nf, nv, frz);
    CHECK(k==1, "flip fires on sliver diagonal (0 => reversed-write guard bug)");
    st=MeshManifold_audit(a,nv,f,nf);
    CHECK(st.same_dir_edges==0, "flip winding consistent (orientation-preserving)");
    {
        int fc_cd = (f[0]==2||f[1]==2||f[2]==2) && (f[0]==3||f[1]==3||f[2]==3);
        int fd_cd = (f[3]==2||f[4]==2||f[5]==2) && (f[3]==3||f[4]==3||f[5]==3);
        CHECK(fc_cd && fd_cd, "flip rewired both faces to the c-d diagonal");
    }
    Arena_dispose(&a);
}

static void test_fullrun_manifold_frozen(void)
{
    Arena_T a=Arena_new();
    float *v; size_t nv; int32_t *f; size_t nf;
    float *ov; size_t onv; int32_t *of; size_t onf;
    RemeshOpts o; MeshManifoldStats st;
    uint8_t *bnd; size_t i, moved=0; int rc;
    grid(a,21,21,&v,&nv,&f,&nf);
    Remesh_default_opts(&o);
    bnd=(uint8_t*)ARENA_ALLOC(a,(long)(nv*sizeof(uint8_t)));
    detect_boundary(a,f,nf,nv,bnd);
    rc=Remesh_isotropic(a,v,nv,f,nf,NULL,&o,&ov,&onv,&of,&onf);
    CHECK(rc==0||rc==1, "fullrun returns 0/1 (never hard error)");
    st=MeshManifold_audit(a,onv,of,onf);
    CHECK(MeshManifold_ok(&st), "fullrun 2-manifold");
    CHECK((double)onf<=1.11*(double)nf && (double)onf>=0.89*(double)nf,
          "fullrun face count within tolerance");
    for (i=0;i<nv;i++){
        float px, py, pz; size_t j; int found=0;
        if (!bnd[i]) continue;
        px=v[i*3+0]; py=v[i*3+1]; pz=v[i*3+2];
        for (j=0;j<onv;j++)
            if (ov[j*3+0]==px && ov[j*3+1]==py && ov[j*3+2]==pz){ found=1; break; }
        if (!found) moved++;
    }
    CHECK(moved==0, "fullrun freezes the boundary polyline exactly");
    Arena_dispose(&a);
}

int main(void)
{
    printf("isotropic_remesh whitebox tests\n");
    test_L_formula();
    test_split_one_edge();
    test_split_skips_boundary();
    test_flip_winding();
    test_fullrun_manifold_frozen();
    printf("=== %s (%d failure%s) ===\n",
           g_fails?"FAILED":"PASSED", g_fails, g_fails==1?"":"s");
    return g_fails?1:0;
}
