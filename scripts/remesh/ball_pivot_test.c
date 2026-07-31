/* ball_pivot_test.c — unit tests for Ball-Pivoting reconstruction
 * (src/remesh/ball_pivot.c). Whitebox: it #includes the .c so it can reach the
 * static pivot_push / glue_folds_back helpers in addition to the public
 * BallPivot_reconstruct.
 *
 * Covers the 2026-06-02 BPA fixes:
 *   - pivot_push keeps the K smallest-theta candidates, ascending, with an
 *     index tie-break (the Part-1 candidate-retry mechanism);
 *   - glue_folds_back flags a same-side (folded) apex pair and passes a flat one;
 *   - a flat sheet reconstructs into a 2-manifold disk with NO interior slits
 *     (the over-rejection regression the retry fixes);
 *   - two anti-parallel sheets within 2*rho are NOT bridged and stay 2-manifold
 *     (the under-rejection regression the anti-parallel filter fixes).
 *
 * Exit 0 = all pass, 1 = a failure.
 */
#include "remesh/ball_pivot.c"   /* whitebox include (pulls in ball_pivot.h, arena.h) */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static int g_fails = 0;
#define CHECK(cond, msg) do { \
    if (!(cond)) { printf("  FAIL: %s\n", (msg)); g_fails++; } \
    else         { printf("  ok:   %s\n", (msg)); } \
} while (0)

/* ---------- mesh analysis on (nv, faces) -------------------------------- */
typedef struct {
    int max_edge_mult;     /* >2 => non-manifold edge                       */
    int boundary_edges;    /* edges used by exactly 1 face                   */
    int boundary_loops;    /* simple boundary loops (valid iff bad_deg==0)   */
    int boundary_bad_deg;  /* boundary verts with degree != 2               */
    int components;        /* connected components of the face graph         */
    int cross_set_edges;   /* edges straddling the index `split` boundary    */
} MeshStats;

static int uf_find(int *p, int x){ while(p[x]!=x){p[x]=p[p[x]];x=p[x];} return x; }
static void uf_union(int *p,int a,int b){ a=uf_find(p,a); b=uf_find(p,b); if(a!=b)p[a]=b; }

static MeshStats analyze(int nv, const int32_t *F, size_t nf, int split)
{
    MeshStats s; memset(&s, 0, sizeof s);
    int *mult = (int *)calloc((size_t)nv*(size_t)nv, sizeof(int));
    int *uf   = (int *)malloc((size_t)nv*sizeof(int));
    int *bdeg = (int *)calloc((size_t)nv, sizeof(int));
    int *bn1  = (int *)malloc((size_t)nv*sizeof(int));
    int *bn2  = (int *)malloc((size_t)nv*sizeof(int));
    for (int i=0;i<nv;i++){ uf[i]=i; bn1[i]=-1; bn2[i]=-1; }

    for (size_t f=0; f<nf; f++){
        int v[3] = { F[f*3+0], F[f*3+1], F[f*3+2] };
        int e[3][2] = {{v[0],v[1]},{v[1],v[2]},{v[0],v[2]}};
        for (int k=0;k<3;k++){
            int u=e[k][0], w=e[k][1]; if(u>w){int t=u;u=w;w=t;}
            mult[(size_t)u*nv+w]++;
            uf_union(uf,u,w);
        }
    }
    for (int u=0;u<nv;u++) for (int w=u+1;w<nv;w++){
        int m = mult[(size_t)u*nv+w];
        if (m==0) continue;
        if (m > s.max_edge_mult) s.max_edge_mult = m;
        if ((u<split) != (w<split)) s.cross_set_edges++;
        if (m==1){
            s.boundary_edges++;
            if (bdeg[u]==0) bn1[u]=w; else if (bdeg[u]==1) bn2[u]=w; bdeg[u]++;
            if (bdeg[w]==0) bn1[w]=u; else if (bdeg[w]==1) bn2[w]=u; bdeg[w]++;
        }
    }
    {   /* components: distinct roots referenced by faces */
        int *seen=(int*)calloc((size_t)nv,sizeof(int));
        for (size_t f=0; f<nf; f++) for(int k=0;k<3;k++){
            int r=uf_find(uf,F[f*3+k]); if(!seen[r]){seen[r]=1;s.components++;}
        }
        free(seen);
    }
    for (int i=0;i<nv;i++) if (bdeg[i]!=0 && bdeg[i]!=2) s.boundary_bad_deg++;
    if (s.boundary_bad_deg==0){
        int *vis=(int*)calloc((size_t)nv,sizeof(int));
        for (int i=0;i<nv;i++) if (bdeg[i]==2 && !vis[i]){
            int cur=i, prev=-1, guard=0;
            while (cur!=-1 && !vis[cur] && guard++ < nv+2){
                vis[cur]=1;
                int nx = (bn1[cur]!=prev)?bn1[cur]:bn2[cur];
                prev=cur; cur=nx;
            }
            s.boundary_loops++;
        }
        free(vis);
    }
    free(mult); free(uf); free(bdeg); free(bn1); free(bn2);
    return s;
}

/* Hex lattice (nearest-neighbour spacing 1.0) flat patch at z=zc, normal (nz,0,0).
 * Coords are (z,y,x). Writes R*C verts starting at vertex `off`; returns count. */
static int build_hex(float *V, float *N, int off, int R, int C, float zc, float nz)
{
    int n=0;
    for (int r=0;r<R;r++) for (int c=0;c<C;c++){
        int i = off + n;
        V[i*3+0]=zc;
        V[i*3+1]=(float)(r*0.86602540378);
        V[i*3+2]=(float)c + ((r&1)?0.5f:0.0f);
        N[i*3+0]=nz; N[i*3+1]=0.0f; N[i*3+2]=0.0f;
        n++;
    }
    return n;
}

/* ---------- test 1: pivot_push top-K ------------------------------------ */
static void test_pivot_push(void)
{
    printf("[pivot_push]\n");
    double O0[3] = {0,0,0};
    PivotCtx p; p.n_k = 0;
    double th[] = {2.0, 0.5, 3.0, 1.0, 0.1, 2.5, 0.7};
    int    id[] = {10,  11,  12,  13,  14,  15,  16 };
    for (int i=0;i<7;i++) pivot_push(&p, id[i], th[i], O0);
    CHECK(p.n_k == BPA_PIVOT_K, "caps at BPA_PIVOT_K");
    int asc=1; for (int i=1;i<p.n_k;i++) if (p.k_theta[i] < p.k_theta[i-1]) asc=0;
    CHECK(asc, "kept in ascending theta order");
    CHECK(fabs(p.k_theta[0]-0.1)<1e-12 && fabs(p.k_theta[BPA_PIVOT_K-1]-1.0)<1e-12,
          "kept the K smallest thetas (0.1..1.0)");

    PivotCtx q; q.n_k=0;
    pivot_push(&q, 20, 1.0, O0);
    pivot_push(&q,  5, 1.0, O0);
    CHECK(q.k_v[0]==5 && q.k_v[1]==20, "equal-theta tie-break by smaller index");
}

/* ---------- test 2: glue_folds_back ------------------------------------- */
static void test_glue_folds_back(void)
{
    printf("[glue_folds_back]\n");
    /* shared edge a-b along +x in the z=0 plane (coords z,y,x) */
    Vec3 V[4];
    V[0].z=0; V[0].y=0; V[0].x=0;          /* a            */
    V[1].z=0; V[1].y=0; V[1].x=1;          /* b            */
    V[2].z=0; V[2].y=1;    V[2].x=0.5f;    /* new apex +y  */
    V[3].z=0; V[3].y=-1;   V[3].x=0.5f;    /* old apex -y (opposite => flat) */
    CHECK(glue_folds_back(V,0,1,2,3)==0, "opposite-side apexes => not a fold");
    V[3].z=0; V[3].y=0.9f; V[3].x=0.5f;    /* old apex now +y (same side => fold) */
    CHECK(glue_folds_back(V,0,1,2,3)==1, "same-side apexes => fold detected");
}

/* ---------- test 2b: anti-parallel filter decision + threshold --------- */
static void test_antiparallel_decision(void)
{
    printf("[cand_antiparallel]\n");
    PivotCtx p;
    p.n_front[0]=1;   p.n_front[1]=0;   p.n_front[2]=0;     /* front faces +z      */
    p.n_front_v[0]=1; p.n_front_v[1]=0; p.n_front_v[2]=0;
    g_antiparallel_cos = BPA_ANTIPARALLEL_COS;              /* -0.3 (~107deg)      */

    double up[3]    = { 1.0, 0.0, 0.0};   /* parallel        => accept */
    double down[3]  = {-1.0, 0.0, 0.0};   /* anti-parallel   => reject */
    double side[3]  = { 0.0, 1.0, 0.0};   /* 90deg (dot 0)   => accept */
    double tilt[3]  = { 0.6, 0.8, 0.0};   /* ~53deg          => accept */
    double steep[3] = {-0.5, 0.866, 0.0}; /* 120deg (dot-.5) => reject */
    CHECK(cand_antiparallel(up,   &p)==0, "parallel normal accepted");
    CHECK(cand_antiparallel(down, &p)==1, "anti-parallel normal rejected");
    CHECK(cand_antiparallel(side, &p)==0, "perpendicular (dot 0 > -0.3) accepted");
    CHECK(cand_antiparallel(tilt, &p)==0, "moderate curvature accepted");
    CHECK(cand_antiparallel(steep,&p)==1, "beyond ~107deg rejected");

    /* The geometric/MLS disagreement guard: if EITHER reference normal is
     * anti-parallel to the candidate, reject (use n_front_v as the dissenter). */
    p.n_front_v[0]=-1; p.n_front_v[1]=0; p.n_front_v[2]=0;
    CHECK(cand_antiparallel(up, &p)==1, "either-normal-disagrees rejects");
}

/* ---------- test 3: flat sheet -> manifold disk, no interior slits ------ */
static void test_flat_sheet(void)
{
    printf("[flat sheet]\n");
    int R=10,C=10, nv=R*C;
    float *V=(float*)malloc((size_t)nv*3*sizeof(float));
    float *N=(float*)malloc((size_t)nv*3*sizeof(float));
    build_hex(V,N,0,R,C,0.0f,1.0f);
    Arena_T a=Arena_new();
    int32_t *F=NULL; size_t nf=0;
    int rc=BallPivot_reconstruct(a,V,N,(size_t)nv,1.2f,&F,&nf);
    CHECK(rc==0 && nf>0, "BPA produced faces");
    if (rc==0 && nf>0){
        MeshStats s=analyze(nv,F,nf,nv);
        printf("    nf=%zu maxmult=%d bnd_edges=%d loops=%d baddeg=%d comp=%d\n",
               nf,s.max_edge_mult,s.boundary_edges,s.boundary_loops,
               s.boundary_bad_deg,s.components);
        CHECK(s.max_edge_mult<=2, "edge-manifold (no >2-face edge)");
        CHECK(s.components==1, "single connected component");
        CHECK(s.boundary_bad_deg==0, "boundary verts all degree 2 (no bowtie boundary)");
        CHECK(s.boundary_loops==1, "exactly one boundary loop (no interior slits)");
    }
    Arena_dispose(&a); free(V); free(N);
}

/* ---------- test 4: anti-parallel sheets are not bridged ---------------- */
static void test_antiparallel(void)
{
    printf("[anti-parallel sheets]\n");
    /* Two back-to-back sheets 0.8 vox apart with ANTI-parallel normals (A faces
     * -z, B faces +z, so they face away from each other). Each seeds cleanly on
     * its own outward side; the question is whether BPA stitches a bridge across
     * the 0.8-vox gap. It must not -- the per-vertex normal gate AND the new
     * front-coherence filter both forbid an anti-parallel join. */
    int R=8,C=8, half=R*C, nv=2*half;
    float *V=(float*)malloc((size_t)nv*3*sizeof(float));
    float *N=(float*)malloc((size_t)nv*3*sizeof(float));
    build_hex(V,N,0,    R,C, 0.0f, -1.0f);   /* sheet A: z=0,   normal -z */
    build_hex(V,N,half, R,C, 0.8f,  1.0f);   /* sheet B: z=0.8, normal +z */
    Arena_T a=Arena_new();
    int32_t *F=NULL; size_t nf=0;
    int rc=BallPivot_reconstruct(a,V,N,(size_t)nv,1.2f,&F,&nf);   /* filter ON by default */
    CHECK(rc==0 && nf>0, "BPA produced faces");
    if (rc==0 && nf>0){
        MeshStats s=analyze(nv,F,nf,half);
        printf("    nf=%zu maxmult=%d cross=%d comp=%d antip_rej=%ld\n",
               nf,s.max_edge_mult,s.cross_set_edges,s.components,g_dbg_antiparallel);
        CHECK(s.max_edge_mult<=2, "edge-manifold");
        CHECK(s.cross_set_edges==0, "NO bridge between the anti-parallel sheets");
        CHECK(s.components>=2, "sheets remain separate components");
    }
    Arena_dispose(&a); free(V); free(N);
}

/* ---------- test 2c: wall-guard decision (glue_is_wall) ----------------- */
static void test_wall_guard_decision(void)
{
    printf("[glue_is_wall]\n");
    /* sheet normal +z (coords z,y,x). Triangle is (t, v_new, h) = (0, 2, 1). */
    Vec3 V[3], N[3];
    N[0].z=1; N[0].y=0; N[0].x=0;
    N[1].z=1; N[1].y=0; N[1].x=0;
    N[2].z=1; N[2].y=0; N[2].x=0;
    /* (a) a flat surface triangle in z=0: face normal ~parallel to +z, short edges */
    V[0].z=0; V[0].y=0; V[0].x=0;
    V[1].z=0; V[1].y=0; V[1].x=1;
    V[2].z=0; V[2].y=1; V[2].x=0;
    CHECK(glue_is_wall(V,N,0,1,2, 0.34,2.0)==0, "flat surface triangle is not a wall");
    /* (b) a steep, long bridge: v_new lifted 2.5 vox along the normal */
    V[2].z=2.5f; V[2].y=0; V[2].x=0.5f;
    CHECK(glue_is_wall(V,N,0,1,2, 0.34,2.0)==1, "steep + long bridge IS a wall");
    /* (c) steep but SHORT (a fold/rim, v_new near): not a wall */
    V[2].z=0.5f; V[2].y=0; V[2].x=0.5f;
    CHECK(glue_is_wall(V,N,0,1,2, 0.34,2.0)==0, "steep but short edge is not a wall");
    /* (d) NULL normals -> never a wall (guard no-op) */
    V[2].z=2.5f;
    CHECK(glue_is_wall(V,NULL,0,1,2, 0.34,2.0)==0, "no normals -> not a wall");
}

/* ---------- test 5: parallel stacked sheets are not fused --------------- */
static void test_wall_guard(void)
{
    printf("[wall-guard: parallel stacked sheets]\n");
    /* Two PARALLEL sheets (both normal +z) 2.2 vox apart -- within BPA's 2*rho=2.4
     * reach, so BPA would otherwise weld them into one fused stack. The
     * anti-parallel filter does NOT catch this (the normals agree); the wall-guard
     * must, by rejecting the steep long "wall" bridge triangles. */
    int R=8,C=8, half=R*C, nv=2*half;
    float *V=(float*)malloc((size_t)nv*3*sizeof(float));
    float *N=(float*)malloc((size_t)nv*3*sizeof(float));
    build_hex(V,N,0,    R,C, 0.0f, 1.0f);   /* sheet A: z=0,   normal +z */
    build_hex(V,N,half, R,C, 2.2f, 1.0f);   /* sheet B: z=2.2, normal +z (PARALLEL) */
    Arena_T a=Arena_new();
    int32_t *F=NULL; size_t nf=0;
    int rc=BallPivot_reconstruct(a,V,N,(size_t)nv,1.2f,&F,&nf);   /* guard ON by default */
    CHECK(rc==0 && nf>0, "BPA produced faces");
    if (rc==0 && nf>0){
        MeshStats s=analyze(nv,F,nf,half);
        printf("    nf=%zu maxmult=%d cross=%d comp=%d wall_rej=%ld\n",
               nf,s.max_edge_mult,s.cross_set_edges,s.components,g_dbg_wall);
        CHECK(s.cross_set_edges==0, "NO bridge between the parallel stacked sheets");
        CHECK(s.components>=2, "sheets remain separate components");
        /* wall_rej may be 0 in this CLEAN flat synthetic: at 2.2 vox the earlier
         * per-vertex normal / empty-ball tests already refuse the bridge. The
         * wall-guard is the BACKSTOP for CURVED/compacted stacks those miss (the
         * real monster cube -- see the A/B); glue_is_wall above unit-tests its
         * decision directly, and the flat-sheet test confirms it harms nothing. */
        printf("    (wall_rej=%ld -- backstop; earlier guards may suffice here)\n", g_dbg_wall);
    }
    Arena_dispose(&a); free(V); free(N);
}

/* ---------- test 6: accumulated phase catches a gradual staircase -------- */
static void test_grow_phase_anchor(void)
{
    printf("[growth winding anchor]\n");
    BpaBuild b;
    memset(&b, 0, sizeof b);
    double phase[6];
    for (int i=0;i<6;i++) phase[i]=NAN;
    b.grow_phase=phase;
    b.grow_pitch=9.5;
    b.grow_tol=0.45;

    Vec3 V[6]; memset(V,0,sizeof V);
    /* One true Archimedean wrap crossing atan2's +pi/-pi branch. */
    double th0=M_PI-0.02, th1=M_PI+0.02;
    double r0=100.0+9.5*th0/(2.0*M_PI);
    double r1=100.0+9.5*th1/(2.0*M_PI);
    V[0].y=(float)(r0*sin(th0)); V[0].x=(float)(r0*cos(th0));
    V[1].y=(float)(r1*sin(th1)); V[1].x=(float)(r1*cos(th1));
    double ds=bpa_grow_phase_step(&b,V,0,1);
    CHECK(fabs(ds)<1e-5, "same spiral stays phase-constant across atan2 branch");

    /* A wall climbed in 0.2-turn radial increments. Every local step is below
     * tolerance, but the seed-relative phase reaches 0.6 and must be rejected. */
    for (int i=2;i<6;i++) {
        double r=100.0+9.5*0.2*(double)(i-2);
        V[i].y=0.0f; V[i].x=(float)r;
    }
    phase[2]=0.0; phase[3]=0.2;
    double pc=0.0;
    CHECK(bpa_grow_phase_candidate(&b,V,2,3,4,&pc)==1 && fabs(pc-0.4)<1e-5,
          "gradual staircase remains allowed inside seed tolerance");
    phase[4]=pc;
    CHECK(bpa_grow_phase_candidate(&b,V,3,4,5,&pc)==0,
          "gradual staircase is rejected after cumulative phase exceeds tolerance");

    /* A seed triangle sampled from one ideal spiral is centered near zero. */
    double th2=M_PI+0.04;
    double r2=100.0+9.5*th2/(2.0*M_PI);
    V[2].y=(float)(r2*sin(th2)); V[2].x=(float)(r2*cos(th2));
    phase[0]=phase[1]=phase[2]=NAN;
    bpa_grow_phase_seed(&b,V,0,1,2);
    CHECK(fabs(phase[0])<1e-5 && fabs(phase[1])<1e-5 && fabs(phase[2])<1e-5,
          "same-wrap seed receives one coherent phase anchor");
}

int main(void)
{
    setvbuf(stdout, NULL, _IONBF, 0);   /* unbuffered: survive an abort() mid-test */
    printf("=== ball_pivot_test ===\n");
    test_pivot_push();
    test_glue_folds_back();
    test_antiparallel_decision();
    test_wall_guard_decision();
    test_flat_sheet();
    test_antiparallel();
    test_wall_guard();
    test_grow_phase_anchor();
    printf("=== %s (%d failure%s) ===\n",
           g_fails==0 ? "PASS" : "FAIL", g_fails, g_fails==1?"":"s");
    return g_fails==0 ? 0 : 1;
}
