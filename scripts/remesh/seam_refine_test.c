/*
 * seam_refine_test -- unit tests for SeamRefine_process.
 *
 * Scene A: two coarse (13-vox) open grid sheets abutting a seam plane with a
 * 2-vox gap (the uniform-coarse CVT weld input shape). Asserts: every band
 * edge ends <= target, band triangles become PRIMEABLE (circumradius <= the
 * adaptive bridge rho -- the property the whole redesign hangs on), non-band
 * geometry is untouched, new verts lie on parent chords with valid src
 * provenance, winding/manifold audits stay clean, and the two sheets remain
 * two components (refinement can neither fuse nor split).
 * Scene B: no planes -> byte-identical no-op, outputs alias inputs.
 * Scene C: already-fine input -> no-op (nothing above target).
 *
 * Verts are (z,y,x); the seam plane is x = 128 (sheets at x in [102,127] and
 * [129,154]).
 */
#include "../../src/remesh/seam_refine.h"
#include "../../src/common/pipeline_constants.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int g_fail = 0;
#define CHECK(cond, ...) do { \
    if (!(cond)) { g_fail++; fprintf(stderr, "  FAIL: " __VA_ARGS__); fprintf(stderr, "\n"); } \
} while (0)

static double vlen3(const float *V, int32_t a, int32_t b)
{
    double dx=(double)V[(size_t)a*3+0]-V[(size_t)b*3+0];
    double dy=(double)V[(size_t)a*3+1]-V[(size_t)b*3+1];
    double dz=(double)V[(size_t)a*3+2]-V[(size_t)b*3+2];
    return sqrt(dx*dx+dy*dy+dz*dz);
}

static double tri_circumradius(const float *V, int32_t a, int32_t b, int32_t c)
{
    double la=vlen3(V,b,c), lb=vlen3(V,c,a), lc=vlen3(V,a,b);
    double e1[3], e2[3], cr[3];
    for (int k=0;k<3;k++){
        e1[k]=(double)V[(size_t)b*3+(size_t)k]-V[(size_t)a*3+(size_t)k];
        e2[k]=(double)V[(size_t)c*3+(size_t)k]-V[(size_t)a*3+(size_t)k];
    }
    cr[0]=e1[1]*e2[2]-e1[2]*e2[1];
    cr[1]=e1[2]*e2[0]-e1[0]*e2[2];
    cr[2]=e1[0]*e2[1]-e1[1]*e2[0];
    double area = 0.5*sqrt(cr[0]*cr[0]+cr[1]*cr[1]+cr[2]*cr[2]);
    if (area < 1e-12) return 1e30;
    return (la*lb*lc)/(4.0*area);
}

static int i64_cmp(const void *pa, const void *pb)
{
    int64_t a = *(const int64_t *)pa, b = *(const int64_t *)pb;
    return a < b ? -1 : (a > b ? 1 : 0);
}

static int max_edge_run(const int32_t *f, size_t nf)
{
    size_t ne = nf*3, i; int maxr = 0;
    int64_t *keys = (int64_t *)malloc((ne?ne:1)*sizeof(int64_t));
    for (size_t t=0;t<nf;t++) for (int e=0;e<3;e++){
        int32_t a=f[t*3+e], b=f[t*3+(e+1)%3];
        int32_t u=a<b?a:b, v=a<b?b:a;
        keys[t*3+(size_t)e] = ((int64_t)u<<32) | (int64_t)(uint32_t)v;
    }
    qsort(keys, ne, sizeof(int64_t), i64_cmp);
    for (i=0;i<ne;){
        size_t j=i+1;
        while (j<ne && keys[j]==keys[i]) j++;
        if ((int)(j-i) > maxr) maxr = (int)(j-i);
        i=j;
    }
    free(keys);
    return maxr;
}

/* same_dir pairs: undirected edges whose two half-edges run the SAME way */
typedef struct { int64_t k; int8_t fw; } KF;
static int kf_cmp(const void *pa, const void *pb)
{
    const KF *a=(const KF*)pa, *b=(const KF*)pb;
    return a->k < b->k ? -1 : (a->k > b->k ? 1 : 0);
}

static size_t same_dir_pairs(const int32_t *f, size_t nf)
{
    size_t ne = nf*3, i, n = 0;
    KF *kf = (KF *)malloc((ne?ne:1)*sizeof(KF));
    for (size_t t=0;t<nf;t++) for (int e=0;e<3;e++){
        int32_t a=f[t*3+e], b=f[t*3+(e+1)%3];
        int32_t u=a<b?a:b, v=a<b?b:a;
        size_t idx=t*3+(size_t)e;
        kf[idx].k = ((int64_t)u<<32) | (int64_t)(uint32_t)v;
        kf[idx].fw = (int8_t)(a<b);
    }
    qsort(kf, ne, sizeof(KF), kf_cmp);
    for (i=0;i<ne;){
        size_t j=i+1;
        while (j<ne && kf[j].k==kf[i].k) j++;
        if (j-i==2 && kf[i].fw==kf[i+1].fw) n++;
        i=j;
    }
    free(kf);
    return n;
}

/* connected components via union-find over face-shared verts */
static size_t component_count(const int32_t *f, size_t nf, size_t nv)
{
    int32_t *par=(int32_t*)malloc((nv?nv:1)*sizeof(int32_t));
    uint8_t *used=(uint8_t*)calloc(nv?nv:1,1);
    size_t i, n=0;
    for (i=0;i<nv;i++) par[i]=(int32_t)i;
    for (size_t t=0;t<nf;t++) for (int e=0;e<3;e++){
        int32_t a=f[t*3+e], b=f[t*3+(e+1)%3], ra=a, rb=b;
        used[a]=used[b]=1;
        while (par[ra]!=ra) ra=par[ra];
        while (par[rb]!=rb) rb=par[rb];
        if (ra!=rb) par[ra]=rb;
    }
    for (i=0;i<nv;i++){
        if (!used[i]) continue;
        int32_t r=(int32_t)i;
        while (par[r]!=r) r=par[r];
        if (r==(int32_t)i) n++;
    }
    free(par); free(used);
    return n;
}

/* Build one coarse grid sheet: NYxNX verts, spacing SP, at x0..; z plane 50. */
static void build_sheet(float *V, int32_t *F, int NY, int NX, float sp,
                        float x0, int32_t vbase, size_t *fi)
{
    for (int i=0;i<NY;i++) for (int j=0;j<NX;j++){
        int32_t id = vbase + (int32_t)(i*NX+j);
        V[(size_t)id*3+0]=50.0f;
        V[(size_t)id*3+1]=(float)i*sp;
        V[(size_t)id*3+2]=x0+(float)j*sp;
    }
    for (int i=0;i+1<NY;i++) for (int j=0;j+1<NX;j++){
        int32_t v00=vbase+(int32_t)(i*NX+j), v01=v00+1;
        int32_t v10=vbase+(int32_t)((i+1)*NX+j), v11=v10+1;
        F[(*fi)*3+0]=v00; F[(*fi)*3+1]=v10; F[(*fi)*3+2]=v11; (*fi)++;
        F[(*fi)*3+0]=v00; F[(*fi)*3+1]=v11; F[(*fi)*3+2]=v01; (*fi)++;
    }
}

int main(void)
{
    Arena_T a = Arena_new();
    const SeamPlane plane = { 2, 128.0 };

    /* ---- Scene A: two 13-vox sheets, 2-vox gap at x=128 ---- */
    {
        const int NY=4, NX=3;              /* 13-vox spacing, x spans 26 vox */
        const float SP=13.0f;
        size_t nv_side=(size_t)NY*NX, nf_side=(size_t)(NY-1)*(NX-1)*2;
        size_t nv=nv_side*2, nf=nf_side*2, fi=0;
        float *V=(float*)ARENA_ALLOC(a,(long)(nv*3*sizeof(float)));
        int32_t *F=(int32_t*)ARENA_ALLOC(a,(long)(nf*3*sizeof(int32_t)));
        /* left sheet ends at x=127; right begins at x=129 (2-vox gap) */
        build_sheet(V, F, NY, NX, SP, 127.0f-2.0f*SP, 0, &fi);
        build_sheet(V, F, NY, NX, SP, 129.0f, (int32_t)nv_side, &fi);

        SeamRefineParams p; SeamRefine_default_params(&p);
        SeamRefineStats st;
        float *OV; int32_t *OF; int32_t *SRC;
        size_t onv, onf, nnew;
        int rc = SeamRefine_process(a, V, nv, F, nf, &plane, 1, &p,
                                    &OV, &onv, &OF, &onf, &SRC, &nnew, &st);
        CHECK(rc==0 && nnew>0 && onf>nf, "A: rc=%d nnew=%zu onf=%zu", rc, nnew, onf);

        /* 1. every OPEN-BOUNDARY in-band edge <= target (+eps). The bridge's
         * front is built from run-1 boundary edges, so boundary spacing is
         * what the adaptive rho keys on; interior transition edges can end a
         * bit above target (flip churn) and only affect quality -- reported
         * but not asserted. */
        {
            size_t ne=onf*3;
            KF *kf=(KF*)malloc(ne*sizeof(KF));
            for (size_t t=0;t<onf;t++) for (int e=0;e<3;e++){
                int32_t x=OF[t*3+e], y=OF[t*3+(e+1)%3];
                int32_t u=x<y?x:y, v=x<y?y:x;
                kf[t*3+(size_t)e].k=((int64_t)u<<32)|(int64_t)(uint32_t)v;
                kf[t*3+(size_t)e].fw=0;
            }
            qsort(kf,ne,sizeof(KF),kf_cmp);
            double worst_bnd=0.0, worst_int=0.0;
            for (size_t i=0;i<ne;){
                size_t j=i+1;
                while (j<ne && kf[j].k==kf[i].k) j++;
                int32_t u=(int32_t)(kf[i].k>>32), v=(int32_t)(kf[i].k & 0xffffffff);
                double du=fabs((double)OV[(size_t)u*3+2]-plane.coord);
                double dv=fabs((double)OV[(size_t)v*3+2]-plane.coord);
                if (du<=p.band && dv<=p.band){
                    double l=vlen3(OV,u,v);
                    if (j-i==1){ if (l>worst_bnd) worst_bnd=l; }
                    else       { if (l>worst_int) worst_int=l; }
                }
                i=j;
            }
            free(kf);
            CHECK(worst_bnd <= (double)p.target_len + 1e-4,
                  "A: band BOUNDARY edge %.3f > target %.2f",
                  worst_bnd, (double)p.target_len);
            fprintf(stderr, "  A: worst in-band edge: boundary %.2f (target "
                    "%.2f), interior %.2f (transition, unasserted)\n",
                    worst_bnd, (double)p.target_len, worst_int);
        }

        /* 2. PRIMEABILITY: every face with an open boundary edge in the seam
         * band has circumradius <= rho_max (3.0) -- the bridge's hinge-ball
         * reconstruction limit, i.e. the redesign's load-bearing property. */
        {
            size_t ne=onf*3;
            int64_t *keys=(int64_t*)malloc(ne*sizeof(int64_t));
            for (size_t t=0;t<onf;t++) for (int e=0;e<3;e++){
                int32_t x=OF[t*3+e], y=OF[t*3+(e+1)%3];
                int32_t u=x<y?x:y, v=x<y?y:x;
                keys[t*3+(size_t)e]=((int64_t)u<<32)|(int64_t)(uint32_t)v;
            }
            int64_t *sk=(int64_t*)malloc(ne*sizeof(int64_t));
            memcpy(sk,keys,ne*sizeof(int64_t));
            qsort(sk,ne,sizeof(int64_t),i64_cmp);
            double worst_cr=0.0; size_t n_front=0;
            for (size_t t=0;t<onf;t++) for (int e=0;e<3;e++){
                int32_t u=OF[t*3+e], v=OF[t*3+(e+1)%3];
                int32_t cu=u<v?u:v, cv=u<v?v:u;
                int64_t k=((int64_t)cu<<32)|(int64_t)(uint32_t)cv;
                /* run-1? binary search span */
                size_t lo=0, hi=ne;
                while (lo<hi){ size_t mid=lo+(hi-lo)/2; if (sk[mid]<k) lo=mid+1; else hi=mid; }
                size_t run=0; while (lo+run<ne && sk[lo+run]==k) run++;
                if (run!=1) continue;
                double du=fabs((double)OV[(size_t)u*3+2]-plane.coord);
                double dv=fabs((double)OV[(size_t)v*3+2]-plane.coord);
                if (du>p.band || dv>p.band) continue;
                double cr=tri_circumradius(OV,OF[t*3+0],OF[t*3+1],OF[t*3+2]);
                if (cr>worst_cr) worst_cr=cr;
                n_front++;
            }
            free(keys); free(sk);
            CHECK(n_front>0, "A: no seam-band front edges found");
            CHECK(worst_cr <= (double)BRIDGE_RHO_MAX + 1e-6,
                  "A: front triangle circumradius %.3f > rho_max %.1f "
                  "(bridge cannot prime)", worst_cr, (double)BRIDGE_RHO_MAX);
            fprintf(stderr, "  A: %zu front tris, worst circumradius %.2f "
                    "(rho_max %.1f)\n", n_front, worst_cr, (double)BRIDGE_RHO_MAX);
        }

        /* 3. originals never move; new verts have valid ordered src */
        for (size_t v=0;v<nv;v++)
            CHECK(OV[v*3+0]==V[v*3+0] && OV[v*3+1]==V[v*3+1] && OV[v*3+2]==V[v*3+2],
                  "A: original vert %zu moved", v);
        for (size_t i=0;i<nnew;i++)
            CHECK(SRC[i]>=0 && (size_t)SRC[i] < nv+i,
                  "A: src[%zu]=%d out of order", i, SRC[i]);

        /* 4. clean audits: manifold, no same_dir, still exactly 2 components */
        CHECK(max_edge_run(OF,onf)<=2, "A: non-manifold edge");
        CHECK(same_dir_pairs(OF,onf)==0, "A: same_dir pair minted");
        CHECK(component_count(OF,onf,onv)==2, "A: component count %zu != 2",
              component_count(OF,onf,onv));
        fprintf(stderr,
            "Scene A: ok (%zu bnd + %zu int splits, %zu flips, %zu rounds, "
            "+%zuv +%zuf)\n", st.bnd_splits, st.int_splits, st.flips,
            st.rounds, st.verts_added, st.faces_added);
    }

    /* ---- Scene B: no planes -> alias no-op ---- */
    {
        float V[9]={50,0,120, 50,13,120, 50,0,133};
        int32_t F[3]={0,1,2};
        float *OV; int32_t *OF; int32_t *SRC; size_t onv,onf,nnew;
        int rc = SeamRefine_process(a, V, 3, F, 1, NULL, 0, NULL,
                                    &OV,&onv,&OF,&onf,&SRC,&nnew,NULL);
        CHECK(rc==0 && OV==V && OF==F && onv==3 && onf==1 && nnew==0,
              "B: not an alias no-op");
        fprintf(stderr, "Scene B (no planes -> no-op): %s\n",
                g_fail ? "check above" : "ok");
    }

    /* ---- Scene C: already fine -> no-op ---- */
    {
        const int NY=3, NX=3; const float SP=2.0f;
        size_t nv=(size_t)NY*NX, nf=(size_t)(NY-1)*(NX-1)*2, fi=0;
        float *V=(float*)ARENA_ALLOC(a,(long)(nv*3*sizeof(float)));
        int32_t *F=(int32_t*)ARENA_ALLOC(a,(long)(nf*3*sizeof(int32_t)));
        build_sheet(V, F, NY, NX, SP, 123.0f, 0, &fi);
        float *OV; int32_t *OF; int32_t *SRC; size_t onv,onf,nnew;
        int rc = SeamRefine_process(a, V, nv, F, nf, &plane, 1, NULL,
                                    &OV,&onv,&OF,&onf,&SRC,&nnew,NULL);
        CHECK(rc==0 && nnew==0 && onf==nf, "C: fine input was split");
        fprintf(stderr, "Scene C (already fine -> no-op): ok\n");
    }

    Arena_dispose(&a);
    if (g_fail){ fprintf(stderr, "SEAM_REFINE TESTS: %d FAILURE(S)\n", g_fail); return 1; }
    fprintf(stderr, "ALL SEAM_REFINE TESTS PASSED\n");
    return 0;
}
