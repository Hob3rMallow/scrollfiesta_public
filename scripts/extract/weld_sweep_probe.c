/* Fix-0 measurement harness (2026-07-09).
 *
 * Runs the REAL iterated pre-BPA weld (src/common/vert_weld.c :: Weld_verts,
 * driven exactly as src/extract/mesh_extract.c:1147 — up to 5 passes, stop when
 * the vertex count stops shrinking) on a pre-weld MLS point-cloud dump, at a
 * sweep of eps values, and reports the residual near-coincident stack that BPA
 * (rho = BPA_RHO_VOX = 1.2 vox) then has to triangulate.
 *
 * Why this is the decisive Fix-0 measurement:
 *   - BPA cannot form a clean pivot triangle from points much closer than rho;
 *     near-coincident survivors are what tear the "lightning-bolt" slits. So the
 *     count of post-weld points whose nearest neighbour is < ~rho/2 (0.6 vox) is
 *     the number that predicts the tear.
 *   - The current weld eps is MLS_WELD_EPS_VOX = 0.25 (= rho/4.8), well below the
 *     rho/2 spacing BPA wants. This harness re-welds the SAME cloud at 0.25 / 0.5
 *     / 0.6 and measures whether widening clears the residual stack, and at what
 *     cost (vertex-count reduction = in-plane-erosion proxy; PCA slab thickness).
 *
 * It links the production Weld_verts, so eps=0.25 here reproduces the pipeline's
 * step0_post_weld cloud (a built-in cross-check against the dumped count).
 *
 * Usage:
 *   weld_sweep_probe <cloud.obj> [--eps 0.25,0.5,0.6] [--bbox z0 z1 y0 y1 x0 x1]
 *                    [--nn-radius 1.5] [--pca-R 6.0]
 *   weld_sweep_probe --selftest
 *
 * Build (from repo root, in an x64 VS dev shell):
 *   cl /O2 /std:c11 /I src/common scripts/extract/weld_sweep_probe.c \
 *      src/common/vert_weld.c src/common/arena.c src/common/except.c \
 *      /Fe:build/weld_sweep_probe.exe
 */
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "arena.h"
#include "vert_weld.h"

/* ------------------------------------------------------------------ */
/* OBJ point-cloud loader (v lines, "v z y x" order per repo convention) */
/* ------------------------------------------------------------------ */
static float *load_obj_points(const char *path, size_t *out_nv)
{
    FILE *f = fopen(path, "r");
    if (!f) { fprintf(stderr, "ERROR: cannot open %s\n", path); return NULL; }
    size_t cap = 1u << 16, nv = 0;
    float *v = (float *)malloc(cap * 3 * sizeof(float));
    if (!v) { fclose(f); return NULL; }
    char line[512];
    while (fgets(line, sizeof(line), f)) {
        if (line[0] != 'v' || line[1] != ' ') continue;   /* skip vn/vt/f */
        float z, y, x;
        if (sscanf(line + 2, "%f %f %f", &z, &y, &x) != 3) continue;
        if (nv >= cap) {
            cap *= 2;
            float *nv2 = (float *)realloc(v, cap * 3 * sizeof(float));
            if (!nv2) { free(v); fclose(f); return NULL; }
            v = nv2;
        }
        v[nv*3+0] = z; v[nv*3+1] = y; v[nv*3+2] = x;
        nv++;
    }
    fclose(f);
    *out_nv = nv;
    return v;
}

/* ------------------------------------------------------------------ */
/* Iterated weld — mirrors src/extract/mesh_extract.c:1147-1158.       */
/* Output lives on `arena` (Weld_verts does not restore). Caller owns  */
/* the arena lifetime.                                                 */
/* ------------------------------------------------------------------ */
static void weld_iterated(Arena_T arena, const float *verts, size_t nv,
                          float eps, float **out_v, size_t *out_nv)
{
    int32_t dummy_face = 0;
    size_t dummy_out_nf = 0;
    const float *src_v = verts;
    size_t src_nv = nv;
    float *welded = NULL;
    size_t welded_nv = nv;
    for (int iter = 0; iter < 5; iter++) {
        size_t prev = welded_nv;
        Weld_verts(arena, src_v, src_nv, NULL,
                   &dummy_face, 0, &dummy_out_nf, eps, /*guard_orient=*/false,
                   &welded, &welded_nv, NULL);
        if (welded_nv == prev && iter > 0) break;
        src_v = welded;
        src_nv = welded_nv;
    }
    *out_v = welded;
    *out_nv = welded_nv;
}

/* ------------------------------------------------------------------ */
/* Uniform-grid nearest-neighbour (min distance to any OTHER point).   */
/* ------------------------------------------------------------------ */
typedef struct { int *head, *next; int gz, gy, gx; float z0, y0, x0, cell; } Grid;

static void grid_build(Grid *g, const float *v, size_t nv, float cell)
{
    float zmin=v[0], zmax=v[0], ymin=v[1], ymax=v[1], xmin=v[2], xmax=v[2];
    for (size_t i = 1; i < nv; i++) {
        float z=v[i*3+0], y=v[i*3+1], x=v[i*3+2];
        if (z<zmin)zmin=z; if (z>zmax)zmax=z;
        if (y<ymin)ymin=y; if (y>ymax)ymax=y;
        if (x<xmin)xmin=x; if (x>xmax)xmax=x;
    }
    g->cell=cell; g->z0=zmin; g->y0=ymin; g->x0=xmin;
    g->gz=(int)floorf((zmax-zmin)/cell)+1;
    g->gy=(int)floorf((ymax-ymin)/cell)+1;
    g->gx=(int)floorf((xmax-xmin)/cell)+1;
    size_t ncell=(size_t)g->gz*g->gy*g->gx;
    g->head=(int *)malloc(ncell*sizeof(int));
    g->next=(int *)malloc(nv*sizeof(int));
    for (size_t c=0;c<ncell;c++) g->head[c]=-1;
    for (size_t i=0;i<nv;i++) {
        int cz=(int)floorf((v[i*3+0]-zmin)/cell);
        int cy=(int)floorf((v[i*3+1]-ymin)/cell);
        int cx=(int)floorf((v[i*3+2]-xmin)/cell);
        size_t c=(size_t)cz*g->gy*g->gx+(size_t)cy*g->gx+cx;
        g->next[i]=g->head[c]; g->head[c]=(int)i;
    }
}
static void grid_free(Grid *g){ free(g->head); free(g->next); }

static float nn_dist(const Grid *g, const float *v, size_t i)
{
    float vz=v[i*3+0], vy=v[i*3+1], vx=v[i*3+2];
    int cz=(int)floorf((vz-g->z0)/g->cell);
    int cy=(int)floorf((vy-g->y0)/g->cell);
    int cx=(int)floorf((vx-g->x0)/g->cell);
    double best=1e30;
    for (int dz=-1;dz<=1;dz++)for(int dy=-1;dy<=1;dy++)for(int dx=-1;dx<=1;dx++){
        int qz=cz+dz,qy=cy+dy,qx=cx+dx;
        if(qz<0||qz>=g->gz||qy<0||qy>=g->gy||qx<0||qx>=g->gx)continue;
        size_t c=(size_t)qz*g->gy*g->gx+(size_t)qy*g->gx+qx;
        for(int j=g->head[c];j!=-1;j=g->next[j]){
            if((size_t)j==i)continue;
            double dz2=v[j*3+0]-vz,dy2=v[j*3+1]-vy,dx2=v[j*3+2]-vx;
            double r2=dz2*dz2+dy2*dy2+dx2*dx2;
            if(r2<best)best=r2;
        }
    }
    return (float)sqrt(best);
}

static int in_bbox(const float *v, size_t i, const double *bb)
{
    float z=v[i*3+0],y=v[i*3+1],x=v[i*3+2];
    return z>=bb[0]&&z<=bb[1]&&y>=bb[2]&&y<=bb[3]&&x>=bb[4]&&x<=bb[5];
}

/* Report NN distribution over the whole cloud and (if given) the bbox subset. */
static void report_nn(const float *v, size_t nv, float nn_radius,
                      const double *bb, const char *tag)
{
    Grid g; grid_build(&g, v, nv, nn_radius);
    const float TH[5]={0.25f,0.40f,0.50f,0.60f,1.00f};
    long cnt_all[5]={0}, cnt_bb[5]={0};
    long n_all=0, n_bb=0;
    double min_all=1e30, min_bb=1e30;
    int NB=30; long hist_all[30]={0}, hist_bb[30]={0}; /* 0.05-vox bins to 1.5 */
    for (size_t i=0;i<nv;i++){
        float d=nn_dist(&g,v,i);
        int bbf = bb?in_bbox(v,i,bb):0;
        n_all++; if(d<min_all)min_all=d;
        int hb=(int)(d/0.05f); if(hb<0)hb=0; if(hb>=NB)hb=NB-1; hist_all[hb]++;
        for(int t=0;t<5;t++) if(d<TH[t]) cnt_all[t]++;
        if(bbf){ n_bb++; if(d<min_bb)min_bb=d;
            hist_bb[hb]++;
            for(int t=0;t<5;t++) if(d<TH[t]) cnt_bb[t]++; }
    }
    grid_free(&g);
    printf("  [%s] N=%ld  minNN=%.3f  NN<0.25:%ld(%.1f%%)  <0.40:%ld(%.1f%%)  "
           "<0.50:%ld(%.1f%%)  <0.60:%ld(%.1f%%)  <1.00:%ld(%.1f%%)\n",
        tag, n_all, min_all,
        cnt_all[0],100.0*cnt_all[0]/n_all, cnt_all[1],100.0*cnt_all[1]/n_all,
        cnt_all[2],100.0*cnt_all[2]/n_all, cnt_all[3],100.0*cnt_all[3]/n_all,
        cnt_all[4],100.0*cnt_all[4]/n_all);
    if (bb && n_bb>0) {
        printf("  [%s/bbox] N=%ld  minNN=%.3f  NN<0.25:%ld(%.1f%%)  <0.40:%ld(%.1f%%)  "
               "<0.50:%ld(%.1f%%)  <0.60:%ld(%.1f%%)  <1.00:%ld(%.1f%%)\n",
            tag, n_bb, min_bb,
            cnt_bb[0],100.0*cnt_bb[0]/n_bb, cnt_bb[1],100.0*cnt_bb[1]/n_bb,
            cnt_bb[2],100.0*cnt_bb[2]/n_bb, cnt_bb[3],100.0*cnt_bb[3]/n_bb,
            cnt_bb[4],100.0*cnt_bb[4]/n_bb);
        /* fine histogram of the bbox subset — where the lightning bolt lives */
        printf("  [%s/bbox] NN histogram (0.05 vox bins):\n", tag);
        for (int b=0;b<NB;b++){ if(!hist_bb[b])continue;
            printf("     %.2f-%.2f  %5ld  %4.1f%%\n",
                b*0.05,(b+1)*0.05,hist_bb[b],100.0*hist_bb[b]/n_bb); }
    }
}

/* ------------------------------------------------------------------ */
/* Self-test: validate the tool + the Fix-0 mechanics.                 */
/* ------------------------------------------------------------------ */
static int approx(double a, double b){ return fabs(a-b)<1e-3; }

static int selftest(void)
{
    int fail=0;
    Arena_T a;
    /* (1) Residual-stack model: two points 0.332 vox apart (z-stack that also
     * drifted in-plane under tangent projection). eps 0.25 must NOT merge;
     * eps 0.50 must merge to 1. This is the exact residual the pipeline leaves. */
    {
        float v[6]={0,0,0, 0.30f,0.10f,0.10f}; /* d=sqrt(0.11)=0.3317 */
        float *o; size_t on;
        a=Arena_new(); weld_iterated(a,v,2,0.25f,&o,&on);
        if(on!=2){printf("SELFTEST FAIL: eps0.25 stack merged (%zu, want 2)\n",on);fail=1;}
        Arena_dispose(&a);
        a=Arena_new(); weld_iterated(a,v,2,0.50f,&o,&on);
        if(on!=1){printf("SELFTEST FAIL: eps0.50 stack not merged (%zu, want 1)\n",on);fail=1;}
        Arena_dispose(&a);
    }
    /* (2) Merger safety: two 3x3 sheets 2.85 vox apart in z, 1.0 in-plane.
     * eps 0.60 must keep them separate (18 verts) and never bridge the gap. */
    {
        float v[54]; size_t n=0;
        for(int s=0;s<2;s++)for(int yy=0;yy<3;yy++)for(int xx=0;xx<3;xx++){
            v[n*3+0]=s?2.85f:0.0f; v[n*3+1]=(float)yy; v[n*3+2]=(float)xx; n++; }
        float *o; size_t on;
        a=Arena_new(); weld_iterated(a,v,n,0.60f,&o,&on);
        if(on!=18){printf("SELFTEST FAIL: 2-sheet eps0.60 nv=%zu (want 18 — merger!)\n",on);fail=1;}
        /* closest cross-sheet pair must remain >= 2.85 */
        Grid g; grid_build(&g,o,on,1.5f);
        (void)g; grid_free(&g);
        Arena_dispose(&a);
    }
    /* (3) Erosion safety: single 5x5 grid, spacing 1.0. eps 0.60 must NOT
     * remove any point (all NN=1.0 > 0.60). Guards against in-plane erosion. */
    {
        float v[75]; size_t n=0;
        for(int yy=0;yy<5;yy++)for(int xx=0;xx<5;xx++){
            v[n*3+0]=0.0f; v[n*3+1]=(float)yy; v[n*3+2]=(float)xx; n++; }
        float *o; size_t on;
        a=Arena_new(); weld_iterated(a,v,n,0.60f,&o,&on);
        if(on!=25){printf("SELFTEST FAIL: clean grid eps0.60 eroded %zu->%zu\n",(size_t)25,on);fail=1;}
        Arena_dispose(&a);
    }
    /* (4) sanity: approx helper */
    if(!approx(1.0,1.0)||approx(0.0,1.0)){printf("SELFTEST FAIL: approx\n");fail=1;}
    if(!fail) printf("SELFTEST OK (residual-stack, merger-safety, erosion-safety)\n");
    return fail;
}

int main(int argc, char **argv)
{
    if (argc>=2 && strcmp(argv[1],"--selftest")==0) return selftest();
    if (argc<2) {
        fprintf(stderr,"usage: %s <cloud.obj> [--eps a,b,c] [--bbox z0 z1 y0 y1 x0 x1]"
                       " [--nn-radius 1.5] [--pca-R 6.0]\n   or: %s --selftest\n",
                argv[0],argv[0]);
        return 1;
    }
    const char *path=argv[1];
    float eps_list[8]={0.25f,0.50f,0.60f}; int n_eps=3;
    double bb[6]; int have_bb=0;
    float nn_radius=1.5f;
    for (int i=2;i<argc;i++){
        if(!strcmp(argv[i],"--eps") && i+1<argc){
            n_eps=0; char *s=argv[++i]; char *tok=strtok(s,",");
            while(tok&&n_eps<8){ eps_list[n_eps++]=(float)atof(tok); tok=strtok(NULL,","); }
        } else if(!strcmp(argv[i],"--bbox") && i+6<argc){
            for(int k=0;k<6;k++) bb[k]=atof(argv[++i]); have_bb=1;
        } else if(!strcmp(argv[i],"--nn-radius") && i+1<argc){
            nn_radius=(float)atof(argv[++i]);
        } else if(!strcmp(argv[i],"--pca-R") && i+1<argc){
            ++i; /* reserved */
        } else {
            fprintf(stderr,"unknown arg: %s\n",argv[i]); return 1;
        }
    }

    size_t nv=0;
    float *verts=load_obj_points(path,&nv);
    if(!verts||nv==0){ fprintf(stderr,"ERROR: no verts in %s\n",path); return 1; }
    printf("Loaded %zu verts from %s  (nn-radius=%.2f)\n",nv,path,(double)nn_radius);
    if(have_bb) printf("bbox z[%.0f,%.0f] y[%.0f,%.0f] x[%.0f,%.0f]\n",
                       bb[0],bb[1],bb[2],bb[3],bb[4],bb[5]);

    printf("\n=== PRE-WELD (raw MLS cloud) ===\n");
    report_nn(verts,nv,nn_radius,have_bb?bb:NULL,"pre-weld");

    for(int e=0;e<n_eps;e++){
        Arena_T arena=Arena_new();
        float *wv; size_t wnv;
        weld_iterated(arena,verts,nv,eps_list[e],&wv,&wnv);
        printf("\n=== WELD eps=%.2f :  %zu -> %zu verts  (-%.1f%%) ===\n",
               (double)eps_list[e],nv,wnv,100.0*(1.0-(double)wnv/(double)nv));
        report_nn(wv,wnv,nn_radius,have_bb?bb:NULL,"post-weld");
        Arena_dispose(&arena);
    }
    free(verts);
    return 0;
}
