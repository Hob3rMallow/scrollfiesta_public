/* atlas_pack.c -- pack many 2D flat-unroll OBJs (from flatten_obj's *_flat.obj)
 * into ONE atlas OBJ. Each input is a planar chart (layout in the first two
 * vertex coords, third ~0). Pieces are shelf-packed into non-overlapping cells
 * (tallest-first rows), translated, and merged with face renumbering. Each piece
 * gets a distinct palette colour so the sheets are tellable apart.
 *
 * Usage:
 *   atlas_pack <out.obj> [--pad P] <flat1.obj> [flat2.obj ...]
 *   atlas_pack --selftest
 */
#include "../common/ves_platform.h"
#include "../common/arena.h"
#include "../common/obj_io.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef struct {
    float   *V;            /* [nv*3], coords 0,1 = layout plane, 2 ~ 0 */
    int32_t *F;            /* [nf*3] */
    size_t   nv, nf;
    double   minx, miny, w, h;   /* layout bbox (coords 0,1) */
    double   ox, oy;             /* placement translation */
} Piece;

static void compute_bbox(Piece *p){
    double mnx=1e300, mny=1e300, mxx=-1e300, mxy=-1e300;
    for(size_t i=0;i<p->nv;i++){
        double x=p->V[i*3+0], y=p->V[i*3+1];
        if(x<mnx)mnx=x; if(x>mxx)mxx=x; if(y<mny)mny=y; if(y>mxy)mxy=y;
    }
    if(p->nv==0){ mnx=mny=0; mxx=mxy=0; }
    p->minx=mnx; p->miny=mny; p->w=mxx-mnx; p->h=mxy-mny;
}

/* tallest-first for tidy shelves */
static int cmp_h_desc(const void *a, const void *b){
    double ha=((const Piece*)a)->h, hb=((const Piece*)b)->h;
    return (ha<hb)?1:((ha>hb)?-1:0);
}

/* shelf-pack: set ox,oy for each piece so cells don't overlap. */
static void pack_pieces(Piece *pieces, size_t n, double pad){
    if(pad<=0) pad=4.0;
    double area=0, maxw=0;
    for(size_t i=0;i<n;i++){ area += (pieces[i].w+pad)*(pieces[i].h+pad);
        if(pieces[i].w>maxw) maxw=pieces[i].w; }
    double target_w = sqrt(area)*1.3; if(target_w<maxw) target_w=maxw;
    qsort(pieces, n, sizeof(Piece), cmp_h_desc);
    double cur_x=0, cur_y=0, row_h=0;
    for(size_t i=0;i<n;i++){
        Piece *p=&pieces[i];
        if(cur_x>0 && cur_x + p->w > target_w){ cur_x=0; cur_y += row_h + pad; row_h=0; }
        p->ox = cur_x - p->minx;
        p->oy = cur_y - p->miny;
        cur_x += p->w + pad;
        if(p->h>row_h) row_h=p->h;
    }
}

static const float PAL[8][3] = {
    {0.90f,0.30f,0.30f},{0.30f,0.65f,0.90f},{0.40f,0.80f,0.40f},{0.90f,0.75f,0.25f},
    {0.70f,0.45f,0.85f},{0.95f,0.55f,0.25f},{0.35f,0.80f,0.75f},{0.85f,0.45f,0.65f},
};

static int write_atlas(const char *out, Piece *pieces, size_t n){
    size_t nv=0, nf=0;
    for(size_t i=0;i<n;i++){ nv+=pieces[i].nv; nf+=pieces[i].nf; }
    float   *V=(float*)malloc((nv?nv:1)*3*sizeof(float));
    int32_t *F=(int32_t*)malloc((nf?nf:1)*3*sizeof(int32_t));
    float   *C=(float*)malloc((nv?nv:1)*3*sizeof(float));
    size_t vo=0, fo=0;
    for(size_t i=0;i<n;i++){
        Piece *p=&pieces[i];
        const float *col=PAL[i%8];
        for(size_t k=0;k<p->nv;k++){
            V[(vo+k)*3+0]=(float)(p->V[k*3+0]+p->ox);
            V[(vo+k)*3+1]=(float)(p->V[k*3+1]+p->oy);
            V[(vo+k)*3+2]=0.0f;
            C[(vo+k)*3+0]=col[0]; C[(vo+k)*3+1]=col[1]; C[(vo+k)*3+2]=col[2];
        }
        for(size_t k=0;k<p->nf;k++){
            F[(fo+k)*3+0]=p->F[k*3+0]+(int32_t)vo;
            F[(fo+k)*3+1]=p->F[k*3+1]+(int32_t)vo;
            F[(fo+k)*3+2]=p->F[k*3+2]+(int32_t)vo;
        }
        vo+=p->nv; fo+=p->nf;
    }
    ves_ensure_parent_dir(out);
    int rc=ObjIO_write_per_vertex_color(out, V, nv, F, nf, C);
    free(V); free(F); free(C);
    return rc;
}

static int run_selftest(void){
    /* three unit squares, all initially at the origin (overlapping). After pack
     * their cells must be disjoint -> combined bbox area >= sum of cell areas. */
    Piece p[3];
    for(int i=0;i<3;i++){
        p[i].nv=4; p[i].nf=2;
        p[i].V=(float*)malloc(4*3*sizeof(float));
        float sq[12]={0,0,0, 1,0,0, 1,1,0, 0,1,0};
        memcpy(p[i].V, sq, sizeof sq);
        p[i].F=(int32_t*)malloc(2*3*sizeof(int32_t));
        int32_t f[6]={0,1,2, 0,2,3}; memcpy(p[i].F,f,sizeof f);
        compute_bbox(&p[i]);
    }
    pack_pieces(p, 3, 1.0);
    /* verify pairwise non-overlap of placed cells */
    int overlap=0;
    for(int i=0;i<3;i++) for(int j=i+1;j<3;j++){
        double ax0=p[i].ox+p[i].minx, ay0=p[i].oy+p[i].miny, ax1=ax0+p[i].w, ay1=ay0+p[i].h;
        double bx0=p[j].ox+p[j].minx, by0=p[j].oy+p[j].miny, bx1=bx0+p[j].w, by1=by0+p[j].h;
        if(ax0<bx1 && bx0<ax1 && ay0<by1 && by0<ay1) overlap=1;
    }
    const char *out="output/_atlas_pack_test.obj";
    int rc=write_atlas(out, p, 3);
    Arena_T ar=Arena_new(); float*V=NULL;int32_t*F=NULL;size_t nv=0,nf=0;
    int rr=ObjIO_read(ar,out,&V,&nv,&F,&nf);
    int fail = overlap || rc || rr || nv!=12 || nf!=6;
    fprintf(stderr,"[selftest] pack 3 squares: overlap=%d nv=%zu(want12) nf=%zu(want6) -> %s\n",
            overlap, nv, nf, fail?"FAIL":"ok");
    Arena_dispose(&ar); remove(out);
    for(int i=0;i<3;i++){ free(p[i].V); free(p[i].F); }
    fprintf(stderr,"=== atlas_pack selftest %s ===\n", fail?"FAILED":"PASSED");
    return fail?1:0;
}

int main(int argc, char **argv){
    if(argc>=2 && !strcmp(argv[1],"--selftest")) return run_selftest();
    if(argc<3){
        fprintf(stderr,"Usage: %s <out.obj> [--pad P] <flat1.obj> [flat2.obj ...]\n"
                       "       %s --selftest\n", argv[0], argv[0]);
        return 1;
    }
    const char *out=argv[1];
    double pad=4.0;
    int a=2;
    if(a<argc && !strcmp(argv[a],"--pad")){ pad=atof(argv[a+1]); a+=2; }
    int nin=argc-a;
    if(nin<1){ fprintf(stderr,"no input OBJs\n"); return 1; }
    Piece *pieces=(Piece*)calloc((size_t)nin, sizeof(Piece));
    size_t np=0;
    for(int i=a;i<argc;i++){
        Arena_T ar=Arena_new();
        float *v=NULL; int32_t *f=NULL; size_t lnv=0,lnf=0;
        if(ObjIO_read(ar,argv[i],&v,&lnv,&f,&lnf)!=0 || lnv==0){
            fprintf(stderr,"atlas_pack: skip %s\n",argv[i]); Arena_dispose(&ar); continue;
        }
        Piece *p=&pieces[np];
        p->nv=lnv; p->nf=lnf;
        p->V=(float*)malloc(lnv*3*sizeof(float)); memcpy(p->V, v, lnv*3*sizeof(float));
        p->F=(int32_t*)malloc((lnf?lnf:1)*3*sizeof(int32_t));
        if(lnf) memcpy(p->F, f, lnf*3*sizeof(int32_t));
        compute_bbox(p);
        np++;
        Arena_dispose(&ar);
    }
    pack_pieces(pieces, np, pad);
    int rc=write_atlas(out, pieces, np);
    double aw=0, ah=0;
    for(size_t i=0;i<np;i++){ double r=pieces[i].ox+pieces[i].minx+pieces[i].w; if(r>aw)aw=r;
        double t=pieces[i].oy+pieces[i].miny+pieces[i].h; if(t>ah)ah=t; }
    fprintf(stderr,"atlas_pack: %zu/%d charts -> %.0f x %.0f atlas -> %s (%s)\n",
            np, nin, aw, ah, out, rc?"FAIL":"ok");
    for(size_t i=0;i<np;i++){ free(pieces[i].V); free(pieces[i].F); }
    free(pieces);
    return rc?1:0;
}
