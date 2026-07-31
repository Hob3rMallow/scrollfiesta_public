/* ============================================================================
 * surface_snap.c -- DETECTION stage. See surface_snap.h. Moves nothing; only
 * classifies dark vertices (off-surface FIXABLE vs real CRACK) and groups them
 * into small regions, for inspection before a snap solver is built.
 * ==========================================================================*/

#include "surface_snap.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../common/csr.h"
#include "../common/kdtree.h"
#include "../common/raw_sample.h"
#include "../common/tiff_io.h"
#include "../common/ves_platform.h"
#include "snap_darkcut.h"
#include "snap_quilt.h"

/* ------------------------------------------------------------------ util */
static void *xmalloc(size_t n){void*p=malloc(n?n:1);if(!p){fprintf(stderr,"snap OOM\n");exit(1);}return p;}
static void *xcalloc(size_t c,size_t s){void*p=calloc(c?c:1,s);if(!p){fprintf(stderr,"snap OOM\n");exit(1);}return p;}

typedef struct { int32_t *p; } UF;
static void uf_init(UF*u,size_t n){u->p=(int32_t*)xmalloc(n*sizeof(int32_t));for(size_t i=0;i<n;i++)u->p[i]=(int32_t)i;}
static int32_t uf_find(UF*u,int32_t x){while(u->p[x]!=x){u->p[x]=u->p[u->p[x]];x=u->p[x];}return x;}
static void uf_union(UF*u,int32_t a,int32_t b){int32_t ra=uf_find(u,a),rb=uf_find(u,b);if(ra!=rb)u->p[rb]=ra;}
static void uf_free(UF*u){free(u->p);u->p=NULL;}

/* intensity window: 256-bin histogram over the sampled population, percentile
 * cut points (mirrors obj_bake_raw stretch_window). */
static void stretch_window(const double *val, const uint8_t *has, size_t nv,
                           double pct_lo, double pct_hi, double *out_lo, double *out_hi){
    size_t hist[256]; memset(hist,0,sizeof hist);
    size_t n=0,cum=0,tlo,thi; int b,lo=0,hi=255;
    for(size_t i=0;i<nv;i++){ if(!has[i])continue; int bb=(int)(val[i]+0.5); if(bb<0)bb=0; if(bb>255)bb=255; hist[bb]++; n++; }
    *out_lo=0; *out_hi=255; if(n==0)return;
    tlo=(size_t)(pct_lo/100.0*(double)n); thi=(size_t)(pct_hi/100.0*(double)n);
    cum=0; for(b=0;b<256;b++){cum+=hist[b]; if(cum>tlo){lo=b;break;}}
    cum=0; for(b=0;b<256;b++){cum+=hist[b]; if(cum>=thi){hi=b;break;}}
    if(hi<=lo){lo=0;hi=255;} *out_lo=lo; *out_hi=hi;
}

void SnapOpts_default(SnapOpts *o){
    memset(o,0,sizeof *o);
    o->axis_point[1]=3405.0f; o->axis_point[2]=2878.0f; o->axis_dir[0]=1.0f;
    o->chunk=128; o->reach=8.0; o->occ_thresh=0.9; o->local_r=2.0; o->gain_bias=1.0;
    o->tensor_weight=0.0; o->tensor_radius=2.0;
    o->step=0.25; o->band_frac=0.30; o->min_gain=25.0;
    o->anchor_frac=0.25; o->region_cap=4000; o->min_region=8; o->pct_lo=1.0; o->pct_hi=99.0;
    o->dark_thresh=64; o->dark_lambda=20.0; o->dark_sigma=30.0;   /* inclusive: catch dim gray, FP>FN */
    o->normal_range=2.0; o->normal_samples=5; o->verbose=0;
    o->target_mode=SNAP_TARGET_QUILT_MRF; o->depth_bins=33;
    o->w_match=30.0; o->w_close=200.0;
    o->smooth_mu=2000.0; o->smooth_tau=4.0;
}

int SnapDetect_run(Arena_T arena,
                   const float *verts, size_t nv,
                   const int32_t *faces, size_t nf,
                   const float *uv,
                   const SnapOpts *opts, SnapResult *out){
    memset(out,0,sizeof *out);
    if(nv<3||nf<1||uv==NULL) return -1;
    SnapOpts o=*opts;
    double reach=o.reach>0?o.reach:8.0;                 /* march ceiling (safety cap) */
    CubeTable *ct=(CubeTable*)ARENA_ALLOC(arena,(long)sizeof(CubeTable));
    if(o.raw_dir==NULL || cubetable_init(ct,arena,o.raw_dir,o.chunk,verts,nv,o.normal_range+reach+2.0)!=0)
        return -1;
    float *nrm=vertex_normals(verts,nv,faces,nf);
    /* occupancy index over the EXISTING mesh -- "have we entered another sheet?" */
    KDTree_T occ=KDTree_new(arena,verts,nv);

    /* Geometric truth must be the intensity AT the vertex. A normal-max sample
     * lets an off-sheet vertex in a dark gap borrow brightness from papyrus up
     * to normal_range voxels away, hiding exactly the defects snap must find. */
    double *cv=(double*)xmalloc(nv*sizeof(double));
    uint8_t *has=(uint8_t*)xcalloc(nv,1);
    for(size_t i=0;i<nv;i++){
        double s=sample_trilinear(ct,verts[i*3+0],verts[i*3+1],verts[i*3+2]);
        cv[i]=s>=0?s:0.0; has[i]=s>=0?1:0;
    }
    double lo,hi; stretch_window(cv,has,nv,o.pct_lo,o.pct_hi,&lo,&hi);
    double dk=lo+((double)o.dark_thresh/255.0)*(hi-lo);
    double band=lo+o.band_frac*(hi-lo);

    /* adjacency (built once; also used to skip ORPHAN verts -- degree 0, i.e.
     * referenced by no kept face, e.g. left behind by overlap-drop -- which are
     * not part of the surface and must not count as bad). */
    CSR_T adj=CSR_from_faces(arena,faces,nf,nv);
    const int32_t *aoff=CSR_offset(adj), *atgt=CSR_target(adj);

    /* DARK segmentation: a binary Boykov-Jolly graph cut (GCO) grows dark cores
     * to their true extent instead of a hard per-vertex threshold. dark_lambda=0
     * (or a GCO failure) falls back to the exact hard cutoff cv<dk. */
    uint8_t *dark=(uint8_t*)xmalloc(nv);
    double band_hi=band+0.5*(hi-band);   /* bright anchor: never flag clearly-bright as dark */
    { SnapDarkOpts dopt; dopt.dk=dk; dopt.lambda=o.dark_lambda; dopt.sigma=o.dark_sigma; dopt.band_hi=band_hi;
      if(o.dark_lambda<=0.0 || SnapDark_segment(cv,has,nv,adj,&dopt,dark)!=0)
          for(size_t i=0;i<nv;i++) dark[i]=(has[i]&&cv[i]<dk)?(uint8_t)1:(uint8_t)0; }
    for(size_t i=0;i<nv;i++) if(dark[i]&&cv[i]>band) out->n_overgrow++;  /* diag: bright-ish in dark set */

    /* outputs */
    out->vclass=(uint8_t*)ARENA_CALLOC(arena,(long)nv,1);
    out->voff=(float*)ARENA_CALLOC(arena,(long)nv,sizeof(float));
    out->vdir=(float*)ARENA_CALLOC(arena,(long)nv,3*sizeof(float));
    out->vgain=(float*)ARENA_CALLOC(arena,(long)nv,sizeof(float));
    out->vtensor=(float*)ARENA_CALLOC(arena,(long)nv,sizeof(float));
    out->vblock=(float*)ARENA_CALLOC(arena,(long)nv,sizeof(float));
    out->vcur=(float*)ARENA_ALLOC(arena,(long)(nv*sizeof(float)));
    out->vkdefect=(float*)ARENA_CALLOC(arena,(long)nv,sizeof(float));
    out->vregion=(int32_t*)ARENA_ALLOC(arena,(long)(nv*sizeof(int32_t)));
    for(size_t i=0;i<nv;i++){ out->vcur[i]=(float)cv[i]; out->vregion[i]=-1; }

    /* Pass 1 is deliberately a REPAIR pass, not final surface placement.  It
     * moves only obvious dark whiffs.  Candidate brightness is an eligibility
     * gate; selection minimizes clean-boundary quilting mismatch + snap distance
     * and jointly regularizes signed depth over the dark 1-ring. */
    {
        SnapQuiltOpts qopt;
        SnapQuiltStats qst;
        memset(&qopt,0,sizeof qopt);
        memcpy(qopt.axis_point,o.axis_point,sizeof qopt.axis_point);
        memcpy(qopt.axis_dir,o.axis_dir,sizeof qopt.axis_dir);
        qopt.reach=reach; qopt.step=o.step;
        qopt.occ_thresh=o.occ_thresh; qopt.local_r=o.local_r;
        qopt.band=band; qopt.min_gain=o.min_gain;
        qopt.target_mode=o.target_mode; qopt.depth_bins=o.depth_bins;
        qopt.w_match=o.w_match; qopt.w_close=o.w_close;
        qopt.smooth_mu=o.smooth_mu; qopt.smooth_tau=o.smooth_tau;
        qopt.verbose=o.verbose;
        if(SnapQuilt_select(arena,ct,occ,verts,0,verts,nrm,nv,adj,cv,has,dark,
                            &qopt,out->vclass,out->voff,out->vdir,out->vgain,
                            out->vblock,&qst)!=0){
            free(nrm);free(cv);free(has);free(dark);
            return -1;
        }
        out->n_dark=qst.n_dark; out->n_fixable=qst.n_fixable;
        out->n_crack=qst.n_crack; out->n_bidir=qst.n_bidir;
        out->n_blocked=qst.n_blocked;
        out->quilt_cost_mean=qst.mean_quilt_cost;
        out->target_dist_mean=qst.mean_distance;
        out->n_quilt_fallback=qst.n_gco_fallback;
        out->raw_table=ct;
    }
    out->n_good=nv-out->n_dark;

    /* angle-defect K_v (interior approx: 2pi - sum incident angles; boundary
     * verts read high, noted) -- geometry-anomaly cross-signal */
    { double *asum=(double*)xcalloc(nv,sizeof(double));
      for(size_t f=0;f<nf;f++){ for(int k=0;k<3;k++){
          size_t a=(size_t)faces[f*3+k],b=(size_t)faces[f*3+(k+1)%3],c=(size_t)faces[f*3+(k+2)%3];
          double e1v[3],e2v[3],l1=0,l2=0,dot=0;
          for(int m=0;m<3;m++){e1v[m]=verts[b*3+m]-verts[a*3+m];e2v[m]=verts[c*3+m]-verts[a*3+m];l1+=e1v[m]*e1v[m];l2+=e2v[m]*e2v[m];dot+=e1v[m]*e2v[m];}
          l1=sqrt(l1);l2=sqrt(l2); if(l1>1e-9&&l2>1e-9){double cc=dot/(l1*l2); if(cc>1)cc=1; if(cc<-1)cc=-1; asum[a]+=acos(cc);}
      }}
      for(size_t i=0;i<nv;i++) out->vkdefect[i]=(float)(6.283185307-asum[i]);
      free(asum);
    }

    /* group dark verts into regions (dev_cut pattern: union dark-dark edges) */
    UF uf; uf_init(&uf,nv);
    for(size_t i=0;i<nv;i++) if(out->vclass[i]!=SNAP_GOOD)
        for(int32_t e=aoff[i];e<aoff[i+1];e++){ int32_t j=atgt[e];
            if(out->vclass[j]!=SNAP_GOOD) uf_union(&uf,(int32_t)i,j); }
    int32_t nreg=0; int32_t *root2rid=(int32_t*)xmalloc(nv*sizeof(int32_t));
    for(size_t i=0;i<nv;i++) root2rid[i]=-1;
    for(size_t i=0;i<nv;i++){ if(out->vclass[i]==SNAP_GOOD)continue;
        int32_t r=uf_find(&uf,(int32_t)i); if(root2rid[r]<0)root2rid[r]=nreg++; out->vregion[i]=root2rid[r]; }

    /* per-region tallies + good-boundary fraction */
    int32_t *rsize=(int32_t*)xcalloc((size_t)nreg,sizeof(int32_t));
    int32_t *rfix=(int32_t*)xcalloc((size_t)nreg,sizeof(int32_t));
    int32_t *rcrk=(int32_t*)xcalloc((size_t)nreg,sizeof(int32_t));
    int32_t *rbnd=(int32_t*)xcalloc((size_t)nreg,sizeof(int32_t));
    int32_t *rgood=(int32_t*)xcalloc((size_t)nreg,sizeof(int32_t));
    for(size_t i=0;i<nv;i++){ int32_t r=out->vregion[i]; if(r<0)continue;
        rsize[r]++; if(out->vclass[i]==SNAP_FIXABLE)rfix[r]++; else rcrk[r]++;
        for(int32_t e=aoff[i];e<aoff[i+1];e++){ int32_t j=atgt[e];
            if(out->vregion[j]!=r){ rbnd[r]++; if(out->vclass[j]==SNAP_GOOD)rgood[r]++; } } }

    out->rclass=(uint8_t*)ARENA_ALLOC(arena,(long)(nreg>0?nreg:1));
    out->rsize=(int32_t*)ARENA_ALLOC(arena,(long)((nreg>0?nreg:1)*sizeof(int32_t)));
    out->nreg=(size_t)nreg;
    for(int32_t r=0;r<nreg;r++){
        out->rsize[r]=rsize[r]; if((size_t)rsize[r]>out->max_region_size)out->max_region_size=(size_t)rsize[r];
        double gfrac=rbnd[r]>0?(double)rgood[r]/(double)rbnd[r]:0.0;
        double ffrac=rsize[r]>0?(double)rfix[r]/(double)rsize[r]:0.0;
        uint8_t rc;
        if(o.min_region>1 && rsize[r]<o.min_region) rc=SNAPREG_REJECT;   /* speckle -> drop */
        else if(gfrac<o.anchor_frac || rsize[r]>o.region_cap) rc=SNAPREG_ANCHORLESS;
        else if(ffrac>=0.5) rc=SNAPREG_FIXABLE;
        else if((double)rcrk[r]/(double)rsize[r]>=0.5) rc=SNAPREG_CRACK;
        else rc=SNAPREG_MIXED;
        out->rclass[r]=rc;
        switch(rc){case SNAPREG_FIXABLE:out->n_reg_fixable++;break;case SNAPREG_CRACK:out->n_reg_crack++;break;
                   case SNAPREG_ANCHORLESS:out->n_reg_anchorless++;break;
                   case SNAPREG_REJECT:out->n_reg_reject++;break;default:out->n_reg_mixed++;}
    }
    /* min-region reject: verts in dropped regions return to GOOD (burn-artifact
     * speckle; see feedback). vregion kept so the reject is still inspectable. */
    for(size_t i=0;i<nv;i++){ int32_t r=out->vregion[i];
        if(r>=0 && out->rclass[r]==SNAPREG_REJECT && out->vclass[i]!=SNAP_GOOD){
            if(out->vclass[i]==SNAP_FIXABLE)out->n_fixable--; else out->n_crack--;
            out->vclass[i]=SNAP_GOOD; out->n_dark--; out->n_good++; out->n_rejected++;
        } }
    out->win_lo=lo; out->win_hi=hi;

    if(o.verbose)
        fprintf(stderr,"[snap] verts=%zu dark=%zu (fixable=%zu [bidir=%zu] crack=%zu [blocked=%zu]) "
                "regions=%d (fix-island=%zu crack=%zu anchorless=%zu mixed=%zu reject=%zu[-%zu v]) "
                "max_region=%zu quilt=%.1f target_dist=%.2f qfb=%zu "
                "window[%.0f,%.0f] dk=%.1f band=%.1f reach=%.1f\n",
                nv,out->n_dark,out->n_fixable,out->n_bidir,out->n_crack,out->n_blocked,nreg,
                out->n_reg_fixable,out->n_reg_crack,out->n_reg_anchorless,out->n_reg_mixed,
                out->n_reg_reject,out->n_rejected,
                out->max_region_size,out->quilt_cost_mean,out->target_dist_mean,
                out->n_quilt_fallback,lo,hi,dk,band,reach);

    free(nrm);free(cv);free(has);free(dark);uf_free(&uf);free(root2rid);
    free(rsize);free(rfix);free(rcrk);free(rbnd);free(rgood);
    return 0;
}

/* ============================================================================
 * Self-test (injected synthetic CubeTable, like obj_bake_raw)
 * ==========================================================================*/
#define CK(c,m) do{ if(!(c)){fprintf(stderr,"  FAIL: %s\n",(m));nf++;} else fprintf(stderr,"  ok: %s\n",(m)); }while(0)

/* Build a flat grid patch: verts at (z=z0+j, y=py, x=x0+i), a nu x nv_ grid.
 * radial (about Z through origin) at these verts ~ +y. Appends. */
static void st_patch(double z0,double py,double x0,int nu,int nvv,
                     float **V,int32_t **F,size_t *nv,size_t *nfc,size_t *cv,size_t *cf){
    size_t base=*nv;
    for(int j=0;j<=nvv;j++)for(int i=0;i<=nu;i++){
        if(*nv+1>*cv){*cv=(*cv?*cv*2:1024);*V=(float*)realloc(*V,*cv*3*sizeof(float));}
        /* tiny deterministic jitter (<0.03 vox) so no two verts share an exact
         * coordinate -- a perfect integer lattice is a KD-tree build degeneracy
         * (many equal split-axis values); real LOP-projected meshes are generic
         * floats, so this keeps the synthetic test representative. */
        unsigned hh=(unsigned)(*nv)*2654435761u;
        double jz=((double)((hh    )&255)/255.0-0.5)*0.05;
        double jy=((double)((hh>> 8)&255)/255.0-0.5)*0.05;
        double jx=((double)((hh>>16)&255)/255.0-0.5)*0.05;
        (*V)[*nv*3+0]=(float)(z0+j+jz);(*V)[*nv*3+1]=(float)(py+jy);(*V)[*nv*3+2]=(float)(x0+i+jx);(*nv)++; }
    for(int j=0;j<nvv;j++)for(int i=0;i<nu;i++){
        int32_t a=(int32_t)(base+(size_t)j*(nu+1)+i),b=a+1,c=a+nu+1,d=c+1;
        if(*nfc+2>*cf){*cf=(*cf?*cf*2:1024);*F=(int32_t*)realloc(*F,*cf*3*sizeof(int32_t));}
        (*F)[*nfc*3+0]=a;(*F)[*nfc*3+1]=b;(*F)[*nfc*3+2]=c;(*nfc)++;
        (*F)[*nfc*3+0]=b;(*F)[*nfc*3+1]=d;(*F)[*nfc*3+2]=c;(*nfc)++; }
}
static void st_uv(int nu,int nvv,double u0,float **T,size_t *nvt,size_t *ct){
    for(int j=0;j<=nvv;j++)for(int i=0;i<=nu;i++){
        if(*nvt+1>*ct){*ct=(*ct?*ct*2:1024);*T=(float*)realloc(*T,*ct*2*sizeof(float));}
        (*T)[*nvt*2+0]=(float)(u0+i);(*T)[*nvt*2+1]=(float)j;(*nvt)++; }
}

int Snap_selftest(void){
    int nf=0; Arena_T arena=Arena_new();
    fprintf(stderr,"[selftest] surface_snap\n");
    const char *dir="output/surface_snap_selftest";

    /* one shared 32^3 cube: dark (20) except a bright slab (200) at y in [8,10];
     * written to disk so cubetable_init can load it. */
    long chunk=32; uint8_t *buf=(uint8_t*)malloc((size_t)chunk*chunk*chunk);
    for(long z=0;z<chunk;z++)for(long y=0;y<chunk;y++)for(long x=0;x<chunk;x++)
        buf[((z*chunk+y)*chunk+x)]=(uint8_t)((y>=8&&y<=10)?200:20);
    { char path[512]; snprintf(path,sizeof path,"%s/z00000_y00000_x00000.tif",dir);
      ves_ensure_parent_dir(path); TiffIO_save(path,buf,(int)chunk,(int)chunk,(int)chunk); }

    SnapOpts o; SnapOpts_default(&o);
    o.axis_point[0]=0;o.axis_point[1]=0;o.axis_point[2]=8;   /* axis Z through x=8 -> radial ~ +y */
    o.axis_dir[0]=1;o.axis_dir[1]=0;o.axis_dir[2]=0;
    o.raw_dir=dir; o.chunk=chunk; o.reach=8.0; o.min_gain=25;
    o.normal_range=2.0; o.normal_samples=5; o.verbose=0;

    /* t1: BRIGHT context patch (on the slab, y=9) + a DARK patch near the slab
     * (y=7, 1-3 vox below). The window is set by the bright context; the dark
     * patch's radial (~+y) probe reaches the slab -> FIXABLE, offset ~ +1..+3. */
    {
        float *V=NULL,*T=NULL;int32_t *F=NULL;size_t nv=0,nvt=0,nfc=0,cvv=0,ctt=0,cff=0;
        st_patch(2.0, 9.0, 4.0, 8,4, &V,&F,&nv,&nfc,&cvv,&cff); st_uv(8,4, 0,&T,&nvt,&ctt); /* bright */
        st_patch(20.0,7.0, 4.0, 8,4, &V,&F,&nv,&nfc,&cvv,&cff); st_uv(8,4,40,&T,&nvt,&ctt); /* dark  */
        SnapResult R; int rc=SnapDetect_run(arena,V,nv,F,nfc,T,&o,&R);
        CK(rc==0,"t1 run ok");
        CK(R.n_dark>0,"t1 dark patch flagged (relative to bright context)");
        CK(R.n_fixable>0,"t1 fixable found (slab reachable along radial)");
        int found=0; for(size_t i=0;i<nv;i++) if(R.vclass[i]==SNAP_FIXABLE && R.voff[i]>=0.5f && R.voff[i]<=3.5f) found=1;
        CK(found,"t1 fixable offset positive toward slab (<= reach)");
        int honest=0,near_target=0; for(size_t i=0;i<nv;i++) if(R.vclass[i]==SNAP_FIXABLE){
            if(R.vcur[i]<40.0f)honest=1;
            if(R.voff[i]<=1.5f)near_target=1;
        }
        CK(honest,"t1 current intensity is sampled at the vertex (no borrowed normal-max brightness)");
        CK(near_target,"t1 repair balances quilting against distance (no core-seeking overshoot)");
        free(V);free(T);free(F);
    }
    /* t2: bright context + DARK patch far from the slab (y=22, ~12 vox away) ->
     * no target within reach -> CRACK. */
    {
        float *V=NULL,*T=NULL;int32_t *F=NULL;size_t nv=0,nvt=0,nfc=0,cvv=0,ctt=0,cff=0;
        st_patch(2.0, 9.0, 4.0, 8,4, &V,&F,&nv,&nfc,&cvv,&cff); st_uv(8,4, 0,&T,&nvt,&ctt); /* bright */
        st_patch(20.0,22.0,4.0, 8,4, &V,&F,&nv,&nfc,&cvv,&cff); st_uv(8,4,40,&T,&nvt,&ctt); /* dark far */
        SnapResult R; int rc=SnapDetect_run(arena,V,nv,F,nfc,T,&o,&R);
        CK(rc==0,"t2 run ok");
        CK(R.n_dark>0 && R.n_fixable==0 && R.n_crack>0,"t2 slab out of reach -> CRACK (no false snap)");
        free(V);free(T);free(F);
    }
    /* t3: bright context + two SEPARATE dark patches -> two dark regions */
    {
        float *V=NULL,*T=NULL;int32_t *F=NULL;size_t nv=0,nvt=0,nfc=0,cvv=0,ctt=0,cff=0;
        st_patch(2.0, 9.0, 4.0, 8,4, &V,&F,&nv,&nfc,&cvv,&cff); st_uv(8,4, 0,&T,&nvt,&ctt);
        st_patch(18.0,22.0,4.0, 4,3, &V,&F,&nv,&nfc,&cvv,&cff); st_uv(4,3,40,&T,&nvt,&ctt);
        st_patch(26.0,22.0,4.0, 4,3, &V,&F,&nv,&nfc,&cvv,&cff); st_uv(4,3,60,&T,&nvt,&ctt);
        SnapResult R; int rc=SnapDetect_run(arena,V,nv,F,nfc,T,&o,&R);
        CK(rc==0 && R.nreg>=2,"t3 two separate dark patches -> >=2 regions");
        free(V);free(T);free(F);
    }
    /* t5: FAR-but-FREE. Dark source ~5 vox below the slab, nothing meshed in
     * between -> reachable in free space -> FIXABLE even though > 3 vox (the
     * whole point: distance is not the guard, occupancy is). */
    {
        float *V=NULL,*T=NULL;int32_t *F=NULL;size_t nv=0,nvt=0,nfc=0,cvv=0,ctt=0,cff=0;
        st_patch(2.0, 9.0, 4.0, 8,4, &V,&F,&nv,&nfc,&cvv,&cff); st_uv(8,4, 0,&T,&nvt,&ctt); /* bright ctx (z2-6) */
        size_t base_src=nv;
        st_patch(20.0,3.0, 4.0, 8,4, &V,&F,&nv,&nfc,&cvv,&cff); st_uv(8,4,40,&T,&nvt,&ctt); /* dark src y3 (z20-28) */
        SnapResult R; int rc=SnapDetect_run(arena,V,nv,F,nfc,T,&o,&R);
        int far_fix=0; for(size_t i=base_src;i<nv;i++)
            if(R.vclass[i]==SNAP_FIXABLE && R.voff[i]>=4.0f
               && R.voff[i]<=(float)(o.reach+0.1)) far_fix=1;
        CK(rc==0,"t5 run ok");
        CK(far_fix,"t5 far (~5 vox) free-space target -> FIXABLE (distance is not the guard)");
        free(V);free(T);free(F);
    }
    /* t6: OCCUPANCY BLOCK. Same far slab, but another mesh SHEET sits between the
     * dark source and the slab -> the probe halts at that sheet -> the shadowed
     * source verts are CRACK, never snapped through an existing sheet. */
    {
        float *V=NULL,*T=NULL;int32_t *F=NULL;size_t nv=0,nvt=0,nfc=0,cvv=0,ctt=0,cff=0;
        st_patch(2.0, 9.0, 4.0, 8,4, &V,&F,&nv,&nfc,&cvv,&cff); st_uv(8,4, 0,&T,&nvt,&ctt); /* bright ctx */
        size_t base_src=nv;
        st_patch(20.0,3.0, 6.0, 4,4, &V,&F,&nv,&nfc,&cvv,&cff); st_uv(4,4,40,&T,&nvt,&ctt); /* dark src y3 x6-10 (narrow) */
        size_t end_src=nv;
        st_patch(20.0,6.0, 2.0, 12,4, &V,&F,&nv,&nfc,&cvv,&cff); st_uv(12,4,80,&T,&nvt,&ctt); /* WIDE blocker sheet y6 x2-14 */
        SnapResult R; int rc=SnapDetect_run(arena,V,nv,F,nfc,T,&o,&R);
        int all_crack=1; for(size_t i=base_src;i<end_src;i++) if(R.vclass[i]==SNAP_FIXABLE) all_crack=0;
        CK(rc==0,"t6 run ok");
        CK(all_crack && R.n_blocked>0,"t6 slab behind another sheet -> shadowed src CRACK (occupancy block)");
        free(V);free(T);free(F);
    }
    /* t7/t8: DARK graph cut. Concentric 11x11 patch -- dark core (20) + dim ring
     * (55, just ABOVE dk so the hard threshold misses it) + bright surround (200).
     * lambda>0 must GROW the ring into DARK without bleeding into the surround;
     * lambda=0 must reduce EXACTLY to cv<dk (only the core). */
    {
        float *V=NULL,*T=NULL;int32_t *F=NULL;size_t nvv=0,nvt=0,nfc=0,cvc=0,ctt=0,cff=0;
        int N=10; (void)T;(void)nvt;(void)ctt;
        st_patch(0.0,0.0,0.0,N,N,&V,&F,&nvv,&nfc,&cvc,&cff);
        double *cvi=(double*)malloc(nvv*sizeof(double)); uint8_t *hasv=(uint8_t*)malloc(nvv);
        for(size_t k=0;k<nvv;k++){ int i=(int)(k%(size_t)(N+1)),j=(int)(k/(size_t)(N+1));
            int cheb=(abs(i-5)>abs(j-5))?abs(i-5):abs(j-5);
            cvi[k]=(cheb<=1)?20.0:(cheb==2)?55.0:200.0; hasv[k]=1; }
        CSR_T adj2=CSR_from_faces(arena,F,nfc,nvv);
        uint8_t *dk2=(uint8_t*)malloc(nvv);
        SnapDarkOpts dop; dop.dk=50.0; dop.lambda=12.0; dop.sigma=40.0; dop.band_hi=0.0;
        int rc7=SnapDark_segment(cvi,hasv,nvv,adj2,&dop,dk2);
        int core_d=0,core_n=0,ring_d=0,ring_n=0,surr_d=0;
        for(size_t k=0;k<nvv;k++){ int i=(int)(k%(size_t)(N+1)),j=(int)(k/(size_t)(N+1));
            int cheb=(abs(i-5)>abs(j-5))?abs(i-5):abs(j-5);
            if(cheb<=1){core_n++;core_d+=dk2[k];} else if(cheb==2){ring_n++;ring_d+=dk2[k];} else surr_d+=dk2[k]; }
        CK(rc7==0,"t7 darkcut run ok");
        CK(core_d==core_n,"t7 dark core all DARK");
        CK(ring_d>=ring_n*3/4,"t7 dim ring GROWN into DARK (cut beats threshold)");
        CK(surr_d==0,"t7 bright surround NOT flagged (no bleed)");
        dop.lambda=0.0; int rc8=SnapDark_segment(cvi,hasv,nvv,adj2,&dop,dk2);
        int par=1; for(size_t k=0;k<nvv;k++){ int want=(cvi[k]<50.0)?1:0; if(dk2[k]!=want)par=0; }
        CK(rc8==0 && par,"t8 lambda=0 == hard threshold cv<dk (parity)");
        free(cvi);free(hasv);free(dk2);free(V);free(F);
    }
    /* t9: two reachable bright slabs.  The old intensity/tensor selector could
     * prefer the farther response and jump to the verso.  Quilting mode treats
     * brightness only as eligibility, balances clean-boundary match vs distance,
     * and is invariant to the retired tensor knob. */
    {
        static const uint8_t fiber[8]={160,181,190,181,160,139,130,139};
        for(long z=0;z<chunk;z++)for(long y=0;y<chunk;y++)for(long x=0;x<chunk;x++){
            uint8_t v=20;
            if(y>=6&&y<=8)v=200;              /* brighter, tangent-flat */
            if(y>=10&&y<=12)v=fiber[z&7];     /* dimmer, coherent axial texture */
            buf[((z*chunk+y)*chunk+x)]=v;
        }
        { char path[512]; snprintf(path,sizeof path,"%s/z00000_y00000_x00000.tif",dir);
          TiffIO_save(path,buf,(int)chunk,(int)chunk,(int)chunk); }

        float *V=NULL,*T=NULL;int32_t *F=NULL;size_t nv=0,nvt=0,nfc=0,cvv=0,ctt=0,cff=0;
        o.axis_point[2]=7.5f;
        st_patch(2.0,7.0,7.0,1,8,&V,&F,&nv,&nfc,&cvv,&cff); st_uv(1,8,0,&T,&nvt,&ctt);
        size_t base_src=nv;
        st_patch(20.0,3.0,7.0,1,8,&V,&F,&nv,&nfc,&cvv,&cff); st_uv(1,8,40,&T,&nvt,&ctt);
        SnapResult R0,R1;
        o.reach=9.0; o.min_region=0; o.tensor_radius=2.0;
        o.tensor_weight=0.0;
        int rc0=SnapDetect_run(arena,V,nv,F,nfc,T,&o,&R0);
        o.tensor_weight=256.0;
        int rc1=SnapDetect_run(arena,V,nv,F,nfc,T,&o,&R1);
        int paired=0,same=1,near_choice=1; double d0=0.0,d1=0.0;
        for(size_t i=base_src;i<nv;i++){
            if(R0.vclass[i]!=SNAP_FIXABLE||R1.vclass[i]!=SNAP_FIXABLE)continue;
            paired++; d0+=R0.voff[i]; d1+=R1.voff[i];
            if(fabs(R1.voff[i]-R0.voff[i])>1e-5f)same=0;
            if(R1.voff[i]>4.5f)near_choice=0;
        }
        fprintf(stderr,"  t9 diag: paired=%d same=%d near=%d mean_d=%.2f->%.2f\n",
                paired,same,near_choice,paired?d0/paired:0.0,paired?d1/paired:0.0);
        CK(rc0==0&&rc1==0&&paired>0,"t9 quilting selectors both run");
        CK(same,"t9 retired tensor knob cannot pull repair targets to another ply");
        CK(near_choice,"t9 quilting+distance chooses the nearer coherent slab");
        CK(R1.quilt_cost_mean>=0.0,"t9 quilting cost is reported");
        free(V);free(T);free(F);
        o.axis_point[2]=8.0f; o.reach=8.0; o.min_region=8; o.tensor_weight=0.0;
    }
    /* t4: degenerate */
    { SnapResult R; CK(SnapDetect_run(arena,NULL,0,NULL,0,NULL,&o,&R)!=0,"t4 empty -> error"); }

    free(buf);
    Arena_dispose(&arena);
    fprintf(stderr,"[selftest] %s (%d failures)\n",nf==0?"ALL PASS":"FAILURES",nf);
    return nf;
}
