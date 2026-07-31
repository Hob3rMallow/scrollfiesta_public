/* ============================================================================
 * snap_solve.c -- see snap_solve.h. Displace FIXABLE-ISLAND verts onto the
 * detector's bright target via one data-weighted Laplacian solve (SnapCG_solve),
 * with an occupancy revert so no move crosses into another sheet.
 * ==========================================================================*/

#include "snap_solve.h"

#include "../common/csr.h"
#include "../common/kdtree.h"
#include "../common/snap_cg.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static void *xmalloc(size_t n){void*p=malloc(n?n:1);if(!p){fprintf(stderr,"snapsolve OOM\n");exit(1);}return p;}

void SnapSolveOpts_default(SnapSolveOpts *o){
    memset(o,0,sizeof *o);
    o->alpha=0.85; o->w_anchor=1e4; o->max_iter=96; o->tol=1e-4;
    o->occ_thresh=0.9; o->local_r=2.0; o->max_disp=8.0; o->verbose=0;
}

/* is vertex i a member of a FIXABLE-ISLAND region? (the only movable class) */
static int is_island(const SnapResult *det, size_t i){
    int32_t r=det->vregion[i];
    return det->vclass[i]==SNAP_FIXABLE && r>=0 && det->rclass[r]==SNAPREG_FIXABLE;
}

int SnapSolve_run(Arena_T arena, float *verts, size_t nv,
                  const int32_t *faces, size_t nf,
                  const SnapResult *det, const SnapSolveOpts *opts,
                  SnapSolveResult *out){
    memset(out,0,sizeof *out);
    if(nv<3||nf<1||det==NULL||det->vclass==NULL||det->vregion==NULL||
       det->rclass==NULL||det->voff==NULL||det->vdir==NULL) return -1;
    SnapSolveOpts o=*opts;
    double occ2=(o.occ_thresh>0?o.occ_thresh:0.9); occ2*=occ2;
    double local2=(o.local_r>0?o.local_r:2.0); local2*=local2;
    double aw=(o.alpha>0.0&&o.alpha<1.0)?o.alpha/(1.0-o.alpha):1.0;
    double maxd_cap=o.max_disp>0?o.max_disp:8.0;

    CSR_T adj=CSR_from_faces(arena,faces,nf,nv);
    const int32_t *aoff=CSR_offset(adj), *atgt=CSR_target(adj);

    /* occupancy over ORIGINAL positions (KDTree copies its points, so it keeps
     * them even after `verts` is mutated by the solve). */
    KDTree_T occ=KDTree_new(arena,verts,nv);

    float *v0=(float*)xmalloc(nv*3*sizeof(float)); memcpy(v0,verts,nv*3*sizeof(float));
    float *w=(float*)xmalloc(nv*sizeof(float));
    float *target=(float*)xmalloc(nv*3*sizeof(float));

    /* fixable-island verts pulled toward target; every other vertex pinned */
    size_t n_fix=0;
    for(size_t i=0;i<nv;i++){
        if(is_island(det,i)){
            double deg=(double)(aoff[i+1]-aoff[i]); if(deg<1.0)deg=1.0;
            w[i]=(float)(aw*deg);
            target[i*3+0]=(float)(verts[i*3+0]+(double)det->voff[i]*det->vdir[i*3+0]);
            target[i*3+1]=(float)(verts[i*3+1]+(double)det->voff[i]*det->vdir[i*3+1]);
            target[i*3+2]=(float)(verts[i*3+2]+(double)det->voff[i]*det->vdir[i*3+2]);
            n_fix++;
        } else {
            w[i]=(float)o.w_anchor;
            target[i*3+0]=verts[i*3+0]; target[i*3+1]=verts[i*3+1]; target[i*3+2]=verts[i*3+2];
        }
    }
    if(n_fix==0){ if(o.verbose)fprintf(stderr,"[snapsolve] no fixable-island verts\n");
                  free(v0);free(w);free(target); return 0; }

    SnapCG_solve(arena,adj,w,target,verts,nv,o.max_iter,o.tol);   /* verts in place */

    /* per moved vert: cap overshoot + revert any that entered OTHER mesh */
    double sumd=0.0,maxd=0.0; int32_t hit[16];
    for(size_t i=0;i<nv;i++){
        if(!is_island(det,i)) continue;
        double dz=(double)verts[i*3+0]-v0[i*3+0], dy=(double)verts[i*3+1]-v0[i*3+1], dx=(double)verts[i*3+2]-v0[i*3+2];
        double d=sqrt(dz*dz+dy*dy+dx*dx);
        int revert=(d>maxd_cap);
        if(!revert){
            float q[3]={verts[i*3+0],verts[i*3+1],verts[i*3+2]};
            size_t nh=KDTree_ball_query(occ,q,(float)occ2,hit,16),h;
            if(nh>16)nh=16;
            for(h=0;h<nh;h++){ int32_t bi=hit[h];
                double ez=(double)v0[bi*3+0]-v0[i*3+0], ey=(double)v0[bi*3+1]-v0[i*3+1], ex=(double)v0[bi*3+2]-v0[i*3+2];
                int own=(bi==(int32_t)i);
                for(int32_t e=aoff[i];!own&&e<aoff[i+1];e++)own=(atgt[e]==bi);
                if(ez*ez+ey*ey+ex*ex>local2||!own){revert=1;break;} /* OTHER sheet under the new pos */
            }
        }
        if(revert){ verts[i*3+0]=v0[i*3+0]; verts[i*3+1]=v0[i*3+1]; verts[i*3+2]=v0[i*3+2]; out->n_reverted++; }
        else if(d>1e-6){ out->n_moved++; sumd+=d; if(d>maxd)maxd=d; }
    }
    out->mean_disp = out->n_moved>0 ? sumd/(double)out->n_moved : 0.0;
    out->max_disp = maxd;
    if(o.verbose)
        fprintf(stderr,"[snapsolve] fixable-island=%zu moved=%zu reverted=%zu mean_disp=%.2f max_disp=%.2f\n",
                n_fix,out->n_moved,out->n_reverted,out->mean_disp,out->max_disp);

    free(v0);free(w);free(target);
    return 0;
}

/* ============================================================================
 * Self-test
 * ==========================================================================*/
#define CK(c,m) do{ if(!(c)){fprintf(stderr,"  FAIL: %s\n",(m));nf++;} else fprintf(stderr,"  ok: %s\n",(m)); }while(0)

/* append an ny x nx grid patch in the plane z=z0 (verts (z0,row,col)) to growable
 * V/F; returns the base vertex index of this patch. */
static size_t sg_patch(double z0,int ny,int nx,float **V,int32_t **F,
                       size_t *nv,size_t *nfc,size_t *cv,size_t *cf){
    size_t base=*nv; int r,c;
    for(r=0;r<ny;r++)for(c=0;c<nx;c++){
        if(*nv+1>*cv){*cv=(*cv?*cv*2:256);*V=(float*)realloc(*V,*cv*3*sizeof(float));}
        /* tiny deterministic jitter (<0.02 vox): a perfect integer lattice is a
         * KD-tree build degeneracy (ball_query mis-prunes); real meshes are
         * generic floats. Keeps the plane essentially flat for the assertions. */
        unsigned hh=(unsigned)(*nv)*2654435761u;
        double jz=((double)((hh)&255)/255.0-0.5)*0.03;
        double jy=((double)((hh>>8)&255)/255.0-0.5)*0.03;
        double jx=((double)((hh>>16)&255)/255.0-0.5)*0.03;
        (*V)[*nv*3+0]=(float)(z0+jz);(*V)[*nv*3+1]=(float)(r+jy);(*V)[*nv*3+2]=(float)(c+jx);(*nv)++; }
    for(r=0;r<ny-1;r++)for(c=0;c<nx-1;c++){
        int32_t a=(int32_t)(base+(size_t)r*(size_t)nx+c),b=a+1,d=a+nx,e=d+1;
        if(*nfc+2>*cf){*cf=(*cf?*cf*2:256);*F=(int32_t*)realloc(*F,*cf*3*sizeof(int32_t));}
        (*F)[*nfc*3+0]=a;(*F)[*nfc*3+1]=b;(*F)[*nfc*3+2]=d;(*nfc)++;
        (*F)[*nfc*3+0]=b;(*F)[*nfc*3+1]=e;(*F)[*nfc*3+2]=d;(*nfc)++; }
    return base;
}

/* build a SnapResult over nv verts with a single region 0; caller fills vclass/
 * vregion/rclass/voff/vdir. Everything malloc'd; freed by sr_free. */
static void sr_alloc(SnapResult *R,size_t nv,size_t nreg){
    memset(R,0,sizeof *R);
    R->vclass=(uint8_t*)calloc(nv,1);
    R->vregion=(int32_t*)malloc(nv*sizeof(int32_t));
    R->voff=(float*)calloc(nv,sizeof(float));
    R->vdir=(float*)calloc(nv,3*sizeof(float));
    R->rclass=(uint8_t*)calloc(nreg?nreg:1,1);
    R->rsize=(int32_t*)calloc(nreg?nreg:1,sizeof(int32_t));
    R->nreg=nreg;
    for(size_t i=0;i<nv;i++)R->vregion[i]=-1;
}
static void sr_free(SnapResult *R){ free(R->vclass);free(R->vregion);free(R->voff);free(R->vdir);
                                    free(R->rclass);free(R->rsize); }

int SnapSolve_selftest(void){
    int nf=0; Arena_T arena=Arena_new();
    fprintf(stderr,"[selftest] snap_solve\n");
    SnapSolveOpts o; SnapSolveOpts_default(&o); o.verbose=0;

    /* t1: 7x7 patch in plane z=4; central 3x3 island pulled to a FREE-SPACE
     * target at z=2 (nothing meshed there). Island dimples toward z=2; the GOOD
     * boundary stays put; nothing reverts. */
    {
        float *V=NULL;int32_t *F=NULL;size_t nv=0,nfc=0,cv=0,cf=0;
        sg_patch(4.0,7,7,&V,&F,&nv,&nfc,&cv,&cf);
        SnapResult R; sr_alloc(&R,nv,1); R.rclass[0]=SNAPREG_FIXABLE;
        size_t ctr=0;
        for(size_t k=0;k<nv;k++){ int r=(int)(k/7),c=(int)(k%7);
            if(r>=2&&r<=4&&c>=2&&c<=4){ R.vclass[k]=SNAP_FIXABLE; R.vregion[k]=0;
                R.voff[k]=2.0f; R.vdir[k*3+0]=-1.0f; if(r==3&&c==3)ctr=k; }
            else R.vclass[k]=SNAP_GOOD; }
        float z_ctr0=V[ctr*3+0], z_corner0=V[0*3+0];
        SnapSolveResult sr; int rc=SnapSolve_run(arena,V,nv,F,nfc,&R,&o,&sr);
        CK(rc==0,"t1 run ok");
        CK(sr.n_moved>0 && sr.n_reverted==0,"t1 island moved, none reverted (free-space target)");
        CK(V[ctr*3+0] < z_ctr0-1.4f,"t1 island center lands close to target (aggressive pull)");
        CK(fabs(V[0*3+0]-z_corner0)<1e-3f,"t1 GOOD boundary stays pinned");
        sr_free(&R); free(V);free(F);
    }
    /* t2: idempotent -- all-good patch => nothing moves. */
    {
        float *V=NULL;int32_t *F=NULL;size_t nv=0,nfc=0,cv=0,cf=0;
        sg_patch(0.0,6,6,&V,&F,&nv,&nfc,&cv,&cf);
        SnapResult R; sr_alloc(&R,nv,1); R.rclass[0]=SNAPREG_CRACK;
        for(size_t k=0;k<nv;k++)R.vclass[k]=SNAP_GOOD;
        float *V0=(float*)malloc(nv*3*sizeof(float)); memcpy(V0,V,nv*3*sizeof(float));
        SnapSolveResult sr; int rc=SnapSolve_run(arena,V,nv,F,nfc,&R,&o,&sr);
        int same=1; for(size_t k=0;k<nv*3;k++) if(fabs(V[k]-V0[k])>1e-4f)same=0;
        CK(rc==0 && sr.n_moved==0 && same,"t2 all-good -> no movement (idempotent)");
        free(V0); sr_free(&R); free(V);free(F);
    }
    /* t3: occupancy revert -- island patch at z=4 whose target (z=0) lands ON a
     * SEPARATE blocker sheet at z=0 => every moved vert re-enters other mesh and
     * is rolled back (no snapping THROUGH a sheet). */
    {
        float *V=NULL;int32_t *F=NULL;size_t nv=0,nfc=0,cv=0,cf=0;
        size_t bblk=sg_patch(0.0,5,5,&V,&F,&nv,&nfc,&cv,&cf);      /* blocker plane z=0 */
        size_t bisl=sg_patch(4.0,3,3,&V,&F,&nv,&nfc,&cv,&cf);      /* island z=4 (cols 0-2 within blocker) */
        (void)bblk;
        SnapResult R; sr_alloc(&R,nv,1); R.rclass[0]=SNAPREG_FIXABLE;
        for(size_t k=bisl;k<nv;k++){ R.vclass[k]=SNAP_FIXABLE; R.vregion[k]=0;
            R.voff[k]=4.0f; R.vdir[k*3+0]=-1.0f; }
        float *V0=(float*)malloc(nv*3*sizeof(float)); memcpy(V0,V,nv*3*sizeof(float));
        SnapSolveResult sr; int rc=SnapSolve_run(arena,V,nv,F,nfc,&R,&o,&sr);
        int island_put=1; for(size_t k=bisl;k<nv;k++) for(int m=0;m<3;m++)
            if(fabs(V[k*3+m]-V0[k*3+m])>1e-3f)island_put=0;
        CK(rc==0,"t3 run ok");
        CK(sr.n_reverted>0 && sr.n_moved==0 && island_put,
           "t3 target on another sheet -> all moves reverted (no cross-wrap)");
        free(V0); sr_free(&R); free(V);free(F);
    }
    /* t4: degenerate */
    { SnapSolveResult sr; CK(SnapSolve_run(arena,NULL,0,NULL,0,NULL,&o,&sr)!=0,"t4 empty -> error"); }

    Arena_dispose(&arena);
    fprintf(stderr,"[selftest] %s (%d failures)\n",nf==0?"ALL PASS":"FAILURES",nf);
    return nf;
}
