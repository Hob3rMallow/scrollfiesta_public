/* ============================================================================
 * scroll_snap -- inspect or run the two-pass snap stage.  Detection labels
 * dark whiffs; --solve first repairs those islands with coherent quilting
 * targets, then globally refines the mesh onto the oriented recto edge.
 * ==========================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../common/arena.h"
#include "../common/ves_platform.h"
#include "../common/tiff_io.h"
#include "../common/obj_io.h"
#include "../flatten/surface_snap.h"
#include "../flatten/snap_solve.h"
#include "../flatten/recto_refine.h"

typedef struct { float *v,*t; int32_t *f; size_t nv,nvt,nf; } Obj;
static void *grow(void*p,size_t need,size_t*cap,size_t el){size_t c=*cap;if(need<=c)return p;c=c?c:4096;while(c<need)c*=2;p=realloc(p,c*el);if(!p){fprintf(stderr,"OOM\n");exit(1);}*cap=c;return p;}

static int obj_load(const char*path,Obj*m){
    FILE*fp=fopen(path,"rb");char line[1024];size_t cv=0,ctv=0,cf=0;memset(m,0,sizeof*m);
    if(!fp){fprintf(stderr,"cannot open %s\n",path);return -1;}
    while(fgets(line,sizeof line,fp)){
        if(line[0]=='v'&&line[1]==' '){char*s=line+2,*e;double a=strtod(s,&e);s=e;double b=strtod(s,&e);s=e;double c=strtod(s,&e);
            m->v=(float*)grow(m->v,(m->nv+1)*3,&cv,sizeof(float));m->v[m->nv*3+0]=(float)a;m->v[m->nv*3+1]=(float)b;m->v[m->nv*3+2]=(float)c;m->nv++;}
        else if(line[0]=='v'&&line[1]=='t'&&line[2]==' '){char*s=line+3,*e;double u=strtod(s,&e);s=e;double v=strtod(s,&e);
            m->t=(float*)grow(m->t,(m->nvt+1)*2,&ctv,sizeof(float));m->t[m->nvt*2+0]=(float)u;m->t[m->nvt*2+1]=(float)v;m->nvt++;}
        else if(line[0]=='f'&&(line[1]==' '||line[1]=='\t')){long idx[16];int ni=0;char*s=line+2;
            while(ni<16){char*e;while(*s==' '||*s=='\t')s++;if(!*s||*s=='\n'||*s=='\r')break;long v=strtol(s,&e,10);if(e==s)break;idx[ni++]=v;s=e;while(*s&&*s!=' '&&*s!='\t'&&*s!='\n'&&*s!='\r')s++;}
            for(int k=2;k<ni;k++){m->f=(int32_t*)grow(m->f,(m->nf+1)*3,&cf,sizeof(int32_t));m->f[m->nf*3+0]=(int32_t)(idx[0]-1);m->f[m->nf*3+1]=(int32_t)(idx[k-1]-1);m->f[m->nf*3+2]=(int32_t)(idx[k]-1);m->nf++;}}
    }
    fclose(fp);
    if(m->nv==0||m->nvt!=m->nv){fprintf(stderr,"%s: nv=%zu nvt=%zu (need vt==v)\n",path,m->nv,m->nvt);return -1;}
    return 0;
}

/* colour per vertex: GOOD grey; dark verts by REGION class */
static void vcolor(const SnapResult*R,size_t i,double*r,double*g,double*b){
    if(R->vclass[i]==SNAP_GOOD){*r=*g=*b=0.55;return;}
    int32_t reg=R->vregion[i]; uint8_t rc=(reg>=0)?R->rclass[reg]:SNAPREG_MIXED;
    switch(rc){
        case SNAPREG_FIXABLE:   *r=0.13;*g=0.80;*b=0.25; break;  /* green  */
        case SNAPREG_CRACK:     *r=0.88;*g=0.12;*b=0.12; break;  /* red    */
        case SNAPREG_ANCHORLESS:*r=0.96;*g=0.55;*b=0.10; break;  /* orange */
        default:                *r=0.90;*g=0.85;*b=0.15; break;  /* yellow */
    }
}
static int write_colored(const char*path,const float*v,size_t nv,const int32_t*f,size_t nf,
                         const SnapResult*R,int flat,const float*t){
    FILE*fp=fopen(path,"wb");if(!fp)return -1;size_t i;
    for(i=0;i<nv;i++){double r,g,b;vcolor(R,i,&r,&g,&b);
        if(flat) fprintf(fp,"v %.4f %.4f 0 %.3f %.3f %.3f\n",(double)t[i*2+0],(double)t[i*2+1],r,g,b);
        else     fprintf(fp,"v %.4f %.4f %.4f %.3f %.3f %.3f\n",(double)v[i*3+0],(double)v[i*3+1],(double)v[i*3+2],r,g,b);}
    for(i=0;i<nf;i++) fprintf(fp,"f %d %d %d\n",f[i*3+0]+1,f[i*3+1]+1,f[i*3+2]+1);
    fclose(fp);return 0;
}

/* per-vertex display class for the raster/OBJ: 1 good, 2 fix-island, 3 crack,
 * 4 anchorless, 5 mixed (matches the diag_render "snap" colormap). */
static int vdisp(const SnapResult*R,size_t i){
    if(R->vclass[i]==SNAP_GOOD)return 1;
    int32_t reg=R->vregion[i]; uint8_t rc=(reg>=0)?R->rclass[reg]:SNAPREG_MIXED;
    switch(rc){case SNAPREG_FIXABLE:return 2;case SNAPREG_CRACK:return 3;
               case SNAPREG_ANCHORLESS:return 4;default:return 5;}
}
/* paint priority so rare/important classes win a shared pixel:
 * crack > anchorless > mixed > fixable > good */
static int drank(int d){switch(d){case 3:return 4;case 4:return 3;case 5:return 2;case 2:return 1;default:return 0;}}

/* Rasterize the per-vertex class onto the UV grid (du=dv=1, same W/H as
 * obj_bake_raw so it aligns with the rawtex ribbons), each pixel = its face's
 * max-priority vertex class. Written as a uint8 TIFF for diag_render (snap mode). */
static int write_class_tif(const char*path,const float*uv,size_t nv,
                           const int32_t*f,size_t nf,const SnapResult*R){
    double umax=0,vmax=0; size_t i;
    for(i=0;i<nv;i++){ if((double)uv[i*2]>umax)umax=uv[i*2]; if((double)uv[i*2+1]>vmax)vmax=uv[i*2+1]; }
    size_t W=(size_t)floor(umax+0.5)+1, H=(size_t)floor(vmax+0.5)+1;
    if(W<1||H<1||W*H>((size_t)1<<28)){fprintf(stderr,"class raster %zux%zu unreasonable\n",W,H);return -1;}
    uint8_t*img=(uint8_t*)calloc(W*H,1); if(!img){fprintf(stderr,"OOM\n");return -1;}
    for(size_t fi=0;fi<nf;fi++){
        size_t a=(size_t)f[fi*3],b=(size_t)f[fi*3+1],c=(size_t)f[fi*3+2];
        double ua=uv[a*2],va=uv[a*2+1],ub=uv[b*2],vb=uv[b*2+1],uc=uv[c*2],vc=uv[c*2+1];
        double A2=(ub-ua)*(vc-va)-(vb-va)*(uc-ua); if(A2>-1e-12&&A2<1e-12)continue;
        int da=vdisp(R,a),db=vdisp(R,b),dc=vdisp(R,c);
        int ra=drank(da),rb=drank(db),rc2=drank(dc);
        int cls=(ra>=rb&&ra>=rc2)?da:(rb>=rc2?db:dc);
        double lox=ua<ub?(ua<uc?ua:uc):(ub<uc?ub:uc), hix=ua>ub?(ua>uc?ua:uc):(ub>uc?ub:uc);
        double loy=va<vb?(va<vc?va:vc):(vb<vc?vb:vc), hiy=va>vb?(va>vc?va:vc):(vb>vc?vb:vc);
        long x0=(long)ceil(lox-0.001),x1=(long)floor(hix+0.001);
        long y0=(long)ceil(loy-0.001),y1=(long)floor(hiy+0.001),xx,yy;
        if(x0<0)x0=0; if(y0<0)y0=0; if(x1>=(long)W)x1=(long)W-1; if(y1>=(long)H)y1=(long)H-1;
        for(yy=y0;yy<=y1;yy++)for(xx=x0;xx<=x1;xx++){
            double l1=(((double)xx-ua)*(vc-va)-((double)yy-va)*(uc-ua))/A2;
            double l2=((ub-ua)*((double)yy-va)-(vb-va)*((double)xx-ua))/A2;
            double l0=1.0-l1-l2;
            if(l0<-1e-6||l1<-1e-6||l2<-1e-6)continue;
            img[(size_t)yy*W+(size_t)xx]=(uint8_t)cls;
        }
    }
    int rc=TiffIO_save(path,img,1,(int)H,(int)W); free(img);
    if(rc==0)fprintf(stderr,"  wrote %s (%zux%zu class raster)\n",path,W,H);
    return rc;
}

int main(int argc,char**argv){
    if(argc>=2&&!strcmp(argv[1],"--selftest"))
        return (Snap_selftest()+SnapSolve_selftest()+RectoRefine_selftest())==0?0:1;
    if(argc<4){fprintf(stderr,"usage: %s <uv.obj> <raw_cube_dir> <out_dir> [--id S]\n"
        "   [--axis-point z y x] [--axis-dir z y x] [--reach F] [--min-gain F] [--region-cap N] [--min-region N]\n"
        "   [--occ-thresh F] [--local-r F] [--gain-bias F]\n"
        "   [--tensor-weight F] [--tensor-radius F]  (legacy; ignored by quilting)\n"
        "   [--dark-lambda F] [--dark-sigma F] [--dark-thresh N]  (lambda 0 = hard threshold)\n"
        "   [--solve] [--no-recto] [--recto-iters N] [--recto-range F]\n"
        "       --solve writes pass-1 <id>_repaired.obj, then recto <id>_snapped.obj\n"
        "       %s --selftest\n",argv[0],argv[0]);return 1;}
    const char*in=argv[1],*raw=argv[2],*outdir=argv[3],*id="snap";
    int do_solve=0,do_recto=1,recto_iters=4; double recto_range=3.0;
    SnapOpts o;SnapOpts_default(&o);o.raw_dir=raw;o.verbose=1;
    for(int i=4;i<argc;i++){
        if(!strcmp(argv[i],"--id")&&i+1<argc)id=argv[++i];
        else if(!strcmp(argv[i],"--axis-point")&&i+3<argc){o.axis_point[0]=(float)atof(argv[++i]);o.axis_point[1]=(float)atof(argv[++i]);o.axis_point[2]=(float)atof(argv[++i]);}
        else if(!strcmp(argv[i],"--axis-dir")&&i+3<argc){o.axis_dir[0]=(float)atof(argv[++i]);o.axis_dir[1]=(float)atof(argv[++i]);o.axis_dir[2]=(float)atof(argv[++i]);}
        else if(!strcmp(argv[i],"--reach")&&i+1<argc)o.reach=atof(argv[++i]);
        else if(!strcmp(argv[i],"--min-gain")&&i+1<argc)o.min_gain=atof(argv[++i]);
        else if(!strcmp(argv[i],"--region-cap")&&i+1<argc)o.region_cap=atoi(argv[++i]);
        else if(!strcmp(argv[i],"--min-region")&&i+1<argc)o.min_region=atoi(argv[++i]);
        else if(!strcmp(argv[i],"--occ-thresh")&&i+1<argc)o.occ_thresh=atof(argv[++i]);
        else if(!strcmp(argv[i],"--local-r")&&i+1<argc)o.local_r=atof(argv[++i]);
        else if(!strcmp(argv[i],"--gain-bias")&&i+1<argc)o.gain_bias=atof(argv[++i]);
        else if(!strcmp(argv[i],"--tensor-weight")&&i+1<argc)o.tensor_weight=atof(argv[++i]);
        else if(!strcmp(argv[i],"--tensor-radius")&&i+1<argc)o.tensor_radius=atof(argv[++i]);
        else if(!strcmp(argv[i],"--dark-lambda")&&i+1<argc)o.dark_lambda=atof(argv[++i]);
        else if(!strcmp(argv[i],"--dark-sigma")&&i+1<argc)o.dark_sigma=atof(argv[++i]);
        else if(!strcmp(argv[i],"--dark-thresh")&&i+1<argc)o.dark_thresh=atoi(argv[++i]);
        else if(!strcmp(argv[i],"--solve"))do_solve=1;
        else if(!strcmp(argv[i],"--no-recto"))do_recto=0;
        else if(!strcmp(argv[i],"--recto-iters")&&i+1<argc)recto_iters=atoi(argv[++i]);
        else if(!strcmp(argv[i],"--recto-range")&&i+1<argc)recto_range=atof(argv[++i]);
        else{fprintf(stderr,"unknown %s\n",argv[i]);return 1;}
    }
    fprintf(stderr,"[scroll_snap] in=%s raw=%s out=%s id=%s\n",in,raw,outdir,id);
    double t0=ves_clock_sec(); Obj m; if(obj_load(in,&m)!=0)return 1;
    fprintf(stderr,"  loaded %zu v, %zu f (%.1fs)\n",m.nv,m.nf,ves_clock_sec()-t0);

    Arena_T arena=Arena_new(); SnapResult R; t0=ves_clock_sec();
    if(SnapDetect_run(arena,m.v,m.nv,m.f,m.nf,m.t,&o,&R)!=0){fprintf(stderr,"SnapDetect failed\n");return 1;}
    fprintf(stderr,"  detect: dark=%zu (fixable=%zu [bidir=%zu] crack=%zu [blocked=%zu]) regions=%zu "
            "(fix-island=%zu crack=%zu anchorless=%zu mixed=%zu reject=%zu[-%zu v]) max_region=%zu (%.1fs)\n",
            R.n_dark,R.n_fixable,R.n_bidir,R.n_crack,R.n_blocked,R.nreg,R.n_reg_fixable,R.n_reg_crack,
            R.n_reg_anchorless,R.n_reg_mixed,R.n_reg_reject,R.n_rejected,R.max_region_size,ves_clock_sec()-t0);
    fprintf(stderr,"  white-region diag: %zu dark verts have cv>band (residual over-grow into bright)\n",R.n_overgrow);

    /* histograms */
    { long rsz[7]={0,0,0,0,0,0,0}; /* 1 2-5 6-20 21-100 101-500 501-4000 >4000 */
      for(size_t r=0;r<R.nreg;r++){int32_t s=R.rsize[r];int b=s<=1?0:s<=5?1:s<=20?2:s<=100?3:s<=500?4:s<=4000?5:6;rsz[b]++;}
      fprintf(stderr,"  region-size hist [1 2-5 6-20 21-100 101-500 501-4000 >4000]: ");
      for(int b=0;b<7;b++)fprintf(stderr,"%ld ",rsz[b]); fprintf(stderr,"\n");
      long gz[4]={0,0,0,0}; /* gain <50 50-100 100-150 >150 (fixable only) */
      for(size_t i=0;i<m.nv;i++) if(R.vclass[i]==SNAP_FIXABLE){
          double g=R.vgain[i]; int gb=g<50?0:g<100?1:g<150?2:3; gz[gb]++; }
      fprintf(stderr,"  fixable gain hist [<50 50-100 100-150 >150]: %ld %ld %ld %ld\n",gz[0],gz[1],gz[2],gz[3]);
      fprintf(stderr,"  pass1 quilting mean=%.1f target_dist=%.2f gco_fallback=%zu\n",
              R.quilt_cost_mean,R.target_dist_mean,R.n_quilt_fallback);
      /* chosen free-space target distance over FIXABLE verts (how far the snap
       * would move) -- with occupancy gating a spread out past 3 vox is legit */
      { long df[8]={0,0,0,0,0,0,0,0}; /* |off|: <1 1-2 2-3 3-4 4-5 5-6 6-7 >=7 */
        for(size_t i=0;i<m.nv;i++) if(R.vclass[i]==SNAP_FIXABLE){
            double a=R.voff[i]; int b=a<1?0:a<2?1:a<3?2:a<4?3:a<5?4:a<6?5:a<7?6:7; df[b]++; }
        fprintf(stderr,"  fixable target-dist hist [<1 1-2 2-3 3-4 4-5 5-6 6-7 >=7]: "
                "%ld %ld %ld %ld %ld %ld %ld %ld\n",df[0],df[1],df[2],df[3],df[4],df[5],df[6],df[7]); }
      /* how close OTHER mesh was to each dark vertex's probe (blocker distance):
       * small => wedged tightly between sheets (real crack); ==reach => open void */
      { long bz[5]={0,0,0,0,0}; /* vblock: <1 1-2 2-4 4-8 ==reach(open) */
        for(size_t i=0;i<m.nv;i++){ if(R.vclass[i]==SNAP_GOOD)continue;
            double a=R.vblock[i]; int b=a<1?0:a<2?1:a<4?2:a<7.999?3:4; bz[b]++; }
        fprintf(stderr,"  dark blocker-dist hist [<1 1-2 2-4 4-8 open]: %ld %ld %ld %ld %ld\n",
                bz[0],bz[1],bz[2],bz[3],bz[4]); }
    }

    char path[4096];
    snprintf(path,sizeof path,"%s/%s_snapdetect.obj",outdir,id); ves_ensure_parent_dir(path);
    if(write_colored(path,m.v,m.nv,m.f,m.nf,&R,0,m.t)==0)fprintf(stderr,"  wrote %s\n",path);
    snprintf(path,sizeof path,"%s/%s_snapdetect_flat.obj",outdir,id);
    if(write_colored(path,m.v,m.nv,m.f,m.nf,&R,1,m.t)==0)fprintf(stderr,"  wrote %s\n",path);
    snprintf(path,sizeof path,"%s/%s_snapclass.tif",outdir,id);
    write_class_tif(path,m.t,m.nv,m.f,m.nf,&R);
    snprintf(path,sizeof path,"%s/%s_snapdetect_stats.json",outdir,id);
    {FILE*fp=fopen(path,"wb");if(fp){fprintf(fp,
        "{ \"id\":\"%s\", \"verts\":%zu, \"dark\":%zu, \"fixable\":%zu, \"bidir\":%zu, \"crack\":%zu, \"blocked\":%zu,\n"
        "  \"regions\":%zu, \"reg_fixable\":%zu, \"reg_crack\":%zu, \"reg_anchorless\":%zu, \"reg_mixed\":%zu,\n"
        "  \"max_region\":%zu, \"quilt_cost_mean\":%.4f, \"target_dist_mean\":%.4f, \"quilt_fallback\":%zu,\n"
        "  \"window_lo\":%.1f, \"window_hi\":%.1f }\n",
        id,m.nv,R.n_dark,R.n_fixable,R.n_bidir,R.n_crack,R.n_blocked,R.nreg,R.n_reg_fixable,R.n_reg_crack,
        R.n_reg_anchorless,R.n_reg_mixed,R.max_region_size,R.quilt_cost_mean,
        R.target_dist_mean,R.n_quilt_fallback,R.win_lo,R.win_hi);
        fclose(fp);fprintf(stderr,"  wrote %s\n",path);}}

    /* --solve: pass 1 repairs only fixable islands.  Pass 2 then uses the same
     * RAW table for a global, outward-oriented recto-edge optimization. */
    if(do_solve){
        SnapSolveOpts so; SnapSolveOpts_default(&so); so.verbose=1;
        SnapSolveResult sres; t0=ves_clock_sec();
        if(SnapSolve_run(arena,m.v,m.nv,m.f,m.nf,&R,&so,&sres)!=0){
            fprintf(stderr,"SnapSolve failed\n"); return 1;
        }
        fprintf(stderr,"  pass1 repair: moved=%zu reverted=%zu mean_disp=%.2f max_disp=%.2f (%.1fs)\n",
                sres.n_moved,sres.n_reverted,sres.mean_disp,sres.max_disp,ves_clock_sec()-t0);
        snprintf(path,sizeof path,"%s/%s_repaired.obj",outdir,id);
        if(ObjIO_write_uv(path,m.v,m.nv,m.f,m.nf,m.t)==0)fprintf(stderr,"  wrote %s\n",path);

        RectoRefineStats rr; memset(&rr,0,sizeof rr);
        if(do_recto){
            RectoRefineOpts ro; RectoRefineOpts_default(&ro);
            memcpy(ro.axis_point,o.axis_point,sizeof ro.axis_point);
            memcpy(ro.axis_dir,o.axis_dir,sizeof ro.axis_dir);
            ro.window_lo=R.win_lo; ro.window_hi=R.win_hi;
            ro.occ_thresh=o.occ_thresh; ro.local_r=o.local_r;
            ro.outer_iters=recto_iters>0?recto_iters:1;
            ro.range=recto_range>0.0?recto_range:3.0; ro.verbose=1;
            if(RectoRefine_run(arena,R.raw_table,m.v,m.nv,m.f,m.nf,&ro,&rr)!=0){
                fprintf(stderr,"RectoRefine failed\n"); return 1;
            }
        }
        snprintf(path,sizeof path,"%s/%s_snapped.obj",outdir,id);
        if(ObjIO_write_uv(path,m.v,m.nv,m.f,m.nf,m.t)==0)fprintf(stderr,"  wrote %s\n",path);
        snprintf(path,sizeof path,"%s/%s_snapped_stats.json",outdir,id);
        {FILE*fp=fopen(path,"wb");if(fp){fprintf(fp,
            "{ \"pass1_moved\":%zu, \"pass1_reverted\":%zu, \"pass1_mean_disp\":%.5f,"
            " \"recto_enabled\":%s, \"recto_iterations\":%d, \"recto_supported\":%zu,"
            " \"recto_moved\":%zu, \"recto_reverted\":%zu, \"recto_slope_limited\":%zu,"
            " \"recto_mean_disp\":%.5f, \"recto_max_disp\":%.5f }\n",
            sres.n_moved,sres.n_reverted,sres.mean_disp,do_recto?"true":"false",
            rr.iterations,rr.n_supported,rr.n_moved,rr.n_reverted,
            rr.n_slope_limited,rr.mean_disp,rr.max_disp);fclose(fp);fprintf(stderr,"  wrote %s\n",path);}}
    }

    free(m.v);free(m.t);free(m.f); Arena_dispose(&arena); return 0;
}
