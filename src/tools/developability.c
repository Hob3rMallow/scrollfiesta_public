/* developability.c -- measure how far a triangle mesh is from DEVELOPABLE
 * (zero Gaussian curvature), via discrete angle-defect Gaussian curvature.
 *
 * Why: a papyrus wrap is a developable surface (it unrolls flat -> K=0
 * everywhere). Where the mesh BRIDGES two wraps at a tight fold/turn, the
 * surface acquires cone/saddle points -> concentrated |K|. So |angle defect|
 * is a direct detector for the "bad slices" the overlap/bridge splitters leave
 * behind: a clean sheet reads ~0, a tangle lights up.
 *
 * Per interior vertex v:  K_v = 2*pi - sum(incident triangle angles at v).
 * Per boundary vertex:    K_v =   pi - sum(...)         (flat-boundary datum).
 * A developable patch has K_v ~ 0 at every interior vertex; total absolute
 * curvature  S|K_int| = sum over interior verts of |K_v|  is the global
 * non-developability score (0 = perfectly developable).
 *
 * Outputs stats, and optionally:
 *   --bad  PATH : the non-developable faces (>= 1 interior vert over --thresh)
 *   --keep PATH : the developable remainder           ("detect a bad slice and
 *                                                       just remove it")
 *   --heatmap PATH : per-vertex |K| heatmap OBJ (blue=flat .. red=high), for
 *                    MeshLab/Blender.
 *
 * Usage:
 *   developability <in.obj> [--thresh RAD] [--bad o.obj] [--keep o.obj]
 *                           [--heatmap o.obj] [--heatmax RAD]
 *   developability --selftest
 */
#include "../common/ves_platform.h"
#include "../common/arena.h"
#include "../common/obj_io.h"
#include "../common/dev_metric.h"        /* DevMetric_compute -- the shared K metric */
#include "../flatten/developability.h"   /* Develop_vertex_energy -- Crane lambda */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define DEV_PI 3.14159265358979323846

static double clampd(double v, double lo, double hi){ return v<lo?lo:(v>hi?hi:v); }

/* The angle-defect K metric moved to src/common/dev_metric.c (DevMetric_compute)
 * so the QC tool, the end-of-pipeline gate, and the pipeline share one
 * definition. This file keeps only the CLI: stats reporting, heatmap colouring,
 * and bad/keep extraction. */

static void heat_color(double t, float *rgb){           /* t in [0,1] -> blue..red */
    t=clampd(t,0.0,1.0);
    /* 4-stop: blue(0,0,1)->cyan(0,1,1)->yellow(1,1,0)->red(1,0,0) */
    if(t<1.0/3.0){ double u=t*3.0;            rgb[0]=0.0f;          rgb[1]=(float)u;        rgb[2]=1.0f; }
    else if(t<2.0/3.0){ double u=(t-1.0/3.0)*3.0; rgb[0]=(float)u;  rgb[1]=1.0f;            rgb[2]=(float)(1.0-u); }
    else { double u=(t-2.0/3.0)*3.0;          rgb[0]=1.0f;          rgb[1]=(float)(1.0-u);  rgb[2]=0.0f; }
}

/* write the subset of faces whose bad[f]==want, remapping verts. */
static int write_submesh(const char *path, const float *V, size_t nv,
                         const int32_t *F, size_t nf,
                         const unsigned char *bad, int want)
{
    int32_t *remap=(int32_t*)malloc(nv*sizeof(int32_t));
    for(size_t i=0;i<nv;i++) remap[i]=-1;
    size_t onv=0, onf=0;
    for(size_t f=0; f<nf; f++) if((bad[f]?1:0)==want) onf++;
    float *ov=(float*)malloc((onf?onf*3*3:3)*sizeof(float));
    int32_t *of=(int32_t*)malloc((onf?onf*3:3)*sizeof(int32_t));
    size_t fi=0;
    for(size_t f=0; f<nf; f++){
        if((bad[f]?1:0)!=want) continue;
        for(int e=0;e<3;e++){
            int32_t vi=F[f*3+e];
            if(remap[vi]<0){ remap[vi]=(int32_t)onv; ov[onv*3]=V[(size_t)vi*3]; ov[onv*3+1]=V[(size_t)vi*3+1]; ov[onv*3+2]=V[(size_t)vi*3+2]; onv++; }
            of[fi*3+e]=remap[vi];
        }
        fi++;
    }
    ves_ensure_parent_dir(path);
    int rc=ObjIO_write(path, ov, onv, of, onf);
    free(remap); free(ov); free(of);
    return rc;
}

/* ---- core analysis + reporting; returns S|K_int| ---- */
static double analyze(const float *V, size_t nv, const int32_t *F, size_t nf,
                      double thresh, int verbose,
                      size_t *out_nbad_faces)
{
    double *defect=(double*)malloc((nv?nv:1)*sizeof(double));
    unsigned char *boundary=(unsigned char*)malloc(nv?nv:1);
    unsigned char *used=(unsigned char*)malloc(nv?nv:1);
    DevMetric_compute(V,nv,F,nf,defect,boundary,used);

    size_t n_int=0,n_bnd=0,n_over=0; double sumabs=0,maxabs=0,signed_sum=0;
    /* histogram of interior |defect|, bins of 0.1 rad up to 1.0+ */
    long hist[12]; memset(hist,0,sizeof hist);
    for(size_t i=0;i<nv;i++){
        if(!used[i]) continue;
        double ad=fabs(defect[i]);
        if(boundary[i]){ n_bnd++; continue; }
        n_int++; sumabs+=ad; signed_sum+=defect[i];
        if(ad>maxabs) maxabs=ad;
        if(ad>thresh) n_over++;
        int b=(int)(ad/0.1); if(b>11)b=11; hist[b]++;
    }
    /* face badness: any interior vertex over thresh */
    size_t nbad=0;
    if(out_nbad_faces){
        for(size_t f=0; f<nf; f++){
            int bad=0;
            for(int e=0;e<3;e++){ int32_t vi=F[f*3+e];
                if(!boundary[vi] && used[vi] && fabs(defect[vi])>thresh){ bad=1; break; } }
            if(bad) nbad++;
        }
        *out_nbad_faces=nbad;
    }
    if(verbose){
        printf("  verts=%zu  faces=%zu  interior=%zu  boundary=%zu\n", nv,nf,n_int,n_bnd);
        printf("  S|K_int| (total abs Gaussian curvature) = %.2f rad\n", sumabs);
        printf("  mean|K_int| = %.4f rad   max|K| = %.3f rad   signed_sum = %.2f\n",
               n_int? sumabs/(double)n_int : 0.0, maxabs, signed_sum);
        printf("  interior verts over thresh %.2f rad: %zu (%.2f%%)\n",
               thresh, n_over, n_int? 100.0*(double)n_over/(double)n_int : 0.0);
        if(out_nbad_faces)
            printf("  non-developable faces: %zu (%.2f%%)\n",
                   nbad, nf? 100.0*(double)nbad/(double)nf : 0.0);
        printf("  |K_int| histogram (0.1 rad bins): ");
        for(int b=0;b<12;b++) printf("%ld%s", hist[b], b<11?",":(b==11?"+\n":"\n"));
    }
    free(defect); free(boundary); free(used);
    return sumabs;
}

static int run_selftest(void){
    int fail=0;
    /* 5x5 flat grid in z=0 plane: every interior vertex defect ~ 0 */
    int N=5; size_t nv=(size_t)N*N;
    float *V=(float*)malloc(nv*3*sizeof(float));
    for(int j=0;j<N;j++) for(int i=0;i<N;i++){ size_t idx=(size_t)j*N+i; V[idx*3+0]=0.0f; V[idx*3+1]=(float)j; V[idx*3+2]=(float)i; }
    size_t nf=(size_t)(N-1)*(N-1)*2;
    int32_t *F=(int32_t*)malloc(nf*3*sizeof(int32_t));
    size_t k=0;
    for(int j=0;j<N-1;j++) for(int i=0;i<N-1;i++){
        int32_t a=(int32_t)(j*N+i), b=(int32_t)(j*N+i+1), c=(int32_t)((j+1)*N+i), d=(int32_t)((j+1)*N+i+1);
        F[k*3]=a;F[k*3+1]=b;F[k*3+2]=d;k++; F[k*3]=a;F[k*3+1]=d;F[k*3+2]=c;k++;
    }
    double *defect=(double*)malloc(nv*sizeof(double));
    unsigned char *bnd=(unsigned char*)malloc(nv), *used=(unsigned char*)malloc(nv);
    DevMetric_compute(V,nv,F,nf,defect,bnd,used);
    size_t ctr=(size_t)2*N+2;   /* center vertex (2,2) */
    double flat_def=fabs(defect[ctr]);
    int f1=(flat_def>1e-3); fail|=f1;
    fprintf(stderr,"[selftest] flat-grid center |K|=%.2e -> %s\n", flat_def, f1?"FAIL":"ok");

    /* raise center vertex -> cone point, |K| must become clearly nonzero */
    V[ctr*3+0]=0.8f;
    DevMetric_compute(V,nv,F,nf,defect,bnd,used);
    double cone_def=fabs(defect[ctr]);
    int f2=(cone_def<0.1);  /* must detect curvature */
    fail|=f2;
    fprintf(stderr,"[selftest] bumped center |K|=%.3f rad -> %s\n", cone_def, f2?"FAIL":"ok");

    /* the bumped grid's total abs curvature must exceed the flat grid's */
    size_t nb=0; double s_bump=analyze(V,nv,F,nf,0.3,0,&nb);
    V[ctr*3+0]=0.0f; size_t nb0=0; double s_flat=analyze(V,nv,F,nf,0.3,0,&nb0);
    int f3=!(s_bump>s_flat+0.3); fail|=f3;
    fprintf(stderr,"[selftest] S|K| bumped %.3f > flat %.3f, bad-faces bump=%zu flat=%zu -> %s\n",
            s_bump, s_flat, nb, nb0, f3?"FAIL":"ok");

    /* lambda (Crane) sanity: the flat grid is developable AND ruled -> lambda~0
     * everywhere (exercises Develop_vertex_energy, the dev-cut detector). V's
     * centre was reset to flat above. */
    {
        Arena_T a2=Arena_new();
        double *lam=(double*)malloc(nv*sizeof(double));
        Develop_vertex_energy(a2, V,nv,F,nf, lam);
        double lmax=0; for(size_t i=0;i<nv;i++) if(lam[i]>lmax) lmax=lam[i];
        int f4=(lmax>1e-6); fail|=f4;
        fprintf(stderr,"[selftest] flat-grid max_lambda=%.2e -> %s\n", lmax, f4?"FAIL":"ok");
        free(lam); Arena_dispose(&a2);
    }

    free(V);free(F);free(defect);free(bnd);free(used);
    fprintf(stderr,"=== developability selftest %s ===\n", fail?"FAILED":"PASSED");
    return fail?1:0;
}

int main(int argc, char **argv){
    if(argc>=2 && !strcmp(argv[1],"--selftest")) return run_selftest();
    if(argc<2){
        fprintf(stderr,"Usage: %s <in.obj> [--lambda] [--thresh V] [--bad o.obj] [--keep o.obj] [--heatmap o.obj] [--heatmax V]\n"
                       "       --lambda: colour Crane covariance energy (dev-cut metric) instead of angle-defect K\n"
                       "       %s --selftest\n", argv[0], argv[0]);
        return 1;
    }
    const char *in=argv[1], *bad_path=NULL, *keep_path=NULL, *heat_path=NULL;
    double thresh=0.5, heatmax=1.0;
    int use_lambda=0, thresh_set=0, heatmax_set=0;
    for(int i=2;i<argc;i++){
        if(!strcmp(argv[i],"--thresh")&&i+1<argc){ thresh=atof(argv[++i]); thresh_set=1; }
        else if(!strcmp(argv[i],"--bad")&&i+1<argc) bad_path=argv[++i];
        else if(!strcmp(argv[i],"--keep")&&i+1<argc) keep_path=argv[++i];
        else if(!strcmp(argv[i],"--heatmap")&&i+1<argc) heat_path=argv[++i];
        else if(!strcmp(argv[i],"--heatmax")&&i+1<argc){ heatmax=atof(argv[++i]); heatmax_set=1; }
        else if(!strcmp(argv[i],"--lambda")) use_lambda=1;
        else { fprintf(stderr,"unknown arg %s\n",argv[i]); return 1; }
    }
    /* Crane lambda lives on a much smaller scale than angle-defect K (rad), so
     * give --lambda its own defaults: threshold = the dev-cut DEV_CUT_EPS, and a
     * heatmap full-scale a few x that so real seams read warm. */
    if(use_lambda){
        if(!thresh_set)  thresh  = 0.05;   /* DEV_CUT_EPS */
        if(!heatmax_set) heatmax = 0.15;
    }

    Arena_T arena=Arena_new();
    float *V=NULL; int32_t *F=NULL; size_t nv=0,nf=0;
    if(ObjIO_read(arena,in,&V,&nv,&F,&nf)!=0){
        fprintf(stderr,"developability: cannot read %s\n",in); Arena_dispose(&arena); return 1;
    }
    printf("== %s ==\n", in);

    if(use_lambda){
        /* Crane covariance energy lambda_min(A_v), A_v = sum_corner theta*N N^T --
         * the metric the dev-cut splitter thresholds. Unlike K (angle defect) it
         * reads ~0 on a cone/cylinder (developable AND ruled) but lights up on a
         * corrugated or doubly-curved (merged) junction. */
        double *lam =(double*)malloc((nv?nv:1)*sizeof(double));
        double *kdef=(double*)malloc((nv?nv:1)*sizeof(double));
        unsigned char *bnd=(unsigned char*)malloc(nv?nv:1), *used=(unsigned char*)malloc(nv?nv:1);
        Develop_vertex_energy(arena, V,nv,F,nf, lam);
        DevMetric_compute(V,nv,F,nf, kdef, bnd, used);   /* for the boundary/used mask */
        size_t n_int=0,n_over=0; double sum=0,mx=0;
        long hist[12]; memset(hist,0,sizeof hist);
        double binw = heatmax>0.0 ? heatmax/11.0 : 1.0;
        for(size_t i=0;i<nv;i++){
            if(!used[i]||bnd[i]) continue;
            double l=lam[i]; n_int++; sum+=l; if(l>mx)mx=l; if(l>thresh)n_over++;
            int b=(int)(l/binw); if(b<0)b=0; if(b>11)b=11; hist[b]++;
        }
        printf("  [lambda] verts=%zu faces=%zu interior=%zu\n", nv,nf,n_int);
        printf("  mean_lambda=%.5f  max_lambda=%.5f\n", n_int?sum/(double)n_int:0.0, mx);
        printf("  interior over eps %.3f: %zu (%.2f%%)\n",
               thresh, n_over, n_int?100.0*(double)n_over/(double)n_int:0.0);
        printf("  lambda histogram (%.4f bins to %.3f): ", binw, heatmax);
        for(int b=0;b<12;b++) printf("%ld%s", hist[b], b<11?",":"+\n");
        if(heat_path){
            float *colors=(float*)malloc((nv?nv:1)*3*sizeof(float));
            for(size_t i=0;i<nv;i++) heat_color(lam[i]/heatmax, &colors[i*3]);
            ves_ensure_parent_dir(heat_path);
            int rc=ObjIO_write_per_vertex_color(heat_path, V,nv,F,nf,colors);
            printf("  wrote lambda heatmap -> %s (%s)\n", heat_path, rc?"FAIL":"ok");
            free(colors);
        }
        if(bad_path||keep_path){
            unsigned char *badf=(unsigned char*)calloc(nf?nf:1,1);
            for(size_t f=0;f<nf;f++) for(int e=0;e<3;e++){ int32_t vi=F[f*3+e];
                if(!bnd[vi]&&used[vi]&&lam[vi]>thresh){ badf[f]=1; break; } }
            if(bad_path){  int rc=write_submesh(bad_path, V,nv,F,nf,badf,1);  printf("  wrote bad  -> %s (%s)\n", bad_path, rc?"FAIL":"ok"); }
            if(keep_path){ int rc=write_submesh(keep_path,V,nv,F,nf,badf,0); printf("  wrote keep -> %s (%s)\n", keep_path, rc?"FAIL":"ok"); }
            free(badf);
        }
        free(lam);free(kdef);free(bnd);free(used);
        Arena_dispose(&arena);
        return 0;
    }

    size_t nbad=0;
    analyze(V,nv,F,nf,thresh,1,&nbad);

    if(bad_path || keep_path || heat_path){
        double *defect=(double*)malloc((nv?nv:1)*sizeof(double));
        unsigned char *boundary=(unsigned char*)malloc(nv?nv:1), *used=(unsigned char*)malloc(nv?nv:1);
        DevMetric_compute(V,nv,F,nf,defect,boundary,used);
        if(bad_path || keep_path){
            unsigned char *badf=(unsigned char*)calloc(nf?nf:1,1);
            for(size_t f=0; f<nf; f++) for(int e=0;e<3;e++){ int32_t vi=F[f*3+e];
                if(!boundary[vi] && used[vi] && fabs(defect[vi])>thresh){ badf[f]=1; break; } }
            if(bad_path){  int rc=write_submesh(bad_path, V,nv,F,nf,badf,1);  printf("  wrote bad  -> %s (%s)\n", bad_path, rc?"FAIL":"ok"); }
            if(keep_path){ int rc=write_submesh(keep_path,V,nv,F,nf,badf,0); printf("  wrote keep -> %s (%s)\n", keep_path, rc?"FAIL":"ok"); }
            free(badf);
        }
        if(heat_path){
            float *colors=(float*)malloc((nv?nv:1)*3*sizeof(float));
            for(size_t i=0;i<nv;i++) heat_color(fabs(defect[i])/heatmax, &colors[i*3]);
            ves_ensure_parent_dir(heat_path);
            int rc=ObjIO_write_per_vertex_color(heat_path, V,nv,F,nf,colors);
            printf("  wrote heatmap -> %s (%s)\n", heat_path, rc?"FAIL":"ok");
            free(colors);
        }
        free(defect); free(boundary); free(used);
    }
    Arena_dispose(&arena);
    return 0;
}
