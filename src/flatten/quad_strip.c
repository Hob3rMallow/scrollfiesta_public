/* ============================================================================
 * quad_strip.c -- see quad_strip.h.  One quad strip per scroll component.
 * Production keeps the component rectangle coverage; the conservative bounded
 * topology remains available for diagnostics and deliberately sparse exports.
 * ==========================================================================*/

#include "quad_strip.h"

#include "../common/csr.h"
#include "../common/snap_cg.h"
#include "../common/eig3.h"
#include "../common/ves_platform.h"
#include "../common/pipeline_constants.h"
#ifndef QUAD_STRIP_NO_TAUCS
#include "sparse_solve.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#define QUAD_PI     3.1415926535897932384626433832795
#define QUAD_TWO_PI 6.283185307179586476925286766559

/* Harmonic-fill solver: 1 = TAUCS supernodal Cholesky, factored once per
 * component and solved for every channel (the CLAUDE.md solver doctrine);
 * 0 = the legacy Jacobi-PCG, kept for A/B and for executables that do not
 * link TAUCS (QUAD_STRIP_NO_TAUCS forces 0). */
#ifdef QUAD_STRIP_NO_TAUCS
#define QUAD_STRIP_HARMONIC_SOLVER_DEFAULT 0
#endif
#ifndef QUAD_STRIP_HARMONIC_SOLVER_DEFAULT
#define QUAD_STRIP_HARMONIC_SOLVER_DEFAULT 1
#endif

#ifdef _OPENMP
#endif

void QuadStrip_defaults(QuadStripOpts *o) {
    memset(o, 0, sizeof *o);
    o->mode = QUAD_STRIP_RECT;
    o->max_hole_distance = 16;
    o->max_row_gap = 1;
    o->fit_stiffness = 1e4;
    o->exact_fitted = 1;
    o->pde_max_iter = 1024;
    o->pde_tol = 1e-6;
    o->cylindrical_fill = 0;
    o->axis_y = 0.0;
    o->axis_x = 0.0;
    o->relax_rounds = 0;
    o->relax_cg_iters = 96;
    o->trend_per_row = 0;
    o->harmonic_solver = QUAD_STRIP_HARMONIC_SOLVER_DEFAULT;
}

void QuadStripArap_defaults(QuadStripArapOpts *o) {
    memset(o, 0, sizeof *o);
    o->max_iterations = 1000;
    o->movement_tolerance = 1e-5;
    o->movement_tolerance_rel = 1e-3;
    o->movement_stall_window = 50;
    o->movement_stall_fraction = 0.01;
    o->cg_max_iterations = 256;
    o->cg_relative_tolerance = 1e-8;
    o->multigrid_levels = 16;
    o->multigrid_cycles = 1;
    o->multigrid_pre_sweeps = 2;
    o->multigrid_post_sweeps = 2;
    o->multigrid_coarse_sweeps = 32;
    o->nonlinear_multigrid_levels = 4;
    o->nonlinear_coarse_iterations = 128;
    o->nonlinear_min_vertices = 4000000;
    o->source_weight = 1.0;
    o->fill_weight = 0.01;
    o->huber_delta = 2.0;
    o->minimum_source_scale = 0.05;
    o->minimum_area_ratio = 0.05;
    o->maximum_vertex_step = 2.0;
    o->preserve_axial = 0;
    o->fixed_vertices = NULL;
    o->repair_weight = NULL;
    o->repair_source_scale = 0.0001;
    o->repair_frame_blend = 0.85;
    o->repair_frame_screen = 0.002;
    o->repair_frame_sweeps = 128;
    o->axis_y = 0.0; o->axis_x = 0.0;
    o->axis_normal_weight = 0.0;
    o->axis_normal_radius = 0.0;
    o->axis_fold_limit = 0.0;
    o->verbose = 0;
}

/* ------------------------------------------------------------------------- */
/* Full structured metric ARAP                                               */

/* The legacy fill-only ARAP below estimates a rank-deficient 3x3 polar
 * factor.  A surface has a 3x2 differential, so the full-ribbon solve uses the
 * actual closest orthonormal 3x2 frame.  S is stored row-major (three world
 * channels by the two flat-rest channels). */
static int metric_closest_frame(const double S[6], const double previous[6],
                                double F[6]) {
    double a=0.0,b=0.0,c=0.0;
    for(int k=0;k<3;k++){
        a+=S[k*2]*S[k*2];
        b+=S[k*2]*S[k*2+1];
        c+=S[k*2+1]*S[k*2+1];
    }
    double tr=a+c,disc=hypot(a-c,2.0*b);
    double l0=0.5*(tr+disc),l1=0.5*(tr-disc);
    if(l0>1e-24&&l1>1e-24){
        double angle=0.5*atan2(2.0*b,a-c);
        double q00=cos(angle),q10=sin(angle);
        double q01=-q10,q11=q00;
        double d0=1.0/sqrt(l0),d1=1.0/sqrt(l1);
        double m00=d0*q00*q00+d1*q01*q01;
        double m01=d0*q00*q10+d1*q01*q11;
        double m11=d0*q10*q10+d1*q11*q11;
        for(int k=0;k<3;k++){
            F[k*2]=S[k*2]*m00+S[k*2+1]*m01;
            F[k*2+1]=S[k*2]*m01+S[k*2+1]*m11;
        }
        return 0;
    }

    /* A one-row/one-column or temporarily collapsed neighbourhood has rank
     * one.  Complete its reliable tangent against the preceding frame rather
     * than inventing a discontinuous normal. */
    double u[3]={S[0],S[2],S[4]},v[3]={S[1],S[3],S[5]};
    double un=sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);
    if(un<1e-12){u[0]=previous[0];u[1]=previous[2];u[2]=previous[4];
        un=sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);}
    if(un<1e-12)return -1;
    for(int k=0;k<3;k++)u[k]/=un;
    double dot=v[0]*u[0]+v[1]*u[1]+v[2]*u[2];
    for(int k=0;k<3;k++)v[k]-=dot*u[k];
    double vn=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
    if(vn<1e-12){
        v[0]=previous[1];v[1]=previous[3];v[2]=previous[5];
        dot=v[0]*u[0]+v[1]*u[1]+v[2]*u[2];
        for(int k=0;k<3;k++)v[k]-=dot*u[k];
        vn=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
    }
    if(vn<1e-12)return -1;
    for(int k=0;k<3;k++){F[k*2]=u[k];F[k*2+1]=v[k]/vn;}
    return 0;
}

static double metric_stencil_sum(int r,int c,int H,int W,
                                 double wh,double wv,double wd){
    double s=0.0;
    if(c>0)s+=wh;if(c+1<W)s+=wh;
    if(r>0)s+=wv;if(r+1<H)s+=wv;
    if(r>0&&c>0)s+=wd;if(r+1<H&&c+1<W)s+=wd;
    return s;
}

static void metric_matvec(const double *x,double *y,const double *diag,
                          int H,int W,double wh,double wv,double wd){
    size_t n=(size_t)H*(size_t)W;
    int kk;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for(kk=0;kk<(int)n;kk++){
        size_t k=(size_t)kk;int r=(int)(k/(size_t)W),c=(int)(k%(size_t)W);
        double s=diag[k]*x[k];
        if(c>0)s-=wh*x[k-1];if(c+1<W)s-=wh*x[k+1];
        if(r>0)s-=wv*x[k-(size_t)W];if(r+1<H)s-=wv*x[k+(size_t)W];
        if(r>0&&c>0)s-=wd*x[k-(size_t)W-1];
        if(r+1<H&&c+1<W)s-=wd*x[k+(size_t)W+1];
        y[k]=s;
    }
}

/* ------------------------------------------------------------------------- */
/* Structured Galerkin multigrid for the metric-ARAP global step.            */

/* Piecewise-constant 2x2 aggregation is a good match for these very long
 * regular ribbons.  P is constant within an aggregate and restriction is
 * P^T, so each stored coarse operator is the exact Galerkin P^T A P of the
 * preceding level.  The NW-SE triangle diagonal is retained: four colors
 * make every smoother set independent, and reversing the color order in the
 * post-sweep gives a symmetric V-cycle.  The fine level stays implicit to
 * avoid duplicating another whole-ribbon operator. */
#define METRIC_MG_MAX_LEVELS 16

typedef struct MetricMGLevel {
    int H,W;
    size_t n;
    float *diag,*right,*down,*diag_dr;
    float *x,*rhs,*residual;
} MetricMGLevel;

typedef struct MetricMG {
    int nlevel; /* stored coarse levels; the fine level is implicit */
    int pre_sweeps,post_sweeps,coarse_sweeps;
    MetricMGLevel level[METRIC_MG_MAX_LEVELS-1];
} MetricMG;

static int metric_mg_level_alloc(MetricMGLevel *L,int H,int W){
    memset(L,0,sizeof *L);
    if(H<1||W<1||(size_t)H>SIZE_MAX/(size_t)W)return -1;
    size_t n=(size_t)H*(size_t)W;
    if(n>SIZE_MAX/sizeof(float))return -1;
    L->H=H;L->W=W;L->n=n;
    L->diag=(float *)calloc(n,sizeof *L->diag);
    L->right=(float *)calloc(n,sizeof *L->right);
    L->down=(float *)calloc(n,sizeof *L->down);
    L->diag_dr=(float *)calloc(n,sizeof *L->diag_dr);
    L->x=(float *)calloc(n,sizeof *L->x);
    L->rhs=(float *)calloc(n,sizeof *L->rhs);
    L->residual=(float *)calloc(n,sizeof *L->residual);
    if(!L->diag||!L->right||!L->down||!L->diag_dr||!L->x||!L->rhs||!L->residual){
        free(L->diag);free(L->right);free(L->down);free(L->diag_dr);
        free(L->x);free(L->rhs);free(L->residual);memset(L,0,sizeof *L);return -1;
    }
    return 0;
}

static void metric_mg_level_free(MetricMGLevel *L){
    free(L->diag);free(L->right);free(L->down);free(L->diag_dr);
    free(L->x);free(L->rhs);free(L->residual);memset(L,0,sizeof *L);
}

static void metric_mg_free(MetricMG *M){
    if(!M)return;
    for(int l=0;l<M->nlevel;l++)metric_mg_level_free(&M->level[l]);
    memset(M,0,sizeof *M);
}

static int metric_mg_init(MetricMG *M,int fine_h,int fine_w,
                          int requested_levels,int pre,int post,int coarse){
    memset(M,0,sizeof *M);M->pre_sweeps=pre;M->post_sweeps=post;M->coarse_sweeps=coarse;
    if(requested_levels<=1)return 0;
    int H=fine_h,W=fine_w;
    size_t n=(size_t)H*(size_t)W;
    while(M->nlevel<requested_levels-1&&M->nlevel<METRIC_MG_MAX_LEVELS-1&&n>32){
        H=(H+1)/2;W=(W+1)/2;n=(size_t)H*(size_t)W;
        if(metric_mg_level_alloc(&M->level[M->nlevel],H,W)!=0){metric_mg_free(M);return -1;}
        M->nlevel++;
        if(H==1&&W==1)break;
    }
    return 0;
}

/* Accumulate a preceding-level edge into a 2x2 aggregate operator.  The
 * preceding diagonal already contains both endpoint conductances; an edge
 * internal to one aggregate therefore cancels twice. */
static int metric_mg_add_edge(MetricMGLevel *coarse,
                              int ar,int ac,int br,int bc,double weight){
    int par=ar>>1,pac=ac>>1,pbr=br>>1,pbc=bc>>1;
    size_t a=(size_t)par*(size_t)coarse->W+(size_t)pac;
    size_t b=(size_t)pbr*(size_t)coarse->W+(size_t)pbc;
    if(a==b){coarse->diag[a]=(float)((double)coarse->diag[a]-2.0*weight);return 0;}
    if(par==pbr&&pbc==pac+1){coarse->right[a]=(float)((double)coarse->right[a]+weight);return 0;}
    if(pac==pbc&&pbr==par+1){coarse->down[a]=(float)((double)coarse->down[a]+weight);return 0;}
    if(pbr==par+1&&pbc==pac+1){coarse->diag_dr[a]=(float)((double)coarse->diag_dr[a]+weight);return 0;}
    return -1;
}

static int metric_mg_build_first(MetricMGLevel *coarse,const double *fine_diag,
                                 int H,int W,double wh,double wv,double wd){
    memset(coarse->diag,0,coarse->n*sizeof *coarse->diag);
    memset(coarse->right,0,coarse->n*sizeof *coarse->right);
    memset(coarse->down,0,coarse->n*sizeof *coarse->down);
    memset(coarse->diag_dr,0,coarse->n*sizeof *coarse->diag_dr);
    for(int r=0;r<H;r++)for(int c=0;c<W;c++){
        size_t v=(size_t)r*(size_t)W+(size_t)c;
        size_t p=(size_t)(r>>1)*(size_t)coarse->W+(size_t)(c>>1);
        coarse->diag[p]=(float)((double)coarse->diag[p]+fine_diag[v]);
    }
    for(int r=0;r<H;r++)for(int c=0;c<W;c++){
        if(c+1<W&&metric_mg_add_edge(coarse,r,c,r,c+1,wh)!=0)return -1;
        if(r+1<H&&metric_mg_add_edge(coarse,r,c,r+1,c,wv)!=0)return -1;
        if(r+1<H&&c+1<W&&metric_mg_add_edge(coarse,r,c,r+1,c+1,wd)!=0)return -1;
    }
    return 0;
}

static int metric_mg_coarsen(const MetricMGLevel *fine,MetricMGLevel *coarse){
    memset(coarse->diag,0,coarse->n*sizeof *coarse->diag);
    memset(coarse->right,0,coarse->n*sizeof *coarse->right);
    memset(coarse->down,0,coarse->n*sizeof *coarse->down);
    memset(coarse->diag_dr,0,coarse->n*sizeof *coarse->diag_dr);
    for(int r=0;r<fine->H;r++)for(int c=0;c<fine->W;c++){
        size_t v=(size_t)r*(size_t)fine->W+(size_t)c;
        size_t p=(size_t)(r>>1)*(size_t)coarse->W+(size_t)(c>>1);
        coarse->diag[p]=(float)((double)coarse->diag[p]+fine->diag[v]);
    }
    for(int r=0;r<fine->H;r++)for(int c=0;c<fine->W;c++){
        size_t v=(size_t)r*(size_t)fine->W+(size_t)c;
        if(c+1<fine->W&&fine->right[v]>0.0f&&
           metric_mg_add_edge(coarse,r,c,r,c+1,fine->right[v])!=0)return -1;
        if(r+1<fine->H&&fine->down[v]>0.0f&&
           metric_mg_add_edge(coarse,r,c,r+1,c,fine->down[v])!=0)return -1;
        if(r+1<fine->H&&c+1<fine->W&&fine->diag_dr[v]>0.0f&&
           metric_mg_add_edge(coarse,r,c,r+1,c+1,fine->diag_dr[v])!=0)return -1;
    }
    return 0;
}

static int metric_mg_rebuild(MetricMG *M,const double *fine_diag,
                             int H,int W,double wh,double wv,double wd){
    if(M->nlevel==0)return 0;
    if(metric_mg_build_first(&M->level[0],fine_diag,H,W,wh,wv,wd)!=0)return -1;
    for(int l=1;l<M->nlevel;l++)if(metric_mg_coarsen(&M->level[l-1],&M->level[l])!=0)return -1;
    return 0;
}

static void metric_mg_level_smooth(MetricMGLevel *L,int sweeps,int reverse){
    for(int sweep=0;sweep<sweeps;sweep++)for(int order=0;order<4;order++){
        int color=reverse?3-order:order;
        int row_parity=color>>1,column_parity=color&1;
        int r;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for(r=row_parity;r<L->H;r+=2)for(int c=column_parity;c<L->W;c+=2){
            size_t v=(size_t)r*(size_t)L->W+(size_t)c;
            if(L->diag[v]<=1e-12f)continue;
            double num=L->rhs[v];
            if(c>0)num+=(double)L->right[v-1]*L->x[v-1];
            if(c+1<L->W)num+=(double)L->right[v]*L->x[v+1];
            if(r>0)num+=(double)L->down[v-(size_t)L->W]*L->x[v-(size_t)L->W];
            if(r+1<L->H)num+=(double)L->down[v]*L->x[v+(size_t)L->W];
            if(r>0&&c>0)num+=(double)L->diag_dr[v-(size_t)L->W-1]*L->x[v-(size_t)L->W-1];
            if(r+1<L->H&&c+1<L->W)num+=(double)L->diag_dr[v]*L->x[v+(size_t)L->W+1];
            L->x[v]=(float)(num/(double)L->diag[v]);
        }
    }
}

static void metric_mg_level_residual(MetricMGLevel *L){
    int kk;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for(kk=0;kk<(int)L->n;kk++){
        size_t v=(size_t)kk;int r=(int)(v/(size_t)L->W),c=(int)(v%(size_t)L->W);
        double ax=(double)L->diag[v]*L->x[v];
        if(c>0)ax-=(double)L->right[v-1]*L->x[v-1];
        if(c+1<L->W)ax-=(double)L->right[v]*L->x[v+1];
        if(r>0)ax-=(double)L->down[v-(size_t)L->W]*L->x[v-(size_t)L->W];
        if(r+1<L->H)ax-=(double)L->down[v]*L->x[v+(size_t)L->W];
        if(r>0&&c>0)ax-=(double)L->diag_dr[v-(size_t)L->W-1]*L->x[v-(size_t)L->W-1];
        if(r+1<L->H&&c+1<L->W)ax-=(double)L->diag_dr[v]*L->x[v+(size_t)L->W+1];
        L->residual[v]=(float)((double)L->rhs[v]-ax);
    }
}

static void metric_mg_restrict_float(const MetricMGLevel *fine,MetricMGLevel *coarse){
    int kk;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for(kk=0;kk<(int)coarse->n;kk++){
        size_t p=(size_t)kk;int pr=(int)(p/(size_t)coarse->W),pc=(int)(p%(size_t)coarse->W);
        double sum=0.0;
        for(int dr=0;dr<2;dr++)for(int dc=0;dc<2;dc++){
            int r=2*pr+dr,c=2*pc+dc;
            if(r<fine->H&&c<fine->W)sum+=fine->residual[(size_t)r*(size_t)fine->W+(size_t)c];
        }
        coarse->rhs[p]=(float)sum;coarse->x[p]=0.0f;
    }
}

static void metric_mg_prolong_float(MetricMGLevel *fine,const MetricMGLevel *coarse){
    int kk;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for(kk=0;kk<(int)fine->n;kk++){
        size_t v=(size_t)kk;int r=(int)(v/(size_t)fine->W),c=(int)(v%(size_t)fine->W);
        size_t p=(size_t)(r>>1)*(size_t)coarse->W+(size_t)(c>>1);
        fine->x[v]+=coarse->x[p];
    }
}

static void metric_mg_vcycle(MetricMG *M,int level){
    MetricMGLevel *L=&M->level[level];
    metric_mg_level_smooth(L,M->pre_sweeps,0);
    if(level+1==M->nlevel){
        metric_mg_level_smooth(L,M->coarse_sweeps,0);
        metric_mg_level_smooth(L,M->coarse_sweeps,1);
        return;
    }
    metric_mg_level_residual(L);
    metric_mg_restrict_float(L,&M->level[level+1]);
    metric_mg_vcycle(M,level+1);
    metric_mg_prolong_float(L,&M->level[level+1]);
    metric_mg_level_smooth(L,M->post_sweeps,1);
}

static void metric_mg_fine_smooth(double *x,const double *rhs,const double *diag,
                                  int H,int W,double wh,double wv,double wd,
                                  int sweeps,int reverse){
    for(int sweep=0;sweep<sweeps;sweep++)for(int order=0;order<4;order++){
        int color=reverse?3-order:order;
        int row_parity=color>>1,column_parity=color&1;
        int r;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for(r=row_parity;r<H;r+=2)for(int c=column_parity;c<W;c+=2){
            size_t v=(size_t)r*(size_t)W+(size_t)c;
            double num=rhs[v];
            if(c>0)num+=wh*x[v-1];if(c+1<W)num+=wh*x[v+1];
            if(r>0)num+=wv*x[v-(size_t)W];if(r+1<H)num+=wv*x[v+(size_t)W];
            if(r>0&&c>0)num+=wd*x[v-(size_t)W-1];
            if(r+1<H&&c+1<W)num+=wd*x[v+(size_t)W+1];
            x[v]=num/diag[v];
        }
    }
}

static double metric_mg_fine_residual(const double *x,const double *rhs,
                                      const double *diag,int H,int W,
                                      double wh,double wv,double wd,double *residual){
    size_t n=(size_t)H*(size_t)W;double r2=0.0,b2=0.0;int kk;
    metric_matvec(x,residual,diag,H,W,wh,wv,wd);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:r2,b2)
#endif
    for(kk=0;kk<(int)n;kk++){
        size_t v=(size_t)kk;residual[v]=rhs[v]-residual[v];
        r2+=residual[v]*residual[v];b2+=rhs[v]*rhs[v];
    }
    return sqrt(r2)/fmax(sqrt(b2),1.0);
}

static void metric_mg_restrict_double(const double *fine_residual,int H,int W,
                                      MetricMGLevel *coarse){
    int kk;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for(kk=0;kk<(int)coarse->n;kk++){
        size_t p=(size_t)kk;int pr=(int)(p/(size_t)coarse->W),pc=(int)(p%(size_t)coarse->W);
        double sum=0.0;
        for(int dr=0;dr<2;dr++)for(int dc=0;dc<2;dc++){
            int r=2*pr+dr,c=2*pc+dc;
            if(r<H&&c<W)sum+=fine_residual[(size_t)r*(size_t)W+(size_t)c];
        }
        coarse->rhs[p]=(float)sum;coarse->x[p]=0.0f;
    }
}

static void metric_mg_prolong_double(double *fine_x,int H,int W,
                                     const MetricMGLevel *coarse){
    size_t n=(size_t)H*(size_t)W;int kk;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for(kk=0;kk<(int)n;kk++){
        size_t v=(size_t)kk;int r=(int)(v/(size_t)W),c=(int)(v%(size_t)W);
        size_t p=(size_t)(r>>1)*(size_t)coarse->W+(size_t)(c>>1);
        fine_x[v]+=(double)coarse->x[p];
    }
}

static int metric_mg_solve(MetricMG *M,double *x,const double *rhs,const double *diag,
                           int H,int W,double wh,double wv,double wd,int cycles,
                           double *scratch,double *out_before,double *out_after){
    if(M->nlevel==0||cycles<=0){if(out_before)*out_before=0.0;if(out_after)*out_after=0.0;return 0;}
    double before=metric_mg_fine_residual(x,rhs,diag,H,W,wh,wv,wd,scratch),after=before;
    for(int cycle=0;cycle<cycles&&after>1e-14;cycle++){
        metric_mg_fine_smooth(x,rhs,diag,H,W,wh,wv,wd,M->pre_sweeps,0);
        metric_mg_fine_residual(x,rhs,diag,H,W,wh,wv,wd,scratch);
        metric_mg_restrict_double(scratch,H,W,&M->level[0]);
        metric_mg_vcycle(M,0);
        metric_mg_prolong_double(x,H,W,&M->level[0]);
        metric_mg_fine_smooth(x,rhs,diag,H,W,wh,wv,wd,M->post_sweeps,1);
        after=metric_mg_fine_residual(x,rhs,diag,H,W,wh,wv,wd,scratch);
        if(!isfinite(after))return -1;
    }
    if(out_before)*out_before=before;if(out_after)*out_after=after;
    return 0;
}

static int metric_precondition(MetricMG *M,const double *residual,double *z,
                               const double *diag,int H,int W,double wh,double wv,
                               double wd,int cycles,double *scratch,
                               double *out_before,double *out_after){
    size_t n=(size_t)H*(size_t)W;
    if(M&&M->nlevel>0&&cycles>0){
        memset(z,0,n*sizeof *z);
        return metric_mg_solve(M,z,residual,diag,H,W,wh,wv,wd,cycles,scratch,
                               out_before,out_after);
    }
    for(size_t i=0;i<n;i++)z[i]=residual[i]/diag[i];
    if(out_before)*out_before=0.0;if(out_after)*out_after=0.0;
    return 0;
}

static int metric_pcg(double *x,const double *rhs,const double *diag,
                      int H,int W,double wh,double wv,double wd,
                      int max_iter,double tol,double *r,double *z,double *d,
                      double *Ad,MetricMG *mg,int mg_cycles,
                      int *out_iter,double *out_relres,
                      double *out_mg_before,double *out_mg_after){
    size_t n=(size_t)H*(size_t)W;
    metric_matvec(x,Ad,diag,H,W,wh,wv,wd);
    double r2=0.0,rz=0.0;
    int kk;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:r2)
#endif
    for(kk=0;kk<(int)n;kk++){
        size_t k=(size_t)kk;
        r[k]=rhs[k]-Ad[k];r2+=r[k]*r[k];
    }
    /* Normalizing by ||rhs|| makes a translated scroll (coordinates around
     * 3000) appear solved while its correction is still visibly moving.  Use
     * the initial residual instead, with an absolute floor below one voxel. */
    double norm=fmax(sqrt(r2),1.0),rel=sqrt(r2)/norm;
    double mg_before=0.0,mg_after=0.0;
    if(rel>tol){
        if(metric_precondition(mg,r,z,diag,H,W,wh,wv,wd,mg_cycles,Ad,
                               &mg_before,&mg_after)!=0)return -1;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:rz)
#endif
        for(kk=0;kk<(int)n;kk++){size_t k=(size_t)kk;d[k]=z[k];rz+=r[k]*z[k];}
    }
    int it=0;
    for(it=0;it<max_iter&&rel>tol;it++){
        metric_matvec(d,Ad,diag,H,W,wh,wv,wd);
        double dAd=0.0;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:dAd)
#endif
        for(kk=0;kk<(int)n;kk++){
            size_t k=(size_t)kk;dAd+=d[k]*Ad[k];
        }
        if(!(dAd>0.0)||!isfinite(dAd)||!isfinite(rz))return -1;
        double alpha=rz/dAd;
        r2=0.0;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:r2)
#endif
        for(kk=0;kk<(int)n;kk++){
            size_t k=(size_t)kk;x[k]+=alpha*d[k];r[k]-=alpha*Ad[k];r2+=r[k]*r[k];
        }
        rel=sqrt(r2)/norm;
        if(rel<=tol){it++;break;}
        double rz2=0.0;
        if(metric_precondition(mg,r,z,diag,H,W,wh,wv,wd,mg_cycles,Ad,
                               &mg_before,&mg_after)!=0)return -1;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:rz2)
#endif
        for(kk=0;kk<(int)n;kk++){
            size_t k=(size_t)kk;rz2+=r[k]*z[k];
        }
        if(!isfinite(rz2)||!(rz>0.0))return -1;
        double beta=rz2/rz;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for(kk=0;kk<(int)n;kk++){
            size_t k=(size_t)kk;d[k]=z[k]+beta*d[k];
        }
        rz=rz2;
    }
    if(out_iter)*out_iter=it;if(out_relres)*out_relres=rel;
    if(out_mg_before)*out_mg_before=mg_before;if(out_mg_after)*out_mg_after=mg_after;
    return isfinite(rel)?0:-1;
}

static int metric_trial_valid(const double *p,const double *candidate,
                              int H,int W,double alpha,double minimum_ratio){
    int bad=0;
    int r;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(|:bad)
#endif
    for(r=0;r<H-1;r++)for(int c=0;c+1<W;c++){
        size_t a=(size_t)r*(size_t)W+(size_t)c;
        size_t b=a+1,cc=a+(size_t)W,d=cc+1;
        size_t tri[2][3]={{a,b,d},{a,d,cc}};
        for(int t=0;t<2;t++){
            double oldp[3][3],newp[3][3];
            for(int q=0;q<3;q++)for(int ch=0;ch<3;ch++){
                size_t k=tri[t][q]*3+(size_t)ch;
                oldp[q][ch]=p[k];
                newp[q][ch]=p[k]+alpha*(candidate[k]-p[k]);
            }
            double e10[3],e20[3],e11[3],e21[3],n0[3],n1[3];
            for(int ch=0;ch<3;ch++){
                e10[ch]=oldp[1][ch]-oldp[0][ch];e20[ch]=oldp[2][ch]-oldp[0][ch];
                e11[ch]=newp[1][ch]-newp[0][ch];e21[ch]=newp[2][ch]-newp[0][ch];
            }
            n0[0]=e10[1]*e20[2]-e10[2]*e20[1];
            n0[1]=e10[2]*e20[0]-e10[0]*e20[2];
            n0[2]=e10[0]*e20[1]-e10[1]*e20[0];
            n1[0]=e11[1]*e21[2]-e11[2]*e21[1];
            n1[1]=e11[2]*e21[0]-e11[0]*e21[2];
            n1[2]=e11[0]*e21[1]-e11[1]*e21[0];
            double a0=n0[0]*n0[0]+n0[1]*n0[1]+n0[2]*n0[2];
            double a1=n1[0]*n1[0]+n1[1]*n1[1]+n1[2]*n1[2];
            double dot=n0[0]*n1[0]+n0[1]*n1[1]+n0[2]*n1[2];
            if(!isfinite(a1)||a1<1e-24||(a0>=1e-24&&dot<minimum_ratio*a0))bad=1;
        }
    }
    return !bad;
}

/* A single bad triangle must not throttle an entire multi-million-vertex
 * ribbon.  Barrier scales live on a much coarser regular control grid and are
 * bilinearly interpolated at vertices.  A bad triangle halves a one-control-
 * ring neighbourhood, giving a smooth ~64-cell transition instead of a cut
 * where adjacent per-vertex line-search factors disagree. */
#define METRIC_BARRIER_TILE 64
#define METRIC_BARRIER_HALVING_ROUNDS 24
#define METRIC_BARRIER_ZERO_ROUNDS 8
/* Below this radius the circumferential direction is ill-conditioned, so the
 * radial-normal prior applies only its V=z-hat target there. */
#define METRIC_AXIS_MIN_RADIUS 1.0
/* A candidate triangle must get axially worse by more than this before the
 * no-fold guard trims the step (guarantees a feasible step at scale 0). */
#define METRIC_AXIS_FOLD_EPS 1e-6
static double metric_barrier_control_at(const float *control,int TH,int TW,
                                        int H,int W,int r,int c){
    int gr=r/METRIC_BARRIER_TILE,gc=c/METRIC_BARRIER_TILE;
    if(gr>=TH-1)gr=TH-2;if(gc>=TW-1)gc=TW-2;
    int r0=gr*METRIC_BARRIER_TILE,c0=gc*METRIC_BARRIER_TILE;
    int r1=(gr+1)*METRIC_BARRIER_TILE,c1=(gc+1)*METRIC_BARRIER_TILE;
    if(r1>=H)r1=H-1;if(c1>=W)c1=W-1;
    double tr=r1>r0?(double)(r-r0)/(double)(r1-r0):0.0;
    double tc=c1>c0?(double)(c-c0)/(double)(c1-c0):0.0;
    size_t a=(size_t)gr*(size_t)TW+(size_t)gc,b=a+1;
    size_t q=a+(size_t)TW,d=q+1;
    double top=(1.0-tc)*control[a]+tc*control[b];
    double bottom=(1.0-tc)*control[q]+tc*control[d];
    return (1.0-tr)*top+tr*bottom;
}

static int metric_local_barrier_scales(const double *p,const double *displacement,
                                       int H,int W,double minimum_ratio,double maximum_step,
                                       double axis_fold_limit,
                                       float *scale,int *bad,int *dilated,
                                       int *out_rounds,double *out_min_scale,
                                       int *out_zero_fallbacks){
    size_t n=(size_t)H*(size_t)W;int kk;
    int TH=((H-1)+METRIC_BARRIER_TILE-1)/METRIC_BARRIER_TILE+1;
    int TW=((W-1)+METRIC_BARRIER_TILE-1)/METRIC_BARRIER_TILE+1;
    size_t tn=(size_t)TH*(size_t)TW;
    float *control=(float *)malloc(tn*sizeof *control);
    if(!control)return -1;
    for(size_t i=0;i<tn;i++)control[i]=1.0f;
    int rounds=0,zero_fallbacks=0;
    for(;;){
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for(kk=0;kk<(int)n;kk++){
            size_t i=(size_t)kk;int r=(int)(i/(size_t)W),c=(int)(i%(size_t)W);
            double m2=0.0;
            for(int ch=0;ch<3;ch++){double d=displacement[i*3+(size_t)ch];m2+=d*d;}
            double length=sqrt(m2),s=metric_barrier_control_at(control,TH,TW,H,W,r,c);
            if(length>maximum_step)s=fmin(s,maximum_step/length);
            scale[i]=(float)s;
        }
        memset(bad,0,tn*sizeof *bad);
        int any=0;int r;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(|:any)
#endif
        for(r=0;r<H-1;r++)for(int c=0;c+1<W;c++){
            size_t a=(size_t)r*(size_t)W+(size_t)c,b=a+1,q=a+(size_t)W,d=q+1;
            size_t tri[2][3]={{a,b,d},{a,d,q}};
            for(int t=0;t<2;t++){
                double oldp[3][3],newp[3][3];
                for(int v=0;v<3;v++)for(int ch=0;ch<3;ch++){
                    size_t vertex=tri[t][v],k=vertex*3+(size_t)ch;
                    oldp[v][ch]=p[k];
                    newp[v][ch]=p[k]+(double)scale[vertex]*displacement[k];
                }
                double e10[3],e20[3],e11[3],e21[3],n0[3],n1[3];
                for(int ch=0;ch<3;ch++){
                    e10[ch]=oldp[1][ch]-oldp[0][ch];e20[ch]=oldp[2][ch]-oldp[0][ch];
                    e11[ch]=newp[1][ch]-newp[0][ch];e21[ch]=newp[2][ch]-newp[0][ch];
                }
                n0[0]=e10[1]*e20[2]-e10[2]*e20[1];
                n0[1]=e10[2]*e20[0]-e10[0]*e20[2];
                n0[2]=e10[0]*e20[1]-e10[1]*e20[0];
                n1[0]=e11[1]*e21[2]-e11[2]*e21[1];
                n1[1]=e11[2]*e21[0]-e11[0]*e21[2];
                n1[2]=e11[0]*e21[1]-e11[1]*e21[0];
                double a0=n0[0]*n0[0]+n0[1]*n0[1]+n0[2]*n0[2];
                double a1=n1[0]*n1[0]+n1[1]*n1[1]+n1[2]*n1[2];
                double dot=n0[0]*n1[0]+n0[1]*n1[1]+n0[2]*n1[2];
                /* No-fold guard: forbid a step that makes a triangle's normal
                 * MORE axial (a bucket lid) past the limit.  Ratchet only -- at
                 * scale 0, axi1->axi0 so this never removes the feasible step. */
                int axial_fold=0;
                if(axis_fold_limit>0.0&&axis_fold_limit<1.0&&a1>1e-24){
                    double axi0=a0>1e-24?fabs(n0[0])/sqrt(a0):0.0;
                    double axi1=fabs(n1[0])/sqrt(a1);
                    axial_fold=axi1>axis_fold_limit&&axi1>axi0+METRIC_AXIS_FOLD_EPS;
                }
                if(!isfinite(a1)||a1<1e-24||(a0>=1e-24&&dot<minimum_ratio*a0)||axial_fold){
                    any=1;
#ifdef _OPENMP
#pragma omp critical(metric_barrier_mark_double)
#endif
                    {
                        for(int v=0;v<3;v++){
                            int vr=(int)(tri[t][v]/(size_t)W),vc=(int)(tri[t][v]%(size_t)W);
                            int gr=vr/METRIC_BARRIER_TILE,gc=vc/METRIC_BARRIER_TILE;
                            if(gr>=TH-1)gr=TH-2;if(gc>=TW-1)gc=TW-2;
                            size_t node=(size_t)gr*(size_t)TW+(size_t)gc;
                            bad[node]=bad[node+1]=bad[node+(size_t)TW]=bad[node+(size_t)TW+1]=1;
                        }
                    }
                }
            }
        }
        if(!any)break;
        if(rounds>=METRIC_BARRIER_HALVING_ROUNDS+METRIC_BARRIER_ZERO_ROUNDS){
            /* A line search always contains alpha=0.  Retain the current valid
             * iterate rather than turning an exhausted local search into a
             * fatal unwrap error.  Return 1 so the caller records a genuine
             * barrier stall instead of declaring movement convergence. */
            for(size_t i=0;i<n;i++)scale[i]=0.0f;
            free(control);
            if(out_rounds)*out_rounds=rounds;
            if(out_min_scale)*out_min_scale=0.0;
            if(out_zero_fallbacks)*out_zero_fallbacks=zero_fallbacks+1;
            return 1;
        }
        for(size_t k=0;k<tn;k++){
            int row=(int)(k/(size_t)TW),c=(int)(k%(size_t)TW);
            int marked=0;
            for(int dr=-1;dr<=1&&!marked;dr++)for(int dc=-1;dc<=1;dc++){
                int rr=row+dr,cc=c+dc;
                if(rr>=0&&rr<TH&&cc>=0&&cc<TW&&bad[(size_t)rr*(size_t)TW+(size_t)cc]){
                    marked=1;break;
                }
            }
            dilated[k]=marked;
        }
        for(size_t k=0;k<tn;k++)if(dilated[k]){
            if(rounds>=METRIC_BARRIER_HALVING_ROUNDS)control[k]=0.0f;
            else control[k]*=0.5f;
        }
        if(rounds>=METRIC_BARRIER_HALVING_ROUNDS)zero_fallbacks++;
        rounds++;
    }
    double minimum=1.0;
    for(size_t i=0;i<n;i++)if((double)scale[i]<minimum)minimum=scale[i];
    free(control);
    if(out_rounds)*out_rounds=rounds;if(out_min_scale)*out_min_scale=minimum;
    if(out_zero_fallbacks)*out_zero_fallbacks=zero_fallbacks;
    return 0;
}

static void metric_edge_error(const double *p,int H,int W,double du,double dv,
                              double *out_rms,double *out_max){
    double ss=0.0,mx=0.0;size_t count=0;
    double ld=hypot(du,dv);
    for(int r=0;r<H;r++)for(int c=0;c<W;c++){
        size_t i=(size_t)r*(size_t)W+(size_t)c;
        size_t js[3];double rest[3];int n=0;
        if(c+1<W){js[n]=i+1;rest[n++]=du;}
        if(r+1<H){js[n]=i+(size_t)W;rest[n++]=dv;}
        if(r+1<H&&c+1<W){js[n]=i+(size_t)W+1;rest[n++]=ld;}
        for(int q=0;q<n;q++){
            size_t j=js[q];double d0=p[i*3]-p[j*3],d1=p[i*3+1]-p[j*3+1],
                d2=p[i*3+2]-p[j*3+2];
            double length=sqrt(d0*d0+d1*d1+d2*d2);
            if(length>0.0&&rest[q]>0.0){double e=fabs(log(length/rest[q]));
                ss+=e*e;if(e>mx)mx=e;count++;}
        }
    }
    *out_rms=count?sqrt(ss/(double)count):0.0;*out_max=mx;
}

static double metric_repair_weight_at(const QuadStripArapOpts *o,size_t i){
    if(!o->repair_weight)return 0.0;
    double w=(double)o->repair_weight[i];
    if(!isfinite(w)||w<=0.0)return 0.0;
    return w>=1.0?1.0:w;
}

/* Estimate the UV->XYZ orthonormal frame of the current iterate.  Keeping this
 * in one routine is important for marble repair: the fair reference and every
 * ordinary ARAP local step must use exactly the same rest-edge convention. */
static int metric_compute_frames(const double *p,const float *uv,int H,int W,
                                 double wh,double wv,double wd,double *frame){
    size_t n=(size_t)H*(size_t)W;int frame_fail=0;int kk;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:frame_fail)
#endif
    for(kk=0;kk<(int)n;kk++){
        size_t i=(size_t)kk;int rr=(int)(i/(size_t)W),cc=(int)(i%(size_t)W);
        double S[6]={0,0,0,0,0,0};
#define METRIC_FRAME_ACCUM(J,WEIGHT) do { \
            size_t _j=(J);double _du=(double)uv[i*2]-(double)uv[_j*2]; \
            double _dv=(double)uv[i*2+1]-(double)uv[_j*2+1]; \
            for(int _ch=0;_ch<3;_ch++){ \
                double _ec=p[i*3+(size_t)_ch]-p[_j*3+(size_t)_ch]; \
                S[_ch*2]+=(WEIGHT)*_ec*_du; \
                S[_ch*2+1]+=(WEIGHT)*_ec*_dv; \
            } \
        } while(0)
        if(cc>0)METRIC_FRAME_ACCUM(i-1,wh);
        if(cc+1<W)METRIC_FRAME_ACCUM(i+1,wh);
        if(rr>0)METRIC_FRAME_ACCUM(i-(size_t)W,wv);
        if(rr+1<H)METRIC_FRAME_ACCUM(i+(size_t)W,wv);
        if(rr>0&&cc>0)METRIC_FRAME_ACCUM(i-(size_t)W-1,wd);
        if(rr+1<H&&cc+1<W)METRIC_FRAME_ACCUM(i+(size_t)W+1,wd);
#undef METRIC_FRAME_ACCUM
        double F[6];
        if(metric_closest_frame(S,&frame[i*6],F)!=0)frame_fail++;
        else memcpy(&frame[i*6],F,sizeof F);
    }
    return frame_fail? -1:0;
}

/* Build a low-frequency, screened tangent-frame target in UV.  This is the
 * missing bending prior in ordinary metric ARAP: edge lengths alone cannot
 * distinguish a smooth scroll from an isometric crumple.  The repair mask is
 * soft.  Clean vertices remain Dirichlet data; confidence below one retains a
 * proportional amount of the original local frame, making false positives
 * conservative.
 *
 * A single fine-grid relaxation only communicates clean orientation O(sqrt N)
 * cells in N sweeps. Real marble bands are hundreds of vertices wide, so that
 * merely polishes the crumple. The hierarchy below restricts both confidence
 * and clean-frame data, solves coarse-to-fine, then performs the same projected
 * red/black relaxation at every level. Thus distant clean sides influence the
 * middle in one setup pass while the terminal ARAP global step remains the
 * metric authority. */
#define METRIC_REPAIR_FRAME_LEVELS 12
typedef struct MetricRepairFrameLevel {
    int H,W;
    float *initial,*target,*mask;
} MetricRepairFrameLevel;

static void metric_repair_frame_level_free(MetricRepairFrameLevel *L){
    if(!L)return;free(L->initial);free(L->target);free(L->mask);memset(L,0,sizeof *L);
}

static int metric_repair_frame_level_alloc(MetricRepairFrameLevel *L,int H,int W){
    memset(L,0,sizeof *L);size_t n=(size_t)H*(size_t)W;
    if(H<2||W<2||n>SIZE_MAX/(6*sizeof(float)))return -1;
    L->H=H;L->W=W;L->initial=(float *)malloc(n*6*sizeof *L->initial);
    L->target=(float *)malloc(n*6*sizeof *L->target);
    L->mask=(float *)malloc(n*sizeof *L->mask);
    if(!L->initial||!L->target||!L->mask){metric_repair_frame_level_free(L);return -1;}
    return 0;
}

static int metric_repair_frame_project(const double S[6],const double previous[6],
                                       float *out){
    double F[6];if(metric_closest_frame(S,previous,F)!=0)return -1;
    for(int q=0;q<6;q++)out[q]=(float)F[q];return 0;
}

static int metric_repair_frame_coarsen(const MetricRepairFrameLevel *fine,
                                       MetricRepairFrameLevel *coarse){
    for(int r=0;r<coarse->H;r++)for(int c=0;c<coarse->W;c++){
        double S[6]={0,0,0,0,0,0},fallback[6]={0,0,0,0,0,0};
        double data_sum=0.0,mask_sum=0.0;int count=0;double previous[6];
        int nr=2*r;if(nr>=fine->H)nr=fine->H-1;
        int nc=2*c;if(nc>=fine->W)nc=fine->W-1;
        size_t nearest=(size_t)nr*(size_t)fine->W+(size_t)nc;
        for(int q=0;q<6;q++)previous[q]=fine->initial[nearest*6+(size_t)q];
        for(int dr=0;dr<2;dr++)for(int dc=0;dc<2;dc++){
            int fr=2*r+dr,fc=2*c+dc;if(fr>=fine->H||fc>=fine->W)continue;
            size_t i=(size_t)fr*(size_t)fine->W+(size_t)fc;
            double m=fmin(fmax((double)fine->mask[i],0.0),1.0),data=1.0-m;
            mask_sum+=m;data_sum+=data;count++;
            for(int q=0;q<6;q++){
                double v=fine->initial[i*6+(size_t)q];S[q]+=data*v;fallback[q]+=v;
            }
        }
        if(count<1)return -1;
        for(int q=0;q<6;q++)S[q]=data_sum>1e-8?S[q]/data_sum:fallback[q]/(double)count;
        size_t k=(size_t)r*(size_t)coarse->W+(size_t)c;
        if(metric_repair_frame_project(S,previous,&coarse->initial[k*6])!=0)return -1;
        memcpy(&coarse->target[k*6],&coarse->initial[k*6],6*sizeof(float));
        coarse->mask[k]=(float)(mask_sum/(double)count);
    }
    return 0;
}

static void metric_repair_frame_prolong(const MetricRepairFrameLevel *coarse,
                                        MetricRepairFrameLevel *fine){
    for(int r=0;r<fine->H;r++)for(int c=0;c<fine->W;c++){
        size_t i=(size_t)r*(size_t)fine->W+(size_t)c;
        double m=fmin(fmax((double)fine->mask[i],0.0),1.0);
        if(m<=0.0){memcpy(&fine->target[i*6],&fine->initial[i*6],6*sizeof(float));continue;}
        double fr=(double)r*(double)(coarse->H-1)/(double)(fine->H-1);
        double fc=(double)c*(double)(coarse->W-1)/(double)(fine->W-1);
        int r0=(int)floor(fr),c0=(int)floor(fc),r1=r0+1,c1=c0+1;
        if(r1>=coarse->H)r1=coarse->H-1;if(c1>=coarse->W)c1=coarse->W-1;
        double tr=fr-r0,tc=fc-c0,S[6],previous[6];
        size_t a=(size_t)r0*(size_t)coarse->W+(size_t)c0;
        size_t b=(size_t)r0*(size_t)coarse->W+(size_t)c1;
        size_t q=(size_t)r1*(size_t)coarse->W+(size_t)c0;
        size_t d=(size_t)r1*(size_t)coarse->W+(size_t)c1;
        for(int ch=0;ch<6;ch++){
            double top=(1.0-tc)*coarse->target[a*6+(size_t)ch]+tc*coarse->target[b*6+(size_t)ch];
            double bottom=(1.0-tc)*coarse->target[q*6+(size_t)ch]+tc*coarse->target[d*6+(size_t)ch];
            double cv=(1.0-tr)*top+tr*bottom;
            previous[ch]=fine->initial[i*6+(size_t)ch];
            S[ch]=(1.0-m)*previous[ch]+m*cv;
        }
        if(metric_repair_frame_project(S,previous,&fine->target[i*6])!=0)
            memcpy(&fine->target[i*6],&fine->initial[i*6],6*sizeof(float));
    }
}

static int metric_repair_frame_smooth(const QuadStripArapOpts *o,
                                      MetricRepairFrameLevel *L){
    size_t n=(size_t)L->H*(size_t)L->W;
    for(int sweep=0;sweep<o->repair_frame_sweeps;sweep++)for(int color=0;color<2;color++){
        int frame_fail=0,kk;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:frame_fail)
#endif
        for(kk=0;kk<(int)n;kk++){
            size_t i=(size_t)kk;int rr=(int)(i/(size_t)L->W),cc=(int)(i%(size_t)L->W);
            double m=fmin(fmax((double)L->mask[i],0.0),1.0);
            if(m<=0.0||((rr+cc)&1)!=color)continue;
            size_t neighbour[4];int nn=0;
            if(cc>0)neighbour[nn++]=i-1;if(cc+1<L->W)neighbour[nn++]=i+1;
            if(rr>0)neighbour[nn++]=i-(size_t)L->W;if(rr+1<L->H)neighbour[nn++]=i+(size_t)L->W;
            if(nn==0)continue;
            double data=(1.0-m)+m*o->repair_frame_screen;
            double smooth=m,den=data+smooth*(double)nn,S[6],previous[6];
            for(int q=0;q<6;q++){
                double sum=data*L->initial[i*6+(size_t)q];
                for(int e=0;e<nn;e++)sum+=smooth*L->target[neighbour[e]*6+(size_t)q];
                S[q]=sum/den;previous[q]=L->target[i*6+(size_t)q];
            }
            if(metric_repair_frame_project(S,previous,&L->target[i*6])!=0)frame_fail++;
        }
        if(frame_fail)return -1;
    }
    return 0;
}

static int metric_build_repair_frame(const QuadStripArapOpts *o,
                                     const double *initial_frame,
                                     int H,int W,QuadStripArapStats *stats,
                                     float **out_target){
    *out_target=NULL;
    if(!o->repair_weight||o->repair_frame_sweeps<=0||o->repair_frame_blend<=0.0)
        return 0;
    size_t n=(size_t)H*(size_t)W;
    MetricRepairFrameLevel L[METRIC_REPAIR_FRAME_LEVELS];memset(L,0,sizeof L);
    int levels=1;if(metric_repair_frame_level_alloc(&L[0],H,W)!=0)return -1;
    size_t active=0;double weight_sum=0.0;
    for(size_t i=0;i<n;i++){
        double w=metric_repair_weight_at(o,i);
        if(w>0.0){active++;weight_sum+=w;}
        L[0].mask[i]=(float)w;
        for(int q=0;q<6;q++)L[0].initial[i*6+(size_t)q]
            =L[0].target[i*6+(size_t)q]=(float)initial_frame[i*6+(size_t)q];
    }
    if(stats){stats->repair_vertices=active;stats->repair_weight_sum=weight_sum;}
    if(active==0){metric_repair_frame_level_free(&L[0]);return 0;}
    while(levels<METRIC_REPAIR_FRAME_LEVELS){
        int CH=(L[levels-1].H+1)/2,CW=(L[levels-1].W+1)/2;
        if(CH<2||CW<2||((size_t)CH*(size_t)CW)<16)break;
        if(metric_repair_frame_level_alloc(&L[levels],CH,CW)!=0)goto fail;
        if(metric_repair_frame_coarsen(&L[levels-1],&L[levels])!=0)goto fail;
        levels++;
    }
    for(int l=levels-1;l>=0;l--){
        if(l+1<levels)metric_repair_frame_prolong(&L[l+1],&L[l]);
        if(metric_repair_frame_smooth(o,&L[l])!=0)goto fail;
    }
    *out_target=L[0].target;L[0].target=NULL;
    for(int l=0;l<levels;l++)metric_repair_frame_level_free(&L[l]);
    return 0;
fail:
    for(int l=0;l<METRIC_REPAIR_FRAME_LEVELS;l++)metric_repair_frame_level_free(&L[l]);
    return -1;
}

static int metric_arap_core(float *verts,const float *anchor_verts,
                            const float *uv,const uint8_t *filled,
                            int H,int W,const QuadStripArapOpts *opts,
                            QuadStripArapStats *stats){
    QuadStripArapOpts o;if(opts)o=*opts;else QuadStripArap_defaults(&o);
    if(stats)memset(stats,0,sizeof *stats);
    if(stats)stats->minimum_step_scale=1.0;
    if(!verts||!anchor_verts||!uv||!filled||H<2||W<2||o.max_iterations<1||
       !(o.movement_tolerance>0.0)||o.cg_max_iterations<1||
       !(o.cg_relative_tolerance>0.0)||o.source_weight<0.0||o.fill_weight<0.0||
       o.multigrid_levels<1||o.multigrid_levels>METRIC_MG_MAX_LEVELS||
       o.multigrid_cycles<0||o.multigrid_pre_sweeps<0||o.multigrid_post_sweeps<0||
       o.multigrid_coarse_sweeps<1||
       o.nonlinear_multigrid_levels<1||o.nonlinear_multigrid_levels>8||
       o.nonlinear_coarse_iterations<1||o.nonlinear_min_vertices<0||
       !(o.source_weight>0.0||o.fill_weight>0.0)||!(o.huber_delta>0.0)||
       o.minimum_source_scale<0.0||o.minimum_source_scale>1.0||
       o.minimum_area_ratio<=0.0||o.minimum_area_ratio>1.0||
       !(o.maximum_vertex_step>0.0)||!isfinite(o.maximum_vertex_step)||
       (o.preserve_axial!=0&&o.preserve_axial!=1)||
       o.repair_source_scale<0.0||o.repair_source_scale>1.0||
       o.repair_frame_blend<0.0||o.repair_frame_blend>1.0||
       !(o.repair_frame_screen>0.0)||!isfinite(o.repair_frame_screen)||
       o.repair_frame_sweeps<0||o.repair_frame_sweeps>100000||
       !isfinite(o.axis_normal_weight)||o.axis_normal_weight<0.0||o.axis_normal_weight>1.0||
       !isfinite(o.axis_normal_radius)||o.axis_normal_radius<0.0||
       !isfinite(o.axis_fold_limit)||o.axis_fold_limit<0.0||o.axis_fold_limit>=1.0||
       (o.axis_normal_weight>0.0&&(!isfinite(o.axis_y)||!isfinite(o.axis_x))))return -1;
    size_t n=(size_t)H*(size_t)W;
    if(n>(size_t)INT_MAX)return -1; /* MSVC OpenMP 2.0 requires int loop indices. */
    double du=fabs((double)uv[2]-(double)uv[0]);
    double dv=fabs((double)uv[(size_t)W*2+1]-(double)uv[1]);
    if(!(du>0.0)||!(dv>0.0)||!isfinite(du)||!isfinite(dv))return -1;
    double wh=1.0/(du*du),wv=1.0/(dv*dv),wd=1.0/(du*du+dv*dv);
    double *original=(double *)malloc(n*3*sizeof *original);
    double *p=(double *)malloc(n*3*sizeof *p);
    double *candidate=(double *)malloc(n*3*sizeof *candidate);
    double *frame=(double *)malloc(n*6*sizeof *frame);
    double *anchor=(double *)malloc(n*sizeof *anchor);
    double *diag=(double *)malloc(n*sizeof *diag);
    double *rhs=(double *)malloc(n*sizeof *rhs);
    double *x=(double *)malloc(n*sizeof *x);
    double *r=(double *)malloc(n*sizeof *r),*z=(double *)malloc(n*sizeof *z);
    double *d=(double *)malloc(n*sizeof *d),*Ad=(double *)malloc(n*sizeof *Ad);
    float *step_scale=(float *)malloc(n*sizeof *step_scale);
    size_t barrier_nodes=(size_t)(((H-1)+METRIC_BARRIER_TILE-1)/METRIC_BARRIER_TILE+1)*
                         (size_t)(((W-1)+METRIC_BARRIER_TILE-1)/METRIC_BARRIER_TILE+1);
    int *barrier_bad=(int *)malloc(barrier_nodes*sizeof *barrier_bad);
    int *barrier_dilated=(int *)malloc(barrier_nodes*sizeof *barrier_dilated);
    float *repair_frame=NULL;
    MetricMG mg;
    if(!original||!p||!candidate||!frame||!anchor||!diag||!rhs||!x||!r||!z||!d||!Ad||
       !step_scale||!barrier_bad||!barrier_dilated){
        free(original);free(p);free(candidate);free(frame);free(anchor);free(diag);
        free(rhs);free(x);free(r);free(z);free(d);free(Ad);free(step_scale);
        free(barrier_bad);free(barrier_dilated);return -1;
    }
    if(metric_mg_init(&mg,H,W,o.multigrid_levels,o.multigrid_pre_sweeps,
                      o.multigrid_post_sweeps,o.multigrid_coarse_sweeps)!=0){
        free(original);free(p);free(candidate);free(frame);free(anchor);free(diag);
        free(rhs);free(x);free(r);free(z);free(d);free(Ad);free(step_scale);
        free(barrier_bad);free(barrier_dilated);return -1;
    }
    if(stats)stats->multigrid_levels_used=mg.nlevel+1;
    for(size_t i=0;i<n;i++){
        for(int ch=0;ch<3;ch++){
            original[i*3+(size_t)ch]=anchor_verts[i*3+(size_t)ch];
            p[i*3+(size_t)ch]=verts[i*3+(size_t)ch];
        }
        if(o.preserve_axial)p[i*3]=original[i*3];
        /* Stable fallback: rest-U initially points to x and rest-V to z in z/y/x order. */
        frame[i*6]=0.0;frame[i*6+1]=1.0;
        frame[i*6+2]=0.0;frame[i*6+3]=0.0;
        frame[i*6+4]=1.0;frame[i*6+5]=0.0;
    }
    if(o.repair_weight&&o.repair_frame_sweeps>0&&o.repair_frame_blend>0.0){
        if(metric_compute_frames(p,uv,H,W,wh,wv,wd,frame)!=0||
           metric_build_repair_frame(&o,frame,H,W,stats,&repair_frame)!=0){
            free(original);free(p);free(candidate);free(frame);free(anchor);free(diag);
            free(rhs);free(x);free(r);free(z);free(d);free(Ad);free(step_scale);
            free(barrier_bad);free(barrier_dilated);metric_mg_free(&mg);return -1;
        }
        if(o.verbose&&stats&&stats->repair_vertices>0)fprintf(stderr,
            "      metric ARAP marble repair: %zu vertices, weight sum %.3f, "
            "frame sweeps %d blend %.3g screen %.3g\n",
            stats->repair_vertices,stats->repair_weight_sum,o.repair_frame_sweeps,
            o.repair_frame_blend,o.repair_frame_screen);
    }
    double initial_rms=0.0,initial_max=0.0;
    metric_edge_error(p,H,W,du,dv,&initial_rms,&initial_max);
    if(stats){stats->initial_edge_log_rms=initial_rms;stats->initial_edge_log_max=initial_max;}
    int result=0;
    double previous_maxmove=INFINITY;
    double effective_tolerance=o.movement_tolerance;
    double first_maxmove=-1.0,stall_best=INFINITY;
    int stall_count=0;
    if(stats)stats->stop_reason=QUAD_STRIP_ARAP_STOP_MAX_ITERATIONS;
    for(int outer=0;outer<o.max_iterations;outer++){
        /* Robust soft attachment to the measured field.  Low-confidence movement
         * never removes a vertex; it only lets metric consistency outrank one
         * locally compressed observation. */
        int kk;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for(kk=0;kk<(int)n;kk++){
            size_t i=(size_t)kk;double aweight=o.fill_weight;
            if(!filled[i]){
                double q0=p[i*3]-original[i*3],q1=p[i*3+1]-original[i*3+1],
                    q2=p[i*3+2]-original[i*3+2];double disp=sqrt(q0*q0+q1*q1+q2*q2);
                double scale=disp>o.huber_delta?o.huber_delta/disp:1.0;
                if(scale<o.minimum_source_scale)scale=o.minimum_source_scale;
                aweight=o.source_weight*scale;
            }
            {
                double repair=metric_repair_weight_at(&o,i);
                aweight*=1.0-repair+repair*o.repair_source_scale;
            }
            /* Keep Dirichlet rows strongly attached without an enormous
             * penalty that would dominate the relative residual norm. The
             * proposal is stamped exactly below, before the orientation
             * barrier, so no finite-penalty error reaches output. */
            if(o.fixed_vertices&&o.fixed_vertices[i])
                aweight=fmax(o.source_weight,100.0);
            anchor[i]=aweight;
            int rr=(int)(i/(size_t)W),cc=(int)(i%(size_t)W);
            diag[i]=metric_stencil_sum(rr,cc,H,W,wh,wv,wd)+aweight;
        }
        if(metric_mg_rebuild(&mg,diag,H,W,wh,wv,wd)!=0){result=-1;break;}

        /* Local step: closest 3x2 orthonormal frame per ribbon vertex.  A
         * marble repair then rotates only the selected frames toward their
         * screened low-frequency UV target; ordinary ARAP is bit-for-bit on
         * the null-mask path. */
        if(metric_compute_frames(p,uv,H,W,wh,wv,wd,frame)!=0){result=-1;break;}
        if(repair_frame){
            int repair_fail=0;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:repair_fail)
#endif
            for(kk=0;kk<(int)n;kk++){
                size_t i=(size_t)kk;double m=metric_repair_weight_at(&o,i);
                double beta=o.repair_frame_blend*m;
                if(beta<=0.0)continue;
                double S[6],F[6];
                for(int q=0;q<6;q++)S[q]=(1.0-beta)*frame[i*6+(size_t)q]
                    +beta*(double)repair_frame[i*6+(size_t)q];
                if(metric_closest_frame(S,&frame[i*6],F)!=0)repair_fail++;
                else memcpy(&frame[i*6],F,sizeof F);
            }
            if(repair_fail){result=-1;break;}
        }
        /* Radial-normal prior: rotate each vertex frame toward its local cylinder
         * frame (V = axial z-hat, U = circumferential) so an axial "bucket-lid"
         * normal becomes high-energy.  Self-gating (target == current where the
         * wrap is already cylindrical) and bit-for-bit no-op when weight == 0. */
        if(o.axis_normal_weight>0.0){
            int axis_fail=0;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:axis_fail)
#endif
            for(kk=0;kk<(int)n;kk++){
                size_t i=(size_t)kk;
                double dy=p[i*3+1]-o.axis_y,dx=p[i*3+2]-o.axis_x;
                double R=hypot(dy,dx);
                double beta=o.axis_normal_weight;
                if(o.axis_normal_radius>0.0){
                    double t=R/o.axis_normal_radius;beta*=t>=1.0?0.0:(1.0-t);
                }
                if(!(beta>0.0))continue;
                double tgt[6];
                tgt[1]=frame[i*6+1]<0.0?-1.0:1.0;tgt[3]=0.0;tgt[5]=0.0; /* V -> +/- z-hat */
                if(R>METRIC_AXIS_MIN_RADIUS){
                    double c1=-dx/R,c2=dy/R;                           /* circumferential (0,-dx,dy)/R */
                    if(frame[i*6+2]*c1+frame[i*6+4]*c2<0.0){c1=-c1;c2=-c2;}
                    tgt[0]=0.0;tgt[2]=c1;tgt[4]=c2;                     /* U -> +/- circumferential */
                }else{
                    tgt[0]=frame[i*6];tgt[2]=frame[i*6+2];tgt[4]=frame[i*6+4]; /* keep U near axis */
                }
                double S[6],F[6];
                for(int q=0;q<6;q++)S[q]=(1.0-beta)*frame[i*6+(size_t)q]+beta*tgt[q];
                if(metric_closest_frame(S,&frame[i*6],F)!=0)axis_fail++;
                else memcpy(&frame[i*6],F,sizeof F);
            }
            if(axis_fail){result=-1;break;}
        }

        /* Global step: fixed-frame Poisson solve for z, y, and x. */
        int maximum_cg=0;double last_rel=0.0,last_mg_before=0.0,last_mg_after=0.0;
        /* An exact global solve is wasteful while the nonlinear local frames
         * are still changing by whole voxels.  This Eisenstat-Walker-style
         * forcing schedule keeps early solves inexact, then reaches the user's
         * requested tolerance before the physical stopping threshold. */
        double movement_ratio=previous_maxmove/effective_tolerance;
        double linear_tol=o.cg_relative_tolerance;
        if(movement_ratio>1e5)linear_tol=fmax(linear_tol,1e-2);
        else if(movement_ratio>1e4)linear_tol=fmax(linear_tol,3e-3);
        else if(movement_ratio>1e3)linear_tol=fmax(linear_tol,1e-3);
        else if(movement_ratio>1e2)linear_tol=fmax(linear_tol,1e-4);
        else if(movement_ratio>10.0)linear_tol=fmax(linear_tol,1e-5);
        for(int ch=0;ch<3;ch++){
            if(o.preserve_axial&&ch==0){
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
                for(kk=0;kk<(int)n;kk++)
                    candidate[(size_t)kk*3]=original[(size_t)kk*3];
                continue;
            }
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
            for(kk=0;kk<(int)n;kk++){
                size_t i=(size_t)kk;int rr=(int)(i/(size_t)W),cc=(int)(i%(size_t)W);
                double bi=anchor[i]*original[i*3+(size_t)ch];
#define METRIC_RHS(J,WEIGHT) do { \
                    size_t _j=(J);double _du=(double)uv[i*2]-(double)uv[_j*2]; \
                    double _dv=(double)uv[i*2+1]-(double)uv[_j*2+1]; \
                    double _ri=frame[i*6+(size_t)ch*2]*_du+frame[i*6+(size_t)ch*2+1]*_dv; \
                    double _rj=frame[_j*6+(size_t)ch*2]*_du+frame[_j*6+(size_t)ch*2+1]*_dv; \
                    bi+=(WEIGHT)*0.5*(_ri+_rj); \
                } while(0)
                if(cc>0)METRIC_RHS(i-1,wh);if(cc+1<W)METRIC_RHS(i+1,wh);
                if(rr>0)METRIC_RHS(i-(size_t)W,wv);if(rr+1<H)METRIC_RHS(i+(size_t)W,wv);
                if(rr>0&&cc>0)METRIC_RHS(i-(size_t)W-1,wd);
                if(rr+1<H&&cc+1<W)METRIC_RHS(i+(size_t)W+1,wd);
#undef METRIC_RHS
                rhs[i]=bi;x[i]=p[i*3+(size_t)ch];
            }
            double mg_before=0.0,mg_after=0.0;
            int cg_it=0;double rel=0.0;
            if(metric_pcg(x,rhs,diag,H,W,wh,wv,wd,o.cg_max_iterations,
                          linear_tol,r,z,d,Ad,&mg,o.multigrid_cycles,
                          &cg_it,&rel,&mg_before,&mg_after)!=0){result=-1;break;}
            last_mg_before=mg_before;last_mg_after=mg_after;
            if(stats&&mg.nlevel>0)stats->multigrid_cycles+=(cg_it+1)*o.multigrid_cycles;
            if(cg_it>maximum_cg)maximum_cg=cg_it;last_rel=rel;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
            for(kk=0;kk<(int)n;kk++)candidate[(size_t)kk*3+(size_t)ch]=x[(size_t)kk];
        }
        if(result)break;
        if(o.fixed_vertices){
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
            for(kk=0;kk<(int)n;kk++)if(o.fixed_vertices[(size_t)kk])
                for(int ch=0;ch<3;ch++)candidate[(size_t)kk*3+(size_t)ch]
                    =original[(size_t)kk*3+(size_t)ch];
        }
        if(stats){if(maximum_cg>stats->maximum_cg_iterations)stats->maximum_cg_iterations=maximum_cg;
            stats->final_cg_relative_residual=last_rel;
            stats->final_multigrid_residual_before=last_mg_before;
            stats->final_multigrid_residual_after=last_mg_after;}

        /* Keep the proposal as a displacement so a localized barrier can
         * scale each neighbourhood independently without another n*3 array. */
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for(kk=0;kk<(int)n;kk++)for(int ch=0;ch<3;ch++){
            size_t k=(size_t)kk*3+(size_t)ch;candidate[k]-=p[k];
        }
        int barrier_rounds=0,barrier_zero_fallbacks=0;
        double minimum_step_scale=1.0;
        int barrier_status=metric_local_barrier_scales(
            p,candidate,H,W,o.minimum_area_ratio,o.maximum_vertex_step,
            o.axis_fold_limit,
            step_scale,barrier_bad,barrier_dilated,&barrier_rounds,
            &minimum_step_scale,&barrier_zero_fallbacks);
        if(barrier_status<0){result=-2;break;}
        if(stats){stats->barrier_rejections+=barrier_rounds;
            stats->barrier_zero_fallbacks+=barrier_zero_fallbacks;
            if(minimum_step_scale<stats->minimum_step_scale)
                stats->minimum_step_scale=minimum_step_scale;}
        if(barrier_status>0){
            if(stats){stats->iterations=outer+1;stats->barrier_stalled=1;
                stats->stop_reason=QUAD_STRIP_ARAP_STOP_BARRIER;
                stats->maximum_vertex_movement=0.0;stats->rms_vertex_movement=0.0;}
            if(o.verbose)fprintf(stderr,
                "      metric ARAP iter %d: local barrier retained the current "
                "valid iterate after %d rounds (not convergence)\n",
                outer+1,barrier_rounds);
            break;
        }
        double maxmove=0.0,sumsq=0.0;
        /* `max` reductions require MSVC's LLVM OpenMP runtime, so calculate
         * the convergence norm once, portably, before the parallel update. */
        for(size_t i=0;i<n;i++){
            double m2=0.0;
            for(int ch=0;ch<3;ch++){
                size_t k=i*3+(size_t)ch;
                double dd=(double)step_scale[i]*candidate[k];
                m2+=dd*dd;
            }
            double mv=sqrt(m2);if(mv>maxmove)maxmove=mv;
        }
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:sumsq)
#endif
        for(kk=0;kk<(int)n;kk++){
            size_t i=(size_t)kk;double m2=0.0;
            for(int ch=0;ch<3;ch++){size_t k=i*3+(size_t)ch;
                double dd=(double)step_scale[i]*candidate[k];p[k]+=dd;m2+=dd*dd;}
            sumsq+=m2;
        }
        double rms=sqrt(sumsq/(double)n);
        previous_maxmove=maxmove;
        if(first_maxmove<0.0){
            first_maxmove=maxmove;
            if(o.movement_tolerance_rel>0.0){
                double rel=o.movement_tolerance_rel*first_maxmove;
                if(rel>effective_tolerance)effective_tolerance=rel;
            }
        }
        if(stats){stats->iterations=outer+1;stats->maximum_vertex_movement=maxmove;
            stats->rms_vertex_movement=rms;
            stats->effective_movement_tolerance=effective_tolerance;}
        if(o.verbose&&(outer<10||((outer+1)%100)==0||maxmove<=effective_tolerance))
            fprintf(stderr,"      metric ARAP iter %d: max/rms move %.9g/%.9g, "
                           "MG L%d %.3g->%.3g, CG <=%d, relres %.3g (target %.1g), "
                           "local barrier rounds %d min-scale %.3g\n",
                    outer+1,maxmove,rms,mg.nlevel+1,last_mg_before,last_mg_after,
                    maximum_cg,last_rel,linear_tol,barrier_rounds,minimum_step_scale);
        if(maxmove<=effective_tolerance&&linear_tol<=o.cg_relative_tolerance){
            if(stats){stats->converged=1;
                stats->stop_reason=maxmove<=o.movement_tolerance
                    ?QUAD_STRIP_ARAP_STOP_TOLERANCE
                    :QUAD_STRIP_ARAP_STOP_RELATIVE_TOLERANCE;}
            break;
        }
        /* Plateau detection: further iterations are provably useless when the
         * running-min max move has not improved by the stall fraction for a
         * whole window.  Recorded as stalled, NOT converged, and distinct from
         * the max-iterations cap so report counters stay honest. */
        if(o.movement_stall_window>0){
            if(maxmove<stall_best*(1.0-o.movement_stall_fraction)){
                stall_best=maxmove;stall_count=0;
            }else{
                if(maxmove<stall_best)stall_best=maxmove;
                stall_count++;
                if(stall_count>=o.movement_stall_window){
                    if(stats)stats->stop_reason=QUAD_STRIP_ARAP_STOP_STALLED;
                    if(o.verbose)fprintf(stderr,
                        "      metric ARAP iter %d: stalled (max move %.9g, "
                        "best %.9g, window %d)\n",
                        outer+1,maxmove,stall_best,o.movement_stall_window);
                    break;
                }
            }
        }
    }
    if(result==0){
        for(size_t i=0;i<n;i++)for(int ch=0;ch<3;ch++)verts[i*3+(size_t)ch]=(float)p[i*3+(size_t)ch];
        if(stats)metric_edge_error(p,H,W,du,dv,&stats->final_edge_log_rms,
                                   &stats->final_edge_log_max);
    }
    free(original);free(p);free(candidate);free(frame);free(anchor);free(diag);
    free(rhs);free(x);free(r);free(z);free(d);free(Ad);free(step_scale);
    free(barrier_bad);free(barrier_dilated);free(repair_frame);
    metric_mg_free(&mg);
    return result;
}

/* Nonlinear multigrid ------------------------------------------------------ */

typedef struct MetricArapCoarse {
    int H,W;
    float *anchor,*work,*uv;
    float *repair;
    uint8_t *filled;
} MetricArapCoarse;

static void metric_arap_coarse_free(MetricArapCoarse *C){
    if(!C)return;
    free(C->anchor);free(C->work);free(C->uv);free(C->repair);free(C->filled);
    memset(C,0,sizeof *C);
}

static void metric_bilinear_sample(const float *field,int H,int W,int channels,
                                   double row,double column,float *out){
    int r0=(int)floor(row),c0=(int)floor(column);
    int r1=r0+1<H?r0+1:r0,c1=c0+1<W?c0+1:c0;
    double tr=row-r0,tc=column-c0;
    size_t a=(size_t)r0*(size_t)W+(size_t)c0;
    size_t b=(size_t)r0*(size_t)W+(size_t)c1;
    size_t c=(size_t)r1*(size_t)W+(size_t)c0;
    size_t d=(size_t)r1*(size_t)W+(size_t)c1;
    for(int ch=0;ch<channels;ch++){
        double top=(1.0-tc)*field[a*(size_t)channels+(size_t)ch]
                  +tc*field[b*(size_t)channels+(size_t)ch];
        double bottom=(1.0-tc)*field[c*(size_t)channels+(size_t)ch]
                     +tc*field[d*(size_t)channels+(size_t)ch];
        out[ch]=(float)((1.0-tr)*top+tr*bottom);
    }
}

static int metric_arap_build_coarse(const float *work,const float *anchor,
                                    const float *uv,const uint8_t *filled,
                                    const float *repair,
                                    int H,int W,MetricArapCoarse *C){
    memset(C,0,sizeof *C);
    if(H<3||W<3)return 1;
    int CH=(H+1)/2,CW=(W+1)/2;
    if(CH<2||CW<2||CH==H||CW==W)return 1;
    size_t cn=(size_t)CH*(size_t)CW;
    if(cn>SIZE_MAX/(3*sizeof(float))||cn>SIZE_MAX/(2*sizeof(float)))return -1;
    C->H=CH;C->W=CW;
    C->anchor=(float *)malloc(cn*3*sizeof *C->anchor);
    C->work=(float *)malloc(cn*3*sizeof *C->work);
    C->uv=(float *)malloc(cn*2*sizeof *C->uv);
    if(repair)C->repair=(float *)malloc(cn*sizeof *C->repair);
    C->filled=(uint8_t *)malloc(cn*sizeof *C->filled);
    if(!C->anchor||!C->work||!C->uv||(repair&&!C->repair)||!C->filled){
        metric_arap_coarse_free(C);return -1;
    }
    for(int r=0;r<CH;r++)for(int c=0;c<CW;c++){
        double fr=(double)r*(double)(H-1)/(double)(CH-1);
        double fc=(double)c*(double)(W-1)/(double)(CW-1);
        size_t k=(size_t)r*(size_t)CW+(size_t)c;
        metric_bilinear_sample(anchor,H,W,3,fr,fc,&C->anchor[k*3]);
        metric_bilinear_sample(work,H,W,3,fr,fc,&C->work[k*3]);
        metric_bilinear_sample(uv,H,W,2,fr,fc,&C->uv[k*2]);
        if(repair)metric_bilinear_sample(repair,H,W,1,fr,fc,&C->repair[k]);
        int nr=(int)floor(fr+0.5),nc=(int)floor(fc+0.5);
        if(nr>=H)nr=H-1;if(nc>=W)nc=W-1;
        C->filled[k]=filled[(size_t)nr*(size_t)W+(size_t)nc];
    }
    return 0;
}

static void metric_prolong_candidate(const float *base,float *candidate,
                                     int H,int W,const MetricArapCoarse *C,
                                     double alpha){
    size_t n=(size_t)H*(size_t)W;int kk;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for(kk=0;kk<(int)n;kk++){
        size_t k=(size_t)kk;int r=(int)(k/(size_t)W),c=(int)(k%(size_t)W);
        double cr=(double)r*(double)(C->H-1)/(double)(H-1);
        double cc=(double)c*(double)(C->W-1)/(double)(W-1);
        int r0=(int)floor(cr),c0=(int)floor(cc);
        int r1=r0+1<C->H?r0+1:r0,c1=c0+1<C->W?c0+1:c0;
        double tr=cr-r0,tc=cc-c0;
        size_t a=(size_t)r0*(size_t)C->W+(size_t)c0;
        size_t b=(size_t)r0*(size_t)C->W+(size_t)c1;
        size_t q=(size_t)r1*(size_t)C->W+(size_t)c0;
        size_t d=(size_t)r1*(size_t)C->W+(size_t)c1;
        for(int ch=0;ch<3;ch++){
#define METRIC_DISP(I) ((double)C->work[(I)*3+(size_t)ch]-(double)C->anchor[(I)*3+(size_t)ch])
            double top=(1.0-tc)*METRIC_DISP(a)+tc*METRIC_DISP(b);
            double bottom=(1.0-tc)*METRIC_DISP(q)+tc*METRIC_DISP(d);
#undef METRIC_DISP
            candidate[k*3+(size_t)ch]=(float)((double)base[k*3+(size_t)ch]
                                                +alpha*((1.0-tr)*top+tr*bottom));
        }
    }
}

static int metric_float_mesh_valid(const float *base,const float *candidate,
                                   int H,int W,double minimum_ratio){
    int bad=0;int r;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(|:bad)
#endif
    for(r=0;r<H-1;r++)for(int c=0;c+1<W;c++){
        size_t a=(size_t)r*(size_t)W+(size_t)c,b=a+1,q=a+(size_t)W,d=q+1;
        size_t tri[2][3]={{a,b,d},{a,d,q}};
        for(int t=0;t<2;t++){
            double oldp[3][3],newp[3][3];
            for(int v=0;v<3;v++)for(int ch=0;ch<3;ch++){
                size_t k=tri[t][v]*3+(size_t)ch;
                oldp[v][ch]=base[k];newp[v][ch]=candidate[k];
            }
            double e10[3],e20[3],e11[3],e21[3],n0[3],n1[3];
            for(int ch=0;ch<3;ch++){
                e10[ch]=oldp[1][ch]-oldp[0][ch];e20[ch]=oldp[2][ch]-oldp[0][ch];
                e11[ch]=newp[1][ch]-newp[0][ch];e21[ch]=newp[2][ch]-newp[0][ch];
            }
            n0[0]=e10[1]*e20[2]-e10[2]*e20[1];
            n0[1]=e10[2]*e20[0]-e10[0]*e20[2];
            n0[2]=e10[0]*e20[1]-e10[1]*e20[0];
            n1[0]=e11[1]*e21[2]-e11[2]*e21[1];
            n1[1]=e11[2]*e21[0]-e11[0]*e21[2];
            n1[2]=e11[0]*e21[1]-e11[1]*e21[0];
            double a0=n0[0]*n0[0]+n0[1]*n0[1]+n0[2]*n0[2];
            double a1=n1[0]*n1[0]+n1[1]*n1[1]+n1[2]*n1[2];
            double dot=n0[0]*n1[0]+n0[1]*n1[1]+n0[2]*n1[2];
            if(!isfinite(a1)||a1<1e-24||(a0>=1e-24&&dot<minimum_ratio*a0))bad=1;
        }
    }
    return !bad;
}

static int metric_float_local_barrier(const float *base,float *proposal,
                                      int H,int W,double minimum_ratio,double maximum_step,
                                      double axis_fold_limit,
                                      int *out_rounds,double *out_min_scale,
                                      int *out_zero_fallbacks){
    size_t n=(size_t)H*(size_t)W;int kk;
    int TH=((H-1)+METRIC_BARRIER_TILE-1)/METRIC_BARRIER_TILE+1;
    int TW=((W-1)+METRIC_BARRIER_TILE-1)/METRIC_BARRIER_TILE+1;
    size_t tn=(size_t)TH*(size_t)TW;
    float *scale=(float *)malloc(n*sizeof *scale);
    float *control=(float *)malloc(tn*sizeof *control);
    int *bad=(int *)malloc(tn*sizeof *bad),*dilated=(int *)malloc(tn*sizeof *dilated);
    if(!scale||!control||!bad||!dilated){free(scale);free(control);free(bad);free(dilated);return -1;}
    for(size_t i=0;i<tn;i++)control[i]=1.0f;
    int rounds=0,zero_fallbacks=0;
    for(;;){
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for(kk=0;kk<(int)n;kk++){
            size_t i=(size_t)kk;int row=(int)(i/(size_t)W),column=(int)(i%(size_t)W);
            double m2=0.0;
            for(int ch=0;ch<3;ch++){double d=(double)proposal[i*3+(size_t)ch]-
                (double)base[i*3+(size_t)ch];m2+=d*d;}
            double length=sqrt(m2),s=metric_barrier_control_at(control,TH,TW,H,W,row,column);
            if(length>maximum_step)s=fmin(s,maximum_step/length);
            scale[i]=(float)s;
        }
        memset(bad,0,tn*sizeof *bad);int any=0;int r;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(|:any)
#endif
        for(r=0;r<H-1;r++)for(int c=0;c+1<W;c++){
            size_t a=(size_t)r*(size_t)W+(size_t)c,b=a+1,q=a+(size_t)W,d=q+1;
            size_t tri[2][3]={{a,b,d},{a,d,q}};
            for(int t=0;t<2;t++){
                double oldp[3][3],newp[3][3];
                for(int v=0;v<3;v++)for(int ch=0;ch<3;ch++){
                    size_t vertex=tri[t][v],k=vertex*3+(size_t)ch;
                    oldp[v][ch]=base[k];
                    newp[v][ch]=(double)base[k]+(double)scale[vertex]*
                        ((double)proposal[k]-(double)base[k]);
                }
                double e10[3],e20[3],e11[3],e21[3],n0[3],n1[3];
                for(int ch=0;ch<3;ch++){
                    e10[ch]=oldp[1][ch]-oldp[0][ch];e20[ch]=oldp[2][ch]-oldp[0][ch];
                    e11[ch]=newp[1][ch]-newp[0][ch];e21[ch]=newp[2][ch]-newp[0][ch];
                }
                n0[0]=e10[1]*e20[2]-e10[2]*e20[1];
                n0[1]=e10[2]*e20[0]-e10[0]*e20[2];
                n0[2]=e10[0]*e20[1]-e10[1]*e20[0];
                n1[0]=e11[1]*e21[2]-e11[2]*e21[1];
                n1[1]=e11[2]*e21[0]-e11[0]*e21[2];
                n1[2]=e11[0]*e21[1]-e11[1]*e21[0];
                double a0=n0[0]*n0[0]+n0[1]*n0[1]+n0[2]*n0[2];
                double a1=n1[0]*n1[0]+n1[1]*n1[1]+n1[2]*n1[2];
                double dot=n0[0]*n1[0]+n0[1]*n1[1]+n0[2]*n1[2];
                int axial_fold=0;
                if(axis_fold_limit>0.0&&axis_fold_limit<1.0&&a1>1e-24){
                    double axi0=a0>1e-24?fabs(n0[0])/sqrt(a0):0.0;
                    double axi1=fabs(n1[0])/sqrt(a1);
                    axial_fold=axi1>axis_fold_limit&&axi1>axi0+METRIC_AXIS_FOLD_EPS;
                }
                if(!isfinite(a1)||a1<1e-24||(a0>=1e-24&&dot<minimum_ratio*a0)||axial_fold){
                    any=1;
#ifdef _OPENMP
#pragma omp critical(metric_barrier_mark_float)
#endif
                    {
                        for(int v=0;v<3;v++){
                            int vr=(int)(tri[t][v]/(size_t)W),vc=(int)(tri[t][v]%(size_t)W);
                            int gr=vr/METRIC_BARRIER_TILE,gc=vc/METRIC_BARRIER_TILE;
                            if(gr>=TH-1)gr=TH-2;if(gc>=TW-1)gc=TW-2;
                            size_t node=(size_t)gr*(size_t)TW+(size_t)gc;
                            bad[node]=bad[node+1]=bad[node+(size_t)TW]=bad[node+(size_t)TW+1]=1;
                        }
                    }
                }
            }
        }
        if(!any)break;
        if(rounds>=METRIC_BARRIER_HALVING_ROUNDS+METRIC_BARRIER_ZERO_ROUNDS){
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
            for(kk=0;kk<(int)n;kk++)for(int ch=0;ch<3;ch++){
                size_t k=(size_t)kk*3+(size_t)ch;proposal[k]=base[k];
            }
            free(scale);free(control);free(bad);free(dilated);
            if(out_rounds)*out_rounds=rounds;
            if(out_min_scale)*out_min_scale=0.0;
            if(out_zero_fallbacks)*out_zero_fallbacks=zero_fallbacks+1;
            return 1;
        }
        for(size_t k=0;k<tn;k++){
            int row=(int)(k/(size_t)TW),c=(int)(k%(size_t)TW);
            int marked=0;
            for(int dr=-1;dr<=1&&!marked;dr++)for(int dc=-1;dc<=1;dc++){
                int rr=row+dr,cc=c+dc;
                if(rr>=0&&rr<TH&&cc>=0&&cc<TW&&bad[(size_t)rr*(size_t)TW+(size_t)cc]){
                    marked=1;break;
                }
            }
            dilated[k]=marked;
        }
        for(size_t k=0;k<tn;k++)if(dilated[k]){
            if(rounds>=METRIC_BARRIER_HALVING_ROUNDS)control[k]=0.0f;
            else control[k]*=0.5f;
        }
        if(rounds>=METRIC_BARRIER_HALVING_ROUNDS)zero_fallbacks++;
        rounds++;
    }
    double minimum=1.0;
    for(size_t i=0;i<n;i++)if(scale[i]<minimum)minimum=scale[i];
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for(kk=0;kk<(int)n;kk++){
        size_t i=(size_t)kk;
        for(int ch=0;ch<3;ch++){size_t k=i*3+(size_t)ch;
            proposal[k]=(float)((double)base[k]+(double)scale[i]*
                                ((double)proposal[k]-(double)base[k]));}
    }
    free(scale);free(control);free(bad);free(dilated);
    if(out_rounds)*out_rounds=rounds;if(out_min_scale)*out_min_scale=minimum;
    if(out_zero_fallbacks)*out_zero_fallbacks=zero_fallbacks;
    return 0;
}

static void metric_edge_error_float(const float *p,int H,int W,double du,double dv,
                                    double *out_rms,double *out_max){
    double ss=0.0,mx=0.0;size_t count=0;double ld=hypot(du,dv);
    for(int r=0;r<H;r++)for(int c=0;c<W;c++){
        size_t i=(size_t)r*(size_t)W+(size_t)c,js[3];double rest[3];int ne=0;
        if(c+1<W){js[ne]=i+1;rest[ne++]=du;}
        if(r+1<H){js[ne]=i+(size_t)W;rest[ne++]=dv;}
        if(r+1<H&&c+1<W){js[ne]=i+(size_t)W+1;rest[ne++]=ld;}
        for(int e=0;e<ne;e++){
            size_t j=js[e];double z=(double)p[i*3]-p[j*3];
            double y=(double)p[i*3+1]-p[j*3+1],x=(double)p[i*3+2]-p[j*3+2];
            double length=sqrt(z*z+y*y+x*x);
            if(length>0.0&&rest[e]>0.0){double err=fabs(log(length/rest[e]));
                ss+=err*err;if(err>mx)mx=err;count++;}
        }
    }
    *out_rms=count?sqrt(ss/(double)count):0.0;*out_max=mx;
}

static int metric_arap_hierarchy(float *work,const float *anchor,const float *uv,
                                 const uint8_t *filled,int H,int W,
                                 const QuadStripArapOpts *opts,
                                 QuadStripArapStats *stats){
    size_t n=(size_t)H*(size_t)W;
    double du=fabs((double)uv[2]-(double)uv[0]);
    double dv=fabs((double)uv[(size_t)W*2+1]-(double)uv[1]);
    double initial_rms=0.0,initial_max=0.0;
    metric_edge_error_float(anchor,H,W,du,dv,&initial_rms,&initial_max);
    int activate=!opts->fixed_vertices&&opts->nonlinear_multigrid_levels>1&&H>2&&W>2&&
        (opts->nonlinear_min_vertices==0||n>=(size_t)opts->nonlinear_min_vertices);
    if(!activate){
        int rc=metric_arap_core(work,anchor,uv,filled,H,W,opts,stats);
        if(rc==0&&stats){
            stats->initial_edge_log_rms=initial_rms;stats->initial_edge_log_max=initial_max;
            stats->nonlinear_min_prolongation_scale=1.0;
        }
        return rc;
    }

    MetricArapCoarse C;int brc=metric_arap_build_coarse(
        work,anchor,uv,filled,opts->repair_weight,H,W,&C);
    if(brc!=0){
        if(brc<0)return -1;
        return metric_arap_core(work,anchor,uv,filled,H,W,opts,stats);
    }
    QuadStripArapOpts co=*opts;
    co.max_iterations=opts->nonlinear_coarse_iterations;
    co.nonlinear_multigrid_levels=opts->nonlinear_multigrid_levels-1;
    co.nonlinear_min_vertices=0;
    co.repair_weight=C.repair;
    if(opts->verbose)fprintf(stderr,
        "      nonlinear MG descend %dx%d -> %dx%d (%zu vertices, %d coarse rounds)\n",
        H,W,C.H,C.W,(size_t)C.H*(size_t)C.W,co.max_iterations);
    QuadStripArapStats coarse_stats;
    int rc=metric_arap_hierarchy(C.work,C.anchor,C.uv,C.filled,C.H,C.W,&co,&coarse_stats);
    if(rc!=0){metric_arap_coarse_free(&C);return rc;}

    if(n>SIZE_MAX/(3*sizeof(float))){metric_arap_coarse_free(&C);return -1;}
    float *candidate=(float *)malloc(n*3*sizeof *candidate);
    if(!candidate){metric_arap_coarse_free(&C);return -1;}
    metric_prolong_candidate(work,candidate,H,W,&C,1.0);
    if(opts->preserve_axial)
        for(size_t i=0;i<n;i++)candidate[i*3]=anchor[i*3];
    int barrier_rounds=0,barrier_zero_fallbacks=0;double minimum_scale=1.0;
    int prolong_status=metric_float_local_barrier(
        work,candidate,H,W,opts->minimum_area_ratio,
        8.0*opts->maximum_vertex_step,opts->axis_fold_limit,&barrier_rounds,&minimum_scale,
        &barrier_zero_fallbacks);
    if(prolong_status<0){
        free(candidate);metric_arap_coarse_free(&C);return -2;
    }
    memcpy(work,candidate,n*3*sizeof *work);
    free(candidate);metric_arap_coarse_free(&C);
    if(opts->verbose)fprintf(stderr,
        "      nonlinear MG prolong %dx%d: local barrier rounds %d min-scale %.9g\n",
        H,W,barrier_rounds,minimum_scale);

    QuadStripArapStats fine_stats;
    rc=metric_arap_core(work,anchor,uv,filled,H,W,opts,&fine_stats);
    if(rc!=0)return rc;
    fine_stats.initial_edge_log_rms=initial_rms;
    fine_stats.initial_edge_log_max=initial_max;
    fine_stats.nonlinear_coarse_levels=coarse_stats.nonlinear_coarse_levels+1;
    fine_stats.nonlinear_coarse_iterations=
        coarse_stats.nonlinear_coarse_iterations+coarse_stats.iterations;
    fine_stats.nonlinear_coarse_caps=coarse_stats.nonlinear_coarse_caps+
        (coarse_stats.converged?0:1);
    fine_stats.nonlinear_min_prolongation_scale=fmin(
        minimum_scale,coarse_stats.nonlinear_min_prolongation_scale);
    fine_stats.minimum_step_scale=fmin(
        fine_stats.minimum_step_scale,fmin(minimum_scale,
                                          coarse_stats.minimum_step_scale));
    fine_stats.barrier_rejections+=coarse_stats.barrier_rejections+barrier_rounds;
    fine_stats.barrier_zero_fallbacks+=coarse_stats.barrier_zero_fallbacks+
        barrier_zero_fallbacks;
    fine_stats.barrier_stalled|=coarse_stats.barrier_stalled;
    fine_stats.multigrid_cycles+=coarse_stats.multigrid_cycles;
    if(coarse_stats.maximum_cg_iterations>fine_stats.maximum_cg_iterations)
        fine_stats.maximum_cg_iterations=coarse_stats.maximum_cg_iterations;
    if(coarse_stats.multigrid_levels_used>fine_stats.multigrid_levels_used)
        fine_stats.multigrid_levels_used=coarse_stats.multigrid_levels_used;
    if(stats)*stats=fine_stats;
    return 0;
}

int QuadStrip_metric_arap(float *verts,const float *uv,const uint8_t *filled,
                          int H,int W,const QuadStripArapOpts *opts,
                          QuadStripArapStats *stats){
    QuadStripArapOpts o;if(opts)o=*opts;else QuadStripArap_defaults(&o);
    if(!verts||!uv||!filled||H<2||W<2)return -1;
    size_t n=(size_t)H*(size_t)W;
    int activate=o.nonlinear_multigrid_levels>1&&H>2&&W>2&&
        (o.nonlinear_min_vertices==0||n>=(size_t)o.nonlinear_min_vertices);
    if(!activate)return metric_arap_core(verts,verts,uv,filled,H,W,&o,stats);
    if(n>SIZE_MAX/(3*sizeof(float)))return -1;
    float *anchor=(float *)malloc(n*3*sizeof *anchor);
    if(!anchor)return -1;
    memcpy(anchor,verts,n*3*sizeof *anchor);
    int rc=metric_arap_hierarchy(verts,anchor,uv,filled,H,W,&o,stats);
    free(anchor);return rc;
}

/* Closest rotation to a 3x3 matrix S (polar factor), via the symmetric eigen-
 * decomposition of M = S^T S:  R = S * (M^-1/2), with a determinant flip so the
 * result is a proper rotation.  Reuses the shared Eig3 solver. */
static void closest_rotation(const double S[9], double R[9]) {
    double M[3][3];
    for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++) {
            double s = 0.0;
            for (int k = 0; k < 3; k++) s += S[k * 3 + a] * S[k * 3 + b];
            M[a][b] = s;
        }
    double ev[3], Q[3][3];
    Eig3_sym(M, ev, Q);            /* ev ascending; Q columns are eigenvectors */
    double is[3];
    for (int k = 0; k < 3; k++) is[k] = ev[k] > 1e-18 ? 1.0 / sqrt(ev[k]) : 0.0;
    /* Ms = Q diag(is) Q^T */
    double Ms[3][3];
    for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++) {
            double s = 0.0;
            for (int k = 0; k < 3; k++) s += Q[a][k] * is[k] * Q[b][k];
            Ms[a][b] = s;
        }
    for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++) {
            double s = 0.0;
            for (int k = 0; k < 3; k++) s += S[a * 3 + k] * Ms[k][b];
            R[a * 3 + b] = s;
        }
    double det = R[0] * (R[4] * R[8] - R[5] * R[7])
               - R[1] * (R[3] * R[8] - R[5] * R[6])
               + R[2] * (R[3] * R[7] - R[4] * R[6]);
    if (det < 0.0) {  /* reflection -> flip the smallest-singular-value axis (k=0) */
        for (int a = 0; a < 3; a++)
            for (int b = 0; b < 3; b++) Ms[a][b] -= 2.0 * is[0] * Q[a][0] * Q[b][0];
        for (int a = 0; a < 3; a++)
            for (int b = 0; b < 3; b++) {
                double s = 0.0;
                for (int k = 0; k < 3; k++) s += S[a * 3 + k] * Ms[k][b];
                R[a * 3 + b] = s;
            }
    }
}

/* Preconditioned CG for one channel of (L + diag(w)) x = rhs, where L = D - A is
 * the graph Laplacian of `adj`.  matvec: y_i = (deg_i + w_i) x_i - sum_j x_j.
 * Jacobi preconditioner 1/(deg_i + w_i).  x carries the initial guess in/out. */
static void laplacian_cg(const int32_t *off, const int32_t *tgt,
                         const double *w, const double *rhs,
                         double *x, size_t nv, int iters) {
    double *r = (double *)malloc(nv * sizeof *r);
    double *z = (double *)malloc(nv * sizeof *z);
    double *p = (double *)malloc(nv * sizeof *p);
    double *Ap = (double *)malloc(nv * sizeof *Ap);
    double *diag = (double *)malloc(nv * sizeof *diag);
    if (!r || !z || !p || !Ap || !diag) { free(r); free(z); free(p); free(Ap); free(diag); return; }
    for (size_t i = 0; i < nv; i++) diag[i] = (double)(off[i + 1] - off[i]) + w[i];
    /* r = rhs - A_op(x) */
    for (size_t i = 0; i < nv; i++) {
        double s = diag[i] * x[i];
        for (int32_t e = off[i]; e < off[i + 1]; e++) s -= x[tgt[e]];
        r[i] = rhs[i] - s;
    }
    double rz = 0.0;
    for (size_t i = 0; i < nv; i++) { z[i] = r[i] / diag[i]; p[i] = z[i]; rz += r[i] * z[i]; }
    for (int it = 0; it < iters && rz > 1e-20; it++) {
        double pAp = 0.0;
        for (size_t i = 0; i < nv; i++) {
            double s = diag[i] * p[i];
            for (int32_t e = off[i]; e < off[i + 1]; e++) s -= p[tgt[e]];
            Ap[i] = s; pAp += p[i] * s;
        }
        if (pAp <= 0.0) break;
        double alpha = rz / pAp;
        for (size_t i = 0; i < nv; i++) { x[i] += alpha * p[i]; r[i] -= alpha * Ap[i]; }
        double rz2 = 0.0;
        for (size_t i = 0; i < nv; i++) { z[i] = r[i] / diag[i]; rz2 += r[i] * z[i]; }
        double beta = rz2 / rz;
        for (size_t i = 0; i < nv; i++) p[i] = z[i] + beta * p[i];
        rz = rz2;
    }
    free(r); free(z); free(p); free(Ap); free(diag);
}

/* ARAP developable relaxation of the FILLED cells.  The rest shape is the flat
 * UV lattice embedded as (u, v, 0) -- at true arc-length spacing it is what the
 * papyrus would be if unrolled -- so driving each cell's local frame to be
 * isometric to it unfolds the harmonic fill developably instead of letting it
 * chord across a gap through a neighbouring wrap.  The fitted rim is held fixed
 * (large anchor weight); `verts` carries the harmonic warm start in/out. */
static void arap_developable(const CSR_T adj, const float *uv, float *verts,
                             const uint8_t *filled, size_t nv,
                             int rounds, int cg_iters) {
    const int32_t *off = CSR_offset(adj), *tgt = CSR_target(adj);
    const double ANCHOR = 1e6;
    double *w  = (double *)malloc(nv * sizeof *w);
    double *pf = (double *)malloc(nv * 3 * sizeof *pf);   /* fixed fitted positions */
    double *p  = (double *)malloc(nv * 3 * sizeof *p);    /* working positions (double) */
    double *R  = (double *)malloc(nv * 9 * sizeof *R);    /* per-vertex rotation */
    double *b  = (double *)malloc(nv * 3 * sizeof *b);
    double *chan = (double *)malloc(nv * sizeof *chan);
    if (!w || !pf || !p || !R || !b || !chan) { free(w); free(pf); free(p); free(R); free(b); free(chan); return; }
    for (size_t i = 0; i < nv; i++) {
        w[i] = filled[i] ? 0.0 : ANCHOR;
        for (int c = 0; c < 3; c++) { p[i * 3 + c] = verts[i * 3 + c]; pf[i * 3 + c] = verts[i * 3 + c]; }
    }
    for (int round = 0; round < rounds; round++) {
        /* local: closest rotation of current edges to the flat-UV rest edges */
        for (size_t i = 0; i < nv; i++) {
            double S[9] = {0,0,0,0,0,0,0,0,0};
            double ui = (double)uv[i * 2 + 0], vi = (double)uv[i * 2 + 1];
            for (int32_t e = off[i]; e < off[i + 1]; e++) {
                size_t j = (size_t)tgt[e];
                double ec[3] = { p[i*3+0]-p[j*3+0], p[i*3+1]-p[j*3+1], p[i*3+2]-p[j*3+2] };
                double er[3] = { ui-(double)uv[j*2+0], vi-(double)uv[j*2+1], 0.0 };
                for (int a = 0; a < 3; a++) for (int c = 0; c < 3; c++) S[a*3+c] += ec[a]*er[c];
            }
            closest_rotation(S, &R[i * 9]);
        }
        /* global RHS: b_i = sum_j (Ri+Rj)/2 (rest_i - rest_j) */
        for (size_t i = 0; i < nv; i++) {
            double ui = (double)uv[i * 2 + 0], vi = (double)uv[i * 2 + 1];
            double bi[3] = {0,0,0};
            for (int32_t e = off[i]; e < off[i + 1]; e++) {
                size_t j = (size_t)tgt[e];
                double er[3] = { ui-(double)uv[j*2+0], vi-(double)uv[j*2+1], 0.0 };
                for (int a = 0; a < 3; a++) {
                    double rr = 0.0;
                    for (int c = 0; c < 3; c++) rr += 0.5 * (R[i*9+a*3+c] + R[j*9+a*3+c]) * er[c];
                    bi[a] += rr;
                }
            }
            b[i*3+0] = bi[0]; b[i*3+1] = bi[1]; b[i*3+2] = bi[2];
        }
        /* global solve per channel: (L + diag(w)) p = b + w * pf */
        for (int c = 0; c < 3; c++) {
            for (size_t i = 0; i < nv; i++) { chan[i] = p[i*3+c]; }
            double *rhs = (double *)malloc(nv * sizeof *rhs);
            if (!rhs) break;
            for (size_t i = 0; i < nv; i++) rhs[i] = b[i*3+c] + w[i] * pf[i*3+c];
            laplacian_cg(off, tgt, w, rhs, chan, nv, cg_iters);
            for (size_t i = 0; i < nv; i++) p[i*3+c] = chan[i];
            free(rhs);
        }
    }
    for (size_t i = 0; i < nv; i++)
        for (int c = 0; c < 3; c++) verts[i*3+c] = (float)p[i*3+c];
    free(w); free(pf); free(p); free(R); free(b); free(chan);
}

/* Proper Dirichlet harmonic fill: solve L x = 0 on the FILLED cells with the
 * fitted cells as hard boundary data, one coordinate at a time.  This reduces to
 * the fill-only system, so the residual (and its tolerance) measure the fill's
 * own convergence -- unlike a penalty/anchor solve, whose relative residual is
 * swamped by the huge anchor RHS (||w*target|| ~ 1e4 * position) and therefore
 * "converges" after one iteration, leaving every fill cell at its centroid warm
 * start (a blob in the inter-wrap void).  Jacobi(diagonal)-preconditioned CG. */
static double harmonic_uv_edge_weight(const float *uv,size_t a,size_t b) {
    double du=fabs((double)uv[a*2]-(double)uv[b*2]);
    double dv=fabs((double)uv[a*2+1]-(double)uv[b*2+1]);
    /* The triangulation contributes a bookkeeping diagonal to every quad.  It
     * is not a sample direction and has zero cotangent weight on the flat UV
     * rectangle, so exclude it from the Poisson stencil. */
    if(du>1e-6&&dv>1e-6)return 0.0;
    double d=du>dv?du:dv;
    return d>1e-12?1.0/(d*d):0.0;
}

#ifndef QUAD_STRIP_NO_TAUCS
/* TAUCS branch of the Dirichlet harmonic fill: the reduced fill-only system
 * (deg on the diagonal, -w on fill/fill edges, the observed rim folded into
 * the right-hand side) is assembled ONCE as lower-triangle COO, factored once
 * (supernodal multifrontal Cholesky) and solved for every channel in a single
 * multi-RHS call.  Returns 0 on success. */
static int harmonic_fill_taucs(const int32_t *off, const int32_t *tgt,
                               float *verts, const float *uv, size_t nfill,
                               const int32_t *g2f, const size_t *f2g,
                               const double *deg, int first_ch, int end_ch,
                               const char *const channel_name[3]) {
    size_t cap = nfill * 5u + 1u, nt = 0;
    int nch = end_ch - first_ch;
    int *rows = NULL, *cols = NULL;
    double *vals = NULL, *B = NULL, *X = NULL;
    SparseFactor_T F = NULL;
    double t0 = ves_clock_sec(), t_factor = 0.0, t_solve = 0.0;
    int rc = -1;
    if (nfill == 0) return 0;
    if (nfill > (size_t)INT_MAX / 8u || cap > (size_t)INT_MAX || nch <= 0) return -1;
    rows = (int *)malloc(cap * sizeof *rows);
    cols = (int *)malloc(cap * sizeof *cols);
    vals = (double *)malloc(cap * sizeof *vals);
    B = (double *)malloc(nfill * (size_t)nch * sizeof *B);
    X = (double *)malloc(nfill * (size_t)nch * sizeof *X);
    if (!rows || !cols || !vals || !B || !X) goto done;
    for (size_t k = 0; k < nfill; k++) {
        size_t i = f2g[k];
        rows[nt] = (int)k; cols[nt] = (int)k; vals[nt] = deg[k]; nt++;
        for (int32_t e = off[i]; e < off[i + 1]; e++) {
            int32_t j = tgt[e], kj = g2f[j];
            double w = 0.0;
            if (kj < 0 || (size_t)kj >= k) continue;   /* lower triangle, once */
            w = harmonic_uv_edge_weight(uv, i, (size_t)j);
            if (w <= 0.0) continue;                     /* quad diagonal */
            if (nt >= cap) goto done;
            rows[nt] = (int)k; cols[nt] = kj; vals[nt] = -w; nt++;
        }
    }
    for (int ch = first_ch; ch < end_ch; ch++) {
        double *b = B + (size_t)(ch - first_ch) * nfill;
        for (size_t k = 0; k < nfill; k++) {
            size_t i = f2g[k];
            double sum = 0.0;
            for (int32_t e = off[i]; e < off[i + 1]; e++) {
                int32_t j = tgt[e];
                if (g2f[j] < 0)
                    sum += harmonic_uv_edge_weight(uv, i, (size_t)j) *
                           (double)verts[(size_t)j * 3 + (size_t)ch];
            }
            b[k] = sum;
        }
    }
    memset(X, 0, nfill * (size_t)nch * sizeof *X);
    if (Sparse_factor_spd((int)nfill, (int)nt, rows, cols, vals, &F) != 0) {
        fprintf(stderr, "      harmonic fill: TAUCS factorization failed "
                        "(n=%zu nnz=%zu)\n", nfill, nt);
        goto done;
    }
    t_factor = ves_clock_sec() - t0;
    if (Sparse_factor_solve_multi(F, B, X, nch) != 0) {
        fprintf(stderr, "      harmonic fill: TAUCS solve failed\n");
        goto done;
    }
    t_solve = ves_clock_sec() - t0 - t_factor;
    for (int ch = first_ch; ch < end_ch; ch++) {
        const double *b = B + (size_t)(ch - first_ch) * nfill;
        const double *x = X + (size_t)(ch - first_ch) * nfill;
        double rn = 0.0, bn = 0.0;
        for (size_t k = 0; k < nfill; k++) {
            size_t i = f2g[k];
            double ax = deg[k] * x[k];
            for (int32_t e = off[i]; e < off[i + 1]; e++) {
                int32_t j = tgt[e];
                if (g2f[j] >= 0)
                    ax -= harmonic_uv_edge_weight(uv, i, (size_t)j) * x[g2f[j]];
            }
            rn += (b[k] - ax) * (b[k] - ax);
            bn += b[k] * b[k];
            verts[i * 3 + (size_t)ch] = (float)x[k];
        }
        fprintf(stderr,
                "      harmonic fill %s: TAUCS n=%zu nnz=%zu factor=%.2fs solve=%.2fs "
                "rel_res=%.2e\n",
                channel_name[ch], nfill, nt, t_factor, t_solve,
                bn > 1e-30 ? sqrt(rn / bn) : sqrt(rn));
    }
    rc = 0;
done:
    if (F != NULL) Sparse_factor_free(&F);
    free(rows); free(cols); free(vals); free(B); free(X);
    return rc;
}
#endif

static int harmonic_fill_channels(const int32_t *off, const int32_t *tgt,
                                  const uint8_t *filled, float *verts,
                                  const float *uv,size_t nv,
                                  int first_ch,int end_ch,
                                  const char *const channel_name[3],
                                  int max_iter, double tol, int solver) {
    size_t nfill = 0;
    for (size_t i = 0; i < nv; i++) nfill += filled[i] ? 1u : 0u;
    if (nfill == 0) return 0;
    if(first_ch<0||end_ch>3||first_ch>=end_ch)return -1;
    int32_t *g2f = (int32_t *)malloc(nv * sizeof *g2f);
    size_t  *f2g = (size_t *)malloc(nfill * sizeof *f2g);
    double  *deg = (double *)malloc(nfill * sizeof *deg);
    double  *x   = (double *)malloc(nfill * sizeof *x);
    double  *r   = (double *)malloc(nfill * sizeof *r);
    double  *z   = (double *)malloc(nfill * sizeof *z);
    double  *pp  = (double *)malloc(nfill * sizeof *pp);
    double  *Ap  = (double *)malloc(nfill * sizeof *Ap);
    double  *rhs = (double *)malloc(nfill * sizeof *rhs);
    if (!g2f||!f2g||!deg||!x||!r||!z||!pp||!Ap||!rhs) {
        free(g2f);free(f2g);free(deg);free(x);free(r);free(z);free(pp);free(Ap);free(rhs);
        return -1;
    }
    { size_t k = 0;
      for (size_t i = 0; i < nv; i++) {
          if (filled[i]) { g2f[i] = (int32_t)k; f2g[k] = i; k++; }
          else g2f[i] = -1;
      } }
    for (size_t k = 0; k < nfill; k++) {
        size_t i = f2g[k];
        deg[k]=0.0;
        for(int32_t e=off[i];e<off[i+1];e++)
            deg[k]+=harmonic_uv_edge_weight(uv,i,(size_t)tgt[e]);
        if(deg[k]<=0.0){
            fprintf(stderr,"      harmonic fill: unknown vertex has no UV-axis neighbour\n");
            free(g2f);free(f2g);free(deg);free(x);free(r);free(z);free(pp);free(Ap);free(rhs);
            return -1;
        }
    }
#ifndef QUAD_STRIP_NO_TAUCS
    if (solver == 1) {
        int trc = harmonic_fill_taucs(off, tgt, verts, uv, nfill, g2f, f2g, deg,
                                      first_ch, end_ch, channel_name);
        free(g2f);free(f2g);free(deg);free(x);free(r);free(z);free(pp);free(Ap);free(rhs);
        return trc;
    }
#else
    (void)solver;
#endif
    for (int ch = first_ch; ch < end_ch; ch++) {
        double channel_started = ves_clock_sec();
        double last_progress = channel_started;
        fprintf(stderr,
                "      harmonic fill %s: start %zu unknown cells, cap=%d tol=%.2e\n",
                channel_name[ch], nfill, max_iter, tol);
        fflush(stderr);
        double bnorm = 0.0;
        for (size_t k = 0; k < nfill; k++) {
            size_t i = f2g[k];
            double s = 0.0;
            for (int32_t e = off[i]; e < off[i + 1]; e++) {
                int32_t j = tgt[e];
                double ew=harmonic_uv_edge_weight(uv,i,(size_t)j);
                if (g2f[j] < 0) s += ew*(double)verts[(size_t)j * 3 + (size_t)ch];
            }
            rhs[k] = s;
            x[k] = (double)verts[i * 3 + (size_t)ch];  /* centroid warm start */
            bnorm += s * s;
        }
        bnorm = sqrt(bnorm); if (bnorm < 1e-30) bnorm = 1.0;
        double rz = 0.0;
        for (size_t k = 0; k < nfill; k++) {
            size_t i = f2g[k];
            double ax = deg[k] * x[k];
            for (int32_t e = off[i]; e < off[i + 1]; e++) {
                int32_t j = tgt[e];
                if (g2f[j] >= 0)
                    ax-=harmonic_uv_edge_weight(uv,i,(size_t)j)*x[g2f[j]];
            }
            r[k] = rhs[k] - ax; z[k] = r[k] / deg[k]; pp[k] = z[k]; rz += r[k] * z[k];
        }
        int it = 0; double relres = 0.0;
        for (it = 0; it < max_iter; it++) {
            double pAp = 0.0;
            for (size_t k = 0; k < nfill; k++) {
                size_t i = f2g[k];
                double ax = deg[k] * pp[k];
                for (int32_t e = off[i]; e < off[i + 1]; e++) {
                    int32_t j = tgt[e];
                    if (g2f[j] >= 0)
                        ax-=harmonic_uv_edge_weight(uv,i,(size_t)j)*pp[g2f[j]];
                }
                Ap[k] = ax; pAp += pp[k] * ax;
            }
            if (pAp <= 1e-30) break;
            double alpha = rz / pAp, rn = 0.0;
            for (size_t k = 0; k < nfill; k++) { x[k] += alpha * pp[k]; r[k] -= alpha * Ap[k]; rn += r[k] * r[k]; }
            relres = sqrt(rn) / bnorm;
            if (relres < tol) { it++; break; }
            double rz2 = 0.0;
            for (size_t k = 0; k < nfill; k++) { z[k] = r[k] / deg[k]; rz2 += r[k] * z[k]; }
            double beta = (rz != 0.0) ? rz2 / rz : 0.0;
            for (size_t k = 0; k < nfill; k++) pp[k] = z[k] + beta * pp[k];
            rz = rz2;
            {
                double now = ves_clock_sec();
                if (now - last_progress >= 10.0) {
                    fprintf(stderr,
                            "      harmonic fill %s: iteration=%d rel_res=%.2e elapsed=%.1fs\n",
                            channel_name[ch], it + 1, relres,
                            now - channel_started);
                    fflush(stderr);
                    last_progress = now;
                }
            }
        }
        fprintf(stderr,
                "      harmonic fill %s: done %d iterations, rel_res=%.2e, elapsed=%.1fs\n",
                channel_name[ch], it, relres,
                ves_clock_sec() - channel_started);
        fflush(stderr);
        for (size_t k = 0; k < nfill; k++) verts[f2g[k] * 3 + (size_t)ch] = (float)x[k];
    }
    free(g2f);free(f2g);free(deg);free(x);free(r);free(z);free(pp);free(Ap);free(rhs);
    return 0;
}

/* SPAN mode: solid = union of each row's [first..last domain] and each column's
 * [first..last domain] -- the axis-aligned span hull, which fills any hole
 * spanned in a row or a column while staying inside the component's extent. */
static void span_fill(const uint8_t *domain, int H, int W, uint8_t *solid) {
    size_t HW = (size_t)H * (size_t)W;
    memset(solid, 0, HW);
    for (int j = 0; j < H; j++) {
        int lo = -1, hi = -1;
        for (int i = 0; i < W; i++) {
            if (domain[(size_t)j * (size_t)W + (size_t)i]) {
                if (lo < 0) lo = i;
                hi = i;
            }
        }
        for (int i = lo; i >= 0 && i <= hi; i++) {
            solid[(size_t)j * (size_t)W + (size_t)i] = 1;
        }
    }
    for (int i = 0; i < W; i++) {
        int lo = -1, hi = -1;
        for (int j = 0; j < H; j++) {
            if (domain[(size_t)j * (size_t)W + (size_t)i]) {
                if (lo < 0) lo = j;
                hi = j;
            }
        }
        for (int j = lo; j >= 0 && j <= hi; j++) {
            solid[(size_t)j * (size_t)W + (size_t)i] = 1;
        }
    }
}

/* Production fill mask.  A zero-domain component is loftable only when:
 *
 *   1. it does not touch the lattice boundary (so it is a hole, not exterior),
 *   2. every fitted rim cell belongs to the same 4-connected fitted island, and
 *   3. its deepest cell is within max_hole_distance of that rim.
 *
 * This is the topological distinction RECT/SPAN omitted: an open absence is not
 * evidence for a physical surface.  Keeping it out of `solid` makes the raw bake
 * transparent there instead of sampling a harmonic extrapolation through air. */
static void bounded_fill(Arena_T arena, const uint8_t *domain, int H, int W,
                         int max_hole_distance, uint8_t *solid) {
    size_t HW = (size_t)H * (size_t)W;
    int32_t *fit_cc = ARENA_ALLOC(arena, HW * sizeof *fit_cc);
    int32_t *seen_zero = ARENA_ALLOC(arena, HW * sizeof *seen_zero);
    int32_t *dist = ARENA_ALLOC(arena, HW * sizeof *dist);
    int32_t *queue = ARENA_ALLOC(arena, HW * sizeof *queue);
    int32_t *frontier = ARENA_ALLOC(arena, HW * sizeof *frontier);
    int fit_components = 0;
    size_t accepted_components = 0, accepted_cells = 0;
    size_t open_components = 0, mixed_components = 0, wide_components = 0;

    memcpy(solid, domain, HW);
    for (size_t k = 0; k < HW; k++) {
        fit_cc[k] = -1;
        seen_zero[k] = 0;
        dist[k] = -1;
    }

    /* Label fitted islands first so a geometrically enclosed void between two
     * unrelated charts cannot masquerade as a same-sheet hole. */
    for (size_t seed = 0; seed < HW; seed++) {
        if (!domain[seed] || fit_cc[seed] >= 0) continue;
        size_t qh = 0, qt = 0;
        fit_cc[seed] = fit_components;
        queue[qt++] = (int32_t)seed;
        while (qh < qt) {
            size_t k = (size_t)queue[qh++];
            int r = (int)(k / (size_t)W), c = (int)(k % (size_t)W);
            const int dr[4] = {-1, 1, 0, 0};
            const int dc[4] = {0, 0, -1, 1};
            for (int d = 0; d < 4; d++) {
                int rr = r + dr[d], cc = c + dc[d];
                if (rr < 0 || rr >= H || cc < 0 || cc >= W) continue;
                size_t j = (size_t)rr * (size_t)W + (size_t)cc;
                if (domain[j] && fit_cc[j] < 0) {
                    fit_cc[j] = fit_components;
                    queue[qt++] = (int32_t)j;
                }
            }
        }
        fit_components++;
    }

    for (size_t seed = 0; seed < HW; seed++) {
        if (domain[seed] || seen_zero[seed]) continue;
        size_t qh = 0, qt = 0;
        int touches_outer = 0, rim_component = -1, mixed_rim = 0;
        size_t rim_edges = 0;
        seen_zero[seed] = 1;
        queue[qt++] = (int32_t)seed;
        while (qh < qt) {
            size_t k = (size_t)queue[qh++];
            int r = (int)(k / (size_t)W), c = (int)(k % (size_t)W);
            if (r == 0 || r + 1 == H || c == 0 || c + 1 == W) touches_outer = 1;
            const int dr[4] = {-1, 1, 0, 0};
            const int dc[4] = {0, 0, -1, 1};
            for (int d = 0; d < 4; d++) {
                int rr = r + dr[d], cc = c + dc[d];
                if (rr < 0 || rr >= H || cc < 0 || cc >= W) continue;
                size_t j = (size_t)rr * (size_t)W + (size_t)cc;
                if (domain[j]) {
                    int fc = fit_cc[j];
                    rim_edges++;
                    if (rim_component < 0) rim_component = fc;
                    else if (fc != rim_component) mixed_rim = 1;
                } else if (!seen_zero[j]) {
                    seen_zero[j] = 1;
                    queue[qt++] = (int32_t)j;
                }
            }
        }

        if (touches_outer || rim_edges < 4 || rim_component < 0) {
            open_components++;
            continue;
        }
        if (mixed_rim) {
            mixed_components++;
            continue;
        }

        /* Multi-source distance from the fitted rim through this zero component. */
        size_t fh = 0, ft = 0;
        for (size_t q = 0; q < qt; q++) {
            size_t k = (size_t)queue[q];
            dist[k] = -1;
            int r = (int)(k / (size_t)W), c = (int)(k % (size_t)W);
            int on_rim = 0;
            if (r > 0 && domain[k-(size_t)W]) on_rim = 1;
            if (r + 1 < H && domain[k+(size_t)W]) on_rim = 1;
            if (c > 0 && domain[k-1]) on_rim = 1;
            if (c + 1 < W && domain[k+1]) on_rim = 1;
            if (on_rim) { dist[k] = 1; frontier[ft++] = (int32_t)k; }
        }
        int deepest = 0;
        while (fh < ft) {
            size_t k = (size_t)frontier[fh++];
            int dk = dist[k];
            if (dk > deepest) deepest = dk;
            int r = (int)(k / (size_t)W), c = (int)(k % (size_t)W);
            const int dr[4] = {-1, 1, 0, 0};
            const int dc[4] = {0, 0, -1, 1};
            for (int d = 0; d < 4; d++) {
                int rr = r + dr[d], cc = c + dc[d];
                if (rr < 0 || rr >= H || cc < 0 || cc >= W) continue;
                size_t j = (size_t)rr * (size_t)W + (size_t)cc;
                if (!domain[j] && seen_zero[j] && dist[j] < 0) {
                    dist[j] = dk + 1;
                    frontier[ft++] = (int32_t)j;
                }
            }
        }
        if (max_hole_distance <= 0 || deepest > max_hole_distance) {
            wide_components++;
            continue;
        }
        for (size_t q = 0; q < qt; q++) solid[(size_t)queue[q]] = 1;
        accepted_components++;
        accepted_cells += qt;
    }

    fprintf(stderr,
            "      bounded topology: fitted islands=%d, closed holes=%zu/%zu cells; "
            "alpha components open/mixed/wide=%zu/%zu/%zu (max depth %d)\n",
            fit_components, accepted_components, accepted_cells,
            open_components, mixed_components, wide_components,
            max_hole_distance);
}

/* Vertices which participate in no emitted quad cannot affect the bake and make
 * the face-CSR relaxer singular.  Remove them from sparse production topology. */
static void retain_renderable_vertices(Arena_T arena, int H, int W,
                                       int max_row_gap, uint8_t *solid) {
    size_t HW = (size_t)H * (size_t)W;
    uint8_t *used = ARENA_CALLOC(arena, HW, 1);
    uint8_t *covered = ARENA_CALLOC(arena, HW, 1);
    uint8_t *candidate = ARENA_ALLOC(arena, W > 1 ? (size_t)(W - 1) : 1u);
    for (int gap = 1; gap <= max_row_gap; gap++) {
        for (int r = 0; r + gap < H; r++) {
            /* Decide the whole row at once.  This matches NumPy's vectorized
             * tracks_to_obj implementation: adjacent quads at the same gap
             * share vertices and must not suppress each other while marking. */
            for (int c = 0; c + 1 < W; c++) {
                size_t a = (size_t)r*(size_t)W+(size_t)c;
                size_t b = a+1;
                size_t cc = (size_t)(r+gap)*(size_t)W+(size_t)c;
                size_t d = cc+1;
                candidate[c] = (uint8_t)(solid[a] && solid[b] && solid[cc] &&
                                          solid[d] && !covered[a] && !covered[b]);
            }
            for (int c = 0; c + 1 < W; c++) {
                if (!candidate[c]) continue;
                size_t a = (size_t)r*(size_t)W+(size_t)c;
                size_t b = a+1;
                size_t cc = (size_t)(r+gap)*(size_t)W+(size_t)c;
                size_t d = cc+1;
                used[a] = used[b] = used[cc] = used[d] = 1;
                covered[a] = covered[b] = 1;
            }
        }
    }
    memcpy(solid, used, HW);
}

/* Match tracks_to_obj's sparse-row topology: a horizontal vertex pair joins
 * the nearest later row containing the matching pair.  Once a nearer quad
 * owns either top vertex, a larger row gap may not bridge across it. */
static size_t emit_ragged_faces(Arena_T arena, const int32_t *idx,
                                int H, int W, int max_row_gap,
                                int32_t **out_faces) {
    size_t max_quads = H > 1 && W > 1
                     ? (size_t)(H - 1) * (size_t)(W - 1) : 0;
    int32_t *faces = ARENA_ALLOC(
        arena, (max_quads ? max_quads * 2u : 1u) * 3u * sizeof *faces);
    uint8_t *covered = ARENA_CALLOC(arena, (size_t)H * (size_t)W, 1);
    uint8_t *candidate = ARENA_ALLOC(arena, W > 1 ? (size_t)(W - 1) : 1u);
    size_t nf = 0;
    for (int gap = 1; gap <= max_row_gap; gap++) {
        for (int r = 0; r + gap < H; r++) {
            for (int c = 0; c + 1 < W; c++) {
                size_t ka = (size_t)r*(size_t)W+(size_t)c;
                size_t kb = ka+1;
                size_t kc = (size_t)(r+gap)*(size_t)W+(size_t)c;
                size_t kd = kc+1;
                candidate[c] = (uint8_t)(idx[ka] >= 0 && idx[kb] >= 0 &&
                                          idx[kc] >= 0 && idx[kd] >= 0 &&
                                          !covered[ka] && !covered[kb]);
            }
            for (int c = 0; c + 1 < W; c++) {
                if (!candidate[c]) continue;
                size_t ka = (size_t)r*(size_t)W+(size_t)c;
                size_t kb = ka+1;
                size_t kc = (size_t)(r+gap)*(size_t)W+(size_t)c;
                size_t kd = kc+1;
                int32_t a=idx[ka],b=idx[kb],cc=idx[kc],d=idx[kd];
                faces[nf*3]=a;faces[nf*3+1]=b;faces[nf*3+2]=d;nf++;
                faces[nf*3]=a;faces[nf*3+1]=d;faces[nf*3+2]=cc;nf++;
                covered[ka]=covered[kb]=1;
            }
        }
    }
    *out_faces=faces;
    return nf;
}

typedef struct PhaseIsland {
    int id;
    size_t n;
    double mean_col, mean_radius, mean_theta;
} PhaseIsland;

static int phase_island_col_cmp(const void *aa, const void *bb) {
    const PhaseIsland *a = (const PhaseIsland *)aa;
    const PhaseIsland *b = (const PhaseIsland *)bb;
    return a->mean_col < b->mean_col ? -1 : a->mean_col > b->mean_col ? 1 : 0;
}

/* Unwrap the fitted scroll phase over each 4-connected domain island.  Adjacent
 * fitted cells are samples of one sheet, so their physical phase difference is
 * far below pi even across the atan2 branch cut.  Disconnected islands are
 * ordered in developed-u and shifted by whole turns toward
 *
 *     delta(theta) ~= sign(d theta / d u) * delta(u) / mean(radius).
 *
 * This keeps the winding gauge across missing stretches and is essential for
 * the experimental giant ribbon: merely shifting every island to the same mean
 * aliases different wraps onto one polar angle. */
static double *unwrap_cylindrical_phase(Arena_T arena, const uint8_t *domain,
                                         const float *field, int H, int W,
                                         double axis_y, double axis_x,
                                         double grid_du,
                                         int *out_islands) {
    const double TWO_PI = QUAD_TWO_PI;
    size_t HW = (size_t)H * (size_t)W;
    double *theta = ARENA_ALLOC(arena, HW * sizeof *theta);
    int32_t *island_id = ARENA_ALLOC(arena, HW * sizeof *island_id);
    int32_t *queue = ARENA_ALLOC(arena, HW * sizeof *queue);
    for (size_t c = 0; c < HW; c++) island_id[c] = -1;
    PhaseIsland *stats = NULL;
    size_t stats_cap = 0;
    int islands = 0;

    for (size_t seed = 0; seed < HW; seed++) {
        if (!domain[seed] || island_id[seed] >= 0) continue;
        size_t qh = 0, qt = 0;
        double dy = (double)field[seed*3+1] - axis_y;
        double dx = (double)field[seed*3+2] - axis_x;
        theta[seed] = atan2(dx, dy);
        island_id[seed] = islands;
        queue[qt++] = (int32_t)seed;
        while (qh < qt) {
            size_t c = (size_t)queue[qh++];
            int r = (int)(c / (size_t)W), col = (int)(c % (size_t)W);
            size_t nb[4]; int nn = 0;
            if (r > 0) nb[nn++] = c - (size_t)W;
            if (r + 1 < H) nb[nn++] = c + (size_t)W;
            if (col > 0) nb[nn++] = c - 1;
            if (col + 1 < W) nb[nn++] = c + 1;
            for (int k = 0; k < nn; k++) {
                size_t j = nb[k];
                if (!domain[j] || island_id[j] >= 0) continue;
                double jy = (double)field[j*3+1] - axis_y;
                double jx = (double)field[j*3+2] - axis_x;
                double raw = atan2(jx, jy);
                raw += nearbyint((theta[c] - raw) / TWO_PI) * TWO_PI;
                theta[j] = raw;
                island_id[j] = islands;
                queue[qt++] = (int32_t)j;
            }
        }

        if ((size_t)islands == stats_cap) {
            stats_cap = stats_cap ? stats_cap * 2 : 16;
            PhaseIsland *grown = (PhaseIsland *)realloc(stats, stats_cap * sizeof *stats);
            if (!grown) { free(stats); return theta; }
            stats = grown;
        }
        double sum_col = 0.0, sum_r = 0.0, sum_theta = 0.0;
        for (size_t k = 0; k < qt; k++) {
            size_t c = (size_t)queue[k];
            double cy = (double)field[c*3+1] - axis_y;
            double cx = (double)field[c*3+2] - axis_x;
            sum_col += (double)(c % (size_t)W);
            sum_r += sqrt(cy*cy + cx*cx);
            sum_theta += theta[c];
        }
        stats[islands].id = islands;
        stats[islands].n = qt;
        stats[islands].mean_col = sum_col / (double)qt;
        stats[islands].mean_radius = sum_r / (double)qt;
        stats[islands].mean_theta = sum_theta / (double)qt;
        islands++;
    }

    if (islands > 1 && stats) {
        /* Infer winding orientation only from within-island horizontal edges;
         * this is immune to the unknown whole-turn shifts between islands. */
        double slope_sum = 0.0;
        size_t slope_n = 0;
        for (int r = 0; r < H; r++) for (int c = 0; c + 1 < W; c++) {
            size_t a = (size_t)r*(size_t)W+(size_t)c, b = a+1;
            if (domain[a] && domain[b] && island_id[a] == island_id[b]) {
                slope_sum += theta[b] - theta[a];
                slope_n++;
            }
        }
        double orient = (slope_n == 0 || slope_sum >= 0.0) ? 1.0 : -1.0;
        qsort(stats, (size_t)islands, sizeof *stats, phase_island_col_cmp);
        double *shift = (double *)calloc((size_t)islands, sizeof *shift);
        if (shift) {
            double prev_theta = stats[0].mean_theta;
            for (int k = 1; k < islands; k++) {
                double dc = stats[k].mean_col - stats[k-1].mean_col;
                double mean_r = 0.5*(stats[k-1].mean_radius + stats[k].mean_radius);
                if (mean_r < 1.0) mean_r = 1.0;
                double expected = prev_theta + orient*dc*grid_du/mean_r;
                double sh = nearbyint((expected - stats[k].mean_theta) / TWO_PI) * TWO_PI;
                shift[stats[k].id] = sh;
                prev_theta = stats[k].mean_theta + sh;
            }
            for (size_t c = 0; c < HW; c++)
                if (domain[c]) theta[c] += shift[island_id[c]];
            free(shift);
        }
    }
    free(stats);
    if (out_islands) *out_islands = islands;
    return theta;
}

/* Preserve StrokeStrip's universal-cover decision without treating its phase
 * estimate as a geometric target.  atan2(current XYZ) is the observation; the
 * reference contributes only the nearest integer multiple of 2*pi.  Thus this
 * operation cannot move a measured point, while disconnected support islands
 * retain the upstream cycle-consistent winding gauge. */
static int lift_cylindrical_phase_from_reference(
        Arena_T arena,const uint8_t *domain,const float *field,
        const float *reference_phase,int H,int W,double axis_y,double axis_x,
        double **out_theta) {
    const double TWO_PI=QUAD_TWO_PI;
    size_t HW=(size_t)H*(size_t)W,n=0,large_edges=0,edges=0;
    double *theta=ARENA_ALLOC(arena,HW*sizeof *theta);
    double sum2=0.0,max_mismatch=0.0,min_branch=1e300,max_branch=-1e300;
    if(!reference_phase)return -1;
    for(size_t c=0;c<HW;c++) {
        theta[c]=NAN;
        if(!domain[c])continue;
        double ref=(double)reference_phase[c];
        double dy=(double)field[c*3+1]-axis_y;
        double dx=(double)field[c*3+2]-axis_x;
        if(!isfinite(ref)||!isfinite(dy)||!isfinite(dx)) {
            fprintf(stderr,"      cylindrical phase: non-finite reference/XYZ at cell %zu\n",c);
            return -1;
        }
        double raw=atan2(dx,dy);
        double branch=nearbyint((ref-raw)/TWO_PI);
        double lifted=raw+branch*TWO_PI;
        double mismatch=fabs(lifted-ref);
        theta[c]=lifted;
        sum2+=mismatch*mismatch;
        if(mismatch>max_mismatch)max_mismatch=mismatch;
        if(branch<min_branch)min_branch=branch;
        if(branch>max_branch)max_branch=branch;
        n++;
    }
    if(n==0)return -1;
    for(int r=0;r<H;r++)for(int c=0;c<W;c++) {
        size_t a=(size_t)r*(size_t)W+(size_t)c;
        if(!domain[a])continue;
        if(c+1<W&&domain[a+1]) {
            edges++;
            if(fabs(theta[a+1]-theta[a])>QUAD_PI)large_edges++;
        }
        if(r+1<H&&domain[a+(size_t)W]) {
            edges++;
            if(fabs(theta[a+(size_t)W]-theta[a])>QUAD_PI)large_edges++;
        }
    }
    fprintf(stderr,
            "      cylindrical phase branch: upstream sidecar, samples=%zu "
            "integer turns=[%.0f,%.0f], reference mismatch rms/max=%.4f/%.4f rad, "
            "local edges >pi=%zu/%zu\n",
            n,min_branch,max_branch,sqrt(sum2/(double)n),max_mismatch,
            large_edges,edges);
    if(large_edges)
        fprintf(stderr,
                "      warning: lifted phase has %zu locally discontinuous measured "
                "edges; initializer quarantine should remove their geometry\n",
                large_edges);
    *out_theta=theta;
    return 0;
}

static int compare_double(const void *aa,const void *bb) {
    double a=*(const double *)aa,b=*(const double *)bb;
    return a<b?-1:a>b?1:0;
}

typedef struct CylPitchPoint {
    double phase;
    double radius;
} CylPitchPoint;

typedef struct CylPitchSample {
    double row;
    double phase;
    double pitch;
} CylPitchSample;

/* Positive inter-turn pitch sampled on a small (row,lifted-phase) lattice.
 * `pitch` is voxels per complete turn.  `integral` stores
 * integral pitch/(2*pi) dphi from phase_min at every phase node, so evaluating
 * the spiral trend is cheap and its phase derivative is positive everywhere. */
typedef struct CylPitchModel {
    int valid;
    int nr, np;
    double row_step;
    double phase_min, phase_step;
    double p05, p50, p95;
    double *pitch;
    double *integral;
} CylPitchModel;

static int cyl_pitch_point_cmp(const void *aa,const void *bb) {
    const CylPitchPoint *a=(const CylPitchPoint *)aa;
    const CylPitchPoint *b=(const CylPitchPoint *)bb;
    return a->phase<b->phase?-1:a->phase>b->phase?1:0;
}

static size_t cyl_pitch_lower_bound(const CylPitchPoint *p,size_t n,double value) {
    size_t lo=0,hi=n;
    while(lo<hi) {
        size_t mid=lo+(hi-lo)/2;
        if(p[mid].phase<value)lo=mid+1;else hi=mid;
    }
    return lo;
}

static double cyl_pitch_quantile(const double *v,size_t n,double q) {
    if(n==0)return NAN;
    double x=q*(double)(n-1);
    size_t a=(size_t)floor(x),b=a+1<n?a+1:a;
    double t=x-(double)a;
    return v[a]*(1.0-t)+v[b]*t;
}

static double cyl_pitch_row_integral(const CylPitchModel *m,int ir,double phase) {
    const double inv_two_pi=1.0/QUAD_TWO_PI;
    const double *p=&m->pitch[(size_t)ir*(size_t)m->np];
    const double *a=&m->integral[(size_t)ir*(size_t)m->np];
    double phase_max=m->phase_min+m->phase_step*(double)(m->np-1);
    if(phase<=m->phase_min)
        return (phase-m->phase_min)*p[0]*inv_two_pi;
    if(phase>=phase_max)
        return a[m->np-1]+(phase-phase_max)*p[m->np-1]*inv_two_pi;
    double x=(phase-m->phase_min)/m->phase_step;
    int j=(int)floor(x);
    if(j<0)j=0;if(j>=m->np-1)j=m->np-2;
    double t=x-(double)j;
    /* Exact integral of the linearly interpolated positive pitch. */
    return a[j]+m->phase_step*inv_two_pi*
        (p[j]*t+0.5*(p[j+1]-p[j])*t*t);
}

static double cyl_pitch_trend(const CylPitchModel *m,double row,double phase) {
    if(!m||!m->valid||m->nr<1||m->np<2)return 0.0;
    if(m->nr==1)return cyl_pitch_row_integral(m,0,phase);
    double x=row/m->row_step;
    if(x<0.0)x=0.0;
    if(x>(double)(m->nr-1))x=(double)(m->nr-1);
    int i=(int)floor(x);
    if(i>=m->nr-1)i=m->nr-2;
    double t=x-(double)i;
    return (1.0-t)*cyl_pitch_row_integral(m,i,phase)
          +t*cyl_pitch_row_integral(m,i+1,phase);
}

/* Remove the scroll's LOCAL inter-turn radial pitch before the Poisson solve. Raw
 * radius is not close to periodic: advancing lifted phase by 2*pi advances to
 * the next physical wrap.  Its large missing regions can therefore harmonically
 * collapse two turns onto one another.  The reduced coordinate
 *
 *     rho = radius - integral P(z,phi)/(2*pi) dphi
 *
 * describes the repeating cross-section instead. P is a smooth positive field
 * fitted from exact same-row phi -> phi+2*pi observations. Reintroducing its
 * integral through the already-filled phase restores the spiral while allowing
 * pitch to vary both axially and from turn to turn.
 *
 * Local estimates are windowed medians, clamped only to the observed robust
 * p05..p95 range and smoothed in log space.  Positivity is therefore structural,
 * not a post-hoc collision projection.  The remaining reduced-radius trend is
 * only affine in z; a column term would put a second winding trend back in. */
static int fit_local_pitch_model(
        Arena_T arena,
        const uint8_t *domain,const float *field,const double *theta,
        int H,int W,double axis_y,double axis_x,
        CylPitchModel *model,double coef[3]) {
    const double MATCH_TOL=0.05;
    const size_t MIN_PAIRS=32;
    size_t nfit=0;
    memset(model,0,sizeof *model);
    for(size_t c=0;c<(size_t)H*(size_t)W;c++)nfit+=domain[c]!=0;
    CylPitchPoint *row=(CylPitchPoint *)malloc((size_t)W*sizeof *row);
    CylPitchSample *sample=(CylPitchSample *)malloc((nfit?nfit:1)*sizeof *sample);
    double *value=(double *)malloc((nfit?nfit:1)*sizeof *value);
    double *local=(double *)malloc((nfit?nfit:1)*sizeof *local);
    if(!row||!sample||!value||!local){
        free(row);free(sample);free(value);free(local);return -1;
    }
    size_t ns=0;
    for(int r=0;r<H;r++) {
        size_t n=0;
        for(int c=0;c<W;c++) {
            size_t k=(size_t)r*(size_t)W+(size_t)c;
            if(!domain[k])continue;
            double dy=(double)field[k*3+1]-axis_y;
            double dx=(double)field[k*3+2]-axis_x;
            row[n].phase=theta[k];
            row[n].radius=sqrt(dy*dy+dx*dx);
            n++;
        }
        if(n<2)continue;
        qsort(row,n,sizeof *row,cyl_pitch_point_cmp);
        for(size_t i=0;i<n;i++) {
            double target=row[i].phase+QUAD_TWO_PI;
            size_t j=cyl_pitch_lower_bound(row,n,target),best=n;
            double error=1e300;
            if(j<n){best=j;error=fabs(row[j].phase-target);}
            if(j>0&&fabs(row[j-1].phase-target)<error) {
                best=j-1;error=fabs(row[j-1].phase-target);
            }
            if(best<n&&best!=i&&error<=MATCH_TOL) {
                sample[ns].row=(double)r;
                sample[ns].phase=0.5*(row[i].phase+row[best].phase);
                sample[ns].pitch=row[best].radius-row[i].radius;
                ns++;
            }
        }
    }
    free(row);
    if(ns<MIN_PAIRS) {
        fprintf(stderr,
                "      cylindrical reduced-radius: only %zu trusted phi+2pi "
                "pair%s; retaining legacy radius trend\n",ns,ns==1?"":"s");
        free(sample);free(value);free(local);return 0;
    }
    for(size_t i=0;i<ns;i++)value[i]=sample[i].pitch;
    qsort(value,ns,sizeof *value,compare_double);
    double signed_median=cyl_pitch_quantile(value,ns,0.50);
    if(!isfinite(signed_median)||fabs(signed_median)<1e-3) {
        fprintf(stderr,
                "      cylindrical reduced-radius: no stable signed turn pitch "
                "(median %.6g); retaining legacy radius trend\n",signed_median);
        free(sample);free(value);free(local);return 0;
    }
    double sign=signed_median>0.0?1.0:-1.0;
    size_t agree=0;
    for(size_t i=0;i<ns;i++)if(sign*sample[i].pitch>0.0) {
        sample[agree]=sample[i];sample[agree].pitch=sign*sample[i].pitch;agree++;
    }
    if(agree*10<ns*9) {
        fprintf(stderr,
                "      cylindrical reduced-radius: refusing mixed radial order "
                "(%zu/%zu trusted pairs agree)\n",agree,ns);
        free(sample);free(value);free(local);return -1;
    }
    double phase_min=1e300,phase_max=-1e300;
    for(size_t i=0;i<agree;i++) {
        value[i]=sample[i].pitch;
        if(sample[i].phase<phase_min)phase_min=sample[i].phase;
        if(sample[i].phase>phase_max)phase_max=sample[i].phase;
    }
    qsort(value,agree,sizeof *value,compare_double);
    double p05=cyl_pitch_quantile(value,agree,0.05);
    double p50=cyl_pitch_quantile(value,agree,0.50);
    double p95=cyl_pitch_quantile(value,agree,0.95);
    if(!(p05>0.0)||!isfinite(p95)||p95<p05) {
        free(sample);free(value);free(local);return -1;
    }

    int nr=H>1?1+(H-1+31)/32:1;
    if(nr<2&&H>1)nr=2;if(nr>17)nr=17;
    double row_step=nr>1?(double)(H-1)/(double)(nr-1):1.0;
    const double desired_phase_step=0.5*QUAD_PI;
    double phase_span=phase_max-phase_min;
    if(phase_span<desired_phase_step)phase_span=desired_phase_step;
    phase_min=floor(phase_min/desired_phase_step)*desired_phase_step;
    phase_max=ceil(phase_max/desired_phase_step)*desired_phase_step;
    if(phase_max<=phase_min)phase_max=phase_min+desired_phase_step;
    int np=1+(int)ceil((phase_max-phase_min)/desired_phase_step);
    if(np<2)np=2;if(np>65)np=65;
    double phase_step=(phase_max-phase_min)/(double)(np-1);
    size_t nodes=(size_t)nr*(size_t)np;
    model->pitch=ARENA_ALLOC(arena,nodes*sizeof *model->pitch);
    model->integral=ARENA_ALLOC(arena,nodes*sizeof *model->integral);
    double *tmp=(double *)malloc(nodes*sizeof *tmp);
    double *support=(double *)malloc(nodes*sizeof *support);
    if(!tmp||!support){
        free(sample);free(value);free(local);free(tmp);free(support);return -1;
    }
    model->nr=nr;model->np=np;model->row_step=row_step;
    model->phase_min=phase_min;model->phase_step=phase_step;
    model->p05=p05;model->p50=p50;model->p95=p95;

    for(int ir=0;ir<nr;ir++)for(int ip=0;ip<np;ip++) {
        double zr=row_step*(double)ir;
        double ph=phase_min+phase_step*(double)ip;
        double zr_window=fmax(48.0,1.75*row_step);
        double ph_window=fmax(QUAD_PI,1.75*phase_step);
        size_t nl=0;
        for(int widen=0;widen<3&&nl<MIN_PAIRS;widen++) {
            nl=0;
            double zw=zr_window*(double)(1<<widen);
            double pw=ph_window*(double)(1<<widen);
            for(size_t i=0;i<agree;i++)
                if(fabs(sample[i].row-zr)<=zw&&fabs(sample[i].phase-ph)<=pw)
                    local[nl++]=sample[i].pitch;
        }
        double estimate=p50;
        if(nl) {
            qsort(local,nl,sizeof *local,compare_double);
            estimate=cyl_pitch_quantile(local,nl,0.50);
        }
        if(estimate<p05)estimate=p05;if(estimate>p95)estimate=p95;
        size_t k=(size_t)ir*(size_t)np+(size_t)ip;
        model->pitch[k]=log(estimate);
        support[k]=(double)(nl<256?nl:256)/64.0;
        if(support[k]<0.25)support[k]=0.25;
    }

    /* A few screened Jacobi rounds remove bin boundaries without flattening
     * genuine axial/turn variation.  Work in log pitch so every iterate stays
     * positive after exponentiation. */
    for(int sweep=0;sweep<8;sweep++) {
        for(int ir=0;ir<nr;ir++)for(int ip=0;ip<np;ip++) {
            size_t k=(size_t)ir*(size_t)np+(size_t)ip;
            double sum=support[k]*model->pitch[k],den=support[k];
            if(ir>0){sum+=model->pitch[k-(size_t)np];den+=1.0;}
            if(ir+1<nr){sum+=model->pitch[k+(size_t)np];den+=1.0;}
            if(ip>0){sum+=model->pitch[k-1];den+=1.0;}
            if(ip+1<np){sum+=model->pitch[k+1];den+=1.0;}
            tmp[k]=sum/den;
        }
        memcpy(model->pitch,tmp,nodes*sizeof *tmp);
    }
    double node_min=1e300,node_max=0.0,node_sum=0.0;
    for(size_t k=0;k<nodes;k++) {
        double p=exp(model->pitch[k]);
        if(p<p05)p=p05;if(p>p95)p=p95;
        model->pitch[k]=p;
        if(p<node_min)node_min=p;if(p>node_max)node_max=p;node_sum+=p;
    }
    for(int ir=0;ir<nr;ir++) {
        size_t base=(size_t)ir*(size_t)np;
        model->integral[base]=0.0;
        for(int ip=1;ip<np;ip++)
            model->integral[base+(size_t)ip]=model->integral[base+(size_t)ip-1]
              +0.5*(model->pitch[base+(size_t)ip-1]+model->pitch[base+(size_t)ip])
                    *phase_step/QUAD_TWO_PI;
    }
    model->valid=1;

    double mean_r=0.0,mean_q=0.0;
    size_t n=0;
    for(int r=0;r<H;r++)for(int c=0;c<W;c++) {
        size_t k=(size_t)r*(size_t)W+(size_t)c;
        if(!domain[k])continue;
        double dy=(double)field[k*3+1]-axis_y;
        double dx=(double)field[k*3+2]-axis_x;
        double q=sqrt(dy*dy+dx*dx)
                -sign*cyl_pitch_trend(model,(double)r,theta[k]);
        mean_r+=(double)r;mean_q+=q;n++;
    }
    mean_r/=(double)n;mean_q/=(double)n;
    double srr=0.0,srq=0.0;
    for(int r=0;r<H;r++)for(int c=0;c<W;c++) {
        size_t k=(size_t)r*(size_t)W+(size_t)c;
        if(!domain[k])continue;
        double dy=(double)field[k*3+1]-axis_y;
        double dx=(double)field[k*3+2]-axis_x;
        double q=sqrt(dy*dy+dx*dx)
                -sign*cyl_pitch_trend(model,(double)r,theta[k]);
        double dr=(double)r-mean_r;
        srr+=dr*dr;srq+=dr*(q-mean_q);
    }
    double br=srr>1e-20?srq/srr:0.0;
    coef[0]=mean_q-br*mean_r;
    coef[1]=br;
    coef[2]=sign;
    fprintf(stderr,
            "      cylindrical local-pitch: trusted pairs=%zu agree=%zu "
            "pitch[p05/p50/p95]=%.3f/%.3f/%.3f vox; field=%dx%d "
            "node[min/mean/max]=%.3f/%.3f/%.3f, radial-order=%+.0f, "
            "rho z-slope=%+.6f\n",
            ns,agree,p05,p50,p95,nr,np,node_min,node_sum/(double)nodes,node_max,
            sign,br);
    free(sample);free(value);free(local);free(tmp);free(support);return 1;
}

/* Cylindrical mode is explicitly axis-aligned: v determines world z exactly.
 * Infer the one additive world-z offset robustly, then reject a malformed
 * handoff instead of diffusing z or trusting a damaged fitted value. */
static int infer_axial_origin(Arena_T arena,const uint8_t *domain,
                              const float *field,int H,int W,
                              double *out_z0,double *out_max_error) {
    size_t HW=(size_t)H*(size_t)W,n=0;
    for(size_t c=0;c<HW;c++)n+=domain[c]!=0;
    if(n==0)return -1;
    double *value=ARENA_ALLOC(arena,n*sizeof *value);
    size_t q=0;
    for(int r=0;r<H;r++)for(int c=0;c<W;c++) {
        size_t k=(size_t)r*(size_t)W+(size_t)c;
        if(!domain[k])continue;
        double z=(double)field[k*3];
        if(!isfinite(z))return -1;
        value[q++]=z-(double)r;
    }
    qsort(value,n,sizeof *value,compare_double);
    double z0=n&1u?value[n/2]:0.5*(value[n/2-1]+value[n/2]);
    double max_error=0.0;
    for(int r=0;r<H;r++)for(int c=0;c<W;c++) {
        size_t k=(size_t)r*(size_t)W+(size_t)c;
        if(!domain[k])continue;
        double error=fabs((double)field[k*3]-(z0+(double)r));
        if(error>max_error)max_error=error;
    }
    fprintf(stderr,
            "      cylindrical exact-z: z(row)=%.6f+row, fitted max error=%.3e\n",
            z0,max_error);
    if(max_error>1e-3) {
        fprintf(stderr,
                "      cylindrical exact-z: refusing inconsistent fixed-v handoff "
                "(limit 1e-3 vox)\n");
        return -1;
    }
    *out_z0=z0;
    if(out_max_error)*out_max_error=max_error;
    return 0;
}

/* Fit q(row,col) = a + br*row + bc*col over the fitted cells, independently
 * for q=(z,r,theta).  A Dirichlet harmonic extension reproduces an affine field
 * inside an enclosed hole, but at a one-sided/trailing fill its natural outer
 * boundary condition collapses toward the fitted-boundary mean.  Removing the
 * fitted affine trend before the solve gives the residual a sensible zero-flux
 * condition; adding it back afterwards continues the ribbon's axial and winding
 * directions beyond the last observation instead of making a polar fan. */
static void fit_cylindrical_affine(const uint8_t *domain, const float *field,
                                   const double *theta, int H, int W,
                                   double axis_y, double axis_x,
                                   double coef[3][3]) {
    double mr = 0.0, mc = 0.0, mq[3] = {0.0, 0.0, 0.0};
    size_t n = 0;
    for (int r = 0; r < H; r++) for (int c = 0; c < W; c++) {
        size_t k = (size_t)r * (size_t)W + (size_t)c;
        if (!domain[k]) continue;
        double dy = (double)field[k*3+1] - axis_y;
        double dx = (double)field[k*3+2] - axis_x;
        double q[3] = { (double)field[k*3+0], sqrt(dy*dy + dx*dx), theta[k] };
        mr += (double)r; mc += (double)c;
        for (int ch = 0; ch < 3; ch++) mq[ch] += q[ch];
        n++;
    }
    memset(coef, 0, 9 * sizeof(double));
    if (n == 0) return;
    mr /= (double)n; mc /= (double)n;
    for (int ch = 0; ch < 3; ch++) mq[ch] /= (double)n;

    double srr = 0.0, scc = 0.0, src = 0.0;
    double srq[3] = {0.0,0.0,0.0}, scq[3] = {0.0,0.0,0.0};
    for (int r = 0; r < H; r++) for (int c = 0; c < W; c++) {
        size_t k = (size_t)r * (size_t)W + (size_t)c;
        if (!domain[k]) continue;
        double rr = (double)r - mr, cc = (double)c - mc;
        double dy = (double)field[k*3+1] - axis_y;
        double dx = (double)field[k*3+2] - axis_x;
        double q[3] = { (double)field[k*3+0], sqrt(dy*dy + dx*dx), theta[k] };
        srr += rr*rr; scc += cc*cc; src += rr*cc;
        for (int ch = 0; ch < 3; ch++) {
            srq[ch] += rr * (q[ch] - mq[ch]);
            scq[ch] += cc * (q[ch] - mq[ch]);
        }
    }
    double det = srr*scc - src*src;
    double scale = srr*scc + src*src + 1.0;
    for (int ch = 0; ch < 3; ch++) {
        double br = 0.0, bc = 0.0;
        if (fabs(det) > 1e-12 * scale) {
            br = (srq[ch]*scc - scq[ch]*src) / det;
            bc = (scq[ch]*srr - srq[ch]*src) / det;
        } else if (srr >= scc && srr > 1e-20) {
            br = srq[ch] / srr;
        } else if (scc > 1e-20) {
            bc = scq[ch] / scc;
        }
        coef[ch][0] = mq[ch] - br*mr - bc*mc;
        coef[ch][1] = br;
        coef[ch][2] = bc;
    }
}

/* Per-row trend for the cylindrical channels (1 = radius, 2 = lifted phase):
 * trend(row, col) = row_a[ch][row] + col_c[ch] * col.  The column slope is the
 * MEDIAN of the differences between adjacent observed cells in the same row,
 * so a run whose fragments sit out of angular order cannot flip it the way a
 * least-squares fit over the whole rectangle does; each row's offset comes
 * from that row's own observations, and rows without any are interpolated
 * between the nearest observed rows (held constant beyond the ends).  A
 * channel with no adjacent pair falls back to the affine slope. */
typedef struct CylRowTrend {
    double *row_a[3];     /* [H] per channel (channel 0 unused) */
    double  col_c[3];
} CylRowTrend;

static int fit_cylindrical_rows(Arena_T arena, const uint8_t *domain,
                                const float *field, const double *theta,
                                int H, int W, double axis_y, double axis_x,
                                const double affine[3][3], CylRowTrend *t) {
    size_t HW = (size_t)H * (size_t)W;
    double *diff = ARENA_ALLOC(arena, (HW > 0 ? HW : 1) * sizeof *diff);
    uint8_t *has = ARENA_CALLOC(arena, (size_t)H, 1);
    memset(t, 0, sizeof *t);
    for (int ch = 1; ch < 3; ch++) {
        size_t nd = 0;
        t->row_a[ch] = ARENA_ALLOC(arena, (size_t)H * sizeof(double));
        for (int r = 0; r < H; r++) for (int c = 0; c + 1 < W; c++) {
            size_t k = (size_t)r * (size_t)W + (size_t)c;
            double qa = 0.0, qb = 0.0, dy = 0.0, dx = 0.0;
            if (!domain[k] || !domain[k + 1]) continue;
            if (ch == 2) { qa = theta[k]; qb = theta[k + 1]; }
            else {
                dy = (double)field[k*3+1] - axis_y; dx = (double)field[k*3+2] - axis_x;
                qa = sqrt(dy*dy + dx*dx);
                dy = (double)field[(k+1)*3+1] - axis_y; dx = (double)field[(k+1)*3+2] - axis_x;
                qb = sqrt(dy*dy + dx*dx);
            }
            diff[nd++] = qb - qa;
        }
        if (ch == 1) {
            /* radius is not linear in u over several turns (eccentricity and
             * ellipticity dominate) and the spiral pitch is negligible over a
             * fill reach: continue one-sided fills at the rim radius */
            t->col_c[ch] = 0.0;
        } else if (nd > 0) {
            qsort(diff, nd, sizeof *diff, compare_double);
            t->col_c[ch] = diff[nd / 2];
        } else {
            t->col_c[ch] = affine[ch][2];
        }
        for (int r = 0; r < H; r++) {
            size_t n = 0;
            for (int c = 0; c < W; c++) {
                size_t k = (size_t)r * (size_t)W + (size_t)c;
                double q = 0.0, dy = 0.0, dx = 0.0;
                if (!domain[k]) continue;
                if (ch == 2) q = theta[k];
                else { dy = (double)field[k*3+1] - axis_y; dx = (double)field[k*3+2] - axis_x; q = sqrt(dy*dy + dx*dx); }
                diff[n++] = q - t->col_c[ch] * (double)c;   /* scratch reuse */
            }
            if (n > 0) {
                qsort(diff, n, sizeof *diff, compare_double);
                t->row_a[ch][r] = n & 1u ? diff[n / 2] : 0.5 * (diff[n / 2 - 1] + diff[n / 2]);
                if (ch == 1) has[r] = 1;
            } else {
                t->row_a[ch][r] = NAN;
            }
        }
        /* rows without observations: interpolate between neighbours */
        {
            int prev = -1;
            for (int r = 0; r < H; r++) {
                if (isfinite(t->row_a[ch][r])) { prev = r; continue; }
                {
                    int next = r + 1;
                    while (next < H && !isfinite(t->row_a[ch][next])) next++;
                    if (prev >= 0 && next < H) {
                        double w = (double)(r - prev) / (double)(next - prev);
                        t->row_a[ch][r] = t->row_a[ch][prev] * (1.0 - w) + t->row_a[ch][next] * w;
                    } else if (prev >= 0) {
                        t->row_a[ch][r] = t->row_a[ch][prev];
                    } else if (next < H) {
                        t->row_a[ch][r] = t->row_a[ch][next];
                    } else {
                        t->row_a[ch][r] = affine[ch][0] + affine[ch][1] * (double)r;
                    }
                }
            }
        }
        /* Smooth the RADIUS baseline along the rows.  Each row's value is the
         * median of that row's own observations, so it carries the sampling
         * noise of which columns happened to be observed there; the fills that
         * ride on it then step between rows.  A wrap's radius is smooth in z,
         * so a short running median removes that without moving any
         * observation -- the harmonic residual still reproduces those exactly. */
        if (ch == 1 && QS_ROW_TREND_SMOOTH > 1 && H > QS_ROW_TREND_SMOOTH) {
            int hwv = QS_ROW_TREND_SMOOTH / 2;
            double *sm = ARENA_ALLOC(arena, (size_t)H * sizeof(double));
            if (sm != NULL) {
                for (int r = 0; r < H; r++) {
                    double buf[QS_ROW_TREND_SMOOTH > 1 ? QS_ROW_TREND_SMOOTH : 2];
                    int n = 0, k = 0;
                    for (k = -hwv; k <= hwv; k++) {
                        int rr = r + k;
                        if (rr < 0 || rr >= H) continue;
                        if (!isfinite(t->row_a[ch][rr])) continue;
                        buf[n++] = t->row_a[ch][rr];
                    }
                    if (n == 0) { sm[r] = t->row_a[ch][r]; continue; }
                    qsort(buf, (size_t)n, sizeof *buf, compare_double);
                    sm[r] = n & 1 ? buf[n / 2] : 0.5 * (buf[n / 2 - 1] + buf[n / 2]);
                }
                for (int r = 0; r < H; r++) t->row_a[ch][r] = sm[r];
            }
        }
        fprintf(stderr,
                "      cylindrical per-row trend %s: column slope %.5f (%s; %zu "
                "adjacent pairs; affine gave %.5f), row offsets = medians\n",
                ch == 1 ? "radius" : "phase", t->col_c[ch],
                ch == 1 ? "held at 0" : "median adjacent step", nd, affine[ch][2]);
    }
    (void)has;
    return 0;
}

static double cyl_trend_at(const CylRowTrend *t, const double affine[3][3],
                           int ch, int j, int i) {
    if (t->row_a[ch] != NULL)
        return t->row_a[ch][j] + t->col_c[ch] * (double)i;
    return affine[ch][0] + affine[ch][1] * (double)j + affine[ch][2] * (double)i;
}

int QuadStrip_build_topology_with_phase_ex(Arena_T arena,
                    const uint8_t *topology,
                    const uint8_t *domain, const float *field,
                    const float *lifted_phase,
                    int H, int W, int gr0, int c0, double grid_du,
                    const QuadStripOpts *opts,
                    float **out_verts, size_t *out_nv,
                    int32_t **out_faces, size_t *out_nf,
                    float **out_uv, uint8_t **out_filled,
                    float **out_phase) {
    float *phase_out = NULL;
    if (domain == NULL || field == NULL || H < 1 || W < 1) return -1;
    QuadStripOpts o = *opts;
    if (o.max_row_gap < 1) return -1;
    size_t HW = (size_t)H * (size_t)W;
    if (HW > (size_t)INT32_MAX) return -1;  /* faces are int32-indexed */
    const uint8_t *topology_seed = topology ? topology : domain;

    uint8_t *solid = ARENA_ALLOC(arena, HW);
    if (o.mode == QUAD_STRIP_SPAN) {
        span_fill(topology_seed, H, W, solid);
    } else if (o.mode == QUAD_STRIP_BOUNDED) {
        bounded_fill(arena, topology_seed, H, W, o.max_hole_distance, solid);
    } else {
        memset(solid, 1, HW);  /* RECT: the full rectangle */
    }
    retain_renderable_vertices(arena, H, W, o.max_row_gap, solid);

    /* Compact vertex indices over solid cells (row-major emit order). */
    int32_t *idx = ARENA_ALLOC(arena, HW * sizeof *idx);
    size_t nv = 0;
    for (size_t c = 0; c < HW; c++) idx[c] = solid[c] ? (int32_t)nv++ : -1;
    if (nv == 0) return -1;

    float   *verts  = ARENA_ALLOC(arena, nv * 3 * sizeof *verts);
    float   *uv     = ARENA_ALLOC(arena, nv * 2 * sizeof *uv);
    uint8_t *filled = ARENA_ALLOC(arena, nv);
    if (out_phase) {
        phase_out = ARENA_ALLOC(arena, nv * sizeof *phase_out);
        for (size_t v = 0; v < nv; v++) phase_out[v] = NAN;
    }

    /* Optional winding-aware coordinates.  Cartesian harmonic interpolation
     * draws a chord through a rolled sheet; (z,r,unwrapped theta) harmonic
     * interpolation follows the same wrap around the known scroll axis. */
    double *theta = NULL;
    int phase_islands = 0;
    int phase_from_reference=0;
    int reduced_radius=0;
    CylPitchModel pitch_model={0};
    double axial_z0=0.0;
    if(o.cylindrical_fill) {
        if(lifted_phase) {
            if(lift_cylindrical_phase_from_reference(
                    arena,domain,field,lifted_phase,H,W,o.axis_y,o.axis_x,
                    &theta)!=0)return -1;
            phase_from_reference=1;
        } else {
            fprintf(stderr,
                    "      warning: no lifted-phase sidecar; using legacy local "
                    "atan2 island reconstruction\n");
            theta=unwrap_cylindrical_phase(arena,domain,field,H,W,
                                           o.axis_y,o.axis_x,grid_du,
                                           &phase_islands);
        }
        if(infer_axial_origin(arena,domain,field,H,W,&axial_z0,NULL)!=0)
            return -1;
    }

    /* Centroid of fitted solve coordinates -> a cheap warm start for ordinary
     * Cartesian fill.  Cylindrical fill instead solves the residual after an
     * affine lattice trend, which is near zero and extrapolates one-sided tails. */
    double cz = 0.0, cy = 0.0, cx = 0.0;
    double cyl_affine[3][3] = {{0}};
    CylRowTrend row_trend;
    memset(&row_trend, 0, sizeof row_trend);
    if (theta) {
        fit_cylindrical_affine(domain, field, theta, H, W,
                               o.axis_y, o.axis_x, cyl_affine);
        if (o.trend_per_row) {
            if (fit_cylindrical_rows(arena, domain, field, theta, H, W,
                                     o.axis_y, o.axis_x, cyl_affine,
                                     &row_trend) != 0)
                return -1;
            reduced_radius = 0;
        } else {
            reduced_radius=fit_local_pitch_model(
                    arena,domain,field,theta,H,W,o.axis_y,o.axis_x,
                    &pitch_model,cyl_affine[1]);
            if(reduced_radius<0)return -1;
        }
    }
    if(theta) {
        cyl_affine[0][0]=axial_z0;
        cyl_affine[0][1]=1.0;
        cyl_affine[0][2]=0.0;
    }
    size_t nfit = 0;
    for (size_t c = 0; c < HW; c++) {
        if (solid[c] && domain[c]) {
            cz += field[c * 3 + 0];
            if (theta) {
                double dy = (double)field[c*3+1] - o.axis_y;
                double dx = (double)field[c*3+2] - o.axis_x;
                cy += sqrt(dy*dy + dx*dx);  /* temporary channel 1 = radius */
                cx += theta[c];             /* temporary channel 2 = phase */
            } else {
                cy += field[c * 3 + 1];
                cx += field[c * 3 + 2];
            }
            nfit++;
        }
    }
    if (nfit > 0) { cz /= (double)nfit; cy /= (double)nfit; cx /= (double)nfit; }
    size_t nfilled = nv - nfit;

    for (int j = 0; j < H; j++) {
        for (int i = 0; i < W; i++) {
            size_t c = (size_t)j * (size_t)W + (size_t)i;
            if (!solid[c]) continue;
            size_t vi = (size_t)idx[c];
            uv[vi * 2 + 0] = (float)(((double)c0 + (double)i) * grid_du);
            uv[vi * 2 + 1] = (float)((double)gr0 + (double)j);
            if (domain[c]) {
                if (theta) {
                    double dy = (double)field[c*3+1] - o.axis_y;
                    double dx = (double)field[c*3+2] - o.axis_x;
                    double q[3] = { axial_z0+(double)j,
                                    sqrt(dy*dy + dx*dx), theta[c] };
                    for (int ch = 0; ch < 3; ch++) {
                        double trend = reduced_radius&&ch==1
                                     ? cyl_affine[1][0]
                                      +cyl_affine[1][1]*(double)j
                                      +cyl_affine[1][2]
                                       *cyl_pitch_trend(&pitch_model,(double)j,q[2])
                                     : (ch == 0
                                        ? cyl_affine[0][0]+cyl_affine[0][1]*(double)j
                                        : cyl_trend_at(&row_trend, cyl_affine, ch, j, i));
                        verts[vi*3+(size_t)ch] = (float)(q[ch] - trend);
                    }
                } else {
                    verts[vi * 3 + 0] = field[c * 3 + 0];
                    verts[vi * 3 + 1] = field[c * 3 + 1];
                    verts[vi * 3 + 2] = field[c * 3 + 2];
                }
                filled[vi] = 0;
            } else {
                verts[vi * 3 + 0] = theta ? 0.0f : (float)cz;
                verts[vi * 3 + 1] = theta ? 0.0f : (float)cy;
                verts[vi * 3 + 2] = theta ? 0.0f : (float)cx;
                filled[vi] = 1;
            }
        }
    }

    /* Two triangles (a,b,d)/(a,d,c), with sparse alternating rows joined only
     * to their nearest later live row. */
    int32_t *faces = NULL;
    size_t nf = emit_ragged_faces(
        arena, idx, H, W, o.max_row_gap, &faces);

    /* UV-metric Dirichlet Poisson fill.  Fitted cells are hard boundary data;
     * filled cells become the smoothest scalar continuation matching the rim. */
    if (nf >= 1 && nfilled > 0) {
        CSR_T adj = CSR_from_faces(arena, faces, nf, nv);
        const int32_t *off = CSR_offset(adj), *tgt = CSR_target(adj);
        static const char *cartesian_name[3]={"z","y","x"};
        static const char *cylindrical_name[3]={
            "exact-z (not solved)","radius residual","lifted-phase residual"};
        /* Proper Dirichlet UV-Poisson fill (fill-only unknowns, fitted rim
         * fixed).  Cylindrical z is prescribed analytically, so only r and the
         * already-lifted phase are unknown. */
        if(harmonic_fill_channels(off,tgt,filled,verts,uv,nv,
                theta?1:0,3,theta?cylindrical_name:cartesian_name,
                o.pde_max_iter,o.pde_tol,o.harmonic_solver)!=0)return -1;

        /* Phase B: unfold the harmonic fill developably (fitted rim fixed). */
        if (!theta && o.relax_rounds > 0)
            arap_developable(adj, uv, verts, filled, nv,
                             o.relax_rounds, o.relax_cg_iters);
    }

    if (theta) {
        /* Restore the extrapolated trend and reconstruct Cartesian y/x before
         * any later geometric/CT relaxation.  This is deliberately outside the
         * nfilled/nfaces branch so a dense or one-row cylindrical strip cannot
         * leak temporary (z,r,theta) coordinates to its caller. */
        for (int j = 0; j < H; j++) for (int i = 0; i < W; i++) {
            size_t c = (size_t)j * (size_t)W + (size_t)i;
            if (!solid[c]) continue;
            size_t v = (size_t)idx[c];
            double phase=(double)verts[v*3+2]
                        +cyl_trend_at(&row_trend, cyl_affine, 2, j, i);
            double radius=(double)verts[v*3+1]
                          +(reduced_radius
                            ?cyl_affine[1][0]+cyl_affine[1][1]*(double)j
                             +cyl_affine[1][2]
                              *cyl_pitch_trend(&pitch_model,(double)j,phase)
                            :cyl_trend_at(&row_trend, cyl_affine, 1, j, i));
            if(!isfinite(radius)||!isfinite(phase)||radius<=0.0) {
                fprintf(stderr,
                        "      cylindrical fill: invalid radius/phase at row=%d "
                        "col=%d (r=%.6g phi=%.6g)\n",j,i,radius,phase);
                return -1;
            }
            if (phase_out) phase_out[v] = (float)phase;
            verts[v*3+0]=(float)(axial_z0+(double)j);
            verts[v*3+1]=(float)(o.axis_y+radius*cos(remainder(phase,QUAD_TWO_PI)));
            verts[v*3+2]=(float)(o.axis_x+radius*sin(remainder(phase,QUAD_TWO_PI)));
        }
        if(phase_from_reference)
            fprintf(stderr,
                    "      cylindrical fill: axis=(%.1f,%.1f), upstream integer "
                    "phase branches, exact z, %s UV-Poisson r/phi\n",
                    o.axis_y,o.axis_x,reduced_radius?
                    "local-pitch-reduced-radius":"affine-residual");
        else
            fprintf(stderr,
                    "      cylindrical fill: axis=(%.1f,%.1f), %d legacy fitted "
                    "phase island%s, exact z, affine-residual UV-Poisson r/phi\n",
                    o.axis_y,o.axis_x,phase_islands,
                    phase_islands==1?"":"s");
        if (o.relax_rounds > 0)
            fprintf(stderr, "      warning: --relax-rounds ignored with cylindrical fill (strict ARAP folds large fills)\n");
    }

    /* Re-stamp fitted cells with the exact sample: the anchor weight is finite,
     * so the solve leaves them within ~1e-4; the deliverable samples the real
     * faces exactly where they exist. */
    if (o.exact_fitted) {
        for (size_t c = 0; c < HW; c++) {
            if (solid[c] && domain[c]) {
                size_t vi = (size_t)idx[c];
                verts[vi * 3 + 0] = theta
                                  ? (float)(axial_z0+(double)(c/(size_t)W))
                                  : field[c * 3 + 0];
                verts[vi * 3 + 1] = field[c * 3 + 1];
                verts[vi * 3 + 2] = field[c * 3 + 2];
            }
        }
    }

    if (out_verts)  *out_verts = verts;
    if (out_nv)     *out_nv = nv;
    if (out_faces)  *out_faces = faces;
    if (out_nf)     *out_nf = nf;
    if (out_uv)     *out_uv = uv;
    if (out_filled) *out_filled = filled;
    if (out_phase)  *out_phase = phase_out;
    return 0;
}

int QuadStrip_build_topology_with_phase(Arena_T arena,
                    const uint8_t *topology,
                    const uint8_t *domain, const float *field,
                    const float *lifted_phase,
                    int H, int W, int gr0, int c0, double grid_du,
                    const QuadStripOpts *opts,
                    float **out_verts, size_t *out_nv,
                    int32_t **out_faces, size_t *out_nf,
                    float **out_uv, uint8_t **out_filled) {
    return QuadStrip_build_topology_with_phase_ex(
        arena,topology,domain,field,lifted_phase,H,W,gr0,c0,grid_du,opts,
        out_verts,out_nv,out_faces,out_nf,out_uv,out_filled,NULL);
}

int QuadStrip_build_with_phase(Arena_T arena,
                    const uint8_t *domain, const float *field,
                    const float *lifted_phase,
                    int H, int W, int gr0, int c0, double grid_du,
                    const QuadStripOpts *opts,
                    float **out_verts, size_t *out_nv,
                    int32_t **out_faces, size_t *out_nf,
                    float **out_uv, uint8_t **out_filled) {
    return QuadStrip_build_topology_with_phase(
        arena,NULL,domain,field,lifted_phase,H,W,gr0,c0,grid_du,opts,
        out_verts,out_nv,out_faces,out_nf,out_uv,out_filled);
}

int QuadStrip_build(Arena_T arena,
                    const uint8_t *domain, const float *field,
                    int H, int W, int gr0, int c0, double grid_du,
                    const QuadStripOpts *opts,
                    float **out_verts, size_t *out_nv,
                    int32_t **out_faces, size_t *out_nf,
                    float **out_uv, uint8_t **out_filled) {
    return QuadStrip_build_with_phase(arena,domain,field,NULL,H,W,gr0,c0,
                    grid_du,opts,out_verts,out_nv,out_faces,out_nf,out_uv,
                    out_filled);
}

/* ============================================================================
 * Self-test
 * ==========================================================================*/
#define CK(c, m) do { if (!(c)) { fprintf(stderr, "  FAIL: %s\n", (m)); fails++; } \
                      else fprintf(stderr, "  ok: %s\n", (m)); } while (0)

int QuadStrip_selftest(void) {
    int fails = 0;
    Arena_T arena = Arena_new();
    fprintf(stderr, "[selftest] quad_strip\n");
    QuadStripOpts o; QuadStrip_defaults(&o);

    /* A plane field pos(j,i) = (10, j, 2i).  Harmonic (Laplacian) fill of a
     * linear field reproduces it exactly, so every filled cell must land on the
     * plane -- a clean numeric check that the PDE and the anchoring are right. */
    const int H = 7, W = 8;
    Arena_Mark mk = Arena_save(arena);
    uint8_t *domain = ARENA_ALLOC(arena, (size_t)H * W);
    float *field = ARENA_ALLOC(arena, (size_t)H * W * 3 * sizeof *field);
    for (int j = 0; j < H; j++) for (int i = 0; i < W; i++) {
        size_t c = (size_t)j * W + (size_t)i;
        domain[c] = 1;
        field[c * 3 + 0] = 10.0f;
        field[c * 3 + 1] = (float)j;
        field[c * 3 + 2] = 2.0f * (float)i;
    }
    /* Punch two interior holes (missing fitted data). */
    domain[(size_t)2 * W + 3] = 0;
    domain[(size_t)4 * W + 5] = 0;
    domain[(size_t)4 * W + 6] = 0;  /* a 2-wide hole */

    /* t1: production rectangle build with exact planar harmonic fill. */
    {
        float *verts = NULL, *uv = NULL; int32_t *faces = NULL; uint8_t *fl = NULL;
        size_t nv = 0, nf = 0;
        int rc = QuadStrip_build(arena, domain, field, H, W, 100, 5, 2.0, &o,
                                &verts, &nv, &faces, &nf, &uv, &fl);
        CK(rc == 0, "t1 build ok");
        CK(nv == (size_t)H * W, "t1 full grid: nv == H*W");
        CK(nf == (size_t)(H - 1) * (W - 1) * 2, "t1 full grid: nf == 2(H-1)(W-1)");
        /* filled count == number of punched holes */
        size_t nfl = 0; for (size_t k = 0; k < nv; k++) nfl += fl[k];
        CK(nfl == 3, "t1 three cells flagged filled");
        /* every face index in range */
        int inrange = 1;
        for (size_t k = 0; k < nf * 3; k++) if (faces[k] < 0 || (size_t)faces[k] >= nv) inrange = 0;
        CK(inrange, "t1 face indices in range");
        /* Dense output remains row-major here because every cell is retained. */
        size_t hv = (size_t)2 * W + 3;
        float hz = verts[hv * 3 + 0], hy = verts[hv * 3 + 1], hx = verts[hv * 3 + 2];
        CK(fabsf(hz - 10.0f) < 0.05f && fabsf(hy - 2.0f) < 0.05f &&
           fabsf(hx - 6.0f) < 0.05f, "t1 hole (2,3) filled onto the plane");
        /* the 2-wide hole cell (4,6) also reproduced */
        size_t hv2 = (size_t)4 * W + 6;
        CK(fabsf(verts[hv2 * 3 + 1] - 4.0f) < 0.05f &&
           fabsf(verts[hv2 * 3 + 2] - 12.0f) < 0.05f, "t1 wide hole (4,6) on the plane");
        /* a fitted cell stays exact; uv is (u=(c0+i)*du, v=gr0+j) */
        size_t fv = (size_t)3 * W + 2;
        CK(verts[fv * 3 + 2] == 4.0f, "t1 fitted cell exact");
        CK(uv[fv * 2 + 0] == (float)((5 + 2) * 2.0) && uv[fv * 2 + 1] == 103.0f,
           "t1 uv = ((c0+i)*du, gr0+j)");
    }

    /* t2: SPAN mode still solid & consistent (hull >= domain). */
    {
        QuadStripOpts os = o; os.mode = QUAD_STRIP_SPAN;
        float *verts = NULL; int32_t *faces = NULL; size_t nv = 0, nf = 0;
        int rc = QuadStrip_build(arena, domain, field, H, W, 0, 0, 2.0, &os,
                                &verts, &nv, &faces, &nf, NULL, NULL);
        CK(rc == 0 && nv > 0 && nf > 0, "t2 span build ok");
        CK(nv == (size_t)H * W, "t2 dense domain -> span hull == full grid");
    }

    /* t3: winding-aware harmonic coordinates reproduce a cylindrical strip.
     * Cartesian interpolation would chord inside radius 10 at the punched cell.
     * More importantly, affine-residual extension must continue a one-sided
     * missing tail; plain harmonic values collapse that tail to the rim mean. */
    {
        QuadStripOpts oc = o;
        oc.mode = QUAD_STRIP_RECT; /* deliberately exercise unsafe one-sided
                                      extrapolation as a diagnostic */
        oc.cylindrical_fill = 1; oc.axis_y = 0.0; oc.axis_x = 0.0;
        float *phase=ARENA_ALLOC(arena,(size_t)H*W*sizeof *phase);
        for (int j = 0; j < H; j++) for (int i = 0; i < W; i++) {
            size_t c = (size_t)j * W + (size_t)i;
            double a = -0.7 + 0.2 * (double)i;
            field[c*3+0] = (float)j;
            field[c*3+1] = (float)(10.0*cos(a));
            field[c*3+2] = (float)(10.0*sin(a));
            /* Deliberate nonzero winding gauge and fractional disagreement:
             * only the integer branch may affect the solve. */
            phase[c]=(float)(a+4.0*QUAD_TWO_PI+0.35);
        }
        float *verts = NULL; size_t nv = 0;
        int rc = QuadStrip_build_with_phase(
                arena,domain,field,phase,H,W,0,0,2.0,&oc,
                &verts,&nv,NULL,NULL,NULL,NULL);
        size_t hv = (size_t)2 * W + 3;
        double rr = sqrt((double)verts[hv*3+1]*verts[hv*3+1] +
                         (double)verts[hv*3+2]*verts[hv*3+2]);
        double aa = atan2((double)verts[hv*3+2], (double)verts[hv*3+1]);
        CK(rc == 0 && nv == (size_t)H*W, "t3 cylindrical build ok");
        CK(fabs(rr - 10.0) < 0.05, "t3 filled cell stays on cylinder radius");
        CK(fabs(aa - (-0.7 + 0.2*3.0)) < 0.02, "t3 filled phase follows winding");
        CK(fabs((double)verts[hv*3]-2.0)<1e-6,
           "t3 filled z comes exactly from v");

        phase[0]=NAN;
        CK(QuadStrip_build_with_phase(
                arena,domain,field,phase,H,W,0,0,2.0,&oc,
                NULL,NULL,NULL,NULL,NULL,NULL)==-1,
           "t3 non-finite measured phase is rejected");
        phase[0]=(float)(-0.7+4.0*QUAD_TWO_PI+0.35);

        uint8_t *tail_domain = ARENA_ALLOC(arena, (size_t)H * W);
        for (int j = 0; j < H; j++) for (int i = 0; i < W; i++)
            tail_domain[(size_t)j*W+(size_t)i] = (uint8_t)(i >= 3);
        field[((size_t)3*W)*3]=12345.0f; /* ignored: this cell is not observed */
        float *tail_verts = NULL; size_t tail_nv = 0;
        rc = QuadStrip_build_with_phase(
                arena,tail_domain,field,phase,H,W,0,0,2.0,&oc,
                &tail_verts,&tail_nv,NULL,NULL,NULL,NULL);
        size_t tv = (size_t)3 * W;  /* row 3, missing leftmost column */
        double tr = sqrt((double)tail_verts[tv*3+1]*tail_verts[tv*3+1] +
                         (double)tail_verts[tv*3+2]*tail_verts[tv*3+2]);
        double ta = atan2((double)tail_verts[tv*3+2], (double)tail_verts[tv*3+1]);
        CK(rc == 0 && tail_nv == (size_t)H*W, "t3 one-sided cylindrical tail build ok");
        CK(fabs(tr - 10.0) < 0.05, "t3 one-sided tail keeps cylinder radius");
        CK(fabs(ta - (-0.7)) < 0.02, "t3 one-sided tail extrapolates winding phase");
        CK(fabs((double)tail_verts[tv*3]-3.0)<1e-6,
           "t3 one-sided tail ignores missing field z and uses v");
    }

    /* t4: disconnected fitted islands separated by more than one full turn
     * retain their whole-turn gauge from developed-u. */
    {
        const int GH = 3, GW = 40;
        uint8_t *gd = ARENA_ALLOC(arena, (size_t)GH * GW);
        float *gf = ARENA_ALLOC(arena, (size_t)GH * GW * 3 * sizeof *gf);
        float *gp = ARENA_ALLOC(arena, (size_t)GH * GW * sizeof *gp);
        for (int j = 0; j < GH; j++) for (int i = 0; i < GW; i++) {
            size_t c = (size_t)j*GW+(size_t)i;
            double a = -1.0 + 0.2*(double)i;
            gd[c] = (uint8_t)(i < 4 || i >= 36);
            gf[c*3+0] = (float)j;
            gf[c*3+1] = (float)(10.0*cos(a));
            gf[c*3+2] = (float)(10.0*sin(a));
            gp[c]=(float)(a+8.0*QUAD_TWO_PI);
        }
        QuadStripOpts og = o; og.mode = QUAD_STRIP_RECT; og.cylindrical_fill = 1;
        float *gv = NULL; size_t gnv = 0;
        int rc = QuadStrip_build_with_phase(
                arena,gd,gf,gp,GH,GW,0,0,2.0,&og,
                &gv,&gnv,NULL,NULL,NULL,NULL);
        size_t mv = (size_t)1*GW+20;
        double ma = atan2((double)gv[mv*3+2], (double)gv[mv*3+1]);
        CK(rc == 0 && gnv == (size_t)GH*GW, "t4 disconnected winding build ok");
        CK(fabs(ma - 3.0) < 0.05, "t4 disconnected islands preserve whole-turn phase");
    }

    /* t4b: a two-turn missing region is continued in pitch-reduced radius.
     * Both queried samples are generated; advancing approximately 2*pi must
     * still advance to the next radial wrap rather than collapsing onto it. */
    {
        const int PH=8,PW=100;
        const double expected_pitch=8.0;
        uint8_t *pd=ARENA_ALLOC(arena,(size_t)PH*PW);
        float *pf=ARENA_ALLOC(arena,(size_t)PH*PW*3*sizeof *pf);
        float *pp=ARENA_ALLOC(arena,(size_t)PH*PW*sizeof *pp);
        for(int r=0;r<PH;r++)for(int c=0;c<PW;c++) {
            size_t k=(size_t)r*PW+(size_t)c;
            double a=0.16*(double)c;
            double radius=30.0+expected_pitch*a/QUAD_TWO_PI+2.0*cos(a);
            pd[k]=(uint8_t)!((c>=20&&c<=29)||(c>=59&&c<=68));
            pf[k*3]=(float)r;
            pf[k*3+1]=(float)(radius*cos(a));
            pf[k*3+2]=(float)(radius*sin(a));
            pp[k]=(float)(a+4.0*QUAD_TWO_PI);
        }
        QuadStripOpts po=o;po.mode=QUAD_STRIP_RECT;po.cylindrical_fill=1;
        float *pv=NULL;size_t pnv=0;
        int rc=QuadStrip_build_with_phase(
                arena,pd,pf,pp,PH,PW,0,0,2.0,&po,
                &pv,&pnv,NULL,NULL,NULL,NULL);
        size_t inner=(size_t)3*PW+25,outer=(size_t)3*PW+64;
        double ri=hypot((double)pv[inner*3+1],(double)pv[inner*3+2]);
        double ro=hypot((double)pv[outer*3+1],(double)pv[outer*3+2]);
        CK(rc==0&&pnv==(size_t)PH*PW,"t4b reduced-radius build ok");
        CK(ro-ri>6.0,"t4b generated adjacent turns retain radial pitch");
    }

    /* t4c: pitch is a positive LOCAL field, not one scroll-wide scalar.  This
     * synthetic spiral thickens both down the axis and on later turns.  Three
     * repeated missing bands force the continuation (rather than the observed
     * samples) to demonstrate both variations. */
    {
        const int PH=96,PW=220;
        const double phase_step=0.10,phase_k=0.12;
        uint8_t *pd=ARENA_ALLOC(arena,(size_t)PH*PW);
        float *pf=ARENA_ALLOC(arena,(size_t)PH*PW*3*sizeof *pf);
        float *pp=ARENA_ALLOC(arena,(size_t)PH*PW*sizeof *pp);
        for(int r=0;r<PH;r++)for(int c=0;c<PW;c++) {
            size_t q=(size_t)r*PW+(size_t)c;
            double a=phase_step*(double)c;
            double base=7.0+0.12*(double)r;
            double trend=base*a/QUAD_TWO_PI
                        +phase_k*a*a/(2.0*QUAD_TWO_PI);
            double radius=45.0+trend+1.2*cos(a);
            int hole=(c>=30&&c<=40)||(c>=93&&c<=103)||
                     (c>=156&&c<=166);
            pd[q]=(uint8_t)!hole;
            pf[q*3]=(float)r;
            pf[q*3+1]=(float)(radius*cos(a));
            pf[q*3+2]=(float)(radius*sin(a));
            pp[q]=(float)(a+6.0*QUAD_TWO_PI);
        }
        QuadStripOpts po=o;po.mode=QUAD_STRIP_RECT;po.cylindrical_fill=1;
        float *pv=NULL;size_t pnv=0;
        int rc=QuadStrip_build_with_phase(
                arena,pd,pf,pp,PH,PW,0,0,2.0,&po,
                &pv,&pnv,NULL,NULL,NULL,NULL);
        size_t l0=(size_t)2*PW+35,l1=(size_t)2*PW+98,l2=(size_t)2*PW+161;
        size_t h0=(size_t)87*PW+35,h1=(size_t)87*PW+98;
        double lr0=hypot((double)pv[l0*3+1],(double)pv[l0*3+2]);
        double lr1=hypot((double)pv[l1*3+1],(double)pv[l1*3+2]);
        double lr2=hypot((double)pv[l2*3+1],(double)pv[l2*3+2]);
        double hr0=hypot((double)pv[h0*3+1],(double)pv[h0*3+2]);
        double hr1=hypot((double)pv[h1*3+1],(double)pv[h1*3+2]);
        CK(rc==0&&pnv==(size_t)PH*PW,"t4c local-pitch build ok");
        CK((lr2-lr1)>(lr1-lr0)+0.25,
           "t4c later generated turn retains larger local pitch");
        CK((hr1-hr0)>(lr1-lr0)+0.75,
           "t4c lower generated row retains larger local pitch");
        CK(fabs((double)pv[l1*3]-2.0)<1e-6&&
           fabs((double)pv[h1*3]-87.0)<1e-6,
           "t4c local-pitch continuation preserves exact axial rows");
    }

    /* t5: production topology fills an enclosed same-island hole but leaves an
     * open exterior notch absent. */
    {
        const int BH = 9, BW = 12;
        uint8_t *bd = ARENA_ALLOC(arena, (size_t)BH*BW);
        float *bf = ARENA_ALLOC(arena, (size_t)BH*BW*3*sizeof *bf);
        memset(bd, 1, (size_t)BH*BW);
        for (int r = 0; r < BH; r++) for (int c = 0; c < BW; c++) {
            size_t k=(size_t)r*BW+(size_t)c;
            bf[k*3]=3.0f; bf[k*3+1]=(float)r; bf[k*3+2]=(float)c;
        }
        bd[(size_t)3*BW+3]=0;                       /* closed hole */
        for (int r=0;r<=3;r++) bd[(size_t)r*BW+8]=0; /* open notch */
        QuadStripOpts ob=o; ob.mode=QUAD_STRIP_BOUNDED; ob.max_hole_distance=2;
        float *bv=NULL,*bu=NULL;int32_t *bt=NULL;uint8_t *bl=NULL;
        size_t bnv=0,bnf=0;
        int rc=QuadStrip_build(arena,bd,bf,BH,BW,0,0,1.0,&ob,
                               &bv,&bnv,&bt,&bnf,&bu,&bl);
        int found_hole=0,found_notch=0;
        for(size_t v=0;v<bnv;v++) {
            int c=(int)lround((double)bu[v*2]);
            int r=(int)lround((double)bu[v*2+1]);
            if(r==3&&c==3&&bl[v])found_hole=1;
            if(r==2&&c==8)found_notch=1;
        }
        CK(rc==0&&bnf>0,"t5 bounded topology build ok");
        CK(found_hole,"t5 closed same-sheet hole is lofted");
        CK(!found_notch,"t5 open exterior notch remains alpha");
    }

    /* t5b: the output footprint and the initializer evidence are independent.
     * Removing an observation must PDE-fill it without punching the surface. */
    {
        enum {TH=3,TW=3};
        uint8_t topology[TH*TW],initializer[TH*TW];
        float tf[TH*TW*3];
        memset(topology,1,sizeof topology);
        memset(initializer,1,sizeof initializer);
        initializer[TW+1]=0;
        for(int r=0;r<TH;r++)for(int c=0;c<TW;c++) {
            size_t k=(size_t)r*TW+(size_t)c;
            tf[k*3]=11.0f;tf[k*3+1]=(float)r;tf[k*3+2]=2.0f*(float)c;
        }
        QuadStripOpts ot=o;ot.mode=QUAD_STRIP_BOUNDED;
        float *tv=NULL,*tu=NULL;int32_t *tt=NULL;uint8_t *tl=NULL;
        size_t tnv=0,tnf=0;
        int rc=QuadStrip_build_topology_with_phase(
            arena,topology,initializer,tf,NULL,TH,TW,0,0,1.0,&ot,
            &tv,&tnv,&tt,&tnf,&tu,&tl);
        int center=-1;
        for(size_t v=0;v<tnv;v++)if(lround((double)tu[v*2])==1&&
                                           lround((double)tu[v*2+1])==1)
            center=(int)v;
        CK(rc==0&&tnv==9&&tnf==8,"t5b explicit topology keeps the full surface");
        CK(center>=0&&tl[center]&&fabs((double)tv[(size_t)center*3+2]-2.0)<0.05,
           "t5b missing initializer is harmonically filled inside topology");
    }

    /* t5c: alternating live rows form a continuous ragged strip when the
     * nearest-row bridge is explicitly permitted. */
    {
        enum {RH=3,RW=4};
        uint8_t rt[RH*RW],rd[RH*RW];float rf[RH*RW*3];
        memset(rt,0,sizeof rt);memset(rd,0,sizeof rd);
        for(int r=0;r<RH;r++)for(int c=0;c<RW;c++) {
            size_t k=(size_t)r*RW+(size_t)c;
            if(r==0||r==2)rt[k]=rd[k]=1;
            rf[k*3]=(float)r;rf[k*3+1]=0.0f;rf[k*3+2]=(float)c;
        }
        QuadStripOpts ro=o;ro.mode=QUAD_STRIP_BOUNDED;ro.max_row_gap=2;
        float *rv=NULL,*ru=NULL;int32_t *rfaces=NULL;uint8_t *rl=NULL;
        size_t rnv=0,rnf=0;
        int rc=QuadStrip_build_topology_with_phase(
            arena,rt,rd,rf,NULL,RH,RW,0,0,1.0,&ro,
            &rv,&rnv,&rfaces,&rnf,&ru,&rl);
        CK(rc==0&&rnv==8&&rnf==6,
           "t5c gap-two rows emit all three nearest-row quads");
        int valid=1;for(size_t f=0;f<rnf*3;f++)
            if(rfaces[f]<0||(size_t)rfaces[f]>=rnv)valid=0;
        CK(valid,"t5c ragged face indices are valid");
    }

    /* t6: degenerate inputs */
    {
        size_t nv = 0;
        CK(QuadStrip_build(arena, NULL, field, H, W, 0, 0, 2.0, &o,
                           NULL, &nv, NULL, NULL, NULL, NULL) != 0, "t6 null domain -> err");
        CK(QuadStrip_build(arena, domain, field, 0, W, 0, 0, 2.0, &o,
                           NULL, &nv, NULL, NULL, NULL, NULL) != 0, "t6 H=0 -> err");
    }

    /* t7: full-ribbon metric ARAP must use a real coarse hierarchy and
     * uncompress a deliberately half-width planar lattice. */
    {
        const int AH=17,AW=19;size_t an=(size_t)AH*(size_t)AW;
        float *av=ARENA_ALLOC(arena,an*3*sizeof *av);
        float *au=ARENA_ALLOC(arena,an*2*sizeof *au);
        uint8_t *af=ARENA_ALLOC(arena,an);
        for(int r=0;r<AH;r++)for(int c=0;c<AW;c++){
            size_t v=(size_t)r*(size_t)AW+(size_t)c;
            av[v*3]=(float)r;av[v*3+1]=0.0f;av[v*3+2]=0.5f*(float)c;
            au[v*2]=(float)c;au[v*2+1]=(float)r;af[v]=1;
        }
        QuadStripArapOpts ao;QuadStripArap_defaults(&ao);
        ao.max_iterations=80;ao.movement_tolerance=1e-6;
        ao.source_weight=0.0;ao.fill_weight=1e-3;ao.verbose=0;
        QuadStripArapStats as;int rc=QuadStrip_metric_arap(av,au,af,AH,AW,&ao,&as);
        CK(rc==0,"t7 metric ARAP executes safely");
        CK(as.multigrid_levels_used>1&&as.multigrid_cycles>0,
           "t7 metric ARAP uses a coarse Galerkin hierarchy");
        CK(as.final_edge_log_rms<0.5*as.initial_edge_log_rms,
           "t7 metric ARAP uncompresses the rest metric");
        CK(isfinite(av[(an-1)*3+2]),"t7 metric ARAP output is finite");
    }

    /* t8: force the nonlinear hierarchy on a small fixture so resampling,
     * coarse ARAP, barriered prolongation, and the original-anchor fine solve
     * remain covered without allocating a production-scale ribbon. */
    {
        const int AH=33,AW=35;size_t an=(size_t)AH*(size_t)AW;
        float *av=ARENA_ALLOC(arena,an*3*sizeof *av);
        float *au=ARENA_ALLOC(arena,an*2*sizeof *au);
        uint8_t *af=ARENA_ALLOC(arena,an);
        for(int r=0;r<AH;r++)for(int c=0;c<AW;c++){
            size_t v=(size_t)r*(size_t)AW+(size_t)c;
            av[v*3]=(float)r;av[v*3+1]=0.02f*(float)r;av[v*3+2]=0.6f*(float)c;
            au[v*2]=(float)c;au[v*2+1]=(float)r;af[v]=1;
        }
        QuadStripArapOpts ao;QuadStripArap_defaults(&ao);
        ao.max_iterations=24;ao.movement_tolerance=1e-5;
        ao.source_weight=0.0;ao.fill_weight=1e-3;ao.verbose=0;
        ao.nonlinear_multigrid_levels=3;ao.nonlinear_coarse_iterations=12;
        ao.nonlinear_min_vertices=0;
        QuadStripArapStats as;int rc=QuadStrip_metric_arap(av,au,af,AH,AW,&ao,&as);
        CK(rc==0,"t8 nonlinear metric ARAP executes safely");
        CK(as.nonlinear_coarse_levels==2&&as.nonlinear_coarse_iterations>0,
           "t8 nonlinear metric ARAP visits both coarse levels");
        CK(as.nonlinear_min_prolongation_scale>0.0,
           "t8 nonlinear coarse displacement passes the orientation barrier");
        CK(as.final_edge_log_rms<as.initial_edge_log_rms,
           "t8 nonlinear metric ARAP improves the original fine rest metric");
    }

    /* t9: a nearly collapsed but still valid triangle can require a step below
     * 2^-24.  The smooth local barrier must hard-freeze its implicated control
     * tile and retain a valid iterate, never abort the whole ribbon. */
    {
        const int BH=2,BW=2;const double eps=1e-9;
        double bp[12]={0.0,0.0,0.0, 0.0,0.0,1.0,
                       0.0,eps,0.0, 0.0,eps,1.0};
        double bd[12]={0};float bs[4];int bad[4],dilated[4];
        bd[2*3+1]=-2.0;bd[3*3+1]=-2.0;
        int rounds=0,zeros=0;double minimum=1.0;
        int rc=metric_local_barrier_scales(bp,bd,BH,BW,0.05,2.0,0.0,
                                            bs,bad,dilated,&rounds,&minimum,&zeros);
        double proposal[12];
        for(int i=0;i<12;i++)proposal[i]=bp[i]+(double)bs[i/3]*bd[i];
        CK(rc==0&&rounds>METRIC_BARRIER_HALVING_ROUNDS&&zeros>0,
           "t9 exhausted halving uses a local zero-step fallback");
        CK(metric_trial_valid(bp,proposal,BH,BW,1.0,0.05),
           "t9 local zero-step fallback preserves triangle orientation");
    }

    /* t10: metric ARAP alone has no bending preference.  The repair path must
     * fair a selected high-frequency tangent crumple while the unselected rim
     * remains an effectively fixed Dirichlet boundary. */
    {
        const int RH=17,RW=33;size_t rn=(size_t)RH*(size_t)RW;
        float *rv=ARENA_ALLOC(arena,rn*3*sizeof *rv);
        float *ru=ARENA_ALLOC(arena,rn*2*sizeof *ru);
        float *rw=ARENA_ALLOC(arena,rn*sizeof *rw);
        uint8_t *rf=ARENA_ALLOC(arena,rn);
        uint8_t *rfix=ARENA_ALLOC(arena,rn);
        double rough_before=0.0,rough_after=0.0;
        for(int r=0;r<RH;r++)for(int c=0;c<RW;c++){
            size_t v=(size_t)r*(size_t)RW+(size_t)c;
            int active=r>=4&&r<=12&&c>=8&&c<=24;
            rv[v*3]=(float)r;
            rv[v*3+1]=active?(c&1?1.5f:-1.5f):0.0f;
            rv[v*3+2]=(float)c;
            ru[v*2]=(float)c;ru[v*2+1]=(float)r;
            rw[v]=active?1.0f:0.0f;rf[v]=0;
            rfix[v]=(uint8_t)(r==0||r==RH-1||c==0||c==RW-1);
        }
        for(int r=4;r<=12;r++)for(int c=9;c<24;c++){
            size_t v=(size_t)r*(size_t)RW+(size_t)c;
            rough_before+=fabs((double)rv[(v-1)*3+1]-2.0*(double)rv[v*3+1]
                               +(double)rv[(v+1)*3+1]);
        }
        QuadStripArapOpts ro;QuadStripArap_defaults(&ro);
        ro.max_iterations=96;ro.movement_tolerance=1e-5;
        ro.source_weight=100.0;ro.fill_weight=0.01;
        ro.fixed_vertices=rfix;
        ro.repair_weight=rw;ro.repair_source_scale=0.0001;
        ro.repair_frame_blend=0.9;ro.repair_frame_screen=0.002;
        ro.repair_frame_sweeps=64;ro.nonlinear_multigrid_levels=1;
        QuadStripArapStats rs;
        int rc=QuadStrip_metric_arap(rv,ru,rf,RH,RW,&ro,&rs);
        for(int r=4;r<=12;r++)for(int c=9;c<24;c++){
            size_t v=(size_t)r*(size_t)RW+(size_t)c;
            rough_after+=fabs((double)rv[(v-1)*3+1]-2.0*(double)rv[v*3+1]
                              +(double)rv[(v+1)*3+1]);
        }
        CK(rc==0&&rs.repair_vertices==9u*17u,
           "t10 marble-aware metric ARAP executes on only the repair mask");
        CK(rough_after<0.35*rough_before,
           "t10 fair UV tangent target removes the selected crumple");
        CK(rv[0]==0.0f&&rv[1]==0.0f&&rv[2]==0.0f&&
           fabs((double)rv[((size_t)(RH-1)*RW+(RW-1))*3]-(RH-1))<1e-7&&
           fabs((double)rv[((size_t)(RH-1)*RW+(RW-1))*3+1])<1e-7&&
           fabs((double)rv[((size_t)(RH-1)*RW+(RW-1))*3+2]-(RW-1))<1e-7,
           "t10 exact Dirichlet rim remains bit-exact");
    }

    /* t11: the radial-normal prior removes axial "bucket-lid" normals (n ~ z)
     * that ordinary ARAP is metric-blind to; it is wired end to end and a
     * bit-for-bit no-op at weight 0. */
    {
        const int LH=21,LW=41;size_t ln=(size_t)LH*(size_t)LW;
        float *lv=ARENA_ALLOC(arena,ln*3*sizeof *lv);
        float *lvA=ARENA_ALLOC(arena,ln*3*sizeof *lvA);
        float *lvB=ARENA_ALLOC(arena,ln*3*sizeof *lvB);
        float *lvC=ARENA_ALLOC(arena,ln*3*sizeof *lvC);
        float *lu=ARENA_ALLOC(arena,ln*2*sizeof *lu);
        uint8_t *lf=ARENA_ALLOC(arena,ln);
        const double AY=64.0,AX=64.0,R0=20.0;
        for(int r=0;r<LH;r++)for(int c=0;c<LW;c++){
            size_t v=(size_t)r*(size_t)LW+(size_t)c;
            double th=0.10*(double)c;
            /* cylinder (radial normals) for r<10; then a near-isometric inward
             * curl (radius shrinks) that ARAP cannot see but turns the normal
             * axial -- a bucket lid. */
            double zc=(r<10)?(double)r:10.0+0.2*(double)(r-10);
            double rad=(r<10)?R0:R0-1.5*(double)(r-10);
            lv[v*3]=(float)zc;lv[v*3+1]=(float)(AY+rad*cos(th));
            lv[v*3+2]=(float)(AX+rad*sin(th));
            lu[v*2]=(float)c;lu[v*2+1]=(float)r;lf[v]=1;
        }
        memcpy(lvA,lv,ln*3*sizeof *lv);memcpy(lvB,lv,ln*3*sizeof *lv);
        memcpy(lvC,lv,ln*3*sizeof *lv);
#define METRIC_AXIAL_FRAC(V,OUT) do { size_t _ax=0,_tot=0; \
        for(int _r=0;_r+1<LH;_r++)for(int _c=0;_c+1<LW;_c++){ \
            size_t _a=(size_t)_r*(size_t)LW+(size_t)_c,_b=_a+1,_q=_a+(size_t)LW,_d=_q+1; \
            size_t _tri[2][3]={{_a,_b,_d},{_a,_d,_q}}; \
            for(int _t=0;_t<2;_t++){double _e1[3],_e2[3],_nn[3]; \
                for(int _k=0;_k<3;_k++){_e1[_k]=(double)(V)[_tri[_t][1]*3+_k]-(double)(V)[_tri[_t][0]*3+_k]; \
                    _e2[_k]=(double)(V)[_tri[_t][2]*3+_k]-(double)(V)[_tri[_t][0]*3+_k];} \
                _nn[0]=_e1[1]*_e2[2]-_e1[2]*_e2[1];_nn[1]=_e1[2]*_e2[0]-_e1[0]*_e2[2]; \
                _nn[2]=_e1[0]*_e2[1]-_e1[1]*_e2[0]; \
                double _L=sqrt(_nn[0]*_nn[0]+_nn[1]*_nn[1]+_nn[2]*_nn[2]); \
                if(_L>1e-9){_tot++;if(fabs(_nn[0])/_L>0.7)_ax++;}}} \
        (OUT)=_tot?(double)_ax/(double)_tot:0.0; } while(0)
        double before=0.0;METRIC_AXIAL_FRAC(lv,before);
        QuadStripArapOpts la;QuadStripArap_defaults(&la);
        la.max_iterations=300;la.movement_tolerance=1e-6;
        la.source_weight=0.0;la.fill_weight=1e-3;la.nonlinear_multigrid_levels=1;
        la.axis_y=AY;la.axis_x=AX;              /* axis set, weight 0 -> no-op */
        QuadStripArapStats sA;int rcA=QuadStrip_metric_arap(lvA,lu,lf,LH,LW,&la,&sA);
        QuadStripArapOpts lb=la;lb.axis_normal_weight=0.5;lb.axis_fold_limit=0.9;
        QuadStripArapStats sB;int rcB=QuadStrip_metric_arap(lvB,lu,lf,LH,LW,&lb,&sB);
        QuadStripArapOpts lc;QuadStripArap_defaults(&lc); /* no axis fields at all */
        lc.max_iterations=300;lc.movement_tolerance=1e-6;
        lc.source_weight=0.0;lc.fill_weight=1e-3;lc.nonlinear_multigrid_levels=1;
        QuadStripArapStats sC;int rcC=QuadStrip_metric_arap(lvC,lu,lf,LH,LW,&lc,&sC);
        double afterA=0.0,afterB=0.0;METRIC_AXIAL_FRAC(lvA,afterA);METRIC_AXIAL_FRAC(lvB,afterB);
#undef METRIC_AXIAL_FRAC
        CK(rcA==0&&rcB==0&&rcC==0,"t11 radial-normal prior executes safely");
        CK(before>0.15,"t11 lid fixture starts with axial normals");
        CK(memcmp(lvA,lvC,ln*3*sizeof *lvA)==0,
           "t11 weight-0 with axis set is bit-for-bit identical to plain ARAP");
        CK(memcmp(lvA,lvB,ln*3*sizeof *lvA)!=0,
           "t11 radial-normal prior changes the solve");
        CK(afterB<=afterA+0.02,
           "t11 radial-normal prior does not increase axial-normal lids");
    }

    /* t12: cylindrical post-snap fairing may change radius/phase, but lattice v
     * already determines world z.  preserve_axial is therefore an exact
     * per-vertex constraint, not another finite penalty. */
    {
        const int PH=11,PW=15;size_t pn=(size_t)PH*(size_t)PW;
        float *pv=ARENA_ALLOC(arena,pn*3*sizeof *pv);
        float *pu=ARENA_ALLOC(arena,pn*2*sizeof *pu);
        float *pz=ARENA_ALLOC(arena,pn*sizeof *pz);
        uint8_t *pf=ARENA_ALLOC(arena,pn);
        for(int r=0;r<PH;r++)for(int c=0;c<PW;c++){
            size_t v=(size_t)r*(size_t)PW+(size_t)c;
            pv[v*3]=(float)(100.0+r+0.07*sin(0.4*c));pz[v]=pv[v*3];
            pv[v*3+1]=(float)(20.0+0.03*r);pv[v*3+2]=0.55f*(float)c;
            pu[v*2]=(float)c;pu[v*2+1]=(float)r;pf[v]=1;
        }
        QuadStripArapOpts po;QuadStripArap_defaults(&po);
        po.max_iterations=48;po.movement_tolerance=1e-6;
        po.source_weight=0.0;po.fill_weight=1e-3;
        po.nonlinear_multigrid_levels=1;po.preserve_axial=1;
        QuadStripArapStats ps;
        int rc=QuadStrip_metric_arap(pv,pu,pf,PH,PW,&po,&ps),zexact=1;
        for(size_t v=0;v<pn;v++)if(pv[v*3]!=pz[v]){zexact=0;break;}
        CK(rc==0,"t12 axial-preserving metric ARAP executes safely");
        CK(zexact,"t12 axial-preserving metric ARAP keeps every z bit-exact");
        CK(ps.final_edge_log_rms<ps.initial_edge_log_rms,
           "t12 transverse relaxation still improves the rest metric");
    }

    Arena_restore(arena, mk);
    Arena_dispose(&arena);
    fprintf(stderr, "[selftest] %s (%d failures)\n", fails == 0 ? "ALL PASS" : "FAILURES", fails);
    return fails;
}
