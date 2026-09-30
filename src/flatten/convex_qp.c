#include "convex_qp.h"
#include "osqp.h"
#include "../common/ves_platform.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef struct CqpEntry {int row,col;double value;} CqpEntry;
static int cqp_order(const void *a,const void *b)
{
    const CqpEntry *x=a,*y=b;if(x->col!=y->col)return(x->col>y->col)-(x->col<y->col);
    return(x->row>y->row)-(x->row<y->row);
}
static int cqp_csc(OSQPCscMatrix *out,int m,int n,int count,const int *row,const int *col,
                   const double *value,const double *scale,int hessian)
{
    memset(out,0,sizeof *out);out->m=m;out->n=n;out->nz=-1;
    CqpEntry *entries=malloc((size_t)(count?count:1)*sizeof *entries);if(!entries)return 0;
    int ok=0;
    for(int i=0;i<count;i++){
        if(row[i]<0 || row[i]>=m || col[i]<0 || col[i]>=n || !isfinite(value[i]) || (hessian && row[i]<col[i]))goto done;
        /* C's local/global Hessian and target are twice Python's normal. */
        double v=value[i]*scale[col[i]]*(hessian?.5*scale[row[i]]:1);
        if(!isfinite(v))goto done;
        entries[i]=(CqpEntry){hessian?col[i]:row[i],hessian?row[i]:col[i],v};
    }
    qsort(entries,(size_t)count,sizeof *entries,cqp_order);
    out->p=calloc((size_t)n+1,sizeof(OSQPInt));out->i=malloc((size_t)(count?count:1)*sizeof(OSQPInt));out->x=malloc((size_t)(count?count:1)*sizeof(OSQPFloat));
    if(!out->p || !out->i || !out->x)goto done;
    int used=0;
    for(int at=0;at<count;){int end=at+1;double sum=entries[at].value;
        while(end<count && entries[end].col==entries[at].col && entries[end].row==entries[at].row)sum+=entries[end++].value;
        if(!isfinite(sum))goto done;
        if(sum){out->i[used]=entries[at].row;out->x[used++]=sum;out->p[entries[at].col+1]++;}at=end;
    }
    for(int c=1;c<=n;c++)out->p[c]+=out->p[c-1];out->nzmax=used;ok=1;
done:
    free(entries);return ok;
}
static void cqp_free(OSQPCscMatrix *m){free(m->p);free(m->i);free(m->x);}
/* Per thread: concurrent solves each keep their own. */
static VES_THREAD_LOCAL double cqp_time_limit=0;
static VES_THREAD_LOCAL double cqp_deadline=0;
void ConvexQp_deadline(double deadline){cqp_deadline=deadline>0 && isfinite(deadline)?deadline:0;}
void ConvexQp_time_limit(double seconds){cqp_time_limit=seconds>0 && isfinite(seconds)?seconds:0;}
static int cqp_solve(int n,int hn,const int *hr,const int *hc,const double *hv,
    const double *target,int m,int an,const int *ar,const int *ac,const double *av,
    const double *upper,double *result,int initialization)
{
    if(cqp_deadline>0 && ves_clock_sec()>=cqp_deadline)return 0;
    if(n<=0 || hn<0 || m<0 || an<0 || !hr || !hc || !hv || !target || !result || (m && !upper) || (an && (!ar || !ac || !av)))return 0;
    double *scale=calloc((size_t)n,sizeof(double)),*q=malloc((size_t)n*sizeof(double));
    double *l=malloc((size_t)(m?m:1)*sizeof(double));OSQPCscMatrix h={0},a={0};OSQPSolver *solver=NULL;int ok=0;
    if(!scale || !q || !l)goto done;
    for(int k=0;k<hn;k++){
        if(hr[k]<0 || hr[k]>=n || hc[k]<0 || hc[k]>hr[k] || !isfinite(hv[k]))goto done;
        if(hr[k]==hc[k])scale[hr[k]]+=.5*hv[k];
    }
    for(int i=0;i<n;i++){if(!(scale[i]>0) || !isfinite(scale[i]) || !isfinite(target[i]))goto done;scale[i]=initialization?1:1/sqrt(scale[i]);q[i]=-.5*scale[i]*target[i];}
    for(int i=0;i<m;i++){if(!isfinite(upper[i]))goto done;l[i]=-OSQP_INFTY;}
    if(!cqp_csc(&h,n,n,hn,hr,hc,hv,scale,1) || !cqp_csc(&a,m,n,an,ar,ac,av,scale,0))goto done;
    OSQPSettings settings;osqp_set_default_settings(&settings);
    settings.verbose=0;settings.polishing=1;settings.time_limit=cqp_time_limit>0?cqp_time_limit:OSQP_TIME_LIMIT;  /* OSQP validates time_limit > 0; its default 1e10 means none */
    if(cqp_deadline>0){double left=cqp_deadline-ves_clock_sec();if(left<=0)goto done;settings.time_limit=fmin(settings.time_limit,left);}
    if(initialization){settings.eps_abs=1e-7;settings.eps_rel=1e-8;settings.max_iter=20000;}
    else{
        settings.eps_abs=1e-6;settings.eps_rel=1e-7;settings.eps_prim_inf=1e-9;settings.eps_dual_inf=1e-9;
        settings.max_iter=4000;settings.adaptive_rho_interval=25;
    }
    if(osqp_setup(&solver,&h,q,&a,l,upper,m,n,&settings) || osqp_solve(solver) || !solver->solution || !solver->solution->x)goto done;
    int status=solver->info->status_val;
    if(status!=OSQP_SOLVED && status!=OSQP_SOLVED_INACCURATE && (initialization || status!=OSQP_MAX_ITER_REACHED))goto done;
    for(int i=0;i<n;i++){double x=scale[i]*solver->solution->x[i];if(!isfinite(x))goto done;result[i]=x;}
    ok=1;
done:
    if(solver)osqp_cleanup(solver);cqp_free(&h);cqp_free(&a);free(scale);free(q);free(l);return ok;
}
int ConvexQp_solve(int n,int hn,const int *hr,const int *hc,const double *hv,
    const double *target,int m,int an,const int *ar,const int *ac,const double *av,const double *upper,double *result)
{return cqp_solve(n,hn,hr,hc,hv,target,m,an,ar,ac,av,upper,result,0);}
int ConvexQp_translation_seed(int n,int hn,const int *hr,const int *hc,const double *hv,
    const double *target,int m,int an,const int *ar,const int *ac,const double *av,const double *upper,double *result)
{return cqp_solve(n,hn,hr,hc,hv,target,m,an,ar,ac,av,upper,result,1);}
int ConvexQp_selftest(void)
{
    int hr[]={0,1,1},hc[]={0,0,1},ar[]={0,1},ac[]={0,1};
    double hv[]={4,1,2},rhs[]={4,3},av[]={1,1},upper[]={.5,.75},x[2]={0};
    int ok=ConvexQp_solve(2,3,hr,hc,hv,rhs,2,2,ar,ac,av,upper,x) && fabs(x[0]-.5)<1e-6 && fabs(x[1]-.75)<1e-6;
    ConvexQp_deadline(ves_clock_sec()-1);double saved[2]={x[0],x[1]};
    ok=ok && !ConvexQp_solve(2,3,hr,hc,hv,rhs,2,2,ar,ac,av,upper,x) && !memcmp(saved,x,sizeof saved);
    ConvexQp_deadline(0);
    ok=ok && ConvexQp_translation_seed(2,3,hr,hc,hv,rhs,2,2,ar,ac,av,upper,x) && fabs(x[0]-.5)<1e-7 && fabs(x[1]-.75)<1e-7;
    upper[1]=-1;ac[1]=0;av[1]=-1;
    /* x <= 0.5 and -x <= -1 are contradictory. */
    ok=ok && !ConvexQp_solve(2,3,hr,hc,hv,rhs,2,2,ar,ac,av,upper,x);
    ok=ok && !ConvexQp_translation_seed(2,3,hr,hc,hv,rhs,2,2,ar,ac,av,upper,x);
    fprintf(stderr,"  native OSQP equilibrated bounds and infeasibility: %s\n",ok?"ok":"FAIL");return ok?0:1;
}
