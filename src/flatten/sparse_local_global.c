/* Local convex projections and global weighted least squares. Splitting
 * C*x=z, z>=lower gives independent scalar projections and the SPD system
 * (H+rho*C'C)x=b+rho*C'(z-u). The scaled correction u retains the hard
 * constraint in the limit; this is not a fixed soft-penalty replacement.
 * Boyd et al., Distributed Optimization ... ADMM (2011), sections 3, 4.2.
 * No indefinite system or primal-dual Newton solve is formed. */
#include "sparse_local_global.h"
#include "sparse_solve.h"
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>

static double slg_dot(const SparseLocalRow *r, const double *x)
{
    double q=0;
    for (size_t k=0;k<6;k++) q+=r->value[k]*x[r->index[k]];
    return q;
}

static void slg_multiply(size_t n,size_t nt,const int *row,const int *col,
    const double *h,const double *x,double *out)
{
    memset(out,0,n*sizeof *out);
    for(size_t k=0;k<nt;k++) {
        out[row[k]]+=h[k]*x[col[k]];
        if(row[k]!=col[k])out[col[k]]+=h[k]*x[row[k]];
    }
}

static int slg_certificate(size_t n,size_t nt,const int *row,const int *col,
    const double *h,const double *b,size_t nr,const SparseLocalRow *c,
    const double *x,const double *lambda,double *work,SparseLocalGlobalStats *s)
{
    double q=0,bx=0,scale=1,primal=0,stationarity=0,gap=0;
    slg_multiply(n,nt,row,col,h,x,work);
    for(size_t k=0;k<n;k++) {
        if(!isfinite(x[k])||!isfinite(work[k]))return -1;
        q+=x[k]*work[k];bx+=b[k]*x[k];work[k]-=b[k];scale=fmax(scale,fabs(b[k]));
    }
    for(size_t r=0;r<nr;r++) {
        double slack=slg_dot(c+r,x)-c[r].lower;
        if(!isfinite(slack)||!isfinite(lambda[r])||lambda[r]<0)return -1;
        primal=fmax(primal,-slack);gap+=fabs(lambda[r]*slack);
        for(size_t k=0;k<6;k++)work[c[r].index[k]]-=lambda[r]*c[r].value[k];
    }
    for(size_t k=0;k<n;k++)stationarity=fmax(stationarity,fabs(work[k]));
    s->primal_residual=primal;s->stationarity_residual=stationarity/scale;
    s->complementarity=gap/fmax(1,fmax(fabs(q),fabs(bx)));s->objective=.5*q-bx;
    if(!isfinite(s->objective)||!isfinite(s->complementarity)||!isfinite(s->stationarity_residual))return -1;
    return primal<=1e-9&&s->stationarity_residual<=1e-9&&s->complementarity<=1e-10;
}

#include "sparse_local_global_refine.inc"

int Sparse_local_global(Arena_T arena,size_t n,size_t nt,
    const int *row,const int *col,const double *h,const double *b,
    size_t nr,const SparseLocalRow *original,const double *initial,
    const double *initial_multipliers,int limit,double *out,double *out_multipliers,
    SparseLocalGlobalStats *stats)
{
    if(!stats)return -1;
    memset(stats,0,sizeof *stats);
    if(!arena||!n||n>INT_MAX||!nt||nt>INT_MAX||!row||!col||!h||!b||!out||
        (nr&&(!original||!out_multipliers))||limit<0||nr>((size_t)INT_MAX-nt)/21)return -1;
    size_t capacity=nt+21*nr;
    /* Allocate before owning an external factor; arena exceptions cannot
     * leak it. Scratch lifetime is controlled by the caller's arena mark. */
    int *ar=ARENA_ALLOC(arena,capacity*sizeof *ar),*ac=ARENA_ALLOC(arena,capacity*sizeof *ac);
    double *av=ARENA_ALLOC(arena,capacity*sizeof *av),*scaled=ARENA_ALLOC(arena,capacity*sizeof *scaled);
    double *diagonal=ARENA_CALLOC(arena,n,sizeof *diagonal),*scale=ARENA_ALLOC(arena,n*sizeof *scale);
    double *x=ARENA_CALLOC(arena,n,sizeof *x),*work=ARENA_ALLOC(arena,n*sizeof *work);
    double *rhs=ARENA_ALLOC(arena,n*sizeof *rhs),*solution=ARENA_ALLOC(arena,n*sizeof *solution);
    double *correction=ARENA_ALLOC(arena,n*sizeof *correction),*residual=ARENA_ALLOC(arena,n*sizeof *residual);
    SparseLocalRow *c=ARENA_CALLOC(arena,nr?nr:1,sizeof *c);
    double *row_scale=ARENA_ALLOC(arena,(nr?nr:1)*sizeof *row_scale);
    double *z=ARENA_ALLOC(arena,(nr?nr:1)*sizeof *z),*u=ARENA_ALLOC(arena,(nr?nr:1)*sizeof *u);
    double *lambda=ARENA_CALLOC(arena,nr?nr:1,sizeof *lambda);
    for(size_t k=0;k<nt;k++) {
        if(col[k]<0||row[k]<col[k]||(size_t)row[k]>=n||!isfinite(h[k]))return -1;
        if(row[k]==col[k])diagonal[row[k]]+=h[k];
    }
    for(size_t k=0;k<n;k++) {
        if(!(diagonal[k]>0)||!isfinite(diagonal[k])||!isfinite(b[k]))return -1;
        if(initial){if(!isfinite(initial[k]))return -1;x[k]=initial[k];}
    }
    for(size_t r=0;r<nr;r++) {
        size_t used=0;
        if(!isfinite(original[r].lower))return -1;
        for(size_t k=0;k<6;k++) {
            int i=original[r].index[k];double g=original[r].value[k];
            if(i<0||(size_t)i>=n||!isfinite(g))return -1;
            size_t j=0;for(;j<used;j++)if(c[r].index[j]==i)break;
            if(j==used)c[r].index[used++]=i;
            c[r].value[j]+=g;
        }
        double norm=0;
        for(size_t k=0;k<used;k++)norm+=c[r].value[k]*c[r].value[k]/diagonal[c[r].index[k]];
        if(!isfinite(norm))return -1;
        if(norm==0) {
            if(original[r].lower>0)return -1;
            row_scale[r]=1;
        } else row_scale[r]=1/sqrt(norm);
        for(size_t k=0;k<used;k++)c[r].value[k]*=row_scale[r];
        c[r].lower=original[r].lower*row_scale[r];
        if(initial_multipliers) {
            if(!isfinite(initial_multipliers[r])||initial_multipliers[r]<0)return -1;
            lambda[r]=initial_multipliers[r];
        }
    }
    stats->warm=initial!=NULL;
    int certified=slg_certificate(n,nt,row,col,h,b,nr,original,x,lambda,work,stats);
    if(certified<0)return -1;
    if(certified) {
        memcpy(out,x,n*sizeof *out);if(nr)memcpy(out_multipliers,lambda,nr*sizeof *lambda);
        stats->converged=1;return 0;
    }
    if(!limit)return 1;
    double rho=1,last_refinement_merit=HUGE_VAL;
    for(size_t r=0;r<nr;r++) {
        z[r]=fmax(c[r].lower,slg_dot(c+r,x));u[r]=-lambda[r]/(rho*row_scale[r]);
    }
    SparseFactor_T factor=NULL;
    int status=1;size_t at=0;
    for(int iteration=0;iteration<limit&&stats->iterations<(size_t)limit;iteration++) {
        if(!factor) {
            memcpy(ar,row,nt*sizeof *ar);memcpy(ac,col,nt*sizeof *ac);memcpy(av,h,nt*sizeof *av);at=nt;
            for(size_t r=0;r<nr;r++)for(size_t k=0;k<6;k++)if(c[r].value[k]!=0)
                for(size_t j=0;j<=k;j++)if(c[r].value[j]!=0) {
                    int p=c[r].index[k],q=c[r].index[j];
                    ar[at]=p>q?p:q;ac[at]=p>q?q:p;av[at++]=rho*c[r].value[k]*c[r].value[j];
                }
            memset(diagonal,0,n*sizeof *diagonal);
            for(size_t k=0;k<at;k++)if(ar[k]==ac[k])diagonal[ar[k]]+=av[k];
            for(size_t k=0;k<n;k++) {
                if(!(diagonal[k]>0)||!isfinite(diagonal[k])){status=-1;goto done;}
                scale[k]=1/sqrt(diagonal[k]);
            }
            for(size_t k=0;k<at;k++)scaled[k]=av[k]*scale[ar[k]]*scale[ac[k]];
            if(Sparse_factor_spd((int)n,(int)at,ar,ac,scaled,&factor)){status=-1;goto done;}
            stats->factors++;
        }
        memcpy(rhs,b,n*sizeof *rhs);
        for(size_t r=0;r<nr;r++)for(size_t k=0;k<6;k++)rhs[c[r].index[k]]+=rho*c[r].value[k]*(z[r]-u[r]);
        for(size_t k=0;k<n;k++)rhs[k]*=scale[k];
        if(Sparse_factor_solve(factor,rhs,solution)){status=-1;goto done;}
        stats->solves++;
        /* Residual refinement uses the SAME scaled SPD factor. */
        for(int refine=0;refine<3;refine++) {
            double error=0,norm=0;
            slg_multiply(n,at,ar,ac,scaled,solution,residual);
            for(size_t k=0;k<n;k++) {
                residual[k]=rhs[k]-residual[k];error=fmax(error,fabs(residual[k]));norm=fmax(norm,fabs(rhs[k]));
            }
            double relative=error/(norm>0?norm:1);
            stats->linear_residual=fmax(stats->linear_residual,relative);
            if(!isfinite(relative)){status=-1;goto done;}
            if(relative<=1e-12)break;
            if(refine==2||Sparse_factor_solve(factor,residual,correction)){status=-1;goto done;}
            stats->solves++;
            for(size_t k=0;k<n;k++)solution[k]+=correction[k];
        }
        for(size_t k=0;k<n;k++)x[k]=solution[k]*scale[k];
        double primal=0,primal_scale=1;
        for(size_t r=0;r<nr;r++) {
            double cx=slg_dot(c+r,x),relaxed=1.6*cx-.6*z[r];
            z[r]=fmax(c[r].lower,relaxed+u[r]);u[r]+=relaxed-z[r];
            /* Projection implies u<=0, exactly apart from roundoff. */
            u[r]=fmin(0,u[r]);lambda[r]=-rho*u[r]*row_scale[r];
            primal=fmax(primal,fabs(cx-z[r]));primal_scale=fmax(primal_scale,fmax(fabs(cx),fabs(z[r])));
        }
        stats->iterations++;stats->penalty=rho;
        certified=slg_certificate(n,nt,row,col,h,b,nr,original,x,lambda,work,stats);
        if(certified<0){status=-1;goto done;}
        if(certified) {
            memcpy(out,x,n*sizeof *out);if(nr)memcpy(out_multipliers,lambda,nr*sizeof *lambda);
            stats->converged=1;status=0;goto done;
        }
        /* A validated pair needs less settling than a cold map. An unfinished
         * refinement may become useful again after ADMM improves its start.
         * Require a factor-two reduction of the ORIGINAL certificate merit
         * since the last start, so repeated phases cannot hide stagnation. */
        double merit=fmax(stats->primal_residual/1e-9,fmax(stats->stationarity_residual/1e-9,stats->complementarity/1e-10));
        int first_refinement=iteration==(initial&&initial_multipliers?4:24);
        int improved_start=(iteration+1)%50==0&&merit<.5*last_refinement_merit;
        if((first_refinement||improved_start)&&nr&&stats->iterations<(size_t)limit) {
            last_refinement_merit=merit;
            fprintf(stderr,"[local global refinement start] iteration=%d merit=%.9g first=%d\n",iteration+1,merit,first_refinement);
            Sparse_factor_free(&factor);
            Arena_Mark refinement_mark=Arena_save(arena);
            int refined=slg_refine(arena,n,nt,row,col,h,b,nr,original,c,row_scale,x,lambda,limit,stats);
            Arena_restore(arena,refinement_mark);
            if(refined==0) {
                memcpy(out,x,n*sizeof *out);memcpy(out_multipliers,lambda,nr*sizeof *lambda);
                status=0;goto done;
            }
            if(refined<0){status=-1;goto done;}
            /* An unfinished convex refinement remains a useful warm start.
             * Reinitialize local copies before continuing the SPD iteration. */
            for(size_t r=0;r<nr;r++) {
                z[r]=fmax(c[r].lower,slg_dot(c+r,x));u[r]=-lambda[r]/(rho*row_scale[r]);
            }
            /* The residuals preceding refinement describe its old start.
             * Recompute them before any ADMM penalty adjustment. */
            continue;
        }
        if((iteration+1)%50==0&&nr) {
            double dual=0,dual_scale=1;
            for(size_t k=0;k<n;k++) {
                dual=fmax(dual,fabs(work[k])*scale[k]);dual_scale=fmax(dual_scale,fabs(b[k])*scale[k]);
            }
            double ratio=sqrt((primal/primal_scale)/fmax(1e-30,dual/dual_scale));
            double next=fmin(1e6,fmax(1e-6,rho*fmin(10,fmax(.1,ratio))));
            if(next>rho*5||next<rho/5) {
                for(size_t r=0;r<nr;r++)u[r]*=rho/next;
                rho=next;stats->penalty_changes++;Sparse_factor_free(&factor);
            }
        }
        if((iteration+1)%1000==0)fprintf(stderr,"[local global] iteration=%d rho=%.6g primal=%.6g stationarity=%.6g complementarity=%.6g\n",
            iteration+1,rho,stats->primal_residual,stats->stationarity_residual,stats->complementarity);
    }
done:
    Sparse_factor_free(&factor);
    return status;
}

#include "sparse_local_global_test.inc"
