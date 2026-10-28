/* Refine a useful local/global start by minimizing a convex augmented
 * least-squares objective. With multipliers fixed, the local squared hinge
 * determines which rows participate in the next SPD global solve. A line
 * search handles changes in that local classification. Updating the local
 * multipliers then enforces agreement without taking a penalty to infinity.
 * This is a method-of-multipliers refinement, not an indefinite polish. */
static int slg_refine(Arena_T arena,size_t n,size_t nt,const int *row,const int *col,
    const double *h,const double *b,size_t nr,const SparseLocalRow *original,
    const SparseLocalRow *c,const double *row_scale,double *x,double *lambda,
    int limit,SparseLocalGlobalStats *stats)
{
    size_t capacity=nt+21*nr;
    int *ar=ARENA_ALLOC(arena,capacity*sizeof *ar),*ac=ARENA_ALLOC(arena,capacity*sizeof *ac);
    double *av=ARENA_ALLOC(arena,capacity*sizeof *av),*scaled=ARENA_ALLOC(arena,capacity*sizeof *scaled);
    double *diag=ARENA_ALLOC(arena,n*sizeof *diag),*scale=ARENA_ALLOC(arena,n*sizeof *scale);
    double *rhs=ARENA_ALLOC(arena,n*sizeof *rhs),*solution=ARENA_ALLOC(arena,n*sizeof *solution);
    double *gradient=ARENA_ALLOC(arena,n*sizeof *gradient),*hx=ARENA_ALLOC(arena,n*sizeof *hx);
    double *p=ARENA_ALLOC(arena,n*sizeof *p),*hp=ARENA_ALLOC(arena,n*sizeof *hp);
    double *residual=ARENA_ALLOC(arena,n*sizeof *residual),*correction=ARENA_ALLOC(arena,n*sizeof *correction);
    double *local=ARENA_ALLOC(arena,nr*sizeof *local),*shift=ARENA_ALLOC(arena,nr*sizeof *shift);
    uint8_t *active=ARENA_ALLOC(arena,nr),*factored=ARENA_CALLOC(arena,nr,1);
    for(size_t r=0;r<nr;r++)local[r]=lambda[r]/row_scale[r];
    double rho=100,previous_violation=HUGE_VAL;
    /* Early multiplier updates do not require a fully minimized penalty
     * problem. Tighten the inner forcing with the outer feasibility/gap
     * residuals; only the original certificate may finish the whole solve. */
    double inner_tolerance=fmin(1e-3,fmax(5e-10,.1*fmax(stats->primal_residual,stats->complementarity)));
    SparseFactor_T factor=NULL;size_t at=0;int status=1;
    for(int outer=0;outer<64;outer++) {
        int minimized=0;
        for(int inner=0;inner<64;inner++) {
            slg_multiply(n,nt,row,col,h,x,hx);
            double residual_max=0,bscale=1;
            for(size_t k=0;k<n;k++){hx[k]-=b[k];gradient[k]=hx[k];bscale=fmax(bscale,fabs(b[k]));}
            for(size_t r=0;r<nr;r++) {
                shift[r]=local[r]+rho*(c[r].lower-slg_dot(c+r,x));active[r]=(uint8_t)(shift[r]>0);
                double q=fmax(0,shift[r]);
                for(size_t k=0;k<6;k++)gradient[c[r].index[k]]-=c[r].value[k]*q;
            }
            for(size_t k=0;k<n;k++)residual_max=fmax(residual_max,fabs(gradient[k]));
            /* The projected local multipliers correspond to this x, even
             * before the stricter inner minimization target is reached. */
            for(size_t r=0;r<nr;r++)lambda[r]=fmax(0,shift[r])*row_scale[r];
            stats->penalty=rho;
            int certified=slg_certificate(n,nt,row,col,h,b,nr,original,x,lambda,hp,stats);
            if(certified<0){status=-1;goto finish;}
            if(certified){stats->converged=1;status=0;goto finish;}
            if(residual_max/bscale<=inner_tolerance){minimized=1;break;}
            if(stats->iterations>=(size_t)limit)goto finish;
            if(!factor||memcmp(active,factored,nr)) {
                Sparse_factor_free(&factor);
                memcpy(ar,row,nt*sizeof *ar);memcpy(ac,col,nt*sizeof *ac);memcpy(av,h,nt*sizeof *av);at=nt;
                for(size_t r=0;r<nr;r++)if(active[r])for(size_t k=0;k<6;k++)if(c[r].value[k]!=0)
                    for(size_t j=0;j<=k;j++)if(c[r].value[j]!=0) {
                        int a=c[r].index[k],d=c[r].index[j];
                        ar[at]=a>d?a:d;ac[at]=a>d?d:a;av[at++]=rho*c[r].value[k]*c[r].value[j];
                    }
                memset(diag,0,n*sizeof *diag);
                for(size_t k=0;k<at;k++)if(ar[k]==ac[k])diag[ar[k]]+=av[k];
                for(size_t k=0;k<n;k++) {
                    if(!(diag[k]>0)||!isfinite(diag[k])){status=-1;goto finish;}
                    scale[k]=1/sqrt(diag[k]);
                }
                for(size_t k=0;k<at;k++)scaled[k]=av[k]*scale[ar[k]]*scale[ac[k]];
                if(Sparse_factor_spd((int)n,(int)at,ar,ac,scaled,&factor)){status=-1;goto finish;}
                stats->factors++;memcpy(factored,active,nr);
            }
            /* Solve for a correction, preserving a useful starting point and
             * avoiding cancellation in the absolute least-squares solution. */
            for(size_t k=0;k<n;k++)rhs[k]=-gradient[k]*scale[k];
            if(Sparse_factor_solve(factor,rhs,solution)){status=-1;goto finish;}
            stats->solves++;
            for(int refine=0;refine<3;refine++) {
                double error=0,norm=0;
                slg_multiply(n,at,ar,ac,scaled,solution,residual);
                for(size_t k=0;k<n;k++) {
                    residual[k]=rhs[k]-residual[k];error=fmax(error,fabs(residual[k]));norm=fmax(norm,fabs(rhs[k]));
                }
                /* This is an inexact descent direction, not the final
                 * optimization certificate. Recompute the original gradient
                 * after the update instead of oversolving a conditioned proxy. */
                double relative=error/(norm>0?norm:1);
                stats->linear_residual=fmax(stats->linear_residual,relative);
                if(!isfinite(relative)){status=-1;goto finish;}
                double forcing=fmin(.1,fmax(1e-6,sqrt(residual_max/bscale)));
                if(relative<=forcing)break;
                if(refine==2) {
                    fprintf(stderr,"[local global linear floor] rho=%.6g relative=%.6g outer=%d inner=%d\n",rho,relative,outer,inner);
                    goto finish;
                }
                if(Sparse_factor_solve(factor,residual,correction)){status=-1;goto finish;}
                stats->solves++;
                for(size_t k=0;k<n;k++)solution[k]+=correction[k];
            }
            double slope=0,linear=0,quadratic=0;
            for(size_t k=0;k<n;k++){p[k]=solution[k]*scale[k];slope+=gradient[k]*p[k];linear+=hx[k]*p[k];}
            if(!isfinite(slope)||!(slope<0))goto finish;
            slg_multiply(n,nt,row,col,h,p,hp);
            for(size_t k=0;k<n;k++)quadratic+=p[k]*hp[k];
            double alpha=1;int accepted=0;
            for(int trial=0;trial<32;trial++,alpha*=.5) {
                double change=alpha*linear+.5*alpha*alpha*quadratic;
                for(size_t r=0;r<nr;r++) {
                    double q=fmax(0,shift[r]),delta=-rho*alpha*slg_dot(c+r,p);
                    double next=fmax(0,shift[r]+delta);
                    double difference=shift[r]>0&&shift[r]+delta>0?delta:next-q;
                    change+=difference*(next+q)/(2*rho);
                }
                if(isfinite(change)&&change<=1e-4*alpha*slope){accepted=1;break;}
            }
            if(!accepted)goto finish;
            for(size_t k=0;k<n;k++)x[k]+=alpha*p[k];
            stats->iterations++;stats->local_refinements++;
        }
        if(!minimized)goto finish;
        double violation=0;
        for(size_t r=0;r<nr;r++) {
            double slack=slg_dot(c+r,x)-c[r].lower;
            local[r]=fmax(0,local[r]-rho*slack);lambda[r]=local[r]*row_scale[r];
            violation=fmax(violation,fabs(slack)*(local[r]>0?1:(slack<0?1:0)));
        }
        int certified=slg_certificate(n,nt,row,col,h,b,nr,original,x,lambda,gradient,stats);
        stats->penalty=rho;
        if(certified<0){status=-1;goto finish;}
        if(certified){stats->converged=1;status=0;goto finish;}
        inner_tolerance=fmin(inner_tolerance,fmax(5e-10,.1*fmax(stats->primal_residual,stats->complementarity)));
        if(violation>.25*previous_violation&&rho<1e8) {
            rho*=10;stats->penalty_changes++;Sparse_factor_free(&factor);
        }
        previous_violation=violation;
    }
finish:
    Sparse_factor_free(&factor);
    fprintf(stderr,"[local global refinement] steps=%zu factors=%zu rho=%.6g primal=%.6g stationarity=%.6g complementarity=%.6g status=%d\n",
        stats->local_refinements,stats->factors,rho,stats->primal_residual,stats->stationarity_residual,stats->complementarity,status);
    return status;
}
