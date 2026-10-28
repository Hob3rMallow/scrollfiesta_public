/* A constrained fit when releasing Dirichlet samples cannot preserve metric.
 * The unknowns are the SAME continuous bilinear map, not independent chart
 * coordinates. Local projections constrain every original triangle; a shared
 * TAUCS factor couples those constraints back to the requested fibre map. */
typedef struct FmStencil {
    int node[12], count;
    double x[12], y[12], base[4], dual[4];
    double weight, pressure;
} FmStencil;

static void fm_stencil(const FmSample *sample,const int32_t face[3],const FmMetric *metric,
                       size_t cols,const double *base,FmStencil *s)
{
    memset(s,0,sizeof *s);
    double gx[3]={-1/metric->length,1/metric->length,0};
    double gy[3]={(metric->along/metric->length-1)/metric->height,
                 -metric->along/(metric->length*metric->height),1/metric->height};
    for(int k=0;k<3;k++){
        const FmSample *p=sample+(size_t)face[k];
        size_t node[4]={p->node,p->node+1,p->node+cols,p->node+cols+1};
        double weight[4]={(1-p->x)*(1-p->y),p->x*(1-p->y),(1-p->x)*p->y,p->x*p->y};
        for(int j=0;j<4;j++){
            if(weight[j]==0)continue;
            int q=0;while(q<s->count && s->node[q]!=(int)node[j])q++;
            if(q==s->count)s->node[s->count++]=(int)node[j];
            s->x[q]+=weight[j]*gx[k];s->y[q]+=weight[j]*gy[k];
        }
    }
    /* Use differences from one node, avoiding cancellation in large canvases. */
    s->x[0]=s->y[0]=0;
    for(int k=1;k<s->count;k++){s->x[0]-=s->x[k];s->y[0]-=s->y[k];}
    size_t origin=2*(size_t)s->node[0];
    for(int k=1;k<s->count;k++){
        size_t q=2*(size_t)s->node[k];double u=base[q]-base[origin],v=base[q+1]-base[origin+1];
        s->base[0]+=s->x[k]*u;s->base[1]+=s->y[k]*u;
        s->base[2]+=s->x[k]*v;s->base[3]+=s->y[k]*v;
    }
    s->weight=.5*metric->length*metric->height;s->pressure=1;
}

/* Closest signed SVD projection. The permitted orientation is the source's,
 * and the interval is the same original-face contract used by the final audit.
 * Use the closed source bounds: a feasible face is a fixed point, including
 * a pre-existing distortion that the shared map must merely not worsen. */
static int fm_project(const double input[4],const FmMetric *metric,int reversed,double out[4])
{
    double sign=reversed?-1:1,a=input[0],b=input[1],c=sign*input[2],d=sign*input[3];
    if(!isfinite(a) || !isfinite(b) || !isfinite(c) || !isfinite(d))return -1;
    double angle=.5*atan2(2*(a*b+c*d),a*a+c*c-b*b-d*d),vx=cos(angle),vy=sin(angle);
    double ux=a*vx+b*vy,uy=c*vx+d*vy,hi=hypot(ux,uy),lo=0;
    if(hi>0){ux/=hi;uy/=hi;lo=(-uy)*(-a*vy+b*vx)+ux*(-c*vy+d*vx);}
    else{ux=1;uy=0;vx=1;vy=0;}
    double margin=FIBER_MAP_PROJECT_INTERIOR*(metric->upper-metric->lower);
    double lower=metric->lower+margin,upper=metric->upper-margin;
    hi=fmax(lower,fmin(upper,hi));lo=fmax(lower,fmin(upper,lo));
    out[0]=hi*ux*vx+lo*uy*vy;out[1]=hi*ux*vy-lo*uy*vx;
    out[2]=sign*(hi*uy*vx-lo*ux*vy);out[3]=sign*(hi*uy*vy+lo*ux*vx);
    return 0;
}

/* Convex inner part of the physical interval, in the baseline polar frame Q:
 * sym(Q^T J) >= lower I and ||Q^T J||_2 <= upper. The first inequality
 * implies sigma_min(J) >= lower, and fixes orientation along every segment
 * from the feasible baseline. Unlike signed-SVD clipping, this local set is
 * convex. It is conservative about rotation; the final physical audit remains
 * authoritative. upper == 0 requests only the orientation half of this set.
 *
 * Write Q^T J = [p+r, s-q; s+q, p-r]. Its operator norm is
 * hypot(p,q)+hypot(r,s), and its symmetric minimum eigenvalue is p-hypot(r,s).
 * Orthogonal projection preserves the direction of (r,s). First try each
 * individual cone projection. If neither is feasible in the intersection,
 * both boundaries are active: p=lower+t and q^2=(upper+lower)(upper-lower-2t).
 * The remaining squared-distance derivative is monotone in t, so bisection
 * finds the unique minimizer without a nested matrix solver. */
static int fm_project_inner(const double input[4],const double base[4],
                            double lower,double upper,int reversed,double out[4])
{
    double sign=reversed?-1:1,cs=base[0]+sign*base[3],sn=sign*base[2]-base[1];
    double norm=hypot(cs,sn);
    if(!(lower>0) || !isfinite(lower) || !(norm>0) || !isfinite(norm) ||
       !isfinite(upper) || (upper!=0 && upper<lower))return -1;
    for(int k=0;k<4;k++)if(!isfinite(input[k]))return -1;
    cs/=norm;sn/=norm;
    double a=cs*input[0]+sn*sign*input[2],b=cs*input[1]+sn*sign*input[3];
    double c=-sn*input[0]+cs*sign*input[2],d=-sn*input[1]+cs*sign*input[3];
    double p=.5*(a+d),q=.5*(c-b),r=.5*(a-d),s=.5*(b+c),t=hypot(r,s);
    double pp=p,qq=q,rr=r,ss=s,z=p-lower;
    if(z<t){
        double next=z<=-t?0:.5*(z+t),scale=t>0?next/t:0;
        pp=lower+next;rr=scale*r;ss=scale*s;
    }
    if(upper==0 || hypot(pp,qq)+hypot(rr,ss)<=upper)goto rotate;

    double h=hypot(p,q),hh=h,tt=t;
    if(h+t>upper){
        if(h-t>=upper){hh=upper;tt=0;}
        else if(t-h>=upper){hh=0;tt=upper;}
        else{double shift=.5*(h+t-upper);hh=h-shift;tt=t-shift;}
    }
    pp=h>0?p*(hh/h):0;qq=h>0?q*(hh/h):0;
    rr=t>0?r*(tt/t):0;ss=t>0?s*(tt/t):0;
    if(pp>=lower+tt)goto rotate;

    double low=0,high=.5*(upper-lower),sum=upper+lower,aq=fabs(q);
    if(t==0 || high==0)tt=0;
    else if(aq==0)tt=fmax(0,fmin(high,.5*(p+t+upper)));
    else{
        for(int step=0;step<DBL_MANT_DIG;step++){
            double middle=.5*(low+high),height=sqrt(fmax(0,sum*(upper-lower-2*middle)));
            double derivative=height>0?2*middle-p-t-upper+sum*aq/height:INFINITY;
            if(derivative>0)high=middle;else low=middle;
        }
        tt=.5*(low+high);
    }
    pp=lower+tt;qq=copysign(sqrt(fmax(0,sum*(upper-lower-2*tt))),q);
    rr=t>0?r*(tt/t):0;ss=t>0?s*(tt/t):0;
rotate:
    a=pp+rr;b=ss-qq;c=ss+qq;d=pp-rr;
    out[0]=cs*a-sn*c;out[1]=cs*b-sn*d;
    out[2]=sign*(sn*a+cs*c);out[3]=sign*(sn*b+cs*d);
    for(int k=0;k<4;k++)if(!isfinite(out[k]))return -1;
    return 0;
}

static void fm_stencil_j(const FmStencil *s,const double *delta,size_t nn,double j[4])
{
    memcpy(j,s->base,sizeof s->base);
    for(int k=0;k<s->count;k++){
        size_t q=(size_t)s->node[k];
        j[0]+=s->x[k]*delta[q];j[1]+=s->y[k]*delta[q];
        j[2]+=s->x[k]*delta[nn+q];j[3]+=s->y[k]*delta[nn+q];
    }
}

/* The four corner Jacobians certify a bilinear cell, including cells with no
 * source triangles. Each corner uses the same map coefficients as the faces.
 * Quadrature weights sum to one cell's area; these are orientation constraints,
 * not additional physical strain limits on empty canvas. */
static void fm_cell_stencil(size_t corner,const FiberMapResult *map,const double *base,FmStencil *s)
{
    size_t cells_x=(size_t)map->cols-1,cell=corner/4;
    size_t a=(cell/cells_x)*(size_t)map->cols+cell%cells_x,b=a+1,c=a+(size_t)map->cols,d=c+1;
    size_t nodes[4][3]={{a,b,c},{a,b,d},{a,c,d},{b,c,d}};
    const double gx[4][3]={{-1,1,0},{-1,1,0},{0,-1,1},{0,-1,1}};
    const double gy[4][3]={{-1,0,1},{0,-1,1},{-1,1,0},{-1,0,1}};
    size_t t=corner%4;memset(s,0,sizeof *s);s->count=3;
    for(int k=0;k<3;k++){
        s->node[k]=(int)nodes[t][k];s->x[k]=gx[t][k]/map->cell;s->y[k]=gy[t][k]/map->cell;
    }
    for(int k=1;k<3;k++){
        size_t q=2*(size_t)s->node[k],origin=2*(size_t)s->node[0];
        double u=base[q]-base[origin],v=base[q+1]-base[origin+1];
        s->base[0]+=s->x[k]*u;s->base[1]+=s->y[k]*u;
        s->base[2]+=s->x[k]*v;s->base[3]+=s->y[k]*v;
    }
    s->weight=.25*map->cell*map->cell;s->pressure=1;
}

static int fm_cell_path(const double base[4],const double j[4],double alpha)
{
    const int32_t triangle[3]={0,1,2};
    double from[6]={0,0,base[0],base[2],base[1],base[3]};
    double to[6]={0,0,j[0],j[2],j[1],j[3]};
    return UvGuard_interval(from,to,triangle,alpha);
}

/* The same convex orientation cone at all four corners protects a whole
 * bilinear cell and its path from the affine baseline. Empty canvas receives
 * no upper strain bound. The independent whole-path guard still checks output. */
static int fm_project_cell(const double input[4],const FmStencil *s,double out[4])
{
    double scale=hypot(hypot(s->base[0],s->base[1]),hypot(s->base[2],s->base[3]));
    /* The path guard bounds the sum of all determinant-polynomial terms,
     * which can be four times the baseline determinant near collapse. Leave
     * arithmetic room for that check, not just a positive endpoint. */
    return fm_project_inner(input,s->base,sqrt(512*DBL_EPSILON)*scale,0,0,out);
}

static double fm_target_error_values(const FiberMapGrid *grid,const FiberMapResult *map,const double *uv)
{
    double mean[2]={0,0},sum=0;size_t count=0;int pad=FIBER_MAP_PAD_CELLS;
    for(int y=0;y<grid->rows;y++)for(int x=0;x<grid->cols;x++){
        size_t p=(size_t)y*(size_t)grid->cols+(size_t)x;
        if(!isfinite(grid->u[p]) || !isfinite(grid->v[p]))continue;
        size_t q=(size_t)(y+pad)*(size_t)map->cols+(size_t)(x+pad);
        mean[0]+=uv[2*q]-grid->u[p];mean[1]+=uv[2*q+1]-grid->v[p];count++;
    }
    if(!count)return 0;
    mean[0]/=(double)count;mean[1]/=(double)count;
    /* Translation fixes the coordinate gauge, but does not align fibres.
     * Two passes avoid subtracting nearly equal squared norms. */
    for(int y=0;y<grid->rows;y++)for(int x=0;x<grid->cols;x++){
        size_t p=(size_t)y*(size_t)grid->cols+(size_t)x;
        if(!isfinite(grid->u[p]) || !isfinite(grid->v[p]))continue;
        size_t q=(size_t)(y+pad)*(size_t)map->cols+(size_t)(x+pad);
        double u=uv[2*q]-grid->u[p]-mean[0],v=uv[2*q+1]-grid->v[p]-mean[1];sum+=u*u+v*v;
    }
    return sum;
}

static double fm_target_error(const FiberMapGrid *grid,const FiberMapResult *map)
{
    return fm_target_error_values(grid,map,map->uv);
}

/* Releasing all but one target can produce a rigid map plus translation.
 * Physical feasibility alone must not call that successful nonrigid alignment.
 * Apply the same useful-fit contract as the constrained fallback. An exact
 * affine target needs no improvement over its already exact baseline. */
static int fm_harmonic_useful(const FiberMapGrid *grid,const FiberMapResult *map,
                              const double *base,FiberMapStats *stats)
{
    stats->metric_fit_before=fm_target_error_values(grid,map,base);
    stats->metric_fit_after=fm_target_error(grid,map);
    if(!isfinite(stats->metric_fit_before) || !isfinite(stats->metric_fit_after))return 0;
    if(!(stats->metric_fit_before>0))return !(stats->metric_fit_after>0);
    return 1-stats->metric_fit_after/stats->metric_fit_before>=FIBER_MAP_PROJECT_MIN_GAIN;
}

/* Independent residual of the supplied lower-triangle COO system, including
 * duplicate entries. Compensated accumulation avoids mistaking cancellation
 * in high-pressure rows for a failed solve. Two RHS use the factor API's
 * column-major layout. This audit never changes the fitted coefficients. */
static void fm_linear_residual(size_t nn,size_t nt,const int *row,const int *col,
                                const double *value,const double *value_v,const double *rhs,const double *delta,
                                double *residual,double *compensation)
{
    memcpy(residual,rhs,2*nn*sizeof *residual);memset(compensation,0,2*nn*sizeof *compensation);
    for(size_t k=0;k<nt;k++)for(size_t d=0;d<2;d++){
        size_t a=d*nn+(size_t)row[k],b=d*nn+(size_t)col[k];
        double coefficient=d && value_v?value_v[k]:value[k];
        double add=-coefficient*delta[b]-compensation[a],next=residual[a]+add;
        compensation[a]=(next-residual[a])-add;residual[a]=next;
        if(a!=b){
            add=-coefficient*delta[a]-compensation[b];next=residual[b]+add;
            compensation[b]=(next-residual[b])-add;residual[b]=next;
        }
    }
}

static int fm_check_linear_solve(Arena_T arena,SparseFactor_T factor,SparseFactor_T factor_v,size_t nn,size_t nt,
                                  const int *row,const int *col,const double *value,const double *value_v,
                                  const double *rhs,const double *delta,
                                  const FmStencil *stencil,size_t nf,int iteration)
{
    Arena_Mark mark=Arena_save(arena);int rc=-1;
    double *residual=ARENA_ALLOC(arena,2*nn*sizeof *residual);
    double *correction=ARENA_ALLOC(arena,2*nn*sizeof *correction);
    fm_linear_residual(nn,nt,row,col,value,value_v,rhs,delta,residual,correction);
    double maximum=0,rhs_maximum=0,coefficient=0,jacobian=0;size_t worst_face=0;
    for(size_t k=0;k<2*nn;k++){
        maximum=fmax(maximum,fabs(residual[k]));rhs_maximum=fmax(rhs_maximum,fabs(rhs[k]));
    }
    if(factor_v){
        if(Sparse_factor_solve(factor,residual,correction) ||
           Sparse_factor_solve(factor_v,residual+nn,correction+nn))goto done;
    }else if(Sparse_factor_solve_multi(factor,residual,correction,2))goto done;
    for(size_t k=0;k<2*nn;k++){
        if(!isfinite(correction[k]))goto done;
        coefficient=fmax(coefficient,fabs(correction[k]));
    }
    for(size_t f=0;f<nf;f++){
        const FmStencil *s=stencil+f;double j[4]={0};size_t origin=(size_t)s->node[0];
        for(int k=1;k<s->count;k++){
            size_t q=(size_t)s->node[k];
            double u=correction[q]-correction[origin],v=correction[nn+q]-correction[nn+origin];
            j[0]+=s->x[k]*u;j[1]+=s->y[k]*u;j[2]+=s->x[k]*v;j[3]+=s->y[k]*v;
        }
        double norm=hypot(hypot(j[0],j[1]),hypot(j[2],j[3]));
        if(norm>jacobian){jacobian=norm;worst_face=f;}
    }
    fprintf(stderr,"[fiber_map linear audit] iteration %d: residual_inf %.12g relative %.12g; discarded correction coefficient %.12g Jacobian %.12g at face %zu\n",
            iteration,maximum,maximum/fmax(DBL_MIN,rhs_maximum),coefficient,jacobian,worst_face);fflush(stderr);
    rc=0;
done:
    Arena_restore(arena,mark);return rc;
}

#include "fiber_map_direction.inc"

/* Center the original quadratic DATA metric at the current iterate. This is
 * a matrix-vector product, not a new solve: the existing TAUCS factor and all
 * physical constraints stay in force. It removes the distant target's pull
 * during the final feasibility projection. Final quality is still measured
 * against the ORIGINAL observations, never this temporary center. */
static void fm_polish_center(size_t nn,size_t entries,const int *row,const int *col,
                              const double *value,const double *value_v,const double *delta,
                              double *prior,FmStencil *stencil,size_t constraints)
{
    memset(prior,0,2*nn*sizeof *prior);
    for(size_t k=0;k<entries;k++)for(size_t d=0;d<2;d++){
        size_t a=d*nn+(size_t)row[k],b=d*nn+(size_t)col[k];
        double weight=d && value_v?value_v[k]:value[k];
        prior[a]+=weight*delta[b];
        if(a!=b)prior[b]+=weight*delta[a];
    }
    for(size_t f=0;f<constraints;f++)memset(stencil[f].dual,0,sizeof stencil[f].dual);
}

static int fm_metric_fit(Arena_T arena,const FiberMapGrid *grid,FiberMapResult *map,
                         const double *base,const double *safe,const FmSample *sample,
                         const FmMetric *metric,const uint8_t *reversed,const int32_t *faces,
                         size_t nv,size_t nf,double *output,FiberMapStats *stats)
{
    Arena_Mark mark=Arena_save(arena);SparseFactor_T factor=NULL,factor_v=NULL;int rc=-1;
    double seconds=grid->max_seconds>0?grid->max_seconds:FIBER_MAP_PROJECT_SECONDS;
    double started=ves_clock_sec(),deadline=started+seconds;
    size_t nn=(size_t)map->rows*(size_t)map->cols,nt=nn;
    size_t nc=4*(size_t)(map->rows-1)*(size_t)(map->cols-1),constraints=0;
    if(!nf || nn>(size_t)INT_MAX/3 || nc>SIZE_MAX-nf)goto done;
    FmDirectionStencil *direction=NULL;size_t nd=0;
    if(fm_direction_stencils(arena,grid->directions,map,base,&direction,&nd))goto done;
    stats->direction_fit=nd>0;stats->direction_samples=grid->directions?grid->directions->count:0;
    for(size_t i=0;i<nd;i++){
        size_t n=(size_t)direction[i].count,add=n*(n+1)/2;
        if(add>(size_t)INT_MAX-nt)goto done;nt+=add;
    }
    constraints=nf+nc;
    if(constraints>SIZE_MAX/sizeof(FmStencil))goto done;
    FmStencil *stencil=ARENA_ALLOC(arena,constraints*sizeof *stencil);
    /* Direction residuals and face Jacobians are both dimensionless, with
     * physical-area quadrature. The legacy nodal-position objective instead
     * needs cell^2 to balance its distance-squared units. Retaining that factor
     * for directions changes their relative pressure under a canvas rescale. */
    double rho=FIBER_MAP_PROJECT_PENALTY*(nd?1:map->cell*map->cell);
    for(size_t f=0;f<constraints;f++){
        if(f<nf)fm_stencil(sample,faces+3*f,metric+f,(size_t)map->cols,base,stencil+f);
        else fm_cell_stencil(f-nf,map,base,stencil+f);
        size_t n=(size_t)stencil[f].count,add=n*(n+1)/2;
        if(add>(size_t)INT_MAX-nt)goto done;nt+=add;
    }
    size_t edges=(size_t)map->rows*(size_t)(map->cols-1)+(size_t)(map->rows-1)*(size_t)map->cols;
    if(edges>((size_t)INT_MAX-nt)/3)goto done;nt+=3*edges;
    fprintf(stderr,"[fiber_map metric] %zu nodes, %zu original faces, %zu cell corners, %zu matrix entries; %s physical projection; rho %.9g; TAUCS factor\n",
            nn,nf,nc,nt,nd?"signed SVD":"convex polar",rho);fflush(stderr);
    int *row=ARENA_ALLOC(arena,nt*sizeof *row),*col=ARENA_ALLOC(arena,nt*sizeof *col);
    double *value=ARENA_ALLOC(arena,nt*sizeof *value),*prior=ARENA_CALLOC(arena,2*nn,sizeof *prior);
    double *value_v=value;
    double *rhs=ARENA_ALLOC(arena,2*nn*sizeof *rhs),*delta=ARENA_CALLOC(arena,2*nn,sizeof *delta);
    double *best=ARENA_ALLOC(arena,2*nn*sizeof *best),best_error=INFINITY;
    size_t at=0;int pad=FIBER_MAP_PAD_CELLS;
    for(int y=0;y<map->rows;y++)for(int x=0;x<map->cols;x++){
        size_t q=(size_t)y*(size_t)map->cols+(size_t)x;double weight=FIBER_MAP_PROJECT_PRIOR;
        if(!nd && y>=pad && y<map->rows-pad && x>=pad && x<map->cols-pad){
            size_t p=(size_t)(y-pad)*(size_t)grid->cols+(size_t)(x-pad);
            if(isfinite(grid->u[p]) && isfinite(grid->v[p])){
                weight+=1;prior[q]=grid->u[p]-base[2*q];prior[nn+q]=grid->v[p]-base[2*q+1];
            }
        }
        row[at]=col[at]=(int)q;value[at++]=weight;
        size_t neighbor[2]={q+1,q+(size_t)map->cols};int have[2]={x+1<map->cols,y+1<map->rows};
        for(int d=0;d<2;d++)if(have[d]){
            int a=(int)q,b=(int)neighbor[d];double w=FIBER_MAP_PROJECT_SMOOTH;
            row[at]=col[at]=a;value[at++]=w;row[at]=col[at]=b;value[at++]=w;
            row[at]=b;col[at]=a;value[at++]=-w;
        }
    }
    if(nd){
        double alpha=FIBER_MAP_DIRECTION_LENGTH_WEIGHT,beta=FIBER_MAP_DIRECTION_PENALTY;
        if(!(alpha>0) || !isfinite(alpha) || !(beta>0) || !isfinite(beta))goto done;
        for(size_t k=0;k<nd;k++){
            const FmDirectionStencil *s=direction+k;
            for(int a=0;a<s->count;a++){
                for(int b=0;b<=a;b++){
                    int i=s->node[a],j=s->node[b];row[at]=i>j?i:j;col[at]=i>j?j:i;
                    value[at++]=beta*s->weight*(s->h[a]*s->h[b]+s->v[a]*s->v[b]);
                }
            }
        }
        fprintf(stderr,"[fiber_map direction] %zu original observations, %zu incident-cell derivatives; quarter-turn cross; soft length weight %.6g, proximal penalty %.6g\n",
                stats->direction_samples,nd,alpha,beta);fflush(stderr);
    }
    size_t face_entries=at;
    for(size_t f=0;f<constraints;f++){
        const FmStencil *s=stencil+f;
        for(int a=0;a<s->count;a++)for(int b=0;b<=a;b++){
            int i=s->node[a],j=s->node[b];row[at]=i>j?i:j;col[at]=i>j?j:i;
            value[at]=rho*s->weight*s->pressure*(s->x[a]*s->x[b]+s->y[a]*s->y[b]);
            value_v[at]=value[at];at++;
        }
    }
    if(at!=nt || Sparse_factor_spd((int)nn,(int)nt,row,col,value,&factor))goto done;
    memcpy(map->uv,base,2*nn*sizeof *base);
    stats->metric_fit_before=fm_target_error(grid,map);
    double error_before=nd?fm_direction_error(direction,nd,map->uv):stats->metric_fit_before;
    if(!isfinite(error_before))goto done;
    if(nd)stats->direction_angle_before=error_before;
    size_t baseline_bad=0;double baseline_worst=0;
    for(size_t v=0;v<nv;v++)fm_value(map,sample+v,output+2*v);
    for(size_t f=0;f<nf;f++){
        double lo,hi;fm_metric(metric+f,faces+3*f,output,1,1,&lo,&hi);
        baseline_bad+=(size_t)!fm_metric_preserved(metric+f,lo,hi);
        double violation=fmax(metric[f].lower-lo,hi-metric[f].upper);
        baseline_worst=fmax(baseline_worst,violation);
    }
    fprintf(stderr,"[fiber_map metric] sampled baseline: %zu failures, maximum excess %.12g\n",baseline_bad,baseline_worst);fflush(stderr);
    uint8_t *seed=ARENA_ALLOC(arena,nn);
    double previous=INFINITY,next_log=started,gain=0,worst=0,sum=0;int stalled=0,polished=0;
    size_t bad=0,orientation=0,bad_old=0,worst_face=0,cells=0;
    while(ves_clock_sec()<deadline){
        memcpy(rhs,prior,2*nn*sizeof *rhs);
        if(nd && !polished)fm_direction_rhs(direction,nd,nn,delta,rhs);
        for(size_t f=0;f<constraints;f++){
            FmStencil *s=stencil+f;double j[4],z[4];fm_stencil_j(s,delta,nn,j);
            for(int k=0;k<4;k++)j[k]+=s->dual[k];
            if(f<nf){
                /* The physical band is invariant under output rotation. A
                 * fixed polar frame would penalize a pure rotation at an
                 * active source bound. Direct alignment needs that freedom;
                 * independently audit the coupled result before accepting it. */
                int projected=nd?fm_project(j,metric+f,reversed[f],z):
                    fm_project_inner(j,s->base,metric[f].lower,metric[f].upper,reversed[f],z);
                if(projected)goto done;
            }
            else if(fm_project_cell(j,s,z))goto done;
            for(int k=0;k<4;k++){s->dual[k]=j[k]-z[k];z[k]-=s->dual[k]+s->base[k];}
            for(int k=0;k<s->count;k++){
                size_t q=(size_t)s->node[k];double w=rho*s->weight*s->pressure;
                rhs[q]+=w*(s->x[k]*z[0]+s->y[k]*z[1]);
                rhs[nn+q]+=w*(s->x[k]*z[2]+s->y[k]*z[3]);
            }
        }
        if(Sparse_factor_solve_multi(factor,rhs,delta,2))goto done;
        if((stats->metric_iterations==0 || (stats->metric_iterations+1)%128==0) &&
           fm_check_linear_solve(arena,factor,factor_v,nn,nt,row,col,value,value_v,rhs,delta,stencil,nf,stats->metric_iterations+1))goto done;
        for(size_t q=0;q<nn;q++){
            if(!isfinite(delta[q]) || !isfinite(delta[nn+q]))goto done;
            map->uv[2*q]=base[2*q]+delta[q];map->uv[2*q+1]=base[2*q+1]+delta[nn+q];
        }
        for(size_t v=0;v<nv;v++)fm_value(map,sample+v,output+2*v);
        size_t pressured=0;bad=orientation=bad_old=worst_face=0;worst=sum=0;
        for(size_t f=0;f<nf;f++){
            double lo,hi;fm_metric(metric+f,faces+3*f,output,1,1,&lo,&hi);
            if(!fm_metric_preserved(metric+f,lo,hi)){
                bad++;bad_old+=(size_t)(metric[f].lower<FIBER_MAP_METRIC_MIN-ASM_METRIC_ROUNDOFF ||
                                      metric[f].upper>FIBER_MAP_METRIC_MAX+ASM_METRIC_ROUNDOFF);
                double violation=fmax(metric[f].lower-lo,hi-metric[f].upper);
                sum+=violation*violation;
                if(violation>worst){worst=violation;worst_face=f;}
                if((stats->metric_iterations+1)%FIBER_MAP_PROJECT_REWEIGHT_STEPS==0 &&
                   stencil[f].pressure<FIBER_MAP_PROJECT_PRESSURE_MAX){
                    FmStencil *s=stencil+f;
                    double next=fmin(FIBER_MAP_PROJECT_PRESSURE_MAX,s->pressure*FIBER_MAP_PROJECT_PRESSURE_GAIN);
                    /* Preserve the unscaled multiplier when changing rho. */
                    for(int k=0;k<4;k++)s->dual[k]*=s->pressure/next;
                    s->pressure=next;pressured++;
                }
            }
            int32_t idx[3]={faces[3*f],faces[3*f+(reversed[f]?2:1)],faces[3*f+(reversed[f]?1:2)]};
            orientation+=(size_t)!UvGuard_interval(safe,output,idx,1);
        }
        size_t cell_pressured=0;
        if((stats->metric_iterations+1)%FIBER_MAP_PROJECT_REWEIGHT_STEPS==0){
            for(size_t f=nf;f<constraints;f++){
                FmStencil *s=stencil+f;double j[4];fm_stencil_j(s,delta,nn,j);
                if(fm_cell_path(s->base,j,1) || s->pressure>=FIBER_MAP_PROJECT_PRESSURE_MAX)continue;
                double next=fmin(FIBER_MAP_PROJECT_PRESSURE_MAX,s->pressure*FIBER_MAP_PROJECT_PRESSURE_GAIN);
                for(int k=0;k<4;k++)s->dual[k]*=s->pressure/next;
                s->pressure=next;cell_pressured++;
            }
        }
        memset(seed,0,nn);cells=fm_cells(map,base,seed);
        stats->metric_iterations++;stats->metric_faces_after=bad;
        stats->metric_fit_after=fm_target_error(grid,map);
        double error_after=nd?fm_direction_error(direction,nd,map->uv):stats->metric_fit_after;
        if(nd)stats->direction_angle_after=error_after;
        gain=fm_fit_gain(error_before,error_after,nd>0);
        if(ves_clock_sec()>=next_log){
            fprintf(stderr,"[fiber_map metric] iteration %d: %zu metric failures, %zu triangle paths, %zu cells; fit gain %.6f; %.3f s\n",
                    stats->metric_iterations,bad,orientation,cells,gain,ves_clock_sec()-started);
            fprintf(stderr,"[fiber_map metric] excess: %zu originally bad, max %.12g at face %zu, rms %.12g\n",
                    bad_old,worst,worst_face,bad?sqrt(sum/(double)bad):0);fflush(stderr);
            if(nd)fprintf(stderr,"[fiber_map direction] angle RMS degrees %.8g -> %.8g (both families, weighted incident cells)\n",
                          sqrt(error_before)*180/3.141592653589793,sqrt(error_after)*180/3.141592653589793);
            next_log=ves_clock_sec()+FIBER_MAP_PROJECT_LOG_SECONDS;
        }
        if(!bad && !orientation && !cells && gain>=FIBER_MAP_PROJECT_MIN_GAIN && error_after<best_error){
            best_error=error_after;memcpy(best,map->uv,2*nn*sizeof *best);
            if(grid->checkpoint) {
                for(size_t v=0;v<nv;v++)fm_value(map,sample+v,output+2*v);
                int decision=grid->checkpoint(grid->checkpoint_context,output,nv,stats->metric_iterations);
                if(decision<0)goto done;
                if(decision>0)break;
            }
        }
        if((pressured || cell_pressured) && ves_clock_sec()<deadline){
            at=face_entries;
            for(size_t f=0;f<constraints;f++){
                const FmStencil *s=stencil+f;
                for(int a=0;a<s->count;a++)for(int b=0;b<=a;b++){
                    value[at]=rho*s->weight*s->pressure*(s->x[a]*s->x[b]+s->y[a]*s->y[b]);
                    value_v[at]=value[at];at++;
                }
            }
            Sparse_factor_free(&factor);
            Sparse_factor_free(&factor_v);
            fprintf(stderr,"[fiber_map metric] increasing constraint pressure on %zu failing original faces and %zu cell corners\n",pressured,cell_pressured);fflush(stderr);
            if(at!=nt || Sparse_factor_spd((int)nn,(int)nt,row,col,value,&factor))goto done;
        }
        /* Measure angular stagnation relative to the initial error, rather
         * than demanding a vanishing tolerance as alignment approaches zero.
         * This changes only stopping AFTER a physically certified useful map
         * exists. Arithmetic-zero angles must also be able to stop. */
        double stall_scale=nd?fmax(error_before,previous):fmax(1,previous);
        double stall_tolerance=FIBER_MAP_PROJECT_STALL_REL*stall_scale;
        if(nd)stall_tolerance=fmax(256*DBL_EPSILON*DBL_EPSILON,stall_tolerance);
        stalled=isfinite(previous) && fabs(previous-error_after)<=stall_tolerance?stalled+1:0;
        previous=error_after;
        if(!polished && !isfinite(best_error) && !orientation && !cells &&
           gain>=FIBER_MAP_PROJECT_MIN_GAIN &&
           ves_clock_sec()>=deadline-FIBER_MAP_PROJECT_POLISH_SHARE*seconds){
            fm_polish_center(nn,face_entries,row,col,value,value_v,delta,prior,stencil,constraints);
            polished=1;stalled=0;previous=INFINITY;
            fprintf(stderr,"[fiber_map metric] feasibility polish centered on useful iterate: original gain %.9g, %zu metric failures; hard bounds unchanged\n",
                    gain,bad);fflush(stderr);
        }
        /* A stationary target objective does not imply primal feasibility.
         * In particular, nearly active original-face bounds can keep improving
         * after the objective has settled. Keep the existing time budget unless
         * a fully certified, useful map is already available to return. */
        if(!grid->max_iterations&&stalled>=FIBER_MAP_PROJECT_STALL_STEPS && isfinite(best_error))break;
        if(grid->max_iterations>0&&stats->metric_iterations>=grid->max_iterations)break;
        if(stats->metric_iterations==INT_MAX)break;
    }
    fprintf(stderr,"[fiber_map metric] iteration %d: %zu metric failures, %zu triangle paths, %zu cells; fit gain %.6f; %.3f s\n",
            stats->metric_iterations,bad,orientation,cells,gain,ves_clock_sec()-started);
    fprintf(stderr,"[fiber_map metric] excess: %zu originally bad, max %.12g at face %zu, rms %.12g\n",
            bad_old,worst,worst_face,bad?sqrt(sum/(double)bad):0);fflush(stderr);
    if(isfinite(best_error)){
        memcpy(map->uv,best,2*nn*sizeof *best);
        for(size_t v=0;v<nv;v++)fm_value(map,sample+v,output+2*v);
        stats->metric_fit_after=fm_target_error(grid,map);stats->metric_faces_after=0;rc=0;
        if(nd)stats->direction_angle_after=best_error;
        fprintf(stderr,"[fiber_map metric] selected feasible fit: gain %.6f after %d iterations\n",
                fm_fit_gain(error_before,best_error,nd>0),stats->metric_iterations);fflush(stderr);
    }
done:
    Sparse_factor_free(&factor);Sparse_factor_free(&factor_v);Arena_restore(arena,mark);return rc;
}
