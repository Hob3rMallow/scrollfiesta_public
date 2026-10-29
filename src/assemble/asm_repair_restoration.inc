/* Full original-vertex restoration from c3_seam_restoration.restore_support:
 * original-area symmetric Dirichlet, physical Huber seams and a source-scaled
 * boundary barrier. Incremental variables eliminate the immutable source gauge.
 * Original-face/contact guards decide acceptance independently of the QP. */
typedef struct ArsBoundary ArsBoundary;
typedef struct ArsCoordinates ArsCoordinates;
typedef struct ArsStretch ArsStretch;
typedef struct ArsHash {uint64_t key;int index;} ArsHash;
typedef struct ArsProblem {
    Arena_T arena;AprSolver *s;ArsBoundary *boundary;
    ArsCoordinates *coordinates;ArsHash *hslots;size_t hcapacity;
    ArsStretch *stretch;
    AlqMatrix h,a;int *id,n,m,capacity,derivatives;
    double *gradient,*diagonal,*upper,*answer;
    double bound,energy,violation,stretch_energy;
} ArsProblem;

#include "asm_repair_restoration_coordinates.inc"

static double ars_cross(const double a[2],const double b[2])
{return a[0]*b[1]-a[1]*b[0];}
static void ars_roots(double a,double b,double c,double roots[2])
{
    roots[0]=roots[1]=INFINITY;double scale=fmax(1,fmax(fabs(a),fmax(fabs(b),fabs(c))));
    if(fabs(a)<=16*DBL_EPSILON*scale){if(fabs(b)>16*DBL_EPSILON*scale)roots[0]=-c/b;return;}
    double disc=b*b-4*a*c;if(disc<0)return;
    double q=-.5*(b+(b>=0?1:-1)*sqrt(disc));
    if(q){roots[0]=q/a;roots[1]=c/q;}else roots[0]=-b/(2*a);
}
static void ars_psd4(double h[16])
{
    double v[16]={1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1};
    for(int iteration=0;iteration<96;iteration++){
        int p=0,q=1;double largest=fabs(h[1]),scale=1;
        for(int i=0;i<4;i++){scale=fmax(scale,fabs(h[5*i]));for(int j=i+1;j<4;j++)if(fabs(h[4*i+j])>largest){p=i;q=j;largest=fabs(h[4*i+j]);}}
        if(largest<=8*DBL_EPSILON*scale)break;
        double tau=(h[5*q]-h[5*p])/(2*h[4*p+q]);
        double t=(tau>=0?1:-1)/(fabs(tau)+hypot(1,tau)),c=1/hypot(1,t),s=t*c;
        double app=h[5*p],aqq=h[5*q],apq=h[4*p+q];
        h[5*p]=app-t*apq;h[5*q]=aqq+t*apq;h[4*p+q]=h[4*q+p]=0;
        for(int k=0;k<4;k++){
            if(k!=p && k!=q){double x=h[4*k+p],y=h[4*k+q];h[4*k+p]=h[4*p+k]=c*x-s*y;h[4*k+q]=h[4*q+k]=s*x+c*y;}
            double x=v[4*k+p],y=v[4*k+q];v[4*k+p]=c*x-s*y;v[4*k+q]=s*x+c*y;
        }
    }
    double result[16]={0};
    for(int k=0;k<4;k++)for(int i=0;i<4;i++)for(int j=0;j<4;j++)
        result[4*i+j]+=v[4*i+k]*fmax(h[5*k],1e-7)*v[4*j+k];
    memcpy(h,result,sizeof result);
}
static double ars_sd(const double j[4],double area,double gradient[4],double hessian[16])
{
    double det=j[0]*j[3]-j[1]*j[2],square=0;
    if(!(det>1e-10))return INFINITY;
    for(int k=0;k<4;k++)square+=j[k]*j[k];
    double inv=1/det,d2=inv*inv,d3=d2*inv,d4=d3*inv,c[]={j[3],-j[2],-j[1],j[0]};
    if(gradient)for(int k=0;k<4;k++)gradient[k]=area*(2*j[k]*(1+d2)-2*c[k]*square*d3);
    if(hessian){
        const double hd[16]={0,0,0,1,0,0,-1,0,0,-1,0,0,1,0,0,0};
        for(int i=0;i<4;i++)for(int k=0;k<4;k++)
            hessian[4*i+k]=2*(1+d2)*(i==k)-4*d3*(j[i]*c[k]+c[i]*j[k])+6*square*d4*c[i]*c[k]-2*square*d3*hd[4*i+k];
        ars_psd4(hessian);for(int k=0;k<16;k++)hessian[k]*=area;
    }
    return area*(square*(1+d2)-4);
}

#include "asm_repair_restoration_stretch.inc"
static void ars_h(ArsProblem *p,int row,int col,double value)
{
    if(row<0 || col<0 || row<col || !value)return;
    if(row==col)p->diagonal[row]+=value;
    /* ConvexQp_solve consumes the native doubled normal convention. */
    if(p->coordinates)ars_hashed_entry(p,row,col,2*value);
    else alq_entry(&p->h,row,col,2*value);
}
static int ars_plane(ArsProblem *p,int count,const int *ids,const double *values,double upper)
{
    if(p->m==INT_MAX || !isfinite(upper))return 0;
    if(p->m==p->capacity){
        if(p->capacity>INT_MAX/2)return 0;int capacity=p->capacity?2*p->capacity:1024;
        double *next=realloc(p->upper,(size_t)capacity*sizeof(double));if(!next)return 0;
        p->upper=next;p->capacity=capacity;
    }
    for(int k=0;k<count;k++)if(ids[k]>=0)alq_entry(&p->a,p->m,ids[k],values[k]);
    /* A failed seam needs a nonzero correction: negative upper bounds stay. */
    p->upper[p->m++]=upper;return !p->a.failed;
}
static int ars_incident(const AprSolver *s,const AsmFieldSeam *seam)
{
    return seam->required && seam->complete && (!s->opts.fixed_charts ||
        !s->opts.fixed_charts[seam->a] || !s->opts.fixed_charts[seam->b]);
}

#include "asm_repair_restoration_boundary.inc"

static double ars_energy(ArsProblem *p,const double *uv,int derivatives)
{
    AprSolver *s=p->s;const AsmField *f=s->field;p->derivatives=derivatives;p->energy=p->stretch_energy=0;
    if(derivatives){p->h.n=0;p->h.failed=0;if(p->hslots)memset(p->hslots,0,p->hcapacity*sizeof(ArsHash));memset(p->gradient,0,(size_t)p->n*sizeof(double));memset(p->diagonal,0,(size_t)p->n*sizeof(double));}
    for(size_t at=0;at<s->n_active;at++){
        size_t t=s->active_faces[at];const int32_t *face=f->faces+3*t;double j[4],lo,hi,det,gradient[4],hessian[16],g[6];
        if(!UvGuard_positive(uv,face))return INFINITY;
        AsmField_face(f,t,uv,j,&lo,&hi,&det);
        double energy=ars_sd(j,f->metric[t].area,derivatives?gradient:NULL,derivatives?hessian:NULL);
        if(!isfinite(energy))return INFINITY;p->energy+=energy;
        if(p->stretch){
            double sg[4],sh[16],value=ars_stretch_face(p->stretch+at,j,f->metric[t].area,derivatives?sg:NULL,derivatives?sh:NULL);
            if(!isfinite(value))return INFINITY;p->stretch_energy+=value;
            if(derivatives){for(int k=0;k<4;k++)gradient[k]+=sg[k];for(int k=0;k<16;k++)hessian[k]+=sh[k];}
        }
        if(!derivatives)continue;
        apr_gradient_basis(f->metric+t,g);ArsStencil stencil={0};
        for(int a=0;a<3;a++)if(!ars_stencil_add(p,&stencil,face[a],g[2*a],g[2*a+1]))return INFINITY;
        for(int a=0;a<stencil.n;a++)for(int d=0;d<2;d++){
            int row=stencil.id[a]+d,k=2*d;p->gradient[row]+=stencil.x[a]*gradient[k]+stencil.y[a]*gradient[k+1];
            for(int b=0;b<stencil.n;b++)for(int e=0;e<2;e++){
                int col=stencil.id[b]+e,q=2*e;if(row<col)continue;
                double h=stencil.x[a]*(hessian[4*k+q]*stencil.x[b]+hessian[4*k+q+1]*stencil.y[b])+
                         stencil.y[a]*(hessian[4*(k+1)+q]*stencil.x[b]+hessian[4*(k+1)+q+1]*stencil.y[b]);
                ars_h(p,row,col,h);
            }
        }
    }
    for(size_t i=0;i<f->n_seams;i++)if(ars_incident(s,f->seams+i)){
        const AsmFieldSeam *seam=f->seams+i;
        for(size_t at=0;at<seam->count;at++){
            const AsmWitness *w=f->witness+seam->first+at;
            if(!w->valid || !(w->weight>0))return INFINITY;
            for(int side=0;side<2;side++){
                double residual[2],influence;AsmField_residual(w,side,uv,residual);
                p->energy+=16*w->weight*apr_seam_penalty(residual,&influence);
                if(!derivatives)continue;double weight=16*w->weight*influence;ArsStencil stencil={0};
                for(int a=0;a<4;a++)if(!ars_stencil_add(p,&stencil,w->vertex[side][a],w->coefficient[side][a],0))return INFINITY;
                for(int a=0;a<stencil.n;a++)for(int d=0;d<2;d++){
                    int row=stencil.id[a]+d;double ca=stencil.x[a];
                    p->gradient[row]+=weight*ca*residual[d];
                    for(int b=0;b<stencil.n;b++)ars_h(p,row,stencil.id[b]+d,weight*ca*stencil.x[b]);
                }
            }
        }
    }
    /* Match the reference ordering: regularize the material/seam normal before
     * adding the positive boundary Hessian. */
    if(derivatives)for(int i=0;i<p->n;i++)ars_h(p,i,i,1e-8*fmax(p->diagonal[i],1));
    if(!ars_boundary_energy(p,uv) || (derivatives && p->h.failed))return INFINITY;
    return p->energy+p->stretch_energy;
}
static double ars_violation(const ArsProblem *p,const double *uv)
{
    const AprSolver *s=p->s;const AsmField *f=s->field;double violation=0;
    for(size_t i=0;i<f->n_seams;i++)if(ars_incident(s,f->seams+i)){
        const AsmFieldSeam *seam=f->seams+i;
        for(size_t k=0;k<seam->count;k++)for(int side=0;side<2;side++){
            double r[2];AsmField_residual(f->witness+seam->first+k,side,uv,r);
            for(int q=0;q<16;q++){double angle=q*(6.2831853071795864769/16);violation=fmax(violation,cos(angle)*r[0]+sin(angle)*r[1]-p->bound);}
        }
    }
    return violation;
}
static int ars_constraints(ArsProblem *p)
{
    const AprSolver *s=p->s;const AsmField *f=s->field;double maximum=0;
    for(size_t i=0;i<2*f->nv;i++)maximum=fmax(maximum,fabs(f->uv[i]));
    double margin=fmax(1e-4,8*(nextafter(maximum,INFINITY)-maximum));
    if(!(margin<ASM_SEAM_TOLERANCE_VOX))return 0;p->bound=(ASM_SEAM_TOLERANCE_VOX-margin)*cos(3.14159265358979323846/16);
    p->m=0;p->a.n=0;p->a.failed=0;
    for(size_t i=0;i<f->n_seams;i++)if(ars_incident(s,f->seams+i)){
        const AsmFieldSeam *seam=f->seams+i;
        for(size_t k=0;k<seam->count;k++)for(int side=0;side<2;side++){
            const AsmWitness *w=f->witness+seam->first+k;double r[2];AsmField_residual(w,side,f->uv,r);
            ArsStencil stencil={0};
            for(int a=0;a<4;a++)if(!ars_stencil_add(p,&stencil,w->vertex[side][a],w->coefficient[side][a],0))return 0;
            for(int q=0;q<16;q++){
                double angle=q*(6.2831853071795864769/16),normal[]={cos(angle),sin(angle)},value[2*ARS_STENCIL_MAX];int ids[2*ARS_STENCIL_MAX];
                for(int a=0;a<stencil.n;a++)for(int d=0;d<2;d++){
                    ids[2*a+d]=stencil.id[a]+d;value[2*a+d]=stencil.x[a]*normal[d];
                }
                if(!ars_plane(p,2*stencil.n,ids,value,p->bound-normal[0]*r[0]-normal[1]*r[1]))return 0;
            }
        }
    }
    return 1;
}
static double ars_jacobian_limit(const AprSolver *s)
{
    const AsmField *f=s->field;double limit=INFINITY;
    for(size_t at=0;at<s->n_active;at++){
        size_t t=s->active_faces[at];double j[4],dj[4],lo,hi,det,roots[2];
        AsmField_face(f,t,f->uv,j,&lo,&hi,&det);if(!(det>1e-10))return 0;
        AsmField_face(f,t,s->direction,dj,&lo,&hi,&det);
        double a=dj[0]*dj[3]-dj[1]*dj[2],b=dj[0]*j[3]+j[0]*dj[3]-dj[1]*j[2]-j[1]*dj[2];
        ars_roots(a,b,j[0]*j[3]-j[1]*j[2]-1e-10,roots);
        for(int k=0;k<2;k++)if(roots[k]>0 && isfinite(roots[k]))limit=fmin(limit,roots[k]);
    }
    return limit;
}
typedef struct ArsRank {int provisional;size_t bad_faces,bad_seams;double metric,seam,stretch,energy;} ArsRank;
static ArsRank ars_rank(AprSolver *s,double energy)
{
    const AsmField *f=s->field;ArsRank rank={0};rank.energy=energy;apr_checkpoint(s);
    for(size_t at=0;at<s->n_active;at++){
        size_t t=s->active_faces[at];rank.bad_faces+=!AsmMetric_within25(s->lo[t],s->hi[t]);
        rank.stretch=fmax(rank.stretch,fmax(.75/s->lo[t]-1,s->hi[t]/1.25-1));
    }
    for(size_t c=0;c<f->run->n_charts;c++)if(s->area[c]>0 && (!s->opts.fixed_charts || !s->opts.fixed_charts[c])){
        double failure=fmax(0,.95-s->within10[c]/s->area[c])+fmax(0,.99-s->within25[c]/s->area[c]);
        rank.metric+=failure;rank.provisional|=failure>0;
    }
    for(size_t i=0;i<f->n_seams;i++)if(ars_incident(s,f->seams+i)){
        const AsmSeamMeasure *m=s->seam+i;rank.provisional|=!m->pass;rank.bad_seams+=!m->pass;
        rank.seam+=m->mass>0?fmax(0,m->rms-ASM_SEAM_TOLERANCE_VOX)+fmax(0,.9-m->support/m->mass):INFINITY;
    }
    return rank;
}
static int ars_better(ArsRank a,ArsRank b,int stretch)
{
    if(a.provisional!=b.provisional)return a.provisional<b.provisional;
    if(a.metric!=b.metric)return a.metric<b.metric;
    if(a.seam!=b.seam)return a.seam<b.seam;
    if(stretch && a.stretch!=b.stretch)return a.stretch<b.stretch;
    return a.energy<b.energy;
}
static int ars_archive_safe(ArsProblem *p,const double *initial)
{
    /* Python restoration explores valid, separated fields before ranking its
     * archive. Applying the local/global metric guard to each Newton step
     * changes that solve and can prevent stretch recovery. The normal C
     * pipeline still retains only an archive that preserves its entry gates.
     * Each accepted segment already passed the continuous source/boundary
     * checks; the straight chord from the entry field is not its trajectory. */
    AprSolver *s=p->s;AsmField *f=s->field;
    double *working=f->uv,*trial=s->trial;f->uv=(double *)initial;s->trial=working;
    apr_checkpoint(s);int safe=apr_guard_mode(s,NULL,0,0);
    f->uv=working;s->trial=trial;return safe;
}
#include "asm_repair_initialization.inc"

typedef struct ArsArchive {
    double *initial,*best;ArsRank best_rank,encoded_rank;
    int initialized,has_encoded;
} ArsArchive;
static void ars_remember(ArsProblem *p,ArsArchive *archive,ArsRank rank,ArsRank *level_best,int stretch)
{
    if(ars_better(rank,*level_best,stretch))*level_best=rank;
    /* Preserve Python's area-feasible archive for level/phase decisions. It
     * never bypasses the independent guard on the field C will install. */
    if(!rank.provisional && (!archive->has_encoded || ars_better(rank,archive->encoded_rank,stretch))){
        archive->encoded_rank=rank;archive->has_encoded=1;
    }
    if(ars_better(rank,archive->best_rank,stretch) && ars_archive_safe(p,archive->initial)){
        archive->best_rank=rank;memcpy(archive->best,p->s->field->uv,2*p->s->field->nv*sizeof(double));
    }
}
static int ars_level(AprSolver *s,double grid,int budget,int stretch,ArsArchive *archive,ArsRank *level_best)
{
    AsmField *f=s->field;ArsProblem p={0};p.arena=Arena_new();p.s=s;int rc=0;size_t bytes=2*f->nv*sizeof(double);
    if(!ars_coordinates_new(&p,grid)){rc=-1;goto done;}
    if(!p.n)goto done;
    p.gradient=ARENA_ALLOC(p.arena,(size_t)p.n*sizeof(double));p.diagonal=ARENA_ALLOC(p.arena,(size_t)p.n*sizeof(double));
    p.answer=ARENA_ALLOC(p.arena,(size_t)p.n*sizeof(double));double *target=ARENA_ALLOC(p.arena,(size_t)p.n*sizeof(double));
    if(!ars_boundary_new(&p)){rc=-1;goto done;}
    double energy=ars_energy(&p,f->uv,0);
    if(!isfinite(energy)){fprintf(stderr,"[assemble restoration] deferred: initial boundary is not strictly separated\n");rc=1;goto done;}
    *level_best=ars_rank(s,energy);
    if(!archive->initialized){archive->best_rank=*level_best;archive->initialized=1;}
    ars_remember(&p,archive,*level_best,level_best,stretch);
    if(archive->has_encoded && ars_better(archive->encoded_rank,*level_best,stretch))*level_best=archive->encoded_rank;
    if(stretch)p.stretch=ARENA_CALLOC(p.arena,s->n_active,sizeof(ArsStretch));
    ArsPressure pressure={0};pressure.weight=1;
    s->stats->restoration_levels++;
    if(stretch)s->stats->stretch_levels++;
    if((size_t)p.n>s->stats->restoration_max_coordinates)s->stats->restoration_max_coordinates=(size_t)p.n;
    fprintf(stderr,"[assemble restoration] %d coordinates, grid %.9g, %zu original moving faces, phase %s\n",p.n,grid,s->n_active,stretch?"stretch":"closure");
    /* The binding rows of one closure iteration seed the next: the row layout
     * is fixed within a level (seams, witnesses, sides, sixteen angles) and
     * consecutive solves bind nearly the same set. */
    enum { ARS_WARM_CAP=8192 };
    int *warm=ARENA_ALLOC(p.arena,ARS_WARM_CAP*sizeof(int)),nwarm=0,warm_m=-1;
    for(int iteration=0;iteration<budget && !apr_expired(&s->opts);iteration++){
        double stretch_violation=0;
        if(stretch && !ars_stretch_begin(&p,&pressure,&stretch_violation)){rc=-1;break;}
        double before=ars_energy(&p,f->uv,1);
        double base_before=p.energy;
        if(!isfinite(before) || !ars_constraints(&p)){rc=-1;break;}
        for(int i=0;i<p.n;i++)target[i]=-2*p.gradient[i];
        s->stats->restoration_trials++;
        int solved=0;double qp_t0=ves_clock_sec();
#if ASM_RESTORATION_ACTIVE_SET
        {
            /* The step is the exact optimum of the linearised model on ONE
             * Cholesky factor of H; the cost scales with the binding rows. */
            ActiveSetQpStats qs;ActiveSetQpFactor_T fac=NULL;double model=0;
            int factored=ActiveSetQp_factor(p.n,p.h.n,p.h.row,p.h.col,p.h.value,&fac,&qs);
            int modelled=factored && ActiveSetQp_model(fac,target,&model);
            if(modelled && model< -ASM_RESTORATION_MODEL_LIMIT*fmax(1.0,fabs(base_before))){
                fprintf(stderr,"[assemble restoration] iteration %d: linearised model ill-posed (value %.3g against energy %.6g, ridge %d); level ends\n",iteration+1,model,base_before,qs.ridge_retries);
                ActiveSetQp_free(&fac);break;
            }
            if(factored)ActiveSetQp_limit(fac,ASM_RESTORATION_MODEL_LIMIT*fmax(1.0,fabs(base_before)));
            if(p.m!=warm_m){nwarm=0;warm_m=p.m;}
            solved=factored && ActiveSetQp_solve_warm(fac,target,p.m,p.a.n,p.a.row,p.a.col,p.a.value,p.upper,p.answer,&qs,warm,&nwarm,ARS_WARM_CAP);
            if(!solved)nwarm=0;
            ActiveSetQp_free(&fac);
            if(!solved && qs.reason==ASQP_REASON_ILLPOSED){
                fprintf(stderr,"[assemble restoration] iteration %d: linearised model ill-posed (objective passed %.3g x energy %.6g after %d steps, %.2f s; n %d, rows %d, ridge %d, unconstrained value %.3g); level ends\n",
                    iteration+1,ASM_RESTORATION_MODEL_LIMIT,base_before,qs.iterations,ves_clock_sec()-qp_t0,p.n,p.m,qs.ridge_retries,modelled?model:0.0);
                break;
            }
            double qp_seconds=ves_clock_sec()-qp_t0;
            fprintf(stderr,"[assemble restoration] active-set QP: n %d, nnz %d, rows %d, active %d, warm %d, steps %d, dropped %d, ridge %d, model %.3g, factor %.2f s, solves %.2f s, total %.2f s%s%s\n",
                p.n,p.h.n,p.m,qs.active,qs.warm,qs.iterations,qs.dropped,qs.ridge_retries,modelled?model:0.0,qs.factor_seconds,qs.solve_seconds,qp_seconds,solved?"":" -> REFUSED: ",solved?"":ActiveSetQp_reason(qs.reason));
            if(!solved)solved=ConvexQp_solve(p.n,p.h.n,p.h.row,p.h.col,p.h.value,target,p.m,p.a.n,p.a.row,p.a.col,p.a.value,p.upper,p.answer);
        }
#else
        solved=ConvexQp_solve(p.n,p.h.n,p.h.row,p.h.col,p.h.value,target,p.m,p.a.n,p.a.row,p.a.col,p.a.value,p.upper,p.answer);
        fprintf(stderr,"[assemble restoration] OSQP: n %d, nnz %d, rows %d, %.2f s\n",p.n,p.h.n,p.m,ves_clock_sec()-qp_t0);
#endif
        if(!solved){
            fprintf(stderr,"[assemble restoration] iteration %d: QP refused\n",iteration+1);break;
        }
        ars_coordinate_direction(&p,s->direction);double slope=0;
        for(int i=0;i<p.n;i++)slope+=p.gradient[i]*p.answer[i];
        double old_violation=ars_violation(&p,f->uv),limit=fmin(ars_jacobian_limit(s),ars_boundary_limit(&p,f->uv,s->direction));
        double alpha=fmin(1,.8*limit),violation=old_violation,after=before,base_after=base_before;int accepted=0;s->rigid=0;
        for(int backtrack=0;backtrack<30 && alpha>=1e-12;backtrack++,alpha*=.5){
            memcpy(s->trial,f->uv,bytes);
            for(size_t v=0;v<f->nv;v++)if(p.id[v]>=0)for(int d=0;d<2;d++)s->trial[2*v+d]=f->uv[2*v+d]+alpha*s->direction[2*v+d];
            if(!memcmp(s->trial,f->uv,bytes))break;
            violation=ars_violation(&p,s->trial);
            if(old_violation>1e-6?!(violation<old_violation*(1-1e-4*alpha)):violation>1e-6)continue;
            AsmContactStats contacts={0};
            if(AsmContacts_measure(f->contacts,s->trial,NULL,NULL,NULL,&contacts) || contacts.pairs)continue;
            after=ars_energy(&p,s->trial,0);
            if(!isfinite(after) || (old_violation<=1e-6 && !(slope<0 && after<=before+1e-4*alpha*slope)))continue;
            base_after=p.energy;accepted=1;break;
        }
        fprintf(stderr,"[assemble restoration] iteration %d: violation %.9g -> %.9g, objective %.9g -> %.9g, alpha %.9g, accepted %d\n",
            iteration+1,old_violation,accepted?violation:old_violation,base_before,accepted?base_after:base_before,alpha,accepted);
        if(stretch)fprintf(stderr,"[assemble stretch] pressure %.9g, violation %.9g, merit %.9g -> %.9g\n",pressure.weight,stretch_violation,before,accepted?after:before);
        if(!accepted){
            if(stretch && old_violation<=1e-6 && ars_pressure_increase(&pressure,stretch_violation)){
                s->stats->stretch_pressure_increases++;continue;
            }
            break;
        }
        memcpy(f->uv,s->trial,bytes);s->stats->restoration_steps++;if(stretch)s->stats->stretch_steps++;
        ArsRank rank=ars_rank(s,base_after);ars_remember(&p,archive,rank,level_best,stretch);
        if(alpha<ASM_RESTORATION_STALL_ALPHA){
            fprintf(stderr,"[assemble restoration] iteration %d: alpha %.3g is a stall; level ends\n",iteration+1,alpha);break;
        }
        if(!stretch && old_violation>1e-6 && old_violation-violation<ASM_RESTORATION_STALL_PROGRESS*old_violation && !(base_after<base_before)){
            fprintf(stderr,"[assemble restoration] iteration %d: violation %.9g -> %.9g and energy not reduced is a stall; level ends\n",iteration+1,old_violation,violation);break;
        }
    }
    if(archive->has_encoded)*level_best=archive->encoded_rank;
    fprintf(stderr,"[assemble restoration] level retained: area %s, bad faces %zu, failed seams %zu; continuing from last valid iterate\n",
        level_best->provisional?"provisional":"passed",level_best->bad_faces,level_best->bad_seams);
done:
    if(p.boundary)free(p.boundary->edges);if(p.coordinates)free(p.coordinates->weights);
    free(p.hslots);alq_free(&p.h);alq_free(&p.a);free(p.upper);Arena_dispose(&p.arena);
    return rc;
}
static int apr_restore_source(AprSolver *s)
{
    if(!s->opts.iterations || !s->opts.closure_iterations || !s->n_active || apr_expired(&s->opts))return 0;
    int mode=s->opts.restoration_only;
    if(s->contacts.pairs && mode!=0 && mode!=2 && mode!=6){s->stats->restoration_skipped_contacts++;return 0;}
    AsmField *f=s->field;if(f->nv>(size_t)INT_MAX/2)return -1;
    Arena_T arena=Arena_new();size_t bytes=2*f->nv*sizeof(double);ArsArchive archive={0};
    archive.initial=ARENA_ALLOC(arena,bytes);archive.best=ARENA_ALLOC(arena,bytes);
    memcpy(archive.initial,f->uv,bytes);memcpy(archive.best,f->uv,bytes);
    int rc=0,initialized=0;ArsRank level_best={0};
    if(s->contacts.pairs || mode==6){
        /* Discrete separation may temporarily open a protected seam. Keep
         * the original field and its gates outside that initialization. */
        archive.best_rank=ars_rank(s,INFINITY);archive.initialized=1;
        int separated=apr_initialize_retry(s,archive.initial,mode==6);
        if(separated<=0){s->stats->restoration_skipped_contacts++;rc=separated;goto done;}
        initialized=1;
    }
    if(mode==1 || mode==3)rc=ars_level(s,s->opts.restoration_grid_step,s->opts.closure_iterations,mode==3,&archive,&level_best);
    else{
        const double grids[]={8,4,0};
        for(int level=0;level<3 && !apr_expired(&s->opts);level++){
            rc=ars_level(s,grids[level],s->opts.closure_iterations,0,&archive,&level_best);
            /* Zero-area border contact can also make the boundary barrier
             * singular. The touching query supplies those initial regions. */
            if(rc==1 && !level && !initialized){
                archive.best_rank=ars_rank(s,INFINITY);archive.initialized=1;
                int separated=apr_initialize_retry(s,archive.initial,0);
                if(separated<=0){rc=separated;goto done;}
                initialized=1;rc=ars_level(s,grids[level],s->opts.closure_iterations,0,&archive,&level_best);
            }
            if(rc || (!level_best.provisional && !level_best.bad_faces))break;
        }
        if(!rc && archive.initialized && !level_best.bad_seams && (level_best.provisional || level_best.bad_faces))
            rc=ars_level(s,0,40,1,&archive,&level_best);
    }
    /* The finer solve consumes the last valid iterate, while the retained
     * archive is always checked against the immutable outer entry field. */
done:
    memcpy(f->uv,archive.best,bytes);Arena_dispose(&arena);
    /* The retained archive can be far from the temporary separated field.
     * Restore the spatial partition as well as its exact coordinates. */
    if(AsmContacts_rebuild(f->contacts,f->uv))return -1;
    s->energy=apr_objective(s,f->uv,1,&s->contacts);return isfinite(s->energy) && rc>=0?0:-1;
}
