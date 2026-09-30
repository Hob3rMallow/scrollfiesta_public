/* Native counterparts of c3_seam_constrained_pose, c3_metric_restoring_step
 * and c3_metric_ray_reallocation. Original rows are never changed by a solve. */
enum { APC_RAY=1, APC_DEFORM=2, APC_POSE=3, APC_AREA=4, APC_AFFINE=5, APC_CONTACT_POSE=6 };
enum { APC_WIDTH=4, APC_PAIR_WIDTH=2*APC_WIDTH };
typedef struct ApcPlane {int count,ids[APC_PAIR_WIDTH];double values[APC_PAIR_WIDTH],rhs;} ApcPlane;
typedef struct ApcOrder { size_t id; double error; } ApcOrder;
typedef struct ApcSpace {
    Arena_T arena;
    AprSolver *solver;
    int kind,width,n,*offset,*count;
    double *basis,*center,*pose_scale,*lower,*upper,*rotation,*support_radius,*rms_radius,*caps;
    uint8_t *supported;
    int *rows,*cols,nt;
    double *values,*target,*delta,*unscaled;
    ApcPlane *planes;
    size_t nplanes,capacity;
    double closure_fraction,closure_violation;
} ApcSpace;
static int apc_energy_orthonormalize(ApcSpace *p);

static int apc_contact_kind(int kind)
{ return kind==APC_AREA || kind==APC_CONTACT_POSE; }

static int apc_rotating_kind(int kind)
{ return kind==APC_POSE || kind==APC_CONTACT_POSE; }

static int apc_order(const void *a,const void *b)
{
    const ApcOrder *x=a,*y=b;
    if(x->error!=y->error)return x->error<y->error?-1:1;
    return (x->id>y->id)-(x->id<y->id);
}
static double *apc_basis(ApcSpace *p,size_t v,int q)
{ return p->basis+2*(v*(size_t)p->width+q); }
static void apc_derivative(ApcSpace *p,size_t t,int q,double j[4])
{
    const AsmField *f=p->solver->field;double g[6];apr_gradient_basis(f->metric+t,g);
    memset(j,0,4*sizeof(double));
    for(int k=0;k<3;k++){
        double *b=apc_basis(p,(size_t)f->faces[3*t+k],q);
        j[0]+=g[2*k]*b[0];j[1]+=g[2*k+1]*b[0];j[2]+=g[2*k]*b[1];j[3]+=g[2*k+1]*b[1];
    }
}
static void apc_singular(const double j[4],double *lo,double *hi)
{
    *hi=.5*(hypot(j[0]+j[3],j[2]-j[1])+hypot(j[0]-j[3],j[2]+j[1]));
    *lo=*hi>0?fabs(j[0]*j[3]-j[1]*j[2])/ *hi:0;
}
static double apc_symmetric_min(const double j[4],const double r[4],double v[2])
{
    double a=r[0]*j[0]+r[2]*j[2],c=r[1]*j[1]+r[3]*j[3];
    double b=.5*(r[0]*j[1]+r[2]*j[3]+r[1]*j[0]+r[3]*j[2]);
    double lambda=.5*(a+c)-hypot(.5*(a-c),b);
    if(v){
        if(fabs(b)>fabs(lambda-a)){v[0]=b;v[1]=lambda-a;}
        else{v[0]=lambda-c;v[1]=b;}
        double length=hypot(v[0],v[1]);
        if(length>0){v[0]/=length;v[1]/=length;}else{v[0]=a<=c?1:0;v[1]=a<=c?0:1;}
    }
    return lambda;
}

static int apc_envelope(ApcSpace *p)
{
    AprSolver *s=p->solver;AsmField *f=s->field;size_t nc=f->run->n_charts;
    p->lower=ARENA_ALLOC(p->arena,f->nf*sizeof(double));p->upper=ARENA_ALLOC(p->arena,f->nf*sizeof(double));
    p->rotation=ARENA_ALLOC(p->arena,4*f->nf*sizeof(double));
    ApcOrder *order=ARENA_ALLOC(p->arena,f->nf*sizeof(ApcOrder));
    for(size_t c=0;c<nc;c++)if(p->count[c]){
        size_t start=f->face_offset[c],nf=f->run->charts[c].nf,n=0;
        double needed=fmin(.95*s->area[c],s->within10[c]),selected=0,threshold=-1;
        for(size_t k=0;k<nf;k++){
            size_t t=start+k;
            if(AsmMetric_within10(s->lo[t],s->hi[t]))order[n++]=(ApcOrder){t,fmax(1-s->lo[t],s->hi[t]-1)};
        }
        qsort(order,n,sizeof *order,apc_order);
        for(size_t k=0;k<n && selected<needed;k++){selected+=f->metric[order[k].id].area;threshold=order[k].error;}
        for(size_t k=0;k<nf;k++){
            size_t t=start+k;double j[4],lo,hi,det,g[6],magnitude[4]={0};
            AsmField_face(f,t,f->uv,j,&lo,&hi,&det);apr_rotation(j,p->rotation+4*t);apr_gradient_basis(f->metric+t,g);
            int tight=AsmMetric_within10(lo,hi) && fmax(1-lo,hi-1)<=threshold;
            double lower=tight?.9:fmin(.75,s->lo[t]),upper=tight?1.1:fmax(1.25,s->hi[t]);
            for(int q=0;q<3;q++){
                size_t v=(size_t)f->faces[3*t+q];
                for(int d=0;d<2;d++)for(int a=0;a<2;a++)magnitude[2*d+a]+=fabs(g[2*q+a]*f->uv[2*v+d]);
            }
            double margin=8*DBL_EPSILON*hypot(hypot(magnitude[0],magnitude[1]),hypot(magnitude[2],magnitude[3]));
            p->lower[t]=fmin(lower+margin,apc_symmetric_min(j,p->rotation+4*t,NULL));
            p->upper[t]=fmax(upper-margin,hi);
        }
    }
    p->supported=ARENA_CALLOC(p->arena,f->n_witness?f->n_witness:1,1);
    p->support_radius=ARENA_CALLOC(p->arena,f->n_witness?f->n_witness:1,sizeof(double));
    p->rms_radius=ARENA_CALLOC(p->arena,f->n_seams?f->n_seams:1,sizeof(double));
    /* The order buffer for faces need not be large enough for seam samples. */
    size_t largest=1;for(size_t i=0;i<f->n_seams;i++)if(f->seams[i].count>largest)largest=f->seams[i].count;
    ApcOrder *samples=ARENA_ALLOC(p->arena,largest*sizeof(ApcOrder));
    if(p->closure_fraction>0)for(size_t i=0;i<f->n_seams;i++){
        const AsmFieldSeam *seam=f->seams+i;const AsmSeamMeasure *m=s->seam+i;
        if(!seam->required || m->pass || !seam->complete || !seam->parity || !seam->source_run ||
           m->invalid || !(m->mass>0) || (!p->count[seam->a] && !p->count[seam->b]))continue;
        for(size_t k=0;k<seam->count;k++)for(int d=0;d<2;d++){
            double r[2];AsmField_residual(f->witness+seam->first+k,d,f->uv,r);
            p->closure_violation=fmax(p->closure_violation,hypot(r[0],r[1])-ASM_SEAM_TOLERANCE_VOX*(1-1e-4));
        }
    }
    for(size_t i=0;i<f->n_seams;i++){
        const AsmFieldSeam *seam=f->seams+i;const AsmSeamMeasure *m=s->seam+i;
        if(!seam->required)continue;
        if(!m->pass){
            if(p->closure_fraction>0 && seam->complete && seam->parity && seam->source_run &&
               !m->invalid && m->mass>0 && (p->count[seam->a] || p->count[seam->b])){
                /* AllBoundarySamples in Python restoration retains both
                 * directions of every failed physical observation. A common
                 * intermediate radius supplies a bounded closure step. */
                double radius=ASM_SEAM_TOLERANCE_VOX*(1-1e-4)+(1-p->closure_fraction)*p->closure_violation;
                for(size_t k=0;k<seam->count;k++){
                    p->supported[seam->first+k]=1;p->support_radius[seam->first+k]=radius;
                }
                p->rms_radius[i]=radius;
            }
            continue;
        }
        for(size_t k=0;k<seam->count;k++){
            size_t at=seam->first+k;double a[2],b[2];
            AsmField_residual(f->witness+at,0,f->uv,a);AsmField_residual(f->witness+at,1,f->uv,b);
            samples[k]=(ApcOrder){at,fmax(hypot(a[0],a[1]),hypot(b[0],b[1]))};
        }
        qsort(samples,seam->count,sizeof *samples,apc_order);double mass=0,threshold=0;
        for(size_t k=0;k<seam->count && mass<.9*m->mass;k++){mass+=f->witness[samples[k].id].weight;threshold=samples[k].error;}
        if(threshold>ASM_SEAM_TOLERANCE_VOX)return 0;
        for(size_t k=0;k<seam->count;k++)if(samples[k].error<=threshold){
            p->supported[samples[k].id]=1;p->support_radius[samples[k].id]=fmax(ASM_SEAM_TOLERANCE_VOX*(1-1e-4),samples[k].error);
        }
        p->rms_radius[i]=fmax(ASM_SEAM_TOLERANCE_VOX*(1-1e-4),m->rms);
    }
    return 1;
}

static void apc_intrinsic_direction(AprSolver *s,Arena_T arena,double *direction)
{
    /* A qualified original intrinsic chart supplies a shape-recovery
     * direction even when the local ARAP gradient is nearly stationary.
     * Fit only its proper rigid pose. Never fit to an accepted sheet. */
    AsmField *f=s->field;size_t bytes=2*f->nv*sizeof(double);
    double *target=ARENA_ALLOC(arena,bytes);memcpy(target,f->uv,bytes);memset(direction,0,bytes);
    for(size_t c=0;c<f->run->n_charts;c++)if(f->vertex_offset[c]!=SIZE_MAX){
        const AsmChart *chart=f->run->charts+c;size_t first=f->vertex_offset[c];
        if(!chart->uv || (s->opts.fixed_charts && s->opts.fixed_charts[c]))continue;
        double source[2]={0},world[2]={0},mass=0,cosine=0,sine=0,sign=(chart->flags&ASM_CHART_MIRROR)?-1:1;
        int finite=1;
        for(size_t k=0;k<chart->nv;k++){
            size_t v=first+k;double x=sign*chart->uv[2*k],y=chart->uv[2*k+1],w=f->mass[v];
            if(!isfinite(x) || !isfinite(y)){finite=0;break;}
            mass+=w;source[0]+=w*x;source[1]+=w*y;world[0]+=w*f->uv[2*v];world[1]+=w*f->uv[2*v+1];
        }
        if(!finite || !(mass>0))continue;
        for(int d=0;d<2;d++){source[d]/=mass;world[d]/=mass;}
        for(size_t k=0;k<chart->nv;k++){
            size_t v=first+k;double x=sign*chart->uv[2*k]-source[0],y=chart->uv[2*k+1]-source[1];
            double u=f->uv[2*v]-world[0],w=f->uv[2*v+1]-world[1];
            cosine+=f->mass[v]*(x*u+y*w);sine+=f->mass[v]*(x*w-y*u);
        }
        double length=hypot(cosine,sine);if(!(length>0) || !isfinite(length))continue;cosine/=length;sine/=length;
        for(size_t k=0;k<chart->nv;k++){
            size_t v=first+k;if(s->fixed[v])continue;
            double x=sign*chart->uv[2*k]-source[0],y=chart->uv[2*k+1]-source[1];
            target[2*v]=world[0]+cosine*x-sine*y;target[2*v+1]=world[1]+sine*x+cosine*y;
        }
        double area=0,good=0;int valid=1;
        for(size_t k=0;k<chart->nf;k++){
            size_t t=f->face_offset[c]+k;double j[4],lo,hi,det;AsmField_face(f,t,target,j,&lo,&hi,&det);
            if(!UvGuard_positive(target,f->faces+3*t) || !AsmMetric_within25(lo,hi)){valid=0;break;}
            area+=f->metric[t].area;if(AsmMetric_within10(lo,hi))good+=f->metric[t].area;
        }
        if(valid && area>0 && good>=.95*area-32*DBL_EPSILON*area)
            for(size_t k=0;k<2*chart->nv;k++)direction[2*first+k]=target[2*first+k]-f->uv[2*first+k];
    }
}

static int apc_space(ApcSpace *p,AprSolver *s,Arena_T arena,int kind,const double *ray)
{
    memset(p,0,sizeof *p);p->arena=arena;p->solver=s;p->kind=kind;
    p->width=kind==APC_AREA?2:kind==APC_AFFINE?4:kind==APC_CONTACT_POSE?3:kind;
    p->closure_fraction=s->closure_fraction;
    AsmField *f=s->field;size_t nc=f->run->n_charts,nv=f->nv;
    p->offset=ARENA_ALLOC(arena,nc*sizeof(int));p->count=ARENA_CALLOC(arena,nc,sizeof(int));
    p->basis=ARENA_CALLOC(arena,2*nv*(size_t)p->width,sizeof(double));p->center=ARENA_CALLOC(arena,2*nc,sizeof(double));
    p->pose_scale=ARENA_ALLOC(arena,nc*sizeof(double));
    double *first=NULL,*gradient=NULL;
    if(kind==APC_DEFORM){
        first=ARENA_ALLOC(arena,2*nv*sizeof(double));memcpy(first,s->direction,2*nv*sizeof(double));
        gradient=ARENA_ALLOC(arena,2*nv*sizeof(double));memcpy(gradient,s->gradient,2*nv*sizeof(double));
        memset(s->gradient,0,2*nv*sizeof(double));
        for(size_t at=0;at<(s->active_faces?s->n_active:f->nf);at++){
            size_t t=s->active_faces?s->active_faces[at]:at;
            double j[4],r[4],g[6],lo,hi,det;AsmField_face(f,t,f->uv,j,&lo,&hi,&det);apr_rotation(j,r);apr_gradient_basis(f->metric+t,g);
            for(int k=0;k<3;k++){size_t v=(size_t)f->faces[3*t+k];if(s->fixed[v])continue;
                for(int d=0;d<2;d++)s->gradient[2*v+d]+=2*f->metric[t].area*((j[2*d]-r[2*d])*g[2*k]+(j[2*d+1]-r[2*d+1])*g[2*k+1]);}
        }
        int recovered=apr_linear(s);memcpy(s->gradient,gradient,2*nv*sizeof(double));
        if(!recovered)memset(s->direction,0,2*nv*sizeof(double));
    }
    if(kind==APC_AFFINE){
        first=ARENA_ALLOC(arena,2*nv*sizeof(double));memcpy(first,s->direction,2*nv*sizeof(double));
        apc_intrinsic_direction(s,arena,s->direction);
    }
    for(size_t c=0;c<nc;c++){
        p->offset[c]=-1;size_t start=f->vertex_offset[c];if(start==SIZE_MAX)continue;
        if((apc_rotating_kind(kind) || kind==APC_AREA) && s->rigid_fixed[c])continue;
        size_t count=f->run->charts[c].nv;double mass=0;
        for(size_t k=0;k<count;k++){size_t v=start+k;mass+=f->mass[v];for(int d=0;d<2;d++)p->center[2*c+d]+=f->mass[v]*f->uv[2*v+d];}
        if(!(mass>0))continue;for(int d=0;d<2;d++)p->center[2*c+d]/=mass;
        if(apc_rotating_kind(kind) || kind==APC_AREA){
            double radius=0;
            for(size_t k=0;k<count;k++){size_t v=start+k;double x=f->uv[2*v]-p->center[2*c],y=f->uv[2*v+1]-p->center[2*c+1];radius+=f->mass[v]*(x*x+y*y);}
            p->pose_scale[c]=fmax(1,sqrt(radius/mass));
            p->count[c]=p->width;
            for(size_t k=0;k<count;k++){size_t v=start+k;apc_basis(p,v,0)[0]=1;apc_basis(p,v,1)[1]=1;
                if(apc_rotating_kind(kind)){apc_basis(p,v,2)[0]=-(f->uv[2*v+1]-p->center[2*c+1])/p->pose_scale[c];apc_basis(p,v,2)[1]=(f->uv[2*v]-p->center[2*c])/p->pose_scale[c];}}
        }else for(int q=0;q<p->width;q++){
            double norm=0;int slot=p->count[c];
            for(size_t k=0;k<count;k++){
                size_t v=start+k;double *b=apc_basis(p,v,slot);
                if(s->fixed[v])continue;
                if(kind==APC_AFFINE && q<3){
                    double x=f->uv[2*v]-p->center[2*c],y=f->uv[2*v+1]-p->center[2*c+1];
                    b[0]=q==0?x:q==2?y:0;b[1]=q==1?y:q==2?x:0;
                }else for(int d=0;d<2;d++)b[d]=kind==APC_RAY?ray[2*v+d]:q? s->direction[2*v+d]:first[2*v+d];
                for(int d=0;d<2;d++)norm+=b[d]*b[d];
            }
            if(!(norm>0)){
                for(size_t k=0;k<count;k++)memset(apc_basis(p,start+k,slot),0,2*sizeof(double));continue;
            }
            p->count[c]++;
        }
        if(p->count[c]){p->offset[c]=p->n;p->n+=p->count[c];}
    }
    if(first)memcpy(s->direction,first,2*nv*sizeof(double));
    if((kind==APC_DEFORM || kind==APC_AFFINE) && !apc_energy_orthonormalize(p))return 0;
    return p->n>0 && apc_envelope(p);
}

static void apc_entry(ApcSpace *p,int row,int col,double value)
{
    if(!value)return;
    if(row<col){int temp=row;row=col;col=temp;}
    p->rows[p->nt]=row;p->cols[p->nt]=col;p->values[p->nt++]=value;
}
static int apc_witness(ApcSpace *p,const AsmWitness *w,int side,int ids[APC_PAIR_WIDTH],double bx[APC_PAIR_WIDTH],double by[APC_PAIR_WIDTH])
{
    const AsmField *f=p->solver->field;int n=0;
    for(int k=0;k<4;k++){
        size_t v=(size_t)w->vertex[side][k],chart=(size_t)f->vertex_chart[v];
        for(int q=0;q<p->count[chart];q++){
            int id=p->offset[chart]+q,at=0;while(at<n && ids[at]!=id)at++;
            if(at==n){ids[n]=id;bx[n]=by[n]=0;n++;}
            double *b=apc_basis(p,v,q);bx[at]+=w->coefficient[side][k]*b[0];by[at]+=w->coefficient[side][k]*b[1];
        }
    }
    return n;
}
static int apc_normal(ApcSpace *p)
{
    AprSolver *s=p->solver;AsmField *f=s->field;size_t nc=f->run->n_charts;
    size_t capacity=10*nc+72*f->n_witness+36*s->local_count+(size_t)p->n;
    if(capacity>INT_MAX || capacity>SIZE_MAX/sizeof(double))return 0;
    p->rows=ARENA_ALLOC(p->arena,capacity*sizeof(int));p->cols=ARENA_ALLOC(p->arena,capacity*sizeof(int));
    p->values=ARENA_ALLOC(p->arena,capacity*sizeof(double));p->target=ARENA_CALLOC(p->arena,p->n,sizeof(double));
    p->delta=ARENA_CALLOC(p->arena,p->n,sizeof(double));
    double *block=ARENA_CALLOC(p->arena,APC_WIDTH*APC_WIDTH*nc,sizeof(double));
    for(size_t c=0;c<nc;c++)if(p->count[c]){
        for(size_t k=0;k<f->run->charts[c].nv;k++){size_t v=f->vertex_offset[c]+k;
            for(int q=0;q<p->count[c];q++){double *b=apc_basis(p,v,q);
                const double *gradient=apc_contact_kind(p->kind)?s->contact_gradient:s->gradient;
                p->target[p->offset[c]+q]-=b[0]*gradient[2*v]+b[1]*gradient[2*v+1];}}
    }
    if(apc_contact_kind(p->kind)){
        double mean=0;size_t charts=0;for(size_t c=0;c<nc;c++)if(p->count[c]){mean+=s->area[c];charts++;}
        mean/=charts;double energy=0,largest=0;
        for(size_t c=0;c<nc;c++)if(p->count[c]){
            double mass=s->area[c]/mean,*force=p->target+p->offset[c],square=0;
            for(int q=0;q<p->count[c];q++)square+=force[q]*force[q];
            largest=fmax(largest,sqrt(square)/mass);energy+=square/mass;
        }
        double scale=fmax(largest/2,energy/fmax(1.05*s->contacts.area,DBL_MIN));
        if(!(scale>0) || !isfinite(scale))return 0;
        for(size_t c=0;c<nc;c++)if(p->count[c])for(int q=0;q<p->count[c];q++)apc_entry(p,p->offset[c]+q,p->offset[c]+q,s->area[c]/mean*scale);
        return 1;
    }
    for(size_t t=0;t<f->nf;t++){
        size_t c=(size_t)f->face_chart[t];int n=p->count[c];if(!n)continue;double j[APC_WIDTH][4];
        for(int q=0;q<n;q++)apc_derivative(p,t,q,j[q]);
        for(int a=0;a<n;a++)for(int b=0;b<=a;b++){
            /* Proper rotations preserve source strain exactly. Python's
             * pose system excludes this otherwise artificial stiffness. */
            double value=0;if(p->kind!=APC_POSE)for(int d=0;d<4;d++)value+=s->face_weight[t]*j[a][d]*j[b][d];
            double mass=s->opts.proximal/(s->opts.proximal_length*s->opts.proximal_length)*f->metric[t].area/6;
            for(int k=0;k<3;k++)for(int l=0;l<3;l++){
                double *ba=apc_basis(p,(size_t)f->faces[3*t+k],a),*bb=apc_basis(p,(size_t)f->faces[3*t+l],b);
                value+=mass*(k==l?2:1)*(ba[0]*bb[0]+ba[1]*bb[1]);
            }
            block[APC_WIDTH*APC_WIDTH*c+APC_WIDTH*a+b]+=value;
        }
    }
    for(size_t c=0;c<nc;c++)for(int a=0;a<p->count[c];a++)for(int b=0;b<=a;b++)apc_entry(p,p->offset[c]+a,p->offset[c]+b,block[APC_WIDTH*APC_WIDTH*c+APC_WIDTH*a+b]);
    for(size_t i=0;i<f->n_witness;i++)if(f->witness[i].valid && f->witness[i].required)for(int side=0;side<2;side++){
        int ids[APC_PAIR_WIDTH];double bx[APC_PAIR_WIDTH],by[APC_PAIR_WIDTH];int n=apc_witness(p,f->witness+i,side,ids,bx,by);
        for(int a=0;a<n;a++)for(int b=0;b<=a;b++)apc_entry(p,ids[a],ids[b],s->seam_scale[2*i+side]*(bx[a]*bx[b]+by[a]*by[b]));
    }
    for(size_t i=0;i<s->local_count;i++){
        const AprContactCopy *contact=s->local_contacts+i;int ids[APC_PAIR_WIDTH],n=0;double bx[6][APC_PAIR_WIDTH]={{0}},by[6][APC_PAIR_WIDTH]={{0}},sx[APC_PAIR_WIDTH]={0},sy[APC_PAIR_WIDTH]={0};
        for(int k=0;k<contact->count;k++){
            size_t v=(size_t)contact->ids[k],c=(size_t)f->vertex_chart[v];
            for(int q=0;q<p->count[c];q++){
                int id=p->offset[c]+q,at=0;while(at<n && ids[at]!=id)at++;if(at==n)ids[n++]=id;
                double *b=apc_basis(p,v,q);bx[k][at]=b[0];by[k][at]=b[1];sx[at]+=b[0];sy[at]+=b[1];
            }
        }
        for(int a=0;a<n;a++)for(int b=0;b<=a;b++){
            /* Absolute copies for proper poses, centered copies for the
             * deformation solve: same merit, different tangent quadratics. */
            double value=p->kind==APC_POSE?0:-(sx[a]*sx[b]+sy[a]*sy[b])/contact->count;
            for(int k=0;k<contact->count;k++)value+=bx[k][a]*bx[k][b]+by[k][a]*by[k][b];
            apc_entry(p,ids[a],ids[b],2*contact->weight*value);
        }
    }
    return 1;
}

static int apc_energy_orthonormalize(ApcSpace *p)
{
    AsmField *f=p->solver->field;size_t nc=f->run->n_charts;
    int width=p->width;size_t stride=(size_t)width*width;
    double *blocks=ARENA_CALLOC(p->arena,stride*nc,sizeof(double));
    Arena_Mark mark=Arena_save(p->arena);
    if(!p->n || !apc_normal(p)){Arena_restore(p->arena,mark);return 0;}
    int *chart=ARENA_ALLOC(p->arena,(size_t)p->n*sizeof(int));
    for(size_t c=0;c<nc;c++)for(int q=0;q<p->count[c];q++)chart[p->offset[c]+q]=(int)c;
    for(int k=0;k<p->nt;k++){
        int c=chart[p->rows[k]];if(c!=chart[p->cols[k]])continue;
        int a=p->rows[k]-p->offset[c],b=p->cols[k]-p->offset[c];
        blocks[stride*c+width*a+b]+=.5*p->values[k];if(a!=b)blocks[stride*c+width*b+a]+=.5*p->values[k];
    }
    Arena_restore(p->arena,mark);p->rows=p->cols=NULL;p->values=p->target=p->delta=NULL;p->nt=0;p->n=0;
    for(size_t c=0;c<nc;c++){
        int old=p->count[c],count=0;double basis[APC_WIDTH][APC_WIDTH]={{0}},*h=blocks+stride*c;
        for(int kind=0;kind<old;kind++){
            double d[APC_WIDTH]={0},original=h[(width+1)*kind];d[kind]=1;if(!(original>0) || !isfinite(original))continue;
            for(int pass=0;pass<2;pass++)for(int q=0;q<count;q++){
                double dot=0;
                for(int a=0;a<width;a++){double hd=0;for(int b=0;b<width;b++)hd+=h[width*a+b]*d[b];dot+=basis[q][a]*hd;}
                for(int j=0;j<width;j++)d[j]-=basis[q][j]*dot;
            }
            double energy=0;
            for(int a=0;a<width;a++){double hd=0;for(int b=0;b<width;b++)hd+=h[width*a+b]*d[b];energy+=d[a]*hd;}
            if(!(energy>1e-12*original))continue;
            for(int j=0;j<width;j++)basis[count][j]=d[j]/sqrt(energy);count++;
        }
        if(old)for(size_t k=0;k<f->run->charts[c].nv;k++){
            size_t v=f->vertex_offset[c]+k;double raw[APC_WIDTH][2]={{0}};
            for(int q=0;q<old;q++)memcpy(raw[q],apc_basis(p,v,q),2*sizeof(double));
            for(int q=0;q<width;q++)for(int d=0;d<2;d++){
                double value=0;for(int component=0;component<width;component++)value+=basis[q][component]*raw[component][d];apc_basis(p,v,q)[d]=value;
            }
        }
        p->count[c]=count;p->offset[c]=count?p->n:-1;p->n+=count;
    }
    return p->n>0;
}

static void apc_plane(ApcSpace *p,int n,const int ids[APC_PAIR_WIDTH],const double gradient[APC_PAIR_WIDTH],double upper)
{
    if(!n || !isfinite(upper))return;
    if(p->nplanes==p->capacity){
        size_t capacity=p->capacity?2*p->capacity:256;
        ApcPlane *rows=ARENA_ALLOC(p->arena,capacity*sizeof *rows);
        if(p->nplanes)memcpy(rows,p->planes,p->nplanes*sizeof *rows);
        p->planes=rows;p->capacity=capacity;
    }
    ApcPlane row={0};double norm=0;
    for(int k=0;k<n;k++)if(gradient[k]){
        row.ids[row.count]=ids[k];row.values[row.count++]=-gradient[k];norm+=gradient[k]*gradient[k];
    }
    if(!row.count || !(norm>0))return;
    /* Preservation cones include zero; closure cones deliberately may not.
     * Only preservation gets the subtraction-roundoff allowance. */
    if(p->closure_fraction>0)row.rhs=-upper;
    else{if(upper < -1e-8)return;row.rhs=-fmax(0,upper);}
    p->planes[p->nplanes++]=row;
}

static double apc_metric(ApcSpace *p,size_t t,const double *delta,double j[4])
{
    AsmField *f=p->solver->field;size_t c=(size_t)f->face_chart[t];double lo,hi,det;
    AsmField_face(f,t,f->uv,j,&lo,&hi,&det);
    for(int q=0;q<p->count[c];q++){
        double d[4];apc_derivative(p,t,q,d);
        for(int k=0;k<4;k++)j[k]+=delta[p->offset[c]+q]*d[k];
    }
    apc_singular(j,&lo,&hi);
    return fmax(hi-p->upper[t],p->lower[t]-apc_symmetric_min(j,p->rotation+4*t,NULL));
}

static void apc_metric_plane(ApcSpace *p,size_t t,int lower)
{
    AsmField *f=p->solver->field;size_t c=(size_t)f->face_chart[t];double j[4],gradient[4],v[2],lo,hi;
    apc_metric(p,t,p->delta,j);apc_singular(j,&lo,&hi);
    if(lower){
        const double *r=p->rotation+4*t;apc_symmetric_min(j,r,v);
        double x=r[0]*v[0]+r[1]*v[1],y=r[2]*v[0]+r[3]*v[1];
        gradient[0]=-x*v[0];gradient[1]=-x*v[1];gradient[2]=-y*v[0];gradient[3]=-y*v[1];
    }else{
        double a=j[0]*j[0]+j[2]*j[2],b=j[0]*j[1]+j[2]*j[3],c2=j[1]*j[1]+j[3]*j[3];
        double value=.5*(a+c2)+hypot(.5*(a-c2),b);
        if(fabs(b)>fabs(value-a)){v[0]=b;v[1]=value-a;}else{v[0]=value-c2;v[1]=b;}
        double length=hypot(v[0],v[1]);if(length>0){v[0]/=length;v[1]/=length;}else{v[0]=a>=c2?1:0;v[1]=a>=c2?0:1;}
        double x=(j[0]*v[0]+j[1]*v[1])/hi,y=(j[2]*v[0]+j[3]*v[1])/hi;
        gradient[0]=x*v[0];gradient[1]=x*v[1];gradient[2]=y*v[0];gradient[3]=y*v[1];
    }
    double base[4],det;AsmField_face(f,t,f->uv,base,&lo,&hi,&det);
    double upper=lower?-p->lower[t]:p->upper[t],row[APC_PAIR_WIDTH]={0};int ids[APC_PAIR_WIDTH];
    for(int k=0;k<4;k++)upper-=gradient[k]*base[k];
    for(int q=0;q<p->count[c];q++){double d[4];apc_derivative(p,t,q,d);ids[q]=p->offset[c]+q;for(int k=0;k<4;k++)row[q]+=gradient[k]*d[k];}
    apc_plane(p,p->count[c],ids,row,upper);
}

static double apc_constraints(ApcSpace *p,int planes)
{
    AprSolver *s=p->solver;AsmField *f=s->field;double ratio=0;
    for(int q=0;q<p->n;q++)if(!isfinite(p->delta[q]))return INFINITY;
    for(size_t i=0;i<f->n_seams;i++)if(p->rms_radius[i]>0){
        const AsmFieldSeam *seam=f->seams+i;double rms2=0,gradient[APC_PAIR_WIDTH]={0};int group_ids[APC_PAIR_WIDTH],ng=0;
        for(size_t k=0;k<seam->count;k++){
            size_t at=seam->first+k;const AsmWitness *w=f->witness+at;
            double residual[2][2],bx[2][APC_PAIR_WIDTH],by[2][APC_PAIR_WIDTH],length[2];int ids[2][APC_PAIR_WIDTH],n[2];
            for(int side=0;side<2;side++){
                double base[2];AsmField_residual(w,side,f->uv,base);residual[side][0]=base[0];residual[side][1]=base[1];
                n[side]=apc_witness(p,w,side,ids[side],bx[side],by[side]);
                for(int q=0;q<n[side];q++){residual[side][0]+=bx[side][q]*p->delta[ids[side][q]];residual[side][1]+=by[side][q]*p->delta[ids[side][q]];}
                length[side]=hypot(residual[side][0],residual[side][1]);
                if(p->supported[at]){
                    double excess=length[side]-p->support_radius[at];ratio=fmax(ratio,excess/2e-7);
                    if(planes && excess>2e-7){
                        double ux=residual[side][0]/length[side],uy=residual[side][1]/length[side],row[APC_PAIR_WIDTH];
                        for(int q=0;q<n[side];q++)row[q]=ux*bx[side][q]+uy*by[side][q];
                        apc_plane(p,n[side],ids[side],row,p->support_radius[at]-ux*base[0]-uy*base[1]);
                    }
                }
            }
            int side=length[1]>=length[0]?1:0;double weight=w->weight/s->seam[i].mass;
            rms2+=weight*length[side]*length[side];
            for(int q=0;q<n[side];q++){
                int pos=0;while(pos<ng && group_ids[pos]!=ids[side][q])pos++;if(pos==ng)group_ids[ng++]=ids[side][q];
                gradient[pos]+=weight*(residual[side][0]*bx[side][q]+residual[side][1]*by[side][q]);
            }
        }
        double rms=sqrt(rms2),excess=rms-p->rms_radius[i];ratio=fmax(ratio,excess/2e-7);
        if(planes && excess>2e-7){
            double upper=p->rms_radius[i]-rms;
            for(int q=0;q<ng;q++){gradient[q]/=rms;upper+=gradient[q]*p->delta[group_ids[q]];}
            apc_plane(p,ng,group_ids,gradient,upper);
        }
    }
    if(p->kind==APC_DEFORM || p->kind==APC_RAY || p->kind==APC_AFFINE)for(size_t c=0;c<f->run->n_charts;c++)if(p->count[c]){
        double upper=0,lower=0;
        for(size_t k=0;k<f->run->charts[c].nf;k++){
            size_t t=f->face_offset[c]+k;double j[4],lo,hi;apc_metric(p,t,p->delta,j);apc_singular(j,&lo,&hi);
            double a=hi-p->upper[t],b=p->lower[t]-apc_symmetric_min(j,p->rotation+4*t,NULL);
            upper=fmax(a,upper);lower=fmax(b,lower);
            if(planes && a>1e-9)apc_metric_plane(p,t,0);
            if(planes && b>1e-9)apc_metric_plane(p,t,1);
        }
        ratio=fmax(ratio,fmax(upper,lower)/1e-9);
        if(p->caps){double x=p->delta[p->offset[c]];ratio=fmax(ratio,fmax(-x,x-p->caps[c])/2e-7);}
    }
    return ratio;
}

static double apc_model(const ApcSpace *p)
{
    double result=0;for(int q=0;q<p->n;q++)result-=p->target[q]*p->delta[q];
    for(int k=0;k<p->nt;k++)result+=(p->rows[k]==p->cols[k]?.5:1)*p->values[k]*p->delta[p->rows[k]]*p->delta[p->cols[k]];
    return result;
}

static void apc_ray_caps(ApcSpace *p)
{
    AsmField *f=p->solver->field;p->caps=ARENA_CALLOC(p->arena,f->run->n_charts,sizeof(double));
    for(size_t c=0;c<f->run->n_charts;c++)if(p->count[c]){
        int id=p->offset[c];double cap=1;
        for(size_t k=0;k<f->run->charts[c].nf;k++){
            size_t t=f->face_offset[c]+k;double j[4];p->delta[id]=cap;
            if(apc_metric(p,t,p->delta,j)<=0)continue;
            double lo=0,hi=cap;
            for(int b=0;b<36;b++){double mid=.5*(lo+hi);p->delta[id]=mid;if(apc_metric(p,t,p->delta,j)<=0)lo=mid;else hi=mid;}
            cap=.995*lo;
        }
        p->delta[id]=0;p->caps[c]=cap;
        double a[APC_PAIR_WIDTH]={1};int ids[APC_PAIR_WIDTH]={id};apc_plane(p,1,ids,a,cap);a[0]=-1;apc_plane(p,1,ids,a,0);
    }
}

static int apc_project(ApcSpace *p)
{
    if(!p->nplanes && !p->closure_fraction)return AsmContinuity_project_step(p->arena,p->n,p->nt,p->rows,p->cols,p->values,p->target,p->delta,NULL,0);
    if(p->nplanes>INT_MAX/APC_PAIR_WIDTH)return 0;
    Arena_T arena=Arena_new();size_t max=APC_PAIR_WIDTH*p->nplanes;
    int *rows=ARENA_ALLOC(arena,max*sizeof(int)),*cols=ARENA_ALLOC(arena,max*sizeof(int)),n=0;
    double *values=ARENA_ALLOC(arena,max*sizeof(double)),*upper=ARENA_ALLOC(arena,p->nplanes*sizeof(double));
    for(size_t r=0;r<p->nplanes;r++){
        const ApcPlane *row=p->planes+r;upper[r]=-row->rhs;
        for(int k=0;k<row->count;k++){rows[n]=(int)r;cols[n]=row->ids[k];values[n++]=-row->values[k];}
    }
    int ok=ConvexQp_solve(p->n,p->nt,p->rows,p->cols,p->values,p->target,(int)p->nplanes,n,rows,cols,values,upper,p->delta);
    Arena_dispose(&arena);return ok;
}

static int apc_solve_space(ApcSpace *p)
{
    if(p->kind==APC_RAY)apc_ray_caps(p);
    if(!apc_project(p))return 0;
    for(int cut=0;cut<16;cut++){
        size_t before=p->nplanes;double violation=apc_constraints(p,1);
        if(violation<=1)break;
        if(before==p->nplanes || !apc_project(p))break;
    }
    /* A closure cone need not contain zero. Scaling toward zero would undo
     * its purpose, and energy descent is not a closure certificate. Retain
     * only a finite proposal; the caller measures actual seam progress and
     * all original geometry before accepting it. */
    if(p->closure_fraction>0)return isfinite(apc_model(p));
    if(apc_constraints(p,0)>1){
        double *raw=ARENA_ALLOC(p->arena,p->n*sizeof(double));memcpy(raw,p->delta,p->n*sizeof(double));
        p->unscaled=raw;
        for(int q=0;q<p->n;q++)if(!isfinite(raw[q]))return 0;
        double lo=0,hi=1;
        for(int k=0;k<48;k++){
            double mid=.5*(lo+hi);for(int q=0;q<p->n;q++)p->delta[q]=mid*raw[q];
            if(apc_constraints(p,0)<=1)lo=mid;else hi=mid;
        }
        for(int q=0;q<p->n;q++)p->delta[q]=lo*raw[q];
    }
    return apc_constraints(p,0)<=1 && isfinite(apc_model(p)) && apc_model(p)<0;
}

static void apc_field(ApcSpace *p,double *direction)
{
    AsmField *f=p->solver->field;memset(direction,0,2*f->nv*sizeof(double));
    for(size_t c=0;c<f->run->n_charts;c++)if(p->count[c])for(size_t k=0;k<f->run->charts[c].nv;k++){
        size_t v=f->vertex_offset[c]+k;
        for(int q=0;q<p->count[c];q++){double *b=apc_basis(p,v,q);for(int d=0;d<2;d++)direction[2*v+d]+=b[d]*p->delta[p->offset[c]+q];}
    }
}

static double apr_direction_model(AprSolver *s,const double *direction,double *work)
{
    apr_matrix(s,direction,work);return .5*apr_dot(direction,work,2*s->field->nv)+apr_dot(direction,s->gradient,2*s->field->nv);
}

static int apr_constrained(AprSolver *s,int kind)
{
    if(apr_expired(&s->opts))return 0;
    Arena_T arena=Arena_new();ApcSpace p;int ok=0;AsmField *f=s->field;
    if(!apc_space(&p,s,arena,kind,NULL) || !apc_normal(&p) || !apc_solve_space(&p))goto done;
    apc_field(&p,s->direction);
    if(kind==APC_DEFORM || kind==APC_AFFINE){
        if(kind==APC_DEFORM && !s->closure_fraction){
        /* As in metric_ray_reallocation: keep the common ray unless a
         * coupled chart allocation gives a lower value of the SAME model. */
        double *common=ARENA_ALLOC(arena,2*f->nv*sizeof(double));memcpy(common,s->direction,2*f->nv*sizeof(double));
        double *raw=common;
        if(p.unscaled){
            raw=ARENA_ALLOC(arena,2*f->nv*sizeof(double));double *saved=p.delta;p.delta=p.unscaled;apc_field(&p,raw);p.delta=saved;
        }
        double *work=ARENA_ALLOC(arena,2*f->nv*sizeof(double));double before=apr_direction_model(s,common,work);
        Arena_T ray_arena=Arena_new();ApcSpace ray;
        if(apc_space(&ray,s,ray_arena,APC_RAY,raw) && apc_normal(&ray) && apc_solve_space(&ray)){
            apc_field(&ray,s->direction);double after=apr_direction_model(s,s->direction,work);
            if(!isfinite(after) || after>=before)memcpy(s->direction,common,2*f->nv*sizeof(double));
        }
        Arena_dispose(&ray_arena);
        }
        s->rigid=0;
    }else{
        double cap=1;
        for(size_t c=0;c<f->run->n_charts;c++)if(p.count[c]){
            const double *x=p.delta+p.offset[c];double length=hypot(x[0],x[1]);
            if(kind==APC_CONTACT_POSE)length=hypot(length,x[2]);
            if(apc_contact_kind(kind) && length>2)cap=fmin(cap,2/length);
            if(apc_rotating_kind(kind)){double angle=fabs(x[2])/p.pose_scale[c];if(angle>1.5707963267948966)cap=fmin(cap,1.5707963267948966/angle);}
        }
        memset(s->rigid_motion,0,5*f->run->n_charts*sizeof(double));s->rigid=1;
        for(size_t c=0;c<f->run->n_charts;c++)if(p.count[c]){
            const double *x=p.delta+p.offset[c];double *m=s->rigid_motion+5*c;
            m[0]=p.center[2*c];m[1]=p.center[2*c+1];m[2]=cap*x[0];m[3]=cap*x[1];m[4]=apc_rotating_kind(kind)?cap*x[2]/p.pose_scale[c]:0;
        }
        apr_trial(s,1);for(size_t v=0;v<2*f->nv;v++)s->direction[v]=s->trial[v]-f->uv[v];
    }
    s->stats->constrained_steps++;ok=1;
done:
    Arena_dispose(&arena);return ok;
}
