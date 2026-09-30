/* Explicit local-vertex / chart-translation / contact-slack QP. The lift is
 * implicit (at most two coordinates per UV axis); seam RMS rows may involve
 * arbitrarily many local vertices, so their storage is genuinely sparse. */
typedef struct AlqMatrix { int n,capacity,*row,*col;double *value;int failed; } AlqMatrix;
typedef struct AlqProblem {
    Arena_T arena;AprSolver *s;ApcSpace envelope;AlqMatrix h,a;
    int n,m,bcapacity,base_m,base_n,nt,*translation,*vertex;
    double *upper,*target,*delta,*scratch;
    int *touched,*mark;int ntouched,generation;int local;
} AlqProblem;
static void alq_entry(AlqMatrix *a,int r,int c,double value)
{
    if(!value || a->failed)return;
    if(!isfinite(value)){a->failed=1;return;}
    if(a->n==a->capacity){
        if(a->capacity>INT_MAX/2){a->failed=1;return;}
        int cap=a->capacity?2*a->capacity:1024;
        int *rows=malloc((size_t)cap*sizeof(int)),*cols=malloc((size_t)cap*sizeof(int));
        double *values=malloc((size_t)cap*sizeof(double));
        if(!rows || !cols || !values){free(rows);free(cols);free(values);a->failed=1;return;}
        if(a->n){memcpy(rows,a->row,(size_t)a->n*sizeof(int));memcpy(cols,a->col,(size_t)a->n*sizeof(int));memcpy(values,a->value,(size_t)a->n*sizeof(double));}
        free(a->row);free(a->col);free(a->value);a->row=rows;a->col=cols;a->value=values;a->capacity=cap;
    }
    a->row[a->n]=r;a->col[a->n]=c;a->value[a->n++]=value;
}
static void alq_free(AlqMatrix *a){free(a->row);free(a->col);free(a->value);}
static void alq_h(AlqProblem *p,int r,int c,double value)
{ if(r<c){int t=r;r=c;c=t;}alq_entry(&p->h,r,c,2*value); }
static void alq_plane(AlqProblem *p,int n,const int *ids,const double *values,double upper)
{
    if(p->a.failed)return;
    if(!isfinite(upper) || upper < -1e-8){p->a.failed=1;return;}
    if(p->m==p->bcapacity){
        if(p->bcapacity>INT_MAX/2){p->a.failed=1;return;}
        int capacity=p->bcapacity?2*p->bcapacity:1024;
        double *b=realloc(p->upper,(size_t)capacity*sizeof(double));if(!b){p->a.failed=1;return;}
        p->upper=b;p->bcapacity=capacity;
    }
    for(int k=0;k<n;k++)alq_entry(&p->a,p->m,ids[k],values[k]);
    p->upper[p->m++]=fmax(0,upper);
}
static int alq_columns(const AlqProblem *p,size_t v,int d,int ids[2])
{
    int n=0,c=p->s->field->vertex_chart[v];
    if(p->translation[c]>=0)ids[n++]=p->translation[c]+d;
    if(p->vertex[v]>=0)ids[n++]=p->vertex[v]+d;
    return n;
}
static void alq_accumulate(AlqProblem *p,int id,double value)
{
    if(p->mark[id]!=p->generation){p->mark[id]=p->generation;p->touched[p->ntouched++]=id;p->scratch[id]=0;}
    p->scratch[id]+=value;
}
static int alq_witness(AlqProblem *p,const AsmWitness *w,int side,int ids[16],double bx[16],double by[16])
{
    int n=0;
    for(int k=0;k<4;k++)for(int d=0;d<2;d++){
        int columns[2],count=alq_columns(p,(size_t)w->vertex[side][k],d,columns);
        for(int q=0;q<count;q++){
            int at=0;while(at<n && ids[at]!=columns[q])at++;
            if(at==n){ids[n]=columns[q];bx[n]=by[n]=0;n++;}
            (d?by:bx)[at]+=w->coefficient[side][k];
        }
    }
    return n;
}
static void alq_metric_gradient(const double j[4],const double r[4],int lower,double gradient[4])
{
    double v[2],x,y;
    if(lower){apc_symmetric_min(j,r,v);x=-(r[0]*v[0]+r[1]*v[1]);y=-(r[2]*v[0]+r[3]*v[1]);}
    else{
        double a=j[0]*j[0]+j[2]*j[2],b=j[0]*j[1]+j[2]*j[3],c=j[1]*j[1]+j[3]*j[3];
        double lambda=.5*(a+c)+hypot(.5*(a-c),b);
        if(fabs(b)>fabs(lambda-a)){v[0]=b;v[1]=lambda-a;}else{v[0]=lambda-c;v[1]=b;}
        double length=hypot(v[0],v[1]);if(length>0){v[0]/=length;v[1]/=length;}else{v[0]=a>=c?1:0;v[1]=a>=c?0:1;}
        double lo,hi;apc_singular(j,&lo,&hi);x=(j[0]*v[0]+j[1]*v[1])/hi;y=(j[2]*v[0]+j[3]*v[1])/hi;
    }
    gradient[0]=x*v[0];gradient[1]=x*v[1];gradient[2]=y*v[0];gradient[3]=y*v[1];
}
static double alq_constraints(AlqProblem *p,int planes)
{
    AprSolver *s=p->s;AsmField *f=s->field;ApcSpace *e=&p->envelope;double violation=0;
    for(int i=0;i<p->n;i++)if(!isfinite(p->delta[i]))return INFINITY;
    double *residual=ARENA_ALLOC(p->arena,(size_t)(p->base_m?p->base_m:1)*sizeof(double));
    for(int r=0;r<p->base_m;r++)residual[r]=-p->upper[r];
    for(int k=0;k<p->base_n;k++)residual[p->a.row[k]]+=p->a.value[k]*p->delta[p->a.col[k]];
    for(int r=0;r<p->base_m;r++)violation=fmax(violation,residual[r]/2e-7);
    for(size_t si=0;si<f->n_seams;si++)if(e->rms_radius[si]>0){
        const AsmFieldSeam *seam=f->seams+si;double rms2=0;p->ntouched=0;
        if(++p->generation==INT_MAX){memset(p->mark,0,(size_t)p->n*sizeof(int));p->generation=1;}
        for(size_t k=0;k<seam->count;k++){
            size_t at=seam->first+k;const AsmWitness *w=f->witness+at;
            int ids[2][16],n[2];double bx[2][16],by[2][16],value[2][2],length[2];
            for(int side=0;side<2;side++){
                double base[2];AsmField_residual(w,side,f->uv,base);value[side][0]=base[0];value[side][1]=base[1];
                n[side]=alq_witness(p,w,side,ids[side],bx[side],by[side]);
                for(int q=0;q<n[side];q++){value[side][0]+=bx[side][q]*p->delta[ids[side][q]];value[side][1]+=by[side][q]*p->delta[ids[side][q]];}
                length[side]=hypot(value[side][0],value[side][1]);
                if(e->supported[at]){
                    double excess=length[side]-e->support_radius[at];violation=fmax(violation,excess/2e-7);
                    if(planes && excess>2e-7){
                        double ux=value[side][0]/length[side],uy=value[side][1]/length[side],row[16];
                        for(int q=0;q<n[side];q++)row[q]=ux*bx[side][q]+uy*by[side][q];
                        alq_plane(p,n[side],ids[side],row,e->support_radius[at]-ux*base[0]-uy*base[1]);
                    }
                }
            }
            int side=length[1]>=length[0]?1:0;double weight=w->weight/s->seam[si].mass;
            rms2+=weight*length[side]*length[side];
            if(planes)for(int q=0;q<n[side];q++)alq_accumulate(p,ids[side][q],weight*(value[side][0]*bx[side][q]+value[side][1]*by[side][q]));
        }
        double rms=sqrt(rms2),excess=rms-e->rms_radius[si];violation=fmax(violation,excess/2e-7);
        if(planes && excess>2e-7){
            double upper=e->rms_radius[si]-rms;
            double *row=ARENA_ALLOC(p->arena,(size_t)(p->ntouched?p->ntouched:1)*sizeof(double));
            for(int q=0;q<p->ntouched;q++){int id=p->touched[q];row[q]=p->scratch[id]/rms;upper+=row[q]*p->delta[id];}
            alq_plane(p,p->ntouched,p->touched,row,upper);
        }
    }
    if(p->local)for(size_t t=0;t<f->nf;t++){
        int ids[6],n=0;double derivative[6][4]={{0}},g[6];apr_gradient_basis(f->metric+t,g);
        for(int k=0;k<3;k++){
            int col=p->vertex[f->faces[3*t+k]];if(col<0)continue;
            for(int d=0;d<2;d++){ids[n]=col+d;derivative[n][2*d]=g[2*k];derivative[n++][2*d+1]=g[2*k+1];}
        }
        if(!n)continue;
        double j[4],base[4],lo,hi,det;AsmField_face(f,t,f->uv,base,&lo,&hi,&det);memcpy(j,base,sizeof j);
        for(int q=0;q<n;q++)for(int d=0;d<4;d++)j[d]+=p->delta[ids[q]]*derivative[q][d];
        apc_singular(j,&lo,&hi);
        for(int lower=0;lower<2;lower++){
            double excess=lower?e->lower[t]-apc_symmetric_min(j,e->rotation+4*t,NULL):hi-e->upper[t];
            violation=fmax(violation,excess/1e-9);
            if(planes && excess>1e-9){
                double gradient[4],row[6]={0},upper=lower?-e->lower[t]:e->upper[t];
                alq_metric_gradient(j,e->rotation+4*t,lower,gradient);
                for(int d=0;d<4;d++)upper-=gradient[d]*base[d];
                for(int q=0;q<n;q++)for(int d=0;d<4;d++)row[q]+=gradient[d]*derivative[q][d];
                alq_plane(p,n,ids,row,upper);
            }
        }
    }
    return p->a.failed?INFINITY:violation;
}
/* OSQP, deliberately: the clearance QP carries a slack per contact pair and
 * every slack bound binds, so its active set is as large as the problem.
 * A dual active set on it (tried 2026-09-15, factor-once) crawled at a
 * two-thousand-row working set where ADMM takes milliseconds. */
static int alq_solve_bounded(AlqProblem *p);
static int alq_solve(AlqProblem *p)
{
    ConvexQp_time_limit(ASM_CLEARANCE_QP_TIME_LIMIT);
    int ok=alq_solve_bounded(p);
    ConvexQp_time_limit(0);
    return ok;
}
static int alq_solve_bounded(AlqProblem *p)
{
    if(p->h.failed || p->a.failed)return 0;
    double start=ves_clock_sec();
    if(!ConvexQp_solve(p->n,p->h.n,p->h.row,p->h.col,p->h.value,p->target,0,0,NULL,NULL,NULL,NULL,p->delta))return 0;
    for(int cut=0;cut<16;cut++){
        Arena_Mark mark=Arena_save(p->arena);double violation=alq_constraints(p,1);Arena_restore(p->arena,mark);
        if(violation<=1)break;
        if(ves_clock_sec()-start>ASM_CLEARANCE_QP_BUDGET){
            fprintf(stderr,"[assemble clearance] step refused: %d cutting-plane rounds exceeded the %.0f s budget (%d coordinates, %d slacks)\n",cut,(double)ASM_CLEARANCE_QP_BUDGET,p->n,p->m);
            return 0;
        }
        if(p->a.failed || !ConvexQp_solve(p->n,p->h.n,p->h.row,p->h.col,p->h.value,p->target,p->m,p->a.n,p->a.row,p->a.col,p->a.value,p->upper,p->delta))return 0;
    }
    Arena_Mark mark=Arena_save(p->arena);double violation=alq_constraints(p,0);Arena_restore(p->arena,mark);
    if(violation>1){
        double *raw=ARENA_ALLOC(p->arena,(size_t)p->n*sizeof(double));memcpy(raw,p->delta,(size_t)p->n*sizeof(double));
        double lo=0,hi=1;
        for(int k=0;k<48;k++){
            double mid=.5*(lo+hi);for(int i=0;i<p->n;i++)p->delta[i]=mid*raw[i];
            mark=Arena_save(p->arena);violation=alq_constraints(p,0);Arena_restore(p->arena,mark);
            if(violation<=1)lo=mid;else hi=mid;
        }
        for(int i=0;i<p->n;i++)p->delta[i]=lo*raw[i];
    }
    double model=0;for(int i=0;i<p->n;i++)model-=p->target[i]*p->delta[i];
    for(int k=0;k<p->h.n;k++)model+=(p->h.row[k]==p->h.col[k]?.5:1)*p->h.value[k]*p->delta[p->h.row[k]]*p->delta[p->h.col[k]];
    return isfinite(model) && model<0;
}
static int apr_clearance_qp(AprSolver *s,const double *distance,double radius,double penalty)
{
    AlqProblem p={0};p.s=s;p.local=distance!=NULL;p.arena=Arena_new();AsmField *f=s->field;int accepted=0;size_t nc=f->run->n_charts;
    p.translation=ARENA_ALLOC(p.arena,nc*sizeof(int));p.vertex=ARENA_ALLOC(p.arena,f->nv*sizeof(int));
    p.envelope.solver=s;p.envelope.arena=p.arena;p.envelope.count=ARENA_CALLOC(p.arena,nc,sizeof(int));
    for(size_t v=0;v<f->nv;v++)p.vertex[v]=-1;
    for(size_t c=0;c<nc;c++){
        p.translation[c]=-1;if(f->vertex_offset[c]==SIZE_MAX)continue;
        size_t local=0;for(size_t k=0;k<f->run->charts[c].nv;k++){size_t v=f->vertex_offset[c]+k;local+=p.local && !s->fixed[v] && distance[v]<=radius;}
        if(!s->rigid_fixed[c] && (!p.local || local<f->run->charts[c].nv)){p.translation[c]=p.n;p.n+=2;}
        p.envelope.count[c]=(int)(local || p.translation[c]>=0);
    }
    p.nt=p.n;
    if(p.local)for(size_t v=0;v<f->nv;v++)if(!s->fixed[v] && distance[v]<=radius){p.vertex[v]=p.n;p.n+=2;}
    int motion=p.n;if(!motion)goto done;
    size_t workspace=apr_workspace_estimate(f->nv,f->nf,f->n_witness,nc);
    /* Contact-plane storage depends on actual active contacts, not the
     * tessellation size. Include 9 rows with up to 9 coefficients per pair,
     * plus slacks and native sparse triplets before constructing them. */
    apr_bytes_add(&workspace,s->contacts.pairs,9*9*(2*sizeof(int)+sizeof(double))+16*sizeof(double));
    if(s->opts.compact && workspace>s->opts.patch_bytes){
        fprintf(stderr,"[assemble clearance] sparse contact workspace deferred: %zu estimated bytes, %zu-byte budget\n",workspace,s->opts.patch_bytes);
        goto done;
    }
    /* Contact slacks are unknowns too. Refuse their conservative upper bound
     * before allocating contact rows or entering an uninterruptible factor.
     * The caller can still try its bounded material-area direction. */
    if(s->opts.compact && ((size_t)motion>s->opts.patch_coordinates ||
       s->contacts.pairs>s->opts.patch_coordinates-(size_t)motion)){
        fprintf(stderr,"[assemble clearance] deferred: %d motion coordinates plus %zu possible contact slacks exceed the %zu-coordinate patch cap\n",
            motion,s->contacts.pairs,s->opts.patch_coordinates);
        goto done;
    }
    if(!apc_envelope(&p.envelope))goto done;
    AclQuery q={s,ARENA_ALLOC(p.arena,s->contacts.pairs*sizeof(AclPair)),0,p.local};AsmContactStats measured;
    if(AsmContacts_measure(f->contacts,f->uv,NULL,acl_pair,&q,&measured) || !q.count || q.count>(size_t)(INT_MAX-p.n))goto done;
    p.n+=(int)q.count;qsort(q.pairs,q.count,sizeof(AclPair),acl_order);
    for(size_t at=0;at<q.count;){size_t end=at+1;while(end<q.count && q.pairs[end].a==q.pairs[at].a && q.pairs[end].b==q.pairs[at].b)end++;
        for(size_t k=at;k<end;k++)q.pairs[k].weight=penalty/(end-at);at=end;}
    qsort(q.pairs,q.count,sizeof(AclPair),acl_face_order);
    p.target=ARENA_CALLOC(p.arena,p.n,sizeof(double));p.delta=ARENA_CALLOC(p.arena,p.n,sizeof(double));
    p.scratch=ARENA_ALLOC(p.arena,(size_t)p.n*sizeof(double));p.touched=ARENA_ALLOC(p.arena,(size_t)p.n*sizeof(int));p.mark=ARENA_CALLOC(p.arena,p.n,sizeof(int));
    if(p.local){
        for(size_t at=0;at<s->n_active;at++){
            size_t t=s->active_faces[at];double g[6],area=f->metric[t].area;apr_gradient_basis(f->metric+t,g);
            for(int k=0;k<3;k++){
                size_t v=(size_t)f->faces[3*t+k];
                for(int d=0;d<2;d++){
                    int ids[2],n=alq_columns(&p,v,d,ids);
                    for(int a=0;a<n;a++)for(int b=0;b<=a;b++)alq_h(&p,ids[a],ids[b],.01*area/3);
                }
                if(p.vertex[v]<0)continue;
                for(int l=0;l<=k;l++){
                    size_t w=(size_t)f->faces[3*t+l];if(p.vertex[w]<0)continue;
                    double value=area*(g[2*k]*g[2*l]+g[2*k+1]*g[2*l+1]);
                    for(int d=0;d<2;d++)alq_h(&p,p.vertex[v]+d,p.vertex[w]+d,value);
                }
            }
        }
    }else{
        double total=0;size_t count=0;
        for(size_t c=0;c<nc;c++)if(p.translation[c]>=0){total+=s->area[c];count++;}
        if(!(total>0))goto done;
        for(size_t c=0;c<nc;c++)if(p.translation[c]>=0)for(int d=0;d<2;d++)alq_h(&p,p.translation[c]+d,p.translation[c]+d,s->area[c]*count/total);
    }
    for(size_t at=0;at<q.count;at++){
        const AclPair *pair=q.pairs+at;int slack=motion+(int)at;alq_h(&p,slack,slack,pair->weight);p.target[slack]=-2*pair->weight*pair->slack;
        for(int a=0;a<(p.local?3:1);a++)for(int b=0;b<(p.local?3:1);b++){
            int ids[9],n=0;double values[9],bound=0;
            size_t va=(size_t)f->faces[3*pair->face[0]+a],vb=(size_t)f->faces[3*pair->face[1]+b];
            for(int side=0;side<2;side++)for(int d=0;d<2;d++){
                int columns[2],count=alq_columns(&p,side?vb:va,d,columns);
                for(int k=0;k<count;k++){ids[n]=columns[k];values[n++]=(side?-1:1)*pair->normal[d];}
            }
            ids[n]=slack;values[n++]=-1;
            if(p.local){for(int d=0;d<2;d++)bound+=pair->normal[d]*(f->uv[2*vb+d]-f->uv[2*va+d]);bound+=pair->slack-1e-4;}
            /* For translations the extreme-vertex inequality is exact. */
            alq_plane(&p,n,ids,values,fmax(0,bound));
        }
        int id=slack;double value=-1;alq_plane(&p,1,&id,&value,pair->slack);
    }
    p.base_m=p.m;p.base_n=p.a.n;
    fprintf(stderr,"[assemble clearance] %s: %d motion coordinates, %zu contact slacks\n",p.local?"original vertices":"chart translations",motion,q.count);
    if(!alq_solve(&p))goto done;
    memset(s->direction,0,2*f->nv*sizeof(double));s->rigid=0;
    for(size_t v=0;v<f->nv;v++)for(int d=0;d<2;d++){
        int ids[2],n=alq_columns(&p,v,d,ids);for(int k=0;k<n;k++)s->direction[2*v+d]+=p.delta[ids[k]];
    }
    s->stats->constrained_steps++;accepted=apr_search(s);
done:
    alq_free(&p.h);alq_free(&p.a);free(p.upper);Arena_dispose(&p.arena);return accepted;
}
