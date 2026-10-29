/* Factor the complete original-vertex local/global system, as Python does.
 * Both UV axes and the independent ARAP recovery share this factorization. */
typedef struct AprCoo { int *row,*col,n;double *value; } AprCoo;
static void apr_coo(AprCoo *a,int row,int col,double value)
{
    if(!value || row<0 || col<0)return;
    if(row<col){int t=row;row=col;col=t;}
    a->row[a->n]=row;a->col[a->n]=col;a->value[a->n++]=value;
}
static int apr_linear(AprSolver *s)
{
    if(apr_expired(&s->opts))return 0;
    AsmField *f=s->field;s->rigid=0;memset(s->direction,0,2*f->nv*sizeof(double));
    if(!s->factor){
        apr_linear_clear(s);s->linear_arena=Arena_new();
        s->free_id=ARENA_ALLOC(s->linear_arena,f->nv*sizeof(int));
        for(size_t v=0;v<f->nv;v++)s->free_id[v]=s->fixed[v]?-1:s->n_free++;
        if(2*(size_t)s->n_free>s->stats->active_coordinates)s->stats->active_coordinates=2*(size_t)s->n_free;
        if(!s->n_free)return 0;
        size_t faces=s->active_faces?s->n_active:f->nf;
        size_t capacity=6*faces+20*f->n_witness+21*s->local_count;
        if(capacity>INT_MAX)return 0;
        Arena_T arena=Arena_new();AprCoo a={0};
        a.row=ARENA_ALLOC(arena,capacity*sizeof(int));a.col=ARENA_ALLOC(arena,capacity*sizeof(int));a.value=ARENA_ALLOC(arena,capacity*sizeof(double));
        for(size_t at=0;at<faces;at++){
            size_t t=s->active_faces?s->active_faces[at]:at;double g[6];apr_gradient_basis(f->metric+t,g);
            double mass=s->opts.proximal/(s->opts.proximal_length*s->opts.proximal_length)*f->metric[t].area/6;
            for(int k=0;k<3;k++)for(int l=0;l<=k;l++){
                double value=s->face_weight[t]*(g[2*k]*g[2*l]+g[2*k+1]*g[2*l+1])+mass*(k==l?2:1);
                apr_coo(&a,s->free_id[f->faces[3*t+k]],s->free_id[f->faces[3*t+l]],value);
            }
        }
        for(size_t i=0;i<f->n_witness;i++){
            const AsmWitness *w=f->witness+i;if(!w->valid || !w->required)continue;
            for(int side=0;side<2;side++){
                int ids[4],n=0;double values[4]={0};
                for(int k=0;k<4;k++){
                    int id=s->free_id[w->vertex[side][k]];if(id<0)continue;
                    int at=0;while(at<n && ids[at]!=id)at++;if(at==n)ids[n++]=id;
                    values[at]+=w->coefficient[side][k];
                }
                for(int k=0;k<n;k++)for(int l=0;l<=k;l++)apr_coo(&a,ids[k],ids[l],s->seam_scale[2*i+side]*values[k]*values[l]);
            }
        }
        for(size_t i=0;i<s->local_count;i++){
            const AprContactCopy *c=s->local_contacts+i;
            for(int k=0;k<c->count;k++)for(int l=0;l<=k;l++)apr_coo(&a,s->free_id[c->ids[k]],s->free_id[c->ids[l]],2*c->weight*((k==l?1.:0.)-1./c->count));
        }
        double factor_start=ves_clock_sec();
        int rc=Sparse_factor_spd(s->n_free,a.n,a.row,a.col,a.value,&s->factor);
        s->stats->factor_seconds+=ves_clock_sec()-factor_start;Arena_dispose(&arena);
        if(rc){fprintf(stderr,"[assemble repair] original-vertex factorization failed\n");return 0;}
    }
    Arena_Mark mark=Arena_save(s->linear_arena);
    double *rhs=ARENA_ALLOC(s->linear_arena,2*(size_t)s->n_free*sizeof(double)),*answer=ARENA_ALLOC(s->linear_arena,2*(size_t)s->n_free*sizeof(double));
    for(size_t v=0;v<f->nv;v++)if(s->free_id[v]>=0)for(int d=0;d<2;d++)rhs[(size_t)d*s->n_free+s->free_id[v]]=-s->gradient[2*v+d];
    int rc=Sparse_factor_solve_multi(s->factor,rhs,answer,2);
    if(!rc)for(size_t v=0;v<f->nv;v++)if(s->free_id[v]>=0)for(int d=0;d<2;d++)s->direction[2*v+d]=answer[(size_t)d*s->n_free+s->free_id[v]];
    Arena_restore(s->linear_arena,mark);
    if(rc)return 0;
    apr_matrix(s,s->direction,s->ap);double error=0,scale=0;
    for(size_t v=0;v<f->nv;v++)if(!s->fixed[v])for(int d=0;d<2;d++){
        double value=s->ap[2*v+d]+s->gradient[2*v+d];error+=value*value;scale+=s->gradient[2*v+d]*s->gradient[2*v+d];
    }
    if(!isfinite(error) || error>1e-12*fmax(scale,DBL_MIN)){
        fprintf(stderr,"[assemble repair] original-vertex solve residual %.9g\n",sqrt(error/fmax(scale,DBL_MIN)));return 0;
    }
    return apr_dot(s->direction,s->gradient,2*f->nv)<0;
}
