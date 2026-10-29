/* Explicit source closure, separate from local/global energy descent.
 * Python's restoration retains every directional sample of a failed seam
 * and permits an intermediate radius. Here the existing native pose and
 * metric-recovery span, affine strain and intrinsic shape recovery propose
 * that motion; all original-face and collision
 * checks still decide acceptance. This does not implement the historical
 * full-coordinate SD/coarse-to-fine restoration. */
static int apr_closure_seam(const AprSolver *s,size_t i)
{
    const AsmFieldSeam *seam=s->field->seams+i;const AsmSeamMeasure *m=s->seam+i;
    return seam->required && !m->pass && seam->complete && seam->parity && seam->source_run &&
        !m->invalid && m->mass>0 && (!s->opts.fixed_charts ||
        !s->opts.fixed_charts[seam->a] || !s->opts.fixed_charts[seam->b]);
}
static double apr_closure_violation(const AprSolver *s,const double *uv)
{
    const AsmField *f=s->field;double violation=0;
    for(size_t i=0;i<f->n_seams;i++)if(apr_closure_seam(s,i)){
        const AsmFieldSeam *seam=f->seams+i;
        for(size_t k=0;k<seam->count;k++)for(int d=0;d<2;d++){
            double r[2];AsmField_residual(f->witness+seam->first+k,d,uv,r);
            if(!isfinite(r[0]) || !isfinite(r[1]))return INFINITY;
            violation=fmax(violation,hypot(r[0],r[1])-ASM_SEAM_TOLERANCE_VOX*(1-1e-4));
        }
    }
    return violation;
}
static double apr_closure_failure(const AprSolver *s)
{
    double failure=0;
    for(size_t i=0;i<s->field->n_seams;i++)if(s->field->seams[i].required){
        const AsmSeamMeasure *m=s->seam+i;
        if(m->mass>0 && isfinite(m->rms))failure+=fmax(0,m->rms-ASM_SEAM_TOLERANCE_VOX)+fmax(0,.9-m->support/m->mass);
    }
    return failure;
}
static int apr_closure_search(AprSolver *s,double before)
{
    AsmField *f=s->field;size_t n=2*f->nv;
    for(size_t i=0;i<n;i++)s->p[i]=f->uv[i]+s->direction[i];
    for(int ls=0;ls<12;ls++){
        double alpha=ldexp(1.,-ls);s->stats->trials++;apr_trial(s,alpha);
        if(!memcmp(s->trial,f->uv,n*sizeof(double)))break;
        double violation=apr_closure_violation(s,s->trial);
        if(!(violation<before*(1-1e-4*alpha)))continue;
        if(!apr_guard(s,s->p,alpha))continue;
        AsmContactStats contacts={0};double energy=apr_objective(s,s->trial,0,&contacts);
        if(!isfinite(energy)){s->stats->rejected_energy++;continue;}
        /* Closure cannot buy continuity by adding overlap. In particular,
         * a previously separated field must remain exactly separated. */
        if(s->self_contact_regressions || contacts.area>s->contacts.area+1e-10 ||
           (!s->contacts.pairs && contacts.pairs)){
            s->stats->rejected_contacts++;continue;
        }
        memcpy(f->uv,s->trial,n*sizeof(double));s->energy=energy;s->contacts=contacts;
        return 1;
    }
    return 0;
}
static int apr_closure_better(AprRank a,double failure,AprRank b,double old_failure)
{
    if(a.bad_faces!=b.bad_faces)return a.bad_faces<b.bad_faces;
    if(a.metric_failure!=b.metric_failure)return a.metric_failure<b.metric_failure;
    if(a.bad_seams!=b.bad_seams)return a.bad_seams<b.bad_seams;
    if(failure!=old_failure)return failure<old_failure;
    return apr_better(a,b);
}
static int apr_close_seams(AprSolver *s)
{
    if(!s->opts.iterations || !s->opts.closure_iterations || apr_expired(&s->opts))return 0;
    AsmField *f=s->field;Arena_T arena=Arena_new();size_t bytes=2*f->nv*sizeof(double);int rc=0;
    double *best=ARENA_ALLOC(arena,bytes);memcpy(best,f->uv,bytes);
    AprRank best_rank=apr_rank(s);double best_failure=apr_closure_failure(s);
    for(int iteration=0;iteration<s->opts.closure_iterations && !apr_expired(&s->opts);iteration++){
        apr_checkpoint(s);double before=apr_closure_violation(s,f->uv);if(!(before>0))break;
        int accepted=0;
        const int kinds[]={APC_POSE,APC_AFFINE,APC_DEFORM};
        for(int phase=0;phase<3 && !accepted && !apr_expired(&s->opts);phase++){
            int kind=kinds[phase];
            for(int fraction=0;fraction<4 && !accepted;fraction++){
                s->closure_fraction=ldexp(1.,-fraction);
                if(kind==APC_POSE)apr_rigid(s);else if(kind==APC_DEFORM && !apr_linear(s))break;
                if(apr_constrained(s,kind))accepted=apr_closure_search(s,before);
            }
        }
        s->closure_fraction=0;
        fprintf(stderr,"[assemble closure] iteration %d: violation %.9g -> %.9g, overlap %.9g (%zu pairs), accepted %d\n",
            iteration+1,before,apr_closure_violation(s,f->uv),s->contacts.area,s->contacts.pairs,accepted);
        if(!accepted)fprintf(stderr,"[assemble closure] rejected trials: face metric %zu, chart metric %zu, seams %zu, contacts %zu\n",
            s->stats->rejected_metric,s->stats->rejected_chart,s->stats->rejected_seam,s->stats->rejected_contacts);
        if(!accepted)break;
        s->stats->closure_steps++;
        s->energy=apr_objective(s,f->uv,1,&s->contacts);if(!isfinite(s->energy)){rc=-1;break;}
        AprRank rank=apr_rank(s);double failure=apr_closure_failure(s);
        if(apr_closure_better(rank,failure,best_rank,best_failure)){
            best_rank=rank;best_failure=failure;memcpy(best,f->uv,bytes);
        }
    }
    s->closure_fraction=0;memcpy(f->uv,best,bytes);
    s->energy=apr_objective(s,f->uv,1,&s->contacts);if(!isfinite(s->energy))rc=-1;
    Arena_dispose(&arena);return rc;
}

static int apr_closure_plane_control(void)
{
    /* min x^2 subject to x>=1 has an uphill, nonzero solution. A closure
     * constraint cannot be dropped because the incumbent x=0 is infeasible. */
    Arena_T arena=Arena_new();ApcSpace p={0};int hrow=0,hcol=0,id=0;
    double h=4,target=0,answer=0,row=-1;
    p.arena=arena;p.n=1;p.nt=1;p.rows=&hrow;p.cols=&hcol;p.values=&h;
    p.target=&target;p.delta=&answer;p.closure_fraction=1;
    apc_plane(&p,1,&id,&row,-1);
    int fails=p.nplanes!=1 || !apc_project(&p) || fabs(answer-1)>1e-6;
    row=1;apc_plane(&p,1,&id,&row,0);
    if(apc_project(&p))fails++;
    Arena_dispose(&arena);
    fprintf(stderr,"  native nonzero source-closure constraints: %s (%d failures)\n",fails?"FAIL":"ok",fails);
    return fails;
}

static int apr_closure_intrinsic_control(void)
{
    int fails=0;
    for(int mirror=0;mirror<2;mirror++){
        Arena_T arena=Arena_new();float xyz[]={0,0,0,4,0,0,4,4,0,0,4,0},uv[]={0,0,4,0,4,4,0,4};
        int32_t faces[]={0,1,2,0,2,3};double placed[8],direction[8];uint8_t fixed[4]={0};
        double ct=cos(.43),sn=sin(.43);
        for(int k=0;k<4;k++){
            double x=(mirror?-1:1)*1.04*uv[2*k],y=uv[2*k+1];
            placed[2*k]=10000+ct*x-sn*y;placed[2*k+1]=-3000+sn*x+ct*y;
        }
        AsmChart chart={0};chart.id=0;chart.parent=-1;chart.xyz=xyz;chart.uv=uv;chart.faces=faces;
        chart.nv=4;chart.nf=2;chart.placed=1;chart.placement_state=ASM_PLACE_ROOT;chart.placed_uv=placed;
        if(mirror)chart.flags=ASM_CHART_MIRROR;
        AsmRun run={0};run.arena=arena;run.charts=&chart;run.n_charts=1;
        AsmField *field=AsmField_new(arena,&run);if(!field){fails++;Arena_dispose(&arena);continue;}
        AprSolver s={0};s.field=field;s.fixed=fixed;
        apc_intrinsic_direction(&s,arena,direction);
        double target[8];for(int k=0;k<8;k++)target[k]=field->uv[k]+direction[k];
        for(size_t t=0;t<2;t++){
            double j[4],lo,hi,det;AsmField_face(field,t,target,j,&lo,&hi,&det);
            if(!UvGuard_positive(target,field->faces+3*t) || fabs(lo-1)>1e-10 || fabs(hi-1)>1e-10)fails++;
        }
        if(memcmp(field->uv,placed,sizeof placed) || memcmp(field->xyz,xyz,sizeof xyz))fails++;
        fixed[0]=1;apc_intrinsic_direction(&s,arena,direction);if(direction[0] || direction[1])fails++;
        uv[0]=NAN;apc_intrinsic_direction(&s,arena,direction);
        for(int k=0;k<8;k++)if(direction[k]!=0)fails++;
        Arena_dispose(&arena);
    }
    fprintf(stderr,"  native intrinsic closure parity, fixed vertices and invalid seeds: %s (%d failures)\n",fails?"FAIL":"ok",fails);
    return fails;
}
