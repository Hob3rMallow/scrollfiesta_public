/* Compact transactions keep original vertices, source faces and witnesses.
 * Only moving charts are copied back. The incumbent BVH is the authority for
 * every face outside the field; an omitted obstacle is never free space. */
typedef struct AppResult {
    size_t moving_charts,moving_faces,coordinates,field_faces,accepted;
    double prepare_seconds,solve_seconds,check_seconds,seconds;
    AsmContactStats before,after;
    int kept,oversize,expired,expansions;
    const char *reason;
} AppResult;
typedef struct AppProbe {
    const AsmField *parent;
    const size_t *original_face;
    uint8_t *member;
    size_t regressions;
    int compare;
} AppProbe;
static void app_probe_visit(void *context,size_t a,size_t b,double area)
{
    AppProbe *p=context;
    if(p->member)p->member[p->parent->face_chart[b]]=1;
    if(p->compare){
        double old=AsmContacts_pair_faces(p->parent->contacts,p->parent->uv,p->original_face[a],b,NULL,NULL);
        if(!(old>0) || area>old+64*DBL_EPSILON*fmax(1,old))p->regressions++;
    }
}
static void app_faces(const AsmField *f,const uint8_t *core,size_t *faces)
{
    size_t at=0;
    for(size_t c=0;c<f->run->n_charts;c++)if(core[c] && f->face_offset[c]!=SIZE_MAX)
        for(size_t k=0;k<f->run->charts[c].nf;k++){
            faces[at++]=f->face_offset[c]+k;
        }
}
static void app_entry_init(AprSolver *s,AsmField *f,const double *before,const uint8_t *fixed,
                           AsmRepairStats *stats,Arena_T arena)
{
    memset(s,0,sizeof *s);size_t nc=f->run->n_charts;double *candidate=f->uv;
    s->field=f;s->stats=stats;s->fixed=ARENA_CALLOC(arena,f->nv,1);
    s->lo=ARENA_ALLOC(arena,f->nf*sizeof(double));s->hi=ARENA_ALLOC(arena,f->nf*sizeof(double));
    s->within10=ARENA_ALLOC(arena,nc*sizeof(double));s->within25=ARENA_ALLOC(arena,nc*sizeof(double));s->area=ARENA_ALLOC(arena,nc*sizeof(double));
    s->trial10=ARENA_ALLOC(arena,nc*sizeof(double));s->trial25=ARENA_ALLOC(arena,nc*sizeof(double));
    s->seam=ARENA_ALLOC(arena,(f->n_seams?f->n_seams:1)*sizeof(AsmSeamMeasure));
    for(size_t c=0;c<nc;c++)if(fixed[c] && f->vertex_offset[c]!=SIZE_MAX)memset(s->fixed+f->vertex_offset[c],1,f->run->charts[c].nv);
    f->uv=(double *)before;apr_checkpoint(s);f->uv=candidate;
}
typedef struct AppExterior {
    const AsmField *parent;const uint8_t *core,*exclude;size_t *faces;size_t count;
    double *encountered;AprSolver *entry;
} AppExterior;
static void app_exterior_visit(void *context,size_t a,size_t b,double area)
{
    AppExterior *e=context;(void)a;
    if(e->encountered)e->encountered[e->parent->face_chart[b]]+=area;
}
static int app_exterior_guard(const AsmField *f,const double *uv,void *context)
{
    AppExterior *e=context;AsmContactStats contacts={0};
    /* A sequence of rigid steps can accumulate arithmetic drift on inherited
     * thin/bad faces. Use the final transaction's unchanged entry limits now,
     * so a later refused step does not erase all earlier safe progress. */
    if(e->entry){e->entry->trial=(double *)uv;if(!apr_guard_mode(e->entry,uv,1,0))return 0;}
    app_faces(f,e->core,e->faces);
    /* All entry contacts are already in the compact field. Omitted material
     * therefore has a zero-contact budget. Reject unsafe steps early enough
     * for line search to retain a shorter useful move, not just at commit. */
    return !AsmContacts_probe_mesh_prepared(e->parent->contacts,f->contacts,uv,e->count,e->faces,e->exclude,app_exterior_visit,e,&contacts) && !contacts.pairs;
}
/* Check against the immutable entry field, in addition to per-step guards.
 * This catches cumulative metric drift and any temporary restoration state
 * accidentally escaping its archive. */
static int app_guard(AsmField *f,const double *before,const uint8_t *fixed,
                     AsmRepairStats *stats,Arena_T arena)
{
    AprSolver s={0};double *candidate=f->uv;
    AsmRepairStats check={0};check.first_metric_chart=-1;
    app_entry_init(&s,f,before,fixed,&check,arena);s.trial=candidate;
    f->uv=(double *)before;int ok=apr_guard_mode(&s,candidate,1,0);
    if(!ok)fprintf(stderr,"[assemble patch] entry guard rejected: fixed %zu orientation %zu metric %zu chart %zu seam %zu; first metric chart %d face %zu\n",check.rejected_fixed,check.rejected_orientation,check.rejected_metric,check.rejected_chart,check.rejected_seam,check.first_metric_chart,check.first_metric_face);
    for(size_t i=0;i<f->n_seams && ok;i++)if(s.seam[i].pass){AsmSeamMeasure m;AsmField_seam(f,i,candidate,&m);ok=m.pass;if(!ok)fprintf(stderr,"[assemble patch] passing relation %zu regressed\n",f->seams[i].relation);}
    AsmContactStats contacts={0};
    if(ok && (AsmContacts_measure(f->contacts,candidate,NULL,apr_preserve_self_contact,&s,&contacts) || s.self_contact_regressions))ok=0;
    stats->rejected_fixed+=check.rejected_fixed;stats->rejected_orientation+=check.rejected_orientation;
    stats->rejected_metric+=check.rejected_metric;stats->rejected_chart+=check.rejected_chart;stats->rejected_seam+=check.rejected_seam;
    if(stats->first_metric_chart<0 && check.first_metric_chart>=0){stats->first_metric_chart=check.first_metric_chart;stats->first_metric_face=check.first_metric_face;}
    f->uv=candidate;return ok;
}

/* pending[] contains parent seam indices temporarily required by this solve.
 * compact=0 is the diagnostic whole-field arm with the SAME moving mask.
 * No caller must save or copy the full incumbent for the compact arm. */
static int app_attempt_once(AsmField *parent,const uint8_t *core,const size_t *pending,size_t np,
                            const uint8_t *extra,double *encountered,int step_guard,
                            const AsmRepairOpts *opts,AppResult *out,AsmRepairStats *stats)
{
    double start=ves_clock_sec();memset(out,0,sizeof *out);out->reason="deferred";
    Arena_T arena=Arena_new();size_t nc=parent->run->n_charts;int rc=0;
    uint8_t *member=ARENA_ALLOC(arena,nc),*fixed=ARENA_ALLOC(arena,nc);
    memcpy(member,core,nc);
    for(size_t c=0;c<nc;c++){
        if(extra && extra[c])member[c]=1;
        fixed[c]=(uint8_t)!core[c];
        if(core[c] && parent->face_offset[c]!=SIZE_MAX){out->moving_charts++;out->moving_faces+=parent->run->charts[c].nf;out->coordinates+=2*parent->run->charts[c].nv;}
    }
    if(out->moving_faces>opts->patch_faces || out->coordinates>opts->patch_coordinates){
        out->oversize=1;out->reason="moving_resource_cap";goto done;
    }
    if(apr_expired(opts)){out->expired=1;out->reason="deadline";goto done;}
    /* Every complete incident seam is represented. Passing non-required
     * source seams receive the same final protection as required seams. */
    for(size_t i=0;i<parent->n_seams;i++){
        const AsmFieldSeam *s=parent->seams+i;
        if(s->complete && (core[s->a] || core[s->b]))member[s->a]=member[s->b]=1;
    }
    size_t *original=ARENA_ALLOC(arena,(out->moving_faces?out->moving_faces:1)*sizeof(size_t));
    size_t *query=ARENA_ALLOC(arena,(out->moving_faces?out->moving_faces:1)*sizeof(size_t)),at=0;
    for(size_t c=0;c<nc;c++)if(core[c] && parent->face_offset[c]!=SIZE_MAX)
        for(size_t k=0;k<parent->run->charts[c].nf;k++)original[at++]=parent->face_offset[c]+k;
    uint8_t *exclude=ARENA_CALLOC(arena,parent->nf,1);
    for(size_t c=0;c<nc;c++)if(core[c] && parent->face_offset[c]!=SIZE_MAX)memset(exclude+parent->face_offset[c],1,parent->run->charts[c].nf);
    app_faces(parent,core,query);
    AppProbe probe={parent,original,member,0,0};AsmContactStats external_before={0},external_after={0};
    if(AsmContacts_probe_mesh(parent->contacts,parent->uv,parent->contacts,parent->uv,out->moving_faces,query,exclude,app_probe_visit,&probe,&external_before)){
        out->reason="external_walk_budget";goto done;
    }
    for(size_t c=0;c<nc;c++)if(member[c] && parent->face_offset[c]!=SIZE_MAX)out->field_faces+=parent->run->charts[c].nf;
    if(opts->compact && out->field_faces>opts->patch_field_faces){out->oversize=1;out->reason="halo_resource_cap";goto done;}
    if(!opts->compact)memset(member,1,nc);
    AsmField *f=AsmField_new_local(arena,parent->run,member,parent);
    if(!f || AsmField_build_contacts(f)){rc=-1;goto done;}
    out->field_faces=f->nf;AsmContacts_threads(f->contacts,opts->threads);
    uint8_t *active=ARENA_CALLOC(arena,f->nf,1);
    for(size_t c=0;c<nc;c++)if(core[c] && f->face_offset[c]!=SIZE_MAX)memset(active+f->face_offset[c],1,f->run->charts[c].nf);
    AsmContacts_active_faces(f->contacts,active);
    size_t *by_relation=ARENA_ALLOC(arena,parent->run->n_rels*sizeof(size_t));
    for(size_t i=0;i<parent->run->n_rels;i++)by_relation[i]=SIZE_MAX;
    for(size_t i=0;i<f->n_seams;i++)by_relation[f->seams[i].relation]=i;
    for(size_t i=0;i<parent->n_seams;i++)if(by_relation[parent->seams[i].relation]!=SIZE_MAX)
        apa_require(f,by_relation[parent->seams[i].relation],parent->seams[i].required);
    for(size_t k=0;k<np;k++){
        size_t si=by_relation[parent->seams[pending[k]].relation];
        if(si==SIZE_MAX || !f->seams[si].complete){out->reason="missing_source_halo";goto done;}
        apa_require(f,si,1);
    }
    for(size_t i=0;i<f->n_seams;i++)if(f->seams[i].complete && !f->seams[i].required){
        AsmSeamMeasure m;AsmField_seam(f,i,f->uv,&m);if(m.pass)apa_require(f,i,1);
    }
    double *before=ARENA_ALLOC(arena,2*f->nv*sizeof(double));memcpy(before,f->uv,2*f->nv*sizeof(double));
    if(AsmContacts_measure(f->contacts,f->uv,NULL,NULL,NULL,&out->before)){out->reason="local_walk_budget";goto done;}
    memset(exclude,0,parent->nf);
    for(size_t c=0;c<nc;c++)if(member[c] && parent->face_offset[c]!=SIZE_MAX)memset(exclude+parent->face_offset[c],1,parent->run->charts[c].nf);
    /* The complete first probe inserted EVERY initially intersecting chart
     * into member. Thus no entry contact remains outside this field; another
     * identical global probe would only prove the same zero again. */
    memset(&external_before,0,sizeof external_before);external_before.complete=1;
    if(step_guard)AsmContacts_prepare_exclusion(parent->contacts,exclude);
    AsmRepairOpts local=*opts;local.fixed_charts=fixed;local.admission=local.regional=local.source_trials=0;local.checkpoint=NULL;
    AprSolver entry={0};if(!np && step_guard)app_entry_init(&entry,f,before,fixed,stats,arena);
    AppExterior exterior={parent,core,exclude,query,out->moving_faces,encountered,!np && step_guard?&entry:NULL};
    local.external_guard=step_guard?app_exterior_guard:NULL;local.external_context=step_guard?&exterior:NULL;
    double stop=start+opts->patch_seconds;
    if(stop>0 && (local.deadline<=0 || stop<local.deadline))local.deadline=stop;
    out->prepare_seconds=ves_clock_sec()-start;
    fprintf(stderr,"[assemble patch] %zu charts, %zu moving faces, %zu coordinates, %zu field faces, %zu target seams; preparation %.3f s\n",
        out->moving_charts,out->moving_faces,out->coordinates,f->nf,np,out->prepare_seconds);
    double solve_start=ves_clock_sec();int solved=apr_solve(f,&local,stats);
    out->solve_seconds=ves_clock_sec()-solve_start;out->expired=apr_expired(&local);
    if(solved){out->reason=solved==1?"local_walk_budget":"solve_failed";goto done;}
    double check_start=ves_clock_sec();AsmContacts_active_faces(f->contacts,active);
    if(!app_guard(f,before,fixed,stats,arena)){out->reason="entry_guard";goto done;}
    if(AsmContacts_measure(f->contacts,f->uv,NULL,NULL,NULL,&out->after)){out->reason="local_walk_budget";goto done;}
    app_faces(f,core,query);probe.member=NULL;probe.compare=1;probe.regressions=0;
    /* The candidate has lived in its own field; the incumbent UVs and boxes
     * have not changed since the preparation probe. */
    if(AsmContacts_probe_mesh_prepared(parent->contacts,f->contacts,f->uv,out->moving_faces,query,exclude,app_probe_visit,&probe,&external_after)){
        out->reason="external_walk_budget";goto done;
    }
    out->after.pairs+=external_after.pairs;out->after.area+=external_after.area;
    if(probe.regressions){out->reason="external_contact_regression";goto done;}
    if(out->after.pairs>out->before.pairs || out->after.area>out->before.area+1e-9*fmax(1,out->before.area)){
        out->reason="affected_contact_regression";goto done;
    }
    int changed=memcmp(before,f->uv,2*f->nv*sizeof(double))!=0;
    for(size_t k=0;k<np;k++){
        size_t si=by_relation[parent->seams[pending[k]].relation];AsmSeamMeasure m;AsmField_seam(f,si,f->uv,&m);
        if(m.pass)out->accepted++;
    }
    if(!changed && !out->accepted){out->reason="unchanged";goto done;}
    for(size_t c=0;c<nc;c++)if(core[c] && f->vertex_offset[c]!=SIZE_MAX)
        memcpy(parent->uv+2*parent->vertex_offset[c],f->uv+2*f->vertex_offset[c],2*parent->run->charts[c].nv*sizeof(double));
    for(size_t k=0;k<np;k++){
        AsmSeamMeasure m;AsmField_seam(parent,pending[k],parent->uv,&m);
        if(m.pass){apa_require(parent,pending[k],1);parent->run->rels[parent->seams[pending[k]].relation].continuity|=ASM_CONT_REQUIRED;}
    }
    out->check_seconds=ves_clock_sec()-check_start;out->kept=1;out->reason="kept";
done:
    out->seconds=ves_clock_sec()-start;
    AsmContacts_prepare_exclusion(parent->contacts,NULL);
    if(!out->prepare_seconds)out->prepare_seconds=out->seconds;
    if(!out->check_seconds && out->solve_seconds)out->check_seconds=out->seconds-out->prepare_seconds-out->solve_seconds;
    fprintf(stderr,"[assemble patch] %s: %zu targets passed, contacts %zu -> %zu, area %.9g -> %.9g, %.3f s (solve %.3f, check %.3f)%s\n",
        out->reason,out->accepted,out->before.pairs,out->after.pairs,out->before.area,out->after.area,out->seconds,out->solve_seconds,out->check_seconds,out->expired?", deadline reached":"");
    Arena_dispose(&arena);return rc;
}

static int app_attempt(AsmField *parent,const uint8_t *core,const size_t *pending,size_t np,
                       const AsmRepairOpts *opts,AppResult *out,AsmRepairStats *stats)
{
    /* Regional separation can encounter an obstacle absent at entry. Bring
     * that original chart into the next fixed context, within the same caps
     * and deadline. Every subtransaction is guarded before any copyback, so
     * a refused expansion cannot erase the preceding safe improvement. */
    double start=ves_clock_sec();AsmRepairOpts bounded=*opts;
    if(bounded.deadline<=0 || start+opts->patch_seconds<bounded.deadline)bounded.deadline=start+opts->patch_seconds;
    if(np || !opts->compact){
        /* Most source closures never encounter omitted material. Their full
         * commit probe is sufficient; pay for per-step exterior queries only
         * after that probe proves the first proposal needs them. */
        int rc=app_attempt_once(parent,core,pending,np,NULL,NULL,0,&bounded,out,stats);
        if(!rc && opts->compact && !strcmp(out->reason,"external_contact_regression") && bounded.deadline-ves_clock_sec()>3){
            AppResult first=*out;
            fprintf(stderr,"[assemble patch] retrying source closure with exterior line-search feedback\n");
            rc=app_attempt_once(parent,core,pending,np,NULL,NULL,1,&bounded,out,stats);
            out->before=first.before;out->prepare_seconds+=first.prepare_seconds;out->solve_seconds+=first.solve_seconds;out->check_seconds+=first.check_seconds;
            if(first.field_faces>out->field_faces)out->field_faces=first.field_faces;
            out->seconds=ves_clock_sec()-start;
        }
        return rc;
    }
    Arena_T arena=Arena_new();size_t nc=parent->run->n_charts;
    uint8_t *extra=ARENA_CALLOC(arena,nc,1);
    double *encountered=ARENA_CALLOC(arena,nc,sizeof(double));
    ApcOrder *order=ARENA_ALLOC(arena,nc*sizeof(ApcOrder));
    int rc=0,kept=0;AppResult total={0};
    for(int attempt=0;attempt<3;attempt++){
        memset(encountered,0,nc*sizeof(double));AppResult next={0};
        rc=app_attempt_once(parent,core,pending,np,extra,encountered,1,&bounded,&next,stats);
        if(!attempt)total=next;
        else{
            total.expansions++;
            total.prepare_seconds+=next.prepare_seconds;total.solve_seconds+=next.solve_seconds;total.check_seconds+=next.check_seconds;
            if(next.field_faces>total.field_faces)total.field_faces=next.field_faces;
            total.expired|=next.expired;
            if(next.kept){total.after=next.after;total.accepted=next.accepted;}
        }
        kept|=next.kept;
        if(rc || next.oversize || next.expired || bounded.deadline-ves_clock_sec()<3)break;
        size_t count=0,added=0,faces=next.field_faces;
        for(size_t c=0;c<nc;c++)if(encountered[c]>0 && !extra[c])order[count++]=(ApcOrder){c,-encountered[c]};
        qsort(order,count,sizeof *order,apc_order);
        for(size_t i=0;i<count;i++){
            size_t c=order[i].id,nf=parent->run->charts[c].nf;
            if(nf>opts->patch_field_faces-faces)continue;
            extra[c]=1;faces+=nf;added++;
        }
        if(!added)break;
        if(attempt+1<3)fprintf(stderr,"[assemble patch] adding %zu encountered fixed obstacle charts; bounded field estimate %zu faces\n",added,faces);
    }
    total.kept=kept;if(kept){total.reason="kept";total.oversize=0;}
    total.seconds=ves_clock_sec()-start;*out=total;
    Arena_dispose(&arena);return rc;
}

typedef struct AppContactsByChart { const AsmField *field;double *area; } AppContactsByChart;
/* A regional repair needs tree branches as well as cyclic blocks. A causal
 * core from a chart to itself contains only that chart when it is a leaf,
 * pinning its sole source neighbor and preventing component motion. */
static void app_neighborhood(const AsmField *f,size_t seed,uint8_t *core,int hops,uint8_t *scratch)
{
    size_t nc=f->run->n_charts;memset(core,0,nc);core[seed]=1;
    for(int hop=0;hop<hops;hop++){
        memcpy(scratch,core,nc);
        for(size_t i=0;i<f->n_seams;i++){
            const AsmFieldSeam *s=f->seams+i;
            if(s->required && s->complete && (scratch[s->a] || scratch[s->b]))core[s->a]=core[s->b]=1;
        }
    }
}
static void app_contacts_by_chart(void *context,size_t a,size_t b,double area)
{
    AppContactsByChart *p=context;size_t ca=p->field->face_chart[a],cb=p->field->face_chart[b];
    p->area[ca]+=area;if(ca!=cb)p->area[cb]+=area;
}
static uint64_t app_version(const AsmField *f,const uint8_t *core,const size_t *version)
{
    uint64_t hash=UINT64_C(1469598103934665603);
    for(size_t c=0;c<f->run->n_charts;c++)if(core[c]){hash^=(uint64_t)c+UINT64_C(1000003)*version[c];hash*=UINT64_C(1099511628211);}
    /* A fixed source neighbor can unlock this patch without its core moving. */
    for(size_t i=0;i<f->n_seams;i++){const AsmFieldSeam *s=f->seams+i;if(s->complete && (core[s->a] || core[s->b])){
        hash^=(uint64_t)version[s->a]+UINT64_C(1000033)*version[s->b];hash*=UINT64_C(1099511628211);
    }}return hash;
}
/* Only already passing required seams supply connected area. A deferred
 * hypothesis cannot gain priority by pretending to connect material. The
 * queue changes, but each attempted patch keeps the same complete checks. */
static size_t app_source_order(const AsmField *f,Arena_T arena,ApcOrder *order,int connectivity_first)
{
    Arena_Mark mark=Arena_save(arena);size_t nc=f->run->n_charts,count=0;
    size_t *parent=NULL;double *area=NULL;
    if(connectivity_first){
        parent=ARENA_ALLOC(arena,nc*sizeof *parent);area=ARENA_CALLOC(arena,nc,sizeof *area);
        for(size_t c=0;c<nc;c++)parent[c]=c;
        for(size_t i=0;i<f->n_seams;i++){
            const AsmFieldSeam *s=f->seams+i;if(!s->required)continue;
            AsmSeamMeasure m;AsmField_seam(f,i,f->uv,&m);
            if(m.pass)parent[apr_root(parent,s->a)]=apr_root(parent,s->b);
        }
        for(size_t t=0;t<f->nf;t++)area[apr_root(parent,(size_t)f->face_chart[t])]+=f->metric[t].area;
    }
    for(size_t i=0;i<f->n_seams;i++){
        const AsmFieldSeam *s=f->seams+i;if(!s->complete || !s->parity || !apa_source(f->run->rels+s->relation))continue;
        AsmSeamMeasure m;AsmField_seam(f,i,f->uv,&m);
        if(m.pass || m.invalid || m.observations<6 || !(m.mass>=8) || !(m.rms<=ASM_SOURCE_TRIAL_MAX_RMS))continue;
        double score=m.rms+(s->required?64:0);
        if(connectivity_first && s->required){
            size_t a=apr_root(parent,s->a),b=apr_root(parent,s->b);
            if(a!=b)score=-fmin(area[a],area[b]);
        }
        order[count++]=(ApcOrder){i,score};
    }
    qsort(order,count,sizeof *order,apc_order);
    Arena_restore(arena,mark);return count;
}

static int app_schedule(AsmField *f,const AsmRepairOpts *opts,const char *dir,AsmRepairStats *stats)
{
    Arena_T arena=Arena_new();size_t nc=f->run->n_charts;int rc=0;FILE *log=NULL;
    uint8_t *core=ARENA_ALLOC(arena,nc),*scratch=ARENA_ALLOC(arena,nc),*tried=ARENA_CALLOC(arena,f->n_seams,1);
    size_t *version=ARENA_CALLOC(arena,nc,sizeof(size_t)),*pending=ARENA_ALLOC(arena,f->n_seams*sizeof(size_t));
    uint64_t *signature=ARENA_CALLOC(arena,f->n_seams,sizeof(uint64_t));
    ApcOrder *order=ARENA_ALLOC(arena,f->n_seams*sizeof(ApcOrder));
    size_t count=app_source_order(f,arena,order,opts->connectivity_first);
    if(opts->connectivity_first && count)fprintf(stderr,"[assemble patch] connectivity priority: first relation %zu, connected-area score %.9g\n",
                                               f->seams[order[0].id].relation,-order[0].error);
    stats->source_candidates=count;
    if(dir){char path[2048];snprintf(path,sizeof path,"%s/stage6_patches.csv",dir);log=fopen(path,"wb");if(!log){rc=-1;goto done;}
        fprintf(log,"phase,seed,relation,hops,moving_charts,moving_faces,coordinates,field_faces,targets,passed,kept,reason,pairs_before,pairs_after,area_before,area_after,prepare_seconds,solve_seconds,check_seconds,seconds,deadline,obstacle_expansions\n");}
    /* Source patches get the first 75%% of available work time. The remainder
     * clears affected contacts with the same transaction and exterior gates. */
    double source_deadline=opts->deadline>0?ves_clock_sec()+.75*fmax(0,opts->deadline-ves_clock_sec()):0;
    for(int round=0;round<2 && opts->source_trials && !apr_expired(opts);round++){
        size_t kept=0;
        for(size_t k=0;k<count && stats->source_trials<(size_t)opts->source_trials;k++){
            if(apr_expired(opts) || (source_deadline>0 && ves_clock_sec()>=source_deadline))break;
            size_t si=order[k].id;const AsmFieldSeam *seed=f->seams+si;AsmSeamMeasure m;AsmField_seam(f,si,f->uv,&m);if(m.pass)continue;
            int hops=opts->context_hops>0?opts->context_hops:2;
            apr_causal_core(f,seed->a,seed->b,core,hops);
            uint64_t current=app_version(f,core,version);if(tried[si] && current==signature[si])continue;
            size_t np=0;
            for(size_t j=0;j<count;j++){
                size_t id=order[j].id;const AsmFieldSeam *s=f->seams+id;
                if(core[s->a] && core[s->b]){AsmSeamMeasure sm;AsmField_seam(f,id,f->uv,&sm);if(!sm.pass)pending[np++]=id;}
            }
            AsmRepairOpts local=*opts;if(source_deadline>0)local.deadline=source_deadline;
            AppResult result;AsmRepairStats part={0};part.first_metric_chart=-1;
            rc=app_attempt(f,core,pending,np,&local,&result,&part);stats->source_trials++;stats->patches++;
            if(result.oversize && hops>1 && !apr_expired(&local)){
                stats->patches_oversize++;stats->patches_deferred++;
                hops=1;apr_causal_core(f,seed->a,seed->b,core,hops);np=1;pending[0]=si;
                rc=app_attempt(f,core,pending,np,&local,&result,&part);stats->patches++;
            }
            if(rc)goto done;
            apr_accumulate(stats,&part);stats->patches_kept+=result.kept;stats->patches_oversize+=result.oversize;
            stats->patches_deferred+=!result.kept;stats->recovered_connections+=result.accepted;
            if(result.kept){kept++;for(size_t c=0;c<nc;c++)if(core[c])version[c]++;}
            for(size_t j=0;j<np;j++){
                const AsmFieldSeam *s=f->seams+pending[j];apr_causal_core(f,s->a,s->b,scratch,opts->context_hops>0?opts->context_hops:2);
                tried[pending[j]]=1;signature[pending[j]]=app_version(f,scratch,version);
            }
            if(log)fprintf(log,"source,%zu,%zu,%d,%zu,%zu,%zu,%zu,%zu,%zu,%d,%s,%zu,%zu,%.17g,%.17g,%.6f,%.6f,%.6f,%.6f,%d,%d\n",
                k,seed->relation,hops,result.moving_charts,result.moving_faces,result.coordinates,result.field_faces,np,result.accepted,result.kept,result.reason,result.before.pairs,result.after.pairs,result.before.area,result.after.area,result.prepare_seconds,result.solve_seconds,result.check_seconds,result.seconds,result.expired,result.expansions);
            if(log && fflush(log)){rc=-1;goto done;}
        }
        if(!kept)break;
    }
    if(opts->regional && opts->iterations && !apr_expired(opts)){
        double *area=ARENA_CALLOC(arena,nc,sizeof(double));AppContactsByChart by={f,area};AsmContactStats contacts={0};
        AsmContacts_active_faces(f->contacts,NULL);
        if(AsmContacts_measure(f->contacts,f->uv,NULL,app_contacts_by_chart,&by,&contacts)){stats->regions_refused_budget++;goto done;}
        ApcOrder *regions=ARENA_ALLOC(arena,nc*sizeof(ApcOrder));size_t nr=0;
        for(size_t c=0;c<nc;c++)if(area[c]>0)regions[nr++]=(ApcOrder){c,-area[c]};qsort(regions,nr,sizeof *regions,apc_order);
        uint8_t *visited=ARENA_CALLOC(arena,nc,1);
        for(size_t k=0;k<nr && !apr_expired(opts);k++){
            size_t c=regions[k].id;if(visited[c])continue;int hops=2;app_neighborhood(f,c,core,hops,scratch);
            /* Spend this queue's iteration budget on a coarse area solve.
             * Its two coordinates per chart avoid both a vertex solve and
             * thousands of per-contact slack variables. Original geometry,
             * passing seams and omitted obstacles still decide every step. */
            AsmRepairOpts local=*opts;local.clearance_iterations=opts->iterations;local.iterations=0;local.closure_iterations=0;
            local.clearance_method=3; /* two coordinates per chart; no variable per contact */
            AppResult result;AsmRepairStats part={0};part.first_metric_chart=-1;
            rc=app_attempt(f,core,NULL,0,&local,&result,&part);stats->patches++;
            if(result.oversize && !apr_expired(opts)){
                stats->patches_oversize++;stats->patches_deferred++;
                hops=1;app_neighborhood(f,c,core,hops,scratch);rc=app_attempt(f,core,NULL,0,&local,&result,&part);stats->patches++;
            }
            if(rc)goto done;
            apr_accumulate(stats,&part);stats->patches_kept+=result.kept;stats->patches_oversize+=result.oversize;stats->patches_deferred+=!result.kept;
            visited[c]=1;if(result.kept)for(size_t j=0;j<nc;j++)if(core[j])visited[j]=1;
            if(log)fprintf(log,"contacts,%zu,-1,%d,%zu,%zu,%zu,%zu,0,0,%d,%s,%zu,%zu,%.17g,%.17g,%.6f,%.6f,%.6f,%.6f,%d,%d\n",
                c,hops,result.moving_charts,result.moving_faces,result.coordinates,result.field_faces,result.kept,result.reason,result.before.pairs,result.after.pairs,result.before.area,result.after.area,result.prepare_seconds,result.solve_seconds,result.check_seconds,result.seconds,result.expired,result.expansions);
            if(log && fflush(log)){rc=-1;goto done;}
        }
    }
done:
    stats->deadline_reached=apr_expired(opts);
    for(size_t k=0;k<count;k++){
        AsmSeamMeasure m;AsmField_seam(f,order[k].id,f->uv,&m);
        if(!m.pass){stats->source_unresolved++;stats->source_unattempted+=!tried[order[k].id];}
    }
    if(log && fclose(log))rc=-1;Arena_dispose(&arena);return rc;
}
