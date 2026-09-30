/* Native material admission: the entire missing SOURCE component, all of its
 * incident hypotheses, and an immutable copy of every emitted original face.
 * An unsuccessful trial changes neither membership nor the incumbent field. */
static int apa_source(const AsmRelation *r)
{return !(r->flags&ASM_REL_CONTACT) && (r->continuity&(ASM_CONT_SOURCE|ASM_CONT_CUT));}

static int apa_required_source(const AsmRelation *r)
{return apa_source(r) && ((r->continuity&(ASM_CONT_REQUIRED|ASM_CONT_CUT)) || AsmRel_is_join(r));}

/* An ordered source run is a hypothesis, not by itself a selected join.
 * Retain every original obligation. If none anchors this single-chart
 * diagnostic, propose the strongest eligible run and certify it below.
 * Cross-wrap and cleaning vetoes cannot create a new anchor here. */
static int apa_select_frontier_anchor(AsmRun *trial,const uint8_t *core,size_t requested)
{
    size_t count=0,best=SIZE_MAX;int required=0;
    for(size_t c=0;c<trial->n_charts;c++)count+=core[c]!=0;
    if(count!=1)return 0;
    for(size_t i=0;i<trial->n_rels;i++){
        AsmRelation *r=trial->rels+i;
        if(!apa_source(r) || core[r->a]==core[r->b] ||
           !AsmChart_registered(trial->charts+r->a) || !AsmChart_registered(trial->charts+r->b))continue;
        if(apa_required_source(r)){required=1;continue;}
        if(r->flags&(ASM_REL_DROPPED|ASM_REL_CONTACT|ASM_REL_CROSSWRAP))continue;
        if(best==SIZE_MAX || r->seam_len>trial->rels[best].seam_len ||
           (r->seam_len==trial->rels[best].seam_len && r->corr_count>trial->rels[best].corr_count))best=i;
    }
    if(required)return 1;
    if(requested!=SIZE_MAX){
        if(requested>=trial->n_rels)return 0;
        AsmRelation *r=trial->rels+requested;
        if(!apa_source(r) || core[r->a]==core[r->b] ||
           !AsmChart_registered(trial->charts+r->a) || !AsmChart_registered(trial->charts+r->b) ||
           (r->flags&(ASM_REL_DROPPED|ASM_REL_CONTACT|ASM_REL_CROSSWRAP)))return 0;
        best=requested;
    }
    if(best==SIZE_MAX)return 0;
    trial->rels[best].continuity|=ASM_CONT_REQUIRED;
    return 1;
}

static int apa_chart_metric(const AsmField *f,size_t c)
{
    size_t first=f->face_offset[c];if(first==SIZE_MAX)return 0;
    double area=0,good=0;
    for(size_t k=0;k<f->run->charts[c].nf;k++){
        size_t t=first+k;double j[4],lo,hi,det;AsmField_face(f,t,f->uv,j,&lo,&hi,&det);
        if(!UvGuard_positive(f->uv,f->faces+3*t) || !AsmMetric_within25(lo,hi))return 0;
        area+=f->metric[t].area;if(AsmMetric_within10(lo,hi))good+=f->metric[t].area;
    }
    return area>0 && good>=.95*area-32*DBL_EPSILON*area;
}

/* Rigid packing is exact in real arithmetic; stored coordinates round.
 * Opt-in refinements reserve room inside the unchanged final face bounds. */
static int apa_chart_metric_margin(const AsmField *f,size_t c,double margin)
{
    if(!isfinite(margin) || margin<0 || margin>=.25)return 0;
    size_t first=f->face_offset[c];if(first==SIZE_MAX)return 0;
    for(size_t k=0;k<f->run->charts[c].nf;k++){
        double j[4],lo,hi,det;AsmField_face(f,first+k,f->uv,j,&lo,&hi,&det);
        /* An interior reserve must remain inside the nominal band; adding the
         * comparison allowance here would cancel a reserve of the same size. */
        if(margin>0){if(!isfinite(lo) || !isfinite(hi) || lo<.75+margin || hi>1.25-margin)return 0;}
        else if(!AsmMetric_within25(lo,hi))return 0;
    }
    return 1;
}

/* Every directional row is an affine source point on each endpoint. Fitting
 * both directions transports the measured physical gap instead of welding. */
static int apa_fit(AsmField *f,size_t c,const uint8_t *known,size_t selected)
{
    double center[2]={0},target[2]={0},mass=0,cosine=0,sine=0;
    const AsmFieldSeam *pair=f->seams+selected;
    for(int pass=0;pass<2;pass++)for(size_t si=0;si<f->n_seams;si++){
        const AsmFieldSeam *seam=f->seams+si;
        if(!((seam->a==pair->a && seam->b==pair->b) || (seam->a==pair->b && seam->b==pair->a)))continue;
        if(!seam->required || !seam->complete || (seam->a!=(int32_t)c && seam->b!=(int32_t)c))continue;
        int other=seam->a==(int32_t)c?seam->b:seam->a;if(!known[other])continue;
        for(size_t k=0;k<seam->count;k++)for(int d=0;d<2;d++){
            const AsmWitness *w=f->witness+seam->first+k;if(!w->valid)return 0;
            double u[2]={0},v[2]={0},sum=0;
            for(int q=0;q<4;q++){
                size_t id=(size_t)w->vertex[d][q];double coefficient=w->coefficient[d][q];
                int own=f->vertex_chart[id]==(int32_t)c;if(own)sum+=coefficient;
                for(int a=0;a<2;a++)if(own)u[a]+=coefficient*f->uv[2*id+a];else v[a]-=coefficient*f->uv[2*id+a];
            }
            if(fabs(fabs(sum)-1)>1e-8)return 0;
            for(int a=0;a<2;a++){u[a]/=sum;v[a]/=sum;}
            double weight=.5*w->weight;
            if(!pass){mass+=weight;for(int a=0;a<2;a++){center[a]+=weight*u[a];target[a]+=weight*v[a];}}
            else{
                for(int a=0;a<2;a++){u[a]-=center[a]/mass;v[a]-=target[a]/mass;}
                cosine+=weight*(u[0]*v[0]+u[1]*v[1]);sine+=weight*(u[0]*v[1]-u[1]*v[0]);
            }
        }
    }
    if(!(mass>0) || !isfinite(mass) || !isfinite(cosine) || !isfinite(sine))return 0;
    for(int a=0;a<2;a++){center[a]/=mass;target[a]/=mass;}
    double angle=atan2(sine,cosine),ct=cos(angle),sn=sin(angle);size_t first=f->vertex_offset[c];
    for(size_t k=0;k<f->run->charts[c].nv;k++){
        size_t v=first+k;double x=f->uv[2*v]-center[0],y=f->uv[2*v+1]-center[1];
        f->uv[2*v]=ct*x-sn*y+target[0];f->uv[2*v+1]=sn*x+ct*y+target[1];
    }
    return 1;
}

typedef struct ApaEdge {size_t seam;int32_t a,b;double mass;} ApaEdge;
typedef struct ApaNeighbor {size_t seam;int32_t chart;} ApaNeighbor;
static int apa_edge_order(const void *a,const void *b)
{
    const ApaEdge *x=a,*y=b;
    if(x->mass!=y->mass)return x->mass>y->mass?-1:1;
    if(x->a!=y->a)return x->a<y->a?-1:1;
    if(x->b!=y->b)return x->b<y->b?-1:1;
    return x->seam<y->seam?-1:x->seam>y->seam;
}
static int apa_neighbor_order(const void *a,const void *b)
{
    const ApaNeighbor *x=a,*y=b;return x->chart<y->chart?-1:x->chart>y->chart;
}
static int apa_pair_order(const void *a,const void *b)
{
    const ApaEdge *x=a,*y=b;
    if(x->a!=y->a)return x->a<y->a?-1:1;
    if(x->b!=y->b)return x->b<y->b?-1:1;
    return x->seam<y->seam?-1:x->seam>y->seam;
}
static int apa_initialize(AsmField *f,const uint8_t *core,uint8_t *known)
{
    /* Python rigid_initialize uses a maximum-support spanning forest of
     * source hypotheses, for missing material and existing groups alike.
     * Fitting all current anchors at once is
     * a different seed and can put a whole chart across a neighbouring edge.
     * The omitted forest cycles remain required in the solve and final audit. */
    Arena_T arena=Arena_new();size_t nc=f->run->n_charts,count=0,edges=0;int ok=0;
    ApaEdge *order=ARENA_ALLOC(arena,f->n_seams*sizeof(ApaEdge));
    size_t *parent=ARENA_ALLOC(arena,nc*sizeof(size_t)),*offset=ARENA_CALLOC(arena,nc+1,sizeof(size_t));
    for(size_t c=0;c<nc;c++)parent[c]=c;
    for(size_t i=0;i<f->n_seams;i++){
        const AsmFieldSeam *seam=f->seams+i;
        if(!seam->required || !seam->complete || (!core[seam->a] && !core[seam->b]))continue;
        double mass=0;for(size_t k=0;k<seam->count;k++)mass+=f->witness[seam->first+k].weight;
        if(!(mass>0) || !isfinite(mass))goto done;
        int32_t a=seam->a<seam->b?seam->a:seam->b,b=seam->a<seam->b?seam->b:seam->a;
        order[count++]=(ApaEdge){i,a,b,mass};
    }
    /* Python groups all physical rows for a chart pair before choosing the
     * forest. Splitting that evidence into relation records must not change
     * the selected path or omit a row from its proper rigid fit. */
    qsort(order,count,sizeof *order,apa_pair_order);size_t pairs=0;
    for(size_t i=0;i<count;i++){
        if(pairs && order[i].a==order[pairs-1].a && order[i].b==order[pairs-1].b)
            order[pairs-1].mass+=order[i].mass;
        else order[pairs++]=order[i];
    }
    count=pairs;qsort(order,count,sizeof *order,apa_edge_order);
    for(size_t i=0;i<count;i++){
        size_t a=apr_root(parent,order[i].a),b=apr_root(parent,order[i].b);if(a==b)continue;
        parent[a>b?a:b]=a<b?a:b;order[edges++]=order[i];
        offset[order[i].a+1]++;offset[order[i].b+1]++;
    }
    for(size_t c=1;c<=nc;c++)offset[c]+=offset[c-1];
    ApaNeighbor *neighbor=ARENA_ALLOC(arena,2*edges*sizeof(ApaNeighbor));
    size_t *fill=ARENA_ALLOC(arena,nc*sizeof(size_t));memcpy(fill,offset,nc*sizeof(size_t));
    for(size_t i=0;i<edges;i++){
        neighbor[fill[order[i].a]++]=(ApaNeighbor){order[i].seam,order[i].b};
        neighbor[fill[order[i].b]++]=(ApaNeighbor){order[i].seam,order[i].a};
    }
    size_t *queue=ARENA_ALLOC(arena,nc*sizeof(size_t)),head=0,tail=0;
    for(size_t c=0;c<nc;c++){
        qsort(neighbor+offset[c],offset[c+1]-offset[c],sizeof *neighbor,apa_neighbor_order);
        if(known[c])queue[tail++]=c;
    }
    while(head<tail){size_t a=queue[head++];
        for(size_t i=offset[a];i<offset[a+1];i++){
            size_t b=(size_t)neighbor[i].chart;if(known[b])continue;
            if(!core[b] || !apa_fit(f,b,known,neighbor[i].seam))goto done;
            known[b]=1;queue[tail++]=b;
        }
    }
    for(size_t c=0;c<nc;c++)if(core[c] && !known[c])goto done;
    ok=1;
done:
    Arena_dispose(&arena);return ok;
}

/* Take one more seam into the source-repair batch, or refuse it because the batch's MOVING
 * CHARTS -- the union of its members' causal cores -- would pass
 * ASM_SOURCE_TRIALS_BATCH_MAX_CHARTS.  The first member is always taken however large its own
 * core is, so no seam becomes unreachable.  On success the core mask and *moving are updated. */
static int apa_batch_take(uint8_t *core, const uint8_t *core_k, size_t nc, const size_t *vertex_offset,
                          size_t *moving, size_t members)
{
    size_t add=0,c=0;
    for(c=0;c<nc;c++)add+=core_k[c] && !core[c] && vertex_offset[c]!=SIZE_MAX;
    if(members>0 && ASM_SOURCE_TRIALS_BATCH_MAX_CHARTS>0 &&
       *moving+add>(size_t)ASM_SOURCE_TRIALS_BATCH_MAX_CHARTS)return 0;
    for(c=0;c<nc;c++)core[c]|=core_k[c];
    *moving+=add;
    return 1;
}

static int apa_intrinsic(AsmRun *trial,const uint8_t *core,size_t *oversize,const AsmRepairOpts *opts,
                         const char **reason,FILE *trace)
{
    /* Preserve an already qualified intrinsic seed. Cold recovery uses only
     * the candidate's original XYZ/F, never an accepted/repaired target UV. */
    if(oversize)*oversize=0;
    for(size_t c=0;c<trial->n_charts;c++)if(core[c]){
        if(apr_expired(opts)){*reason="intrinsic_recovery_deadline";return 0;}
        const AsmChart *source=trial->charts+c;
        if(!apr_workspace_check(opts,apr_workspace_estimate(source->nv,source->nf,0,1))){
            *reason="workspace_memory_budget";return 0;
        }
        Arena_T scratch=Arena_new();AsmChart local=trial->charts[c];local.id=0;
        AsmRun one={0};one.arena=scratch;one.charts=&local;one.n_charts=1;
        AsmField *seed=local.uv?AsmField_new_local(scratch,&one,NULL,NULL):NULL;
        if(seed && seed->vertex_offset[0]==SIZE_MAX)seed=NULL;
        int valid=seed && apa_chart_metric(seed,0);
        if(!valid){
            AsmChart *chart=trial->charts+c;
            /* Compact recovery corrects the cached original metric within
             * this transaction's deadline and workspace. A certified closed
             * source fan needs a UV cut, not more iterations on this topology. */
            if(opts->compact){
                int recovered=api_recover(seed,chart,trial->arena,opts,reason,trace);
                Arena_dispose(&scratch);
                if(!recovered && !strcmp(*reason,"intrinsic_source_needs_cut"))
                    recovered=api_cut_recover(trial,c,opts,reason,trace);
                if(!recovered)return 0;continue;
            }
            if(apr_expired(opts)){
                *reason="intrinsic_recovery_deadline";Arena_dispose(&scratch);return 0;
            }
            /* A cold re-flattening of a big chart is the admission's whole cost on the big rungs
             * and it does not succeed: measured 2026-09-17, the 21x5x5 and pherc343 spend 100-190 s
             * per candidate running 4,000 SLIM steps over 45-136k-face charts that end "proposal
             * refused", and the candidate is then refused anyway.  Above the cap the candidate is
             * refused on the spot.  The 4x5x5's admissions cold-recover nothing, so this is free
             * on the champion pile. */
            if(ASM_ADMISSION_RECOVER_MAX_FACES>0 && chart->nf>(size_t)ASM_ADMISSION_RECOVER_MAX_FACES){
                if(oversize)*oversize=chart->nf;
                *reason="candidate_metric_unrecoverable";
                Arena_dispose(&scratch);return 0;
            }
            float *uv=ARENA_ALLOC(trial->arena,2*chart->nv*sizeof(float));AsmFlattenStats st;
            if(AsmFlatten_chart(scratch,chart->xyz,chart->nv,chart->faces,chart->nf,.1,100,uv,&st)){
                *reason="intrinsic_recovery_failed";Arena_dispose(&scratch);return 0;
            }
            chart->uv=uv;chart->placed_uv=NULL;chart->pose_theta=chart->pose_x=chart->pose_y=0;
        }
        Arena_dispose(&scratch);
    }
    return 1;
}

static void apa_trace(const AsmField *f,const uint8_t *core,const char *phase,FILE *log)
{
    if(!log)return;
    for(size_t c=0;c<f->run->n_charts;c++)if(core[c] && f->vertex_offset[c]!=SIZE_MAX){
        double lo[2]={DBL_MAX,DBL_MAX},hi[2]={-DBL_MAX,-DBL_MAX},minor=DBL_MAX,major=0,area=0,good=0;size_t bad=0;
        for(size_t k=0;k<f->run->charts[c].nv;k++)for(int d=0;d<2;d++){double x=f->uv[2*(f->vertex_offset[c]+k)+d];lo[d]=fmin(lo[d],x);hi[d]=fmax(hi[d],x);}
        for(size_t k=0;k<f->run->charts[c].nf;k++){
            size_t t=f->face_offset[c]+k;double j[4],a,b,det;AsmField_face(f,t,f->uv,j,&a,&b,&det);
            minor=fmin(minor,a);major=fmax(major,b);area+=f->metric[t].area;
            bad+=!UvGuard_positive(f->uv,f->faces+3*t) || !AsmMetric_within25(a,b);
            if(a>=.9 && b<=1.1)good+=f->metric[t].area;
        }
        fprintf(log,"{\"phase\":\"%s\",\"chart\":%zu,\"uv_box\":[%.17g,%.17g,%.17g,%.17g],\"faces\":%zu,\"strict_bad\":%zu,\"sigma_min\":%.17g,\"sigma_max\":%.17g,\"within10\":%.17g}\n",
            phase,c,lo[0],lo[1],hi[0],hi[1],f->run->charts[c].nf,bad,minor,major,area?good/area:0);
    }
    for(size_t i=0;i<f->n_seams;i++)if(f->seams[i].required && (core[f->seams[i].a] || core[f->seams[i].b])){
        AsmSeamMeasure m;AsmField_seam(f,i,f->uv,&m);
        fprintf(log,"{\"phase\":\"%s\",\"relation\":%zu,\"a\":%d,\"b\":%d,\"pass\":%s,\"rms\":%.17g,\"supported_fraction\":%.17g,\"invalid\":%zu}\n",
            phase,f->seams[i].relation,f->seams[i].a,f->seams[i].b,m.pass?"true":"false",isfinite(m.rms)?m.rms:-1,m.mass?m.support/m.mass:0,m.invalid);
    }
    fflush(log);
}

static int apa_present(const AsmRun *run,int32_t a,int32_t b)
{return AsmChart_registered(run->charts+a) && AsmChart_registered(run->charts+b);}

#include "asm_repair_admission_contacts.inc"

/* The worst RMS over the required incident seams the qualification checks
 * (both charts present on the frontier); a seam without valid witnesses is
 * infinitely far.  0 when there is none. */
static double apa_worst_start_seam(const AsmField *f,const uint8_t *moving,int frontier)
{
    double worst=0;
    for(size_t i=0;i<f->n_seams;i++)if(f->seams[i].required && (moving[f->seams[i].a] || moving[f->seams[i].b]) &&
        (!frontier || apa_present(f->run,f->seams[i].a,f->seams[i].b))){
        AsmSeamMeasure m;AsmField_seam(f,i,f->uv,&m);
        double rms=isfinite(m.rms)?m.rms:INFINITY;
        if(rms>worst)worst=rms;
    }
    return worst;
}

static const char *apa_qualification(AsmField *f,const uint8_t *moving,const uint8_t *core,const double *reference,
                                     const uint8_t *active,FILE *trace,int frontier,int *complete)
{
    if(complete)*complete=1;
    /* Admission never inherits relaxed non-regression bands from a bad
     * seed. Every new original face and every incident hypothesis must pass. */
    for(size_t c=0;c<f->run->n_charts;c++)if(f->vertex_offset[c]!=SIZE_MAX){
        if(moving[c]){if(!apa_chart_metric(f,c))return "new_material_metric_failed";}
        else if(memcmp(reference+2*f->vertex_offset[c],f->uv+2*f->vertex_offset[c],2*f->run->charts[c].nv*sizeof(double)))return "fixed_context_changed";
    }
    for(size_t i=0;i<f->n_seams;i++)if(f->seams[i].required && (moving[f->seams[i].a] || moving[f->seams[i].b]) &&
        (!frontier || apa_present(f->run,f->seams[i].a,f->seams[i].b))){
        AsmSeamMeasure m;AsmField_seam(f,i,f->uv,&m);if(!m.pass)return "incident_source_seam_failed";
    }
    return apa_contact_certificate(f,moving,core,reference,active,trace,complete);
}

/* Contact-only correction for an otherwise qualified chart. This proposes
 * UVs only: the complete admission certificate and
 * immutable exterior query still decide whether any result is retained. */
static size_t apa_contact_polish(AsmField *f,const uint8_t *core,const AsmRepairOpts *opts,FILE *trace,int *complete)
{
    if(complete)*complete=0;
    if(apr_expired(opts))return 0;
    size_t chart=SIZE_MAX;
    for(size_t c=0;c<f->run->n_charts;c++)if(core[c]){
        if(chart!=SIZE_MAX)return 0;
        chart=c;
    }
    if(chart==SIZE_MAX)return 0;
    double started=ves_clock_sec(),deadline=started+.25;
    /* A bounded solve should spend time on useful contact
     * progress before refitting already qualified metrics and seams. The
     * vertex search returns early when it stalls, leaving joint work its
     * remaining allowance. Reserve the other share for moving neighbors:
     * a tiny core-only improvement does not prove it can finish alone.
     * Zero-iteration endpoint probes remain short. */
    int work_budget=opts->deadline>0 && (opts->iterations>0 || opts->clearance_iterations>0);
    if(work_budget)deadline=started+ASM_REPAIR_PRE_FIT_CLEARANCE_SHARE*(opts->deadline-started);
    else if(opts->deadline>0)deadline=fmin(deadline,opts->deadline);
    /* Only the new core changes in this correction. Stationary inherited
     * contacts must not dominate its physical progress measurement. */
    Arena_T arena=Arena_new();uint8_t *active=ARENA_CALLOC(arena,f->nf,1);
    memset(active+f->face_offset[chart],1,f->run->charts[chart].nf);
    AsmContacts_active_faces(f->contacts,active);
    AmrPolishStats stats={0};size_t steps=amr_vertex_polish_report(f,chart,deadline,&stats);
    AsmContacts_active_faces(f->contacts,NULL);Arena_dispose(&arena);
    if(complete)*complete=ves_clock_sec()<deadline;
    if(trace)fprintf(trace,"{\"phase\":\"contact_polish\",\"chart\":%zu,\"steps\":%zu,\"seconds\":%.17g,"
        "\"remaining_work_budget\":%s,\"allocated_seconds\":%.17g,"
        "\"stop\":\"%s\",\"rounds\":%zu,\"vertex_visits\":%zu,\"trials\":%zu,\"metric_refused\":%zu,"
        "\"seam_refused\":%zu,\"contact_queries\":%zu,\"contact_refused\":%zu,\"query_failed\":%zu,"
        "\"zero_gradient_visits\":%zu,\"candidate_vertices\":%zu,\"pair_increasing_steps\":%zu,"
        "\"focused_rounds\":%zu,\"boundary_candidates\":%zu,\"plateau_steps\":%zu,"
        "\"pairs_before\":%zu,\"pairs_after\":%zu,\"area_before\":%.17g,\"area_after\":%.17g}\n",
        chart,steps,ves_clock_sec()-started,work_budget?"true":"false",deadline-started,
        stats.stop,stats.rounds,stats.vertices,stats.trials,stats.metric_refused,
        stats.seam_refused,stats.contact_queries,stats.contact_refused,stats.query_failed,stats.zero_gradient,stats.candidates,stats.pair_increasing_steps,
        stats.focused_rounds,stats.boundary_candidates,stats.plateau_steps,
        stats.pairs_before,stats.pairs_after,stats.area_before,stats.area_after);
    return steps;
}

typedef struct ApaPolishOwners {
    const AsmField *field;
    const uint8_t *core,*moving;
    double *area;
} ApaPolishOwners;

static void apa_polish_owner(void *context,size_t a,size_t b,double area)
{
    ApaPolishOwners *p=context;
    size_t ca=(size_t)p->field->face_chart[a],cb=(size_t)p->field->face_chart[b];
    if(p->core[ca] && !p->core[cb] && p->moving[cb])p->area[cb]+=area;
    if(p->core[cb] && !p->core[ca] && p->moving[ca])p->area[ca]+=area;
}

/* A new chart can reach its metric limit while its movable collider still
 * has room. Apply the same certified vertex cleanup on measured collision
 * owners already authorized by the joint moving set. Each chart's active
 * contact mask includes all its local neighbors, not only the new material;
 * the caller rechecks the complete transaction and exterior afterward. */
static size_t apa_neighbor_polish(AsmField *f,const uint8_t *moving,const uint8_t *core,
                                   const AsmRepairOpts *opts,FILE *trace,int *complete)
{
    if(complete)*complete=1;
    if(!opts->admission_moving || apr_expired(opts))return 0;
    Arena_T arena=Arena_new();size_t nc=f->run->n_charts,steps=0;
    double *area=ARENA_CALLOC(arena,nc,sizeof *area);
    uint8_t *active=ARENA_CALLOC(arena,f->nf,1),*selected=ARENA_CALLOC(arena,nc,1);
    for(size_t c=0;c<nc;c++)if(core[c] && f->face_offset[c]!=SIZE_MAX)
        memset(active+f->face_offset[c],1,f->run->charts[c].nf);
    ApaPolishOwners owners={f,core,moving,area};AsmContactStats contacts={0};
    AsmContacts_active_faces(f->contacts,active);
    int rc=AsmContacts_measure(f->contacts,f->uv,NULL,apa_polish_owner,&owners,&contacts);
    AsmContacts_active_faces(f->contacts,NULL);
    if(rc || !contacts.complete){if(complete)*complete=0;Arena_dispose(&arena);return 0;}
    while(!apr_expired(opts)){
        size_t best=SIZE_MAX;
        for(size_t c=0;c<nc;c++)if(area[c]>0 && !(opts->fixed_charts && opts->fixed_charts[c]) &&
            (best==SIZE_MAX || area[c]>area[best]))best=c;
        if(best==SIZE_MAX)break;
        double blocked=area[best];area[best]=0;selected[best]=1;int finished=1;
        size_t taken=apa_contact_polish(f,selected,opts,trace,&finished);steps+=taken;selected[best]=0;
        if(complete && !finished)*complete=0;
        if(trace)fprintf(trace,"{\"phase\":\"neighbor_contact_polish\",\"chart\":%zu,\"blocking_new_material_area\":%.17g,\"steps\":%zu}\n",
                         best,blocked,taken);
    }
    Arena_dispose(&arena);return steps;
}

typedef struct ApaResidualPull {
    const AsmField *field;
    const uint8_t *core;
    double direction[2];
    size_t pairs,other;
} ApaResidualPull;

/* Area-weighted separation of the new chart from each residual collider:
 * from the collider triangle's centroid to the new triangle's. */
static void apa_residual_pull(void *context,size_t a,size_t b,double area)
{
    ApaResidualPull *p=context;const AsmField *f=p->field;
    int ca=p->core[f->face_chart[a]]!=0,cb=p->core[f->face_chart[b]]!=0;
    if(ca==cb){p->other++;return;}
    size_t mine=ca?a:b,theirs=ca?b:a;double d[2]={0,0};
    for(int k=0;k<3;k++)for(int q=0;q<2;q++)
        d[q]+=(f->uv[2*(size_t)f->faces[3*mine+(size_t)k]+(size_t)q]-f->uv[2*(size_t)f->faces[3*theirs+(size_t)k]+(size_t)q])/3;
    double length=hypot(d[0],d[1]);
    if(length>0 && isfinite(length))for(int q=0;q<2;q++)p->direction[q]+=area*d[q]/length;
    p->pairs++;
}

static const int apa_residual_nudge_on=ASM_REPAIR_RESIDUAL_NUDGE;
/* 1: refinement endpoints (frontier 4); 2: growth's rigid and repair attempts too */
static int apa_nudge_frontier(int frontier)
{return apa_residual_nudge_on>=2?frontier>=2:apa_residual_nudge_on==1 && frontier==4;}
/* Rigid residual clearance (see ASM_REPAIR_RESIDUAL_NUDGE). A solve that ends a
 * few new-material pairs short usually stalls there: its contact separation is
 * a finite penalty against metric and seam forces, and a contained triangle has
 * no area gradient. A rigid translation of the single new chart away from its
 * residual colliders cannot change its metric; each trial re-checks every seam
 * and the complete contact certificate, and a failed search restores the UVs. */
static const char *apa_residual_nudge(AsmField *f,const uint8_t *moving,const uint8_t *core,const double *reference,
                                      const uint8_t *active,const AsmRepairOpts *opts,FILE *trace,int frontier,
                                      int *complete,const char *refused)
{
    /* A solve leaves this search ASM_REPAIR_RESIDUAL_NUDGE_SEC of its quantum
     * (apa_attempt_anchor); past the deadline the endpoint changes nothing. */
    double until=ves_clock_sec()+ASM_REPAIR_RESIDUAL_NUDGE_SEC;
    if(opts->deadline>0)until=fmin(until,opts->deadline);
    size_t chart=SIZE_MAX;
    for(size_t c=0;c<f->run->n_charts;c++)if(core[c]){if(chart!=SIZE_MAX)return refused;chart=c;}
    if(chart==SIZE_MAX || f->vertex_offset[chart]==SIZE_MAX || apr_expired(opts))return refused;
    ApaResidualPull pull={f,core,{0,0},0,0};AsmContactStats stats={0};
    AsmContacts_active_faces(f->contacts,active);
    int rc=AsmContacts_measure(f->contacts,f->uv,NULL,apa_residual_pull,&pull,&stats);
    double length=hypot(pull.direction[0],pull.direction[1]);
    if(rc || !stats.complete || !pull.pairs || pull.pairs>APA_RESIDUAL_PAIRS || !(length>0) || !isfinite(length))return refused;
    double direction[2]={pull.direction[0]/length,pull.direction[1]/length};
    size_t first=2*f->vertex_offset[chart],count=2*f->run->charts[chart].nv;
    Arena_T arena=Arena_new();double *saved=ARENA_ALLOC(arena,count*sizeof(double));memcpy(saved,f->uv+first,count*sizeof(double));
    const double steps[]={.0625,.125,.25,.5,1,2};const char *reason=refused;size_t trials=0;double taken=0;
    for(size_t k=0;k<sizeof steps/sizeof *steps && ves_clock_sec()<until;k++){
        for(size_t i=0;i<count;i++)f->uv[first+i]=saved[i]+steps[k]*direction[i&1];
        trials++;int finished=0;
        reason=apa_qualification(f,moving,core,reference,active,NULL,frontier,&finished);
        if(!reason){taken=steps[k];if(complete)*complete=finished;break;}
        /* Farther only helps a contact refusal; a seam or an incomplete walk ends the search. */
        if(strcmp(reason,"new_material_contact_failed"))break;
    }
    if(reason){memcpy(f->uv+first,saved,count*sizeof(double));reason=refused;}
    if(trace)fprintf(trace,"{\"phase\":\"residual_nudge\",\"chart\":%zu,\"pairs\":%zu,\"direction\":[%.17g,%.17g],"
        "\"trials\":%zu,\"step\":%.17g,\"result\":\"%s\"}\n",chart,pull.pairs,direction[0],direction[1],trials,taken,
        reason?"refused":"qualified");
    Arena_dispose(&arena);return reason;
}

static const char *apa_qualified_endpoint(AsmField *f,const uint8_t *moving,const uint8_t *core,
                                         const double *reference,const uint8_t *active,const AsmRepairOpts *opts,
                                         FILE *trace,int frontier,int revisit,int *complete,int *polish_complete)
{
    *polish_complete=1;
    const char *refused=apa_qualification(f,moving,core,reference,active,trace,frontier,complete);
    /* A collider's accepted displacement can release the core from its
     * previous vertex/metric active set. Revisit both sides while a joint
     * transaction has time and makes progress. A single sweep leaves that
     * newly available direction unexplored. Every sweep retains the same
     * complete endpoint qualification; the exterior is checked at commit. */
    for(size_t pass=0;frontier==4 && refused && !strcmp(refused,"new_material_contact_failed");pass++){
        int finished=1;size_t core_steps=apa_contact_polish(f,core,opts,trace,&finished),neighbor_steps=0;
        *polish_complete&=finished;
        if(core_steps){
            apa_trace(f,moving,"after_polish",trace);
            refused=apa_qualification(f,moving,core,reference,active,trace,frontier,complete);
        }
        if(refused && !strcmp(refused,"new_material_contact_failed")){
            finished=1;neighbor_steps=apa_neighbor_polish(f,moving,core,opts,trace,&finished);
            if(neighbor_steps){
                apa_trace(f,moving,"after_neighbor_polish",trace);
                refused=apa_qualification(f,moving,core,reference,active,trace,frontier,complete);
            }
            *polish_complete&=finished;
        }
        if(trace)fprintf(trace,"{\"phase\":\"joint_contact_pass\",\"pass\":%zu,\"core_steps\":%zu,\"neighbor_steps\":%zu,\"result\":\"%s\"}\n",
                         pass,core_steps,neighbor_steps,refused?refused:"qualified");
        if(!revisit || !opts->admission_moving || !(opts->deadline>0) || apr_expired(opts) ||
           (!core_steps && !neighbor_steps))break;
    }
    if(apa_nudge_frontier(frontier) && refused && !strcmp(refused,"new_material_contact_failed")){
        refused=apa_residual_nudge(f,moving,core,reference,active,opts,trace,frontier,complete,refused);
        if(!refused)apa_trace(f,moving,"after_nudge",trace);
    }
    if(!refused && opts->admission_metric_margin>0)for(size_t c=0;c<f->run->n_charts;c++)
        if(moving[c] && !apa_chart_metric_margin(f,c,opts->admission_metric_margin))return "packing_metric_margin";
    return refused;
}

/* ---- the LOCAL admission field (2026-09-17) ---------------------------------------------------
 * A candidate used to copy the run and build AsmField_new over ALL registered charts before
 * solving with everything but its own component fixed: 21-33M vertices of metric, witnesses and
 * BVH per candidate, so the initial pass attempted 1 of 7,070 candidates on pherc343 and 4 of
 * 1,098 on the 21x5x5 inside its 600-s budget.  The field is now built over the candidate's
 * CONTEXT -- its own charts, every partner of a source relation incident to them, and the
 * registered charts whose uv box comes within reach of those partners -- and the material the
 * field leaves out is checked by probing the candidate's faces against the incumbent whole
 * field's tree (minus the charts the context pulled in) and against an overlay of the charts this
 * pass admitted earlier, which the incumbent does not contain.  Same predicate, same tolerance,
 * so `pairs == 0` still means what it meant. */

typedef struct ApaBox { double lo[2], hi[2]; int have; } ApaBox;

static void apa_chart_box(const AsmChart *ch,ApaBox *box)
{
    box->have=0;box->lo[0]=box->lo[1]=1e300;box->hi[0]=box->hi[1]=-1e300;
    if(!AsmChart_registered(ch) || !ch->nv)return;
    for(size_t v=0;v<ch->nv;v++){
        double p[2];AsmChart_point(ch,v,p);
        if(!isfinite(p[0]) || !isfinite(p[1]))continue;
        for(int d=0;d<2;d++){if(p[d]<box->lo[d])box->lo[d]=p[d];if(p[d]>box->hi[d])box->hi[d]=p[d];}
    }
    box->have=box->hi[0]>=box->lo[0];
}

static void apa_chart_boxes(const AsmRun *run,ApaBox *box)
{for(size_t c=0;c<run->n_charts;c++)apa_chart_box(run->charts+c,box+c);}

static size_t apa_context(Arena_T arena,const AsmRun *run,const uint8_t *core,const ApaBox *box,uint8_t *member,
                          double *reach_out,int *anchored)
{
    size_t nc=run->n_charts,count=0,partners=0;
    double extent=0.0;
    memset(member,0,nc);
    *anchored=0;
    for(size_t c=0;c<nc;c++)if(core[c]){
        member[c]=1;count++;
        const AsmChart *ch=run->charts+c;
        if(!ch->uv || !ch->nv)continue;
        double a[2]={1e300,1e300},b[2]={-1e300,-1e300};
        for(size_t v=0;v<ch->nv;v++)for(int d=0;d<2;d++){
            double x=ch->uv[2*v+d];
            if(x<a[d])a[d]=x;if(x>b[d])b[d]=x;
        }
        if(b[0]>=a[0])extent=fmax(extent,hypot(b[0]-a[0],b[1]-a[1]));
    }
    for(size_t i=0;i<run->n_rels;i++){
        const AsmRelation *r=run->rels+i;
        if(!apa_source(r) || r->a<0 || r->b<0 || (size_t)r->a>=nc || (size_t)r->b>=nc)continue;
        if(core[r->a]==core[r->b])continue;
        size_t other=(size_t)(core[r->a]?r->b:r->a);
        if(member[other])continue;
        member[other]=1;count++;partners++;
    }
    /* the fit puts the core beside its partners, so the reach is EACH partner's box grown by the
     * core's own extent plus the motion a bounded repair can add.  The union of the partners'
     * boxes would be the whole sheet whenever they are scattered, which is how the first arm
     * pulled 931-1,077 charts into a "local" field (2026-09-17). */
    double reach=extent+ASM_ADMISSION_CONTEXT_VOX;
    if(reach_out)*reach_out=reach;
    if(!partners)return count;
    /* the partners with a box, listed once so the halo test is O(partners) per chart */
    Arena_Mark mark=Arena_save(arena);
    size_t *anchor=ARENA_ALLOC(arena,(partners?partners:1)*sizeof(size_t)),np=0;
    for(size_t c=0;c<nc;c++)if(member[c] && !core[c] && box[c].have)anchor[np++]=c;
    if(!np){Arena_restore(arena,mark);return count;}
    *anchored=1;
    for(size_t c=0;c<nc;c++){
        if(member[c] || !box[c].have)continue;
        for(size_t q=0;q<np;q++){
            const ApaBox *b=box+anchor[q];
            if(box[c].lo[0]>b->hi[0]+reach || box[c].hi[0]<b->lo[0]-reach)continue;
            if(box[c].lo[1]>b->hi[1]+reach || box[c].hi[1]<b->lo[1]-reach)continue;
            member[c]=1;count++;break;
        }
    }
    Arena_restore(arena,mark);
    return count;
}

enum {APA_OVERLAY_LEVELS=32};
typedef struct ApaOverlayBatch {
    Arena_T arena;
    size_t nv,nf;
    double *uv;int32_t *faces,*face_chart,*material;
    AsmContacts *tree;
} ApaOverlayBatch;

typedef struct ApaOverlay {
    Arena_T arena;
    size_t nv,nf,cap_nv,cap_nf;
    double *uv;int32_t *faces,*face_chart,*material;
    AsmContacts *tree;
    double build_seconds;
    size_t rebuilds, rebuilt_faces;
    int batched,query_scaled_budget;
    /* Workers for the block trees' builds (<= 1: serial).  The parallel
     * build is the serial tree node for node (ac_parallel_build_control), so
     * this changes only how fast a replacement block is built. */
    int threads;
    ApaOverlayBatch *level[APA_OVERLAY_LEVELS];
    /* Non-null for transactions that may replace registered chart UVs.
     * The immutable incumbent excludes these charts; this tree owns their
     * current versions, including charts moved more than once. */
    uint8_t *changed_charts;
    size_t n_charts;
    /* Incremental replacement (ASM_REPAIR_OVERLAY_INCREMENTAL): chart_block[c]
     * is 1 + the level of the one block holding chart c's current version, 0
     * for none. A face of another block owned by c is a stale copy, skipped by
     * every query. Without it one rebuilt block holds every current version. */
    int incremental;
    uint8_t *chart_block,*scratch;
    size_t compactions;
} ApaOverlay;

static void apa_overlay_clear(ApaOverlay *ov)
{
    for(int i=0;i<APA_OVERLAY_LEVELS;i++)if(ov->level[i]){
        Arena_dispose(&ov->level[i]->arena);ov->level[i]=NULL;
    }
    if(ov->chart_block)memset(ov->chart_block,0,ov->n_charts);
}

/* A probe visits the stored tree once per external triangle. Its workload
 * includes both sides of that query, unlike the tree's self-contact walk.
 * Keep one finite allowance for the entire overlay, never one per block. */
static size_t apa_overlay_probe_budget(const ApaOverlay *ov,size_t ntri)
{
    size_t limit=(size_t)ASM_CONTACT_MAX_LEAF_PAIRS,budget=0;
    if(!ov || !ov->batched)return 0;
    for(int i=0;i<APA_OVERLAY_LEVELS;i++)if(ov->level[i]){
        size_t part=AsmContacts_budget(ov->level[i]->tree);
        if(!part || part>=limit-budget)budget=limit;
        else budget+=part;
    }
    if(ov->query_scaled_budget){
        size_t per=(size_t)ASM_CONTACT_MAX_LEAF_PAIRS_PER_LEAF;
        size_t query=ntri>limit/per?limit:ntri*per;
        if(query>budget)budget=query;
    }
    return budget;
}

/* Binary carries merge immutable blocks of earlier admissions. A face is
 * recopied only when its block doubles; unchanged larger blocks keep their
 * complete contact trees. Each block owns its arena so merged trees are
 * released, rather than retained for the lifetime of the admission pass. */
static int apa_overlay_batch_add(ApaOverlay *ov,const AsmRun *run,const uint8_t *core)
{
    size_t add_nv=0,add_nf=0;
    for(size_t c=0;c<run->n_charts;c++)if(core[c]){add_nv+=run->charts[c].nv;add_nf+=run->charts[c].nf;}
    if(!add_nf)return 0;
    if(ov->nv>(size_t)INT_MAX || add_nv>(size_t)INT_MAX-ov->nv || add_nf>SIZE_MAX-ov->nf)return -1;
    int level=0;while(level<APA_OVERLAY_LEVELS && ov->level[level])level++;
    if(level==APA_OVERLAY_LEVELS)return -1;
    size_t nv=add_nv,nf=add_nf;
    for(int i=0;i<level;i++){nv+=ov->level[i]->nv;nf+=ov->level[i]->nf;}
    if(nv>(size_t)INT_MAX || nf>SIZE_MAX/(3*sizeof(int32_t)))return -1;
    /* 4 MB chunks: a block's large arrays get chunks of their own size, and a
     * small block no longer commits a whole 64 MB default chunk. */
    ApaOverlayBatch *batch=ARENA_CALLOC(ov->arena,1,sizeof *batch);batch->arena=Arena_new_sized((size_t)4<<20);
    batch->nv=nv;batch->nf=nf;batch->uv=ARENA_ALLOC(batch->arena,2*nv*sizeof(double));
    batch->faces=ARENA_ALLOC(batch->arena,3*nf*sizeof(int32_t));batch->face_chart=ARENA_ALLOC(batch->arena,nf*sizeof(int32_t));
    size_t *offset=ARENA_ALLOC(batch->arena,run->n_charts*sizeof(size_t));
    for(size_t c=0;c<run->n_charts;c++)offset[c]=SIZE_MAX;
    size_t v0=0,t0=0;
    for(int i=0;i<level;i++){
        const ApaOverlayBatch *old=ov->level[i];
        memcpy(batch->uv+2*v0,old->uv,2*old->nv*sizeof(double));
        memcpy(batch->face_chart+t0,old->face_chart,old->nf*sizeof(int32_t));
        for(size_t k=0;k<3*old->nf;k++)batch->faces[3*t0+k]=(int32_t)(v0+(size_t)old->faces[k]);
        v0+=old->nv;t0+=old->nf;
    }
    for(size_t c=0;c<run->n_charts;c++)if(core[c]){
        const AsmChart *ch=run->charts+c;
        offset[c]=v0;
        for(size_t v=0;v<ch->nv;v++)AsmChart_point(ch,v,batch->uv+2*(v0+v));
        for(size_t t=0;t<ch->nf;t++){
            batch->face_chart[t0+t]=(int32_t)c;
            for(int k=0;k<3;k++){
                int source=(ch->flags&ASM_CHART_MIRROR)&&k?3-k:k;
                batch->faces[3*(t0+t)+k]=(int32_t)(v0+(size_t)ch->faces[3*t+source]);
            }
        }
        v0+=ch->nv;t0+=ch->nf;
    }
    if(AsmField_material_map(batch->arena,run,offset,nv,&batch->material)){Arena_dispose(&batch->arena);return -1;}
    v0=0;
    for(int i=0;i<level;i++){
        const ApaOverlayBatch *old=ov->level[i];
        if(old->material){
            if(!batch->material){batch->material=ARENA_ALLOC(batch->arena,nv*sizeof(int32_t));for(size_t v=0;v<nv;v++)batch->material[v]=(int32_t)v;}
            for(size_t v=0;v<old->nv;v++)batch->material[v0+v]=(int32_t)(v0+(size_t)old->material[v]);
        }
        v0+=old->nv;
    }
    batch->tree=AsmContacts_new_source(batch->arena,nv,batch->faces,nf,batch->uv,AsmContacts_material_tolerance(),
                                       ov->threads>1?ov->threads:1,batch->material);
    if(!batch->tree || AsmContacts_refit(batch->tree,batch->uv)){Arena_dispose(&batch->arena);return -1;}
    for(int i=0;i<level;i++){Arena_dispose(&ov->level[i]->arena);ov->level[i]=NULL;}
    ov->level[level]=batch;ov->nv+=add_nv;ov->nf+=add_nf;ov->rebuilds++;ov->rebuilt_faces+=nf;
    return 0;
}

static int apa_overlay_replace(ApaOverlay *ov,const AsmRun *run,const uint8_t *moving)
{
    if(!ov->batched || ov->n_charts!=run->n_charts)return -1;
    uint8_t *all=ARENA_ALLOC(ov->arena,run->n_charts);
    for(size_t c=0;c<run->n_charts;c++)all[c]=ov->changed_charts[c] || moving[c];
    ApaOverlay replacement={0};replacement.arena=ov->arena;replacement.batched=1;replacement.threads=ov->threads;
    /* Build before publishing. A refused replacement leaves both the old
     * delta tree and its incumbent exclusion mask intact. */
    if(apa_overlay_batch_add(&replacement,run,all))return -1;
    apa_overlay_clear(ov);ov->level[0]=replacement.level[0];
    ov->nv=replacement.nv;ov->nf=replacement.nf;
    ov->rebuilds+=replacement.rebuilds;ov->rebuilt_faces+=replacement.rebuilt_faces;
    memcpy(ov->changed_charts,all,run->n_charts);return 0;
}

/* Rebuild block `level` from the committed run over the charts in `take`:
 * build before publishing, so a refused build leaves every block intact. */
static int apa_overlay_rebuild_block(ApaOverlay *ov,const AsmRun *run,const uint8_t *take,ApaOverlayBatch **out)
{
    ApaOverlay replacement={0};replacement.arena=ov->arena;replacement.batched=1;replacement.threads=ov->threads;
    if(apa_overlay_batch_add(&replacement,run,take))return -1;
    ov->rebuilds+=replacement.rebuilds;ov->rebuilt_faces+=replacement.rebuilt_faces;
    *out=replacement.level[0];return 0;
}

/* Current (not stale) vertex and face totals, and the current faces of each block. */
static void apa_overlay_count(ApaOverlay *ov,const AsmRun *run,size_t current[APA_OVERLAY_LEVELS])
{
    ov->nv=ov->nf=0;memset(current,0,APA_OVERLAY_LEVELS*sizeof *current);
    for(size_t c=0;c<run->n_charts;c++)if(ov->chart_block[c]){
        ov->nv+=run->charts[c].nv;ov->nf+=run->charts[c].nf;current[ov->chart_block[c]-1]+=run->charts[c].nf;
    }
}

/* Binary carries over current versions only. The new block takes the moved
 * charts and the current versions held by every lower block it replaces, all
 * from the committed run (run->charts always holds each chart's current
 * placement). Copies left in higher blocks by charts moved again turn stale;
 * queries skip them until their block merges, or until more than half of a
 * block is stale and it is rebuilt from its current charts alone. */
static int apa_overlay_replace_incremental(ApaOverlay *ov,const AsmRun *run,const uint8_t *moving)
{
    size_t nc=run->n_charts,current[APA_OVERLAY_LEVELS];
    if(!ov->batched || ov->n_charts!=nc || !nc)return -1;
    if(!ov->chart_block){ov->chart_block=ARENA_CALLOC(ov->arena,nc,1);ov->scratch=ARENA_ALLOC(ov->arena,nc);}
    int level=0;while(level<APA_OVERLAY_LEVELS && ov->level[level])level++;
    if(level==APA_OVERLAY_LEVELS)return -1;
    uint8_t *take=ov->scratch;
    for(size_t c=0;c<nc;c++){
        int block=ov->chart_block[c];
        take[c]=(uint8_t)(moving[c] || (ov->changed_charts[c] && (!block || block<=level)));
    }
    ApaOverlayBatch *built=NULL;
    if(apa_overlay_rebuild_block(ov,run,take,&built))return -1;
    for(int i=0;i<level;i++)if(ov->level[i]){Arena_dispose(&ov->level[i]->arena);ov->level[i]=NULL;}
    ov->level[level]=built;
    for(size_t c=0;c<nc;c++)if(take[c]){ov->chart_block[c]=(uint8_t)(built?level+1:0);ov->changed_charts[c]=1;}
    apa_overlay_count(ov,run,current);
    for(int i=level+1;i<APA_OVERLAY_LEVELS;i++){
        ApaOverlayBatch *old=ov->level[i];
        if(!old || 2*current[i]>=old->nf)continue;
        for(size_t c=0;c<nc;c++)take[c]=(uint8_t)(ov->chart_block[c]==i+1);
        ApaOverlayBatch *compact=NULL;
        /* A refused compaction keeps the complete (partly stale) block. */
        if(apa_overlay_rebuild_block(ov,run,take,&compact))continue;
        Arena_dispose(&old->arena);ov->level[i]=compact;ov->compactions++;
        if(!compact)for(size_t c=0;c<nc;c++)if(take[c])ov->chart_block[c]=0;
    }
    apa_overlay_count(ov,run,current);
    return 0;
}

static int apa_overlay_add(ApaOverlay *ov,const AsmRun *run,const uint8_t *core)
{
    double started=ves_clock_sec();
    if(ov->changed_charts){
        int rc=ov->incremental?apa_overlay_replace_incremental(ov,run,core):apa_overlay_replace(ov,run,core);
        ov->build_seconds+=ves_clock_sec()-started;return rc;
    }
    if(ov->batched){int rc=apa_overlay_batch_add(ov,run,core);ov->build_seconds+=ves_clock_sec()-started;return rc;}
    size_t add_nv=0,add_nf=0;
    for(size_t c=0;c<run->n_charts;c++)if(core[c]){add_nv+=run->charts[c].nv;add_nf+=run->charts[c].nf;}
    if(!add_nf)return 0;
    if(ov->nv+add_nv>ov->cap_nv || ov->nf+add_nf>ov->cap_nf){
        size_t cv=ov->cap_nv?ov->cap_nv:1024,cf=ov->cap_nf?ov->cap_nf:1024;
        while(cv<ov->nv+add_nv)cv*=2;
        while(cf<ov->nf+add_nf)cf*=2;
        double *uv=ARENA_ALLOC(ov->arena,2*cv*sizeof(double));
        int32_t *fa=ARENA_ALLOC(ov->arena,3*cf*sizeof(int32_t)),*fc=ARENA_ALLOC(ov->arena,cf*sizeof(int32_t));
        if(ov->nv)memcpy(uv,ov->uv,2*ov->nv*sizeof(double));
        if(ov->nf){memcpy(fa,ov->faces,3*ov->nf*sizeof(int32_t));memcpy(fc,ov->face_chart,ov->nf*sizeof(int32_t));}
        ov->uv=uv;ov->faces=fa;ov->face_chart=fc;ov->cap_nv=cv;ov->cap_nf=cf;
    }
    size_t old_nv=ov->nv;int32_t *old_material=ov->material;
    size_t *offset=ARENA_ALLOC(ov->arena,run->n_charts*sizeof(size_t));
    for(size_t c=0;c<run->n_charts;c++)offset[c]=SIZE_MAX;
    for(size_t c=0;c<run->n_charts;c++)if(core[c]){
        const AsmChart *ch=run->charts+c;size_t v0=ov->nv;
        offset[c]=v0;
        for(size_t v=0;v<ch->nv;v++)AsmChart_point(ch,v,ov->uv+2*(v0+v));
        for(size_t t=0;t<ch->nf;t++){
            ov->face_chart[ov->nf+t]=(int32_t)c;
            for(int k=0;k<3;k++){
                int source=(ch->flags&ASM_CHART_MIRROR)&&k?3-k:k;
                ov->faces[3*(ov->nf+t)+k]=(int32_t)(v0+(size_t)ch->faces[3*t+source]);
            }
        }
        ov->nv+=ch->nv;ov->nf+=ch->nf;
    }
    if(AsmField_material_map(ov->arena,run,offset,ov->nv,&ov->material))return -1;
    if(old_material){
        if(!ov->material){ov->material=ARENA_ALLOC(ov->arena,ov->nv*sizeof(int32_t));for(size_t v=0;v<ov->nv;v++)ov->material[v]=(int32_t)v;}
        memcpy(ov->material,old_material,old_nv*sizeof(int32_t));
    }
    ov->tree=AsmContacts_new_source(ov->arena,ov->nv,ov->faces,ov->nf,ov->uv,AsmContacts_material_tolerance(),1,ov->material);
    ov->build_seconds+=ves_clock_sec()-started;ov->rebuilds++;ov->rebuilt_faces+=ov->nf;
    return ov->tree?0:-1;
}

typedef struct ApaObstacleProbe {
    const int32_t *owner;
    uint8_t *colliders;
} ApaObstacleProbe;

static void apa_obstacle_visit(void *context,size_t a,size_t b,double area)
{
    ApaObstacleProbe *probe=context;(void)a;(void)area;
    probe->colliders[probe->owner[b]]=1;
}

typedef struct ApaInitialContacts {
    const AsmField *field;
    const uint8_t *core,*moving;
    size_t blocking_pairs;
    double blocking_area;
    size_t far_pairs;   /* blocking pairs whose triangles lie apart in 3-D */
    double far_area;
} ApaInitialContacts;
/* The original 3-D point of field face t under uv point p: p's barycentric
 * coordinates in the face's current uv triangle, clamped to the triangle,
 * applied to its xyz.  A degenerate uv triangle gives its centroid. */
static void apa_face_point(const AsmField *f,size_t t,const double p[2],double out[3])
{
    const int32_t *v=f->faces+3*t;double l[3]={1.0/3,1.0/3,1.0/3};
    double x0=f->uv[2*(size_t)v[0]],y0=f->uv[2*(size_t)v[0]+1];
    double x1=f->uv[2*(size_t)v[1]],y1=f->uv[2*(size_t)v[1]+1];
    double x2=f->uv[2*(size_t)v[2]],y2=f->uv[2*(size_t)v[2]+1];
    double det=(x1-x0)*(y2-y0)-(x2-x0)*(y1-y0);
    if(fabs(det)>1e-12 && isfinite(det)){
        l[1]=((p[0]-x0)*(y2-y0)-(x2-x0)*(p[1]-y0))/det;
        l[2]=((x1-x0)*(p[1]-y0)-(p[0]-x0)*(y1-y0))/det;
        l[0]=1-l[1]-l[2];
        double sum=0;for(int k=0;k<3;k++){l[k]=fmin(1.0,fmax(0.0,l[k]));sum+=l[k];}
        if(sum>0)for(int k=0;k<3;k++)l[k]/=sum;else l[0]=l[1]=l[2]=1.0/3;
    }
    for(int d=0;d<3;d++){
        out[d]=0;for(int k=0;k<3;k++)out[d]+=l[k]*f->xyz[3*(size_t)v[k]+(size_t)d];
    }
}
/* 3-D gap between two faces' materials at the same uv spot (the midpoint of
 * their uv centroids): independent of triangle size, unlike a centroid
 * distance, which is a full triangle apart on a coarse same-place overlap. */
static double apa_face_gap(const AsmField *f,size_t a,size_t b)
{
    double p[2]={0,0},qa[3],qb[3],squared=0;
    for(int j=0;j<3;j++)for(int d=0;d<2;d++)
        p[d]+=f->uv[2*(size_t)f->faces[3*a+(size_t)j]+(size_t)d]+f->uv[2*(size_t)f->faces[3*b+(size_t)j]+(size_t)d];
    p[0]/=6;p[1]/=6;
    apa_face_point(f,a,p,qa);apa_face_point(f,b,p,qb);
    for(int d=0;d<3;d++)squared+=(qa[d]-qb[d])*(qa[d]-qb[d]);
    return sqrt(squared);
}
static void apa_initial_contact(void *context,size_t a,size_t b,double area)
{
    ApaInitialContacts *s=context;int32_t ca=s->field->face_chart[a],cb=s->field->face_chart[b];
    /* The saturation shortcut concerns new material against immovable
     * obstacles (or against itself). Existing neighbors explicitly selected
     * for joint repair must not inflate a ratio divided by candidate faces. */
    if((s->core[ca] && (!s->moving[cb] || s->core[cb])) ||
       (s->core[cb] && (!s->moving[ca] || s->core[ca]))){
        s->blocking_pairs++;s->blocking_area+=area;
        /* The same place on the sheet overlaps within the discovery pass's
         * same-surface distance; anything farther is material elsewhere. */
        if(apa_face_gap(s->field,a,b)>ASM_DISCOVER_SAME_SURFACE_VOX){s->far_pairs++;s->far_area+=area;}
    }
}
static const int apa_area_3d_test=ASM_ADMISSION_AREA_3D_TEST;
/* The pre-solve shortcut against fixed material (see ASM_ADMISSION_AREA_3D_TEST):
 * 1 = the candidate starts where no bounded local solve can clear it. */
static int apa_start_hopeless(const ApaInitialContacts *start,double candidate_area,size_t candidate_faces)
{
    if(!candidate_faces)return 0;
    if(!apa_area_3d_test || !(candidate_area>0))
        return (double)start->blocking_pairs>ASM_ADMISSION_MAX_CONTACT_RATIO*(double)candidate_faces;
    return start->far_area>ASM_ADMISSION_MAX_FAR_OVERLAP*candidate_area ||
           start->blocking_area>ASM_ADMISSION_MAX_PLACED_OVERLAP*candidate_area;
}

static int apa_external_query(const AsmField *f,const uint8_t *core,const uint8_t *member,
                              const AsmField *incumbent,const ApaOverlay *ov,Arena_T arena,
                              size_t *pairs,double *area,int *complete,int prepared,uint8_t *colliders)
{
    *pairs=0;*area=0.0;*complete=1;
    size_t ntri=0;
    for(size_t c=0;c<f->run->n_charts;c++)if(core[c] && f->face_offset[c]!=SIZE_MAX)ntri+=f->run->charts[c].nf;
    if(!ntri)return 0;
    Arena_Mark mark=Arena_save(arena);
    size_t *query=ARENA_ALLOC(arena,ntri*sizeof(size_t)),at=0;
    for(size_t c=0;c<f->run->n_charts;c++)if(core[c] && f->face_offset[c]!=SIZE_MAX){
        const AsmChart *ch=f->run->charts+c;
        for(size_t t=0;t<ch->nf;t++){
            query[at++]=f->face_offset[c]+t;
        }
    }
    int rc=0;
    if(incumbent && incumbent->nf && incumbent->contacts){
        uint8_t *ex=ARENA_CALLOC(arena,incumbent->nf,1);
        for(size_t c=0;c<incumbent->run->n_charts;c++)
            if((member[c] || (ov && ov->changed_charts && ov->changed_charts[c])) && incumbent->face_offset[c]!=SIZE_MAX)
                memset(ex+incumbent->face_offset[c],1,incumbent->run->charts[c].nf);
        AsmContactStats st={0};
        ApaObstacleProbe probe={incumbent->face_chart,colliders};
        int r=prepared?AsmContacts_probe_mesh_prepared(incumbent->contacts,f->contacts,f->uv,at,query,ex,colliders?apa_obstacle_visit:NULL,&probe,&st):
            AsmContacts_probe_mesh(incumbent->contacts,incumbent->uv,f->contacts,f->uv,at,query,ex,colliders?apa_obstacle_visit:NULL,&probe,&st);
        if(r==-2)*complete=0;else if(r)rc=-1;
        *pairs+=st.pairs;*area+=st.area;
    }
    /* A small immutable block must not impose its own tiny walk allowance
     * on a query against the complete overlay. Share the collection's finite
     * allowance across blocks, without multiplying the global resource cap. */
    size_t remaining=apa_overlay_probe_budget(ov,ntri);
    if(!rc && ov && ov->nf)for(int i=0;i<(ov->batched?APA_OVERLAY_LEVELS:1) && !rc;i++){
        const ApaOverlayBatch *batch=ov->batched?ov->level[i]:NULL;
        if(ov->batched && !batch)continue;
        size_t nf=batch?batch->nf:ov->nf;
        const int32_t *labels=batch?batch->face_chart:ov->face_chart;
        AsmContacts *tree=batch?batch->tree:ov->tree;
        if(!tree){rc=-1;break;}
        if(batch && !remaining){*complete=0;break;}
        size_t saved_budget=AsmContacts_budget(tree);
        if(batch)AsmContacts_leaf_pair_budget(tree,remaining);
        Arena_Mark block=Arena_save(arena);
        uint8_t *ex=ARENA_CALLOC(arena,nf,1);
        /* An incremental block also holds stale copies of charts whose
         * current version has moved to another block. */
        const uint8_t *owner=batch?ov->chart_block:NULL;
        for(size_t t=0;t<nf;t++){
            size_t c=(size_t)labels[t];
            if(member[c] || (owner && owner[c]!=(uint8_t)(i+1)))ex[t]=1;
        }
        AsmContactStats st={0};
        ApaObstacleProbe probe={labels,colliders};
        int r=batch?AsmContacts_probe_mesh_prepared(tree,f->contacts,f->uv,at,query,ex,colliders?apa_obstacle_visit:NULL,&probe,&st):
            AsmContacts_probe_mesh(tree,ov->uv,f->contacts,f->uv,at,query,ex,colliders?apa_obstacle_visit:NULL,&probe,&st);
        if(batch){
            AsmContacts_leaf_pair_budget(tree,saved_budget);
            remaining=st.leaf_pairs>=remaining?0:remaining-st.leaf_pairs;
        }
        if(r==-2){
            *complete=0;
            fprintf(stderr,"[admission external] overlay probe stopped: %zu query triangles, %zu stored faces, %zu visited leaf pairs, query-scaled %d\n",
                    ntri,nf,st.leaf_pairs,ov->query_scaled_budget);
        }else if(r)rc=-1;
        *pairs+=st.pairs;*area+=st.area;
        Arena_restore(arena,block);
        if(!*complete)break;
    }
    Arena_restore(arena,mark);
    return rc;
}

static int apa_external_pairs(const AsmField *f,const uint8_t *core,const uint8_t *member,
                              const AsmField *incumbent,const ApaOverlay *ov,Arena_T arena,
                              size_t *pairs,double *area,int *complete,int prepared)
{return apa_external_query(f,core,member,incumbent,ov,arena,pairs,area,complete,prepared,NULL);}

/* Add only owners of measured entry collisions. Their original geometry and
 * current UVs remain fixed; all other incumbent faces still receive the full
 * exterior certificate. Partial walks and oversized fields are refusals. */
static AsmField *apa_obstacle_field(AsmField *f,const uint8_t *core,uint8_t *member,
                                    const AsmField *incumbent,const ApaOverlay *overlay,
                                    const AsmRepairOpts *opts,const char **reason,FILE *trace)
{
    Arena_T arena=f->arena;size_t nc=f->run->n_charts,pairs=0,added=0,faces=f->nf;
    uint8_t *colliders=ARENA_CALLOC(arena,nc,1);double area=0;int complete=0;
    if(apa_external_query(f,core,member,incumbent,overlay,arena,&pairs,&area,&complete,1,colliders)){
        *reason="external_probe_failed";return NULL;
    }
    if(!complete){*reason="contact_walk_over_budget";return NULL;}
    for(size_t c=0;c<nc;c++)if(colliders[c] && !member[c]){
        faces+=f->run->charts[c].nf;added++;
        if(trace)fprintf(trace,"{\"phase\":\"obstacle_owner\",\"chart\":%zu,\"faces\":%zu}\n",c,f->run->charts[c].nf);
    }
    if(trace)fprintf(trace,"{\"phase\":\"obstacle_intake\",\"charts\":%zu,\"field_faces\":%zu,\"pairs\":%zu,\"area\":%.17g}\n",added,faces,pairs,area);
    if(!added)return f;
    if(faces>opts->patch_field_faces){*reason="obstacle_resource_cap";return NULL;}
    for(size_t c=0;c<nc;c++)if(colliders[c])member[c]=1;
    size_t workspace=apr_member_workspace(f->run,member);
    if(trace)fprintf(trace,"{\"phase\":\"workspace_preflight\",\"context\":\"obstacles\",\"estimated_bytes\":%zu,\"budget_bytes\":%zu,\"field_faces\":%zu}\n",workspace,opts->patch_bytes,faces);
    if(!apr_workspace_check(opts,workspace)){*reason="workspace_memory_budget";return NULL;}
    AsmField *expanded=AsmField_new_local(arena,f->run,member,f);
    if(!expanded || AsmField_build_contacts(expanded)){*reason="obstacle_field_failed";return NULL;}
    return expanded;
}

typedef struct ApaReport {
    size_t context_charts,field_vertices,field_faces,local_pairs,external_pairs;
    int    local_field,endpoint_export_failed,cache_hit;
    double sec,solver_seconds,solver_wait_seconds;
} ApaReport;

/* Feedback observes a refused transaction before its scratch field is
 * released. It never changes geometry, obligations or the acceptance rule. */
typedef struct ApaFeedback {
    size_t count;
    double *contact_area,*seam_ratio;
    uint8_t *context;
    double deadline;
    int complete,observe_only;
    double violation,material_area,blocking_material_area;
} ApaFeedback;

/* Retaining a refused candidate is not publication. Every retry still checks
 * the original incumbent, all incident seams and the complete exterior. */
typedef struct ApaContinuation {
    Arena_T arena;
    size_t count;
    double **uv;
    size_t *vertices;
    uint8_t *mirror;
    int changed;
} ApaContinuation;

static void apa_continuation_clear(ApaContinuation *state)
{
    if(state->arena)Arena_dispose(&state->arena);
    memset(state,0,sizeof *state);
}

static int apa_continuation_apply(const ApaContinuation *state,AsmField *f,const uint8_t *moving)
{
    if(!state || !state->arena || state->count!=f->run->n_charts)return 0;
    /* Validate every entry before copying any of them. A changed source
     * context invalidates the whole candidate, never just one endpoint. */
    for(size_t c=0;c<state->count;c++)if(state->uv[c] &&
       (!moving[c] || f->vertex_offset[c]==SIZE_MAX || state->vertices[c]!=f->run->charts[c].nv ||
        state->mirror[c]!=!!(f->run->charts[c].flags&ASM_CHART_MIRROR)))return 0;
    for(size_t c=0;c<state->count;c++)if(state->uv[c])
        memcpy(f->uv+2*f->vertex_offset[c],state->uv[c],2*state->vertices[c]*sizeof(double));
    return 1;
}

static void apa_continuation_save(ApaContinuation *state,const AsmField *f,const uint8_t *moving,
                                   const double *reference)
{
    state->changed=0;
    for(size_t c=0;c<f->run->n_charts;c++)if(moving[c] &&
       (f->vertex_offset[c]==SIZE_MAX || !apa_chart_metric(f,c)))return;
    int changed=0;
    for(size_t c=0;c<f->run->n_charts;c++)if(moving[c]){
        const double *old=state->arena && c<state->count && state->uv[c]?state->uv[c]:reference+2*f->vertex_offset[c];
        if(memcmp(old,f->uv+2*f->vertex_offset[c],2*f->run->charts[c].nv*sizeof(double)))changed=1;
    }
    if(!changed)return;
    /* A continuation lives across the pass, hundreds at once in a large region:
     * size its arena to its content (a default arena commits a 64 MB chunk). */
    size_t bytes_needed=f->run->n_charts*(sizeof(double *)+sizeof(size_t)+1)+64;
    for(size_t c=0;c<f->run->n_charts;c++)if(moving[c])bytes_needed+=2*f->run->charts[c].nv*sizeof(double)+32;
    ApaContinuation next={0};next.arena=Arena_new_sized(bytes_needed);next.count=f->run->n_charts;next.changed=1;
    next.uv=ARENA_CALLOC(next.arena,next.count,sizeof *next.uv);
    next.vertices=ARENA_CALLOC(next.arena,next.count,sizeof *next.vertices);
    next.mirror=ARENA_CALLOC(next.arena,next.count,1);
    for(size_t c=0;c<next.count;c++)if(moving[c]){
        size_t bytes=2*f->run->charts[c].nv*sizeof(double);
        next.uv[c]=ARENA_ALLOC(next.arena,bytes);next.vertices[c]=f->run->charts[c].nv;
        next.mirror[c]=(uint8_t)!!(f->run->charts[c].flags&ASM_CHART_MIRROR);
        memcpy(next.uv[c],f->uv+2*f->vertex_offset[c],bytes);
    }
    apa_continuation_clear(state);*state=next;
}
typedef struct ApaFeedbackContacts {
    const AsmField *field;
    const uint8_t *moving;
    ApaFeedback *feedback;
    ApaContactCertificate *certificate;
} ApaFeedbackContacts;

static void apa_feedback_contact(void *context,size_t a,size_t b,double area)
{
    ApaFeedbackContacts *p=context;
    apa_contact_visit(p->certificate,a,b,area);
    if(p->certificate->before)return;
    size_t ca=(size_t)p->field->face_chart[a],cb=(size_t)p->field->face_chart[b];
    if(p->moving[ca] || p->moving[cb])p->feedback->material_area+=area;
    if(p->moving[ca] && !p->moving[cb])p->feedback->contact_area[cb]+=area;
    if(p->moving[cb] && !p->moving[ca])p->feedback->contact_area[ca]+=area;
}

static void apa_feedback(const AsmField *f,const uint8_t *core,const uint8_t *moving,const double *reference,const AsmField *incumbent,
                         const ApaOverlay *overlay,Arena_T arena,ApaFeedback *out)
{
    if(!out || !f || !reference || !f->contacts || out->count!=f->run->n_charts)return;
    size_t nc=out->count,largest=1;double source_area=0;
    /* A wider context is not geometric progress. Keep the normalization
     * tied to the original material requested by this transaction. */
    for(size_t c=0;c<nc;c++)if(core[c])source_area+=f->run->charts[c].area3d;
    for(size_t i=0;i<f->n_seams;i++)if(f->seams[i].count>largest)largest=f->seams[i].count;
    ApcOrder *samples=ARENA_ALLOC(arena,largest*sizeof *samples);
    for(size_t c=0;c<nc;c++)if(f->vertex_offset[c]!=SIZE_MAX)out->context[c]=1;
    for(size_t i=0;i<f->n_seams;i++){
        const AsmFieldSeam *s=f->seams+i;
        if(!s->required || (!moving[s->a] && !moving[s->b]) || !apa_present(f->run,s->a,s->b))continue;
        AsmSeamMeasure m;
        AsmField_seam(f,i,f->uv,&m);
        if(!s->complete || !s->source_run || !s->parity || m.invalid || !(m.mass>0) || !isfinite(m.rms) ||
           (!s->original_cut && (m.observations<6 || m.mass<8))){out->violation=INFINITY;continue;}
        double threshold=0;
        for(size_t k=0;k<s->count;k++){
            double a[2],b[2];AsmField_residual(f->witness+s->first+k,0,f->uv,a);AsmField_residual(f->witness+s->first+k,1,f->uv,b);
            samples[k]=(ApcOrder){s->first+k,fmax(hypot(a[0],a[1]),hypot(b[0],b[1]))};
        }
        qsort(samples,s->count,sizeof *samples,apc_order);double mass=0,needed=(s->original_cut?1:.9)*m.mass;
        for(size_t k=0;k<s->count && mass<needed;k++){
            mass+=f->witness[samples[k].id].weight;threshold=samples[k].error;
        }
        double ratio=fmax(m.rms,threshold)/ASM_SEAM_TOLERANCE_VOX;
        if(!m.pass)out->violation+=fmax(0,ratio-1);
        if(!!moving[s->a]!=!!moving[s->b]){
            size_t fixed=(size_t)(moving[s->a]?s->b:s->a);
            out->seam_ratio[fixed]=fmax(out->seam_ratio[fixed],ratio);
        }
    }
    if(out->deadline>0 && ves_clock_sec()>=out->deadline)return;
    AsmContactStats contacts={0};ApaContactCertificate certificate={0};
    certificate.arena=arena;certificate.field=f;certificate.core=core;
    certificate.reference=reference;certificate.okay=1;
    ApaFeedbackContacts probe={f,moving,out,&certificate};
    uint8_t *active=ARENA_CALLOC(arena,f->nf,1);int neighbors=0;
    for(size_t c=0;c<nc;c++)if(moving[c] && f->face_offset[c]!=SIZE_MAX){
        memset(active+f->face_offset[c],1,f->run->charts[c].nf);neighbors|=!core[c];
    }
    AsmContacts_active_faces(f->contacts,active);
    int measured=0;
    if(neighbors){
        certificate.before=1;
        measured=AsmContacts_measure(f->contacts,reference,NULL,apa_feedback_contact,&probe,&contacts);
        if(!contacts.complete || !certificate.okay)measured=-1;
        certificate.before=0;
    }
    if(!measured)measured=AsmContacts_measure(f->contacts,f->uv,NULL,apa_feedback_contact,&probe,&contacts);
    AsmContacts_active_faces(f->contacts,NULL);
    if(measured || !contacts.complete || !certificate.okay)return;
    size_t pairs=0;double area=0;int complete=1;
    uint8_t *outside=ARENA_CALLOC(arena,nc,1);
    if(apa_external_query(f,moving,out->context,incumbent,overlay,arena,&pairs,&area,&complete,1,outside) || !complete)return;
    out->material_area+=area;
    out->blocking_material_area=apa_blocking_contact_area(&certificate)+area;
    out->violation+=sqrt(out->blocking_material_area/fmax(source_area,DBL_MIN));
    for(size_t c=0;c<nc;c++)if(outside[c]){
        /* The external query supplies a set, not a per-chart area. Positive
         * membership is sufficient to add a real obstacle before seam-only
         * expansion; do not invent an overlap-area measurement. */
        out->context[c]=1;
        if(!(out->contact_area[c]>0))out->contact_area[c]=DBL_MIN;
    }
    out->complete=1;
}

/* A region owns its field, certificates and immutable exterior. Only the
 * legacy sparse/QP solver has shared scratch state. Charge queue waiting to
 * the region's wall budget without consuming its candidate work allowance. */
static int apa_solve_owned(AsmField *field,AsmRepairOpts *local,AsmRepairOpts *request,
                           AsmRepairStats *stats,ApaReport *report)
{
    double waiting=ves_clock_sec(),started=waiting,waited=0;int rc=0;
    if(request->admission_solver_deadline>0){
#pragma omp critical(asm_region_deformation)
        {
            started=ves_clock_sec();waited=started-waiting;
            request->deadline=fmin(request->admission_solver_deadline,request->deadline+waited);
            local->deadline=request->deadline;
            rc=apr_solve(field,local,stats);
            if(report)report->solver_seconds+=ves_clock_sec()-started;
        }
    }else{
        rc=apr_solve(field,local,stats);
        if(report)report->solver_seconds+=ves_clock_sec()-started;
    }
    if(report)report->solver_wait_seconds+=waited;
    return rc;
}

/* A proof of infeasibility for a SINGLE moving chart, not a new acceptance
 * threshold. Every admitted triangle has sigma_max <= 1.25, so a source-edge
 * path of length d bounds its UV endpoint distance by 1.25*d. A seam's
 * weighted mean residual is at most the shared certified seam RMS. Thus two
 * fixed target means must fit inside their path-radius sum plus twice that
 * bound.
 * Extrapolation uses absolute affine coefficients, and disconnected or
 * numerically indeterminate paths are inconclusive. No candidate is accepted
 * by this bound; it only avoids solving a provably incompatible fixed frame. */
static int apa_incompatible_fixed_anchors(const AsmField *f,const uint8_t *core,FILE *trace)
{
    size_t chart=SIZE_MAX,ncore=0,nseams=0;
    for(size_t c=0;c<f->run->n_charts;c++)if(core[c]){chart=c;ncore++;}
    if(ncore!=1 || f->vertex_offset[chart]==SIZE_MAX)return 0;
    for(size_t i=0;i<f->n_seams;i++)if(f->seams[i].required && f->seams[i].complete &&
        (f->seams[i].a==(int32_t)chart || f->seams[i].b==(int32_t)chart))nseams++;
    if(nseams<2)return 0;
    Arena_T arena=Arena_new();const AsmChart *ch=f->run->charts+chart;int incompatible=0;
    size_t nv=ch->nv,base=f->vertex_offset[chart],root=0;
    /* A boundary root keeps the bound useful without claiming a Euclidean
     * diameter is an upper bound on a folded source surface's path length. */
    for(size_t v=0;v<nv;v++)if(ch->boundary && ch->boundary[v]){root=v;break;}
    size_t *offset=ARENA_CALLOC(arena,nv+1,sizeof(size_t));
    for(size_t t=0;t<ch->nf;t++)for(int k=0;k<3;k++)offset[(size_t)ch->faces[3*t+k]+1]++;
    for(size_t v=0;v<nv;v++)offset[v+1]+=offset[v];
    size_t *incident=ARENA_ALLOC(arena,offset[nv]*sizeof(size_t)),*fill=ARENA_CALLOC(arena,nv,sizeof(size_t));
    for(size_t t=0;t<ch->nf;t++)for(int k=0;k<3;k++){
        size_t v=(size_t)ch->faces[3*t+k];incident[offset[v]+fill[v]++]=t;
    }
    double *distance=ARENA_ALLOC(arena,nv*sizeof(double));size_t *heap=ARENA_ALLOC(arena,nv*sizeof(size_t)),count=1;
    int32_t *position=ARENA_ALLOC(arena,nv*sizeof(int32_t));
    for(size_t v=0;v<nv;v++){distance[v]=DBL_MAX;position[v]=-1;}
    distance[root]=0;heap[0]=root;position[root]=0;
    while(count){
        size_t v=acl_heap_pop(heap,position,&count,distance);
        for(size_t t=offset[v];t<offset[v+1];t++)for(int k=0;k<3;k++){
            size_t w=(size_t)ch->faces[3*incident[t]+k];if(position[w]==-2)continue;
            double square=0;
            for(int d=0;d<3;d++){double delta=(double)ch->xyz[3*w+d]-ch->xyz[3*v+d];square+=delta*delta;}
            double next=distance[v]+sqrt(square);if(!isfinite(next) || next>=distance[w])continue;
            distance[w]=next;
            if(position[w]<0){position[w]=(int32_t)count;heap[count++]=w;}
            acl_heap_up(heap,position,(size_t)position[w],distance);
        }
    }
    double origin[2]={f->uv[2*base],f->uv[2*base+1]};
    double *center=ARENA_CALLOC(arena,2*nseams,sizeof(double)),*radius=ARENA_CALLOC(arena,nseams,sizeof(double));
    size_t *relation=ARENA_ALLOC(arena,nseams*sizeof(size_t)),at=0;
    for(size_t i=0;i<f->n_seams;i++){
        const AsmFieldSeam *s=f->seams+i;
        if(!s->required || !s->complete || (s->a!=(int32_t)chart && s->b!=(int32_t)chart))continue;
        double mass=0;
        for(size_t k=0;k<s->count;k++)for(int d=0;d<2;d++){
            const AsmWitness *w=f->witness+s->first+k;double own=0,bound=0,target[2]={0},sum=0;
            if(!w->valid || !(w->weight>0))goto done;
            for(int q=0;q<3;q++)sum+=w->coefficient[d][q];
            for(int q=0;q<4;q++){
                size_t v=(size_t)w->vertex[d][q];double coefficient=q==3?-sum:w->coefficient[d][q];
                if(f->vertex_chart[v]==(int32_t)chart){
                    if(distance[v-base]==DBL_MAX)goto done;
                    own+=coefficient;bound+=fabs(coefficient)*distance[v-base];
                }else for(int axis=0;axis<2;axis++)target[axis]-=coefficient*(f->uv[2*v+axis]-origin[axis]);
            }
            if(fabs(fabs(own)-1)>1e-8)goto done;
            double weight=.5*w->weight;mass+=weight;radius[at]+=weight*bound/fabs(own);
            for(int axis=0;axis<2;axis++)center[2*at+axis]+=weight*target[axis]/own;
        }
        if(!(mass>0) || !isfinite(mass))goto done;
        radius[at]*=1.25/mass;relation[at]=s->relation;
        for(int axis=0;axis<2;axis++)center[2*at+axis]/=mass;
        if(!isfinite(radius[at]) || !isfinite(center[2*at]) || !isfinite(center[2*at+1]))goto done;
        for(size_t j=0;j<at;j++){
            double separation=hypot(center[2*at]-center[2*j],center[2*at+1]-center[2*j+1]);
            double limit=radius[at]+radius[j]+2*ASM_SEAM_TOLERANCE_VOX/(1.0-1e-8);
            if(isfinite(separation) && isfinite(limit) && separation>limit+1e-6*fmax(1,limit)){
                if(trace)fprintf(trace,"{\"phase\":\"fixed_anchor_infeasible\",\"chart\":%zu,\"relation_a\":%zu,\"relation_b\":%zu,\"target_distance\":%.17g,\"maximum_distance\":%.17g}\n",
                    chart,relation[j],relation[at],separation,limit);
                incompatible=1;goto done;
            }
        }
        at++;
    }
done:
    Arena_dispose(&arena);return incompatible;
}

static int apa_export_endpoint(AsmField *f,const double *before,const char *dir)
{
    if(!dir)return 1;
    ves_mkdir(dir);
    for(size_t c=0;c<f->run->n_charts;c++)if(f->vertex_offset[c]!=SIZE_MAX){
        char path[2048];snprintf(path,sizeof path,"%s/chart_%zu",dir,c);
        if(!amr_trace(f,c,before,path))return 0;
    }
    return 1;
}

#include "asm_repair_failure_cache.inc"

static int apa_attempt_anchor(AsmRun *run,const AsmRepairOpts *opts,const uint8_t *core,
                       size_t *anchors,AsmRepairStats *stats,const char **reason,FILE *trace,
                       const ApaBox *boxes,const AsmField *incumbent,ApaOverlay *overlay,ApaReport *rep,int frontier,size_t requested_anchor)
{
    double t_attempt=ves_clock_sec();
    AsmRepairOpts request=*opts;opts=&request;
    AprWorkspace workspace={0};
    if(opts->compact && !opts->workspace){workspace.measured=apr_process_bytes(&workspace.entry);request.workspace=&workspace;}
    if(rep){memset(rep,0,sizeof *rep);}
    Arena_T arena=Arena_new();size_t nc=run->n_charts,oversize=0;int accepted=0,exported=0,solved=0;*anchors=0;
    AsmField *f=NULL;double *reference=NULL;const uint8_t *moving=core;
    ApaFailureCache *cache=opts->admission_failure_cache;char cache_key[65]={0};int cache_ready=0;
    if(frontier!=4 || opts->admission_rounds!=1 || opts->admission_moving || opts->admission_trace_dir || opts->external_guard ||
       opts->external_context || opts->checkpoint || opts->checkpoint_context || opts->layer_place || opts->layer_conflict || opts->restoration_only || opts->admission_metric_margin || opts->admission_incident_sources)cache=NULL;
    if(apr_expired(opts)){*reason="shared_deadline";goto done;}
    if(opts->admission_moving){
        if(frontier<3 || (overlay && !overlay->changed_charts) || !incumbent){*reason="joint_context_unavailable";goto done;}
        uint8_t *joint=ARENA_ALLOC(arena,nc);memcpy(joint,core,nc);moving=joint;
        for(size_t c=0;c<nc;c++)if(opts->admission_moving[c] && !core[c]){
            if(!AsmChart_registered(run->charts+c)){*reason="joint_chart_not_represented";goto done;}
            joint[c]=1;
        }
    }
    if(opts->compact){size_t faces=0,coordinates=0;
        for(size_t c=0;c<nc;c++)if(moving[c]){faces+=run->charts[c].nf;coordinates+=2*run->charts[c].nv;}
        if(faces>opts->patch_faces || coordinates>opts->patch_coordinates){*reason="moving_resource_cap";goto done;}
    }
    AsmRun trial=*run;trial.arena=arena;
    trial.charts=ARENA_ALLOC(arena,nc*sizeof(AsmChart));memcpy(trial.charts,run->charts,nc*sizeof(AsmChart));
    trial.rels=ARENA_ALLOC(arena,run->n_rels*sizeof(AsmRelation));memcpy(trial.rels,run->rels,run->n_rels*sizeof(AsmRelation));
    uint8_t *fixed=ARENA_CALLOC(arena,nc,1),*known=ARENA_CALLOC(arena,nc,1),*anchor=ARENA_CALLOC(arena,nc,1);
    int8_t *parity=ARENA_ALLOC(arena,nc);memset(parity,-1,nc);
    *reason="no_fixed_source_anchor";
    for(size_t c=0;c<nc;c++){
        fixed[c]=(uint8_t)!moving[c];
        if(AsmChart_registered(run->charts+c)){parity[c]=(int8_t)!!(run->charts[c].flags&ASM_CHART_MIRROR);known[c]=1;}
        if(core[c]){
            AsmChart *ch=trial.charts+c;
            if(!ch->xyz || !ch->faces || !ch->boundary || !ch->nv || !ch->nf ||
               (ch->flags&(ASM_CHART_BLOB|ASM_CHART_NONMANIFOLD|ASM_CHART_TINY))){*reason="unusable_original_component";goto done;}
            for(size_t k=0;k<nc;k++)if(trial.charts[k].parent==(int32_t)c){*reason="component_contains_replaced_parent";goto done;}
            ch->flags&=~(ASM_CHART_EXTRAS|ASM_CHART_FLAT_FAILED);ch->placed=1;ch->placement_state=ASM_PLACE_SEAM;
            ch->placed_uv=NULL;ch->pose_x=ch->pose_y=ch->pose_theta=0;
        }
    }
    if(frontier>=2 && !apa_select_frontier_anchor(&trial,core,requested_anchor)){*reason="no_eligible_source_anchor";goto done;}
    /* Reconcile the observed cavity using eligible deferred source runs.
     * Add constraints only in this private transaction. Existing anchors,
     * observations and vetoes remain intact, and every added run must pass
     * the ordinary parity, source, seam and full material certificate. */
    if(frontier==4 && opts->admission_incident_sources)for(size_t i=0;i<trial.n_rels;i++){
        AsmRelation *r=trial.rels+i;
        if(!apa_source(r) || apa_required_source(r) || core[r->a]==core[r->b] ||
           !apa_present(&trial,r->a,r->b) || (r->flags&(ASM_REL_DROPPED|ASM_REL_CONTACT|ASM_REL_CROSSWRAP)))continue;
        r->continuity|=ASM_CONT_REQUIRED;
        if(trace)fprintf(trace,"{\"phase\":\"proposed_incident_source\",\"relation\":%zu,\"a\":%d,\"b\":%d}\n",i,r->a,r->b);
    }
    for(size_t i=0;i<trial.n_rels;i++)if(apa_source(trial.rels+i) && (core[trial.rels[i].a] || core[trial.rels[i].b])){
        AsmRelation *r=trial.rels+i;
        /* The diagnostic grows through measured source boundaries. Missing
         * endpoints stay unresolved; every incident hypothesis whose two
         * endpoints are now represented is still required and checked. */
        if(frontier && !apa_present(&trial,r->a,r->b))continue;
        if(frontier>=2 && !apa_required_source(r))continue;
        r->continuity|=ASM_CONT_REQUIRED;
        for(int d=0;d<2;d++){size_t c=(size_t)(d?r->b:r->a);if(!core[c]){
            if(!known[c]){*reason="incomplete_missing_source_component";goto done;}anchor[c]=1;
        }}
    }
    for(size_t c=0;c<nc;c++)*anchors+=anchor[c];
    int layer_seed=!*anchors;
    if(layer_seed){
        if(opts->compact){*reason="layer_seed_deferred";goto done;}
        if(!opts->layer_place || !opts->layer_conflict)goto done;
        /* Retain the component's stage-five intrinsic registration. Its
         * global placement is reconsidered using only current observations. */
        for(size_t c=0;c<nc;c++)if(core[c]){
            trial.charts[c]=run->charts[c];
            trial.charts[c].flags&=~(ASM_CHART_EXTRAS|ASM_CHART_FLAT_FAILED);
            trial.charts[c].placed=1;trial.charts[c].placement_state=ASM_PLACE_SEAM;
        }
        if(!apa_intrinsic(&trial,core,&oversize,opts,reason,trace))goto done;
        for(size_t c=0;c<nc;c++)if(core[c])trial.charts[c].placement_state=ASM_PLACE_NONE;
        AsmPlaceStats placement={0};
        int placed=AsmPlace_retry_component(&trial,core,opts->layer_place,opts->layer_conflict,&placement);
        if(placed<=0){*reason=placed<0?"layer_retry_failed":"layer_placement_unsupported";goto done;}
    }
    for(size_t round=0;!layer_seed && round<nc;round++){
        int changed=0;
        for(size_t i=0;i<trial.n_rels;i++){
            const AsmRelation *r=trial.rels+i;if(!apa_source(r) || (!core[r->a] && !core[r->b]))continue;
            if(frontier && !apa_present(&trial,r->a,r->b))continue;
            if(frontier>=2 && !apa_required_source(r))continue;
            int p=!!(r->flags&ASM_REL_PARITY);
            if(parity[r->a]>=0 && parity[r->b]>=0 && (parity[r->a]^parity[r->b])!=p){*reason="contradictory_source_parity";goto done;}
            if(parity[r->a]>=0 && parity[r->b]<0){parity[r->b]=(int8_t)(parity[r->a]^p);changed=1;}
            if(parity[r->b]>=0 && parity[r->a]<0){parity[r->a]=(int8_t)(parity[r->b]^p);changed=1;}
        }
        if(!changed)break;
    }
    for(size_t c=0;!layer_seed && c<nc;c++)if(core[c]){
        if(parity[c]<0){*reason="unanchored_source_parity";goto done;}
        trial.charts[c].flags=(trial.charts[c].flags&~ASM_CHART_MIRROR)|(parity[c]?ASM_CHART_MIRROR:0);
    }
    if(!layer_seed && !apa_intrinsic(&trial,core,&oversize,opts,reason,trace))goto done;
    /* the candidate's own context, or the whole field when nothing places it (a layer seed, or a
     * component with no registered partner to localise around) */
    uint8_t *member=NULL;
    if(ASM_ADMISSION_LOCAL_FIELD && !layer_seed && boxes){
        double reach=0.0;int anchored=0;
        uint8_t *want=ARENA_ALLOC(arena,nc);
        size_t context=apa_context(arena,&trial,core,boxes,want,&reach,&anchored);
        /* The frontier's incumbent is immutable and its exterior is queried
         * in full below. Only actual source partners need solve variables;
         * nearby unrelated charts remain explicit exterior obstacles. */
        if(frontier)for(size_t c=0;c<nc;c++)if(want[c] && !core[c] && !anchor[c]){want[c]=0;context--;}
        if(opts->admission_moving){
            for(size_t c=0;c<nc;c++)if(moving[c] && !want[c]){want[c]=1;context++;}
            /* A freed neighbor brings every represented source partner into
             * the fixed halo, including deferred relations that may pass. */
            for(size_t i=0;i<trial.n_rels;i++){
                const AsmRelation *r=trial.rels+i;
                if(!apa_source(r) || (!moving[r->a] && !moving[r->b]))continue;
                for(int d=0;d<2;d++){size_t c=(size_t)(d?r->b:r->a);
                    if(AsmChart_registered(trial.charts+c) && !want[c]){want[c]=1;context++;}}
            }
        }
        if(anchored){
            member=want;
            if(opts->compact){size_t faces=0;for(size_t c=0;c<nc;c++)if(want[c])faces+=trial.charts[c].nf;
                if(faces>opts->patch_field_faces){*reason="halo_resource_cap";goto done;}}
            size_t estimated=apr_member_workspace(&trial,member);
            if(trace)fprintf(trace,"{\"phase\":\"workspace_preflight\",\"context\":\"source_halo\",\"estimated_bytes\":%zu,\"budget_bytes\":%zu}\n",estimated,opts->patch_bytes);
            if(!apr_workspace_check(opts,estimated)){*reason="workspace_memory_budget";goto done;}
            f=AsmField_new_local(arena,&trial,member,NULL);
            if(rep){rep->context_charts=context;rep->local_field=1;}
            if(!f)member=NULL;
        }
    }
    if(!f && opts->compact){*reason="unbounded_context_deferred";goto done;}
    if(!f)f=AsmField_new(arena,&trial);
    if(!f){*reason="incomplete_source_field";goto done;}
    if(rep){rep->field_vertices=f->nv;rep->field_faces=f->nf;if(!member){rep->local_field=0;rep->context_charts=f->charts;}}
    reference=ARENA_ALLOC(arena,2*f->nv*sizeof(double));memcpy(reference,f->uv,2*f->nv*sizeof(double));
    if(opts->admission_moving)for(size_t c=0;c<nc;c++)if(moving[c] && !core[c] && !apa_chart_metric(f,c)){
        *reason="joint_incumbent_metric_failed";goto done;
    }
    for(size_t i=0;i<f->n_seams;i++)if(moving[f->seams[i].a] || moving[f->seams[i].b]){
        if(frontier && !apa_present(&trial,f->seams[i].a,f->seams[i].b))continue;
        if(f->seams[i].required && (!f->seams[i].complete || !f->seams[i].source_run)){
            *reason="unresolved_incident_source_witnesses";goto done;
        }
    }
    if(!layer_seed && !apa_initialize(f,core,known)){*reason="source_pose_initialization_failed";goto done;}
    /* The tree must be built on the TRANSPORTED candidate.  A local field has none yet, by
     * design; a whole field has one from before the fit, and refitting it would stretch every
     * leaf that mixed candidate and origin material across the move (2026-09-17: 2^34 leaf pairs
     * and a 256 GiB request on the 21x5x5).  A layer seed was placed before the field was built. */
    if(member){if(AsmField_build_contacts(f)){*reason="contact_tree_failed";goto done;}}
    else if(ASM_ADMISSION_REBUILD_AFTER_FIT && !layer_seed && AsmContacts_rebuild(f->contacts,f->uv)){*reason="contact_rebuild_failed";goto done;}
    apa_trace(f,core,layer_seed?"initial_layer_fit":"initial_source_fit",trace);
    if(frontier && opts->iterations && !opts->admission_moving && apa_incompatible_fixed_anchors(f,core,trace)){
        *reason="incompatible_fixed_anchors";goto done;
    }
    if(frontier>=3 && member){
        f=apa_obstacle_field(f,moving,member,incumbent,overlay,opts,reason,trace);
        if(!f)goto done;
        reference=ARENA_ALLOC(arena,2*f->nv*sizeof(double));memcpy(reference,f->uv,2*f->nv*sizeof(double));
        if(rep){rep->context_charts=f->charts;rep->field_vertices=f->nv;rep->field_faces=f->nf;}
    }
    if(opts->admission_moving){
        for(size_t i=0;i<f->n_seams;i++){
            AsmFieldSeam *s=f->seams+i;
            if(core[s->a] || core[s->b] || (!moving[s->a] && !moving[s->b]))continue;
            AsmSeamMeasure old;AsmField_seam(f,i,f->uv,&old);
            if(old.pass)s->required=1; /* Scratch requirement only; retain the source ledger. */
        }
        apa_trace(f,moving,"joint_initial",trace);
    }
    /* Keep reference and the incumbent's passing-seam requirements from
     * before the continuation. Scratch progress cannot redefine either. */
    if(apa_continuation_apply(opts->admission_continuation,f,moving)){
        if(AsmContacts_rebuild(f->contacts,f->uv)){*reason="continuation_contact_tree_failed";goto done;}
        apa_trace(f,moving,"resumed_candidate",trace);
    }
    AsmRepairOpts local=*opts;local.admission=local.regional=0;local.fixed_charts=fixed;
    ApaSolverContacts contact_merit={0};
    if(frontier==4){
        contact_merit.certificate.arena=arena;contact_merit.certificate.field=f;
        contact_merit.certificate.core=core;contact_merit.certificate.reference=reference;
        contact_merit.certificate.okay=1;contact_merit.trace=trace;
        local.contact_merit=apa_solver_contact_merit;local.contact_merit_context=&contact_merit;
        local.contact_objective=apa_solver_contact_objective;local.contact_copy_active=apa_solver_copy_active;
    }
    uint8_t *active=ARENA_CALLOC(arena,f->nf,1);size_t candidate_faces=0;double candidate_area=0;
    for(size_t c=0;c<nc;c++)if(moving[c]){
        memset(active+f->face_offset[c],1,trial.charts[c].nf);
        if(core[c]){
            candidate_faces+=trial.charts[c].nf;
            for(size_t t=0;t<trial.charts[c].nf;t++)candidate_area+=f->metric[f->face_offset[c]+t].area;
        }
    }
    {
        /* Bound hopeless starts against fixed material. All moving contacts
         * still enter the solver and final certificate; this shortcut must
         * not preclude the neighboring poses that joint repair can change. */
        AsmContactStats initial={0};AsmContacts_active_faces(f->contacts,active);
        ApaInitialContacts start={f,core,moving,0,0,0,0};
        int measured=AsmContacts_measure(f->contacts,f->uv,NULL,apa_initial_contact,&start,&initial);
        AsmContacts_active_faces(f->contacts,NULL);
        if(measured==-2){
            fprintf(stderr,"[assemble admission] candidate of %zu faces: the contact walk exceeded its budget (%zu leaf pairs); refused before the solve\n",candidate_faces,initial.leaf_pairs);
            *reason="contact_walk_over_budget";goto done;
        }
        if(measured || !initial.complete){*reason="contact_query_incomplete";goto done;}
        size_t pairs=initial.pairs;double contact_area=initial.area;
        if(!measured && member){
            /* the placed material this local field leaves out counts too */
            size_t ep=0;double ea=0.0;int ec=1;
            if(apa_external_pairs(f,moving,member,incumbent,overlay,arena,&ep,&ea,&ec,frontier)){*reason="external_probe_failed";goto done;}
            if(!ec){*reason="contact_walk_over_budget";goto done;}
            pairs+=ep;contact_area+=ea;
            if(opts->admission_moving){
                if(apa_external_pairs(f,core,member,incumbent,overlay,arena,&ep,&ea,&ec,frontier)){*reason="external_probe_failed";goto done;}
                if(!ec){*reason="contact_walk_over_budget";goto done;}
            }
            start.blocking_pairs+=ep;start.blocking_area+=ea;
            start.far_pairs+=ep;start.far_area+=ea;   /* no geometry out there: far */
        }
        if(trace && !measured){
            double source_area=0;
            for(size_t t=0;t<f->nf;t++)if(active[t])source_area+=f->metric[t].area;
            fprintf(trace,"{\"phase\":\"initial_contacts\",\"pairs\":%zu,\"area\":%.17g,\"source_area\":%.17g,\"faces\":%zu,\"fixed_obstacle_pairs\":%zu,\"fixed_obstacle_area\":%.17g,"
                    "\"candidate_area\":%.17g,\"far_obstacle_pairs\":%zu,\"far_obstacle_area\":%.17g,\"complete\":true}\n",
                    pairs,contact_area,source_area,candidate_faces,start.blocking_pairs,start.blocking_area,candidate_area,start.far_pairs,start.far_area);
        }
        if(!measured && apa_start_hopeless(&start,candidate_area,candidate_faces)){
            fprintf(stderr,"[assemble admission] candidate of %zu faces (%.1f vox^2) starts on %zu contact pairs (%.1f vox^2, %.1f vox^2 of it farther than %.1f vox in 3-D) with fixed material; refused before the solve\n",
                    candidate_faces,candidate_area,start.blocking_pairs,start.blocking_area,start.far_area,(double)ASM_DISCOVER_SAME_SURFACE_VOX);
            *reason="candidate_overlaps_placed_material";goto done;
        }
        if(cache && !measured && initial.complete && !apr_expired(opts)){
            double began=ves_clock_sec();afc_key(f,&local,moving,reference,frontier,cache_key);
            cache->hash_seconds+=ves_clock_sec()-began;cache_ready=1;
            const AfcEntry *previous=afc_lookup(cache,cache_key);
            if(trace)fprintf(trace,"{\"phase\":\"local_failure_key\",\"key\":\"%s\",\"hit\":%s}\n",cache_key,previous?"true":"false");
            if(previous){
                if(rep)rep->cache_hit=1;*reason=previous->reason;
                if(trace)fprintf(trace,"{\"phase\":\"cached_local_failure\",\"key\":\"%s\",\"reason\":\"%s\",\"previous_work_seconds\":%.17g}\n",cache_key,previous->reason,previous->seconds);
                goto done;
            }
        }
    }
    if(opts->admission_feedback && opts->admission_feedback->observe_only){*reason="contact_context_survey";goto done;}
    double *round_start=ARENA_ALLOC(arena,2*f->nv*sizeof(double));
    for(int round=0;round<opts->admission_rounds;round++){
        double solve_started=ves_clock_sec();int complete=0,polish_complete=1;
        /* A rigid source fit can already satisfy the complete certificate.
         * Check every starting point before spending a solve quantum. Fresh
         * candidates receive a read-only check; saved candidates retain
         * their bounded polishing route. Exterior checks remain mandatory. */
        int have_warm=opts->admission_continuation && opts->admission_continuation->arena;
        const char *warm_refused=have_warm?
            apa_qualified_endpoint(f,moving,core,reference,active,&local,trace,frontier,0,&complete,&polish_complete):
            apa_qualification(f,moving,core,reference,active,trace,frontier,&complete);
        if(!warm_refused && !have_warm && local.admission_metric_margin>0)for(size_t c=0;c<nc;c++)
            if(moving[c] && !apa_chart_metric_margin(f,c,local.admission_metric_margin))warm_refused="packing_metric_margin";
        int initial_qualified=!warm_refused,warm_qualified=have_warm && initial_qualified;
        /* Bound hopeless starts on seams, as the contact ratio does against
         * fixed material: with only the candidate moving, a required seam
         * that starts ASM_ADMISSION_MAX_START_SEAM_RMS off is not closed by
         * the solve within the metric and contact bounds.  The refusal is the
         * one the solve would have returned. */
        if(round==0 && !have_warm && !opts->admission_moving && warm_refused && !strcmp(warm_refused,"incident_source_seam_failed") &&
           ASM_ADMISSION_MAX_START_SEAM_RMS>0){
            double worst=apa_worst_start_seam(f,moving,frontier);
            if(worst>ASM_ADMISSION_MAX_START_SEAM_RMS){
                if(trace)fprintf(trace,"{\"phase\":\"hopeless_seam_start\",\"worst_rms\":%.17g,\"bound\":%.17g}\n",
                                 isfinite(worst)?worst:-1.0,(double)ASM_ADMISSION_MAX_START_SEAM_RMS);
                *reason="incident_source_seam_failed";goto done;
            }
        }
        local.clearance_first=frontier==4 && have_warm && complete && warm_refused &&
            (!strcmp(warm_refused,"new_material_contact_failed") || !strcmp(warm_refused,"incumbent_contact_regression"));
        memcpy(round_start,f->uv,2*f->nv*sizeof(double));
        AsmRepairStats part={0};part.first_metric_chart=-1;
        /* A pure rigid-fit diagnostic has no requested solver iterations.
         * Its complete admission certificate is still evaluated below.
         * The residual nudge (ASM_REPAIR_RESIDUAL_NUDGE) runs inside the same
         * quantum: the solve stops its allowance early when the quantum is
         * longer than twice that allowance. */
        double reserve=apa_nudge_frontier(frontier) && local.deadline>0 &&
            local.deadline-ves_clock_sec()>2*ASM_REPAIR_RESIDUAL_NUDGE_SEC?ASM_REPAIR_RESIDUAL_NUDGE_SEC:0;
        double request_reserve=request.deadline>0?reserve:0;
        local.deadline-=reserve;request.deadline-=request_reserve;
        int solve_rc=initial_qualified || (frontier && !local.iterations && !local.clearance_iterations && !local.closure_iterations &&
            !local.restoration_only)?0:apa_solve_owned(f,&local,&request,&part,rep);
        /* Include refused private work too: scheduling and diagnostics must
         * see the pose/deformation steps that actually consumed the quantum. */
        apr_accumulate(stats,&part);
        local.deadline+=reserve;request.deadline+=request_reserve;
        if(trace)fprintf(trace,"{\"phase\":\"solver_work\",\"warm_qualified\":%s,\"initial_qualified\":%s,\"iterations\":%zu,"
            "\"clearance_steps\":%zu,\"closure_steps\":%zu,\"metric_refused\":%zu,\"chart_metric_refused\":%zu,"
            "\"seam_refused\":%zu,\"contact_refused\":%zu,\"energy_refused\":%zu,\"deadline_reached\":%s}\n",
            warm_qualified?"true":"false",initial_qualified?"true":"false",part.iterations,part.clearance_steps,part.closure_steps,part.rejected_metric,
            part.rejected_chart,part.rejected_seam,part.rejected_contacts,part.rejected_energy,part.deadline_reached?"true":"false");
        if(trace)fprintf(trace,"{\"phase\":\"solver_cost\",\"active_coordinates\":%zu,\"general_seconds\":%.17g,\"clearance_seconds\":%.17g,\"factor_seconds\":%.17g}\n",
            part.active_coordinates,part.general_seconds,part.clearance_seconds,part.factor_seconds);
        if(trace && local.pose_iterations)fprintf(trace,"{\"phase\":\"joint_pose_work\",\"steps\":%zu,\"seconds\":%.17g}\n",part.pose_phase_steps,part.pose_phase_seconds);
        if(trace)fprintf(trace,"{\"phase\":\"phase_order\",\"clearance_first\":%s,\"pre_fit_clearance_steps\":%zu,\"pre_fit_clearance_seconds\":%.17g,\"fitting_iterations\":%zu}\n",
            local.clearance_first?"true":"false",part.pre_fit_clearance_steps,part.pre_fit_clearance_seconds,part.iterations);
        if(trace)fprintf(trace,"{\"phase\":\"checkpoint_merit\",\"policy\":\"physical_blocking_contact_area_v1\",\"evaluations\":%zu,\"seconds\":%.17g}\n",
            part.contact_merit_evaluations,part.contact_merit_seconds);
        if(solve_rc){*reason=solve_rc==1?"contact_walk_over_budget":"candidate_solve_failed";goto done;}
        if(!apr_workspace_check(opts,0)){*reason="workspace_memory_budget";goto done;}
        solved=local.iterations>0 || local.clearance_iterations>0;
        apa_trace(f,moving,"after_repair",trace);
        const char *refused=initial_qualified?NULL:
            apa_qualified_endpoint(f,moving,core,reference,active,&local,trace,frontier,1,&complete,&polish_complete);
        if(refused && cache_ready && polish_complete && afc_reason(refused,complete,part.deadline_reached || apr_expired(opts))){
            size_t before=cache->count;
            afc_insert(cache,cache_key,refused,ves_clock_sec()-solve_started);
            if(trace && cache->count>before)fprintf(trace,"{\"phase\":\"retained_local_failure\",\"key\":\"%s\",\"reason\":\"%s\"}\n",cache_key,refused);
        }
        if(!refused && member){
            /* the local field proved its own material clear; the placed sheet it leaves out is
             * proved by probing the candidate's faces against the incumbent and this pass's
             * earlier admissions */
            size_t ep=0;double ea=0.0;int ec=1;
            if(apa_external_pairs(f,moving,member,incumbent,overlay,arena,&ep,&ea,&ec,frontier))refused="external_probe_failed";
            else if(!ec)refused="contact_walk_over_budget";
            else if(ep)refused="new_material_contact_failed";
            if(rep)rep->external_pairs=ep;
            if(trace)fprintf(trace,"{\"phase\":\"external_contacts\",\"pairs\":%zu,\"area\":%.17g,\"complete\":%s}\n",ep,ea,ec?"true":"false");
        }
        if(!refused)break;
        *reason=refused;
        if(round+1==opts->admission_rounds || !memcmp(round_start,f->uv,2*f->nv*sizeof(double)))goto done;
        fprintf(stderr,"[assemble admission] continuing candidate, round %d/%d: %s\n",round+2,opts->admission_rounds,refused);
    }
    if(opts->admission_trace_dir){
        exported=1;if(!apa_export_endpoint(f,reference,opts->admission_trace_dir)){
            if(rep)rep->endpoint_export_failed=1;*reason="endpoint_export_failed";goto done;
        }
    }
    /* Export and prepare the exterior copy before committing membership.
     * After publication no fallible operation may split the run and tree. */
    AsmRun staged;
    if(!api_cut_stage(run,&trial,core,&staged)){*reason="cut_commit_source_mismatch";goto done;}
    if(overlay){
        for(size_t c=0;c<nc;c++)if(moving[c])trial.charts[c].placed_uv=f->uv+2*f->vertex_offset[c];
        if(apa_overlay_add(overlay,&trial,moving)){*reason="overlay_failed";goto done;}
    }
    *run=staged;
    for(size_t c=0;c<nc;c++)if(core[c]){
        AsmChart *ch=run->charts+c;size_t bytes=2*ch->nv*sizeof(double);
        double *uv=ARENA_ALLOC(run->arena,bytes);memcpy(uv,f->uv+2*f->vertex_offset[c],bytes);
        /* Cold flattening owns its seed in the trial arena. Retain that
         * seed in the run before releasing the transaction, including the
         * case where the original chart never had intrinsic coordinates. */
        if(ch->uv!=trial.charts[c].uv){
            ch->uv=ARENA_ALLOC(run->arena,2*ch->nv*sizeof(float));
            memcpy(ch->uv,trial.charts[c].uv,2*ch->nv*sizeof(float));
        }
        ch->placed_uv=uv;ch->flags=trial.charts[c].flags;ch->placed=1;
        ch->placement_state=layer_seed?trial.charts[c].placement_state:ASM_PLACE_SEAM;
        ch->placement_parent=layer_seed?trial.charts[c].placement_parent:-1;
        ch->placement_support=layer_seed?trial.charts[c].placement_support:0;
        ch->placement_rms=layer_seed?trial.charts[c].placement_rms:0;
        for(size_t i=0;i<f->n_seams;i++)if(f->seams[i].required && (f->seams[i].a==(int32_t)c || f->seams[i].b==(int32_t)c)){
            if(frontier && !apa_present(&trial,f->seams[i].a,f->seams[i].b))continue;
            AsmSeamMeasure m;AsmField_seam(f,i,f->uv,&m);ch->placement_support+=m.observations;ch->placement_rms=fmax(ch->placement_rms,m.rms);
            int other=f->seams[i].a==(int32_t)c?f->seams[i].b:f->seams[i].a;
            if(ch->placement_parent<0 && !core[other])ch->placement_parent=other;
        }
        stats->admitted_charts++;
    }
    if(opts->admission_moving)for(size_t c=0;c<nc;c++)if(moving[c] && !core[c]){
        size_t bytes=2*run->charts[c].nv*sizeof(double);double *uv=ARENA_ALLOC(run->arena,bytes);
        memcpy(uv,f->uv+2*f->vertex_offset[c],bytes);run->charts[c].placed_uv=uv;
    }
    for(size_t i=0;i<run->n_rels;i++)if(apa_source(run->rels+i) && (core[run->rels[i].a] || core[run->rels[i].b]) &&
        (!frontier || apa_present(&trial,run->rels[i].a,run->rels[i].b)) &&
        (frontier<2 || apa_required_source(trial.rels+i)))run->rels[i].continuity|=ASM_CONT_REQUIRED;
    stats->admissions++;*reason=layer_seed?"accepted_complete_original_layer_component":"accepted_complete_original_component";accepted=1;
done:
    /* UV-only continuation records cannot represent an uncommitted change
     * of topology. Reconstruct that private cut on the next attempt; never
     * serialize an alias-sized UV array under the old source vertex count. */
    if(!accepted && solved && opts->admission_continuation && f && reference && f->run->n_rels==run->n_rels)
        apa_continuation_save(opts->admission_continuation,f,moving,reference);
    if(!accepted && opts->admission_feedback)
        apa_feedback(f,core,moving,reference,incumbent,overlay,arena,opts->admission_feedback);
    if(!exported && f && reference && opts->admission_trace_dir && !apa_export_endpoint(f,reference,opts->admission_trace_dir)){
        if(rep)rep->endpoint_export_failed=1;
        fprintf(stderr,"[joint admission] endpoint export failed after geometry refusal: %s\n",*reason);
    }
    if(rep)rep->sec=ves_clock_sec()-t_attempt;
    apr_workspace_check(opts,0);
    if(trace && opts->workspace)fprintf(trace,"{\"phase\":\"workspace_usage\",\"estimated_bytes\":%zu,\"observed_process_growth_bytes\":%zu,\"budget_bytes\":%zu,\"measured\":%s,\"exceeded\":%s}\n",
        opts->workspace->estimate,opts->workspace->peak,opts->patch_bytes,opts->workspace->measured?"true":"false",opts->workspace->exceeded?"true":"false");
    Arena_dispose(&arena);return accepted;
}

static int apa_attempt(AsmRun *run,const AsmRepairOpts *opts,const uint8_t *core,
                       size_t *anchors,AsmRepairStats *stats,const char **reason,FILE *trace,
                       const ApaBox *boxes,const AsmField *incumbent,ApaOverlay *overlay,ApaReport *rep,int frontier)
{return apa_attempt_anchor(run,opts,core,anchors,stats,reason,trace,boxes,incumbent,overlay,rep,frontier,SIZE_MAX);}

static int apr_save_progress(AsmRun *run,const char *dir);
static int apr_admit_to(AsmRun *run,const AsmRepairOpts *opts,const char *dir,AsmRepairStats *stats,const char *ledger,
                        const AsmField *incumbent)
{
    if(!opts->admission || opts->fixed_charts || apr_expired(opts))return 0;
    Arena_T arena=Arena_new();size_t nc=run->n_charts;
    size_t *parent=ARENA_ALLOC(arena,nc*sizeof(size_t));uint8_t *core=ARENA_ALLOC(arena,nc);
    for(size_t c=0;c<nc;c++)parent[c]=c;
    for(size_t i=0;i<run->n_rels;i++){
        const AsmRelation *r=run->rels+i;
        if(r->a<0 || r->b<0 || (size_t)r->a>=nc || (size_t)r->b>=nc){Arena_dispose(&arena);return -1;}
        if(apa_source(r) && !AsmChart_registered(run->charts+r->a) && !AsmChart_registered(run->charts+r->b)){
            size_t a=apr_root(parent,r->a),b=apr_root(parent,r->b);parent[a>b?a:b]=a<b?a:b;
        }
    }
    FILE *log=NULL;int io=1;
    if(dir){char path[2048];snprintf(path,sizeof path,"%s/%s",dir,ledger);log=fopen(path,"wb");if(!log){Arena_dispose(&arena);return -1;}
        io=fprintf(log,"component,original_charts,fixed_source_anchors,accepted,reason,local_field,context_charts,field_vertices,field_faces,external_pairs,seconds\n")>=0;}
    /* candidates in the order of their evidence: fixed source anchors, then area */
    typedef struct { size_t root, count, anchors; double area; } ApaCandidate;
    ApaCandidate *cands=ARENA_ALLOC(arena,(nc?nc:1)*sizeof(ApaCandidate));size_t ncand=0;
    for(size_t root=0;root<nc;root++)if(!AsmChart_registered(run->charts+root) && apr_root(parent,root)==root){
        ApaCandidate cd={root,0,0,0.0};
        for(size_t c=0;c<nc;c++)if(apr_root(parent,c)==root){cd.count++;cd.area+=run->charts[c].area3d;}
        for(size_t i=0;i<run->n_rels;i++){const AsmRelation *r=run->rels+i;if(!apa_source(r))continue;
            int ia=apr_root(parent,(size_t)r->a)==root,ib=apr_root(parent,(size_t)r->b)==root;
            if(ia!=ib && (AsmChart_registered(run->charts+(ia?r->b:r->a))))cd.anchors++;}
        cands[ncand++]=cd;
    }
    for(size_t i=1;i<ncand;i++){ApaCandidate k=cands[i];size_t j=i;
        while(j>0 && (cands[j-1].anchors<k.anchors || (cands[j-1].anchors==k.anchors && cands[j-1].area<k.area))){cands[j]=cands[j-1];j--;}
        cands[j]=k;}
    /* every registered chart's uv box, for the candidates' contexts, and the overlay this pass
     * fills with what it admits */
    ApaBox *boxes=NULL;
    ApaOverlay overlay;memset(&overlay,0,sizeof overlay);overlay.arena=arena;
    if(ASM_ADMISSION_LOCAL_FIELD && incumbent){
        double t_box=ves_clock_sec();
        boxes=ARENA_ALLOC(arena,(nc?nc:1)*sizeof(ApaBox));
        apa_chart_boxes(run,boxes);
        fprintf(stderr,"[assemble admission] local trial fields: %zu candidates against an incumbent of %zu charts / %zu faces; chart boxes in %.1f s\n",
                ncand,incumbent->charts,incumbent->nf,ves_clock_sec()-t_box);
    }
    double pass_start=ves_clock_sec();size_t skipped=0,local=0,whole=0;double solve_sec=0.0;
    for(size_t ci=0;ci<ncand;ci++){
        size_t root=cands[ci].root,count=cands[ci].count,anchors=0;
        for(size_t c=0;c<nc;c++)core[c]=(uint8_t)(apr_root(parent,c)==root);
        const char *reason;int accepted=0;ApaReport rep;memset(&rep,0,sizeof rep);
        if(apr_expired(opts) || ves_clock_sec()-pass_start>ASM_ADMISSION_BUDGET_SEC){reason="admission_budget_exhausted";anchors=cands[ci].anchors;skipped++;}
        else{
            AsmRepairOpts bounded=*opts;
            if(opts->compact){double stop=ves_clock_sec()+opts->patch_seconds;if(bounded.deadline<=0 || stop<bounded.deadline)bounded.deadline=stop;}
            accepted=apa_attempt(run,&bounded,core,&anchors,stats,&reason,NULL,boxes,incumbent,boxes?&overlay:NULL,&rep,0);
            if(rep.local_field)local++;else whole++;
            solve_sec+=rep.sec;
        }
        if(!skipped || accepted)fprintf(stderr,"[assemble admission] component %zu (%zu charts, %zu anchors): %s%s\n",root,count,anchors,reason,
                                        rep.local_field?" [local field]":"");
        if(log && fprintf(log,"%zu,%zu,%zu,%d,%s,%d,%zu,%zu,%zu,%zu,%.3f\n",root,count,anchors,accepted,reason,
                          rep.local_field,rep.context_charts,rep.field_vertices,rep.field_faces,rep.external_pairs,rep.sec)<0)io=0;
        if(log && fflush(log))io=0;
        if(accepted && dir && !opts->compact && apr_save_progress(run,dir)){io=0;break;}
    }
    fprintf(stderr,"[assemble admission] %zu attempted (%zu on a local field, %zu on the whole field) in %.1f s, %.2f s each; %zu admitted\n",
            local+whole,local,whole,solve_sec,local+whole?solve_sec/(double)(local+whole):0.0,stats->admissions);
    if(skipped)fprintf(stderr,"[assemble admission] budget %.0f s reached: %zu of %zu candidates not attempted (ledger %s)\n",(double)ASM_ADMISSION_BUDGET_SEC,skipped,ncand,ledger);
    if(log && fclose(log))io=0;Arena_dispose(&arena);return io?0:-1;
}
static int apr_admit(AsmRun *run,const AsmRepairOpts *opts,const char *dir,AsmRepairStats *stats,const AsmField *incumbent)
{return apr_admit_to(run,opts,dir,stats,"stage6_admission.csv",incumbent);}

static void apa_require(AsmField *f,size_t seam,int required)
{
    f->seams[seam].required=required;
    for(size_t k=0;k<f->seams[seam].count;k++)f->witness[f->seams[seam].first+k].required=required;
}

/* Reconstructed/deferred hypotheses become obligations only after a complete
 * transaction. Existing placement veto flags survive; this is a geometric
 * repair hypothesis, not a newly asserted physical material identity. */
static int apr_reconnect(AsmField *f,AsmRun *run,const AsmRepairOpts *opts,const char *dir,AsmRepairStats *stats)
{
    if(!opts->source_trials || opts->fixed_charts)return 0;
    Arena_T arena=Arena_new();size_t nc=run->n_charts;
    size_t *parent=ARENA_ALLOC(arena,nc*sizeof(size_t));double *area=ARENA_ALLOC(arena,nc*sizeof(double));
    uint8_t *fixed=ARENA_ALLOC(arena,nc),*core=ARENA_ALLOC(arena,nc),*passing=ARENA_ALLOC(arena,f->n_seams);
    double *saved=ARENA_ALLOC(arena,2*f->nv*sizeof(double));ApcOrder *order=ARENA_ALLOC(arena,f->n_seams*sizeof(ApcOrder));size_t count=0;
    FILE *log=NULL;int io=1;
    if(dir){char path[2048];snprintf(path,sizeof path,"%s/stage6_connections.csv",dir);log=fopen(path,"wb");if(!log){Arena_dispose(&arena);return -1;}
        io=fprintf(log,"relation,chart_a,chart_b,initial_rms,final_rms,accepted,reason\n")>=0;if(fflush(log))io=0;}
    for(size_t i=0;i<f->n_seams;i++){
        const AsmFieldSeam *seam=f->seams+i;
        if(seam->required || !seam->complete || !seam->parity || !apa_source(run->rels+seam->relation))continue;
        AsmSeamMeasure m;AsmField_seam(f,i,f->uv,&m);
        if(m.observations>=6 && m.mass>=8 && !m.invalid)order[count++]=(ApcOrder){i,m.rms};
    }
    qsort(order,count,sizeof *order,apc_order);
    uint8_t *done=ARENA_CALLOC(arena,count?count:1,1);
#if ASM_SOURCE_TRIALS_BATCHED
    /* One restoration per large component: every queued failing seam whose
     * causal core is that component is required at once, the ladder runs
     * once, and each seam is judged on its own afterwards.  The invariants of
     * the sequential path hold for the batch as a whole -- no previously
     * passing seam regresses, no fixed context moves -- or the batch is
     * reverted and its seams fall through to the sequential loop below.
     * Accepted seams become required and merge components, so the pass
     * repeats until no batch of two forms; a seam that fails twice in a
     * batch is left as it is (2026-09-16: 15 of 21 primary seams closed in
     * one ladder where the sequential path ran one ladder per seam). */
    {
        size_t *root_of=ARENA_ALLOC(arena,(count?count:1)*sizeof(size_t)),*core_size=ARENA_ALLOC(arena,(count?count:1)*sizeof(size_t));
        size_t *batch=ARENA_ALLOC(arena,(count?count:1)*sizeof(size_t));double *rms0=ARENA_ALLOC(arena,(count?count:1)*sizeof(double));
        uint8_t *core_k=ARENA_ALLOC(arena,nc),*fail=ARENA_CALLOC(arena,count?count:1,1),*attempts=ARENA_CALLOC(arena,count?count:1,1);
        for(;;){
            if(stats->source_trials>=(size_t)opts->source_trials || apr_expired(opts))break;
            for(size_t c=0;c<nc;c++)parent[c]=c;
            for(size_t i=0;i<f->n_seams;i++)if(f->seams[i].required && f->seams[i].complete){
                size_t a=apr_root(parent,f->seams[i].a),b=apr_root(parent,f->seams[i].b);parent[a>b?a:b]=a<b?a:b;
            }
            for(size_t k=0;k<count;k++){
                fail[k]=0;core_size[k]=0;root_of[k]=SIZE_MAX;
                if(done[k] || attempts[k]>=ASM_SOURCE_TRIALS_BATCH_ROUNDS)continue;
                AsmSeamMeasure m;AsmField_seam(f,order[k].id,f->uv,&m);
                if(m.pass || !(m.rms<=ASM_SOURCE_TRIAL_MAX_RMS))continue;
                fail[k]=1;rms0[k]=m.rms;const AsmFieldSeam *s=f->seams+order[k].id;
                core_size[k]=apr_causal_core(f,(size_t)s->a,(size_t)s->b,core_k,opts->context_hops);root_of[k]=apr_root(parent,(size_t)s->a);
            }
            /* The batch takes every failing seam of one component, and nothing bounded what
             * that costs: the moving set is the UNION of their causal cores, so the healthier
             * the sheet the bigger the batch.  The 21x5x5's components are rubble and its
             * batches run 70-150 charts; the 10x10x10's primary IS the sheet (63% of it in one
             * component) and its first batch came to 1,975 moving charts, which the solve then
             * sat on for over an hour (2026-09-19).  Better placement must not mean a slower
             * repair.  Cap the moving charts and let the rest form the next batch; a single
             * seam whose own core is over the cap still goes alone, so nothing is unreachable. */
            size_t members=0,root=SIZE_MAX,moving=0;memset(core,0,nc);
            for(size_t k=0;k<count;k++){
                if(!fail[k] || core_size[k]<ASM_SOURCE_TRIALS_BATCH_MIN_CHARTS)continue;
                if(root==SIZE_MAX)root=root_of[k];
                if(root_of[k]!=root)continue;
                const AsmFieldSeam *s=f->seams+order[k].id;apr_causal_core(f,(size_t)s->a,(size_t)s->b,core_k,opts->context_hops);
                if(!apa_batch_take(core,core_k,nc,f->vertex_offset,&moving,members))break;
                batch[members++]=k;
            }
            if(members<1)break;   /* a batch of one keeps the pass alive: its acceptance merges components and makes the next batch */
            stats->source_trials++;
            for(size_t c=0;c<nc;c++)fixed[c]=(uint8_t)!core[c];
            memcpy(saved,f->uv,2*f->nv*sizeof(double));
            for(size_t i=0;i<f->n_seams;i++){AsmSeamMeasure m;AsmField_seam(f,i,f->uv,&m);passing[i]=(uint8_t)(f->seams[i].required && m.pass);}
            for(size_t t=0;t<members;t++)apa_require(f,order[batch[t]].id,1);
            AsmRepairOpts local=*opts;local.fixed_charts=fixed;local.admission=local.regional=0;
            size_t moving_charts=0;for(size_t c=0;c<nc;c++)moving_charts+=!fixed[c] && f->vertex_offset[c]!=SIZE_MAX;
            fprintf(stderr,"[assemble source repair] batch trial %zu: %zu seams of component %zu, %zu moving charts\n",stats->source_trials,members,root,moving_charts);
            AsmRepairStats part={0};part.first_metric_chart=-1;
            int rc=apr_solve(f,&local,&part);int keep=!rc;
            if(rc==1)stats->trials_refused_budget++;
            for(size_t i=0;i<f->n_seams && keep;i++)if(passing[i]){AsmSeamMeasure m;AsmField_seam(f,i,f->uv,&m);if(!m.pass)keep=0;}
            for(size_t c=0;c<nc && keep;c++)if(fixed[c] && f->vertex_offset[c]!=SIZE_MAX &&
                memcmp(saved+2*f->vertex_offset[c],f->uv+2*f->vertex_offset[c],2*run->charts[c].nv*sizeof(double)))keep=0;
            if(!keep){
                memcpy(f->uv,saved,2*f->nv*sizeof(double));
                for(size_t t=0;t<members;t++){apa_require(f,order[batch[t]].id,0);attempts[batch[t]]=ASM_SOURCE_TRIALS_BATCH_ROUNDS;}
                fprintf(stderr,"[assemble source repair] batch trial %zu: %s; %zu seams take the sequential path\n",stats->source_trials,rc==1?"contact_walk_over_budget":rc?"candidate_solve_failed":"existing_seam_regression_or_fixed_context_changed",members);
                continue;
            }
            size_t accepted=0;
            /* A kept batch that left the field exactly as it found it will decide
             * exactly the same way next round: same charts, same required seams,
             * same deterministic solve.  Deferring those seams buys a byte-identical
             * repeat -- the 21x5x5's component 22 reran 73 moving charts through the
             * same refused initialization (19.33 / 288.97 / 1006.01 / 1001.26 vox)
             * for nothing, 2026-09-17.  A batch that DID move the field may decide
             * differently and keeps its rounds. */
            int field_unchanged=!memcmp(saved,f->uv,2*f->nv*sizeof(double));
            for(size_t t=0;t<members;t++){
                size_t k=batch[t],sj=order[k].id;const AsmFieldSeam *s=f->seams+sj;AsmSeamMeasure m;AsmField_seam(f,sj,f->uv,&m);
                const char *reason;
                if(m.pass){run->rels[s->relation].continuity|=ASM_CONT_REQUIRED;stats->recovered_connections++;accepted++;done[k]=1;reason="accepted_in_batch";}
                else if(field_unchanged){apa_require(f,sj,0);attempts[k]=ASM_SOURCE_TRIALS_BATCH_ROUNDS;done[k]=1;reason="failed_in_batch_field_unchanged";}
                else{apa_require(f,sj,0);attempts[k]++;if(attempts[k]>=ASM_SOURCE_TRIALS_BATCH_ROUNDS){done[k]=1;reason="failed_in_batch";}else reason="deferred_to_next_batch";}
                fprintf(stderr,"[assemble source repair] relation %zu (%d,%d): %s\n",s->relation,s->a,s->b,reason);
                if(done[k] && log && fprintf(log,"%zu,%d,%d,%.17g,%.17g,%d,%s\n",s->relation,s->a,s->b,rms0[k],m.rms,m.pass?1:0,reason)<0)io=0;
            }
            if(log && fflush(log))io=0;
            fprintf(stderr,"[assemble source repair] batch trial %zu: %zu of %zu seams accepted; field kept\n",stats->source_trials,accepted,members);
            if(accepted && opts->checkpoint && opts->checkpoint(f,opts->checkpoint_context)){if(log)fclose(log);Arena_dispose(&arena);return -1;}
        }
    }
#endif
    for(size_t k=0;k<count && !apr_expired(opts);k++){
        if(done[k])continue;
        size_t si=order[k].id;AsmFieldSeam *seam=f->seams+si;AsmSeamMeasure before,after;AsmField_seam(f,si,f->uv,&before);
        const char *reason="already_satisfies_complete_source_geometry";int accept=before.pass;
        if(!accept && !(before.rms<=ASM_SOURCE_TRIAL_MAX_RMS)){reason="beyond_source_repair_reach";after=before;}
        else if(!accept && stats->source_trials>=(size_t)opts->source_trials){reason="bounded_trial_budget";after=before;}
        else if(!accept){
            stats->source_trials++;
            for(size_t c=0;c<nc;c++){parent[c]=c;area[c]=0;}
            for(size_t i=0;i<f->n_seams;i++)if(f->seams[i].required && f->seams[i].complete){
                size_t a=apr_root(parent,f->seams[i].a),b=apr_root(parent,f->seams[i].b);parent[a>b?a:b]=a<b?a:b;
            }
            for(size_t t=0;t<f->nf;t++)area[apr_root(parent,f->face_chart[t])]+=f->metric[t].area;
            size_t a=apr_root(parent,seam->a),b=apr_root(parent,seam->b),moving=area[a]<area[b]?a:b;
            /* Start with the causal source path, retaining both sides of
             * every traversed cycle. Full context remains in the field. */
            apr_causal_core(f,(size_t)seam->a,(size_t)seam->b,core,opts->context_hops);
            for(size_t c=0;c<nc;c++)fixed[c]=(uint8_t)!core[c];
            memcpy(saved,f->uv,2*f->nv*sizeof(double));
            for(size_t i=0;i<f->n_seams;i++){AsmSeamMeasure m;AsmField_seam(f,i,f->uv,&m);passing[i]=(uint8_t)(f->seams[i].required && m.pass);}
            apa_require(f,si,1);AsmRepairOpts local=*opts;local.fixed_charts=fixed;local.admission=local.regional=0;
            size_t moving_charts=0;for(size_t c=0;c<nc;c++)moving_charts+=!fixed[c] && f->vertex_offset[c]!=SIZE_MAX;
            fprintf(stderr,"[assemble source repair] trial %zu: relation %zu (%d,%d), %zu moving charts, initial RMS %.9g\n",stats->source_trials,seam->relation,seam->a,seam->b,moving_charts,before.rms);
            AsmRepairStats part={0};part.first_metric_chart=-1;
            int rc=apr_solve(f,&local,&part);AsmField_seam(f,si,f->uv,&after);accept=!rc && after.pass;reason=rc==1?"contact_walk_over_budget":rc?"candidate_solve_failed":"source_seam_failed";
            if(rc==1)stats->trials_refused_budget++;
            if(!rc && !accept){
                /* A pinned causal neighborhood may be too small. Retry the
                 * complete smaller source component from the same seed. */
                int different=0;size_t component_charts=0;
                for(size_t c=0;c<nc;c++){different|=fixed[c]!=(apr_root(parent,c)!=moving);component_charts+=apr_root(parent,c)==moving && f->vertex_offset[c]!=SIZE_MAX;}
                if(different && component_charts>(size_t)ASM_SOURCE_TRIAL_EXPAND_MAX_CHARTS){
                    fprintf(stderr,"[assemble source repair] relation %zu: complete source component has %zu charts, beyond the %d-chart retry cap; not expanded\n",seam->relation,component_charts,ASM_SOURCE_TRIAL_EXPAND_MAX_CHARTS);
                    different=0;
                }
                if(different){
                    memcpy(f->uv,saved,2*f->nv*sizeof(double));
                    for(size_t c=0;c<nc;c++)fixed[c]=(uint8_t)(apr_root(parent,c)!=moving);
                    fprintf(stderr,"[assemble source repair] relation %zu: expanding causal neighborhood to complete source component\n",seam->relation);
                    memset(&part,0,sizeof part);part.first_metric_chart=-1;rc=apr_solve(f,&local,&part);
                    AsmField_seam(f,si,f->uv,&after);accept=!rc && after.pass;reason=rc==1?"contact_walk_over_budget":rc?"candidate_solve_failed":"source_seam_failed";
                    if(rc==1)stats->trials_refused_budget++;
                }
            }
            for(size_t i=0;i<f->n_seams && accept;i++)if(passing[i]){
                AsmSeamMeasure m;AsmField_seam(f,i,f->uv,&m);if(!m.pass){accept=0;reason="existing_seam_regression";}
            }
            for(size_t c=0;c<nc && accept;c++)if(fixed[c] && f->vertex_offset[c]!=SIZE_MAX &&
                memcmp(saved+2*f->vertex_offset[c],f->uv+2*f->vertex_offset[c],2*run->charts[c].nv*sizeof(double))){accept=0;reason="fixed_context_changed";}
            if(!accept){memcpy(f->uv,saved,2*f->nv*sizeof(double));apa_require(f,si,0);}
            else reason="accepted_complete_source_transaction";
        }else after=before;
        if(accept){apa_require(f,si,1);run->rels[seam->relation].continuity|=ASM_CONT_REQUIRED;stats->recovered_connections++;}
        if(accept && !before.pass && opts->checkpoint && opts->checkpoint(f,opts->checkpoint_context)){if(log)fclose(log);Arena_dispose(&arena);return -1;}
        fprintf(stderr,"[assemble source repair] relation %zu (%d,%d): %s\n",seam->relation,seam->a,seam->b,reason);
        if(log && fprintf(log,"%zu,%d,%d,%.17g,%.17g,%d,%s\n",seam->relation,seam->a,seam->b,before.rms,after.rms,accept,reason)<0)io=0;
        if(log && fflush(log))io=0;
    }
    if(log && fclose(log))io=0;Arena_dispose(&arena);return io?0:-1;
}
