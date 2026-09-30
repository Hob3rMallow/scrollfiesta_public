/* Spend remaining time on missing branches of already represented regions.
 * The source graph ranks work; only the existing native admission transaction
 * can accept geometry. No new root, optional join or exclusion is introduced. */
typedef struct ArrStats {
    size_t targets,attempts,admitted,grown,expansions,index_builds,reactivated,continued,context_probes;
    size_t resume_loaded,resume_stale,resume_invalid,resume_saved;
    double added_area,seconds,index_seconds,context_seconds,resume_seconds;
    int expired;
} ArrStats;

static size_t arr_neighbors(const AsmRun *run,size_t target,uint8_t *moving);

/* A scheduling estimate only. The complete native transaction remains the
 * acceptance test, and stalled targets retain the same finite work cap. */
static double arr_remaining(double initial,double current,double spent,double limit)
{
    double left=limit;
    if(isfinite(initial) && isfinite(current) && initial>current && current>=0 && spent>0)
        left=spent*current/(initial-current);
    return fmax(ASM_REGION_REFINE_ESTIMATE_FLOOR_SEC,left);
}

/* A joint quantum that cut the measured violation below the ratio of its own
 * warm starting value is still converging (see ASM_REGION_REFINE_CONVERGING_SEC). */
static int arr_converging(int joint,double before,double after)
{
    return joint && isfinite(before) && isfinite(after) && before>0 && after<ASM_REGION_REFINE_CONVERGING_RATIO*before;
}

/* Seconds a target may spend in total: the patch budget, extended while it converges. */
static double arr_target_limit(int converging,double patch_seconds)
{
    return converging?fmax(patch_seconds,ASM_REGION_REFINE_CONVERGING_SEC):patch_seconds;
}

static int arr_cost(const AsmRun *run,size_t target,const uint8_t *moving,const AsmRepairOpts *opts,
                    size_t *faces,size_t *coordinates)
{
    *faces=*coordinates=0;
    for(size_t c=0;c<run->n_charts;c++)if(c==target || moving[c]){
        size_t nf=run->charts[c].nf,nv=run->charts[c].nv;
        if(nf>opts->patch_faces-*faces || nv>(opts->patch_coordinates-*coordinates)/2)return 0;
        *faces+=nf;*coordinates+=2*nv;
    }
    return *faces>0 && apr_workspace_estimate(*coordinates/2,*faces,0,run->n_charts)<=opts->patch_bytes;
}

static const int arr_hole_share=ASM_REGION_REFINE_HOLE_SHARE;
/* A hole waiting for its joint rounds (see ASM_REGION_REFINE_HOLE_SHARE):
 * enclosed by at least `sides` placed partners, surveyed and requeued. */
static int arr_is_waiting_hole(size_t anchors,uint8_t tried,size_t rounds,size_t sides)
{
    return anchors>=sides && !tried && rounds>=1;
}

static const int arr_gateway_share=ASM_REGION_REFINE_GATEWAY_SHARE;
/* A requeued target's value as a gateway (see ASM_REGION_REFINE_GATEWAY_SHARE):
 * the missing branch it opens, per forecast second, or 0 when nothing lies
 * behind it (its branch is itself). */
static double arr_gateway_rate(double branch_area,double own_area,double forecast)
{
    if(!(branch_area>own_area*(1+1e-9)))return 0;
    return branch_area/fmax(ASM_REGION_REFINE_ESTIMATE_FLOOR_SEC,isfinite(forecast)?forecast:ASM_REGION_REFINE_ESTIMATE_FLOOR_SEC);
}

static size_t arr_next(const AsmRun *run,const uint8_t *tried,size_t *group,
                       size_t *anchors,double *area,double minimum,double *potential,
                       const AsmRepairOpts *opts,uint8_t *scratch,const double *spent,uint8_t *const *pending,
                       const size_t *rounds,const double *represented,const double *forecast,int holes,int gateways)
{
    size_t nc=run->n_charts,best=SIZE_MAX,hole=SIZE_MAX,gate=SIZE_MAX;double best_score=0,hole_score=0,gate_rate=0;*potential=0;
    for(size_t c=0;c<nc;c++){group[c]=c;anchors[c]=0;area[c]=0;}
    for(size_t i=0;i<run->n_rels;i++){
        const AsmRelation *r=run->rels+i;if(!apa_required_source(r))continue;
        int a=AsmChart_registered(run->charts+r->a),b=AsmChart_registered(run->charts+r->b);
        if(a!=b){
            size_t c=(size_t)(a?r->b:r->a),known=(size_t)(a?r->a:r->b);
            if(AsmChart_in_layout(run->charts+c) && (!opts->admission_targets || opts->admission_targets[known]))anchors[c]++;
        }
        else if(!a && AsmChart_in_layout(run->charts+r->a) && AsmChart_in_layout(run->charts+r->b)){
            size_t x=apr_root(group,r->a),y=apr_root(group,r->b);group[x>y?x:y]=x<y?x:y;
        }
    }
    for(size_t c=0;c<nc;c++){
        group[c]=apr_root(group,c);
        if(AsmChart_in_layout(run->charts+c) && !AsmChart_registered(run->charts+c))area[group[c]]+=run->charts[c].area3d;
    }
    for(size_t c=0;c<nc;c++)if(anchors[c] && !tried[c] && area[group[c]]>=minimum &&
                              (!opts->admission_targets || opts->admission_targets[c])){
        size_t faces=0,coordinates=0;
        if(pending && pending[c])memcpy(scratch,pending[c],nc);else arr_neighbors(run,c,scratch);
        if(opts->fixed_charts)for(size_t k=0;k<nc;k++)if(opts->fixed_charts[k])scratch[k]=0;
        if(!arr_cost(run,c,scratch,opts,&faces,&coordinates))continue;
        /* Face count is a resource bound, not evidence that dense original
         * material is less useful. Survey every frontier first, then share
         * bounded rounds within the largest represented source region. */
        double score=run->charts[c].area3d;
        if(forecast && rounds && rounds[c]>1)score/=fmax(ASM_REGION_REFINE_ESTIMATE_FLOOR_SEC,forecast[c]);
        else if(spent)score/=1+spent[c]/ASM_REGION_REFINE_ATTEMPT_SEC;
        if(holes && rounds && arr_is_waiting_hole(anchors[c],tried[c],rounds[c],ASM_REGION_REFINE_HOLE_SIDES) &&
           (hole==SIZE_MAX || score>hole_score)){hole=c;hole_score=score;}
        if(gateways){   /* surveyed or not: a gateway's first quantum is its survey */
            double rate=arr_gateway_rate(area[group[c]],run->charts[c].area3d,forecast?forecast[c]:INFINITY);
            if(rate>0 && (gate==SIZE_MAX || rate>gate_rate)){gate=c;gate_rate=rate;}
        }
        int better=best==SIZE_MAX;
        if(best!=SIZE_MAX){
            if(rounds && (!rounds[c])!=(!rounds[best]))better=!rounds[c];
            else if(represented && represented[c]!=represented[best])better=represented[c]>represented[best];
            else if(rounds && (rounds[c]==1)!=(rounds[best]==1))better=rounds[c]==1;
            else better=score>best_score || (score==best_score &&
                (area[group[c]]>area[group[best]] || (area[group[c]]==area[group[best]] && anchors[c]>anchors[best])));
        }
        if(better){
            best=c;best_score=score;
        }
    }
    if(gate!=SIZE_MAX){*potential=area[group[gate]];return gate;}
    if(hole!=SIZE_MAX){*potential=area[group[hole]];return hole;}
    if(best!=SIZE_MAX)*potential=area[group[best]];return best;
}

static size_t arr_neighbors(const AsmRun *run,size_t target,uint8_t *moving)
{
    memset(moving,0,run->n_charts);size_t count=0;
    for(size_t i=0;i<run->n_rels;i++){
        const AsmRelation *r=run->rels+i;if(!apa_required_source(r))continue;
        size_t c=r->a==(int32_t)target?(size_t)r->b:r->b==(int32_t)target?(size_t)r->a:SIZE_MAX;
        if(c!=SIZE_MAX && AsmChart_registered(run->charts+c) && !moving[c]){moving[c]=1;count++;}
    }
    return count;
}

/* Expand only toward a measured fixed collider or a constrained boundary.
 * All physical obligations retain their existing interpretation. */
static size_t arr_expand_mode(const AsmRun *run,size_t target,uint8_t *moving,const AsmRepairOpts *opts,
                               const ApaFeedback *feedback,int boundary)
{
    if(!feedback->complete)return 0;
    size_t faces=0,coordinates=0,added=0;
    if(!arr_cost(run,target,moving,opts,&faces,&coordinates))return 0;
    for(int phase=0;phase<(boundary?2:1);phase++)for(;;){
        size_t best=SIZE_MAX;double best_score=0;
        for(size_t c=0;c<run->n_charts;c++){
            if(c==target || moving[c] || !AsmChart_registered(run->charts+c) ||
               (opts->fixed_charts && opts->fixed_charts[c]))continue;
            size_t nf=run->charts[c].nf,nv=run->charts[c].nv;
            if(nf>opts->patch_faces-faces || nv>(opts->patch_coordinates-coordinates)/2)continue;
            if(apr_workspace_estimate(coordinates/2+nv,faces+nf,0,run->n_charts)>opts->patch_bytes)continue;
            double score=phase?feedback->seam_ratio[c]:feedback->contact_area[c];
            if(score>0 && (best==SIZE_MAX || score>best_score)){best=c;best_score=score;}
        }
        if(best==SIZE_MAX)break;
        moving[best]=1;faces+=run->charts[best].nf;coordinates+=2*run->charts[best].nv;added++;
        if(phase)break; /* One boundary step; measure again before going farther. */
    }
    return added;
}

static size_t arr_expand(const AsmRun *run,size_t target,uint8_t *moving,const AsmRepairOpts *opts,
                          const ApaFeedback *feedback)
{
    return arr_expand_mode(run,target,moving,opts,feedback,1);
}

static void arr_feedback_reset(ApaFeedback *feedback)
{
    memset(feedback->contact_area,0,feedback->count*sizeof(double));
    memset(feedback->seam_ratio,0,feedback->count*sizeof(double));
    memset(feedback->context,0,feedback->count);feedback->complete=0;
    feedback->violation=feedback->material_area=feedback->blocking_material_area=0;
}

#include "asm_repair_resume.inc"

/* A transaction in this pass changes only qualified placements and source
 * obligations, never original geometry or parity. Its neighbors' saved UVs
 * remain useful guesses. Every retry rebuilds current fixed context and
 * rechecks all seams, metrics and local/exterior contacts before publishing. */
static void arr_reactivate(ApaContinuation *state,uint8_t **pending,size_t *round,
                           uint8_t *tried,double *spent,double *forecast,double limit)
{
    *tried=0;*spent=0;*forecast=limit;
    if(!state->arena){*round=0;*pending=NULL;}
}

static int arr_work(AsmRun *run,const char *dir,double deadline,double minimum,
                    const AsmRepairOpts *policy,const AsmField *prepared,ArrStats *out)
{
    if(!run || !dir || !out || !policy || !policy->patch_faces || !policy->patch_coordinates ||
       !isfinite(deadline) || !(minimum>0) || !isfinite(minimum))return -1;
    /* This pass has one owner in each private OS process. */
    AsmField_threads(policy->threads);
    memset(out,0,sizeof *out);double started=ves_clock_sec();int rc=-1;
    Arena_T arena=Arena_new(),field_arena=Arena_new();FILE *ledger=NULL;char path[2048];
    size_t nc=run->n_charts,initial_count=0;double initial_area=0;
    ApaOverlay overlay={0};overlay.arena=arena;overlay.batched=overlay.query_scaled_budget=1;
    overlay.n_charts=nc;overlay.changed_charts=ARENA_CALLOC(arena,nc,1);overlay.threads=policy->threads;
    overlay.incremental=ASM_REPAIR_OVERLAY_INCREMENTAL;
    size_t *group=ARENA_ALLOC(arena,nc*sizeof *group),*anchors=ARENA_ALLOC(arena,nc*sizeof *anchors);
    double *area=ARENA_ALLOC(arena,nc*sizeof *area);
    uint8_t *tried=ARENA_CALLOC(arena,nc,1),*core=ARENA_CALLOC(arena,nc,1),*moving=ARENA_CALLOC(arena,nc,1);
    uint8_t *converging=ARENA_CALLOC(arena,nc,1);   /* last joint quantum still cut the violation */
    const char **last_reason=ARENA_CALLOC(arena,nc,sizeof *last_reason);
    uint8_t **dependencies=ARENA_CALLOC(arena,nc,sizeof *dependencies);
    uint8_t **pending=ARENA_CALLOC(arena,nc,sizeof *pending),*seen=ARENA_CALLOC(arena,nc,1);
    double *spent=ARENA_CALLOC(arena,nc,sizeof *spent);
    double *forecast=ARENA_ALLOC(arena,nc*sizeof *forecast);
    for(size_t c=0;c<nc;c++)forecast[c]=policy->patch_seconds;
    size_t *rounds=ARENA_CALLOC(arena,nc,sizeof *rounds),*region=ARENA_ALLOC(arena,nc*sizeof *region);
    double *represented=ARENA_CALLOC(arena,nc,sizeof *represented),*region_area=ARENA_CALLOC(arena,nc,sizeof *region_area);
    ApaContinuation *continuations=ARENA_CALLOC(arena,nc,sizeof *continuations);
    double resume_started=ves_clock_sec();
    if(arr_resume_read(run,policy,arena,continuations,pending,dependencies,rounds,forecast,out))
        fprintf(stderr,"[region refine] malformed or incompatible private continuation ignored\n");
    out->resume_seconds+=ves_clock_sec()-resume_started;
    for(size_t c=0;c<nc;c++)region[c]=c;
    for(size_t i=0;i<run->n_rels;i++)if(apa_required_source(run->rels+i)){
        const AsmRelation *r=run->rels+i;
        if(AsmChart_in_layout(run->charts+r->a) && AsmChart_in_layout(run->charts+r->b)){
            size_t a=apr_root(region,r->a),b=apr_root(region,r->b);region[a>b?a:b]=a<b?a:b;
        }
    }
    for(size_t c=0;c<nc;c++)if(AsmChart_registered(run->charts+c))region_area[apr_root(region,c)]+=run->charts[c].area3d;
    for(size_t c=0;c<nc;c++)represented[c]=region_area[apr_root(region,c)];
    ApaFeedback feedback={0};feedback.count=nc;
    feedback.contact_area=ARENA_CALLOC(arena,nc,sizeof(double));feedback.seam_ratio=ARENA_CALLOC(arena,nc,sizeof(double));
    feedback.context=ARENA_CALLOC(arena,nc,1);
    ApaBox *boxes=ARENA_ALLOC(arena,nc*sizeof *boxes);const AsmField *incumbent=prepared;int indexed=0;
    for(size_t c=0;c<nc;c++)if(AsmChart_registered(run->charts+c)){initial_count++;initial_area+=run->charts[c].area3d;}
    ves_mkdir(dir);if(snprintf(path,sizeof path,"%s/attempts.csv",dir)>=(int)sizeof path)goto done;
    ledger=fopen(path,"wb");if(!ledger)goto done;
    if(fputs("chart,potential_area,moving_charts,accepted,reason,seconds,solver_seconds,phase,moving_faces,moving_coordinates,feedback_complete,violation_before,violation_after,estimated_remaining_seconds\n",ledger)<0)goto done;
    unsigned hole_turn=0,gate_turn=0;   /* every other pick may take a waiting hole / the best gateway */
    while(ves_clock_sec()<deadline){
        double potential=0;size_t target=arr_next(run,tried,group,anchors,area,minimum,&potential,policy,moving,spent,pending,rounds,represented,forecast,
                                                  arr_hole_share && (hole_turn++&1),arr_gateway_share && (gate_turn++&1));
        if(target==SIZE_MAX)break;tried[target]=1;if(!seen[target]){seen[target]=1;out->targets++;}core[target]=1;
        if(!indexed){
            double index_start=ves_clock_sec();
            if(!incumbent){incumbent=AsmField_new(field_arena,run);out->index_builds++;}
            if(!incumbent || AsmContacts_refit(incumbent->contacts,incumbent->uv))goto done;
            AsmContacts_threads(incumbent->contacts,policy->threads);apa_chart_boxes(run,boxes);
            indexed=1;out->index_seconds+=ves_clock_sec()-index_start;
        }
        size_t neighbors=0;int accepted=0;
        if(pending[target]){memcpy(moving,pending[target],nc);for(size_t c=0;c<nc;c++)neighbors+=moving[c]!=0;}
        else neighbors=arr_neighbors(run,target,moving);
        if(policy->fixed_charts){neighbors=0;for(size_t c=0;c<nc;c++){
            if(policy->fixed_charts[c])moving[c]=0;neighbors+=moving[c]!=0;
        }}
        feedback.deadline=fmin(deadline,ves_clock_sec()+fmax(0,arr_target_limit(converging[target],policy->patch_seconds)-spent[target]));
        for(size_t phase=rounds[target];ves_clock_sec()<feedback.deadline;){
            int joint=phase>0;if(joint && !neighbors)break;
            double before=INFINITY;
            if(joint){
                /* Close over measured colliders of the whole movable set
                 * before paying for deformation. A newly included neighbor
                 * can itself touch another fixed chart. No solve or live
                 * publication occurs in this observation-only transaction. */
                AsmRepairOpts probe=*policy;probe.compact=1;
                probe.iterations=probe.clearance_iterations=probe.closure_iterations=0;
                probe.admission_rounds=1;probe.regional=probe.source_trials=0;
                probe.admission_moving=moving;probe.admission_continuation=continuations+target;
                probe.admission_feedback=&feedback;probe.deadline=feedback.deadline;feedback.observe_only=1;
                while(ves_clock_sec()<feedback.deadline){
                    arr_feedback_reset(&feedback);size_t fixed=0;const char *reason="";
                    AsmRepairStats stats={0};ApaReport report={0};
                    if(apa_attempt_anchor(run,&probe,core,&fixed,&stats,&reason,NULL,boxes,incumbent,&overlay,&report,4,SIZE_MAX))goto done;
                    out->context_probes++;out->context_seconds+=report.sec;spent[target]+=report.sec;
                    if(!dependencies[target])dependencies[target]=ARENA_CALLOC(arena,nc,1);
                    for(size_t c=0;c<nc;c++)dependencies[target][c]|=feedback.context[c] || moving[c] || c==target;
                    size_t added=arr_expand_mode(run,target,moving,&probe,&feedback,0);neighbors+=added;
                    if(!added){if(feedback.complete)before=feedback.violation;break;}out->expansions++;
                    if(!pending[target])pending[target]=ARENA_ALLOC(arena,nc);
                    memcpy(pending[target],moving,nc);tried[target]=0;
                }
                feedback.observe_only=0;
                if(ves_clock_sec()>=feedback.deadline)break;tried[target]=1;
            }
            AsmRepairOpts opts=*policy;opts.compact=1;
            opts.iterations=opts.clearance_iterations=joint?ASM_REGION_REFINE_ITERATIONS:0;
            opts.closure_iterations=0;opts.clearance_method=0;opts.admission_rounds=1;
            opts.regional=opts.source_trials=0;opts.admission_moving=joint?moving:NULL;
            double quantum=phase==1?ASM_REGION_REFINE_FIRST_SEC:ASM_REGION_REFINE_ATTEMPT_SEC;
            opts.deadline=fmin(feedback.deadline,ves_clock_sec()+quantum);
            opts.admission_feedback=&feedback;arr_feedback_reset(&feedback);
            opts.admission_continuation=joint?continuations+target:NULL;
            continuations[target].changed=0;
            size_t moving_faces=0,moving_coordinates=0;
            if(!arr_cost(run,target,joint?moving:core,&opts,&moving_faces,&moving_coordinates))break;
            AsmRepairStats stats={0};ApaReport report={0};size_t fixed=0;const char *reason="";
            accepted=apa_attempt_anchor(run,&opts,core,&fixed,&stats,&reason,NULL,boxes,incumbent,&overlay,&report,4,SIZE_MAX);
            last_reason[target]=reason;
            out->attempts++;spent[target]+=report.sec;rounds[target]++;
            /* Compare the latest quantum against its own warm starting
             * state in the same context. Earlier progress cannot hide a
             * subsequent stall, nor can context expansion fake progress. */
            double after=accepted?0:feedback.complete?feedback.violation:INFINITY;
            forecast[target]=arr_remaining(before,after,report.sec,policy->patch_seconds);
            converging[target]=(uint8_t)arr_converging(joint,before,after);
            if(fprintf(ledger,"%zu,%.17g,",target,potential)<0)goto done;
            if(joint){int first=1;for(size_t c=0;c<nc;c++)if(moving[c]){
                if(fprintf(ledger,"%s%zu",first?"":";",c)<0)goto done;first=0;
            }}
            if(fprintf(ledger,",%d,%s,%.17g,%.17g,%zu,%zu,%zu,%d,%.17g,%.17g,%.17g\n",accepted,reason,report.sec,report.solver_seconds,
                       phase,moving_faces,moving_coordinates,feedback.complete,before,after,forecast[target])<0 || fflush(ledger))goto done;
            fprintf(stderr,"[region refine] chart %zu, %s: %s, %.3f s\n",target,joint?"neighbors movable":"fixed neighbors",reason,report.sec);
            if(accepted){out->admitted++;break;}
            if(!dependencies[target])dependencies[target]=ARENA_CALLOC(arena,nc,1);
            for(size_t c=0;c<nc;c++)dependencies[target][c]|=feedback.context[c] || moving[c] || c==target;
            {
                size_t added=arr_expand(run,target,moving,&opts,&feedback);
                int retry=added || continuations[target].changed || (!joint && feedback.complete);
                if(retry && spent[target]<arr_target_limit(converging[target],policy->patch_seconds)){
                    if(!pending[target])pending[target]=ARENA_ALLOC(arena,nc);
                    memcpy(pending[target],moving,nc);tried[target]=0;if(added)out->expansions++;
                    if(continuations[target].changed)out->continued++;
                    fprintf(stderr,"[region refine] requeued with %zu additional measured context charts\n",added);
                }
                break;
            }
        }
        core[target]=0;
        if(accepted){
            apa_continuation_clear(continuations+target);
            for(size_t c=0;c<nc;c++)if(dependencies[c] && !AsmChart_registered(run->charts+c))
                for(size_t k=0;k<nc;k++)if((k==target || moving[k]) && dependencies[c][k]){
                    arr_reactivate(continuations+c,pending+c,rounds+c,tried+c,spent+c,forecast+c,policy->patch_seconds);
                    out->reactivated++;break;
                }
            /* The replacement overlay now owns current versions of moved
             * charts; the much larger immutable index remains reusable. */
            for(size_t c=0;c<nc;c++)if(c==target || moving[c])apa_chart_box(run->charts+c,boxes+c);
            /* The next detection pass discovers exposed successors itself;
             * a second whole-field growth driver would duplicate setup. */
        }
    }
    {
        /* Report actual deferral separately from a historical geometry
         * refusal. A chart the queue never reached has not failed a solve. */
        double potential=0;
        (void)arr_next(run,tried,group,anchors,area,minimum,&potential,policy,moving,spent,pending,rounds,represented,forecast,0,0);
        snprintf(path,sizeof path,"%s/queue.csv",dir);FILE *queue=fopen(path,"wb");if(!queue)goto done;
        int okay=fputs("chart,state,attempts,spent_seconds,moving_faces,moving_coordinates,represented_region_area,estimated_remaining_seconds\n",queue)>=0;
        for(size_t c=0;c<nc && okay;c++)if(anchors[c] && !AsmChart_registered(run->charts+c)){
            if(pending[c])memcpy(moving,pending[c],nc);else arr_neighbors(run,c,moving);
            if(policy->fixed_charts)for(size_t k=0;k<nc;k++)if(policy->fixed_charts[k])moving[k]=0;
            size_t faces=0,coordinates=0;int fits=arr_cost(run,c,moving,policy,&faces,&coordinates);
            const char *state=!fits?"workspace_or_index_budget":!rounds[c]?"not_surveyed_budget":
                spent[c]>=arr_target_limit(converging[c],policy->patch_seconds)?"target_budget_spent":tried[c] && last_reason[c] &&
                    (strstr(last_reason[c],"resource_cap") || !strcmp(last_reason[c],"workspace_memory_budget"))?last_reason[c]:tried[c]?"no_further_progress":
                rounds[c]==1?"surveyed_waiting_for_joint":"continuation_waiting";
            okay=fprintf(queue,"%zu,%s,%zu,%.17g,%zu,%zu,%.17g,%.17g\n",c,state,rounds[c],spent[c],faces,coordinates,represented[c],forecast[c])>=0;
        }
        if(fclose(queue))okay=0;if(!okay)goto done;
    }
    {
        size_t final_count=0;double final_area=0;
        for(size_t c=0;c<nc;c++)if(AsmChart_registered(run->charts+c)){final_count++;final_area+=run->charts[c].area3d;}
        if(final_count<initial_count+out->admitted || final_area<initial_area)goto done;
        out->grown=final_count-initial_count-out->admitted;out->added_area=final_area-initial_area;
    }
    resume_started=ves_clock_sec();
    if(arr_resume_write(run,continuations,pending,dependencies,rounds,forecast,out))goto done;
    out->resume_seconds+=ves_clock_sec()-resume_started;
    out->seconds=ves_clock_sec()-started;out->expired=ves_clock_sec()>=deadline;
    if(fclose(ledger)){ledger=NULL;goto done;}ledger=NULL;
    snprintf(path,sizeof path,"%s/region_refine.json",dir);FILE *report=fopen(path,"wb");if(!report)goto done;
    int okay=fprintf(report,"{\"schema\":\"native-required-region-refinement-v1\",\"targets\":%zu,\"attempts\":%zu,"
        "\"admitted_charts\":%zu,\"grown_charts\":%zu,\"added_area\":%.17g,\"seconds\":%.17g,\"deadline_reached\":%s,"
        "\"minimum_branch_potential\":%.17g,\"iterations\":%d,\"candidate_seconds_limit\":%.17g,"
        "\"contact_threads\":%d,\"metric_margin\":%.17g,\"clearance_method\":0,\"expansions\":%zu,"
        "\"reactivated_targets\":%zu,\"index_builds\":%zu,\"index_seconds\":%.17g,"
        "\"reused_incumbent\":%s,\"replacement_index_builds\":%zu,\"replacement_index_seconds\":%.17g,\"replacement_faces\":%zu,"
        "\"replacement_incremental\":%s,\"replacement_rebuilt_faces\":%zu,\"replacement_compactions\":%zu,"
        "\"selection\":\"survey_then_blocking_progress_per_remaining_second\",\"continued_candidates\":%zu,"
        "\"context_probes\":%zu,\"context_seconds\":%.17g,\"resume_loaded\":%zu,\"resume_stale\":%zu,"
        "\"resume_invalid\":%zu,\"resume_saved\":%zu,\"resume_bytes\":%zu,\"resume_seconds\":%.17g,\"audit_performed\":false}\n",
        out->targets,out->attempts,out->admitted,out->grown,out->added_area,out->seconds,out->expired?"true":"false",
        minimum,ASM_REGION_REFINE_ITERATIONS,(double)ASM_REGION_REFINE_ATTEMPT_SEC,
        policy->threads,policy->admission_metric_margin,out->expansions,out->reactivated,out->index_builds,out->index_seconds,
        prepared?"true":"false",overlay.rebuilds,overlay.build_seconds,overlay.nf,
        overlay.incremental?"true":"false",overlay.rebuilt_faces,overlay.compactions,out->continued,out->context_probes,out->context_seconds,
        out->resume_loaded,out->resume_stale,out->resume_invalid,out->resume_saved,run->repair_resume_size,out->resume_seconds)>=0;
    if(fclose(report))okay=0;if(okay)rc=0;
done:
    for(size_t c=0;c<nc;c++)apa_continuation_clear(continuations+c);
    if(ledger)fclose(ledger);apa_overlay_clear(&overlay);Arena_dispose(&field_arena);Arena_dispose(&arena);return rc;
}

static int arr_run(AsmRun *run,const char *dir,double deadline,double minimum,ArrStats *out)
{
    AsmRepairOpts policy;AsmRepair_defaults(&policy);policy.compact=1;
    policy.threads=ASM_REGION_REFINE_CONTACT_THREADS;policy.admission_metric_margin=ASM_REGION_REFINE_METRIC_MARGIN;
    int rc=arr_work(run,dir,deadline,minimum,&policy,NULL,out);AsmField_threads(1);return rc;
}
