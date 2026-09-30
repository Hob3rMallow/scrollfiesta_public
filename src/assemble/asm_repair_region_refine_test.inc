#include "asm_repair_admission_contacts_test.inc"

static int arr_margin_control(void)
{
    /* The reserve mechanism itself, at an explicit test reserve: the configured
     * refinement reserve may be 0 when the solver keeps faces in band. */
    const double reserve=1e-8;
    Arena_T arena=Arena_new();AsmRun run={0};int failures=0;
    if(ags_fixture(arena,&run)){Arena_dispose(&arena);return 1;}
    AsmField *field=AsmField_new(arena,&run);
    if(!field){Arena_dispose(&arena);return 1;}
    if(!apa_chart_metric_margin(field,0,reserve) ||
       apa_chart_metric_margin(field,1,0) || apa_chart_metric_margin(field,0,-1) ||
       apa_chart_metric_margin(field,0,NAN) || apa_chart_metric_margin(field,0,.25)){failures++;fprintf(stderr,"  metric reserve initial control failed\n");}
    size_t first=field->vertex_offset[0];
    for(size_t v=0;v<run.charts[0].nv;v++){
        field->uv[2*(first+v)]=1.25*(double)(v%9);
        field->uv[2*(first+v)+1]=(double)(v/9);
    }
    if(!apa_chart_metric_margin(field,0,0) || apa_chart_metric_margin(field,0,reserve)){failures++;fprintf(stderr,"  metric reserve upper boundary failed\n");}
    for(size_t v=0;v<run.charts[0].nv;v++)field->uv[2*(first+v)]=.75*(double)(v%9);
    if(!apa_chart_metric_margin(field,0,0) || apa_chart_metric_margin(field,0,reserve)){failures++;fprintf(stderr,"  metric reserve lower boundary failed\n");}
    for(size_t v=0;v<run.charts[0].nv;v++)field->uv[2*(first+v)]=1.125*(double)(v%9);
    if(!apa_chart_metric_margin(field,0,reserve))failures++;
    /* A nearly rigid new chart passes all old metric/contact/seam checks,
     * but a deliberately tight reserve must refuse without publishing it. */
    run.charts[2].flags|=ASM_CHART_EXTRAS;run.charts[3].flags|=ASM_CHART_EXTRAS;
    run.charts[1].uv[2*(4*9+4)]+=.05f;
    AsmField *incumbent=AsmField_new(arena,&run);uint8_t core[5]={0,1,0,0,0};ApaBox boxes[5];
    if(!incumbent || AsmContacts_refit(incumbent->contacts,incumbent->uv))failures++;
    else{
        apa_chart_boxes(&run,boxes);AsmChart before[5];memcpy(before,run.charts,sizeof before);
        AsmRepairOpts opts;AsmRepair_defaults(&opts);opts.iterations=opts.clearance_iterations=opts.closure_iterations=0;
        opts.admission_rounds=1;opts.admission_metric_margin=.24;
        AsmRepairStats stats={0};ApaReport report={0};size_t anchors=0;const char *reason="";
        if(apa_attempt_anchor(&run,&opts,core,&anchors,&stats,&reason,NULL,boxes,incumbent,NULL,&report,4,SIZE_MAX) ||
           strcmp(reason,"packing_metric_margin") || memcmp(before,run.charts,sizeof before)){failures++;fprintf(stderr,"  metric reserve rejection failed: %s\n",reason);}
        opts.admission_metric_margin=0;
        if(!apa_attempt_anchor(&run,&opts,core,&anchors,&stats,&reason,NULL,boxes,incumbent,NULL,&report,4,SIZE_MAX) ||
           !AsmChart_registered(run.charts+1)){failures++;fprintf(stderr,"  metric reserve release failed: %s\n",reason);}
    }
    Arena_dispose(&arena);return failures;
}

static int arr_policy_control(void)
{
    Arena_T arena=Arena_new();AsmRun run={0};int failures=0;
    if(ags_fixture(arena,&run)){Arena_dispose(&arena);return 1;}
    run.charts[3].placed=0;run.charts[3].placement_state=ASM_PLACE_NONE;
    AsmRepairOpts policy;AsmRepair_defaults(&policy);
    size_t group[5],anchors[5];double area[5],potential=0,spent[5]={0};uint8_t tried[5]={0},moving[5]={0};
    if(arr_next(&run,tried,group,anchors,area,1,&potential,&policy,moving,spent,NULL,NULL,NULL,NULL,0,0)!=1)failures++;
    size_t original_faces=run.charts[1].nf;run.charts[1].nf*=3000;
    if(arr_next(&run,tried,group,anchors,area,1,&potential,&policy,moving,spent,NULL,NULL,NULL,NULL,0,0)!=1)failures++;
    run.charts[1].nf=original_faces;
    spent[1]=ASM_REGION_REFINE_ATTEMPT_SEC;
    if(arr_next(&run,tried,group,anchors,area,1,&potential,&policy,moving,spent,NULL,NULL,NULL,NULL,0,0)!=3)failures++;
    size_t rounds[5]={0,1,0,0,0};double represented[5]={1000,1000,64,64,1000};
    if(arr_next(&run,tried,group,anchors,area,1,&potential,&policy,moving,spent,NULL,rounds,represented,NULL,0,0)!=3)failures++;
    rounds[3]=1;
    if(arr_next(&run,tried,group,anchors,area,1,&potential,&policy,moving,spent,NULL,rounds,represented,NULL,0,0)!=1)failures++;
    represented[3]=1000;rounds[3]=2;
    if(arr_next(&run,tried,group,anchors,area,1,&potential,&policy,moving,spent,NULL,rounds,represented,NULL,0,0)!=1)failures++;
    double forecast[5]={0,20,0,2,0};rounds[1]=2;rounds[3]=4;
    if(arr_next(&run,tried,group,anchors,area,1,&potential,&policy,moving,spent,NULL,rounds,represented,forecast,0,0)!=3)failures++;
    rounds[1]=1;
    if(arr_next(&run,tried,group,anchors,area,1,&potential,&policy,moving,spent,NULL,rounds,represented,forecast,0,0)!=1)failures++;
    if(arr_remaining(3,2,8,45)!=16 || arr_remaining(3,.25,8,45)!=ASM_REGION_REFINE_ESTIMATE_FLOOR_SEC ||
       arr_remaining(3,3,8,45)!=45 || arr_remaining(INFINITY,1,8,45)!=45)failures++;
    policy.patch_faces=2*run.charts[0].nf-1;
    if(arr_next(&run,tried,group,anchors,area,1,&potential,&policy,moving,spent,NULL,NULL,NULL,NULL,0,0)!=SIZE_MAX)failures++;
    AsmRepair_defaults(&policy);
    policy.patch_bytes=1;
    if(arr_next(&run,tried,group,anchors,area,1,&potential,&policy,moving,spent,NULL,NULL,NULL,NULL,0,0)!=SIZE_MAX)failures++;
    AsmRepair_defaults(&policy);
    /* Relabel the same source graph. Selection must follow the geometry
     * and effort ledger rather than a particular chart number. */
    const int permutation[5]={1,3,0,2,4};AsmChart remapped[5];AsmRelation relations[3];double reordered[5]={0};
    for(int c=0;c<5;c++){remapped[permutation[c]]=run.charts[c];remapped[permutation[c]].id=permutation[c];reordered[permutation[c]]=spent[c];}
    for(int i=0;i<3;i++){relations[i]=run.rels[i];relations[i].a=permutation[relations[i].a];relations[i].b=permutation[relations[i].b];}
    AsmRun renamed=run;renamed.charts=remapped;renamed.rels=relations;
    if(arr_next(&renamed,tried,group,anchors,area,1,&potential,&policy,moving,reordered,NULL,NULL,NULL,NULL,0,0)!=(size_t)permutation[3])failures++;
    fprintf(stderr,"  automatic neighborhood ranking, spent-time fairness and relabeling: %s\n",failures?"FAIL":"ok");
    Arena_dispose(&arena);return failures;
}

static int arr_obstacle_feedback_control(void)
{
    Arena_T arena=Arena_new();AsmRun run={0};int failures=0;
    if(ags_fixture(arena,&run)){Arena_dispose(&arena);return 1;}
    run.charts[4].flags|=ASM_CHART_EXTRAS;apr_test_scale(&run,32);
    AsmField *incumbent=AsmField_new(arena,&run);ApaBox boxes[5];apa_chart_boxes(&run,boxes);
    if(!incumbent || AsmContacts_refit(incumbent->contacts,incumbent->uv)){Arena_dispose(&arena);return 1;}
    uint8_t core[5]={0,1,0,0,0},moving[5]={0},context[5]={0};double contact[5]={0},seam[5]={0};
    ApaFeedback feedback={5,contact,seam,context,ves_clock_sec()+10,0};
    AsmRepairOpts opts;AsmRepair_defaults(&opts);opts.compact=1;
    opts.iterations=opts.clearance_iterations=opts.closure_iterations=0;opts.admission_rounds=1;
    opts.deadline=feedback.deadline;opts.admission_feedback=&feedback;
    AsmChart before[5];memcpy(before,run.charts,sizeof before);
    AsmRepairStats stats={0};ApaReport report={0};size_t anchors=0;const char *reason="";
    int accepted=apa_attempt_anchor(&run,&opts,core,&anchors,&stats,&reason,NULL,boxes,incumbent,NULL,&report,4,SIZE_MAX);
    if(accepted || !feedback.complete || !(contact[2]>0) || !context[2] || memcmp(before,run.charts,sizeof before))failures++;
    arr_neighbors(&run,1,moving);
    if(moving[2] || arr_expand(&run,1,moving,&opts,&feedback)!=1 || !moving[2])failures++;
    arr_neighbors(&run,1,moving);opts.patch_faces=2*run.charts[0].nf;
    if(arr_expand(&run,1,moving,&opts,&feedback) || moving[2])failures++;
    opts.patch_faces=ASM_REPAIR_PATCH_FACES;feedback.complete=0;
    if(arr_expand(&run,1,moving,&opts,&feedback))failures++;
    fprintf(stderr,"  measured non-neighbor blocker, rollback, cap and incomplete-feedback controls: %s (%s, area %.9g)\n",
            failures?"FAIL":"ok",reason,contact[2]);
    Arena_dispose(&arena);return failures;
}

static int arr_continuation_control(void)
{
    Arena_T arena=Arena_new();AsmRun run={0};int failures=0;ApaContinuation state={0};
    if(ags_fixture(arena,&run)){Arena_dispose(&arena);return 1;}
    AsmField *field=AsmField_new(arena,&run);if(!field){Arena_dispose(&arena);return 1;}
    size_t bytes=2*field->nv*sizeof(double);double *original=ARENA_ALLOC(arena,bytes),*shifted=ARENA_ALLOC(arena,bytes);
    memcpy(original,field->uv,bytes);uint8_t moving[5]={1,0,0,0,0};
    for(size_t v=0;v<run.charts[0].nv;v++)field->uv[2*(field->vertex_offset[0]+v)]-=32;
    memcpy(shifted,field->uv,bytes);apa_continuation_save(&state,field,moving,original);
    if(!state.changed || !state.arena)failures++;
    memcpy(field->uv,original,bytes);
    if(!apa_continuation_apply(&state,field,moving) || memcmp(shifted,field->uv,bytes))failures++;
    apa_continuation_save(&state,field,moving,original);if(state.changed)failures++;
    memcpy(field->uv,original,bytes);moving[0]=0;moving[2]=1;
    if(apa_continuation_apply(&state,field,moving) || memcmp(original,field->uv,bytes))failures++;
    moving[0]=1;
    run.charts[0].flags^=ASM_CHART_MIRROR;
    if(apa_continuation_apply(&state,field,moving) || memcmp(original,field->uv,bytes))failures++;
    run.charts[0].flags^=ASM_CHART_MIRROR;
    if(!apa_continuation_apply(&state,field,moving) || memcmp(shifted,field->uv,bytes))failures++;
    /* A retained pose with a failed incident seam must still be refused,
     * leaving the live run and original fixed frame byte-for-byte intact. */
    memcpy(field->uv,original,bytes);AsmChart before[5];memcpy(before,run.charts,sizeof before);
    uint8_t core[5]={0,1,0,0,0},neighbors[5]={1,0,0,0,0};ApaBox boxes[5];apa_chart_boxes(&run,boxes);
    uint8_t tried=1,*pending=neighbors;size_t round=7;double spent=45,forecast=10000;
    Arena_T retained=state.arena;
    arr_reactivate(&state,&pending,&round,&tried,&spent,&forecast,45);
    if(tried || spent || forecast!=45 || round!=7 || pending!=neighbors || state.arena!=retained)failures++;
    AsmRepairOpts opts;AsmRepair_defaults(&opts);opts.compact=1;
    opts.iterations=opts.clearance_iterations=opts.closure_iterations=0;opts.admission_rounds=1;
    opts.admission_moving=neighbors;opts.admission_continuation=&state;
    AsmRepairStats stats={0};ApaReport report={0};size_t anchors=0;const char *reason="";
    if(apa_attempt_anchor(&run,&opts,core,&anchors,&stats,&reason,NULL,boxes,field,NULL,&report,4,SIZE_MAX) ||
       strcmp(reason,"incident_source_seam_failed") || memcmp(before,run.charts,sizeof before))failures++;
    apa_continuation_clear(&state);if(state.arena || state.uv)failures++;
    tried=1;spent=12;arr_reactivate(&state,&pending,&round,&tried,&spent,&forecast,45);
    if(tried || spent || round || pending)failures++;
    fprintf(stderr,"  uncommitted continuation, wider context, unchanged seed, parity and seam rollback: %s\n",failures?"FAIL":"ok");
    Arena_dispose(&arena);return failures;
}

static int arr_context_survey_control(void)
{
    Arena_T arena=Arena_new();AsmRun run={0};int failures=0;
    if(ags_fixture(arena,&run)){Arena_dispose(&arena);return 1;}
    run.charts[2].flags|=ASM_CHART_EXTRAS;run.charts[3].flags|=ASM_CHART_EXTRAS;
    AsmField *field=AsmField_new(arena,&run);if(!field){Arena_dispose(&arena);return 1;}
    uint8_t core[5]={0,1,0,0,0},moving[5]={1,0,0,0,0},context[5]={0};double contact[5]={0},seam[5]={0};
    ApaFeedback feedback={5,contact,seam,context,ves_clock_sec()+5,0,1};
    ApaBox boxes[5];apa_chart_boxes(&run,boxes);AsmChart original[5];memcpy(original,run.charts,sizeof original);
    AsmRepairOpts opts;AsmRepair_defaults(&opts);opts.compact=1;
    opts.iterations=opts.clearance_iterations=opts.closure_iterations=0;opts.admission_rounds=1;
    opts.admission_feedback=&feedback;opts.admission_moving=moving;opts.deadline=feedback.deadline;
    AsmRepairStats stats={0};ApaReport report={0};size_t anchors=0;const char *reason="";
    if(apa_attempt_anchor(&run,&opts,core,&anchors,&stats,&reason,NULL,boxes,field,NULL,&report,4,SIZE_MAX) ||
       strcmp(reason,"contact_context_survey") || !feedback.complete || memcmp(original,run.charts,sizeof original))failures++;
    feedback.observe_only=0;
    if(!apa_attempt_anchor(&run,&opts,core,&anchors,&stats,&reason,NULL,boxes,field,NULL,&report,4,SIZE_MAX))failures++;
    /* A source boundary with no measured collision does not become a
     * contact-closure neighbor merely because it has seam pressure. */
    arr_feedback_reset(&feedback);feedback.complete=1;seam[0]=100;memset(moving,0,sizeof moving);
    if(arr_expand_mode(&run,1,moving,&opts,&feedback,0) || moving[0])failures++;
    if(arr_expand(&run,1,moving,&opts,&feedback)!=1 || !moving[0])failures++;
    fprintf(stderr,"  contact-only context closure observes without publication or seam-only growth: %s\n",failures?"FAIL":"ok");
    Arena_dispose(&arena);return failures;
}

static int arr_anchor_bound_control(void)
{
    Arena_T arena=Arena_new();AsmRun run={0};int failures=0;
    if(ags_fixture(arena,&run)){Arena_dispose(&arena);return 1;}
    AsmField *field=AsmField_new(arena,&run);if(!field){Arena_dispose(&arena);return 1;}
    AsmField local=*field;AsmFieldSeam seams[2]={{0}};AsmWitness witness[2]={{0}};
    local.seams=seams;local.n_seams=2;local.witness=witness;local.n_witness=2;
    uint8_t core[5]={1,0,0,0,0};size_t base=field->vertex_offset[0];
    for(size_t i=0;i<2;i++){
        seams[i].a=0;seams[i].b=(int32_t)i+2;seams[i].first=i;seams[i].count=1;
        seams[i].required=seams[i].complete=seams[i].source_run=1;
        witness[i].valid=1;witness[i].weight=1;
        size_t fixed=field->vertex_offset[i+2];
        field->uv[2*fixed]=i?6:0;field->uv[2*fixed+1]=0;
        for(int d=0;d<2;d++){
            for(int k=0;k<3;k++)witness[i].vertex[d][k]=(int32_t)base;
            witness[i].vertex[d][3]=(int32_t)fixed;witness[i].coefficient[d][0]=1;
        }
    }
    /* Both anchors observe the same original point. Six voxels of target
     * separation permits two residuals of three; the old 2+2 bound wrongly
     * claimed infeasibility. Nine cannot fit two four-voxel residuals. */
    if(apa_incompatible_fixed_anchors(&local,core,NULL))failures++;
    field->uv[2*field->vertex_offset[3]]=9;
    if(!apa_incompatible_fixed_anchors(&local,core,NULL))failures++;
    fprintf(stderr,"  fixed-anchor infeasibility uses both shared four-voxel seam allowances: %s\n",failures?"FAIL":"ok");
    Arena_dispose(&arena);return failures;
}

static int arr_progress_context_control(void)
{
    Arena_T arena=Arena_new();AsmRun run={0};int failures=0;
    if(ags_fixture(arena,&run)){Arena_dispose(&arena);return 1;}
    apr_test_scale(&run,32);run.n_rels=0;
    run.charts[2].pose_x=run.charts[0].pose_x;
    run.charts[3].pose_x=30000;run.charts[3].area3d*=100;
    AsmField *field=AsmField_new(arena,&run);
    if(!field){Arena_dispose(&arena);return 1;}
    uint8_t core[5]={1,0,0,0,0},moving[5]={1,0,0,0,0},context[5]={0};
    double contact[5]={0},seam[5]={0};ApaFeedback feedback={5,contact,seam,context,ves_clock_sec()+10,0};
    apa_feedback(field,core,moving,field->uv,field,NULL,arena,&feedback);
    double initial=feedback.violation,overlap=feedback.material_area;
    if(!feedback.complete || !(overlap>0) || !(initial>0) || !isfinite(initial))failures++;
    /* Widening the movable set to a large, distant chart changes no UV,
     * contact or seam. It cannot reduce the measured remaining defect. */
    moving[3]=1;arr_feedback_reset(&feedback);
    apa_feedback(field,core,moving,field->uv,field,NULL,arena,&feedback);
    if(!feedback.complete || fabs(feedback.violation-initial)>1e-12 ||
       fabs(feedback.material_area-overlap)>1e-9 || arr_remaining(initial,feedback.violation,8,45)!=45)failures++;
    fprintf(stderr,"  unchanged geometry in wider context cannot imitate repair progress: %s\n",failures?"FAIL":"ok");
    Arena_dispose(&arena);return failures;
}

static int arr_resume_control(void)
{
    Arena_T arena=Arena_new();AsmRun run={0};int failures=0;
    if(ags_fixture(arena,&run)){Arena_dispose(&arena);return 1;}
    AsmRun trial=run;AsmChart placed[5];memcpy(placed,run.charts,sizeof placed);trial.charts=placed;
    placed[1].placed=1;placed[1].placement_state=ASM_PLACE_SEAM;
    AsmField *field=AsmField_new(arena,&trial);if(!field){Arena_dispose(&arena);return 1;}
    size_t bytes=2*field->nv*sizeof(double);double *reference=ARENA_ALLOC(arena,bytes);memcpy(reference,field->uv,bytes);
    uint8_t moving[5]={1,1,0,0,0},neighbors[5]={1,0,0,0,0},context[5]={1,1,0,0,0};
    for(size_t v=0;v<run.charts[1].nv;v++)field->uv[2*(field->vertex_offset[1]+v)]+=.25;
    ApaContinuation states[5]={{0}},loaded[5]={{0}};
    apa_continuation_save(states+1,field,moving,reference);
    uint8_t *pending[5]={0},*dependencies[5]={0},*loaded_pending[5]={0},*loaded_dependencies[5]={0};
    pending[1]=neighbors;dependencies[1]=context;size_t rounds[5]={0,3,0,0,0},loaded_rounds[5]={0};
    double forecast[5]={0,2.5,0,0,0},loaded_forecast[5]={0};ArrStats stats={0};
    AsmRepairOpts policy;AsmRepair_defaults(&policy);AsmChart before[5];memcpy(before,run.charts,sizeof before);
    if(arr_resume_write(&run,states,pending,dependencies,rounds,forecast,&stats) || stats.resume_saved!=1 ||
       memcmp(before,run.charts,sizeof before))failures++;
    char path[512];snprintf(path,sizeof path,"output/asm_resume_%d.asr",ves_getpid());
    AsmRun copy={0};copy.arena=arena;
    if(AsmStore_write_run(path,&run) || AsmStore_read_run(path,&copy) ||
       arr_resume_read(&copy,&policy,arena,loaded,loaded_pending,loaded_dependencies,loaded_rounds,loaded_forecast,&stats) ||
       stats.resume_loaded!=1 || !loaded_pending[1] || !loaded_pending[1][0] || loaded_rounds[1]!=3 ||
       loaded_forecast[1]!=2.5 || AsmChart_registered(copy.charts+1) ||
       !loaded[1].uv[1] || memcmp(loaded[1].uv[1],states[1].uv[1],2*run.charts[1].nv*sizeof(double)))failures++;
    if(loaded[1].arena){
        /* Source-compatible cached UVs are still speculative. This rigid
         * but badly separated seed must fail its incident seam, without
         * changing the live checkpoint or gaining placement evidence. */
        AsmChart original[5];memcpy(original,copy.charts,sizeof original);
        AsmField *incumbent=AsmField_new(arena,&copy);ApaBox boxes[5];apa_chart_boxes(&copy,boxes);
        AsmRepairOpts opts=policy;opts.compact=1;
        opts.iterations=opts.clearance_iterations=opts.closure_iterations=0;opts.admission_rounds=1;
        opts.admission_continuation=loaded+1;opts.admission_moving=neighbors;
        uint8_t core[5]={0,1,0,0,0};AsmRepairStats part={0};ApaReport report={0};
        size_t anchors=0;const char *reason="";
        if(!incumbent || apa_attempt_anchor(&copy,&opts,core,&anchors,&part,&reason,NULL,boxes,incumbent,NULL,&report,4,SIZE_MAX) ||
           strcmp(reason,"incident_source_seam_failed") || memcmp(original,copy.charts,sizeof original))failures++;
    }
    for(size_t c=0;c<5;c++)apa_continuation_clear(loaded+c);
    /* Preserve a guess across a known field commit, including a new source
     * obligation. A foreign source/neighbor change or corrupt scratch blob
     * before capture must still discard it. Neither route places the core. */
    for(int mode=0;mode<5;mode++){
        int before_failures=failures;
        AsmRun next={0};next.arena=arena;
        if(AsmStore_read_run(path,&next)){failures++;continue;}
        if(mode==1){
            /* Begin with this source hypothesis not yet promoted, then
             * model its validated promotion after the guess is held. */
            next.rels[0].continuity&=~ASM_CONT_REQUIRED;ArrStats encoded={0};
            if(arr_resume_write(&next,states,pending,dependencies,rounds,forecast,&encoded)){failures++;continue;}
        }
        if(mode==2)next.charts[1].xyz[0]+=1;
        if(mode==3)next.charts[0].pose_x+=1;
        if(mode==4)next.repair_resume_size--;
        ArrHeld held={0};arr_hold(&next,&policy,&held);
        AsmField *live=AsmField_new(arena,&next);
        if(!live){failures++;arr_held_clear(&held);continue;}
        if(mode==1){
            for(size_t v=0;v<next.charts[0].nv;v++)live->uv[2*(live->vertex_offset[0]+v)]+=.125;
            next.rels[0].continuity|=ASM_CONT_REQUIRED;
        }
        if(AsmField_commit(live,&next) || arr_refresh(&next,&policy,&held))failures++;
        arr_held_clear(&held);memset(&stats,0,sizeof stats);
        if(arr_resume_read(&next,&policy,arena,loaded,loaded_pending,loaded_dependencies,loaded_rounds,loaded_forecast,&stats) ||
           stats.resume_loaded!=(size_t)(mode<2) || AsmChart_registered(next.charts+1))failures++;
        if(mode<2 && (!loaded[1].uv[1] || loaded_rounds[1]!=3 ||
           loaded_forecast[1]!=(mode?policy.patch_seconds:2.5) ||
           memcmp(loaded[1].uv[1],states[1].uv[1],2*run.charts[1].nv*sizeof(double))))failures++;
        if(mode==1 && loaded[1].arena){
            AsmChart original[5];memcpy(original,next.charts,sizeof original);ApaBox boxes[5];apa_chart_boxes(&next,boxes);
            AsmRepairOpts opts=policy;opts.compact=1;opts.admission_moving=neighbors;opts.admission_continuation=loaded+1;
            opts.iterations=opts.clearance_iterations=opts.closure_iterations=0;opts.admission_rounds=1;
            uint8_t core[5]={0,1,0,0,0};AsmRepairStats part={0};ApaReport report={0};size_t anchors=0;const char *reason="";
            if(apa_attempt_anchor(&next,&opts,core,&anchors,&part,&reason,NULL,boxes,live,NULL,&report,4,SIZE_MAX) ||
               strcmp(reason,"incident_source_seam_failed") || memcmp(original,next.charts,sizeof original))failures++;
        }
        for(size_t c=0;c<5;c++)apa_continuation_clear(loaded+c);
        fprintf(stderr,"  private guess finalization %d: %s\n",mode,failures==before_failures?"ok":"FAIL");
    }
    /* Independently change every kind of dependency. Unrelated placement
     * does not invalidate a candidate; source, parity, witnesses and fixed
     * neighbor motion do. Too-small current caps also refuse reuse. */
    for(int mode=0;mode<7;mode++){
        float xyz=run.charts[1].xyz[0];double pose=run.charts[0].pose_x;uint32_t flags=run.charts[1].flags;
        float gap=run.corr[0].gap_a[0];double distant=run.charts[2].pose_x;size_t cap=policy.patch_coordinates;
        double *encoded=run.charts[0].placed_uv;
        if(mode==0)run.charts[1].xyz[0]+=1;
        if(mode==1)run.charts[0].pose_x+=1;
        if(mode==2)run.charts[1].flags^=ASM_CHART_MIRROR;
        if(mode==3)run.corr[0].gap_a[0]+=1;
        if(mode==4)policy.patch_coordinates=2;
        if(mode==5)run.charts[2].pose_x+=1;
        if(mode==6){
            double *uv=ARENA_ALLOC(arena,2*run.charts[0].nv*sizeof(double));
            for(size_t v=0;v<run.charts[0].nv;v++)AsmChart_point(run.charts,v,uv+2*v);
            run.charts[0].placed_uv=uv;
        }
        memset(&stats,0,sizeof stats);
        if(arr_resume_read(&run,&policy,arena,loaded,loaded_pending,loaded_dependencies,loaded_rounds,loaded_forecast,&stats) ||
           stats.resume_loaded!=(size_t)(mode>=5) || stats.resume_stale!=(size_t)(mode<5))failures++;
        for(size_t c=0;c<5;c++)apa_continuation_clear(loaded+c);
        run.charts[1].xyz[0]=xyz;run.charts[0].pose_x=pose;run.charts[1].flags=flags;
        run.corr[0].gap_a[0]=gap;run.charts[2].pose_x=distant;policy.patch_coordinates=cap;
        run.charts[0].placed_uv=encoded;
    }
    /* The first frozen resume experiment used representation-sensitive
     * fingerprints. Read those records with their recorded schema. */
    memset(&stats,0,sizeof stats);forecast[1]=1000000;
    if(arr_resume_encode(&run,states,pending,dependencies,rounds,forecast,&stats,1) ||
       arr_resume_read(&run,&policy,arena,loaded,loaded_pending,loaded_dependencies,loaded_rounds,loaded_forecast,&stats) ||
       stats.resume_loaded!=1 || loaded_forecast[1]!=policy.patch_seconds)failures++;
    for(size_t c=0;c<5;c++)apa_continuation_clear(loaded+c);
    if(arr_resume_write(&run,states,pending,dependencies,rounds,forecast,&stats))failures++;
    memcpy(run.repair_resume,"ARRUV002",8);memset(&stats,0,sizeof stats);
    if(arr_resume_read(&run,&policy,arena,loaded,loaded_pending,loaded_dependencies,loaded_rounds,loaded_forecast,&stats) ||
       stats.resume_loaded!=1 || loaded_forecast[1]!=policy.patch_seconds)failures++;
    for(size_t c=0;c<5;c++)apa_continuation_clear(loaded+c);
    memcpy(run.repair_resume,"ARRUV003",8);memset(&stats,0,sizeof stats);
    if(arr_resume_read(&run,&policy,arena,loaded,loaded_pending,loaded_dependencies,loaded_rounds,loaded_forecast,&stats) ||
       stats.resume_loaded!=1 || loaded_forecast[1]!=policy.patch_seconds)failures++;
    for(size_t c=0;c<5;c++)apa_continuation_clear(loaded+c);
    memcpy(run.repair_resume,"ARRUV004",8);memset(&stats,0,sizeof stats);
    if(arr_resume_read(&run,&policy,arena,loaded,loaded_pending,loaded_dependencies,loaded_rounds,loaded_forecast,&stats) ||
       stats.resume_loaded!=1 || loaded_forecast[1]!=1000000)failures++;
    for(size_t c=0;c<5;c++)apa_continuation_clear(loaded+c);
    forecast[1]=2.5;memset(&stats,0,sizeof stats);
    if(arr_resume_encode(&run,states,pending,dependencies,rounds,forecast,&stats,1) ||
       arr_resume_read(&run,&policy,arena,loaded,loaded_pending,loaded_dependencies,loaded_rounds,loaded_forecast,&stats) ||
       stats.resume_loaded!=1 || loaded_forecast[1]!=2.5)failures++;
    for(size_t c=0;c<5;c++)apa_continuation_clear(loaded+c);
    size_t length=run.repair_resume_size;run.repair_resume_size--;
    memset(&stats,0,sizeof stats);
    if(!arr_resume_read(&run,&policy,arena,loaded,loaded_pending,loaded_dependencies,loaded_rounds,loaded_forecast,&stats) ||
       stats.resume_loaded || !stats.resume_invalid || loaded[1].arena)failures++;
    run.repair_resume_size=length;uint8_t *blob=run.repair_resume;
    uint64_t bad=UINT64_MAX;memcpy(blob+24,&bad,sizeof bad);memset(&stats,0,sizeof stats);
    if(!arr_resume_read(&run,&policy,arena,loaded,loaded_pending,loaded_dependencies,loaded_rounds,loaded_forecast,&stats) ||
       stats.resume_loaded || loaded[1].arena)failures++;
    for(size_t c=0;c<5;c++){apa_continuation_clear(states+c);apa_continuation_clear(loaded+c);}
    remove(path);fprintf(stderr,"  checkpointed private UVs, current-context validation, caps and corrupt-data rollback: %s\n",failures?"FAIL":"ok");
    Arena_dispose(&arena);return failures;
}

/* Holes before edges, half the time: a hole is enclosed, surveyed and requeued. */
static int arr_hole_control(void)
{
    int failures=0;
    failures+=!arr_is_waiting_hole(3,0,1,3);    /* surveyed, requeued, three placed partners */
    failures+=!arr_is_waiting_hole(6,0,4,3);    /* a continuation */
    failures+=arr_is_waiting_hole(2,0,1,3);     /* an edge target */
    failures+=arr_is_waiting_hole(3,0,0,3);     /* not surveyed yet: the ordinary order surveys it */
    failures+=arr_is_waiting_hole(3,1,2,3);     /* tried and not requeued */
    if(failures)fprintf(stderr,"  region refine hole-turn control: %d failures\n",failures);
    return failures;
}

/* Gateway turn (ASM_REGION_REFINE_GATEWAY_SHARE): a requeued target is worth
 * the missing branch it opens per forecast second, nothing when that branch is
 * itself; the turn takes the best such target while the ordinary order still
 * surveys first, and a single-chart hole never wins it. */
static int arr_gateway_control(void)
{
    int failures=0;
    failures+=!(arr_gateway_rate(192,64,45)>0);                                   /* area behind it */
    failures+=arr_gateway_rate(64,64,45)!=0;                                      /* its branch is itself */
    failures+=!(arr_gateway_rate(192,64,10)>arr_gateway_rate(192,64,45));        /* sooner is worth more */
    failures+=!(arr_gateway_rate(192,64,INFINITY)>0);                             /* no forecast yet */
    /* 0 and 2 placed; 1 opens 5 and 6; 3 is a requeued single chart; 4 is unsurveyed */
    enum {NC=7};
    Arena_T arena=Arena_new();AsmRun run={0};run.arena=arena;run.n_charts=NC;
    run.charts=ARENA_CALLOC(arena,NC,sizeof *run.charts);float uv[8]={0};
    for(int c=0;c<NC;c++){
        AsmChart *ch=run.charts+c;ch->id=c;ch->parent=-1;ch->nv=4;ch->nf=2;ch->area3d=64;ch->uv=uv;
        if(c==0 || c==2){ch->placed=1;ch->placement_state=ASM_PLACE_ROOT;}
    }
    const int pairs[5][2]={{0,1},{2,3},{2,4},{1,5},{5,6}};
    for(int i=0;i<5;i++){
        AsmRelation r={0};r.a=pairs[i][0];r.b=pairs[i][1];r.continuity=ASM_CONT_SOURCE|ASM_CONT_REQUIRED;
        AsmRun_push_rel(&run,&r);
    }
    size_t group[NC],anchors[NC],rounds[NC]={0,1,0,1,0,0,0};
    double area[NC],spent[NC]={0},forecast[NC],potential=0;uint8_t tried[NC]={0},scratch[NC];
    for(int c=0;c<NC;c++)forecast[c]=45;
    AsmRepairOpts policy;AsmRepair_defaults(&policy);
    failures+=arr_next(&run,tried,group,anchors,area,1,&potential,&policy,scratch,spent,NULL,rounds,NULL,forecast,0,0)!=4;
    failures+=arr_next(&run,tried,group,anchors,area,1,&potential,&policy,scratch,spent,NULL,rounds,NULL,forecast,0,1)!=1 ||
              fabs(potential-192)>1e-9;
    /* An unsurveyed gateway wins the turn too (its survey is its first
     * quantum), while the ordinary order takes the larger unsurveyed chart. */
    rounds[1]=0;run.charts[4].area3d=100;
    failures+=arr_next(&run,tried,group,anchors,area,1,&potential,&policy,scratch,spent,NULL,rounds,NULL,forecast,0,0)!=4;
    failures+=arr_next(&run,tried,group,anchors,area,1,&potential,&policy,scratch,spent,NULL,rounds,NULL,forecast,0,1)!=1;
    tried[1]=1;   /* the gateway was taken: a turn without a gateway keeps the ordinary pick */
    failures+=arr_next(&run,tried,group,anchors,area,1,&potential,&policy,scratch,spent,NULL,rounds,NULL,forecast,0,1)!=4;
    Arena_dispose(&arena);
    if(failures)fprintf(stderr,"  region refine gateway-turn control: %d failures\n",failures);
    return failures;
}

/* A target outlives the patch budget only while its last joint quantum still
 * cuts the violation; survey quanta, stalls and unmeasured states never do. */
static int arr_converging_control(void)
{
    int failures=0;
    failures+=!arr_converging(1,10.0,5.0);                          /* halved: converging */
    failures+=arr_converging(1,10.0,9.8);                           /* 2% cut: stalled */
    failures+=arr_converging(1,10.0,12.0);                          /* rising */
    failures+=arr_converging(0,10.0,1.0);                           /* a survey quantum never counts */
    failures+=arr_converging(1,INFINITY,1.0) || arr_converging(1,10.0,INFINITY) || arr_converging(1,0.0,0.0);
    failures+=fabs(arr_target_limit(0,45)-45)>1e-12;
    failures+=fabs(arr_target_limit(1,45)-fmax(45.0,ASM_REGION_REFINE_CONVERGING_SEC))>1e-12;
    failures+=fabs(arr_target_limit(1,400)-fmax(400.0,ASM_REGION_REFINE_CONVERGING_SEC))>1e-12;   /* never below the budget */
    if(failures)fprintf(stderr,"  region refine converging-budget control: %d failures\n",failures);
    return failures;
}

/* Refinement may add a placement but never edits a chart's source inventory. */
static int arr_inventory_same(const AsmChart *a,const AsmChart *b)
{
    return a->id==b->id && a->cube==b->cube && a->comp==b->comp && a->parent==b->parent &&
        a->nv==b->nv && a->nf==b->nf && a->xyz==b->xyz && a->faces==b->faces &&
        a->vid==b->vid && a->uv==b->uv && a->nrm==b->nrm && a->boundary==b->boundary &&
        a->repaired==b->repaired && a->area3d==b->area3d;
}

static int arr_selftest(void)
{
    Arena_T arena=Arena_new();AsmRun run={0};int failures=0;char root[256],path[512];
    if(ags_fixture(arena,&run)){Arena_dispose(&arena);return 1;}
    size_t group[5],anchors[5];double area[5],potential=0;uint8_t tried[5]={0},moving[5]={0};
    AsmRepairOpts policy;AsmRepair_defaults(&policy);
    if(arr_next(&run,tried,group,anchors,area,100,&potential,&policy,moving,NULL,NULL,NULL,NULL,NULL,0,0)!=1 || fabs(potential-128)>1e-12)failures++;
    /* A hole turn without a waiting hole keeps the ordinary pick. */
    {size_t survey_rounds[5]={0};
     if(arr_next(&run,tried,group,anchors,area,100,&potential,&policy,moving,NULL,NULL,survey_rounds,NULL,NULL,1,0)!=1)failures++;}
    tried[1]=1;if(arr_next(&run,tried,group,anchors,area,100,&potential,&policy,moving,NULL,NULL,NULL,NULL,NULL,0,0)!=SIZE_MAX)failures++;tried[1]=0;
    run.charts[4].flags|=ASM_CHART_EXTRAS;
    if(arr_next(&run,tried,group,anchors,area,100,&potential,&policy,moving,NULL,NULL,NULL,NULL,NULL,0,0)!=SIZE_MAX)failures++;
    /* A 64-vox² attached chart must remain repairable even when an output
     * piece has a 20,000-vox² floor. Excluded successors stay excluded. */
    if(arr_next(&run,tried,group,anchors,area,ASM_REGION_REFINE_MIN_BRANCH_AREA,&potential,&policy,moving,NULL,NULL,NULL,NULL,NULL,0,0)!=1 || fabs(potential-64)>1e-12)failures++;
    run.charts[4].flags&=~ASM_CHART_EXTRAS;
    if(arr_neighbors(&run,1,moving)!=1 || !moving[0] || moving[1] || moving[2] || moving[3] || moving[4])failures++;
    AsmRelation saved=run.rels[0];run.rels[0].flags|=ASM_REL_WEAK;
    run.rels[0].continuity&=~(ASM_CONT_REQUIRED|ASM_CONT_CUT);
    if(arr_next(&run,tried,group,anchors,area,32,&potential,&policy,moving,NULL,NULL,NULL,NULL,NULL,0,0)!=SIZE_MAX)failures++;
    run.rels[0]=saved;
    AsmChart original[5];memcpy(original,run.charts,sizeof original);
    snprintf(root,sizeof root,"output/asm_region_refine_%d",ves_getpid());ves_mkdir("output");ves_mkdir(root);
    snprintf(path,sizeof path,"%s/expired",root);ArrStats stats={0};
    if(arr_run(&run,path,ves_clock_sec()-1,32,&stats) || stats.attempts || !stats.expired ||
       memcmp(original,run.charts,sizeof original))failures++;
    if(arr_run(&run,path,ves_clock_sec()+1,0,&stats)!=-1 || arr_run(NULL,path,ves_clock_sec()+1,32,&stats)!=-1)failures++;
    AsmField_threads(1);snprintf(path,sizeof path,"%s/active",root);
    if(arr_run(&run,path,ves_clock_sec()+15,32,&stats) || !stats.admitted || !AsmChart_registered(run.charts+1))failures++;
    for(size_t c=0;c<5;c++)if(!arr_inventory_same(original+c,run.charts+c) ||
        (AsmChart_registered(original+c) && !AsmChart_registered(run.charts+c)))failures++;
    AsmField *field=AsmField_new(arena,&run);AsmContactStats contacts={0};
    if(!field || AsmContacts_measure(field->contacts,field->uv,NULL,NULL,NULL,&contacts) || contacts.pairs)failures++;
    else{
        for(size_t c=0;c<5;c++)if(AsmChart_registered(run.charts+c) && !apa_chart_metric(field,c))failures++;
        for(size_t i=0;i<field->n_seams;i++)if(field->seams[i].required && apa_present(&run,field->seams[i].a,field->seams[i].b)){
            AsmSeamMeasure measure;AsmField_seam(field,i,field->uv,&measure);if(!measure.pass)failures++;
        }
    }
    snprintf(path,sizeof path,"%s/active/proposal.asr",root);FILE *unexpected=fopen(path,"rb");if(unexpected){fclose(unexpected);failures++;}
    failures+=arr_converging_control()+arr_hole_control()+arr_gateway_control()+arr_margin_control()+arr_policy_control()+arr_obstacle_feedback_control()+arr_continuation_control()+
        arr_context_survey_control()+arr_anchor_bound_control()+arr_progress_context_control()+arr_resume_control()+
        apa_contact_certificate_control()+apa_residual_trace_control()+apa_residual_nudge_control()+apa_contact_density_control()+apa_inherited_self_contact_control()+apa_joint_core_polish_control()+
        apa_polish_density_control()+apa_polish_interface_control();AsmField_threads(24);Arena_dispose(&arena);
    fprintf(stderr,"  native required-region refinement: %s (%d failures)\n",failures?"FAIL":"ok",failures);return failures;
}
