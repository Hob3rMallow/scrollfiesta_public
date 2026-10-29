/* A missing third chart is an unresolved audit obligation, not an infinite
 * contribution to the progress of a different, measurable local seam. */
static int apr_incident_stop_control(void)
{
    int failures=0;
    for(int defect=0;defect<3;defect++){
        Arena_T arena=Arena_new();AsmRun run={0};AsmRepairStats stats={0};
        if(ags_grid_fixture(arena,&run)){Arena_dispose(&arena);return 1;}
        for(size_t c=0;c<run.n_charts;c++){
            run.charts[c].placed=1;run.charts[c].placement_state=ASM_PLACE_ROOT;
            run.charts[c].pose_x=20*(double)(c%3);run.charts[c].pose_y=20*(double)(c/3);
        }
        if(defect==1)run.charts[4].pose_x+=6;
        if(defect==2)for(size_t v=0;v<run.charts[4].nv;v++)run.charts[4].uv[2*v]*=1.4f;
        AsmField *field=AsmField_new(arena,&run);
        if(!field){failures++;Arena_dispose(&arena);continue;}
        size_t bytes=2*field->nv*sizeof(double);double *before=ARENA_ALLOC(arena,bytes);memcpy(before,field->uv,bytes);
        AsmRepairOpts opts;AsmRepair_defaults(&opts);opts.admission_incident_sources=1;
        opts.iterations=2;opts.clearance_iterations=opts.closure_iterations=0;
        if(apr_solve(field,&opts,&stats))failures++;
        if(!defect){
            AsmContactStats contacts={0};
            if(stats.iterations || memcmp(before,field->uv,bytes) ||
               AsmContacts_measure(field->contacts,field->uv,NULL,NULL,NULL,&contacts) || !contacts.complete || contacts.pairs)failures++;
            for(size_t c=0;c<run.n_charts;c++)if(!apa_chart_metric(field,c))failures++;
            for(size_t i=0;i<field->n_seams;i++){
                AsmSeamMeasure seam;AsmField_seam(field,i,field->uv,&seam);if(!seam.pass)failures++;
            }
        }else if(!stats.iterations)failures++;
        Arena_dispose(&arena);
    }
    fprintf(stderr,"  incident-source feasible stop and defective-field continuation: %s\n",failures?"FAIL":"ok");
    return failures;
}

static int apr_rank_control(void)
{
    Arena_T arena=Arena_new();AsmRun run={0};int failures=apr_incident_stop_control();
    if(ags_fixture(arena,&run)){Arena_dispose(&arena);return 1;}
    const double x[5]={0,30,23,53,1000};
    for(size_t c=0;c<5;c++){
        run.charts[c].placed=c!=4;run.charts[c].placement_state=c!=4?ASM_PLACE_ROOT:ASM_PLACE_NONE;
        run.charts[c].pose_x=x[c];run.charts[c].pose_y=0;
    }
    for(size_t i=0;i<run.n_rels;i++)run.rels[i].continuity|=ASM_CONT_REQUIRED;
    uint8_t member[5]={1,1,1,1,0},fixed[5]={1,0,1,1,1},moving[5]={0,1,0,0,0};
    AsmField *f=AsmField_new_local(arena,&run,member,NULL);
    if(!f || AsmField_build_contacts(f)){Arena_dispose(&arena);return 1;}
    size_t bytes=2*f->nv*sizeof(double),nc=run.n_charts;
    double *reference=ARENA_ALLOC(arena,bytes);memcpy(reference,f->uv,bytes);
    AprSolver s={0};AsmRepairStats stats={0};s.field=f;s.stats=&stats;AsmRepair_defaults(&s.opts);s.opts.fixed_charts=fixed;
    s.lo=ARENA_ALLOC(arena,f->nf*sizeof(double));s.hi=ARENA_ALLOC(arena,f->nf*sizeof(double));
    s.within10=ARENA_ALLOC(arena,nc*sizeof(double));s.within25=ARENA_ALLOC(arena,nc*sizeof(double));s.area=ARENA_ALLOC(arena,nc*sizeof(double));
    s.fixed=ARENA_CALLOC(arena,f->nv,1);s.seam=ARENA_ALLOC(arena,f->n_seams*sizeof *s.seam);
    for(size_t c=0;c<nc;c++)if(fixed[c] && f->vertex_offset[c]!=SIZE_MAX)memset(s.fixed+f->vertex_offset[c],1,run.charts[c].nv);
    s.best_uv=ARENA_ALLOC(arena,bytes);memcpy(s.best_uv,reference,bytes);
    if(AsmContacts_measure(f->contacts,f->uv,NULL,NULL,NULL,&s.contacts) || !s.contacts.complete){failures++;goto done;}
    AprRank before=apr_rank(&s);s.best_rank=before;
    size_t active=SIZE_MAX,absent=SIZE_MAX,stationary=SIZE_MAX;
    for(size_t i=0;i<f->n_seams;i++){
        if(f->seams[i].a==0)active=i;
        else if(f->seams[i].a==2)stationary=i;
        else absent=i;
    }
    int okay=f->n_seams==3 && active!=SIZE_MAX && absent!=SIZE_MAX && stationary!=SIZE_MAX &&
        before.bad_seams==1 && before.unavailable_seams==1 && before.fixed_seams==1 && isfinite(before.seam_failure);
    if(!okay){failures++;goto done;}
    for(size_t v=0;v<run.charts[1].nv;v++)f->uv[2*(f->vertex_offset[1]+v)]-=4;
    if(AsmContacts_measure(f->contacts,f->uv,NULL,NULL,NULL,&s.contacts) || !s.contacts.complete){failures++;goto done;}
    AprRank after=apr_rank(&s);
    okay=after.bad_seams==1 && after.seam_failure<before.seam_failure && after.area>before.area && apr_better(after,before);
    apr_keep_best(&s);okay=okay && !memcmp(s.best_uv,f->uv,bytes) && memcmp(s.best_uv,reference,bytes);
    AsmSeamMeasure missing;AsmField_seam(f,absent,f->uv,&missing);
    okay=okay && f->seams[absent].required && missing.invalid && !missing.pass && f->n_seams==3;
    int complete=0;const char *reason=apa_qualification(f,moving,moving,reference,NULL,NULL,4,&complete);
    okay=okay && reason && !strcmp(reason,"incident_source_seam_failed");
    /* A corrupt moving observation remains a real failure. The same defect
     * in wholly fixed context cannot erase another seam's finite progress. */
    size_t wa=f->seams[active].first,wf=f->seams[stationary].first;
    int saved=f->witness[wf].valid;f->witness[wf].valid=0;
    AprRank fixed_invalid=apr_rank(&s);okay=okay && isfinite(fixed_invalid.seam_failure) && fixed_invalid.bad_seams==1;
    f->witness[wf].valid=saved;saved=f->witness[wa].valid;f->witness[wa].valid=0;
    AprRank moving_invalid=apr_rank(&s);okay=okay && !isfinite(moving_invalid.seam_failure) && moving_invalid.bad_seams==1;
    f->witness[wa].valid=saved;s.opts.fixed_charts=NULL;
    AprRank whole=apr_rank(&s);okay=okay && whole.bad_seams==2 && !whole.fixed_seams && whole.unavailable_seams==1;
    /* The hopeless-start seam bound measures only present required seams:
     * the near start is below it, the same chart 60 vox off is above it. */
    double near_start=apa_worst_start_seam(f,moving,4);
    for(size_t v=0;v<run.charts[1].nv;v++)f->uv[2*(f->vertex_offset[1]+v)]+=60;
    double far_start=apa_worst_start_seam(f,moving,4);
    for(size_t v=0;v<run.charts[1].nv;v++)f->uv[2*(f->vertex_offset[1]+v)]-=60;
    okay=okay && isfinite(near_start) && near_start<ASM_ADMISSION_MAX_START_SEAM_RMS && far_start>ASM_ADMISSION_MAX_START_SEAM_RMS;
    if(!okay)fprintf(stderr,"  hopeless seam start: near %.9g, far %.9g vox\n",near_start,far_start);
    if(!okay)fprintf(stderr,"  seam rank: bad %zu -> %zu, deficit %.9g -> %.9g, contact %.9g -> %.9g, qualification %s\n",
        before.bad_seams,after.bad_seams,before.seam_failure,after.seam_failure,before.area,after.area,reason?reason:"accepted");
    failures+=!okay;
done:
    fprintf(stderr,"  local seam progress with absent/fixed obligations and unchanged admission: %s\n",failures?"FAIL":"ok");
    Arena_dispose(&arena);return failures;
}
