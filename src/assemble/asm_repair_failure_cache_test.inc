static int afc_selftest(void)
{
    int failures=0;Arena_T arena=Arena_new();AsmRun run={0};
    if(ags_fixture(arena,&run)){failures++;goto done;}
    run.n_charts=3;run.n_rels=1;run.n_corr=9;
    AsmChart *obstacle=run.charts+2;obstacle->nv=4;obstacle->nf=2;obstacle->pose_x=116;
    int32_t quad[]={0,1,2,0,2,3};float uv[]={0,0,8,0,8,8,0,8};
    obstacle->faces=ARENA_ALLOC(arena,sizeof quad);memcpy(obstacle->faces,quad,sizeof quad);
    memcpy(obstacle->uv,uv,sizeof uv);
    /* The rigid source fit keeps the 2-unit 3-D gap to chart 0, so the
     * candidate lands on uv [110,118] (z = u - 100) and overlaps the obstacle
     * on [116,118].  The obstacle repeats the candidate's sheet there (z =
     * u - 100 as well): a same-place start contact that the pre-solve
     * shortcut leaves to the solve under either contact rule, so the cache
     * sees the solver's outcome. */
    for(size_t v=0;v<4;v++){obstacle->xyz[3*v]=16+uv[2*v];obstacle->xyz[3*v+1]=uv[2*v+1];obstacle->xyz[3*v+2]=0;}
    apr_test_scale(&run,8);
    AsmField *parent=AsmField_new(arena,&run);ApaBox boxes[3];apa_chart_boxes(&run,boxes);
    if(!parent || AsmContacts_refit(parent->contacts,parent->uv)){failures++;goto done;}
    uint8_t core[]={0,1,0};AsmRepairOpts opts;AsmRepair_defaults(&opts);
    opts.compact=1;opts.iterations=opts.clearance_iterations=opts.closure_iterations=0;opts.admission_rounds=1;
    opts.regional=opts.source_trials=0;opts.deadline=ves_clock_sec()+10;
    ApaFailureCache cache;afc_init(&cache,arena,16);opts.admission_failure_cache=&cache;
    AsmChart before[3];memcpy(before,run.charts,sizeof before);
    for(int repeat=0;repeat<2;repeat++){
        size_t anchors=0;AsmRepairStats stats={0};ApaReport report={0};const char *reason="";
        int accepted=apa_attempt(&run,&opts,core,&anchors,&stats,&reason,NULL,boxes,parent,NULL,&report,4);
        /* The physical polisher now works until its short deadline; that
         * unfinished work must remain retryable instead of being cached. */
        int okay=!accepted && !strcmp(reason,"new_material_contact_failed") && !report.cache_hit &&
            !memcmp(before,run.charts,sizeof before) && !cache.count && !cache.hits;
        if(!okay)fprintf(stderr,"  cache repeat %d: accepted %d reason %s hit %d entries %zu hits %zu\n",repeat,accepted,reason,report.cache_hit,cache.count,cache.hits);
        failures+=!okay;
    }
    /* A complete failed seam still caches. Conflicting recorded boundary
     * observations cannot all fit in this zero-iteration pose trial.  The
     * perturbed fit starts the candidate more than the same-surface distance
     * off, where the obstacle would be refused before the solve (the area and
     * 3-D shortcut) and never reach the seam certificate this case caches, so
     * the obstacle stands aside for it and returns afterwards. */
    AsmCorr saved_corr=run.corr[4];run.corr[4].gap_a[1]+=32;
    double obstacle_x=obstacle->pose_x;obstacle->pose_x=8000;apa_chart_boxes(&run,boxes);
    parent=AsmField_new(arena,&run);if(!parent || AsmContacts_refit(parent->contacts,parent->uv)){failures++;goto done;}
    AsmChart aside[3];memcpy(aside,run.charts,sizeof aside);
    for(int repeat=0;repeat<2;repeat++){
        size_t anchors=0;AsmRepairStats stats={0};ApaReport report={0};const char *reason="";
        int accepted=apa_attempt(&run,&opts,core,&anchors,&stats,&reason,NULL,boxes,parent,NULL,&report,4);
        int okay=!accepted && !strcmp(reason,"incident_source_seam_failed") && report.cache_hit==repeat &&
            cache.count==1 && cache.hits==(size_t)repeat && !memcmp(aside,run.charts,sizeof aside);
        if(!okay)fprintf(stderr,"  seam cache repeat %d: accepted %d reason %s hit %d entries %zu hits %zu\n",repeat,accepted,reason,report.cache_hit,cache.count,cache.hits);
        failures+=!okay;
    }
    run.corr[4]=saved_corr;obstacle->pose_x=obstacle_x;apa_chart_boxes(&run,boxes);
    parent=AsmField_new(arena,&run);if(!parent || AsmContacts_refit(parent->contacts,parent->uv)){failures++;goto done;}
    size_t probes=cache.probes;opts.deadline=ves_clock_sec()-1;
    size_t anchors=0;AsmRepairStats stats={0};ApaReport report={0};const char *reason="";
    if(apa_attempt(&run,&opts,core,&anchors,&stats,&reason,NULL,boxes,parent,NULL,&report,4) || cache.probes!=probes || strcmp(reason,"shared_deadline"))failures++;
    opts.deadline=ves_clock_sec()+10;
    /* Same candidate after its obstacle moves away must be solved afresh. */
    obstacle->pose_x=8000;apa_chart_boxes(&run,boxes);parent=AsmField_new(arena,&run);
    if(!parent || AsmContacts_refit(parent->contacts,parent->uv) ||
       !apa_attempt(&run,&opts,core,&anchors,&stats,&reason,NULL,boxes,parent,NULL,&report,4) ||
       report.cache_hit || cache.hits!=1 || !AsmChart_registered(run.charts+1))failures++;
    AsmField *field=AsmField_new(arena,&run);
    if(!field){failures++;goto done;}
    size_t bytes=2*field->nv*sizeof(double);double *reference=ARENA_ALLOC(arena,bytes);memcpy(reference,field->uv,bytes);
    char key[65],changed[65];afc_key(field,&opts,core,reference,4,key);
    for(int variant=0;variant<13;variant++){
        AsmRepairOpts changed_opts=opts;AsmFieldSeam *seam=field->seams;AsmWitness *w=field->witness;
        if(variant==0)field->uv[0]+=.125;
        if(variant==1)reference[0]+=.125;
        if(variant==2)field->xyz[0]+=.125f;
        if(variant==3){int32_t t=field->faces[0];field->faces[0]=field->faces[1];field->faces[1]=t;}
        if(variant==4)run.charts[0].uv[0]+=.125f;
        if(variant==5)run.charts[0].boundary[0]^=1;
        if(variant==6)run.charts[0].flags^=ASM_CHART_MIRROR;
        if(variant==7)seam->required^=1;
        if(variant==8)w->coefficient[0][0]+=.125;
        if(variant==9)changed_opts.iterations++;
        if(variant==10)changed_opts.clearance_method=(opts.clearance_method+1)%4;
        if(variant==11)changed_opts.seam_weight*=2;
        if(variant==12)core[0]^=1;
        afc_key(field,&changed_opts,core,reference,4,changed);failures+=!strcmp(key,changed);
        if(variant==0)field->uv[0]-=.125;
        if(variant==1)reference[0]-=.125;
        if(variant==2)field->xyz[0]-=.125f;
        if(variant==3){int32_t t=field->faces[0];field->faces[0]=field->faces[1];field->faces[1]=t;}
        if(variant==4)run.charts[0].uv[0]-=.125f;
        if(variant==5)run.charts[0].boundary[0]^=1;
        if(variant==6)run.charts[0].flags^=ASM_CHART_MIRROR;
        if(variant==7)seam->required^=1;
        if(variant==8)w->coefficient[0][0]-=.125;
        if(variant==12)core[0]^=1;
        afc_key(field,&opts,core,reference,4,changed);failures+=strcmp(key,changed)!=0;
    }
    failures+=afc_reason("new_material_contact_failed",0,0) || afc_reason("new_material_contact_failed",1,1) ||
        afc_reason("contact_walk_over_budget",1,0) || afc_reason("external_probe_failed",1,0);
    ApaFailureCache small;afc_init(&small,arena,1);
    afc_insert(&small,key,"incident_source_seam_failed",.25);afc_insert(&small,key,"incident_source_seam_failed",1);
    changed[0]=key[0]=='a'?'b':'a';afc_insert(&small,changed,"incident_source_seam_failed",1);
    if(small.count!=1 || small.full!=1 || !afc_lookup(&small,key) || afc_lookup(&small,changed) || small.previous_work_seconds!=.25)failures++;
done:
    Arena_dispose(&arena);fprintf(stderr,"  native admission failure cache: %s (%d failures)\n",failures?"FAIL":"ok",failures);return failures;
}
