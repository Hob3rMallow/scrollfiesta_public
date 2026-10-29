/* Equal material at two tessellations: translating an aligned overlap by a
 * fractional cell creates more triangle intersections while reducing its
 * physical area. Both the certificate and solver merit must accept this. */
static int apa_contact_density_control(void)
{
    int failures=0;
    for(size_t n=5;n<=9;n+=4){
        Arena_T arena=Arena_new();AsmRun run={0};run.arena=arena;run.n_charts=3;
        run.charts=ARENA_CALLOC(arena,3,sizeof *run.charts);double step=128.0/(double)(n-1);
        for(size_t c=0;c<3;c++){
            AsmChart *ch=run.charts+c;ch->id=(int32_t)c;ch->parent=-1;
            ch->nv=n*n;ch->nf=2*(n-1)*(n-1);ch->area3d=16384;
            ch->placed=1;ch->placement_state=ASM_PLACE_ROOT;ch->pose_x=c==1?32:c==2?1024:0;
            ch->xyz=ARENA_ALLOC(arena,3*ch->nv*sizeof(float));ch->uv=ARENA_ALLOC(arena,2*ch->nv*sizeof(float));
            ch->faces=ARENA_ALLOC(arena,3*ch->nf*sizeof(int32_t));ch->boundary=ARENA_ALLOC(arena,ch->nv);
            for(size_t y=0;y<n;y++)for(size_t x=0;x<n;x++){
                size_t v=y*n+x;
                ch->xyz[3*v]=ch->uv[2*v]=(float)(step*(double)x);
                ch->xyz[3*v+1]=ch->uv[2*v+1]=(float)(step*(double)y);ch->xyz[3*v+2]=(float)(512*c);
                ch->boundary[v]=(uint8_t)(!x || !y || x+1==n || y+1==n);
                if(x+1<n && y+1<n){
                    int32_t a=(int32_t)v,w=(int32_t)n,tri[6]={a,a+1,a+w+1,a,a+w+1,a+w};
                    memcpy(ch->faces+6*(y*(n-1)+x),tri,sizeof tri);
                }
            }
        }
        AsmField *f=AsmField_new(arena,&run);if(!f){failures++;Arena_dispose(&arena);continue;}
        size_t bytes=2*f->nv*sizeof(double);double *before=ARENA_ALLOC(arena,bytes);memcpy(before,f->uv,bytes);
        uint8_t moving[3]={0,1,1},core[3]={0,0,1};uint8_t *active=ARENA_CALLOC(arena,f->nf,1);
        for(size_t c=1;c<3;c++)memset(active+f->face_offset[c],1,run.charts[c].nf);
        AsmContacts_active_faces(f->contacts,active);AsmContactStats old={0},now={0};
        if(AsmContacts_measure(f->contacts,f->uv,NULL,NULL,NULL,&old))failures++;
        ApaSolverContacts state={0};state.certificate.arena=arena;state.certificate.field=f;
        state.certificate.core=core;state.certificate.reference=before;state.certificate.okay=1;
        double merit=-1;if(apa_solver_contact_merit(f,f->uv,&merit,&state) || merit!=0)failures++;
        for(size_t v=0;v<run.charts[1].nv;v++){
            size_t k=2*(f->vertex_offset[1]+v);f->uv[k]+=16.25;f->uv[k+1]+=.25;
        }
        int complete=0;const char *reason=apa_qualification(f,moving,core,before,active,NULL,4,&complete);
        int okay=!reason && complete && !AsmContacts_measure(f->contacts,f->uv,NULL,NULL,NULL,&now) &&
            now.pairs>old.pairs && now.area<old.area && !apa_solver_contact_merit(f,f->uv,&merit,&state) && merit==0;
        if(!okay)fprintf(stderr,"  contact density %zu: pairs %zu -> %zu, area %.9g -> %.9g, %s\n",
            n,old.pairs,now.pairs,old.area,now.area,reason?reason:"accepted");
        failures+=!okay;
        memcpy(f->uv,before,bytes);
        for(size_t v=0;v<run.charts[1].nv;v++)f->uv[2*(f->vertex_offset[1]+v)]-=16;
        reason=apa_qualification(f,moving,core,before,active,NULL,4,&complete);
        if(!reason || strcmp(reason,"incumbent_contact_regression") || !complete ||
           apa_solver_contact_merit(f,f->uv,&merit,&state) || !(merit>0))failures++;
        memcpy(f->uv,before,bytes);
        for(size_t v=0;v<run.charts[2].nv;v++)f->uv[2*(f->vertex_offset[2]+v)]-=960;
        reason=apa_qualification(f,moving,core,before,active,NULL,4,&complete);
        if(!reason || strcmp(reason,"new_material_contact_failed") || !complete ||
           apa_solver_contact_merit(f,f->uv,&merit,&state) || !(merit>0))failures++;
        AsmContacts_leaf_pair_budget(f->contacts,1);
        if(!apa_solver_contact_merit(f,f->uv,&merit,&state))failures++;
        Arena_dispose(&arena);
    }
    AprRank a={0},b={0};a.pairs=1000;a.area=10000;a.blocking_area=1;
    b.pairs=1;b.area=2;b.blocking_area=2;
    if(!apr_better(a,b) || apr_better(b,a))failures++;
    b=a;b.pairs=1;if(apr_better(a,b) || apr_better(b,a))failures++;
    b=a;b.bad_seams=1;if(!apr_better(a,b))failures++;
    fprintf(stderr,"  physical contact merit across tessellations, regressions and failed walks: %s\n",failures?"FAIL":"ok");
    return failures;
}

static int apa_contact_certificate_control(void)
{
    Arena_T arena=Arena_new();AsmRun run={0};int failures=0;
    if(ags_fixture(arena,&run)){Arena_dispose(&arena);return 1;}
    apr_test_scale(&run,16);run.n_rels=run.n_corr=0;
    for(size_t c=0;c<run.n_charts;c++){
        run.charts[c].placed=1;run.charts[c].placement_state=ASM_PLACE_ROOT;
        run.charts[c].pose_x=c==2?32:512*(double)c;run.charts[c].pose_y=0;
    }
    AsmField *f=AsmField_new(arena,&run);if(!f){Arena_dispose(&arena);return 1;}
    size_t bytes=2*f->nv*sizeof(double);double *reference=ARENA_ALLOC(arena,bytes);
    memcpy(reference,f->uv,bytes);
    uint8_t moving[5]={0,1,1,0,0},core[5]={0,1,0,0,0};
    uint8_t *active=ARENA_CALLOC(arena,f->nf,1);
    for(size_t c=0;c<5;c++)if(moving[c])memset(active+f->face_offset[c],1,run.charts[c].nf);
    AsmContactStats original={0};AsmContacts_active_faces(f->contacts,active);
    if(AsmContacts_measure(f->contacts,f->uv,NULL,NULL,NULL,&original) || !original.pairs)failures++;
    for(int test=0;test<7;test++){
        memcpy(f->uv,reference,bytes);size_t chart=test==4?1:test==5?0:2;
        double shift=test==1?16:test==2?-16:test==3?1536:test==4?-396:test==5?1:0;
        for(size_t v=0;v<run.charts[chart].nv;v++)f->uv[2*(f->vertex_offset[chart]+v)]+=shift;
        if(test==4)for(size_t v=0;v<run.charts[2].nv;v++)f->uv[2*(f->vertex_offset[2]+v)]+=256;
        if(test==6)AsmContacts_leaf_pair_budget(f->contacts,1);
        int complete=0;const char *reason=apa_qualification(f,moving,core,reference,active,NULL,4,&complete);
        const char *expected=test<2?NULL:test<4?"incumbent_contact_regression":test==4?"new_material_contact_failed":
            test==5?"fixed_context_changed":"contact_walk_over_budget";
        int okay=expected?reason && !strcmp(reason,expected):!reason && complete;
        if(test==6 && complete)okay=0;
        if(test<5){
            uint8_t context[5]={0};double contact[5]={0},seam[5]={0};
            ApaFeedback feedback={5,contact,seam,context,ves_clock_sec()+5,0};
            apa_feedback(f,core,moving,reference,f,NULL,arena,&feedback);
            okay=okay && feedback.complete && isfinite(feedback.violation) &&
                (test<2?feedback.blocking_material_area==0 && feedback.violation==0:
                        feedback.blocking_material_area>0 && feedback.violation>0);
            if(test==4)okay=okay && feedback.material_area<original.area;
        }
        if(!okay)fprintf(stderr,"  inherited contact control %d: %s\n",test,reason?reason:"accepted");
        failures+=!okay;
    }
    /* Growing the sparse chart-pair ledger cannot alias another pair or
     * discard its original allowance when the table is rehashed. */
    ApaContactCertificate book={0};book.arena=arena;book.okay=1;
    for(size_t k=0;k<96;k++){
        uint64_t key=((uint64_t)(INT32_MAX-k)<<32)|(uint64_t)(k+1);
        ApaContactGroup *g=apa_contact_group(&book,key);if(!g){failures++;break;}g->before_pairs=k+1;
    }
    for(size_t k=0;k<96;k++){
        uint64_t key=((uint64_t)(INT32_MAX-k)<<32)|(uint64_t)(k+1);
        ApaContactGroup *g=apa_contact_group(&book,key);if(!g || g->before_pairs!=k+1)failures++;
    }
    if(book.count!=96)failures++;
    fprintf(stderr,"  inherited contact preservation, reduction, new-pair/core refusal, fixed context and walk budget: %s\n",failures?"FAIL":"ok");
    Arena_dispose(&arena);return failures;
}

/* A refusal by a few new-material pairs writes one residual_contact record
 * per pair (charts, faces, area, 3-D gap, containment, gradients) with or
 * without moving neighbors; the trace never changes the verdict, and a
 * refusal by more pairs than APA_RESIDUAL_PAIRS writes none. */
static int apa_residual_trace_control(void)
{
    int failures=0;
    for(int test=0;test<3;test++){
        enum {N=5};Arena_T arena=Arena_new();AsmRun run={0};run.arena=arena;run.n_charts=2;
        run.charts=ARENA_CALLOC(arena,2,sizeof *run.charts);double step=128.0/(N-1);
        for(size_t c=0;c<2;c++){
            AsmChart *ch=run.charts+c;ch->id=(int32_t)c;ch->parent=-1;ch->nv=N*N;ch->nf=2*(N-1)*(N-1);ch->area3d=16384;
            ch->placed=1;ch->placement_state=ASM_PLACE_ROOT;
            /* a 12 x 12 corner overlap, deeper than the band: a few pairs; test 2 a 40-vox strip: many */
            ch->pose_x=c?(test==2?88:116):0;ch->pose_y=c?(test==2?0:116):0;
            ch->xyz=ARENA_ALLOC(arena,3*ch->nv*sizeof(float));ch->uv=ARENA_ALLOC(arena,2*ch->nv*sizeof(float));
            ch->faces=ARENA_ALLOC(arena,3*ch->nf*sizeof(int32_t));ch->boundary=ARENA_ALLOC(arena,ch->nv);
            for(size_t y=0;y<N;y++)for(size_t x=0;x<N;x++){
                size_t v=y*N+x;
                ch->xyz[3*v]=ch->uv[2*v]=(float)(step*(double)x);
                ch->xyz[3*v+1]=ch->uv[2*v+1]=(float)(step*(double)y);ch->xyz[3*v+2]=(float)(40*c);
                ch->boundary[v]=(uint8_t)(!x || !y || x+1==N || y+1==N);
                if(x+1<N && y+1<N){
                    int32_t a=(int32_t)v,w=(int32_t)N,tri[6]={a,a+1,a+w+1,a,a+w+1,a+w};
                    memcpy(ch->faces+6*(y*(N-1)+x),tri,sizeof tri);
                }
            }
        }
        AsmField *f=AsmField_new(arena,&run);if(!f){failures++;Arena_dispose(&arena);continue;}
        size_t bytes=2*f->nv*sizeof(double);double *reference=ARENA_ALLOC(arena,bytes);memcpy(reference,f->uv,bytes);
        /* test 0: only the new chart moves; tests 1-2: chart 0 is a moving neighbor */
        uint8_t core[2]={0,1},moving[2]={(uint8_t)(test>0),1};uint8_t *active=ARENA_CALLOC(arena,f->nf,1);
        for(size_t c=0;c<2;c++)if(moving[c])memset(active+f->face_offset[c],1,run.charts[c].nf);
        int complete[2]={0,0};const char *reason[2]={NULL,NULL};
        reason[0]=apa_contact_certificate(f,moving,core,reference,active,NULL,complete);
        FILE *trace=tmpfile();if(!trace){failures++;Arena_dispose(&arena);continue;}
        reason[1]=apa_contact_certificate(f,moving,core,reference,active,trace,complete+1);
        AsmContactStats stats={0};AsmContacts_active_faces(f->contacts,active);
        int okay=!AsmContacts_measure(f->contacts,f->uv,NULL,NULL,NULL,&stats) && stats.pairs>0 &&
            reason[0] && reason[1] && !strcmp(reason[0],"new_material_contact_failed") && !strcmp(reason[1],reason[0]) &&
            complete[0] && complete[1];
        size_t records=0,charts_ok=1;char line[4096];rewind(trace);
        while(fgets(line,sizeof line,trace))if(strstr(line,"\"phase\":\"residual_contact\"")){
            records++;
            charts_ok=charts_ok && strstr(line,"\"chart\":0,") && strstr(line,"\"chart\":1,") && strstr(line,"\"contained\":");
        }
        fclose(trace);
        size_t want=test==2?0:stats.pairs;
        okay=okay && (test==2?stats.pairs>APA_RESIDUAL_PAIRS:stats.pairs<=APA_RESIDUAL_PAIRS) && records==want && charts_ok;
        if(!okay)fprintf(stderr,"  residual trace %d: %s / %s, %zu pairs, %zu records\n",test,
                         reason[0]?reason[0]:"accepted",reason[1]?reason[1]:"accepted",stats.pairs,records);
        failures+=!okay;Arena_dispose(&arena);
    }
    fprintf(stderr,"  residual new-material pairs traced without changing the verdict: %s\n",failures?"FAIL":"ok");
    return failures;
}

/* The rigid residual nudge (ASM_REPAIR_RESIDUAL_NUDGE) clears a small corner
 * overlap by translating only the new chart away from its collider, leaving
 * its triangles congruent; a deep overlap it cannot clear within 2 vox, a
 * refusal that is not a new-material contact, or an expired deadline leaves
 * every UV bit-identical. */
static const char *apa_residual_nudge(AsmField *f,const uint8_t *moving,const uint8_t *core,const double *reference,
                                      const uint8_t *active,const AsmRepairOpts *opts,FILE *trace,int frontier,
                                      int *complete,const char *refused);
static int apa_residual_nudge_control(void)
{
    int failures=0;
    for(int test=0;test<4;test++){
        enum {N=9};Arena_T arena=Arena_new();AsmRun run={0};run.arena=arena;run.n_charts=2;
        run.charts=ARENA_CALLOC(arena,2,sizeof *run.charts);double step=64.0/(N-1);
        /* test 0: a 5 x 5 corner overlap (centre 2.5 vox inside both boundaries);
         * test 1: 7 x 7, a few pairs deeper than any 2-vox translation clears; test 2: no overlap;
         * test 3: test 0 past its deadline */
        double overlap=test==0 || test==3?5:test==1?7:-8;
        for(size_t c=0;c<2;c++){
            AsmChart *ch=run.charts+c;ch->id=(int32_t)c;ch->parent=-1;ch->nv=N*N;ch->nf=2*(N-1)*(N-1);ch->area3d=4096;
            ch->placed=1;ch->placement_state=ASM_PLACE_ROOT;ch->pose_x=ch->pose_y=c?64-overlap:0;
            ch->xyz=ARENA_ALLOC(arena,3*ch->nv*sizeof(float));ch->uv=ARENA_ALLOC(arena,2*ch->nv*sizeof(float));
            ch->faces=ARENA_ALLOC(arena,3*ch->nf*sizeof(int32_t));ch->boundary=ARENA_ALLOC(arena,ch->nv);
            for(size_t y=0;y<N;y++)for(size_t x=0;x<N;x++){
                size_t v=y*N+x;
                ch->xyz[3*v]=ch->uv[2*v]=(float)(step*(double)x);
                ch->xyz[3*v+1]=ch->uv[2*v+1]=(float)(step*(double)y);ch->xyz[3*v+2]=(float)(40*c);
                ch->boundary[v]=(uint8_t)(!x || !y || x+1==N || y+1==N);
                if(x+1<N && y+1<N){
                    int32_t a=(int32_t)v,w=(int32_t)N,tri[6]={a,a+1,a+w+1,a,a+w+1,a+w};
                    memcpy(ch->faces+6*(y*(N-1)+x),tri,sizeof tri);
                }
            }
        }
        AsmField *f=AsmField_new(arena,&run);if(!f){failures++;Arena_dispose(&arena);continue;}
        size_t bytes=2*f->nv*sizeof(double);double *reference=ARENA_ALLOC(arena,bytes),*before=ARENA_ALLOC(arena,bytes);
        memcpy(reference,f->uv,bytes);memcpy(before,f->uv,bytes);
        uint8_t core[2]={0,1},moving[2]={0,1};uint8_t *active=ARENA_CALLOC(arena,f->nf,1);
        memset(active+f->face_offset[1],1,run.charts[1].nf);
        AsmRepairOpts opts={0};opts.deadline=ves_clock_sec()+(test==3?-1:30);int complete=0;
        const char *initial=apa_qualification(f,moving,core,reference,active,NULL,4,&complete);
        const char *reason=initial?apa_residual_nudge(f,moving,core,reference,active,&opts,NULL,4,&complete,initial):NULL;
        size_t first=2*f->vertex_offset[1];double dx=f->uv[first]-before[first],dy=f->uv[first+1]-before[first+1];
        int rigid=1;
        for(size_t v=0;v<run.charts[1].nv;v++)
            rigid=rigid && fabs(f->uv[first+2*v]-before[first+2*v]-dx)<1e-12 && fabs(f->uv[first+2*v+1]-before[first+2*v+1]-dy)<1e-12;
        int fixed_same=!memcmp(f->uv,before,2*run.charts[0].nv*sizeof(double));
        int okay=fixed_same && rigid;
        if(test==0)okay=okay && initial && !strcmp(initial,"new_material_contact_failed") && !reason && complete &&
            dx>0 && dy>0 && hypot(dx,dy)<=2+1e-12;
        if(test==1 || test==3)okay=okay && initial && !strcmp(initial,"new_material_contact_failed") && reason &&
            !strcmp(reason,initial) && !memcmp(f->uv,before,bytes);
        if(test==2)okay=okay && !initial && !memcmp(f->uv,before,bytes);
        if(!okay)fprintf(stderr,"  residual nudge %d: %s -> %s, shift (%.4g, %.4g)\n",test,initial?initial:"accepted",
                         reason?reason:"accepted",dx,dy);
        failures+=!okay;Arena_dispose(&arena);
    }
    fprintf(stderr,"  rigid residual nudge clears a corner sliver and restores an uncleared search: %s\n",failures?"FAIL":"ok");
    return failures;
}

static int apa_contact_objective_control(AsmField *f,const uint8_t *core,const double *reference)
{
    Arena_T arena=Arena_new();int failures=0;size_t bytes=2*f->nv*sizeof(double);
    double *gradient=ARENA_ALLOC(arena,bytes);ApaSolverContacts state={0};
    state.certificate.arena=arena;state.certificate.field=f;state.certificate.core=core;
    state.certificate.reference=reference;state.certificate.okay=1;
    size_t first=f->vertex_offset[3],count=f->run->charts[3].nv;
    for(int direction=0;direction<3;direction++){
        double shift=direction==0?0:direction==1?1:-1,energy=0,force=0,new_force=0;
        memcpy(f->uv,reference,bytes);
        for(size_t v=first;v<first+count;v++)f->uv[2*v]+=shift;
        AsmContactStats all={0},focused={0};
        int okay=!AsmContacts_measure(f->contacts,f->uv,NULL,NULL,NULL,&all) && all.complete &&
            !apa_solver_contact_objective(f,f->uv,gradient,NULL,NULL,&focused,&energy,&state) && focused.complete &&
            focused.pairs==all.pairs && focused.area==all.area;
        for(size_t v=first;v<first+count;v++)force+=gradient[2*v];
        for(size_t v=f->vertex_offset[2];v<f->vertex_offset[2]+f->run->charts[2].nv;v++)new_force+=gradient[2*v];
        double expected=4.01*4.01+fmax(0,-128*shift),expected_force=shift<0?-128:0;
        okay=okay && fabs(energy-expected)<1e-7 && fabs(force-expected_force)<1e-7 && fabs(new_force+4.01)<1e-7;
        /* Central differences avoid the inherited allowance's actual kink. */
        if(direction){
            double values[2]={0},h=1e-5;
            for(int side=0;side<2;side++){
                for(size_t v=first;v<first+count;v++)f->uv[2*v]=reference[2*v]+shift+(side?1:-1)*h;
                okay=okay && !apa_solver_contact_objective(f,f->uv,NULL,NULL,NULL,&focused,values+side,&state);
            }
            okay=okay && fabs((values[1]-values[0])/(2*h)-force)<1e-4;
        }
        if(!okay)fprintf(stderr,"  contact objective variant %d: energy %.17g expected %.17g, inherited force %.17g expected %.17g, new force %.17g\n",
            direction,energy,expected,force,expected_force,new_force);
        failures+=!okay;
    }
    memcpy(f->uv,reference,bytes);AsmContacts_leaf_pair_budget(f->contacts,1);
    AsmContactStats incomplete={0};double energy=0;
    if(!apa_solver_contact_objective(f,f->uv,gradient,NULL,NULL,&incomplete,&energy,&state) ||
       memcmp(f->uv,reference,bytes))failures++;
    AsmContacts_leaf_pair_budget(f->contacts,0);
    fprintf(stderr,"  admission raw-excess objective, fixed allowances, exact forces and complete-query requirement: %s\n",failures?"FAIL":"ok");
    Arena_dispose(&arena);return failures;
}

/* A clearance step is ranked by unfinished admission work. Its total old
 * overlap may move in the opposite direction without changing that verdict. */
static int apa_clearance_merit_control(AsmField *f,const uint8_t *core,const double *reference)
{
    Arena_T arena=Arena_new();size_t bytes=2*f->nv*sizeof(double);int failures=0;
    double *before=ARENA_ALLOC(arena,bytes);
    for(int scenario=0;scenario<3;scenario++){
        memcpy(f->uv,reference,bytes);
        if(scenario==1)for(size_t v=0;v<f->run->charts[1].nv;v++)f->uv[2*(f->vertex_offset[1]+v)]-=16;
        memcpy(before,f->uv,bytes);
        AprSolver s={0};AsmRepairStats stats={0};s.field=f;s.stats=&stats;AsmRepair_defaults(&s.opts);
        ApaSolverContacts state={0};state.certificate.arena=arena;state.certificate.field=f;
        state.certificate.core=core;state.certificate.reference=reference;state.certificate.okay=1;
        s.opts.contact_merit=apa_solver_contact_merit;s.opts.contact_merit_context=&state;s.area_clearance=1;
        s.trial=ARENA_ALLOC(arena,bytes);s.direction=ARENA_CALLOC(arena,2*f->nv,sizeof(double));s.p=ARENA_ALLOC(arena,bytes);
        s.lo=ARENA_ALLOC(arena,f->nf*sizeof(double));s.hi=ARENA_ALLOC(arena,f->nf*sizeof(double));
        size_t nc=f->run->n_charts;
        s.within10=ARENA_ALLOC(arena,nc*sizeof(double));s.within25=ARENA_ALLOC(arena,nc*sizeof(double));
        s.area=ARENA_ALLOC(arena,nc*sizeof(double));s.trial10=ARENA_ALLOC(arena,nc*sizeof(double));s.trial25=ARENA_ALLOC(arena,nc*sizeof(double));
        s.fixed=ARENA_CALLOC(arena,f->nv,1);
        for(size_t c=0;c<nc;c++)for(size_t v=0;v<f->run->charts[c].nv;v++){
            size_t at=f->vertex_offset[c]+v;
            if(c>1)s.fixed[at]=1;
            else s.direction[2*at]=c==0?(scenario==1?-.02:.05):(scenario==1?8:-1);
        }
        apr_checkpoint(&s);s.energy=apr_objective(&s,f->uv,0,&s.contacts);
        double initial_merit=0,next_merit=0,initial_total=s.contacts.area;AsmContactStats proposed={0};
        apr_trial(&s,1);
        int okay=isfinite(s.energy) && apr_guard(&s,s.trial,1) &&
            !apa_solver_contact_merit(f,f->uv,&initial_merit,&state) &&
            !apa_solver_contact_merit(f,s.trial,&next_merit,&state) &&
            !AsmContacts_measure(f->contacts,s.trial,NULL,NULL,NULL,&proposed);
        if(scenario==1)okay=okay && next_merit<initial_merit && proposed.area>s.contacts.area;
        else okay=okay && next_merit>initial_merit && proposed.area<s.contacts.area;
        if(scenario==2)AsmContacts_leaf_pair_budget(f->contacts,1);
        int accepted=apr_search(&s);
        if(scenario==2)AsmContacts_leaf_pair_budget(f->contacts,0);
        okay=okay && accepted==(scenario==1) && (accepted || !memcmp(f->uv,before,bytes));
        for(size_t c=2;c<nc;c++)okay=okay && !memcmp(f->uv+2*f->vertex_offset[c],reference+2*f->vertex_offset[c],2*f->run->charts[c].nv*sizeof(double));
        if(accepted){
            uint8_t moving[4]={1,1,0,0};int complete=0;
            const char *reason=apa_qualification(f,moving,core,reference,NULL,NULL,4,&complete);
            okay=okay && !reason && complete;
            if(reason)fprintf(stderr,"  clearance merit final qualification: %s\n",reason);
        }
        if(!okay)fprintf(stderr,"  clearance merit control %d: accepted %d, blocking %.9g -> %.9g, total %.9g -> %.9g\n",
            scenario,accepted,initial_merit,next_merit,initial_total,proposed.area);
        failures+=!okay;
    }
    memcpy(f->uv,reference,bytes);
    fprintf(stderr,"  admission clearance physical merit, inherited allowance and incomplete-query rollback: %s\n",failures?"FAIL":"ok");
    Arena_dispose(&arena);return failures;
}

static int apa_joint_core_polish_control(void)
{
    Arena_T arena=Arena_new();AsmRun run={0};run.arena=arena;run.n_charts=4;int failures=0;
    run.charts=ARENA_CALLOC(arena,4,sizeof *run.charts);
    for(size_t c=0;c<4;c++){
        AsmChart *ch=run.charts+c;size_t n=c%2?9:2;double step=c%2?16:64;
        ch->id=(int32_t)c;ch->parent=-1;ch->nv=n*n;ch->nf=2*(n-1)*(n-1);
        ch->area3d=step*step*(double)((n-1)*(n-1));ch->placed=1;ch->placement_state=ASM_PLACE_ROOT;
        ch->pose_x=c==1?1024:c==2?59.99:c==3?1056:0;ch->pose_y=c==2?59.99:0;
        ch->xyz=ARENA_ALLOC(arena,3*ch->nv*sizeof(float));ch->uv=ARENA_ALLOC(arena,2*ch->nv*sizeof(float));
        ch->faces=ARENA_ALLOC(arena,3*ch->nf*sizeof(int32_t));ch->boundary=ARENA_ALLOC(arena,ch->nv);
        for(size_t y=0;y<n;y++)for(size_t x=0;x<n;x++){
            size_t v=y*n+x;ch->xyz[3*v]=(float)(step*(double)x);ch->xyz[3*v+1]=(float)(step*(double)y);ch->xyz[3*v+2]=(float)(256*c);
            ch->uv[2*v]=(float)(step*(double)x);ch->uv[2*v+1]=(float)(step*(double)y);
            ch->boundary[v]=(uint8_t)(!x || !y || x==n-1 || y==n-1);
            if(x+1<n && y+1<n){
                int32_t at=(int32_t)v,stride=(int32_t)n,tri[6]={at,at+1,at+stride+1,at,at+stride+1,at+stride};
                memcpy(ch->faces+6*(y*(n-1)+x),tri,sizeof tri);
            }
        }
    }
    AsmField *f=AsmField_new(arena,&run);if(!f){Arena_dispose(&arena);return 1;}
    size_t bytes=2*f->nv*sizeof(double);double *reference=ARENA_ALLOC(arena,bytes);memcpy(reference,f->uv,bytes);
    uint8_t moving[4]={1,1,0,0},core[4]={1,0,0,0};uint8_t *active=ARENA_CALLOC(arena,f->nf,1);
    failures+=apa_contact_objective_control(f,core,reference);
    failures+=apa_clearance_merit_control(f,core,reference);
    for(size_t c=0;c<4;c++)if(moving[c])memset(active+f->face_offset[c],1,run.charts[c].nf);
    AsmContactStats initial={0};
    if(AsmContacts_measure(f->contacts,f->uv,NULL,NULL,NULL,&initial) || initial.pairs<=16)failures++;
    AsmRepairOpts opts;AsmRepair_defaults(&opts);int complete=0,polish_complete=0;
    opts.deadline=ves_clock_sec()-1;
    const char *reason=apa_qualified_endpoint(f,moving,core,reference,active,&opts,NULL,4,1,&complete,&polish_complete);
    if(!reason || strcmp(reason,"new_material_contact_failed") || memcmp(reference,f->uv,bytes))failures++;
    opts.deadline=ves_clock_sec()+5;
    reason=apa_qualified_endpoint(f,moving,core,reference,active,&opts,NULL,4,1,&complete,&polish_complete);
    if(reason || !complete || !memcmp(reference,f->uv,2*run.charts[0].nv*sizeof(double)))failures++;
    for(size_t c=1;c<4;c++)if(memcmp(reference+2*f->vertex_offset[c],f->uv+2*f->vertex_offset[c],2*run.charts[c].nv*sizeof(double)))failures++;
    /* Rejection evidence distinguishes a bounded refusal from a completed
     * search and reports the actual retained endpoint, without changing it. */
    memcpy(f->uv,reference,bytes);AsmContacts_active_faces(f->contacts,NULL);
    AmrPolishStats observed={0};size_t steps=amr_vertex_polish_report(f,0,ves_clock_sec()-1,&observed);
    if(steps || observed.rounds || observed.trials || strcmp(observed.stop,"deadline") || memcmp(f->uv,reference,bytes))failures++;
    steps=amr_vertex_polish_report(f,0,ves_clock_sec()+.25,&observed);
    if(!steps || observed.pairs_before<=16 || !observed.trials || !(observed.area_after<observed.area_before))failures++;
    for(size_t c=1;c<4;c++)if(memcmp(reference+2*f->vertex_offset[c],f->uv+2*f->vertex_offset[c],2*run.charts[c].nv*sizeof(double)))failures++;
    memcpy(f->uv,reference,bytes);
    memset(active,0,f->nf);memset(active+f->face_offset[0],1,run.charts[0].nf);AsmContacts_active_faces(f->contacts,active);
    steps=amr_vertex_polish_report(f,0,ves_clock_sec()+5,&observed);
    AsmContactStats polished={0};
    if(!steps || !observed.trials || !observed.pairs_before || observed.pairs_after || observed.area_after ||
       observed.trials!=observed.metric_refused+observed.seam_refused+observed.contact_queries || observed.query_failed ||
       AsmContacts_measure(f->contacts,f->uv,NULL,NULL,NULL,&polished) || polished.pairs || !apa_chart_metric(f,0))failures++;
    for(size_t c=1;c<4;c++)if(memcmp(reference+2*f->vertex_offset[c],f->uv+2*f->vertex_offset[c],2*run.charts[c].nv*sizeof(double)))failures++;
    fprintf(stderr,"  joint core polish beside large inherited contact, fixed neighbors and deadline: %s (%s)\n",
        failures?"FAIL":"ok",reason?reason:"accepted");
    Arena_dispose(&arena);return failures;
}

/* The same material overlap remains polishable after subdivision, even
 * when more intersecting triangle pairs accompany a smaller physical area. */
/* On fine triangles the certified interior has no area gradient: the force
 * belongs to original boundary vertices outside those retained triangles.
 * A third chart has tolerated overlap only and must not bias the fallback
 * objective, although it remains in every physical contact measurement. */
static int apa_polish_interface_control(void)
{
    Arena_T arena=Arena_new();AsmRun run={0};run.arena=arena;run.n_charts=3;
    run.charts=ARENA_CALLOC(arena,3,sizeof *run.charts);size_t n=33;
    for(size_t c=0;c<3;c++){
        AsmChart *ch=run.charts+c;ch->id=(int32_t)c;ch->parent=-1;ch->nv=n*n;ch->nf=2*(n-1)*(n-1);
        ch->area3d=1024;ch->placed=1;ch->placement_state=ASM_PLACE_ROOT;ch->pose_x=c==1?16:c==2?31:0;
        ch->xyz=ARENA_ALLOC(arena,3*ch->nv*sizeof(float));ch->uv=ARENA_ALLOC(arena,2*ch->nv*sizeof(float));
        ch->faces=ARENA_ALLOC(arena,3*ch->nf*sizeof(int32_t));ch->boundary=ARENA_ALLOC(arena,ch->nv);
        for(size_t y=0;y<n;y++){
            for(size_t x=0;x<n;x++){
                size_t v=y*n+x;ch->xyz[3*v]=ch->uv[2*v]=(float)x;
                ch->xyz[3*v+1]=ch->uv[2*v+1]=(float)y;ch->xyz[3*v+2]=(float)(512*c);
                ch->boundary[v]=(uint8_t)(!x || !y || x+1==n || y+1==n);
                if(x+1<n && y+1<n){
                    int32_t a=(int32_t)v,w=(int32_t)n,tri[]={a,a+1,a+w+1,a,a+w+1,a+w};
                    memcpy(ch->faces+6*(y*(n-1)+x),tri,sizeof tri);
                }
            }
        }
    }
    AsmField *f=AsmField_new(arena,&run);
    if(!f){Arena_dispose(&arena);return 1;}
    size_t bytes=2*f->nv*sizeof(double);double *before=ARENA_ALLOC(arena,bytes),*gradient=ARENA_ALLOC(arena,bytes);
    memcpy(before,f->uv,bytes);
    uint8_t *active=ARENA_CALLOC(arena,f->nf,1),*vertices=ARENA_CALLOC(arena,f->nv,1),opponents[3]={0};
    memset(active+f->face_offset[0],1,run.charts[0].nf);AsmContacts_active_faces(f->contacts,active);
    AmrPolishInterface interface={{f,0,vertices},opponents};AsmContactStats initial={0},weighted={0},after={0};
    double raw=0,after_raw=0,translation=0;
    int okay=!AsmContacts_measure(f->contacts,f->uv,NULL,amr_polish_interface,&interface,&initial) &&
        initial.complete && initial.pairs && opponents[1] && !opponents[2];
    okay=okay && !AsmContacts_weighted_objective(f->contacts,f->uv,gradient,NULL,NULL,
        amr_polish_weight,&interface,&weighted,&raw) && weighted.complete && weighted.pairs==initial.pairs &&
        fabs(weighted.area-initial.area)<1e-9 && fabs(raw-512)<1e-8;
    for(size_t v=0;v<run.charts[0].nv;v++){translation+=gradient[2*v];}
    okay=okay && fabs(translation-32)<1e-8;
    AmrPolishStats stats={0};size_t steps=amr_vertex_polish_report(f,0,ves_clock_sec()+1,&stats);
    okay=okay && steps && stats.focused_rounds && stats.boundary_candidates && stats.plateau_steps &&
        !AsmContacts_weighted_objective(f->contacts,f->uv,NULL,NULL,NULL,
            amr_polish_weight,&interface,&after,&after_raw) && after.complete && after.area<=initial.area &&
        after_raw<raw && apa_chart_metric(f,0);
    for(size_t c=1;c<3;c++){
        okay=okay && !memcmp(before+2*f->vertex_offset[c],f->uv+2*f->vertex_offset[c],2*run.charts[c].nv*sizeof(double));
    }
    fprintf(stderr,"  whole blocking-interface polish: %s (%zu steps, %zu focused rounds, %zu boundary candidates, %zu plateau steps)\n",
        okay?"ok":"FAIL",steps,stats.focused_rounds,stats.boundary_candidates,stats.plateau_steps);
    Arena_dispose(&arena);return !okay;
}

static int apa_polish_density_control(void)
{
    int failures=0;
    for(size_t n=5;n<=17;n+=12){
        Arena_T arena=Arena_new();AsmRun run={0};run.arena=arena;run.n_charts=2;
        run.charts=ARENA_CALLOC(arena,2,sizeof *run.charts);double step=128.0/(double)(n-1);
        for(size_t c=0;c<2;c++){
            AsmChart *ch=run.charts+c;ch->id=(int32_t)c;ch->parent=-1;ch->nv=n*n;ch->nf=2*(n-1)*(n-1);
            ch->area3d=16384;ch->placed=1;ch->placement_state=ASM_PLACE_ROOT;ch->pose_x=c?32:0;
            ch->xyz=ARENA_ALLOC(arena,3*ch->nv*sizeof(float));ch->uv=ARENA_ALLOC(arena,2*ch->nv*sizeof(float));
            ch->faces=ARENA_ALLOC(arena,3*ch->nf*sizeof(int32_t));ch->boundary=ARENA_ALLOC(arena,ch->nv);
            for(size_t y=0;y<n;y++)for(size_t x=0;x<n;x++){
                size_t v=y*n+x;ch->xyz[3*v]=ch->uv[2*v]=(float)(step*(double)x);
                ch->xyz[3*v+1]=ch->uv[2*v+1]=(float)(step*(double)y);ch->xyz[3*v+2]=(float)(512*c);
                ch->boundary[v]=(uint8_t)(!x || !y || x+1==n || y+1==n);
                if(x+1<n && y+1<n){int32_t a=(int32_t)v,w=(int32_t)n,tri[]={a,a+1,a+w+1,a,a+w+1,a+w};memcpy(ch->faces+6*(y*(n-1)+x),tri,sizeof tri);}
            }
        }
        AsmField *f=AsmField_new(arena,&run);if(!f){failures++;Arena_dispose(&arena);continue;}
        size_t bytes=2*f->nv*sizeof(double);double *before=ARENA_ALLOC(arena,bytes);memcpy(before,f->uv,bytes);
        uint8_t *vertices=ARENA_CALLOC(arena,f->nv,1);AmbContactSeeds seeds={f,0,vertices};AsmContactStats initial={0},after={0};
        int okay=!AsmContacts_measure(f->contacts,f->uv,NULL,amb_contact_seed,&seeds,&initial) && initial.complete && initial.pairs>16;
        size_t involved=0;for(size_t v=0;v<f->nv;v++)involved+=vertices[v]!=0;
        if(n==17 && involved<=96)okay=0;
        AmrPolishStats stats={0};size_t steps=amr_vertex_polish_report(f,0,ves_clock_sec()+.5,&stats);
        okay=okay && steps && stats.trials && stats.pair_increasing_steps &&
            !AsmContacts_measure(f->contacts,f->uv,NULL,NULL,NULL,&after) && after.complete && after.area<initial.area &&
            after.pairs==stats.pairs_after && fabs(after.area-stats.area_after)<1e-9 && apa_chart_metric(f,0) &&
            !memcmp(before+2*f->vertex_offset[1],f->uv+2*f->vertex_offset[1],2*run.charts[1].nv*sizeof(double));
        if(n==17){
            /* Useful contact-only work continues past the historical .25 s
             * probe, honors the transaction deadline, and remains private. */
            memcpy(f->uv,before,bytes);uint8_t core[2]={1,0};int complete=1;
            AsmRepairOpts opts;AsmRepair_defaults(&opts);opts.iterations=1;
            double started=ves_clock_sec();opts.deadline=started+.8;
            steps=apa_contact_polish(f,core,&opts,NULL,&complete);
            AsmContactStats retained={0};
            okay=okay && steps && !complete &&
                ves_clock_sec()>=started+.8*ASM_REPAIR_PRE_FIT_CLEARANCE_SHARE &&
                !AsmContacts_measure(f->contacts,f->uv,NULL,NULL,NULL,&retained) && retained.complete &&
                retained.area<initial.area && apa_chart_metric(f,0) &&
                !memcmp(before+2*f->vertex_offset[1],f->uv+2*f->vertex_offset[1],2*run.charts[1].nv*sizeof(double));
        }
        if(n==5){
            /* The collider must already belong to the joint movable set.
             * Its cleanup preserves the new chart and every caller-fixed
             * chart, while retaining the same metric/contact measurements. */
            memcpy(f->uv,before,bytes);
            uint8_t core[2]={1,0},moving[2]={1,0},fixed[2]={0,1};int complete=1;
            AsmRepairOpts opts;AsmRepair_defaults(&opts);opts.iterations=1;
            opts.admission_moving=moving;opts.deadline=ves_clock_sec()+.5;
            okay=okay && !apa_neighbor_polish(f,moving,core,&opts,NULL,&complete) && !memcmp(before,f->uv,bytes);
            moving[1]=1;opts.fixed_charts=fixed;
            okay=okay && !apa_neighbor_polish(f,moving,core,&opts,NULL,&complete) && !memcmp(before,f->uv,bytes);
            opts.fixed_charts=NULL;opts.deadline=ves_clock_sec()+.5;
            steps=apa_neighbor_polish(f,moving,core,&opts,NULL,&complete);
            AsmContactStats retained={0};
            okay=okay && steps && !AsmContacts_measure(f->contacts,f->uv,NULL,NULL,NULL,&retained) &&
                retained.complete && retained.area<initial.area && apa_chart_metric(f,1) &&
                !memcmp(before+2*f->vertex_offset[0],f->uv+2*f->vertex_offset[0],2*run.charts[0].nv*sizeof(double));
            /* Returning to the core after a neighbor moves is real bounded
             * geometry work, with the same physical endpoint certificate.
             * The short probe shares permit both directions to be revisited
             * even when neither overlapping chart clears in its first sweep. */
            memcpy(f->uv,before,bytes);opts.iterations=opts.clearance_iterations=0;
            opts.deadline=ves_clock_sec()+1.5;
            uint8_t *joint_active=ARENA_ALLOC(arena,f->nf);memset(joint_active,1,f->nf);
            char path[256];snprintf(path,sizeof path,"output/asm_joint_polish_%d.jsonl",ves_getpid());
            FILE *trace=fopen(path,"wb+");int certified=0,polished=0;size_t passes=0;
            const char *reason=apa_qualified_endpoint(f,moving,core,before,joint_active,&opts,trace,4,1,&certified,&polished);
            if(trace){
                char line[4096];rewind(trace);
                while(fgets(line,sizeof line,trace))passes+=strstr(line,"\"phase\":\"joint_contact_pass\"")!=NULL;
                fclose(trace);remove(path);
            }
            AsmContacts_active_faces(f->contacts,NULL);
            okay=okay && trace && (passes>=2 || !reason) && certified &&
                (!reason || !strcmp(reason,"new_material_contact_failed")) &&
                !AsmContacts_measure(f->contacts,f->uv,NULL,NULL,NULL,&retained) && retained.complete &&
                retained.area<initial.area && apa_chart_metric(f,0) && apa_chart_metric(f,1);
            fprintf(stderr,"  alternating core/collider polish: %zu passes, residual %.9g, %s\n",
                    passes,retained.area,reason?reason:"qualified");
        }
        memcpy(f->uv,before,bytes);AsmContacts_leaf_pair_budget(f->contacts,1);
        steps=amr_vertex_polish_report(f,0,ves_clock_sec()+.5,&stats);
        okay=okay && !steps && stats.query_failed && !strcmp(stats.stop,"contact_query_failed") && !memcmp(before,f->uv,bytes);
        fprintf(stderr,"  physical polish density %zu: %zu involved vertices, pairs %zu -> %zu, area %.9g -> %.9g: %s\n",
            n,involved,initial.pairs,after.pairs,initial.area,after.area,okay?"ok":"FAIL");
        failures+=!okay;Arena_dispose(&arena);
    }
    return failures;
}

static int apa_inherited_self_contact_control(void)
{
    Arena_T arena=Arena_new();AsmRun run={0};run.arena=arena;run.n_charts=2;int failures=0;
    run.charts=ARENA_CALLOC(arena,2,sizeof *run.charts);
    for(size_t c=0;c<2;c++){
        AsmChart *ch=run.charts+c;size_t pieces=c?1:3;ch->id=(int32_t)c;ch->parent=-1;
        ch->nv=4*pieces;ch->nf=2*pieces;ch->area3d=4096*(double)pieces;
        ch->placed=1;ch->placement_state=ASM_PLACE_ROOT;
        ch->xyz=ARENA_ALLOC(arena,3*ch->nv*sizeof(float));ch->uv=ARENA_ALLOC(arena,2*ch->nv*sizeof(float));
        ch->faces=ARENA_ALLOC(arena,3*ch->nf*sizeof(int32_t));ch->boundary=ARENA_ALLOC(arena,ch->nv);
        memset(ch->boundary,1,ch->nv);
        for(size_t p=0;p<pieces;p++){
            const float xy[8]={0,0,64,0,64,64,0,64};
            float position=c?1000:p==0?0:p==1?16:256;
            for(size_t v=0;v<4;v++){
                size_t k=4*p+v;ch->xyz[3*k]=(float)(256*p)+xy[2*v];ch->xyz[3*k+1]=xy[2*v+1];ch->xyz[3*k+2]=(float)(32*c);
                ch->uv[2*k]=position+xy[2*v];ch->uv[2*k+1]=xy[2*v+1];
            }
            int32_t first=(int32_t)(4*p),tri[6]={first,first+1,first+2,first,first+2,first+3};
            memcpy(ch->faces+6*p,tri,sizeof tri);
        }
    }
    AsmField *f=AsmField_new(arena,&run);if(!f){Arena_dispose(&arena);return 1;}
    double *before=ARENA_ALLOC(arena,2*f->nv*sizeof(double));memcpy(before,f->uv,2*f->nv*sizeof(double));
    uint8_t moving[2]={1,1},core[2]={0,1};uint8_t *active=ARENA_ALLOC(arena,f->nf);memset(active,1,f->nf);
    /* Replace an old A/B self-overlap with a smaller B/C self-overlap in
     * the same chart. Chart totals improve, but new self intersections are
     * still forbidden at the original-face-pair level. */
    for(size_t v=4;v<8;v++)f->uv[2*v]+=272;
    char path[256];snprintf(path,sizeof path,"output/asm_inherited_self_%d.jsonl",ves_getpid());
    FILE *trace=fopen(path,"wb+");int complete=0;
    const char *reason=apa_qualification(f,moving,core,before,active,trace,4,&complete);
    if(!reason || strcmp(reason,"incumbent_contact_regression") || !complete || !trace)failures++;
    if(trace){
        char line[2048];rewind(trace);
        if(!fgets(line,sizeof line,trace) || !strstr(line,"\"neighbor_chart_pair_regressions\":0") ||
           strstr(line,"\"self_contact_regressions\":0"))failures++;
        fclose(trace);remove(path);
    }
    fprintf(stderr,"  reduced total overlap cannot exchange old self intersections for new ones: %s\n",failures?"FAIL":"ok");
    Arena_dispose(&arena);return failures;
}
