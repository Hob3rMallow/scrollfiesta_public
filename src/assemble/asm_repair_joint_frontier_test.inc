static int ajf_incident_sources_control(void)
{
    int failures=0;
    for(int variant=0;variant<4;variant++){
        Arena_T arena=Arena_new();AsmRun run={0};
        if(ags_grid_fixture(arena,&run)){Arena_dispose(&arena);return 1;}
        for(size_t c=0;c<run.n_charts;c++){run.charts[c].pose_x=20*(double)(c%3);run.charts[c].pose_y=20*(double)(c/3);}
        size_t first=SIZE_MAX,optional=SIZE_MAX,incident=0;
        for(size_t i=0;i<run.n_rels;i++)if(run.rels[i].a==4 || run.rels[i].b==4){
            incident++;
            if(first==SIZE_MAX)first=i;
            else{run.rels[i].continuity&=~ASM_CONT_REQUIRED;if(optional==SIZE_MAX)optional=i;}
        }
        if(first==SIZE_MAX || optional==SIZE_MAX){failures++;Arena_dispose(&arena);continue;}
        if(variant)run.rels[optional].flags|=variant==1?ASM_REL_CROSSWRAP:variant==2?ASM_REL_DROPPED:ASM_REL_PARITY;
        AsmChart before[9];AsmRelation relations[12];memcpy(before,run.charts,sizeof before);memcpy(relations,run.rels,sizeof relations);
        uint8_t core[9]={0},moving[9]={0};core[4]=1;for(size_t c=1;c<9;c++)moving[c]=c!=4;
        AsmField *parent=AsmField_new(arena,&run);ApaBox boxes[9];apa_chart_boxes(&run,boxes);
        AsmRepairOpts opts;AsmRepair_defaults(&opts);opts.compact=1;opts.admission_moving=moving;opts.admission_incident_sources=1;
        opts.iterations=opts.clearance_iterations=opts.closure_iterations=0;opts.admission_rounds=1;opts.deadline=ves_clock_sec()+5;
        AsmRepairStats stats={0};ApaReport report={0};size_t anchors=0;const char *reason="";
        int accepted=parent && !AsmContacts_refit(parent->contacts,parent->uv) &&
            apa_attempt(&run,&opts,core,&anchors,&stats,&reason,NULL,boxes,parent,NULL,&report,4);
        int okay=accepted==(variant!=3);
        if(accepted){
            size_t selected=0;
            for(size_t i=0;i<run.n_rels;i++)if(run.rels[i].a==4 || run.rels[i].b==4)selected+=apa_required_source(run.rels+i)!=0;
            okay=okay && selected==incident-(variant==1 || variant==2) && (run.rels[first].continuity&ASM_CONT_REQUIRED);
            for(size_t i=0;i<run.n_rels;i++){
                AsmRelation next=run.rels[i];next.continuity=relations[i].continuity;
                okay=okay && !memcmp(&next,relations+i,sizeof next);
            }
        }else okay=okay && !strcmp(reason,"contradictory_source_parity") &&
            !memcmp(before,run.charts,sizeof before) && !memcmp(relations,run.rels,sizeof relations);
        fprintf(stderr,"  joint incident source closure %d: %s (%s)\n",variant,okay?"ok":"FAIL",reason);
        failures+=!okay;Arena_dispose(&arena);
    }
    return failures;
}

static int ajf_selftest(void)
{
    enum {N=9,NV=N*N,NF=2*(N-1)*(N-1),NC=5};int failures=ajf_incident_sources_control();
    for(int test=0;test<22;test++){
        Arena_T arena=Arena_new();AsmRun run={0};run.arena=arena;
        AsmChart charts[NC]={{0}};float xyz[NC][3*NV],uv[NC][2*NV];int32_t faces[3*NF];uint8_t boundary[NV];
        for(int y=0;y<N;y++)for(int x=0;x<N;x++){
            int v=y*N+x;boundary[v]=(uint8_t)(!x || !y || x==N-1 || y==N-1);
            for(int c=0;c<NC;c++){
                xyz[c][3*v]=(float)(10*c+x);xyz[c][3*v+1]=(float)y;xyz[c][3*v+2]=0;
                uv[c][2*v]=(float)x;uv[c][2*v+1]=(float)y;
            }
            if(x+1<N && y+1<N){int32_t tri[]={v,v+1,v+N+1,v,v+N+1,v+N};memcpy(faces+6*(y*(N-1)+x),tri,sizeof tri);}
        }
        if(test==2)for(int v=0;v<NV;v++)uv[2][2*v]*=.5f;
        for(int c=0;c<NC;c++){
            AsmChart *ch=charts+c;ch->id=c;ch->parent=-1;ch->nv=NV;ch->nf=NF;ch->area3d=64;
            ch->xyz=xyz[c];ch->uv=uv[c];ch->faces=faces;ch->boundary=boundary;
            if(c!=1 && !(test==4 && c==4)){ch->placed=1;ch->placement_state=ASM_PLACE_ROOT;}
            ch->pose_x=c==0?100:c==2?117.8:c==3?127.8:1000;ch->pose_y=200;
            /* A wide exterior overlap must be refused before admission,
             * including when the diagnostic endpoint path is unwritable. */
            if((test==10 || test==13) && c==4)ch->pose_x=110.3;
            if(test==9){ch->pose_x+=1e9;ch->pose_y+=1e9;}
        }
        if(test==18){
            charts[2].pose_x=118;charts[3].pose_x=128;
            charts[4].pose_x=120;charts[4].pose_y=202;
        }
        run.charts=charts;run.n_charts=NC;
        for(int i=0;i<2;i++){
            AsmRelation r={0};r.a=i?2:0;r.b=i?3:1;if(test==7 && i)r.flags=ASM_REL_WEAK;
            r.corr_first=(int32_t)run.n_corr;r.corr_count=N;r.n_corr=N;
            for(int y=0;y<N;y++){AsmCorr w={0};w.va=y*N+N-1;w.vb=y*N;AsmRun_push_corr(&run,&w);}
            AsmRun_push_rel(&run,&r);
        }
        if(AsmContinuity_prepare(&run)){failures++;Arena_dispose(&arena);continue;}
        apr_test_scale(&run,64);
        AsmChart before[NC];AsmRelation relation=run.rels[1];AsmCorr corr[2*N];
        float old_xyz[NC][3*NV],old_uv[NC][2*NV];int32_t old_faces[3*NF];
        memcpy(before,charts,sizeof charts);memcpy(corr,run.corr,sizeof corr);memcpy(old_xyz,xyz,sizeof xyz);memcpy(old_uv,uv,sizeof uv);memcpy(old_faces,faces,sizeof faces);
        AsmField *parent=AsmField_new(arena,&run);ApaBox boxes[NC];apa_chart_boxes(&run,boxes);
        if(!parent || AsmContacts_refit(parent->contacts,parent->uv)){failures++;Arena_dispose(&arena);continue;}
        uint8_t core[NC]={0,1,0,0,0},moving[NC]={0};moving[test==4?4:2]=1;
        AsmRepairOpts opts;AsmRepair_defaults(&opts);opts.compact=1;opts.iterations=opts.clearance_iterations=8;
        opts.closure_iterations=0;opts.clearance_method=2;opts.admission_rounds=1;opts.admission_moving=test==1?NULL:moving;
        opts.deadline=ves_clock_sec()+(test==3?-1:5);
        if(test==18)opts.iterations=opts.clearance_iterations=0;
        ApaContinuation guess={0};FILE *trace=NULL;char trace_path[512]={0};
        if(test>=20){
            /* This warm field passes metrics and seams but cannot finish
             * with the short clearance pass alone. Require the fitting
             * fallback and the full certificate, including rollback when
             * its shared deadline has already expired. */
            AsmRun trial=run;AsmChart candidate[NC];memcpy(candidate,charts,sizeof candidate);trial.charts=candidate;
            candidate[1].placed=1;candidate[1].placement_state=ASM_PLACE_SEAM;
            AsmField *warm=AsmField_new_local(arena,&trial,NULL,NULL);
            if(!warm){failures++;Arena_dispose(&arena);continue;}
            size_t bytes=2*warm->nv*sizeof(double);double *entry=ARENA_ALLOC(arena,bytes);memcpy(entry,warm->uv,bytes);
            for(size_t v=0;v<candidate[1].nv;v++)warm->uv[2*(warm->vertex_offset[1]+v)]+=(110.-1000.)*64.;
            uint8_t active[NC]={0,1,1,0,0};apa_continuation_save(&guess,warm,active,entry);
            opts.admission_continuation=&guess;opts.clearance_method=0;opts.clearance_iterations=1;
            if(test==21)opts.deadline=ves_clock_sec()-1;
            ves_mkdir("output");snprintf(trace_path,sizeof trace_path,"output/asm_joint_warm_%d_%d.jsonl",ves_getpid(),test);
            trace=fopen(trace_path,"wb+");
            if(!guess.arena || !trace){failures++;if(trace)fclose(trace);apa_continuation_clear(&guess);Arena_dispose(&arena);continue;}
        }
        char endpoint[1024];
        if((test>=11 && test<=13) || test==16){
            ves_mkdir("output");snprintf(endpoint,sizeof endpoint,"output/asm_joint_endpoint_%d_%d",ves_getpid(),test);
            if(test>=12){FILE *blocked=fopen(endpoint,"wb");if(!blocked){failures++;Arena_dispose(&arena);continue;}fputs("file blocks directory",blocked);fclose(blocked);}
            opts.admission_trace_dir=endpoint;
        }
        if(test==5)opts.patch_faces=2*NF-1;if(test==6)opts.patch_field_faces=2*NF;
        if(test==14)opts.patch_coordinates=4*NV-1;
        if(test==19)opts.patch_bytes=1;
        ApaOverlay overlay={0};overlay.arena=arena;uint8_t changed[NC]={0};
        if(test>=15){overlay.changed_charts=changed;overlay.n_charts=test==17?NC-1:NC;overlay.batched=overlay.query_scaled_budget=1;}
        AsmRepairStats stats={0};ApaReport report={0};size_t anchors=0;const char *reason="";
        int accepted=apa_attempt(&run,&opts,core,&anchors,&stats,&reason,trace,boxes,parent,test==8 || test>=15?&overlay:NULL,&report,4);
        int expected=test==0 || test==1 || test==7 || test==9 || test==11 || test==15 || test==18 || test==20;
        int okay=accepted==expected && !memcmp(old_xyz,xyz,sizeof xyz) && !memcmp(old_uv,uv,sizeof uv) && !memcmp(old_faces,faces,sizeof faces) &&
            !memcmp(corr,run.corr,sizeof corr) && !memcmp(&relation,run.rels+1,sizeof relation);
        if(trace){
            int observed=0;char line[2048];rewind(trace);
            while(fgets(line,sizeof line,trace))if(strstr(line,"\"phase\":\"phase_order\"")){
                fprintf(stderr,"  warm phase order: %s",line);
                if(strstr(line,"\"clearance_first\":true") && !strstr(line,"\"pre_fit_clearance_steps\":0") &&
                   !strstr(line,"\"fitting_iterations\":0"))observed=1;
            }
            okay=okay && (test==20?observed:!observed);fclose(trace);remove(trace_path);
        }
        for(int c=0;c<NC;c++)if(c!=1 && c!=2)okay=okay && !memcmp(before+c,charts+c,sizeof(AsmChart));
        if(!accepted)okay=okay && !memcmp(before,charts,sizeof charts) && !stats.admitted_charts;
        else{
            okay=okay && stats.admitted_charts==1 && AsmChart_registered(charts+1);
            AsmField *after=AsmField_new(arena,&run);AsmContactStats contacts={0};
            if(!after || AsmContacts_measure(after->contacts,after->uv,NULL,NULL,NULL,&contacts) ||
               (test==18?!contacts.pairs:contacts.pairs!=0)){okay=0;fprintf(stderr,"  joint frontier %d final field: %zu pairs, complete %d\n",test,contacts.pairs,contacts.complete);}
            else{
                for(size_t c=0;c<NC;c++)if(AsmChart_registered(charts+c))okay=okay && apa_chart_metric(after,c);
                for(size_t i=0;i<after->n_seams;i++){AsmSeamMeasure m;AsmField_seam(after,i,after->uv,&m);okay=okay && m.pass;if(!m.pass)fprintf(stderr,"  joint frontier %d final seam %zu rms %.17g\n",test,i,m.rms);}
            }
            if(test==1)okay=okay && !memcmp(before+2,charts+2,sizeof(AsmChart));
            else if(test==18){double old[2],next[2];AsmChart_point(before+2,0,old);AsmChart_point(charts+2,0,next);okay=okay && old[0]==next[0] && old[1]==next[1];}
            else{double old[2],next[2];AsmChart_point(before+2,0,old);AsmChart_point(charts+2,0,next);okay=okay && hypot(next[0]-old[0],next[1]-old[1])>1e-5;}
        }
        if(test==2)okay=okay && !strcmp(reason,"joint_incumbent_metric_failed");
        if(test==3 || test==21)okay=okay && !strcmp(reason,"shared_deadline");
        if(test==4)okay=okay && !strcmp(reason,"joint_chart_not_represented");
        if(test==5 || test==14)okay=okay && !strcmp(reason,"moving_resource_cap");
        if(test==6)okay=okay && !strcmp(reason,"halo_resource_cap");
        if(test==19)okay=okay && !strcmp(reason,"workspace_memory_budget");
        if(test==8)okay=okay && !strcmp(reason,"joint_context_unavailable");
        if(test==10 || test==13)okay=okay && !strcmp(reason,"candidate_overlaps_placed_material");
        if(test==11){
            char path[1200],line[2048];snprintf(path,sizeof path,"%s/chart_0/candidate_vertices.csv",endpoint);FILE *file=fopen(path,"rb");
            size_t lines=0;if(file){while(fgets(line,sizeof line,file))lines++;fclose(file);}
            okay=okay && lines==NV+1 && !report.endpoint_export_failed;
        }
        if((test>=12 && test<=13) || test==16)okay=okay && report.endpoint_export_failed;
        if(test==12 || test==16)okay=okay && !strcmp(reason,"endpoint_export_failed");
        if(test==17)okay=okay && !strcmp(reason,"overlay_failed");
        if(test>=15){
            okay=okay && overlay.nf==(accepted?2*NF:0);
            for(int c=0;c<NC;c++)okay=okay && changed[c]==(accepted && (c==1 || c==2));
        }
        fprintf(stderr,"  native joint frontier %d: accepted %d (%s), %s\n",test,accepted,reason,okay?"ok":"FAIL");
        failures+=!okay;apa_overlay_clear(&overlay);apa_continuation_clear(&guess);Arena_dispose(&arena);
    }
    return failures;
}

/* The pre-solve shortcut against fixed material decides on overlap area and
 * the 3-D truth: the candidate's triangle, placed material at the same place
 * (1 vox away) and placed material a wrap away (9.5 vox along the normal).
 * Tests 10 and 13 above keep the whole-transaction refusal (71% overlap by
 * material 30 vox away in 3-D). */
static int apa_start_area_control(void)
{
    /* Face 0 (the candidate) and face 1 (placed) are coarse triangles of one
     * sheet overlapping in uv: xyz = (u, v, 0) on both, so the materials
     * coincide under every shared uv point although their centroids lie 9.5
     * vox apart.  Face 2 is placed material a wrap away: the same uv
     * triangle as face 1 at xyz = (u, v, 9.5). */
    double uv[18]={0,0, 24,0, 0,24,  9,3, 33,3, 9,27,  9,3, 33,3, 9,27};
    float xyz[27];
    for(int v=0;v<9;v++){xyz[3*v]=(float)uv[2*v];xyz[3*v+1]=(float)uv[2*v+1];xyz[3*v+2]=v>=6?9.5f:0.0f;}
    int32_t faces[9]={0,1,2, 3,4,5, 6,7,8},face_chart[3]={0,1,2};
    AsmField f;memset(&f,0,sizeof f);f.xyz=xyz;f.uv=uv;f.faces=faces;f.face_chart=face_chart;f.nf=3;f.nv=9;
    uint8_t core[3]={1,0,0},moving[3]={1,0,0};int failures=0;
    failures+=apa_face_gap(&f,0,1)>1e-6 || fabs(apa_face_gap(&f,0,2)-9.5)>1e-6 || fabs(apa_face_gap(&f,2,0)-9.5)>1e-6;
    ApaInitialContacts s={&f,core,moving,0,0,0,0};
    apa_initial_contact(&s,0,1,2.0);   /* the same place on the sheet */
    apa_initial_contact(&s,2,0,3.0);   /* a wrap away, in either order */
    apa_initial_contact(&s,1,2,5.0);   /* two fixed charts: not the candidate's contact */
    failures+=s.blocking_pairs!=2 || fabs(s.blocking_area-5.0)>1e-12 || s.far_pairs!=1 || fabs(s.far_area-3.0)>1e-12;
    /* A 100-vox^2 candidate of 100 faces starting on 400 pairs. */
    ApaInitialContacts d={&f,core,moving,400,40.0,0,10.0};
    int near_ok=apa_start_hopeless(&d,100,100);                  /* far 10%, all 40%: solve it */
    d.far_area=16.0;int far_refused=apa_start_hopeless(&d,100,100);
    d.far_area=0.0;d.blocking_area=51.0;int covered_refused=apa_start_hopeless(&d,100,100);
    d.blocking_area=40.0;int no_area=apa_start_hopeless(&d,0,100);   /* degenerate area: the pair ratio */
    int no_faces=apa_start_hopeless(&d,100,0);
    if(apa_area_3d_test)failures+=near_ok || !far_refused || !covered_refused || !no_area || no_faces;
    else failures+=!near_ok || !far_refused || !covered_refused || !no_area || no_faces;
    d.blocking_pairs=50;d.blocking_area=90.0;d.far_area=90.0;         /* few pairs, fully on another wrap */
    int fine_refused=apa_start_hopeless(&d,100,100);
    failures+=apa_area_3d_test?!fine_refused:fine_refused;
    fprintf(stderr,"  admission start-area control (%s): %s\n",apa_area_3d_test?"area and 3-D":"pair ratio",failures?"FAIL":"ok");
    return failures;
}
