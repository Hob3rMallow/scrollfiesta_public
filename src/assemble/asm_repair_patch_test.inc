static int apr_contact_pose_control(void)
{
    /* Equal overlaps on opposite ends of a bar have zero translation force
     * but nonzero torque. Translation-only clearance must stall; a proper
     * rotation must reduce real material contact without straining a face. */
    Arena_T arena=Arena_new();AsmRun run={0};run.arena=arena;run.n_charts=3;
    AsmChart charts[3]={{0}};float xyz[3][12],uv[3][8];int32_t faces[]={0,1,2,0,2,3};
    const double box[3][4]={{-20,-5,20,5},{10,0,30,10},{-30,-10,-10,0}};
    run.charts=charts;
    for(int c=0;c<3;c++){
        AsmChart *ch=charts+c;ch->id=c;ch->parent=-1;ch->nv=4;ch->nf=2;ch->xyz=xyz[c];ch->uv=uv[c];ch->faces=faces;
        ch->placed=1;ch->placement_state=ASM_PLACE_ROOT;
        for(int k=0;k<4;k++){
            double x=box[c][k==1 || k==2?2:0],y=box[c][k>=2?3:1];
            uv[c][2*k]=(float)x;uv[c][2*k+1]=(float)y;
            xyz[c][3*k]=(float)x;xyz[c][3*k+1]=(float)y;xyz[c][3*k+2]=(float)(20*c);
        }
    }
    AsmField *f=AsmField_new(arena,&run);if(!f){Arena_dispose(&arena);return 1;}
    double saved[24];memcpy(saved,f->uv,sizeof saved);uint8_t fixed[]={0,1,1};
    AsmRepairOpts opts;AsmRepair_defaults(&opts);opts.compact=1;opts.fixed_charts=fixed;
    opts.iterations=opts.closure_iterations=0;opts.clearance_iterations=1;
    AsmContactStats before={0},after={0};int failures=0;
    if(AsmContacts_measure(f->contacts,f->uv,NULL,NULL,NULL,&before) || !before.complete || !before.pairs)failures++;
    for(int method=3;method>=0;method-=3){
        memcpy(f->uv,saved,sizeof saved);opts.clearance_method=method;opts.deadline=ves_clock_sec()+2;
        AsmRepairStats stats={0};stats.first_metric_chart=-1;
        if(apr_solve(f,&opts,&stats) || AsmContacts_measure(f->contacts,f->uv,NULL,NULL,NULL,&after) || !after.complete ||
           memcmp(saved+8,f->uv+8,16*sizeof(double)))failures++;
        if(method==3){if(stats.clearance_steps || fabs(after.area-before.area)>1e-10)failures++;}
        else{
            if(!stats.clearance_steps || !(after.area<before.area) || fabs(f->uv[3]-f->uv[1])<1e-6)failures++;
            for(size_t t=0;t<2;t++){
                double j[4],lo,hi,det;AsmField_face(f,t,f->uv,j,&lo,&hi,&det);
                if(fabs(lo-1)>1e-10 || fabs(hi-1)>1e-10 || !(det>0))failures++;
            }
        }
    }
    fprintf(stderr,"  contact torque, translation cancellation, proper rotation and fixed exterior: %s (%d failures, %.9g -> %.9g area)\n",
            failures?"FAIL":"ok",failures,before.area,after.area);
    Arena_dispose(&arena);return failures;
}

static int app_selftest(void)
{
    enum {N=9,NV=N*N,NF=2*(N-1)*(N-1)};
    int failures=apr_contact_pose_control();
    for(int obstructed=0;obstructed<2;obstructed++){
        Arena_T arena=Arena_new();AsmRun run={0};run.arena=arena;
        AsmChart charts[3]={{0}};float xyz[3][3*NV],uv[3][2*NV];int32_t faces[3*NF];
        for(int y=0;y<N;y++)for(int x=0;x<N;x++){
            int v=y*N+x;
            for(int c=0;c<3;c++){xyz[c][3*v]=(float)(x+10*c);xyz[c][3*v+1]=(float)y;xyz[c][3*v+2]=0;uv[c][2*v]=(float)x;uv[c][2*v+1]=(float)y;}
            if(x+1<N && y+1<N){int32_t tri[]={v,v+1,v+N+1,v,v+N+1,v+N};memcpy(faces+6*(y*(N-1)+x),tri,sizeof tri);}
        }
        run.charts=charts;run.n_charts=3;
        for(int c=0;c<3;c++){
            charts[c].id=c;charts[c].parent=-1;charts[c].nv=NV;charts[c].nf=NF;charts[c].xyz=xyz[c];charts[c].uv=uv[c];charts[c].faces=faces;
            charts[c].placed=1;charts[c].placement_state=ASM_PLACE_ROOT;charts[c].pose_x=c==1?11.5:c==2?(obstructed?2.5:1000):0;
        }
        AsmRelation relation={0};relation.a=0;relation.b=1;relation.corr_count=N;relation.n_corr=N;relation.continuity=ASM_CONT_SOURCE;
        AsmCorr corr[N]={{0}};run.rels=&relation;run.n_rels=1;run.corr=corr;run.n_corr=N;
        for(int y=0;y<N;y++){
            AsmCorr *c=corr+y;c->va=y*N+N-1;c->vb=y*N;c->valid=3;c->gap_a[0]=c->gap_b[0]=2;
            c->trim_a=apr_test_trim(charts,c->va);c->trim_b=apr_test_trim(charts+1,c->vb);
        }
        apr_test_scale(&run,16);
        AsmField *f=AsmField_new(arena,&run);if(!f){failures++;Arena_dispose(&arena);continue;}
        double saved[6*NV];memcpy(saved,f->uv,sizeof saved);uint8_t core[]={0,1,0};size_t pending=0;
        /* A real failed source seam connects two 16384-square-voxel pieces.
         * The opt-in queue recognizes this without changing any geometry;
         * an unselected hypothesis must retain its ordinary priority. */
        ApcOrder source_order[1];
        if(app_source_order(f,arena,source_order,0)!=1 || fabs(source_order[0].error-88)>1e-8 ||
           app_source_order(f,arena,source_order,1)!=1 || fabs(source_order[0].error+16384)>1e-8 ||
           memcmp(saved,f->uv,sizeof saved))failures++;
        f->seams[0].required=0;
        if(app_source_order(f,arena,source_order,1)!=1 || fabs(source_order[0].error-24)>1e-8)failures++;
        f->seams[0].required=1;
        AsmRepairOpts opts;AsmRepair_defaults(&opts);opts.compact=1;opts.iterations=8;opts.clearance_iterations=0;opts.closure_iterations=3;
        AsmRepairStats stats={0};AppResult result;
        opts.patch_faces=NF-1;
        if(app_attempt(f,core,&pending,1,&opts,&result,&stats) || !result.oversize || result.kept || memcmp(saved,f->uv,sizeof saved))failures++;
        opts.patch_faces=NF;opts.patch_field_faces=2*NF;opts.deadline=ves_clock_sec()-1;
        if(app_attempt(f,core,&pending,1,&opts,&result,&stats) || !result.expired || result.kept || memcmp(saved,f->uv,sizeof saved))failures++;
        opts.deadline=0;
        if(app_attempt(f,core,&pending,1,&opts,&result,&stats))failures++;
        if(obstructed){
            if(!stats.rejected_contacts || memcmp(saved,f->uv,2*NV*sizeof(double)) || memcmp(saved+4*NV,f->uv+4*NV,2*NV*sizeof(double)))failures++;
            if(result.kept){
                /* A shortened step can retain shallow overlap. For these
                 * aligned right triangles the SAT band allows tolerance*sqrt(2). */
                for(size_t k=0;k<NV;k++)if(f->uv[2*(NV+k)]<168-AsmContacts_material_tolerance()*sqrt(2.0)-1e-7)failures++;
            }else if(memcmp(saved,f->uv,sizeof saved) || relation.continuity!=ASM_CONT_SOURCE)failures++;
        }else{
            AsmSeamMeasure seam;AsmField_seam(f,0,f->uv,&seam);
            if(!result.kept || !seam.pass || result.field_faces!=2*NF || memcmp(saved,f->uv,2*NV*sizeof(double)) ||
               memcmp(saved+4*NV,f->uv+4*NV,2*NV*sizeof(double)) || !(relation.continuity&ASM_CONT_REQUIRED))failures++;
            if(app_source_order(f,arena,source_order,1))failures++;
        }
        Arena_dispose(&arena);
    }
    for(int area_only=0;area_only<2;area_only++){
        /* A fully included component has no physical fixed source anchor.
         * Numerical gauge pins must not stop it moving away from an obstacle.
         * The regional neighborhood also has to include leaf-tree partners. */
        Arena_T arena=Arena_new();AsmRun run={0};run.arena=arena;
        float xyz[]={0,0,0,2,0,0,0,2,0},uv[]={0,0,2,0,0,2};int32_t faces[]={0,1,2};AsmChart charts[3]={{0}};
        for(int c=0;c<3;c++){charts[c].id=c;charts[c].parent=-1;charts[c].nv=3;charts[c].nf=1;charts[c].xyz=xyz;charts[c].uv=uv;charts[c].faces=faces;
            charts[c].placed=1;charts[c].placement_state=ASM_PLACE_ROOT;charts[c].pose_x=c==2?-2.1:.5*c;}
        run.charts=charts;run.n_charts=3;apr_test_scale(&run,16);AsmField *f=AsmField_new(arena,&run);
        if(!f)failures++;else{
            uint8_t core[]={1,0,0};double fixed[12];memcpy(fixed,f->uv+6,sizeof fixed);
            AsmRepairOpts opts;AsmRepair_defaults(&opts);opts.compact=1;opts.iterations=area_only?0:8;
            opts.clearance_iterations=area_only?12:0;opts.closure_iterations=0;opts.clearance_method=area_only?3:0;
            if(area_only){
                /* A two-coordinate translation plus its contact slack must
                 * not enter a factorization with a two-unknown cap. */
                uint8_t fixed_charts[]={0,1,1};double saved[18];memcpy(saved,f->uv,sizeof saved);
                AsmRepairOpts capped=opts;capped.fixed_charts=fixed_charts;capped.clearance_method=2;capped.clearance_iterations=1;capped.patch_coordinates=2;
                AsmRepairStats check={0};check.first_metric_chart=-1;
                if(apr_solve(f,&capped,&check) || check.clearance_steps || memcmp(saved,f->uv,sizeof saved))failures++;
            }
            AsmRepairStats stats={0};AppResult result;
            /* The area proposal can approach the new obstacle before its
             * expansion would add a second pair. Retain the earlier useful
             * state when that later expansion is correctly refused. */
            if(app_attempt(f,core,NULL,0,&opts,&result,&stats) || !result.kept || result.field_faces>3 || !result.expansions ||
               result.after.pairs>result.before.pairs || !(result.after.area<(area_only?1:.5)*result.before.area) ||
               memcmp(fixed,f->uv+6,sizeof fixed))failures++;
            AsmContactStats complete={0};
            if(AsmContacts_measure(f->contacts,f->uv,NULL,NULL,NULL,&complete) || !complete.complete ||
               complete.pairs!=result.after.pairs || fabs(complete.area-result.after.area)>1e-12)failures++;
            AsmFieldSeam seam={0};seam.a=0;seam.b=1;seam.required=seam.complete=1;f->seams=&seam;f->n_seams=1;
            uint8_t scratch[3];app_neighborhood(f,0,core,1,scratch);if(!core[0] || !core[1] || core[2])failures++;
        }
        Arena_dispose(&arena);
    }
    fprintf(stderr,"  compact repair: source closure, exact exterior, new obstacle guard, free component pose, resource/deadline deferral: %s (%d failures)\n",failures?"FAIL":"ok",failures);
    return failures;
}
