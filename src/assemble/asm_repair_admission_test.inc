static int apr_overlay_contact_control(void)
{
    int fails=0;
    for(int batched=0;batched<2;batched++){
    enum {NC=9};Arena_T arena=Arena_new();AsmRun run={0};run.arena=arena;run.continuity_ready=1;
    float uv[]={0,0,2,0,2,2,0,2},xyz[]={0,0,0,2,0,0,2,2,0,0,2,0};
    int32_t faces[]={0,1,2,0,2,3};uint8_t core[NC]={0},member[NC]={0};
    for(int c=0;c<NC;c++){
        AsmChart chart={0};chart.id=c;chart.parent=-1;chart.nv=4;chart.nf=2;
        chart.xyz=xyz;chart.uv=uv;chart.faces=faces;chart.area3d=4;
        chart.placed=1;chart.placement_state=ASM_PLACE_ROOT;chart.pose_x=10*c;
        if(c%2)chart.flags|=ASM_CHART_MIRROR;
        AsmRun_push_chart(&run,&chart);
    }
    apr_test_scale(&run,8);
    AsmField *field=AsmField_new(arena,&run);int okay=field!=NULL;
    ApaOverlay overlay={0};overlay.arena=arena;overlay.batched=batched;
    for(int c=0;c<NC && okay;c++){
        core[c]=1;okay=!apa_overlay_add(&overlay,&run,core);core[c]=0;
    }
    /* Probe every block, including mirrored faces. Geometrically identical
     * external triangles must be seen; excluding their owner must remove
     * precisely those contacts, never unrelated stored material. */
    for(int c=0;c<NC && okay;c++){
        core[c]=1;size_t pairs=0;double area=0;int complete=0;
        uint8_t colliders[NC]={0};
        okay=!apa_external_query(field,core,member,NULL,&overlay,arena,&pairs,&area,&complete,0,colliders) &&
            complete && pairs==2 && fabs(area-256)<1e-12;
        for(int k=0;k<NC;k++)okay=okay && colliders[k]==(k==c);
        member[c]=1;
        memset(colliders,0,sizeof colliders);
        okay=okay && !apa_external_query(field,core,member,NULL,&overlay,arena,&pairs,&area,&complete,0,colliders) &&
            complete && !pairs && !area;
        for(int k=0;k<NC;k++)okay=okay && !colliders[k];
        member[c]=0;core[c]=0;
    }
    if(okay && batched){
        AsmContacts *small=overlay.level[0]->tree;
        size_t saved=AsmContacts_budget(small);AsmContacts_leaf_pair_budget(small,1);
        core[NC-1]=1;size_t pairs=0;double area=0;int complete=0;
        okay=!apa_external_pairs(field,core,member,NULL,&overlay,arena,&pairs,&area,&complete,0) &&
            complete && pairs==2 && fabs(area-256)<1e-12 && AsmContacts_budget(small)==1;
        core[NC-1]=0;AsmContacts_leaf_pair_budget(small,saved);
    }
    /* Exhaust a real external query. Partial contacts must remain visibly
     * incomplete, including when the shared allowance spans several blocks. */
    if(okay){
        memset(core,1,sizeof core);
        if(batched){for(int i=0;i<APA_OVERLAY_LEVELS;i++)if(overlay.level[i])
            AsmContacts_leaf_pair_budget(overlay.level[i]->tree,1);}
        else AsmContacts_leaf_pair_budget(overlay.tree,1);
        size_t pairs=0;double area=0;int complete=1;
        okay=!apa_external_pairs(field,core,member,NULL,&overlay,arena,&pairs,&area,&complete,0) && !complete;
        if(batched){for(int i=0;i<APA_OVERLAY_LEVELS;i++)if(overlay.level[i])
            okay=okay && AsmContacts_budget(overlay.level[i]->tree)==1;}
    }
    fprintf(stderr,"  native admission overlay contacts and exclusions, batched %d: %s\n",batched,okay?"ok":"FAIL");
    fails+=!okay;apa_overlay_clear(&overlay);Arena_dispose(&arena);
    }
    return fails;
}

static int apr_overlay_replacement_control(void)
{
    enum {NC=4};Arena_T arena=Arena_new();AsmRun run={0};run.arena=arena;run.continuity_ready=1;
    float uv[]={0,0,20,0,20,20,0,20},xyz[]={0,0,0,20,0,0,20,20,0,0,20,0};
    int32_t faces[]={0,1,2,0,2,3};uint8_t changed[NC]={0},moving[NC]={0},core[NC]={0,0,0,1},member[NC]={0,0,0,1};
    for(int c=0;c<NC;c++){
        AsmChart chart={0};chart.id=c;chart.parent=-1;chart.nv=4;chart.nf=2;
        chart.xyz=xyz;chart.uv=uv;chart.faces=faces;chart.area3d=400;
        chart.placed=c!=2;chart.placement_state=c==2?ASM_PLACE_NONE:ASM_PLACE_ROOT;chart.pose_x=100*c;
        if(c==1)chart.flags|=ASM_CHART_MIRROR;
        AsmRun_push_chart(&run,&chart);
    }
    AsmField *baseline=AsmField_new(arena,&run),*probe=AsmField_new(arena,&run);
    int okay=baseline && probe && !AsmContacts_refit(baseline->contacts,baseline->uv);
    ApaOverlay overlay={0};overlay.arena=arena;overlay.batched=overlay.query_scaled_budget=1;
    overlay.changed_charts=changed;overlay.n_charts=NC;
    for(int round=0;round<3 && okay;round++){
        memset(moving,0,sizeof moving);int c=round==1?1:0;moving[c]=1;
        run.charts[c].pose_x=round==0?50:round==1?150:400;
        if(round==1){run.charts[2].placed=1;run.charts[2].placement_state=ASM_PLACE_ROOT;moving[2]=1;}
        okay=!apa_overlay_add(&overlay,&run,moving) && overlay.nf==(round?6:2);
        AsmField *fresh=AsmField_new(arena,&run);
        okay=okay && fresh && !AsmContacts_refit(fresh->contacts,fresh->uv);
        /* Compare both old and new locations against a separately rebuilt
         * full index, including mirror parity, a newly admitted chart, and
         * the second replacement of a chart already in the delta tree. */
        for(int x=0;x<=400 && okay;x+=25){
            size_t offset=probe->vertex_offset[3];
            for(size_t v=0;v<4;v++){probe->uv[2*(offset+v)]=x+uv[2*v];probe->uv[2*(offset+v)+1]=uv[2*v+1];}
            size_t a=0,b=0;double aa=0,ba=0;int ac=0,bc=0;uint8_t ca[NC]={0},cb[NC]={0};
            okay=!AsmContacts_refit(probe->contacts,probe->uv) &&
                !apa_external_query(probe,core,member,baseline,&overlay,arena,&a,&aa,&ac,1,ca) &&
                !apa_external_query(probe,core,member,fresh,NULL,arena,&b,&ba,&bc,1,cb) &&
                ac && bc && a==b && fabs(aa-ba)<1e-9 && !memcmp(ca,cb,sizeof ca);
        }
    }
    if(okay){
        ApaOverlayBatch *saved=overlay.level[0];size_t old_nf=overlay.nf,nv=run.charts[0].nv;
        uint8_t old_changed[NC];memcpy(old_changed,changed,sizeof changed);
        run.charts[0].nv=(size_t)INT_MAX+1;moving[0]=1;
        okay=apa_overlay_add(&overlay,&run,moving)==-1 && overlay.level[0]==saved &&
            overlay.nf==old_nf && !memcmp(old_changed,changed,sizeof changed);
        run.charts[0].nv=nv;
    }
    fprintf(stderr,"  replacement contact index versus full rebuild, repeated motion and rollback: %s\n",okay?"ok":"FAIL");
    apa_overlay_clear(&overlay);Arena_dispose(&arena);return !okay;
}

/* The incremental replacement index (ASM_REPAIR_OVERLAY_INCREMENTAL) answers
 * every probe exactly as the full rebuild and as a freshly built full index,
 * through 40 rounds of repeated motion, admissions, mirrored charts, stale
 * copies in older blocks and their compaction, while rebuilding fewer faces;
 * a refused build changes nothing. */
static int apr_overlay_incremental_control(void)
{
    enum {NC=13,PROBE=NC-1,ROUNDS=40};Arena_T arena=Arena_new();AsmRun run={0};run.arena=arena;run.continuity_ready=1;
    float uv[]={0,0,20,0,20,20,0,20},xyz[]={0,0,0,20,0,0,20,20,0,0,20,0};
    int32_t faces[]={0,1,2,0,2,3};uint8_t changed[2][NC]={{0}},moving[NC]={0},core[NC]={0},member[NC]={0};
    core[PROBE]=member[PROBE]=1;
    for(int c=0;c<NC;c++){
        AsmChart chart={0};chart.id=c;chart.parent=-1;chart.nv=4;chart.nf=2;
        chart.xyz=xyz;chart.uv=uv;chart.faces=faces;chart.area3d=400;
        chart.placed=c<9 || c==PROBE;chart.placement_state=chart.placed?ASM_PLACE_ROOT:ASM_PLACE_NONE;
        chart.pose_x=c==PROBE?5000:100*c;
        if(c%3==1)chart.flags|=ASM_CHART_MIRROR;
        AsmRun_push_chart(&run,&chart);
    }
    AsmField *baseline=AsmField_new(arena,&run),*probe=AsmField_new(arena,&run);
    int okay=baseline && probe && !AsmContacts_refit(baseline->contacts,baseline->uv);
    ApaOverlay ov[2];
    for(int m=0;m<2;m++){
        memset(ov+m,0,sizeof ov[m]);ov[m].arena=arena;ov[m].batched=ov[m].query_scaled_budget=1;
        ov[m].changed_charts=changed[m];ov[m].n_charts=NC;ov[m].incremental=m;
    }
    uint32_t seed=12345u;size_t contacts=0,stale_seen=0;
    for(int round=0;round<ROUNDS && okay;round++){
        memset(moving,0,sizeof moving);
        int count=1+(int)((seed=seed*1664525u+1013904223u)>>29)%3;
        for(int k=0;k<count;k++){
            int c=(int)((seed=seed*1664525u+1013904223u)>>16)%PROBE;
            if(!run.charts[c].placed){run.charts[c].placed=1;run.charts[c].placement_state=ASM_PLACE_ROOT;}
            run.charts[c].pose_x=(double)(((seed=seed*1664525u+1013904223u)>>16)%1200);moving[c]=1;
        }
        for(int m=0;m<2 && okay;m++)okay=!apa_overlay_add(ov+m,&run,moving);
        okay=okay && ov[0].nf==ov[1].nf && !memcmp(changed[0],changed[1],NC);
        size_t stored=0;for(int i=0;i<APA_OVERLAY_LEVELS;i++)if(ov[1].level[i])stored+=ov[1].level[i]->nf;
        if(stored>ov[1].nf)stale_seen++;
        AsmField *fresh=AsmField_new(arena,&run);
        okay=okay && fresh && !AsmContacts_refit(fresh->contacts,fresh->uv);
        for(int x=0;x<=1225 && okay;x+=25){
            size_t offset=probe->vertex_offset[PROBE];
            for(size_t v=0;v<4;v++){probe->uv[2*(offset+v)]=x+uv[2*v];probe->uv[2*(offset+v)+1]=uv[2*v+1];}
            size_t pairs[3]={0};double area[3]={0};int complete[3]={0};uint8_t colliders[3][NC]={{0}};
            okay=!AsmContacts_refit(probe->contacts,probe->uv);
            for(int m=0;m<3 && okay;m++)
                okay=!apa_external_query(probe,core,member,m<2?baseline:fresh,m<2?ov+m:NULL,arena,pairs+m,area+m,complete+m,1,colliders[m]);
            okay=okay && complete[0] && complete[1] && complete[2] && pairs[0]==pairs[2] && pairs[1]==pairs[2] &&
                fabs(area[0]-area[2])<1e-9 && fabs(area[1]-area[2])<1e-9 &&
                !memcmp(colliders[0],colliders[2],NC) && !memcmp(colliders[1],colliders[2],NC);
            contacts+=pairs[2];
        }
    }
    okay=okay && contacts>0 && stale_seen>0 && ov[1].compactions>0 && ov[1].rebuilt_faces<ov[0].rebuilt_faces;
    if(okay){
        ApaOverlayBatch *saved[APA_OVERLAY_LEVELS];memcpy(saved,ov[1].level,sizeof saved);
        uint8_t blocks[NC],before[NC];memcpy(blocks,ov[1].chart_block,NC);memcpy(before,changed[1],NC);
        size_t nf=ov[1].nf,nv=run.charts[0].nv;
        run.charts[0].nv=(size_t)INT_MAX+1;memset(moving,0,sizeof moving);moving[0]=1;
        okay=apa_overlay_add(ov+1,&run,moving)==-1 && !memcmp(saved,ov[1].level,sizeof saved) &&
            !memcmp(blocks,ov[1].chart_block,NC) && !memcmp(before,changed[1],NC) && ov[1].nf==nf;
        run.charts[0].nv=nv;
    }
    fprintf(stderr,"  incremental replacement index versus full rebuild and fresh index, %d rounds: %s "
            "(%zu probe contacts, %zu rounds with stale copies, %zu compactions, rebuilt faces %zu vs %zu)\n",
            ROUNDS,okay?"ok":"FAIL",contacts,stale_seen,ov[1].compactions,ov[1].rebuilt_faces,ov[0].rebuilt_faces);
    for(int m=0;m<2;m++)apa_overlay_clear(ov+m);
    Arena_dispose(&arena);return !okay;
}

/* A replacement block built on four workers is the serial block: the same
 * tree size and the same external answer, bit for bit, for a candidate over
 * an 8,192-face overlay block (past the parallel build's threshold). */
static int apr_overlay_threads_control(void)
{
    enum {N=65,NV=N*N,NF=2*(N-1)*(N-1)};
    Arena_T arena=Arena_new();AsmRun run={0};run.arena=arena;run.continuity_ready=1;
    float coarse_uv[]={0,0,1,0,1,1,0,1},coarse_xyz[]={0,0,0,1,0,0,1,1,0,0,1,0};
    int32_t coarse_faces[]={0,1,2,0,2,3};
    float *uv=ARENA_ALLOC(arena,2*NV*sizeof(float)),*xyz=ARENA_ALLOC(arena,3*NV*sizeof(float));
    int32_t *faces=ARENA_ALLOC(arena,3*NF*sizeof(int32_t));
    for(int y=0;y<N;y++)for(int x=0;x<N;x++){
        int v=y*N+x;uv[2*v]=xyz[3*v]=(float)x/(N-1);uv[2*v+1]=xyz[3*v+1]=(float)y/(N-1);xyz[3*v+2]=0;
        if(x+1<N && y+1<N){int32_t t[]={v,v+1,v+N+1,v,v+N+1,v+N};memcpy(faces+6*(y*(N-1)+x),t,sizeof t);}
    }
    for(int c=0;c<3;c++){
        AsmChart ch={0};ch.id=c;ch.parent=-1;ch.nv=c==2?NV:4;ch.nf=c==2?NF:2;
        ch.xyz=c==2?xyz:coarse_xyz;ch.uv=c==2?uv:coarse_uv;ch.faces=c==2?faces:coarse_faces;
        ch.area3d=1;ch.placed=1;ch.placement_state=ASM_PLACE_ROOT;AsmRun_push_chart(&run,&ch);
    }
    apr_test_scale(&run,1024);
    uint8_t core[3]={1,0,0},member[3]={1,0,0},add[3]={0,1,1};
    AsmField *field=AsmField_new_local(arena,&run,member,NULL);
    int okay=field && !AsmField_build_contacts(field);
    size_t pairs[2]={0},leaves[2]={0};double area[2]={0};int complete[2]={0};
    for(int k=0;k<2 && okay;k++){
        ApaOverlay ov={0};ov.arena=arena;ov.batched=ov.query_scaled_budget=1;ov.threads=k?4:1;
        okay=!apa_overlay_add(&ov,&run,add) && ov.level[0] && ov.level[0]->nf==(size_t)NF+2 &&
            !apa_external_pairs(field,core,member,NULL,&ov,arena,pairs+k,area+k,complete+k,1);
        if(okay)leaves[k]=AsmContacts_leaves(ov.level[0]->tree);
        apa_overlay_clear(&ov);
    }
    okay=okay && complete[0] && complete[1] && pairs[0]>0 && pairs[0]==pairs[1] &&
         !memcmp(area,area+1,sizeof *area) && leaves[0]==leaves[1];
    fprintf(stderr,"  replacement block built on 4 workers equals the serial block: %s (%zu pairs, %zu leaves)\n",
            okay?"ok":"FAIL",pairs[0],leaves[0]);
    Arena_dispose(&arena);return !okay;
}

static int apr_overlay_query_budget_control(void)
{
    enum {N=65,NV=N*N,NF=2*(N-1)*(N-1)};int failures=0;
    for(int blocks=1;blocks<=3;blocks++){
        Arena_T arena=Arena_new();AsmRun run={0};run.arena=arena;run.continuity_ready=1;
        float coarse_uv[]={0,0,1,0,1,1,0,1},coarse_xyz[]={0,0,0,1,0,0,1,1,0,0,1,0};
        int32_t coarse_faces[]={0,1,2,0,2,3};
        float *uv=ARENA_ALLOC(arena,2*NV*sizeof(float)),*xyz=ARENA_ALLOC(arena,3*NV*sizeof(float));
        int32_t *faces=ARENA_ALLOC(arena,3*NF*sizeof(int32_t));
        for(int y=0;y<N;y++)for(int x=0;x<N;x++){
            int v=y*N+x;uv[2*v]=xyz[3*v]=(float)x/(N-1);uv[2*v+1]=xyz[3*v+1]=(float)y/(N-1);xyz[3*v+2]=0;
            if(x+1<N && y+1<N){int32_t t[]={v,v+1,v+N+1,v,v+N+1,v+N};memcpy(faces+6*(y*(N-1)+x),t,sizeof t);}
        }
        for(int c=0;c<=blocks;c++){
            AsmChart ch={0};ch.id=c;ch.parent=-1;ch.nv=c==blocks?NV:4;ch.nf=c==blocks?NF:2;
            ch.xyz=c==blocks?xyz:coarse_xyz;ch.uv=c==blocks?uv:coarse_uv;ch.faces=c==blocks?faces:coarse_faces;
            ch.area3d=1;ch.placed=1;ch.placement_state=ASM_PLACE_ROOT;AsmRun_push_chart(&run,&ch);
        }
        apr_test_scale(&run,1024);
        uint8_t core[4]={0},member[4]={0};core[blocks]=member[blocks]=1;
        AsmField *field=AsmField_new_local(arena,&run,member,NULL);
        ApaOverlay ov={0};ov.arena=arena;ov.batched=1;int okay=field && !AsmField_build_contacts(field);
        uint8_t add[4]={0};
        for(int c=0;c<blocks && okay;c++){add[c]=1;okay=!apa_overlay_add(&ov,&run,add);add[c]=0;}
        size_t pairs=0;double area=0;int complete=1;
        if(okay)okay=!apa_external_pairs(field,core,member,NULL,&ov,arena,&pairs,&area,&complete,1) && !complete;
        ov.query_scaled_budget=1;
        size_t saved[APA_OVERLAY_LEVELS]={0};
        for(int i=0;i<APA_OVERLAY_LEVELS;i++)if(ov.level[i])saved[i]=AsmContacts_budget(ov.level[i]->tree);
        for(int repeat=0;repeat<2 && okay;repeat++){
            okay=!apa_external_pairs(field,core,member,NULL,&ov,arena,&pairs,&area,&complete,1) &&
                complete && pairs==(size_t)blocks*NF && fabs(area-1048576.0*blocks)<1e-10;
            for(int i=0;i<APA_OVERLAY_LEVELS;i++)if(ov.level[i])okay=okay && saved[i]==AsmContacts_budget(ov.level[i]->tree);
        }
        /* Exclusion remains exact under the larger query allowance. */
        member[0]=1;
        if(okay)okay=!apa_external_pairs(field,core,member,NULL,&ov,arena,&pairs,&area,&complete,1) &&
            complete && pairs==(size_t)(blocks-1)*NF && fabs(area-1048576.0*(blocks-1))<1e-10;
        size_t cap=(size_t)ASM_CONTACT_MAX_LEAF_PAIRS;
        okay=okay && apa_overlay_probe_budget(&ov,SIZE_MAX)==cap && apa_overlay_probe_budget(&ov,NF)<=cap;
        fprintf(stderr,"  native asymmetric overlay query, %d charts, %d incoming triangles: %s\n",blocks,NF,okay?"ok":"FAIL");
        failures+=!okay;apa_overlay_clear(&ov);Arena_dispose(&arena);
    }
    return failures;
}

static int apr_frontier_obligation_control(void)
{
    enum {N=9,NV=N*N,NF=2*(N-1)*(N-1)};int failures=0;
    for(int test=0;test<21;test++){
        int incompatible_required=test==1 || test==8,compatible_required=test==7 || test==9;
        int alternative=test>=10 && test<=13,obstacle=test>=14;
        Arena_T arena=Arena_new();AsmRun run={0};run.arena=arena;
        float xyz[4][3*NV],uv[4][2*NV];int32_t faces[3*NF];uint8_t boundary[NV];
        for(int y=0;y<N;y++)for(int x=0;x<N;x++){
            int v=y*N+x;boundary[v]=(uint8_t)(!x || !y || x==N-1 || y==N-1);
            for(int c=0;c<4;c++){
                xyz[c][3*v]=(float)(10*c+x);xyz[c][3*v+1]=(float)y;xyz[c][3*v+2]=(float)(c==3?30:0);
                uv[c][2*v]=(float)x;uv[c][2*v+1]=(float)y;
            }
            if(x+1<N && y+1<N){int32_t tri[]={v,v+1,v+N+1,v,v+N+1,v+N};memcpy(faces+6*(y*(N-1)+x),tri,sizeof tri);}
        }
        for(int c=0;c<4;c++){
            AsmChart ch={0};ch.id=c;ch.parent=-1;ch.nv=NV;ch.nf=NF;ch.xyz=xyz[c];ch.faces=faces;
            ch.boundary=boundary;ch.uv=uv[c];ch.area3d=64;
            if(c!=1){ch.placed=1;ch.placement_state=ASM_PLACE_ROOT;ch.pose_y=200;
                ch.pose_x=c==0?100:c==2?(compatible_required?120:1000):obstacle?(test==15?110:117.625):alternative?990:test==6?110:10000;
                if(c==3 && (test==18 || test==19)){ch.pose_x=117.9999;ch.pose_y=207.9999;}
                if(test==8 || test==9){ch.pose_x+=1e12;ch.pose_y+=1e12;}}
            AsmRun_push_chart(&run,&ch);
        }
        for(int i=0;i<2;i++){
            AsmRelation r={0};r.a=i;r.b=i+1;r.flags=i?ASM_REL_WEAK:0;
            if(!i && (test==2 || test==3 || alternative))r.flags=ASM_REL_WEAK;
            if(!i && test==3)r.flags|=ASM_REL_CROSSWRAP;
            if(!i && test==12)r.flags|=ASM_REL_CROSSWRAP;
            if(!i && test==13)r.flags|=ASM_REL_DROPPED;
            if(i && (test==2 || test==3))r.flags|=ASM_REL_DROPPED;
            if(i && test==5)r.flags|=ASM_REL_PARITY;
            if(alternative)r.seam_len=(float)(i+1);
            r.corr_first=(int32_t)run.n_corr;r.corr_count=r.n_corr=N;
            for(int y=0;y<N;y++){AsmCorr corr={0};corr.va=N*y+N-1;corr.vb=N*y;AsmRun_push_corr(&run,&corr);}
            AsmRun_push_rel(&run,&r);
        }
        if(AsmContinuity_prepare(&run)){Arena_dispose(&arena);failures++;continue;}
        if(incompatible_required || compatible_required || test==11)run.rels[1].continuity|=ASM_CONT_REQUIRED;
        if(test==4)run.corr[0].valid=0;
        if(test==6 || (test>=10 && test<=17) || test==20)apr_test_scale(&run,16);
        AsmChart before[4];memcpy(before,run.charts,sizeof before);AsmRelation deferred=run.rels[1];
        AsmCorr observations[2*N];memcpy(observations,run.corr,sizeof observations);
        ApaBox boxes[4];apa_chart_boxes(&run,boxes);AsmField *incumbent=AsmField_new(arena,&run);
        if(!incumbent || AsmContacts_refit(incumbent->contacts,incumbent->uv)){Arena_dispose(&arena);failures++;continue;}
        ApaOverlay overlay={0};overlay.arena=arena;overlay.batched=1;uint8_t core[]={0,1,0,0};
        AsmRepairOpts opts;AsmRepair_defaults(&opts);opts.compact=1;
        opts.iterations=opts.clearance_iterations=opts.closure_iterations=0;opts.admission_rounds=1;
        AsmRepairStats stats={0};size_t anchors=0;const char *reason="";ApaReport report={0};
        int all=compatible_required?0:apa_attempt(&run,&opts,core,&anchors,&stats,&reason,NULL,boxes,incumbent,&overlay,&report,1);
        int okay=!all && !memcmp(before,run.charts,sizeof before) && !memcmp(&deferred,run.rels+1,sizeof deferred);
        if(incompatible_required || compatible_required)opts.iterations=1;
        if(test==19)opts.deadline=ves_clock_sec()-1;
        int selected=apa_attempt(&run,&opts,core,&anchors,&stats,&reason,NULL,boxes,incumbent,&overlay,&report,2);
        if(alternative){
            /* The strongest deferred anchor lands on an external obstacle.
             * A weaker clear anchor can succeed, but cannot replace an
             * existing obligation or bypass a source-relation veto. */
            okay=okay && !selected && !memcmp(before,run.charts,sizeof before) && !stats.admitted_charts;
            selected=apa_attempt_anchor(&run,&opts,core,&anchors,&stats,&reason,NULL,boxes,incumbent,&overlay,&report,2,0);
        }
        if(obstacle && test<18){
            okay=okay && !selected && !memcmp(before,run.charts,sizeof before) && !stats.admitted_charts;
            opts.iterations=opts.clearance_iterations=2;opts.clearance_method=test==17?0:3;
            if(test==16)opts.patch_field_faces=2*NF;
            selected=apa_attempt(&run,&opts,core,&anchors,&stats,&reason,NULL,boxes,incumbent,&overlay,&report,3);
            if(test==16)okay=okay && !strcmp(reason,"obstacle_resource_cap");
        }
        /* Shallow contact now succeeds directly. Deadline refusal still
         * precedes admission, and a wider unresolved edge stays blocked. */
        if(test==18)okay=okay && selected && !report.external_pairs;
        if(test==19)okay=okay && !strcmp(reason,"shared_deadline");
        if(test==20){
            okay=okay && !selected && !memcmp(before,run.charts,sizeof before);
            selected=apa_attempt(&run,&opts,core,&anchors,&stats,&reason,NULL,boxes,incumbent,&overlay,&report,4);
        }
        int expect=test==0 || test==2 || test==5 || compatible_required || test==10 || test==14 || test==17 || test==18;
        okay=okay && selected==expect && !memcmp(before,run.charts,sizeof(AsmChart)) &&
            !memcmp(before+2,run.charts+2,2*sizeof(AsmChart)) && !memcmp(&deferred,run.rels+1,sizeof deferred) &&
            !memcmp(observations,run.corr,sizeof observations);
        if(expect)okay=okay && AsmChart_registered(run.charts+1) && (run.rels[0].continuity&ASM_CONT_REQUIRED);
        else okay=okay && !memcmp(before,run.charts,sizeof before) && !stats.admitted_charts;
        if(incompatible_required)okay=okay && !strcmp(reason,"incompatible_fixed_anchors") && !stats.iterations;
        fprintf(stderr,"  native retained frontier obligations %d: accepted %d (%s), %s\n",test,selected,reason,okay?"ok":"FAIL");
        failures+=!okay;apa_overlay_clear(&overlay);Arena_dispose(&arena);
    }
    return failures;
}

static int apr_frontier_admission_control(void)
{
    enum {N=9,NV=N*N,NF=2*(N-1)*(N-1)};int failures=0;
    for(int test=0;test<5;test++){
        Arena_T arena=Arena_new();AsmRun run={0};run.arena=arena;
        float xyz[4][3*NV],uv[4][2*NV];int32_t faces[3*NF];uint8_t boundary[NV];
        for(int y=0;y<N;y++)for(int x=0;x<N;x++){
            int v=y*N+x;boundary[v]=(uint8_t)(!x || !y || x==N-1 || y==N-1);
            for(int c=0;c<4;c++){
                xyz[c][3*v]=(float)(10*c+x);xyz[c][3*v+1]=(float)y;xyz[c][3*v+2]=(float)(c==3?20:0);
                uv[c][2*v]=(float)x;uv[c][2*v+1]=(float)y;
            }
            if(x+1<N && y+1<N){int32_t tri[]={v,v+1,v+N+1,v,v+N+1,v+N};memcpy(faces+6*(y*(N-1)+x),tri,sizeof tri);}
        }
        for(int c=0;c<4;c++){
            AsmChart ch={0};ch.id=c;ch.parent=-1;ch.nv=NV;ch.nf=NF;ch.xyz=xyz[c];ch.faces=faces;ch.boundary=boundary;ch.uv=uv[c];ch.area3d=64;
            if(c==0 || c==3){ch.placed=1;ch.placement_state=ASM_PLACE_ROOT;ch.pose_x=c?(test==1?120:10000):100;ch.pose_y=200;}
            AsmRun_push_chart(&run,&ch);
        }
        for(int i=0;i<2;i++){
            AsmRelation r={0};r.a=i;r.b=i+1;r.flags=i?ASM_REL_DROPPED|ASM_REL_WEAK:0;
            r.corr_first=(int32_t)run.n_corr;r.corr_count=r.n_corr=N;
            for(int y=0;y<N;y++){AsmCorr corr={0};corr.va=N*y+N-1;corr.vb=N*y;AsmRun_push_corr(&run,&corr);}
            AsmRun_push_rel(&run,&r);
        }
        if(AsmContinuity_prepare(&run)){Arena_dispose(&arena);failures++;continue;}
        if(test==3)run.rels[1].continuity|=ASM_CONT_REQUIRED;
        if(test==2)run.corr[N].valid=0;
        if(test==1)apr_test_scale(&run,8);
        AsmChart before[4];memcpy(before,run.charts,sizeof before);
        AsmRelation pending=run.rels[1],first_relation=run.rels[0];AsmCorr observations[2*N];memcpy(observations,run.corr,sizeof observations);
        ApaBox boxes[4];apa_chart_boxes(&run,boxes);AsmField *incumbent=AsmField_new(arena,&run);
        if(!incumbent || AsmContacts_refit(incumbent->contacts,incumbent->uv)){Arena_dispose(&arena);failures++;continue;}
        ApaOverlay overlay={0};overlay.arena=arena;overlay.batched=1;uint8_t core[]={0,1,0,0};
        if(test==4)overlay.nv=INT_MAX; /* force a bounded overlay-capacity failure */
        AsmRepairOpts opts;AsmRepair_defaults(&opts);opts.compact=1;opts.iterations=opts.clearance_iterations=opts.closure_iterations=0;opts.admission_rounds=1;
        AsmRepairStats stats={0};size_t anchors=0;const char *reason="";ApaReport report={0};
        int normal=apa_attempt(&run,&opts,core,&anchors,&stats,&reason,NULL,boxes,incumbent,&overlay,&report,0);
        int okay=!normal && !memcmp(before,run.charts,sizeof before);
        int first=apa_attempt(&run,&opts,core,&anchors,&stats,&reason,NULL,boxes,incumbent,&overlay,&report,1);
        if(test==4){
            okay=okay && !first && !strcmp(reason,"overlay_failed") && !stats.admitted_charts &&
                !memcmp(before,run.charts,sizeof before) && !memcmp(&first_relation,run.rels,sizeof first_relation) &&
                !memcmp(&pending,run.rels+1,sizeof pending) && !memcmp(observations,run.corr,sizeof observations);
            fprintf(stderr,"  native source-front failed-overlay rollback: %s\n",okay?"ok":"FAIL");
            failures+=!okay;apa_overlay_clear(&overlay);Arena_dispose(&arena);continue;
        }
        okay=okay && first && stats.admitted_charts==1 && AsmChart_registered(run.charts+1) &&
            !AsmChart_registered(run.charts+2) && !memcmp(&pending,run.rels+1,sizeof pending) &&
            !memcmp(before,run.charts,sizeof(AsmChart)) && !memcmp(before+3,run.charts+3,sizeof(AsmChart));
        AsmField *partial=AsmField_new(arena,&run);
        if(!partial)okay=0;
        else{
            int found=0;
            for(size_t i=0;i<partial->n_seams;i++)if(partial->seams[i].relation==1){
                AsmSeamMeasure measure;AsmField_seam(partial,i,partial->uv,&measure);
                found=1;if(partial->seams[i].complete || measure.pass)okay=0;
            }
            if(!found)okay=0;
        }
        apa_chart_box(run.charts+1,boxes+1);AsmChart retained=run.charts[1];
        core[1]=0;core[2]=1;
        int second=apa_attempt(&run,&opts,core,&anchors,&stats,&reason,NULL,boxes,incumbent,&overlay,&report,1);
        int expect=test==0 || test==3;
        okay=okay && second==expect && !memcmp(&retained,run.charts+1,sizeof retained) &&
            !memcmp(before,run.charts,sizeof(AsmChart)) && !memcmp(before+3,run.charts+3,sizeof(AsmChart)) &&
            !memcmp(observations,run.corr,sizeof observations);
        if(expect)okay=okay && AsmChart_registered(run.charts+2) && (run.rels[1].continuity&ASM_CONT_REQUIRED);
        else okay=okay && !memcmp(before+2,run.charts+2,sizeof(AsmChart)) && !memcmp(&pending,run.rels+1,sizeof pending);
        fprintf(stderr,"  native source-front admission %d: first %d second %d (%s), %s\n",test,first,second,reason,okay?"ok":"FAIL");
        failures+=!okay;apa_overlay_clear(&overlay);Arena_dispose(&arena);
    }
    return failures;
}

static int apr_factor_failure_control(void)
{
    /* A late nonpositive pivot must reject a factorization without reading
     * outside its supernode table or attempting a solve through a null factor. */
    int failures=0;const int sizes[]={8,65,128,260};
    for(size_t test=0;test<sizeof sizes/sizeof *sizes;test++){
        Arena_T arena=Arena_new();int n=sizes[test],nt=n*(n+1)/2,at=0;
        int *rows=ARENA_ALLOC(arena,nt*sizeof(int)),*cols=ARENA_ALLOC(arena,nt*sizeof(int));
        double *values=ARENA_ALLOC(arena,nt*sizeof(double)),*rhs=ARENA_CALLOC(arena,n,sizeof(double)),*x=ARENA_CALLOC(arena,n,sizeof(double));
        for(int i=0;i<n;i++)for(int j=0;j<=i;j++){
            rows[at]=i;cols[at]=j;values[at++]=i==j?(i==n-1?-1:1):.001;
        }
        AsmContinuityHalfspace plane={0};plane.count=1;plane.ids[0]=0;plane.values[0]=1;plane.rhs=-1;
        fprintf(stderr,"  native factor failure control: n=%d, late nonpositive pivot\n",n);
        if(AsmContinuity_project_step(arena,n,nt,rows,cols,values,rhs,x,&plane,1))failures++;
        Arena_dispose(&arena);
    }
    return failures;
}

static int apr_admission_forest_control(void)
{
    /* A weaker inconsistent cycle must not pull the initial forest away
     * from two stronger source paths. It remains an audit obligation. */
    AsmRun run={0};AsmChart charts[3]={{0}};run.charts=charts;run.n_charts=3;
    AsmField field={0};field.run=&run;field.n_seams=3;field.n_witness=6;
    double uv[]={100,200,104,200,100,204, 0,0,4,0,0,4, 0,0,4,0,0,4};
    int32_t labels[]={0,0,0,1,1,1,2,2,2};size_t offsets[]={0,3,6};
    AsmFieldSeam seams[3]={{0}};AsmWitness witness[6]={{0}};
    field.uv=uv;field.vertex_chart=labels;field.vertex_offset=offsets;field.seams=seams;field.witness=witness;
    const int pair[][2]={{0,1},{0,2},{1,2}};const double gap[]={13,20,10};
    for(int c=0;c<3;c++)charts[c].nv=3;
    for(int i=0;i<3;i++){
        AsmFieldSeam *s=seams+i;s->a=pair[i][0];s->b=pair[i][1];s->first=2*i;s->count=2;s->complete=s->required=1;
        for(int k=0;k<2;k++){
            AsmWitness *w=witness+2*i+k;w->valid=w->required=1;w->weight=i?1:.25;
            for(int side=0;side<2;side++){
                int a=pair[i][side],b=pair[i][1-side];double along=(side?-1:1)*gap[i]/4;
                for(int q=0;q<3;q++)w->vertex[side][q]=3*a+q;
                w->vertex[side][3]=3*b+(k?2:0);
                w->coefficient[side][0]=(k?0:1)-along;w->coefficient[side][1]=along;
                w->coefficient[side][2]=k?1:0;w->coefficient[side][3]=-1;
            }
        }
    }
    uint8_t core[]={0,1,1},known[]={1,0,0};int failures=!apa_initialize(&field,core,known);
    for(int c=0;c<3;c++)for(int v=0;v<3;v++){
        double x=100+10*c+(v==1?4:0),y=200+(v==2?4:0);
        if(hypot(uv[6*c+2*v]-x,uv[6*c+2*v+1]-y)>1e-10)failures++;
    }
    for(int i=0;i<3;i++)if(!seams[i].required)failures++;
    fprintf(stderr,"  native maximum-support admission forest and retained cycles: %s\n",failures?"FAIL":"ok");return failures;
}

static int apr_admission_control(void)
{
    enum {N=9,NV=N*N,NF=2*(N-1)*(N-1)};int failures=0;
    for(int test=0;test<8;test++){
        Arena_T arena=Arena_new();AsmRun run={0};run.arena=arena;
        float xyz[4][3*NV],uv[4][2*NV];double fixed_uv[2*NV];int32_t faces[3*NF];uint8_t boundary[NV];
        for(int y=0;y<N;y++)for(int x=0;x<N;x++){
            int v=y*N+x;boundary[v]=(uint8_t)(!x || !y || x==N-1 || y==N-1);
            for(int c=0;c<4;c++){
                xyz[c][3*v]=(float)(10*c+x);xyz[c][3*v+1]=(float)y;xyz[c][3*v+2]=(float)(c==3?20:0);
                uv[c][2*v]=(float)x;uv[c][2*v+1]=(float)y;
            }
            fixed_uv[2*v]=x+100.;fixed_uv[2*v+1]=y+200.;
            if(x+1<N && y+1<N){int32_t tri[]={v,v+1,v+N+1,v,v+N+1,v+N};memcpy(faces+6*(y*(N-1)+x),tri,sizeof tri);}
        }
        for(int c=0;c<(test==1?4:3);c++){
            AsmChart ch={0};ch.id=c;ch.parent=-1;ch.nv=NV;ch.nf=NF;ch.xyz=xyz[c];ch.faces=faces;ch.boundary=boundary;ch.uv=uv[c];ch.area3d=64;
            if(c==0 || c==3){ch.placed=1;ch.placement_state=ASM_PLACE_ROOT;ch.pose_x=c?112:100;ch.pose_y=200;}
            if(c==0)ch.placed_uv=fixed_uv;
            AsmRun_push_chart(&run,&ch);
        }
        for(int i=0;i<(test==2?3:2);i++){
            AsmRelation r={0};r.a=i==2?0:i;r.b=i==2?2:i+1;r.flags=i?ASM_REL_DROPPED|ASM_REL_WEAK:0;
            if(i==2)r.flags|=ASM_REL_PARITY;
            r.corr_first=(int32_t)run.n_corr;r.corr_count=N;r.n_corr=N;
            for(int y=0;y<N;y++){AsmCorr corr={0};corr.va=N*y+N-1;corr.vb=N*y;AsmRun_push_corr(&run,&corr);}
            AsmRun_push_rel(&run,&r);
        }
        if(AsmContinuity_prepare(&run)){Arena_dispose(&arena);failures++;continue;}
        if(test==3){run.corr[N-1].gap_a[0]+=1000;run.corr[N-1].gap_b[0]+=1000;}
        /* Source witnesses are already frozen. A missing, nonfinite or
         * folded intrinsic seed must be recovered from the original mesh,
         * without modifying its stored source observations or fixed anchor. */
        if(test==4){run.charts[1].uv=NULL;run.charts[1].flags|=ASM_CHART_FLAT_FAILED;}
        if(test==5)uv[1][0]=NAN;
        if(test==6)uv[1][2*(N+1)]=-3;
        if(test==1)apr_test_scale(&run,16);
        AsmChart before[4];memcpy(before,run.charts,run.n_charts*sizeof(AsmChart));
        AsmCorr samples[3*N];memcpy(samples,run.corr,run.n_corr*sizeof(AsmCorr));
        float original_xyz[4][3*NV],original_uv[4][2*NV];int32_t original_faces[3*NF];
        memcpy(original_xyz,xyz,sizeof xyz);memcpy(original_uv,uv,sizeof uv);memcpy(original_faces,faces,sizeof faces);
        AsmRepairOpts opts;AsmRepair_defaults(&opts);opts.iterations=2;opts.clearance_iterations=1;
        AsmRepairStats stats={0};int rc;
        if(test==7){opts.iterations=opts.clearance_iterations=opts.source_trials=0;rc=AsmRepair_run(&run,&opts,NULL,&stats);}
        else {AsmField *incumbent=AsmField_new(arena,&run);rc=incumbent?apr_admit(&run,&opts,NULL,&stats,incumbent):-1;}
        int pass=test==0 || test>=4;
        int anchor_changed=test==7?memcmp(run.charts[0].placed_uv,fixed_uv,sizeof fixed_uv):memcmp(run.charts,before,sizeof(AsmChart));
        if(rc || stats.admitted_charts!=(size_t)(pass?2:0) || stats.admissions!=(size_t)pass ||
           anchor_changed || memcmp(samples,run.corr,run.n_corr*sizeof(AsmCorr)) ||
           memcmp(original_xyz,xyz,sizeof xyz) || memcmp(original_uv,uv,sizeof uv) || memcmp(original_faces,faces,sizeof faces))failures++;
        for(int c=1;c<=2;c++){
            if(!pass && memcmp(run.charts+c,before+c,sizeof(AsmChart)))failures++;
            if(pass && !AsmChart_registered(run.charts+c))failures++;
        }
        if(pass && (!(run.rels[1].flags&ASM_REL_DROPPED) || !(run.rels[1].continuity&ASM_CONT_REQUIRED)))failures++;
        if(pass){
            AsmField *field=AsmField_new(arena,&run);
            if(!field || field->charts!=3 || !apa_chart_metric(field,1) || !apa_chart_metric(field,2))failures++;
        }
        if(test==1 && memcmp(run.charts+3,before+3,sizeof(AsmChart)))failures++;
        fprintf(stderr,"  native whole-component admission control %d: accepted %zu, charts %zu, cumulative failures %d\n",test,stats.admissions,stats.admitted_charts,failures);
        Arena_dispose(&arena);
    }
    return failures;
}

/* THE SOURCE-REPAIR BATCH CAP (2026-09-19).  A batch's cost is its MOVING CHARTS, the union of
 * its members' causal cores, and nothing bounded it: the 21x5x5's rubble components batch
 * 70-150 charts while the 10x10x10, whose primary holds 63% of its sheet, produced a first
 * batch of 1,975 and sat in one solve for over an hour.  Three properties are pinned here:
 * the union is counted once (overlapping cores are not charged twice), a member that would
 * pass the cap is refused so it forms the next batch, and the FIRST member is always taken
 * however large its own core is, so no seam becomes unreachable. */
static int apa_batch_cap_control(void)
{
    enum { NC = 64 };
    int fails = 0;
    size_t vertex_offset[NC], absent[NC], moving = 0, i = 0, cap = 0;
    uint8_t core[NC], core_k[NC];
    for (i = 0; i < NC; i++) { vertex_offset[i] = i; absent[i] = SIZE_MAX; }

    /* (a) the union is counted once: the same core twice adds nothing the second time */
    memset(core, 0, sizeof core); memset(core_k, 0, sizeof core_k);
    for (i = 0; i < 10; i++) core_k[i] = 1;
    if (!apa_batch_take(core, core_k, NC, vertex_offset, &moving, 0) || moving != 10)
        { fprintf(stderr, "  batch cap: first member should move 10 charts, moved %zu\n", moving); fails++; }
    if (!apa_batch_take(core, core_k, NC, vertex_offset, &moving, 1) || moving != 10)
        { fprintf(stderr, "  batch cap: an overlapping core must add nothing, moving %zu\n", moving); fails++; }

    /* (b) a chart absent from the field is not a moving chart */
    memset(core, 0, sizeof core); moving = 0;
    if (!apa_batch_take(core, core_k, NC, absent, &moving, 0) || moving != 0)
        { fprintf(stderr, "  batch cap: absent charts must not count, moving %zu\n", moving); fails++; }

    /* (c) the cap refuses the member that would pass it; the FIRST is always taken */
    if (ASM_SOURCE_TRIALS_BATCH_MAX_CHARTS > 0 && ASM_SOURCE_TRIALS_BATCH_MAX_CHARTS < NC) {
        cap = (size_t)ASM_SOURCE_TRIALS_BATCH_MAX_CHARTS;
        memset(core, 0, sizeof core); memset(core_k, 0, sizeof core_k);
        for (i = 0; i < cap; i++) core_k[i] = 1;
        moving = 0;
        if (!apa_batch_take(core, core_k, NC, vertex_offset, &moving, 0) || moving != cap)
            { fprintf(stderr, "  batch cap: a core of exactly the cap must be taken\n"); fails++; }
        memset(core_k, 0, sizeof core_k); core_k[cap] = 1;
        if (apa_batch_take(core, core_k, NC, vertex_offset, &moving, 1) || moving != cap)
            { fprintf(stderr, "  batch cap: one chart over the cap must be refused\n"); fails++; }
    }
    memset(core, 0, sizeof core); memset(core_k, 0, sizeof core_k);
    for (i = 0; i < NC; i++) core_k[i] = 1;
    moving = 0;
    if (!apa_batch_take(core, core_k, NC, vertex_offset, &moving, 0) || moving != NC)
        { fprintf(stderr, "  batch cap: a lone oversized member must still run\n"); fails++; }

    fprintf(stderr, "  native source-repair batch cap (%d moving charts): %s (%d failures)\n",
            ASM_SOURCE_TRIALS_BATCH_MAX_CHARTS, fails ? "FAIL" : "ok", fails);
    return fails;
}

/* DETERMINISM (2026-09-17).  The batched source-repair pass stops retrying a trial
 * that left the field byte-identical, on the ground that the same field with the
 * same required seams and the same options decides the same way.  That ground is a
 * property of apr_solve, and this control is it: two solves from one saved state
 * produce the same status and the same coordinates byte for byte.  (The 21x5x5 spent
 * two full trials on component 22 -- 73 moving charts, the same refused
 * initialization both times -- for exactly this reason.) */
static int apr_solve_determinism_control(void)
{
    enum {N=9,NV=N*N,NF=2*(N-1)*(N-1)};
    int failures=0;
    Arena_T arena=Arena_new();AsmRun run={0};run.arena=arena;
    float xyz[3][3*NV],uv[3][2*NV];int32_t faces[3*NF];uint8_t boundary[NV];
    for(int y=0;y<N;y++)for(int x=0;x<N;x++){
        int v=y*N+x;boundary[v]=(uint8_t)(!x || !y || x==N-1 || y==N-1);
        for(int c=0;c<3;c++){
            xyz[c][3*v]=(float)(10*c+x);xyz[c][3*v+1]=(float)y;xyz[c][3*v+2]=0.0f;
            /* the middle chart is sheared, so the solve has real work to do */
            uv[c][2*v]=(float)(x+(c==1?0.25*y:0));uv[c][2*v+1]=(float)y;
        }
        if(x+1<N && y+1<N){int32_t tri[]={v,v+1,v+N+1,v,v+N+1,v+N};memcpy(faces+6*(y*(N-1)+x),tri,sizeof tri);}
    }
    for(int c=0;c<3;c++){
        AsmChart ch={0};ch.id=c;ch.parent=-1;ch.nv=NV;ch.nf=NF;ch.xyz=xyz[c];ch.faces=faces;ch.boundary=boundary;ch.uv=uv[c];ch.area3d=64;
        ch.placed=1;ch.placement_state=ASM_PLACE_ROOT;ch.pose_x=12.0*c;ch.pose_y=0;
        AsmRun_push_chart(&run,&ch);
    }
    for(int i=0;i<2;i++){
        AsmRelation r={0};r.a=i;r.b=i+1;
        r.corr_first=(int32_t)run.n_corr;r.corr_count=N;r.n_corr=N;
        for(int y=0;y<N;y++){AsmCorr corr={0};corr.va=N*y+N-1;corr.vb=N*y;AsmRun_push_corr(&run,&corr);}
        AsmRun_push_rel(&run,&r);
    }
    if(AsmContinuity_prepare(&run)){Arena_dispose(&arena);return 1;}
    AsmField *f=AsmField_new(arena,&run);
    if(!f){Arena_dispose(&arena);return 1;}
    size_t bytes=2*f->nv*sizeof(double);
    double *saved=ARENA_ALLOC(arena,bytes),*first=ARENA_ALLOC(arena,bytes);
    memcpy(saved,f->uv,bytes);
    AsmRepairOpts o;AsmRepair_defaults(&o);o.admission=o.regional=0;o.source_trials=0;o.iterations=3;o.clearance_iterations=1;
    AsmRepairStats s1={0},s2={0};s1.first_metric_chart=s2.first_metric_chart=-1;
    int rc1=apr_solve(f,&o,&s1);
    memcpy(first,f->uv,bytes);
    memcpy(f->uv,saved,bytes);
    int rc2=apr_solve(f,&o,&s2);
    if(rc1!=rc2 || memcmp(first,f->uv,bytes)){
        fprintf(stderr,"  native repair solve determinism: FAIL (status %d then %d, coordinates %s)\n",
                rc1,rc2,memcmp(first,f->uv,bytes)?"differ":"match");
        failures++;
    }
    fprintf(stderr,"  native repair solve determinism from one saved field: %s (%d failures)\n",failures?"FAIL":"ok",failures);
    Arena_dispose(&arena);
    return failures;
}

/* LOCALITY (2026-09-17).  A local trial field must decide exactly what the whole field decides,
 * and the material it leaves out must still be seen.  Part one runs the same candidate twice,
 * once with its context and once over everything, with a registered chart 10,000 vox away: the
 * far chart is outside the context, and the two paths must accept the same coordinates byte for
 * byte.  Part two proves the external probe is what sees material the field omits: a registered
 * chart sitting exactly where the candidate lands, deliberately left out of the member set. */
static int apr_admission_locality_control(void)
{
    enum {N=9,NV=N*N,NF=2*(N-1)*(N-1)};
    int failures=0;
    double kept[2*NV];int have_kept=0,accepted_local=0,accepted_whole=0;
    for(int localise=1;localise>=0;localise--){
        Arena_T arena=Arena_new();AsmRun run={0};run.arena=arena;
        float xyz[4][3*NV],uv[4][2*NV];double fixed_uv[2*NV];int32_t faces[3*NF];uint8_t boundary[NV];
        for(int y=0;y<N;y++)for(int x=0;x<N;x++){
            int v=y*N+x;boundary[v]=(uint8_t)(!x || !y || x==N-1 || y==N-1);
            for(int c=0;c<4;c++){
                xyz[c][3*v]=(float)(10*c+x);xyz[c][3*v+1]=(float)y;xyz[c][3*v+2]=(float)(c==3?20:0);
                uv[c][2*v]=(float)x;uv[c][2*v+1]=(float)y;
            }
            fixed_uv[2*v]=x+100.;fixed_uv[2*v+1]=y+200.;
            if(x+1<N && y+1<N){int32_t tri[]={v,v+1,v+N+1,v,v+N+1,v+N};memcpy(faces+6*(y*(N-1)+x),tri,sizeof tri);}
        }
        for(int c=0;c<4;c++){
            AsmChart ch={0};ch.id=c;ch.parent=-1;ch.nv=NV;ch.nf=NF;ch.xyz=xyz[c];ch.faces=faces;ch.boundary=boundary;ch.uv=uv[c];ch.area3d=64;
            if(c==0){ch.placed=1;ch.placement_state=ASM_PLACE_ROOT;ch.pose_x=100;ch.pose_y=200;ch.placed_uv=fixed_uv;}
            if(c==3){ch.placed=1;ch.placement_state=ASM_PLACE_ROOT;ch.pose_x=10000;ch.pose_y=200;}   /* far away */
            AsmRun_push_chart(&run,&ch);
        }
        for(int i=0;i<2;i++){
            AsmRelation r={0};r.a=i;r.b=i+1;r.flags=i?ASM_REL_DROPPED|ASM_REL_WEAK:0;
            r.corr_first=(int32_t)run.n_corr;r.corr_count=N;r.n_corr=N;
            for(int y=0;y<N;y++){AsmCorr corr={0};corr.va=N*y+N-1;corr.vb=N*y;AsmRun_push_corr(&run,&corr);}
            AsmRun_push_rel(&run,&r);
        }
        if(AsmContinuity_prepare(&run)){Arena_dispose(&arena);failures++;continue;}
        apr_test_scale(&run,8);
        ApaBox *boxes=ARENA_ALLOC(arena,run.n_charts*sizeof(ApaBox));apa_chart_boxes(&run,boxes);
        AsmField *incumbent=AsmField_new(arena,&run);
        uint8_t core[4]={0,1,1,0};
        AsmRepairOpts opts;AsmRepair_defaults(&opts);opts.iterations=2;opts.clearance_iterations=1;
        AsmRepairStats stats={0};size_t anchors=0;const char *reason="";
        ApaOverlay overlay;memset(&overlay,0,sizeof overlay);overlay.arena=arena;
        ApaReport rep;memset(&rep,0,sizeof rep);
        int acc=incumbent?apa_attempt(&run,&opts,core,&anchors,&stats,&reason,NULL,localise?boxes:NULL,incumbent,localise?&overlay:NULL,&rep,0):0;
        if(localise){
            accepted_local=acc;
            if(!rep.local_field || rep.context_charts!=3)failures++;   /* charts 0, 1, 2; never the far one */
            if(acc && run.charts[1].placed_uv){memcpy(kept,run.charts[1].placed_uv,sizeof kept);have_kept=1;}
            /* part two: the probe must see a chart the member set leaves out.  Put one exactly
             * where the candidate landed and omit it from the members. */
            if(acc){
                AsmChart *obstruct=run.charts+3;
                obstruct->placed_uv=ARENA_ALLOC(arena,2*NV*sizeof(double));
                memcpy(obstruct->placed_uv,run.charts[1].placed_uv,2*NV*sizeof(double));
                uint8_t member[4]={1,1,1,0},all[4]={1,1,1,1};
                AsmField *local=AsmField_new_local(arena,&run,member,NULL);
                AsmField *wide=AsmField_new_local(arena,&run,all,NULL);
                AsmField *fresh=AsmField_new(arena,&run);
                size_t pairs=0;double area=0.0;int complete=0;
                if(!local || !wide || !fresh || AsmField_build_contacts(local) ||
                   apa_external_pairs(local,core,member,fresh,NULL,arena,&pairs,&area,&complete,0) ||
                   !complete || !pairs)failures++;
                size_t inside=0;double area2=0.0;int complete2=0;
                if(!wide || AsmField_build_contacts(wide) || apa_external_pairs(wide,core,all,fresh,NULL,arena,&inside,&area2,&complete2,0) ||
                   !complete2 || inside)failures++;   /* in the member set, the local field owns it */
                fprintf(stderr,"  admission external probe sees omitted material: %zu pairs outside the members, %zu inside\n",pairs,inside);
            }
        } else {
            accepted_whole=acc;
            /* the decision AND the coordinates must be the whole field's, byte for byte */
            if(acc && have_kept && (!run.charts[1].placed_uv || memcmp(kept,run.charts[1].placed_uv,sizeof kept)))failures++;
        }
        Arena_dispose(&arena);
    }
    if(accepted_local!=accepted_whole || !accepted_local || !have_kept)failures++;
    fprintf(stderr,"  native local-field admission locality: %s (%d failures; local %d, whole %d)\n",
            failures?"FAIL":"ok",failures,accepted_local,accepted_whole);
    return failures;
}

static int apr_causal_control(void)
{
    const int edges[][2]={{0,1},{1,2},{2,0},{2,3},{3,4},{4,5},{5,3},{1,6}};
    AsmRun run={0};run.n_charts=8;AsmField field={0};field.run=&run;AsmFieldSeam seams[8]={{0}};
    field.seams=seams;field.n_seams=8;for(int i=0;i<8;i++){seams[i].a=edges[i][0];seams[i].b=edges[i][1];seams[i].complete=seams[i].required=1;}
    uint8_t core[8];int failures=0;
    if(apr_causal_core(&field,0,5,core,0)!=6)failures++;
    for(int i=0;i<8;i++)if(core[i]!=(i<6))failures++;
    /* the hop cap: from {0,5} one hop reaches 1, 2, 3, 4 (all six); from {0,0} one hop reaches 1 and 2 */
    if(apr_causal_core(&field,0,5,core,1)!=6)failures++;
    if(apr_causal_core(&field,0,0,core,1)!=3)failures++;
    if(apr_causal_core(&field,0,0,core,0)!=3 || !core[0] || !core[1] || !core[2])failures++;
    if(apr_causal_core(&field,6,6,core,0)!=1 || !core[6])failures++;
    if(apr_causal_core(&field,6,7,core,0)!=2 || !core[6] || !core[7])failures++;
    fprintf(stderr,"  native causal paths retain complete cyclic blocks: %s\n",failures?"FAIL":"ok");return failures;
}

static int apr_pose_model_control(void)
{
    /* Rigid strain invariance and absolute-copy translation curvature are
     * independent properties of Python's proper-pose quadratic. */
    Arena_T arena=Arena_new();float xyz[]={0,0,0,4,0,0,0,4,0},uv[]={0,0,3.2f,0,0,3.8f};int32_t face[]={0,1,2};
    AsmChart chart={0};chart.id=0;chart.nv=3;chart.nf=1;chart.xyz=xyz;chart.faces=face;chart.uv=uv;
    chart.placed=1;chart.placement_state=ASM_PLACE_ROOT;
    AsmRun run={0};run.arena=arena;run.charts=&chart;run.n_charts=1;AsmField *f=AsmField_new(arena,&run);
    if(!f){Arena_dispose(&arena);return 1;}
    AprSolver s={0};s.field=f;AsmRepair_defaults(&s.opts);
    double gradient[6]={0},weight[1]={2},lo[1],hi[1],within10[1],within25[1],area[1];uint8_t fixed[3]={0},rigid[1]={0};
    s.gradient=gradient;s.face_weight=weight;s.lo=lo;s.hi=hi;s.within10=within10;s.within25=within25;s.area=area;
    s.fixed=fixed;s.rigid_fixed=rigid;AprContactCopy contact={{0,1,2},3,2};s.local_contacts=&contact;s.local_count=1;
    apr_checkpoint(&s);double normal[2][9]={{0}};int failures=0;
    for(int test=0;test<2;test++){
        ApcSpace p;weight[0]=test?2e6:2;
        if(!apc_space(&p,&s,arena,APC_POSE,NULL) || !apc_normal(&p) || p.n!=3){failures++;continue;}
        for(int k=0;k<p.nt;k++){normal[test][3*p.rows[k]+p.cols[k]]+=p.values[k];if(p.rows[k]!=p.cols[k])normal[test][3*p.cols[k]+p.rows[k]]+=p.values[k];}
    }
    double translation=12+2*s.opts.proximal/(s.opts.proximal_length*s.opts.proximal_length)*f->metric[0].area;
    if(fabs(normal[0][0]-translation)>1e-10 || fabs(normal[0][4]-translation)>1e-10)failures++;
    for(int k=0;k<9;k++)if(normal[0][k]!=normal[1][k])failures++;
    Arena_dispose(&arena);fprintf(stderr,"  native proper-pose strain invariance and absolute contacts: %s\n",failures?"FAIL":"ok");return failures;
}

static int apr_source_reconstruction_control(void)
{
    enum {N=9,NV=N*N,NF=2*(N-1)*(N-1)};
    int total_failures=0;
    for(int variant=0;variant<15;variant++){
    int internal=variant>=5,reject=variant>=10,scenario=reject?0:variant%5;
    int mirror=scenario>=3;double angle=(scenario==2 || scenario==4)?.43:0;
    double ct=cos(angle),sn=sin(angle),shift[2]={scenario>=2?128:0,scenario>=2?-256:0};
    Arena_T arena=Arena_new();AsmRun run={0};run.arena=arena;run.cube_size=internal?32:10;
    MeshPileEntry cubes[2]={{0}};cubes[0].has_id=cubes[1].has_id=1;cubes[1].oz=10;run.pile=cubes;run.n_cubes=2;
    float xyz[3][3*NV],uv[3][2*NV];int32_t faces[3*NF];uint8_t boundary[NV];
    for(int y=0;y<N;y++)for(int x=0;x<N;x++){
        int v=y*N+x;boundary[v]=(uint8_t)(!x || !y || x==N-1 || y==N-1);
        for(int c=0;c<3;c++){
            double offset=variant==10 || variant==11?0:variant==12?1:variant==14?13:10;
            double height=variant==11?.25:variant==13?1.5:0;
            xyz[c][3*v]=(float)(x+(c==1?offset:c==2?9:0));xyz[c][3*v+1]=(float)(y*(internal?2:1));xyz[c][3*v+2]=(float)(c==1?height:c==2?.1:0);
            double u=x*(c==0 && scenario?1.0625:1),v0=y*(internal?2:1);
            uv[c][2*v]=(float)((c==1 && mirror?-1:1)*(ct*u-sn*v0+shift[0]));
            uv[c][2*v+1]=(float)(sn*u+ct*v0+shift[1]);
        }
        if(x+1<N && y+1<N){int32_t tri[]={v,v+1,v+N+1,v,v+N+1,v+N};memcpy(faces+6*(y*(N-1)+x),tri,sizeof tri);}
    }
    for(int c=0;c<3;c++){
        AsmChart ch={0};ch.id=c;ch.cube=internal?0:c?1:0;ch.parent=-1;ch.nv=NV;ch.nf=NF;ch.xyz=xyz[c];ch.faces=faces;ch.boundary=boundary;ch.uv=uv[c];
        for(int d=0;d<3;d++){ch.bbox_lo[d]=xyz[c][d];ch.bbox_hi[d]=xyz[c][3*(NV-1)+d];}
        AsmRun_push_chart(&run,&ch);
    }
    float original[3][3*NV],original_uv[3][2*NV];memcpy(original,xyz,sizeof xyz);memcpy(original_uv,uv,sizeof uv);
    int failures=AsmContinuity_reconstruct(&run)!=0,found=0;
    for(size_t i=0;i<run.n_rels;i++)if(run.rels[i].a==0 && run.rels[i].b==1){
        AsmRelation *r=run.rels+i;found++;if(r->corr_count<6 || !(r->continuity&ASM_CONT_SOURCE) || !(r->continuity&ASM_CONT_RECONSTRUCTED) || AsmRel_is_join(r))failures++;
        if(reject){failures++;continue;}
        for(int k=0;k<r->corr_count;k++){const AsmCorr *c=run.corr+r->corr_first+k;if(c->valid!=3 || !c->trim_a.valid || !c->trim_b.valid || c->run<0)failures++;}
        /* A six-percent intrinsic extension carries the physical two-voxel
         * gap to 2.125 on A and 2 on B. The two directed targets for B are
         * 10.625 and 10.5: their least-squares translation is 10.5625 and
         * their RMS is .0625. Averaged midpoints conceal that residual. */
        double expected=scenario?10.5625:10,expected_rms=scenario?.0625:0;
        /* Compare the transform at the common origin: float UV rounding
         * can produce a tiny rotation whose translation depends on gauge. */
        double fx=cos(r->theta)*shift[0]-sin(r->theta)*shift[1]+r->tx-shift[0];
        double fy=sin(r->theta)*shift[0]+cos(r->theta)*shift[1]+r->ty-shift[1];
        if(fabs(r->theta)>2e-5 || fabs(fx-ct*expected)>5e-5 || fabs(fy-sn*expected)>5e-5 ||
           fabs(r->rms-expected_rms)>5e-5 || !!(r->flags&ASM_REL_PARITY)!=mirror){
            fprintf(stderr,"  native source fit variant %d: theta %.9g, translation %.9g %.9g, RMS %.9g; expected %.9g %.9g, RMS %.9g\n",
                variant,r->theta,r->tx,r->ty,r->rms,ct*expected,sn*expected,expected_rms);failures++;
        }
    }
    size_t count=run.n_rels,ncorr=run.n_corr;
    /* Negative internal cases: coincident copies, parallel layers,
     * overlapping coplanar patches, excessive normal gap, excessive distance. */
    if(AsmContinuity_reconstruct(&run) || run.n_rels!=count || run.n_corr!=ncorr || found!=(reject?0:1) ||
       memcmp(original,xyz,sizeof xyz) || memcmp(original_uv,uv,sizeof uv))failures++;
    fprintf(stderr,"  native per-chart source reconstruction variant %d with closer third chart: %s (%d failures)\n",variant,failures?"FAIL":"ok",failures);
    Arena_dispose(&arena);total_failures+=failures;
    }
    return total_failures;
}
