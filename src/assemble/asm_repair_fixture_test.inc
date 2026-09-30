static int ags_fixture(Arena_T arena,AsmRun *run)
{
    enum {N=9,NV=N*N,NF=2*(N-1)*(N-1),NC=5};
    memset(run,0,sizeof *run);run->arena=arena;run->n_charts=NC;
    run->charts=ARENA_CALLOC(arena,NC,sizeof *run->charts);
    int32_t *faces=ARENA_ALLOC(arena,3*NF*sizeof *faces);uint8_t *boundary=ARENA_ALLOC(arena,NV);
    for(int y=0;y<N;y++)for(int x=0;x<N;x++){
        int v=y*N+x;boundary[v]=(uint8_t)(!x || !y || x==N-1 || y==N-1);
        if(x+1<N && y+1<N){int32_t tri[]={v,v+1,v+N+1,v,v+N+1,v+N};memcpy(faces+6*(y*(N-1)+x),tri,sizeof tri);}
    }
    for(int c=0;c<NC;c++){
        AsmChart *ch=run->charts+c;ch->id=c;ch->parent=-1;ch->nv=NV;ch->nf=NF;ch->area3d=64;
        ch->xyz=ARENA_ALLOC(arena,3*NV*sizeof *ch->xyz);ch->uv=ARENA_ALLOC(arena,2*NV*sizeof *ch->uv);
        ch->vid=ARENA_ALLOC(arena,NV*sizeof *ch->vid);ch->faces=faces;ch->boundary=boundary;
        for(int y=0;y<N;y++)for(int x=0;x<N;x++){
            int v=y*N+x;ch->vid[v]=v;ch->xyz[3*v]=(float)(10*c+x);ch->xyz[3*v+1]=(float)y;ch->xyz[3*v+2]=0;
            ch->uv[2*v]=(float)x;ch->uv[2*v+1]=(float)y;
        }
        if(c!=1 && c!=4){ch->placed=1;ch->placement_state=ASM_PLACE_ROOT;}
        ch->pose_x=c==0?100:c==2?117.8:c==3?127.8:1000;ch->pose_y=200;
    }
    for(int i=0;i<3;i++){
        AsmRelation r={0};r.a=i==0?0:i==1?2:1;r.b=i==0?1:i==1?3:4;
        r.corr_first=(int32_t)run->n_corr;r.corr_count=N;r.n_corr=N;
        for(int y=0;y<N;y++){AsmCorr w={0};w.va=y*N+N-1;w.vb=y*N;AsmRun_push_corr(run,&w);}
        AsmRun_push_rel(run,&r);
    }
    return AsmContinuity_prepare(run);
}

/* A 3x3 grid of 9x9 charts joined by required source seams, the centre
 * chart unplaced. Shared by the admission, rank and review controls. */
static int ags_grid_fixture(Arena_T arena,AsmRun *run)
{
    enum {N=9,NV=N*N,NF=2*(N-1)*(N-1),NC=9};
    memset(run,0,sizeof *run);run->arena=arena;run->n_charts=NC;
    run->charts=ARENA_CALLOC(arena,NC,sizeof *run->charts);
    const double shift[NC]={-4,0,4,-8,0,8,-4,0,4};
    for(int c=0;c<NC;c++){
        AsmChart *ch=run->charts+c;ch->id=c;ch->parent=-1;ch->nv=NV;ch->nf=NF;ch->area3d=400;
        ch->xyz=ARENA_ALLOC(arena,3*NV*sizeof(float));ch->uv=ARENA_ALLOC(arena,2*NV*sizeof(float));
        ch->vid=ARENA_ALLOC(arena,NV*sizeof(int32_t));ch->faces=ARENA_ALLOC(arena,3*NF*sizeof(int32_t));
        ch->boundary=ARENA_ALLOC(arena,NV);ch->placed=c!=4;ch->placement_state=c==4?ASM_PLACE_NONE:ASM_PLACE_ROOT;
        ch->pose_x=20*(c%3)+shift[c];ch->pose_y=20*(c/3);
        for(int y=0;y<N;y++)for(int x=0;x<N;x++){
            int v=y*N+x;ch->vid[v]=v;ch->boundary[v]=(uint8_t)(!x || !y || x==N-1 || y==N-1);
            ch->xyz[3*v]=(float)(20*(c%3)+2.5*x);ch->xyz[3*v+1]=(float)(20*(c/3)+2.5*y);ch->xyz[3*v+2]=0;
            ch->uv[2*v]=(float)(2.5*x);ch->uv[2*v+1]=(float)(2.5*y);
            if(x+1<N && y+1<N){int32_t t[]={v,v+1,v+N+1,v,v+N+1,v+N};memcpy(ch->faces+6*(y*(N-1)+x),t,sizeof t);}
        }
    }
    for(int c=0;c<NC;c++)for(int d=0;d<2;d++)if(d?c<6:c%3<2){
        AsmRelation r={0};r.a=c;r.b=c+(d?3:1);r.corr_first=(int32_t)run->n_corr;r.corr_count=N;r.n_corr=N;
        r.flags=ASM_REL_WEAK;r.continuity=ASM_CONT_SOURCE|ASM_CONT_REQUIRED;
        for(int k=0;k<N;k++){
            AsmCorr w={0};w.va=d?(N-1)*N+k:k*N+N-1;w.vb=d?k:k*N;w.valid=3;
            w.trim_a=apr_test_trim(run->charts+r.a,w.va);w.trim_b=apr_test_trim(run->charts+r.b,w.vb);
            AsmRun_push_corr(run,&w);
        }
        AsmRun_push_rel(run,&r);
    }
    run->continuity_ready=1;
    return 0;
}

/* A private repair solve whose shared deadline has already passed returns
 * without touching the field, and the admission solver's deadline never
 * extends past the request's. */
static int ags_solver_owner_control(void)
{
    Arena_T arena=Arena_new();AsmRun run={0};int failures=0;
    if(ags_fixture(arena,&run)){Arena_dispose(&arena);return 1;}
    apr_test_scale(&run,32);
    AsmField_threads(1);
    {
        AsmField *field=AsmField_new(arena,&run);AsmRepairOpts request,local;AsmRepair_defaults(&request);
        request.iterations=2;request.deadline=ves_clock_sec()-1;request.admission_solver_deadline=request.deadline;local=request;
        AsmRepairStats stats={0};ApaReport report={0};
        if(!field)failures++;else{
            size_t bytes=2*field->nv*sizeof(double);double *saved=ARENA_ALLOC(arena,bytes);memcpy(saved,field->uv,bytes);
            if(apa_solve_owned(field,&local,&request,&stats,&report) || !stats.deadline_reached ||
               request.deadline>request.admission_solver_deadline || memcmp(saved,field->uv,bytes))failures++;
        }
    }
    AsmField_threads(24);Arena_dispose(&arena);
    fprintf(stderr,"  private repair solver deadline ownership: %s (%d failures)\n",failures?"FAIL":"ok",failures);return failures;
}
