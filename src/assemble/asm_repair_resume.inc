/* Private candidates travel inside the checksummed native checkpoint. They
 * are starting guesses, never placement evidence or reusable certificates. */
typedef struct ArrBytes {const uint8_t *at;size_t left;} ArrBytes;
typedef struct ArrSaved {
    size_t target,rounds;
    double forecast;
    uint8_t *moving,*dependencies;
    ApaContinuation state;
} ArrSaved;

static int arr_take(ArrBytes *in,void *out,size_t bytes)
{
    if(bytes>in->left)return 0;
    if(out && bytes)memcpy(out,in->at,bytes);
    in->at+=bytes;in->left-=bytes;return 1;
}

static void arr_chart_key(const AsmRun *run,size_t chart,char key[65],int legacy)
{
    const AsmChart *c=run->charts+chart;VesSha256 hash;VesSha256_init(&hash);
    uint64_t shape[]={chart,c->nv,c->nf,c->flags,AsmChart_registered(c),
        c->vid!=NULL,c->nrm!=NULL,c->uv!=NULL,c->boundary!=NULL,c->repaired!=NULL,c->placed_uv!=NULL};
    VesSha256_update(&hash,shape,legacy?sizeof shape:sizeof shape-sizeof shape[0]);
    VesSha256_update(&hash,&c->area3d,sizeof c->area3d);
#define ARR_HASH(name,count) do{if(c->name)VesSha256_update(&hash,c->name,(count)*sizeof *c->name);}while(0)
    ARR_HASH(xyz,3*c->nv);ARR_HASH(faces,3*c->nf);ARR_HASH(vid,c->nv);ARR_HASH(nrm,3*c->nv);
    ARR_HASH(uv,2*c->nv);ARR_HASH(boundary,c->nv);ARR_HASH(repaired,c->nv);
    if(AsmChart_registered(c)){
        if(!legacy){
            /* Field commit may encode an unchanged pose as float64 vertex
             * UVs. Hash the actual coordinates, not their representation. */
            for(size_t v=0;v<c->nv;v++){double p[2];AsmChart_point(c,v,p);VesSha256_update(&hash,p,sizeof p);}
        }else if(c->placed_uv){ARR_HASH(placed_uv,2*c->nv);}
        else{double pose[]={c->pose_x,c->pose_y,c->pose_theta};VesSha256_update(&hash,pose,sizeof pose);}
    }
#undef ARR_HASH
    /* The ASR format preserves these exact source rows, including rejected
     * relations. A changed witness or newly required boundary invalidates
     * the candidate just as changed source faces or neighbor UVs do. */
    for(size_t i=0;i<run->n_rels;i++){
        const AsmRelation *r=run->rels+i;if(r->a!=(int32_t)chart && r->b!=(int32_t)chart)continue;
        VesSha256_update(&hash,&i,sizeof i);VesSha256_update(&hash,r,sizeof *r);
        if(r->corr_count)VesSha256_update(&hash,run->corr+r->corr_first,(size_t)r->corr_count*sizeof *run->corr);
    }
    VesSha256_hex(&hash,key);
}

static const char *arr_cached_key(const AsmRun *run,size_t chart,char *keys,int legacy)
{
    char *key=keys+65*chart;if(!*key)arr_chart_key(run,chart,key,legacy);return key;
}

static int arr_resume_read(const AsmRun *run,const AsmRepairOpts *policy,Arena_T arena,
                           ApaContinuation *states,uint8_t **pending,uint8_t **dependencies,
                           size_t *rounds,double *forecast,ArrStats *stats)
{
    if(!run->repair_resume_size)return 0;
    Arena_T scratch=Arena_new();ArrBytes in={run->repair_resume,run->repair_resume_size};
    size_t nc=run->n_charts,count=0;int okay=0,legacy=0,current_progress=0;char magic[8];uint64_t header[2];
    ArrSaved *entries=NULL;char *keys=ARENA_CALLOC(scratch,nc,65);uint8_t *seen=ARENA_CALLOC(scratch,nc,1);
    if(!arr_take(&in,magic,8))goto done;legacy=!memcmp(magic,"ARRUV001",8);
    current_progress=!memcmp(magic,"ARRUV004",8);
    if((!legacy && !current_progress && memcmp(magic,"ARRUV002",8) && memcmp(magic,"ARRUV003",8)) || !arr_take(&in,header,sizeof header) ||
       header[0]!=nc || header[1]>nc)goto done;
    count=(size_t)header[1];entries=ARENA_CALLOC(scratch,count,sizeof *entries);
    for(size_t e=0;e<count;e++){
        ArrSaved *s=entries+e;uint64_t meta[5];double prediction=0;
        if(!arr_take(&in,meta,sizeof meta) || !arr_take(&in,&prediction,sizeof prediction) ||
           meta[0]>=nc || seen[meta[0]] || !meta[1] || meta[1]>UINT32_MAX ||
           !meta[2] || meta[2]>nc || meta[3]>=nc || !meta[4] || meta[4]>nc ||
           !isfinite(prediction) || prediction<0)goto done;
        s->target=(size_t)meta[0];s->rounds=(size_t)meta[1];
        /* Older forecasts used a different contact merit and checkpoint
         * policy. Bound their pessimism by one current target allowance, while
         * retaining useful short estimates until the next fresh measurement. */
        s->forecast=current_progress?prediction:fmin(prediction,policy->patch_seconds);seen[s->target]=1;
        s->moving=ARENA_CALLOC(scratch,nc,1);s->dependencies=ARENA_CALLOC(scratch,nc,1);
        int stale=AsmChart_registered(run->charts+s->target) || !AsmChart_in_layout(run->charts+s->target);
        size_t changed_chart=SIZE_MAX;
        for(uint64_t k=0;k<meta[2];k++){
            uint64_t chart=0;char stored[64];
            if(!arr_take(&in,&chart,sizeof chart) || chart>=nc || s->dependencies[chart] ||
               !arr_take(&in,stored,sizeof stored))goto done;
            s->dependencies[chart]=1;
            if(memcmp(stored,arr_cached_key(run,(size_t)chart,keys,legacy),64)){
                stale=1;if(changed_chart==SIZE_MAX)changed_chart=(size_t)chart;
            }
        }
        if(!s->dependencies[s->target])goto done;
        for(uint64_t k=0;k<meta[3];k++){
            uint64_t chart=0;
            if(!arr_take(&in,&chart,sizeof chart) || chart>=nc || chart==s->target || s->moving[chart] ||
               !s->dependencies[chart])goto done;
            s->moving[chart]=1;if(!AsmChart_registered(run->charts+chart))stale=1;
        }
        size_t faces=0,coordinates=0;
        if(!arr_cost(run,s->target,s->moving,policy,&faces,&coordinates))stale=1;
        /* sized like a saved continuation's: its chart tables, then chunks as the UVs need */
        s->state.arena=Arena_new_sized(nc*(sizeof(double *)+sizeof(size_t)+1)+((size_t)1<<20));s->state.count=nc;
        s->state.uv=ARENA_CALLOC(s->state.arena,nc,sizeof *s->state.uv);
        s->state.vertices=ARENA_CALLOC(s->state.arena,nc,sizeof *s->state.vertices);
        s->state.mirror=ARENA_CALLOC(s->state.arena,nc,1);
        uint8_t *uv_seen=ARENA_CALLOC(scratch,nc,1);size_t total=0;
        for(uint64_t k=0;k<meta[4];k++){
            uint64_t shape[3];
            if(!arr_take(&in,shape,sizeof shape) || shape[0]>=nc || shape[1]>SIZE_MAX/(2*sizeof(double)) ||
               shape[2]>1 || uv_seen[shape[0]] || (shape[0]!=s->target && !s->moving[shape[0]]))goto done;
            size_t chart=(size_t)shape[0],nv=(size_t)shape[1],bytes=2*nv*sizeof(double);uv_seen[chart]=1;
            if(bytes>in.left || nv>SIZE_MAX-total)goto done;total+=nv;
            if(nv!=run->charts[chart].nv || total>policy->patch_coordinates/2 ||
               shape[2]!=(uint64_t)!!(run->charts[chart].flags&ASM_CHART_MIRROR))stale=1;
            if(!stale){
                double *uv=ARENA_ALLOC(s->state.arena,bytes);
                if(!arr_take(&in,uv,bytes))goto done;
                for(size_t v=0;v<2*nv;v++)if(!isfinite(uv[v]))goto done;
                s->state.uv[chart]=uv;s->state.vertices[chart]=nv;s->state.mirror[chart]=(uint8_t)shape[2];
            }else if(!arr_take(&in,NULL,bytes))goto done;
        }
        if(!uv_seen[s->target])goto done;
        if(stale){
            if(changed_chart!=SIZE_MAX)fprintf(stderr,"[region resume] target %zu discarded: dependency chart %zu changed\n",s->target,changed_chart);
            apa_continuation_clear(&s->state);stats->resume_stale++;
        }
    }
    if(in.left)goto done;
    /* No decoded entry reaches the scheduler before the complete payload
     * has passed structural validation. Live chart arrays are never touched. */
    for(size_t e=0;e<count;e++)if(entries[e].state.arena){
        ArrSaved *s=entries+e;size_t c=s->target;
        states[c]=s->state;memset(&s->state,0,sizeof s->state);
        pending[c]=ARENA_ALLOC(arena,nc);memcpy(pending[c],s->moving,nc);
        dependencies[c]=ARENA_ALLOC(arena,nc);memcpy(dependencies[c],s->dependencies,nc);
        rounds[c]=s->rounds;forecast[c]=s->forecast;stats->resume_loaded++;
        /* The per-target work allowance is for this run. Previous effort
         * affects the forecast, but does not permanently ban a candidate. */
    }
    okay=1;
done:
    if(entries)for(size_t e=0;e<count;e++)apa_continuation_clear(&entries[e].state);
    if(!okay)stats->resume_invalid++;
    Arena_dispose(&scratch);return okay?0:-1;
}

static int arr_resume_encode(AsmRun *run,const ApaContinuation *states,uint8_t *const *pending,
                             uint8_t *const *dependencies,const size_t *rounds,const double *forecast,ArrStats *stats,int legacy)
{
    Arena_T scratch=Arena_new();size_t nc=run->n_charts,bytes=24,count=0;
    uint8_t *keep=ARENA_CALLOC(scratch,nc,1);char *keys=ARENA_CALLOC(scratch,nc,65);
    uint64_t *meta=ARENA_CALLOC(scratch,nc,5*sizeof *meta);int rc=-1;
    for(size_t c=0;c<nc;c++)if(states[c].arena && pending[c] && dependencies[c] &&
       rounds[c] && !AsmChart_registered(run->charts+c)){
        size_t need=48;uint64_t *m=meta+5*c;m[0]=c;m[1]=rounds[c];
        for(size_t k=0;k<nc;k++){
            m[2]+=dependencies[c][k] || pending[c][k] || k==c;
            m[3]+=pending[c][k]!=0;
            if(states[c].uv[k]){
                size_t nv=states[c].vertices[k];if(nv>(SIZE_MAX-24-need)/(2*sizeof(double)))goto done;
                need+=24+2*nv*sizeof(double);m[4]++;
            }
        }
        if(m[2]>(SIZE_MAX-need)/72)goto done;need+=(size_t)m[2]*72;
        if(m[3]>(SIZE_MAX-need)/8)goto done;need+=(size_t)m[3]*8;
        if(need>SIZE_MAX-bytes)goto done;bytes+=need;keep[c]=1;count++;
    }
    if(!count){run->repair_resume=NULL;run->repair_resume_size=0;rc=0;goto done;}
    uint8_t *blob=ARENA_ALLOC(run->arena,bytes),*at=blob;uint64_t header[2]={nc,count};
#define ARR_PUT(p,n) do{size_t amount=(n);memcpy(at,(p),amount);at+=amount;}while(0)
    ARR_PUT(legacy?"ARRUV001":"ARRUV004",8);ARR_PUT(header,sizeof header);
    for(size_t c=0;c<nc;c++)if(keep[c]){
        double prediction=isfinite(forecast[c])?forecast[c]:ASM_REGION_REFINE_ATTEMPT_SEC;
        ARR_PUT(meta+5*c,5*sizeof(uint64_t));ARR_PUT(&prediction,sizeof prediction);
        for(size_t k=0;k<nc;k++)if(dependencies[c][k] || pending[c][k] || k==c){
            uint64_t id=k;ARR_PUT(&id,sizeof id);ARR_PUT(arr_cached_key(run,k,keys,legacy),64);
        }
        for(size_t k=0;k<nc;k++)if(pending[c][k]){uint64_t id=k;ARR_PUT(&id,sizeof id);}
        for(size_t k=0;k<nc;k++)if(states[c].uv[k]){
            uint64_t shape[3]={k,states[c].vertices[k],states[c].mirror[k]};
            ARR_PUT(shape,sizeof shape);ARR_PUT(states[c].uv[k],2*states[c].vertices[k]*sizeof(double));
        }
    }
    if((size_t)(at-blob)!=bytes)goto done;
    run->repair_resume=blob;run->repair_resume_size=bytes;stats->resume_saved=count;rc=0;
#undef ARR_PUT
done:
    Arena_dispose(&scratch);return rc;
}

static int arr_resume_write(AsmRun *run,const ApaContinuation *states,uint8_t *const *pending,
                            uint8_t *const *dependencies,const size_t *rounds,const double *forecast,ArrStats *stats)
{
    return arr_resume_encode(run,states,pending,dependencies,rounds,forecast,stats,0);
}

/* Keep untrusted guesses across the native stage's validated seam patches
 * and admissions. Capture before those operations change required relations
 * or incumbent UVs; refresh only after their normal certificates succeed.
 * The decoder remains strict for externally changed checkpoint inputs. */
typedef struct ArrHeld {
    Arena_T arena;
    size_t count;
    ApaContinuation *states;
    uint8_t **pending,**dependencies;
    size_t *rounds;
    double *forecast;
    double seconds;
    char *keys;
    ArrStats stats;
} ArrHeld;

static void arr_held_clear(ArrHeld *held)
{
    if(held->arena)for(size_t c=0;c<held->count;c++)apa_continuation_clear(held->states+c);
    if(held->arena)Arena_dispose(&held->arena);memset(held,0,sizeof *held);
}

static void arr_hold(AsmRun *run,const AsmRepairOpts *policy,ArrHeld *held)
{
    memset(held,0,sizeof *held);if(!run->repair_resume_size)return;
    double started=ves_clock_sec();
    Arena_T arena=Arena_new();held->arena=arena;held->count=run->n_charts;size_t nc=held->count;
    held->states=ARENA_CALLOC(arena,nc,sizeof *held->states);
    held->pending=ARENA_CALLOC(arena,nc,sizeof *held->pending);held->dependencies=ARENA_CALLOC(arena,nc,sizeof *held->dependencies);
    held->rounds=ARENA_CALLOC(arena,nc,sizeof *held->rounds);held->forecast=ARENA_CALLOC(arena,nc,sizeof *held->forecast);
    held->keys=ARENA_CALLOC(arena,nc,65);
    if(arr_resume_read(run,policy,arena,held->states,held->pending,held->dependencies,held->rounds,held->forecast,&held->stats)){
        /* Corrupt optional scratch data must not prevent valid geometry from
         * being retained or acquire fresh fingerprints by being rewritten. */
        run->repair_resume=NULL;run->repair_resume_size=0;arr_held_clear(held);return;
    }
    for(size_t c=0;c<nc;c++)if(held->states[c].arena)
        for(size_t k=0;k<nc;k++)if(held->dependencies[c][k])arr_cached_key(run,k,held->keys,0);
    held->seconds=ves_clock_sec()-started;
}

static int arr_refresh(AsmRun *run,const AsmRepairOpts *policy,ArrHeld *held)
{
    if(!held->arena)return 0;
    if(held->count!=run->n_charts)return -1;
    double started=ves_clock_sec();
    size_t nc=held->count,refreshed=0;char *after=ARENA_CALLOC(held->arena,nc,65);
    for(size_t c=0;c<nc;c++)if(held->states[c].arena && !AsmChart_registered(run->charts+c))
        for(size_t k=0;k<nc;k++)if(held->dependencies[c][k] &&
            memcmp(held->keys+65*k,arr_cached_key(run,k,after,0),64)){
            held->forecast[c]=policy->patch_seconds;refreshed++;break;
        }
    int rc=arr_resume_write(run,held->states,held->pending,held->dependencies,held->rounds,held->forecast,&held->stats);
    fprintf(stderr,"[region resume] after validated repairs: %zu guesses retained, %zu context forecasts refreshed, %zu stale discarded, %.3f seconds\n",
        held->stats.resume_saved,refreshed,held->stats.resume_stale,held->seconds+ves_clock_sec()-started);
    return rc;
}
