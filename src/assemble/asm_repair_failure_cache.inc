/* Cache a completed bounded local attempt, not geometric infeasibility.
 * The cache is private to one process/executable. Exact local bytes include
 * all source geometry, current/intrinsic UVs, witnesses and solver settings.
 * Exterior-only failures, timeouts and incomplete walks are never retained. */
typedef struct AfcEntry { char key[65];const char *reason;double seconds; } AfcEntry;
typedef struct ApaFailureCache {
    AfcEntry *entries;size_t capacity,count,probes,hits,full;
    double hash_seconds,previous_work_seconds;
} ApaFailureCache;

static void afc_init(ApaFailureCache *cache,Arena_T arena,size_t capacity)
{
    memset(cache,0,sizeof *cache);cache->capacity=capacity;
    cache->entries=ARENA_CALLOC(arena,capacity,sizeof *cache->entries);
}

static AfcEntry *afc_slot(ApaFailureCache *cache,const char key[65])
{
    if(!cache || !cache->capacity)return NULL;
    size_t code=0;for(size_t i=0;i<64;i++)code=code*33+(unsigned char)key[i];
    for(size_t n=0;n<cache->capacity;n++){
        AfcEntry *entry=cache->entries+(code+n)%cache->capacity;
        if(!entry->reason || !memcmp(entry->key,key,65))return entry;
    }
    return NULL;
}

static const AfcEntry *afc_lookup(ApaFailureCache *cache,const char key[65])
{
    cache->probes++;AfcEntry *entry=afc_slot(cache,key);
    if(!entry || !entry->reason)return NULL;
    cache->hits++;cache->previous_work_seconds+=entry->seconds;return entry;
}

static int afc_reason(const char *reason,int complete,int expired)
{
    return complete && !expired && (!strcmp(reason,"new_material_metric_failed") ||
        !strcmp(reason,"incident_source_seam_failed") || !strcmp(reason,"new_material_contact_failed"));
}

static void afc_insert(ApaFailureCache *cache,const char key[65],const char *reason,double seconds)
{
    AfcEntry *entry=afc_slot(cache,key);if(!entry){cache->full++;return;}
    if(entry->reason)return;
    memcpy(entry->key,key,65);entry->reason=reason;entry->seconds=seconds;cache->count++;
}

static void afc_key(const AsmField *f,const AsmRepairOpts *opts,const uint8_t *moving,
                    const double *reference,int frontier,char key[65])
{
    VesSha256 hash;VesSha256_init(&hash);const char schema[]="local-admission-failure-v1";
    VesSha256_update(&hash,schema,sizeof schema);VesSha256_update(&hash,&frontier,sizeof frontier);
    /* Scalar options retain their exact bits. Pointer identities and absolute
     * deadlines are not solver inputs; their payloads are bound below. Cache
     * eligibility forbids callbacks, joint motion, restoration and exports. */
    AsmRepairOpts settings=*opts;settings.deadline=0;settings.fixed_charts=NULL;
    settings.admission_moving=NULL;settings.admission_trace_dir=NULL;settings.admission_failure_cache=NULL;
    settings.external_guard=NULL;settings.external_context=NULL;settings.checkpoint=NULL;settings.checkpoint_context=NULL;
    settings.layer_place=NULL;settings.layer_conflict=NULL;
    VesSha256_update(&hash,&settings,sizeof settings);
    size_t sizes[]={f->run->n_charts,f->nv,f->nf,f->n_seams,f->n_witness,AsmContacts_budget(f->contacts)};
    VesSha256_update(&hash,sizes,sizeof sizes);VesSha256_update(&hash,moving,f->run->n_charts);
    VesSha256_update(&hash,f->vertex_offset,f->run->n_charts*sizeof(size_t));
    VesSha256_update(&hash,f->face_offset,f->run->n_charts*sizeof(size_t));
    VesSha256_update(&hash,f->xyz,3*f->nv*sizeof(float));
    VesSha256_update(&hash,f->faces,3*f->nf*sizeof(int32_t));
    VesSha256_update(&hash,f->uv,2*f->nv*sizeof(double));
    VesSha256_update(&hash,reference,2*f->nv*sizeof(double));
    for(size_t c=0;c<f->run->n_charts;c++)if(f->vertex_offset[c]!=SIZE_MAX){
        const AsmChart *ch=f->run->charts+c;
        size_t shape[]={c,ch->nv,ch->nf};int present[]={!!ch->uv,!!ch->boundary,!!ch->nrm,!!ch->repaired};
        VesSha256_update(&hash,shape,sizeof shape);VesSha256_update(&hash,present,sizeof present);
        VesSha256_update(&hash,&ch->flags,sizeof ch->flags);
        if(ch->uv)VesSha256_update(&hash,ch->uv,2*ch->nv*sizeof(float));
        if(ch->boundary)VesSha256_update(&hash,ch->boundary,ch->nv);
        if(ch->nrm)VesSha256_update(&hash,ch->nrm,3*ch->nv*sizeof(float));
        if(ch->repaired)VesSha256_update(&hash,ch->repaired,ch->nv);
    }
    for(size_t i=0;i<f->n_seams;i++){
        const AsmFieldSeam *s=f->seams+i;
        size_t span[]={s->relation,s->first,s->count};
        int state[]={s->a,s->b,s->complete,s->parity,s->source_run,s->required,s->original_cut,
            AsmChart_registered(f->run->charts+s->a),AsmChart_registered(f->run->charts+s->b)};
        VesSha256_update(&hash,span,sizeof span);VesSha256_update(&hash,state,sizeof state);
    }
    /* AsmWitness contains no padding: eight int32s, nine doubles, two ints. */
    VesSha256_update(&hash,f->witness,f->n_witness*sizeof(AsmWitness));
    VesSha256_hex(&hash,key);
}
