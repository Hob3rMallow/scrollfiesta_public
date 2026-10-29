/* A native checkpoint, not a replacement clean cache. Pointer values are
 * cleared before serialization; every optional array has an explicit bit. */
typedef struct AsRunIo { FILE *fp; VesSha256 sha; uint64_t remaining; } AsRunIo;
static int asr_put(AsRunIo *io,const void *p,size_t n)
{
    if(n && (!p || fwrite(p,1,n,io->fp)!=n))return -1;
    VesSha256_update(&io->sha,p,n);return 0;
}
static int asr_get(AsRunIo *io,void *p,size_t n)
{
    if(n>io->remaining || (n && fread(p,1,n,io->fp)!=n))return -1;
    io->remaining-=n;VesSha256_update(&io->sha,p,n);return 0;
}
static void asr_abi(uint64_t sizes[6])
{
    sizes[0]=sizeof(AsmChart);sizes[1]=sizeof(AsmRelation);sizes[2]=sizeof(AsmCorr);
    sizes[3]=sizeof(AsmLayerPair);sizes[4]=sizeof(AsmLayerHit);sizes[5]=sizeof(MeshPileEntry);
}
int AsmStore_write_run(const char *path,const AsmRun *run)
{
    if(!path || !run || !run->n_charts)return -1;
    char tmp[2048];if(snprintf(tmp,sizeof tmp,"%s.tmp",path)>=(int)sizeof tmp)return -1;
    ves_ensure_parent_dir(path);AsRunIo io={0};io.fp=fopen(tmp,"wb");if(!io.fp)return -1;
    VesSha256_init(&io.sha);int rc=-1;uint64_t abi[6];asr_abi(abi);
    uint64_t counts[7]={run->n_cubes,run->n_charts,run->n_rels,run->n_corr,run->n_layers,run->n_layer_hits,run->chart_layer_d!=NULL};
    double scalar[2]={run->cube_size,run->layer_d_global};int32_t ready=run->continuity_ready;
#define ASR_PUT(p,n) do{if(asr_put(&io,(p),(n)))goto done;}while(0)
    ASR_PUT(run->repair_resume_size?"ASMRUN02":"ASMRUN01",8);ASR_PUT(abi,sizeof abi);ASR_PUT(counts,sizeof counts);ASR_PUT(scalar,sizeof scalar);ASR_PUT(&ready,sizeof ready);
    ASR_PUT(run->pile,run->n_cubes*sizeof(MeshPileEntry));
    for(size_t c=0;c<run->n_charts;c++){
        const AsmChart *source=run->charts+c;AsmChart row=*source;uint32_t present=0;
#define ASR_PTR(name,bit) do{if(row.name)present|=1u<<(bit);row.name=NULL;}while(0)
        ASR_PTR(vid,0);ASR_PTR(xyz,1);ASR_PTR(nrm,2);ASR_PTR(uv,3);
        ASR_PTR(placed_uv,4);ASR_PTR(faces,5);ASR_PTR(boundary,6);ASR_PTR(repaired,7);
#undef ASR_PTR
        ASR_PUT(&row,sizeof row);ASR_PUT(&present,sizeof present);
#define ASR_ARRAY(name,count) do{if(source->name)ASR_PUT(source->name,(count)*sizeof *source->name);}while(0)
        ASR_ARRAY(vid,row.nv);ASR_ARRAY(xyz,3*row.nv);ASR_ARRAY(nrm,3*row.nv);ASR_ARRAY(uv,2*row.nv);
        ASR_ARRAY(placed_uv,2*row.nv);ASR_ARRAY(faces,3*row.nf);ASR_ARRAY(boundary,row.nv);ASR_ARRAY(repaired,row.nv);
#undef ASR_ARRAY
    }
    ASR_PUT(run->rels,run->n_rels*sizeof(AsmRelation));ASR_PUT(run->corr,run->n_corr*sizeof(AsmCorr));
    ASR_PUT(run->layers,run->n_layers*sizeof(AsmLayerPair));ASR_PUT(run->layer_hits,run->n_layer_hits*sizeof(AsmLayerHit));
    if(run->chart_layer_d)ASR_PUT(run->chart_layer_d,run->n_charts*sizeof(double));
    if(run->repair_resume_size){
        uint64_t bytes=run->repair_resume_size;ASR_PUT(&bytes,sizeof bytes);ASR_PUT(run->repair_resume,run->repair_resume_size);
    }
    {char digest[65];VesSha256_hex(&io.sha,digest);if(fwrite(digest,1,64,io.fp)!=64)goto done;}
    rc=0;
done:
    if(fclose(io.fp))rc=-1;
    if(!rc){
#ifdef _WIN32
        if(!MoveFileExA(tmp,path,MOVEFILE_REPLACE_EXISTING|MOVEFILE_WRITE_THROUGH))rc=-1;
#else
        if(rename(tmp,path))rc=-1;
#endif
    }
    if(rc)remove(tmp);return rc;
#undef ASR_PUT
}
/* Skip n payload bytes. Hashed, and so still verified, unless `seek`, which
 * leaves the digest incomplete and is used only by the graph read. */
static int asr_skip(AsRunIo *io,size_t n,int seek)
{
    if(n>io->remaining)return -1;
    if(seek){
#ifdef _WIN32
        if(n && _fseeki64(io->fp,(__int64)n,SEEK_CUR))return -1;
#else
        if(n && fseeko(io->fp,(off_t)n,SEEK_CUR))return -1;
#endif
        io->remaining-=n;return 0;
    }
    unsigned char buffer[65536];
    while(n){size_t k=n<sizeof buffer?n:sizeof buffer;if(asr_get(io,buffer,k))return -1;n-=k;}
    return 0;
}

/* mode 0: everything; 1: no chart arrays, digest not verified; 2: chart
 * arrays only where keep[c] (NULL keeps all), digest verified. */
static int asr_read(const char *path,AsmRun *run,int mode,const uint8_t *keep,uint8_t **has_uv);
int AsmStore_read_run(const char *path,AsmRun *run){return asr_read(path,run,0,NULL,NULL);}

static int asr_read(const char *path,AsmRun *run,int mode,const uint8_t *keep,uint8_t **has_uv)
{
    if(!path || !run || !run->arena || run->n_charts || mode<0 || mode>2)return -1;
    if(has_uv)*has_uv=NULL;
    uint8_t *uv_present=NULL;
    AsRunIo io={0};io.fp=fopen(path,"rb");if(!io.fp)return -1;
    int rc=-1;AsmRun candidate={0};candidate.arena=run->arena;
#ifdef _WIN32
    if(_fseeki64(io.fp,0,SEEK_END))goto done;__int64 length=_ftelli64(io.fp);
    if(length<196 || _fseeki64(io.fp,0,SEEK_SET))goto done;
#else
    if(fseeko(io.fp,0,SEEK_END))goto done;off_t length=ftello(io.fp);
    if(length<196 || fseeko(io.fp,0,SEEK_SET))goto done;
#endif
    io.remaining=(uint64_t)length-64;VesSha256_init(&io.sha);
    char magic[8];uint64_t abi[6],expected[6],counts[7];double scalar[2];int32_t ready;
#define ASR_GET(p,n) do{if(asr_get(&io,(p),(n)))goto done;}while(0)
    ASR_GET(magic,8);ASR_GET(abi,sizeof abi);asr_abi(expected);
    int resume=!memcmp(magic,"ASMRUN02",8);
    if((!resume && memcmp(magic,"ASMRUN01",8)) || memcmp(abi,expected,sizeof abi))goto done;
    ASR_GET(counts,sizeof counts);ASR_GET(scalar,sizeof scalar);ASR_GET(&ready,sizeof ready);
    if(!counts[1] || counts[6]>1 || ready<0 || ready>1 || !isfinite(scalar[0]) || !isfinite(scalar[1]))goto done;
    for(int i=0;i<6;i++)if(counts[i]>INT32_MAX)goto done;
    candidate.n_cubes=(size_t)counts[0];candidate.n_charts=candidate.cap_charts=(size_t)counts[1];
    candidate.n_rels=candidate.cap_rels=(size_t)counts[2];candidate.n_corr=candidate.cap_corr=(size_t)counts[3];
    candidate.n_layers=candidate.cap_layers=(size_t)counts[4];candidate.n_layer_hits=candidate.cap_layer_hits=(size_t)counts[5];
    candidate.cube_size=scalar[0];candidate.layer_d_global=scalar[1];candidate.continuity_ready=ready;
#define ASR_ALLOC(name,count) do{size_t n=(size_t)(count);if(n>io.remaining/sizeof *candidate.name)goto done;candidate.name=ARENA_ALLOC(run->arena,n*sizeof *candidate.name);ASR_GET(candidate.name,n*sizeof *candidate.name);}while(0)
    ASR_ALLOC(pile,candidate.n_cubes);
    for(size_t i=0;i<candidate.n_cubes;i++)if(!memchr(candidate.pile[i].path,0,sizeof candidate.pile[i].path) || !memchr(candidate.pile[i].cube_id,0,sizeof candidate.pile[i].cube_id))goto done;
    if(candidate.n_charts>io.remaining/(sizeof(AsmChart)+4))goto done;
    candidate.charts=ARENA_CALLOC(run->arena,candidate.n_charts,sizeof(AsmChart));
    if(has_uv)uv_present=ARENA_CALLOC(run->arena,candidate.n_charts,1);
    for(size_t c=0;c<candidate.n_charts;c++){
        AsmChart *row=candidate.charts+c;uint32_t present;ASR_GET(row,sizeof *row);ASR_GET(&present,sizeof present);
        if(row->id!=(int32_t)c || row->nv>INT32_MAX || row->nf>INT32_MAX || present>255 ||
            row->vid || row->xyz || row->nrm || row->uv || row->placed_uv || row->faces || row->boundary || row->repaired)goto done;
        int store=mode==0 || (mode==2 && (!keep || keep[c]));
        if(uv_present)uv_present[c]=(present&(1u<<3))!=0;
#define ASR_ARRAY(name,count,bit) do{if(present&(1u<<(bit))){size_t n=(size_t)(count);if(n>io.remaining/sizeof *row->name)goto done; \
            if(store){row->name=ARENA_ALLOC(run->arena,n*sizeof *row->name);ASR_GET(row->name,n*sizeof *row->name);} \
            else if(asr_skip(&io,n*sizeof *row->name,mode==1))goto done;}}while(0)
        ASR_ARRAY(vid,row->nv,0);ASR_ARRAY(xyz,3*row->nv,1);ASR_ARRAY(nrm,3*row->nv,2);ASR_ARRAY(uv,2*row->nv,3);
        ASR_ARRAY(placed_uv,2*row->nv,4);ASR_ARRAY(faces,3*row->nf,5);ASR_ARRAY(boundary,row->nv,6);ASR_ARRAY(repaired,row->nv,7);
#undef ASR_ARRAY
        if(!store){
            /* the record must still carry the geometry a stored chart would need */
            if((row->nv && !(present&(1u<<1))) || (row->nf && !(present&(1u<<5))))goto done;
            continue;
        }
        if((row->nv && !row->xyz) || (row->nf && !row->faces))goto done;
        for(size_t k=0;k<3*row->nf;k++)if(row->faces[k]<0 || (size_t)row->faces[k]>=row->nv)goto done;
    }
    ASR_ALLOC(rels,candidate.n_rels);ASR_ALLOC(corr,candidate.n_corr);
    ASR_ALLOC(layers,candidate.n_layers);ASR_ALLOC(layer_hits,candidate.n_layer_hits);
    if(counts[6])ASR_ALLOC(chart_layer_d,candidate.n_charts);
    if(resume){
        uint64_t bytes=0;ASR_GET(&bytes,sizeof bytes);if(bytes>io.remaining || bytes>SIZE_MAX)goto done;
        candidate.repair_resume_size=(size_t)bytes;ASR_ALLOC(repair_resume,candidate.repair_resume_size);
    }
    {char stored[64],digest[65];if(io.remaining || fread(stored,1,64,io.fp)!=64 || fgetc(io.fp)!=EOF || ferror(io.fp))goto done;
        VesSha256_hex(&io.sha,digest);if(mode!=1 && memcmp(stored,digest,64))goto done;}
    for(size_t i=0;i<candidate.n_rels;i++){
        const AsmRelation *r=candidate.rels+i;
        if(r->a<0 || r->b<0 || (size_t)r->a>=candidate.n_charts || (size_t)r->b>=candidate.n_charts || r->corr_count<0 ||
           (r->corr_count && (r->corr_first<0 || (size_t)r->corr_first>candidate.n_corr || (size_t)r->corr_count>candidate.n_corr-(size_t)r->corr_first)))goto done;
    }
    *run=candidate;if(has_uv)*has_uv=uv_present;rc=0;
done:
    if(fclose(io.fp))rc=-1;return rc;
#undef ASR_GET
#undef ASR_ALLOC
}
