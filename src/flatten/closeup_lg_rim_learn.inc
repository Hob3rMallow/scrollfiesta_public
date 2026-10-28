/* Persist only actual native rim-refusal pairs. Each later solve recomputes
 * their current separating direction; no vertex or material boundary is fixed. */
struct SmRimLearn {
    Arena_T arena;
    int32_t *pairs;
    size_t count,capacity,learned,calls,rows,near_duplicates,constant_rows;
    double seconds;
};

static int smrl_compare(const void *aa,const void *bb)
{
    const int32_t *a=aa,*b=bb;
    for(size_t k=0;k<4;k++)if(a[k]!=b[k])return a[k]<b[k]?-1:1;
    return 0;
}

static size_t smrl_lower(const SmRimLearn *s,const int32_t key[4])
{
    size_t lo=0,hi=s->count;
    while(lo<hi){size_t mid=lo+(hi-lo)/2;if(smrl_compare(s->pairs+4*mid,key)<0)lo=mid+1;else hi=mid;}
    return lo;
}

static int smrl_has_edge(const SmRim *rim,int32_t a,int32_t b)
{
    size_t lo=0,hi=rim->nb;
    while(lo<hi) {
        size_t mid=lo+(hi-lo)/2;int32_t u=rim->edges[2*mid],v=rim->edges[2*mid+1];
        if(u<a||(u==a&&v<b))lo=mid+1;else hi=mid;
    }
    return lo<rim->nb&&rim->edges[2*lo]==a&&rim->edges[2*lo+1]==b;
}

static int smrl_key(const SmRim *rim,const int32_t *input,int32_t key[4])
{
    for(size_t k=0;k<4;k++) {
        if(input[k]<0||(size_t)input[k]>=rim->nv)return -1;
        for(size_t j=0;j<k;j++)if(input[j]==input[k])return -1;
        key[k]=input[k];
    }
    for(size_t k=0;k<4;k+=2)if(key[k]>key[k+1]){int32_t swap=key[k];key[k]=key[k+1];key[k+1]=swap;}
    if(key[0]>key[2]||(key[0]==key[2]&&key[1]>key[3]))
        for(size_t k=0;k<2;k++){int32_t swap=key[k];key[k]=key[k+2];key[k+2]=swap;}
    return smrl_has_edge(rim,key[0],key[1])&&smrl_has_edge(rim,key[2],key[3])?0:-1;
}

static void smrl_boxes(const double *uv,const int32_t key[4],SmRimBox *a,SmRimBox *b)
{
    SmRimBox *boxes[2]={a,b};
    for(size_t k=0;k<2;k++) {
        SmRimBox *box=boxes[k];box->a=key[2*k];box->b=key[2*k+1];
        for(size_t xy=0;xy<2;xy++) {
            box->lo[xy]=fmin(uv[2*(size_t)box->a+xy],uv[2*(size_t)box->b+xy]);
            box->hi[xy]=fmax(uv[2*(size_t)box->a+xy],uv[2*(size_t)box->b+xy]);
        }
    }
    if(smr_compare_box(a,b)>0){SmRimBox swap=*a;*a=*b;*b=swap;}
}

/* Validate the entire observation before changing persistent membership.
 * New pairs returned to the caller live in its scratch arena, not the cache. */
static int smrl_learn(Arena_T scratch,SmRim *rim,const double *uv,const int32_t *pairs,size_t count,
    int32_t **added_pairs,size_t *added_count)
{
    if(!scratch||!rim||!rim->learning||!rim->learning->arena||!uv||!added_pairs||!added_count||
       (count&&!pairs)||count>SIZE_MAX/4/sizeof(int32_t))return -1;
    SmRimLearn *s=rim->learning;double started=ves_clock_sec();
    *added_pairs=NULL;*added_count=0;if(!count)return 0;
    int32_t *keys=ARENA_ALLOC(scratch,4*count*sizeof *keys);
    for(size_t k=0;k<count;k++) {
        int32_t *key=keys+4*k;if(smrl_key(rim,pairs+4*k,key))return -1;
        for(size_t j=0;j<4;j++)for(size_t xy=0;xy<2;xy++)if(!isfinite(uv[2*(size_t)key[j]+xy]))return -1;
        SmRimBox a,b;double normal[2],gap;smrl_boxes(uv,key,&a,&b);
        if(!smr_axis(uv,&a,&b,normal,&gap)||!isfinite(gap))return -1;
    }
    qsort(keys,count,4*sizeof *keys,smrl_compare);
    size_t unique=0,new_count=0;
    for(size_t k=0;k<count;k++)if(!unique||smrl_compare(keys+4*k,keys+4*(unique-1))) {
        if(unique!=k)memcpy(keys+4*unique,keys+4*k,4*sizeof *keys);unique++;
    }
    int32_t *fresh=ARENA_ALLOC(scratch,4*unique*sizeof *fresh);
    for(size_t k=0;k<unique;k++) {
        size_t at=smrl_lower(s,keys+4*k);
        if(at==s->count||smrl_compare(s->pairs+4*at,keys+4*k))memcpy(fresh+4*new_count++,keys+4*k,4*sizeof *fresh);
    }
    if(new_count>SIZE_MAX-s->count||(s->count+new_count)>SIZE_MAX/4/sizeof *s->pairs)return -1;
    size_t required=s->count+new_count;
    if(required>s->capacity) {
        size_t capacity=s->capacity?s->capacity:8;
        while(capacity<required) {
            if(capacity>SIZE_MAX/8/sizeof *s->pairs){capacity=required;break;}
            capacity*=2;
        }
        int32_t *next=ARENA_ALLOC(s->arena,4*capacity*sizeof *next);
        if(s->count)memcpy(next,s->pairs,4*s->count*sizeof *next);
        s->pairs=next;s->capacity=capacity;
    }
    for(size_t k=0;k<new_count;k++) {
        size_t at=smrl_lower(s,fresh+4*k);
        memmove(s->pairs+4*(at+1),s->pairs+4*at,4*(s->count-at)*sizeof *s->pairs);
        memcpy(s->pairs+4*at,fresh+4*k,4*sizeof *s->pairs);s->count++;
    }
    s->learned+=new_count;s->seconds+=ves_clock_sec()-started;
    *added_pairs=fresh;*added_count=new_count;return 0;
}

static int smrl_bounds(Arena_T arena,SmRim *rim,const int32_t *map,const double *uv,SmBound **rows,size_t *count)
{
    if(!arena||!rim||!rim->learning||!map||!uv||!rows||!count)return -1;
    SmRimLearn *s=rim->learning;double started=ves_clock_sec();
    if(s->count>(SIZE_MAX-*count)/4||(*count+4*s->count)>SIZE_MAX/sizeof(SmBound))return -1;
    SmBound *result=ARENA_ALLOC(arena,(*count+4*s->count)*sizeof *result);size_t at=*count,duplicates=0,constant=0;
    if(*count)memcpy(result,*rows,*count*sizeof *result);
    for(size_t k=0;k<s->count;k++) {
        SmRimBox a,b;SmBound pair[4];double gap=DBL_MAX;size_t folded=0;
        smrl_boxes(uv,s->pairs+4*k,&a,&b);
        int n=smr_pair(uv,map,&a,&b,DBL_MAX,pair,&folded,&gap);if(n<0)return -1;
        int box_near=b.lo[0]<=a.hi[0]+1e-4&&b.hi[1]>=a.lo[1]-1e-4&&b.lo[1]<=a.hi[1]+1e-4;
        /* The existing sweep already emitted every near-pair row. */
        if(box_near&&gap<=1e-4){duplicates++;continue;}
        memcpy(result+at,pair,(size_t)n*sizeof *pair);at+=(size_t)n;constant+=folded;
    }
    s->calls++;s->rows+=at-*count;s->near_duplicates+=duplicates;s->constant_rows+=constant;
    s->seconds+=ves_clock_sec()-started;*rows=result;*count=at;return 0;
}

static void smrl_close(SmRimLearn *s)
{
    if(!s)return;if(s->arena)Arena_dispose(&s->arena);s->pairs=NULL;s->capacity=0;
}
