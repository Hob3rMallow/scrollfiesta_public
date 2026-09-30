/* assembly_cases.causal_core: a shortest source path and every cyclic block
 * traversed by it. Other emitted faces remain immutable collision context. */
typedef struct AprBlocks {
    const size_t *offset;const int32_t *neighbor;const size_t *edge;
    size_t *visit,*low,*stack_a,*stack_b,clock,top;
    const uint8_t *path;uint8_t *core,*block;
    size_t charts;int same,first;
} AprBlocks;
static void apr_blocks_visit(AprBlocks *b,size_t vertex,size_t parent_edge)
{
    b->visit[vertex]=b->low[vertex]=++b->clock;
    for(size_t at=b->offset[vertex];at<b->offset[vertex+1];at++){
        size_t next=(size_t)b->neighbor[at],edge=b->edge[at];if(edge==parent_edge)continue;
        if(!b->visit[next]){
            size_t base=b->top;b->stack_a[b->top]=vertex;b->stack_b[b->top++]=next;
            apr_blocks_visit(b,next,edge);if(b->low[next]<b->low[vertex])b->low[vertex]=b->low[next];
            if(b->low[next]>=b->visit[vertex]){
                memset(b->block,0,b->charts);size_t count=0,hits=0;
                while(b->top>base){b->top--;b->block[b->stack_a[b->top]]=b->block[b->stack_b[b->top]]=1;}
                for(size_t c=0;c<b->charts;c++)if(b->block[c]){count++;hits+=b->path[c]!=0;}
                if(count>2 && (b->same?b->block[b->first]:hits>=2))for(size_t c=0;c<b->charts;c++)b->core[c]|=b->block[c];
            }
        }else if(b->visit[next]<b->visit[vertex]){
            b->stack_a[b->top]=vertex;b->stack_b[b->top++]=next;
            if(b->visit[next]<b->low[vertex])b->low[vertex]=b->visit[next];
        }
    }
}
static size_t apr_causal_core(const AsmField *f,size_t first,size_t last,uint8_t *core,int hops)
{
    size_t nc=f->run->n_charts;memset(core,0,nc);core[first]=core[last]=1;
    Arena_T arena=Arena_new();size_t *offset=ARENA_CALLOC(arena,nc+1,sizeof(size_t));
    for(size_t i=0;i<f->n_seams;i++){const AsmFieldSeam *s=f->seams+i;if(s->required && s->complete){offset[s->a+1]++;offset[s->b+1]++;}}
    for(size_t c=1;c<=nc;c++)offset[c]+=offset[c-1];
    int32_t *neighbor=ARENA_ALLOC(arena,offset[nc]*sizeof(int32_t));size_t *edges=ARENA_ALLOC(arena,offset[nc]*sizeof(size_t));
    size_t *fill=ARENA_ALLOC(arena,nc*sizeof(size_t));memcpy(fill,offset,nc*sizeof(size_t));
    for(size_t i=0;i<f->n_seams;i++){const AsmFieldSeam *s=f->seams+i;if(s->required && s->complete){
        size_t a=fill[s->a]++,b=fill[s->b]++;neighbor[a]=s->b;neighbor[b]=s->a;edges[a]=edges[b]=i;
    }}
    size_t *parent=ARENA_ALLOC(arena,nc*sizeof(size_t)),*queue=ARENA_ALLOC(arena,nc*sizeof(size_t));
    for(size_t c=0;c<nc;c++)parent[c]=SIZE_MAX;
    size_t head=0,tail=0;parent[first]=first;queue[tail++]=first;
    while(head<tail && parent[last]==SIZE_MAX){size_t c=queue[head++];for(size_t at=offset[c];at<offset[c+1];at++){
        size_t d=(size_t)neighbor[at];if(parent[d]==SIZE_MAX){parent[d]=c;queue[tail++]=d;}
    }}
    if(parent[last]!=SIZE_MAX){
        uint8_t *path=ARENA_CALLOC(arena,nc,1);size_t c=last;
        for(;;){path[c]=core[c]=1;if(c==first)break;c=parent[c];}
        AprBlocks blocks={0};blocks.offset=offset;blocks.neighbor=neighbor;blocks.edge=edges;blocks.charts=nc;
        blocks.visit=ARENA_CALLOC(arena,nc,sizeof(size_t));blocks.low=ARENA_CALLOC(arena,nc,sizeof(size_t));
        blocks.stack_a=ARENA_ALLOC(arena,offset[nc]*sizeof(size_t));blocks.stack_b=ARENA_ALLOC(arena,offset[nc]*sizeof(size_t));
        blocks.path=path;blocks.core=core;blocks.block=ARENA_CALLOC(arena,nc,1);blocks.same=first==last;blocks.first=(int)first;
        apr_blocks_visit(&blocks,first,SIZE_MAX);
    }
    /* keep the part of the core within reach of the seam: charts at most
     * `hops` required seams from either endpoint (0 = the whole core) */
    if(hops>0){
        size_t *dist=ARENA_ALLOC(arena,nc*sizeof(size_t));for(size_t c=0;c<nc;c++)dist[c]=SIZE_MAX;
        head=tail=0;dist[first]=0;queue[tail++]=first;if(last!=first){dist[last]=0;queue[tail++]=last;}
        while(head<tail){size_t c=queue[head++];if(dist[c]>=(size_t)hops)continue;
            for(size_t at=offset[c];at<offset[c+1];at++){size_t d=(size_t)neighbor[at];if(dist[d]==SIZE_MAX){dist[d]=dist[c]+1;queue[tail++]=d;}}}
        for(size_t c=0;c<nc;c++)if(core[c] && dist[c]==SIZE_MAX)core[c]=0;
    }
    size_t count=0;for(size_t c=0;c<nc;c++)count+=core[c]!=0;
    Arena_dispose(&arena);return count;
}
