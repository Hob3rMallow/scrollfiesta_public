/* Original edge-connected material topology and a boundary BVH per patch.
 * The face BVH still enumerates ALL contacts. Material labels only select
 * which true boundary defines the distance band; they never filter pairs. */
typedef struct AcMaterialEdge { int32_t a,b;size_t face; } AcMaterialEdge;
typedef struct AcMaterialNode {
    double lo[2],hi[2];
    size_t begin,end;
    int32_t left,right;
    int active;
} AcMaterialNode;
typedef struct AcMaterial {
    AcMaterialEdge *edge;
    AcMaterialNode *node;
    int32_t *component,*root;
    uint8_t *strict;
    size_t *order,*offset,edges,nodes,capacity,components;
    const double *uv;
} AcMaterial;

static int acm_key_order(const void *a,const void *b)
{
    uint64_t x=*(const uint64_t *)a,y=*(const uint64_t *)b;return(x>y)-(x<y);
}

static int acm_node_box(AcMaterial *m,int32_t id,const double *uv)
{
    AcMaterialNode *node=m->node+id;
    for(int d=0;d<2;d++){node->lo[d]=DBL_MAX;node->hi[d]=-DBL_MAX;}
    for(size_t i=node->begin;i<node->end;i++){
        const AcMaterialEdge *edge=m->edge+m->order[i];
        for(int side=0;side<2;side++)for(int d=0;d<2;d++){
            double x=uv[2*(size_t)(side?edge->b:edge->a)+d];if(!isfinite(x))return 0;
            node->lo[d]=ac_min(node->lo[d],x);node->hi[d]=ac_max(node->hi[d],x);
        }
    }
    return 1;
}

static double acm_center(const AcMaterial *m,size_t id,int axis)
{
    const AcMaterialEdge *edge=m->edge+id;
    return .5*m->uv[2*(size_t)edge->a+axis]+.5*m->uv[2*(size_t)edge->b+axis];
}

static int acm_less(const AcMaterial *m,size_t a,size_t b,int axis)
{
    double x=acm_center(m,a,axis),y=acm_center(m,b,axis);return x<y || (x==y && a<b);
}

static void acm_select(AcMaterial *m,size_t begin,size_t end,size_t mid,int axis)
{
    while(end-begin>1){
        size_t pivot=m->order[begin+(end-begin)/2],i=begin,j=end-1;
        for(;;){
            while(acm_less(m,m->order[i],pivot,axis))i++;
            while(acm_less(m,pivot,m->order[j],axis))j--;
            if(i>=j)break;size_t swap=m->order[i];m->order[i++]=m->order[j];m->order[j--]=swap;
        }
        size_t split=i;if(split==begin)split++;
        if(mid<split)end=split;else begin=split;
    }
}

static int32_t acm_build(AcMaterial *m,size_t begin,size_t end)
{
    if(m->nodes>=m->capacity)return -1;int32_t id=(int32_t)m->nodes++;AcMaterialNode *node=m->node+id;
    node->begin=begin;node->end=end;node->left=node->right=-1;node->active=1;
    if(!acm_node_box(m,id,m->uv))return -1;
    if(end-begin>8){
        size_t mid=begin+(end-begin)/2;int axis=node->hi[1]-node->lo[1]>node->hi[0]-node->lo[0];
        acm_select(m,begin,end,mid,axis);node->left=acm_build(m,begin,mid);node->right=acm_build(m,mid,end);
        if(node->left<0 || node->right<0)return -1;
    }
    return id;
}

static int acm_refit(AcMaterial *m,const double *uv,const uint8_t *active)
{
    if(!m || !uv)return 0;m->uv=uv;
    for(size_t i=m->nodes;i-->0;){
        AcMaterialNode *node=m->node+i;node->active=0;
        if(node->left<0){
            for(size_t k=node->begin;k<node->end && !node->active;k++)node->active=!active || active[m->edge[m->order[k]].face];
            if(node->active && !acm_node_box(m,(int32_t)i,uv))return 0;
        }else{
            const AcMaterialNode *a=m->node+node->left,*b=m->node+node->right;node->active=a->active || b->active;
            if(node->active)for(int d=0;d<2;d++){node->lo[d]=ac_min(a->lo[d],b->lo[d]);node->hi[d]=ac_max(a->hi[d],b->hi[d]);}
        }
    }
    return 1;
}

static AcMaterial *acm_new(Arena_T arena,size_t nv,const int32_t *faces,size_t nf,const double *uv,const int32_t *material_vertex)
{
    if(!arena || !faces || !uv || !nf || nf>INT32_MAX || nf>SIZE_MAX/3/sizeof(uint64_t) || nv>INT32_MAX)return NULL;
    Arena_T scratch=Arena_new();AcMaterial *result=NULL;
    size_t *offset=ARENA_CALLOC(scratch,nv+1,sizeof *offset),*fill=ARENA_CALLOC(scratch,nv,sizeof *fill);
    uint64_t *key=ARENA_ALLOC(scratch,3*nf*sizeof *key);UnionFind uf=UF_new(scratch,(int32_t)nf);
    for(size_t f=0;f<nf;f++)for(int k=0;k<3;k++){
        int32_t a=faces[3*f+k],b=faces[3*f+(k+1)%3];
        if(a<0 || b<0 || (size_t)a>=nv || (size_t)b>=nv || a==b)goto done;
        if(material_vertex){a=material_vertex[a];b=material_vertex[b];if(a==b)goto done;}
        offset[(size_t)(a<b?a:b)+1]++;
    }
    for(size_t v=0;v<nv;v++)offset[v+1]+=offset[v];
    for(size_t f=0;f<nf;f++)for(int k=0;k<3;k++){
        int32_t a=faces[3*f+k],b=faces[3*f+(k+1)%3];
        if(material_vertex){a=material_vertex[a];b=material_vertex[b];}
        int32_t lo=a<b?a:b,hi=a<b?b:a;
        key[offset[lo]+fill[lo]++]=((uint64_t)(uint32_t)hi<<32)|((uint64_t)f<<1)|(uint64_t)(a<b);
    }
    size_t boundary=0;
    for(size_t v=0;v<nv;v++){
        uint64_t *keys=key+offset[v];size_t count=offset[v+1]-offset[v];
        /* Normal mesh valence is small; only exceptional vertices need a
         * generic sort. This avoids sorting all 3*nf edges globally. */
        if(count<32){for(size_t i=1;i<count;i++){uint64_t value=keys[i];size_t j=i;while(j && keys[j-1]>value){keys[j]=keys[j-1];j--;}keys[j]=value;}}
        else qsort(keys,count,sizeof *keys,acm_key_order);
        for(size_t i=0;i<count;){
            size_t j=i+1;while(j<count && (keys[j]>>32)==(keys[i]>>32))j++;
            if(j-i==1)boundary++;
            for(size_t k=i+1;k<j;k++)uf_union(&uf,(int32_t)((uint32_t)keys[i]>>1),(int32_t)((uint32_t)keys[k]>>1));
            i=j;
        }
    }
    if(boundary>INT32_MAX/2 || boundary>SIZE_MAX/2/sizeof(AcMaterialNode))goto done;
    AcMaterial *m=ARENA_CALLOC(arena,1,sizeof *m);m->uv=uv;m->edges=boundary;
    m->component=ARENA_ALLOC(arena,nf*sizeof *m->component);
    for(size_t i=0;i<nf;i++)uf.rank[i]=-1;
    for(size_t i=0;i<nf;i++){
        int32_t root=uf_find(&uf,(int32_t)i);if(uf.rank[root]<0)uf.rank[root]=(int32_t)m->components++;
        m->component[i]=uf.rank[root];
    }
    m->strict=ARENA_CALLOC(arena,m->components,1);m->offset=ARENA_CALLOC(arena,m->components+1,sizeof *m->offset);
    for(size_t v=0;v<nv;v++)for(size_t i=offset[v];i<offset[v+1];){
        size_t j=i+1;while(j<offset[v+1] && (key[j]>>32)==(key[i]>>32))j++;
        int32_t c=m->component[(uint32_t)key[i]>>1];
        if(j-i==1)m->offset[c+1]++;
        else if(j-i!=2 || (key[i]&1)==(key[i+1]&1))m->strict[c]=1;
        i=j;
    }
    for(size_t c=0;c<m->components;c++)m->offset[c+1]+=m->offset[c];
    m->edge=ARENA_ALLOC(arena,(boundary?boundary:1)*sizeof *m->edge);
    m->order=ARENA_ALLOC(arena,(boundary?boundary:1)*sizeof *m->order);
    size_t *next=ARENA_CALLOC(scratch,m->components,sizeof *next);
    for(size_t v=0;v<nv;v++)for(size_t i=offset[v];i<offset[v+1];){
        size_t j=i+1;while(j<offset[v+1] && (key[j]>>32)==(key[i]>>32))j++;
        if(j-i==1){
            size_t face=(uint32_t)key[i]>>1;int32_t c=m->component[face];size_t at=m->offset[c]+next[c]++;
            int32_t a=(int32_t)v,b=(int32_t)(key[i]>>32);
            if(material_vertex){
                int found=0;
                for(int k=0;k<3;k++){
                    int32_t x=faces[3*face+k],y=faces[3*face+(k+1)%3];
                    int32_t u=material_vertex[x],w=material_vertex[y];
                    if((u==a && w==b) || (u==b && w==a)){m->edge[at]=(AcMaterialEdge){x,y,face};found++;}
                }
                if(found!=1)goto done;
            }else m->edge[at]=(AcMaterialEdge){a,b,face};
            m->order[at]=at;
        }
        i=j;
    }
    m->capacity=boundary?2*boundary:1;m->node=ARENA_ALLOC(arena,m->capacity*sizeof *m->node);
    m->root=ARENA_ALLOC(arena,m->components*sizeof *m->root);
    for(size_t c=0;c<m->components;c++){
        m->root[c]=-1;if(m->offset[c]==m->offset[c+1]){m->strict[c]=1;continue;}
        m->root[c]=acm_build(m,m->offset[c],m->offset[c+1]);if(m->root[c]<0)goto done;
    }
    result=m;
done:
    Arena_dispose(&scratch);return result;
}

typedef struct AcMaterialGather {
    double *segments;
    double box[4],origin[2];
    const double *uv;
    size_t count;
    int complete;
} AcMaterialGather;

static void acm_gather_edge(const AcMaterial *m,size_t edge,AcMaterialGather *g)
{
    const AcMaterialEdge *e=m->edge+edge;double xy[4];
    for(int d=0;d<2;d++){
        xy[d]=g->uv[2*(size_t)e->a+d]-g->origin[d];xy[2+d]=g->uv[2*(size_t)e->b+d]-g->origin[d];
        if(ac_min(xy[d],xy[2+d])>g->box[2+d] || ac_max(xy[d],xy[2+d])<g->box[d])return;
    }
    if(g->count==ACB_SEGMENTS){g->complete=0;return;}
    memcpy(g->segments+4*g->count++,xy,sizeof xy);
}

static void acm_gather_node(const AcMaterial *m,int32_t id,AcMaterialGather *g)
{
    if(!g->complete || id<0)return;const AcMaterialNode *node=m->node+id;
    for(int d=0;d<2;d++)if(node->lo[d]-g->origin[d]>g->box[2+d] || node->hi[d]-g->origin[d]<g->box[d])return;
    if(node->left<0){for(size_t i=node->begin;i<node->end && g->complete;i++)acm_gather_edge(m,m->order[i],g);}
    else{acm_gather_node(m,node->left,g);acm_gather_node(m,node->right,g);}
}

static void acm_gather(const AcMaterial *m,size_t face,const double *uv,AcMaterialGather *g)
{
    int32_t c=m->component[face];g->uv=uv;
    if(uv==m->uv)acm_gather_node(m,m->root[c],g);
    else for(size_t i=m->offset[c];i<m->offset[c+1] && g->complete;i++)acm_gather_edge(m,i,g);
}

static int acm_selftest(void)
{
    Arena_T arena=Arena_new();int failures=0;
    const int32_t faces[12]={0,1,2,0,2,3,4,5,6,4,6,7};
    double uv[16]={0,0,20,0,20,20,0,20,14,0,34,0,34,20,14,20};
    AcMaterial *m=acm_new(arena,8,faces,4,uv,NULL);
    if(!m || m->components!=2 || m->edges!=8 || m->strict[0] || m->strict[1])failures++;
    if(m){
        double segments[4*ACB_SEGMENTS];AcMaterialGather g={segments,{13,5,21,15},{0,0},uv,0,1};
        acm_gather(m,0,uv,&g);acm_gather(m,2,uv,&g);
        if(!g.complete || g.count!=2)failures++;
        uint8_t moving[4]={0,0,1,1};for(int i=4;i<8;i++)uv[2*i]+=10;
        if(!acm_refit(m,uv,moving))failures++;
        g.count=0;acm_gather(m,0,uv,&g);acm_gather(m,2,uv,&g);if(g.count!=1)failures++;
        double previous[16]={0,0,20,0,20,20,0,20,14,0,34,0,34,20,14,20};
        g.count=0;acm_gather(m,0,previous,&g);acm_gather(m,2,previous,&g);if(g.count!=2)failures++;
    }
    /* Two patches meeting at a vertex do not acquire each other's interior
     * edges; contradictory face winding cannot certify a tolerance band. */
    {
        const int32_t point_touch[6]={0,1,2,0,4,3},duplicate[6]={0,1,2,0,1,2};
        double xy[10]={0,0,2,0,0,2,-2,0,0,-2};
        AcMaterial *point=acm_new(arena,5,point_touch,2,xy,NULL),*bad=acm_new(arena,5,duplicate,2,xy,NULL);
        if(!point || point->components!=2 || point->edges!=6 || !bad || !bad->strict[0])failures++;
    }
    fprintf(stderr,"  original material boundary topology/index: %s (%d failures)\n",failures?"FAIL":"ok",failures);
    Arena_dispose(&arena);return failures;
}
