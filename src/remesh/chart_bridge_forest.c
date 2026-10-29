/* chart_bridge_forest.c -- chart-level transactional BPA bridge selection. */
#include "chart_bridge_forest.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct { int32_t lo, hi, face; int8_t dir; } CbfEdgeFace;
typedef struct { int32_t lo, hi; uint32_t count; } CbfInputEdge;
typedef struct { int32_t patch, value; } CbfPatchValue;
typedef struct { int32_t patch, side, vertex; } CbfAttachVertex;
typedef struct { int32_t patch, center, a, b; } CbfLinkEdge;
typedef struct {
    int32_t root, chart_a, chart_b;
    size_t faces;
    double area;
} CbfCandidate;

typedef struct {
    int32_t chart_a, chart_b;
    uint8_t parity;
    size_t faces;
    double area;
} CbfOrientConstraint;

typedef struct {
    int32_t root;
    size_t charts, input_faces, bridge_faces;
    double bbox_min[3], bbox_max[3];
} CbfGroupDiag;

typedef struct {
    int32_t lo, hi;
    uint32_t count;
    int32_t orient_sum;
    uint8_t used;
} CbfEdgeSlot;

typedef struct {
    CbfEdgeSlot *slot;
    size_t cap;
} CbfEdgeMap;

enum {
    CBF_REJECT_NONE = 0,
    CBF_REJECT_ATTACHMENT = 1,
    CBF_REJECT_MULTIPATH = 2,
    CBF_REJECT_EDGE = 3,
    CBF_REJECT_ORIENTATION = 4
};

static int32_t cbf_find(int32_t *parent, int32_t x)
{
    while (parent[x] != x) {
        parent[x] = parent[parent[x]];
        x = parent[x];
    }
    return x;
}

static void cbf_union(int32_t *parent, uint8_t *rank, int32_t a, int32_t b)
{
    a = cbf_find(parent, a); b = cbf_find(parent, b);
    if (a == b) return;
    if (rank[a] < rank[b]) { int32_t t = a; a = b; b = t; }
    parent[b] = a;
    if (rank[a] == rank[b] && rank[a] < UINT8_MAX) rank[a]++;
}

static int cbf_edge_face_cmp(const void *pa, const void *pb)
{
    const CbfEdgeFace *a = (const CbfEdgeFace *)pa;
    const CbfEdgeFace *b = (const CbfEdgeFace *)pb;
    if (a->lo != b->lo) return a->lo < b->lo ? -1 : 1;
    if (a->hi != b->hi) return a->hi < b->hi ? -1 : 1;
    return a->face < b->face ? -1 : (a->face > b->face ? 1 : 0);
}

static int cbf_patch_edge_cmp(const void *pa, const void *pb)
{
    const CbfEdgeFace *a = (const CbfEdgeFace *)pa;
    const CbfEdgeFace *b = (const CbfEdgeFace *)pb;
    if (a->face != b->face) return a->face < b->face ? -1 : 1;
    if (a->lo != b->lo) return a->lo < b->lo ? -1 : 1;
    if (a->hi != b->hi) return a->hi < b->hi ? -1 : 1;
    return 0;
}

static uint64_t cbf_edge_hash(int32_t lo, int32_t hi)
{
    uint64_t x = ((uint64_t)(uint32_t)lo << 32) | (uint32_t)hi;
    x ^= x >> 30; x *= UINT64_C(0xbf58476d1ce4e5b9);
    x ^= x >> 27; x *= UINT64_C(0x94d049bb133111eb);
    return x ^ (x >> 31);
}

static int cbf_edge_map_init(CbfEdgeMap *map, size_t edge_hint)
{
    size_t cap = 16;
    if (!map || edge_hint > SIZE_MAX / 2) return -1;
    while (cap < edge_hint * 2) {
        if (cap > SIZE_MAX / 2) return -1;
        cap *= 2;
    }
    map->slot = (CbfEdgeSlot *)calloc(cap, sizeof(CbfEdgeSlot));
    if (!map->slot) return -1;
    map->cap = cap;
    return 0;
}

static CbfEdgeSlot *cbf_edge_map_slot(CbfEdgeMap *map,
                                      int32_t lo, int32_t hi)
{
    size_t mask = map->cap - 1;
    size_t at = (size_t)cbf_edge_hash(lo, hi) & mask;
    while (map->slot[at].used &&
           (map->slot[at].lo != lo || map->slot[at].hi != hi))
        at = (at + 1) & mask;
    return &map->slot[at];
}

static uint32_t cbf_edge_map_get(CbfEdgeMap *map, int32_t lo, int32_t hi)
{
    CbfEdgeSlot *slot = cbf_edge_map_slot(map, lo, hi);
    return slot->used ? slot->count : 0;
}

static int cbf_edge_map_add(CbfEdgeMap *map, int32_t lo, int32_t hi,
                            uint32_t count, int32_t orient_sum)
{
    CbfEdgeSlot *slot = cbf_edge_map_slot(map, lo, hi);
    if (!slot->used) {
        slot->used = 1;
        slot->lo = lo;
        slot->hi = hi;
        slot->count = 0;
        slot->orient_sum = 0;
    }
    if (count > 2 || slot->count > 2 - count) return -1;
    slot->count += count;
    slot->orient_sum += orient_sum;
    if ((slot->count == 1 &&
         !(slot->orient_sum == 1 || slot->orient_sum == -1)) ||
        (slot->count == 2 && slot->orient_sum != 0))
        return -1;
    return 0;
}

static int cbf_boundary_add(int32_t *n0, int32_t *n1,
                            int32_t v, int32_t other)
{
    if (n0[v] == other || n1[v] == other) return 0;
    if (n0[v] < 0) { n0[v] = other; return 0; }
    if (n1[v] < 0) { n1[v] = other; return 0; }
    return -1;
}

static int cbf_boundary_remove(int32_t *n0, int32_t *n1,
                               int32_t v, int32_t other)
{
    if (n0[v] == other) { n0[v] = n1[v]; n1[v] = -1; return 0; }
    if (n1[v] == other) { n1[v] = -1; return 0; }
    return -1;
}

static int cbf_boundary_same_loop(const int32_t *n0, const int32_t *n1,
                                  size_t nv, int32_t start, int32_t target)
{
    int32_t previous = -1, current = start;
    if (start == target) return 1;
    for (size_t step = 0; step <= nv; step++) {
        int32_t next;
        if (current < 0 || (size_t)current >= nv) return 0;
        next = n0[current] != previous ? n0[current] : n1[current];
        if (next < 0) return 0;
        if (next == target) return 1;
        previous = current;
        current = next;
        if (current == start) return 0;
    }
    return 0;
}

static int cbf_input_edge_cmp(const void *pa, const void *pb)
{
    const CbfInputEdge *a = (const CbfInputEdge *)pa;
    const CbfInputEdge *b = (const CbfInputEdge *)pb;
    if (a->lo != b->lo) return a->lo < b->lo ? -1 : 1;
    if (a->hi != b->hi) return a->hi < b->hi ? -1 : 1;
    return 0;
}

static int cbf_patch_value_cmp(const void *pa, const void *pb)
{
    const CbfPatchValue *a = (const CbfPatchValue *)pa;
    const CbfPatchValue *b = (const CbfPatchValue *)pb;
    if (a->patch != b->patch) return a->patch < b->patch ? -1 : 1;
    if (a->value != b->value) return a->value < b->value ? -1 : 1;
    return 0;
}

static int cbf_attach_vertex_cmp(const void *pa, const void *pb)
{
    const CbfAttachVertex *a = (const CbfAttachVertex *)pa;
    const CbfAttachVertex *b = (const CbfAttachVertex *)pb;
    if (a->patch != b->patch) return a->patch < b->patch ? -1 : 1;
    if (a->side != b->side) return a->side < b->side ? -1 : 1;
    if (a->vertex != b->vertex) return a->vertex < b->vertex ? -1 : 1;
    return 0;
}

static int cbf_link_edge_cmp(const void *pa, const void *pb)
{
    const CbfLinkEdge *a = (const CbfLinkEdge *)pa;
    const CbfLinkEdge *b = (const CbfLinkEdge *)pb;
    if (a->patch != b->patch) return a->patch < b->patch ? -1 : 1;
    if (a->center != b->center) return a->center < b->center ? -1 : 1;
    if (a->a != b->a) return a->a < b->a ? -1 : 1;
    if (a->b != b->b) return a->b < b->b ? -1 : 1;
    return 0;
}

static int cbf_candidate_cmp(const void *pa, const void *pb)
{
    const CbfCandidate *a = (const CbfCandidate *)pa;
    const CbfCandidate *b = (const CbfCandidate *)pb;
    if (a->area != b->area) return a->area > b->area ? -1 : 1;
    if (a->faces != b->faces) return a->faces > b->faces ? -1 : 1;
    return a->root < b->root ? -1 : (a->root > b->root ? 1 : 0);
}

static double cbf_face_area(const float *v,const int32_t *f);

static int cbf_orient_constraint_cmp(const void *pa, const void *pb)
{
    const CbfOrientConstraint *a=(const CbfOrientConstraint *)pa;
    const CbfOrientConstraint *b=(const CbfOrientConstraint *)pb;
    if(a->area!=b->area)return a->area>b->area?-1:1;
    if(a->faces!=b->faces)return a->faces>b->faces?-1:1;
    if(a->chart_a!=b->chart_a)return a->chart_a<b->chart_a?-1:1;
    return a->chart_b<b->chart_b?-1:(a->chart_b>b->chart_b?1:0);
}

static int32_t cbf_parity_find(int32_t *parent,uint8_t *parity,int32_t x)
{
    if(parent[x]!=x){
        int32_t old=parent[x];
        parent[x]=cbf_parity_find(parent,parity,old);
        parity[x]^=parity[old];
    }
    return parent[x];
}

/* Add x[a] xor x[b] == relation.  Return 1 for a new tree edge, 0 for
 * an already-implied relation, and -1 for a conflicting cycle. */
static int cbf_parity_union(int32_t *parent,uint8_t *rank,uint8_t *parity,
                            int32_t a,int32_t b,uint8_t relation)
{
    int32_t ra=cbf_parity_find(parent,parity,a);
    int32_t rb=cbf_parity_find(parent,parity,b);
    uint8_t link=(uint8_t)(parity[a]^parity[b]^relation);
    if(ra==rb)return (parity[a]^parity[b])==relation?0:-1;
    if(rank[ra]<rank[rb]){
        parent[ra]=rb;parity[ra]=link;
    }else{
        parent[rb]=ra;parity[rb]=link;
        if(rank[ra]==rank[rb]&&rank[ra]<UINT8_MAX)rank[ra]++;
    }
    return 1;
}

int ChartBridgeForest_orient_source(Arena_T arena,
                                    const float *verts,size_t nv,
                                    int32_t *faces,size_t nf,
                                    size_t input_nf,
                                    const int32_t *vertex_chart,
                                    size_t n_charts,
                                    ChartBridgeOrientationStats *stats)
{
    size_t nbf,nedges,nconstraint=0;
    int32_t *face_parent=NULL,*pair_a=NULL,*pair_b=NULL;
    uint8_t *face_rank=NULL,*face_valid=NULL,*patch_bad=NULL;
    size_t *patch_faces=NULL,*patch_vertices=NULL,*patch_edges=NULL;
    size_t *patch_boundary=NULL,*attach_a=NULL,*attach_b=NULL;
    double *patch_area=NULL;
    int8_t *req_a=NULL,*req_b=NULL;
    CbfEdgeFace *edge_face=NULL,*patch_edge=NULL;
    CbfPatchValue *patch_vertex=NULL;
    CbfOrientConstraint *constraint=NULL;
    int32_t *parent=NULL;
    uint8_t *rank=NULL,*parity=NULL,*root_toggle=NULL;
    size_t *chart_faces=NULL,*weight0=NULL,*weight1=NULL;
    CbfEdgeMap source_edges;
    ChartBridgeOrientationStats st;
    int rc=-1;

    (void)arena;
    memset(&st,0,sizeof st);
    memset(&source_edges,0,sizeof source_edges);
    if(stats)*stats=st;
    if(!verts||!faces||!vertex_chart||input_nf>nf||nv==0||n_charts==0||
       n_charts>(size_t)INT32_MAX)return -1;
    nbf=nf-input_nf;
    if(nbf==0)return 0;
    if(nbf>(size_t)INT32_MAX||nbf>SIZE_MAX/3)return -1;

    chart_faces=(size_t*)calloc(n_charts,sizeof(*chart_faces));
    if(!chart_faces)goto cleanup_orient;
    if(cbf_edge_map_init(&source_edges,input_nf*3+1)!=0)
        goto cleanup_orient;
    for(size_t f=0;f<input_nf;f++){
        const int32_t *tri=&faces[f*3];
        int32_t c;
        for(int k=0;k<3;k++)if(tri[k]<0||(size_t)tri[k]>=nv)
            goto cleanup_orient;
        c=vertex_chart[tri[0]];
        if(c<0||(size_t)c>=n_charts||vertex_chart[tri[1]]!=c||
           vertex_chart[tri[2]]!=c)goto cleanup_orient;
        chart_faces[c]++;
        for(int k=0;k<3;k++){
            int32_t a=tri[k],b=tri[(k+1)%3];
            int32_t lo=a<b?a:b,hi=a<b?b:a;
            if(cbf_edge_map_add(&source_edges,lo,hi,1,a<b?1:-1)!=0)
                goto cleanup_orient;
        }
    }

    face_parent=(int32_t*)malloc(nbf*sizeof(*face_parent));
    face_rank=(uint8_t*)calloc(nbf,1);
    face_valid=(uint8_t*)calloc(nbf,1);
    pair_a=(int32_t*)malloc(nbf*sizeof(*pair_a));
    pair_b=(int32_t*)malloc(nbf*sizeof(*pair_b));
    edge_face=(CbfEdgeFace*)malloc(nbf*3*sizeof(*edge_face));
    if(!face_parent||!face_rank||!face_valid||!pair_a||!pair_b||!edge_face)
        goto cleanup_orient;
    for(size_t i=0;i<nbf;i++){
        const int32_t *tri=&faces[(input_nf+i)*3];
        int32_t chart[3];int nc=0;
        face_parent[i]=(int32_t)i;pair_a[i]=pair_b[i]=-1;
        for(int k=0;k<3;k++){
            int32_t v=tri[k],c;int seen=0;
            if(v<0||(size_t)v>=nv){nc=4;break;}
            c=vertex_chart[v];
            if(c<0||(size_t)c>=n_charts){nc=4;break;}
            for(int q=0;q<nc;q++)if(chart[q]==c)seen=1;
            if(!seen&&nc<3)chart[nc++]=c;
        }
        if(nc==2){
            if(chart[0]>chart[1]){int32_t t=chart[0];chart[0]=chart[1];chart[1]=t;}
            pair_a[i]=chart[0];pair_b[i]=chart[1];face_valid[i]=1;
        }
        for(int k=0;k<3;k++){
            int32_t a=tri[k],b=tri[(k+1)%3];
            edge_face[i*3+(size_t)k].lo=a<b?a:b;
            edge_face[i*3+(size_t)k].hi=a<b?b:a;
            edge_face[i*3+(size_t)k].face=(int32_t)i;
            edge_face[i*3+(size_t)k].dir=a<b?1:-1;
        }
    }
    nedges=nbf*3;
    qsort(edge_face,nedges,sizeof(*edge_face),cbf_edge_face_cmp);
    for(size_t i=0;i<nedges;){
        size_t j=i+1;
        while(j<nedges&&edge_face[j].lo==edge_face[i].lo&&
              edge_face[j].hi==edge_face[i].hi)j++;
        for(size_t a=i;a<j;a++)for(size_t b=a+1;b<j;b++){
            int32_t fa=edge_face[a].face,fb=edge_face[b].face;
            if(face_valid[fa]&&face_valid[fb]&&pair_a[fa]==pair_a[fb]&&
               pair_b[fa]==pair_b[fb])
                cbf_union(face_parent,face_rank,fa,fb);
        }
        i=j;
    }
    for(size_t i=0;i<nbf;i++)face_parent[i]=cbf_find(face_parent,(int32_t)i);

    patch_faces=(size_t*)calloc(nbf,sizeof(*patch_faces));
    patch_vertices=(size_t*)calloc(nbf,sizeof(*patch_vertices));
    patch_edges=(size_t*)calloc(nbf,sizeof(*patch_edges));
    patch_boundary=(size_t*)calloc(nbf,sizeof(*patch_boundary));
    attach_a=(size_t*)calloc(nbf,sizeof(*attach_a));
    attach_b=(size_t*)calloc(nbf,sizeof(*attach_b));
    patch_area=(double*)calloc(nbf,sizeof(*patch_area));
    patch_bad=(uint8_t*)calloc(nbf,1);
    req_a=(int8_t*)calloc(nbf,1);req_b=(int8_t*)calloc(nbf,1);
    patch_edge=(CbfEdgeFace*)malloc(nbf*3*sizeof(*patch_edge));
    patch_vertex=(CbfPatchValue*)malloc(nbf*3*sizeof(*patch_vertex));
    constraint=(CbfOrientConstraint*)malloc(nbf*sizeof(*constraint));
    if(!patch_faces||!patch_vertices||!patch_edges||!patch_boundary||
       !attach_a||!attach_b||!patch_area||!patch_bad||!req_a||!req_b||
       !patch_edge||!patch_vertex||!constraint)goto cleanup_orient;
    for(size_t i=0;i<nbf;i++){
        const int32_t *tri=&faces[(input_nf+i)*3];
        int32_t r=face_valid[i]?face_parent[i]:-1;
        if(r>=0){patch_faces[r]++;patch_area[r]+=cbf_face_area(verts,tri);}
        for(int k=0;k<3;k++){
            int32_t a=tri[k],b=tri[(k+1)%3];
            patch_vertex[i*3+(size_t)k].patch=r;
            patch_vertex[i*3+(size_t)k].value=r>=0?tri[k]:-1;
            patch_edge[i*3+(size_t)k].face=r;
            patch_edge[i*3+(size_t)k].lo=r>=0?(a<b?a:b):-1;
            patch_edge[i*3+(size_t)k].hi=r>=0?(a<b?b:a):-1;
            patch_edge[i*3+(size_t)k].dir=r>=0?(a<b?1:-1):0;
        }
    }
    qsort(patch_vertex,nbf*3,sizeof(*patch_vertex),cbf_patch_value_cmp);
    for(size_t i=0;i<nbf*3;){
        size_t j=i+1;
        while(j<nbf*3&&patch_vertex[j].patch==patch_vertex[i].patch&&
              patch_vertex[j].value==patch_vertex[i].value)j++;
        if(patch_vertex[i].patch>=0)patch_vertices[patch_vertex[i].patch]++;
        i=j;
    }
    qsort(patch_edge,nbf*3,sizeof(*patch_edge),cbf_patch_edge_cmp);
    for(size_t i=0;i<nbf*3;){
        size_t j=i+1;int dirsum=patch_edge[i].dir;
        while(j<nbf*3&&patch_edge[j].face==patch_edge[i].face&&
              patch_edge[j].lo==patch_edge[i].lo&&
              patch_edge[j].hi==patch_edge[i].hi){
            dirsum+=patch_edge[j].dir;j++;
        }
        if(patch_edge[i].face>=0){
            int32_t r=patch_edge[i].face;
            size_t count=j-i;
            patch_edges[r]++;
            if(count==1)patch_boundary[r]++;
            if(count>2||(count==2&&dirsum!=0))patch_bad[r]=1;
            if(count==1){
                CbfEdgeSlot *slot=cbf_edge_map_slot(
                    &source_edges,patch_edge[i].lo,patch_edge[i].hi);
                if(slot->used&&slot->count==1){
                    int32_t c0=vertex_chart[patch_edge[i].lo];
                    int32_t c1=vertex_chart[patch_edge[i].hi];
                    int8_t req=(int8_t)(dirsum==-slot->orient_sum?1:-1);
                    int8_t *dst=NULL;
                    size_t *nattach=NULL;
                    if(c0==c1&&c0==pair_a[r]){dst=&req_a[r];nattach=&attach_a[r];}
                    else if(c0==c1&&c0==pair_b[r]){dst=&req_b[r];nattach=&attach_b[r];}
                    if(dst){
                        (*nattach)++;
                        if(*dst==0)*dst=req;
                        else if(*dst!=req)patch_bad[r]=1;
                    }
                }
            }
        }
        i=j;
    }
    for(size_t r=0;r<nbf;r++)if(face_valid[r]&&face_parent[r]==(int32_t)r){
        int topology_ok=!patch_bad[r]&&patch_boundary[r]>0&&
            (int64_t)patch_vertices[r]-(int64_t)patch_edges[r]+
            (int64_t)patch_faces[r]==1;
        if(!topology_ok||attach_a[r]==0||attach_b[r]==0||
           req_a[r]==0||req_b[r]==0)continue;
        constraint[nconstraint++]=(CbfOrientConstraint){
            pair_a[r],pair_b[r],(uint8_t)(req_a[r]!=req_b[r]),
            patch_faces[r],patch_area[r]};
    }
    st.candidate_patches=nconstraint;
    qsort(constraint,nconstraint,sizeof(*constraint),cbf_orient_constraint_cmp);

    parent=(int32_t*)malloc(n_charts*sizeof(*parent));
    rank=(uint8_t*)calloc(n_charts,1);
    parity=(uint8_t*)calloc(n_charts,1);
    root_toggle=(uint8_t*)calloc(n_charts,1);
    weight0=(size_t*)calloc(n_charts,sizeof(*weight0));
    weight1=(size_t*)calloc(n_charts,sizeof(*weight1));
    if(!parent||!rank||!parity||!root_toggle||!weight0||!weight1)
        goto cleanup_orient;
    for(size_t c=0;c<n_charts;c++)parent[c]=(int32_t)c;
    for(size_t i=0;i<nconstraint;i++){
        int u=cbf_parity_union(parent,rank,parity,constraint[i].chart_a,
                               constraint[i].chart_b,constraint[i].parity);
        if(u<0)st.constraints_conflicted++;
        else st.constraints_used++;
    }
    for(size_t c=0;c<n_charts;c++){
        int32_t root=cbf_parity_find(parent,parity,(int32_t)c);
        if(parity[c])weight1[root]+=chart_faces[c];
        else weight0[root]+=chart_faces[c];
    }
    for(size_t c=0;c<n_charts;c++)if(parent[c]==(int32_t)c)
        root_toggle[c]=(uint8_t)(weight0[c]<weight1[c]);
    for(size_t c=0;c<n_charts;c++){
        int32_t root=cbf_parity_find(parent,parity,(int32_t)c);
        if((parity[c]^root_toggle[root])&&chart_faces[c]>0)
            st.charts_flipped++;
    }
    for(size_t f=0;f<input_nf;f++){
        int32_t c=vertex_chart[faces[f*3]];
        int32_t root=cbf_parity_find(parent,parity,c);
        if(parity[c]^root_toggle[root]){
            int32_t t=faces[f*3+1];faces[f*3+1]=faces[f*3+2];faces[f*3+2]=t;
            st.faces_flipped++;
        }
    }
    if(stats)*stats=st;
    rc=0;

cleanup_orient:
    free(weight1);free(weight0);free(root_toggle);free(parity);free(rank);
    free(parent);free(constraint);free(patch_vertex);free(patch_edge);
    free(req_b);free(req_a);free(patch_bad);free(patch_area);
    free(attach_b);free(attach_a);free(patch_boundary);free(patch_edges);
    free(patch_vertices);free(patch_faces);free(edge_face);free(pair_b);
    free(pair_a);free(face_valid);free(face_rank);free(face_parent);
    free(source_edges.slot);free(chart_faces);
    return rc;
}

static int cbf_group_diag_cmp(const void *pa, const void *pb)
{
    const CbfGroupDiag *a = (const CbfGroupDiag *)pa;
    const CbfGroupDiag *b = (const CbfGroupDiag *)pb;
    size_t af = a->input_faces + a->bridge_faces;
    size_t bf = b->input_faces + b->bridge_faces;
    if (af != bf) return af > bf ? -1 : 1;
    return a->root < b->root ? -1 : (a->root > b->root ? 1 : 0);
}

static const char *cbf_reject_name(uint8_t reason)
{
    switch (reason) {
    case CBF_REJECT_MULTIPATH: return "multipath";
    case CBF_REJECT_EDGE: return "edge";
    case CBF_REJECT_ORIENTATION: return "orientation";
    default: return "attachment";
    }
}

static const CbfInputEdge *cbf_input_edge_find(const CbfInputEdge *edge,
                                                size_t n,
                                                int32_t lo, int32_t hi)
{
    size_t left = 0, right = n;
    while (left < right) {
        size_t mid = left + (right - left) / 2;
        if (edge[mid].lo < lo || (edge[mid].lo == lo && edge[mid].hi < hi))
            left = mid + 1;
        else right = mid;
    }
    if (left < n && edge[left].lo == lo && edge[left].hi == hi)
        return &edge[left];
    return NULL;
}

static double cbf_face_area(const float *v, const int32_t *f)
{
    const float *a = &v[(size_t)f[0]*3];
    const float *b = &v[(size_t)f[1]*3];
    const float *c = &v[(size_t)f[2]*3];
    double u0=b[0]-a[0], u1=b[1]-a[1], u2=b[2]-a[2];
    double w0=c[0]-a[0], w1=c[1]-a[1], w2=c[2]-a[2];
    double n0=u1*w2-u2*w1, n1=u2*w0-u0*w2, n2=u0*w1-u1*w0;
    return 0.5*sqrt(n0*n0+n1*n1+n2*n2);
}

static int cbf_filter_round(Arena_T arena,
                            const float *verts, size_t nv,
                            const int32_t *faces, size_t nf,
                            size_t input_nf,
                            const int32_t *vertex_chart,
                            size_t n_charts,
                            int32_t **out_faces, size_t *out_nf,
                            uint8_t *out_bridge_kept,
                            int promote_support,
                            size_t support_begin,
                            size_t support_end,
                            ChartBridgeForestStats *stats)
{
    size_t nbf, nedges, nvalid = 0, ninput_edges;
    int32_t *face_parent = NULL, *face_chart_a = NULL, *face_chart_b = NULL;
    int32_t *face_chart_c = NULL;
    uint8_t *face_rank = NULL, *face_valid = NULL, *patch_valid = NULL;
    uint8_t *face_chart_count = NULL;
    uint8_t *patch_topology = NULL, *patch_reject = NULL;
    uint8_t *patch_link_bad = NULL;
    uint8_t *patch_flip = NULL;
    uint8_t *patch_accept = NULL, *chart_rank = NULL;
    int32_t *chart_parent = NULL;
    size_t *patch_faces = NULL, *patch_vertices = NULL, *patch_edges = NULL;
    size_t *patch_boundary = NULL, *attach_edges = NULL;
    size_t *attach_vertices = NULL, *attach_endpoints = NULL;
    size_t *patch_edge_begin = NULL, *patch_edge_end = NULL;
    size_t *patch_vertex_begin = NULL, *patch_vertex_end = NULL;
    double *patch_area = NULL;
    CbfEdgeFace *edge_face = NULL, *patch_edge = NULL;
    CbfInputEdge *input_edge = NULL;
    CbfPatchValue *patch_vertex = NULL;
    CbfAttachVertex *attach_vertex = NULL;
    CbfLinkEdge *link_edge = NULL;
    CbfCandidate *candidate = NULL, *growth = NULL;
    size_t ncandidate = 0, n_attach_vertex = 0;
    size_t ngrowth = 0;
    uint8_t *current_used = NULL;
    uint8_t *path_seen = NULL;
    int32_t *current_owner = NULL;
    int32_t *boundary_n0 = NULL, *boundary_n1 = NULL;
    int32_t *shared_a = NULL, *shared_b = NULL;
    int32_t *roots = NULL;
    uint32_t *intersection_degree = NULL;
    CbfEdgeMap current_edges;
    ChartBridgeForestStats st;
    int rc = -1;

    memset(&st, 0, sizeof st);
    memset(&current_edges, 0, sizeof current_edges);
    st.input_charts = n_charts;
    if (!arena || !verts || !faces || !vertex_chart || !out_faces || !out_nf ||
        input_nf > nf || nv == 0 || n_charts == 0 ||
        n_charts > (size_t)INT32_MAX)
        return -1;
    nbf = nf - input_nf;
    if(support_begin>support_end||support_end>nbf)return -1;
    st.bridge_faces_in = nbf;
    if (nbf == 0) {
        int32_t *copy = (int32_t *)ARENA_ALLOC(
            arena, ((nf ? nf : 1) * 3 * sizeof(int32_t)));
        memcpy(copy, faces, nf*3*sizeof(int32_t));
        *out_faces = copy; *out_nf = nf;
        if (stats) *stats = st;
        return 0;
    }
    if (nbf > (size_t)INT32_MAX || nbf > SIZE_MAX/3) return -1;

    face_parent = (int32_t *)malloc(nbf*sizeof(int32_t));
    face_rank = (uint8_t *)calloc(nbf, 1);
    face_valid = (uint8_t *)calloc(nbf, 1);
    face_chart_a = (int32_t *)malloc(nbf*sizeof(int32_t));
    face_chart_b = (int32_t *)malloc(nbf*sizeof(int32_t));
    face_chart_c = (int32_t *)malloc(nbf*sizeof(int32_t));
    face_chart_count = (uint8_t *)calloc(nbf, 1);
    edge_face = (CbfEdgeFace *)malloc(nbf*3*sizeof(CbfEdgeFace));
    if (!face_parent || !face_rank || !face_valid || !face_chart_a ||
        !face_chart_b || !face_chart_c || !face_chart_count || !edge_face)
        goto cleanup;
    for (size_t i = 0; i < nbf; i++) {
        const int32_t *tri = &faces[(input_nf+i)*3];
        int32_t chart[3]; int nc = 0;
        face_parent[i] = (int32_t)i;
        for (int k = 0; k < 3; k++) {
            int32_t v = tri[k], c;
            if (v < 0 || (size_t)v >= nv) { nc = 4; break; }
            c = vertex_chart[v];
            if (c < 0 || (size_t)c >= n_charts) { nc = 4; break; }
            int seen = 0;
            for (int j = 0; j < nc; j++) if (chart[j] == c) seen = 1;
            if (!seen && nc < 3) chart[nc++] = c;
        }
        if (nc >= 2 && chart[0] > chart[1]) {
            int32_t t=chart[0]; chart[0]=chart[1]; chart[1]=t;
        }
        if (nc == 3) {
            if (chart[1] > chart[2]) {
                int32_t t=chart[1]; chart[1]=chart[2]; chart[2]=t;
            }
            if (chart[0] > chart[1]) {
                int32_t t=chart[0]; chart[0]=chart[1]; chart[1]=t;
            }
        }
        face_chart_count[i]=(uint8_t)(nc >= 1 && nc <= 3 ? nc : 0);
        face_chart_c[i]=nc==3?chart[2]:-1;
        if (nc == 2) {
            face_chart_a[i]=chart[0]; face_chart_b[i]=chart[1];
            face_valid[i]=1; nvalid++;
        } else {
            face_chart_a[i]=nc>=1&&nc<=3?chart[0]:-1;
            face_chart_b[i]=nc>=2&&nc<=3?chart[1]:-1;
        }
        for (int k = 0; k < 3; k++) {
            int32_t a=tri[k], b=tri[(k+1)%3];
            edge_face[i*3+(size_t)k].lo = a < b ? a : b;
            edge_face[i*3+(size_t)k].hi = a < b ? b : a;
            edge_face[i*3+(size_t)k].face = (int32_t)i;
            edge_face[i*3+(size_t)k].dir = a < b ? 1 : -1;
        }
    }
    (void)nvalid;
    nedges = nbf*3;
    qsort(edge_face, nedges, sizeof(*edge_face), cbf_edge_face_cmp);
    if(promote_support||support_begin<support_end){
        /* Once earlier rounds have safely joined two source charts, some raw
         * pair faces collapse to a one-chart label.  They are not disposable:
         * an edge-connected run of those faces can be the geometric shoulder
         * that carries a neighbouring two-chart strip all the way onto the
         * coarsened chart boundary.  Recover only unambiguous shoulders.
         *
         * First form edge-connected components of one-chart faces (never
         * vertex-only).  A component is promoted only if every adjacent
         * two-chart seed requests the same pair and the component's chart is
         * one side of that pair.  Competing pair claims leave the whole
         * shoulder rejected. */
        int32_t *support_parent=(int32_t*)malloc(nbf*sizeof(*support_parent));
        uint8_t *support_rank=(uint8_t*)calloc(nbf,1);
        int32_t *support_a=(int32_t*)malloc(nbf*sizeof(*support_a));
        int32_t *support_b=(int32_t*)malloc(nbf*sizeof(*support_b));
        uint8_t *support_ambiguous=(uint8_t*)calloc(nbf,1);
        if(!support_parent||!support_rank||!support_a||!support_b||
           !support_ambiguous){
            free(support_ambiguous);free(support_b);free(support_a);
            free(support_rank);free(support_parent);goto cleanup;
        }
        for(size_t i=0;i<nbf;i++){
            support_parent[i]=(int32_t)i;
            support_a[i]=support_b[i]=-1;
        }
        for(size_t i=0;i<nedges;){
            size_t j=i+1;
            while(j<nedges&&edge_face[j].lo==edge_face[i].lo&&
                  edge_face[j].hi==edge_face[i].hi)j++;
            for(size_t a=i;a<j;a++)for(size_t b=a+1;b<j;b++){
                int32_t fa=edge_face[a].face,fb=edge_face[b].face;
                int eligible_a=promote_support||
                    ((size_t)fa>=support_begin&&(size_t)fa<support_end);
                int eligible_b=promote_support||
                    ((size_t)fb>=support_begin&&(size_t)fb<support_end);
                if(eligible_a&&eligible_b&&
                   face_chart_count[fa]==1&&face_chart_count[fb]==1&&
                   face_chart_a[fa]==face_chart_a[fb])
                    cbf_union(support_parent,support_rank,fa,fb);
            }
            i=j;
        }
        for(size_t i=0;i<nedges;){
            size_t j=i+1;
            while(j<nedges&&edge_face[j].lo==edge_face[i].lo&&
                  edge_face[j].hi==edge_face[i].hi)j++;
            for(size_t a=i;a<j;a++){
                int32_t self=edge_face[a].face;
                int eligible=promote_support||
                    ((size_t)self>=support_begin&&(size_t)self<support_end);
                if(!eligible||face_chart_count[self]!=1)continue;
                for(size_t q=i;q<j;q++){
                    int32_t seed=edge_face[q].face;
                    int32_t root,pa,pb;
                    if(seed==self||!face_valid[seed]||
                       face_chart_count[seed]!=2)continue;
                    if(face_chart_a[self]!=face_chart_a[seed]&&
                       face_chart_a[self]!=face_chart_b[seed])continue;
                    root=cbf_find(support_parent,self);
                    pa=face_chart_a[seed];pb=face_chart_b[seed];
                    if(pa>pb){int32_t t=pa;pa=pb;pb=t;}
                    if(support_a[root]<0){
                        support_a[root]=pa;support_b[root]=pb;
                    }else if(support_a[root]!=pa||support_b[root]!=pb){
                        support_ambiguous[root]=1;
                    }
                }
            }
            i=j;
        }
        for(size_t i=0;i<nbf;i++)if(face_chart_count[i]==1&&
           (promote_support||(i>=support_begin&&i<support_end))){
            int32_t root=cbf_find(support_parent,(int32_t)i);
            if(support_a[root]>=0&&!support_ambiguous[root]){
                face_chart_a[i]=support_a[root];
                face_chart_b[i]=support_b[root];
                face_valid[i]=1;
                nvalid++;
            }
        }
        free(support_ambiguous);free(support_b);free(support_a);
        free(support_rank);free(support_parent);
    }
    for (size_t i = 0; i < nedges; ) {
        size_t j=i+1;
        while (j<nedges && edge_face[j].lo==edge_face[i].lo &&
               edge_face[j].hi==edge_face[i].hi) j++;
        for (size_t a=i; a<j; a++) {
            int32_t fa=edge_face[a].face;
            if (!face_valid[fa]) continue;
            for (size_t b=a+1; b<j; b++) {
                int32_t fb=edge_face[b].face;
                if (face_valid[fb] &&
                    face_chart_a[fa]==face_chart_a[fb] &&
                    face_chart_b[fa]==face_chart_b[fb])
                    cbf_union(face_parent, face_rank, fa, fb);
            }
        }
        i=j;
    }
    for (size_t i=0;i<nbf;i++) face_parent[i]=cbf_find(face_parent,(int32_t)i);

    patch_faces=(size_t*)calloc(nbf,sizeof(size_t));
    patch_vertices=(size_t*)calloc(nbf,sizeof(size_t));
    patch_edges=(size_t*)calloc(nbf,sizeof(size_t));
    patch_boundary=(size_t*)calloc(nbf,sizeof(size_t));
    patch_area=(double*)calloc(nbf,sizeof(double));
    patch_valid=(uint8_t*)calloc(nbf,1);
    patch_topology=(uint8_t*)calloc(nbf,1);
    patch_reject=(uint8_t*)calloc(nbf,1);
    patch_link_bad=(uint8_t*)calloc(nbf,1);
    patch_flip=(uint8_t*)calloc(nbf,1);
    patch_accept=(uint8_t*)calloc(nbf,1);
    patch_vertex=(CbfPatchValue*)malloc(nbf*3*sizeof(CbfPatchValue));
    patch_edge=(CbfEdgeFace*)malloc(nbf*3*sizeof(CbfEdgeFace));
    patch_edge_begin=(size_t*)malloc(nbf*sizeof(size_t));
    patch_edge_end=(size_t*)malloc(nbf*sizeof(size_t));
    patch_vertex_begin=(size_t*)malloc(nbf*sizeof(size_t));
    patch_vertex_end=(size_t*)malloc(nbf*sizeof(size_t));
    if(!patch_faces||!patch_vertices||!patch_edges||!patch_boundary||
       !patch_area||!patch_valid||!patch_topology||!patch_reject||
       !patch_link_bad||!patch_flip||
       !patch_accept||!patch_vertex||!patch_edge||!patch_edge_begin||
       !patch_edge_end||!patch_vertex_begin||!patch_vertex_end)
        goto cleanup;
    for(size_t i=0;i<nbf;i++){
        patch_edge_begin[i]=patch_edge_end[i]=SIZE_MAX;
        patch_vertex_begin[i]=patch_vertex_end[i]=SIZE_MAX;
    }
    for(size_t i=0;i<nbf;i++) if(face_valid[i]){
        int32_t r=face_parent[i];
        patch_faces[r]++;
        patch_area[r]+=cbf_face_area(verts,&faces[(input_nf+i)*3]);
        for(int k=0;k<3;k++){
            patch_vertex[i*3+(size_t)k].patch=r;
            patch_vertex[i*3+(size_t)k].value=faces[(input_nf+i)*3+(size_t)k];
            patch_edge[i*3+(size_t)k]=edge_face[0]; /* fields overwritten below */
            int32_t a=faces[(input_nf+i)*3+(size_t)k];
            int32_t b=faces[(input_nf+i)*3+(size_t)((k+1)%3)];
            patch_edge[i*3+(size_t)k].lo=a<b?a:b;
            patch_edge[i*3+(size_t)k].hi=a<b?b:a;
            patch_edge[i*3+(size_t)k].face=r; /* patch id */
            patch_edge[i*3+(size_t)k].dir=a<b?1:-1;
        }
    } else {
        for(int k=0;k<3;k++){
            patch_vertex[i*3+(size_t)k].patch=-1;
            patch_vertex[i*3+(size_t)k].value=-1;
            patch_edge[i*3+(size_t)k].face=-1;
            patch_edge[i*3+(size_t)k].lo=patch_edge[i*3+(size_t)k].hi=-1;
            patch_edge[i*3+(size_t)k].dir=0;
        }
    }
    qsort(patch_vertex,nbf*3,sizeof(*patch_vertex),cbf_patch_value_cmp);
    for(size_t i=0;i<nbf*3;){
        size_t j=i+1;
        while(j<nbf*3&&patch_vertex[j].patch==patch_vertex[i].patch&&
              patch_vertex[j].value==patch_vertex[i].value)j++;
        if(patch_vertex[i].patch>=0){
            int32_t r=patch_vertex[i].patch;
            patch_vertices[r]++;
            if(patch_vertex_begin[r]==SIZE_MAX)patch_vertex_begin[r]=i;
            patch_vertex_end[r]=j;
        }
        i=j;
    }
    qsort(patch_edge,nbf*3,sizeof(*patch_edge),cbf_patch_edge_cmp);

    /* A bridge patch must be vertex-manifold in its own right.  Edge counts
     * and chi=1 do not detect a bow-tie vertex: the link of every patch vertex
     * must be one connected path (boundary vertex) or cycle (interior vertex).
     * Build one link edge between the two opposite triangle vertices for every
     * face corner, then validate each (patch,center) link graph. */
    link_edge=(CbfLinkEdge*)malloc(nbf*3*sizeof(*link_edge));
    if(!link_edge)goto cleanup;
    {
        size_t nl=0;
        for(size_t i=0;i<nbf;i++)if(face_valid[i]){
            const int32_t *tri=&faces[(input_nf+i)*3];
            int32_t r=face_parent[i];
            for(int k=0;k<3;k++){
                int32_t a=tri[(k+1)%3],b=tri[(k+2)%3];
                if(a>b){int32_t t=a;a=b;b=t;}
                link_edge[nl++]=(CbfLinkEdge){r,tri[k],a,b};
            }
        }
        qsort(link_edge,nl,sizeof(*link_edge),cbf_link_edge_cmp);
        for(size_t i=0;i<nl;){
            size_t j=i+1,m,nnode=0,endpoints=0,nroot=0;
            int bad=0;
            while(j<nl&&link_edge[j].patch==link_edge[i].patch&&
                  link_edge[j].center==link_edge[i].center)j++;
            m=j-i;
            int32_t *node=(int32_t*)malloc(2*m*sizeof(int32_t));
            int32_t *lparent=(int32_t*)malloc(2*m*sizeof(int32_t));
            uint8_t *lrank=(uint8_t*)calloc(2*m,1);
            uint32_t *degree=(uint32_t*)calloc(2*m,sizeof(uint32_t));
            if(!node||!lparent||!lrank||!degree){
                free(degree);free(lrank);free(lparent);free(node);
                goto cleanup;
            }
            for(size_t q=i;q<j;q++){
                int32_t value[2]={link_edge[q].a,link_edge[q].b};
                size_t idx[2];
                if(value[0]==value[1]){bad=1;break;}
                for(int side=0;side<2;side++){
                    size_t p=0;
                    while(p<nnode&&node[p]!=value[side])p++;
                    if(p==nnode){
                        node[nnode]=value[side];
                        lparent[nnode]=(int32_t)nnode;
                        nnode++;
                    }
                    idx[side]=p;
                    degree[p]++;
                    if(degree[p]>2)bad=1;
                }
                cbf_union(lparent,lrank,(int32_t)idx[0],(int32_t)idx[1]);
            }
            if(!bad){
                for(size_t p=0;p<nnode;p++){
                    if(degree[p]==1)endpoints++;
                    if(cbf_find(lparent,(int32_t)p)==(int32_t)p)nroot++;
                }
                if(nroot!=1||!(endpoints==0||endpoints==2))bad=1;
            }
            if(bad)patch_link_bad[link_edge[i].patch]=1;
            free(degree);free(lrank);free(lparent);free(node);
            i=j;
        }
    }

    /* Input edge multiplicity tells which patch edges truly attach to an
     * original chart boundary. */
    ninput_edges=input_nf*3;
    input_edge=(CbfInputEdge*)malloc((ninput_edges?ninput_edges:1)*sizeof(*input_edge));
    if(!input_edge)goto cleanup;
    for(size_t f=0;f<input_nf;f++)for(int k=0;k<3;k++){
        int32_t a=faces[f*3+(size_t)k],b=faces[f*3+(size_t)((k+1)%3)];
        input_edge[f*3+(size_t)k].lo=a<b?a:b;
        input_edge[f*3+(size_t)k].hi=a<b?b:a;
        input_edge[f*3+(size_t)k].count=1;
    }
    qsort(input_edge,ninput_edges,sizeof(*input_edge),cbf_input_edge_cmp);
    {
        size_t w=0;
        for(size_t i=0;i<ninput_edges;){
            size_t j=i+1;
            while(j<ninput_edges&&input_edge[j].lo==input_edge[i].lo&&
                  input_edge[j].hi==input_edge[i].hi)j++;
            input_edge[w]=input_edge[i]; input_edge[w].count=(uint32_t)(j-i); w++;
            i=j;
        }
        ninput_edges=w;
    }
    attach_edges=(size_t*)calloc(nbf*2,sizeof(size_t));
    attach_vertices=(size_t*)calloc(nbf*2,sizeof(size_t));
    attach_endpoints=(size_t*)calloc(nbf*2,sizeof(size_t));
    attach_vertex=(CbfAttachVertex*)malloc(nbf*6*sizeof(*attach_vertex));
    if(!attach_edges||!attach_vertices||!attach_endpoints||!attach_vertex)
        goto cleanup;
    for(size_t i=0;i<nbf*3;){
        size_t j=i+1;
        while(j<nbf*3&&patch_edge[j].face==patch_edge[i].face&&
              patch_edge[j].lo==patch_edge[i].lo&&patch_edge[j].hi==patch_edge[i].hi)j++;
        if(patch_edge[i].face>=0){
            int32_t r=patch_edge[i].face;
            size_t count=j-i;
            int dirsum=0;
            for(size_t q=i;q<j;q++)dirsum+=patch_edge[q].dir;
            if(patch_edge_begin[r]==SIZE_MAX)patch_edge_begin[r]=i;
            patch_edge_end[r]=j;
            patch_edges[r]++;
            if(count==1)patch_boundary[r]++;
            if(count>2||(count==2&&dirsum!=0))
                patch_valid[r]=2; /* invalid/nonorientable sentinel */
            const CbfInputEdge *ie=cbf_input_edge_find(input_edge,ninput_edges,
                                                       patch_edge[i].lo,
                                                       patch_edge[i].hi);
            if(ie&&ie->count==1){
                int32_t c0=vertex_chart[patch_edge[i].lo];
                int32_t c1=vertex_chart[patch_edge[i].hi];
                if(c0==c1&&(c0==face_chart_a[r]||c0==face_chart_b[r])){
                    int side=c0==face_chart_a[r]?0:1;
                    attach_edges[(size_t)r*2+(size_t)side]++;
                    attach_vertex[n_attach_vertex++] =
                        (CbfAttachVertex){r,side,patch_edge[i].lo};
                    attach_vertex[n_attach_vertex++] =
                        (CbfAttachVertex){r,side,patch_edge[i].hi};
                }
            }
        }
        i=j;
    }
    qsort(attach_vertex,n_attach_vertex,sizeof(*attach_vertex),
          cbf_attach_vertex_cmp);
    for(size_t i=0;i<n_attach_vertex;){
        size_t j=i+1;
        while(j<n_attach_vertex&&attach_vertex[j].patch==attach_vertex[i].patch&&
              attach_vertex[j].side==attach_vertex[i].side&&
              attach_vertex[j].vertex==attach_vertex[i].vertex)j++;
        size_t idx=(size_t)attach_vertex[i].patch*2+(size_t)attach_vertex[i].side;
        attach_vertices[idx]++;
        if(j-i==1)attach_endpoints[idx]++;
        else if(j-i!=2)patch_valid[attach_vertex[i].patch]=2;
        i=j;
    }

    candidate=(CbfCandidate*)malloc(nbf*sizeof(*candidate));
    chart_parent=(int32_t*)malloc(n_charts*sizeof(int32_t));
    chart_rank=(uint8_t*)calloc(n_charts,1);
    if(!candidate||!chart_parent||!chart_rank)goto cleanup;
    for(size_t c=0;c<n_charts;c++)chart_parent[c]=(int32_t)c;
    for(size_t r=0;r<nbf;r++)if(face_valid[r]&&face_parent[r]==(int32_t)r){
        int topology_ok=patch_valid[r]!=2&&!patch_link_bad[r]&&
            patch_boundary[r]>0&&
            (int64_t)patch_vertices[r]-(int64_t)patch_edges[r]+
            (int64_t)patch_faces[r]==1;
        int attach_ok=1;
        st.patches++;
        for(int side=0;side<2;side++){
            size_t idx=r*2+(size_t)side;
            if(attach_edges[idx]==0||
               attach_vertices[idx]!=attach_edges[idx]+1||
               attach_endpoints[idx]!=2)attach_ok=0;
        }
        if(!topology_ok){st.rejected_patch_topology++;continue;}
        patch_topology[r]=1;
        if(!attach_ok)continue;
        patch_valid[r]=1;
        candidate[ncandidate++]=(CbfCandidate){(int32_t)r,
            face_chart_a[r],face_chart_b[r],patch_faces[r],patch_area[r]};
    }
    st.rejected_self_or_mixed=nbf-nvalid;
    qsort(candidate,ncandidate,sizeof(*candidate),cbf_candidate_cmp);
    /* Do not pre-commit even the apparently ideal two-sided joins.  Two such
     * patches can be individually valid yet collide at one already-used
     * boundary vertex.  Every patch, including the first chart join, is
     * admitted below against the live surface edge/link state. */

    /* Grow the chart forest maximally.  A disk patch can be appended without
     * changing disk topology iff its intersection with the current surface is
     * one boundary path, or two boundary paths belonging to two distinct disk
     * components.  The latter is another safe join; the former expands seam
     * coverage around an already-accepted core.  Isolated point contacts and
     * multiple paths on one component are precisely the pinch/handle cases. */
    current_used=(uint8_t*)calloc(nv,1);
    path_seen=(uint8_t*)calloc(nv,1);
    current_owner=(int32_t*)malloc(nv*sizeof(int32_t));
    boundary_n0=(int32_t*)malloc(nv*sizeof(int32_t));
    boundary_n1=(int32_t*)malloc(nv*sizeof(int32_t));
    shared_a=(int32_t*)malloc(nbf*3*sizeof(int32_t));
    shared_b=(int32_t*)malloc(nbf*3*sizeof(int32_t));
    intersection_degree=(uint32_t*)calloc(nv,sizeof(uint32_t));
    growth=(CbfCandidate*)malloc(nbf*sizeof(*growth));
    roots=(int32_t*)malloc(n_charts*sizeof(*roots));
    if(!current_used||!path_seen||!current_owner||!boundary_n0||
       !boundary_n1||!shared_a||!shared_b||!intersection_degree||!growth||
       !roots)
        goto cleanup;
    for(size_t v=0;v<nv;v++){
        current_owner[v]=-1;
        boundary_n0[v]=boundary_n1[v]=-1;
    }
    if(cbf_edge_map_init(&current_edges,(input_nf+nbf)*3)!=0)goto cleanup;
    for(size_t f=0;f<input_nf;f++){
        for(int k=0;k<3;k++){
            int32_t v=faces[f*3+(size_t)k];
            int32_t a=faces[f*3+(size_t)k];
            int32_t b=faces[f*3+(size_t)((k+1)%3)];
            int32_t lo=a<b?a:b,hi=a<b?b:a;
            if(vertex_chart[v]<0||(size_t)vertex_chart[v]>=n_charts||
               cbf_edge_map_add(&current_edges,lo,hi,1,a<b?1:-1)!=0)
                goto cleanup;
            current_used[v]=1;
            current_owner[v]=vertex_chart[v];
        }
    }
    for(size_t s=0;s<current_edges.cap;s++)if(current_edges.slot[s].used&&
       current_edges.slot[s].count==1){
        int32_t a=current_edges.slot[s].lo,b=current_edges.slot[s].hi;
        if(cbf_boundary_add(boundary_n0,boundary_n1,a,b)!=0||
           cbf_boundary_add(boundary_n0,boundary_n1,b,a)!=0)goto cleanup;
    }
    for(size_t i=0;i<nbf;i++)if(face_valid[i]&&patch_accept[face_parent[i]]){
        int32_t r=face_parent[i];
        for(int k=0;k<3;k++){
            int32_t v=faces[(input_nf+i)*3+(size_t)k];
            int32_t a=faces[(input_nf+i)*3+(size_t)k];
            int32_t b=faces[(input_nf+i)*3+(size_t)((k+1)%3)];
            int32_t lo=a<b?a:b,hi=a<b?b:a;
            if(cbf_edge_map_add(&current_edges,lo,hi,1,a<b?1:-1)!=0)
                goto cleanup;
            current_used[v]=1;
            if(current_owner[v]<0)current_owner[v]=face_chart_a[r];
        }
    }
    for(size_t r=0;r<nbf;r++)if(face_valid[r]&&
       face_parent[r]==(int32_t)r&&patch_topology[r])
        growth[ngrowth++]=(CbfCandidate){(int32_t)r,face_chart_a[r],
            face_chart_b[r],patch_faces[r],patch_area[r]};
    qsort(growth,ngrowth,sizeof(*growth),cbf_candidate_cmp);
    /* Admission is lexicographic: first maximize genuine cross-component chart
     * joins to a fixed point, only then spend boundary arcs on redundant
     * same-component coverage strips.  Mixing both in one area-sorted pass let
     * a large annular strip consume an arc before a small neighbouring chart's
     * otherwise-safe join was visited, stranding visible micro-charts. */
    for(int admission_phase=0;admission_phase<2;admission_phase++){
      for(;;){
        size_t accepted_this_round=0;
        for(size_t ci=0;ci<ngrowth;ci++){
            int32_t r=growth[ci].root;
            size_t shared_vertices=0,shared_edges=0,endpoints=0;
            size_t paths=0;
            size_t nshared=0;
            int desired_flip=0;
            size_t nroots=0;
            int reason=CBF_REJECT_NONE;
            if(patch_accept[r])continue;
            for(size_t i=patch_vertex_begin[r];i<patch_vertex_end[r];){
                size_t j=i+1;
                while(j<patch_vertex_end[r]&&
                      patch_vertex[j].value==patch_vertex[i].value)j++;
                intersection_degree[patch_vertex[i].value]=0;
                path_seen[patch_vertex[i].value]=0;
                i=j;
            }
            for(size_t i=patch_edge_begin[r];i<patch_edge_end[r];){
                size_t j=i+1;
                while(j<patch_edge_end[r]&&patch_edge[j].lo==patch_edge[i].lo&&
                      patch_edge[j].hi==patch_edge[i].hi)j++;
                uint32_t pc=(uint32_t)(j-i);
                uint32_t cc=cbf_edge_map_get(&current_edges,
                                             patch_edge[i].lo,
                                             patch_edge[i].hi);
                int32_t pdir=0;
                for(size_t q=i;q<j;q++)pdir+=patch_edge[q].dir;
                if(pc>2||cc>2||pc+cc>2){reason=CBF_REJECT_EDGE;break;}
                if(cc){
                    int32_t a=patch_edge[i].lo,b=patch_edge[i].hi;
                    CbfEdgeSlot *slot=cbf_edge_map_slot(
                        &current_edges,patch_edge[i].lo,patch_edge[i].hi);
                    int edge_flip;
                    if(!current_used[a]||!current_used[b]||
                       cbf_find(chart_parent,current_owner[a])!=
                       cbf_find(chart_parent,current_owner[b])){
                        reason=CBF_REJECT_EDGE;break;
                    }
                    if(pc!=1||cc!=1||
                       !(pdir==1||pdir==-1)||
                       !(slot->orient_sum==1||slot->orient_sum==-1)){
                        reason=CBF_REJECT_ORIENTATION;break;
                    }
                    edge_flip=(pdir==-slot->orient_sum)?1:-1;
                    if(desired_flip==0)desired_flip=edge_flip;
                    else if(desired_flip!=edge_flip){
                        reason=CBF_REJECT_ORIENTATION;break;
                    }
                    intersection_degree[a]++;
                    intersection_degree[b]++;
                    shared_a[nshared]=a;
                    shared_b[nshared]=b;
                    nshared++;
                    shared_edges++;
                }
                i=j;
            }
            if(reason==CBF_REJECT_NONE){
                for(size_t i=patch_vertex_begin[r];i<patch_vertex_end[r];){
                    size_t j=i+1;
                    int32_t v=patch_vertex[i].value;
                    while(j<patch_vertex_end[r]&&patch_vertex[j].value==v)j++;
                    if(current_used[v]){
                        uint32_t d=intersection_degree[v];
                        int32_t root;
                        if(d==0||d>2||current_owner[v]<0){
                            if(getenv("SEAM_CHART_DIAGNOSTICS")!=NULL)
                                fprintf(stderr,
                                    "  [chart-attachment] patch=%d vertex=%d "
                                    "xyz=(%.3f %.3f %.3f) used=%u degree=%u "
                                    "owner=%d reason=%s\n",
                                    r,v,(double)verts[(size_t)v*3],
                                    (double)verts[(size_t)v*3+1],
                                    (double)verts[(size_t)v*3+2],
                                    (unsigned)current_used[v],(unsigned)d,
                                    current_owner[v],
                                    d==0?"point-contact":
                                    (d>2?"branch":"unowned"));
                            reason=CBF_REJECT_ATTACHMENT;break;
                        }
                        shared_vertices++;
                        if(d==1)endpoints++;
                        root=cbf_find(chart_parent,current_owner[v]);
                        size_t q=0;
                        while(q<nroots&&roots[q]!=root)q++;
                        if(q==nroots)roots[nroots++]=root;
                    }else if(intersection_degree[v]!=0){
                        reason=CBF_REJECT_EDGE;break;
                    }
                    i=j;
                }
            }
            if(reason==CBF_REJECT_NONE){
                if(shared_vertices<=shared_edges){
                    reason=CBF_REJECT_MULTIPATH;
                }else{
                    paths=shared_vertices-shared_edges;
                    if(paths<1||endpoints!=2*paths||
                       !((nroots==1&&paths<=2)||
                          (nroots>=2&&paths==nroots)))
                        reason=CBF_REJECT_MULTIPATH;
                }
            }
            if(reason==CBF_REJECT_NONE&&paths==2&&nroots==1){
                int32_t first=-1,second=-1;
                if(nshared==0){
                    reason=CBF_REJECT_MULTIPATH;
                }else{
                    first=shared_a[0];
                    path_seen[first]=1;
                    for(;;){
                        int changed=0;
                        for(size_t e=0;e<nshared;e++){
                            if(path_seen[shared_a[e]]&&!path_seen[shared_b[e]]){
                                path_seen[shared_b[e]]=1;changed=1;
                            }else if(path_seen[shared_b[e]]&&
                                    !path_seen[shared_a[e]]){
                                path_seen[shared_a[e]]=1;changed=1;
                            }
                        }
                        if(!changed)break;
                    }
                    for(size_t i=patch_vertex_begin[r];
                        i<patch_vertex_end[r];){
                        size_t j=i+1;
                        int32_t v=patch_vertex[i].value;
                        while(j<patch_vertex_end[r]&&
                              patch_vertex[j].value==v)j++;
                        if(current_used[v]&&intersection_degree[v]>0&&
                           !path_seen[v]){second=v;break;}
                        i=j;
                    }
                    if(second<0||!cbf_boundary_same_loop(
                            boundary_n0,boundary_n1,nv,first,second))
                        reason=CBF_REJECT_MULTIPATH;
                }
            }
            patch_reject[r]=(uint8_t)reason;
            if(reason!=CBF_REJECT_NONE)continue;
            if(admission_phase==0&&nroots<2)continue;
            if(desired_flip==0){
                patch_reject[r]=CBF_REJECT_ORIENTATION;
                continue;
            }
            if(nroots>=2){
                for(size_t q=1;q<nroots;q++)
                    cbf_union(chart_parent,chart_rank,roots[0],roots[q]);
                st.accepted_join_patches++;
            }else if(paths==2){
                /* A second strip on one genus-zero component creates another
                 * ordinary boundary loop.  The chart-level hole-fill stage
                 * caps that loop before any vertex-level cleanup. */
                st.accepted_hole_patches++;
            }else{
                st.accepted_growth_patches++;
            }
            {
                int32_t owner=cbf_find(chart_parent,roots[0]);
                patch_flip[r]=(uint8_t)(desired_flip<0);
                /* Remove attachment paths before adding the replacement patch
                 * rim, so a boundary vertex never transiently has degree 3. */
                for(size_t i=patch_edge_begin[r];i<patch_edge_end[r];){
                    size_t j=i+1;
                    while(j<patch_edge_end[r]&&
                          patch_edge[j].lo==patch_edge[i].lo&&
                          patch_edge[j].hi==patch_edge[i].hi)j++;
                    if((uint32_t)(j-i)==1&&cbf_edge_map_get(
                            &current_edges,patch_edge[i].lo,
                            patch_edge[i].hi)==1){
                        if(cbf_boundary_remove(boundary_n0,boundary_n1,
                                               patch_edge[i].lo,
                                               patch_edge[i].hi)!=0||
                           cbf_boundary_remove(boundary_n0,boundary_n1,
                                               patch_edge[i].hi,
                                               patch_edge[i].lo)!=0)
                            goto cleanup;
                    }
                    i=j;
                }
                for(size_t i=patch_edge_begin[r];i<patch_edge_end[r];){
                    size_t j=i+1;
                    int32_t add_dir=0;
                    while(j<patch_edge_end[r]&&
                          patch_edge[j].lo==patch_edge[i].lo&&
                          patch_edge[j].hi==patch_edge[i].hi)j++;
                    for(size_t q=i;q<j;q++)add_dir+=patch_edge[q].dir;
                    {
                        uint32_t pc=(uint32_t)(j-i);
                        uint32_t cc=cbf_edge_map_get(
                            &current_edges,patch_edge[i].lo,patch_edge[i].hi);
                        if(pc==1&&cc==0){
                            if(cbf_boundary_add(boundary_n0,boundary_n1,
                                                patch_edge[i].lo,
                                                patch_edge[i].hi)!=0||
                               cbf_boundary_add(boundary_n0,boundary_n1,
                                                patch_edge[i].hi,
                                                patch_edge[i].lo)!=0)
                                goto cleanup;
                        }
                    }
                    if(cbf_edge_map_add(&current_edges,patch_edge[i].lo,
                                        patch_edge[i].hi,
                                        (uint32_t)(j-i),
                                        add_dir*desired_flip)!=0)goto cleanup;
                    i=j;
                }
                for(size_t i=patch_vertex_begin[r];i<patch_vertex_end[r];){
                    size_t j=i+1;
                    int32_t v=patch_vertex[i].value;
                    while(j<patch_vertex_end[r]&&patch_vertex[j].value==v)j++;
                    if(!current_used[v]){
                        current_used[v]=1;
                        current_owner[v]=owner;
                    }
                    i=j;
                }
            }
            patch_accept[r]=1;
            patch_reject[r]=CBF_REJECT_NONE;
            st.accepted_patches++;
            st.bridge_faces_kept+=patch_faces[r];
            accepted_this_round++;
        }
        if(accepted_this_round==0)break;
      }
    }
    st.rejected_attachment=0;
    st.rejected_cycle=0;
    st.rejected_edge_conflict=0;
    st.rejected_orientation=0;
    for(size_t r=0;r<nbf;r++)if(face_valid[r]&&
       face_parent[r]==(int32_t)r&&patch_topology[r]&&!patch_accept[r]){
        if(patch_reject[r]==CBF_REJECT_EDGE)st.rejected_edge_conflict++;
        else if(patch_reject[r]==CBF_REJECT_ORIENTATION)
            st.rejected_orientation++;
        else if(patch_reject[r]==CBF_REJECT_MULTIPATH)st.rejected_cycle++;
        else st.rejected_attachment++;
    }
    if(getenv("SEAM_CHART_DIAGNOSTICS")!=NULL){
        CbfGroupDiag *by_root=(CbfGroupDiag*)calloc(n_charts,sizeof(*by_root));
        CbfGroupDiag *ordered=(CbfGroupDiag*)malloc(n_charts*sizeof(*ordered));
        int32_t *root_rank=(int32_t*)malloc(n_charts*sizeof(*root_rank));
        size_t ngroup=0,missed=0,mixed_cross=0,self_faces=0,three_faces=0;
        if(!by_root||!ordered||!root_rank){
            free(root_rank);free(ordered);free(by_root);goto cleanup;
        }
        for(size_t c=0;c<n_charts;c++){
            int32_t r=cbf_find(chart_parent,(int32_t)c);
            by_root[r].root=r;by_root[r].charts++;
            root_rank[c]=-1;
        }
        for(size_t f=0;f<input_nf;f++){
            const int32_t *tri=&faces[f*3];
            int32_t chart=vertex_chart[tri[0]];
            int32_t r=cbf_find(chart_parent,chart);
            CbfGroupDiag *g=&by_root[r];
            int first=g->input_faces==0;
            g->input_faces++;
            for(int k=0;k<3;k++)for(int d=0;d<3;d++){
                double x=verts[(size_t)tri[k]*3+(size_t)d];
                if(first&&k==0){g->bbox_min[d]=g->bbox_max[d]=x;}
                else{
                    if(x<g->bbox_min[d])g->bbox_min[d]=x;
                    if(x>g->bbox_max[d])g->bbox_max[d]=x;
                }
            }
        }
        for(size_t i=0;i<nbf;i++)if(face_valid[i]&&
           patch_accept[face_parent[i]]){
            int32_t r=cbf_find(chart_parent,face_chart_a[i]);
            by_root[r].bridge_faces++;
        }
        for(size_t c=0;c<n_charts;c++)if(by_root[c].charts)
            ordered[ngroup++]=by_root[c];
        qsort(ordered,ngroup,sizeof(*ordered),cbf_group_diag_cmp);
        for(size_t g=0;g<ngroup;g++)root_rank[ordered[g].root]=(int32_t)g+1;
        fprintf(stderr,"  [chart-diagnostic] %zu round chart group(s):\n",ngroup);
        for(size_t g=0;g<ngroup;g++){
            const CbfGroupDiag *q=&ordered[g];
            fprintf(stderr,"    group %zu root=%d charts=%zu faces=%zu+%zu "
                    "bbox=[%.2f %.2f %.2f]-[%.2f %.2f %.2f]\n",
                    g+1,q->root,q->charts,q->input_faces,q->bridge_faces,
                    q->bbox_min[0],q->bbox_min[1],q->bbox_min[2],
                    q->bbox_max[0],q->bbox_max[1],q->bbox_max[2]);
        }
        fprintf(stderr,"  [chart-diagnostic] rejected cross-group pair patches:\n");
        for(size_t r=0;r<nbf;r++)if(face_valid[r]&&
           face_parent[r]==(int32_t)r&&patch_topology[r]&&!patch_accept[r]){
            int32_t ra=cbf_find(chart_parent,face_chart_a[r]);
            int32_t rb=cbf_find(chart_parent,face_chart_b[r]);
            if(ra!=rb){
                double lo[3]={DBL_MAX,DBL_MAX,DBL_MAX};
                double hi[3]={-DBL_MAX,-DBL_MAX,-DBL_MAX};
                for(size_t i=patch_vertex_begin[r];i<patch_vertex_end[r];){
                    size_t j=i+1;int32_t v=patch_vertex[i].value;
                    while(j<patch_vertex_end[r]&&patch_vertex[j].value==v)j++;
                    for(int d=0;d<3;d++){
                        double x=verts[(size_t)v*3+(size_t)d];
                        if(x<lo[d])lo[d]=x;if(x>hi[d])hi[d]=x;
                    }
                    i=j;
                }
                fprintf(stderr,"    patch=%zu groups=%d/%d charts=%d/%d "
                        "faces=%zu area=%.3f reject=%s bbox="
                        "[%.2f %.2f %.2f]-[%.2f %.2f %.2f]\n",
                        r,root_rank[ra],root_rank[rb],face_chart_a[r],
                        face_chart_b[r],patch_faces[r],patch_area[r],
                        cbf_reject_name(patch_reject[r]),
                        lo[0],lo[1],lo[2],hi[0],hi[1],hi[2]);
                if(getenv("SEAM_CHART_DIAGNOSTICS_VERBOSE")!=NULL){
                    for(size_t i=0;i<nbf;i++)if(face_valid[i]&&
                       face_parent[i]==(int32_t)r){
                        const int32_t *tri=&faces[(input_nf+i)*3];
                        fprintf(stderr,
                            "      raw_face=%zu tri=%d/%d/%d labels=%d/%d/%d "
                            "xyz0=(%.3f %.3f %.3f) xyz1=(%.3f %.3f %.3f) "
                            "xyz2=(%.3f %.3f %.3f)\n",
                            i,tri[0],tri[1],tri[2],
                            vertex_chart[tri[0]],vertex_chart[tri[1]],
                            vertex_chart[tri[2]],
                            (double)verts[(size_t)tri[0]*3],
                            (double)verts[(size_t)tri[0]*3+1],
                            (double)verts[(size_t)tri[0]*3+2],
                            (double)verts[(size_t)tri[1]*3],
                            (double)verts[(size_t)tri[1]*3+1],
                            (double)verts[(size_t)tri[1]*3+2],
                            (double)verts[(size_t)tri[2]*3],
                            (double)verts[(size_t)tri[2]*3+1],
                            (double)verts[(size_t)tri[2]*3+2]);
                    }
                }
                missed++;
            }
        }
        fprintf(stderr,"  [chart-diagnostic] mixed raw bridge faces spanning "
                "final groups:\n");
        for(size_t i=0;i<nbf;i++){
            if(face_chart_count[i]==1)self_faces++;
            if(face_chart_count[i]==3){
                int32_t diag_roots[3];size_t nr=0;
                int32_t chart[3]={face_chart_a[i],face_chart_b[i],face_chart_c[i]};
                three_faces++;
                for(int k=0;k<3;k++){
                    int32_t r=cbf_find(chart_parent,chart[k]);size_t q=0;
                    while(q<nr&&diag_roots[q]!=r)q++;
                    if(q==nr)diag_roots[nr++]=r;
                }
                if(nr>1){
                    const int32_t *tri=&faces[(input_nf+i)*3];
                    double lo[3]={DBL_MAX,DBL_MAX,DBL_MAX};
                    double hi[3]={-DBL_MAX,-DBL_MAX,-DBL_MAX};
                    for(int k=0;k<3;k++)for(int d=0;d<3;d++){
                        double x=verts[(size_t)tri[k]*3+(size_t)d];
                        if(x<lo[d])lo[d]=x;if(x>hi[d])hi[d]=x;
                    }
                    fprintf(stderr,"    face=%zu charts=%d/%d/%d groups="
                            "%d/%d/%d bbox=[%.2f %.2f %.2f]-"
                            "[%.2f %.2f %.2f]\n",i,chart[0],chart[1],chart[2],
                             root_rank[diag_roots[0]],root_rank[diag_roots[1]],
                             nr>2?root_rank[diag_roots[2]]:
                                  root_rank[diag_roots[1]],
                            lo[0],lo[1],lo[2],hi[0],hi[1],hi[2]);
                    mixed_cross++;
                }
            }
        }
        fprintf(stderr,"  [chart-diagnostic] summary: missed_pair=%zu "
                "self_faces=%zu three_chart_faces=%zu mixed_cross=%zu\n",
                missed,self_faces,three_faces,mixed_cross);
        free(root_rank);free(ordered);free(by_root);
    }
    {
        size_t total=input_nf+st.bridge_faces_kept,w=input_nf;
        int32_t *out=(int32_t*)ARENA_ALLOC(arena,((total?total:1)*3*sizeof(int32_t)));
        memcpy(out,faces,input_nf*3*sizeof(int32_t));
        for(size_t i=0;i<nbf;i++)if(face_valid[i]&&patch_accept[face_parent[i]]){
            const int32_t *src=&faces[(input_nf+i)*3];
            if(out_bridge_kept)out_bridge_kept[i]=1;
            out[w*3]=src[0];
            if(patch_flip[face_parent[i]]){
                out[w*3+1]=src[2];out[w*3+2]=src[1];
            }else{
                out[w*3+1]=src[1];out[w*3+2]=src[2];
            }
            w++;
        }
        *out_faces=out;*out_nf=w;
    }
    if(stats)*stats=st;
    rc=0;

cleanup:
    free(current_edges.slot);
    free(roots);free(shared_b);free(shared_a);free(boundary_n1);free(boundary_n0);
    free(intersection_degree);free(current_owner);free(path_seen);
    free(current_used);free(growth);
    free(chart_rank);free(chart_parent);free(candidate);
    free(link_edge);free(attach_vertex);free(attach_endpoints);
    free(attach_vertices);free(attach_edges);
    free(input_edge);free(patch_vertex);free(patch_edge);
    free(patch_vertex_end);free(patch_vertex_begin);
    free(patch_edge_end);free(patch_edge_begin);
    free(patch_area);free(patch_accept);free(patch_reject);
    free(patch_flip);free(patch_link_bad);free(patch_topology);free(patch_valid);
    free(patch_boundary);
    free(patch_edges);free(patch_vertices);free(patch_faces);
    free(edge_face);free(face_chart_count);free(face_chart_c);
    free(face_chart_b);free(face_chart_a);free(face_valid);
    free(face_rank);free(face_parent);
    return rc;
}

int ChartBridgeForest_filter_with_support(
                             Arena_T arena,
                             const float *verts, size_t nv,
                             const int32_t *faces, size_t nf,
                             size_t input_nf,
                             const int32_t *vertex_chart,
                             size_t n_charts,
                             size_t support_bridge_begin,
                             size_t support_bridge_end,
                             int32_t **out_faces, size_t *out_nf,
                             ChartBridgeForestStats *stats)
{
    const size_t original_bridge_nf = nf >= input_nf ? nf-input_nf : 0;
    int32_t *current_faces = NULL;
    size_t current_nf = input_nf;
    int32_t *remaining = NULL;
    size_t remaining_nf = original_bridge_nf;
    int32_t *label = NULL;
    size_t current_charts = n_charts;
    ChartBridgeForestStats total;
    int rc = -1;

    memset(&total,0,sizeof total);
    total.input_charts=n_charts;
    total.bridge_faces_in=original_bridge_nf;
    if(!arena||!verts||!faces||!vertex_chart||!out_faces||!out_nf||
       input_nf>nf||nv==0||n_charts==0||
       n_charts>(size_t)INT32_MAX||
       support_bridge_begin>support_bridge_end||
       support_bridge_end>original_bridge_nf)
        return -1;
    if(original_bridge_nf==0)
        return cbf_filter_round(arena,verts,nv,faces,nf,input_nf,
                                 vertex_chart,n_charts,out_faces,out_nf,NULL,0,
                                 0,0,stats);

    current_faces=(int32_t*)faces;
    remaining=(int32_t*)malloc(original_bridge_nf*3*sizeof(*remaining));
    label=(int32_t*)malloc(nv*sizeof(*label));
    if(!remaining||!label)goto cleanup_promote;
    memcpy(remaining,faces+input_nf*3,
           original_bridge_nf*3*sizeof(*remaining));
    memcpy(label,vertex_chart,nv*sizeof(*label));

    for(size_t round=0;remaining_nf>0;round++){
        size_t round_nf=current_nf+remaining_nf;
        int32_t *round_input=(int32_t*)malloc(round_nf*3*sizeof(*round_input));
        uint8_t *kept=(uint8_t*)calloc(remaining_nf,1);
        int32_t *round_out=NULL;
        size_t round_out_nf=0,accepted_faces;
        ChartBridgeForestStats one;
        if(!round_input||!kept){
            free(kept);free(round_input);goto cleanup_promote;
        }
        memcpy(round_input,current_faces,current_nf*3*sizeof(*round_input));
        memcpy(round_input+current_nf*3,remaining,
               remaining_nf*3*sizeof(*round_input));
        if(getenv("SEAM_CHART_DIAGNOSTICS")&&round>0)
            fprintf(stderr,
                    "  [chart-promotion] round %zu: %zu accepted face(s) "
                    "fixed, retrying %zu residual raw face(s) against %zu "
                    "coarsened chart(s)\n",
                    round+1,current_nf-input_nf,remaining_nf,current_charts);
        memset(&one,0,sizeof one);
        if(cbf_filter_round(arena,verts,nv,round_input,round_nf,current_nf,
                             label,current_charts,&round_out,&round_out_nf,
                             kept,round>0,
                             round==0?support_bridge_begin:0,
                             round==0?support_bridge_end:0,&one)!=0){
            free(kept);free(round_input);goto cleanup_promote;
        }
        if(round_out_nf<current_nf){
            free(kept);free(round_input);goto cleanup_promote;
        }
        accepted_faces=round_out_nf-current_nf;
        total.promotion_rounds++;
        total.patches+=one.patches;
        total.accepted_patches+=one.accepted_patches;
        total.accepted_join_patches+=one.accepted_join_patches;
        total.accepted_growth_patches+=one.accepted_growth_patches;
        total.accepted_hole_patches+=one.accepted_hole_patches;
        total.rejected_self_or_mixed+=one.rejected_self_or_mixed;
        total.rejected_patch_topology+=one.rejected_patch_topology;
        total.rejected_attachment+=one.rejected_attachment;
        total.rejected_cycle+=one.rejected_cycle;
        total.rejected_edge_conflict+=one.rejected_edge_conflict;
        total.rejected_orientation+=one.rejected_orientation;
        total.bridge_faces_kept+=accepted_faces;
        if(round==0&&support_bridge_begin<support_bridge_end){
            for(size_t i=support_bridge_begin;i<support_bridge_end;i++)
                if(kept[i]){
                    const int32_t *tri=&remaining[i*3];
                    if(label[tri[0]]==label[tri[1]]&&
                       label[tri[0]]==label[tri[2]])
                        total.accepted_source_shoulder_faces++;
                }
        }
        if(round>0){
            total.promoted_patches+=one.accepted_patches;
            total.promoted_faces+=accepted_faces;
        }
        current_faces=round_out;
        current_nf=round_out_nf;
        free(round_input);

        if(accepted_faces==0||getenv("SEAM_NO_CHART_PROMOTION")){
            free(kept);
            break;
        }
        {
            int32_t *parent=(int32_t*)malloc(current_charts*sizeof(*parent));
            uint8_t *rank=(uint8_t*)calloc(current_charts,1);
            int32_t *remap=(int32_t*)malloc(current_charts*sizeof(*remap));
            int32_t *next_remaining=NULL;
            size_t next_nf=remaining_nf-accepted_faces,w=0,next_charts=0;
            if(next_nf)
                next_remaining=(int32_t*)malloc(next_nf*3*
                                                sizeof(*next_remaining));
            if(!parent||!rank||!remap||(next_nf&&!next_remaining)){
                free(next_remaining);free(remap);free(rank);free(parent);
                free(kept);goto cleanup_promote;
            }
            for(size_t c=0;c<current_charts;c++){
                parent[c]=(int32_t)c;remap[c]=-1;
            }
            for(size_t i=0;i<remaining_nf;i++){
                const int32_t *tri=&remaining[i*3];
                if(kept[i]){
                    int32_t a=label[tri[0]],b=label[tri[1]],c=label[tri[2]];
                    cbf_union(parent,rank,a,b);
                    cbf_union(parent,rank,a,c);
                }else{
                    memcpy(&next_remaining[w*3],tri,3*sizeof(*tri));
                    w++;
                }
            }
            for(size_t v=0;v<nv;v++){
                int32_t root=cbf_find(parent,label[v]);
                if(remap[root]<0)remap[root]=(int32_t)next_charts++;
                label[v]=remap[root];
            }
            free(remaining);
            remaining=next_remaining;
            remaining_nf=next_nf;
            current_charts=next_charts;
            free(remap);free(rank);free(parent);
        }
        free(kept);
    }
    if(getenv("SEAM_CHART_DIAGNOSTICS")&&total.promotion_rounds>1)
        fprintf(stderr,
                "  [chart-promotion] fixed point after %zu round(s): "
                "%zu promoted patch(es), %zu promoted face(s), "
                "%zu residual raw face(s) rejected\n",
                total.promotion_rounds,total.promoted_patches,
                total.promoted_faces,remaining_nf);
    *out_faces=current_faces;
    *out_nf=current_nf;
    if(stats)*stats=total;
    rc=0;

cleanup_promote:
    free(label);
    free(remaining);
    return rc;
}

int ChartBridgeForest_filter(Arena_T arena,
                             const float *verts, size_t nv,
                             const int32_t *faces, size_t nf,
                             size_t input_nf,
                             const int32_t *vertex_chart,
                             size_t n_charts,
                             int32_t **out_faces, size_t *out_nf,
                             ChartBridgeForestStats *stats)
{
    return ChartBridgeForest_filter_with_support(
        arena,verts,nv,faces,nf,input_nf,vertex_chart,n_charts,0,0,
        out_faces,out_nf,stats);
}

int ChartBridgeForest_selftest(void)
{
    static const float v[]={
        0,0,0, 0,1,0, 0,0,1, 0,1,1,
        0,2,0, 0,3,0, 0,2,1, 0,3,1
    };
    static const int32_t f[]={
        0,1,2, 1,3,2, 4,5,6, 5,7,6,
        1,4,3, 3,4,6,
        0,2,5, 2,7,5
    };
    static const int32_t chart[]={0,0,0,0,1,1,1,1};
    Arena_T arena=Arena_new();
    int32_t *out=NULL;size_t nf=0;ChartBridgeForestStats st;
    int fail=ChartBridgeForest_filter(arena,v,8,f,8,4,chart,2,&out,&nf,&st)!=0||
        st.accepted_patches!=2||st.accepted_join_patches!=1||
        st.accepted_hole_patches!=1||st.rejected_cycle!=0||nf!=8;
    fprintf(stderr,"[selftest] chart bridge forest -> %s\n",fail?"FAIL":"ok");
    {
        static const float ov[]={
            -1,0,0, 0,0,0, 0,1,0,
             1,0,0, 1,1,0, 2,0,0
        };
        int32_t of[]={
            0,1,2, 5,3,4,
            1,3,2, 2,3,4
        };
        static const int32_t oc[]={0,0,0,1,1,1};
        ChartBridgeOrientationStats ost;
        ChartBridgeForestStats fst;
        int32_t *oout=NULL;size_t onf=0;
        int orientation_fail=
            ChartBridgeForest_orient_source(arena,ov,6,of,4,2,oc,2,&ost)!=0||
            ost.candidate_patches!=1||ost.constraints_used!=1||
            ost.constraints_conflicted!=0||ost.charts_flipped!=1||
            ost.faces_flipped!=1||of[3]!=5||of[4]!=4||of[5]!=3||
            ChartBridgeForest_filter(arena,ov,6,of,4,2,oc,2,
                                     &oout,&onf,&fst)!=0||
            fst.accepted_join_patches!=1||onf!=4;
        fprintf(stderr,"[selftest] chart bridge source orientation -> %s\n",
                orientation_fail?"FAIL":"ok");
        fail|=orientation_fail;
    }
    {
        /* The first bridge face is a source-chart shoulder: all of its
         * vertices carry chart 1, but it is edge-adjacent to exactly one 0/1
         * bridge strip.  Without the shoulder that strip only point-touches
         * chart 1 at vertices 4 and 6; with explicit provenance the complete
         * three-face patch attaches along one path on each disk. */
        static const float sv[]={
            0,0,0, 0,1,0, 0,0,1,
            0,2,0, 0,3,0, 0,3,1, 0,2,1
        };
        static const int32_t sf[]={
            0,1,2, 3,4,5, 3,5,6,
            4,3,6, 4,6,0, 1,0,6
        };
        static const int32_t sc[]={0,0,0,1,1,1,1};
        int32_t *plain_out=NULL,*support_out=NULL;
        size_t plain_nf=0,support_nf=0;
        ChartBridgeForestStats plain_st,support_st;
        int shoulder_fail=
            ChartBridgeForest_filter(arena,sv,7,sf,6,3,sc,2,
                                     &plain_out,&plain_nf,&plain_st)!=0||
            plain_st.accepted_join_patches!=0||plain_nf!=3||
            ChartBridgeForest_filter_with_support(
                arena,sv,7,sf,6,3,sc,2,0,1,
                &support_out,&support_nf,&support_st)!=0||
            support_st.accepted_join_patches!=1||support_nf!=6||
            support_st.accepted_source_shoulder_faces!=1;
        fprintf(stderr,"[selftest] chart bridge source shoulder -> %s\n",
                shoulder_fail?"FAIL":"ok");
        fail|=shoulder_fail;
    }
    Arena_dispose(&arena);return fail;
}
