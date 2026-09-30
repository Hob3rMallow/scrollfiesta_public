#include "quad_field.h"
#include "sparse_solve.h"
#include "ribbon_lod.h"
#include "../common/pipeline_constants.h"
#include "../common/ves_platform.h"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif

/* Hanging nodes are eliminated algebraically, C^T A C and C^T b, before
 * factorization. This is the usual constrained Q1 element, not a triangle
 * fan used as a surrogate finite element. Each cell stores only its 4x4
 * sufficient statistics, independent of the number of source observations. */
typedef struct { int32_t chart,x,y,id; } QfKey;
typedef struct { int32_t chart,axis,line,start,end,cell; } QfEdge;
typedef struct { int n; const int32_t *id; const double *w; } QfBasis;
typedef struct {
    int32_t chart,x,y,master,parent[2];
    int visiting;
    QfBasis basis;
} QfNode;
typedef struct {
    QuadFieldCell key;
    int32_t node[4],mid[4];
    uint8_t element;
    double a[16],b[12],sum_square,weight;
    double origin[3],centered_b[12],centered_square;
    size_t samples;
} QfCell;
struct QuadField_T {
    Arena_T arena;
    double step;
    size_t nc,nn,ndof,nsamples;
    size_t node_bytes;
    QfCell *cell;
    QfNode *node;
    QfKey *keys;
    QfEdge *edges;
    int32_t *directory; /* sparse leaf hash, index+1; zero is empty */
    size_t directory_size;
    uint32_t levels;
    double *value; /* node XYZ, including constrained nodes */
    uint8_t *fixed; /* optional original-boundary Dirichlet data, master order */
    double *fixed_value;
    int solved;
};

static uint64_t qf_hash(int32_t chart,int32_t x,int32_t y,int32_t size)
{
    uint64_t h=(uint32_t)chart;
    h=(h^(uint32_t)x)*UINT64_C(1099511628211);
    h=(h^(uint32_t)y)*UINT64_C(1099511628211);
    h=(h^(uint32_t)size)*UINT64_C(1099511628211);
    h^=h>>30; h*=UINT64_C(0xbf58476d1ce4e5b9); h^=h>>27;
    return h;
}
static int qf_cell_find(QuadField_T f,int32_t chart,int32_t x,int32_t y,int32_t size,int insert)
{
    size_t slot=(size_t)qf_hash(chart,x,y,size)&(f->directory_size-1);
    for (size_t i=0;i<f->directory_size;i++) {
        int32_t id=f->directory[slot];
        if (id==0) {
            if (insert>=0) f->directory[slot]=insert+1;
            return -1;
        }
        {
            const QuadFieldCell *c=&f->cell[id-1].key;
            if (c->chart==chart && c->x==x && c->y==y && c->size==size) return id-1;
        }
        slot=(slot+1)&(f->directory_size-1);
    }
    return -2;
}
static int32_t qf_floor_cell(int32_t x,int32_t size)
{
    int32_t q=x/size;
    return (x%size<0 ? q-1 : q)*size;
}

static int qf_key_cmp(const void *pa,const void *pb)
{
    const QfKey *a=(const QfKey *)pa,*b=(const QfKey *)pb;
    if (a->chart!=b->chart) return a->chart<b->chart ? -1 : 1;
    if (a->y!=b->y) return a->y<b->y ? -1 : 1;
    if (a->x!=b->x) return a->x<b->x ? -1 : 1;
    return 0;
}
static int qf_edge_cmp(const void *pa,const void *pb)
{
    const QfEdge *a=(const QfEdge *)pa,*b=(const QfEdge *)pb;
    if (a->chart!=b->chart) return a->chart<b->chart ? -1 : 1;
    if (a->axis!=b->axis) return a->axis<b->axis ? -1 : 1;
    if (a->line!=b->line) return a->line<b->line ? -1 : 1;
    if (a->start!=b->start) return a->start<b->start ? -1 : 1;
    if (a->end!=b->end) return a->end>b->end ? -1 : 1;
    return a->cell<b->cell ? -1 : a->cell>b->cell;
}
static int qf_find(QuadField_T f,int32_t chart,int32_t x,int32_t y)
{
    QfKey key={chart,x,y,0};
    size_t lo=0,hi=f->nn;
    while (lo<hi) { size_t mid=lo+(hi-lo)/2; if (qf_key_cmp(f->keys+mid,&key)<0) lo=mid+1; else hi=mid; }
    return lo<f->nn && qf_key_cmp(f->keys+lo,&key)==0 ? (int)lo : -1;
}
static int qf_expand(QuadField_T f,int node,int depth)
{
    static const double unit_weight=1;
    QfNode *n=f->node+node;
    if (n->basis.n) return 0;
    if (n->visiting || depth>32) return -1;
    n->visiting=1;
    if (n->master>=0) n->basis=(QfBasis){1,&n->master,&unit_weight};
    else {
        int32_t ids[32]={0};
        double weights[32]={0};
        int count=0;
        size_t offset=0,bytes=0;
        uint8_t *storage=NULL;
        for (int p=0;p<2;p++) {
            const QfBasis *b=NULL;
            if (n->parent[p]<0 || qf_expand(f,n->parent[p],depth+1)!=0) return -1;
            b=&f->node[n->parent[p]].basis;
            for (int j=0;j<b->n;j++) {
                int at=0;
                while (at<count && ids[at]!=b->id[j]) at++;
                if (at==32) return -1;
                if (at==count) ids[count++]=b->id[j];
                weights[at]+=.5*b->w[j];
            }
        }
        /* Masters need no stencil allocation. A hanging node stores only
         * its nonzeros, in exactly the old parent/accumulation order. Charge
         * the padded allocation against a byte cap, not the free-DOF cap. */
        offset=((size_t)count*sizeof(int32_t)+sizeof(double)-1)/sizeof(double)*sizeof(double);
        bytes=(offset+(size_t)count*sizeof(double)+15)&~(size_t)15;
        if (bytes>QUAD_FIELD_MAX_NODE_BYTES-f->node_bytes) return -1;
        storage=(uint8_t *)ARENA_ALLOC(f->arena,bytes);
        memcpy(storage,ids,(size_t)count*sizeof(int32_t));
        memcpy(storage+offset,weights,(size_t)count*sizeof(double));
        n->basis=(QfBasis){count,(const int32_t *)storage,(const double *)(storage+offset)};
        f->node_bytes+=bytes;
    }
    n->visiting=0;
    return 0;
}

int QuadField_new(Arena_T arena,const QuadFieldCell *cells,size_t count,
                   double step,QuadField_T *out)
{
    QuadField_T f=NULL;
    QfKey *all=NULL;
    QfEdge *edges=NULL;
    static const int dx[4]={0,1,1,0},dy[4]={0,0,1,1};
    if (arena==NULL || out==NULL || (count && cells==NULL) ||
        !isfinite(step) || step<=0) return -1;
    *out=NULL;
    if (count>QUAD_FIELD_MAX_CELLS) {
        fprintf(stderr,"[quad field] initial domain has %zu cells, above %u; refused before allocation\n",
            count,(unsigned)QUAD_FIELD_MAX_CELLS);
        return -1;
    }
    f=(QuadField_T)ARENA_CALLOC(arena,1,sizeof *f);
    f->arena=arena; f->step=step; f->nc=count;
    f->cell=(QfCell *)ARENA_CALLOC(arena,count,sizeof *f->cell);
    f->directory_size=1;
    while (f->directory_size<2*count) f->directory_size*=2;
    f->directory=(int32_t *)ARENA_CALLOC(arena,f->directory_size,sizeof(int32_t));
    all=(QfKey *)ARENA_ALLOC(arena,4*count*sizeof *all);
    edges=(QfEdge *)ARENA_ALLOC(arena,4*count*sizeof *edges);
    for (size_t i=0;i<count;i++) {
        const QuadFieldCell *c=cells+i;
        if (c->chart<0 || c->size<=0 || c->size>(1<<28) ||
            (c->size&(c->size-1)) || c->x%c->size || c->y%c->size ||
            c->x>INT32_MAX-c->size || c->y>INT32_MAX-c->size) return -1;
        f->cell[i].key=*c;
        if (qf_cell_find(f,c->chart,c->x,c->y,c->size,(int)i)!=-1) return -1;
        f->levels|=(uint32_t)c->size;
        for (int j=0;j<4;j++) {
            all[4*i+j]=(QfKey){c->chart,c->x+dx[j]*c->size,c->y+dy[j]*c->size,(int32_t)(4*i+j)};
            f->cell[i].mid[j]=-1;
        }
        edges[4*i]=(QfEdge){c->chart,0,c->y,c->x,c->x+c->size,(int32_t)i};
        edges[4*i+1]=(QfEdge){c->chart,1,c->x+c->size,c->y,c->y+c->size,(int32_t)i};
        edges[4*i+2]=(QfEdge){c->chart,0,c->y+c->size,c->x,c->x+c->size,(int32_t)i};
        edges[4*i+3]=(QfEdge){c->chart,1,c->x,c->y,c->y+c->size,(int32_t)i};
    }
    for (size_t i=0;i<count;i++) {
        const QuadFieldCell *c=cells+i;
        for (int32_t size=2*c->size;size<=(1<<28);size*=2) {
            if (((uint32_t)size&f->levels) &&
                qf_cell_find(f,c->chart,qf_floor_cell(c->x,size),qf_floor_cell(c->y,size),size,-1)>=0) return -1;
        }
    }
    qsort(all,4*count,sizeof *all,qf_key_cmp);
    for (size_t i=0;i<4*count;) {
        size_t j=i+1;
        while (j<4*count && qf_key_cmp(all+i,all+j)==0) j++;
        for (size_t k=i;k<j;k++) f->cell[all[k].id/4].node[all[k].id%4]=(int32_t)f->nn;
        all[f->nn++]=all[i]; i=j;
    }
    f->keys=all;
    if (f->nn>QUAD_FIELD_MAX_NODE_BYTES/sizeof *f->node) {
        fprintf(stderr,"[quad field] %zu geometric nodes exceed the node-storage byte budget; refused\n",f->nn);
        return -1;
    }
    f->node_bytes=(f->nn*sizeof *f->node+15)&~(size_t)15;
    if (f->node_bytes>QUAD_FIELD_MAX_NODE_BYTES) return -1;
    f->node=(QfNode *)ARENA_CALLOC(arena,f->nn,sizeof *f->node);
    f->value=(double *)ARENA_CALLOC(arena,3*f->nn,sizeof(double));
    for (size_t i=0;i<f->nn;i++) {
        f->node[i].chart=all[i].chart; f->node[i].x=all[i].x; f->node[i].y=all[i].y;
        f->node[i].parent[0]=f->node[i].parent[1]=-1;
    }
    qsort(edges,4*count,sizeof *edges,qf_edge_cmp);
    f->edges=edges;
    for (size_t i=0;i<4*count;i++) {
        const QfEdge *a=edges+i;
        int length=a->end-a->start;
        for (size_t j=i+1;j<4*count;j++) {
            const QfEdge *b=edges+j;
            if (a->chart!=b->chart || a->axis!=b->axis || a->line!=b->line || b->start>=a->end) break;
            if (length>2*(b->end-b->start)) return -1; /* non-2:1 interface */
        }
        if (length>1) {
            int mid=(int)(((int64_t)a->start+a->end)/2);
            int id=qf_find(f,a->chart,a->axis ? a->line : mid,a->axis ? mid : a->line);
            if (id>=0) {
                int p0=qf_find(f,a->chart,a->axis ? a->line : a->start,a->axis ? a->start : a->line);
                int p1=qf_find(f,a->chart,a->axis ? a->line : a->end,a->axis ? a->end : a->line);
                QfNode *n=f->node+id;
                if (n->parent[0]>=0 && (n->parent[0]!=p0 || n->parent[1]!=p1)) return -1;
                n->parent[0]=p0; n->parent[1]=p1;
            }
        }
    }
    for (size_t i=0;i<f->nn;i++) f->node[i].master=f->node[i].parent[0]<0 ? (int32_t)f->ndof++ : -1;
    /* Dirichlet data are supplied after construction. The factorization cap
     * applies to free unknowns in solve(), not to known source samples. */
    for (size_t i=0;i<f->nn;i++) if (qf_expand(f,(int)i,0)!=0) return -1;
    for (size_t i=0;i<count;i++) for (int j=0;j<4;j++) {
        int a=f->cell[i].node[j],b=f->cell[i].node[(j+1)%4];
        if (cells[i].size>1) f->cell[i].mid[j]=qf_find(f,cells[i].chart,
            (int32_t)(((int64_t)f->node[a].x+f->node[b].x)/2),
            (int32_t)(((int64_t)f->node[a].y+f->node[b].y)/2));
    }
    *out=f;
    return 0;
}

static int qf_weights(QuadField_T f,size_t index,double u,double v,double w[4])
{
    const QuadFieldCell *c=NULL;
    double x=0,y=0;
    if (f==NULL || index>=f->nc || !isfinite(u) || !isfinite(v)) return -1;
    c=&f->cell[index].key;
    x=(u/f->step-c->x)/c->size; y=(v/f->step-c->y)/c->size;
    if (x<0 || x>1 || y<0 || y>1) return -1;
    if (f->cell[index].element==1) {
        if (x>=y) { w[0]=1-x; w[1]=x-y; w[2]=y; w[3]=0; }
        else { w[0]=1-y; w[1]=0; w[2]=x; w[3]=y-x; }
    } else if (f->cell[index].element==2) {
        if (x+y<=1) { w[0]=1-x-y; w[1]=x; w[2]=0; w[3]=y; }
        else { w[0]=0; w[1]=1-y; w[2]=x+y-1; w[3]=1-x; }
    } else {
        w[0]=(1-x)*(1-y); w[1]=x*(1-y); w[2]=x*y; w[3]=(1-x)*y;
    }
    return 0;
}

int QuadField_set_element(QuadField_T f,size_t cell,int kind)
{
    if (!f || cell>=f->nc || f->solved || f->cell[cell].samples ||
        kind<0 || kind>2 || (kind && f->cell[cell].key.size!=1)) return -1;
    f->cell[cell].element=(uint8_t)kind;
    return 0;
}
int QuadField_element(QuadField_T f,size_t cell)
{ return f && cell<f->nc ? f->cell[cell].element : -1; }

size_t QuadField_observation_count(QuadField_T f,size_t cell)
{ return f && cell<f->nc ? f->cell[cell].samples : 0; }
int QuadField_reuse_observations(QuadField_T f,QuadField_T previous)
{
    if (!f || !previous || f==previous || f->solved || f->nsamples || f->step!=previous->step) return -1;
    for (size_t i=0;i<f->nc;i++) {
        QfCell *c=f->cell+i;
        int id=qf_cell_find(previous,c->key.chart,c->key.x,c->key.y,c->key.size,-1);
        const QfCell *old=id<0 ? NULL : previous->cell+id;
        if (!old || old->element!=c->element) continue;
        memcpy(c->a,old->a,sizeof c->a); memcpy(c->b,old->b,sizeof c->b);
        c->sum_square=old->sum_square; c->weight=old->weight; c->samples=old->samples;
        memcpy(c->origin,old->origin,sizeof c->origin);
        memcpy(c->centered_b,old->centered_b,sizeof c->centered_b);
        c->centered_square=old->centered_square;
        if (SIZE_MAX-f->nsamples<c->samples) return -1;
        f->nsamples+=c->samples;
    }
    return 0;
}

int QuadField_observe(QuadField_T f,size_t index,double u,double v,const double xyz[3],double weight)
{
    double w[4]={0};
    QfCell *c=NULL;
    if (xyz==NULL || !isfinite(weight) || weight<=0 || qf_weights(f,index,u,v,w)!=0 || f->solved) return -1;
    for (int k=0;k<3;k++) if (!isfinite(xyz[k])) return -1;
    c=f->cell+index;
    if (!c->samples) memcpy(c->origin,xyz,sizeof c->origin);
    for (int i=0;i<4;i++) {
        for (int j=0;j<4;j++) c->a[4*i+j]+=weight*w[i]*w[j];
        for (int k=0;k<3;k++) c->b[3*i+k]+=weight*w[i]*xyz[k];
        for (int k=0;k<3;k++) c->centered_b[3*i+k]+=weight*w[i]*(xyz[k]-c->origin[k]);
    }
    for (int k=0;k<3;k++) c->sum_square+=weight*xyz[k]*xyz[k];
    for (int k=0;k<3;k++) c->centered_square+=weight*(xyz[k]-c->origin[k])*(xyz[k]-c->origin[k]);
    c->weight+=weight; c->samples++; f->nsamples++;
    return 0;
}

int QuadField_locate(QuadField_T f,int32_t chart,double u,double v,size_t *out)
{
    double x=0,y=0;
    if (f==NULL || out==NULL || !isfinite(u) || !isfinite(v)) return -1;
    x=u/f->step; y=v/f->step;
    if (x<INT32_MIN || x>INT32_MAX || y<INT32_MIN || y>INT32_MAX) return -1;
    for (int32_t size=1;size<=(1<<28);size*=2) if (f->levels&(uint32_t)size) {
        int64_t ix=(int64_t)floor(x/size)*size,iy=(int64_t)floor(y/size)*size;
        int nx=x==(double)ix ? 2 : 1,ny=y==(double)iy ? 2 : 1;
        for (int j=0;j<ny;j++) for (int i=0;i<nx;i++) {
            int64_t px=ix-i*(int64_t)size,py=iy-j*(int64_t)size;
            int found=-1;
            if (px<INT32_MIN || py<INT32_MIN || px>INT32_MAX || py>INT32_MAX) continue;
            found=qf_cell_find(f,chart,(int32_t)px,(int32_t)py,size,-1);
            if (found>=0) { *out=(size_t)found; return 0; }
        }
    }
    return -1;
}

int QuadField_evaluate(QuadField_T f,size_t index,double u,double v,double xyz[3])
{
    double w[4]={0};
    if (xyz==NULL || qf_weights(f,index,u,v,w)!=0 || !f->solved) return -1;
    for (int k=0;k<3;k++) {
        xyz[k]=0;
        for (int i=0;i<4;i++) xyz[k]+=w[i]*f->value[3*(size_t)f->cell[index].node[i]+k];
    }
    return 0;
}

size_t QuadField_node_count(QuadField_T f) { return f ? f->nn : 0; }

int QuadField_node_key(QuadField_T f,size_t node,QuadFieldNodeKey *key)
{
    if (!f || !key || node>=f->nn) return -1;
    key->chart=f->node[node].chart;
    key->x=f->node[node].x; key->y=f->node[node].y;
    key->constrained=f->node[node].master<0;
    return 0;
}

int QuadField_find_node(QuadField_T f,int32_t chart,int32_t x,int32_t y,size_t *node)
{
    int found=-1;
    if (!f || !node) return -1;
    found=qf_find(f,chart,x,y);
    if (found<0) return 1;
    *node=(size_t)found;
    return 0;
}

int QuadField_fix_node(QuadField_T f,size_t node,const double xyz[3])
{
    size_t id=0;
    if (!f || !xyz || f->solved || node>=f->nn || f->node[node].master<0) return -1;
    for (int k=0;k<3;k++) if (!isfinite(xyz[k])) return -1;
    id=(size_t)f->node[node].master;
    if (!f->fixed) {
        f->fixed=(uint8_t *)ARENA_CALLOC(f->arena,f->ndof,1);
        f->fixed_value=(double *)ARENA_CALLOC(f->arena,3*f->ndof,sizeof(double));
    }
    if (f->fixed[id]) {
        for (int k=0;k<3;k++)
            if (fabs(f->fixed_value[3*id+k]-xyz[k])>1e-10) return -1;
        return 0;
    }
    f->fixed[id]=1;
    memcpy(f->fixed_value+3*id,xyz,3*sizeof(double));
    return 0;
}

int QuadField_solve(QuadField_T f,QuadField_T coarse,QuadFieldReport *report)
{
    Arena_Mark mark;
    size_t nt=0,cap=0,n=0,nfree=0;
    int *row=NULL,*col=NULL,*free_id=NULL;
    double *a=NULL,*rhs=NULL,*x=NULL,*diag=NULL,*seed=NULL,*res=NULL,*mass=NULL;
    SparseFactor_T factor=NULL;
    int rc=-1;
    if (f==NULL || report==NULL || f->solved || (coarse && !coarse->solved)) return -1;
    memset(report,0,sizeof *report);
    report->cells=f->nc; report->nodes=f->nn; report->unknowns=f->ndof;
    report->node_storage_bytes=f->node_bytes;
    report->hanging_nodes=f->nn-f->ndof; report->samples=f->nsamples;
    n=f->ndof;
    if (n==0) { f->solved=1; return 0; }
    for (size_t i=0;i<n;i++) report->fixed_nodes+=f->fixed && f->fixed[i];
    nfree=n-report->fixed_nodes;
    report->unknowns=nfree;
    if (nfree>QUAD_FIELD_MAX_UNKNOWNS) {
        fprintf(stderr,"[quad field] %zu FREE unknowns (%zu exact source values excluded) exceed direct-factor limit %u; requires a bounded multilevel linear solver\n",
            nfree,report->fixed_nodes,(unsigned)QUAD_FIELD_MAX_UNKNOWNS);
        return -1;
    }
    mark=Arena_save(f->arena);
    free_id=(int *)ARENA_ALLOC(f->arena,n*sizeof(int));
    for (size_t i=0,next=0;i<n;i++)
        free_id[i]=f->fixed && f->fixed[i] ? -1 : (int)next++;
    for (size_t c=0;c<f->nc;c++) for (int i=0;i<4;i++) for (int j=0;j<4;j++)
        cap+=(size_t)f->node[f->cell[c].node[i]].basis.n*(size_t)f->node[f->cell[c].node[j]].basis.n;
    cap+=n;
    if (cap>INT_MAX || cap>SIZE_MAX/sizeof(double)) goto done;
    row=(int *)ARENA_ALLOC(f->arena,cap*sizeof(int)); col=(int *)ARENA_ALLOC(f->arena,cap*sizeof(int));
    a=(double *)ARENA_ALLOC(f->arena,cap*sizeof(double));
    rhs=(double *)ARENA_CALLOC(f->arena,3*n,sizeof(double));
    x=(double *)ARENA_CALLOC(f->arena,3*n,sizeof(double));
    seed=(double *)ARENA_CALLOC(f->arena,3*n,sizeof(double));
    diag=(double *)ARENA_CALLOC(f->arena,n,sizeof(double));
    res=(double *)ARENA_CALLOC(f->arena,3*n,sizeof(double));
    mass=(double *)ARENA_CALLOC(f->arena,n,sizeof(double));
    for (size_t c=0;c<f->nc;c++) for (int i=0;i<4;i++) {
        QfBasis *bi=&f->node[f->cell[c].node[i]].basis;
        for (int p=0;p<bi->n;p++) {
            int r=bi->id[p];
            for (int k=0;k<3;k++) rhs[(size_t)k*n+r]+=bi->w[p]*f->cell[c].b[3*i+k];
            for (int j=0;j<4;j++) {
                QfBasis *bj=&f->node[f->cell[c].node[j]].basis;
                for (int q=0;q<bj->n;q++) if (r>=bj->id[q]) {
                    double v=bi->w[p]*bj->w[q]*f->cell[c].a[4*i+j];
                    row[nt]=r; col[nt]=bj->id[q]; a[nt++]=v;
                    if (r==bj->id[q]) diag[r]+=v;
                }
            }
        }
    }
    for (size_t i=0;i<nt;i++) {
        mass[row[i]]+=a[i];
        if (row[i]!=col[i]) mass[col[i]]+=a[i];
    }
    if (!coarse) for (size_t begin=0;begin<f->nn;) {
        size_t end=begin+1;
        double total=0,average[3]={0};
        while (end<f->nn && f->node[end].chart==f->node[begin].chart) end++;
        for (size_t i=begin;i<end;i++) if (f->node[i].master>=0) {
            size_t id=(size_t)f->node[i].master;
            total+=mass[id];
            for (int k=0;k<3;k++) average[k]+=rhs[(size_t)k*n+id];
        }
        for (size_t i=begin;i<end;i++) if (f->node[i].master>=0) {
            size_t id=(size_t)f->node[i].master;
            for (int k=0;k<3;k++) seed[(size_t)k*n+id]=mass[id]>0 ? rhs[(size_t)k*n+id]/mass[id] : total>0 ? average[k]/total : 0;
        }
        begin=end;
    }
    /* A tiny scale-relative regularizer only anchors under-observed modes.
     * On refinement it anchors to the prolongated coarse solution; every
     * master node is still free to move in the same global fine system. */
    for (size_t i=0;i<f->nn;i++) if (f->node[i].master>=0) {
        QfNode *node=f->node+i;
        size_t id=(size_t)node->master,c=0;
        double xyz[3]={seed[id],seed[n+id],seed[2*n+id]},prior=0;
        if (coarse) {
            if (QuadField_locate(coarse,node->chart,node->x*f->step,node->y*f->step,&c)!=0 ||
                QuadField_evaluate(coarse,c,node->x*f->step,node->y*f->step,xyz)!=0) goto done;
        }
        prior=QUAD_FIELD_NULLSPACE_PIN*fmax(diag[id],1.0);
        row[nt]=col[nt]=(int)id; a[nt++]=prior;
        for (int k=0;k<3;k++) { seed[(size_t)k*n+id]=xyz[k]; rhs[(size_t)k*n+id]+=prior*xyz[k]; }
    }
    /* Symmetric Dirichlet elimination, not a large penalty. Original support
     * boundaries remain exact; every unconstrained interior master is free in
     * this GLOBAL system, including at coarse/fine interfaces. */
    if (f->fixed) {
        size_t kept=0;
        for (size_t t=0;t<nt;t++) {
            size_t r=(size_t)row[t],c=(size_t)col[t];
            if (f->fixed[r] || f->fixed[c]) {
                if (r!=c) for (int k=0;k<3;k++) {
                    if (!f->fixed[r]) rhs[(size_t)k*n+r]-=a[t]*f->fixed_value[3*c+k];
                    if (!f->fixed[c]) rhs[(size_t)k*n+c]-=a[t]*f->fixed_value[3*r+k];
                }
            } else {
                row[kept]=row[t]; col[kept]=col[t]; a[kept]=a[t]; kept++;
            }
        }
        nt=kept;
        for (size_t i=0;i<n;i++) if (f->fixed[i]) {
            for (int k=0;k<3;k++)
                rhs[(size_t)k*n+i]=seed[(size_t)k*n+i]=f->fixed_value[3*i+k];
        }
    }
    /* Solve the fine residual about P*x_coarse. This is actual nested
     * iteration, not a dense fine solve followed by mesh decimation. */
    for (size_t t=0;t<nt;t++) for (int k=0;k<3;k++) {
        size_t base=(size_t)k*n;
        rhs[base+row[t]]-=a[t]*seed[base+col[t]];
        if (row[t]!=col[t]) rhs[base+col[t]]-=a[t]*seed[base+row[t]];
    }
    /* Remove the known-value identity block entirely. Monotone compaction
     * preserves the lower-triangle order and in-place RHS reads. The global
     * free system and hanging constraints are unchanged; this is not a set
     * of independent patches or a higher factorization cap. */
    for (int k=0;k<3;k++) for (size_t i=0;i<n;i++) if (free_id[i]>=0)
        rhs[(size_t)k*nfree+(size_t)free_id[i]]=rhs[(size_t)k*n+i];
    for (size_t t=0;t<nt;t++) {
        row[t]=free_id[row[t]]; col[t]=free_id[col[t]];
        if (row[t]<0 || col[t]<0) goto done;
    }
    if (nfree && (Sparse_factor_spd((int)nfree,(int)nt,row,col,a,&factor)!=0 ||
                  Sparse_factor_solve_multi(factor,rhs,x,3)!=0)) goto done;
    {
        double r2=0,b2=0;
        memcpy(res,rhs,3*nfree*sizeof(double));
        for (size_t t=0;t<nt;t++) for (int k=0;k<3;k++) {
            size_t base=(size_t)k*nfree;
            res[base+row[t]]-=a[t]*x[base+col[t]];
            if (row[t]!=col[t]) res[base+col[t]]-=a[t]*x[base+row[t]];
        }
        for (size_t i=0;i<3*nfree;i++) { r2+=res[i]*res[i]; b2+=rhs[i]*rhs[i]; }
        report->relative_residual=sqrt(r2)/fmax(sqrt(b2),1.0);
        if (!isfinite(report->relative_residual) || report->relative_residual>1e-7) goto done;
    }
    /* Expand backwards: each destination is at or above its compact source,
     * so no unread correction is overwritten. Fixed values bypass the solve. */
    for (int k=2;k>=0;k--) for (size_t i=n;i-- >0;)
        x[(size_t)k*n+i]=free_id[i]<0 ? f->fixed_value[3*i+k] :
            x[(size_t)k*nfree+(size_t)free_id[i]]+seed[(size_t)k*n+i];
    for (size_t i=0;i<f->nn;i++) for (int k=0;k<3;k++) {
        double value=0;
        QfBasis *b=&f->node[i].basis;
        for (int j=0;j<b->n;j++) value+=b->w[j]*x[(size_t)k*n+b->id[j]];
        if (!isfinite(value)) goto done;
        f->value[3*i+k]=value;
    }
    {
        double err=0,weight=0;
        for (size_t c=0;c<f->nc;c++) {
            QfCell *cell=f->cell+c;
            /* The same integral as replaying every observation, but in a
             * local frame: world-coordinate cancellation must not turn a
             * sub-voxel residual into zero (or a spurious large error). */
            double sum=cell->centered_square;
            for (int i=0;i<4;i++) for (int k=0;k<3;k++) {
                double v=f->value[3*(size_t)cell->node[i]+k]-cell->origin[k];
                sum-=2*v*cell->centered_b[3*i+k];
                for (int j=0;j<4;j++) sum+=v*cell->a[4*i+j]*(f->value[3*(size_t)cell->node[j]+k]-cell->origin[k]);
            }
            err+=fmax(sum,0); weight+=cell->weight;
        }
        report->rms_error=weight>0 ? sqrt(err/weight) : 0;
        report->observation_weight=weight;
        report->maximum_error=NAN; /* requires a separate source residual pass */
    }
    f->solved=1; rc=0;
done:
    Sparse_factor_free(&factor); Arena_restore(f->arena,mark);
    return rc;
}

int QuadField_mesh(Arena_T arena,QuadField_T f,MeshBinData *out)
{ return QuadField_mesh_selection(arena,f,NULL,out,NULL); }

static void qf_face_normal(const float *xyz,const int32_t *tri,double normal[3])
{
    double a[3]={0},b[3]={0};
    for (int k=0;k<3;k++) {
        a[k]=(double)xyz[3*(size_t)tri[1]+k]-xyz[3*(size_t)tri[0]+k];
        b[k]=(double)xyz[3*(size_t)tri[2]+k]-xyz[3*(size_t)tri[0]+k];
    }
    normal[0]=a[1]*b[2]-a[2]*b[1];
    normal[1]=a[2]*b[0]-a[0]*b[2];
    normal[2]=a[0]*b[1]-a[1]*b[0];
}

int QuadField_project_mesh(Arena_T arena,QuadField_T f,const MeshBinData *in,
                            const int32_t *chart,MeshBinData *out,int32_t **original,
                            QuadFieldProjectionReport *report)
{
    int32_t *map=NULL;
    double square=0;
    if (!arena || !f || !f->solved || !in || !out || in==out || !original || !report ||
        (in->nv && (!in->verts || !in->uv || !chart)) || (in->nf && !in->faces) ||
        in->nv>INT32_MAX || in->nv>SIZE_MAX/(3*sizeof(float)) ||
        in->nf>SIZE_MAX/(3*sizeof(int32_t))) return -1;
    memset(out,0,sizeof *out); memset(report,0,sizeof *report); *original=NULL;
    map=(int32_t *)ARENA_ALLOC(arena,in->nv*sizeof(int32_t));
    for (size_t i=0;i<in->nv;i++) {
        map[i]=-1;
        if (chart[i]>=0) map[i]=(int32_t)out->nv++;
    }
    for (size_t i=0;i<in->nf;i++) {
        const int32_t *t=in->faces+3*i;
        for (int k=0;k<3;k++) if (t[k]<0 || (size_t)t[k]>=in->nv) return -1;
        if (chart[t[0]]!=chart[t[1]] || chart[t[0]]!=chart[t[2]]) return -1;
        if (map[t[0]]>=0) out->nf++;
    }
    out->verts=(float *)ARENA_ALLOC(arena,3*out->nv*sizeof(float));
    out->uv=(float *)ARENA_ALLOC(arena,2*out->nv*sizeof(float));
    out->faces=(int32_t *)ARENA_ALLOC(arena,3*out->nf*sizeof(int32_t));
    *original=(int32_t *)ARENA_ALLOC(arena,out->nv*sizeof(int32_t));
    for (size_t i=0;i<in->nv;i++) if (map[i]>=0) {
        size_t at=(size_t)map[i],cell=0;
        double xyz[3]={0},err=0;
        if (QuadField_locate(f,chart[i],in->uv[2*i],in->uv[2*i+1],&cell)!=0 ||
            QuadField_evaluate(f,cell,in->uv[2*i],in->uv[2*i+1],xyz)!=0) return -1;
        for (int k=0;k<3;k++) {
            float stored=0;
            double d=0;
            if (!isfinite(xyz[k]) || fabs(xyz[k])>FLT_MAX || !isfinite(in->verts[3*i+k])) return -1;
            stored=(float)xyz[k];
            out->verts[3*at+k]=stored;
            d=(double)stored-in->verts[3*i+k]; err+=d*d;
        }
        square+=err;
        if (sqrt(err)>report->maximum_error) report->maximum_error=sqrt(err);
        memcpy(out->uv+2*at,in->uv+2*i,2*sizeof(float));
        (*original)[at]=(int32_t)i;
    }
    for (size_t i=0,at=0;i<in->nf;i++) if (map[in->faces[3*i]]>=0) {
        const int32_t *src=in->faces+3*i;
        int32_t *dst=out->faces+3*at++;
        double a[3]={0},b[3]={0},dot=0,norm=0;
        for (int k=0;k<3;k++) dst[k]=map[src[k]];
        qf_face_normal(in->verts,src,a); qf_face_normal(out->verts,dst,b);
        for (int k=0;k<3;k++) { dot+=a[k]*b[k]; norm+=b[k]*b[k]; }
        if (dot<0) report->reversed_faces++;
        if (norm==0) report->collapsed_faces++;
    }
    report->vertices=out->nv; report->faces=out->nf;
    report->excluded_vertices=in->nv-out->nv; report->excluded_faces=in->nf-out->nf;
    report->rms_error=out->nv ? sqrt(square/out->nv) : 0;
    return 0;
}

int QuadField_mesh_selection(Arena_T arena,QuadField_T f,const uint8_t *selected,
                              MeshBinData *out,int32_t **source_node)
{
    size_t nf=0,nv=0;
    uint8_t *used=NULL;
    int32_t *map=NULL;
    if (arena==NULL || f==NULL || out==NULL || !f->solved) return -1;
    memset(out,0,sizeof *out);
    used=(uint8_t *)ARENA_CALLOC(arena,f->nn,1);
    map=(int32_t *)ARENA_ALLOC(arena,f->nn*sizeof(int32_t));
    /* The native output retains QUADS in the field. Inspection uses the same
     * published 16 restricted-quadtree templates as RibbonLod. */
    for (size_t i=0;i<f->nc;i++) {
        int mask=0,tri[6][3]={{0}},nt=0;
        if (selected && !selected[i]) continue;
        for (int j=0;j<4;j++) { used[f->cell[i].node[j]]=1; if (f->cell[i].mid[j]>=0) used[f->cell[i].mid[j]]=1; }
        for (int j=0;j<4;j++) if (f->cell[i].mid[j]>=0) mask|=1<<j;
        if (RibbonLod_transition(mask,tri,&nt)!=0) return -1;
        nf+=(size_t)nt;
    }
    for (size_t i=0;i<f->nn;i++) { map[i]=-1; if (used[i]) map[i]=(int32_t)nv++; }
    out->nv=nv; out->nf=nf;
    out->verts=(float *)ARENA_ALLOC(arena,3*nv*sizeof(float));
    out->uv=(float *)ARENA_ALLOC(arena,2*nv*sizeof(float));
    out->faces=(int32_t *)ARENA_ALLOC(arena,3*nf*sizeof(int32_t));
    if (source_node) *source_node=(int32_t *)ARENA_ALLOC(arena,nv*sizeof(int32_t));
    for (size_t i=0;i<f->nn;i++) if (used[i]) {
        size_t at=(size_t)map[i];
        if (source_node) (*source_node)[at]=(int32_t)i;
        for (int k=0;k<3;k++) out->verts[3*at+k]=(float)f->value[3*i+k];
        out->uv[2*at]=(float)(f->node[i].x*f->step); out->uv[2*at+1]=(float)(f->node[i].y*f->step);
    }
    nf=0;
    for (size_t i=0;i<f->nc;i++) {
        int mask=0,tri[6][3]={{0}},nt=0,ids[8]={0};
        if (selected && !selected[i]) continue;
        for (int j=0;j<4;j++) {
            ids[j]=map[f->cell[i].node[j]]; ids[j+4]=f->cell[i].mid[j]>=0 ? map[f->cell[i].mid[j]] : -1;
            if (ids[j+4]>=0) mask|=1<<j;
        }
        if (RibbonLod_transition(mask,tri,&nt)!=0) return -1;
        if (f->cell[i].element) {
            static const int source_tri[2][2][3]={{{0,1,2},{0,2,3}},{{0,1,3},{1,2,3}}};
            if (mask || nt!=2) return -1;
            memcpy(tri,source_tri[f->cell[i].element-1],sizeof source_tri[0]);
        }
        for (int t=0;t<nt;t++) { for (int k=0;k<3;k++) out->faces[3*nf+k]=ids[tri[t][k]]; nf++; }
    }
    return 0;
}

size_t QuadField_cell_count(QuadField_T f) { return f ? f->nc : 0; }
const QuadFieldCell *QuadField_cell(QuadField_T f,size_t c)
{ return f && c<f->nc ? &f->cell[c].key : NULL; }

typedef struct {
    uint64_t magic,version,source,nc,nn,checksum;
    double step;
} QfStore;
static uint64_t qf_checksum(uint64_t h,const void *data,size_t bytes)
{
    const uint8_t *p=(const uint8_t *)data;
    for (size_t i=0;i<bytes;i++) { h^=p[i]; h*=UINT64_C(1099511628211); }
    return h;
}
int QuadField_save(QuadField_T f,const char *path,uint64_t source)
{
    QfStore h={0};
    uint64_t hash=UINT64_C(14695981039346656037);
    char temporary[2048];
    FILE *file=NULL;
    int ok=1;
    if (!f || !f->solved || !path || strlen(path)>2000) return -1;
    h.magic=UINT64_C(0x31444c4549465156); h.version=1; h.source=source;
    h.nc=f->nc; h.nn=f->nn; h.step=f->step;
    for (size_t i=0;i<f->nc;i++) if (f->cell[i].element) h.version=2;
    hash=qf_checksum(hash,&h,sizeof h);
    for (size_t i=0;i<f->nc;i++) hash=qf_checksum(hash,&f->cell[i].key,sizeof(QuadFieldCell));
    if (h.version==2) for (size_t i=0;i<f->nc;i++) hash=qf_checksum(hash,&f->cell[i].element,1);
    hash=qf_checksum(hash,f->value,3*f->nn*sizeof(double)); h.checksum=hash;
    snprintf(temporary,sizeof temporary,"%s.tmp",path);
    file=fopen(temporary,"wb");
    if (!file) return -1;
    ok=fwrite(&h,sizeof h,1,file)==1;
    for (size_t i=0;ok && i<f->nc;i++) ok=fwrite(&f->cell[i].key,sizeof(QuadFieldCell),1,file)==1;
    if (h.version==2) for (size_t i=0;ok && i<f->nc;i++) ok=fwrite(&f->cell[i].element,1,1,file)==1;
    if (ok) ok=fwrite(f->value,sizeof(double),3*f->nn,file)==3*f->nn;
    if (fclose(file)!=0) ok=0;
#ifdef _WIN32
    if (ok) ok=MoveFileExA(temporary,path,MOVEFILE_REPLACE_EXISTING|MOVEFILE_WRITE_THROUGH)!=0;
#else
    if (ok) ok=rename(temporary,path)==0;
#endif
    if (!ok) remove(temporary);
    return ok ? 0 : -1;
}
int QuadField_load(Arena_T arena,const char *path,uint64_t source,QuadField_T *out)
{
    QfStore h={0};
    QuadFieldCell *cells=NULL;
    QuadField_T f=NULL;
    uint64_t hash=UINT64_C(14695981039346656037),expected=0;
    FILE *file=NULL;
    int rc=-1;
    if (!arena || !path || !out) return -1;
    *out=NULL;
    file=fopen(path,"rb");
    if (!file) return -1;
    if (fread(&h,sizeof h,1,file)!=1 || h.magic!=UINT64_C(0x31444c4549465156) ||
        (h.version!=1 && h.version!=2) || h.source!=source || h.nc>QUAD_FIELD_MAX_CELLS ||
        h.nn>4*h.nc || h.nn>QUAD_FIELD_MAX_NODE_BYTES/sizeof(QfNode) ||
        !isfinite(h.step) || h.step<=0) goto done;
    /* Loading does not factor a free system. Bound serialized node claims by
     * the actual corner count and storage budget; new() also charges hanging
     * stencils and must reconstruct exactly the claimed number of nodes. */
    expected=h.checksum; h.checksum=0; hash=qf_checksum(hash,&h,sizeof h);
    cells=(QuadFieldCell *)ARENA_ALLOC(arena,(size_t)h.nc*sizeof *cells);
    if (fread(cells,sizeof *cells,(size_t)h.nc,file)!=(size_t)h.nc) goto done;
    hash=qf_checksum(hash,cells,(size_t)h.nc*sizeof *cells);
    if (QuadField_new(arena,cells,(size_t)h.nc,h.step,&f)!=0 || f->nn!=h.nn) goto done;
    if (h.version==2) for (size_t i=0;i<f->nc;i++) {
        uint8_t kind=0;
        if (fread(&kind,1,1,file)!=1 || QuadField_set_element(f,i,kind)!=0) goto done;
        hash=qf_checksum(hash,&kind,1);
    }
    if (fread(f->value,sizeof(double),3*f->nn,file)!=3*f->nn || fgetc(file)!=EOF) goto done;
    hash=qf_checksum(hash,f->value,3*f->nn*sizeof(double));
    if (hash!=expected) goto done;
    for (size_t i=0;i<3*f->nn;i++) if (!isfinite(f->value[i])) goto done;
    f->solved=1; *out=f; rc=0;
done:
    if (fclose(file)!=0) { *out=NULL; rc=-1; }
    return rc;
}

int QuadField_refine(Arena_T arena,QuadField_T f,const uint8_t *marked,QuadField_T *out)
{
    uint8_t *split=NULL;
    QuadFieldCell *cells=NULL;
    size_t count=0,at=0;
    int changes=1;
    if (arena==NULL || f==NULL || out==NULL || (f->nc && marked==NULL)) return -1;
    *out=NULL;
    split=(uint8_t *)ARENA_ALLOC(arena,f->nc);
    for (size_t i=0;i<f->nc;i++) {
        if (marked[i] && f->cell[i].key.size<2) return -1;
        split[i]=(uint8_t)(marked[i]!=0);
    }
    while (changes) {
        changes=0;
        for (size_t i=0;i<4*f->nc;i++) {
            const QfEdge *a=f->edges+i;
            for (size_t j=i+1;j<4*f->nc;j++) {
                const QfEdge *b=f->edges+j;
                int sa=0,sb=0;
                if (a->chart!=b->chart || a->axis!=b->axis || a->line!=b->line || b->start>=a->end) break;
                sa=f->cell[a->cell].key.size/(split[a->cell] ? 2 : 1);
                sb=f->cell[b->cell].key.size/(split[b->cell] ? 2 : 1);
                if (sa>2*sb) { if (split[a->cell]) return -1; split[a->cell]=1; changes++; }
                if (sb>2*sa) { if (split[b->cell]) return -1; split[b->cell]=1; changes++; }
            }
        }
    }
    for (size_t i=0;i<f->nc;i++) count+=split[i] ? 4 : 1;
    if (count>QUAD_FIELD_MAX_CELLS) {
        fprintf(stderr,"[quad field] proposed refinement has %zu cells, above %u; retained the last complete level\n",count,(unsigned)QUAD_FIELD_MAX_CELLS);
        return -1;
    }
    cells=(QuadFieldCell *)ARENA_ALLOC(arena,count*sizeof *cells);
    for (size_t i=0;i<f->nc;i++) {
        const QuadFieldCell *c=&f->cell[i].key;
        if (!split[i]) cells[at++]=*c;
        else for (int y=0;y<2;y++) for (int x=0;x<2;x++)
            cells[at++]=(QuadFieldCell){c->chart,c->x+x*c->size/2,c->y+y*c->size/2,c->size/2};
    }
    if (QuadField_new(arena,cells,count,f->step,out)!=0) return -1;
    at=0;
    for (size_t i=0;i<f->nc;i++) {
        if (!split[i]) (*out)->cell[at++].element=f->cell[i].element;
        else at+=4;
    }
    return 0;
}

typedef struct { int checkpoints; size_t quads[8]; double rms[8]; } QfStreamTest;
static int qf_test_reader(void *context,Arena_T arena,size_t chunk,MeshBinData *m,int32_t **chart)
{
    (void)context;
    if (chunk>=2) return -1;
    memset(m,0,sizeof *m); m->nv=128;
    m->verts=(float *)ARENA_ALLOC(arena,3*m->nv*sizeof(float));
    m->uv=(float *)ARENA_ALLOC(arena,2*m->nv*sizeof(float));
    *chart=(int32_t *)ARENA_CALLOC(arena,m->nv,sizeof(int32_t));
    for (int y=0;y<16;y++) for (int x=0;x<8;x++) {
        size_t i=(size_t)y*8+x;
        double u=.125+.25*(x+8*chunk),v=.125+.25*y;
        m->uv[2*i]=(float)u; m->uv[2*i+1]=(float)v;
        m->verts[3*i]=(float)(4352+v); m->verts[3*i+1]=(float)(3000+u);
        m->verts[3*i+2]=(float)(2000+.3*u*u+.2*v*v);
    }
    return 0;
}
static int qf_test_checkpoint(void *context,int level,QuadField_T f,const QuadFieldReport *r)
{
    QfStreamTest *test=(QfStreamTest *)context;
    Arena_T arena=Arena_new();
    MeshBinData mesh={0};
    int ok=level<8 && level==test->checkpoints && QuadField_mesh(arena,f,&mesh)==0;
    if (ok) {
        test->quads[level]=r->cells; test->rms[level]=r->rms_error; test->checkpoints++;
        ok=mesh.nf>=2*r->cells && mesh.nv==r->nodes && r->relative_residual<1e-7;
    }
    Arena_dispose(&arena); return ok ? 0 : -1;
}
static int qf_test_persistence(Arena_T arena,QuadField_T f)
{
    char path[2048]={0};
    QuadField_T loaded=NULL;
    MeshBinData a={0},b={0};
    FILE *file=NULL;
    int rc=-1;
    (void)ves_mkdir("output");
#ifdef _WIN32
    if (GetTempFileNameA("output","qf",0,path)==0) return -1;
#else
    {
        int fd=-1;
        memcpy(path,"output/quad_field_XXXXXX",sizeof "output/quad_field_XXXXXX");
        fd=mkstemp(path); if (fd<0) return -1; close(fd);
    }
#endif
    if (QuadField_save(f,path,UINT64_C(12345))!=0 ||
        QuadField_load(arena,path,UINT64_C(12345),&loaded)!=0 ||
        QuadField_mesh(arena,f,&a)!=0 || QuadField_mesh(arena,loaded,&b)!=0 ||
        a.nv!=b.nv || a.nf!=b.nf || memcmp(a.verts,b.verts,3*a.nv*sizeof(float))!=0 ||
        memcmp(a.uv,b.uv,2*a.nv*sizeof(float))!=0 ||
        memcmp(a.faces,b.faces,3*a.nf*sizeof(int32_t))!=0 ||
        memcmp(f->value,loaded->value,3*f->nn*sizeof(double))!=0 ||
        QuadField_load(arena,path,UINT64_C(12346),&loaded)==0) goto done;
    file=fopen(path,"r+b");
    if (!file || fseek(file,-1,SEEK_END)!=0) goto done;
    {
        int byte=fgetc(file);
        if (byte==EOF || fseek(file,-1,SEEK_END)!=0 || fputc(byte^1,file)==EOF) goto done;
    }
    if (fclose(file)!=0) { file=NULL; goto done; } file=NULL;
    if (QuadField_load(arena,path,UINT64_C(12345),&loaded)==0) goto done;
    rc=0;
done:
    if (file) fclose(file);
    remove(path); /* Only the unique file created above. */
    return rc;
}

static int qf_test_many_known_values(void)
{
    Arena_T arena=Arena_new();
    const size_t nc=QUAD_FIELD_MAX_UNKNOWNS/2+1;
    QuadFieldCell *cells=(QuadFieldCell *)ARENA_ALLOC(arena,nc*sizeof *cells);
    QuadField_T field=NULL;
    QuadFieldReport report={0};
    int fail=0;
    for (size_t c=0;c<nc;c++) cells[c]=(QuadFieldCell){0,(int32_t)(2*c),0,1};
    if (QuadField_new(arena,cells,nc,1,&field)!=0) { fail++; goto done; }
    /* More than a million geometric nodes and the direct-factor cap, but
     * just one free unknown. Compact master stencils allocate no arrays. */
    for (size_t i=0;i+1<field->nn;i++) {
        double value[3]={field->node[i].x,field->node[i].y,7};
        if (QuadField_fix_node(field,i,value)!=0) { fail++; goto done; }
    }
    if (QuadField_solve(field,NULL,&report)!=0) { fail++; goto done; }
    fail+=report.unknowns!=1 || report.fixed_nodes!=4*nc-1 || report.relative_residual!=0;
    fail+=report.node_storage_bytes!=((field->nn*sizeof(QfNode)+15)&~(size_t)15);
    for (size_t i=0;i<field->nn;i++) {
        double expected[3]={field->node[i].x,field->node[i].y,7};
        if (i+1==field->nn) expected[0]=expected[1]=expected[2]=0;
        fail+=memcmp(field->value+3*i,expected,sizeof expected)!=0;
    }
    /* A constructible/solved field must also survive the native reader.
     * Its geometry is not limited by the number of free factor unknowns. */
    fail+=qf_test_persistence(arena,field)!=0;
done:
    fprintf(stderr,"[selftest] quad field known-value reduction %s (%d failures)\n",fail ? "FAIL" : "PASS",fail);
    Arena_dispose(&arena); return fail;
}

int QuadField_selftest(void)
{
    Arena_T arena=Arena_new();
    QuadFieldCell coarse_cells[2]={{0,0,0,2},{0,2,0,2}};
    QuadFieldCell fine_cells[5]={{0,0,0,2},{0,2,0,1},{0,3,0,1},{0,2,1,1},{0,3,1,1}};
    QuadField_T coarse=NULL,fine=NULL;
    QuadFieldReport cr={0},fr={0};
    MeshBinData mesh={0};
    int fail=QuadField_new(arena,coarse_cells,2,1,&coarse)!=0;
    fail+=QuadField_new(arena,fine_cells,5,1,&fine)!=0;
    for (int l=0;l<2 && !fail;l++) {
        QuadField_T f=l ? fine : coarse;
        const QuadFieldCell *cells=l ? fine_cells : coarse_cells;
        size_t nc=l ? 5 : 2;
        for (size_t c=0;c<nc;c++) for (int j=0;j<5;j++) for (int i=0;i<5;i++) {
            double u=cells[c].x+cells[c].size*i*.25,v=cells[c].y+cells[c].size*j*.25;
            double xyz[3]={4000+v,3000+u+2*v,2000+u*v};
            fail+=QuadField_observe(f,c,u,v,xyz,1)!=0;
        }
        fail+=QuadField_solve(f,l ? coarse : NULL,l ? &fr : &cr)!=0;
    }
    if (!fail) {
        double a[3]={0},b[3]={0};
        fail+=!(cr.cells==2 && fr.cells==5 && fr.hanging_nodes==1 && fr.unknowns<fr.nodes);
        for (int j=0;j<=20;j++) {
            double v=j*.1;
            fail+=QuadField_evaluate(fine,0,2,v,a)!=0;
            fail+=QuadField_evaluate(fine,v<=1 ? 1 : 3,2,v,b)!=0;
            for (int k=0;k<3;k++) fail+=fabs(a[k]-b[k])>1e-9;
            fail+=fabs(a[0]-(4000+v))>1e-5 || fabs(a[2]-(2000+2*v))>1e-5;
        }
        fail+=QuadField_mesh(arena,fine,&mesh)!=0;
        fail+=mesh.nv!=11 || mesh.nf!=11;
    }
    if (!fail) {
        /* A small source patch occupies only part of the control domain.
         * Keep its hole, two disconnected islands and an isolated point;
         * exclude a core triangle explicitly rather than inventing UVs. */
        MeshBinData source={0},projected={0},repeated={0};
        QuadFieldProjectionReport report={0},again={0};
        float uv[30]={.2f,.2f, 1.8f,.2f, 1.8f,1.8f, .2f,1.8f,
                       .6f,.6f, 1.4f,.6f, 1.4f,1.4f, .6f,1.4f,
                       2.2f,.2f, 3.8f,.2f, 2.2f,1.8f, 3.5f,1.8f,
                       0,0, 0,0, 0,0};
        float xyz[45]={0};
        int32_t tri[30]={0,1,5, 0,5,4, 1,2,6, 1,6,5,
                         2,3,7, 2,7,6, 3,0,4, 3,4,7, 8,9,10, 12,13,14};
        int32_t chart[15]={0},*original=NULL,*original_again=NULL;
        for (int i=0;i<15;i++) {
            double u=uv[2*i],v=uv[2*i+1];
            xyz[3*i]=(float)(4000+v); xyz[3*i+1]=(float)(3000+u+2*v);
            xyz[3*i+2]=(float)(2000+u*v);
            if (i>=12) chart[i]=-1;
        }
        source.nv=15; source.nf=10; source.verts=xyz; source.uv=uv; source.faces=tri;
        fail+=QuadField_project_mesh(arena,fine,&source,chart,&projected,&original,&report)!=0;
        fail+=QuadField_project_mesh(arena,fine,&source,chart,&repeated,&original_again,&again)!=0;
        fail+=!(projected.nv==12 && projected.nf==9 && report.excluded_vertices==3 && report.excluded_faces==1);
        if (projected.nv==12 && projected.nf==9 && repeated.nv==12 && repeated.nf==9) {
            fail+=memcmp(projected.uv,uv,24*sizeof(float))!=0 || memcmp(projected.faces,tri,27*sizeof(int32_t))!=0;
            fail+=memcmp(projected.verts,repeated.verts,36*sizeof(float))!=0;
            fail+=report.maximum_error>.001 || report.reversed_faces || report.collapsed_faces;
            for (int i=0;i<12;i++) fail+=original[i]!=i || original_again[i]!=i;
        }
        chart[0]=1;
        fail+=QuadField_project_mesh(arena,fine,&source,chart,&repeated,&original_again,&again)==0;
        chart[0]=0; tri[0]=15;
        fail+=QuadField_project_mesh(arena,fine,&source,chart,&repeated,&original_again,&again)==0;
        source=(MeshBinData){0};
        fail+=QuadField_project_mesh(arena,fine,&source,NULL,&repeated,&original_again,&again)!=0;
    }
    {
        QuadFieldCell bad[2]={{0,0,0,4},{0,4,0,1}};
        QuadField_T rejected=NULL;
        fail+=QuadField_new(arena,bad,2,1,&rejected)==0;
    }
    {
        QuadFieldCell duplicate[2]={{0,0,0,2},{0,0,0,2}};
        QuadFieldCell overlap[2]={{0,0,0,2},{0,1,1,1}};
        QuadFieldCell separate[2]={{0,-2,-2,2},{1,-2,-2,2}};
        QuadField_T f=NULL;
        size_t found=0;
        fail+=QuadField_new(arena,duplicate,2,1,&f)==0;
        fail+=QuadField_new(arena,overlap,2,1,&f)==0;
        fail+=QuadField_new(arena,separate,2,1,&f)!=0;
        if (f) {
            fail+=f->nn!=8;
            fail+=QuadField_locate(f,0,-1,-1,&found)!=0 || found!=0;
            fail+=QuadField_locate(f,1,-1,-1,&found)!=0 || found!=1;
            fail+=QuadField_locate(f,2,-1,-1,&found)==0;
        }
    }
    {
        /* Genuine nested solves: all levels see identical observations;
         * finer coefficients reduce a non-bilinear geometric residual. */
        QuadFieldCell cell={0,0,0,4};
        QuadField_T levels[3]={NULL,NULL,NULL};
        QuadFieldReport reports[3]={{0}};
        fail+=QuadField_new(arena,&cell,1,1,levels)!=0;
        for (int l=0;l<3 && levels[l];l++) {
            for (int y=0;y<=32;y++) for (int x=0;x<=32;x++) {
                double u=x*.125,v=y*.125,xyz[3]={4352+v,3072+u,2560+.3*u*u+.2*v*v};
                size_t c=0;
                if (QuadField_locate(levels[l],0,u,v,&c)!=0) { fail++; continue; }
                fail+=QuadField_observe(levels[l],c,u,v,xyz,1)!=0;
            }
            fail+=QuadField_solve(levels[l],l ? levels[l-1] : NULL,reports+l)!=0;
            if (l<2) {
                uint8_t mask[16];
                memset(mask,1,sizeof mask);
                fail+=QuadField_refine(arena,levels[l],mask,levels+l+1)!=0;
            }
        }
        fail+=!(reports[0].cells==1 && reports[1].cells==4 && reports[2].cells==16);
        fail+=!(reports[1].rms_error<.3*reports[0].rms_error && reports[2].rms_error<.3*reports[1].rms_error);
    }
    {
        /* Exact original-boundary data, a globally free interior, and a
         * constrained midpoint that cannot be fixed independently. */
        QuadField_T f=NULL;
        QuadFieldReport report={0};
        size_t node=0;
        double fixed[3]={0,0,0},bad[3]={0,0,1},a[3]={0},b[3]={0};
        fail+=QuadField_new(arena,fine_cells,5,1,&f)!=0;
        if (f) {
            fail+=QuadField_find_node(f,0,0,0,&node)!=0;
            fail+=QuadField_fix_node(f,node,fixed)!=0;
            fail+=QuadField_fix_node(f,node,fixed)!=0;
            fail+=QuadField_fix_node(f,node,bad)==0;
            fail+=QuadField_find_node(f,0,0,2,&node)!=0;
            fixed[1]=2;
            fail+=QuadField_fix_node(f,node,fixed)!=0;
            fail+=QuadField_find_node(f,0,2,1,&node)!=0;
            fail+=QuadField_fix_node(f,node,fixed)==0;
            fail+=QuadField_find_node(f,99,0,0,&node)!=1;
            for (size_t c=0;c<5;c++) for (int y=0;y<5;y++) for (int x=0;x<5;x++) {
                double u=fine_cells[c].x+fine_cells[c].size*x*.25;
                double v=fine_cells[c].y+fine_cells[c].size*y*.25;
                double xyz[3]={u,v,1};
                fail+=QuadField_observe(f,c,u,v,xyz,1)!=0;
            }
            fail+=QuadField_solve(f,NULL,&report)!=0;
            fail+=report.fixed_nodes!=2 || report.unknowns!=8;
            fail+=QuadField_evaluate(f,0,0,1,a)!=0;
            fail+=fabs(a[0])>1e-12 || fabs(a[1]-1)>1e-12 || fabs(a[2])>1e-12;
            fail+=QuadField_evaluate(f,0,2,1,a)!=0;
            fail+=QuadField_evaluate(f,1,2,1,b)!=0;
            fail+=fabs(a[2]-b[2])>1e-10 || a[2]<.5;
        }
    }
    {
        QuadFieldCell cell={0,0,0,1};
        QuadField_T f=NULL;
        QuadFieldReport report={0};
        double sum=0;
        fail+=QuadField_new(arena,&cell,1,1,&f)!=0;
        if (f) {
            for (int y=0;y<=4;y++) for (int x=0;x<=4;x++) {
                double u=x*.25,v=y*.25,target[3]={1e6+u,2e6+v,3e6+u*u};
                fail+=QuadField_observe(f,0,u,v,target,1)!=0;
            }
            fail+=QuadField_solve(f,NULL,&report)!=0;
            for (int y=0;y<=4;y++) for (int x=0;x<=4;x++) {
                double u=x*.25,v=y*.25,xyz[3]={0},target[3]={1e6+u,2e6+v,3e6+u*u};
                fail+=QuadField_evaluate(f,0,u,v,xyz)!=0;
                for (int k=0;k<3;k++) sum+=(xyz[k]-target[k])*(xyz[k]-target[k]);
            }
            fail+=fabs(report.rms_error-sqrt(sum/25))>1e-9;
            fail+=fabs(report.rms_error-sqrt(.0109375))>1e-8;
        }
    }
    for (int kind=1;kind<=2;kind++) {
        /* A diagonal crease is exactly P1, but not Q1, even with identical
         * corner samples. Exercise both diagonals and persisted evaluation. */
        QuadFieldCell cell={0,0,0,1};
        QuadField_T f=NULL;
        QuadFieldReport report={0};
        MeshBinData detail_mesh={0};
        double xyz[3]={0};
        fail+=QuadField_new(arena,&cell,1,1,&f)!=0;
        if (f) {
            fail+=QuadField_set_element(f,0,kind)!=0;
            for (int y=0;y<=4;y++) for (int x=0;x<=4;x++) {
                double u=x*.25,v=y*.25;
                double target[3]={u,v,3*fabs(kind==1 ? u-v : u+v-1)};
                fail+=QuadField_observe(f,0,u,v,target,1)!=0;
            }
            fail+=QuadField_set_element(f,0,0)==0;
            fail+=QuadField_solve(f,NULL,&report)!=0;
            fail+=report.unknowns!=4 || report.fixed_nodes!=0;
            fail+=QuadField_evaluate(f,0,.5,.5,xyz)!=0 || fabs(xyz[2])>1e-8;
            fail+=qf_test_persistence(arena,f)!=0;
            {
                QuadField_T copy=NULL;
                QuadFieldReport copied={0};
                fail+=QuadField_new(arena,&cell,1,1,&copy)!=0;
                if (copy) {
                    fail+=QuadField_set_element(copy,0,kind)!=0;
                    fail+=QuadField_reuse_observations(copy,f)!=0;
                    fail+=QuadField_observation_count(copy,0)!=25;
                    fail+=QuadField_reuse_observations(copy,f)==0;
                    fail+=QuadField_solve(copy,NULL,&copied)!=0;
                    fail+=memcmp(copy->value,f->value,12*sizeof(double))!=0;
                }
            }
            fail+=QuadField_mesh(arena,f,&detail_mesh)!=0 || detail_mesh.nf!=2;
            if (detail_mesh.nf==2) {
                int32_t a=detail_mesh.faces[0],b=detail_mesh.faces[1],c=detail_mesh.faces[2];
                double mean=(detail_mesh.verts[3*(size_t)a+2]+detail_mesh.verts[3*(size_t)b+2]+detail_mesh.verts[3*(size_t)c+2])/3;
                double u=(detail_mesh.uv[2*(size_t)a]+detail_mesh.uv[2*(size_t)b]+detail_mesh.uv[2*(size_t)c])/3.0;
                double v=(detail_mesh.uv[2*(size_t)a+1]+detail_mesh.uv[2*(size_t)b+1]+detail_mesh.uv[2*(size_t)c+1])/3.0;
                fail+=fabs(mean-3*fabs(kind==1 ? u-v : u+v-1))>1e-7;
            }
        }
    }
    {
        QuadField_T a=NULL,b=NULL;
        QfStreamTest ta={0},tb={0};
        MeshBinData ma={0},mb={0};
        fail+=QuadField_fit_stream(arena,2,qf_test_reader,NULL,4,1,.06,qf_test_checkpoint,&ta,&a)!=0;
        fail+=QuadField_fit_stream(arena,2,qf_test_reader,NULL,4,1,.06,qf_test_checkpoint,&tb,&b)!=0;
        fail+=ta.checkpoints!=3 || tb.checkpoints!=3;
        fail+=!(ta.quads[0]==1 && ta.quads[1]==4 && ta.quads[2]==16);
        if (a && b) {
            fail+=qf_test_persistence(arena,a)!=0;
            fail+=QuadField_mesh(arena,a,&ma)!=0 || QuadField_mesh(arena,b,&mb)!=0;
            if (ma.nv==mb.nv && ma.nf==mb.nf) {
                fail+=memcmp(ma.verts,mb.verts,3*ma.nv*sizeof(float))!=0;
                fail+=memcmp(ma.uv,mb.uv,2*ma.nv*sizeof(float))!=0;
                fail+=memcmp(ma.faces,mb.faces,3*ma.nf*sizeof(int32_t))!=0;
            } else fail++;
        }
    }
    fail+=qf_test_many_known_values();
    fprintf(stderr,"[selftest] constrained quad field %s (%d failures)\n",fail ? "FAIL" : "PASS",fail);
    Arena_dispose(&arena);
    return fail ? -1 : 0;
}
