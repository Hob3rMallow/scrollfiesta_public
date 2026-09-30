/* Original-contact supporting axes and source-surface neighborhoods for the
 * explicit translation and local-vertex clearance QPs below. */
typedef struct AclPair { int32_t a,b;double normal[2],slack,weight;size_t face[2]; } AclPair;
typedef struct AclQuery { AprSolver *solver;AclPair *pairs;size_t count;int local; } AclQuery;
static void acl_pair(void *context,size_t a,size_t b,double area)
{
    (void)area;AclQuery *q=context;AprSolver *s=q->solver;const AsmField *f=s->field;
    if(a>b){size_t t=a;a=b;b=t;}
    int32_t ca=f->face_chart[a],cb=f->face_chart[b];
    if(!q->local && (ca==cb || (s->rigid_fixed[ca] && s->rigid_fixed[cb])))return;
    double origin[2]={0},tri[2][3][2],best=-DBL_MAX,normal[2]={0};
    for(int side=0;side<2;side++)for(int k=0;k<3;k++)for(int d=0;d<2;d++)origin[d]+=f->uv[2*(size_t)f->faces[3*(side?b:a)+k]+d];
    origin[0]/=6;origin[1]/=6;
    for(int side=0;side<2;side++)for(int k=0;k<3;k++)for(int d=0;d<2;d++)tri[side][k][d]=f->uv[2*(size_t)f->faces[3*(side?b:a)+k]+d]-origin[d];
    for(int sign=1;sign>=-1;sign-=2)for(int side=0;side<2;side++)for(int k=0;k<3;k++){
        double nx=tri[side][(k+1)%3][1]-tri[side][k][1],ny=tri[side][k][0]-tri[side][(k+1)%3][0],length=hypot(nx,ny);
        if(!(length>0))continue;nx*=sign/length;ny*=sign/length;
        double high=-DBL_MAX,low=DBL_MAX;
        for(int j=0;j<3;j++){high=fmax(high,nx*tri[0][j][0]+ny*tri[0][j][1]);low=fmin(low,nx*tri[1][j][0]+ny*tri[1][j][1]);}
        if(low-high>best){best=low-high;normal[0]=nx;normal[1]=ny;}
    }
    if(!isfinite(best))return;
    AclPair p={ca,cb,{normal[0],normal[1]},fmax(0,1e-4-best),0,{a,b}};
    if(p.a>p.b){int32_t t=p.a;p.a=p.b;p.b=t;p.normal[0]=-p.normal[0];p.normal[1]=-p.normal[1];p.face[0]=b;p.face[1]=a;}
    q->pairs[q->count++]=p;
}
static int acl_order(const void *a,const void *b)
{
    const AclPair *x=a,*y=b;if(x->a!=y->a)return (x->a>y->a)-(x->a<y->a);return (x->b>y->b)-(x->b<y->b);
}
static int acl_face_order(const void *a,const void *b)
{
    const AclPair *x=a,*y=b;
    if(x->face[0]!=y->face[0])return(x->face[0]>y->face[0])-(x->face[0]<y->face[0]);
    return(x->face[1]>y->face[1])-(x->face[1]<y->face[1]);
}
#include "asm_repair_local_qp.inc"

static int apr_supporting_clearance(AprSolver *s,double penalty)
{ return !apr_expired(&s->opts) && apr_clearance_qp(s,NULL,0,penalty); }

typedef struct AclSeeds { const AsmField *field;double *distance;const uint8_t *fixed_charts; } AclSeeds;
static void acl_seed(void *context,size_t a,size_t b,double area)
{
    (void)area;AclSeeds *s=context;for(int side=0;side<2;side++)for(int k=0;k<3;k++){
        size_t v=(size_t)s->field->faces[3*(side?b:a)+k];
        if(!s->fixed_charts || !s->fixed_charts[s->field->vertex_chart[v]])s->distance[v]=0;
    }
}
static void acl_heap_up(size_t *heap,int32_t *position,size_t at,const double *distance)
{
    size_t v=heap[at];while(at){size_t parent=(at-1)/2;if(distance[heap[parent]]<=distance[v])break;heap[at]=heap[parent];position[heap[at]]=(int32_t)at;at=parent;}heap[at]=v;position[v]=(int32_t)at;
}
static size_t acl_heap_pop(size_t *heap,int32_t *position,size_t *count,const double *distance)
{
    size_t result=heap[0],v=heap[--*count],at=0;position[result]=-2;
    if(!*count)return result;
    while(2*at+1<*count){size_t child=2*at+1;if(child+1<*count && distance[heap[child+1]]<distance[heap[child]])child++;
        if(distance[v]<=distance[heap[child]])break;heap[at]=heap[child];position[heap[at]]=(int32_t)at;at=child;}
    heap[at]=v;position[v]=(int32_t)at;return result;
}
static int apr_local_clearance(AprSolver *s,double radius)
{
    if(apr_expired(&s->opts))return 0;
    Arena_T arena=Arena_new();AsmField *f=s->field;int accepted=0;
    double *distance=ARENA_ALLOC(arena,f->nv*sizeof(double));for(size_t v=0;v<f->nv;v++)distance[v]=DBL_MAX;
    AclSeeds seed={f,distance,s->opts.fixed_charts};AsmContactStats contacts;
    if(AsmContacts_measure(f->contacts,f->uv,NULL,acl_seed,&seed,&contacts) || !contacts.pairs)goto done;
    size_t *offset=ARENA_CALLOC(arena,f->nv+1,sizeof(size_t));
    for(size_t k=0;k<s->n_active;k++)for(int q=0;q<3;q++)offset[(size_t)f->faces[3*s->active_faces[k]+q]+1]++;
    for(size_t v=0;v<f->nv;v++)offset[v+1]+=offset[v];
    size_t *incident=ARENA_ALLOC(arena,offset[f->nv]*sizeof(size_t)),*fill=ARENA_CALLOC(arena,f->nv,sizeof(size_t));
    for(size_t k=0;k<s->n_active;k++)for(int q=0;q<3;q++){size_t v=(size_t)f->faces[3*s->active_faces[k]+q];incident[offset[v]+fill[v]++]=s->active_faces[k];}
    size_t *heap=ARENA_ALLOC(arena,f->nv*sizeof(size_t)),count=0;int32_t *position=ARENA_ALLOC(arena,f->nv*sizeof(int32_t));
    for(size_t v=0;v<f->nv;v++){position[v]=-1;if(distance[v]==0){heap[count]=v;position[v]=(int32_t)count++;}}
    while(count){size_t v=acl_heap_pop(heap,position,&count,distance);
        for(size_t k=offset[v];k<offset[v+1];k++)for(int q=0;q<3;q++){
            size_t w=(size_t)f->faces[3*incident[k]+q];if(w==v || position[w]==-2)continue;
            double length=0;for(int d=0;d<3;d++){double x=(double)f->xyz[3*w+d]-f->xyz[3*v+d];length+=x*x;}
            double next=distance[v]+sqrt(length);if(next>radius || next>=distance[w])continue;distance[w]=next;
            if(position[w]<0){position[w]=(int32_t)count;heap[count++]=w;}acl_heap_up(heap,position,(size_t)position[w],distance);
        }
    }
    accepted=apr_clearance_qp(s,distance,radius,1e4);
    s->energy=apr_objective(s,f->uv,1,&s->contacts);
done:
    Arena_dispose(&arena);return accepted;
}
