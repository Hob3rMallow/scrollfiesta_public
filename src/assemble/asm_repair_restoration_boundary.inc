/* All original boundary vertices and edges, with double-precision AABB queries.
 * Streaming candidates avoids a sampled shell or a truncated contact list.
 * Static forces and swept contact use the same source-distance rule as Python. */
typedef struct ArsEdge {int32_t a,b;} ArsEdge;
typedef struct ArsNode {double lo[2],hi[2];size_t begin,end;int left,right,active;} ArsNode;
struct ArsBoundary {
    ArsProblem *problem;ArsEdge *edges;size_t ne,capacity,nv,nodes;
    size_t *parent,*vertices,*order;uint8_t *moving;ArsNode *tree;double *box;
    const double *uv,*direction;double padding,limit;int failed;
};
static int ars_edge_order(const void *a,const void *b)
{
    const ArsEdge *x=a,*y=b;if(x->a!=y->a)return(x->a>y->a)-(x->a<y->a);
    return(x->b>y->b)-(x->b<y->b);
}
static int ars_edge_append(ArsBoundary *b,ArsEdge edge)
{
    if(b->ne==b->capacity){
        if(b->capacity>SIZE_MAX/(2*sizeof(ArsEdge)))return 0;
        size_t capacity=b->capacity?2*b->capacity:1024;
        ArsEdge *next=realloc(b->edges,capacity*sizeof(ArsEdge));if(!next)return 0;
        b->edges=next;b->capacity=capacity;
    }
    b->edges[b->ne++]=edge;return 1;
}
static void ars_bnode_box(ArsBoundary *b,ArsNode *node)
{
    node->lo[0]=node->lo[1]=DBL_MAX;node->hi[0]=node->hi[1]=-DBL_MAX;node->active=0;
    for(size_t at=node->begin;at<node->end;at++){
        size_t edge=b->order[at];const double *box=b->box+4*edge;const ArsEdge *e=b->edges+edge;
        for(int d=0;d<2;d++){node->lo[d]=fmin(node->lo[d],box[d]);node->hi[d]=fmax(node->hi[d],box[2+d]);}
        node->active|=b->moving[e->a] || b->moving[e->b];
    }
}
static int ars_bless(const ArsBoundary *b,size_t a,size_t c,int axis)
{
    double x=b->box[4*a+axis]+.5*(b->box[4*a+2+axis]-b->box[4*a+axis]);
    double y=b->box[4*c+axis]+.5*(b->box[4*c+2+axis]-b->box[4*c+axis]);
    return x<y || (x==y && a<c);
}
static void ars_bselect(ArsBoundary *b,size_t first,size_t last,size_t nth,int axis)
{
    while(last-first>1){
        size_t pivot=b->order[first+(last-first)/2],i=first,j=last-1;
        for(;;){
            while(ars_bless(b,b->order[i],pivot,axis))i++;
            while(ars_bless(b,pivot,b->order[j],axis))j--;
            if(i>=j)break;size_t tmp=b->order[i];b->order[i++]=b->order[j];b->order[j--]=tmp;
        }
        size_t cut=i;if(cut==first)cut++;if(nth<cut)last=cut;else first=cut;
    }
}
static int ars_bbuild(ArsBoundary *b,size_t first,size_t last)
{
    int id=(int)b->nodes++;ArsNode *node=b->tree+id;node->begin=first;node->end=last;node->left=node->right=-1;
    ars_bnode_box(b,node);
    if(last-first>8){
        size_t mid=first+(last-first)/2;int axis=node->hi[1]-node->lo[1]>node->hi[0]-node->lo[0];
        ars_bselect(b,first,last,mid,axis);node->left=ars_bbuild(b,first,mid);node->right=ars_bbuild(b,mid,last);
    }
    return id;
}
static int ars_bbox_update(ArsBoundary *b,const double *uv,const double *direction)
{
    const AsmField *f=b->problem->s->field;b->uv=uv;b->direction=direction;b->padding=1;
    for(size_t v=0;v<2*f->nv;v++){
        if(!isfinite(uv[v]) || (direction && !isfinite(direction[v])))return 0;
        b->padding=fmax(b->padding,fabs(uv[v]));
        if(direction)b->padding=fmax(b->padding,fmax(fabs(direction[v]),fabs(uv[v]+direction[v])));
    }
    b->padding*=128*DBL_EPSILON;
    for(size_t i=0;i<b->ne;i++){
        const ArsEdge *e=b->edges+i;double *box=b->box+4*i;
        if(hypot(uv[2*(size_t)e->b]-uv[2*(size_t)e->a],uv[2*(size_t)e->b+1]-uv[2*(size_t)e->a+1])<=1e-12)return 0;
        for(int d=0;d<2;d++){
            double a=uv[2*(size_t)e->a+d],c=uv[2*(size_t)e->b+d];box[d]=fmin(a,c);box[2+d]=fmax(a,c);
            if(direction){a+=direction[2*(size_t)e->a+d];c+=direction[2*(size_t)e->b+d];box[d]=fmin(box[d],fmin(a,c));box[2+d]=fmax(box[2+d],fmax(a,c));}
        }
    }
    return 1;
}
static int ars_brefit(ArsBoundary *b,const double *uv,const double *direction)
{
    if(!ars_bbox_update(b,uv,direction))return 0;
    for(size_t i=b->nodes;i-->0;){
        ArsNode *node=b->tree+i;
        if(node->left<0)ars_bnode_box(b,node);
        else{
            const ArsNode *a=b->tree+node->left,*c=b->tree+node->right;node->active=a->active || c->active;
            for(int d=0;d<2;d++){node->lo[d]=fmin(a->lo[d],c->lo[d]);node->hi[d]=fmax(a->hi[d],c->hi[d]);}
        }
    }
    return 1;
}
static int ars_boundary_new(ArsProblem *p)
{
    AsmField *f=p->s->field;ArsBoundary *b=ARENA_CALLOC(p->arena,1,sizeof *b);p->boundary=b;b->problem=p;
    b->parent=ARENA_ALLOC(p->arena,f->nv*sizeof(size_t));b->moving=ARENA_ALLOC(p->arena,f->nv);
    uint8_t *boundary=ARENA_CALLOC(p->arena,f->nv,1);size_t maximum=0;
    for(size_t v=0;v<f->nv;v++){b->parent[v]=v;b->moving[v]=(uint8_t)(!p->s->opts.fixed_charts || !p->s->opts.fixed_charts[f->vertex_chart[v]]);}
    for(size_t c=0;c<f->run->n_charts;c++)if(f->face_offset[c]!=SIZE_MAX && f->run->charts[c].nf>maximum)maximum=f->run->charts[c].nf;
    if(maximum>SIZE_MAX/(3*sizeof(ArsEdge)))return 0;
    ArsEdge *work=malloc(3*maximum*sizeof(ArsEdge));if(!work)return 0;int ok=1;
    for(size_t c=0;c<f->run->n_charts && ok;c++)if(f->face_offset[c]!=SIZE_MAX){
        size_t count=0;
        for(size_t t=0;t<f->run->charts[c].nf;t++){
            const int32_t *face=f->faces+3*(f->face_offset[c]+t);
            for(int k=0;k<3;k++){
                int32_t a=face[k],d=face[(k+1)%3];work[count++]=(ArsEdge){a<d?a:d,a<d?d:a};
                size_t r=apr_root(b->parent,(size_t)a),s=apr_root(b->parent,(size_t)d);b->parent[r>s?r:s]=r<s?r:s;
            }
        }
        qsort(work,count,sizeof *work,ars_edge_order);
        for(size_t at=0;at<count;){
            size_t last=at+1;while(last<count && !ars_edge_order(work+at,work+last))last++;
            if(last==at+1){ok=ars_edge_append(b,work[at]);boundary[work[at].a]=boundary[work[at].b]=1;if(!ok)break;}at=last;
        }
    }
    free(work);if(!ok || b->ne>(size_t)INT_MAX/2 || b->ne>SIZE_MAX/(2*sizeof(ArsNode)))return 0;
    b->vertices=ARENA_ALLOC(p->arena,f->nv*sizeof(size_t));
    for(size_t v=0;v<f->nv;v++){b->parent[v]=apr_root(b->parent,v);if(boundary[v])b->vertices[b->nv++]=v;}
    b->order=ARENA_ALLOC(p->arena,(b->ne?b->ne:1)*sizeof(size_t));
    b->box=ARENA_ALLOC(p->arena,4*(b->ne?b->ne:1)*sizeof(double));
    b->tree=ARENA_ALLOC(p->arena,2*(b->ne?b->ne:1)*sizeof(ArsNode));
    for(size_t i=0;i<b->ne;i++)b->order[i]=i;
    if(!ars_bbox_update(b,f->uv,NULL))return 0;
    if(b->ne)ars_bbuild(b,0,b->ne);return 1;
}
static int ars_bgeometry(const ArsBoundary *b,size_t v,const ArsEdge *e,double *t,double q[2],double *d2)
{
    double edge[2],point[2],length=0,dot=0;
    for(int d=0;d<2;d++){edge[d]=b->uv[2*(size_t)e->b+d]-b->uv[2*(size_t)e->a+d];point[d]=b->uv[2*v+d]-b->uv[2*(size_t)e->a+d];length+=edge[d]*edge[d];dot+=point[d]*edge[d];}
    if(!(length>1e-24))return 0;
    *t=fmax(0,fmin(1,dot/length));*d2=0;
    for(int d=0;d<2;d++){q[d]=point[d]-*t*edge[d];*d2+=q[d]*q[d];}
    return isfinite(*d2);
}
static void ars_bstatic(ArsBoundary *b,size_t v,const ArsEdge *e)
{
    ArsProblem *p=b->problem;const AsmField *f=p->s->field;double t,q[2],d2;
    if(!ars_bgeometry(b,v,e,&t,q,&d2) || !(d2>1e-24)){b->failed=1;return;}
    double height2=.25;
    if(b->parent[v]==b->parent[e->a]){
        double edge[3],point[3],length=0,dot=0;
        for(int d=0;d<3;d++){edge[d]=(double)f->xyz[3*(size_t)e->b+d]-f->xyz[3*(size_t)e->a+d];point[d]=(double)f->xyz[3*v+d]-f->xyz[3*(size_t)e->a+d];length+=edge[d]*edge[d];dot+=edge[d]*point[d];}
        if(!(length>0)){b->failed=1;return;}
        double parameter=fmax(0,fmin(1,dot/length)),distance=0;
        for(int d=0;d<3;d++){double value=point[d]-parameter*edge[d];distance+=value*value;}height2=fmin(height2,distance);
    }
    if(d2>=height2)return;
    double r=d2/height2;p->energy+=(1-r)*(1-r)/r;
    if(!p->derivatives)return;
    size_t vertices[]={v,(size_t)e->a,(size_t)e->b};double coefficient[]={1,t-1,-t};
    double factor=2/height2*(1-1/(r*r)),curvature=2/(r*r*r*height2*height2);ArsStencil stencil={0};
    for(int k=0;k<3;k++)if(!ars_stencil_add(p,&stencil,vertices[k],coefficient[k],0)){b->failed=1;return;}
    for(int a=0;a<stencil.n;a++)for(int d=0;d<2;d++){
        int row=stencil.id[a]+d;p->gradient[row]+=factor*stencil.x[a]*q[d];
        for(int c=0;c<stencil.n;c++)for(int e=0;e<2;e++)
            ars_h(p,row,stencil.id[c]+e,curvature*(2*stencil.x[a]*q[d])*(2*stencil.x[c]*q[e]));
    }
}
static void ars_bswept(ArsBoundary *b,size_t v,const ArsEdge *e)
{
    double edge[2],point[2],de[2],dp[2];
    for(int d=0;d<2;d++){
        edge[d]=b->uv[2*(size_t)e->b+d]-b->uv[2*(size_t)e->a+d];point[d]=b->uv[2*v+d]-b->uv[2*(size_t)e->a+d];
        de[d]=b->direction[2*(size_t)e->b+d]-b->direction[2*(size_t)e->a+d];dp[d]=b->direction[2*v+d]-b->direction[2*(size_t)e->a+d];
    }
    double a=ars_cross(de,dp),c=ars_cross(edge,point),middle=ars_cross(de,point)+ars_cross(edge,dp),roots[2];
    ars_roots(a,middle,c,roots);
    for(int k=0;k<2;k++)if(roots[k]>0 && roots[k]<=1 && isfinite(roots[k])){
        double t=roots[k],length=0,dot=0;
        for(int d=0;d<2;d++){double ed=edge[d]+t*de[d],po=point[d]+t*dp[d];length+=ed*ed;dot+=ed*po;}
        if(dot>=-1e-10 && dot<=length+1e-10)b->limit=fmin(b->limit,t);
    }
    double scale=fmax(1,fmax(fabs(a),fmax(fabs(middle),fabs(c))));
    if(fmax(fabs(a),fmax(fabs(middle),fabs(c)))<=32*DBL_EPSILON*scale){
        double t,q[2],d2,speed=fmax(hypot(dp[0],dp[1]),hypot(dp[0]-de[0],dp[1]-de[1]));
        if(!ars_bgeometry(b,v,e,&t,q,&d2)){b->failed=1;return;}
        if(speed>0)b->limit=fmin(b->limit,sqrt(d2)/speed);
    }
}
static void ars_bquery(ArsBoundary *b,int node_id,size_t v,const double lo[2],const double hi[2],int swept)
{
    if(b->failed)return;const ArsNode *node=b->tree+node_id;
    if((!b->moving[v] && !node->active) || lo[0]>node->hi[0] || lo[1]>node->hi[1] || hi[0]<node->lo[0] || hi[1]<node->lo[1])return;
    if(node->left>=0){ars_bquery(b,node->left,v,lo,hi,swept);ars_bquery(b,node->right,v,lo,hi,swept);return;}
    for(size_t at=node->begin;at<node->end;at++){
        size_t index=b->order[at];const ArsEdge *e=b->edges+index;const double *box=b->box+4*index;
        if(v==(size_t)e->a || v==(size_t)e->b || (!b->moving[v] && !b->moving[e->a] && !b->moving[e->b]))continue;
        if(lo[0]>box[2] || lo[1]>box[3] || hi[0]<box[0] || hi[1]<box[1])continue;
        if(swept)ars_bswept(b,v,e);else ars_bstatic(b,v,e);
        if(b->failed)return;
    }
}
static int ars_boundary_energy(ArsProblem *p,const double *uv)
{
    ArsBoundary *b=p->boundary;b->failed=0;if(!ars_brefit(b,uv,NULL))return 0;
    for(size_t at=0;at<b->nv && !b->failed;at++){
        size_t v=b->vertices[at];double lo[2],hi[2];
        for(int d=0;d<2;d++){lo[d]=uv[2*v+d]-.5-b->padding;hi[d]=uv[2*v+d]+.5+b->padding;}
        ars_bquery(b,0,v,lo,hi,0);
    }
    return !b->failed && isfinite(p->energy) && !p->h.failed;
}
static double ars_boundary_limit(ArsProblem *p,const double *uv,const double *direction)
{
    ArsBoundary *b=p->boundary;b->failed=0;b->limit=INFINITY;
    if(!ars_brefit(b,uv,direction))return 0;
    for(size_t at=0;at<b->nv && !b->failed;at++){
        size_t v=b->vertices[at];double lo[2],hi[2];
        for(int d=0;d<2;d++){double first=uv[2*v+d],last=first+direction[2*v+d];lo[d]=fmin(first,last)-b->padding;hi[d]=fmax(first,last)+b->padding;}
        ars_bquery(b,0,v,lo,hi,1);
    }
    return b->failed?0:b->limit;
}
