/* c3_coarse_coordinates.coordinate_basis and source_gauge_nullspace.
 * Only the incremental coordinate space changes; all original vertices,
 * faces, physical witnesses and collision context remain in the problem. */
typedef struct ArsWeight {int column;double value;} ArsWeight;
struct ArsCoordinates {size_t *offset,n,capacity;ArsWeight *weights;};
typedef struct ArsGridNode {int64_t x,y;size_t slot;} ArsGridNode;
#define ARS_STENCIL_MAX 48
typedef struct ArsStencil {int n,id[ARS_STENCIL_MAX];double x[ARS_STENCIL_MAX],y[ARS_STENCIL_MAX];} ArsStencil;

static int ars_stencil_entry(ArsStencil *s,int id,double x,double y)
{
    if(id<0 || (!x && !y))return 1;
    int at=0;while(at<s->n && s->id[at]!=id)at++;
    if(at==s->n){if(s->n==ARS_STENCIL_MAX)return 0;s->id[at]=id;s->x[at]=s->y[at]=0;s->n++;}
    s->x[at]+=x;s->y[at]+=y;return isfinite(s->x[at]) && isfinite(s->y[at]);
}
static int ars_stencil_add(const ArsProblem *p,ArsStencil *s,size_t vertex,double x,double y)
{
    if(!p->coordinates)return ars_stencil_entry(s,p->id[vertex],x,y);
    const ArsCoordinates *c=p->coordinates;
    for(size_t k=c->offset[vertex];k<c->offset[vertex+1];k++)
        if(!ars_stencil_entry(s,c->weights[k].column,c->weights[k].value*x,c->weights[k].value*y))return 0;
    return 1;
}
static int ars_coordinate_append(ArsCoordinates *c,int column,double value)
{
    if(!value)return 1;
    if(!isfinite(value) || column<0)return 0;
    if(c->n==c->capacity){
        if(c->capacity>SIZE_MAX/(2*sizeof(ArsWeight)))return 0;
        size_t capacity=c->capacity?2*c->capacity:1024;
        ArsWeight *next=realloc(c->weights,capacity*sizeof *next);if(!next)return 0;
        c->weights=next;c->capacity=capacity;
    }
    c->weights[c->n++]=(ArsWeight){column,value};return 1;
}
static int ars_grid_order(const void *a,const void *b)
{
    const ArsGridNode *x=a,*y=b;if(x->x!=y->x)return(x->x>y->x)-(x->x<y->x);
    return(x->y>y->y)-(x->y<y->y);
}
static int ars_int_order(const void *a,const void *b)
{int x=*(const int *)a,y=*(const int *)b;return(x>y)-(x<y);}

/* A moving native chart has at most the two hard source-gauge points.
 * Different gauge charts have disjoint control support, so their pivoted
 * nullspace blocks can be eliminated independently without a dense global
 * matrix containing all fixed-context vertices. */
static int ars_gauge_pivots(int pins,int count,const double g[2][8],int pivot[2],double inverse[4])
{
    pivot[0]=pivot[1]=-1;memset(inverse,0,4*sizeof(double));if(!pins)return 1;
    double largest=0;
    for(int c=0;c<count;c++){
        double norm=g[0][c]*g[0][c]+(pins==2?g[1][c]*g[1][c]:0);
        if(norm>largest){largest=norm;pivot[0]=c;}
    }
    if(sqrt(largest)<1e-12)return 0;
    if(pins==1){inverse[0]=1/g[0][pivot[0]];return 1;}
    double best=0;
    for(int c=0;c<count;c++)if(c!=pivot[0]){
        double residual=fabs(g[0][pivot[0]]*g[1][c]-g[1][pivot[0]]*g[0][c])/sqrt(largest);
        if(residual>best){best=residual;pivot[1]=c;}
    }
    if(best<1e-12)return 0;
    double a=g[0][pivot[0]],b=g[0][pivot[1]],c=g[1][pivot[0]],d=g[1][pivot[1]],det=a*d-b*c;
    inverse[0]=d/det;inverse[1]=-b/det;inverse[2]=-c/det;inverse[3]=a/det;return 1;
}
static int ars_coordinate_chart(ArsProblem *p,size_t chart,double step)
{
    const AsmField *f=p->s->field;ArsCoordinates *coordinates=p->coordinates;
    size_t first=f->vertex_offset[chart],nv=f->run->charts[chart].nv;
    if(p->s->opts.fixed_charts && p->s->opts.fixed_charts[chart]){
        for(size_t v=first;v<first+nv;v++)coordinates->offset[v+1]=coordinates->n;
        return 1;
    }
    if(nv>SIZE_MAX/(4*sizeof(ArsGridNode)))return 0;
    Arena_T arena=Arena_new();int ok=0;
    ArsGridNode *nodes=ARENA_ALLOC(arena,4*nv*sizeof *nodes);
    int *columns=ARENA_ALLOC(arena,4*nv*sizeof(int));double *weight=ARENA_ALLOC(arena,4*nv*sizeof(double));
    for(size_t v=0;v<nv;v++){
        double qx=f->uv[2*(first+v)]/step,qy=f->uv[2*(first+v)+1]/step;
        if(!isfinite(qx) || !isfinite(qy) || fabs(qx)>4503599627370495. || fabs(qy)>4503599627370495.)goto done;
        double x=floor(qx),y=floor(qy),tx=qx-x,ty=qy-y;
        weight[4*v]=(1-tx)*(1-ty);weight[4*v+1]=tx*(1-ty);weight[4*v+2]=(1-tx)*ty;weight[4*v+3]=tx*ty;
        double sum=0;for(int k=0;k<4;k++){
            nodes[4*v+k]=(ArsGridNode){(int64_t)x+(k&1),(int64_t)y+(k>>1),4*v+(size_t)k};sum+=weight[4*v+k];
        }
        if(fabs(sum-1)>1e-14)goto done;
    }
    qsort(nodes,4*nv,sizeof *nodes,ars_grid_order);int unique=0;int64_t oldx=0,oldy=0;
    for(size_t k=0;k<4*nv;k++){
        ArsGridNode node=nodes[k];
        if(!k || node.x!=oldx || node.y!=oldy){if(unique==INT_MAX)goto done;unique++;oldx=node.x;oldy=node.y;}
        columns[node.slot]=unique-1;
    }
    if((size_t)unique>=nv){
        for(size_t v=first;v<first+nv;v++){
            if(!p->s->fixed[v]){if(p->n>INT_MAX-2 || !ars_coordinate_append(coordinates,p->n,1))goto done;p->n+=2;}
            coordinates->offset[v+1]=coordinates->n;
        }
        ok=1;goto done;
    }
    size_t pin_vertices[2];int pins=0,support[8],count=0;
    for(size_t v=0;v<nv;v++)if(p->s->fixed[first+v]){if(pins==2)goto done;pin_vertices[pins++]=v;}
    for(int row=0;row<pins;row++)for(int k=0;k<4;k++){
        size_t at=4*pin_vertices[row]+k;if(!weight[at])continue;int column=columns[at],found=0;
        for(int j=0;j<count;j++)found|=support[j]==column;
        if(!found)support[count++]=column;
    }
    qsort(support,(size_t)count,sizeof(int),ars_int_order);double gauge[2][8]={{0}};
    for(int row=0;row<pins;row++)for(int k=0;k<4;k++){
        size_t at=4*pin_vertices[row]+k;
        for(int j=0;j<count;j++)if(support[j]==columns[at])gauge[row][j]+=weight[at];
    }
    int pivots[2];double inverse[4];if(!ars_gauge_pivots(pins,count,gauge,pivots,inverse))goto done;
    int raw_pivots[]={pins?support[pivots[0]]:-1,pins==2?support[pivots[1]]:-1};
    int *free_id=ARENA_ALLOC(arena,(size_t)unique*sizeof(int));
    for(int k=0;k<unique;k++){
        free_id[k]=-1;if(k==raw_pivots[0] || k==raw_pivots[1])continue;
        if(p->n>INT_MAX-2)goto done;free_id[k]=p->n;p->n+=2;
    }
    double correction[2][8]={{0}};
    for(int row=0;row<pins;row++)for(int j=0;j<count;j++)
        for(int k=0;k<pins;k++)correction[row][j]-=inverse[2*row+k]*gauge[k][j];
    for(size_t v=0;v<nv;v++){
        ArsStencil row={0};
        for(int k=0;k<4;k++){
            size_t at=4*v+k;int column=columns[at];double w=weight[at];if(!w)continue;
            if(free_id[column]>=0){if(!ars_stencil_entry(&row,free_id[column],w,0))goto done;}
            else{
                int pivot=column==raw_pivots[0]?0:1;
                for(int j=0;j<count;j++)if(free_id[support[j]]>=0 &&
                    !ars_stencil_entry(&row,free_id[support[j]],w*correction[pivot][j],0))goto done;
            }
        }
        if(p->s->fixed[first+v]){
            for(int k=0;k<row.n;k++)if(fabs(row.x[k])>1e-12)goto done;
        }else{
            for(int k=0;k<row.n;k++)if(!ars_coordinate_append(coordinates,row.id[k],row.x[k]))goto done;
        }
        coordinates->offset[first+v+1]=coordinates->n;
    }
    ok=1;
done:
    Arena_dispose(&arena);return ok;
}
static int ars_coordinates_new(ArsProblem *p,double step)
{
    const AsmField *f=p->s->field;if(!isfinite(step) || step<0)return 0;
    p->id=ARENA_ALLOC(p->arena,f->nv*sizeof(int));
    if(!step){for(size_t v=0;v<f->nv;v++){p->id[v]=p->s->fixed[v]?-1:p->n;if(!p->s->fixed[v])p->n+=2;}return 1;}
    p->coordinates=ARENA_CALLOC(p->arena,1,sizeof(ArsCoordinates));
    p->coordinates->offset=ARENA_CALLOC(p->arena,f->nv+1,sizeof(size_t));
    for(size_t c=0;c<f->run->n_charts;c++)if(f->vertex_offset[c]!=SIZE_MAX && !ars_coordinate_chart(p,c,step))return 0;
    for(size_t v=0;v<f->nv;v++)p->id[v]=p->coordinates->offset[v]<p->coordinates->offset[v+1]?0:-1;
    return 1;
}
static void ars_coordinate_direction(const ArsProblem *p,double *direction)
{
    const AsmField *f=p->s->field;memset(direction,0,2*f->nv*sizeof(double));
    for(size_t v=0;v<f->nv;v++)if(p->id[v]>=0){
        if(!p->coordinates){for(int d=0;d<2;d++)direction[2*v+d]=p->answer[p->id[v]+d];}
        else for(size_t k=p->coordinates->offset[v];k<p->coordinates->offset[v+1];k++)
            for(int d=0;d<2;d++)direction[2*v+d]+=p->coordinates->weights[k].value*p->answer[p->coordinates->weights[k].column+d];
    }
}
static size_t ars_hash_index(uint64_t key,size_t mask)
{
    key^=key>>30;key*=UINT64_C(0xbf58476d1ce4e5b9);key^=key>>27;key*=UINT64_C(0x94d049bb133111eb);key^=key>>31;
    return(size_t)key&mask;
}
static void ars_hashed_entry(ArsProblem *p,int row,int col,double value)
{
    if(p->h.failed || !value)return;
    if(!isfinite(value)){p->h.failed=1;return;}
    if((size_t)p->h.n*2>=p->hcapacity){
        if(p->hcapacity>SIZE_MAX/(2*sizeof(ArsHash))){p->h.failed=1;return;}
        size_t capacity=p->hcapacity?2*p->hcapacity:2048;
        ArsHash *next=calloc(capacity,sizeof *next);if(!next){p->h.failed=1;return;}
        for(size_t i=0;i<p->hcapacity;i++)if(p->hslots[i].key){
            size_t at=ars_hash_index(p->hslots[i].key,capacity-1);while(next[at].key)at=(at+1)&(capacity-1);
            next[at]=p->hslots[i];
        }
        free(p->hslots);p->hslots=next;p->hcapacity=capacity;
    }
    uint64_t key=(((uint64_t)(uint32_t)row<<32)|(uint32_t)col)+1;
    size_t at=ars_hash_index(key,p->hcapacity-1);while(p->hslots[at].key && p->hslots[at].key!=key)at=(at+1)&(p->hcapacity-1);
    if(p->hslots[at].key){int index=p->hslots[at].index;p->h.value[index]+=value;if(!isfinite(p->h.value[index]))p->h.failed=1;}
    else{int index=p->h.n;alq_entry(&p->h,row,col,value);if(!p->h.failed)p->hslots[at]=(ArsHash){key,index};}
}
