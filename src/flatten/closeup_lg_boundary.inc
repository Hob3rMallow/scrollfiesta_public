/* Private relative boundary constraints. Preserve the current separating-axis
 * clearance of near edge pairs through the entire linear displacement. No
 * vertex is pinned and no contact/physical acceptance test is changed. */
typedef struct SmRimBox { int32_t a,b; double lo[2],hi[2]; } SmRimBox;
typedef struct SmRimLearn SmRimLearn;
typedef struct SmRim {
    size_t nv,nb,calls,pairs,rows,constant_rows,max_pairs;
    int32_t *edges;double seconds,minimum_gap;SmRimLearn *learning;
} SmRim;
static int smr_compare_key(const void *aa,const void *bb)
{
    uint64_t a=*(const uint64_t *)aa,b=*(const uint64_t *)bb;return a<b?-1:a>b;
}
static int smr_compare_box(const void *aa,const void *bb)
{
    const SmRimBox *a=aa,*b=bb;
    if(a->lo[0]!=b->lo[0])return a->lo[0]<b->lo[0]?-1:1;
    if(a->a!=b->a)return a->a<b->a?-1:1;return a->b<b->b?-1:a->b>b->b;
}
static int smr_setup(Arena_T arena,size_t nv,const int32_t *faces,size_t nf,SmRim *s)
{
    if(!arena||!nv||nv>(size_t)INT32_MAX||!faces||!nf||!s||nf>SIZE_MAX/3/sizeof(uint64_t))return -1;
    memset(s,0,sizeof *s);s->nv=nv;s->minimum_gap=DBL_MAX;
    uint64_t *keys=ARENA_ALLOC(arena,3*nf*sizeof *keys);
    for(size_t f=0;f<nf;f++)for(size_t k=0;k<3;k++) {
        int32_t a=faces[3*f+k],b=faces[3*f+(k+1)%3];
        if(a<0||b<0||(size_t)a>=nv||(size_t)b>=nv||a==b)return -1;
        keys[3*f+k]=((uint64_t)(uint32_t)(a<b?a:b)<<32)|(uint32_t)(a<b?b:a);
    }
    qsort(keys,3*nf,sizeof *keys,smr_compare_key);
    for(size_t i=0;i<3*nf;) {
        size_t end=i+1;while(end<3*nf&&keys[end]==keys[i])end++;
        if(end-i>2)return -1;if(end-i==1)keys[s->nb++]=keys[i];i=end;
    }
    s->edges=ARENA_ALLOC(arena,2*s->nb*sizeof *s->edges);
    for(size_t e=0;e<s->nb;e++){s->edges[2*e]=(int32_t)(keys[e]>>32);s->edges[2*e+1]=(int32_t)(uint32_t)keys[e];}
    return 0;
}
/* Same relative four-axis separating witness as the unchanged native path
 * certificate, including tangents for collinear separated segments. */
static int smr_axis(const double *uv,const SmRimBox *a,const SmRimBox *b,double normal[2],double *gap)
{
    const double *p=uv+2*(size_t)a->a,*q=uv+2*(size_t)a->b;
    const double *r=uv+2*(size_t)b->a,*t=uv+2*(size_t)b->b;
    double edges[2][2]={{q[0]-p[0],q[1]-p[1]},{t[0]-r[0],t[1]-r[1]}},best=0;
    for(int e=0;e<2;e++)for(int perpendicular=0;perpendicular<2;perpendicular++) {
        double nx=perpendicular?-edges[e][1]:edges[e][0],ny=perpendicular?edges[e][0]:edges[e][1];
        double length=hypot(nx,ny);if(!(length>0))continue;nx/=length;ny/=length;
        double aq=(q[0]-p[0])*nx+(q[1]-p[1])*ny;
        double br=(r[0]-p[0])*nx+(r[1]-p[1])*ny,bt=(t[0]-p[0])*nx+(t[1]-p[1])*ny;
        double forward=fmin(br,bt)-fmax(0,aq),reverse=fmin(0,aq)-fmax(br,bt);
        if(forward>best){best=forward;normal[0]=nx;normal[1]=ny;}
        if(reverse>best){best=reverse;normal[0]=-nx;normal[1]=-ny;}
    }
    *gap=best;return best>0;
}
/* Fold repeated source indices exactly, so fixed cut offsets produce constant
 * feasible rows instead of zero Schur pivots. No geometric pair is exempted. */
static int smr_canonical(SmBound *r)
{
    SmBound result={0};result.rhs=r->rhs;result.face=r->face;result.upper=r->upper;
    size_t n=0;
    for(size_t k=0;k<6;k++)if(r->g[k]!=0) {
        size_t at=0;while(at<n&&result.index[at]!=r->index[k])at++;
        if(at==n)result.index[n++]=r->index[k];result.g[at]+=r->g[k];
    }
    *r=result;double norm=0;for(size_t k=0;k<6;k++)norm+=r->g[k]*r->g[k];
    return norm>0?1:r->rhs<=0?0:-1;
}
/* Each four-endpoint inequality retains the minimum current projected gap.
 * This is affine in displacements, so satisfying it at both ends proves it
 * throughout the straight path. Other pairs still require native line search. */
static int smr_pair(const double *uv,const int32_t *map,const SmRimBox *a,const SmRimBox *b,
                     double activation,SmBound out[4],size_t *constant,double *out_gap)
{
    double normal[2],gap;
    if(!smr_axis(uv,a,b,normal,&gap))return -1;
    if(gap>activation)return 0;
    int32_t aa[2]={a->a,a->b},bb[2]={b->a,b->b};double values[4],minimum=DBL_MAX;
    for(size_t i=0;i<2;i++)for(size_t j=0;j<2;j++) {
        size_t p=(size_t)aa[i],q=(size_t)bb[j];
        values[2*i+j]=(uv[2*q]-uv[2*p])*normal[0]+(uv[2*q+1]-uv[2*p+1])*normal[1];
        minimum=fmin(minimum,values[2*i+j]);
    }
    if(!(minimum>0)||!isfinite(minimum))return -1;
    int count=0;*out_gap=minimum;
    for(size_t i=0;i<2;i++)for(size_t j=0;j<2;j++) {
        SmBound r={0};size_t p=(size_t)aa[i],q=(size_t)bb[j];r.upper=2;r.face=SIZE_MAX;
        r.index[0]=2*map[p];r.index[1]=2*map[p]+1;r.index[2]=2*map[q];r.index[3]=2*map[q]+1;
        r.g[0]=-normal[0];r.g[1]=-normal[1];r.g[2]=normal[0];r.g[3]=normal[1];
        r.rhs=minimum-values[2*i+j];int nonzero=smr_canonical(&r);
        if(nonzero<0)return -1;if(nonzero)out[count++]=r;else (*constant)++;
    }
    return count;
}
static int smr_bounds(Arena_T arena,SmRim *s,const int32_t *map,const double *uv,
                       SmBound **rows,size_t *count)
{
    const double activation=1e-4;double started=ves_clock_sec();
    if(!arena||!s||!map||!uv||!rows||!count)return -1;
    for(size_t v=0;v<s->nv;v++)if(map[v]<0||!isfinite(uv[2*v])||!isfinite(uv[2*v+1]))return -1;
    SmRimBox *boxes=ARENA_ALLOC(arena,s->nb*sizeof *boxes);
    for(size_t k=0;k<s->nb;k++) {
        SmRimBox *b=boxes+k;b->a=s->edges[2*k];b->b=s->edges[2*k+1];
        for(size_t xy=0;xy<2;xy++){b->lo[xy]=fmin(uv[2*(size_t)b->a+xy],uv[2*(size_t)b->b+xy]);b->hi[xy]=fmax(uv[2*(size_t)b->a+xy],uv[2*(size_t)b->b+xy]);}
    }
    qsort(boxes,s->nb,sizeof *boxes,smr_compare_box);
    size_t required=*count;SmBound *result=NULL;
    for(int pass=0;pass<2;pass++) {
        size_t at=*count,np=0,nconstant=0;double minimum=DBL_MAX;
        if(pass) {
            if(required>SIZE_MAX/sizeof *result)return -1;
            result=ARENA_ALLOC(arena,required*sizeof *result);
            if(*count)memcpy(result,*rows,*count*sizeof *result);
        }
        for(size_t i=0;i<s->nb;i++)for(size_t j=i+1;j<s->nb;j++) {
            const SmRimBox *a=boxes+i,*b=boxes+j;
            if(b->lo[0]>a->hi[0]+activation)break;
            if(b->hi[1]<a->lo[1]-activation||b->lo[1]>a->hi[1]+activation)continue;
            if(a->a==b->a||a->a==b->b||a->b==b->a||a->b==b->b)continue;
            SmBound pair[4];double gap=DBL_MAX;size_t constant=0;
            int n=smr_pair(uv,map,a,b,activation,pair,&constant,&gap);if(n<0)return -1;
            if(gap<DBL_MAX){np++;minimum=fmin(minimum,gap);nconstant+=constant;}
            if(at>SIZE_MAX-(size_t)n)return -1;
            if(pass){if(at+(size_t)n>required)return -1;memcpy(result+at,pair,(size_t)n*sizeof *pair);}
            at+=(size_t)n;
        }
        if(!pass)required=at;
        else {
            if(at!=required)return -1;s->calls++;s->pairs+=np;s->rows+=required-*count;s->constant_rows+=nconstant;
            s->max_pairs=s->max_pairs>np?s->max_pairs:np;s->minimum_gap=fmin(s->minimum_gap,minimum);
        }
    }
    *rows=result;*count=required;s->seconds+=ves_clock_sec()-started;return 0;
}
