/* EMBEDDED-CHART PRUNING (E18, 2026-10-05).
 *
 * A field's charts never share vertices, so every pair of faces that share a
 * vertex or an edge belongs to one chart; on PHerc1203 cand1 such pairs were
 * 93% of the pair tests of a solve's walks, and they never contribute: their
 * overlap area is zero and their boundary derivative is roundoff (|g| ~ 1e-15).
 * A chart whose map is an embedding has no self-overlap at all, and the
 * embedding has a cheap certificate (the degree of the map at a point off the
 * boundary image is the boundary's winding number there):
 *   - every face of the chart has positive orientation, and
 *   - its boundary loops are simple and pairwise disjoint, and
 *   - exactly one loop is counterclockwise (the others clockwise),
 * so every winding number is 0 or 1 and every point has at most one preimage.
 * The certificate is taken for every chart with an active face before each
 * walk (positivity over its faces, a sort-and-sweep over its boundary with a
 * margin far above roundoff); the walk then skips node pairs inside one
 * certified chart and face pairs of one certified chart. Charts whose boundary
 * is not a set of manifold loops, and any chart that fails a test or comes
 * within the margin of failing, are walked in full as before. */

typedef struct AcSeg { double lo[2], hi[2]; int32_t id; } AcSeg;
typedef struct AcCharts {
    size_t n;                   /* charts present in the field (local ids) */
    int32_t *face;              /* [nf] local chart of each face */
    size_t *first, *faces;      /* CSR: the faces of each chart */
    size_t *seg_first;          /* [n+1] CSR: the boundary segments of each chart */
    int32_t *seg_a, *seg_b;     /* segment endpoints (field vertices), in face orientation */
    int32_t *seg_next, *seg_prev, *seg_loop;   /* loop chaining (segment ids), loop number within the chart */
    uint8_t *certifiable;       /* manifold boundary loops */
    uint8_t *embedded;          /* this walk's certificate */
    uint8_t *active;            /* scratch: chart has an active face */
    int32_t *node;              /* [capacity] the one chart of a node's faces, or -1 */
    AcSeg *sweep; size_t sweep_cap;
    size_t walks, certified, refused;   /* diagnostics */
} AcCharts;

/* Margins: a face is positive when its doubled area exceeds AC_EMBED_AREA of
 * its two edge lengths' product; boundary segments are disjoint when they are
 * more than AC_EMBED_GAP of the chart boundary's coordinate scale apart. */
#define AC_EMBED_AREA 1e-12
#define AC_EMBED_GAP 1e-9

static void acc_nodes(AsmContacts *q)
{
    AcCharts *c=q->charts;if(!c)return;
    for(size_t i=q->nodes;i-->0;){
        AcNode *node=q->tree+i;int32_t chart=-1;
        if(node->left<0){
            chart=c->face[q->order[node->begin]];
            for(size_t j=node->begin+1;j<node->end && chart>=0;j++)if(c->face[q->order[j]]!=chart)chart=-1;
        }else{
            int32_t l=c->node[node->left],r=c->node[node->right];chart=l>=0 && l==r?l:-1;
        }
        c->node[i]=chart;
    }
}

/* Squared distance from p to segment ab. */
static double acc_point_segment(const double p[2],const double a[2],const double b[2])
{
    double dx=b[0]-a[0],dy=b[1]-a[1],len=dx*dx+dy*dy,t=0;
    if(len>0){t=((p[0]-a[0])*dx+(p[1]-a[1])*dy)/len;t=t<0?0:t>1?1:t;}
    double x=a[0]+t*dx-p[0],y=a[1]+t*dy-p[1];return x*x+y*y;
}

/* Whether segments pq and rs cross or come within sqrt(gap2). */
static int acc_segments_near(const double p[2],const double q2[2],const double r[2],const double s[2],double gap2)
{
    double d1=ac_cross(q2[0]-p[0],q2[1]-p[1],r[0]-p[0],r[1]-p[1]),d2=ac_cross(q2[0]-p[0],q2[1]-p[1],s[0]-p[0],s[1]-p[1]);
    double d3=ac_cross(s[0]-r[0],s[1]-r[1],p[0]-r[0],p[1]-r[1]),d4=ac_cross(s[0]-r[0],s[1]-r[1],q2[0]-r[0],q2[1]-r[1]);
    if(((d1>0 && d2<0) || (d1<0 && d2>0)) && ((d3>0 && d4<0) || (d3<0 && d4>0)))return 1;
    return acc_point_segment(p,r,s)<=gap2 || acc_point_segment(q2,r,s)<=gap2 ||
           acc_point_segment(r,p,q2)<=gap2 || acc_point_segment(s,p,q2)<=gap2;
}

static int acc_sweep_order(const void *x,const void *y)
{
    const AcSeg *a=x,*b=y;
    if(a->lo[0]!=b->lo[0])return a->lo[0]<b->lo[0]?-1:1;
    return a->id<b->id?-1:a->id>b->id;
}

/* The embedding certificate of local chart k at q->uv. */
static int acc_certify(AsmContacts *q,size_t k)
{
    AcCharts *c=q->charts;const double *uv=q->uv;
    if(!c->certifiable[k])return 0;
    for(size_t i=c->first[k];i<c->first[k+1];i++){
        const int32_t *f=q->faces+3*c->faces[i];
        const double *a=uv+2*(size_t)f[0],*b=uv+2*(size_t)f[1],*d=uv+2*(size_t)f[2];
        double ex=b[0]-a[0],ey=b[1]-a[1],fx=d[0]-a[0],fy=d[1]-a[1],cross=ac_cross(ex,ey,fx,fy);
        if(!(cross>AC_EMBED_AREA*sqrt((ex*ex+ey*ey)*(fx*fx+fy*fy))))return 0;
    }
    size_t s0=c->seg_first[k],s1=c->seg_first[k+1],ns=s1-s0;
    if(ns<3)return 0;
    double scale=1;int loops=0,positive=0;
    for(size_t s=s0;s<s1;s++){
        const double *a=uv+2*(size_t)c->seg_a[s],*b=uv+2*(size_t)c->seg_b[s];
        AcSeg *g=c->sweep+(s-s0);
        for(int d=0;d<2;d++){g->lo[d]=ac_min(a[d],b[d]);g->hi[d]=ac_max(a[d],b[d]);scale=ac_max(scale,ac_max(fabs(a[d]),fabs(b[d])));}
        g->id=(int32_t)s;
        if(c->seg_loop[s]+1>loops)loops=c->seg_loop[s]+1;
    }
    /* Exactly one counterclockwise loop (shoelace over each loop's segments). */
    for(int l=0;l<loops;l++){
        double area=0;size_t count=0;
        for(size_t s=s0;s<s1;s++)if(c->seg_loop[s]==l){
            const double *a=uv+2*(size_t)c->seg_a[s],*b=uv+2*(size_t)c->seg_b[s];area+=ac_cross(a[0],a[1],b[0],b[1]);count++;
        }
        if(count<3 || !isfinite(area))return 0;
        if(area>0)positive++;
    }
    if(positive!=1)return 0;
    double gap=AC_EMBED_GAP*scale,gap2=gap*gap;
    /* Adjacent segments share one vertex; neither may fold back onto the other. */
    for(size_t s=s0;s<s1;s++){
        size_t t=(size_t)c->seg_next[s];
        const double *p=uv+2*(size_t)c->seg_a[s],*v=uv+2*(size_t)c->seg_b[s],*w=uv+2*(size_t)c->seg_b[t];
        if(acc_point_segment(p,v,w)<=gap2 || acc_point_segment(w,p,v)<=gap2)return 0;
    }
    qsort(c->sweep,ns,sizeof *c->sweep,acc_sweep_order);
    for(size_t i=0;i<ns;i++){
        const AcSeg *x=c->sweep+i;
        for(size_t j=i+1;j<ns && c->sweep[j].lo[0]<=x->hi[0]+gap;j++){
            const AcSeg *y=c->sweep+j;
            if(y->lo[1]>x->hi[1]+gap || x->lo[1]>y->hi[1]+gap)continue;
            int32_t a=x->id,b=y->id;
            if(c->seg_next[a]==b || c->seg_prev[a]==b)continue;   /* adjacent: tested above */
            if(acc_segments_near(uv+2*(size_t)c->seg_a[a],uv+2*(size_t)c->seg_b[a],uv+2*(size_t)c->seg_a[b],uv+2*(size_t)c->seg_b[b],gap2))return 0;
        }
    }
    return 1;
}

/* Before a walk: certify every chart with an active face. */
static void acc_certify_all(AsmContacts *q)
{
    AcCharts *c=q->charts;if(!c)return;
    memset(c->embedded,0,c->n);
    if(q->active_faces){
        memset(c->active,0,c->n);
        for(size_t t=0;t<q->nf;t++)if(q->active_faces[t])c->active[c->face[t]]=1;
    }else memset(c->active,1,c->n);
    c->walks++;
    for(size_t k=0;k<c->n;k++)if(c->active[k]){
        c->embedded[k]=(uint8_t)acc_certify(q,k);
        if(c->embedded[k])c->certified++;else c->refused++;
    }
}

static inline int acc_skip_nodes(const AsmContacts *q,int32_t a,int32_t b)
{
    const AcCharts *c=q->charts;if(!c)return 0;
    int32_t x=c->node[a];return x>=0 && x==c->node[b] && c->embedded[x];
}

static inline int acc_skip_faces(const AsmContacts *q,size_t a,size_t b)
{
    const AcCharts *c=q->charts;if(!c)return 0;
    int32_t x=c->face[a];return x==c->face[b] && c->embedded[x];
}

/* Chart labels for the walks' pruning (face_chart: any ids, one per face).
 * The boundary loops are extracted here once: directed face edges without a
 * reverse twin, chained by their start vertex; a chart whose boundary is not
 * a set of manifold loops (a vertex starting two boundary edges, an edge
 * repeated in one direction) is never certified. Only final arrays in the
 * query's arena escape; construction scratch is borrowed from the caller. */
static AcCharts *acc_build_charts(AsmContacts *q,const int32_t *face_chart,size_t n_ids,Arena_T scratch)
{
    Arena_T arena=q->arena;size_t nf=q->nf,nv=q->nv;
    AcCharts *c=ARENA_CALLOC(arena,1,sizeof *c);
    int32_t *local=ARENA_ALLOC(scratch,n_ids*sizeof *local);
    for(size_t i=0;i<n_ids;i++)local[i]=-1;
    c->face=ARENA_ALLOC(arena,nf*sizeof *c->face);
    for(size_t t=0;t<nf;t++){
        int32_t g=face_chart[t];if(local[g]<0)local[g]=(int32_t)c->n++;
        c->face[t]=local[g];
    }
    size_t n=c->n;
    c->first=ARENA_CALLOC(arena,n+1,sizeof *c->first);c->faces=ARENA_ALLOC(arena,nf*sizeof *c->faces);
    for(size_t t=0;t<nf;t++)c->first[c->face[t]+1]++;
    for(size_t k=0;k<n;k++)c->first[k+1]+=c->first[k];
    {
        size_t *fill=ARENA_ALLOC(scratch,(n?n:1)*sizeof *fill);memcpy(fill,c->first,n*sizeof *fill);
        for(size_t t=0;t<nf;t++)c->faces[fill[c->face[t]]++]=t;
    }
    /* Directed edges by start vertex (CSR), then the boundary: no reverse twin. */
    size_t *out_first=ARENA_CALLOC(scratch,nv+1,sizeof *out_first);int32_t *out_to=ARENA_ALLOC(scratch,3*nf*sizeof *out_to);
    for(size_t t=0;t<nf;t++)for(int k=0;k<3;k++)out_first[(size_t)q->faces[3*t+k]+1]++;
    for(size_t v=0;v<nv;v++)out_first[v+1]+=out_first[v];
    {
        size_t *fill=ARENA_ALLOC(scratch,nv*sizeof *fill);memcpy(fill,out_first,nv*sizeof *fill);
        for(size_t t=0;t<nf;t++)for(int k=0;k<3;k++)out_to[fill[(size_t)q->faces[3*t+k]]++]=q->faces[3*t+(k+1)%3];
    }
    c->certifiable=ARENA_ALLOC(arena,(n?n:1));memset(c->certifiable,1,n);
    c->seg_first=ARENA_CALLOC(arena,n+1,sizeof *c->seg_first);
    int32_t *bstart=ARENA_ALLOC(scratch,nv*sizeof *bstart);
    for(size_t v=0;v<nv;v++)bstart[v]=-1;
    size_t nseg=0;
    for(size_t t=0;t<nf;t++)for(int k=0;k<3;k++){
        int32_t u=q->faces[3*t+k],w=q->faces[3*t+(k+1)%3];int twin=0,same=0;
        if(u==w){c->certifiable[c->face[t]]=0;continue;}
        for(size_t e=out_first[(size_t)w];e<out_first[(size_t)w+1];e++)twin+=out_to[e]==u;
        for(size_t e=out_first[(size_t)u];e<out_first[(size_t)u+1];e++)same+=out_to[e]==w;
        if(same!=1 || twin>1){c->certifiable[c->face[t]]=0;continue;}
        if(!twin){c->seg_first[c->face[t]+1]++;nseg++;}
    }
    for(size_t k=0;k<n;k++)c->seg_first[k+1]+=c->seg_first[k];
    c->seg_a=ARENA_ALLOC(arena,(nseg?nseg:1)*sizeof *c->seg_a);c->seg_b=ARENA_ALLOC(arena,(nseg?nseg:1)*sizeof *c->seg_b);
    c->seg_next=ARENA_ALLOC(arena,(nseg?nseg:1)*sizeof *c->seg_next);c->seg_prev=ARENA_ALLOC(arena,(nseg?nseg:1)*sizeof *c->seg_prev);
    c->seg_loop=ARENA_ALLOC(arena,(nseg?nseg:1)*sizeof *c->seg_loop);
    {
        size_t *fill=ARENA_ALLOC(scratch,(n?n:1)*sizeof *fill);memcpy(fill,c->seg_first,n*sizeof *fill);
        for(size_t t=0;t<nf;t++)for(int k=0;k<3;k++){
            int32_t u=q->faces[3*t+k],w=q->faces[3*t+(k+1)%3];int twin=0,same=0;
            if(u==w)continue;
            for(size_t e=out_first[(size_t)w];e<out_first[(size_t)w+1];e++)twin+=out_to[e]==u;
            for(size_t e=out_first[(size_t)u];e<out_first[(size_t)u+1];e++)same+=out_to[e]==w;
            if(same!=1 || twin)continue;
            size_t s=fill[c->face[t]]++;c->seg_a[s]=u;c->seg_b[s]=w;
            if(bstart[u]==-1)bstart[u]=(int32_t)s;else{bstart[u]=-2;c->certifiable[c->face[t]]=0;}   /* two boundary edges leave u */
        }
    }
    size_t widest=0;
    for(size_t k=0;k<n;k++){
        size_t s0=c->seg_first[k],s1=c->seg_first[k+1];if(s1-s0>widest)widest=s1-s0;
        for(size_t s=s0;s<s1;s++){c->seg_next[s]=c->seg_prev[s]=-1;c->seg_loop[s]=-1;}
        if(!c->certifiable[k] || s1-s0<3){c->certifiable[k]=0;continue;}
        for(size_t s=s0;s<s1 && c->certifiable[k];s++){
            int32_t t=bstart[c->seg_b[s]];
            if(t<0 || (size_t)t<s0 || (size_t)t>=s1){c->certifiable[k]=0;break;}
            c->seg_next[s]=t;
        }
        if(!c->certifiable[k])continue;
        for(size_t s=s0;s<s1;s++){
            int32_t t=c->seg_next[s];
            if(c->seg_prev[t]!=-1){c->certifiable[k]=0;break;}
            c->seg_prev[t]=(int32_t)s;
        }
        if(!c->certifiable[k])continue;
        int32_t loops=0;
        for(size_t s=s0;s<s1;s++)if(c->seg_loop[s]<0){
            size_t at=s,steps=0;
            while(c->seg_loop[at]<0 && steps<=s1-s0){c->seg_loop[at]=loops;at=(size_t)c->seg_next[at];steps++;}
            if(at!=s){c->certifiable[k]=0;break;}
            loops++;
        }
    }
    c->sweep=ARENA_ALLOC(arena,(widest?widest:1)*sizeof *c->sweep);c->sweep_cap=widest;
    c->embedded=ARENA_CALLOC(arena,n?n:1,1);c->active=ARENA_CALLOC(arena,n?n:1,1);
    c->node=ARENA_ALLOC(arena,q->capacity*sizeof *c->node);
    return c;
}

int AsmContacts_set_charts(AsmContacts *q,const int32_t *face_chart,size_t n_ids)
{
    if(!q || !face_chart || !n_ids)return -1;
    for(size_t t=0;t<q->nf;t++)if(face_chart[t]<0 || (size_t)face_chart[t]>=n_ids)return -1;
    /* Only the final chart/segment arrays belong to the query. The local-id
     * map, CSR cursors, directed edges and boundary-start map are construction
     * scratch (20*nv + 12*nf bytes plus chart maps on 64-bit hosts). Dispose,
     * rather than restore, so these chunks cannot remain on a free list. */
    Arena_T scratch=Arena_new_sized((size_t)1<<20);AcCharts *c=NULL;
    TRY
        c=acc_build_charts(q,face_chart,n_ids,scratch);
    FINALLY
        Arena_dispose(&scratch);   /* also release scratch if any arena allocation raises */
    END_TRY;
    /* Publish only after construction and scratch disposal have succeeded. */
    q->charts=c;acc_nodes(q);
    q->version++;   /* a memo enumerated without pruning must not replay */
    return 0;
}

void AsmContacts_chart_counts(const AsmContacts *q,size_t out[3])
{
    out[0]=q && q->charts?q->charts->walks:0;out[1]=q && q->charts?q->charts->certified:0;out[2]=q && q->charts?q->charts->refused:0;
}

/* Control: grids as charts. The pruned walk must report what the full walk
 * reports (pairs and areas) on separate charts, overlapping charts, a chart
 * folded over itself (never certified) and a chart with a flipped face. */
static void acc_grid(double *uv,int32_t *faces,size_t *nv,size_t *nf,int nx,int ny,double x0,double y0,double h,int fold)
{
    size_t base=*nv;
    for(int j=0;j<=ny;j++)for(int i=0;i<=nx;i++){
        double x=x0+h*i,y=y0+h*j;
        if(fold && i>nx/2){x=x0+h*(nx-i)+.37*h;y+=.25*h;}   /* the right half folds back over the left */
        uv[2*(*nv)]=x;uv[2*(*nv)+1]=y;(*nv)++;
    }
    for(int j=0;j<ny;j++)for(int i=0;i<nx;i++){
        int32_t a=(int32_t)(base+(size_t)(j*(nx+1)+i)),b=a+1,c=a+nx+1,d=c+1;
        faces[3*(*nf)]=a;faces[3*(*nf)+1]=b;faces[3*(*nf)+2]=d;(*nf)++;
        faces[3*(*nf)]=a;faces[3*(*nf)+1]=d;faces[3*(*nf)+2]=c;(*nf)++;
    }
}

static int acc_control(void)
{
    int failures=0;
    for(int scene=0;scene<4;scene++){
        Arena_T arena=Arena_new();
        double *uv=ARENA_ALLOC(arena,2*2000*sizeof(double));int32_t *faces=ARENA_ALLOC(arena,3*4000*sizeof(int32_t));
        int32_t *chart=ARENA_ALLOC(arena,4000*sizeof(int32_t));size_t nv=0,nf=0,cut=0;
        acc_grid(uv,faces,&nv,&nf,12,10,0,0,1,scene==2);cut=nf;
        if(scene==0)acc_grid(uv,faces,&nv,&nf,12,10,20,0,1,0);         /* apart */
        else acc_grid(uv,faces,&nv,&nf,12,10,6.3,3.1,1,0);             /* overlapping the first */
        if(scene==3){int32_t t=faces[3*5+1];faces[3*5+1]=faces[3*5+2];faces[3*5+2]=t;}   /* one flipped face in chart 0 */
        for(size_t t=0;t<nf;t++)chart[t]=t<cut?7:3;
        AsmContacts *full=AsmContacts_new(arena,nv,faces,nf,uv,0),*pruned=AsmContacts_new(arena,nv,faces,nf,uv,0);
        if(!full || !pruned || AsmContacts_set_charts(pruned,chart,8)){failures++;Arena_dispose(&arena);continue;}
        AsmContactStats a,b;double ga=0,gb=0;
        int ra=AsmContacts_measure(full,uv,NULL,NULL,NULL,&a),rb=AsmContacts_measure(pruned,uv,NULL,NULL,NULL,&b);
        if(ra || rb || a.pairs!=b.pairs || fabs(a.area-b.area)>1e-9*fmax(1,a.area))failures++;
        double *g1=ARENA_ALLOC(arena,2*nv*sizeof(double)),*g2=ARENA_ALLOC(arena,2*nv*sizeof(double));
        ra=AsmContacts_objective(full,uv,g1,NULL,NULL,&a,&ga);rb=AsmContacts_objective(pruned,uv,g2,NULL,NULL,&b,&gb);
        double worst=0,size=1;for(size_t i=0;i<2*nv;i++){worst=fmax(worst,fabs(g1[i]-g2[i]));size=fmax(size,fabs(g1[i]));}
        if(ra || rb || fabs(ga-gb)>1e-9*fmax(1,fabs(ga)) || worst>1e-9*size)failures++;
        size_t counts[3];AsmContacts_chart_counts(pruned,counts);
        /* scene 0/1: both charts certified; scene 2: the folded chart refused; scene 3: the flipped chart refused */
        size_t expect_refused=scene>=2?1:0;
        if(counts[2]!=2*expect_refused)failures++;
        if(scene<2 && a.pairs!=(scene==1?a.pairs:0))failures++;
        if(scene==1 && !(a.pairs>0))failures++;
        if(scene==2 && !(a.pairs>0))failures++;
        Arena_dispose(&arena);
    }
    fprintf(stderr,"  contacts embedded-chart pruning (apart, overlapping, folded, flipped): %s (%d failures)\n",failures?"FAIL":"ok",failures);
    return failures;
}
