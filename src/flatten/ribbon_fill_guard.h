#ifndef RIBBON_FILL_GUARD_INCLUDED
#define RIBBON_FILL_GUARD_INCLUDED

/* Read-only broad phase for invented continuation chords.  It indexes the
 * original triangles, not whichever claimant happened to win a raster cell.
 * No proximity or material-ID threshold is an intersection certificate.
 *
 * The balanced BVH has linear storage and logarithmic depth. Queries have no
 * candidate cap: one long triangle cannot silently truncate the evidence or
 * inflate a global nearest-neighbour search radius. This is a chord guard,
 * not a proof that whole generated patches (or two such patches) are disjoint.
 */
typedef struct {
    float lo[3], hi[3];
    int32_t face;
} RibFillRef;

typedef struct {
    float lo[3], hi[3];
    int32_t first, count, left, right;
} RibFillNode;

typedef struct {
    const float *verts;
    const int32_t *faces;
    RibFillRef *refs;
    RibFillNode *nodes;
    size_t nnodes;
} RibFillGuard;

#define RIB_FILL_LEAF 16

static double rib_fill_center(const RibFillRef *ref, int axis)
{
    return (double)ref->lo[axis] + ref->hi[axis];
}

static void rib_fill_partition(RibFillRef *refs, size_t lo, size_t hi,
                                size_t middle, int axis)
{
    /* Three-way quickselect: equal centers terminate immediately, including
     * meshes made of many identical or coplanar triangles. */
    while (hi - lo > 1) {
        double a = rib_fill_center(refs + lo, axis);
        double b = rib_fill_center(refs + lo + (hi-lo)/2, axis);
        double c = rib_fill_center(refs + hi-1, axis);
        double pivot = fmax(fmin(a,b), fmin(fmax(a,b),c));
        size_t lower=lo, scan=lo, upper=hi;
        while (scan < upper) {
            double value = rib_fill_center(refs + scan, axis);
            if (value < pivot) {
                RibFillRef tmp=refs[lower]; refs[lower++]=refs[scan]; refs[scan++]=tmp;
            } else if (value > pivot) {
                RibFillRef tmp=refs[--upper]; refs[upper]=refs[scan]; refs[scan]=tmp;
            } else scan++;
        }
        if (middle < lower) hi=lower;
        else if (middle >= upper) lo=upper;
        else return;
    }
}

static int32_t rib_fill_build_node(RibFillGuard *guard, size_t first, size_t count)
{
    int32_t id=(int32_t)guard->nnodes++;
    RibFillNode *node=guard->nodes+id;
    int axis=0;
    for (int d=0; d<3; d++) {
        node->lo[d]=INFINITY; node->hi[d]=-INFINITY;
        for (size_t i=first; i<first+count; i++) {
            node->lo[d]=fminf(node->lo[d],guard->refs[i].lo[d]);
            node->hi[d]=fmaxf(node->hi[d],guard->refs[i].hi[d]);
        }
        if ((double)node->hi[d]-node->lo[d] >
            (double)node->hi[axis]-node->lo[axis]) axis=d;
    }
    node->first=(int32_t)first; node->count=(int32_t)count;
    node->left=node->right=-1;
    if (count > RIB_FILL_LEAF) {
        size_t half=count/2;
        rib_fill_partition(guard->refs,first,first+count,first+half,axis);
        node->left=rib_fill_build_node(guard,first,half);
        node->right=rib_fill_build_node(guard,first+half,count-half);
        node->count=0;
    }
    return id;
}

static int rib_fill_guard_build(Arena_T arena, const float *verts, size_t nv,
                                const int32_t *faces, size_t nf,
                                RibFillGuard *guard)
{
    memset(guard,0,sizeof(*guard));
    if (!verts || !faces || !nv || !nf || nf > INT32_MAX || nv > INT32_MAX)
        return -1;
    guard->verts=verts; guard->faces=faces;
    /* Halving at every split gives leaf sizes >= LEAF/2 (except the root). */
    size_t capacity=4*(nf/RIB_FILL_LEAF+1);
    guard->refs=(RibFillRef *)ARENA_ALLOC(arena,nf*sizeof(*guard->refs));
    guard->nodes=(RibFillNode *)ARENA_ALLOC(arena,capacity*sizeof(*guard->nodes));
    for (size_t f=0; f<nf; f++) {
        RibFillRef *ref=guard->refs+f;
        ref->face=(int32_t)f;
        for (int v=0; v<3; v++)
            if (faces[3*f+v]<0 || (size_t)faces[3*f+v]>=nv) return -1;
        for (int d=0; d<3; d++) {
            float a=verts[3*(size_t)faces[3*f]+d];
            float b=verts[3*(size_t)faces[3*f+1]+d];
            float c=verts[3*(size_t)faces[3*f+2]+d];
            if (!isfinite(a) || !isfinite(b) || !isfinite(c)) return -1;
            ref->lo[d]=fminf(a,fminf(b,c));
            ref->hi[d]=fmaxf(a,fmaxf(b,c));
        }
    }
    rib_fill_build_node(guard,0,nf);
    assert(guard->nnodes <= capacity);
    return 0;
}

static int rib_fill_segment_box(const double *a, const double *b,
                                const float *lo, const float *hi)
{
    double t0=0, t1=1;
    for (int d=0; d<3; d++) {
        double delta=b[d]-a[d];
        if (delta == 0) {
            if (a[d]<lo[d] || a[d]>hi[d]) return 0;
        } else {
            double x=(lo[d]-a[d])/delta, y=(hi[d]-a[d])/delta;
            t0=fmax(t0,fmin(x,y)); t1=fmin(t1,fmax(x,y));
            if (t0>t1) return 0;
        }
    }
    return 1;
}

static int rib_fill_segment_triangle(const double *a, const double *b,
                                     const float *p, const float *q, const float *r)
{
    double e1[3], e2[3], dir[3], offset[3], h[3], cross[3];
    for (int d=0; d<3; d++) {
        e1[d]=(double)q[d]-p[d]; e2[d]=(double)r[d]-p[d];
        dir[d]=b[d]-a[d]; offset[d]=a[d]-p[d];
    }
    h[0]=dir[1]*e2[2]-dir[2]*e2[1];
    h[1]=dir[2]*e2[0]-dir[0]*e2[2];
    h[2]=dir[0]*e2[1]-dir[1]*e2[0];
    double det=v3dot(e1,h);
    double scale=fabs(e1[0]*h[0])+fabs(e1[1]*h[1])+fabs(e1[2]*h[2]);
    if (fabs(det)<=128*DBL_EPSILON*scale) return 0;
    double u=v3dot(offset,h)/det;
    cross[0]=offset[1]*e1[2]-offset[2]*e1[1];
    cross[1]=offset[2]*e1[0]-offset[0]*e1[2];
    cross[2]=offset[0]*e1[1]-offset[1]*e1[0];
    double v=v3dot(dir,cross)/det;
    double t=v3dot(e2,cross)/det;
    /* Include triangle edges so a triangulation diagonal cannot hide a stab.
     * Do not reject source anchors at the chord endpoints or coplanar chords. */
    return u>=-1e-10 && v>=-1e-10 && u+v<=1+1e-10 &&
           t>1e-7 && t<1-1e-7;
}

static int rib_fill_guard_crosses(const RibFillGuard *guard,
                                  const double *a, const double *b)
{
    int32_t stack[64];
    size_t top=0;
    stack[top++]=0;
    while (top) {
        const RibFillNode *node=guard->nodes+stack[--top];
        if (!rib_fill_segment_box(a,b,node->lo,node->hi)) continue;
        if (!node->count) {
            /* INT32_MAX leaves need fewer than 32 simultaneous siblings. */
            assert(top+2<=sizeof(stack)/sizeof(stack[0]));
            stack[top++]=node->right; stack[top++]=node->left;
        } else for (int32_t i=0; i<node->count; i++) {
            const RibFillRef *ref=guard->refs+node->first+i;
            if (!rib_fill_segment_box(a,b,ref->lo,ref->hi)) continue;
            const int32_t *f=guard->faces+3*(size_t)ref->face;
            if (rib_fill_segment_triangle(a,b,guard->verts+3*(size_t)f[0],
                    guard->verts+3*(size_t)f[1],guard->verts+3*(size_t)f[2])) return 1;
        }
    }
    return 0;
}

#endif
