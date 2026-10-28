/* Conservative broad phase for sequential rim relaxation. Each vertex moves
 * at most one pixel per sweep, so boxes expanded by that amount remain valid
 * until the next rebuild, including edges whose other endpoint already moved.
 * The final component-wide boundary certificate remains mandatory. */
#define SR_RELAX_RIM_STEP 1.0

static void sr_relax_rim_index(const SRSystem *s, size_t component,
                                const double *uv)
{
    size_t begin=s->component_boundary_offset[component];
    size_t end=s->component_boundary_offset[component+1],count=end-begin;
    for (size_t k=0;k<count;k++) {
        const SRBoundaryEdge *e=s->boundary+s->component_boundary[begin+k];
        const double *a=uv+2*(size_t)e->a,*b=uv+2*(size_t)e->b;
        SREdgeBox *box=s->relax_boxes+k;
        box->a=e->a;box->b=e->b;box->component=(int32_t)component;
        box->min_u=fmin(a[0],b[0])-SR_RELAX_RIM_STEP;
        box->max_u=fmax(a[0],b[0])+SR_RELAX_RIM_STEP;
        box->min_v=fmin(a[1],b[1])-SR_RELAX_RIM_STEP;
        box->max_v=fmax(a[1],b[1])+SR_RELAX_RIM_STEP;
    }
    qsort(s->relax_boxes,count,sizeof(SREdgeBox),sr_compare_box);
    double hi=-INFINITY;
    for (size_t k=0;k<count;k++) {
        hi=fmax(hi,s->relax_boxes[k].max_u);
        s->relax_prefix[k]=hi;
    }
}

/* Feasible alpha growth stays on each vertex's current-to-destination segment.
 * Index the whole segment hull, including endpoints that have already moved,
 * so its swept boundary tests need no arbitrary displacement cap. */
static void sr_growth_rim_index(const SRSystem *s,size_t component,
                                 const double *current,const double *destination)
{
    size_t begin=s->component_boundary_offset[component];
    size_t count=s->component_boundary_offset[component+1]-begin;
    for (size_t k=0;k<count;k++) {
        const SRBoundaryEdge *e=s->boundary+s->component_boundary[begin+k];
        const double *a=current+2*(size_t)e->a,*b=current+2*(size_t)e->b;
        const double *c=destination+2*(size_t)e->a,*d=destination+2*(size_t)e->b;
        SREdgeBox *box=s->relax_boxes+k;
        box->a=e->a;box->b=e->b;box->component=(int32_t)component;
        box->min_u=fmin(fmin(a[0],b[0]),fmin(c[0],d[0]));
        box->max_u=fmax(fmax(a[0],b[0]),fmax(c[0],d[0]));
        box->min_v=fmin(fmin(a[1],b[1]),fmin(c[1],d[1]));
        box->max_v=fmax(fmax(a[1],b[1]),fmax(c[1],d[1]));
        double magnitude=fmax(fmax(fabs(box->min_u),fabs(box->max_u)),
                              fmax(fabs(box->min_v),fabs(box->max_v)));
        double roundoff=8*DBL_EPSILON*fmax(1,magnitude);
        box->min_u-=roundoff;box->max_u+=roundoff;
        box->min_v-=roundoff;box->max_v+=roundoff;
    }
    qsort(s->relax_boxes,count,sizeof(SREdgeBox),sr_compare_box);
    double hi=-INFINITY;
    for (size_t k=0;k<count;k++) {
        hi=fmax(hi,s->relax_boxes[k].max_u);
        s->relax_prefix[k]=hi;
    }
}

static int sr_swept_triangle_contains(const double *a,const double *b,
                                      const double *c,const double *p)
{
    int orientation=sr_orientation(a,b,c);
    if (!orientation) return 0;
    return sr_orientation(a,b,p)*orientation>=0 &&
           sr_orientation(b,c,p)*orientation>=0 &&
           sr_orientation(c,a,p)*orientation>=0;
}

/* uv already contains the trial vertex. Testing only its new edges would let
 * a complete boundary segment tunnel through the moving triangle. Test all
 * three swept sides and containment of the other segment's endpoints. */
static int sr_relax_rim_safe(const SRSystem *s, size_t component,
                              const double *uv,size_t v,const double old[2])
{
    size_t count=s->component_boundary_offset[component+1]-
                 s->component_boundary_offset[component];
    const double *now=uv+2*v;
    for (size_t p=s->rim_offset[v];p<s->rim_offset[v+1];p++) {
        const SRBoundaryEdge *edge=s->boundary+s->rim_edge[p];
        size_t w=(size_t)edge->a==v?(size_t)edge->b:(size_t)edge->a;
        const double *fixed=uv+2*w;
        double lo_u=fmin(old[0],fmin(now[0],fixed[0]));
        double hi_u=fmax(old[0],fmax(now[0],fixed[0]));
        double lo_v=fmin(old[1],fmin(now[1],fixed[1]));
        double hi_v=fmax(old[1],fmax(now[1],fixed[1]));
        double eps=1e-10*(hypot(hi_u-lo_u,hi_v-lo_v)+1);
        size_t lo=0,hi=count;
        while (lo<hi) {
            size_t mid=lo+(hi-lo)/2;
            if (s->relax_prefix[mid]<lo_u-eps) lo=mid+1;
            else hi=mid;
        }
        for (size_t k=lo;k<count;k++) {
            const SREdgeBox *box=s->relax_boxes+k;
            if (box->min_u>hi_u+eps) break;
            if (box->max_u<lo_u-eps || box->max_v<lo_v-eps || box->min_v>hi_v+eps) continue;
            if ((size_t)box->a==v || (size_t)box->b==v ||
                (size_t)box->a==w || (size_t)box->b==w) continue;
            const double *a=uv+2*(size_t)box->a,*b=uv+2*(size_t)box->b;
            /* The certificate evaluates stored endpoint order. Its relative
             * collinearity tolerance depends on the first endpoint, so the
             * swept side now->fixed alone can miss a contact accepted as
             * collinear by fixed->now. Match the certificate exactly, while
             * retaining the conservative tests of the whole swept triangle. */
            if (sr_segments_intersect(uv+2*(size_t)edge->a,
                                      uv+2*(size_t)edge->b,a,b) ||
                sr_segments_intersect(old,now,a,b) ||
                sr_segments_intersect(now,fixed,a,b) ||
                sr_segments_intersect(fixed,old,a,b) ||
                sr_swept_triangle_contains(old,now,fixed,a) ||
                sr_swept_triangle_contains(old,now,fixed,b)) return 0;
        }
    }
    return 1;
}

/* Keep enough evidence to distinguish a narrow-phase discrepancy from a stale
 * broad-phase box when the mandatory whole-component certificate rejects a
 * relaxation proposal. This runs only for a rejected, verbose proposal. */
static void sr_relax_rim_diagnose(const SRSystem *s, size_t component,
                                  const double *before, const double *uv,
                                  const SREdgeBox *boxes)
{
    size_t count=s->component_boundary_offset[component+1]-
                 s->component_boundary_offset[component],reported=0;
    for (size_t i=0;i<count;i++) {
        const SREdgeBox *a=boxes+i;
        double eps=1e-10*(hypot(a->max_u-a->min_u,a->max_v-a->min_v)+1);
        for (size_t k=i+1;k<count;k++) {
            const SREdgeBox *b=boxes+k;
            if (b->min_u>a->max_u+eps) break;
            if (b->max_v<a->min_v-eps || b->min_v>a->max_v+eps ||
                a->a==b->a || a->a==b->b || a->b==b->a || a->b==b->b) continue;
            const double *p=uv+2*(size_t)a->a,*q=uv+2*(size_t)a->b;
            const double *r=uv+2*(size_t)b->a,*t=uv+2*(size_t)b->b;
            if (!sr_segments_intersect(p,q,r,t)) continue;
            fprintf(stderr,"      [slim rim] c%zu edges=%d,%d/%d,%d directions=%d/%d/%d/%d\n",
                    component,a->a,a->b,b->a,b->b,
                    sr_segments_intersect(p,q,r,t),sr_segments_intersect(q,p,r,t),
                    sr_segments_intersect(p,q,t,r),sr_segments_intersect(q,p,t,r));
            int32_t vertices[4]={a->a,a->b,b->a,b->b};
            for (size_t j=0;j<4;j++) {
                size_t v=(size_t)vertices[j];
                fprintf(stderr,"      [slim rim] v%zu %.17g %.17g -> %.17g %.17g\n",
                        v,before[2*v],before[2*v+1],uv[2*v],uv[2*v+1]);
            }
            if (++reported>=4) return;
        }
    }
}
