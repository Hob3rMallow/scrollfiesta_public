#include "intersection_cleanup.h"
#include "../common/arena.h"
#include "../common/pipeline_constants.h"
#include "../common/ves_platform.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct {
    int32_t face;
    uint32_t morton;
    double lo[3];
    double hi[3];
} FaceBox;

typedef struct {
    double lo[3];
    double hi[3];
    size_t begin, end;
    int32_t left, right;
} FaceBoxNode;

typedef struct {
    uint32_t index;
    uint32_t morton;
} FaceBoxOrder;

typedef int (*FaceBoxPairVisitor)(const FaceBox *a, const FaceBox *b,
                                  void *context);

typedef struct {
    int32_t fa;
    int32_t fb;
    uint8_t kind;
    uint8_t shared;
} Conflict;

#define FACEBOX_BVH_LEAF 16

static uint32_t facebox_part1by2(uint32_t x)
{
    x &= 0x000003ffu;
    x = (x | (x << 16)) & 0x030000ffu;
    x = (x | (x << 8))  & 0x0300f00fu;
    x = (x | (x << 4))  & 0x030c30c3u;
    x = (x | (x << 2))  & 0x09249249u;
    return x;
}

static uint32_t facebox_quantize(double value, double lo, double hi)
{
    double t;
    if (!(hi > lo)) return 0;
    t = (value - lo) / (hi - lo);
    if (t <= 0.0) return 0;
    if (t >= 1.0) return 1023;
    return (uint32_t)(t * 1023.0 + 0.5);
}

/* Stable two-pass radix ordering by the 32-bit Morton key.  Face boxes are
 * created in ascending face order, so stability supplies the comparator's
 * secondary face-id ordering exactly.  Keeping an index order avoids qsort's
 * O(n log n) random swaps of 56-byte FaceBox records on full-scroll meshes. */
static int facebox_morton_order(const FaceBox *boxes, size_t nbox,
                                FaceBoxOrder **out_order)
{
    const size_t nbucket = 65536;
    FaceBoxOrder *order = NULL, *scratch = NULL, *src, *dst;
    size_t *bucket = NULL;
    if (!boxes || !out_order || nbox == 0 ||
        nbox > UINT32_MAX || nbox > SIZE_MAX / sizeof(*order))
        return -1;
    *out_order = NULL;
    order = (FaceBoxOrder *)malloc(nbox * sizeof(*order));
    scratch = (FaceBoxOrder *)malloc(nbox * sizeof(*scratch));
    bucket = (size_t *)malloc(nbucket * sizeof(*bucket));
    if (!order || !scratch || !bucket) goto fail;
    for (size_t i = 0; i < nbox; i++) {
        order[i].index = (uint32_t)i;
        order[i].morton = boxes[i].morton;
    }
    src = order;
    dst = scratch;
    for (unsigned shift = 0; shift < 32; shift += 16) {
        size_t sum = 0;
        memset(bucket, 0, nbucket * sizeof(*bucket));
        for (size_t i = 0; i < nbox; i++)
            bucket[(src[i].morton >> shift) & UINT32_C(0xffff)]++;
        for (size_t b = 0; b < nbucket; b++) {
            size_t count = bucket[b];
            bucket[b] = sum;
            sum += count;
        }
        if (sum != nbox) goto fail;
        for (size_t i = 0; i < nbox; i++) {
            uint32_t key = (src[i].morton >> shift) & UINT32_C(0xffff);
            dst[bucket[key]++] = src[i];
        }
        {
            FaceBoxOrder *swap = src;
            src = dst;
            dst = swap;
        }
    }
    if (src != order) memcpy(order, src, nbox * sizeof(*order));
    free(scratch);
    free(bucket);
    *out_order = order;
    return 0;
fail:
    free(order);
    free(scratch);
    free(bucket);
    return -1;
}

static int32_t facebox_bvh_build(FaceBoxNode *nodes, size_t node_cap,
                                 size_t *next_node, const FaceBox *boxes,
                                 const FaceBoxOrder *order,
                                 size_t begin, size_t end)
{
    size_t index;
    FaceBoxNode *node;
    if (!nodes || !next_node || !boxes || !order || begin >= end ||
        *next_node >= node_cap || *next_node > (size_t)INT32_MAX)
        return -1;
    index = (*next_node)++;
    node = &nodes[index];
    node->begin = begin;
    node->end = end;
    node->left = node->right = -1;
    if (end - begin <= FACEBOX_BVH_LEAF) {
        size_t i;
        int k;
        for (k = 0; k < 3; k++) {
            node->lo[k] = DBL_MAX;
            node->hi[k] = -DBL_MAX;
        }
        for (i = begin; i < end; i++) {
            const FaceBox *box = &boxes[order[i].index];
            for (k = 0; k < 3; k++) {
                if (box->lo[k] < node->lo[k])
                    node->lo[k] = box->lo[k];
                if (box->hi[k] > node->hi[k])
                    node->hi[k] = box->hi[k];
            }
        }
    } else {
        size_t mid = begin + (end - begin) / 2;
        int32_t left = facebox_bvh_build(
            nodes, node_cap, next_node, boxes, order, begin, mid);
        int32_t right = facebox_bvh_build(
            nodes, node_cap, next_node, boxes, order, mid, end);
        int k;
        if (left < 0 || right < 0) return -1;
        node = &nodes[index];
        node->left = left;
        node->right = right;
        for (k = 0; k < 3; k++) {
            node->lo[k] = nodes[left].lo[k] < nodes[right].lo[k] ?
                          nodes[left].lo[k] : nodes[right].lo[k];
            node->hi[k] = nodes[left].hi[k] > nodes[right].hi[k] ?
                          nodes[left].hi[k] : nodes[right].hi[k];
        }
    }
    return (int32_t)index;
}

static int facebox_bounds_overlap(const double alo[3], const double ahi[3],
                                  const double blo[3], const double bhi[3],
                                  double gap)
{
    int k;
    for (k = 0; k < 3; k++)
        if (ahi[k] + gap < blo[k] || bhi[k] + gap < alo[k])
            return 0;
    return 1;
}

static int facebox_bvh_visit_pair(const FaceBoxNode *nodes,
                                  const FaceBox *boxes,
                                  const FaceBoxOrder *order,
                                  int32_t ia, int32_t ib, double gap,
                                  FaceBoxPairVisitor visitor, void *context)
{
    const FaceBoxNode *a = &nodes[ia], *b = &nodes[ib];
    int a_leaf = a->left < 0, b_leaf = b->left < 0;
    if (ia != ib &&
        !facebox_bounds_overlap(a->lo, a->hi, b->lo, b->hi, gap))
        return 0;
    if (ia == ib) {
        if (a_leaf) {
            size_t i, j;
            for (i = a->begin; i < a->end; i++)
                for (j = i + 1; j < a->end; j++)
                    if (facebox_bounds_overlap(
                            boxes[order[i].index].lo,
                            boxes[order[i].index].hi,
                            boxes[order[j].index].lo,
                            boxes[order[j].index].hi, gap) &&
                        visitor(&boxes[order[i].index],
                                &boxes[order[j].index], context) != 0)
                        return -1;
            return 0;
        }
        if (facebox_bvh_visit_pair(nodes, boxes, order, a->left, a->left,
                                   gap, visitor, context) != 0 ||
            facebox_bvh_visit_pair(nodes, boxes, order, a->left, a->right,
                                   gap, visitor, context) != 0 ||
            facebox_bvh_visit_pair(nodes, boxes, order, a->right, a->right,
                                   gap, visitor, context) != 0)
            return -1;
        return 0;
    }
    if (a_leaf && b_leaf) {
        size_t i, j;
        for (i = a->begin; i < a->end; i++)
            for (j = b->begin; j < b->end; j++)
                if (facebox_bounds_overlap(
                        boxes[order[i].index].lo,
                        boxes[order[i].index].hi,
                        boxes[order[j].index].lo,
                        boxes[order[j].index].hi, gap) &&
                    visitor(&boxes[order[i].index],
                            &boxes[order[j].index], context) != 0)
                    return -1;
        return 0;
    }
    if (a_leaf || (!b_leaf && b->end - b->begin > a->end - a->begin)) {
        if (facebox_bvh_visit_pair(nodes, boxes, order, ia, b->left,
                                   gap, visitor, context) != 0 ||
            facebox_bvh_visit_pair(nodes, boxes, order, ia, b->right,
                                   gap, visitor, context) != 0)
            return -1;
    } else {
        if (facebox_bvh_visit_pair(nodes, boxes, order, a->left, ib,
                                   gap, visitor, context) != 0 ||
            facebox_bvh_visit_pair(nodes, boxes, order, a->right, ib,
                                   gap, visitor, context) != 0)
            return -1;
    }
    return 0;
}

/*
 * Morton ordering supplies deterministic 3-D locality; the balanced AABB tree
 * then enumerates every and only potentially overlapping leaf pair.  The exact
 * triangle predicate remains the sole arbiter.  This broad phase is invariant
 * to scroll axis, cube size, meshing density, and region dimensions.
 */
static int facebox_visit_overlaps(FaceBox *boxes, size_t nbox, double gap,
                                  FaceBoxPairVisitor visitor, void *context)
{
    FaceBoxNode *nodes = NULL;
    FaceBoxOrder *order = NULL;
    double center_lo[3] = {DBL_MAX, DBL_MAX, DBL_MAX};
    double center_hi[3] = {-DBL_MAX, -DBL_MAX, -DBL_MAX};
    size_t nleaf, leaf_cap = 1, node_cap, next_node = 0, i;
    int32_t root;
    int rc = -1, k;
    if (!boxes || !visitor || nbox < 2 || gap < 0.0) return nbox < 2 ? 0 : -1;
    for (i = 0; i < nbox; i++) {
        for (k = 0; k < 3; k++) {
            double center = boxes[i].lo[k] +
                            0.5 * (boxes[i].hi[k] - boxes[i].lo[k]);
            if (!isfinite(center)) goto done;
            if (center < center_lo[k]) center_lo[k] = center;
            if (center > center_hi[k]) center_hi[k] = center;
        }
    }
    for (i = 0; i < nbox; i++) {
        uint32_t q[3];
        for (k = 0; k < 3; k++) {
            double center = boxes[i].lo[k] +
                            0.5 * (boxes[i].hi[k] - boxes[i].lo[k]);
            q[k] = facebox_quantize(center, center_lo[k], center_hi[k]);
        }
        boxes[i].morton = facebox_part1by2(q[0]) |
                          (facebox_part1by2(q[1]) << 1) |
                          (facebox_part1by2(q[2]) << 2);
    }
    if (facebox_morton_order(boxes, nbox, &order) != 0) goto done;
    nleaf = (nbox + FACEBOX_BVH_LEAF - 1) / FACEBOX_BVH_LEAF;
    while (leaf_cap < nleaf) {
        if (leaf_cap > SIZE_MAX / 2) goto done;
        leaf_cap *= 2;
    }
    node_cap = leaf_cap + (leaf_cap - 1);
    if (node_cap > (size_t)INT32_MAX ||
        node_cap > SIZE_MAX / sizeof(*nodes))
        goto done;
    nodes = (FaceBoxNode *)malloc(node_cap * sizeof(*nodes));
    if (!nodes) goto done;
    root = facebox_bvh_build(
        nodes, node_cap, &next_node, boxes, order, 0, nbox);
    if (root < 0 || next_node > node_cap) goto done;
    if (facebox_bvh_visit_pair(nodes, boxes, order, root, root, gap,
                               visitor, context) != 0)
        goto done;
    rc = 0;
done:
    free(nodes);
    free(order);
    return rc;
}

static void vsub(double out[3], const double a[3], const double b[3])
{
    out[0] = a[0] - b[0];
    out[1] = a[1] - b[1];
    out[2] = a[2] - b[2];
}

static void vcross(double out[3], const double a[3], const double b[3])
{
    out[0] = a[1] * b[2] - a[2] * b[1];
    out[1] = a[2] * b[0] - a[0] * b[2];
    out[2] = a[0] * b[1] - a[1] * b[0];
}

static double vdot(const double a[3], const double b[3])
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

typedef struct {
    FaceBox *boxes;
    FaceBoxNode *nodes;
    FaceBoxOrder *order;
} ClearanceTree;

static int clearance_tree(Arena_T arena, const float *verts, size_t nv,
                            const int32_t *faces, size_t nf, ClearanceTree *tree)
{
    double lo[3]={DBL_MAX,DBL_MAX,DBL_MAX},hi[3]={-DBL_MAX,-DBL_MAX,-DBL_MAX};
    size_t leaves=1,nodes=0,next=0;
    if (!verts || !faces || !nf || nf>INT32_MAX || nf>SIZE_MAX/sizeof(FaceBox)) return -1;
    tree->boxes=(FaceBox *)ARENA_ALLOC(arena,nf*sizeof(FaceBox));
    for (size_t f=0;f<nf;f++) {
        FaceBox *box=tree->boxes+f;
        box->face=(int32_t)f; box->morton=0;
        for (int j=0;j<3;j++) {
            int32_t id=faces[3*f+j];
            if (id<0 || (size_t)id>=nv) return -1;
            for (int k=0;k<3;k++) {
                double value=verts[3*(size_t)id+k];
                if (!isfinite(value)) return -1;
                if (!j) box->lo[k]=box->hi[k]=value;
                else { box->lo[k]=fmin(box->lo[k],value); box->hi[k]=fmax(box->hi[k],value); }
            }
        }
        for (int k=0;k<3;k++) {
            double mid=.5*(box->lo[k]+box->hi[k]);
            lo[k]=fmin(lo[k],mid); hi[k]=fmax(hi[k],mid);
        }
    }
    for (size_t f=0;f<nf;f++) {
        uint32_t q[3]={0};
        for (int k=0;k<3;k++)
            q[k]=facebox_quantize(.5*(tree->boxes[f].lo[k]+tree->boxes[f].hi[k]),lo[k],hi[k]);
        tree->boxes[f].morton=facebox_part1by2(q[0]) |
            (facebox_part1by2(q[1])<<1) | (facebox_part1by2(q[2])<<2);
    }
    if (facebox_morton_order(tree->boxes,nf,&tree->order)!=0) return -1;
    while (leaves<(nf+FACEBOX_BVH_LEAF-1)/FACEBOX_BVH_LEAF) leaves*=2;
    nodes=2*leaves-1;
    if (nodes>INT32_MAX || nodes>SIZE_MAX/sizeof(FaceBoxNode)) return -1;
    tree->nodes=(FaceBoxNode *)ARENA_ALLOC(arena,nodes*sizeof(FaceBoxNode));
    return facebox_bvh_build(tree->nodes,nodes,&next,tree->boxes,tree->order,0,nf)==0 ? 0 : -1;
}

/* Any separating projection bounds Euclidean distance from below. We do not
 * require a complete SAT or treat failure to separate as proof of intersection.
 * Subtract a dot-product roundoff enclosure; close/degenerate cases get zero. */
static double clearance_axis(const double a[3][3],const double b[3][3],
                               const double axis[3])
{
    double length=sqrt(vdot(axis,axis)),scale=0;
    double alo=DBL_MAX,ahi=-DBL_MAX,blo=DBL_MAX,bhi=-DBL_MAX;
    if (!(length>DBL_MIN) || !isfinite(length)) return 0;
    for (int i=0;i<3;i++) {
        double pa=vdot(a[i],axis),pb=vdot(b[i],axis),sa=0,sb=0;
        for (int k=0;k<3;k++) { sa+=fabs(a[i][k]*axis[k]); sb+=fabs(b[i][k]*axis[k]); }
        scale=fmax(scale,fmax(sa,sb));
        alo=fmin(alo,pa); ahi=fmax(ahi,pa); blo=fmin(blo,pb); bhi=fmax(bhi,pb);
    }
    return fmax(0,fmax(blo-ahi,alo-bhi)-64*DBL_EPSILON*scale)/(length*(1+16*DBL_EPSILON));
}

static double clearance_lower_bound(const float *av,const int32_t *af,
                                      const float *bv,const int32_t *bf,double stop)
{
    double a[3][3]={{0}},b[3][3]={{0}},ea[3][3]={{0}},eb[3][3]={{0}};
    double na[3]={0},nb[3]={0},lower=0;
    for (int i=0;i<3;i++) for (int k=0;k<3;k++) {
        a[i][k]=(double)av[3*(size_t)af[i]+k]-av[3*(size_t)af[0]+k];
        b[i][k]=(double)bv[3*(size_t)bf[i]+k]-av[3*(size_t)af[0]+k];
    }
    for (int i=0;i<3;i++) { vsub(ea[i],a[(i+1)%3],a[i]); vsub(eb[i],b[(i+1)%3],b[i]); }
    vcross(na,ea[0],ea[1]); vcross(nb,eb[0],eb[1]);
    lower=fmax(clearance_axis(a,b,na),clearance_axis(a,b,nb));
    if (lower>=stop) return lower;
    for (int i=0;i<3;i++) {
        double axis[3]={0};
        axis[i]=1; lower=fmax(lower,clearance_axis(a,b,axis));
        vcross(axis,na,ea[i]); lower=fmax(lower,clearance_axis(a,b,axis));
        vcross(axis,nb,eb[i]); lower=fmax(lower,clearance_axis(a,b,axis));
        if (lower>=stop) return lower;
        for (int j=0;j<3;j++) {
            vcross(axis,ea[i],eb[j]); lower=fmax(lower,clearance_axis(a,b,axis));
            if (lower>=stop) return lower;
        }
    }
    return lower;
}

typedef struct {
    const float *verts[2];
    const int32_t *faces[2];
    double *budget[2];
    double gap;
    size_t tested;
    const int32_t *owner; /* optional same-mesh whole-domain batch ownership */
    int32_t *candidate_a,*candidate_b;
    double *candidate_budget;
    size_t pending,chunk;
} ClearancePair;

static int clearance_flush(ClearancePair *p)
{
    ptrdiff_t i=0;
    /* Budgets are immutable throughout this parallel region. A snapshot's
     * short-circuit threshold is at least the eventual serial threshold, so
     * extra projections cannot change either folded minimum. No atomics,
     * per-thread face arrays or reordered floating-point reductions. */
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(p->pending>=INTERSECTION_CLEARANCE_PARALLEL_MIN)
#endif
    for (i=0;i<(ptrdiff_t)p->pending;i++) {
        int32_t a=p->candidate_a[i],b=p->candidate_b[i];
        p->candidate_budget[i]=.25*clearance_lower_bound(
            p->verts[0],p->faces[0]+3*(size_t)a,p->verts[1],p->faces[1]+3*(size_t)b,
            4*fmax(p->budget[0][a],p->budget[1][b]));
    }
    for (size_t j=0;j<p->pending;j++) {
        int32_t a=p->candidate_a[j],b=p->candidate_b[j];
        /* Preserve the exact legacy traversal/count semantics: a preceding
         * pending pair may have already made both envelopes zero. */
        if (p->budget[0][a]==0 && p->budget[1][b]==0) continue;
        if (p->tested==SIZE_MAX) return -1;
        p->tested++;
        p->budget[0][a]=fmin(p->budget[0][a],p->candidate_budget[j]);
        p->budget[1][b]=fmin(p->budget[1][b],p->candidate_budget[j]);
    }
    p->pending=0; return 0;
}

static int clearance_visit(const ClearanceTree *ta,const ClearanceTree *tb,
                             int32_t ia,int32_t ib,ClearancePair *p)
{
    const FaceBoxNode *a=ta->nodes+ia,*b=tb->nodes+ib;
    int al=a->left<0,bl=b->left<0;
    if (!facebox_bounds_overlap(a->lo,a->hi,b->lo,b->hi,p->gap)) return 0;
    if (ta==tb && ia==ib && !al) {
        if (clearance_visit(ta,tb,a->left,a->left,p)!=0 ||
            clearance_visit(ta,tb,a->left,a->right,p)!=0 ||
            clearance_visit(ta,tb,a->right,a->right,p)!=0) return -1;
        return 0;
    }
    if (al && bl) {
        for (size_t i=a->begin;i<a->end;i++) for (size_t j=ta==tb && ia==ib ? i+1 : b->begin;j<b->end;j++) {
            const FaceBox *fa=ta->boxes+ta->order[i].index,*fb=tb->boxes+tb->order[j].index;
            double budget=0;
            if (p->owner && p->owner[fa->face]==p->owner[fb->face]) continue;
            if (!facebox_bounds_overlap(fa->lo,fa->hi,fb->lo,fb->hi,p->gap)) continue;
            if (p->budget[0][fa->face]==0 && p->budget[1][fb->face]==0) continue;
            if (p->chunk) {
                p->candidate_a[p->pending]=fa->face; p->candidate_b[p->pending]=fb->face;
                p->pending++;
                if (p->pending==p->chunk && clearance_flush(p)!=0) return -1;
                continue;
            }
            if (p->tested==SIZE_MAX) return -1;
            p->tested++;
            /* Two quarter-distance envelopes leave half the measured lower
             * bound as clearance; this is not a material-evidence weight. */
            /* The old exhaustive result is the maximum of these same
             * projections. Once one is >= four times BOTH current budgets,
             * every remaining projection leaves both fmin results bitwise
             * unchanged. This skips arithmetic, never a candidate pair. */
            budget=.25*clearance_lower_bound(p->verts[0],p->faces[0]+3*(size_t)fa->face,
                p->verts[1],p->faces[1]+3*(size_t)fb->face,
                4*fmax(p->budget[0][fa->face],p->budget[1][fb->face]));
            p->budget[0][fa->face]=fmin(p->budget[0][fa->face],budget);
            p->budget[1][fb->face]=fmin(p->budget[1][fb->face],budget);
        }
        return 0;
    }
    if (al || (!bl && b->end-b->begin>a->end-a->begin)) {
        if (clearance_visit(ta,tb,ia,b->left,p)!=0 || clearance_visit(ta,tb,ia,b->right,p)!=0) return -1;
    } else if (clearance_visit(ta,tb,a->left,ib,p)!=0 || clearance_visit(ta,tb,a->right,ib,p)!=0) return -1;
    return 0;
}

static int clearance_execute(Arena_T arena,const ClearanceTree *a,const ClearanceTree *b,
                                ClearancePair *p,size_t chunk)
{
    p->chunk=chunk;
    if (chunk>PTRDIFF_MAX || chunk>SIZE_MAX/sizeof(double)) return -1;
    if (chunk) {
        p->candidate_a=(int32_t *)ARENA_ALLOC(arena,chunk*sizeof(int32_t));
        p->candidate_b=(int32_t *)ARENA_ALLOC(arena,chunk*sizeof(int32_t));
        p->candidate_budget=(double *)ARENA_ALLOC(arena,chunk*sizeof(double));
    }
    if (clearance_visit(a,b,0,0,p)!=0) return -1;
    return p->pending ? clearance_flush(p) : 0;
}

int IntersectionCleanup_cross_budget(
    const float *av,size_t anv,const int32_t *af,size_t anf,
    const float *bv,size_t bnv,const int32_t *bf,size_t bnf,
    double maximum,double *ab,double *bb,size_t *tested)
{
    Arena_T arena=Arena_new();
    ClearanceTree a={0},b={0};
    ClearancePair p={{av,bv},{af,bf},{ab,bb},0,0,NULL};
    double started=ves_clock_sec(),built=0;
    int rc=-1;
    if (!tested || !isfinite(maximum) || maximum<0 || maximum>DBL_MAX/2 ||
        (anf && !ab) || (bnf && !bb)) goto done;
    *tested=0;
    p.gap=2*maximum;
    for (size_t i=0;i<anf;i++) if (!isfinite(ab[i]) || ab[i]<0 || ab[i]>maximum) goto done;
    for (size_t i=0;i<bnf;i++) if (!isfinite(bb[i]) || bb[i]<0 || bb[i]>maximum) goto done;
    if (!anf || !bnf) { rc=0; goto done; }
    if (clearance_tree(arena,av,anv,af,anf,&a)!=0 || clearance_tree(arena,bv,bnv,bf,bnf,&b)!=0) goto done;
    built=ves_clock_sec();
    rc=clearance_execute(arena,&a,&b,&p,INTERSECTION_CLEARANCE_CHUNK);
    fprintf(stderr,"[clearance] cross trees %.2fs, pair tests/fold %.2fs, bounded chunk %zu pairs\n",
        built-started,ves_clock_sec()-built,p.chunk);
    *tested=p.tested;
done:
    free(a.order); free(b.order); Arena_dispose(&arena); return rc;
}

int IntersectionCleanup_partition_budget(const float *verts,size_t nv,
    const int32_t *faces,size_t nf,const int32_t *owner,double maximum,
    double *budget,size_t *tested)
{
    Arena_T arena=Arena_new();
    ClearanceTree tree={0};
    ClearancePair p={{verts,verts},{faces,faces},{budget,budget},0,0,owner};
    double started=ves_clock_sec(),built=0;
    int rc=-1;
    if (!tested || !isfinite(maximum) || maximum<0 || maximum>DBL_MAX/2 ||
        (nf && (!owner || !budget))) goto done;
    *tested=0; p.gap=2*maximum;
    for (size_t f=0;f<nf;f++)
        if (owner[f]<0 || !isfinite(budget[f]) || budget[f]<0 || budget[f]>maximum) goto done;
    if (nf<2) { rc=0; goto done; }
    if (clearance_tree(arena,verts,nv,faces,nf,&tree)!=0) goto done;
    built=ves_clock_sec();
    rc=clearance_execute(arena,&tree,&tree,&p,INTERSECTION_CLEARANCE_CHUNK); *tested=p.tested;
    fprintf(stderr,"[clearance] partition tree %.2fs, pair tests/fold %.2fs, bounded chunk %zu pairs\n",
        built-started,ves_clock_sec()-built,p.chunk);
done:
    free(tree.order); Arena_dispose(&arena); return rc;
}

static int tri_unit_normal(const double p0[3], const double p1[3],
                           const double p2[3], double n[3])
{
    double e1[3], e2[3], len;
    vsub(e1, p1, p0);
    vsub(e2, p2, p0);
    vcross(n, e1, e2);
    len = sqrt(vdot(n, n));
    if (len < 1e-12) return 0;
    n[0] /= len;
    n[1] /= len;
    n[2] /= len;
    return 1;
}

static double cross2(double ax, double ay, double bx, double by)
{
    return ax * by - ay * bx;
}

static int pt_in_tri2d(double px, double py,
                       double ax, double ay, double bx, double by,
                       double cx, double cy)
{
    double d1 = cross2(bx - ax, by - ay, px - ax, py - ay);
    double d2 = cross2(cx - bx, cy - by, px - bx, py - by);
    double d3 = cross2(ax - cx, ay - cy, px - cx, py - cy);
    const double eps = 1e-9;
    int neg, pos;
    if (fabs(d1) < eps || fabs(d2) < eps || fabs(d3) < eps) return 0;
    neg = (d1 < 0.0) || (d2 < 0.0) || (d3 < 0.0);
    pos = (d1 > 0.0) || (d2 > 0.0) || (d3 > 0.0);
    return !(neg && pos);
}

static int seg_cross2d(double p1x, double p1y, double p2x, double p2y,
                       double q1x, double q1y, double q2x, double q2y)
{
    double d1x = p2x - p1x, d1y = p2y - p1y;
    double d2x = q2x - q1x, d2y = q2y - q1y;
    double denom = cross2(d1x, d1y, d2x, d2y);
    double sx, sy, t, s;
    if (fabs(denom) < 1e-12) return 0;
    sx = q1x - p1x;
    sy = q1y - p1y;
    t = cross2(sx, sy, d2x, d2y) / denom;
    s = cross2(sx, sy, d1x, d1y) / denom;
    return t > 1e-9 && t < 1.0 - 1e-9 &&
           s > 1e-9 && s < 1.0 - 1e-9;
}

static int inplane_overlap(const double n[3],
                           const double v0[3], const double v1[3],
                           const double v2[3], const double u0[3],
                           const double u1[3], const double u2[3])
{
    double ax = fabs(n[0]), ay = fabs(n[1]), az = fabs(n[2]);
    const double *ta[3] = {v0, v1, v2};
    const double *tb[3] = {u0, u1, u2};
    double a[3][2], b[3][2];
    int i0, i1, e, f, k;
    if (ax > ay && ax > az) {
        i0 = 1; i1 = 2;
    } else if (ay > az) {
        i0 = 0; i1 = 2;
    } else {
        i0 = 0; i1 = 1;
    }
    for (k = 0; k < 3; k++) {
        a[k][0] = ta[k][i0]; a[k][1] = ta[k][i1];
        b[k][0] = tb[k][i0]; b[k][1] = tb[k][i1];
    }
    for (e = 0; e < 3; e++)
        for (f = 0; f < 3; f++)
            if (seg_cross2d(a[e][0], a[e][1],
                            a[(e + 1) % 3][0], a[(e + 1) % 3][1],
                            b[f][0], b[f][1],
                            b[(f + 1) % 3][0], b[(f + 1) % 3][1]))
                return 1;
    for (k = 0; k < 3; k++) {
        if (pt_in_tri2d(a[k][0], a[k][1],
                        b[0][0], b[0][1], b[1][0], b[1][1],
                        b[2][0], b[2][1]))
            return 1;
        if (pt_in_tri2d(b[k][0], b[k][1],
                        a[0][0], a[0][1], a[1][0], a[1][1],
                        a[2][0], a[2][1]))
            return 1;
    }
    /*
     * Identical/coincident triangulations can have every vertex on the other
     * triangle's boundary and every edge collinear, so the strict tests above
     * intentionally see only boundary contact.  Their centroids, however, are
     * strictly interior whenever the common region has positive area.  This
     * closes that degeneracy without turning a bare edge/vertex touch into an
     * overlap.
     */
    {
        double acx = (a[0][0] + a[1][0] + a[2][0]) / 3.0;
        double acy = (a[0][1] + a[1][1] + a[2][1]) / 3.0;
        double bcx = (b[0][0] + b[1][0] + b[2][0]) / 3.0;
        double bcy = (b[0][1] + b[1][1] + b[2][1]) / 3.0;
        if (pt_in_tri2d(acx, acy,
                        b[0][0], b[0][1], b[1][0], b[1][1],
                        b[2][0], b[2][1]) ||
            pt_in_tri2d(bcx, bcy,
                        a[0][0], a[0][1], a[1][0], a[1][1],
                        a[2][0], a[2][1]))
            return 1;
    }
    return 0;
}

static int point_in_triangle_strict_3d(const double p[3],
                                       const double a[3],
                                       const double b[3],
                                       const double c[3])
{
    double e0[3], e1[3], q[3];
    double d00, d01, d11, d20, d21, denom, u, v, w;
    const double eps = 1e-10;
    vsub(e0, b, a);
    vsub(e1, c, a);
    vsub(q, p, a);
    d00 = vdot(e0, e0);
    d01 = vdot(e0, e1);
    d11 = vdot(e1, e1);
    d20 = vdot(q, e0);
    d21 = vdot(q, e1);
    denom = d00 * d11 - d01 * d01;
    if (!(denom > 1e-24 * d00 * d11)) return 0;
    u = (d11 * d20 - d01 * d21) / denom;
    v = (d00 * d21 - d01 * d20) / denom;
    w = 1.0 - u - v;
    return u > eps && v > eps && w > eps;
}

static int segment_triangle_proper(const double p0[3], const double p1[3],
                                   const double a[3], const double b[3],
                                   const double c[3])
{
    double n[3], q0[3], q1[3], d0, d1, t, p[3];
    const double plane_eps = 1e-7;
    if (!tri_unit_normal(a, b, c, n)) return 0;
    vsub(q0, p0, a);
    vsub(q1, p1, a);
    d0 = vdot(n, q0);
    d1 = vdot(n, q1);
    if (!((d0 > plane_eps && d1 < -plane_eps) ||
          (d0 < -plane_eps && d1 > plane_eps)))
        return 0;
    t = d0 / (d0 - d1);
    p[0] = p0[0] + t * (p1[0] - p0[0]);
    p[1] = p0[1] + t * (p1[1] - p0[1]);
    p[2] = p0[2] + t * (p1[2] - p0[2]);
    return point_in_triangle_strict_3d(p, a, b, c);
}

static int moller_interpenetrate(const double v0[3], const double v1[3],
                                 const double v2[3], const double u0[3],
                                 const double u1[3], const double u2[3])
{
    /* A proper non-coplanar triangle intersection has an edge that crosses
     * the strict interior of the other triangle.  Testing all six directed
     * edges makes the verdict explicitly invariant to face winding, cyclic
     * rotation, and pair order; plane/edge/vertex-only contact stays benign. */
    return segment_triangle_proper(v0, v1, u0, u1, u2) ||
           segment_triangle_proper(v1, v2, u0, u1, u2) ||
           segment_triangle_proper(v2, v0, u0, u1, u2) ||
           segment_triangle_proper(u0, u1, v0, v1, v2) ||
           segment_triangle_proper(u1, u2, v0, v1, v2) ||
           segment_triangle_proper(u2, u0, v0, v1, v2);
}

static int tri_pair_test(const double v0[3], const double v1[3],
                         const double v2[3], const double u0[3],
                         const double u1[3], const double u2[3],
                         double cos_parallel, double gap_max)
{
    double n1[3], n2[3], plane_d;
    double s0, s1, s2, gmax_u_to_v, gmax_v_to_u;
    if (!tri_unit_normal(v0, v1, v2, n1)) return 0;
    if (!tri_unit_normal(u0, u1, u2, n2)) return 0;
    if (fabs(vdot(n1, n2)) >= cos_parallel) {
        /*
         * Use the smaller of the two directed slab distances.  A one-sided
         * test is order-dependent: a small triangle can lie wholly within the
         * gap of a coarse triangle's plane even though remote vertices of the
         * coarse triangle lie outside the small triangle's plane slab.  That
         * coarse/fine T-junction is precisely a local doubled-surface defect.
         */
        plane_d = -vdot(n1, v0);
        s0 = fabs(vdot(n1, u0) + plane_d);
        s1 = fabs(vdot(n1, u1) + plane_d);
        s2 = fabs(vdot(n1, u2) + plane_d);
        gmax_u_to_v = s0 > s1 ? (s0 > s2 ? s0 : s2)
                              : (s1 > s2 ? s1 : s2);
        plane_d = -vdot(n2, u0);
        s0 = fabs(vdot(n2, v0) + plane_d);
        s1 = fabs(vdot(n2, v1) + plane_d);
        s2 = fabs(vdot(n2, v2) + plane_d);
        gmax_v_to_u = s0 > s1 ? (s0 > s2 ? s0 : s2)
                              : (s1 > s2 ? s1 : s2);
        if ((gmax_u_to_v <= gap_max || gmax_v_to_u <= gap_max) &&
            (inplane_overlap(n1, v0, v1, v2, u0, u1, u2) ||
             inplane_overlap(n2, v0, v1, v2, u0, u1, u2)))
            return INTERSECTION_HIT_OVERLAP;
        /* A shallow-angle but genuine crossing is still an exact stab.  The
         * angular split controls proximity classification, not topology. */
        return moller_interpenetrate(v0, v1, v2, u0, u1, u2)
                   ? INTERSECTION_HIT_STAB : 0;
    }
    return moller_interpenetrate(v0, v1, v2, u0, u1, u2)
               ? INTERSECTION_HIT_STAB : 0;
}

static int shared_vertices(const float *verts,
                           const int32_t *a, const int32_t *b)
{
    int n = 0, i, j;
    const double eps2 = 1e-10; /* cloned pinch vertices are coordinate-identical */
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            double d0, d1, d2;
            if (a[i] == b[j]) {
                n++;
                break;
            }
            d0 = (double)verts[(size_t)a[i] * 3 + 0] -
                 (double)verts[(size_t)b[j] * 3 + 0];
            d1 = (double)verts[(size_t)a[i] * 3 + 1] -
                 (double)verts[(size_t)b[j] * 3 + 1];
            d2 = (double)verts[(size_t)a[i] * 3 + 2] -
                 (double)verts[(size_t)b[j] * 3 + 2];
            if (d0*d0 + d1*d1 + d2*d2 <= eps2) {
                n++;
                break;
            }
        }
    }
    return n;
}

static int32_t dsu_find(int32_t *parent, int32_t x)
{
    int32_t r = x;
    while (parent[r] != r) r = parent[r];
    while (parent[x] != x) {
        int32_t next = parent[x];
        parent[x] = r;
        x = next;
    }
    return r;
}

static void dsu_union(int32_t *parent, uint8_t *rank, int32_t a, int32_t b)
{
    int32_t ra = dsu_find(parent, a), rb = dsu_find(parent, b);
    if (ra == rb) return;
    if (rank[ra] < rank[rb]) {
        parent[ra] = rb;
    } else if (rank[ra] > rank[rb]) {
        parent[rb] = ra;
    } else {
        parent[rb] = ra;
        rank[ra]++;
    }
}

static double face_area(const float *verts, const int32_t *face)
{
    double a[3], b[3], c[3], ab[3], ac[3], cr[3];
    int k;
    for (k = 0; k < 3; k++) {
        a[k] = verts[(size_t)face[0] * 3 + (size_t)k];
        b[k] = verts[(size_t)face[1] * 3 + (size_t)k];
        c[k] = verts[(size_t)face[2] * 3 + (size_t)k];
    }
    vsub(ab, b, a);
    vsub(ac, c, a);
    vcross(cr, ab, ac);
    return 0.5 * sqrt(vdot(cr, cr));
}

static int append_conflict(Conflict **pairs, size_t *count, size_t *capacity,
                           size_t cap, int32_t fa, int32_t fb,
                           int kind, int shared)
{
    Conflict *grown;
    size_t next;
    if (*count >= cap) return -1;
    if (*count == *capacity) {
        next = *capacity ? *capacity * 2 : 1024;
        if (next > cap) next = cap;
        grown = (Conflict *)realloc(*pairs, next * sizeof(**pairs));
        if (!grown) return -1;
        *pairs = grown;
        *capacity = next;
    }
    (*pairs)[*count].fa = fa;
    (*pairs)[*count].fb = fb;
    (*pairs)[*count].kind = (uint8_t)kind;
    (*pairs)[*count].shared = (uint8_t)shared;
    (*count)++;
    return 0;
}

void IntersectionCleanup_default_params(IntersectionCleanupParams *params)
{
    if (!params) return;
    params->gap_max = 1.0;
    params->parallel_angle_deg = 20.0;
    params->max_delete_fraction = 0.01;
    params->max_delete_fraction_hard = 0.02;
    params->min_delete_budget_faces = 128;
    params->max_conflicts = 1000000;
    params->include_hinges = 1;
    params->hit_kind_mask = INTERSECTION_HIT_MASK_ALL;
}

int IntersectionCleanup_pair_test(const float *verts, size_t nv,
                                  const int32_t face_a[3],
                                  const int32_t face_b[3],
                                  const IntersectionCleanupParams *params,
                                  int *shared_out)
{
    double av[3][3], bv[3][3], cos_parallel;
    int shared, hit, k;
    if (shared_out) *shared_out = 0;
    if (!verts || !face_a || !face_b || !params || nv == 0 ||
        params->gap_max < 0.0 || params->parallel_angle_deg < 0.0 ||
        params->parallel_angle_deg >= 90.0 ||
        (params->hit_kind_mask & ~INTERSECTION_HIT_MASK_ALL) != 0)
        return -1;
    for (k = 0; k < 3; k++) {
        if (face_a[k] < 0 || face_b[k] < 0 ||
            (size_t)face_a[k] >= nv || (size_t)face_b[k] >= nv)
            return -1;
        av[0][k] = verts[(size_t)face_a[0] * 3 + (size_t)k];
        av[1][k] = verts[(size_t)face_a[1] * 3 + (size_t)k];
        av[2][k] = verts[(size_t)face_a[2] * 3 + (size_t)k];
        bv[0][k] = verts[(size_t)face_b[0] * 3 + (size_t)k];
        bv[1][k] = verts[(size_t)face_b[1] * 3 + (size_t)k];
        bv[2][k] = verts[(size_t)face_b[2] * 3 + (size_t)k];
    }
    shared = shared_vertices(verts, face_a, face_b);
    if (shared_out) *shared_out = shared;
    if (shared == 1 && !params->include_hinges) return 0;
    if (shared >= 3)
        return params->hit_kind_mask & INTERSECTION_HIT_MASK_OVERLAP
             ? INTERSECTION_HIT_OVERLAP : 0;
    cos_parallel = cos(params->parallel_angle_deg * M_PI / 180.0);
    hit = tri_pair_test(av[0], av[1], av[2], bv[0], bv[1], bv[2],
                        cos_parallel, params->gap_max);
    if (shared > 0 && hit == INTERSECTION_HIT_STAB) return 0;
    hit = shared > 0 && hit ? INTERSECTION_HIT_FOLD : hit;
    if (hit && !(params->hit_kind_mask & (1u << (hit - 1)))) return 0;
    return hit;
}

typedef struct {
    const float *verts;
    size_t nv;
    const int32_t *faces;
    const IntersectionCleanupParams *params;
    size_t *face_conflict_degree;
    IntersectionCleanupConflictVisitor visitor;
    void *visitor_context;
    IntersectionCleanupStats *stats;
} IntersectionAuditPairContext;

static int intersection_audit_facebox_pair(const FaceBox *box_a,
                                           const FaceBox *box_b,
                                           void *context)
{
    IntersectionAuditPairContext *c =
        (IntersectionAuditPairContext *)context;
    int32_t fa = box_a->face, fb = box_b->face;
    int shared = 0, hit;
    if (fa > fb) {
        int32_t t = fa; fa = fb; fb = t;
    }
    c->stats->candidate_pairs++;
    hit = IntersectionCleanup_pair_test(
        c->verts, c->nv, &c->faces[(size_t)fa * 3],
        &c->faces[(size_t)fb * 3], c->params, &shared);
    if (hit < 0) return -1;
    if (!hit) return 0;
    if (!(c->params->hit_kind_mask & (1u << (hit - 1)))) return 0;
    if (c->stats->conflicts >= c->params->max_conflicts) return -1;
    c->stats->conflicts++;
    if (hit == INTERSECTION_HIT_OVERLAP)
        c->stats->overlap_pairs++;
    else if (hit == INTERSECTION_HIT_STAB)
        c->stats->stab_pairs++;
    else
        c->stats->fold_pairs++;
    if (c->visitor &&
        c->visitor((size_t)fa, (size_t)fb, hit,
                   c->visitor_context) != 0)
        return -1;
    if (c->face_conflict_degree) {
        size_t da = ++c->face_conflict_degree[(size_t)fa];
        size_t db = ++c->face_conflict_degree[(size_t)fb];
        if (da > c->stats->max_conflict_degree)
            c->stats->max_conflict_degree = da;
        if (db > c->stats->max_conflict_degree)
            c->stats->max_conflict_degree = db;
    }
    return 0;
}

int IntersectionCleanup_audit_visit(
                              const float *verts, size_t nv,
                              const int32_t *faces, size_t nf,
                              const uint8_t *face_mask,
                              const IntersectionCleanupParams *params,
                              size_t *face_conflict_degree,
                              IntersectionCleanupConflictVisitor visitor,
                              void *visitor_context,
                              IntersectionCleanupStats *stats)
{
    IntersectionCleanupStats local_stats;
    IntersectionAuditPairContext pair_context;
    FaceBox *boxes = NULL;
    size_t nbox = 0, i, k;
    int rc = -1;

    memset(&local_stats, 0, sizeof local_stats);
    if (stats) memset(stats, 0, sizeof *stats);
    if (!verts || !faces || !params || nv == 0 || nf == 0 ||
        nf > (size_t)INT32_MAX || nf > SIZE_MAX / sizeof(*boxes) ||
        (face_conflict_degree && nf > SIZE_MAX /
                                      sizeof(*face_conflict_degree)) ||
        params->gap_max < 0.0 || params->parallel_angle_deg < 0.0 ||
        params->parallel_angle_deg >= 90.0 || params->max_conflicts == 0 ||
        (params->hit_kind_mask & ~INTERSECTION_HIT_MASK_ALL) != 0)
        return -1;
    if (face_conflict_degree)
        memset(face_conflict_degree, 0, nf * sizeof *face_conflict_degree);
    boxes = (FaceBox *)malloc(nf * sizeof *boxes);
    if (!boxes) return -1;
    for (i = 0; i < nf; i++) {
        const int32_t *f = &faces[i * 3];
        FaceBox box;
        if (face_mask && !face_mask[i]) continue;
        if (f[0] < 0 || f[1] < 0 || f[2] < 0 ||
            (size_t)f[0] >= nv || (size_t)f[1] >= nv ||
            (size_t)f[2] >= nv)
            goto cleanup;
        box.face = (int32_t)i;
        box.morton = 0;
        for (k = 0; k < 3; k++) {
            double a = verts[(size_t)f[0] * 3 + k];
            double b = verts[(size_t)f[1] * 3 + k];
            double c = verts[(size_t)f[2] * 3 + k];
            if (!isfinite(a) || !isfinite(b) || !isfinite(c)) goto cleanup;
            box.lo[k] = a < b ? (a < c ? a : c) : (b < c ? b : c);
            box.hi[k] = a > b ? (a > c ? a : c) : (b > c ? b : c);
        }
        boxes[nbox++] = box;
    }
    if (nbox < 2) { rc = 0; goto cleanup; }
    memset(&pair_context, 0, sizeof pair_context);
    pair_context.verts = verts;
    pair_context.nv = nv;
    pair_context.faces = faces;
    pair_context.params = params;
    pair_context.face_conflict_degree = face_conflict_degree;
    pair_context.visitor = visitor;
    pair_context.visitor_context = visitor_context;
    pair_context.stats = &local_stats;
    if (facebox_visit_overlaps(
            boxes, nbox, params->gap_max,
            intersection_audit_facebox_pair, &pair_context) != 0)
        goto cleanup;
    rc = 0;
cleanup:
    if (stats) *stats = local_stats;
    free(boxes);
    return rc;
}

int IntersectionCleanup_audit(const float *verts, size_t nv,
                              const int32_t *faces, size_t nf,
                              const uint8_t *face_mask,
                              const IntersectionCleanupParams *params,
                              size_t *face_conflict_degree,
                              IntersectionCleanupStats *stats)
{
    return IntersectionCleanup_audit_visit(
        verts, nv, faces, nf, face_mask, params, face_conflict_degree,
        NULL, NULL, stats);
}

/* ---- parallel audit ------------------------------------------------------
 * The serial traversal collects candidate pairs into a bounded chunk; the
 * exact triangle tests -- the dominant cost -- run over each chunk with
 * OpenMP, and stats/visitor/degree fold-in happens serially in collection
 * order.  Results are therefore identical to the serial audit and invariant
 * to thread count.  The one behavioural difference is fail-closed timing: a
 * max_conflicts overflow or visitor abort is detected at the next chunk
 * flush rather than at the exact offending pair (rc is -1 either way, and
 * candidate_pairs may be larger on that failure path). */

#define INTERSECTION_PARALLEL_CHUNK ((size_t)1 << 22)

typedef struct {
    const float *verts;
    size_t nv;
    const int32_t *faces;
    const IntersectionCleanupParams *params;
    size_t *face_conflict_degree;
    IntersectionCleanupConflictVisitor visitor;
    void *visitor_context;
    IntersectionCleanupStats *stats;
    int32_t *cand_a;
    int32_t *cand_b;
    int8_t *cand_hit;
    size_t pending;
    size_t chunk;
} IntersectionParallelContext;

static int intersection_parallel_flush(IntersectionParallelContext *c)
{
    ptrdiff_t i;
    size_t j;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (i = 0; i < (ptrdiff_t)c->pending; i++) {
        int shared = 0;
        int hit = IntersectionCleanup_pair_test(
            c->verts, c->nv, &c->faces[(size_t)c->cand_a[i] * 3],
            &c->faces[(size_t)c->cand_b[i] * 3], c->params, &shared);
        c->cand_hit[i] = (int8_t)hit;
    }
    for (j = 0; j < c->pending; j++) {
        int hit = c->cand_hit[j];
        int32_t fa = c->cand_a[j], fb = c->cand_b[j];
        if (hit < 0) return -1;
        if (!hit) continue;
        if (c->stats->conflicts >= c->params->max_conflicts) return -1;
        c->stats->conflicts++;
        if (hit == INTERSECTION_HIT_OVERLAP)
            c->stats->overlap_pairs++;
        else if (hit == INTERSECTION_HIT_STAB)
            c->stats->stab_pairs++;
        else
            c->stats->fold_pairs++;
        if (c->visitor &&
            c->visitor((size_t)fa, (size_t)fb, hit,
                       c->visitor_context) != 0)
            return -1;
        if (c->face_conflict_degree) {
            size_t da = ++c->face_conflict_degree[(size_t)fa];
            size_t db = ++c->face_conflict_degree[(size_t)fb];
            if (da > c->stats->max_conflict_degree)
                c->stats->max_conflict_degree = da;
            if (db > c->stats->max_conflict_degree)
                c->stats->max_conflict_degree = db;
        }
    }
    c->pending = 0;
    return 0;
}

static int intersection_parallel_collect(const FaceBox *box_a,
                                         const FaceBox *box_b,
                                         void *context)
{
    IntersectionParallelContext *c = (IntersectionParallelContext *)context;
    int32_t fa = box_a->face, fb = box_b->face;
    if (fa > fb) { int32_t t = fa; fa = fb; fb = t; }
    c->stats->candidate_pairs++;
    c->cand_a[c->pending] = fa;
    c->cand_b[c->pending] = fb;
    c->pending++;
    if (c->pending == c->chunk) return intersection_parallel_flush(c);
    return 0;
}

int IntersectionCleanup_audit_visit_parallel(
                              const float *verts, size_t nv,
                              const int32_t *faces, size_t nf,
                              const uint8_t *face_mask,
                              const IntersectionCleanupParams *params,
                              size_t *face_conflict_degree,
                              IntersectionCleanupConflictVisitor visitor,
                              void *visitor_context,
                              IntersectionCleanupStats *stats)
{
    IntersectionCleanupStats local_stats;
    IntersectionParallelContext c;
    FaceBox *boxes = NULL;
    size_t nbox = 0, i, k;
    int rc = -1;

    memset(&local_stats, 0, sizeof local_stats);
    memset(&c, 0, sizeof c);
    if (stats) memset(stats, 0, sizeof *stats);
    if (!verts || !faces || !params || nv == 0 || nf == 0 ||
        nf > (size_t)INT32_MAX || nf > SIZE_MAX / sizeof(*boxes) ||
        (face_conflict_degree && nf > SIZE_MAX /
                                      sizeof(*face_conflict_degree)) ||
        params->gap_max < 0.0 || params->parallel_angle_deg < 0.0 ||
        params->parallel_angle_deg >= 90.0 || params->max_conflicts == 0 ||
        (params->hit_kind_mask & ~INTERSECTION_HIT_MASK_ALL) != 0)
        return -1;
    if (face_conflict_degree)
        memset(face_conflict_degree, 0, nf * sizeof *face_conflict_degree);
    boxes = (FaceBox *)malloc(nf * sizeof *boxes);
    if (!boxes) return -1;
    for (i = 0; i < nf; i++) {
        const int32_t *f = &faces[i * 3];
        FaceBox box;
        if (face_mask && !face_mask[i]) continue;
        if (f[0] < 0 || f[1] < 0 || f[2] < 0 ||
            (size_t)f[0] >= nv || (size_t)f[1] >= nv ||
            (size_t)f[2] >= nv)
            goto cleanup;
        box.face = (int32_t)i;
        box.morton = 0;
        for (k = 0; k < 3; k++) {
            double a = verts[(size_t)f[0] * 3 + k];
            double b = verts[(size_t)f[1] * 3 + k];
            double cc = verts[(size_t)f[2] * 3 + k];
            if (!isfinite(a) || !isfinite(b) || !isfinite(cc)) goto cleanup;
            box.lo[k] = a < b ? (a < cc ? a : cc) : (b < cc ? b : cc);
            box.hi[k] = a > b ? (a > cc ? a : cc) : (b > cc ? b : cc);
        }
        boxes[nbox++] = box;
    }
    if (nbox < 2) { rc = 0; goto cleanup; }
    c.verts = verts;
    c.nv = nv;
    c.faces = faces;
    c.params = params;
    c.face_conflict_degree = face_conflict_degree;
    c.visitor = visitor;
    c.visitor_context = visitor_context;
    c.stats = &local_stats;
    c.chunk = INTERSECTION_PARALLEL_CHUNK;
    c.cand_a = (int32_t *)malloc(c.chunk * sizeof *c.cand_a);
    c.cand_b = (int32_t *)malloc(c.chunk * sizeof *c.cand_b);
    c.cand_hit = (int8_t *)malloc(c.chunk * sizeof *c.cand_hit);
    if (!c.cand_a || !c.cand_b || !c.cand_hit) goto cleanup;
    if (facebox_visit_overlaps(
            boxes, nbox, params->gap_max,
            intersection_parallel_collect, &c) != 0)
        goto cleanup;
    if (c.pending && intersection_parallel_flush(&c) != 0) goto cleanup;
    rc = 0;
cleanup:
    if (stats) *stats = local_stats;
    free(c.cand_a);
    free(c.cand_b);
    free(c.cand_hit);
    free(boxes);
    return rc;
}

typedef struct {
    const float *verts;
    size_t nv;
    const int32_t *faces;
    const IntersectionCleanupParams *params;
    double cos_parallel;
    Conflict **pairs;
    size_t *pair_count;
    size_t *pair_capacity;
    IntersectionCleanupStats *stats;
} IntersectionProcessPairContext;

static int intersection_process_facebox_pair(const FaceBox *box_a,
                                             const FaceBox *box_b,
                                             void *context)
{
    IntersectionProcessPairContext *c =
        (IntersectionProcessPairContext *)context;
    int32_t fa = box_a->face, fb = box_b->face;
    const int32_t *af, *bf;
    double av[3][3], bv[3][3];
    int shared, hit, kind, k;
    if (fa > fb) {
        int32_t t = fa; fa = fb; fb = t;
    }
    c->stats->candidate_pairs++;
    af = &c->faces[(size_t)fa * 3];
    bf = &c->faces[(size_t)fb * 3];
    shared = shared_vertices(c->verts, af, bf);
    if (shared == 1 && !c->params->include_hinges) return 0;
    if (shared >= 3) {
        hit = INTERSECTION_HIT_OVERLAP;
    } else {
        for (k = 0; k < 3; k++) {
            av[0][k] = c->verts[(size_t)af[0] * 3 + (size_t)k];
            av[1][k] = c->verts[(size_t)af[1] * 3 + (size_t)k];
            av[2][k] = c->verts[(size_t)af[2] * 3 + (size_t)k];
            bv[0][k] = c->verts[(size_t)bf[0] * 3 + (size_t)k];
            bv[1][k] = c->verts[(size_t)bf[1] * 3 + (size_t)k];
            bv[2][k] = c->verts[(size_t)bf[2] * 3 + (size_t)k];
        }
        hit = tri_pair_test(av[0], av[1], av[2],
                            bv[0], bv[1], bv[2],
                            c->cos_parallel, c->params->gap_max);
    }
    if (!hit) return 0;
    if (shared > 0 && hit == INTERSECTION_HIT_STAB)
        return 0;
    kind = shared > 0 ? INTERSECTION_HIT_FOLD : hit;
    if (!(c->params->hit_kind_mask & (1u << (kind - 1)))) return 0;
    if (append_conflict(c->pairs, c->pair_count, c->pair_capacity,
                        c->params->max_conflicts, fa, fb,
                        kind, shared) != 0)
        return -1;
    if (kind == INTERSECTION_HIT_OVERLAP)
        c->stats->overlap_pairs++;
    else if (kind == INTERSECTION_HIT_STAB)
        c->stats->stab_pairs++;
    else
        c->stats->fold_pairs++;
    return 0;
}

typedef struct {
    uint32_t edge;
    int32_t other;
} ConflictAdj;

typedef struct {
    int32_t *item;
    int32_t *position;
    size_t n;
    const int32_t *active_face;
    const size_t *degree;
    const size_t *component_size;
    const size_t *shared_degree;
    const double *area;
} ConflictHeap;

static int conflict_heap_better(const ConflictHeap *heap,
                                int32_t active_a, int32_t active_b)
{
    int32_t a = heap->active_face[active_a];
    int32_t b = heap->active_face[active_b];
    int shared_a, shared_b;
    if (heap->degree[a] != heap->degree[b])
        return heap->degree[a] > heap->degree[b];
    if (heap->component_size[a] != heap->component_size[b])
        return heap->component_size[a] < heap->component_size[b];
    shared_a = heap->shared_degree[a] != 0;
    shared_b = heap->shared_degree[b] != 0;
    /* A strict total order makes heap maintenance possible.  Fold faces retain
     * the established small-flap preference; pure stab/overlap faces retain
     * the large-area preference.  Mixed ties prioritize resolving the fold. */
    if (shared_a != shared_b) return shared_a > shared_b;
    if (heap->area[a] != heap->area[b])
        return shared_a ? heap->area[a] < heap->area[b]
                        : heap->area[a] > heap->area[b];
    return a > b;
}

static void conflict_heap_swap(ConflictHeap *heap, size_t a, size_t b)
{
    int32_t ia = heap->item[a], ib = heap->item[b];
    heap->item[a] = ib;
    heap->item[b] = ia;
    heap->position[ia] = (int32_t)b;
    heap->position[ib] = (int32_t)a;
}

static void conflict_heap_sift_down(ConflictHeap *heap, size_t at)
{
    for (;;) {
        size_t left = at * 2 + 1, right = left + 1, best = at;
        if (left < heap->n &&
            conflict_heap_better(heap, heap->item[left], heap->item[best]))
            best = left;
        if (right < heap->n &&
            conflict_heap_better(heap, heap->item[right], heap->item[best]))
            best = right;
        if (best == at) break;
        conflict_heap_swap(heap, at, best);
        at = best;
    }
}

static int conflict_heap_init(ConflictHeap *heap, size_t n,
                              const int32_t *active_face,
                              const size_t *degree,
                              const size_t *component_size,
                              const size_t *shared_degree,
                              const double *area)
{
    memset(heap, 0, sizeof(*heap));
    if (n == 0 || n > (size_t)INT32_MAX ||
        n > SIZE_MAX / sizeof(*heap->item))
        return -1;
    heap->item = (int32_t *)malloc(n * sizeof(*heap->item));
    heap->position = (int32_t *)malloc(n * sizeof(*heap->position));
    if (!heap->item || !heap->position) {
        free(heap->item);
        free(heap->position);
        memset(heap, 0, sizeof(*heap));
        return -1;
    }
    heap->n = n;
    heap->active_face = active_face;
    heap->degree = degree;
    heap->component_size = component_size;
    heap->shared_degree = shared_degree;
    heap->area = area;
    for (size_t i = 0; i < n; i++) {
        heap->item[i] = (int32_t)i;
        heap->position[i] = (int32_t)i;
    }
    for (size_t i = n / 2; i > 0; i--)
        conflict_heap_sift_down(heap, i - 1);
    return 0;
}

static void conflict_heap_dispose(ConflictHeap *heap)
{
    free(heap->item);
    free(heap->position);
    memset(heap, 0, sizeof(*heap));
}

static int32_t conflict_heap_pop(ConflictHeap *heap)
{
    int32_t result;
    if (!heap || heap->n == 0) return -1;
    result = heap->item[0];
    heap->position[result] = -1;
    heap->n--;
    if (heap->n > 0) {
        heap->item[0] = heap->item[heap->n];
        heap->position[heap->item[0]] = 0;
        conflict_heap_sift_down(heap, 0);
    }
    return result;
}

static void conflict_heap_decreased(ConflictHeap *heap, int32_t active)
{
    int32_t position;
    if (!heap || active < 0) return;
    position = heap->position[active];
    if (position >= 0) conflict_heap_sift_down(heap, (size_t)position);
}

static int32_t conflict_active_find(const int32_t *face, size_t n,
                                    int32_t query)
{
    size_t lo = 0, hi = n;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (face[mid] < query) lo = mid + 1;
        else hi = mid;
    }
    return lo < n && face[lo] == query ? (int32_t)lo : -1;
}

int IntersectionCleanup_process(const float *verts, size_t nv,
                                int32_t *faces, size_t *pnf,
                                const uint8_t *face_mask,
                                const IntersectionCleanupParams *params,
                                IntersectionCleanupStats *stats)
{
    IntersectionCleanupStats local_stats;
    FaceBox *boxes = NULL;
    Conflict *pairs = NULL;
    size_t nbox = 0, pair_count = 0, pair_capacity = 0;
    size_t nf, i, k;
    double cos_parallel;
    int rc = -1;
    int32_t *parent = NULL, *first_face = NULL;
    int32_t *active_face = NULL;
    ConflictAdj *adj = NULL;
    size_t *adj_off = NULL, *adj_cursor = NULL;
    ConflictHeap heap = {0};
    uint8_t *rank = NULL, *victim = NULL, *resolved = NULL;
    size_t *component_size = NULL, *degree = NULL, *shared_degree = NULL;
    size_t n_active_face = 0;
    double *area = NULL;

    memset(&local_stats, 0, sizeof(local_stats));
    if (stats) memset(stats, 0, sizeof(*stats));
    if (!verts || !faces || !pnf || !params || nv == 0 || *pnf == 0 ||
        params->gap_max < 0.0 ||
        params->parallel_angle_deg < 0.0 ||
        params->parallel_angle_deg >= 90.0 ||
        params->max_delete_fraction < 0.0 ||
        params->max_delete_fraction > 1.0 ||
        params->max_delete_fraction_hard < 0.0 ||
        params->max_delete_fraction_hard > 1.0 ||
        params->max_conflicts == 0 ||
        (params->hit_kind_mask & ~INTERSECTION_HIT_MASK_ALL) != 0)
        return -1;
    nf = *pnf;
    if (nf > (size_t)INT32_MAX || nf > SIZE_MAX / sizeof(*boxes))
        return -1;
    boxes = (FaceBox *)malloc(nf * sizeof(*boxes));
    if (!boxes) goto cleanup;

    for (i = 0; i < nf; i++) {
        const int32_t *f = &faces[i * 3];
        FaceBox box;
        if (face_mask && !face_mask[i]) continue;
        if (f[0] < 0 || f[1] < 0 || f[2] < 0 ||
            (size_t)f[0] >= nv || (size_t)f[1] >= nv || (size_t)f[2] >= nv)
            goto cleanup;
        box.face = (int32_t)i;
        box.morton = 0;
        for (k = 0; k < 3; k++) {
            double a = verts[(size_t)f[0] * 3 + k];
            double b = verts[(size_t)f[1] * 3 + k];
            double c = verts[(size_t)f[2] * 3 + k];
            if (!isfinite(a) || !isfinite(b) || !isfinite(c)) goto cleanup;
            box.lo[k] = a < b ? (a < c ? a : c) : (b < c ? b : c);
            box.hi[k] = a > b ? (a > c ? a : c) : (b > c ? b : c);
        }
        boxes[nbox++] = box;
    }
    if (nbox < 2) {
        rc = 0;
        goto cleanup;
    }
    cos_parallel = cos(params->parallel_angle_deg * M_PI / 180.0);
    {
        IntersectionProcessPairContext pair_context;
        memset(&pair_context, 0, sizeof pair_context);
        pair_context.verts = verts;
        pair_context.nv = nv;
        pair_context.faces = faces;
        pair_context.params = params;
        pair_context.cos_parallel = cos_parallel;
        pair_context.pairs = &pairs;
        pair_context.pair_count = &pair_count;
        pair_context.pair_capacity = &pair_capacity;
        pair_context.stats = &local_stats;
        if (facebox_visit_overlaps(
                boxes, nbox, params->gap_max,
                intersection_process_facebox_pair, &pair_context) != 0)
            goto cleanup;
    }
    local_stats.conflicts = pair_count;
    if (pair_count == 0) {
        rc = 0;
        goto cleanup;
    }

    parent = (int32_t *)malloc(nf * sizeof(*parent));
    rank = (uint8_t *)calloc(nf, sizeof(*rank));
    first_face = (int32_t *)malloc(nv * sizeof(*first_face));
    component_size = (size_t *)calloc(nf, sizeof(*component_size));
    degree = (size_t *)calloc(nf, sizeof(*degree));
    shared_degree = (size_t *)calloc(nf, sizeof(*shared_degree));
    area = (double *)malloc(nf * sizeof(*area));
    victim = (uint8_t *)calloc(nf, sizeof(*victim));
    resolved = (uint8_t *)calloc(pair_count, sizeof(*resolved));
    if (pair_count > SIZE_MAX / 2 ||
        pair_count * 2 > SIZE_MAX / sizeof(*active_face))
        goto cleanup;
    active_face = (int32_t *)malloc(
        (pair_count ? pair_count * 2 : 1) * sizeof(*active_face));
    if (!parent || !rank || !first_face || !component_size || !degree ||
        !shared_degree || !area || !victim || !resolved || !active_face)
        goto cleanup;

    for (i = 0; i < nf; i++) parent[i] = (int32_t)i;
    for (i = 0; i < nv; i++) first_face[i] = -1;
    for (i = 0; i < nf; i++) {
        const int32_t *f = &faces[i * 3];
        area[i] = face_area(verts, f);
        for (k = 0; k < 3; k++) {
            int32_t v = f[k];
            if (first_face[v] < 0)
                first_face[v] = (int32_t)i;
            else
                dsu_union(parent, rank, (int32_t)i, first_face[v]);
        }
    }
    for (i = 0; i < nf; i++)
        component_size[dsu_find(parent, (int32_t)i)]++;
    for (i = 0; i < nf; i++)
        component_size[i] = component_size[dsu_find(parent, (int32_t)i)];
    for (i = 0; i < pair_count; i++) {
        int32_t a = pairs[i].fa, b = pairs[i].fb;
        degree[a]++;
        degree[b]++;
        if (pairs[i].shared) {
            shared_degree[a]++;
            shared_degree[b]++;
        }
    }
    for (i = 0; i < nf; i++) {
        if (degree[i] == 0) continue;
        if (n_active_face >= pair_count * 2) goto cleanup;
        active_face[n_active_face++] = (int32_t)i;
        if (degree[i] > local_stats.max_conflict_degree)
            local_stats.max_conflict_degree = degree[i];
    }

    /* Conflict-local CSR plus a mutable max heap turns the cover from
     * O(victims * (all_faces + all_conflicts)) into
     * O((active_faces + conflicts) log active_faces). */
    if (pair_count > UINT32_MAX || n_active_face == SIZE_MAX ||
        n_active_face + 1 > SIZE_MAX / sizeof(*adj_off) ||
        n_active_face > SIZE_MAX / sizeof(*adj_cursor) ||
        pair_count > SIZE_MAX / 2 ||
        pair_count * 2 > SIZE_MAX / sizeof(*adj))
        goto cleanup;
    adj_off = (size_t *)calloc(n_active_face + 1, sizeof(*adj_off));
    adj_cursor = (size_t *)malloc(
        (n_active_face ? n_active_face : 1) * sizeof(*adj_cursor));
    adj = (ConflictAdj *)malloc(
        (pair_count ? pair_count * 2 : 1) * sizeof(*adj));
    if (!adj_off || !adj_cursor || !adj) goto cleanup;
    for (i = 0; i < pair_count; i++) {
        int32_t a = conflict_active_find(
            active_face,n_active_face,pairs[i].fa);
        int32_t b = conflict_active_find(
            active_face,n_active_face,pairs[i].fb);
        if (a < 0 || b < 0) goto cleanup;
        adj_off[(size_t)a + 1]++;
        adj_off[(size_t)b + 1]++;
    }
    for (i = 0; i < n_active_face; i++) {
        if (adj_off[i+1] > SIZE_MAX - adj_off[i]) goto cleanup;
        adj_off[i+1] += adj_off[i];
    }
    if (adj_off[n_active_face] != pair_count * 2) goto cleanup;
    memcpy(adj_cursor,adj_off,n_active_face*sizeof(*adj_cursor));
    for (i = 0; i < pair_count; i++) {
        int32_t a = conflict_active_find(
            active_face,n_active_face,pairs[i].fa);
        int32_t b = conflict_active_find(
            active_face,n_active_face,pairs[i].fb);
        size_t qa, qb;
        if (a < 0 || b < 0) goto cleanup;
        qa = adj_cursor[a]++;
        qb = adj_cursor[b]++;
        adj[qa].edge = (uint32_t)i; adj[qa].other = b;
        adj[qb].edge = (uint32_t)i; adj[qb].other = a;
    }
    if (conflict_heap_init(
            &heap,n_active_face,active_face,degree,component_size,
            shared_degree,area) != 0)
        goto cleanup;

    {
        size_t unresolved = pair_count, n_victim = 0;
        while (unresolved > 0) {
            int32_t active = conflict_heap_pop(&heap);
            int32_t best;
            if (active < 0) goto cleanup;
            best = active_face[active];
            if (victim[best] || degree[best] == 0) goto cleanup;
            victim[best] = 1;
            n_victim++;
            for (i = adj_off[active]; i < adj_off[(size_t)active+1]; i++) {
                uint32_t edge = adj[i].edge;
                int32_t other_active = adj[i].other;
                int32_t other = active_face[other_active];
                if (resolved[edge]) continue;
                resolved[edge] = 1;
                unresolved--;
                if (!victim[other] && degree[other] > 0) {
                    degree[other]--;
                    conflict_heap_decreased(&heap,other_active);
                }
            }
            degree[best] = 0;
        }
        size_t nominal_budget = (size_t)ceil(
            params->max_delete_fraction * (double)nf);
        double hard_fraction = params->max_delete_fraction_hard;
        if (hard_fraction < params->max_delete_fraction)
            hard_fraction = params->max_delete_fraction;
        size_t hard_budget = (size_t)ceil(hard_fraction * (double)nf);
        size_t delete_budget = nominal_budget;
        if (delete_budget < params->min_delete_budget_faces)
            delete_budget = params->min_delete_budget_faces;
        if (delete_budget > hard_budget) delete_budget = hard_budget;
        local_stats.faces_deleted = n_victim;
        local_stats.delete_budget = delete_budget;
        if (n_victim > delete_budget) {
            local_stats.budget_rejected = 1;
            goto cleanup;
        }
    }

    {
        uint8_t *component_touched = (uint8_t *)calloc(nf, 1);
        size_t out_nf = 0;
        if (!component_touched) goto cleanup;
        for (i = 0; i < nf; i++)
            if (victim[i])
                component_touched[dsu_find(parent, (int32_t)i)] = 1;
        for (i = 0; i < nf; i++)
            if (component_touched[i]) local_stats.components_touched++;
        for (i = 0; i < nf; i++) {
            if (victim[i]) continue;
            faces[out_nf * 3 + 0] = faces[i * 3 + 0];
            faces[out_nf * 3 + 1] = faces[i * 3 + 1];
            faces[out_nf * 3 + 2] = faces[i * 3 + 2];
            out_nf++;
        }
        free(component_touched);
        *pnf = out_nf;
    }
    rc = 0;

cleanup:
    free(boxes);
    free(pairs);
    free(parent);
    free(active_face);
    free(adj);
    free(adj_off);
    free(adj_cursor);
    conflict_heap_dispose(&heap);
    free(rank);
    free(first_face);
    free(component_size);
    free(degree);
    free(shared_degree);
    free(area);
    free(victim);
    free(resolved);
    if (stats) *stats = local_stats;
    return rc;
}

static void selftest_check(int condition, const char *name, int *fails)
{
    if (condition) {
        fprintf(stderr, "  ok: %s\n", name);
    } else {
        fprintf(stderr, "  FAIL: %s\n", name);
        (*fails)++;
    }
}

typedef struct {
    uint8_t *kind;
    size_t n;
} IntersectionSelftestCollector;

static int intersection_selftest_collect(size_t face_a, size_t face_b,
                                         int hit_kind, void *context)
{
    IntersectionSelftestCollector *c =
        (IntersectionSelftestCollector *)context;
    size_t index;
    if (!c || !c->kind || face_a >= c->n || face_b >= c->n ||
        face_a == face_b)
        return -1;
    if (face_a > face_b) {
        size_t t = face_a; face_a = face_b; face_b = t;
    }
    index = face_a * c->n + face_b;
    if (c->kind[index] != 0) return -1;
    c->kind[index] = (uint8_t)hit_kind;
    return 0;
}

int IntersectionCleanup_selftest(void)
{
    int fails = 0;
    IntersectionCleanupParams params;
    IntersectionCleanupStats stats;
    IntersectionCleanup_default_params(&params);
    params.max_delete_fraction = 0.5;
    params.max_delete_fraction_hard = 0.5;

    {
        float av[9]={0,0,0,2,0,0,0,2,0};
        float bv[9]={0,0,.125f,2,0,.125f,0,2,.125f};
        int32_t face[3]={0,1,2};
        double ab=4,bb=4;
        size_t tested=0;
        selftest_check(IntersectionCleanup_cross_budget(av,3,face,1,bv,3,face,1,4,&ab,&bb,&tested)==0 &&
                       tested==1 && ab>0 && ab<=.03125 && fabs(ab-.03125)<1e-12 && ab==bb,
                       "cross-surface approximation budgets leave positive clearance",&fails);
        for (int i=0;i<9;i++) { av[i]+=1000000; bv[i]+=1000000; }
        ab=bb=4;
        selftest_check(IntersectionCleanup_cross_budget(av,3,face,1,bv,3,face,1,4,&ab,&bb,&tested)==0 &&
                       ab>0 && ab<=.03125 && fabs(ab-.03125)<1e-12,
                       "clearance budget remains conservative at large world coordinates",&fails);
        bv[8]=999999.875f; ab=bb=4;
        selftest_check(IntersectionCleanup_cross_budget(av,3,face,1,bv,3,face,1,4,&ab,&bb,&tested)==0 &&
                       ab==0 && bb==0,
                       "original interpenetration requests exact detail, not deletion",&fails);
        ab=bb=4;
        selftest_check(IntersectionCleanup_cross_budget(av,3,face,1,av,3,face,1,4,&ab,&bb,&tested)==0 &&
                       ab==0 && bb==0,
                       "coincident separately indexed surfaces remain explicit contacts",&fails);
    }

    {
        enum { N=64 };
        float verts[9*N]={0};
        int32_t faces[3*N]={0},owner[N]={0};
        double actual[N]={0},expected[N]={0};
        size_t tested=0;
        int ok=1;
        for (int f=0;f<N;f++) {
            float x=(float)((f/2)%8)*4,y=(float)((f/2)/8)*4,z=(float)(f%2)*.125f;
            float tri[9]={x,y,z,x+2,y,z,x,y+2,z};
            memcpy(verts+9*f,tri,sizeof tri);
            for (int k=0;k<3;k++) faces[3*f+k]=3*f+k;
            owner[f]=f%2; actual[f]=expected[f]=4;
        }
        for (int a=0;a<N;a++) for (int b=a+1;b<N;b++) if (owner[a]!=owner[b]) {
            double limit=.25*clearance_lower_bound(verts,faces+3*a,verts,faces+3*b,DBL_MAX);
            expected[a]=fmin(expected[a],limit); expected[b]=fmin(expected[b],limit);
        }
        ok=IntersectionCleanup_partition_budget(verts,3*N,faces,N,owner,4,actual,&tested)==0 && tested>0;
        for (int f=0;f<N;f++) ok=ok && actual[f]>0 && actual[f]<=.03125 && fabs(actual[f]-expected[f])<1e-12;
        selftest_check(ok,"partition tree covers cross-owner pairs with conservative source budgets",&fails);
        for (int f=0;f<N;f++) { owner[f]=0; actual[f]=4; }
        ok=IntersectionCleanup_partition_budget(verts,3*N,faces,N,owner,4,actual,&tested)==0 && tested==0;
        for (int f=0;f<N;f++) ok=ok && actual[f]==4;
        selftest_check(ok,"a single coupled owner needs no partition envelope",&fails);
        owner[0]=-1;
        selftest_check(IntersectionCleanup_partition_budget(verts,3*N,faces,N,owner,4,actual,&tested)!=0,
                       "negative partition ownership is rejected",&fails);
    }

    {
        uint32_t state=UINT32_C(2847519);
        int32_t face[3]={0,1,2};
        size_t shortened=0;
        int ok=1;
        for (int sample=0;sample<2048;sample++) {
            float a[9]={0},b[9]={0};
            double full=0;
            for (int k=0;k<9;k++) {
                state=UINT32_C(1664525)*state+UINT32_C(1013904223);
                a[k]=(float)(state%4096)/1024;
                state=UINT32_C(1664525)*state+UINT32_C(1013904223);
                b[k]=(float)(state%4096)/1024+(k%3==2 ? 6.0f : 0.0f);
                if (sample%3==0) { a[k]+=1000000; b[k]+=1000000; }
            }
            if (sample%7==0) memcpy(b,a,sizeof a);
            full=.25*clearance_lower_bound(a,face,b,face,DBL_MAX);
            for (int limit=0;limit<8;limit++) {
                double before_a=ldexp(4.,-3*limit),before_b=ldexp(3.,-2*limit);
                double bound=.25*clearance_lower_bound(a,face,b,face,4*fmax(before_a,before_b));
                double fast[2]={fmin(before_a,bound),fmin(before_b,bound)};
                double slow[2]={fmin(before_a,full),fmin(before_b,full)};
                shortened+=bound<full;
                ok=ok && memcmp(fast,slow,sizeof fast)==0;
            }
        }
        selftest_check(ok && shortened>1000,
            "clearance projection short-circuit matches exhaustive budget bits (16384 cases)",&fails);
    }

    {
        enum { N=96 };
        Arena_T arena=Arena_new();
        ClearanceTree a={0},b={0};
        float av[9*N]={0},bv[9*N]={0};
        int32_t faces[3*N]={0},owner[N]={0};
        const size_t chunks[4]={1,7,256,INTERSECTION_CLEARANCE_CHUNK};
        double initial[2*N]={0},expected[2*N]={0},actual[2*N]={0};
        int ok=1;
        for (int f=0;f<N;f++) {
            float x=(float)(f%12),y=(float)((f/12)%4),z=(float)(f%9)*.125f;
            float tri[9]={x,y,z,x+2,y,z,x,y+2,z};
            for (int k=0;k<9;k++) {
                av[9*f+k]=tri[k]+1000000;
                bv[9*f+k]=av[9*f+k]+(k%3==2 && f%7 ? .125f : 0);
            }
            for (int k=0;k<3;k++) faces[3*f+k]=3*f+k;
            owner[f]=f%4;
            initial[f]=f%5 ? 4.0/(1+f%8) : 0;
            initial[N+f]=f%3 ? 4.0/(1+f%7) : 0;
        }
        if (clearance_tree(arena,av,3*N,faces,N,&a)!=0 ||
            clearance_tree(arena,bv,3*N,faces,N,&b)!=0) ok=0;
        for (int same=0;same<2 && ok;same++) {
            ClearancePair serial={{av,same ? av : bv},{faces,faces},
                {expected,same ? expected : expected+N},8,0,same ? owner : NULL};
            memcpy(expected,initial,sizeof expected);
            ok=clearance_execute(arena,&a,same ? &a : &b,&serial,0)==0 && serial.tested>0;
            for (size_t c=0;c<4 && ok;c++) {
                ClearancePair batched={{av,same ? av : bv},{faces,faces},
                    {actual,same ? actual : actual+N},8,0,same ? owner : NULL};
                memcpy(actual,initial,sizeof actual);
                ok=clearance_execute(arena,&a,same ? &a : &b,&batched,chunks[c])==0 &&
                    batched.tested==serial.tested && !memcmp(actual,expected,sizeof actual);
            }
        }
        selftest_check(ok,"batched cross/partition clearance matches serial budget bits AND test counts at four chunk sizes",&fails);
        free(a.order); free(b.order); Arena_dispose(&arena);
    }

    /* Directed slab distance is intentionally symmetric.  All vertices of
     * the small tilted triangle are close to the coarse z=0 plane, while the
     * remote coarse vertex is far from the small triangle's plane. */
    {
        double coarse[3][3] = {{0,0,0},{10,0,0},{0,10,0}};
        double fine[3][3] = {{1,1,0.10},{2,1,0.10},{1,2,0.35}};
        double cp = cos(20.0 * M_PI / 180.0);
        int ab = tri_pair_test(coarse[0], coarse[1], coarse[2],
                               fine[0], fine[1], fine[2], cp, 1.0);
        int ba = tri_pair_test(fine[0], fine[1], fine[2],
                               coarse[0], coarse[1], coarse[2], cp, 1.0);
        selftest_check(ab == INTERSECTION_HIT_OVERLAP &&
                       ba == INTERSECTION_HIT_OVERLAP,
                       "coarse/fine slab test is order-independent", &fails);
    }

    /* A boundary-only contact is not a proper interpenetration, and that
     * verdict must not change when a later winding pass rotates or reverses
     * either face.  The old interval implementation called this regression
     * pair a stab for 24 windings and clean for the other 48. */
    {
        static const int perm[6][3] = {
            {0,1,2},{1,2,0},{2,0,1},{0,2,1},{2,1,0},{1,0,2}
        };
        float v[] = {
            4991.59375f,2916.64209f,2345.71240f,
            4989.46973f,2915.90869f,2347.92822f,
            4989.51221f,2916.33398f,2347.09033f,
            4991.50879f,2915.79150f,2347.38818f,
            4992.49902f,2917.72510f,2349.87964f,
            4991.42383f,2914.94092f,2349.06396f
        };
        IntersectionCleanupParams ap = params;
        int baseline = -1, invariant = 1, n_stab = 0, n_none = 0;
        ap.gap_max = 0.0;
        ap.include_hinges = 1;
        for (int swap = 0; swap < 2; swap++) {
            for (int pa = 0; pa < 6; pa++) {
                for (int pb = 0; pb < 6; pb++) {
                    int32_t a[3], b[3];
                    int shared = 0, hit;
                    for (int k = 0; k < 3; k++) {
                        a[k] = (int32_t)(perm[pa][k] + (swap ? 3 : 0));
                        b[k] = (int32_t)(perm[pb][k] + (swap ? 0 : 3));
                    }
                    hit = IntersectionCleanup_pair_test(
                        v, 6, a, b, &ap, &shared);
                    if (baseline < 0) baseline = hit;
                    if (hit != baseline || shared != 0) invariant = 0;
                    if (hit == INTERSECTION_HIT_STAB) n_stab++;
                    if (hit == 0) n_none++;
                }
            }
        }
        if (!invariant || baseline != 0)
            fprintf(stderr,
                    "[selftest] winding-invariant contact detail: baseline=%d "
                    "stab=%d none=%d\n", baseline, n_stab, n_none);
        selftest_check(invariant && baseline == 0,
                       "boundary contact is invariant to face winding/order",
                       &fails);
    }

    /* The same permutation census must retain a genuine strict-interior stab. */
    {
        static const int perm[6][3] = {
            {0,1,2},{1,2,0},{2,0,1},{0,2,1},{2,1,0},{1,0,2}
        };
        float v[] = {
            0,0,0, 4,0,0, 0,4,0,
            1,1,-2, 1,1,2, 3,1,0
        };
        IntersectionCleanupParams ap = params;
        int invariant = 1, n_stab = 0;
        ap.gap_max = 0.0;
        ap.include_hinges = 1;
        for (int swap = 0; swap < 2; swap++) {
            for (int pa = 0; pa < 6; pa++) {
                for (int pb = 0; pb < 6; pb++) {
                    int32_t a[3], b[3];
                    int shared = 0, hit;
                    for (int k = 0; k < 3; k++) {
                        a[k] = (int32_t)(perm[pa][k] + (swap ? 3 : 0));
                        b[k] = (int32_t)(perm[pb][k] + (swap ? 0 : 3));
                    }
                    hit = IntersectionCleanup_pair_test(
                        v, 6, a, b, &ap, &shared);
                    if (hit != INTERSECTION_HIT_STAB || shared != 0)
                        invariant = 0;
                    if (hit == INTERSECTION_HIT_STAB) n_stab++;
                }
            }
        }
        if (!invariant)
            fprintf(stderr,
                    "[selftest] winding-invariant true-stab detail: stab=%d/72\n",
                    n_stab);
        selftest_check(invariant,
                       "proper stab is invariant to face winding/order", &fails);
    }

    /*
     * The production broad phase must be merely an accelerator.  Compare its
     * complete masked conflict set, classifications, degrees, and candidate
     * count against exhaustive AABB + exact-predicate enumeration.
     */
    {
        enum { NFACE = 96, NVERT = NFACE * 3 };
        float v[NVERT * 3];
        int32_t f[NFACE * 3];
        uint8_t mask[NFACE], found[NFACE * NFACE], expected[NFACE * NFACE];
        size_t degree[NFACE], expected_degree[NFACE];
        FaceBox box[NFACE];
        IntersectionCleanupParams ap = params;
        IntersectionCleanupStats got, want;
        IntersectionSelftestCollector collector;
        int ok = 1;
        memset(found, 0, sizeof found);
        memset(expected, 0, sizeof expected);
        memset(degree, 0, sizeof degree);
        memset(expected_degree, 0, sizeof expected_degree);
        memset(&want, 0, sizeof want);
        ap.gap_max = 0.25;
        ap.max_conflicts = 100000;
        for (size_t q = 0; q < NFACE; q++) {
            size_t group = q / 4;
            size_t variant = q % 4;
            float bx = (float)(group % 6) * 3.0f;
            float by = (float)((group / 6) % 4) * 3.0f;
            float bz = (float)(group / 24) * 3.0f;
            int32_t i0 = (int32_t)(q * 3);
            f[q*3+0] = i0; f[q*3+1] = i0 + 1; f[q*3+2] = i0 + 2;
            v[(size_t)i0*3+0] = bx + 0.04f * (float)variant;
            v[(size_t)i0*3+1] = by;
            v[(size_t)i0*3+2] = bz + 0.03f * (float)(variant & 1);
            v[(size_t)(i0+1)*3+0] = bx + 1.8f;
            v[(size_t)(i0+1)*3+1] = by + 0.03f * (float)variant;
            v[(size_t)(i0+1)*3+2] = bz + 0.04f * (float)(variant & 2);
            v[(size_t)(i0+2)*3+0] = bx + 0.03f * (float)variant;
            v[(size_t)(i0+2)*3+1] = by + 1.8f;
            v[(size_t)(i0+2)*3+2] = bz + 0.02f * (float)((variant+1)&1);
            mask[q] = (uint8_t)((q % 11) != 0);
            box[q].face = (int32_t)q;
            box[q].morton = 0;
            for (int k = 0; k < 3; k++) {
                double a = v[(size_t)i0*3+(size_t)k];
                double b = v[(size_t)(i0+1)*3+(size_t)k];
                double c = v[(size_t)(i0+2)*3+(size_t)k];
                box[q].lo[k] = a < b ? (a < c ? a : c) :
                                         (b < c ? b : c);
                box[q].hi[k] = a > b ? (a > c ? a : c) :
                                         (b > c ? b : c);
            }
        }
        for (size_t a = 0; a < NFACE; a++) {
            if (!mask[a]) continue;
            for (size_t b = a + 1; b < NFACE; b++) {
                int shared = 0, hit;
                if (!mask[b] ||
                    !facebox_bounds_overlap(box[a].lo, box[a].hi,
                                            box[b].lo, box[b].hi,
                                            ap.gap_max))
                    continue;
                want.candidate_pairs++;
                hit = IntersectionCleanup_pair_test(
                    v, NVERT, &f[a*3], &f[b*3], &ap, &shared);
                if (hit < 0) { ok = 0; continue; }
                if (!hit) continue;
                expected[a*NFACE+b] = (uint8_t)hit;
                expected_degree[a]++;
                expected_degree[b]++;
                want.conflicts++;
                if (hit == INTERSECTION_HIT_OVERLAP)
                    want.overlap_pairs++;
                else if (hit == INTERSECTION_HIT_STAB)
                    want.stab_pairs++;
                else
                    want.fold_pairs++;
            }
        }
        for (size_t q = 0; q < NFACE; q++)
            if (expected_degree[q] > want.max_conflict_degree)
                want.max_conflict_degree = expected_degree[q];
        collector.kind = found;
        collector.n = NFACE;
        if (IntersectionCleanup_audit_visit(
                v, NVERT, f, NFACE, mask, &ap, degree,
                intersection_selftest_collect, &collector, &got) != 0)
            ok = 0;
        if (got.candidate_pairs != want.candidate_pairs ||
            got.conflicts != want.conflicts ||
            got.overlap_pairs != want.overlap_pairs ||
            got.stab_pairs != want.stab_pairs ||
            got.fold_pairs != want.fold_pairs ||
            got.max_conflict_degree != want.max_conflict_degree ||
            memcmp(degree, expected_degree, sizeof degree) != 0 ||
            memcmp(found, expected, sizeof found) != 0)
            ok = 0;
        selftest_check(ok,
                       "3-D BVH equals exhaustive exact conflict set",
                       &fails);

        /* The parallel audit must reproduce the serial audit exactly:
         * counts, classifications, per-face degrees, and visitor stream. */
        ok = 1;
        memset(found, 0, sizeof found);
        memset(degree, 0, sizeof degree);
        memset(&got, 0, sizeof got);
        if (IntersectionCleanup_audit_visit_parallel(
                v, NVERT, f, NFACE, mask, &ap, degree,
                intersection_selftest_collect, &collector, &got) != 0)
            ok = 0;
        if (got.candidate_pairs != want.candidate_pairs ||
            got.conflicts != want.conflicts ||
            got.overlap_pairs != want.overlap_pairs ||
            got.stab_pairs != want.stab_pairs ||
            got.fold_pairs != want.fold_pairs ||
            got.max_conflict_degree != want.max_conflict_degree ||
            memcmp(degree, expected_degree, sizeof degree) != 0 ||
            memcmp(found, expected, sizeof found) != 0)
            ok = 0;
        selftest_check(ok,
                       "parallel audit equals serial audit exactly",
                       &fails);
    }

    /*
     * One coarse triangle overlaps four fine triangles in a disconnected
     * component.  The high-degree coarse face is the unique one-face cover.
     */
    {
        float v[] = {
            0,0,0,  0,2,0,  0,0,2,
            0.05f,0,0,  0.05f,1,0,  0.05f,0,1,
            0.05f,2,0,  0.05f,1,1,  0.05f,0,2
        };
        int32_t f[] = {
            0,1,2,
            3,4,5, 4,6,7, 5,7,8, 4,7,5
        };
        size_t nf = 5;
        int rc = IntersectionCleanup_process(v, 9, f, &nf, NULL,
                                             &params, &stats);
        selftest_check(rc == 0 && stats.overlap_pairs >= 4,
                       "coarse/fine overlaps detected", &fails);
        selftest_check(stats.faces_deleted == 1 && nf == 4,
                       "greedy cover deletes one coarse face", &fails);
    }

    /* A kind-selective fold repair must not consume an ordinary disconnected
     * layer overlap merely because both triangles are inside the same BVH. */
    {
        float v[] = {
            0,0,0, 0,2,0, 0,0,2,
            0.05f,0,0, 0.05f,2,0, 0.05f,0,2
        };
        int32_t f[] = {0,1,2, 3,4,5};
        size_t nf = 2;
        IntersectionCleanupParams fp = params;
        fp.hit_kind_mask = INTERSECTION_HIT_MASK_FOLD;
        int rc = IntersectionCleanup_process(v, 6, f, &nf, NULL, &fp, &stats);
        selftest_check(rc == 0 && stats.conflicts == 0 &&
                       stats.faces_deleted == 0 && nf == 2,
                       "fold-only cleanup preserves layer overlap", &fails);
    }

    /* A clean two-triangle quad shares an edge but has disjoint interiors. */
    {
        float v[] = {0,0,0, 0,1,0, 0,1,1, 0,0,1};
        int32_t f[] = {0,1,2, 0,2,3};
        size_t nf = 2;
        int rc = IntersectionCleanup_process(v, 4, f, &nf, NULL,
                                             &params, &stats);
        selftest_check(rc == 0 && stats.conflicts == 0 && nf == 2,
                       "clean adjacent triangles unchanged", &fails);
    }

    /* A pinch split duplicates a shared point geometrically while assigning a
     * new index.  The non-parallel Moller path must still treat that point as
     * shared boundary, not invent a post-repair stab.  Coordinates are a
     * translated PHerc1447 regression pair. */
    {
        float v[] = {
            -0.114258f,  1.436768f,  0.373291f,
            -0.083984f, -0.685303f, -1.675537f,
             0.0f,       0.0f,       0.0f,
            -0.406250f,  0.121826f,  0.813232f,
            -0.105469f,  1.388184f,  0.342773f,
             0.0f,       0.0f,       0.0f
        };
        int32_t f[] = {0,1,2, 3,4,5};
        size_t nf = 2;
        int rc = IntersectionCleanup_process(v, 6, f, &nf, NULL,
                                             &params, &stats);
        selftest_check(rc == 0 && stats.conflicts == 0 && nf == 2,
                       "duplicate-index shared point is boundary contact", &fails);
    }

    /* Edge-adjacent fold flap: delete the smaller overlapping triangle. */
    {
        float v[] = {0,0,0, 0,4,0, 0,0,4, 0.05f,0.5f,0.5f};
        int32_t f[] = {0,1,2, 0,1,3};
        size_t nf = 2;
        IntersectionCleanupParams fp = params;
        fp.hit_kind_mask = INTERSECTION_HIT_MASK_FOLD;
        int rc = IntersectionCleanup_process(v, 4, f, &nf, NULL,
                                             &fp, &stats);
        selftest_check(rc == 0 && stats.fold_pairs == 1,
                       "shared-edge fold detected", &fails);
        selftest_check(stats.faces_deleted == 1 && nf == 1 &&
                       f[0] == 0 && f[1] == 1 && f[2] == 2,
                       "smaller fold flap removed", &fails);
    }

    /* Fail closed: the same conflict exceeds a zero deletion budget. */
    {
        float v[] = {
            0,0,0, 0,2,0, 0,0,2,
            0.05f,0,0, 0.05f,2,0, 0.05f,0,2
        };
        int32_t f[] = {0,1,2, 3,4,5};
        int32_t original[6];
        size_t nf = 2;
        memcpy(original, f, sizeof(f));
        params.max_delete_fraction = 0.0;
        params.max_delete_fraction_hard = 0.0;
        params.min_delete_budget_faces = 0;
        selftest_check(IntersectionCleanup_process(v, 6, f, &nf, NULL,
                                                   &params, &stats) != 0 &&
                       stats.budget_rejected && nf == 2 &&
                       memcmp(f, original, sizeof(f)) == 0,
                       "budget rejection leaves mesh unchanged", &fails);
    }

    fprintf(stderr, "[selftest] intersection_cleanup: %s (%d failures)\n",
            fails ? "FAIL" : "ALL PASS", fails);
    return fails;
}
