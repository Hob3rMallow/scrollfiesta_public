#include "cvt_remesh.h"
#include "rvd_clip.h"
#include "../common/closest_tri.h"
#include "../common/kdtree.h"
#include "../common/cycle_timer.h"
#include "../common/pipeline_constants.h"

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "../common/ves_omp.h"

/* Per-thread Lloyd scratch: the projection candidate buffer and instrument
 * counters (summed into CvtPhaseStats after each pass). Padded so threads
 * never share a cache line. */
typedef struct {
    int32_t  cand[1024];
    double   lbs[1024];       /* projection lower bounds (project_to_mesh_kd_lb) */
    uint64_t tk_bnd, tk_proj;
    uint64_t n_bnd_calls, n_seg_scans, n_proj_calls, n_proj_cand, n_proj_cap;
    char     pad[64];
} LloydThread;

void CVT_default_opts(CvtOpts *o) {
    if (!o) return;
    o->n_iters        = 50;
    o->lambda_na      = 1.0f;
    o->lambda_cvt0    = 1.0f;
    o->cvt_decay      = 0.95f;
    o->lambda_cvt_min = 1e-3f;
    o->seed           = 0;
    o->verbose        = 0;
    o->stats          = NULL;
    o->probe          = NULL;
    o->probe_ctx      = NULL;
    o->probe_mask     = 0;
}

void CVT_phase_stats_add(CvtPhaseStats *d, const CvtPhaseStats *s) {
    if (!d || !s) return;
    d->tk_setup += s->tk_setup;  d->tk_bext += s->tk_bext;  d->tk_seed += s->tk_seed;
    d->tk_rvd += s->tk_rvd;      d->tk_bnd += s->tk_bnd;    d->tk_proj += s->tk_proj;
    d->tk_dual += s->tk_dual;
    d->n_calls += s->n_calls;    d->n_sites += s->n_sites;  d->n_iters += s->n_iters;
    d->n_bnd_segs += s->n_bnd_segs;     d->n_seed_tries += s->n_seed_tries;
    d->n_bnd_calls += s->n_bnd_calls;   d->n_seg_scans += s->n_seg_scans;
    d->n_proj_calls += s->n_proj_calls; d->n_proj_cand += s->n_proj_cand;
    d->n_proj_cap += s->n_proj_cap;
    d->rvd.tk_near += s->rvd.tk_near;   d->rvd.tk_cand += s->rvd.tk_cand;
    d->rvd.tk_clip += s->rvd.tk_clip;   d->rvd.n_passes += s->rvd.n_passes;
    d->rvd.n_faces += s->rvd.n_faces;   d->rvd.n_fast += s->rvd.n_fast;
    d->rvd.n_cand += s->rvd.n_cand;     d->rvd.n_sat += s->rvd.n_sat;
    d->rvd.n_wide += s->rvd.n_wide;
    d->rvd.n_cells += s->rvd.n_cells;   d->rvd.n_planes += s->rvd.n_planes;
    d->rvd.n_preds += s->rvd.n_preds;
    d->tk_probe += s->tk_probe;  d->n_probe += s->n_probe;
}

/* Solve a symmetric 3x3 system A x = b, A stored as {a00,a01,a02,a11,a12,a22}.
 * A is positive-definite here (lambda_CVT*m*I regularizes N), so the cofactor
 * inverse is stable; falls back to a diagonal solve if the determinant underflows. */
static void sym3_solve(const double A[6], const double b[3], double x[3]) {
    double a0=A[0],a1=A[1],a2=A[2],a3=A[3],a4=A[4],a5=A[5];
    double c00=a3*a5-a4*a4, c01=a2*a4-a1*a5, c02=a1*a4-a2*a3;
    double det=a0*c00+a1*c01+a2*c02;
    if (fabs(det) < 1e-20) {
        x[0]= a0>0?b[0]/a0:b[0]; x[1]= a3>0?b[1]/a3:b[1]; x[2]= a5>0?b[2]/a5:b[2];
        return;
    }
    double c11=a0*a5-a2*a2, c12=a1*a2-a0*a4, c22=a0*a3-a1*a1, inv=1.0/det;
    x[0]=inv*(c00*b[0]+c01*b[1]+c02*b[2]);
    x[1]=inv*(c01*b[0]+c11*b[1]+c12*b[2]);
    x[2]=inv*(c02*b[0]+c12*b[1]+c22*b[2]);
}

/* xorshift32 -> uniform [0,1) */
static uint32_t xs32(uint32_t *s) { uint32_t x=*s; x^=x<<13; x^=x>>17; x^=x<<5; *s=x; return x; }
static double   rnd01(uint32_t *s) { return (double)(xs32(s) >> 8) * (1.0/16777216.0); }

/* Closest point on the mesh to p (brute force over triangles — M1/M2). */
static void project_to_mesh(const double *V, const int32_t *F, size_t nf,
                            const double p[3], double out[3]) {
    double best = 1e300;
    out[0]=p[0]; out[1]=p[1]; out[2]=p[2];
    for (size_t f = 0; f < nf; f++) {
        const double *A=&V[(size_t)F[f*3+0]*3];
        const double *B=&V[(size_t)F[f*3+1]*3];
        const double *C=&V[(size_t)F[f*3+2]*3];
        double q[3], u, v;
        ClosestTri_point(p, A, B, C, q, &u, &v);
        double dz=p[0]-q[0], dy=p[1]-q[1], dx=p[2]-q[2];
        double d2 = dz*dz+dy*dy+dx*dx;
        if (d2 < best) { best=d2; out[0]=q[0]; out[1]=q[1]; out[2]=q[2]; }
    }
}

/* Closest point on the mesh via a triangle-centroid KD-tree. Any triangle that
 * could be the closest has its centroid within (nearest-centroid distance +
 * 2*max_tri_r) of p, so one ball query returns all candidates. Replaces the O(nf)
 * brute force -- the dominant cost on dense input meshes. */
static size_t project_to_mesh_kd(KDTree_T kd, const double *V, const int32_t *F,
                                 double max_tri_r, const double p[3], double out[3],
                                 int32_t *cand, size_t cand_cap) {
    float pf[3] = {(float)p[0],(float)p[1],(float)p[2]};
    float d0sq = 0.0f;
    size_t i0 = KDTree_nearest(kd, pf, &d0sq);
    double R = sqrt((double)d0sq) + 2.0*max_tri_r + 1e-6;
    size_t nc = KDTree_ball_query(kd, pf, (float)(R*R), cand, cand_cap);
    if (nc == 0) { cand[0] = (int32_t)i0; nc = 1; }
    double best = 1e300;
    out[0]=p[0]; out[1]=p[1]; out[2]=p[2];
    for (size_t k = 0; k < nc; k++) {
        size_t f = (size_t)cand[k];
        const double *A=&V[(size_t)F[f*3+0]*3], *B=&V[(size_t)F[f*3+1]*3], *C=&V[(size_t)F[f*3+2]*3];
        double q[3], u, v;
        ClosestTri_point(p, A, B, C, q, &u, &v);
        double dz=p[0]-q[0], dy=p[1]-q[1], dx=p[2]-q[2], d2=dz*dz+dy*dy+dx*dx;
        if (d2 < best) { best=d2; out[0]=q[0]; out[1]=q[1]; out[2]=q[2]; }
    }
    return nc;
}

/* project_to_mesh_kd with the same candidate list and the same answer -- the
 * first minimum of d2 in candidate order -- evaluating only the candidates
 * that can still tie or beat the best (2026-10-05).  Every point of triangle f
 * lies within tri_cr[f*4+3] of its centroid tri_cr[f*4+0..2], so
 * (|p - c| - r - 1e-6)^2 bounds its d2 from below (the 1e-6 vox margin is far
 * above double rounding at these coordinates).  The most promising candidate
 * seeds the bound; one pass in candidate order then evaluates the rest whose
 * bound does not exceed it, and equal d2 keep the earliest position, as the
 * original strict '<' scan does.  CVT_projection_selftest compares the two
 * bitwise, including truncated (capped) candidate lists. */
static size_t project_to_mesh_kd_lb(KDTree_T kd, const double *V, const int32_t *F,
                                    double max_tri_r, const double *tri_cr,
                                    const double p[3], double out[3],
                                    int32_t *cand, size_t cand_cap, double *lbs) {
    float pf[3] = {(float)p[0],(float)p[1],(float)p[2]};
    float d0sq = 0.0f;
    size_t i0 = KDTree_nearest(kd, pf, &d0sq);
    double R = sqrt((double)d0sq) + 2.0*max_tri_r + 1e-6;
    size_t nc = KDTree_ball_query(kd, pf, (float)(R*R), cand, cand_cap);
    if (nc == 0) { cand[0] = (int32_t)i0; nc = 1; }
    out[0]=p[0]; out[1]=p[1]; out[2]=p[2];
    size_t kmin = 0;
    double lmin = HUGE_VAL;
    for (size_t k = 0; k < nc; k++) {
        const double *cr = &tri_cr[(size_t)cand[k]*4];
        double dz=p[0]-cr[0], dy=p[1]-cr[1], dx=p[2]-cr[2];
        double lb = sqrt(dz*dz+dy*dy+dx*dx) - cr[3] - 1e-6;
        double lb2 = lb > 0.0 ? lb*lb : 0.0;
        lbs[k] = lb2;
        if (lb2 < lmin) { lmin = lb2; kmin = k; }
    }
    double best = 1e300;
    size_t best_k = 0;
    int have = 0;
    for (size_t step = 0; step <= nc; step++) {
        size_t k;
        if (step == 0) k = kmin;                  /* seed the bound */
        else { k = step - 1; if (k == kmin) continue; }
        if (lbs[k] > best) continue;              /* cannot tie or beat (NaN: evaluate) */
        size_t f = (size_t)cand[k];
        const double *A=&V[(size_t)F[f*3+0]*3], *B=&V[(size_t)F[f*3+1]*3], *C=&V[(size_t)F[f*3+2]*3];
        double q[3], u, v;
        ClosestTri_point(p, A, B, C, q, &u, &v);
        double dz=p[0]-q[0], dy=p[1]-q[1], dx=p[2]-q[2], d2=dz*dz+dy*dy+dx*dx;
        if (d2 < best || (have && d2 == best && k < best_k)) {
            best=d2; best_k=k; have=1;
            out[0]=q[0]; out[1]=q[1]; out[2]=q[2];
        }
    }
    return nc;
}

static double closest_on_seg(const double p[3], const double a[3], const double b[3], double out[3]) {
    double ab[3]={b[0]-a[0],b[1]-a[1],b[2]-a[2]};
    double t = (p[0]-a[0])*ab[0]+(p[1]-a[1])*ab[1]+(p[2]-a[2])*ab[2];
    double d = ab[0]*ab[0]+ab[1]*ab[1]+ab[2]*ab[2];
    t = d>0.0 ? t/d : 0.0; if (t<0) t=0; if (t>1) t=1;
    out[0]=a[0]+t*ab[0]; out[1]=a[1]+t*ab[1]; out[2]=a[2]+t*ab[2];
    double dz=p[0]-out[0], dy=p[1]-out[1], dx=p[2]-out[2];
    return dz*dz+dy*dy+dx*dx;
}

/* Closest point on the union of boundary segments (seg[nseg*6] = a.xyz,b.xyz).
 * Reference scan; production queries go through seg_bvh_project below. */
static void project_to_boundary(const double *seg, size_t nseg, const double p[3], double out[3]) {
    double best=1e300;
    out[0]=p[0]; out[1]=p[1]; out[2]=p[2];
    for (size_t s=0; s<nseg; s++) {
        double q[3];
        double d2 = closest_on_seg(p, &seg[s*6], &seg[s*6+3], q);
        if (d2<best) { best=d2; out[0]=q[0]; out[1]=q[1]; out[2]=q[2]; }
    }
}

/* project_to_boundary + the index of the nearest segment (for LOCAL spacing). */
static void project_to_boundary_seg(const double *seg, size_t nseg,
                                    const double p[3], double out[3],
                                    size_t *out_seg) {
    double best = 1e300;
    out[0]=p[0]; out[1]=p[1]; out[2]=p[2];
    *out_seg = 0;
    for (size_t s=0; s<nseg; s++) {
        double q[3];
        double d2 = closest_on_seg(p, &seg[s*6], &seg[s*6+3], q);
        if (d2<best) { best=d2; out[0]=q[0]; out[1]=q[1]; out[2]=q[2]; *out_seg=s; }
    }
}

/* ---- Boundary segment BVH: exact replacement for the brute-force scans ----
 *
 * project_to_boundary(_seg) return the FIRST segment (lowest index) attaining
 * the minimum computed d2 = closest_on_seg(p, a, b, q), and that segment's q.
 * Ragged BPA charts carry ~0.17 boundary segments per face (22,222 on a 136k-face
 * chart), and the brute-force scan ran up to 9x per seed plus once per site per
 * Lloyd iteration -- a measured CVT hot spot (2026-10-04).
 *
 * The BVH returns the identical segment and bitwise-identical q:
 *  - leaves call closest_on_seg with the same pointers, and a segment is taken
 *    when d2 < best || (d2 == best && s < best_s): the lowest index among the
 *    minimum computed d2, exactly what the ascending strict-'<' scan keeps;
 *  - a node is skipped only when a lower bound on the COMPUTED d2 of every
 *    segment under it exceeds the current best, so no winner or tie is lost.
 *    Node boxes are the exact endpoint min/max padded by 2^-30 (1+M); the
 *    computed q of a clamped segment parameter lies within ~5u*M of its box
 *    (u = 2^-53), and the box distance is shrunk by 2^-40 to cover its own
 *    rounding, both far below the margins (M = max |coordinate|);
 *  - non-finite data, or a query beyond 4 (1+M), falls back to the scan. */
#define SEG_BVH_LEAF      4
#define SEG_BVH_STACK     128
#define SEG_BVH_PAD_REL   9.313225746154785e-10   /* 2^-30 */
#define SEG_BVH_LB_SHRINK (1.0 - 9.094947017729282e-13)   /* 1 - 2^-40 */

typedef struct {
    double  lo[3], hi[3];    /* padded box of the segments below */
    int32_t left, right;     /* children; -1 for a leaf */
    int32_t first, count;    /* leaf: range in idx[] */
} SegBvhNode;

typedef struct {
    SegBvhNode   *nodes;
    int32_t      *idx;       /* segment ids, leaf order */
    const double *seg;       /* [nseg*6], not owned */
    size_t        nseg;
    int32_t       n_nodes;
    double        pad, lim;
    int           ok;        /* 0 => always scan */
} SegBvh;

static double seg_cen(const double *seg, int32_t s, int ax) {
    return 0.5 * (seg[(size_t)s*6 + (size_t)ax] + seg[(size_t)s*6 + 3 + (size_t)ax]);
}

/* strict total order: centroid on `ax`, then segment id */
static int seg_less(const double *seg, int ax, int32_t a, int32_t b) {
    double ca = seg_cen(seg, a, ax), cb = seg_cen(seg, b, ax);
    if (ca < cb) return 1;
    if (ca > cb) return 0;
    return a < b;
}

/* Quickselect: afterwards v[k] is the k-th of v[lo,hi) in seg_less order. */
static void seg_select(const double *seg, int ax, int32_t *v,
                       int32_t lo, int32_t hi, int32_t k) {
    while (hi - lo > 1) {
        int32_t mid = lo + (hi - lo) / 2, t;
        int32_t piv = v[mid];
        v[mid] = v[hi-1]; v[hi-1] = piv;
        int32_t store = lo;
        for (int32_t i = lo; i < hi - 1; i++)
            if (seg_less(seg, ax, v[i], piv)) { t = v[i]; v[i] = v[store]; v[store] = t; store++; }
        t = v[store]; v[store] = v[hi-1]; v[hi-1] = t;
        if (k == store) return;
        if (k < store) hi = store; else lo = store + 1;
    }
}

static int32_t seg_bvh_build_rec(SegBvh *b, int32_t lo, int32_t hi) {
    int32_t id = b->n_nodes++;
    SegBvhNode *nd = &b->nodes[id];
    double bl[3] = { 1e300, 1e300, 1e300 }, bh[3] = { -1e300, -1e300, -1e300 };
    double cl[3] = { 1e300, 1e300, 1e300 }, ch[3] = { -1e300, -1e300, -1e300 };
    for (int32_t i = lo; i < hi; i++) {
        const double *g = &b->seg[(size_t)b->idx[i] * 6];
        for (int k = 0; k < 3; k++) {
            double a0 = g[k], a1 = g[3 + k];
            double mn = a0 < a1 ? a0 : a1, mx = a0 < a1 ? a1 : a0;
            double c = 0.5 * (a0 + a1);
            if (mn < bl[k]) bl[k] = mn;
            if (mx > bh[k]) bh[k] = mx;
            if (c < cl[k]) cl[k] = c;
            if (c > ch[k]) ch[k] = c;
        }
    }
    for (int k = 0; k < 3; k++) { nd->lo[k] = bl[k] - b->pad; nd->hi[k] = bh[k] + b->pad; }
    nd->first = lo; nd->count = hi - lo;
    nd->left = nd->right = -1;
    if (hi - lo <= SEG_BVH_LEAF) return id;
    int ax = 0;
    for (int k = 1; k < 3; k++) if (ch[k] - cl[k] > ch[ax] - cl[ax]) ax = k;
    int32_t mid = lo + (hi - lo) / 2;
    seg_select(b->seg, ax, b->idx, lo, hi, mid);
    int32_t l = seg_bvh_build_rec(b, lo, mid);
    int32_t r = seg_bvh_build_rec(b, mid, hi);
    b->nodes[id].left = l;
    b->nodes[id].right = r;
    return id;
}

static void seg_bvh_build(Arena_T arena, const double *seg, size_t nseg, SegBvh *b) {
    memset(b, 0, sizeof(*b));
    b->seg = seg; b->nseg = nseg;
    if (nseg == 0 || nseg > (size_t)INT32_MAX / 4) return;
    double M = 0.0;
    for (size_t i = 0; i < nseg * 6; i++) {
        if (!isfinite(seg[i])) return;
        double a = fabs(seg[i]);
        if (a > M) M = a;
    }
    b->pad = SEG_BVH_PAD_REL * (1.0 + M);
    b->lim = 4.0 * (1.0 + M);
    b->idx = (int32_t *)ARENA_ALLOC(arena, (size_t)(nseg * sizeof(int32_t)));
    for (size_t i = 0; i < nseg; i++) b->idx[i] = (int32_t)i;
    b->nodes = (SegBvhNode *)ARENA_ALLOC(arena, (size_t)((2 * nseg + 1) * sizeof(SegBvhNode)));
    b->n_nodes = 0;
    (void)seg_bvh_build_rec(b, 0, (int32_t)nseg);
    b->ok = 1;
}

/* Lower bound on the computed closest_on_seg d2 of any segment under nd. */
static double seg_box_lb(const SegBvhNode *nd, const double p[3]) {
    double s = 0.0;
    for (int k = 0; k < 3; k++) {
        double g = 0.0;
        if (p[k] < nd->lo[k]) g = nd->lo[k] - p[k];
        else if (p[k] > nd->hi[k]) g = p[k] - nd->hi[k];
        s += g * g;
    }
    return s * SEG_BVH_LB_SHRINK;
}

/* Same result as project_to_boundary_seg(b->seg, b->nseg, p, out, out_seg);
 * returns the number of segments examined (instrumentation). */
static size_t seg_bvh_project(const SegBvh *b, const double p[3], double out[3],
                              size_t *out_seg) {
    size_t s0 = 0;
    if (!b->ok || !(fabs(p[0]) <= b->lim && fabs(p[1]) <= b->lim &&
                    fabs(p[2]) <= b->lim)) {
        project_to_boundary_seg(b->seg, b->nseg, p, out, &s0);
        if (out_seg) *out_seg = s0;
        return b->nseg;
    }
    double best = 1e300;
    int32_t best_s = -1;
    double bq[3] = { p[0], p[1], p[2] };
    int32_t stack[SEG_BVH_STACK];
    int sp = 0;
    size_t evals = 0;
    stack[sp++] = 0;
    while (sp > 0) {
        const SegBvhNode *nd = &b->nodes[stack[--sp]];
        if (seg_box_lb(nd, p) > best) continue;
        if (nd->left < 0) {
            for (int32_t i = nd->first; i < nd->first + nd->count; i++) {
                int32_t s = b->idx[i];
                double q[3];
                double d2 = closest_on_seg(p, &b->seg[(size_t)s*6], &b->seg[(size_t)s*6+3], q);
                evals++;
                if (d2 < best || (d2 == best && s < best_s)) {
                    best = d2; best_s = s; bq[0] = q[0]; bq[1] = q[1]; bq[2] = q[2];
                }
            }
            continue;
        }
        if (sp + 2 > SEG_BVH_STACK) {          /* cannot happen for a median tree */
            project_to_boundary_seg(b->seg, b->nseg, p, out, &s0);
            if (out_seg) *out_seg = s0;
            return evals + b->nseg;
        }
        double ll = seg_box_lb(&b->nodes[nd->left], p);
        double lr = seg_box_lb(&b->nodes[nd->right], p);
        if (ll <= lr) { stack[sp++] = nd->right; stack[sp++] = nd->left; }
        else          { stack[sp++] = nd->left;  stack[sp++] = nd->right; }
    }
    if (best_s < 0) {                          /* nothing below 1e300: as the scan */
        out[0] = p[0]; out[1] = p[1]; out[2] = p[2];
        if (out_seg) *out_seg = 0;
        return evals;
    }
    out[0] = bq[0]; out[1] = bq[1]; out[2] = bq[2];
    if (out_seg) *out_seg = (size_t)best_s;
    return evals;
}


typedef struct {
    int32_t *loopv;    /* boundary vertex indices, concatenated per loop */
    int32_t *loopoff;  /* [nloops+1] offsets into loopv */
    int      nloops;
    double  *seg;      /* [nseg*6] boundary segment endpoint positions */
    size_t   nseg;
} Boundary;

static int cmp_edge(const void *pa, const void *pb) {
    const int32_t *a=(const int32_t*)pa, *b=(const int32_t*)pb;
    if (a[0]!=b[0]) return a[0]<b[0]?-1:1;
    if (a[1]!=b[1]) return a[1]<b[1]?-1:1;
    return 0;
}

/* Extract boundary loops (edges incident to exactly one face) and chain them.
 * Assumes a manifold boundary (each boundary vertex has degree 2). */
static void extract_boundary(Arena_T a, const double *V, const int32_t *F,
                             size_t nf, size_t nv, Boundary *bd) {
    memset(bd, 0, sizeof(*bd));
    size_t ne = nf*3;
    int32_t *E = (int32_t*)ARENA_ALLOC(a, (size_t)(ne*2*sizeof(int32_t)));
    for (size_t f=0; f<nf; f++)
        for (int e=0; e<3; e++) {
            int32_t u=F[f*3+e], w=F[f*3+(e+1)%3];
            int32_t lo=u<w?u:w, hi=u<w?w:u;
            E[(f*3+e)*2]=lo; E[(f*3+e)*2+1]=hi;
        }
    qsort(E, ne, 2*sizeof(int32_t), cmp_edge);

    /* boundary edges = runs of length 1 */
    int32_t *nbr0=(int32_t*)ARENA_ALLOC(a,(size_t)(nv*sizeof(int32_t)));
    int32_t *nbr1=(int32_t*)ARENA_ALLOC(a,(size_t)(nv*sizeof(int32_t)));
    for (size_t i=0;i<nv;i++){ nbr0[i]=-1; nbr1[i]=-1; }
    size_t nbnd=0;
    /* first pass: count boundary edges */
    for (size_t i=0;i<ne;) {
        size_t j=i+1;
        while (j<ne && E[j*2]==E[i*2] && E[j*2+1]==E[i*2+1]) j++;
        if (j-i==1) nbnd++;
        i=j;
    }
    bd->nseg=nbnd;
    bd->seg=(double*)ARENA_ALLOC(a,(size_t)(nbnd*6*sizeof(double)));
    size_t si=0;
    for (size_t i=0;i<ne;) {
        size_t j=i+1;
        while (j<ne && E[j*2]==E[i*2] && E[j*2+1]==E[i*2+1]) j++;
        if (j-i==1) {
            int32_t u=E[i*2], w=E[i*2+1];
            const double *Pu=&V[(size_t)u*3], *Pw=&V[(size_t)w*3];
            bd->seg[si*6+0]=Pu[0]; bd->seg[si*6+1]=Pu[1]; bd->seg[si*6+2]=Pu[2];
            bd->seg[si*6+3]=Pw[0]; bd->seg[si*6+4]=Pw[1]; bd->seg[si*6+5]=Pw[2];
            si++;
            /* adjacency for chaining */
            if (nbr0[u]<0) nbr0[u]=w; else nbr1[u]=w;
            if (nbr0[w]<0) nbr0[w]=u; else nbr1[w]=u;
        }
        i=j;
    }

    /* chain into loops */
    uint8_t *vis=(uint8_t*)ARENA_CALLOC(a,(size_t)nv,1);
    bd->loopv=(int32_t*)ARENA_ALLOC(a,(size_t)((nbnd+1)*sizeof(int32_t)));
    bd->loopoff=(int32_t*)ARENA_ALLOC(a,(size_t)((nbnd+2)*sizeof(int32_t)));
    int nloops=0; int32_t w=0; bd->loopoff[0]=0;
    for (size_t s=0; s<nv; s++) {
        if (nbr0[s]<0 || vis[s]) continue;
        int32_t prev=-1, cur=(int32_t)s;
        while (cur>=0 && !vis[cur]) {
            vis[cur]=1; bd->loopv[w++]=cur;
            int32_t nx = (nbr0[cur]!=prev) ? nbr0[cur] : nbr1[cur];
            prev=cur; cur=nx;
        }
        bd->loopoff[++nloops]=w;
    }
    bd->nloops=nloops;
}

/* Place m points evenly by arc length along a loop of vertices; writes to out[m*3]. */
static void place_on_loop(const double *V, const int32_t *idx, int n, int m, double *out) {
    if (n < 2 || m < 1) return;
    double L=0.0;
    for (int i=0;i<n;i++) {
        const double *a=&V[(size_t)idx[i]*3], *b=&V[(size_t)idx[(i+1)%n]*3];
        double dz=b[0]-a[0],dy=b[1]-a[1],dx=b[2]-a[2];
        L += sqrt(dz*dz+dy*dy+dx*dx);
    }
    if (L<=0.0) { for (int k=0;k<m;k++){ const double*p=&V[(size_t)idx[0]*3]; out[k*3]=p[0];out[k*3+1]=p[1];out[k*3+2]=p[2]; } return; }
    double step=L/m, target=0.0, acc=0.0; int seg=0;
    const double *a=&V[(size_t)idx[0]*3], *b=&V[(size_t)idx[1%n]*3];
    double dz=b[0]-a[0],dy=b[1]-a[1],dx=b[2]-a[2];
    double seglen=sqrt(dz*dz+dy*dy+dx*dx);
    for (int k=0;k<m;k++) {
        while (acc+seglen < target && seg < n-1) {
            acc+=seglen; seg++;
            a=&V[(size_t)idx[seg]*3]; b=&V[(size_t)idx[(seg+1)%n]*3];
            dz=b[0]-a[0];dy=b[1]-a[1];dx=b[2]-a[2]; seglen=sqrt(dz*dz+dy*dy+dx*dx);
        }
        double t = seglen>0 ? (target-acc)/seglen : 0.0; if (t<0)t=0; if (t>1)t=1;
        out[k*3+0]=a[0]+t*dz; out[k*3+1]=a[1]+t*dy; out[k*3+2]=a[2]+t*dx;
        target+=step;
    }
}

/* Per-chart CVT setup: everything a CVT attempt needs that does not depend on
 * the site count, so the certificate ladder (pipeline_cube.c) builds it once
 * per chart instead of once per attempt (2026-10-04).  Arena-allocated; take
 * per-attempt arena marks AFTER it.  No CVT-level Arena_save/restore (avoids
 * nested-mark aliasing with the dual `tris` Rvd_accumulate returns). */
struct CvtPrep_T {
    const int32_t  *faces;    /* caller's faces, not copied */
    size_t          nv, nf;
    const CvtField *field;
    double         *V;        /* [nv*3] double copy of the vertices */
    Rvd_T           r;
    double          surf_area;
    double          max_tri_r;
    KDTree_T        proj_kd;  /* triangle centroids (float) */
    double         *tri_cr;   /* [nf*4] centroid (double) + centroid-to-vertex radius */
    double         *cum;      /* seeding CDF over faces */
    double          tot;
    Boundary        bd;
    SegBvh          bvh;
};

CvtPrep_T CVT_prepare(Arena_T arena,
                      const float *in_verts, size_t in_nv,
                      const int32_t *in_faces, size_t in_nf,
                      const CvtField *field, CvtPhaseStats *st, int *out_rc) {
    if (out_rc) *out_rc = 0;
    if (!arena || !in_verts || !in_faces || in_nv < 3 || in_nf < 1) {
        if (out_rc) *out_rc = -1;
        return NULL;
    }
    uint64_t tk0 = st ? cycle_now() : 0;
    CvtPrep_T p = (CvtPrep_T)ARENA_CALLOC(arena, 1, sizeof(*p));
    p->faces = in_faces; p->nv = in_nv; p->nf = in_nf; p->field = field;

    double *V = (double *)ARENA_ALLOC(arena, (size_t)(in_nv*3*sizeof(double)));
    for (size_t i = 0; i < in_nv*3; i++) V[i] = (double)in_verts[i];
    p->V = V;

    p->r = Rvd_new(arena, V, in_nv, in_faces, in_nf);
    if (!p->r) { if (out_rc) *out_rc = -2; return NULL; }
    p->surf_area = Rvd_mesh_area(p->r);
    if (p->surf_area <= 0.0) { if (out_rc) *out_rc = -3; return NULL; }

    /* Projection acceleration: KD-tree over triangle centroids (+ the max
     * centroid-to-vertex radius, so the ball query is a provable superset). */
    float *tc = (float *)ARENA_ALLOC(arena, (size_t)(in_nf*3*sizeof(float)));
    double *tri_cr = (double *)ARENA_ALLOC(arena, (size_t)(in_nf*4*sizeof(double)));
    double max_tri_r = 0.0;
    for (size_t f = 0; f < in_nf; f++) {
        const double *A=&V[(size_t)in_faces[f*3+0]*3];
        const double *B=&V[(size_t)in_faces[f*3+1]*3];
        const double *C=&V[(size_t)in_faces[f*3+2]*3];
        double c0=(A[0]+B[0]+C[0])/3.0, c1=(A[1]+B[1]+C[1])/3.0, c2=(A[2]+B[2]+C[2])/3.0;
        tc[f*3+0]=(float)c0; tc[f*3+1]=(float)c1; tc[f*3+2]=(float)c2;
        double ra=sqrt((A[0]-c0)*(A[0]-c0)+(A[1]-c1)*(A[1]-c1)+(A[2]-c2)*(A[2]-c2));
        double rb=sqrt((B[0]-c0)*(B[0]-c0)+(B[1]-c1)*(B[1]-c1)+(B[2]-c2)*(B[2]-c2));
        double rc2=sqrt((C[0]-c0)*(C[0]-c0)+(C[1]-c1)*(C[1]-c1)+(C[2]-c2)*(C[2]-c2));
        double rr=ra>rb?(ra>rc2?ra:rc2):(rb>rc2?rb:rc2);
        if (rr>max_tri_r) max_tri_r=rr;
        tri_cr[f*4+0]=c0; tri_cr[f*4+1]=c1; tri_cr[f*4+2]=c2; tri_cr[f*4+3]=rr;
    }
    p->max_tri_r = max_tri_r;
    p->tri_cr = tri_cr;
    p->proj_kd = KDTree_new(arena, tc, in_nf);

    /* Seeding CDF over triangles. Uniform (field==NULL): weight = area, so sampling
     * is area-weighted. Graded (field!=NULL): weight = area / h(centroid)^2, so seed
     * density ~ 1/h^2 (spacing ~ h); the weighted total then DERIVES the site count
     * N = round(2 * sum(area/h^2)) -- the generalization of the uniform N = 2A/h^2. */
    double *cum = (double *)ARENA_ALLOC(arena, (size_t)(in_nf*sizeof(double)));
    double tot = 0.0;
    for (size_t f = 0; f < in_nf; f++) {
        const double *A=&V[(size_t)in_faces[f*3+0]*3];
        const double *B=&V[(size_t)in_faces[f*3+1]*3];
        const double *C=&V[(size_t)in_faces[f*3+2]*3];
        double ab[3]={B[0]-A[0],B[1]-A[1],B[2]-A[2]};
        double ac[3]={C[0]-A[0],C[1]-A[1],C[2]-A[2]};
        double cr[3]={ab[1]*ac[2]-ab[2]*ac[1], ab[2]*ac[0]-ab[0]*ac[2], ab[0]*ac[1]-ab[1]*ac[0]};
        double af = 0.5*sqrt(cr[0]*cr[0]+cr[1]*cr[1]+cr[2]*cr[2]);
        if (field) {
            double c0=(A[0]+B[0]+C[0])/3.0, c1=(A[1]+B[1]+C[1])/3.0, c2=(A[2]+B[2]+C[2])/3.0;
            af *= sizing_seed_rho(field, c0, c1, c2);
        }
        tot += af;
        cum[f] = tot;
    }
    p->cum = cum; p->tot = tot;

    if (st) { uint64_t tk = cycle_now(); st->tk_setup += tk - tk0; tk0 = tk; }

    /* boundary loops + their segment BVH */
    extract_boundary(arena, V, in_faces, in_nf, in_nv, &p->bd);
    seg_bvh_build(arena, p->bd.seg, p->bd.nseg, &p->bvh);
    if (st) st->tk_bext += cycle_now() - tk0;
    return p;
}

int CVT_remesh(Arena_T arena,
               const float *in_verts, size_t in_nv,
               const int32_t *in_faces, size_t in_nf,
               size_t target_nsites,
               const CvtOpts *opts,
               const CvtField *field,
               float **out_verts, size_t *out_nv,
               int32_t **out_faces, size_t *out_nf) {
    if (!arena || !in_verts || !in_faces || in_nv < 3 || in_nf < 1 || target_nsites < 3)
        return -1;
    if (out_verts) *out_verts = NULL;
    if (out_faces) *out_faces = NULL;
    if (out_nv) *out_nv = 0;
    if (out_nf) *out_nf = 0;
    int rc = 0;
    CvtPrep_T p = CVT_prepare(arena, in_verts, in_nv, in_faces, in_nf, field,
                              opts ? opts->stats : NULL, &rc);
    if (!p) return rc;
    return CVT_remesh_prepared(arena, p, target_nsites, opts,
                               out_verts, out_nv, out_faces, out_nf);
}

int CVT_remesh_prepared(Arena_T arena, CvtPrep_T p,
                        size_t target_nsites,
                        const CvtOpts *opts,
                        float **out_verts, size_t *out_nv,
                        int32_t **out_faces, size_t *out_nf) {
    if (!arena || !p || target_nsites < 3)
        return -1;
    if (out_verts) *out_verts = NULL;
    if (out_faces) *out_faces = NULL;
    if (out_nv) *out_nv = 0;
    if (out_nf) *out_nf = 0;

    CvtOpts o; if (opts) o = *opts; else CVT_default_opts(&o);
    uint32_t rng = o.seed ? o.seed : 12345u;
    size_t N = target_nsites;
    CvtPhaseStats *st = o.stats;           /* instrumentation only */
    uint64_t tk0 = st ? cycle_now() : 0;

    const double   *V = p->V;
    const int32_t  *in_faces = p->faces;
    size_t          in_nf = p->nf;
    const CvtField *field = p->field;
    Rvd_T           r = p->r;
    double          surf_area = p->surf_area;
    double          max_tri_r = p->max_tri_r;
    KDTree_T        proj_kd = p->proj_kd;
    const double   *tri_cr = p->tri_cr;
    const double   *cum = p->cum;
    double          tot = p->tot;
    Boundary        bd = p->bd;
    SegBvh          bvh = p->bvh;
    if (st) Rvd_set_stats(r, &st->rvd);

    if (field) {
        long long nn = (long long)(2.0*tot + 0.5);   /* derived site count */
        if (nn < 4) nn = 4;
        if ((size_t)nn > in_nf) nn = (long long)in_nf;
        N = (size_t)nn;
    }
    double spacing = sqrt(2.0*surf_area/(double)N);   /* approx target edge length */
    if (st) { st->n_calls++; st->n_sites += N; st->n_bnd_segs += bd.nseg; }

    double *sites = (double *)ARENA_ALLOC(arena, (size_t)(N*3*sizeof(double)));
    uint8_t *bflag = (uint8_t *)ARENA_CALLOC(arena, (size_t)N, 1);
    size_t ns = 0;

    /* boundary seeds first (constrained to their loop) */
    for (int k = 0; k < bd.nloops && ns < N; k++) {
        int lo = bd.loopoff[k], hi = bd.loopoff[k+1], n = hi - lo;
        if (n < 2) continue;
        double L=0.0, hc=0.0;
        for (int i=0;i<n;i++){
            const double *aa=&V[(size_t)bd.loopv[lo+i]*3], *bb=&V[(size_t)bd.loopv[lo+(i+1)%n]*3];
            double dz=bb[0]-aa[0],dy=bb[1]-aa[1],dx=bb[2]-aa[2];
            double sl=sqrt(dz*dz+dy*dy+dx*dx); L+=sl;
            if (field) {
                double mz=(aa[0]+bb[0])*0.5, my=(aa[1]+bb[1])*0.5, mx=(aa[2]+bb[2])*0.5;
                hc += sl / sizing_h(field, mz, my, mx);   /* local-h arc-length count */
            }
        }
        /* Count by LOCAL edge length: seam rims (h==h_seam) get dense verts the weld
         * bridge can span; interior-hole rims (h==h_interior) stay coarse. */
        int m = field ? (int)(hc + 0.5) : (int)(L/spacing + 0.5);
        if (m < 3) m = 3;
        if ((size_t)m > N - ns) m = (int)(N - ns);
        if (m < 1) break;
        place_on_loop(V, &bd.loopv[lo], n, m, &sites[ns*3]);
        for (int i=0;i<m;i++) bflag[ns+i]=1;
        ns += (size_t)m;
    }

    size_t n_bnd = ns;   /* boundary-constrained seeds precede the interior ones */

    /* interior seeds (area-weighted barycentric), kept clear of the boundary strip
     * so the boundary seeds own the rim cleanly. Without this an interior cell
     * wedges to the edge, breaks two boundary seeds' adjacency, and the dual
     * boundary saw-tooths inward. */
    double clear = 0.85 * spacing;
    for (; ns < N; ns++) {
        double pos[3] = {0,0,0};
        for (int tries = 0; ; tries++) {
            if (st) st->n_seed_tries++;
            double rr = rnd01(&rng) * tot;
            size_t loo=0, hii=in_nf-1;
            while (loo < hii) { size_t mid=(loo+hii)/2; if (cum[mid] < rr) loo=mid+1; else hii=mid; }
            size_t f = loo;
            double u=rnd01(&rng), v=rnd01(&rng);
            if (u+v > 1.0) { u=1.0-u; v=1.0-v; }
            const double *A=&V[(size_t)in_faces[f*3+0]*3];
            const double *B=&V[(size_t)in_faces[f*3+1]*3];
            const double *C=&V[(size_t)in_faces[f*3+2]*3];
            pos[0]=A[0]+u*(B[0]-A[0])+v*(C[0]-A[0]);
            pos[1]=A[1]+u*(B[1]-A[1])+v*(C[1]-A[1]);
            pos[2]=A[2]+u*(B[2]-A[2])+v*(C[2]-A[2]);
            if (bd.nseg == 0 || tries >= 8) break;
            double q[3]; size_t nsc = seg_bvh_project(&bvh, pos, q, NULL);
            if (st) { st->n_bnd_calls++; st->n_seg_scans += nsc; }
            double dz=pos[0]-q[0],dy=pos[1]-q[1],dx=pos[2]-q[2];
            double cl = field ? 0.85*sizing_h(field, pos[0],pos[1],pos[2]) : clear;
            if (dz*dz+dy*dy+dx*dx >= cl*cl) break;   /* far enough from the rim */
        }
        sites[ns*3+0]=pos[0]; sites[ns*3+1]=pos[1]; sites[ns*3+2]=pos[2];
    }

    if (st) st->tk_seed += cycle_now() - tk0;

    double *area = (double *)ARENA_ALLOC(arena, (size_t)(N*sizeof(double)));
    double *cent = (double *)ARENA_ALLOC(arena, (size_t)(N*3*sizeof(double)));
    double *NA   = (double *)ARENA_ALLOC(arena, (size_t)(N*6*sizeof(double)));
    double *bNA  = (double *)ARENA_ALLOC(arena, (size_t)(N*3*sizeof(double)));
    double lna = (double)o.lambda_na;

    /* The site update reads only site i's own cell moments and writes only
     * sites[i] (the trees, BVH and mesh are read-only), so the loop runs in
     * parallel with results identical to the serial order for any thread
     * count.  All scratch is allocated here, outside the parallel region. */
    int nthr = 1;
#ifdef _OPENMP
    nthr = omp_get_max_threads();
    if (nthr < 1) nthr = 1;
    if (nthr > 64) nthr = 64;
#endif
    LloydThread *lt = (LloydThread *)ARENA_CALLOC(arena, (size_t)nthr, sizeof(LloydThread));
    double *mv = o.verbose ? (double *)ARENA_ALLOC(arena, (size_t)(N*sizeof(double))) : NULL;
    int n_sites_i = (int)N;
    int par_ok = (N >= (size_t)CVT_PARALLEL_MIN_SITES && N < (size_t)INT32_MAX) ? 1 : 0;

    for (int it = 0; it < o.n_iters; it++) {
        double lcvt = (double)o.lambda_cvt0 * pow((double)o.cvt_decay, (double)it);
        if (lcvt < (double)o.lambda_cvt_min) lcvt = (double)o.lambda_cvt_min;
        if (st) tk0 = cycle_now();
        if (Rvd_accumulate(r, sites, N, area, cent, NA, bNA, NULL, NULL, field) != 0) return -4;
        if (st) { st->tk_rvd += cycle_now() - tk0; st->n_iters++; }
        int ii;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 64) num_threads(nthr) if (par_ok)
#endif
        for (ii = 0; ii < n_sites_i; ii++) {
            size_t i = (size_t)ii;
            int tid = 0;
#ifdef _OPENMP
            tid = omp_get_thread_num();
#endif
            LloydThread *th = &lt[tid];
            double np[3];
            if (bflag[i] && bd.nseg > 0) {
                /* boundary site: CVT-relax along its loop */
                uint64_t tb = st ? cycle_now() : 0;
                size_t nsc = seg_bvh_project(&bvh, &cent[i*3], np, NULL);
                if (st) {
                    th->tk_bnd += cycle_now() - tb;
                    th->n_bnd_calls++; th->n_seg_scans += nsc;
                }
            } else if (area[i] <= 0.0) {
                np[0]=sites[i*3+0]; np[1]=sites[i*3+1]; np[2]=sites[i*3+2];  /* empty cell */
            } else {
                /* interior CWF solve: (lna*N + lcvt*m*I) x = lna*b + lcvt*m*c */
                double m = area[i];
                const double *Ni=&NA[i*6], *bi=&bNA[i*3], *ci=&cent[i*3];
                double A6[6] = { lna*Ni[0]+lcvt*m, lna*Ni[1],        lna*Ni[2],
                                 lna*Ni[3]+lcvt*m, lna*Ni[4],        lna*Ni[5]+lcvt*m };
                double rhs[3] = { lna*bi[0]+lcvt*m*ci[0],
                                  lna*bi[1]+lcvt*m*ci[1],
                                  lna*bi[2]+lcvt*m*ci[2] };
                double x[3]; sym3_solve(A6, rhs, x);
                uint64_t tp = st ? cycle_now() : 0;
                size_t ncp = project_to_mesh_kd_lb(proj_kd, V, in_faces, max_tri_r, tri_cr, x, np,
                                                   th->cand, 1024, th->lbs);
                if (st) {
                    th->tk_proj += cycle_now() - tp;
                    th->n_proj_calls++; th->n_proj_cand += ncp;
                    if (ncp >= 1024) th->n_proj_cap++;
                }
                /* keep interior sites off the boundary strip every iteration (not
                 * just at seeding), so boundary seeds keep sole ownership of the
                 * rim and the dual boundary stays the clean boundary-seed loop. */
                if (bd.nseg > 0) {
                    uint64_t tb = st ? cycle_now() : 0;
                    double q[3]; size_t nsc = seg_bvh_project(&bvh, np, q, NULL);
                    if (st) {
                        th->tk_bnd += cycle_now() - tb;
                        th->n_bnd_calls++; th->n_seg_scans += nsc;
                    }
                    double dz=np[0]-q[0], dy=np[1]-q[1], dx=np[2]-q[2];
                    double dd=sqrt(dz*dz+dy*dy+dx*dx);
                    double cl = field ? 0.85*sizing_h(field, np[0],np[1],np[2]) : clear;
                    if (dd < cl && dd > 1e-9) {
                        double sc=cl/dd;
                        double push[3]={ q[0]+dz*sc, q[1]+dy*sc, q[2]+dx*sc };
                        tp = st ? cycle_now() : 0;
                        ncp = project_to_mesh_kd_lb(proj_kd, V, in_faces, max_tri_r, tri_cr, push, np,
                                                    th->cand, 1024, th->lbs);
                        if (st) {
                            th->tk_proj += cycle_now() - tp;
                            th->n_proj_calls++; th->n_proj_cand += ncp;
                            if (ncp >= 1024) th->n_proj_cap++;
                        }
                    }
                }
            }
            if (mv) {
                double dz=np[0]-sites[i*3+0], dy=np[1]-sites[i*3+1], dx=np[2]-sites[i*3+2];
                mv[i] = sqrt(dz*dz+dy*dy+dx*dx);
            }
            sites[i*3+0]=np[0]; sites[i*3+1]=np[1]; sites[i*3+2]=np[2];
        }
        if (st) {
            for (int t = 0; t < nthr; t++) {
                LloydThread *th = &lt[t];
                st->tk_bnd += th->tk_bnd;             st->tk_proj += th->tk_proj;
                st->n_bnd_calls += th->n_bnd_calls;   st->n_seg_scans += th->n_seg_scans;
                st->n_proj_calls += th->n_proj_calls; st->n_proj_cand += th->n_proj_cand;
                st->n_proj_cap += th->n_proj_cap;
                th->tk_bnd = th->tk_proj = 0;
                th->n_bnd_calls = th->n_seg_scans = 0;
                th->n_proj_calls = th->n_proj_cand = th->n_proj_cap = 0;
            }
        }
        if (o.verbose) {
            double moved = 0.0;
            for (size_t i = 0; i < N; i++) moved += mv[i];   /* serial order, as before */
            fprintf(stderr, "  cvt iter %2d: move %.5g  lcvt %.4f  (%zu bnd / %zu)\n",
                    it, moved, lcvt, n_bnd, N);
        }
        /* Mid-Lloyd probe: the dual these sites would emit now, built exactly
         * as the final dual below.  It reads the sites only (area/cent are
         * rebuilt by the next pass), so a probe that does not abort leaves the
         * attempt unchanged. */
        if (o.probe && it + 1 < o.n_iters && it < 64 &&
            ((o.probe_mask >> it) & 1u)) {
            uint64_t tpb = st ? cycle_now() : 0;
            int32_t *ptris = NULL; size_t pnt = 0;
            if (st) Rvd_set_stats(r, NULL);
            int prc = Rvd_accumulate(r, sites, N, area, cent, NULL, NULL, &ptris, &pnt, NULL);
            if (st) Rvd_set_stats(r, &st->rvd);
            if (prc != 0) return -4;
            float *pv = (float *)ARENA_ALLOC(arena, (size_t)(N*3*sizeof(float)));
            for (size_t i = 0; i < N*3; i++) pv[i] = (float)sites[i];
            int abort_now = o.probe(o.probe_ctx, it + 1, pv, N, ptris, pnt);
            if (st) { st->tk_probe += cycle_now() - tpb; st->n_probe++; }
            if (abort_now) return CVT_RC_PROBE_ABORT;
        }
    }

    int32_t *tris = NULL; size_t nt = 0;
    /* Dual is the geometric Restricted-Delaunay of the converged sites -- density
     * doesn't affect it, so pass NULL (the area/centroid outputs here are scratch). */
    if (st) tk0 = cycle_now();
    if (Rvd_accumulate(r, sites, N, area, cent, NULL, NULL, &tris, &nt, NULL) != 0) return -5;
    if (st) { st->tk_dual += cycle_now() - tk0; Rvd_set_stats(r, NULL); }

    float *ov = (float *)ARENA_ALLOC(arena, (size_t)(N*3*sizeof(float)));
    for (size_t i = 0; i < N*3; i++) ov[i] = (float)sites[i];

    *out_verts = ov; *out_nv = N;
    *out_faces = tris; *out_nf = nt;
    return 0;
}

/* Pinned-boundary variant: see cvt_remesh.h. Mirrors CVT_remesh with three
 * deltas -- boundary sites are the EXACT input boundary verts (no arc-length
 * resampling, no along-loop relax, bit-exact positions), the interior count is
 * derived from target_h, and the interior clearance uses the LOCAL boundary
 * segment length so coarse junction rings and fine hole rims each get a
 * proportionate berth. */
int CVT_remesh_pinned(Arena_T arena,
                      const float *in_verts, size_t in_nv,
                      const int32_t *in_faces, size_t in_nf,
                      double target_h,
                      const CvtOpts *opts,
                      float **out_verts, size_t *out_nv,
                      int32_t **out_faces, size_t *out_nf,
                      int32_t **out_site_src, size_t *out_n_pin) {
    if (!arena || !in_verts || !in_faces || in_nv < 3 || in_nf < 1 || target_h <= 0.0)
        return -1;
    if (out_verts) *out_verts = NULL;
    if (out_faces) *out_faces = NULL;
    if (out_nv) *out_nv = 0;
    if (out_nf) *out_nf = 0;
    if (out_site_src) *out_site_src = NULL;
    if (out_n_pin) *out_n_pin = 0;

    CvtOpts o; if (opts) o = *opts; else CVT_default_opts(&o);
    uint32_t rng = o.seed ? o.seed : 12345u;

    double *V = (double *)ARENA_ALLOC(arena, (size_t)(in_nv*3*sizeof(double)));
    for (size_t i = 0; i < in_nv*3; i++) V[i] = (double)in_verts[i];

    Rvd_T r = Rvd_new(arena, V, in_nv, in_faces, in_nf);
    if (!r) return -2;
    double surf_area = Rvd_mesh_area(r);
    if (surf_area <= 0.0) return -3;

    /* Boundary loops: every loop vert is a pinned site, in loop order. */
    Boundary bd; extract_boundary(arena, V, in_faces, in_nf, in_nv, &bd);
    if (bd.nloops == 0 || bd.nseg == 0) return -6;   /* closed patch: not our case */
    SegBvh bvh; seg_bvh_build(arena, bd.seg, bd.nseg, &bvh);
    size_t n_pin = (size_t)bd.loopoff[bd.nloops];
    if (n_pin < 3) return -6;

    /* Per-segment length, for the scale-aware clearance. */
    double *seglen = (double *)ARENA_ALLOC(arena, (size_t)(bd.nseg*sizeof(double)));
    for (size_t s = 0; s < bd.nseg; s++) {
        double dz=bd.seg[s*6+3]-bd.seg[s*6+0];
        double dy=bd.seg[s*6+4]-bd.seg[s*6+1];
        double dx=bd.seg[s*6+5]-bd.seg[s*6+2];
        seglen[s] = sqrt(dz*dz+dy*dy+dx*dx);
    }

    /* Interior count from the target spacing (uniform N = 2A/h^2). */
    size_t n_int = (size_t)(2.0*surf_area/(target_h*target_h) + 0.5);
    size_t N = n_pin + n_int;

    /* Projection acceleration (as CVT_remesh). */
    float *tc = (float *)ARENA_ALLOC(arena, (size_t)(in_nf*3*sizeof(float)));
    double max_tri_r = 0.0;
    for (size_t f = 0; f < in_nf; f++) {
        const double *A=&V[(size_t)in_faces[f*3+0]*3];
        const double *B=&V[(size_t)in_faces[f*3+1]*3];
        const double *C=&V[(size_t)in_faces[f*3+2]*3];
        double c0=(A[0]+B[0]+C[0])/3.0, c1=(A[1]+B[1]+C[1])/3.0, c2=(A[2]+B[2]+C[2])/3.0;
        tc[f*3+0]=(float)c0; tc[f*3+1]=(float)c1; tc[f*3+2]=(float)c2;
        double ra=sqrt((A[0]-c0)*(A[0]-c0)+(A[1]-c1)*(A[1]-c1)+(A[2]-c2)*(A[2]-c2));
        double rb=sqrt((B[0]-c0)*(B[0]-c0)+(B[1]-c1)*(B[1]-c1)+(B[2]-c2)*(B[2]-c2));
        double rc2=sqrt((C[0]-c0)*(C[0]-c0)+(C[1]-c1)*(C[1]-c1)+(C[2]-c2)*(C[2]-c2));
        double rr=ra>rb?(ra>rc2?ra:rc2):(rb>rc2?rb:rc2);
        if (rr>max_tri_r) max_tri_r=rr;
    }
    KDTree_T proj_kd = KDTree_new(arena, tc, in_nf);
    int32_t *proj_cand = (int32_t *)ARENA_ALLOC(arena, (size_t)(1024*sizeof(int32_t)));

    /* Area CDF for interior seeding (uniform). */
    double *cum = (double *)ARENA_ALLOC(arena, (size_t)(in_nf*sizeof(double)));
    double tot = 0.0;
    for (size_t f = 0; f < in_nf; f++) {
        const double *A=&V[(size_t)in_faces[f*3+0]*3];
        const double *B=&V[(size_t)in_faces[f*3+1]*3];
        const double *C=&V[(size_t)in_faces[f*3+2]*3];
        double ab[3]={B[0]-A[0],B[1]-A[1],B[2]-A[2]};
        double ac[3]={C[0]-A[0],C[1]-A[1],C[2]-A[2]};
        double cr[3]={ab[1]*ac[2]-ab[2]*ac[1], ab[2]*ac[0]-ab[0]*ac[2], ab[0]*ac[1]-ab[1]*ac[0]};
        tot += 0.5*sqrt(cr[0]*cr[0]+cr[1]*cr[1]+cr[2]*cr[2]);
        cum[f] = tot;
    }

    double *sites = (double *)ARENA_ALLOC(arena, (size_t)(N*3*sizeof(double)));
    int32_t *src  = (int32_t *)ARENA_ALLOC(arena, (size_t)(N*sizeof(int32_t)));

    /* Pinned sites: exact boundary verts, loop order, bit-exact positions. */
    for (size_t i = 0; i < n_pin; i++) {
        int32_t vid = bd.loopv[i];
        sites[i*3+0]=V[(size_t)vid*3+0];
        sites[i*3+1]=V[(size_t)vid*3+1];
        sites[i*3+2]=V[(size_t)vid*3+2];
        src[i] = vid;
    }

    /* Interior seeds, scale-aware clearance from the boundary polyline. */
    for (size_t ns = n_pin; ns < N; ns++) {
        double pos[3] = {0,0,0};
        for (int tries = 0; ; tries++) {
            double rr = rnd01(&rng) * tot;
            size_t loo=0, hii=in_nf-1;
            while (loo < hii) { size_t mid=(loo+hii)/2; if (cum[mid] < rr) loo=mid+1; else hii=mid; }
            size_t f = loo;
            double u=rnd01(&rng), v=rnd01(&rng);
            if (u+v > 1.0) { u=1.0-u; v=1.0-v; }
            const double *A=&V[(size_t)in_faces[f*3+0]*3];
            const double *B=&V[(size_t)in_faces[f*3+1]*3];
            const double *C=&V[(size_t)in_faces[f*3+2]*3];
            pos[0]=A[0]+u*(B[0]-A[0])+v*(C[0]-A[0]);
            pos[1]=A[1]+u*(B[1]-A[1])+v*(C[1]-A[1]);
            pos[2]=A[2]+u*(B[2]-A[2])+v*(C[2]-A[2]);
            if (tries >= 8) break;
            double q[3]; size_t sidx;
            (void)seg_bvh_project(&bvh, pos, q, &sidx);
            double dz=pos[0]-q[0],dy=pos[1]-q[1],dx=pos[2]-q[2];
            /* Wedge rule: an interior site closer to a boundary segment than
             * HALF its length can capture the segment's midpoint and sever
             * the two pinned endpoints' dual adjacency (= a dropped boundary
             * edge = a crack). 0.85 * seglen clears the half-length threshold
             * for ANY segment; deliberately NOT capped by target_h -- long
             * segments (despur stair-step diagonals) need proportionally wide
             * berths, which is geometrically required, not a density knob. */
            double cl = 0.85 * seglen[sidx];
            if (dz*dz+dy*dy+dx*dx >= cl*cl) break;
        }
        sites[ns*3+0]=pos[0]; sites[ns*3+1]=pos[1]; sites[ns*3+2]=pos[2];
        src[ns] = -1;
    }

    double *area = (double *)ARENA_ALLOC(arena, (size_t)(N*sizeof(double)));
    double *cent = (double *)ARENA_ALLOC(arena, (size_t)(N*3*sizeof(double)));
    double *NA   = (double *)ARENA_ALLOC(arena, (size_t)(N*6*sizeof(double)));
    double *bNA  = (double *)ARENA_ALLOC(arena, (size_t)(N*3*sizeof(double)));
    double lna = (double)o.lambda_na;

    for (int it = 0; it < o.n_iters; it++) {
        double lcvt = (double)o.lambda_cvt0 * pow((double)o.cvt_decay, (double)it);
        if (lcvt < (double)o.lambda_cvt_min) lcvt = (double)o.lambda_cvt_min;
        if (Rvd_accumulate(r, sites, N, area, cent, NA, bNA, NULL, NULL, NULL) != 0) return -4;
        for (size_t i = n_pin; i < N; i++) {   /* pinned sites NEVER move */
            double np[3];
            if (area[i] <= 0.0) {
                np[0]=sites[i*3+0]; np[1]=sites[i*3+1]; np[2]=sites[i*3+2];
            } else {
                double m = area[i];
                const double *Ni=&NA[i*6], *bi=&bNA[i*3], *ci=&cent[i*3];
                double A6[6] = { lna*Ni[0]+lcvt*m, lna*Ni[1],        lna*Ni[2],
                                 lna*Ni[3]+lcvt*m, lna*Ni[4],        lna*Ni[5]+lcvt*m };
                double rhs[3] = { lna*bi[0]+lcvt*m*ci[0],
                                  lna*bi[1]+lcvt*m*ci[1],
                                  lna*bi[2]+lcvt*m*ci[2] };
                double x[3]; sym3_solve(A6, rhs, x);
                project_to_mesh_kd(proj_kd, V, in_faces, max_tri_r, x, np, proj_cand, 1024);
                /* scale-aware push off the boundary strip, every iteration
                 * (same uncapped wedge rule as at seeding) */
                {
                    double q[3]; size_t sidx;
                    (void)seg_bvh_project(&bvh, np, q, &sidx);
                    double dz=np[0]-q[0], dy=np[1]-q[1], dx=np[2]-q[2];
                    double dd=sqrt(dz*dz+dy*dy+dx*dx);
                    double cl = 0.85 * seglen[sidx];
                    if (dd < cl && dd > 1e-9) {
                        double sc=cl/dd;
                        double push[3]={ q[0]+dz*sc, q[1]+dy*sc, q[2]+dx*sc };
                        project_to_mesh_kd(proj_kd, V, in_faces, max_tri_r, push, np, proj_cand, 1024);
                    }
                }
            }
            sites[i*3+0]=np[0]; sites[i*3+1]=np[1]; sites[i*3+2]=np[2];
        }
    }

    int32_t *tris = NULL; size_t nt = 0;
    if (Rvd_accumulate(r, sites, N, area, cent, NULL, NULL, &tris, &nt, NULL) != 0) return -5;
    if (nt == 0) return -7;

    float *ov = (float *)ARENA_ALLOC(arena, (size_t)(N*3*sizeof(float)));
    for (size_t i = 0; i < N*3; i++) ov[i] = (float)sites[i];

    *out_verts = ov; *out_nv = N;
    *out_faces = tris; *out_nf = nt;
    if (out_site_src) *out_site_src = src;
    if (out_n_pin) *out_n_pin = n_pin;
    return 0;
}

/* ---- projection self-test: project_to_mesh_kd MUST equal the O(nf) brute force on
 * every query, including on meshes with wildly-varying triangle sizes (a giant
 * triangle whose closest point is far from any centroid) -- that is the case the
 * "nearest-centroid + 2*max_tri_r" ball radius must still capture. ---- */
static uint32_t g_ps = 0xB0Bu;
static double psr(double lo, double hi) { uint32_t x=g_ps; x^=x<<13; x^=x>>17; x^=x<<5; g_ps=x;
                                          return lo + (hi-lo)*((double)(x>>8)*(1.0/16777216.0)); }

/* ---- boundary BVH selftest: seg_bvh_project must return the identical segment
 * index and the bitwise-identical point as the reference scans, on ragged
 * loops at large offsets, duplicate and zero-length segments, collinear runs,
 * exact equidistant ties, far-away and non-finite queries. ---- */
static uint32_t g_bs = 0xC0FFEEu;
static double bsr(double lo, double hi) { uint32_t x=g_bs; x^=x<<13; x^=x>>17; x^=x<<5; g_bs=x;
                                          return lo + (hi-lo)*((double)(x>>8)*(1.0/16777216.0)); }

int CVT_boundary_bvh_selftest(void) {
    Arena_T a = Arena_new();
    int ok = 1;
    size_t checks = 0;
    g_bs = 0xC0FFEEu;
    for (int kind = 0; kind < 4 && ok; kind++) {
        size_t cap = 6000;
        double *seg = (double *)ARENA_ALLOC(a, (size_t)(cap*6*sizeof(double)));
        size_t ns = 0;
        double cz = 1.0e4, cy = 2.5e3, cx = 3.1e3;
        if (kind == 0 || kind == 3) {        /* jagged loops, shared vertices */
            int npt = kind == 0 ? 2500 : 300;
            double pz0 = 0, py0 = 0, px0 = 0, pz = 0, py = 0, px = 0;
            for (int i = 0; i <= npt; i++) {
                double th = 6.283185307179586 * (double)(i % npt) / (double)npt;
                double r = 50.0 + bsr(-2.0, 2.0);
                double qz = cz + bsr(-0.5, 0.5), qy = cy + r*cos(th), qx = cx + r*sin(th);
                if (i == 0) { pz0 = qz; py0 = qy; px0 = qx; }
                if (i == npt) { qz = pz0; qy = py0; qx = px0; }
                if (i > 0) {
                    double *g = &seg[ns*6]; ns++;
                    g[0]=pz; g[1]=py; g[2]=px; g[3]=qz; g[4]=qy; g[5]=qx;
                    if (i % 97 == 0 && ns < cap) {          /* exact duplicate */
                        memcpy(&seg[ns*6], g, 6*sizeof(double)); ns++;
                    }
                    if (i % 101 == 0 && ns < cap) {         /* zero-length */
                        double *z = &seg[ns*6]; ns++;
                        z[0]=qz; z[1]=qy; z[2]=qx; z[3]=qz; z[4]=qy; z[5]=qx;
                    }
                }
                pz = qz; py = qy; px = qx;
            }
        } else if (kind == 1) {              /* collinear run: endpoint ties */
            for (int i = 0; i < 400; i++) {
                double *g = &seg[ns*6]; ns++;
                g[0]=cz; g[1]=cy + 0.75*i; g[2]=cx; g[3]=cz; g[4]=cy + 0.75*(i+1); g[5]=cx;
            }
        } else {                             /* parallel pairs: equidistant ties */
            for (int i = 0; i < 300; i++) {
                double y = cy + 3.0*i;
                double *g = &seg[ns*6]; ns++;
                g[0]=cz-1.0; g[1]=y; g[2]=cx; g[3]=cz-1.0; g[4]=y+2.0; g[5]=cx;
                g = &seg[ns*6]; ns++;
                g[0]=cz+1.0; g[1]=y; g[2]=cx; g[3]=cz+1.0; g[4]=y+2.0; g[5]=cx;
            }
        }
        SegBvh b; seg_bvh_build(a, seg, ns, &b);
        if (!b.ok) { fprintf(stderr, "CVT_boundary_bvh_selftest FAIL: kind=%d not built\n", kind); ok = 0; break; }
        for (int q = 0; q < 4000 && ok; q++) {
            double p[3];
            int mode = q % 8;
            if (mode == 0) {                 /* exactly on a segment endpoint */
                size_t s = (size_t)(bsr(0.0, (double)ns)); if (s >= ns) s = ns - 1;
                memcpy(p, &seg[s*6 + (q & 1 ? 3 : 0)], 3*sizeof(double));
            } else if (mode == 1) {          /* segment midpoint */
                size_t s = (size_t)(bsr(0.0, (double)ns)); if (s >= ns) s = ns - 1;
                for (int k = 0; k < 3; k++) p[k] = 0.5*(seg[s*6+k] + seg[s*6+3+k]);
            } else if (mode == 2 && kind == 2) {   /* exactly between a parallel pair */
                p[0] = cz; p[1] = cy + 3.0*(double)(q % 300) + bsr(0.0, 2.0); p[2] = cx + bsr(-1.0, 1.0);
            } else if (mode == 3) {          /* far but inside the query limit */
                p[0] = cz + bsr(-3e3, 3e3); p[1] = cy + bsr(-3e3, 3e3); p[2] = cx + bsr(-3e3, 3e3);
            } else if (mode == 4 && q % 16 == 4) {   /* beyond the limit: scan fallback */
                p[0] = 1e9; p[1] = -1e9; p[2] = 1e9;
            } else if (mode == 5 && q % 40 == 5) {   /* non-finite query */
                p[0] = NAN; p[1] = cy; p[2] = cx;
            } else {                         /* near the boundary */
                p[0] = cz + bsr(-4.0, 4.0); p[1] = cy + bsr(-60.0, 360.0); p[2] = cx + bsr(-60.0, 60.0);
            }
            double ob[3], ot[3], obs[3];
            size_t sb = 99, st_ = 77;
            project_to_boundary_seg(seg, ns, p, ob, &sb);
            project_to_boundary(seg, ns, p, obs);
            (void)seg_bvh_project(&b, p, ot, &st_);
            checks++;
            if (sb != st_ || memcmp(ob, ot, sizeof(ob)) != 0 || memcmp(ob, obs, sizeof(ob)) != 0) {
                fprintf(stderr, "CVT_boundary_bvh_selftest FAIL: kind=%d q=%d seg %zu vs %zu\n",
                        kind, q, sb, st_);
                ok = 0;
            }
        }
    }
    Arena_dispose(&a);
    fprintf(stderr, "CVT_boundary_bvh_selftest: %s (%zu queries, bitwise vs scan)\n",
            ok ? "OK" : "FAIL", checks);
    return ok ? 0 : -1;
}

typedef struct { int calls; int last_iter; int abort_at; size_t last_nf; } ProbeSelftestCtx;
static int probe_selftest_cb(void *ctx, int iter, const float *verts, size_t nv,
                             const int32_t *faces, size_t nf)
{
    ProbeSelftestCtx *c = (ProbeSelftestCtx *)ctx;
    (void)verts; (void)nv; (void)faces;
    c->calls++;
    c->last_iter = iter;
    c->last_nf = nf;
    return (c->abort_at > 0 && iter == c->abort_at) ? 1 : 0;
}

/* A non-aborting probe after every Lloyd iteration leaves the attempt
 * bitwise unchanged; an aborting probe stops it at the requested iteration. */
int CVT_probe_selftest(void) {
    Arena_T a = Arena_new();
    int fails = 0;
    const int M = 24, side = M + 1;
    size_t nv = (size_t)side * side, nf = (size_t)M * M * 2;
    float *V = (float *)ARENA_ALLOC(a, nv * 3 * sizeof(float));
    int32_t *F = (int32_t *)ARENA_ALLOC(a, nf * 3 * sizeof(int32_t));
    for (int i = 0; i < side; i++) for (int j = 0; j < side; j++) {
        size_t idx = (size_t)i * side + j;
        double u = (double)i / M, v = (double)j / M;
        V[idx*3+0] = (float)(1.5 * sin(5.0 * u) * cos(4.0 * v));
        V[idx*3+1] = (float)(30.0 * u);
        V[idx*3+2] = (float)(30.0 * v);
    }
    size_t fi = 0;
    for (int i = 0; i < M; i++) for (int j = 0; j < M; j++) {
        int32_t v00 = (int32_t)((size_t)i*side + j), v01 = v00 + 1;
        int32_t v10 = (int32_t)((size_t)(i+1)*side + j), v11 = v10 + 1;
        F[fi*3+0] = v00; F[fi*3+1] = v10; F[fi*3+2] = v11; fi++;
        F[fi*3+0] = v00; F[fi*3+1] = v11; F[fi*3+2] = v01; fi++;
    }
    int prc = 0;
    CvtPrep_T prep = CVT_prepare(a, V, nv, F, nf, NULL, NULL, &prc);
    if (!prep) { Arena_dispose(&a); fprintf(stderr, "CVT_probe_selftest: prepare failed\n"); return 1; }
    for (int pass = 0; pass < 3 && !fails; pass++) {
        size_t N = (size_t)(60 + 70 * pass);
        CvtOpts o; CVT_default_opts(&o); o.n_iters = 12;
        float *v0 = NULL, *v1 = NULL; int32_t *f0 = NULL, *f1 = NULL;
        size_t nv0 = 0, nf0 = 0, nv1 = 0, nf1 = 0;
        int r0 = CVT_remesh_prepared(a, prep, N, &o, &v0, &nv0, &f0, &nf0);
        ProbeSelftestCtx pc = { 0, 0, 0, 0 };
        o.probe = probe_selftest_cb; o.probe_ctx = &pc;
        o.probe_mask = (1u << 11) - 1u;            /* after iterations 1..11 */
        int r1 = CVT_remesh_prepared(a, prep, N, &o, &v1, &nv1, &f1, &nf1);
        if (r0 != 0 || r1 != 0 || nv0 != nv1 || nf0 != nf1 || pc.calls != 11 ||
            pc.last_iter != 11 || pc.last_nf == 0 ||
            memcmp(v0, v1, nv0 * 3 * sizeof(float)) != 0 ||
            memcmp(f0, f1, nf0 * 3 * sizeof(int32_t)) != 0) {
            fprintf(stderr, "CVT_probe_selftest FAIL: pass %d (rc %d/%d, nv %zu/%zu, nf %zu/%zu, probes %d)\n",
                    pass, r0, r1, nv0, nv1, nf0, nf1, pc.calls);
            fails++;
        }
        ProbeSelftestCtx pa = { 0, 0, 4, 0 };
        o.probe_ctx = &pa;
        o.probe_mask = 1u << 3;                     /* after iteration 4 only */
        float *v2 = NULL; int32_t *f2 = NULL; size_t nv2 = 0, nf2 = 0;
        int r2 = CVT_remesh_prepared(a, prep, N, &o, &v2, &nv2, &f2, &nf2);
        if (r2 != CVT_RC_PROBE_ABORT || pa.calls != 1 || pa.last_iter != 4 || nf2 != 0) {
            fprintf(stderr, "CVT_probe_selftest FAIL: abort pass %d (rc %d, probes %d, iter %d)\n",
                    pass, r2, pa.calls, pa.last_iter);
            fails++;
        }
    }
    Arena_dispose(&a);
    fprintf(stderr, "CVT_probe_selftest: %s (non-aborting probe bitwise-neutral at 3 site counts; abort at iter 4)\n",
            fails ? "FAIL" : "OK");
    return fails ? 1 : 0;
}

int CVT_projection_selftest(void) {
    Arena_T a = Arena_new();
    int ok = 1;
    g_ps = 0xB0Bu;

    for (int kind = 0; kind < 3; kind++) {
        int M = 16, side = M+1;
        size_t base_nv = (size_t)side*side;
        size_t base_nf = (size_t)M*M*2;
        int giant = (kind == 2);
        size_t nv = base_nv + (giant ? 3 : 0);
        size_t nf = base_nf + (giant ? 1 : 0);
        double  *V = (double *)ARENA_ALLOC(a, (size_t)(nv*3*sizeof(double)));
        int32_t *F = (int32_t *)ARENA_ALLOC(a, (size_t)(nf*3*sizeof(int32_t)));

        for (int i=0;i<side;i++) for (int j=0;j<side;j++) {
            size_t idx=(size_t)i*side+j;
            double u=(double)i/M, v=(double)j/M;
            if (kind==1) { u=pow(u,1.7); v=pow(v,1.7); }          /* graded: mixed triangle sizes */
            V[idx*3+0]= (kind==0||kind==2)?0.0 : 0.4*sin(6.2831*u)*sin(6.2831*v); /* wavy for kind 1 */
            V[idx*3+1]=10.0*u; V[idx*3+2]=10.0*v;
        }
        size_t fi=0;
        for (int i=0;i<M;i++) for (int j=0;j<M;j++) {
            int32_t v00=(int32_t)((size_t)i*side+j), v01=(int32_t)((size_t)i*side+j+1);
            int32_t v10=(int32_t)((size_t)(i+1)*side+j), v11=(int32_t)((size_t)(i+1)*side+j+1);
            F[fi*3+0]=v00;F[fi*3+1]=v10;F[fi*3+2]=v11;fi++;
            F[fi*3+0]=v00;F[fi*3+1]=v11;F[fi*3+2]=v01;fi++;
        }
        if (giant) { /* one huge triangle far to the side */
            int32_t g0=(int32_t)base_nv, g1=g0+1, g2=g0+2;
            V[(size_t)g0*3+0]=0; V[(size_t)g0*3+1]=-40; V[(size_t)g0*3+2]=-40;
            V[(size_t)g1*3+0]=0; V[(size_t)g1*3+1]= 60; V[(size_t)g1*3+2]=-40;
            V[(size_t)g2*3+0]=0; V[(size_t)g2*3+1]=-40; V[(size_t)g2*3+2]= 60;
            F[base_nf*3+0]=g0; F[base_nf*3+1]=g1; F[base_nf*3+2]=g2;
        }

        /* build the projection KD-tree exactly as CVT_remesh does */
        float *tc = (float *)ARENA_ALLOC(a, (size_t)(nf*3*sizeof(float)));
        double max_tri_r = 0.0;
        for (size_t f=0; f<nf; f++) {
            const double *A=&V[(size_t)F[f*3+0]*3], *B=&V[(size_t)F[f*3+1]*3], *C=&V[(size_t)F[f*3+2]*3];
            double c0=(A[0]+B[0]+C[0])/3.0, c1=(A[1]+B[1]+C[1])/3.0, c2=(A[2]+B[2]+C[2])/3.0;
            tc[f*3+0]=(float)c0; tc[f*3+1]=(float)c1; tc[f*3+2]=(float)c2;
            double ra=sqrt((A[0]-c0)*(A[0]-c0)+(A[1]-c1)*(A[1]-c1)+(A[2]-c2)*(A[2]-c2));
            double rb=sqrt((B[0]-c0)*(B[0]-c0)+(B[1]-c1)*(B[1]-c1)+(B[2]-c2)*(B[2]-c2));
            double rc=sqrt((C[0]-c0)*(C[0]-c0)+(C[1]-c1)*(C[1]-c1)+(C[2]-c2)*(C[2]-c2));
            double rr=ra>rb?(ra>rc?ra:rc):(rb>rc?rb:rc);
            if (rr>max_tri_r) max_tri_r=rr;
        }
        KDTree_T kd = KDTree_new(a, tc, nf);
        int32_t *cand = (int32_t *)ARENA_ALLOC(a, (size_t)((nf+1)*sizeof(int32_t)));
        int32_t *cand2 = (int32_t *)ARENA_ALLOC(a, (size_t)((nf+1)*sizeof(int32_t)));
        double *lbs = (double *)ARENA_ALLOC(a, (size_t)((nf+1)*sizeof(double)));
        double *tri_cr = (double *)ARENA_ALLOC(a, (size_t)(nf*4*sizeof(double)));
        for (size_t f=0; f<nf; f++) {
            const double *A=&V[(size_t)F[f*3+0]*3], *B=&V[(size_t)F[f*3+1]*3], *C=&V[(size_t)F[f*3+2]*3];
            double c0=(A[0]+B[0]+C[0])/3.0, c1=(A[1]+B[1]+C[1])/3.0, c2=(A[2]+B[2]+C[2])/3.0;
            double ra=sqrt((A[0]-c0)*(A[0]-c0)+(A[1]-c1)*(A[1]-c1)+(A[2]-c2)*(A[2]-c2));
            double rb=sqrt((B[0]-c0)*(B[0]-c0)+(B[1]-c1)*(B[1]-c1)+(B[2]-c2)*(B[2]-c2));
            double rc=sqrt((C[0]-c0)*(C[0]-c0)+(C[1]-c1)*(C[1]-c1)+(C[2]-c2)*(C[2]-c2));
            tri_cr[f*4+0]=c0; tri_cr[f*4+1]=c1; tri_cr[f*4+2]=c2;
            tri_cr[f*4+3]=ra>rb?(ra>rc?ra:rc):(rb>rc?rb:rc);
        }
        /* lower-bound projection == the original scan, bitwise, at several caps
         * (a small cap truncates the candidate list as production's 1024 does) */
        static const size_t lb_caps[4] = { 7, 64, 1024, (size_t)-1 };
        for (int q=0; q<4000 && ok; q++) {
            double p[3] = { psr(-6.0,6.0), psr(-45.0,65.0), psr(-45.0,65.0) };
            if (q % 7 == 0) {               /* on a vertex / an edge midpoint: exact ties */
                size_t vi = (size_t)(psr(0.0, (double)base_nv - 1.0));
                p[0] = V[vi*3+0]; p[1] = V[vi*3+1]; p[2] = V[vi*3+2];
                if (q % 14 == 0 && vi + 1 < base_nv) {
                    p[1] = 0.5*(V[vi*3+1]+V[(vi+1)*3+1]); p[2] = 0.5*(V[vi*3+2]+V[(vi+1)*3+2]);
                }
            }
            size_t cap = lb_caps[q % 4] > nf + 1 ? nf + 1 : lb_caps[q % 4];
            double o1[3], o2[3];
            size_t n1 = project_to_mesh_kd(kd, V, F, max_tri_r, p, o1, cand, cap);
            size_t n2 = project_to_mesh_kd_lb(kd, V, F, max_tri_r, tri_cr, p, o2, cand2, cap, lbs);
            if (n1 != n2 || memcmp(o1, o2, sizeof(o1)) != 0) {
                fprintf(stderr, "CVT_projection_selftest FAIL: lower-bound projection differs (kind %d q %d cap %zu)\n",
                        kind, q, cap);
                ok = 0;
            }
        }

        for (int q=0; q<3000; q++) {
            double p[3] = { psr(-6.0,6.0), psr(-45.0,65.0), psr(-45.0,65.0) };
            double okd[3], obf[3];
            project_to_mesh_kd(kd, V, F, max_tri_r, p, okd, cand, nf+1);
            project_to_mesh(V, F, nf, p, obf);
            double dk=(p[0]-okd[0])*(p[0]-okd[0])+(p[1]-okd[1])*(p[1]-okd[1])+(p[2]-okd[2])*(p[2]-okd[2]);
            double db=(p[0]-obf[0])*(p[0]-obf[0])+(p[1]-obf[1])*(p[1]-obf[1])+(p[2]-obf[2])*(p[2]-obf[2]);
            if (fabs(sqrt(dk)-sqrt(db)) > 1e-3*(sqrt(db)+1.0)) {
                fprintf(stderr, "CVT_projection_selftest FAIL: kind=%d q=(%.2f,%.2f,%.2f) kd_d=%.5f brute_d=%.5f\n",
                        kind, p[0],p[1],p[2], sqrt(dk), sqrt(db));
                ok = 0; break;
            }
        }
        if (!ok) break;
    }

    Arena_dispose(&a);
    fprintf(stderr, "CVT_projection_selftest: %s\n", ok ? "OK" : "FAIL");
    return ok ? 0 : -1;
}
