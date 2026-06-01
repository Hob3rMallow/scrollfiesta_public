/*
 * ball_pivot.c — Ball-Pivoting Algorithm (Bernardini 1999) for the
 * vesuvius-c pipeline. Multi-seed: re-seeds when the front empties
 * so the whole cloud is covered. Manifold pre-check rejects pivots
 * that would create a 3-fan on an interior edge.
 *
 * Algorithm same as scripts/ball_pivot_standalone.c; the only
 * difference is the output buffer goes through the arena instead of
 * malloc. Internal scratch (grid, edge store, queue, used bitmap,
 * face accumulator) uses malloc/free since it's local to a single
 * BallPivot_reconstruct call.
 */
#define _USE_MATH_DEFINES
#include "ball_pivot.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Debug-only per-pivot trace (set from BPA_TRACE in BallPivot_bridge). */
static int g_bpa_trace = 0;
/* Diagnostic (BPA_NO_NORMAL_GATE): bypass the per-vertex MLS normal gate to test
 * whether flipped normals are what strands the no_pivot holes. Default 0. */
static int  g_no_normal_gate = 0;
/* Pivot knockout tally (diagnostic): how often the theta-best candidate was
 * rejected by the empty-ball test vs by the normal gate. */
static long g_ko_empty  = 0;
static long g_ko_normal = 0;
/* Dead-end classification (BPA_DEBUG): when a front edge yields NO pivot
 * (becomes a hole boundary), WHY? Partitions the no_pivot count:
 *   gap    — no neighbour within 2*rho admitted a ball center (true sampling
 *            gap; a bigger rho could reach, at the cost of inter-wrap safety)
 *   empty  — candidates existed but every one failed the empty-ball test (a
 *            third point intrudes into the ball — local self-proximity/noise)
 *   normal — candidates existed but every one failed the per-vertex normal gate
 *   mixed  — both empty-ball and normal rejected candidates
 *   theta  — candidates had ball centers but all were behind the front (theta
 *            degenerate); geometrically unreachable from this edge */
static long g_dead_gap = 0, g_dead_empty = 0, g_dead_normal = 0,
            g_dead_mixed = 0, g_dead_theta = 0;
/* Mixed-dead-end margin measurement (BPA_DEBUG). For each candidate rejected at
 * a MIXED dead end we record, first-order: how much smaller rho would have to be
 * to clear the empty-ball intruder (rho - dist(O, nearest other point)), and for
 * normal-gate rejections how far past 90deg the (O-c).n angle is. The cross-tab
 * "recoverable" = empty-rejected AND would still pass the normal gate -> the only
 * candidates a smaller-radius re-arm can actually convert. */
static int    g_dead_measure = 0;
enum { DEADCAP = 128 };
static int    g_n_erec = 0, g_n_nrec = 0;        /* per-pivot rejection records */
static double g_erec_O[DEADCAP][3];
static int    g_erec_idx[DEADCAP];               /* candidate vertex (+normal side) */
static float  g_nrec_past90sin[DEADCAP];         /* -(O-c).n/(|O-c||n|) in [0,1]  */
enum { NBD = 8, NBA = 7, NBH = 9 };
/* empty-ball depth (rho shrink to clear intruder), only for +normal-side (pivot-
 * relevant) candidates; split by whether the BLOCKING intruder's normal is anti-
 * aligned with the candidate's (i.e. it belongs to the OTHER surface layer). */
static long   g_hd_all[NBD], g_hd_anti[NBD];
static long   g_ha[NBA];                         /* normal-gate angle past 90deg    */
/* signed height of the blocker above the candidate's tangent plane (vox): ~0 =
 * coplanar Delaunay violation (radius-invariant); large = off-plane point. */
static long   g_hh[NBH];
static long   g_mix_erej = 0, g_mix_anti = 0, g_mix_nrej = 0;
static double g_depth_sum = 0, g_past90_sum = 0, g_absh_sum = 0;
/* Tunable slack on the per-vertex normal gate (BPA_NORMAL_TOL, default 0 =
 * strict). The gate rejects a pivot candidate whose (O - vertex) points to the
 * wrong side of the vertex normal. With tol > 0 the reject threshold becomes
 * dot <= -tol*|O-c|*|n|, i.e. candidates up to ~asin(tol) past 90 deg are
 * accepted. Tests whether the holes are a too-strict tolerance (closes with a
 * little slack, genus stays 0) or genuine normal inconsistency (slack wrecks
 * genus). */
static double g_normal_tol = 0.0;

typedef struct { float z, y, x; } Vec3;

static inline void v_sub(double *r, const double *a, const double *b)
{ r[0]=a[0]-b[0]; r[1]=a[1]-b[1]; r[2]=a[2]-b[2]; }
static inline void v_cross(double *r, const double *a, const double *b)
{ r[0]=a[1]*b[2]-a[2]*b[1]; r[1]=a[2]*b[0]-a[0]*b[2]; r[2]=a[0]*b[1]-a[1]*b[0]; }
static inline double v_dot(const double *a, const double *b)
{ return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
static inline double v_len(const double *a) { return sqrt(v_dot(a,a)); }
static inline void v_scale(double *r, const double *a, double s)
{ r[0]=a[0]*s; r[1]=a[1]*s; r[2]=a[2]*s; }
static inline void vec3_to_d(double *r, const Vec3 *v)
{ r[0]=v->z; r[1]=v->y; r[2]=v->x; }
static inline double dist2(const double *a, const double *b)
{ double dz=a[0]-b[0], dy=a[1]-b[1], dx=a[2]-b[2]; return dz*dz+dy*dy+dx*dx; }

/* ---- Geometric fold-back guard (used in bpa_grow's glue acceptance) --------
 * Where two triangles share an edge on a consistently-wound surface, the two
 * apexes sit on OPPOSITE sides of that edge. Computing BOTH normals from the
 * shared edge's endpoints (sa,sb) and each apex, dot(n_new,n_old) is then
 * strongly NEGATIVE. If the apexes fall on the SAME side instead -- the new
 * triangle folds flat back over the existing one (dihedral -> 0, interiors
 * overlap) -- the dot is strongly POSITIVE. That is a self-intersection the
 * topological winding guards below cannot see: the glue's winding is still
 * correct, only the geometry has folded. Rejecting it never opens a hole, since
 * the existing triangle already surfaces that area. Reject when the fold
 * dihedral is sharper than ~45 deg (dot > BRIDGE_FOLD_COS); a real recto sheet
 * never folds that tight within one wrap (that would be an inter-wrap touch,
 * forbidden), so legitimate curvature (dot < 0) is untouched. */
#define BRIDGE_FOLD_COS 0.7

/* Unit normal of triangle (a,b,c) by vertex index; 0 if degenerate. */
static int tri_unit_normal_i(const Vec3 *V, int a, int b, int c, double n[3])
{
    double pa[3], pb[3], pc[3], ab[3], ac[3];
    vec3_to_d(pa, &V[a]); vec3_to_d(pb, &V[b]); vec3_to_d(pc, &V[c]);
    v_sub(ab, pb, pa); v_sub(ac, pc, pa); v_cross(n, ab, ac);
    double l = v_len(n);
    if (l < 1e-12) return 0;
    n[0] /= l; n[1] /= l; n[2] /= l;
    return 1;
}

/* 1 iff the new triangle (shared edge sa-sb, apex new_apex) folds flat back over
 * the existing triangle on that edge (apex old_apex). Both normals are taken in
 * the SAME (sa,sb,apex) order, so same-side apexes give dot > +BRIDGE_FOLD_COS. */
static int glue_folds_back(const Vec3 *V, int sa, int sb, int new_apex, int old_apex)
{
    if (old_apex < 0 || old_apex == sa || old_apex == sb) return 0;
    double nn[3], no[3];
    if (!tri_unit_normal_i(V, sa, sb, new_apex, nn)) return 0;
    if (!tri_unit_normal_i(V, sa, sb, old_apex, no)) return 0;
    return v_dot(nn, no) > BRIDGE_FOLD_COS;
}

typedef struct {
    int gz, gy, gx;
    double zmin, ymin, xmin, cell_size;
    int *head, *next;
} Grid;

static Grid *grid_build(const Vec3 *V, int nv, double cell_size)
{
    Grid *g = (Grid *)calloc(1, sizeof(*g));
    g->cell_size = cell_size;
    if (nv <= 0) return g;
    double zmin=V[0].z, zmax=V[0].z, ymin=V[0].y, ymax=V[0].y, xmin=V[0].x, xmax=V[0].x;
    for (int i = 1; i < nv; i++) {
        if (V[i].z<zmin)zmin=V[i].z; if (V[i].z>zmax)zmax=V[i].z;
        if (V[i].y<ymin)ymin=V[i].y; if (V[i].y>ymax)ymax=V[i].y;
        if (V[i].x<xmin)xmin=V[i].x; if (V[i].x>xmax)xmax=V[i].x;
    }
    g->zmin=zmin; g->ymin=ymin; g->xmin=xmin;
    g->gz = (int)ceil((zmax-zmin)/cell_size)+1;
    g->gy = (int)ceil((ymax-ymin)/cell_size)+1;
    g->gx = (int)ceil((xmax-xmin)/cell_size)+1;
    size_t nc = (size_t)g->gz*g->gy*g->gx;
    g->head = (int *)malloc(nc * sizeof(int));
    g->next = (int *)malloc((size_t)nv * sizeof(int));
    for (size_t i = 0; i < nc; i++) g->head[i] = -1;
    for (int i = 0; i < nv; i++) {
        int cz = (int)floor((V[i].z-zmin)/cell_size);
        int cy = (int)floor((V[i].y-ymin)/cell_size);
        int cx = (int)floor((V[i].x-xmin)/cell_size);
        size_t c = (size_t)cz*g->gy*g->gx + (size_t)cy*g->gx + cx;
        g->next[i] = g->head[c]; g->head[c] = i;
    }
    return g;
}
static void grid_free(Grid *g) { if (!g) return; free(g->head); free(g->next); free(g); }

typedef int (*neighbor_cb)(int idx, void *user);
static void grid_query(const Grid *g, double cz, double cy, double cx,
                        double radius, neighbor_cb cb, void *user, const Vec3 *V)
{
    int reach = (int)ceil(radius/g->cell_size);
    int qz = (int)floor((cz-g->zmin)/g->cell_size);
    int qy = (int)floor((cy-g->ymin)/g->cell_size);
    int qx = (int)floor((cx-g->xmin)/g->cell_size);
    double r2 = radius*radius;
    for (int dz=-reach; dz<=reach; dz++) {
        int iz = qz+dz; if (iz<0||iz>=g->gz) continue;
        for (int dy=-reach; dy<=reach; dy++) {
            int iy = qy+dy; if (iy<0||iy>=g->gy) continue;
            for (int dx=-reach; dx<=reach; dx++) {
                int ix = qx+dx; if (ix<0||ix>=g->gx) continue;
                int cell = (iz*g->gy + iy)*g->gx + ix;
                for (int j=g->head[cell]; j!=-1; j=g->next[j]) {
                    double dvz=V[j].z-cz, dvy=V[j].y-cy, dvx=V[j].x-cx;
                    double d2=dvz*dvz+dvy*dvy+dvx*dvx;
                    if (d2 > r2) continue;
                    if (cb(j, user)) return;
                }
            }
        }
    }
}

static int sphere_center(const double *a, const double *b, const double *c,
                          double rho, int sign_side, double *O_out)
{
    double ab[3], ac[3], n[3];
    v_sub(ab,b,a); v_sub(ac,c,a); v_cross(n,ab,ac);
    double n_len2 = v_dot(n,n);
    if (n_len2 < 1e-30) return 0;
    double ab_dot = v_dot(ab,ab);
    double ac_dot = v_dot(ac,ac);
    double ab_ac = v_dot(ab,ac);
    double D = ab_dot*ac_dot - ab_ac*ab_ac;
    if (D < 1e-30) return 0;
    double alpha = ac_dot * (ab_dot - ab_ac) / (2.0 * D);
    double beta  = ab_dot * (ac_dot - ab_ac) / (2.0 * D);
    double c_o[3];
    c_o[0] = a[0] + alpha*ab[0] + beta*ac[0];
    c_o[1] = a[1] + alpha*ab[1] + beta*ac[1];
    c_o[2] = a[2] + alpha*ab[2] + beta*ac[2];
    double r_c2 = dist2(c_o, a);
    if (r_c2 > rho*rho) return 0;
    double h = sqrt(rho*rho - r_c2);
    double n_len = sqrt(n_len2);
    double nhat[3] = { n[0]/n_len, n[1]/n_len, n[2]/n_len };
    if (sign_side < 0) { nhat[0]=-nhat[0]; nhat[1]=-nhat[1]; nhat[2]=-nhat[2]; }
    O_out[0] = c_o[0] + nhat[0]*h;
    O_out[1] = c_o[1] + nhat[1]*h;
    O_out[2] = c_o[2] + nhat[2]*h;
    return 1;
}

typedef struct {
    const Vec3 *V;
    double O[3];
    double r2_shrunk;
    int ea, eb, ec;
    int hit;
} EmptyCtx;
static int empty_cb(int idx, void *user)
{
    EmptyCtx *c = (EmptyCtx *)user;
    if (idx==c->ea || idx==c->eb || idx==c->ec) return 0;
    double dz=c->V[idx].z-c->O[0], dy=c->V[idx].y-c->O[1], dx=c->V[idx].x-c->O[2];
    double d2 = dz*dz+dy*dy+dx*dx;
    if (d2 < c->r2_shrunk) { c->hit = 1; return 1; }
    return 0;
}
static int ball_empty(const Grid *g, const Vec3 *V, const double *O,
                       double rho, int ea, int eb, int ec)
{
    EmptyCtx c;
    c.V=V; c.O[0]=O[0]; c.O[1]=O[1]; c.O[2]=O[2];
    double r = rho - 1e-4;
    c.r2_shrunk = r*r;
    c.ea=ea; c.eb=eb; c.ec=ec; c.hit=0;
    grid_query(g, O[0], O[1], O[2], rho, empty_cb, &c, V);
    return !c.hit;
}

/* --- Mixed-dead-end margin diagnostics (BPA_DEBUG only) ------------------- */
typedef struct { const Vec3 *V; double O[3]; int ea, eb, ec; double dmin2; int imin; } NearCtx;
static int near_cb(int idx, void *user)
{
    NearCtx *c = (NearCtx *)user;
    if (idx==c->ea || idx==c->eb || idx==c->ec) return 0;
    double dz=c->V[idx].z-c->O[0], dy=c->V[idx].y-c->O[1], dx=c->V[idx].x-c->O[2];
    double d2 = dz*dz+dy*dy+dx*dx;
    if (d2 < c->dmin2) { c->dmin2 = d2; c->imin = idx; }
    return 0;
}
/* Distance from O to the nearest point other than ea/eb/ec within rho (the
 * intruder failing the empty-ball test); its index in *out_idx. Returns rho if
 * the ball is empty. */
static double closest_other_dist(const Grid *g, const Vec3 *V, const double *O,
                                 double rho, int ea, int eb, int ec, int *out_idx)
{
    NearCtx c; c.V=V; c.O[0]=O[0]; c.O[1]=O[1]; c.O[2]=O[2];
    c.ea=ea; c.eb=eb; c.ec=ec; c.dmin2 = rho*rho; c.imin = -1;
    grid_query(g, O[0], O[1], O[2], rho, near_cb, &c, V);
    if (out_idx) *out_idx = c.imin;
    return sqrt(c.dmin2);
}
static void deadbin_depth(double depth, int anti)
{
    static const double edg[NBD] = {0.02,0.05,0.1,0.2,0.4,0.8,1.6,1e30};
    int b = 0; while (b < NBD-1 && depth >= edg[b]) b++;
    g_hd_all[b]++; if (anti) g_hd_anti[b]++;
}
static void deadbin_angle(double deg)
{
    static const double edg[NBA] = {1,2,5,10,20,45,1e30};
    int b = 0; while (b < NBA-1 && deg >= edg[b]) b++;
    g_ha[b]++;
}
static void deadbin_height(double h)   /* signed height, vox */
{
    static const double edg[NBH] = {-0.4,-0.2,-0.05,0.05,0.2,0.4,0.8,1.6,1e30};
    int b = 0; while (b < NBH-1 && h >= edg[b]) b++;
    g_hh[b]++;
}
/* Fold this pivot's recorded rejections into the global mixed-dead-end stats. */
static void commit_mixed(const Grid *g, const Vec3 *V, const Vec3 *N,
                         double rho, int va, int vb)
{
    for (int k = 0; k < g_n_erec; k++) {
        int blk = -1;
        double d = closest_other_dist(g, V, g_erec_O[k], rho, va, vb, g_erec_idx[k], &blk);
        double depth = rho - d;                 /* how much smaller rho must be */
        int anti = 0;                           /* blocker on the other layer?  */
        if (blk >= 0) {
            double ni[3], nc[3], cc[3], bb[3];
            vec3_to_d(ni, &N[blk]); vec3_to_d(nc, &N[g_erec_idx[k]]);
            if (v_dot(ni, nc) < 0.0) anti = 1;
            /* signed height of blocker above candidate's tangent plane */
            vec3_to_d(cc, &V[g_erec_idx[k]]); vec3_to_d(bb, &V[blk]);
            double nl = v_len(nc);
            double dvec[3]; v_sub(dvec, bb, cc);
            double h = (nl > 0.0) ? v_dot(dvec, nc) / nl : 0.0;
            deadbin_height(h);
            g_absh_sum += (h < 0 ? -h : h);
        }
        deadbin_depth(depth, anti);
        g_mix_erej++; g_depth_sum += depth; if (anti) g_mix_anti++;
    }
    for (int k = 0; k < g_n_nrec; k++) {
        double s = g_nrec_past90sin[k]; if (s < 0) s = 0; if (s > 1) s = 1;
        double deg = asin(s) * 180.0 / M_PI;
        deadbin_angle(deg);
        g_mix_nrej++; g_past90_sum += deg;
    }
}

#define ES_FRONT     0
#define ES_INTERIOR  1
#define ES_BOUNDARY  2

typedef struct EdgeRec {
    int va, vb, v_opp;       /* va<=vb : UNDIRECTED dedup key */
    int dir_tail;            /* winding tail as the BPA created the edge; the
                              * half-edge runs dir_tail -> (the other of va,vb).
                              * Lets glue require an OPPOSITE-direction coincidence
                              * (Bernardini 1999 Fig 7); a SAME-direction one is a
                              * non-orientable connection (foldover) and is rejected. */
    double O_prev[3];
    int state;
    int next_in_bucket;
} EdgeRec;

typedef struct {
    EdgeRec *e;
    int n, cap;
    int *buckets;
    int n_buckets;
} EdgeStore;

static uint64_t edge_hash(int a, int b)
{
    uint64_t h = ((uint64_t)(uint32_t)a) * 0x9E3779B97F4A7C15ULL + (uint64_t)(uint32_t)b;
    h ^= h>>33; h *= 0xff51afd7ed558ccdULL; h ^= h>>33;
    return h;
}

static EdgeStore *edges_new(int cap)
{
    EdgeStore *es = (EdgeStore *)calloc(1, sizeof(*es));
    if (cap < 1024) cap = 1024;
    es->cap = cap;
    es->e = (EdgeRec *)malloc((size_t)cap * sizeof(EdgeRec));
    int nb = 1; while (nb < cap*2) nb <<= 1;
    es->n_buckets = nb;
    es->buckets = (int *)malloc((size_t)nb * sizeof(int));
    for (int i = 0; i < nb; i++) es->buckets[i] = -1;
    return es;
}
static void edges_free(EdgeStore *es) { if (!es) return; free(es->e); free(es->buckets); free(es); }

static int edges_find(const EdgeStore *es, int a, int b)
{
    if (a > b) { int t=a; a=b; b=t; }
    int bk = (int)(edge_hash(a, b) & (uint64_t)(es->n_buckets - 1));
    int idx = es->buckets[bk];
    while (idx >= 0) {
        if (es->e[idx].va == a && es->e[idx].vb == b) return idx;
        idx = es->e[idx].next_in_bucket;
    }
    return -1;
}

/* Directed coincidence test layered on the undirected store. Returns the slot
 * of {tail,head} (or -1). On hit, *is_reverse = 1 iff the stored half-edge runs
 * head->tail (OPPOSITE direction == a legal glue) and 0 iff it runs tail->head
 * (SAME direction == non-orientable; the pivot must be rejected). */
static int edges_find_dir(const EdgeStore *es, int tail, int head, int *is_reverse)
{
    int idx = edges_find(es, tail, head);
    if (idx < 0) { *is_reverse = 0; return -1; }
    *is_reverse = (es->e[idx].dir_tail == head) ? 1 : 0;
    return idx;
}

/* Insert an ES_FRONT half-edge. The caller passes (a,b) in WINDING order; the
 * undirected dedup key sorts them, but dir_tail records the winding tail (a) so
 * the glue step can tell an opposite-direction coincidence (legal) from a
 * same-direction one (non-orientable). vfront[] (may be NULL) counts incident
 * FRONT edges per vertex for the O(1) bowtie test. */
static int edges_insert(EdgeStore *es, int *vfront, int a, int b, int v_opp, const double *O)
{
    int tail = a;                       /* winding tail, BEFORE the sort */
    if (a > b) { int t=a; a=b; b=t; }
    if (es->n >= es->cap) {
        es->cap *= 2;
        es->e = (EdgeRec *)realloc(es->e, (size_t)es->cap * sizeof(EdgeRec));
    }
    int idx = es->n++;
    EdgeRec *r = &es->e[idx];
    r->va=a; r->vb=b; r->v_opp=v_opp; r->dir_tail=tail;
    r->O_prev[0]=O[0]; r->O_prev[1]=O[1]; r->O_prev[2]=O[2];
    r->state = ES_FRONT;
    int bk = (int)(edge_hash(a, b) & (uint64_t)(es->n_buckets - 1));
    r->next_in_bucket = es->buckets[bk]; es->buckets[bk] = idx;
    if (vfront) { vfront[a]++; vfront[b]++; }
    return idx;
}

/* Transition an edge's state, keeping the per-vertex FRONT-edge counter vfront[]
 * (may be NULL) in sync so the bowtie test is O(1). Only ES_FRONT memberships
 * are counted. */
static void edge_set_state(EdgeStore *es, int *vfront, int eid, int new_state)
{
    int old = es->e[eid].state;
    if (old == new_state) return;
    if (vfront) {
        int a = es->e[eid].va, bb = es->e[eid].vb;
        if (old == ES_FRONT && new_state != ES_FRONT) { vfront[a]--; vfront[bb]--; }
        else if (old != ES_FRONT && new_state == ES_FRONT) { vfront[a]++; vfront[bb]++; }
    }
    es->e[eid].state = new_state;
}

static void sort_by_dist(int *arr, int n, const Vec3 *V, int anchor)
{
    for (int i = 1; i < n; i++) {
        int v = arr[i];
        double dz=V[v].z-V[anchor].z, dy=V[v].y-V[anchor].y, dx=V[v].x-V[anchor].x;
        double dv = dz*dz+dy*dy+dx*dx;
        int j = i-1;
        while (j >= 0) {
            int u = arr[j];
            double duz=V[u].z-V[anchor].z, duy=V[u].y-V[anchor].y, dux=V[u].x-V[anchor].x;
            double du = duz*duz+duy*duy+dux*dux;
            if (du < dv || (du == dv && u < v)) break;
            arr[j+1] = arr[j]; j--;
        }
        arr[j+1] = v;
    }
}

typedef struct { int *cands; int n_cands; int cap; } CandList;

static int seed_collect(int idx, void *user)
{
    CandList *cl = (CandList *)user;
    if (cl->n_cands >= cl->cap) {
        cl->cap *= 2;
        cl->cands = (int *)realloc(cl->cands, (size_t)cl->cap * sizeof(int));
    }
    cl->cands[cl->n_cands++] = idx;
    return 0;
}

static int try_seed(const Grid *g, const Vec3 *V, const Vec3 *N,
                    int i, int j, int k, double rho,
                    int *out_face, double *out_O)
{
    if (i==j || j==k || i==k) return 0;
    double a[3], b[3], c[3];
    vec3_to_d(a,&V[i]); vec3_to_d(b,&V[j]); vec3_to_d(c,&V[k]);
    double na[3], nb[3], nc[3], n_avg[3];
    vec3_to_d(na,&N[i]); vec3_to_d(nb,&N[j]); vec3_to_d(nc,&N[k]);
    n_avg[0]=na[0]+nb[0]+nc[0]; n_avg[1]=na[1]+nb[1]+nc[1]; n_avg[2]=na[2]+nb[2]+nc[2];
    double ab[3], ac[3], tnorm[3];
    v_sub(ab,b,a); v_sub(ac,c,a); v_cross(tnorm,ab,ac);
    int sign = (v_dot(tnorm, n_avg) > 0) ? +1 : -1;
    double O[3];
    if (!sphere_center(a, b, c, rho, sign, O)) return 0;
    if (!g_no_normal_gate) {
        double diff[3], dl, nl;
        v_sub(diff, O, a); dl = v_len(diff); nl = v_len(na);
        if (v_dot(diff, na) <= -g_normal_tol * dl * nl) return 0;
        v_sub(diff, O, b); dl = v_len(diff); nl = v_len(nb);
        if (v_dot(diff, nb) <= -g_normal_tol * dl * nl) return 0;
        v_sub(diff, O, c); dl = v_len(diff); nl = v_len(nc);
        if (v_dot(diff, nc) <= -g_normal_tol * dl * nl) return 0;
    }
    if (!ball_empty(g, V, O, rho, i, j, k)) return 0;
    if (sign > 0) { out_face[0]=i; out_face[1]=j; out_face[2]=k; }
    else          { out_face[0]=i; out_face[1]=k; out_face[2]=j; }
    out_O[0]=O[0]; out_O[1]=O[1]; out_O[2]=O[2];
    return 1;
}

typedef struct {
    const Grid *g;
    const Vec3 *V, *N;
    double rho;
    double a[3], b[3], M[3], ab[3], ab_len;
    int va, vb;
    double O_prev_perp[3];
    double O_prev_perp_len;
    int best_v;
    double best_theta;
    double best_O[3];
    /* dead-end classification (per pivot attempt) */
    int n_sphere;     /* candidates that admitted >=1 ball center           */
    int n_rej_empty;  /* theta-eligible candidates killed by empty-ball test */
    int n_rej_normal; /* theta-eligible candidates killed by normal gate     */
    int n_theta_elig; /* candidates that reached the theta gate              */
} PivotCtx;

static int pivot_cb(int idx, void *user)
{
    PivotCtx *p = (PivotCtx *)user;
    if (idx == p->va || idx == p->vb) return 0;
    double c[3]; vec3_to_d(c, &p->V[idx]);
    double O_plus[3], O_minus[3];
    int has_p = sphere_center(p->a, p->b, c, p->rho, +1, O_plus);
    int has_m = sphere_center(p->a, p->b, c, p->rho, -1, O_minus);
    if (!has_p && !has_m) return 0;
    p->n_sphere++;
    double Os[2][3]; int n_O = 0;
    if (has_p) { Os[n_O][0]=O_plus[0]; Os[n_O][1]=O_plus[1]; Os[n_O][2]=O_plus[2]; n_O++; }
    if (has_m) { Os[n_O][0]=O_minus[0]; Os[n_O][1]=O_minus[1]; Os[n_O][2]=O_minus[2]; n_O++; }
    double u[3]; v_scale(u, p->ab, 1.0/p->ab_len);
    for (int q = 0; q < n_O; q++) {
        double *O_c = Os[q];
        double dM[3]; v_sub(dM, O_c, p->M);
        double t = v_dot(dM, p->ab) / (p->ab_len * p->ab_len);
        double tab[3]; v_scale(tab, p->ab, t);
        double perp[3]; v_sub(perp, dM, tab);
        double perp_len = v_len(perp);
        if (perp_len < 1e-12 || p->O_prev_perp_len < 1e-12) continue;
        double cr[3]; v_cross(cr, p->O_prev_perp, perp);
        double sin_t = v_dot(cr, u);
        double cos_t = v_dot(p->O_prev_perp, perp);
        double inv = 1.0 / (perp_len * p->O_prev_perp_len);
        sin_t *= inv; cos_t *= inv;
        double theta = atan2(sin_t, cos_t);
        if (theta < 0) theta += 2.0 * M_PI;
        if (theta < 1e-3) continue;
        if (theta < p->best_theta) {
            p->n_theta_elig++;
            if (!ball_empty(p->g, p->V, O_c, p->rho, p->va, p->vb, idx)) {
                g_ko_empty++; p->n_rej_empty++;
                if (g_dead_measure && g_n_erec < DEADCAP) {
                    double nv_[3]; vec3_to_d(nv_, &p->N[idx]);
                    double diff[3]; v_sub(diff, O_c, c);
                    if (v_dot(diff, nv_) > 0.0) {   /* +normal side: pivot-relevant */
                        int k = g_n_erec++;
                        g_erec_O[k][0]=O_c[0]; g_erec_O[k][1]=O_c[1]; g_erec_O[k][2]=O_c[2];
                        g_erec_idx[k]=idx;
                    }
                }
                continue;
            }
            double nv_[3]; vec3_to_d(nv_, &p->N[idx]);
            double diff[3]; v_sub(diff, O_c, c);
            if (!g_no_normal_gate) {
                double dl = v_len(diff), nl = v_len(nv_);
                if (v_dot(diff, nv_) <= -g_normal_tol * dl * nl) {
                    g_ko_normal++; p->n_rej_normal++;
                    if (g_dead_measure && g_n_nrec < DEADCAP) {
                        double denom = dl * nl;
                        g_nrec_past90sin[g_n_nrec++] =
                            (float)(denom > 0 ? -v_dot(diff, nv_) / denom : 0.0);
                    }
                    continue;
                }
            }
            p->best_v = idx; p->best_theta = theta;
            p->best_O[0]=O_c[0]; p->best_O[1]=O_c[1]; p->best_O[2]=O_c[2];
        }
    }
    return 0;
}

static int do_pivot(const Grid *g, const Vec3 *V, const Vec3 *N, double rho,
                     int va, int vb, const double *O_prev,
                     int *out_v, double *out_O)
{
    PivotCtx p;
    p.g=g; p.V=V; p.N=N; p.rho=rho;
    vec3_to_d(p.a, &V[va]); vec3_to_d(p.b, &V[vb]);
    v_sub(p.ab, p.b, p.a); p.ab_len = v_len(p.ab);
    if (p.ab_len < 1e-9) return 0;
    p.M[0]=(p.a[0]+p.b[0])*0.5; p.M[1]=(p.a[1]+p.b[1])*0.5; p.M[2]=(p.a[2]+p.b[2])*0.5;
    p.va=va; p.vb=vb;
    double dM[3]; v_sub(dM, O_prev, p.M);
    double t = v_dot(dM, p.ab) / (p.ab_len * p.ab_len);
    double tab[3]; v_scale(tab, p.ab, t);
    v_sub(p.O_prev_perp, dM, tab);
    p.O_prev_perp_len = v_len(p.O_prev_perp);
    p.best_v = -1; p.best_theta = 2.0 * M_PI + 1;
    p.n_sphere = p.n_rej_empty = p.n_rej_normal = p.n_theta_elig = 0;
    if (g_dead_measure) { g_n_erec = 0; g_n_nrec = 0; }
    grid_query(g, p.M[0], p.M[1], p.M[2], 2.0*rho, pivot_cb, &p, V);
    if (p.best_v < 0) {
        /* Dead end -> hole boundary. Classify why no candidate won. */
        if (p.n_theta_elig == 0)
            { if (p.n_sphere == 0) g_dead_gap++; else g_dead_theta++; }
        else if (p.n_rej_empty > 0 && p.n_rej_normal == 0)    g_dead_empty++;
        else if (p.n_rej_normal > 0 && p.n_rej_empty == 0)    g_dead_normal++;
        else if (p.n_rej_empty > 0 && p.n_rej_normal > 0) {
            g_dead_mixed++;
            if (g_dead_measure) commit_mixed(g, V, N, rho, va, vb);
        }
        else                                                  g_dead_theta++;
        return 0;
    }
    *out_v = p.best_v;
    out_O[0]=p.best_O[0]; out_O[1]=p.best_O[1]; out_O[2]=p.best_O[2];
    return 1;
}

/* Growth accumulators shared by the seeded and bridge entry points. */
typedef struct {
    int *F;       /* faces, 3 ints each */
    int nf;
    int face_cap;
    int *queue;   /* front-edge indices to process */
    int qh, qt;
    int queue_cap;
    /* Debug-only pivot tally (BPA_DEBUG). No effect on behaviour, except
     * dbg_fold, whose guard DOES reject pivots (see glue_folds_back). */
    int dbg_no_pivot, dbg_bowtie, dbg_orient, dbg_3fan, dbg_bowtie_ok, dbg_fold;
} BpaBuild;

static void bpa_build_init(BpaBuild *b, int n)
{
    b->face_cap = n * 2; if (b->face_cap < 16) b->face_cap = 16;
    b->F = (int *)malloc((size_t)b->face_cap * 3 * sizeof(int));
    b->nf = 0;
    b->queue_cap = b->face_cap * 3;
    b->queue = (int *)malloc((size_t)b->queue_cap * sizeof(int));
    b->qh = 0; b->qt = 0;
    b->dbg_no_pivot = b->dbg_bowtie = b->dbg_orient = b->dbg_3fan = 0;
    b->dbg_bowtie_ok = 0; b->dbg_fold = 0;
}
static void bpa_build_free(BpaBuild *b) { free(b->F); free(b->queue); }

static void bpa_add_face(BpaBuild *b, int a, int v, int c)
{
    if (b->nf >= b->face_cap) {
        b->face_cap *= 2;
        b->F = (int *)realloc(b->F, (size_t)b->face_cap * 3 * sizeof(int));
    }
    b->F[b->nf*3+0] = a; b->F[b->nf*3+1] = v; b->F[b->nf*3+2] = c; b->nf++;
}
static void bpa_enqueue(BpaBuild *b, int eid)
{
    if (b->qt >= b->queue_cap) {
        b->queue_cap *= 2;
        b->queue = (int *)realloc(b->queue, (size_t)b->queue_cap * sizeof(int));
    }
    b->queue[b->qt++] = eid;
}

/* Drain the front queue, advancing every ES_FRONT edge by one pivot until
 * none remain. Maintains TRUE directed half-edges (Bernardini 1999 §IV-D,
 * matching scripts/bpa_reference.cpp): the front edge ei is directed
 * dir_tail->head, the new triangle (t, v_new, h) is wound to traverse the
 * shared edge in reverse, and a coincident edge is glued ONLY when it runs the
 * opposite direction. A same-direction coincidence is a non-orientable (Mobius)
 * connection -- the foldover bug -- and is rejected (edge marked boundary).
 * Endpoints/O_prev/dir are cached into locals each iteration; reads of
 * es->e[...] after an edges_insert use indices (valid across realloc), never a
 * held EdgeRec pointer (the original inline loop dangled a stale `e`). */
static void bpa_grow(const Grid *g, const Vec3 *V, const Vec3 *N, double rho,
                     EdgeStore *es, uint8_t *used, int *vfront, BpaBuild *b,
                     int relax_bowtie)
{
    while (b->qh < b->qt) {
        int ei = b->queue[b->qh++];
        if (es->e[ei].state != ES_FRONT) continue;
        int va = es->e[ei].va, vb = es->e[ei].vb;       /* sorted: pivot geometry */
        int t  = es->e[ei].dir_tail;                    /* directed tail */
        int h  = (t == va) ? vb : va;                   /* directed head */
        double O_prev[3] = { es->e[ei].O_prev[0], es->e[ei].O_prev[1],
                             es->e[ei].O_prev[2] };
        int v_new; double O_new[3];
        if (!do_pivot(g, V, N, rho, va, vb, O_prev, &v_new, O_new)) {
            edge_set_state(es, vfront, ei, ES_BOUNDARY);
            b->dbg_no_pivot++;
            if (g_bpa_trace)
                fprintf(stderr, "  [trace] edge (%.1f,%.1f)-(%.1f,%.1f) NO_PIVOT (rho=%.2f)\n",
                        V[va].y, V[va].x, V[vb].y, V[vb].x, rho);
            continue;
        }
        /* Bowtie guard (Bernardini case 1 / reference onFront): a pivot landing
         * on a USED vertex carrying no incident FRONT edge would make a
         * non-manifold (bowtie) vertex. Reject. vfront[] makes this O(1). */
        if (used[v_new] && vfront[v_new] == 0) {
            if (!relax_bowtie) {
                edge_set_state(es, vfront, ei, ES_BOUNDARY);
                b->dbg_bowtie++;
                if (g_bpa_trace)
                    fprintf(stderr, "  [trace] edge (%.1f,%.1f)-(%.1f,%.1f) -> v(%.1f,%.1f) "
                            "BOWTIE (used,vfront=0)\n",
                            V[va].y, V[va].x, V[vb].y, V[vb].x, V[v_new].y, V[v_new].x);
                continue;
            }
            /* Bridge meeting (relax_bowtie): the opposite front fully consumed
             * this boundary vertex, so it has no FRONT edge -- but attaching to
             * it is exactly how the two seam fronts zip together. Allow it; the
             * transient vertex-touch is resolved as the zip completes (or by the
             * downstream pinch-split). Edge-manifoldness is still enforced by the
             * 3-fan guard below, and the pivot's per-vertex normal test already
             * rejected any opposite-facing wrap, so this cannot fuse wraps. */
            b->dbg_bowtie_ok++;
        }
        /* The new triangle (t, v_new, h) contributes directed half-edges
         * t->v_new and v_new->h. A coincident edge may be glued only if it runs
         * the REVERSE direction (Bernardini Fig 7 / findReverseEdgeOnFront). */
        int rev_tv = 0, rev_vh = 0;
        int e_tv = edges_find_dir(es, t, v_new, &rev_tv);
        int e_vh = edges_find_dir(es, v_new, h, &rev_vh);
        int abort_pivot = 0;
        /* 3-fan guard: a coincident edge already carrying 2 faces (ES_INTERIOR)
         * cannot take a third -- the one true non-manifold-edge case. */
        if (e_tv >= 0 && es->e[e_tv].state == ES_INTERIOR) { abort_pivot = 1; b->dbg_3fan++; }
        if (e_vh >= 0 && es->e[e_vh].state == ES_INTERIOR) { abort_pivot = 1; b->dbg_3fan++; }
        /* Orientation guard: a coincident still-open edge wound the SAME way is
         * a non-orientable connection (foldover). Reject. */
        if (e_tv >= 0 && es->e[e_tv].state != ES_INTERIOR && !rev_tv) { abort_pivot = 1; b->dbg_orient++; }
        if (e_vh >= 0 && es->e[e_vh].state != ES_INTERIOR && !rev_vh) { abort_pivot = 1; b->dbg_orient++; }
        /* Fold-back guard: the reverse-coincident edges above are legal glues by
         * winding, but reject the pivot if the new triangle folds flat back over
         * the existing triangle across that shared edge (see glue_folds_back).
         * This is the seam-bridge self-intersection the orient guard cannot
         * catch -- the front zipping onto itself at a ~0 dihedral. */
        if (!abort_pivot && e_tv >= 0 && rev_tv &&
            glue_folds_back(V, t, v_new, h, es->e[e_tv].v_opp)) {
            abort_pivot = 1; b->dbg_fold++;
        }
        if (!abort_pivot && e_vh >= 0 && rev_vh &&
            glue_folds_back(V, v_new, h, t, es->e[e_vh].v_opp)) {
            abort_pivot = 1; b->dbg_fold++;
        }
        /* Front-edge fold-back: the new triangle is ALSO glued, across the
         * advancing front edge (t,h), to that edge's originating triangle (whose
         * apex is es->e[ei].v_opp). This is the dominant self-fold -- the front
         * rolling straight back over its own parent -- and the e_tv/e_vh checks
         * above never see it (they test only the two NEW edges, not the front
         * edge). Reject if v_new lands on the parent apex's side of (t,h). */
        if (!abort_pivot &&
            glue_folds_back(V, t, h, v_new, es->e[ei].v_opp)) {
            abort_pivot = 1; b->dbg_fold++;
        }
        if (abort_pivot) {
            edge_set_state(es, vfront, ei, ES_BOUNDARY);
            continue;
        }
        bpa_add_face(b, t, v_new, h);
        used[v_new] = 1;
        edge_set_state(es, vfront, ei, ES_INTERIOR);
        /* Only reverse coincidences survive the guards above, so a hit always
         * closes cleanly to a 2-face manifold edge (FRONT or BOUNDARY -> closed). */
        if (e_tv >= 0) {
            edge_set_state(es, vfront, e_tv, ES_INTERIOR);
        } else {
            int new_id = edges_insert(es, vfront, t, v_new, h, O_new);
            bpa_enqueue(b, new_id);
        }
        if (e_vh >= 0) {
            edge_set_state(es, vfront, e_vh, ES_INTERIOR);
        } else {
            int new_id = edges_insert(es, vfront, v_new, h, t, O_new);
            bpa_enqueue(b, new_id);
        }
    }
}

/* Boundary-loop analysis (BPA_DEBUG): trace ES_BOUNDARY edges into loops,
 * histogram loop lengths, and for every genuine 3-edge loop (a single-triangle
 * pinhole) test whether its fill triangle's ball is actually EMPTY (and normal-
 * consistent). A 3-loop with an empty ball was simply MISSED by the front (an
 * ordering/bowtie issue, not Delaunay); a blocked one has a real intruder, whose
 * depth tells us whether it is borderline-cocircular (precision/Shewchuk would
 * matter) or robustly inside (it would not). */
static void analyze_holes(const EdgeStore *es, const Grid *g, const Vec3 *V,
                          const Vec3 *N, double rho, int n)
{
    int nb = 0;
    for (int i = 0; i < es->n; i++) if (es->e[i].state == ES_BOUNDARY) nb++;
    if (nb == 0) { fprintf(stderr, "[bpa_holes] no boundary edges\n"); return; }

    int *ba = (int *)malloc((size_t)nb*sizeof(int));
    int *bb = (int *)malloc((size_t)nb*sizeof(int));
    uint8_t *evis = (uint8_t *)calloc((size_t)nb, 1);
    { int k=0; for (int i=0;i<es->n;i++) if (es->e[i].state==ES_BOUNDARY)
        { ba[k]=es->e[i].va; bb[k]=es->e[i].vb; k++; } }

    int *deg = (int *)calloc((size_t)n, sizeof(int));
    for (int e=0;e<nb;e++){ deg[ba[e]]++; deg[bb[e]]++; }
    int *off = (int *)malloc((size_t)(n+1)*sizeof(int)); off[0]=0;
    for (int i=0;i<n;i++) off[i+1]=off[i]+deg[i];
    int *adj = (int *)malloc((size_t)off[n]*sizeof(int));
    int *cur = (int *)malloc((size_t)n*sizeof(int));
    memcpy(cur, off, (size_t)n*sizeof(int));
    for (int e=0;e<nb;e++){ adj[cur[ba[e]]++]=e; adj[cur[bb[e]]++]=e; }
    free(cur);

    long hist[6]={0,0,0,0,0,0};         /* 3,4,5-8,9-16,17-32,33+ */
    long n_loops=0, n_nonsimple=0;
    long t3=0, t3_fill=0, t3_fill_nok=0, t3_block=0, t3_block_border=0;
    double t3_block_depth_sum=0;

    for (int e0=0;e0<nb;e0++){
        if (evis[e0]) continue;
        int firstv=ba[e0], cur_v=bb[e0], cur_e=e0, len=1, simple=1, guard=0;
        int vv[8]; int nvv=0; vv[nvv++]=firstv;
        evis[e0]=1;
        while (guard++ < nb+4) {
            if (nvv<8) vv[nvv++]=cur_v;
            if (deg[cur_v]!=2) simple=0;
            int nexte=-1;
            for (int t=off[cur_v]; t<off[cur_v+1]; t++)
                { int ee=adj[t]; if (ee==cur_e || evis[ee]) continue; nexte=ee; break; }
            if (nexte<0) break;
            evis[nexte]=1; len++;
            int other = (ba[nexte]==cur_v) ? bb[nexte] : ba[nexte];
            cur_e=nexte; cur_v=other;
            if (cur_v==firstv) break;
        }
        n_loops++; if (!simple) n_nonsimple++;
        int hb = (len==3)?0:(len==4)?1:(len<=8)?2:(len<=16)?3:(len<=32)?4:5;
        hist[hb]++;

        if (len==3 && nvv>=3) {
            t3++;
            int A=vv[0],B=vv[1],C=vv[2];
            double a[3],b3[3],c3[3];
            vec3_to_d(a,&V[A]); vec3_to_d(b3,&V[B]); vec3_to_d(c3,&V[C]);
            int fill=0, fill_nok=0; double blk_depth=-1;
            for (int s=-1;s<=1;s+=2){
                double O[3];
                if (!sphere_center(a,b3,c3,rho,s,O)) continue;
                int bi=-1;
                double d=closest_other_dist(g,V,O,rho,A,B,C,&bi);
                if (d >= rho-1e-4) {            /* empty: BPA could have filled */
                    fill=1;
                    double na[3],nbv[3],ncv[3];
                    vec3_to_d(na,&N[A]); vec3_to_d(nbv,&N[B]); vec3_to_d(ncv,&N[C]);
                    double dA[3],dB[3],dC[3];
                    v_sub(dA,O,a); v_sub(dB,O,b3); v_sub(dC,O,c3);
                    if (v_dot(dA,na)>0 && v_dot(dB,nbv)>0 && v_dot(dC,ncv)>0) fill_nok=1;
                } else {
                    double depth=rho-d; if (depth>blk_depth) blk_depth=depth;
                }
            }
            if (fill) { t3_fill++; if (fill_nok) t3_fill_nok++; }
            else { t3_block++; if (blk_depth>=0){ t3_block_depth_sum+=blk_depth;
                       if (blk_depth<0.02) t3_block_border++; } }
        }
    }

    fprintf(stderr, "[bpa_holes] %ld boundary edges, %ld loops (%ld non-simple)"
            "  lengths: [3]=%ld [4]=%ld [5-8]=%ld [9-16]=%ld [17-32]=%ld [33+]=%ld\n",
            (long)nb, n_loops, n_nonsimple,
            hist[0],hist[1],hist[2],hist[3],hist[4],hist[5]);
    fprintf(stderr, "[bpa_holes] single-triangle (3-loop) pinholes=%ld:"
            "  EMPTY-ball-fillable=%ld (normal-ok=%ld)  BLOCKED=%ld"
            " (borderline<0.02=%ld, mean depth=%.3f vox)\n",
            t3, t3_fill, t3_fill_nok, t3_block, t3_block_border,
            t3_block ? t3_block_depth_sum/(double)t3_block : 0.0);

    free(ba); free(bb); free(evis); free(deg); free(off); free(adj);
}

/* --- Pinhole fill (closes the single-triangle holes BPA's one-shot front leaves) */

/* Is P strictly inside the circumcircle of coplanar A,B,C? Project A,B,C,P onto an
 * ORTHONORMAL in-plane basis (eu,ev) built from the plane normal, translate to A,
 * then run the lifted-paraboloid in-circle determinant. Two correctness points,
 * both verified by pinhole_verdict.exe --selftest against the ground-truth
 * circumcenter-distance test:
 *   - the basis MUST be orthonormal: an axis-aligned coordinate drop (the previous
 *     implementation) shears the circle into an ellipse and tests the wrong conic
 *     on tilted sheets;
 *   - inside corresponds to the winding-normalised determinant being NEGATIVE
 *     (the previous `inside > 0` was sign-inverted, so the borderline branch of
 *     pinhole_fillable filled cocircular intruders and skipped safe ones).
 * Coords are O(1 vox), so plain double is exact past our 1e-3..1e-2 vox margins.
 * nrm is any normal of the ABC plane. */
static int robust_inside_circum(const double *a, const double *b, const double *c,
                                const double *p, const double *nrm)
{
    double n[3] = { nrm[0], nrm[1], nrm[2] };
    double nl = sqrt(v_dot(n, n));
    if (nl < 1e-20) {                              /* degenerate normal: use ABC */
        double e1[3], e2[3]; v_sub(e1, b, a); v_sub(e2, c, a); v_cross(n, e1, e2);
        nl = sqrt(v_dot(n, n)); if (nl < 1e-20) return 0;
    }
    n[0] /= nl; n[1] /= nl; n[2] /= nl;
    double t[3] = { 1.0, 0.0, 0.0 };
    if (fabs(n[0]) > 0.9) { t[0] = 0.0; t[1] = 1.0; }   /* keep t away from n */
    double eu[3], ev[3]; v_cross(eu, t, n);
    double ul = sqrt(v_dot(eu, eu)); eu[0] /= ul; eu[1] /= ul; eu[2] /= ul;
    v_cross(ev, n, eu);                            /* unit: n,eu orthonormal */
    double db[3], dc[3], dp[3]; v_sub(db, b, a); v_sub(dc, c, a); v_sub(dp, p, a);
    double bx = v_dot(db,eu), by = v_dot(db,ev);
    double cx = v_dot(dc,eu), cy = v_dot(dc,ev);
    double px = v_dot(dp,eu), py = v_dot(dp,ev);
    double orient = bx*cy - by*cx;                 /* >0 : A,B,C CCW in (eu,ev) */
    double bl=bx*bx+by*by, cl=cx*cx+cy*cy, pl=px*px+py*py;
    double det = bx*(cy*pl - cl*py) - by*(cx*pl - cl*px) + bl*(cx*py - cy*px);
    double val = (orient >= 0.0) ? det : -det;     /* winding-normalised */
    return val < 0.0;                              /* <0 == strictly inside */
}

/* Should single-triangle hole ABC be filled? On the normal-consistent ball side
 * the rho-ball must be empty; a borderline-cocircular blocker (the float empty-
 * ball test mis-rejects) is resolved exactly by robust incircle. *out_flip set
 * if the face must wind (A,C,B) for orientation consistent with the normals. */
static int pinhole_fillable(const Grid *g, const Vec3 *V, const Vec3 *N, double rho,
                            int A, int B, int C, int *out_flip)
{
    double a[3],b[3],c[3], na[3],nb[3],nc[3];
    vec3_to_d(a,&V[A]); vec3_to_d(b,&V[B]); vec3_to_d(c,&V[C]);
    vec3_to_d(na,&N[A]); vec3_to_d(nb,&N[B]); vec3_to_d(nc,&N[C]);
    double nsum[3] = { na[0]+nb[0]+nc[0], na[1]+nb[1]+nc[1], na[2]+nb[2]+nc[2] };
    double e1[3],e2[3],ng[3];
    v_sub(e1,b,a); v_sub(e2,c,a); v_cross(ng,e1,e2);
    if (out_flip) *out_flip = (v_dot(ng,nsum) < 0.0);
    double ctr[3] = {(a[0]+b[0]+c[0])/3.0,(a[1]+b[1]+c[1])/3.0,(a[2]+b[2]+c[2])/3.0};
    for (int s=-1; s<=1; s+=2) {
        double O[3];
        if (!sphere_center(a,b,c,rho,s,O)) continue;
        double oc[3]; v_sub(oc,O,ctr);
        if (v_dot(oc,nsum) <= 0.0) continue;           /* not the surface side */
        int blk = -1;
        double d = closest_other_dist(g,V,O,rho,A,B,C,&blk);
        if (d >= rho - 1e-4) return 1;                 /* empty: fill */
        if (blk >= 0) {                                /* borderline: robust test */
            double p[3]; vec3_to_d(p,&V[blk]);
            if (!robust_inside_circum(a,b,c,p,nsum)) return 1;
        }
    }
    return 0;
}

/* Trace boundary 3-loops and fill each whose ball is empty + normal-consistent.
 * Single pass: 3-loops are edge-disjoint, so filling one never exposes another. */
static int fill_pinholes(EdgeStore *es, const Grid *g, const Vec3 *V, const Vec3 *N,
                         double rho, BpaBuild *b, int *vfront, int n)
{
    int nb = 0;
    for (int i=0;i<es->n;i++) if (es->e[i].state==ES_BOUNDARY) nb++;
    if (nb==0) return 0;
    int *ba=(int*)malloc((size_t)nb*sizeof(int));
    int *bbv=(int*)malloc((size_t)nb*sizeof(int));
    int *bid=(int*)malloc((size_t)nb*sizeof(int));
    uint8_t *evis=(uint8_t*)calloc((size_t)nb,1);
    { int k=0; for(int i=0;i<es->n;i++) if(es->e[i].state==ES_BOUNDARY)
        { ba[k]=es->e[i].va; bbv[k]=es->e[i].vb; bid[k]=i; k++; } }
    int *deg=(int*)calloc((size_t)n,sizeof(int));
    for(int e=0;e<nb;e++){ deg[ba[e]]++; deg[bbv[e]]++; }
    int *off=(int*)malloc((size_t)(n+1)*sizeof(int)); off[0]=0;
    for(int i=0;i<n;i++) off[i+1]=off[i]+deg[i];
    int *adj=(int*)malloc((size_t)off[n]*sizeof(int));
    int *cur=(int*)malloc((size_t)n*sizeof(int)); memcpy(cur,off,(size_t)n*sizeof(int));
    for(int e=0;e<nb;e++){ adj[cur[ba[e]]++]=e; adj[cur[bbv[e]]++]=e; }
    free(cur);

    int filled=0;
    for(int e0=0;e0<nb;e0++){
        if(evis[e0]) continue;
        int le[4], lv[4], len=1, simple=1, guard=0;
        int cur_e=e0, firstv=ba[e0], cur_v=bbv[e0];
        lv[0]=firstv; le[0]=e0; evis[e0]=1;
        if(deg[firstv]!=2) simple=0;
        while(guard++ < nb+4){
            if(deg[cur_v]!=2) simple=0;
            int nexte=-1;
            for(int t=off[cur_v]; t<off[cur_v+1]; t++){ int ee=adj[t]; if(ee==cur_e||evis[ee]) continue; nexte=ee; break; }
            if(nexte<0) break;
            evis[nexte]=1;
            if(len<4){ lv[len]=cur_v; le[len]=nexte; }
            len++;
            int other=(ba[nexte]==cur_v)?bbv[nexte]:ba[nexte];
            cur_e=nexte; cur_v=other;
            if(cur_v==firstv) break;
        }
        if(len==3 && simple){
            int A=lv[0], B=lv[1], C=lv[2];
            if(A!=B && B!=C && A!=C){
                int flip=0;
                if(pinhole_fillable(g,V,N,rho,A,B,C,&flip)){
                    if(flip) bpa_add_face(b, A, C, B); else bpa_add_face(b, A, B, C);
                    edge_set_state(es, vfront, bid[le[0]], ES_INTERIOR);
                    edge_set_state(es, vfront, bid[le[1]], ES_INTERIOR);
                    edge_set_state(es, vfront, bid[le[2]], ES_INTERIOR);
                    filled++;
                }
            }
        }
    }
    free(ba); free(bbv); free(bid); free(evis); free(deg); free(off); free(adj);
    return filled;
}

int BallPivot_reconstruct(Arena_T arena,
                          const float *verts,
                          const float *normals,
                          size_t nv,
                          float rho_f,
                          int32_t **out_faces,
                          size_t *out_nf)
{
    *out_faces = NULL;
    *out_nf = 0;
    if (nv < 3) return 0;
    double rho = (double)rho_f;
    /* Experimental radius override (BPA_RHO): sweep the pivot radius without
     * recompiling, to test whether a larger ball reduces front-collision
     * pinholes. Falls back to the caller's rho. Caller still owns the
     * inter-wrap-clearance guarantee — do not set this above the safe cap. */
    { const char *e = getenv("BPA_RHO"); if (e) { double r = atof(e); if (r > 0.0) rho = r; } }
    g_no_normal_gate = (getenv("BPA_NO_NORMAL_GATE") != NULL);
    g_normal_tol = 0.0;
    { const char *e = getenv("BPA_NORMAL_TOL"); if (e) { double t = atof(e); if (t >= 0.0 && t < 1.0) g_normal_tol = t; } }
    g_ko_empty = 0; g_ko_normal = 0;
    g_dead_gap = g_dead_empty = g_dead_normal = g_dead_mixed = g_dead_theta = 0;
    g_dead_measure = (getenv("BPA_DEBUG") != NULL);
    g_mix_erej = g_mix_anti = g_mix_nrej = 0;
    g_depth_sum = g_past90_sum = 0; g_absh_sum = 0;
    for (int i = 0; i < NBD; i++) { g_hd_all[i] = 0; g_hd_anti[i] = 0; }
    for (int i = 0; i < NBA; i++) g_ha[i] = 0;
    for (int i = 0; i < NBH; i++) g_hh[i] = 0;

    const Vec3 *V = (const Vec3 *)verts;
    const Vec3 *N = (const Vec3 *)normals;
    int n = (int)nv;

    Grid *g = grid_build(V, n, 2.0*rho);
    EdgeStore *es = edges_new(n*3);
    uint8_t *used = (uint8_t *)calloc((size_t)n, 1);
    int *vfront = (int *)calloc((size_t)n, sizeof(int));
    BpaBuild b; bpa_build_init(&b, n);

    int seed_scan = 0;
    while (seed_scan < n) {
        int seed_face[3]; double seed_O[3];
        int found = 0;
        for (; seed_scan < n && !found; seed_scan++) {
            if (used[seed_scan]) continue;
            CandList cl;
            cl.cap = 64; cl.n_cands = 0;
            cl.cands = (int *)malloc((size_t)cl.cap * sizeof(int));
            grid_query(g, V[seed_scan].z, V[seed_scan].y, V[seed_scan].x,
                       2.0*rho, seed_collect, &cl, V);
            sort_by_dist(cl.cands, cl.n_cands, V, seed_scan);
            for (int aa = 0; aa < cl.n_cands && !found; aa++) {
                int ja = cl.cands[aa];
                if (ja == seed_scan || used[ja]) continue;
                for (int bb = aa+1; bb < cl.n_cands && !found; bb++) {
                    int jb = cl.cands[bb];
                    if (jb == seed_scan || used[jb]) continue;
                    if (try_seed(g, V, N, seed_scan, ja, jb, rho,
                                 seed_face, seed_O)) found = 1;
                }
            }
            free(cl.cands);
        }
        if (!found) break;

        bpa_add_face(&b, seed_face[0], seed_face[1], seed_face[2]);
        used[seed_face[0]] = used[seed_face[1]] = used[seed_face[2]] = 1;

        int e01 = edges_insert(es, vfront, seed_face[0], seed_face[1], seed_face[2], seed_O);
        int e12 = edges_insert(es, vfront, seed_face[1], seed_face[2], seed_face[0], seed_O);
        int e20 = edges_insert(es, vfront, seed_face[2], seed_face[0], seed_face[1], seed_O);
        bpa_enqueue(&b, e01); bpa_enqueue(&b, e12); bpa_enqueue(&b, e20);

        bpa_grow(g, V, N, rho, es, used, vfront, &b, 0 /*strict bowtie*/);
    }

    /* --- Pinhole-fill pass. BPA's front is one-shot: a dead-ended edge is retired
     * as ES_BOUNDARY and never re-closed, even after the surrounding mesh fills in
     * from other directions and leaves a single-triangle hole with a genuinely
     * empty, normal-consistent ball (~91% of the leftover pinholes, measured).
     * Re-running the front can't re-acquire them — the pivot only ADVANCES a front
     * from one committed side, it does not re-close an interior triangle between
     * three already-committed neighbours (verified: re-arm fills 0). So close them
     * directly: trace boundary 3-loops and add the triangle when its ball is empty
     * (robust Shewchuk incircle resolves the borderline-cocircular blockers the
     * float empty-ball test mis-rejects). Topology-safe: empty ball on a single
     * layer => no inter-wrap bridge; each filled edge goes boundary->interior.
     * Disable with BPA_NO_REARM. */
    if (!getenv("BPA_NO_REARM")) {
        int pf = fill_pinholes(es, g, V, N, rho, &b, vfront, n);
        if (getenv("BPA_DEBUG"))
            fprintf(stderr, "[bpa_recon] pinhole-fill: closed %d single-triangle holes\n", pf);
    }

    if (getenv("BPA_DEBUG")) {
        int n_front = 0, n_int = 0, n_bnd = 0;
        for (int i = 0; i < es->n; i++) {
            if (es->e[i].state == ES_FRONT) n_front++;
            else if (es->e[i].state == ES_INTERIOR) n_int++;
            else n_bnd++;
        }
        fprintf(stderr, "[bpa_recon] nv=%d rho=%.3f faces=%d edges: "
                "front=%d interior=%d boundary=%d | pivot-fail: "
                "no_pivot=%d bowtie=%d orient=%d 3fan=%d fold=%d | "
                "pivot-ko: empty=%ld normal=%ld tol=%.2f%s\n",
                n, rho, b.nf, n_front, n_int, n_bnd,
                b.dbg_no_pivot, b.dbg_bowtie, b.dbg_orient, b.dbg_3fan, b.dbg_fold,
                g_ko_empty, g_ko_normal, g_normal_tol,
                g_no_normal_gate ? " [NORMAL GATE OFF]" : "");
        fprintf(stderr, "[bpa_recon] no_pivot dead-ends by cause: "
                "gap=%ld empty-ball=%ld normal-gate=%ld mixed=%ld theta-degenerate=%ld\n",
                g_dead_gap, g_dead_empty, g_dead_normal, g_dead_mixed, g_dead_theta);
        if (g_mix_erej > 0 || g_mix_nrej > 0) {
            static const char *dl[NBD] = {"<0.02","0.02-.05",".05-0.1","0.1-0.2",
                                          "0.2-0.4","0.4-0.8","0.8-1.6"," >1.6 "};
            static const char *al[NBA] = {"<1deg","1-2","2-5","5-10","10-20","20-45"," >45 "};
            fprintf(stderr, "[bpa_recon] MIXED empty-ball (+normal-side cands)=%ld"
                    "  blocker anti-aligned (other layer)=%ld (%.0f%%)  mean depth=%.3f vox\n",
                    g_mix_erej, g_mix_anti,
                    g_mix_erej ? 100.0*(double)g_mix_anti/(double)g_mix_erej : 0.0,
                    g_mix_erej ? g_depth_sum/(double)g_mix_erej : 0.0);
            fprintf(stderr, "[bpa_recon]   depth(rho shrink to clear intruder) | all / anti-aligned blocker:\n");
            for (int i = 0; i < NBD; i++)
                fprintf(stderr, "[bpa_recon]     %-9s %6ld / %-6ld\n", dl[i], g_hd_all[i], g_hd_anti[i]);
            static const char *hl[NBH] = {"<-0.4","-0.4/-0.2","-0.2/-.05","-.05/.05",
                                          ".05/0.2","0.2/0.4","0.4/0.8","0.8/1.6"," >1.6 "};
            fprintf(stderr, "[bpa_recon]   blocker signed height above cand tangent plane"
                    " (mean|h|=%.3f vox):\n",
                    g_mix_erej ? g_absh_sum/(double)g_mix_erej : 0.0);
            for (int i = 0; i < NBH; i++)
                fprintf(stderr, "[bpa_recon]     %-10s %6ld\n", hl[i], g_hh[i]);
            fprintf(stderr, "[bpa_recon] MIXED normal-gate rejects=%ld  mean past-90=%.1f deg | bins:\n",
                    g_mix_nrej, g_mix_nrej ? g_past90_sum/(double)g_mix_nrej : 0.0);
            for (int i = 0; i < NBA; i++)
                fprintf(stderr, "[bpa_recon]     %-7s %6ld\n", al[i], g_ha[i]);
        }
        analyze_holes(es, g, V, N, rho, n);
    }

    int32_t *out = (int32_t *)ARENA_ALLOC(arena, (long)((size_t)b.nf * 3 * sizeof(int32_t)));
    for (int i = 0; i < b.nf*3; i++) out[i] = (int32_t)b.F[i];
    *out_faces = out;
    *out_nf = (size_t)b.nf;

    bpa_build_free(&b); free(used); free(vfront);
    edges_free(es); grid_free(g);
    return 0;
}

int BallPivot_bridge(Arena_T arena,
                     const float *verts,
                     const float *normals,
                     size_t nv,
                     float rho_f,
                     float rho_max_f,
                     const BpaInitEdge *init_edges,
                     size_t n_init,
                     int32_t **out_faces,
                     size_t *out_nf)
{
    *out_faces = NULL;
    *out_nf = 0;
    if (nv < 3 || n_init == 0) return 0;
    double rho = (double)rho_f;
    double rho_max = (double)rho_max_f;
    if (rho_max < rho) rho_max = rho;   /* rho_max==rho => single radius */

    const Vec3 *V = (const Vec3 *)verts;
    const Vec3 *N = (const Vec3 *)normals;
    int n = (int)nv;
    g_bpa_trace = (getenv("BPA_TRACE") != NULL);
    g_no_normal_gate = 0;   /* the seam bridge always honors the normal gate */
    g_normal_tol = 0.0;

    Grid *g = grid_build(V, n, 2.0*rho);
    EdgeStore *es = edges_new((int)n_init * 4);
    /* `used` only gates the seed scan, which the bridge never runs; keep
     * it so bpa_grow's `used[v_new]=1` writes land somewhere valid. */
    uint8_t *used = (uint8_t *)calloc((size_t)n, 1);
    int *vfront = (int *)calloc((size_t)n, sizeof(int));
    BpaBuild b; bpa_build_init(&b, (int)n_init * 2);

    /* Prime the front with every supplied boundary half-edge. The hinge
     * ball center O_prev is reconstructed from the edge's single existing
     * face (va, vb, v_opp), placed on the outward side (same normal-sign
     * convention as try_seed) so the pivot rolls away from the existing
     * surface, into the gap. */
    for (size_t k = 0; k < n_init; k++) {
        int va = init_edges[k].va;
        int vb = init_edges[k].vb;
        int vo = init_edges[k].v_opp;
        if (va < 0 || vb < 0 || vo < 0 ||
            va >= n || vb >= n || vo >= n) continue;
        if (va == vb || vb == vo || va == vo) continue;
        if (edges_find(es, va, vb) >= 0) continue;  /* already primed */
        double a[3], bb_[3], c[3];
        vec3_to_d(a, &V[va]); vec3_to_d(bb_, &V[vb]); vec3_to_d(c, &V[vo]);
        double na[3], nb[3], nc[3], n_avg[3];
        vec3_to_d(na, &N[va]); vec3_to_d(nb, &N[vb]); vec3_to_d(nc, &N[vo]);
        n_avg[0]=na[0]+nb[0]+nc[0];
        n_avg[1]=na[1]+nb[1]+nc[1];
        n_avg[2]=na[2]+nb[2]+nc[2];
        double ab[3], ac[3], tnorm[3];
        v_sub(ab, bb_, a); v_sub(ac, c, a); v_cross(tnorm, ab, ac);
        int sign = (v_dot(tnorm, n_avg) > 0) ? +1 : -1;
        double O_prev[3];
        if (!sphere_center(a, bb_, c, rho, sign, O_prev)) continue;
        int eid = edges_insert(es, vfront, va, vb, vo, O_prev);
        bpa_enqueue(&b, eid);
    }

    size_t primed = (size_t)b.qt;

    /* Grow at escalating ball radii (Bernardini 1999 §IV-F). At each radius:
     * drain the front, then re-arm every boundary edge and re-grow to a
     * fixpoint -- a failed-pivot edge often succeeds once the OPPOSITE front
     * has rolled in and deposited neighbors. Then bump the radius to span gaps
     * the smaller ball fell through. All radii share one edge store, so the
     * directed glue keeps the result manifold/orientable across passes. rho_max
     * caps the largest ball so 2*rho stays below the inter-wrap clearance.
     * Face count only increases and is bounded, so this terminates. */
    double radii[3]; int n_rho = 0;
    {
        double mm[3] = { 1.0, 1.5, 2.0 }, prev = -1.0;
        for (int m = 0; m < 3; m++) {
            double r = rho * mm[m]; if (r > rho_max) r = rho_max;
            if (r > prev + 1e-6) { radii[n_rho++] = r; prev = r; }
        }
    }
    for (int lvl = 0; lvl < n_rho; lvl++) {
        double r = radii[lvl];
        for (int pass = 0; pass < 16; pass++) {
            if (!(lvl == 0 && pass == 0)) {
                int rearmed = 0;
                for (int i = 0; i < es->n; i++) {
                    if (es->e[i].state == ES_BOUNDARY) {
                        edge_set_state(es, vfront, i, ES_FRONT);
                        bpa_enqueue(&b, i);
                        rearmed++;
                    }
                }
                if (rearmed == 0) break;   /* nothing open at this radius */
            }
            int before = b.nf;
            bpa_grow(g, V, N, r, es, used, vfront, &b, 1 /*relax bowtie: zip fronts*/);
            if (!(lvl == 0 && pass == 0) && b.nf == before) break;
        }
    }

    if (getenv("BPA_DEBUG")) {
        int n_front = 0, n_int = 0, n_bnd = 0;
        for (int i = 0; i < es->n; i++) {
            if (es->e[i].state == ES_FRONT) n_front++;
            else if (es->e[i].state == ES_INTERIOR) n_int++;
            else n_bnd++;
        }
        fprintf(stderr, "[bpa_bridge] nv=%d primed_init=%zu faces=%d "
                "edges: front=%d interior=%d boundary=%d | "
                "pivot-fail: no_pivot=%d bowtie=%d orient=%d 3fan=%d fold=%d bowtie_ok=%d\n",
                n, primed, b.nf, n_front, n_int, n_bnd,
                b.dbg_no_pivot, b.dbg_bowtie, b.dbg_orient, b.dbg_3fan, b.dbg_fold, b.dbg_bowtie_ok);
    }

    int32_t *out = (int32_t *)ARENA_ALLOC(arena, (long)((size_t)b.nf * 3 * sizeof(int32_t)));
    for (int i = 0; i < b.nf*3; i++) out[i] = (int32_t)b.F[i];
    *out_faces = out;
    *out_nf = (size_t)b.nf;

    bpa_build_free(&b); free(used); free(vfront);
    edges_free(es); grid_free(g);
    return 0;
}
