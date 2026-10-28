#include "rvd_clip.h"
#include "../common/kdtree.h"
#include "../common/attene_predicates_wrap.h"
#include "../common/cycle_timer.h"

#include <float.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef _OPENMP
#include <omp.h>
#endif

/* Clip-polygon vertex: an AttenePt id + a tag for the edge that FOLLOWS it
 * (vertex i -> i+1). etag 0 = edge lies on a mesh edge with explicit endpoints
 * (me0,me1); etag 1 = edge lies on bisector(owner, bsj). The tag selects the exact
 * construction when the edge is cut. (fi,fj,fk) name the three sites when the
 * vertex itself is a facet-bisector2 Voronoi vertex (fi >= 0), else fi = -1. */
typedef struct {
    int32_t v; int etag; int32_t me0, me1, bsj;
    int32_t fi, fj, fk;
} PV;

#define RVD_MAX_POLY 256
#define RVD_MAX_CAND 1024

struct Rvd_T {
    Arena_T  arena;
    double  *verts;   /* [nv*3] */
    int32_t *faces;   /* [nf*3] */
    size_t   nv, nf;
    double   mesh_area;
    RvdStats *stats;  /* optional instrumentation sink (NULL = off) */
    int32_t *vnear;   /* [nv] nearest site of each vertex in the last pass (-1 = none):
                       * the next pass's exact bounded search starts from it */
};

static inline void vsub3(const double a[3], const double b[3], double o[3]) {
    o[0]=a[0]-b[0]; o[1]=a[1]-b[1]; o[2]=a[2]-b[2];
}
static inline void vcross3(const double a[3], const double b[3], double o[3]) {
    o[0]=a[1]*b[2]-a[2]*b[1]; o[1]=a[2]*b[0]-a[0]*b[2]; o[2]=a[0]*b[1]-a[1]*b[0];
}
static inline double vdot3(const double a[3], const double b[3]) {
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}
static inline double vlen3(const double a[3]) { return sqrt(vdot3(a,a)); }
static inline double tri_area(const double A[3], const double B[3], const double C[3]) {
    double ab[3], ac[3], cr[3];
    vsub3(B,A,ab); vsub3(C,A,ac); vcross3(ab,ac,cr);
    return 0.5*vlen3(cr);
}

Rvd_T Rvd_new(Arena_T arena, const double *verts, size_t nv,
              const int32_t *faces, size_t nf) {
    if (!arena || !verts || !faces || nv < 3 || nf < 1) return NULL;
    Rvd_T r = (Rvd_T)ARENA_ALLOC(arena, sizeof(*r));
    r->arena = arena; r->nv = nv; r->nf = nf;
    r->verts = (double *)ARENA_ALLOC(arena, (size_t)(nv*3*sizeof(double)));
    memcpy(r->verts, verts, nv*3*sizeof(double));
    r->faces = (int32_t *)ARENA_ALLOC(arena, (size_t)(nf*3*sizeof(int32_t)));
    memcpy(r->faces, faces, nf*3*sizeof(int32_t));
    r->mesh_area = 0.0;
    r->stats = NULL;
    r->vnear = (int32_t *)ARENA_ALLOC(arena, (size_t)(nv*sizeof(int32_t)));
    for (size_t i = 0; i < nv; i++) r->vnear[i] = -1;
    for (size_t f=0; f<nf; f++) {
        const double *A=&verts[(size_t)faces[f*3+0]*3];
        const double *B=&verts[(size_t)faces[f*3+1]*3];
        const double *C=&verts[(size_t)faces[f*3+2]*3];
        r->mesh_area += tri_area(A,B,C);
    }
    return r;
}

double Rvd_mesh_area(Rvd_T r) { return r ? r->mesh_area : 0.0; }

void Rvd_set_stats(Rvd_T r, RvdStats *stats) { if (r) r->stats = stats; }

/* Intersection of the edge following `cur` with bisector(si,sk); sets (fi,fj,fk)
 * to the three sites when the result is a facet-bisector2 Voronoi vertex, else -1. */
static int32_t make_cross(AttenePredCtx *ctx, const PV *cur, int32_t si, int32_t sk,
                          int32_t va, int32_t vb, int32_t vc,
                          int32_t *fi, int32_t *fj, int32_t *fk) {
    if (cur->etag == 0) {
        *fi = *fj = *fk = -1;
        return AttenePred_edge_bisector(ctx, cur->me0, cur->me1, si, sk);
    }
    *fi = si; *fj = cur->bsj; *fk = sk;
    return AttenePred_facet_bisector2(ctx, va, vb, vc, si, cur->bsj, sk);
}

/* Whole-triangle bisector filter (2026-10-05).  The side of bisector(si,sk) is
 * the sign of the affine function |P-SI|^2 - |P-SK|^2, and every polygon vertex
 * -- explicit or constructed -- lies on mesh triangle ABC.  When A, B and C are
 * all CERTAINLY strictly closer to SI, every vertex tests "inside" exactly and
 * the plane leaves the polygon untouched, so clip_cell may skip it without any
 * exact predicate.  Evaluated in double: the rounding error of da - db is below
 * 6u(da+db) (u = 2^-53), the margin used is 1e-12 (da+db) (~1500x); NaN or any
 * doubt falls through to the exact path.  0 disables it (Rvd_selftest reference). */
static int g_rvd_tri_filter = 1;

static inline int rvd_tri_inside(const double A[3], const double B[3], const double C[3],
                                 const double SI[3], const double SK[3]) {
    const double *P[3] = { A, B, C };
    for (int k = 0; k < 3; k++) {
        double a0 = P[k][0]-SI[0], a1 = P[k][1]-SI[1], a2 = P[k][2]-SI[2];
        double b0 = P[k][0]-SK[0], b1 = P[k][1]-SK[1], b2 = P[k][2]-SK[2];
        double da = a0*a0 + a1*a1 + a2*a2;
        double db = b0*b0 + b1*b1 + b2*b2;
        if (!(da - db < -1e-12 * (da + db))) return 0;
    }
    return 1;
}

/* Clip mesh triangle (A,B,C) to the Voronoi cell of owner `si`, against the
 * candidate sites. Returns the vertex count of the clipped convex polygon in
 * `poly` (>= 3), or a value < 3 if the cell is empty on this triangle. */
static int clip_cell(AttenePredCtx *ctx, const double *sites,
                     const int32_t *cand, const double *cxyz, size_t ncand, int32_t si,
                     const double A[3], const double B[3], const double C[3],
                     int32_t va, int32_t vb, int32_t vc, PV *poly,
                     uint64_t *n_planes, uint64_t *n_preds) {
    const double *SI = &sites[(size_t)si*3];
    double da[3], db[3], dc[3];
    vsub3(A,SI,da); vsub3(B,SI,db); vsub3(C,SI,dc);
    double la=vlen3(da), lb=vlen3(db), lc=vlen3(dc);
    double Ri = la>lb ? (la>lc?la:lc) : (lb>lc?lb:lc);
    double cut2 = (2.0*Ri)*(2.0*Ri);

    PV tmp[RVD_MAX_POLY];
    int np = 3;
    poly[0] = (PV){ va,0, va,vb,-1, -1,-1,-1 };
    poly[1] = (PV){ vb,0, vb,vc,-1, -1,-1,-1 };
    poly[2] = (PV){ vc,0, vc,va,-1, -1,-1,-1 };

    /* cxyz: the candidates' coordinates gathered once per face, as three rows
     * of RVD_MAX_CAND (exact copies of sites[]): the filter below reads them in
     * order instead of a random sites[] row per candidate per owner. */
    const double *cx = cxyz, *cy = cxyz + RVD_MAX_CAND, *cz = cxyz + 2*RVD_MAX_CAND;
    for (size_t di = 0; di < ncand; di++) {
        int32_t sk = cand[di];
        if (sk == si) continue;
        double d0 = cx[di]-SI[0], d1 = cy[di]-SI[1], d2 = cz[di]-SI[2];
        if (d0*d0 + d1*d1 + d2*d2 > cut2) continue;      /* = vdot3(SK - SI) */
        const double *SK = &sites[(size_t)sk*3];
        if (g_rvd_tri_filter && rvd_tri_inside(A, B, C, SI, SK)) continue;
        if (n_planes) { (*n_planes)++; *n_preds += (uint64_t)np; }

        /* Each vertex's side of this bisector, evaluated once: the loop below
         * reads it as `cur` and again as the previous edge's `nxt`.  The
         * predicate is a pure function of the (immutable) point definition,
         * and the make_cross sequence is unchanged, so this only halves the
         * exact-predicate calls (2026-10-04). */
        int side[RVD_MAX_POLY];
        AttenePt vid[RVD_MAX_POLY];
        for (int i = 0; i < np; i++) vid[i] = poly[i].v;
        AttenePred_side_bisector_batch(ctx, vid, np, si, sk, side);
        for (int i = 0; i < np; i++) side[i] = (side[i] <= 0);

        int no = 0;
        for (int i = 0; i < np; i++) {
            PV cur = poly[i];
            int cin = side[i];
            int nin = side[(i+1)%np];
            if (cin) {
                if (no >= RVD_MAX_POLY) return -1;
                tmp[no++] = cur;
                if (!nin) {
                    int32_t fi, fj, fk;
                    int32_t X = make_cross(ctx, &cur, si, sk, va, vb, vc, &fi, &fj, &fk);
                    if (no >= RVD_MAX_POLY) return -1;
                    tmp[no++] = (PV){ X,1, -1,-1,sk, fi,fj,fk };   /* follow: clip line */
                }
            } else if (nin) {
                int32_t fi, fj, fk;
                int32_t X = make_cross(ctx, &cur, si, sk, va, vb, vc, &fi, &fj, &fk);
                if (no >= RVD_MAX_POLY) return -1;
                tmp[no++] = (PV){ X, cur.etag, cur.me0, cur.me1, cur.bsj, fi,fj,fk }; /* resume */
            }
        }
        memcpy(poly, tmp, (size_t)no*sizeof(PV));
        np = no;
        if (np < 3) return np;
    }
    return np;
}

/* Tight candidate radius (2026-10-06).  A site sk that clip_cell can use for
 * owner si passes |SK - SI| <= 2 Ri (cut2), Ri = si's farthest triangle vertex.
 * An owner's cell holds a triangle point P, so with c the centroid, rc its
 * farthest vertex and d_near the centroid's nearest-site distance:
 *   |SI - P| <= d_near + rc,  |SI - c| <= d_near + 2 rc,  Ri <= d_near + 3 rc,
 *   |SK - c| <= 2 Ri + |SI - c| <= 3 d_near + 8 rc.
 * (The wide radius is the same bound with the longest edge for rc; rc <= 2/3 of
 * it.)  A ball query's list is a depth-first pre-order whose child order depends
 * only on the query point, and a smaller radius prunes only subtrees without a
 * point inside it, so the tight list is the wide list restricted to the tight
 * ball, in the same order: every owner meets the same cutting candidates in the
 * same sequence, and the clipping is bitwise the same.  That holds only while
 * the wide list was complete (under RVD_MAX_CAND), which a site-count box bound
 * certifies; otherwise the wide radius is kept.  0 = always wide (selftest). */
static int g_rvd_tight = 1;

/* Site counts per grid cell with inclusive 3-D prefix sums: an O(1) count of the
 * sites inside an axis box, an upper bound on any ball inside it. */
typedef struct {
    float lo[3], inv_h;
    int n[3];
    uint32_t *ps;            /* (n0+1)(n1+1)(n2+1), zero first plane/row/column */
} RvdSiteGrid;

static inline int rvd_grid_cell(const RvdSiteGrid *g, int a, float x) {
    float t = (x - g->lo[a]) * g->inv_h;
    int i = t <= 0.0f ? 0 : (t >= (float)g->n[a] ? g->n[a] - 1 : (int)t);
    return i < g->n[a] ? i : g->n[a] - 1;
}

static int rvd_grid_build(Arena_T arena, const float *sf, size_t ns, RvdSiteGrid *g) {
    g->ps = NULL;
    if (ns == 0) return -1;
    float lo[3] = { sf[0], sf[1], sf[2] }, hi[3] = { sf[0], sf[1], sf[2] };
    for (size_t i = 1; i < ns; i++)
        for (int a = 0; a < 3; a++) {
            float x = sf[i*3+a];
            if (x < lo[a]) lo[a] = x;
            if (x > hi[a]) hi[a] = x;
        }
    float ext = 0.0f;
    for (int a = 0; a < 3; a++) if (hi[a] - lo[a] > ext) ext = hi[a] - lo[a];
    if (!(ext < 1e7f)) return -1;                       /* NaN or absurd extent: no bound */
    float h = ext / 48.0f;
    if (h < 2.0f) h = 2.0f;
    g->inv_h = 1.0f / h;
    for (int a = 0; a < 3; a++) {
        g->lo[a] = lo[a];
        g->n[a] = (int)((hi[a] - lo[a]) * g->inv_h) + 1;
        if (g->n[a] < 1) g->n[a] = 1;
        if (g->n[a] > 64) g->n[a] = 64;
    }
    size_t s1 = (size_t)g->n[2] + 1, s0 = s1 * ((size_t)g->n[1] + 1);
    size_t np = s0 * ((size_t)g->n[0] + 1);
    g->ps = (uint32_t *)ARENA_CALLOC(arena, np, sizeof(uint32_t));
    for (size_t i = 0; i < ns; i++) {
        size_t c0 = (size_t)rvd_grid_cell(g, 0, sf[i*3+0]) + 1;
        size_t c1 = (size_t)rvd_grid_cell(g, 1, sf[i*3+1]) + 1;
        size_t c2 = (size_t)rvd_grid_cell(g, 2, sf[i*3+2]) + 1;
        g->ps[c0*s0 + c1*s1 + c2]++;
    }
    for (size_t i = 1; i <= (size_t)g->n[0]; i++)          /* prefix along each axis */
        for (size_t j = 0; j <= (size_t)g->n[1]; j++)
            for (size_t k = 0; k <= (size_t)g->n[2]; k++) g->ps[i*s0 + j*s1 + k] += g->ps[(i-1)*s0 + j*s1 + k];
    for (size_t i = 0; i <= (size_t)g->n[0]; i++)
        for (size_t j = 1; j <= (size_t)g->n[1]; j++)
            for (size_t k = 0; k <= (size_t)g->n[2]; k++) g->ps[i*s0 + j*s1 + k] += g->ps[i*s0 + (j-1)*s1 + k];
    for (size_t i = 0; i <= (size_t)g->n[0]; i++)
        for (size_t j = 0; j <= (size_t)g->n[1]; j++)
            for (size_t k = 1; k <= (size_t)g->n[2]; k++) g->ps[i*s0 + j*s1 + k] += g->ps[i*s0 + j*s1 + k - 1];
    return 0;
}

/* Sites in the cells meeting the box c +- half (cells are found by the same
 * monotone float map as the sites, so every site inside the box is counted). */
static size_t rvd_grid_box_count(const RvdSiteGrid *g, const float c[3], float half) {
    int a0[3], a1[3];
    for (int a = 0; a < 3; a++) {
        a0[a] = rvd_grid_cell(g, a, c[a] - half);
        a1[a] = rvd_grid_cell(g, a, c[a] + half);
    }
    size_t s1 = (size_t)g->n[2] + 1, s0 = s1 * ((size_t)g->n[1] + 1);
    size_t i0 = (size_t)a0[0], i1 = (size_t)a1[0] + 1, j0 = (size_t)a0[1], j1 = (size_t)a1[1] + 1;
    size_t k0 = (size_t)a0[2], k1 = (size_t)a1[2] + 1;
    int64_t v = (int64_t)g->ps[i1*s0 + j1*s1 + k1] - (int64_t)g->ps[i0*s0 + j1*s1 + k1]
              - (int64_t)g->ps[i1*s0 + j0*s1 + k1] - (int64_t)g->ps[i1*s0 + j1*s1 + k0]
              + (int64_t)g->ps[i0*s0 + j0*s1 + k1] + (int64_t)g->ps[i0*s0 + j1*s1 + k0]
              + (int64_t)g->ps[i1*s0 + j0*s1 + k0] - (int64_t)g->ps[i0*s0 + j0*s1 + k0];
    return v < 0 ? 0 : (size_t)v;
}

/* Per-triangle candidate set: a provable superset of owners + cutting neighbours.
 * Returns the count written into cand. */
static size_t triangle_candidates(KDTree_T kt, const RvdSiteGrid *grid,
                                  const float *sf, size_t nsites,
                                  const int32_t near3[3], const double A[3], const double B[3],
                                  const double C[3], int32_t *cand, int *wide) {
    double cz[3] = { (A[0]+B[0]+C[0])/3.0, (A[1]+B[1]+C[1])/3.0, (A[2]+B[2]+C[2])/3.0 };
    double e0[3], e1[3], e2[3];
    vsub3(B,A,e0); vsub3(C,B,e1); vsub3(A,C,e2);
    double le0=vlen3(e0), le1=vlen3(e1), le2=vlen3(e2);
    double maxe = le0>le1 ? (le0>le2?le0:le2) : (le1>le2?le1:le2);
    float czf[3] = { (float)cz[0], (float)cz[1], (float)cz[2] };
    /* d_near: the bounded search seeded by the vertices' nearest sites returns
     * KDTree_nearest's distance exactly (see KDTree_nearest_bounded). */
    float dnn_sq = 0.0f, bnd = FLT_MAX;
    for (int t = 0; t < 3; t++) {
        if (near3[t] < 0 || (size_t)near3[t] >= nsites) continue;
        const float *s = &sf[(size_t)near3[t]*3];
        float g0 = czf[0] - s[0], g1 = czf[1] - s[1], g2 = czf[2] - s[2];
        float d = g0 * g0 + g1 * g1 + g2 * g2;
        if (d < bnd) bnd = d;
    }
    if (bnd < FLT_MAX)
        (void)KDTree_nearest_bounded(kt, czf, nextafterf(bnd, FLT_MAX), &dnn_sq);
    else
        (void)KDTree_nearest(kt, czf, &dnn_sq);
    double d_near = sqrt((double)dnn_sq);
    double R = 3.0*d_near + 8.0*maxe;
    R = R*1.01 + 1.0;                        /* float-kdtree safety */
    *wide = 1;
    if (g_rvd_tight && grid && grid->ps) {
        /* the wide ball sits inside the box c +- (R + pad): fewer sites there
         * than the cap means its list was complete (see g_rvd_tight) */
        float half = (float)(R * 1.01 + 0.01);
        if (rvd_grid_box_count(grid, czf, half) < (size_t)RVD_MAX_CAND) {
            double rc = 0.0, q[3];
            const double *P3[3] = { A, B, C };
            for (int t = 0; t < 3; t++) { vsub3(P3[t], cz, q); double l = vlen3(q); if (l > rc) rc = l; }
            double Rt = (3.0*d_near + 8.0*rc)*1.01 + 1.0;
            if (Rt < R) R = Rt;
            *wide = 0;
        }
    }
    return KDTree_ball_query(kt, czf, (float)(R*R), cand, RVD_MAX_CAND);
}

/* qsort helper: compare a "sorted triple" key stored as int32[6] = {key3, tri3}. */
static int cmp_triple(const void *pa, const void *pb) {
    const int32_t *a = (const int32_t *)pa, *b = (const int32_t *)pb;
    for (int i = 0; i < 3; i++) if (a[i] != b[i]) return a[i] < b[i] ? -1 : 1;
    return 0;
}

/* BFS owner-queue push: enqueue site s for the current triangle if it is unvisited
 * (per-triangle generation stamp) and the queue has room. */
static inline void owner_push(int32_t s, int32_t cur_gen, int32_t *gen,
                              int32_t *q, int *qt) {
    if (s >= 0 && *qt < RVD_MAX_CAND && gen[s] != cur_gen) {
        gen[s] = cur_gen;
        q[(*qt)++] = s;
    }
}

/* ---- Rvd_accumulate: clip in parallel, accumulate in the serial order ----
 *
 * Every per-site sum is a floating-point accumulation whose result depends on
 * the order of its terms, and the dual triples feed an unstable qsort, so the
 * output is defined by the serial face order.  Faces are therefore processed in
 * blocks: the expensive exact clipping of a block's non-fast faces runs on all
 * threads, each face writing its contributions ({site, w, centroid}) and dual
 * triples, in their serial order, into its thread's buffers; then one thread
 * replays the block face by face, applying fast faces inline and clipped faces'
 * records, with the original expressions.  The per-site sums, the dual records
 * (including the dual_cap cut-off) and the log lines are those of the serial
 * loop for any thread count (2026-10-04). */

typedef struct { int32_t si; int32_t pad; double w, ct0, ct1, ct2; } RvdSub;
typedef struct { int32_t k[6]; } RvdDual;   /* {sorted key[3], oriented tri[3]} */

typedef struct {
    AttenePredCtx *ctx;      /* all nsites sites registered first: ids match */
    int32_t  *cand;          /* [RVD_MAX_CAND] */
    double   *cxyz;          /* [3 * RVD_MAX_CAND] the candidates' coordinates, by row */
    int32_t  *owner_gen;     /* [nsites] generation stamps */
    int32_t  *owner_q;       /* [RVD_MAX_CAND] */
    int32_t   cur_gen;
    RvdSub   *sub;  size_t nsub, cap_sub;
    RvdDual  *dual; size_t ndual, cap_dual;
    int       oom;           /* a buffer could not grow */
    RvdStats  st;
    char      pad[64];
} RvdWork;

typedef struct {
    int32_t  tid;            /* -2 fast face, -1 to clip at replay, else worker */
    uint32_t ncand;
    size_t   sub0, nsub, dual0, ndual;
} RvdFaceOut;

#define RVD_BLOCK_FACES   8192
#define RVD_PAR_MIN_FACES 2048

static int rvd_push_sub(RvdWork *w, int32_t si, double wt, double c0, double c1, double c2) {
    if (w->nsub == w->cap_sub) {
        size_t nc = w->cap_sub ? w->cap_sub * 2 : 4096;
        RvdSub *ns = (RvdSub *)realloc(w->sub, nc * sizeof(RvdSub));
        if (!ns) { w->oom = 1; return -1; }
        w->sub = ns; w->cap_sub = nc;
    }
    RvdSub *s = &w->sub[w->nsub++];
    s->si = si; s->pad = 0; s->w = wt; s->ct0 = c0; s->ct1 = c1; s->ct2 = c2;
    return 0;
}

static int rvd_push_dual(RvdWork *w, const int32_t key[3], const int32_t t[3]) {
    if (w->ndual == w->cap_dual) {
        size_t nc = w->cap_dual ? w->cap_dual * 2 : 2048;
        RvdDual *nd = (RvdDual *)realloc(w->dual, nc * sizeof(RvdDual));
        if (!nd) { w->oom = 1; return -1; }
        w->dual = nd; w->cap_dual = nc;
    }
    RvdDual *d = &w->dual[w->ndual++];
    d->k[0] = key[0]; d->k[1] = key[1]; d->k[2] = key[2];
    d->k[3] = t[0];   d->k[4] = t[1];   d->k[5] = t[2];
    return 0;
}

/* Clip one non-fast face into w's buffers (records in the serial order).
 * Returns 0, or -1 when a buffer could not grow (the face is then redone at
 * replay).  Thread-safe for distinct w: everything shared is read-only. */
static int rvd_clip_face(const struct Rvd_T *r, const double *sites, size_t nsites,
                         KDTree_T kt, const RvdSiteGrid *grid, const float *sf,
                         const int32_t *vnear, size_t f,
                         const CvtField *dens, int want_dual, RvdWork *w,
                         RvdFaceOut *fo) {
    int timed = (r->stats != NULL);
    uint64_t tk0 = timed ? cycle_now() : 0;
    const double *A=&r->verts[(size_t)r->faces[f*3+0]*3];
    const double *B=&r->verts[(size_t)r->faces[f*3+1]*3];
    const double *C=&r->verts[(size_t)r->faces[f*3+2]*3];
    double nrm[3]; { double ab[3],ac[3]; vsub3(B,A,ab); vsub3(C,A,ac); vcross3(ab,ac,nrm); }
    int32_t na = vnear[r->faces[f*3+0]], nb = vnear[r->faces[f*3+1]], nc2 = vnear[r->faces[f*3+2]];
    size_t sub0 = w->nsub, dual0 = w->ndual;
    PV poly[RVD_MAX_POLY];

    int32_t va=AttenePred_explicit(w->ctx, A[0],A[1],A[2]);
    int32_t vb=AttenePred_explicit(w->ctx, B[0],B[1],B[2]);
    int32_t vc=AttenePred_explicit(w->ctx, C[0],C[1],C[2]);

    const int32_t near3[3] = { na, nb, nc2 };
    int wide = 1;
    size_t ncand = triangle_candidates(kt, grid, sf, nsites, near3, A, B, C, w->cand, &wide);
    w->st.n_wide += (uint64_t)wide;
    for (size_t k = 0; k < ncand; k++) {
        const double *S = &sites[(size_t)w->cand[k]*3];
        w->cxyz[k] = S[0]; w->cxyz[RVD_MAX_CAND + k] = S[1]; w->cxyz[2*RVD_MAX_CAND + k] = S[2];
    }
    if (timed) {
        uint64_t tk = cycle_now(); w->st.tk_cand += tk - tk0; tk0 = tk;
        w->st.n_cand += ncand;
        if (ncand >= RVD_MAX_CAND) w->st.n_sat++;
    }
    fo->ncand = (uint32_t)ncand;

    /* Propagation: clip only ACTUAL owners. The owners of a convex triangle are
     * connected in the restricted-Voronoi adjacency graph, so a BFS seeded at the
     * three vertices' nearest sites reaches them all; each clipped cell's bisector
     * edges (etag==1) name the neighbouring site (bsj) to enqueue. `cand` stays the
     * membership superset clip_cell tests against. The queue is sized RVD_MAX_CAND
     * (>= ncand), so it never overflows -- worst case degrades to the old bound. */
    int ok = 0;
    w->cur_gen++;
    int qh = 0, qt = 0;
    owner_push(na,  w->cur_gen, w->owner_gen, w->owner_q, &qt);
    owner_push(nb,  w->cur_gen, w->owner_gen, w->owner_q, &qt);
    owner_push(nc2, w->cur_gen, w->owner_gen, w->owner_q, &qt);
    while (qh < qt && ok == 0) {
        int32_t si = w->owner_q[qh++];
        int np = clip_cell(w->ctx, sites, w->cand, w->cxyz, ncand, si, A, B, C, va, vb, vc, poly,
                           timed ? &w->st.n_planes : NULL, timed ? &w->st.n_preds : NULL);
        if (timed) w->st.n_cells++;
        if (np < 3) continue;
        for (int e = 0; e < np; e++)
            if (poly[e].etag == 1)
                owner_push(poly[e].bsj, w->cur_gen, w->owner_gen, w->owner_q, &qt);

        double P0[3];
        if (AttenePred_xyz(w->ctx, poly[0].v, &P0[0],&P0[1],&P0[2]) == 0) {
            for (int i = 1; i+1 < np; i++) {
                double Pi[3], Pj[3];
                if (AttenePred_xyz(w->ctx, poly[i].v,   &Pi[0],&Pi[1],&Pi[2]) != 0) continue;
                if (AttenePred_xyz(w->ctx, poly[i+1].v, &Pj[0],&Pj[1],&Pj[2]) != 0) continue;
                double at = tri_area(P0,Pi,Pj);
                double ct0=(P0[0]+Pi[0]+Pj[0])/3.0, ct1=(P0[1]+Pi[1]+Pj[1])/3.0, ct2=(P0[2]+Pi[2]+Pj[2])/3.0;
                double wt = dens ? at*Rvd_density_rho(dens,ct0,ct1,ct2) : at;  /* NULL => w==at (bit-identical) */
                if (rvd_push_sub(w, si, wt, ct0, ct1, ct2) != 0) { ok = -1; break; }
            }
        }
        if (ok != 0) break;

        /* Dual: each facet-bisector2 corner is a Voronoi vertex of 3 sites. */
        if (want_dual) {
            for (int i = 0; i < np; i++) {
                if (poly[i].fi < 0) continue;
                int32_t t[3] = { poly[i].fi, poly[i].fj, poly[i].fk };
                /* orient by local surface normal */
                const double *Pa=&sites[(size_t)t[0]*3], *Pb=&sites[(size_t)t[1]*3], *Pc=&sites[(size_t)t[2]*3];
                double u[3], v[3], cr[3];
                vsub3(Pb,Pa,u); vsub3(Pc,Pa,v); vcross3(u,v,cr);
                if (vdot3(cr,nrm) < 0.0) { int32_t tt=t[1]; t[1]=t[2]; t[2]=tt; }
                int32_t key[3] = { t[0],t[1],t[2] };
                for (int p=0;p<2;p++) for (int q=0;q<2-p;q++)
                    if (key[q]>key[q+1]) { int32_t s=key[q]; key[q]=key[q+1]; key[q+1]=s; }
                if (rvd_push_dual(w, key, t) != 0) { ok = -1; break; }
            }
        }
    }
    AttenePred_truncate(w->ctx, (int32_t)nsites);
    if (timed) w->st.tk_clip += cycle_now() - tk0;
    if (ok != 0) { w->nsub = sub0; w->ndual = dual0; return -1; }
    fo->sub0 = sub0; fo->nsub = w->nsub - sub0;
    fo->dual0 = dual0; fo->ndual = w->ndual - dual0;
    return 0;
}

int Rvd_accumulate(Rvd_T r, const double *sites, size_t nsites,
                   double *out_area, double *out_centroid,
                   double *out_NA, double *out_bNA,
                   int32_t **out_dual_tris, size_t *out_ndual,
                   const CvtField *dens) {
    if (!r || !sites || nsites == 0 || !out_area || !out_centroid) return -1;
    if (out_dual_tris) *out_dual_tris = NULL;
    if (out_ndual) *out_ndual = 0;
    int want_moments = (out_NA != NULL && out_bNA != NULL);
    int want_dual = (out_dual_tris != NULL);

    for (size_t i = 0; i < nsites; i++) {
        out_area[i] = 0.0;
        out_centroid[i*3+0] = out_centroid[i*3+1] = out_centroid[i*3+2] = 0.0;
    }
    if (want_moments) {
        memset(out_NA,  0, nsites*6*sizeof(double));
        memset(out_bNA, 0, nsites*3*sizeof(double));
    }

    /* Raw dual triples ({sortedkey[3], tri[3]}) persist past the scratch restore. */
    size_t dual_cap = out_dual_tris ? (nsites * 8 + 16) : 0;
    int32_t *raw = NULL;
    size_t   nraw = 0;
    if (out_dual_tris)
        raw = (int32_t *)ARENA_ALLOC(r->arena, (size_t)(dual_cap * 6 * sizeof(int32_t)));

    Arena_Mark mark = Arena_save(r->arena);

    float *sf = (float *)ARENA_ALLOC(r->arena, (size_t)(nsites*3*sizeof(float)));
    for (size_t i = 0; i < nsites; i++) {
        sf[i*3+0]=(float)sites[i*3+0]; sf[i*3+1]=(float)sites[i*3+1]; sf[i*3+2]=(float)sites[i*3+2];
    }
    KDTree_T kt = KDTree_new(r->arena, sf, nsites);
    RvdSiteGrid grid;
    if (rvd_grid_build(r->arena, sf, nsites, &grid) != 0) grid.ps = NULL;   /* no bound: wide radius */
    RvdStats *st = r->stats;
    if (st) st->n_passes++;

    int nthr = 1;
#ifdef _OPENMP
    nthr = omp_get_max_threads();
    if (nthr < 1) nthr = 1;
    if (nthr > 64) nthr = 64;
#endif
    int par = (nthr > 1 && r->nf >= (size_t)RVD_PAR_MIN_FACES && r->nf < (size_t)INT32_MAX);
    if (!par) nthr = 1;

    /* Nearest site of every mesh vertex: the same float query the per-face code
     * made, three times per face, now once per vertex (pure, so parallel). */
    uint64_t tkn = st ? cycle_now() : 0;
    int32_t *vnear = r->vnear;   /* in place: last pass's nearest site seeds this pass */
    {
        int nv_i = (int)r->nv, vi;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nthr) if (par)
#endif
        for (vi = 0; vi < nv_i; vi++) {
            const double *P = &r->verts[(size_t)vi*3];
            float pf[3] = { (float)P[0], (float)P[1], (float)P[2] };
            float d0;
            int32_t prev = vnear[vi];
            float bnd = FLT_MAX;
            if (prev >= 0 && (size_t)prev < nsites) {
                const float *s = &sf[(size_t)prev*3];
                float e0 = pf[0] - s[0], e1 = pf[1] - s[1], e2 = pf[2] - s[2];
                bnd = e0 * e0 + e1 * e1 + e2 * e2;
            }
            if (bnd < FLT_MAX)   /* exact: site `prev` lies below the bound */
                vnear[vi] = (int32_t)KDTree_nearest_bounded(kt, pf, nextafterf(bnd, FLT_MAX), &d0);
            else
                vnear[vi] = (int32_t)KDTree_nearest(kt, pf, &d0);
        }
    }
    if (st) st->tk_near += cycle_now() - tkn;

    RvdWork *work = (RvdWork *)ARENA_CALLOC(r->arena, (size_t)nthr, sizeof(RvdWork));
    for (int t = 0; t < nthr; t++) {
        RvdWork *w = &work[t];
        w->ctx = AttenePred_new();
        for (size_t i = 0; i < nsites; i++)
            (void)AttenePred_explicit(w->ctx, sites[i*3+0], sites[i*3+1], sites[i*3+2]);
        w->cand      = (int32_t *)ARENA_ALLOC(r->arena, (size_t)(RVD_MAX_CAND*sizeof(int32_t)));
        w->cxyz      = (double *)ARENA_ALLOC(r->arena, (size_t)(3*RVD_MAX_CAND*sizeof(double)));
        w->owner_gen = (int32_t *)ARENA_CALLOC(r->arena, (size_t)nsites, sizeof(int32_t));
        w->owner_q   = (int32_t *)ARENA_ALLOC(r->arena, (size_t)(RVD_MAX_CAND*sizeof(int32_t)));
    }
    RvdFaceOut *fout = (RvdFaceOut *)ARENA_ALLOC(r->arena, (size_t)(RVD_BLOCK_FACES*sizeof(RvdFaceOut)));
    int rc = 0;

    for (size_t f0 = 0; f0 < r->nf && rc == 0; f0 += RVD_BLOCK_FACES) {
        size_t f1 = f0 + RVD_BLOCK_FACES < r->nf ? f0 + RVD_BLOCK_FACES : r->nf;
        int nb_i = (int)(f1 - f0), k;
        for (int t = 0; t < nthr; t++) { work[t].nsub = 0; work[t].ndual = 0; }

        /* 1. exact clipping of the block's non-fast faces, on all threads */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 16) num_threads(nthr) if (par)
#endif
        for (k = 0; k < nb_i; k++) {
            size_t f = f0 + (size_t)k;
            RvdFaceOut *fo = &fout[k];
            int32_t na = vnear[r->faces[f*3+0]], nb = vnear[r->faces[f*3+1]], nc2 = vnear[r->faces[f*3+2]];
            if (na == nb && nb == nc2) { fo->tid = -2; continue; }
            int tid = 0;
#ifdef _OPENMP
            tid = omp_get_thread_num();
#endif
            fo->tid = tid;
            if (rvd_clip_face(r, sites, nsites, kt, &grid, sf, vnear, f, dens, want_dual, &work[tid], fo) != 0)
                fo->tid = -1;
        }

        /* 2. serial replay in face order */
        for (size_t f = f0; f < f1 && rc == 0; f++) {
            RvdFaceOut *fo = &fout[f - f0];
            const double *A=&r->verts[(size_t)r->faces[f*3+0]*3];
            const double *B=&r->verts[(size_t)r->faces[f*3+1]*3];
            const double *C=&r->verts[(size_t)r->faces[f*3+2]*3];

            double nrm[3]; { double ab[3],ac[3]; vsub3(B,A,ab); vsub3(C,A,ac); vcross3(ab,ac,nrm); }
            double un[3]={nrm[0],nrm[1],nrm[2]}, nnT[6]={0,0,0,0,0,0};
            if (want_moments) {
                double nl=sqrt(un[0]*un[0]+un[1]*un[1]+un[2]*un[2]);
                if (nl>0.0){ un[0]/=nl; un[1]/=nl; un[2]/=nl; }
                nnT[0]=un[0]*un[0]; nnT[1]=un[0]*un[1]; nnT[2]=un[0]*un[2];
                nnT[3]=un[1]*un[1]; nnT[4]=un[1]*un[2]; nnT[5]=un[2]*un[2];
            }
            if (st) st->n_faces++;

            if (fo->tid == -2) {
                /* Fast path: the three vertices share their nearest site, so the
                 * whole convex planar triangle lies in that Voronoi cell. */
                if (st) st->n_fast++;
                int32_t si=vnear[r->faces[f*3+0]];
                double at=tri_area(A,B,C);
                double ct0=(A[0]+B[0]+C[0])/3.0, ct1=(A[1]+B[1]+C[1])/3.0, ct2=(A[2]+B[2]+C[2])/3.0;
                double w = dens ? at*Rvd_density_rho(dens,ct0,ct1,ct2) : at;  /* NULL => w==at (bit-identical) */
                out_area[si]+=w;
                out_centroid[si*3+0]+=w*ct0; out_centroid[si*3+1]+=w*ct1; out_centroid[si*3+2]+=w*ct2;
                if (want_moments) {
                    double *Ni=&out_NA[(size_t)si*6];
                    Ni[0]+=w*nnT[0]; Ni[1]+=w*nnT[1]; Ni[2]+=w*nnT[2];
                    Ni[3]+=w*nnT[3]; Ni[4]+=w*nnT[4]; Ni[5]+=w*nnT[5];
                    double ndotc=un[0]*ct0+un[1]*ct1+un[2]*ct2;
                    double *bi=&out_bNA[(size_t)si*3];
                    bi[0]+=w*ndotc*un[0]; bi[1]+=w*ndotc*un[1]; bi[2]+=w*ndotc*un[2];
                }
                continue;
            }

            RvdWork *w = NULL;
            if (fo->tid >= 0) {
                w = &work[fo->tid];
            } else {
                /* a buffer could not grow in the parallel pass: redo it here */
                w = &work[0];
                w->oom = 0;
                if (rvd_clip_face(r, sites, nsites, kt, &grid, sf, vnear, f, dens, want_dual, w, fo) != 0) {
                    rc = -6;
                    break;
                }
            }
            if (fo->ncand >= RVD_MAX_CAND)
                fprintf(stderr, "  rvd: face %zu saturated candidates (%zu)\n", f, (size_t)fo->ncand);

            for (size_t j = 0; j < fo->nsub; j++) {
                const RvdSub *s = &w->sub[fo->sub0 + j];
                int32_t si = s->si;
                double wt = s->w, ct0 = s->ct0, ct1 = s->ct1, ct2 = s->ct2;
                out_area[si] += wt;
                out_centroid[si*3+0] += wt*ct0;
                out_centroid[si*3+1] += wt*ct1;
                out_centroid[si*3+2] += wt*ct2;
                if (want_moments) {
                    double *Ni=&out_NA[(size_t)si*6];
                    Ni[0]+=wt*nnT[0]; Ni[1]+=wt*nnT[1]; Ni[2]+=wt*nnT[2];
                    Ni[3]+=wt*nnT[3]; Ni[4]+=wt*nnT[4]; Ni[5]+=wt*nnT[5];
                    double ndotc=un[0]*ct0+un[1]*ct1+un[2]*ct2;
                    double *bi=&out_bNA[(size_t)si*3];
                    bi[0]+=wt*ndotc*un[0]; bi[1]+=wt*ndotc*un[1]; bi[2]+=wt*ndotc*un[2];
                }
            }
            for (size_t j = 0; j < fo->ndual; j++) {
                if (nraw < dual_cap) {
                    const int32_t *d = w->dual[fo->dual0 + j].k;
                    int32_t *rec = &raw[nraw*6];
                    rec[0]=d[0]; rec[1]=d[1]; rec[2]=d[2];
                    rec[3]=d[3]; rec[4]=d[4]; rec[5]=d[5];
                    nraw++;
                }
            }
        }
    }

    for (int t = 0; t < nthr; t++) {
        RvdWork *w = &work[t];
        if (st) {
            st->tk_cand += w->st.tk_cand; st->tk_clip += w->st.tk_clip;
            st->n_cand += w->st.n_cand;   st->n_sat += w->st.n_sat;
            st->n_cells += w->st.n_cells; st->n_planes += w->st.n_planes;
            st->n_preds += w->st.n_preds; st->n_wide += w->st.n_wide;
        }
        AttenePred_free(w->ctx);
        free(w->sub);
        free(w->dual);
    }
    Arena_restore(r->arena, mark);
    if (rc != 0) return rc;

    for (size_t i = 0; i < nsites; i++) {
        if (out_area[i] > 0.0) {
            out_centroid[i*3+0] /= out_area[i];
            out_centroid[i*3+1] /= out_area[i];
            out_centroid[i*3+2] /= out_area[i];
        } else {
            out_centroid[i*3+0]=sites[i*3+0]; out_centroid[i*3+1]=sites[i*3+1]; out_centroid[i*3+2]=sites[i*3+2];
        }
    }

    /* Dedup the dual triples (each Voronoi vertex is seen once per incident cell). */
    if (out_dual_tris && nraw > 0) {
        qsort(raw, nraw, 6*sizeof(int32_t), cmp_triple);
        int32_t *tris = (int32_t *)ARENA_ALLOC(r->arena, (size_t)(nraw*3*sizeof(int32_t)));
        size_t nt = 0;
        for (size_t i = 0; i < nraw; i++) {
            if (i > 0 && cmp_triple(&raw[i*6], &raw[(i-1)*6]) == 0) continue;
            /* skip degenerate (repeated site index) */
            int32_t a=raw[i*6+3], b=raw[i*6+4], c=raw[i*6+5];
            if (a==b || b==c || a==c) continue;
            tris[nt*3+0]=a; tris[nt*3+1]=b; tris[nt*3+2]=c; nt++;
        }
        *out_dual_tris = tris;
        *out_ndual = nt;
    }

    return 0;
}

/* ---- selftest: one restricted-Voronoi pass at 1, 2, 3 and 8 threads must give
 * bitwise-identical cell areas, centroids, CWF moments and dual triangles, and
 * the cells must partition the surface. ---- */
static uint32_t rvd_st_x = 0x1234567u;
static double rvd_st_rnd(void) {
    uint32_t x = rvd_st_x; x ^= x << 13; x ^= x >> 17; x ^= x << 5; rvd_st_x = x;
    return (double)(x >> 8) * (1.0 / 16777216.0);
}

int Rvd_selftest(void) {
    Arena_T a = Arena_new();
    int ok = 1;
    const int M = 120;                       /* (M+1)^2 verts, 2 M^2 faces */
    size_t nv = (size_t)(M + 1) * (size_t)(M + 1), nf = (size_t)M * (size_t)M * 2;
    double  *V = (double *)ARENA_ALLOC(a, (size_t)(nv * 3 * sizeof(double)));
    int32_t *F = (int32_t *)ARENA_ALLOC(a, (size_t)(nf * 3 * sizeof(int32_t)));
    for (int i = 0; i <= M; i++)
        for (int j = 0; j <= M; j++) {
            size_t v = (size_t)i * (size_t)(M + 1) + (size_t)j;
            double u = 0.5 * i, w = 0.5 * j;
            V[v*3+0] = 1.0e4 + 1.5 * sin(0.21 * u) * cos(0.17 * w);   /* wavy sheet */
            V[v*3+1] = 2.5e3 + u;
            V[v*3+2] = 3.1e3 + w;
        }
    size_t fi = 0;
    for (int i = 0; i < M; i++)
        for (int j = 0; j < M; j++) {
            int32_t v00 = (int32_t)(i * (M + 1) + j), v01 = v00 + 1;
            int32_t v10 = v00 + (M + 1), v11 = v10 + 1;
            F[fi*3+0] = v00; F[fi*3+1] = v10; F[fi*3+2] = v11; fi++;
            F[fi*3+0] = v00; F[fi*3+1] = v11; F[fi*3+2] = v01; fi++;
        }
    Rvd_T r = Rvd_new(a, V, nv, F, nf);
    const size_t N = 1500;
    double *sites = (double *)ARENA_ALLOC(a, (size_t)(N * 3 * sizeof(double)));
    rvd_st_x = 0x1234567u;
    for (size_t k = 0; k < N; k++) {
        size_t f = (size_t)(rvd_st_rnd() * (double)nf); if (f >= nf) f = nf - 1;
        double u = rvd_st_rnd(), w = rvd_st_rnd();
        if (u + w > 1.0) { u = 1.0 - u; w = 1.0 - w; }
        const double *A = &V[(size_t)F[f*3+0]*3], *B = &V[(size_t)F[f*3+1]*3], *C = &V[(size_t)F[f*3+2]*3];
        for (int c = 0; c < 3; c++) sites[k*3+c] = A[c] + u * (B[c] - A[c]) + w * (C[c] - A[c]);
    }
    double *ar[2], *ce[2], *na[2], *nb[2];
    int32_t *tr[2]; size_t nt[2];
    for (int s2 = 0; s2 < 2; s2++) {
        ar[s2] = (double *)ARENA_ALLOC(a, (size_t)(N * sizeof(double)));
        ce[s2] = (double *)ARENA_ALLOC(a, (size_t)(N * 3 * sizeof(double)));
        na[s2] = (double *)ARENA_ALLOC(a, (size_t)(N * 6 * sizeof(double)));
        nb[s2] = (double *)ARENA_ALLOC(a, (size_t)(N * 3 * sizeof(double)));
    }
#ifdef _OPENMP
    int saved = omp_get_max_threads();
#endif
    /* reference: 1 thread, whole-triangle filter OFF; then filter ON */
    static const int thr[5] = { 1, 1, 2, 3, 8 };
    RvdStats fst; memset(&fst, 0, sizeof(fst));
    uint64_t planes_ref = 0;
    for (int ti = 0; ti < 5 && ok; ti++) {
#ifdef _OPENMP
        omp_set_num_threads(thr[ti]);
#endif
        g_rvd_tri_filter = (ti == 0) ? 0 : 1;
        if (ti <= 1) { memset(&fst, 0, sizeof(fst)); Rvd_set_stats(r, &fst); }
        int slot = ti == 0 ? 0 : 1;
        if (Rvd_accumulate(r, sites, N, ar[slot], ce[slot], na[slot], nb[slot], NULL, NULL, NULL) != 0 ||
            Rvd_accumulate(r, sites, N, ar[slot], ce[slot], NULL, NULL, &tr[slot], &nt[slot], NULL) != 0) {
            fprintf(stderr, "Rvd_selftest FAIL: accumulate error at %d threads\n", thr[ti]);
            ok = 0; break;
        }
        /* moments come from the first call; redo it so ar/ce match that call too */
        if (Rvd_accumulate(r, sites, N, ar[slot], ce[slot], na[slot], nb[slot], NULL, NULL, NULL) != 0) {
            ok = 0; break;
        }
        if (ti == 0) planes_ref = fst.n_planes;
        if (ti <= 1) Rvd_set_stats(r, NULL);
        if (ti == 0) {
            double sum = 0.0;
            for (size_t k = 0; k < N; k++) sum += ar[0][k];
            double ma = Rvd_mesh_area(r);
            if (!(fabs(sum - ma) <= 1e-9 * ma) || nt[0] == 0) {
                fprintf(stderr, "Rvd_selftest FAIL: partition %.12g vs mesh %.12g, %zu dual\n", sum, ma, nt[0]);
                ok = 0;
            }
            continue;
        }
        if (memcmp(ar[0], ar[1], N * sizeof(double)) || memcmp(ce[0], ce[1], N * 3 * sizeof(double)) ||
            memcmp(na[0], na[1], N * 6 * sizeof(double)) || memcmp(nb[0], nb[1], N * 3 * sizeof(double)) ||
            nt[0] != nt[1] || memcmp(tr[0], tr[1], nt[0] * 3 * sizeof(int32_t))) {
            fprintf(stderr, "Rvd_selftest FAIL: filtered at %d threads differs from the unfiltered reference\n",
                    thr[ti]);
            ok = 0;
        }
    }
    /* Tight vs wide candidate radius: bitwise the same on this sheet (the tight
     * radius engaged) and on a coarse mesh of the same sheet whose wide lists hold
     * all 1500 sites (the guard must keep them wide). */
    uint64_t wide_fine = 0, wide_coarse = 0, nonfast_coarse = 0;
    for (int coarse = 0; coarse < 2 && ok; coarse++) {
        Rvd_T rr = r;
        if (coarse) {
            const int M2 = 8;
            size_t nv2 = (size_t)(M2 + 1) * (size_t)(M2 + 1), nf2 = (size_t)M2 * (size_t)M2 * 2;
            double *V2 = (double *)ARENA_ALLOC(a, (size_t)(nv2 * 3 * sizeof(double)));
            int32_t *F2 = (int32_t *)ARENA_ALLOC(a, (size_t)(nf2 * 3 * sizeof(int32_t)));
            for (int i = 0; i <= M2; i++)
                for (int j = 0; j <= M2; j++) {
                    size_t v = (size_t)i * (size_t)(M2 + 1) + (size_t)j;
                    double u = 60.0 * i / M2, w = 60.0 * j / M2;
                    V2[v*3+0] = 1.0e4 + 1.5 * sin(0.21 * u) * cos(0.17 * w);
                    V2[v*3+1] = 2.5e3 + u;
                    V2[v*3+2] = 3.1e3 + w;
                }
            size_t f2 = 0;
            for (int i = 0; i < M2; i++)
                for (int j = 0; j < M2; j++) {
                    int32_t v00 = (int32_t)(i * (M2 + 1) + j), v01 = v00 + 1;
                    int32_t v10 = v00 + (M2 + 1), v11 = v10 + 1;
                    F2[f2*3+0] = v00; F2[f2*3+1] = v10; F2[f2*3+2] = v11; f2++;
                    F2[f2*3+0] = v00; F2[f2*3+1] = v11; F2[f2*3+2] = v01; f2++;
                }
            rr = Rvd_new(a, V2, nv2, F2, nf2);
        }
        for (int tight = 0; tight < 2 && ok; tight++) {
            g_rvd_tight = tight;
            RvdStats ts; memset(&ts, 0, sizeof(ts));
            Rvd_set_stats(rr, &ts);
            int slot = tight;
            if (Rvd_accumulate(rr, sites, N, ar[slot], ce[slot], na[slot], nb[slot], NULL, NULL, NULL) != 0 ||
                Rvd_accumulate(rr, sites, N, ar[slot], ce[slot], NULL, NULL, &tr[slot], &nt[slot], NULL) != 0 ||
                Rvd_accumulate(rr, sites, N, ar[slot], ce[slot], na[slot], nb[slot], NULL, NULL, NULL) != 0) {
                fprintf(stderr, "Rvd_selftest FAIL: accumulate error (tight %d, coarse %d)\n", tight, coarse);
                ok = 0;
            }
            Rvd_set_stats(rr, NULL);
            if (tight && coarse) { wide_coarse = ts.n_wide; nonfast_coarse = ts.n_faces - ts.n_fast; }
            if (tight && !coarse) wide_fine = ts.n_wide;
        }
        if (ok && (memcmp(ar[0], ar[1], N * sizeof(double)) || memcmp(ce[0], ce[1], N * 3 * sizeof(double)) ||
                   memcmp(na[0], na[1], N * 6 * sizeof(double)) || memcmp(nb[0], nb[1], N * 3 * sizeof(double)) ||
                   nt[0] != nt[1] || memcmp(tr[0], tr[1], nt[0] * 3 * sizeof(int32_t)))) {
            fprintf(stderr, "Rvd_selftest FAIL: tight candidate radius differs from the wide one (coarse %d)\n", coarse);
            ok = 0;
        }
    }
    g_rvd_tight = 1;
    if (ok && (wide_fine != 0 || wide_coarse == 0 || wide_coarse != nonfast_coarse)) {
        fprintf(stderr, "Rvd_selftest FAIL: tight-radius guard (fine wide %llu, coarse wide %llu of %llu)\n",
                (unsigned long long)wide_fine, (unsigned long long)wide_coarse, (unsigned long long)nonfast_coarse);
        ok = 0;
    }
#ifdef _OPENMP
    omp_set_num_threads(saved);
#endif
    g_rvd_tri_filter = 1;
    Arena_dispose(&a);
    fprintf(stderr, "Rvd_selftest: %s (%zu faces, %zu sites; triangle filter at 1/2/3/8 threads bitwise "
            "== unfiltered; planes tested %llu -> %llu; tight candidate radius bitwise == wide, "
            "guard kept %llu of %llu coarse faces wide)\n", ok ? "OK" : "FAIL", nf, N,
            (unsigned long long)planes_ref, (unsigned long long)fst.n_planes,
            (unsigned long long)wide_coarse, (unsigned long long)nonfast_coarse);
    return ok ? 0 : -1;
}
