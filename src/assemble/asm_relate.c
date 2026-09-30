/* asm_relate.c -- seam relations and layer neighbours between charts. */
#include "asm_relate.h"
#include "asm_continuity.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/ves_omp.h"

#include "../common/closest_tri.h"
#include "../common/eig3.h"
#include "../common/kdtree.h"
#include "../common/pipeline_constants.h"
#include "../common/union_find.h"
#include "../common/ves_platform.h"

#define AR_PROBE_INITIAL_CAP 2048
#define AR_PROBE_MIN_RADIUS 6.0
#define AR_PI 3.14159265358979323846

void AsmRelate_default_opts(AsmRelateOpts *o)
{
    o->cube_size = 128.0;
    o->seam_max_vox = ASM_SEAM_MAX_VOX;
    o->seam_band_vox = 8.0;
    o->seam_min_len = ASM_SEAM_MIN_LEN;
    o->normal_dot_min = ASM_SEAM_NORMAL_DOT;
    o->tangency_max_vox = ASM_TANGENCY_MAX_VOX;
    o->iso_tol = ASM_ISO_TOL;
    o->seam_noise_vox = ASM_SEAM_NOISE_VOX;
    o->layer_probe_vox = ASM_LAYER_PROBE_VOX;
    o->layer_samples = 48;
    o->diag_csv = NULL;
}

/* ---- rigid 2-D fit ---------------------------------------------------------- */

int AsmRelate_rigid_fit(const double *qa, const double *qb, const double *w, size_t n,
                        double *theta, double *tx, double *ty, double *rms)
{
    if (n < 2) return -1;
    double sw = 0.0, ca0 = 0.0, ca1 = 0.0, cb0 = 0.0, cb1 = 0.0;
    for (size_t i = 0; i < n; i++) {
        double wi = w ? w[i] : 1.0;
        sw += wi;
        ca0 += wi * qa[i*2]; ca1 += wi * qa[i*2+1];
        cb0 += wi * qb[i*2]; cb1 += wi * qb[i*2+1];
    }
    if (!(sw > 0.0)) return -1;
    ca0 /= sw; ca1 /= sw; cb0 /= sw; cb1 /= sw;
    double sxx = 0.0, sxy = 0.0;  /* sum of cross terms for the angle */
    for (size_t i = 0; i < n; i++) {
        double wi = w ? w[i] : 1.0;
        double ax = qa[i*2] - ca0, ay = qa[i*2+1] - ca1;
        double bx = qb[i*2] - cb0, by = qb[i*2+1] - cb1;
        sxx += wi * (bx * ax + by * ay);
        sxy += wi * (bx * ay - by * ax);
    }
    double th = atan2(sxy, sxx);
    double c = cos(th), s = sin(th);
    *theta = th;
    *tx = ca0 - (c * cb0 - s * cb1);
    *ty = ca1 - (s * cb0 + c * cb1);
    double e2 = 0.0;
    for (size_t i = 0; i < n; i++) {
        double wi = w ? w[i] : 1.0;
        double px = c * qb[i*2] - s * qb[i*2+1] + *tx;
        double py = s * qb[i*2] + c * qb[i*2+1] + *ty;
        double dx = px - qa[i*2], dy = py - qa[i*2+1];
        e2 += wi * (dx*dx + dy*dy);
    }
    *rms = sqrt(e2 / sw);
    return 0;
}

/* ---- thread-local growable buffers ------------------------------------------ */

/* SHORT SEAMS (2026-09-10): a seam under the pose graph's minimum length is no pose-graph edge (its
 * rotation is ill-determined), but when its correspondences fit rigidly it is exact PLACEMENT evidence
 * for the cube-sized pieces of the outer wraps, kept like a WEAK relation.  0 = the pre-2026-09-10
 * behaviour (correspondences discarded). */
#define AR_SHORT_SEAMS 1
#define AR_SHORT_MIN_LEN 12.0     /* vox: below this a seam is a point contact */
#define AR_SHORT_MIN_CORR 6       /* a rigid fit on fewer correspondences is unconstrained (640 of the 10x's 1,403 length rejections had n < 6) */
#define AR_RIGID_MAX_RMS 4.0      /* vox: the rigid-fit limit for accepted, WEAK and SHORT relations alike */

typedef struct ArCand {
    int32_t chart_a, chart_b;   /* chart ids, a < b */
    int32_t va, vb;             /* chart-local vertex indices */
    float   dist;               /* 3-D pair distance */
    float   ngap;               /* |displacement . n_a| */
    int8_t  sign;               /* sign of n_a . n_b */
} ArCand;

typedef struct ArTrim {
    float gap_a[2], gap_b[2];
    uint8_t valid;
} ArTrim;

static void ar_append_relation(AsmRun *run, const AsmRelation *relation,
                                const ArCand *candidate, const ArTrim *trim, size_t n)
{
    AsmRelation r=*relation;
    r.corr_first=(int32_t)run->n_corr; r.corr_count=(int32_t)n;
    for (size_t k=0; k<n; k++) {
        AsmCorr c={candidate[k].va,candidate[k].vb};
        /* Stage 3 needs the same midpoint that supplied the fit's rotation
         * information. These measurements do not certify a boundary run;
         * AsmContinuity_prepare still rebuilds source bases and ordering. */
        c.run=-1;
        if (trim) {
            memcpy(c.gap_a,trim[k].gap_a,sizeof c.gap_a);
            memcpy(c.gap_b,trim[k].gap_b,sizeof c.gap_b);
            c.valid=trim[k].valid;
        }
        AsmRun_push_corr(run,&c);
    }
    AsmRun_push_rel(run,&r);
}

/* Measure once per chart, after gathering candidates. A chart can meet many
 * neighbouring cubes; rebuilding its tangent frames per pair wastes work.
 * This neither filters observations nor consults any assembled chart pose. */
static int ar_candidate_trims(const AsmRun *run, const ArCand *cand, size_t n,
                              Arena_T arena, ArTrim *trim)
{
    Arena_Mark mark = Arena_save(arena);
    size_t nc = run->n_charts;
    size_t *off = ARENA_CALLOC(arena,nc+1,sizeof *off);
    memset(trim,0,n*sizeof *trim);
    for (size_t k = 0; k < n; k++) {
        const ArCand *g = &cand[k];
        if (g->chart_a < 0 || g->chart_b < 0 || (size_t)g->chart_a >= nc || (size_t)g->chart_b >= nc ||
            g->va < 0 || g->vb < 0 || (size_t)g->va >= run->charts[g->chart_a].nv ||
            (size_t)g->vb >= run->charts[g->chart_b].nv) { Arena_restore(arena,mark); return -1; }
        off[g->chart_a+1]++; off[g->chart_b+1]++;
    }
    for (size_t i = 0; i < nc; i++) off[i+1] += off[i];
    size_t *adj = ARENA_ALLOC(arena,(off[nc] ? off[nc] : 1)*sizeof *adj);
    size_t *fill = ARENA_CALLOC(arena,nc ? nc : 1,sizeof *fill);
    for (size_t k = 0; k < n; k++) {
        adj[off[cand[k].chart_a]+fill[cand[k].chart_a]++] = k;
        adj[off[cand[k].chart_b]+fill[cand[k].chart_b]++] = k;
    }
    for (size_t i = 0; i < nc; i++) {
        size_t count = off[i+1]-off[i];
        if (!count) continue;
        Arena_Mark chart_mark = Arena_save(arena);
        int32_t *vertices = ARENA_ALLOC(arena,count*sizeof *vertices);
        double *delta = ARENA_ALLOC(arena,count*3*sizeof *delta);
        float *gaps = ARENA_ALLOC(arena,count*2*sizeof *gaps);
        uint8_t *valid = ARENA_ALLOC(arena,count);
        for (size_t j = 0; j < count; j++) {
            const ArCand *g = &cand[adj[off[i]+j]];
            vertices[j] = (size_t)g->chart_a == i ? g->va : g->vb;
            const float *a = run->charts[g->chart_a].xyz+3*(size_t)g->va;
            const float *b = run->charts[g->chart_b].xyz+3*(size_t)g->vb;
            for (int q = 0; q < 3; q++) delta[3*j+q] = (double)b[q]-a[q];
        }
        if (AsmContinuity_chart_gaps(arena,&run->charts[i],vertices,delta,count,gaps,valid)) {
            Arena_restore(arena,mark); return -1;
        }
        for (size_t j = 0; j < count; j++) {
            size_t k = adj[off[i]+j];
            int is_a = (size_t)cand[k].chart_a == i;
            memcpy(is_a ? trim[k].gap_a : trim[k].gap_b,gaps+2*j,2*sizeof(float));
            if (valid[j]) trim[k].valid |= is_a ? 1u : 2u;
        }
        Arena_restore(arena,chart_mark);
    }
    Arena_restore(arena,mark);
    return 0;
}

typedef struct ArHit {
    int32_t chart, other;
    int32_t va, fb;
    float   l0, l1;
    int8_t  side;
    int8_t  order;   /* crossing order along the probe: 1 = the first surface met */
    int8_t  nsign;   /* sign of n_a . n_face_b at the hit */
    float   dist;
} ArHit;

typedef struct ArBuf {
    Arena_T arena;
    ArCand *cand; size_t n_cand, cap_cand;
    ArHit  *hit;  size_t n_hit,  cap_hit;
    size_t gap_hist[ASM_RELATE_HIST];
    size_t gap_pairs;
} ArBuf;

static void ar_push_cand(ArBuf *b, const ArCand *c)
{
    if (b->n_cand == b->cap_cand) {
        size_t ncap = b->cap_cand ? b->cap_cand * 2 : 4096;
        ArCand *na = ARENA_ALLOC(b->arena, ncap * sizeof(ArCand));
        if (b->n_cand) memcpy(na, b->cand, b->n_cand * sizeof(ArCand));
        b->cand = na; b->cap_cand = ncap;
    }
    b->cand[b->n_cand++] = *c;
}

static void ar_push_hit(ArBuf *b, const ArHit *h)
{
    if (b->n_hit == b->cap_hit) {
        size_t ncap = b->cap_hit ? b->cap_hit * 2 : 4096;
        ArHit *na = ARENA_ALLOC(b->arena, ncap * sizeof(ArHit));
        if (b->n_hit) memcpy(na, b->hit, b->n_hit * sizeof(ArHit));
        b->hit = na; b->cap_hit = ncap;
    }
    b->hit[b->n_hit++] = *h;
}

/* ---- cube neighbourhood ------------------------------------------------------ */

typedef struct ArCubeIndex {
    size_t   n_cubes;
    int32_t *chart_first;  /* [n_cubes+1] offsets into chart_ids (in-layout charts per cube) */
    int32_t *chart_ids;
    int32_t *nbr;          /* [n_cubes*26] neighbour cube index or -1 */
} ArCubeIndex;

static int ar_find_cube(const AsmRun *run, long oz, long oy, long ox)
{
    /* pile is sorted by (oz, oy, ox): binary search */
    size_t lo = 0, hi = run->n_cubes;
    while (lo < hi) {
        size_t mid = (lo + hi) / 2;
        const MeshPileEntry *e = &run->pile[mid];
        int c = (e->oz < oz) ? -1 : (e->oz > oz) ? 1 : (e->oy < oy) ? -1 : (e->oy > oy) ? 1 : (e->ox < ox) ? -1 : (e->ox > ox) ? 1 : 0;
        if (c == 0) return (int)mid;
        if (c < 0) lo = mid + 1; else hi = mid;
    }
    return -1;
}

static void ar_build_cube_index(AsmRun *run, double cube_size, ArCubeIndex *ix)
{
    size_t nc = run->n_cubes;
    ix->n_cubes = nc;
    ix->chart_first = ARENA_CALLOC(run->arena, nc + 1, sizeof(int32_t));
    for (size_t i = 0; i < run->n_charts; i++)
        if (AsmChart_in_layout(&run->charts[i])) ix->chart_first[(size_t)run->charts[i].cube + 1]++;
    for (size_t c = 0; c < nc; c++) ix->chart_first[c+1] += ix->chart_first[c];
    ix->chart_ids = ARENA_ALLOC(run->arena, ((size_t)ix->chart_first[nc] + 1) * sizeof(int32_t));
    int32_t *fill = ARENA_CALLOC(run->arena, nc, sizeof(int32_t));
    for (size_t i = 0; i < run->n_charts; i++) {
        if (!AsmChart_in_layout(&run->charts[i])) continue;
        size_t c = (size_t)run->charts[i].cube;
        ix->chart_ids[ix->chart_first[c] + fill[c]++] = (int32_t)i;
    }
    ix->nbr = ARENA_ALLOC(run->arena, nc * 26 * sizeof(int32_t));
    long step = (long)cube_size;
    for (size_t c = 0; c < nc; c++) {
        int k = 0;
        for (int dz = -1; dz <= 1; dz++) for (int dy = -1; dy <= 1; dy++) for (int dx = -1; dx <= 1; dx++) {
            if (!dz && !dy && !dx) continue;
            ix->nbr[c*26 + (size_t)k++] = ar_find_cube(run, run->pile[c].oz + dz*step, run->pile[c].oy + dy*step, run->pile[c].ox + dx*step);
        }
    }
}

/* ---- seam candidates for one cube pair ---------------------------------------- */

typedef struct ArPt { float p[3]; float n[3]; int32_t chart; int32_t v; } ArPt;

static size_t ar_collect_boundary(const AsmRun *run, const ArCubeIndex *ix, size_t cube,
                                  const float lo[3], const float hi[3], double band,
                                  ArPt *out, size_t cap)
{
    size_t n = 0;
    for (int32_t k = ix->chart_first[cube]; k < ix->chart_first[cube+1]; k++) {
        const AsmChart *c = &run->charts[(size_t)ix->chart_ids[k]];
        for (size_t v = 0; v < c->nv; v++) {
            if (!c->boundary[v]) continue;
            const float *p = &c->xyz[v*3];
            if (p[0] < lo[0] - band || p[0] > hi[0] + band ||
                p[1] < lo[1] - band || p[1] > hi[1] + band ||
                p[2] < lo[2] - band || p[2] > hi[2] + band) continue;
            if (n < cap) {
                memcpy(out[n].p, p, 3 * sizeof(float));
                memcpy(out[n].n, &c->nrm[v*3], 3 * sizeof(float));
                out[n].chart = c->id; out[n].v = (int32_t)v;
            }
            n++;
        }
    }
    return n;
}

static size_t ar_count_boundary(const AsmRun *run, const ArCubeIndex *ix, size_t cube)
{
    size_t n = 0;
    for (int32_t k = ix->chart_first[cube]; k < ix->chart_first[cube+1]; k++) {
        const AsmChart *c = &run->charts[(size_t)ix->chart_ids[k]];
        for (size_t v = 0; v < c->nv; v++) n += c->boundary[v];
    }
    return n;
}

/* nearest point of set S to q: brute force under 64 points, else the tree */
static size_t ar_nearest(const ArPt *S, size_t n, KDTree_T tree, const float q[3], float *d2)
{
    if (tree) return KDTree_nearest(tree, q, d2);
    size_t best = 0; float bd = 1e30f;
    for (size_t i = 0; i < n; i++) {
        float dz = S[i].p[0]-q[0], dy = S[i].p[1]-q[1], dx = S[i].p[2]-q[2];
        float d = dz*dz + dy*dy + dx*dx;
        if (d < bd) { bd = d; best = i; }
    }
    *d2 = bd;
    return best;
}

static void ar_cube_box(const AsmRun *run, size_t cube, double cube_size, float lo[3], float hi[3])
{
    lo[0] = (float)run->pile[cube].oz; lo[1] = (float)run->pile[cube].oy; lo[2] = (float)run->pile[cube].ox;
    hi[0] = lo[0] + (float)cube_size; hi[1] = lo[1] + (float)cube_size; hi[2] = lo[2] + (float)cube_size;
}

static void ar_seam_pair(const AsmRun *run, const ArCubeIndex *ix, size_t ca, size_t cb,
                         const AsmRelateOpts *o, Arena_T scratch, ArBuf *buf)
{
    Arena_Mark mark = Arena_save(scratch);
    float loa[3], hia[3], lob[3], hib[3];
    ar_cube_box(run, ca, o->cube_size, loa, hia);
    ar_cube_box(run, cb, o->cube_size, lob, hib);
    size_t na_all = ar_count_boundary(run, ix, ca), nb_all = ar_count_boundary(run, ix, cb);
    if (na_all == 0 || nb_all == 0) { Arena_restore(scratch, mark); return; }
    ArPt *A = ARENA_ALLOC(scratch, na_all * sizeof(ArPt));
    ArPt *B = ARENA_ALLOC(scratch, nb_all * sizeof(ArPt));
    size_t na = ar_collect_boundary(run, ix, ca, lob, hib, o->seam_band_vox, A, na_all);
    size_t nb = ar_collect_boundary(run, ix, cb, loa, hia, o->seam_band_vox, B, nb_all);
    if (na < 3 || nb < 3) { Arena_restore(scratch, mark); return; }
    float *pa = ARENA_ALLOC(scratch, na * 3 * sizeof(float));
    float *pb = ARENA_ALLOC(scratch, nb * 3 * sizeof(float));
    for (size_t i = 0; i < na; i++) memcpy(&pa[i*3], A[i].p, 3 * sizeof(float));
    for (size_t i = 0; i < nb; i++) memcpy(&pb[i*3], B[i].p, 3 * sizeof(float));
    KDTree_T ta = na >= 64 ? KDTree_new(scratch, pa, na) : NULL;
    KDTree_T tb = nb >= 64 ? KDTree_new(scratch, pb, nb) : NULL;
    float max2 = (float)(o->seam_max_vox * o->seam_max_vox);
    float hist2 = (float)(ASM_RELATE_HIST * ASM_RELATE_HIST);
    for (size_t i = 0; i < na; i++) {
        float d2;
        size_t j = ar_nearest(B, nb, tb, A[i].p, &d2);
        if (d2 > hist2) continue;
        float d2b;
        size_t ii = ar_nearest(A, na, ta, B[j].p, &d2b);
        if (ii != i) continue;
        double nd = A[i].n[0]*B[j].n[0] + A[i].n[1]*B[j].n[1] + A[i].n[2]*B[j].n[2];
        if (fabs(nd) >= o->normal_dot_min) {
            int bin = (int)sqrtf(d2);
            if (bin >= ASM_RELATE_HIST) bin = ASM_RELATE_HIST - 1;
            buf->gap_hist[bin]++;
            buf->gap_pairs++;
        }
        if (d2 > max2) continue;
        if (fabs(nd) < o->normal_dot_min) continue;
        double dz = B[j].p[0]-A[i].p[0], dy = B[j].p[1]-A[i].p[1], dx = B[j].p[2]-A[i].p[2];
        double ng = fabs(dz*A[i].n[0] + dy*A[i].n[1] + dx*A[i].n[2]);
        ArCand c;
        if (A[i].chart < B[j].chart) { c.chart_a = A[i].chart; c.chart_b = B[j].chart; c.va = A[i].v; c.vb = B[j].v; }
        else { c.chart_a = B[j].chart; c.chart_b = A[i].chart; c.va = B[j].v; c.vb = A[i].v; }
        c.dist = sqrtf(d2); c.ngap = (float)ng; c.sign = (int8_t)(nd >= 0.0 ? 1 : -1);
        ar_push_cand(buf, &c);
    }
    Arena_restore(scratch, mark);
}

/* ---- relation from a group of candidates ------------------------------------- */

static int ar_cmp_cand(const void *x, const void *y)
{
    const ArCand *a = x, *b = y;
    if (a->chart_a != b->chart_a) return a->chart_a < b->chart_a ? -1 : 1;
    if (a->chart_b != b->chart_b) return a->chart_b < b->chart_b ? -1 : 1;
    if (a->va != b->va) return a->va < b->va ? -1 : 1;
    return (a->vb > b->vb) - (a->vb < b->vb);
}

static int ar_cmp_float(const void *x, const void *y)
{
    float a = *(const float *)x, b = *(const float *)y;
    return (a > b) - (a < b);
}

static double ar_median_f(float *v, size_t n, Arena_T scratch)
{
    Arena_Mark m = Arena_save(scratch);
    float *c = ARENA_ALLOC(scratch, n * sizeof(float));
    memcpy(c, v, n * sizeof(float));
    qsort(c, n, sizeof(float), ar_cmp_float);
    double r = (n % 2) ? c[n/2] : 0.5 * (c[n/2-1] + c[n/2]);
    Arena_restore(scratch, m);
    return r;
}

static uint32_t ar_rng(uint32_t *s)
{
    *s ^= *s << 13; *s ^= *s >> 17; *s ^= *s << 5;
    return *s;
}

typedef struct ArBridge { double uni; size_t corr; int32_t rel; } ArBridge;

static int ar_cmp_bridge(const void *x, const void *y)
{
    const ArBridge *a = (const ArBridge *)x, *b = (const ArBridge *)y;
    if (a->uni > b->uni) return -1;
    if (a->uni < b->uni) return 1;
    if (a->corr > b->corr) return -1;
    if (a->corr < b->corr) return 1;
    return (a->rel > b->rel) - (a->rel < b->rel);
}

/* The first relation joining two otherwise-disconnected blocks is LOAD-BEARING
 * and uncorroborated: nothing else holds those blocks together, so one marginal
 * seam welds them into a single component, and an over-merged block then
 * occupies at placement the space its neighbours need.  A relation that only
 * corroborates an existing connection is checked by its neighbours and keeps
 * the ordinary bar.  Bridges are judged BEST FIRST, so the premium falls on the
 * strongest candidate for each merge; one that misses it stays as placement
 * evidence (WEAK) and never becomes a join.  A seam whose unimodality was not
 * measured (too few spans) keeps the ordinary bar. */
/* Cut edges of the graph the candidate relations form over the charts: an edge
 * whose removal disconnects its endpoints, so nothing else in the accepted
 * evidence connects the two blocks.  Iterative Tarjan lowlink; parallel edges
 * are distinguished by index, so two relations between one chart pair
 * corroborate each other.  `is_cut` is indexed by position in `ord`. */
static void ar_cut_edges(const AsmRun *run, const ArBridge *ord, size_t nk, uint8_t *is_cut, Arena_T scratch)
{
    size_t nc = run->n_charts;
    Arena_Mark mark = Arena_save(scratch);
    size_t *head = ARENA_CALLOC(scratch, nc + 1, sizeof *head);
    for (size_t k = 0; k < nk; k++) {
        const AsmRelation *R = &run->rels[ord[k].rel];
        head[(size_t)R->a + 1]++; head[(size_t)R->b + 1]++;
    }
    for (size_t i = 0; i < nc; i++) head[i+1] += head[i];
    size_t *fill = ARENA_ALLOC(scratch, (nc ? nc : 1) * sizeof *fill);
    memcpy(fill, head, nc * sizeof *fill);
    size_t nadj = 2 * nk ? 2 * nk : 1;
    int32_t *adj_v = ARENA_ALLOC(scratch, nadj * sizeof *adj_v);
    int32_t *adj_e = ARENA_ALLOC(scratch, nadj * sizeof *adj_e);
    for (size_t k = 0; k < nk; k++) {
        const AsmRelation *R = &run->rels[ord[k].rel];
        adj_v[fill[R->a]] = R->b; adj_e[fill[R->a]] = (int32_t)k; fill[R->a]++;
        adj_v[fill[R->b]] = R->a; adj_e[fill[R->b]] = (int32_t)k; fill[R->b]++;
    }
    size_t *disc = ARENA_CALLOC(scratch, nc ? nc : 1, sizeof *disc);
    size_t *low = ARENA_ALLOC(scratch, (nc ? nc : 1) * sizeof *low);
    size_t *cursor = ARENA_ALLOC(scratch, (nc ? nc : 1) * sizeof *cursor);
    int32_t *parent_e = ARENA_ALLOC(scratch, (nc ? nc : 1) * sizeof *parent_e);
    int32_t *stack = ARENA_ALLOC(scratch, (nc ? nc : 1) * sizeof *stack);
    size_t timer = 1;
    for (size_t s = 0; s < nc; s++) {
        if (disc[s]) continue;
        size_t top = 0;
        stack[top] = (int32_t)s; cursor[s] = head[s]; parent_e[s] = -1;
        disc[s] = low[s] = timer++;
        while (1) {
            int32_t v = stack[top];
            if (cursor[(size_t)v] < head[(size_t)v + 1]) {
                size_t slot = cursor[(size_t)v]++;
                int32_t w = adj_v[slot], e = adj_e[slot];
                if (e == parent_e[(size_t)v]) continue;
                if (disc[(size_t)w]) {
                    if (disc[(size_t)w] < low[(size_t)v]) low[(size_t)v] = disc[(size_t)w];
                    continue;
                }
                disc[(size_t)w] = low[(size_t)w] = timer++;
                cursor[(size_t)w] = head[(size_t)w]; parent_e[(size_t)w] = e;
                stack[++top] = w;
                continue;
            }
            if (top == 0) break;
            int32_t p = stack[top - 1];
            if (low[(size_t)v] < low[(size_t)p]) low[(size_t)p] = low[(size_t)v];
            if (low[(size_t)v] > disc[(size_t)p]) is_cut[(size_t)parent_e[(size_t)v]] = 1;
            top--;
        }
    }
    Arena_restore(scratch, mark);
}

/* A block below this area cannot be the weld the premium guards against (see
 * ASM_SEAM_BRIDGE_LEAF_AREA): a seam attaching it keeps the ordinary bar. */
static int ar_bridge_leaf(const double *block, int32_t ra, int32_t rb, double leaf_area)
{
    return block && leaf_area > 0.0 && fmin(block[ra], block[rb]) < leaf_area;
}

static size_t ar_bridge_premium(AsmRun *run, const double *unimodal, double bar, double leaf_area, Arena_T scratch,
                                size_t *corroborated, size_t *leaves)
{
    if (leaves) *leaves = 0;
    if (!(bar > 0.0) || run->n_rels == 0 || run->n_charts == 0) return 0;
    Arena_Mark mark = Arena_save(scratch);
    ArBridge *ord = ARENA_ALLOC(scratch, run->n_rels * sizeof *ord);
    size_t nk = 0;
    for (size_t r = 0; r < run->n_rels; r++) {
        const AsmRelation *R = &run->rels[r];
        if (R->flags & (ASM_REL_DROPPED | ASM_REL_CONTACT | ASM_REL_PLACEMENT_ONLY)) continue;
        ord[nk].uni = unimodal[r]; ord[nk].corr = R->n_corr; ord[nk].rel = (int32_t)r; nk++;
    }
    qsort(ord, nk, sizeof *ord, ar_cmp_bridge);
    uint8_t *is_cut = ARENA_CALLOC(scratch, nk ? nk : 1, 1);
    if (ASM_SEAM_BRIDGE_CORROBORATED) ar_cut_edges(run, ord, nk, is_cut, scratch);
    UnionFind uf = UF_new(scratch, (int32_t)run->n_charts);
    /* the area of each block so far, at its root (leaf exemption only) */
    double *block = leaf_area > 0.0 ? ARENA_ALLOC(scratch, run->n_charts * sizeof *block) : NULL;
    if (block) for (size_t c = 0; c < run->n_charts; c++) block[c] = run->charts ? run->charts[c].area3d : 0.0;
    size_t demoted = 0, spared = 0, leafs = 0;
    for (size_t k = 0; k < nk; k++) {
        AsmRelation *R = &run->rels[ord[k].rel];
        int32_t ra = uf_find(&uf, R->a), rb = uf_find(&uf, R->b);
        if (ra == rb) continue;   /* corroborated: the ordinary bar stands */
        if (ord[k].uni >= 0.0 && ord[k].uni < bar) {
            /* This pass asks "is it a bridge SO FAR", which demotes both of two
             * independent marginal seams joining one pair of blocks: the first is
             * uncorroborated when it is tried, the second because the first went.
             * A cut edge of the whole accepted graph is the honest test. */
            if (ar_bridge_leaf(block, ra, rb, leaf_area)) leafs++;
            else if (!ASM_SEAM_BRIDGE_CORROBORATED || is_cut[k]) { R->flags |= ASM_REL_WEAK; demoted++; continue; }
            else spared++;
        }
        uf_union(&uf, R->a, R->b);
        if (block) { double sum = block[ra] + block[rb]; block[uf_find(&uf, R->a)] = sum; }
    }
    if (corroborated) *corroborated = spared;
    if (leaves) *leaves = leafs;
    Arena_restore(scratch, mark);
    return demoted;
}

/* Returns 1 and fills rel/corr when the group yields a relation. */
typedef struct ArDiag {
    int32_t chart_a, chart_b; size_t n;
    double seam_len, normal_gap, mean_gap, iso_disc, rms, unimodal, area_a, area_b;
    int parity; int outcome;  /* 0 accepted, 1 tangency, 2 length, 3 iso, 4 unimodal, 5 fit */
    int kept;                 /* 0 rejected, 1 accepted, 2 kept as WEAK, 3 kept as SHORT */
    int crosswrap;            /* the unimodality gate rejected it: a sheet switch inside the seam */
    int trim_transport;
} ArDiag;

static FILE *ar_diag_fp = NULL;

static void ar_diag_write(const ArDiag *d)
{
    if (!ar_diag_fp) return;
    fprintf(ar_diag_fp, "%d,%d,%zu,%.3f,%.3f,%.3f,%.4f,%.3f,%.3f,%.1f,%.1f,%d,%d,%d,%d,%d\n",
            d->chart_a, d->chart_b, d->n, d->seam_len, d->normal_gap, d->mean_gap, d->iso_disc, d->rms, d->unimodal,
            d->area_a, d->area_b, d->parity, d->outcome, d->kept,d->trim_transport, d->crosswrap);
}

static int ar_group_relation(const AsmRun *run, const ArCand *g, size_t n, const ArTrim *trim, const AsmRelateOpts *o,
                             Arena_T scratch, AsmRelation *rel, AsmRelateStats *st, double *unimodal_out)
{
    ArDiag dg; memset(&dg, 0, sizeof dg);
    if (unimodal_out) *unimodal_out = -1.0;
    dg.chart_a = g[0].chart_a; dg.chart_b = g[0].chart_b; dg.n = n;
    dg.area_a = run->charts[(size_t)g[0].chart_a].area3d; dg.area_b = run->charts[(size_t)g[0].chart_b].area3d;
    dg.iso_disc = -1.0; dg.rms = -1.0; dg.unimodal = -1.0;
    Arena_Mark mark = Arena_save(scratch);
    const AsmChart *A = &run->charts[(size_t)g[0].chart_a];
    const AsmChart *B = &run->charts[(size_t)g[0].chart_b];
    /* parity by majority sign */
    long sgn = 0;
    for (size_t i = 0; i < n; i++) sgn += g[i].sign;
    int parity = sgn < 0;
    dg.parity = parity;
    int weak = 0, short_seam = 0, crosswrap = 0;
    double theta = 0.0, tx = 0.0, ty = 0.0, rms = 0.0;
    /* tangency: median normal gap */
    float *ng = ARENA_ALLOC(scratch, n * sizeof(float));
    float *dd = ARENA_ALLOC(scratch, n * sizeof(float));
    for (size_t i = 0; i < n; i++) { ng[i] = g[i].ngap; dd[i] = g[i].dist; }
    double normal_gap = ar_median_f(ng, n, scratch);
    double mean_gap = ar_median_f(dd, n, scratch);
    dg.normal_gap = normal_gap; dg.mean_gap = mean_gap;
    if (normal_gap > o->tangency_max_vox) {
        st->rej_tangency++; dg.outcome = 1; ar_diag_write(&dg); Arena_restore(scratch, mark);
        if (!ASM_RELATE_LEDGER_CONTACT) return 0;
        /* the ledger row: two wraps meeting face to face across this cube face -- a veto for layout
         * joins (ASM_JOIN_VETO_FLAGS), never a pose-graph edge or placement evidence, no correspondences */
        memset(rel, 0, sizeof *rel);
        rel->a = g[0].chart_a; rel->b = g[0].chart_b; rel->flags = ASM_REL_CONTACT | (parity ? ASM_REL_PARITY : 0u);
        rel->normal_gap = normal_gap; rel->mean_gap = mean_gap; rel->n_corr = n; rel->corr_first = -1; rel->corr_count = 0; rel->rms = -1.0;
        return 4;
    }
    /* matched length: principal extent of the A-side 3-D points */
    double c[3] = { 0, 0, 0 };
    for (size_t i = 0; i < n; i++) for (int d = 0; d < 3; d++) c[d] += A->xyz[(size_t)g[i].va*3 + (size_t)d];
    for (int d = 0; d < 3; d++) c[d] /= (double)n;
    double cov[6] = { 0, 0, 0, 0, 0, 0 };
    for (size_t i = 0; i < n; i++) {
        double x = A->xyz[(size_t)g[i].va*3] - c[0], y = A->xyz[(size_t)g[i].va*3+1] - c[1], z = A->xyz[(size_t)g[i].va*3+2] - c[2];
        cov[0] += x*x; cov[1] += x*y; cov[2] += x*z; cov[3] += y*y; cov[4] += y*z; cov[5] += z*z;
    }
    /* A fixed power-iteration seed can be perpendicular to the principal
     * seam direction and report zero length for a long diagonal boundary.
     * Use the symmetric eigensolver so this gate is rotation invariant. */
    double matrix[3][3] = {{cov[0],cov[1],cov[2]},
                           {cov[1],cov[3],cov[4]},
                           {cov[2],cov[4],cov[5]}};
    double evals[3], evecs[3][3]; Eig3_sym(matrix,evals,evecs);
    double lam = fmax(0.0,evals[2]);
    double seam_len = 2.0 * sqrt(3.0) * sqrt(lam / (double)n);
    dg.seam_len = seam_len;
    if (seam_len < o->seam_min_len) {
        st->rej_len++; dg.outcome = 2;   /* the gate's verdict stands for the pose graph */
        if (!(AR_SHORT_SEAMS && seam_len >= AR_SHORT_MIN_LEN && n >= (size_t)AR_SHORT_MIN_CORR)) {
            if (AR_SHORT_SEAMS && seam_len >= AR_SHORT_MIN_LEN) st->short_rej_corr++;
            ar_diag_write(&dg); Arena_restore(scratch, mark); return 0;
        }
        short_seam = 1;   /* kept as placement evidence when its correspondences fit rigidly (below) */
    }
    /* Trimmed boundary vertices are different physical points. Compare
     * transported midpoints for the pairwise gates, and fit BOTH directed
     * endpoint observations. Their half weights retain one physical vote.
     * An unresolved tangent keeps the entire original raw-point trial;
     * silently removing that observation could make a bad seam pass. */
    int transported = trim != NULL;
    for (size_t i = 0; i < n && transported; i++) transported = trim[i].valid == 3;
    dg.trim_transport = transported;
    double *qa = ARENA_ALLOC(scratch, n * 2 * sizeof(double));
    double *qb = ARENA_ALLOC(scratch, n * 2 * sizeof(double));
    double *fit_a = transported ? ARENA_ALLOC(scratch,n*4*sizeof(double)) : qa;
    double *fit_b = transported ? ARENA_ALLOC(scratch,n*4*sizeof(double)) : qb;
    size_t fit_n = transported ? 2*n : n;
    for (size_t i = 0; i < n; i++) {
        qa[i*2] = A->uv[(size_t)g[i].va*2]; qa[i*2+1] = A->uv[(size_t)g[i].va*2+1];
        double ub = B->uv[(size_t)g[i].vb*2], vb = B->uv[(size_t)g[i].vb*2+1];
        qb[i*2] = parity ? -ub : ub; qb[i*2+1] = vb;
        if (transported) for (int q = 0; q < 2; q++) {
            double ga = trim[i].gap_a[q], gb = trim[i].gap_b[q];
            if (parity && q == 0) gb = -gb;
            fit_a[2*i+q] = qa[2*i+q]+ga; fit_a[2*(n+i)+q] = qa[2*i+q];
            fit_b[2*i+q] = qb[2*i+q]; fit_b[2*(n+i)+q] = qb[2*i+q]-gb;
        }
    }
    /* TWO fits, two purposes.  gap_a/gap_b are the pair displacement b-a
     * re-expressed in each chart's frame, so a transported pair agrees with
     * itself by construction: one free parameter per observation.  That makes
     * the transported fit an excellent estimate of the frame-to-frame
     * transform and a useless measure of whether the two boundaries describe
     * the same curve -- transporting the midpoints drove the median
     * unimodality of 34 bimodal seams from 0.639 to 1.000 and their median
     * iso discrepancy from 0.0561 to 0.0040 (familiar 4x5x5 A/B 2026-09-15),
     * admitting wrap-crossing seams that over-merge the graph.  The pairwise
     * gates and the seam's published quality therefore read the UNTRANSPORTED
     * midpoints, where a sheet switch inside a seam stays visible as a step. */
    double raw_theta = 0.0, raw_tx = 0.0, raw_ty = 0.0, raw_rms = 0.0;
    int raw_fitted = AsmRelate_rigid_fit(qa, qb, NULL, n, &raw_theta, &raw_tx, &raw_ty, &raw_rms) == 0;
    if (short_seam) {
        /* a SHORT seam skips the isometry and unimodality tests (they need 8 trials over 8-vox spans) */
        if (!raw_fitted || AsmRelate_rigid_fit(fit_a, fit_b, NULL, fit_n, &theta, &tx, &ty, &rms) != 0 ||
            raw_rms > AR_RIGID_MAX_RMS) {
            st->short_rej_rms++; ar_diag_write(&dg); Arena_restore(scratch, mark); return 0;
        }
        rms = raw_rms;
        dg.rms = rms;
        goto fill;
    }
    /* isometry check on random pairs: the median relative discrepancy of
     * the two charts' seam parameterizations is the measurement; the gate is
     * on that median (a sheet switch inside the seam is a step, not a slope) */
    uint32_t seed = (uint32_t)(0x9E3779B9u * (uint32_t)(A->id * 131 + B->id) + 1u);
    size_t trials = n * 4 > 256 ? 256 : n * 4, tried = 0;
    float *disc = ARENA_ALLOC(scratch, (trials ? trials : 1) * sizeof(float));
    for (size_t t = 0; t < trials; t++) {
        size_t i = ar_rng(&seed) % n, j = ar_rng(&seed) % n;
        if (i == j) continue;
        double da = hypot(qa[i*2]-qa[j*2], qa[i*2+1]-qa[j*2+1]);
        double db = hypot(qb[i*2]-qb[j*2], qb[i*2+1]-qb[j*2+1]);
        if (da < 8.0) continue;
        disc[tried++] = (float)(fmax(0.0, fabs(da - db) - o->seam_noise_vox) / da);   /* one noise budget per pair */
    }
    double iso_disc = tried ? ar_median_f(disc, tried, scratch) : 0.0;
    st->iso_disc_sum += iso_disc; st->iso_disc_n++;
    dg.iso_disc = iso_disc;
    if (tried >= 8 && iso_disc > o->iso_tol) { st->rej_iso++; dg.outcome = 3; weak = 1; }
    if (!raw_fitted || AsmRelate_rigid_fit(fit_a, fit_b, NULL, fit_n, &theta, &tx, &ty, &rms) != 0) { if (!weak) { st->rej_fit++; dg.outcome = 5; } ar_diag_write(&dg); Arena_restore(scratch, mark); return 0; }
    rms = raw_rms;   /* the transform is transported; its quality is not */
    dg.rms = rms;
    /* unimodality: pairwise rotations must agree with the fit.  This is the
     * cross-wrap detector, so it is measured for an isometry-rejected pair too
     * (ASM_WRAP_MEASURE_WEAK): that pair is still handed to placement as an
     * anchor, and until 2026-09-17 it went there untested. */
    int iso_weak = weak;
    if (!iso_weak || ASM_WRAP_MEASURE_WEAK) {
        size_t agree = 0, checked = 0;
        for (size_t t = 0; t < trials; t++) {
            size_t i = ar_rng(&seed) % n, j = ar_rng(&seed) % n;
            if (i == j) continue;
            double da = hypot(qa[i*2]-qa[j*2], qa[i*2+1]-qa[j*2+1]);
            if (da < 8.0) continue;
            double ang = atan2(qa[i*2+1]-qa[j*2+1], qa[i*2]-qa[j*2]) - atan2(qb[i*2+1]-qb[j*2+1], qb[i*2]-qb[j*2]);
            double diff = ang - theta;
            while (diff > AR_PI) diff -= 2.0 * AR_PI;
            while (diff < -AR_PI) diff += 2.0 * AR_PI;
            checked++;
            if (fabs(diff) < 0.05 + o->seam_noise_vox / da) agree++;   /* the angle a noise budget subtends over this span */
        }
        dg.unimodal = checked ? (double)agree / (double)checked : -1.0;
        if (checked >= 8 && (double)agree < 0.8 * (double)checked) {
            crosswrap = 1; weak = 1;
            /* the first gate to fire owns the outcome; the wrap verdict is its own flag */
            if (!iso_weak) { st->rej_unimodal++; dg.outcome = 4; }
            else st->wrap_on_iso++;
        }
    }
    dg.crosswrap = crosswrap;
    if (rms > AR_RIGID_MAX_RMS) {
        if (!weak) { st->rej_fit++; dg.outcome = 5; }
        ar_diag_write(&dg); Arena_restore(scratch, mark); return 0;
    }
    if (!weak) { dg.kept = 1; ar_diag_write(&dg); }
fill:
    if (weak || short_seam) { dg.kept = short_seam ? 3 : 2; ar_diag_write(&dg); }
    if (unimodal_out) *unimodal_out = dg.unimodal;
    memset(rel, 0, sizeof *rel);
    rel->a = A->id; rel->b = B->id;
    rel->theta = theta; rel->tx = tx; rel->ty = ty;
    double var = rms * rms; if (var < 0.25) var = 0.25;
    rel->w_xy = (double)n / var;
    double s2 = 0.0, cb0 = 0.0, cb1 = 0.0;
    for (size_t i = 0; i < n; i++) { cb0 += qb[i*2]; cb1 += qb[i*2+1]; }
    cb0 /= (double)n; cb1 /= (double)n;
    for (size_t i = 0; i < n; i++) s2 += (qb[i*2]-cb0)*(qb[i*2]-cb0) + (qb[i*2+1]-cb1)*(qb[i*2+1]-cb1);
    rel->w_theta = s2 / var;
    rel->rms = rms; rel->seam_len = seam_len; rel->normal_gap = normal_gap; rel->mean_gap = mean_gap;
    rel->n_corr = n;
    rel->robust_w = 1.0;
    if (parity) rel->flags |= ASM_REL_PARITY;
    if (weak) rel->flags |= ASM_REL_WEAK;
    if (short_seam) rel->flags |= ASM_REL_SHORT;
    if (crosswrap) rel->flags |= ASM_REL_CROSSWRAP;
    if ((A->flags & ASM_CHART_SUSPECT) || (B->flags & ASM_CHART_SUSPECT)) rel->flags |= ASM_REL_SUSPECT;
    if (A->repaired || B->repaired) {
        size_t nrep = 0;
        for (size_t i = 0; i < n; i++)
            nrep += (A->repaired && A->repaired[g[i].va]) || (B->repaired && B->repaired[g[i].vb]);
        if (nrep * 4 > n) rel->flags |= ASM_REL_REPAIRED;
    }
    Arena_restore(scratch, mark);
    return short_seam ? 3 : (weak ? 2 : 1);
}

/* ---- layer probes for one cube ------------------------------------------------ */

typedef struct ArFace { float c[3]; float n[3]; int32_t chart; int32_t f; } ArFace;

static void ar_face_normal(const AsmChart *c, size_t f, float n[3])
{
    const int32_t *fv = &c->faces[f*3];
    for (int d = 0; d < 3; d++) n[d] = c->nrm[(size_t)fv[0]*3+(size_t)d] + c->nrm[(size_t)fv[1]*3+(size_t)d] + c->nrm[(size_t)fv[2]*3+(size_t)d];
    float l = sqrtf(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
    if (l > 0) { n[0] /= l; n[1] /= l; n[2] /= l; }
}

/* A point within one voxel of a triangle is within (1 + its maximum
 * vertex-to-centroid distance) of the stored centroid. Use the actual
 * geometry, including large/elongated faces. Outward rounding covers the
 * float distance arithmetic in the centroid tree. */
static float ar_probe_radius_sq(double extent_sq)
{
    double radius = fmax(AR_PROBE_MIN_RADIUS, 1.0 + sqrt(extent_sq));
    double squared = radius*radius*(1.0 + 32.0*FLT_EPSILON);
    return squared >= FLT_MAX ? INFINITY : nextafterf((float)squared, INFINITY);
}

/* KDTree_ball_query reports the number written, not the number found. A
 * full buffer is therefore ambiguous: replay with room for every face.
 * Retain the expanded scratch buffer for the rest of this cube. */
static size_t ar_probe_ball(Arena_T arena, KDTree_T tree, size_t nf,
                            const float q[3], float radius_sq,
                            int32_t **ball, size_t *capacity)
{
    size_t count = KDTree_ball_query(tree,q,radius_sq,*ball,*capacity);
    if (count == *capacity && *capacity < nf) {
        *capacity = nf;
        *ball = ARENA_ALLOC(arena,nf*sizeof(int32_t));
        count = KDTree_ball_query(tree,q,radius_sq,*ball,*capacity);
    }
    return count;
}

static void ar_layer_cube(const AsmRun *run, const ArCubeIndex *ix, size_t cube,
                          const AsmRelateOpts *o, Arena_T scratch, ArBuf *buf)
{
    Arena_Mark mark = Arena_save(scratch);
    /* faces of the 27-neighbourhood */
    size_t nf_total = 0;
    int32_t cubes[27]; int ncubes = 0;
    cubes[ncubes++] = (int32_t)cube;
    for (int k = 0; k < 26; k++) if (ix->nbr[cube*26 + (size_t)k] >= 0) cubes[ncubes++] = ix->nbr[cube*26 + (size_t)k];
    for (int q = 0; q < ncubes; q++)
        for (int32_t k = ix->chart_first[cubes[q]]; k < ix->chart_first[cubes[q]+1]; k++)
            nf_total += run->charts[(size_t)ix->chart_ids[k]].nf;
    if (nf_total == 0 || ix->chart_first[cube] == ix->chart_first[cube+1]) { Arena_restore(scratch, mark); return; }
    ArFace *F = ARENA_ALLOC(scratch, nf_total * sizeof(ArFace));
    float *fc = ARENA_ALLOC(scratch, nf_total * 3 * sizeof(float));
    size_t nf = 0;
    double extent_sq = 0.0;
    for (int q = 0; q < ncubes; q++) {
        for (int32_t k = ix->chart_first[cubes[q]]; k < ix->chart_first[cubes[q]+1]; k++) {
            const AsmChart *c = &run->charts[(size_t)ix->chart_ids[k]];
            for (size_t f = 0; f < c->nf; f++) {
                const int32_t *fv = &c->faces[f*3];
                for (int d = 0; d < 3; d++)
                    F[nf].c[d] = (c->xyz[(size_t)fv[0]*3+(size_t)d] + c->xyz[(size_t)fv[1]*3+(size_t)d] + c->xyz[(size_t)fv[2]*3+(size_t)d]) / 3.0f;
                for (int vertex = 0; vertex < 3; vertex++) {
                    double distance_sq = 0.0;
                    for (int d = 0; d < 3; d++) {
                        double delta = (double)c->xyz[(size_t)fv[vertex]*3+(size_t)d] - F[nf].c[d];
                        distance_sq += delta*delta;
                    }
                    extent_sq = fmax(extent_sq,distance_sq);
                }
                ar_face_normal(c, f, F[nf].n);
                F[nf].chart = c->id; F[nf].f = (int32_t)f;
                memcpy(&fc[nf*3], F[nf].c, 3 * sizeof(float));
                nf++;
            }
        }
    }
    KDTree_T tree = KDTree_new(scratch, fc, nf);
    size_t ball_capacity = nf < AR_PROBE_INITIAL_CAP ? nf : AR_PROBE_INITIAL_CAP;
    int32_t *ball = ARENA_ALLOC(scratch, ball_capacity * sizeof(int32_t));
    const float ball_radius_sq = ar_probe_radius_sq(extent_sq);
    for (int32_t k = ix->chart_first[cube]; k < ix->chart_first[cube+1]; k++) {
        const AsmChart *c = &run->charts[(size_t)ix->chart_ids[k]];
        size_t stride = c->nv / (size_t)(o->layer_samples > 0 ? o->layer_samples : 1);
        if (stride < 1) stride = 1;
        for (size_t v = 0; v < c->nv; v += stride) {
            const float *p = &c->xyz[v*3];
            const float *n = &c->nrm[v*3];
            for (int side = -1; side <= 1; side += 2) {
                /* MULTI-CROSSING (2026-09-09): the probe keeps marching after a surface and records
                 * every crossing with its ORDER, so a pair's wrap count is measured per ray instead
                 * of inferred from a distance or a coordinate.  At the core (6-9 vox spacing, holes
                 * in the mesh) the first surface met is often the second wrap; the median order over
                 * a pair's hits is robust to that where most rays see the near wrap. */
                int found = 0, last_chart = -1;
                double skip_until = 0.0;
                for (double s = 3.0; s <= o->layer_probe_vox && found < 3; s += 1.0) {
                    if (s < skip_until) continue;
                    float q[3] = { (float)(p[0] + side*s*n[0]), (float)(p[1] + side*s*n[1]), (float)(p[2] + side*s*n[2]) };
                    size_t m = ar_probe_ball(scratch,tree,nf,q,ball_radius_sq,&ball,&ball_capacity);
                    double best = 1e30; int32_t best_chart = -1, best_face = -1; double best_u = 0.0, best_v = 0.0, best_nd = 0.0;
                    for (size_t b = 0; b < m; b++) {
                        const ArFace *fa = &F[(size_t)ball[b]];
                        if (fa->chart == c->id) continue;
                        double nd = fa->n[0]*n[0] + fa->n[1]*n[1] + fa->n[2]*n[2];
                        if (fabs(nd) < 0.6) continue;
                        const AsmChart *oc = &run->charts[(size_t)fa->chart];
                        const int32_t *fv = &oc->faces[(size_t)fa->f*3];
                        double qd[3] = { q[0], q[1], q[2] };
                        double a[3] = { oc->xyz[(size_t)fv[0]*3], oc->xyz[(size_t)fv[0]*3+1], oc->xyz[(size_t)fv[0]*3+2] };
                        double bb[3] = { oc->xyz[(size_t)fv[1]*3], oc->xyz[(size_t)fv[1]*3+1], oc->xyz[(size_t)fv[1]*3+2] };
                        double cc[3] = { oc->xyz[(size_t)fv[2]*3], oc->xyz[(size_t)fv[2]*3+1], oc->xyz[(size_t)fv[2]*3+2] };
                        double cp[3], cu, cv2;
                        ClosestTri_point(qd, a, bb, cc, cp, &cu, &cv2);
                        double d2 = (cp[0]-qd[0])*(cp[0]-qd[0]) + (cp[1]-qd[1])*(cp[1]-qd[1]) + (cp[2]-qd[2])*(cp[2]-qd[2]);
                        if (d2 < best) { best = d2; best_chart = fa->chart; best_face = fa->f; best_u = cu; best_v = cv2; best_nd = nd; }
                    }
                    if (best_chart >= 0 && best <= 1.0) {
                        if (best_chart == last_chart) continue;   /* still inside the surface just crossed */
                        found++;
                        ArHit h; h.chart = c->id; h.other = best_chart; h.side = (int8_t)side; h.dist = (float)s;
                        h.va = (int32_t)v; h.fb = best_face; h.order = (int8_t)found; h.nsign = (int8_t)(best_nd < 0.0 ? -1 : 1);
                        /* ClosestTri_point: out = a + u*(b-a) + v*(c-a) -> l1 = u, l2 = v, l0 = 1-u-v */
                        h.l0 = (float)(1.0 - best_u - best_v); h.l1 = (float)best_u;
                        ar_push_hit(buf, &h);
                        last_chart = best_chart;
                        skip_until = s + 3.0;   /* past the papyrus thickness before the next surface can count */
                    }
                }
            }
        }
    }
    Arena_restore(scratch, mark);
}

static int ar_cmp_hit(const void *x, const void *y)
{
    const ArHit *a = x, *b = y;
    if (a->chart != b->chart) return a->chart < b->chart ? -1 : 1;
    if (a->other != b->other) return a->other < b->other ? -1 : 1;
    if (a->side != b->side) return a->side < b->side ? -1 : 1;
    if (a->dist != b->dist) return (a->dist > b->dist) - (a->dist < b->dist);
    /* TOTAL order: the hits are merged from per-thread buffers in scheduling order, and the
     * representative 16 around the median are taken from this sort, so a tie in dist made
     * the placement non-deterministic (three identical runs placed 48/55/49 components,
     * growth 119.8 vs 120.2 vox; measured 2026-09-08) */
    if (a->va != b->va) return a->va < b->va ? -1 : 1;
    if (a->fb != b->fb) return a->fb < b->fb ? -1 : 1;
    if (a->l0 != b->l0) return (a->l0 > b->l0) - (a->l0 < b->l0);
    return (a->l1 > b->l1) - (a->l1 < b->l1);
}

static int ar_cmp_double(const void *x, const void *y)
{
    double a = *(const double *)x, b = *(const double *)y;
    return (a > b) - (a < b);
}

/* ---- driver -------------------------------------------------------------------- */

int AsmRelate_run(AsmRun *run, const AsmRelateOpts *o, int threads, AsmRelateStats *st)
{
    memset(st, 0, sizeof *st);
    double t0 = ves_clock_sec();
    run->cube_size = o->cube_size;
    ArCubeIndex ix;
    ar_build_cube_index(run, o->cube_size, &ix);

    /* unordered cube pairs (each once) */
    size_t npairs = 0;
    for (size_t c = 0; c < run->n_cubes; c++)
        for (int k = 0; k < 26; k++)
            if (ix.nbr[c*26 + (size_t)k] > (int32_t)c) npairs++;
    int32_t *pairs = ARENA_ALLOC(run->arena, (npairs ? npairs : 1) * 2 * sizeof(int32_t));
    {
        size_t q = 0;
        for (size_t c = 0; c < run->n_cubes; c++)
            for (int k = 0; k < 26; k++)
                if (ix.nbr[c*26 + (size_t)k] > (int32_t)c) { pairs[q*2] = (int32_t)c; pairs[q*2+1] = ix.nbr[c*26 + (size_t)k]; q++; }
    }
    st->cube_pairs = npairs;
    if (threads < 1) threads = 1;
    Arena_T *scratch = ARENA_ALLOC(run->arena, (size_t)threads * sizeof(Arena_T));
    ArBuf *bufs = ARENA_CALLOC(run->arena, (size_t)threads, sizeof(ArBuf));
    for (int t = 0; t < threads; t++) { scratch[t] = Arena_new(); bufs[t].arena = scratch[t]; }
    /* thread buffers live in dedicated arenas so scratch save/restore inside
     * the workers cannot clobber them */
    Arena_T *work = ARENA_ALLOC(run->arena, (size_t)threads * sizeof(Arena_T));
    for (int t = 0; t < threads; t++) work[t] = Arena_new();

    int pi = 0, npairs_i = (int)npairs;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 4) num_threads(threads)
#endif
    for (pi = 0; pi < npairs_i; pi++) {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        ar_seam_pair(run, &ix, (size_t)pairs[pi*2], (size_t)pairs[pi*2+1], o, work[tid], &bufs[tid]);
    }
    /* merge candidates */
    size_t ncand = 0;
    for (int t = 0; t < threads; t++) ncand += bufs[t].n_cand;
    ArCand *cand = ARENA_ALLOC(run->arena, (ncand ? ncand : 1) * sizeof(ArCand));
    {
        size_t q = 0;
        for (int t = 0; t < threads; t++) { memcpy(&cand[q], bufs[t].cand, bufs[t].n_cand * sizeof(ArCand)); q += bufs[t].n_cand; }
    }
    st->candidates = ncand;
    for (int t = 0; t < threads; t++) {
        st->gap_pairs += bufs[t].gap_pairs;
        for (int h = 0; h < ASM_RELATE_HIST; h++) st->gap_hist[h] += bufs[t].gap_hist[h];
    }
    qsort(cand, ncand, sizeof(ArCand), ar_cmp_cand);
    ArTrim *trim = ARENA_CALLOC(run->arena,ncand ? ncand : 1,sizeof *trim);
    if (ar_candidate_trims(run,cand,ncand,work[0],trim)) {
        for (int t = 0; t < threads; t++) { Arena_dispose(&scratch[t]); Arena_dispose(&work[t]); }
        return -1;
    }
    /* groups -> relations (serial: the run tables are not thread-safe) */
    if (o->diag_csv && o->diag_csv[0]) {
        ar_diag_fp = fopen(o->diag_csv, "wb");
        if (ar_diag_fp) fprintf(ar_diag_fp, "chart_a,chart_b,n,seam_len,normal_gap,mean_gap,iso_disc,rms,unimodal,area_a,area_b,parity,outcome,kept,trim_transport,crosswrap\n");
    }
    double *rms_list = ARENA_ALLOC(run->arena, (ncand ? ncand : 1) * sizeof(double));
    /* Indexed by RELATION index, so it must cover relations that already exist
     * on entry as well as every one this pass can push (at most one per three
     * candidates); a relation this pass did not measure reads -1 = unmeasured
     * and keeps the ordinary bar. */
    size_t n_unimodal = run->n_rels + ncand + 1;
    double *unimodal = ARENA_ALLOC(run->arena, n_unimodal * sizeof(double));
    for (size_t k = 0; k < n_unimodal; k++) unimodal[k] = -1.0;
    size_t nrms = 0, n_weak = 0;
    AsmRelation *contact_rows = NULL; size_t n_contact = 0, cap_contact = 0;
    for (size_t i = 0; i < ncand;) {
        size_t j = i;
        while (j < ncand && cand[j].chart_a == cand[i].chart_a && cand[j].chart_b == cand[i].chart_b) j++;
        size_t n = j - i;
        if (n >= 3) {
            st->chart_pairs++;
            AsmRelation rel;
            double uni = -1.0;
            int gr = ar_group_relation(run, &cand[i], n, &trim[i], o, work[0], &rel, st, &uni);
            if (gr == 2 || gr == 3) {
                unimodal[run->n_rels] = uni;
                ar_append_relation(run,&rel,cand+i,trim+i,n);
                if (gr == 3) st->short_kept++; else n_weak++;
            } else if (gr == 1) {
                unimodal[run->n_rels] = uni;
                ar_append_relation(run,&rel,cand+i,trim+i,n);
                st->accepted++;
                st->corr_total += n;
                if (rel.flags & ASM_REL_PARITY) st->parity_flipped++;
                rms_list[nrms++] = rel.rms;
            } else if (gr == 4) {
                /* held back, not pushed here: a LEDGER must not change a decision, and a relation
                 * appended between two real ones shifts every later relation index, which moves
                 * clean2's drop order (its flag-mass ties break on the index) and with it the whole
                 * layout -- measured 2026-09-17 on pherc343, clean2 dropped 362 relations instead of
                 * 452 and placement ran 3,897 sweeps instead of 9,781.  They go after the real ones. */
                if (n_contact == cap_contact) {
                    size_t ncap = cap_contact ? cap_contact*2 : 256;
                    AsmRelation *grown = ARENA_ALLOC(run->arena, ncap*sizeof(AsmRelation));
                    if (n_contact) memcpy(grown, contact_rows, n_contact*sizeof(AsmRelation));
                    contact_rows = grown; cap_contact = ncap;
                }
                contact_rows[n_contact++] = rel;
            }
        }
        i = j;
    }
    if (ar_diag_fp) { fclose(ar_diag_fp); ar_diag_fp = NULL; }
    for (size_t k = 0; k < n_contact; k++) { AsmRun_push_rel(run, &contact_rows[k]); st->contact_ledgered++; }
    st->bridge_demoted = ar_bridge_premium(run, unimodal, ASM_SEAM_BRIDGE_UNIMODAL, ASM_SEAM_BRIDGE_LEAF_AREA, work[0],
                                           &st->bridge_corroborated, &st->bridge_leaves);
    if (st->bridge_leaves)
        fprintf(stderr, "[assemble relate] bridge premium: %zu marginal seams kept as joins because one block is below %.0f vox^2 (a leaf, not a weld)\n",
                st->bridge_leaves, (double)ASM_SEAM_BRIDGE_LEAF_AREA);
    if (st->bridge_demoted) {
        st->accepted -= st->bridge_demoted;
        n_weak += st->bridge_demoted;
        fprintf(stderr, "[assemble relate] bridge premium (unimodality >= %.2f to join two otherwise-disconnected blocks): %zu load-bearing seams kept as placement evidence%s\n",
                ASM_SEAM_BRIDGE_UNIMODAL, st->bridge_demoted,
                ASM_SEAM_BRIDGE_CORROBORATED ? "" : " (corroboration test OFF: a bridge of the partial pass, not of the accepted graph)");
        if (st->bridge_corroborated)
            fprintf(stderr, "[assemble relate] bridge premium: %zu marginal seams kept as joins because another seam independently connects the same two blocks\n",
                    st->bridge_corroborated);
    }
    st->weak_kept = n_weak;
    for (size_t r = 0; r < run->n_rels; r++)
        if (run->rels[r].flags & ASM_REL_CROSSWRAP) st->wrap_anchors++;
    if (st->wrap_anchors || st->wrap_on_iso)
        fprintf(stderr, "[assemble relate] wrap gate: %zu relations carry a sheet switch (pairwise rotations agreeing with the fit under 80%%) "
                        "and are %s; separately, %zu of the %zu isometry rejections are bimodal too, which the test only sees now that "
                        "an isometry rejection no longer skips it\n",
                st->wrap_anchors,
                ASM_WRAP_VETO_ANCHOR ? "refused as placement anchors, layout joins and discovery hypotheses"
                                  : "still placement anchors (ASM_WRAP_VETO_ANCHOR 0)",
                st->wrap_on_iso, st->rej_iso);
    if (st->contact_ledgered)
        fprintf(stderr, "[assemble relate] tangency ledger: %zu chart pairs meet face to face across a cube face (ASM_REL_CONTACT rows without correspondences: a veto for layout joins, never an edge)\n", st->contact_ledgered);
    fprintf(stderr, "[assemble relate] kept for placement only (gate-rejected, rigid fit rms <= %.0f): weak %zu, short seams %zu (%.0f-%.0f vox, >= %d correspondences; short seams not kept: corr < %d %zu, rms over %zu)\n",
            AR_RIGID_MAX_RMS, n_weak, st->short_kept, AR_SHORT_MIN_LEN, o->seam_min_len, AR_SHORT_MIN_CORR, AR_SHORT_MIN_CORR, st->short_rej_corr, st->short_rej_rms);
    if (nrms) {
        qsort(rms_list, nrms, sizeof(double), ar_cmp_double);
        st->rms_p50 = rms_list[nrms/2];
        st->rms_p95 = rms_list[(size_t)((double)(nrms-1) * 0.95)];
    }
    st->seam_sec = ves_clock_sec() - t0;

    /* layer probes */
    double t1 = ves_clock_sec();
    int ci = 0, ncubes_i = (int)run->n_cubes;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
#endif
    for (ci = 0; ci < ncubes_i; ci++) {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        ar_layer_cube(run, &ix, (size_t)ci, o, work[tid], &bufs[tid]);
    }
    size_t nhit = 0;
    for (int t = 0; t < threads; t++) nhit += bufs[t].n_hit;
    ArHit *hits = ARENA_ALLOC(run->arena, (nhit ? nhit : 1) * sizeof(ArHit));
    {
        size_t q = 0;
        for (int t = 0; t < threads; t++) { memcpy(&hits[q], bufs[t].hit, bufs[t].n_hit * sizeof(ArHit)); q += bufs[t].n_hit; }
    }
    qsort(hits, nhit, sizeof(ArHit), ar_cmp_hit);
    run->chart_layer_d = ARENA_CALLOC(run->arena, run->n_charts, sizeof(double));
    double *chart_acc = ARENA_CALLOC(run->arena, run->n_charts, sizeof(double));
    double *chart_cnt = ARENA_CALLOC(run->arena, run->n_charts, sizeof(double));
    for (size_t i = 0; i < nhit;) {
        size_t j = i;
        while (j < nhit && hits[j].chart == hits[i].chart && hits[j].other == hits[i].other && hits[j].side == hits[i].side) j++;
        size_t n = j - i;
        if (n >= 3) {
            AsmLayerPair lp;
            lp.a = hits[i].chart; lp.b = hits[i].other; lp.side = hits[i].side; lp.count = (int32_t)n;
            lp.d_median = hits[i + n/2].dist;   /* hits are sorted by dist within the group */
            {
                int cnt[4] = { 0, 0, 0, 0 };
                for (size_t k2 = i; k2 < j; k2++) { int od = hits[k2].order; if (od < 1) od = 1; if (od > 3) od = 3; cnt[od]++; }
                size_t half = n / 2, acc = 0; int med = 1;
                for (int od = 1; od <= 3; od++) { acc += (size_t)cnt[od]; if (acc > half) { med = od; break; } }
                lp.k = (int8_t)med;
                long sgn = 0;
                for (size_t k2 = i; k2 < j; k2++) sgn += hits[k2].nsign;
                lp.nsign = (int8_t)(sgn > 0 ? 1 : sgn < 0 ? -1 : 0);
            }
            /* keep up to 16 representative hits around the median distance */
            size_t keep = n < 16 ? n : 16;
            size_t start = n > keep ? (n - keep) / 2 : 0;
            lp.hit_first = (int32_t)run->n_layer_hits;
            lp.hit_count = (int32_t)keep;
            for (size_t k = 0; k < keep; k++) {
                const ArHit *h = &hits[i + start + k];
                AsmLayerHit lh; lh.a = h->chart; lh.va = h->va; lh.b = h->other; lh.fb = h->fb;
                lh.l0 = h->l0; lh.l1 = h->l1; lh.dist = h->dist; lh.side = h->side; lh.order = h->order; lh.nsign = h->nsign;
                AsmRun_push_layer_hit(run, &lh);
            }
            AsmRun_push_layer(run, &lp);
            st->layer_pairs++;
            /* the chart's layer distance is its NEAREST layer neighbour (min over
             * pairs), never a mean that mixes first- and second-layer hits */
            if (chart_cnt[(size_t)lp.a] == 0.0 || lp.d_median < chart_acc[(size_t)lp.a]) chart_acc[(size_t)lp.a] = lp.d_median;
            chart_cnt[(size_t)lp.a] += (double)n;
        }
        i = j;
    }
    double *all_d = ARENA_ALLOC(run->arena, (run->n_charts ? run->n_charts : 1) * sizeof(double));
    size_t nall = 0;
    for (size_t c = 0; c < run->n_charts; c++) {
        if (chart_cnt[c] > 0.0) {
            run->chart_layer_d[c] = chart_acc[c];
            all_d[nall++] = run->chart_layer_d[c];
            st->charts_with_layer++;
        }
    }
    if (nall) { qsort(all_d, nall, sizeof(double), ar_cmp_double); run->layer_d_global = all_d[nall/2]; }
    else run->layer_d_global = ASM_LAYER_PROBE_VOX * 0.4;
    for (size_t h = 0; h < nhit; h++) {
        int bin = (int)(hits[h].dist / 2.0f);
        if (bin < 0) bin = 0;
        if (bin >= ASM_RELATE_HIST) bin = ASM_RELATE_HIST - 1;
        st->layer_hist[bin]++;
    }
    st->layer_hits = nhit;
    st->layer_d_global = run->layer_d_global;
    st->layer_sec = ves_clock_sec() - t1;

    for (int t = 0; t < threads; t++) { Arena_dispose(&scratch[t]); Arena_dispose(&work[t]); }
    return 0;
}

/* ---- selftest -------------------------------------------------------------------- */

static int ar_trim_selftest(void)
{
    enum { N = 17, NV = 2*N, NF = 2*(N-1) };
    int fails = 0;
    /* Two ragged boundaries: the gap alternates 2 and 6 vox (raw rigid rms 2.0, within the seam
     * bar) and 2 and 8 (rms 3.0, beyond it).  BOTH stay placement evidence rather than becoming
     * joins, because their pairwise rotations disagree and the unimodality statistic is the
     * cross-wrap detector: qualifying it on the rigid fit -- letting the 2/6 case through on the
     * argument that a 2-vox fit cannot hide a switch -- was measured DEAD on three piles
     * (2026-09-17; see ASM_SEAM_UNIMODAL_RMS), because a pair of fused wraps fits one rigid
     * transform very well.  What this control is actually about is the trim transport: the
     * relation must be PUBLISHED either way, must carry its RAW residual (the transport recovers
     * the frame-to-frame transform and must not launder the seam's own consistency), and must be
     * reproducible from the published relation and correspondence tables alone. */
    for (int ragged = 0; ragged < 2; ragged++)
    for (int world = 0; world < 2; world++) for (int parity = 0; parity < 2; parity++) for (int reverse = 0; reverse < 2; reverse++) {
        const int amp = ragged ? 6 : 4;
        Arena_T arena = Arena_new();
        float xyz[2][NV*3], uv[2][NV*2], original_xyz[2][NV*3], original_uv[2][NV*2];
        int32_t faces[2][NF*3], original_faces[2][NF*3]; uint8_t boundary[NV];
        AsmChart charts[2] = {{0}}; ArCand cand[N]; ArTrim trim[N];
        double theta = 0.43, ct = cos(theta), sn = sin(theta);
        for (int side = 0; side < 2; side++) {
            for (int j = 0; j < N; j++) for (int col = 0; col < 2; col++) {
                int v = 2*j+col;
                double x = side ? 22+amp*(j%2)+20*col : 20*col, y = 8*j;
                xyz[side][3*v] = (float)(world ? 10000+cos(.37)*x-sin(.37)*cos(.23)*y : x);
                xyz[side][3*v+1] = (float)(world ? 3500+sin(.37)*x+cos(.37)*cos(.23)*y : y);
                xyz[side][3*v+2] = (float)(world ? 4200+sin(.23)*y : 0);
                double u = side ? ct*(x-17)+sn*(y+9) : x;
                double w = side ? -sn*(x-17)+ct*(y+9) : y;
                uv[side][2*v] = (float)(side && parity ? -u : u); uv[side][2*v+1] = (float)w;
                boundary[v] = 1;
            }
            for (int j = 0; j < N-1; j++) {
                int v = 2*j; int32_t f[6] = {v,v+1,v+3,v,v+3,v+2};
                if (side && parity) for (int k = 0; k < 2; k++) { int32_t t = f[3*k+1]; f[3*k+1] = f[3*k+2]; f[3*k+2] = t; }
                memcpy(faces[side]+6*j,f,sizeof f);
            }
            charts[side].id = side; charts[side].nv = NV; charts[side].nf = NF; charts[side].area3d = 2560;
            charts[side].xyz = xyz[side]; charts[side].uv = uv[side]; charts[side].faces = faces[side]; charts[side].boundary = boundary;
        }
        for (int j = 0; j < N; j++) {
            cand[j] = (ArCand){reverse ? 1 : 0,reverse ? 0 : 1,reverse ? 2*j : 2*j+1,reverse ? 2*j+1 : 2*j,
                              (float)(2+amp*(j%2)),0,parity ? -1 : 1};
        }
        memcpy(original_xyz,xyz,sizeof xyz); memcpy(original_uv,uv,sizeof uv); memcpy(original_faces,faces,sizeof faces);
        AsmRun run = {0}; run.arena = arena; run.charts = charts; run.n_charts = 2;
        AsmRelateOpts opts; AsmRelate_default_opts(&opts); AsmRelateStats stats = {0};
        AsmRelation raw = {0}, corrected = {0}, fallback = {0};
        int raw_result = ar_group_relation(&run,cand,N,NULL,&opts,arena,&raw,&stats,NULL);
        int prepared = ar_candidate_trims(&run,cand,N,arena,trim);
        int corrected_result = prepared ? 0 : ar_group_relation(&run,cand,N,trim,&opts,arena,&corrected,&stats,NULL);
        /* CONTRACT (2026-09-15).  The trim transport exists to recover the
         * frame-to-frame TRANSFORM across a trim gap.  It must never launder
         * the seam's own consistency: gap_a/gap_b are the pair displacement
         * b-a re-expressed in each chart's frame, so transporting the midpoints
         * before the pairwise gates registers every correspondence onto its own
         * partner -- one free parameter per observation -- and a sheet switch
         * inside the seam stops being visible.  This control's boundary is
         * deliberately ragged (its gap alternates 2 and 6 vox over an 8-vox
         * spacing), so the relation must still be published carrying its RAW
         * residual, and ragged evidence alone must not promote it to a join. */
        int ok = raw_result != 1 && raw_result != 0 && corrected_result != 0;
        ok = ok && !prepared && corrected.rms > 1.0 && corrected.n_corr == N;
        if (corrected_result != 0) {
            /* Reconstruct the seam midpoint observations using only the
             * published relation/correspondence tables, as stage 3 does.
             * They must describe the same fit after temporary trims expire. */
            ar_append_relation(&run,&corrected,cand,trim,N);
            const AsmRelation *r=&run.rels[0]; const AsmChart *a=&charts[r->a],*b=&charts[r->b];
            double qa[2*N],qb[2*N]; int transported=1;
            for (int j=0; j<N; j++) if (run.corr[j].valid!=3) transported=0;
            for (int j=0; j<N; j++) {
                const AsmCorr *c=&run.corr[j];
                for (int d=0; d<2; d++) {
                    qa[2*j+d]=a->uv[2*c->va+d]+(transported ? .5*c->gap_a[d] : 0);
                    qb[2*j+d]=b->uv[2*c->vb+d]-(transported ? .5*c->gap_b[d] : 0);
                }
                if (r->flags & ASM_REL_PARITY) qb[2*j]=-qb[2*j];
            }
            double th,tx,ty,rms;
            int fitted=!AsmRelate_rigid_fit(qa,qb,NULL,N,&th,&tx,&ty,&rms);
            if (!fitted || rms>.002 || hypot(tx-r->tx,ty-r->ty)>.002 || fabs(th-r->theta)>1e-4) {
                fprintf(stderr,"  asm_relate selftest FAIL: published trim fit world %d parity %d reverse %d, RMS %.9g translation %.9g\n",
                        world,parity,reverse,rms,hypot(tx-r->tx,ty-r->ty));
                ok=0;
            }
        }
        /* An unavailable tangent cannot remove a difficult observation. */
        trim[N/2].valid = 1;
        int fallback_result = ar_group_relation(&run,cand,N,trim,&opts,arena,&fallback,&stats,NULL);
        ok = ok && fallback_result == raw_result && !memcmp(&raw,&fallback,sizeof raw);
        ar_candidate_trims(&run,cand,N,arena,trim);
        for (int j = 0; j < N; j++) cand[j].ngap = 100;
        {
            /* the tangency rejection leaves a CONTACT ledger row (ASM_RELATE_LEDGER_CONTACT):
             * a veto for layout joins, never a join, no correspondences */
            int contact = ar_group_relation(&run,cand,N,trim,&opts,arena,&fallback,&stats,NULL);
            if (ASM_RELATE_LEDGER_CONTACT)
                ok = ok && contact == 4 && (fallback.flags & ASM_REL_CONTACT) && !AsmRel_is_join(&fallback) &&
                     fallback.corr_count == 0 && fallback.corr_first < 0 && fallback.normal_gap == 100.0 && fallback.n_corr == N;
            else ok = ok && contact == 0;
        }
        for (int j = 0; j < N; j++) cand[j].ngap = 0;
        ok = ok && !memcmp(xyz,original_xyz,sizeof xyz) && !memcmp(uv,original_uv,sizeof uv) && !memcmp(faces,original_faces,sizeof faces);
        /* A nonisometric source parameterization still fails the same gates. */
        for (int v = NV/2; v < NV; v++) uv[1][2*v] += 24;
        ar_candidate_trims(&run,cand,N,arena,trim);
        ok = ok && ar_group_relation(&run,cand,N,trim,&opts,arena,&fallback,&stats,NULL) != 1;
        /* Invalid queries are rejected before any output or source mutation. */
        int32_t bad_vertex = NV; double delta[3] = {1,0,0}; float gap[2] = {7,8}; uint8_t valid = 9;
        ok = ok && AsmContinuity_chart_gaps(arena,&charts[0],&bad_vertex,delta,1,gap,&valid) == -1 && gap[0] == 7 && gap[1] == 8 && valid == 9;
        fprintf(stderr,"  asm_relate trim control: ragged %d world %d parity %d reverse %d, raw %d -> transported %d, RMS %.9g: %s\n",
                ragged,world,parity,reverse,raw_result,corrected_result,corrected.rms,ok ? "PASS" : "FAIL");
        fails += !ok; Arena_dispose(&arena);
    }
    return fails;
}

static int ar_layer_search_selftest(void)
{
    int fails = 0;
    /* The actual layer probe must hit the corner of this large triangle.
     * Its centroid lies over nine voxels from every source probe: the old
     * fixed six-voxel broad phase silently omitted the surface. */
    for (int translated = 0; translated < 2; translated++) {
        Arena_T arena = Arena_new(), scratch = Arena_new();
        AsmChart charts[2] = {{0}}; AsmRun run = {0};
        float xyz[2][9] = {{0,0,0, 0,.25f,0, 0,0,.25f},
                            {5,0,0, 5,20,0, 5,0,20}};
        float normals[9] = {1,0,0, 1,0,0, 1,0,0}, uv[6] = {0,0, 1,0, 0,1};
        int32_t triangle[3] = {0,1,2}, offsets[2] = {0,2}, ids[2] = {0,1}, neighbors[26];
        float shift[3] = {5000,-7000,3000};
        for (int c = 0; c < 2; c++) {
            for (int v = 0; v < 3; v++) for (int d = 0; d < 3; d++) xyz[c][3*v+d] += translated*shift[d];
            charts[c].id = c; charts[c].nv = 3; charts[c].nf = 1;
            charts[c].xyz = xyz[c]; charts[c].nrm = normals; charts[c].uv = uv; charts[c].faces = triangle;
        }
        float original[2][9]; memcpy(original,xyz,sizeof xyz);
        for (int i = 0; i < 26; i++) neighbors[i] = -1;
        run.arena = arena; run.charts = charts; run.n_charts = 2; run.n_cubes = 1;
        ArCubeIndex index = {1,offsets,ids,neighbors}; ArBuf buffer = {0}; buffer.arena = arena;
        AsmRelateOpts opts; AsmRelate_default_opts(&opts); opts.layer_probe_vox = 6;
        ar_layer_cube(&run,&index,0,&opts,scratch,&buffer);
        int found = 0;
        for (size_t i = 0; i < buffer.n_hit; i++) {
            const ArHit *h = buffer.hit+i;
            if (h->chart == 0 && h->va == 0 && h->side == 1 && h->other == 1 && h->fb == 0 &&
                h->order == 1 && h->dist == 4 && h->l0 == 1 && h->l1 == 0) found++;
        }
        if (found != 1 || memcmp(original,xyz,sizeof xyz)) {
            fprintf(stderr,"  asm_relate selftest FAIL: large-face probe, translated %d, hits %d\n",translated,found);
            fails++;
        }
        Arena_dispose(&scratch); Arena_dispose(&arena);
    }
    /* A crowded neighborhood cannot truncate the broad phase. Compare the
     * complete result set with a direct distance test, including a replay
     * after expansion and a query that needs no expansion. */
    {
        enum { N = 4097 };
        Arena_T arena = Arena_new();
        float *points = ARENA_ALLOC(arena,3*N*sizeof(float));
        for (int i = 0; i < N; i++) {
            points[3*i] = (float)((i%31)-15)*.07f;
            points[3*i+1] = (float)(((i/31)%29)-14)*.09f;
            points[3*i+2] = (float)(i/899)*.11f;
        }
        KDTree_T tree = KDTree_new(arena,points,N);
        size_t capacity = AR_PROBE_INITIAL_CAP;
        int32_t *ball = ARENA_ALLOC(arena,capacity*sizeof(int32_t));
        unsigned char seen[N];
        float center[3] = {0,0,0}, radius = ar_probe_radius_sq(0);
        int ok = 1;
        for (int query = 0; query < 3; query++) {
            if (query == 2) center[0] = 100;
            size_t count = ar_probe_ball(arena,tree,N,center,radius,&ball,&capacity), expected = 0;
            memset(seen,0,sizeof seen);
            for (size_t i = 0; i < count; i++) {
                if (ball[i] < 0 || ball[i] >= N || seen[ball[i]]) ok = 0;
                else seen[ball[i]] = 1;
            }
            for (int i = 0; i < N; i++) {
                double squared = 0;
                for (int d = 0; d < 3; d++) { double delta = (double)points[3*i+d]-center[d]; squared += delta*delta; }
                int within = squared <= radius; expected += within;
                if (seen[i] != within) ok = 0;
            }
            if (count != expected || capacity != N) ok = 0;
        }
        if (!ok) { fprintf(stderr,"  asm_relate selftest FAIL: incomplete dense layer neighborhood\n"); fails++; }
        Arena_dispose(&arena);
    }
    fprintf(stderr,"  native complete layer search: %s (%d failures)\n",fails ? "FAIL" : "ok",fails);
    return fails;
}

/* The unimodality statistic is the cross-wrap detector, so it must be MEASURED
 * for every pair that gets a rigid fit and its verdict must survive as a flag.
 * Four seams of 16 correspondences, 45 uv-vox long (over the 30-vox length
 * gate), B's frame built from A's:
 *   0  one rigid motion                  -> accepted, no ASM_REL_CROSSWRAP
 *   1  the top half rotated about the
 *      seam's midpoint (a sheet switch)  -> rejected by unimodality, CROSSWRAP,
 *                                           and still WEAK placement evidence
 *   2  a uniform 1.20 scale              -> rejected by ISOMETRY only: the
 *                                           angles agree, so no CROSSWRAP, and
 *                                           the statistic is measured anyway
 *                                           (before 2026-09-17 an isometry
 *                                           rejection skipped the test)
 *   3  both                              -> the isometry gate owns the outcome,
 *                                           the wrap verdict is its own flag */
static int ar_wrap_flag_selftest(void)
{
    int fails = 0;
    Arena_T ta = Arena_new();
    AsmRelateOpts o; AsmRelate_default_opts(&o);
    for (int cs = 0; cs < 4; cs++) {
        const size_t n = 16;
        const double du = 3.0, mid = du * (double)(n - 1) / 2.0;
        const double scale = (cs == 2 || cs == 3) ? 1.20 : 1.0;
        const double phi = (cs == 1 || cs == 3) ? 0.3 : 0.0;
        const double thr = 0.3, trx = 5.0, tr_y = -2.0, cr = cos(thr), sr = sin(thr);
        AsmChart ch[2]; memset(ch, 0, sizeof ch);
        float xyz[2][16*3], uv[2][16*2];
        for (size_t i = 0; i < n; i++) {
            for (int k = 0; k < 2; k++) {
                xyz[k][i*3] = (float)(4.0 * (double)i);
                xyz[k][i*3+1] = (float)(k == 1 ? 1.0 : 0.0);
                xyz[k][i*3+2] = 0.0f;
            }
            double ua = du * (double)i, va = 0.0;
            uv[0][i*2] = (float)ua; uv[0][i*2+1] = (float)va;
            /* the sheet switch: the far half turns about the seam's midpoint,
             * which leaves every distance within a half exact and the longest
             * cross-half distance 1% short -- invisible to the isometry gate */
            double bu = ua, bv = va;
            if (i >= n/2 && phi != 0.0) {
                double ox = ua - mid, oy = va;
                bu = mid + cos(phi)*ox - sin(phi)*oy;
                bv =       sin(phi)*ox + cos(phi)*oy;
            }
            bu *= scale; bv *= scale;
            uv[1][i*2] = (float)(cr * bu - sr * bv + trx);
            uv[1][i*2+1] = (float)(sr * bu + cr * bv + tr_y);
        }
        for (int k = 0; k < 2; k++) { ch[k].id = k; ch[k].nv = n; ch[k].xyz = xyz[k]; ch[k].uv = uv[k]; ch[k].area3d = 1000.0; }
        AsmRun run; memset(&run, 0, sizeof run); run.charts = ch; run.n_charts = 2;
        ArCand g[16];
        for (size_t i = 0; i < n; i++) { g[i].chart_a = 0; g[i].chart_b = 1; g[i].va = (int32_t)i; g[i].vb = (int32_t)i; g[i].dist = 1.0f; g[i].ngap = 0.0f; g[i].sign = 1; }
        AsmRelation rel = {0}; AsmRelateStats st; memset(&st, 0, sizeof st);
        double uni = -2.0;
        int gr = ar_group_relation(&run, g, n, NULL, &o, ta, &rel, &st, &uni);
        int cw = (rel.flags & ASM_REL_CROSSWRAP) != 0;
        int wk = (rel.flags & ASM_REL_WEAK) != 0;
        int want_cw = (cs == 1 || cs == 3), want_wk = (cs != 0);
        int want_iso = (cs == 2 || cs == 3), want_uni = (cs == 1);
        int ok = gr == (want_wk ? 2 : 1) && cw == want_cw && wk == want_wk && uni >= 0.0 &&
                 st.rej_iso == (size_t)want_iso && st.rej_unimodal == (size_t)want_uni &&
                 st.wrap_on_iso == (size_t)(cs == 3);
        if (!ok) {
            fprintf(stderr, "  asm_relate wrap-flag FAIL case %d: return %d flags %u unimodal %.3f "
                            "(iso %zu unimodal %zu wrap-on-iso %zu; want crosswrap %d weak %d)\n",
                    cs, gr, rel.flags, uni, st.rej_iso, st.rej_unimodal, st.wrap_on_iso, want_cw, want_wk);
            fails++;
        }
    }
    /* the veto set and the switch must agree, so an armed build cannot join a flagged pair */
    if (((ASM_JOIN_VETO_FLAGS & ASM_REL_CROSSWRAP) != 0u) != (ASM_WRAP_VETO_JOIN != 0)) {
        fprintf(stderr, "  asm_relate wrap-flag FAIL: ASM_JOIN_VETO_FLAGS does not follow ASM_WRAP_VETO_JOIN\n");
        fails++;
    }
    Arena_dispose(&ta);
    return fails;
}

/* The cut-edge test the bridge premium's corroboration rests on.  Six charts:
 * a triangle 0-1-2 (no cut edges), a pendant edge 2-3 (a cut edge), a second
 * triangle 3-4-5, and a PARALLEL pair of relations between 1 and 2 -- two
 * independent seams joining one pair of charts corroborate each other, so
 * neither is a cut edge even when the triangle is removed. */
static int ar_cut_edge_selftest(void)
{
    int fails = 0;
    Arena_T ta = Arena_new();
    const int ea[8] = { 0,1,2, 2, 3,4,5, 1 };
    const int eb[8] = { 1,2,0, 3, 4,5,3, 2 };
    for (int cs = 0; cs < 2; cs++) {
        size_t ne = cs == 0 ? 7u : 8u;         /* case 1 adds the parallel 1-2 relation */
        size_t use = cs == 0 ? 7u : 8u;
        AsmRelation rels[8]; memset(rels, 0, sizeof rels);
        ArBridge ord[8];
        for (size_t k = 0; k < use; k++) {
            rels[k].a = (int32_t)ea[k]; rels[k].b = (int32_t)eb[k];
            ord[k].uni = 1.0; ord[k].corr = 8; ord[k].rel = (int32_t)k;
        }
        AsmRun run; memset(&run, 0, sizeof run);
        run.rels = rels; run.n_rels = ne; run.n_charts = 6;
        uint8_t cut[8]; memset(cut, 0, sizeof cut);
        ar_cut_edges(&run, ord, ne, cut, ta);
        /* edge 3 (2-3) is the only cut edge in both cases */
        for (size_t k = 0; k < ne; k++) {
            int want = (k == 3);
            if (cut[k] != (uint8_t)want) {
                fprintf(stderr, "  asm_relate cut-edge FAIL case %d: edge %zu (%d-%d) cut %d, want %d\n",
                        cs, k, ea[k], eb[k], cut[k], want);
                fails++;
            }
        }
    }
    /* a bare path 0-1-2: every edge is a cut edge */
    {
        AsmRelation rels[2]; memset(rels, 0, sizeof rels);
        ArBridge ord[2];
        for (int k = 0; k < 2; k++) { rels[k].a = (int32_t)k; rels[k].b = (int32_t)(k+1); ord[k].uni = 1.0; ord[k].corr = 8; ord[k].rel = k; }
        AsmRun run; memset(&run, 0, sizeof run); run.rels = rels; run.n_rels = 2; run.n_charts = 3;
        uint8_t cut[2] = { 0, 0 };
        ar_cut_edges(&run, ord, 2, cut, ta);
        if (!cut[0] || !cut[1]) { fprintf(stderr, "  asm_relate cut-edge FAIL: a path's edges must all be cut edges\n"); fails++; }
    }
    Arena_dispose(&ta);
    return fails;
}

/* The leaf exemption (ASM_SEAM_BRIDGE_LEAF_AREA).  Blocks X = {0,1,2} and
 * Y = {3,4,5} (30,000 vox^2 each, joined inside by unimodality 1.0), chart 6
 * (5,000 vox^2) and chart 7 (5,000 vox^2).  Marginal seams (0.90, 0.88, 0.86):
 * 2-6 attaches the leaf 6 to X; 6-3 would then join X+6 to Y; 5-7 attaches the
 * leaf 7 to Y.  With a 20,000 leaf area both leaf seams stay joins and the
 * chain through the leaf is still demoted, because X+6 is no longer a leaf;
 * with 0 all three marginal seams are demoted, the premium as before. */
static int ar_bridge_leaf_selftest(void)
{
    int fails = 0;
    Arena_T ta = Arena_new();
    const int ea[7] = { 0,1, 3,4, 2,6,5 }, eb[7] = { 1,2, 4,5, 6,3,7 };
    const double ua[7] = { 1,1, 1,1, 0.90,0.88,0.86 };
    for (int cs = 0; cs < 2; cs++) {
        AsmChart charts[8]; memset(charts, 0, sizeof charts);
        for (int c = 0; c < 8; c++) { charts[c].id = c; charts[c].area3d = c < 6 ? 10000.0 : 5000.0; }
        AsmRelation rels[7]; memset(rels, 0, sizeof rels);
        double uni[7];
        for (int k = 0; k < 7; k++) { rels[k].a = ea[k]; rels[k].b = eb[k]; rels[k].n_corr = 16; uni[k] = ua[k]; }
        AsmRun run; memset(&run, 0, sizeof run);
        run.charts = charts; run.n_charts = 8; run.rels = rels; run.n_rels = 7;
        size_t spared = 0, leaves = 0;
        size_t demoted = ar_bridge_premium(&run, uni, 0.95, cs ? 20000.0 : 0.0, ta, &spared, &leaves);
        int weak[7]; for (int k = 0; k < 7; k++) weak[k] = (rels[k].flags & ASM_REL_WEAK) != 0;
        int ok = cs ? demoted == 1 && leaves == 2 && !weak[4] && weak[5] && !weak[6]
                    : demoted == 3 && leaves == 0 && weak[4] && weak[5] && weak[6];
        for (int k = 0; k < 4; k++) ok = ok && !weak[k];
        if (!ok) {
            fprintf(stderr, "  asm_relate bridge-leaf FAIL case %d: demoted %zu leaves %zu weak %d%d%d\n",
                    cs, demoted, leaves, weak[4], weak[5], weak[6]);
            fails++;
        }
    }
    Arena_dispose(&ta);
    return fails;
}

int AsmRelate_selftest(void)
{
    int fails = ar_trim_selftest() + ar_layer_search_selftest() + ar_wrap_flag_selftest() + ar_cut_edge_selftest() +
                ar_bridge_leaf_selftest();
    /* rigid fit recovers a known transform */
    double qb[16] = { 0,0, 10,0, 10,10, 0,10, 5,5, 3,7, 8,2, 1,9 };
    double qa[16];
    double th = 0.4, tx = 12.0, ty = -3.0, c = cos(th), s = sin(th);
    for (int i = 0; i < 8; i++) {
        qa[i*2] = c*qb[i*2] - s*qb[i*2+1] + tx;
        qa[i*2+1] = s*qb[i*2] + c*qb[i*2+1] + ty;
    }
    double th2, tx2, ty2, rms;
    AsmRelate_rigid_fit(qa, qb, NULL, 8, &th2, &tx2, &ty2, &rms);
    if (fabs(th2 - th) > 1e-9 || fabs(tx2 - tx) > 1e-9 || fabs(ty2 - ty) > 1e-9 || rms > 1e-9) {
        fprintf(stderr, "  asm_relate selftest FAIL: rigid fit %g %g %g %g\n", th2, tx2, ty2, rms); fails++;
    }
    /* mirrored input cannot be fitted rigidly: rms must be large */
    double qm[16];
    for (int i = 0; i < 8; i++) { qm[i*2] = -qb[i*2]; qm[i*2+1] = qb[i*2+1]; }
    AsmRelate_rigid_fit(qa, qm, NULL, 8, &th2, &tx2, &ty2, &rms);
    if (rms < 1.0) { fprintf(stderr, "  asm_relate selftest FAIL: mirror rms %g\n", rms); fails++; }
    /* SHORT seams (2026-09-10): two charts sharing a seam of 8 correspondences over 16 vox (under the
     * 30-vox gate), B's uv a rigid transform of A's: kept as SHORT with its correspondences; 4
     * correspondences: rejected (too few); 12-vox noise on half the points: rejected (rms over 4); the
     * same seam over 63 vox: accepted, SHORT unset */
    {
        Arena_T ta = Arena_new();
        AsmRelateOpts o; AsmRelate_default_opts(&o);
        for (int cs = 0; cs < 7; cs++) {
            size_t n = cs == 1 ? 4 : 8;
            double step = cs >= 3 ? 8.0 : (cs == 1 ? 5.0 : 2.0);   /* 4 points over 15 vox: 19 vox long, too few correspondences */
            AsmChart ch[2]; memset(ch, 0, sizeof ch);
            float xyz[2][8*3], uv[2][8*2];
            double thr = 0.3, trx = 5.0, tr_y = -2.0, cr = cos(thr), sr = sin(thr);
            for (size_t i = 0; i < n; i++) {
                double z = step * (double)i;
                for (int k = 0; k < 2; k++) { xyz[k][i*3] = (float)z; xyz[k][i*3+1] = (float)(k == 1 ? 1.0 : 0.0); xyz[k][i*3+2] = 0.0f; }
                if (cs >= 4) {
                    /* These equally long seams are perpendicular to the
                     * old power iteration's (1,1,1) starting direction. */
                    int p = cs == 6 ? 1 : 0, q = cs == 4 ? 1 : 2;
                    for (int k = 0; k < 2; k++) {
                        xyz[k][i*3] = xyz[k][i*3+1] = xyz[k][i*3+2] = 0;
                        xyz[k][i*3+p] = (float)(z/sqrt(2.0));
                        xyz[k][i*3+q] = -xyz[k][i*3+p];
                        xyz[k][i*3+3-p-q] = (float)k;
                    }
                }
                double ua = z, va = 0.0;
                uv[0][i*2] = (float)ua; uv[0][i*2+1] = (float)va;
                double noise = (cs == 2 && (i & 1)) ? 12.0 : 0.0;
                uv[1][i*2] = (float)(cr * ua - sr * va + trx); uv[1][i*2+1] = (float)(sr * ua + cr * va + tr_y + noise);
            }
            for (int k = 0; k < 2; k++) { ch[k].id = k; ch[k].nv = n; ch[k].xyz = xyz[k]; ch[k].uv = uv[k]; ch[k].area3d = 1000.0; }
            AsmRun run; memset(&run, 0, sizeof run); run.charts = ch; run.n_charts = 2;
            ArCand g[8];
            for (size_t i = 0; i < n; i++) { g[i].chart_a = 0; g[i].chart_b = 1; g[i].va = (int32_t)i; g[i].vb = (int32_t)i; g[i].dist = 1.0f; g[i].ngap = 0.0f; g[i].sign = 1; }
            AsmRelation rel = {0}; AsmRelateStats st; memset(&st, 0, sizeof st);
            int gr = ar_group_relation(&run, g, n, NULL, &o, ta, &rel, &st, NULL);
            int want = cs == 0 ? (AR_SHORT_SEAMS ? 3 : 0) : (cs >= 3 ? 1 : 0);   /* an A/B build with the lever off rejects the short seam */
            int flag_ok = cs == 0 ? (!AR_SHORT_SEAMS || ((rel.flags & ASM_REL_SHORT) != 0 && rel.n_corr == 8)) : (cs >= 3 ? (rel.flags & ASM_REL_SHORT) == 0 : 1);
            int cnt_ok = !AR_SHORT_SEAMS ? 1 : (cs == 1 ? st.short_rej_corr == 1 : (cs == 2 ? st.short_rej_rms == 1 : 1));
            int length_ok = cs < 3 || (gr == 1 && fabs(rel.seam_len-step*sqrt((double)(n*n-1))) < 1e-4);
            if (gr != want || !flag_ok || !cnt_ok || !length_ok) {
                fprintf(stderr, "  asm_relate selftest FAIL: short seam case %d: return %d (want %d), flags %u, rej corr %zu rms %zu\n", cs, gr, want, cs == 0 || cs == 3 ? rel.flags : 0u, st.short_rej_corr, st.short_rej_rms);
                fails++;
            }
        }
        Arena_dispose(&ta);
    }
    if (fails == 0) fprintf(stderr, "  asm_relate selftest: all passed\n");
    return fails;
}
