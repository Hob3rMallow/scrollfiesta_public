/* closeup_adopt.c -- see closeup_adopt.h. */
#include "closeup_adopt.h"
#include "../common/arena.h"
#include "../common/kdtree.h"
#include "../common/pipeline_constants.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/ves_omp.h"

#define CA_MAX_SAMPLES 6000     /* vertices tried per chart (strided) */
#define CA_LAYER_SAMPLES 600    /* ...of which the layer gate tests at most this many (its 5% needs no more) */

void CloseupAdoptOpts_default(CloseupAdoptOpts *o)
{
    o->on_vox = CLOSEUP_ADOPT_ON_VOX;
    o->foot_fill_px = CLOSEUP_ADOPT_FOOT_FILL_PX;
    o->min_extent_px = CLOSEUP_ADOPT_MIN_EXTENT_PX;
    o->rms_px = CLOSEUP_ADOPT_RMS_PX;
    o->conflict_vox = CLOSEUP_ADOPT_CONFLICT_VOX;
    o->agree_frac = CLOSEUP_ADOPT_AGREE_FRAC;
    o->normal_deg = CLOSEUP_ADOPT_NORMAL_DEG;
    o->new_frac = CLOSEUP_ADOPT_NEW_FRAC;
    o->min_corr = CLOSEUP_ADOPT_MIN_CORR;
    o->metric_ratio = CLOSEUP_ADOPT_METRIC_RATIO;
    o->metric_frac = CLOSEUP_ADOPT_METRIC_FRAC;
    o->continuity_px = CLOSEUP_ADOPT_CONTINUITY_PX;
    o->continuity_slope = CLOSEUP_ADOPT_CONTINUITY_SLOPE;
    o->layer_vox = CLOSEUP_ADOPT_LAYER_VOX;
    o->layer_frac = CLOSEUP_ADOPT_LAYER_FRAC;
    o->layer_margin_px = CLOSEUP_ADOPT_LAYER_MARGIN_PX;
}

void CloseupPlacement_apply(const CloseupPlacement *p, double u, double v, double *fu, double *fv)
{
    double x = p->reflect ? -u : u;
    *fu = p->c * x - p->s * v + p->tx;
    *fv = p->s * x + p->c * v + p->ty;
}

/* Project q onto the surface near node n: normal offset *d and the foot's frame UV; the foot node may move
 * to a neighbour once.  Returns the foot node, or -1. */
static long ca_project(const CloseupGrid *g, size_t n, const float *q, double *d, double *fu, double *fv)
{
    int pass = 0;
    for (pass = 0; pass < 2; pass++) {
        size_t i = n % g->ni, j = n / g->ni;
        const float *X = &g->pos[n * 3], *tu = &g->tu[n * 3], *tv = &g->tv[n * 3], *N = &g->nrm[n * 3];
        double r[3], dd = 0.0, a = 0.0, b = 0.0, c = 0.0, e = 0.0, f = 0.0, det = 0.0, du = 0.0, dv = 0.0;
        int k = 0;
        for (k = 0; k < 3; k++) r[k] = (double)q[k] - X[k];
        for (k = 0; k < 3; k++) dd += r[k] * N[k];
        for (k = 0; k < 3; k++) r[k] -= dd * N[k];
        for (k = 0; k < 3; k++) {
            a += (double)tu[k] * tu[k]; b += (double)tu[k] * tv[k]; c += (double)tv[k] * tv[k];
            e += (double)tu[k] * r[k]; f += (double)tv[k] * r[k];
        }
        det = a * c - b * b;
        if (!(det > 1e-12)) return -1;
        du = (c * e - b * f) / det;
        dv = (a * f - b * e) / det;
        if (pass == 0 && (fabs(du) > 0.5 * g->h || fabs(dv) > 0.5 * g->h)) {
            long ni = (long)i + (long)floor(du / g->h + 0.5), nj = (long)j + (long)floor(dv / g->h + 0.5);
            if (ni >= 0 && nj >= 0 && (size_t)ni < g->ni && (size_t)nj < g->nj &&
                g->state[(size_t)nj * g->ni + (size_t)ni] != CLOSEUP_NODE_NONE) {
                n = (size_t)nj * g->ni + (size_t)ni;
                continue;
            }
        }
        *d = dd;
        *fu = g->u0 + g->h * (double)i + du;
        *fv = g->v0 + g->h * (double)j + dv;
        return (long)n;
    }
    return -1;
}

/* The grid node nearest frame (u, v), or -1 outside the grid. */
static long ca_node_at(const CloseupGrid *g, double u, double v)
{
    long i = (long)floor((u - g->u0) / g->h + 0.5), j = (long)floor((v - g->v0) / g->h + 0.5);
    if (i < 0 || j < 0 || (size_t)i >= g->ni || (size_t)j >= g->nj) return -1;
    return j * (long)g->ni + i;
}

typedef struct CaEval {
    CloseupPlacement p;
    size_t screen_hits, data_landed, data_agree, fill_landed, fill_agree;
    size_t ignored_fill_landed, ignored_fill_agree;
    double data_max_error, fill_max_error;
    double ignored_fill_max_error;
    double metric_bad_share;    /* 3-D area share of the chart's faces its own flattening distorts beyond metric_ratio */
    size_t cont_checked, cont_bad;   /* vertices on no supported node with one within continuity_px; breaking it */
    size_t layer_checked, layer_bad; /* the layer gate: sampled vertices; those stacked onto the surface */
    int ok;                 /* 1 = passes every test */
    int why;                /* 0 untouched, 1 extent, 2 rms, 3 conflict, 4 nothing new, 5 ok, 7 metric, 8 continuity,
                             * 9 layered */
} CaEval;

/* The supported node nearest each grid node within continuity_px (4-neighbour rings: Manhattan steps up to
 * continuity_px / h x sqrt 2, so the frame disk is covered); -1 where none.  Supported = DATA, or FILL within
 * foot_fill_px of the data (the band contacts and conflicts use).  malloc'd; NULL on no memory. */
static int32_t *ca_near_supported(const CloseupGrid *g, const CloseupAdoptOpts *o)
{
    const size_t nn = g->ni * g->nj;
    const unsigned maxs = (unsigned)ceil(o->continuity_px / g->h * 1.41421356237);
    int32_t *near = (int32_t *)malloc((nn ? nn : 1) * sizeof *near);
    uint16_t *steps = (uint16_t *)malloc((nn ? nn : 1) * sizeof *steps);
    size_t *queue = (size_t *)malloc((nn ? nn : 1) * sizeof *queue), head = 0, tail = 0, k = 0;
    if (near == NULL || steps == NULL || queue == NULL) { free(near); free(steps); free(queue); return NULL; }
    for (k = 0; k < nn; k++) {
        const int sup = g->state[k] == CLOSEUP_NODE_DATA ||
                        (g->state[k] == CLOSEUP_NODE_FILL && g->dist[k] <= (float)o->foot_fill_px);
        near[k] = sup ? (int32_t)k : -1;
        steps[k] = 0;
        if (sup) queue[tail++] = k;
    }
    while (head < tail) {
        const size_t a = queue[head++], i = a % g->ni, j = a / g->ni;
        int d = 0;
        if (steps[a] >= maxs) continue;
        for (d = 0; d < 4; d++) {
            const long ii = (long)i + (d == 0) - (d == 1), jj = (long)j + (d == 2) - (d == 3);
            size_t b = 0;
            if (ii < 0 || jj < 0 || (size_t)ii >= g->ni || (size_t)jj >= g->nj) continue;
            b = (size_t)jj * g->ni + (size_t)ii;
            if (near[b] >= 0) continue;
            near[b] = near[a];
            steps[b] = (uint16_t)(steps[a] + 1);
            queue[tail++] = b;
        }
    }
    free(steps); free(queue);
    return near;
}

/* A vertex (frame u, w; world q) that lands at node nd where the surface has no supported node: the nearest
 * supported node within continuity_px must lie within continuity_slope x their frame distance + conflict_vox. */
static void ca_continuity(const CloseupGrid *g, const int32_t *near_sup, long nd, double u, double w, const float *q,
                          const CloseupAdoptOpts *o, CaEval *ev)
{
    long ns = 0;
    double du = 0.0, dv = 0.0, duv = 0.0, d3 = 0.0;
    int kk = 0;
    if (near_sup == NULL || nd < 0) return;
    ns = near_sup[nd];
    if (ns < 0) return;
    du = g->u0 + g->h * (double)((size_t)ns % g->ni) - u;
    dv = g->v0 + g->h * (double)((size_t)ns / g->ni) - w;
    duv = sqrt(du * du + dv * dv);
    if (duv > o->continuity_px) return;
    for (kk = 0; kk < 3; kk++) {
        const double e = (double)q[kk] - g->pos[(size_t)ns * 3 + (size_t)kk];
        d3 += e * e;
    }
    ev->cont_checked++;
    if (sqrt(d3) > o->continuity_slope * duv + o->conflict_vox) ev->cont_bad++;
}

/* The layer gate's test of one chart vertex (world q, normal nv, placed at frame (u, w)): 1 when a surface DATA node
 * (kd over the data nodes, pmap -> grid node) lies within 3 vox of q + t n for some t = +-3, +-5, ... while t + 3 <=
 * layer_vox, is displaced from q more across than beside it (|d . n| > d / 2) and lies more than 3 x |d| +
 * layer_margin_px from (u, w) along the surface: the vertex would stack a second layer onto it. */
static int ca_layered(const CloseupGrid *g, KDTree_T kd, const int32_t *pmap, const float *q, const float *nv, double u,
                      double w, const CloseupAdoptOpts *o)
{
    double n[3], nl = 0.0, t = 0.0;
    int kk = 0;
    for (kk = 0; kk < 3; kk++) { n[kk] = (double)nv[kk]; nl += n[kk] * n[kk]; }
    nl = sqrt(nl);
    if (!(nl > 1e-9)) return 0;
    for (kk = 0; kk < 3; kk++) n[kk] /= nl;
    for (t = 3.0; t + 3.0 <= o->layer_vox + 1e-9; t += 2.0) {
        int sgn = 0;
        for (sgn = -1; sgn <= 1; sgn += 2) {
            float c[3], d2 = 0.0f;
            size_t hit = 0, nd = 0;
            double d[3], d3 = 0.0, dn = 0.0, du = 0.0, dv = 0.0;
            for (kk = 0; kk < 3; kk++) c[kk] = (float)((double)q[kk] + (double)sgn * t * n[kk]);
            hit = KDTree_nearest_bounded(kd, c, 9.0f, &d2);
            if (hit == (size_t)-1) continue;
            nd = (size_t)pmap[hit];
            for (kk = 0; kk < 3; kk++) {
                d[kk] = (double)g->pos[nd * 3 + (size_t)kk] - (double)q[kk];
                d3 += d[kk] * d[kk];
                dn += d[kk] * n[kk];
            }
            d3 = sqrt(d3);
            if (!(fabs(dn) > 0.5 * d3)) continue;
            du = g->u0 + g->h * (double)(nd % g->ni) - u;
            dv = g->v0 + g->h * (double)(nd / g->ni) - w;
            if (sqrt(du * du + dv * dv) > 3.0 * d3 + o->layer_margin_px) return 1;
        }
    }
    return 0;
}

/* The vertices of the charts accepted so far in one round (the layer gate's same-round half): a hash of 3 vox cells
 * (open addressing, chained entries) over their world points and frame UV. */
typedef struct CaRound {
    int64_t *cell;          /* cap slots: cell key, or INT64_MIN */
    int32_t *head;          /* cap slots: the first entry of that cell */
    int32_t *next;          /* entries */
    float   *pos;           /* [max*3] */
    double  *uv;            /* [max*2] */
    size_t   cap, n, max;
} CaRound;

static int64_t ca_round_key(const double p[3])
{
    int64_t c[3];
    int k = 0;
    for (k = 0; k < 3; k++) c[k] = (int64_t)floor(p[k] / 3.0) + ((int64_t)1 << 20);
    return (c[0] << 42) | (c[1] << 21) | c[2];
}

static size_t ca_round_slot(const CaRound *r, int64_t key)
{
    size_t i = (size_t)(((uint64_t)key * 0x9E3779B97F4A7C15ull) >> 20) & (r->cap - 1);
    while (r->cell[i] != INT64_MIN && r->cell[i] != key) i = (i + 1) & (r->cap - 1);
    return i;
}

static int ca_round_init(CaRound *r, size_t max)
{
    size_t k = 0;
    memset(r, 0, sizeof *r);
    r->max = max;
    r->cap = 64;
    while (r->cap < 2 * max + 2) r->cap *= 2;
    r->cell = (int64_t *)malloc(r->cap * sizeof *r->cell);
    r->head = (int32_t *)malloc(r->cap * sizeof *r->head);
    r->next = (int32_t *)malloc((max ? max : 1) * sizeof *r->next);
    r->pos = (float *)malloc((max ? max : 1) * 3 * sizeof *r->pos);
    r->uv = (double *)malloc((max ? max : 1) * 2 * sizeof *r->uv);
    if (!r->cell || !r->head || !r->next || !r->pos || !r->uv) return -1;
    for (k = 0; k < r->cap; k++) { r->cell[k] = INT64_MIN; r->head[k] = -1; }
    return 0;
}

static void ca_round_free(CaRound *r)
{
    free(r->cell); free(r->head); free(r->next); free(r->pos); free(r->uv);
    memset(r, 0, sizeof *r);
}

static void ca_round_add(CaRound *r, const float *q, double u, double w)
{
    double p[3];
    size_t s = 0;
    int64_t key = 0;
    int k = 0;
    if (r->n >= r->max) return;
    for (k = 0; k < 3; k++) { p[k] = (double)q[k]; r->pos[r->n * 3 + (size_t)k] = q[k]; }
    r->uv[r->n * 2] = u; r->uv[r->n * 2 + 1] = w;
    key = ca_round_key(p);
    s = ca_round_slot(r, key);
    if (r->cell[s] == INT64_MIN) r->cell[s] = key;
    r->next[r->n] = r->head[s];
    r->head[s] = (int32_t)r->n;
    r->n++;
}

/* ca_layered against the round's earlier charts: 1 when one of their vertices lies within 3 vox of q + t n
 * (t = +-3, +-5, ... while t + 3 <= layer_vox), more across than beside q, and farther from (u, w) in the frame than
 * 3 x its 3-D distance + layer_margin_px. */
static int ca_round_layered(const CaRound *r, const float *q, const float *nv, double u, double w,
                            const CloseupAdoptOpts *o)
{
    double n[3], nl = 0.0, t = 0.0;
    int kk = 0;
    if (r->n == 0) return 0;
    for (kk = 0; kk < 3; kk++) { n[kk] = (double)nv[kk]; nl += n[kk] * n[kk]; }
    nl = sqrt(nl);
    if (!(nl > 1e-9)) return 0;
    for (kk = 0; kk < 3; kk++) n[kk] /= nl;
    for (t = 3.0; t + 3.0 <= o->layer_vox + 1e-9; t += 2.0) {
        int sgn = 0;
        for (sgn = -1; sgn <= 1; sgn += 2) {
            double c[3];
            int dz = 0, dy = 0, dx = 0;
            for (kk = 0; kk < 3; kk++) c[kk] = (double)q[kk] + (double)sgn * t * n[kk];
            for (dz = -1; dz <= 1; dz++)
                for (dy = -1; dy <= 1; dy++)
                    for (dx = -1; dx <= 1; dx++) {
                        double cc[3];
                        size_t s = 0;
                        int32_t e = -1;
                        cc[0] = c[0] + dz * 3.0; cc[1] = c[1] + dy * 3.0; cc[2] = c[2] + dx * 3.0;
                        s = ca_round_slot(r, ca_round_key(cc));
                        for (e = r->cell[s] == INT64_MIN ? -1 : r->head[s]; e >= 0; e = r->next[e]) {
                            double d3 = 0.0, dn = 0.0, ds = 0.0, du = 0.0, dv = 0.0;
                            for (kk = 0; kk < 3; kk++) {
                                const double v = (double)r->pos[(size_t)e * 3 + (size_t)kk];
                                d3 += (v - (double)q[kk]) * (v - (double)q[kk]);
                                dn += (v - (double)q[kk]) * n[kk];
                                ds += (v - c[kk]) * (v - c[kk]);
                            }
                            if (ds > 9.0) continue;
                            d3 = sqrt(d3);
                            if (!(fabs(dn) > 0.5 * d3)) continue;
                            du = r->uv[(size_t)e * 2] - u; dv = r->uv[(size_t)e * 2 + 1] - w;
                            if (sqrt(du * du + dv * dv) > 3.0 * d3 + o->layer_margin_px) return 1;
                        }
                    }
        }
    }
    return 0;
}

/* The claimed node nearest frame (u, w) within radius_px, searched in Chebyshev rings around node nd (out to
 * sqrt 2 x the first hit's ring + 1: a nearer node can sit in a later ring); -1 when none.  *duv: its frame
 * distance. */
static long ca_claim_near(const CloseupGrid *g, const uint16_t *ccount, long nd, double u, double w, double radius_px,
                          double *duv)
{
    const long ni = (long)g->ni, nj = (long)g->nj, i0 = nd % ni, j0 = nd / ni;
    const long rmax = (long)ceil(radius_px / g->h) + 1;
    long r = 0, best = -1, hit = -1;
    double bd = 1e300;
    for (r = 0; r <= rmax && (hit < 0 || (double)r <= 1.41421356237 * (double)hit + 1.0); r++) {
        long di = 0, dj = 0;
        for (dj = -r; dj <= r; dj++)
            for (di = -r; di <= r; di++) {
                long ii = 0, jj = 0, b = 0;
                double du = 0.0, dv = 0.0, d2 = 0.0;
                if (labs(di) != r && labs(dj) != r) continue;          /* the ring only */
                ii = i0 + di; jj = j0 + dj;
                if (ii < 0 || jj < 0 || ii >= ni || jj >= nj) continue;
                b = jj * ni + ii;
                if (ccount[b] == 0) continue;
                du = g->u0 + g->h * (double)ii - u; dv = g->v0 + g->h * (double)jj - w;
                d2 = du * du + dv * dv;
                if (d2 < bd) { bd = d2; best = b; if (hit < 0) hit = r; }
            }
    }
    if (best < 0 || sqrt(bd) > radius_px) return -1;
    *duv = sqrt(bd);
    return best;
}

/* 2-D rigid fit b ~ R a + t over n pairs (weights 0/1 via keep); returns the RMS. */
static double ca_procrustes(const double *a, const double *b, const uint8_t *keep, size_t n, double *c, double *s,
                            double *tx, double *ty)
{
    double ma[2] = { 0.0, 0.0 }, mb[2] = { 0.0, 0.0 }, sxx = 0.0, sxy = 0.0, th = 0.0, e = 0.0;
    size_t k = 0, m = 0;
    for (k = 0; k < n; k++) {
        if (keep != NULL && !keep[k]) continue;
        ma[0] += a[k * 2]; ma[1] += a[k * 2 + 1]; mb[0] += b[k * 2]; mb[1] += b[k * 2 + 1];
        m++;
    }
    if (m == 0) return 1e300;
    ma[0] /= (double)m; ma[1] /= (double)m; mb[0] /= (double)m; mb[1] /= (double)m;
    for (k = 0; k < n; k++) {
        double ax = 0.0, ay = 0.0, bx = 0.0, by = 0.0;
        if (keep != NULL && !keep[k]) continue;
        ax = a[k * 2] - ma[0]; ay = a[k * 2 + 1] - ma[1]; bx = b[k * 2] - mb[0]; by = b[k * 2 + 1] - mb[1];
        sxx += ax * bx + ay * by;
        sxy += ax * by - ay * bx;
    }
    th = atan2(sxy, sxx);
    *c = cos(th); *s = sin(th);
    *tx = mb[0] - (*c * ma[0] - *s * ma[1]);
    *ty = mb[1] - (*s * ma[0] + *c * ma[1]);
    for (k = 0; k < n; k++) {
        double px = 0.0, py = 0.0;
        if (keep != NULL && !keep[k]) continue;
        px = *c * a[k * 2] - *s * a[k * 2 + 1] + *tx - b[k * 2];
        py = *s * a[k * 2] + *c * a[k * 2 + 1] + *ty - b[k * 2 + 1];
        e += px * px + py * py;
    }
    return sqrt(e / (double)m);
}

static void ca_eval_chart(const CloseupGrid *g, KDTree_T kd, const int32_t *pmap, const int32_t *near_sup,
                          KDTree_T kd_data, const int32_t *pmap_data, const CloseupPool *pool, size_t ch,
                          const CloseupAdoptOpts *o, double *A, double *B, uint8_t *keep, CaEval *ev)
{
    size_t v0 = pool->vofs[ch], v1 = pool->vofs[ch + 1], nvc = v1 - v0, stride = 1, f = 0, v = 0, n = 0, k = 0;
    double vote = 0.0, cosn = cos(o->normal_deg * 3.14159265358979 / 180.0), mean[2] = { 0.0, 0.0 };
    double sxx = 0.0, syy = 0.0, sxy = 0.0, lam = 0.0, rms = 0.0, area3 = 0.0, area_bad = 0.0;
    int reflect = 0;
    memset(ev, 0, sizeof *ev);
    ev->p.chart = ch;
    if (nvc < 3) return;
    /* the chart's own handedness: sign(det uv) (e1 x e2) . inward, summed */
    for (f = pool->fofs[ch]; f < pool->fofs[ch + 1]; f++) {
        size_t a = (size_t)pool->faces[f * 3], b = (size_t)pool->faces[f * 3 + 1], c = (size_t)pool->faces[f * 3 + 2];
        double det = (pool->uv[b * 2] - pool->uv[a * 2]) * (pool->uv[c * 2 + 1] - pool->uv[a * 2 + 1]) -
                     (pool->uv[b * 2 + 1] - pool->uv[a * 2 + 1]) * (pool->uv[c * 2] - pool->uv[a * 2]);
        double e1[3], e2[3], cr[3];
        int q = 0;
        for (q = 0; q < 3; q++) {
            e1[q] = (double)pool->verts[b * 3 + (size_t)q] - pool->verts[a * 3 + (size_t)q];
            e2[q] = (double)pool->verts[c * 3 + (size_t)q] - pool->verts[a * 3 + (size_t)q];
        }
        cr[0] = e1[1] * e2[2] - e1[2] * e2[1]; cr[1] = e1[2] * e2[0] - e1[0] * e2[2]; cr[2] = e1[0] * e2[1] - e1[1] * e2[0];
        vote += (det > 0.0 ? 1.0 : det < 0.0 ? -1.0 : 0.0) *
                (cr[0] * pool->inward[a * 3] + cr[1] * pool->inward[a * 3 + 1] + cr[2] * pool->inward[a * 3 + 2]);
        {   /* the metric gate: singular values of the chart's own UV -> 3-D map on this face (a placement is rigid,
             * so this is the metric the close-up gets); a degenerate UV face counts as distorted */
            const double a3 = 0.5 * sqrt(cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]);
            int badf = 1;
            if (fabs(det) > 1e-12) {
                const double u1 = pool->uv[b * 2] - pool->uv[a * 2], w1 = pool->uv[b * 2 + 1] - pool->uv[a * 2 + 1];
                const double u2 = pool->uv[c * 2] - pool->uv[a * 2], w2 = pool->uv[c * 2 + 1] - pool->uv[a * 2 + 1];
                /* E2 = [[u1, u2], [w1, w2]] (UV edges as columns), Inv = E2^-1, J^T J = Inv^T (E3^T E3) Inv */
                const double i00 = w2 / det, i01 = -u2 / det, i10 = -w1 / det, i11 = u1 / det;
                const double g00 = e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2];
                const double g01 = e1[0] * e2[0] + e1[1] * e2[1] + e1[2] * e2[2];
                const double g11 = e2[0] * e2[0] + e2[1] * e2[1] + e2[2] * e2[2];
                const double t00 = g00 * i00 + g01 * i10, t01 = g00 * i01 + g01 * i11;
                const double t10 = g01 * i00 + g11 * i10, t11 = g01 * i01 + g11 * i11;
                const double m00 = i00 * t00 + i10 * t10, m01 = i00 * t01 + i10 * t11, m11 = i01 * t01 + i11 * t11;
                const double tr = m00 + m11, dt = m00 * m11 - m01 * m01;
                const double disc = sqrt(fmax(0.25 * tr * tr - dt, 0.0));
                const double smax = sqrt(fmax(0.5 * tr + disc, 0.0)), smin = sqrt(fmax(0.5 * tr - disc, 0.0));
                badf = !(smax <= o->metric_ratio && smin >= 1.0 / o->metric_ratio);
            }
            area3 += a3;
            if (badf) area_bad += a3;
        }
    }
    ev->metric_bad_share = area3 > 0.0 ? area_bad / area3 : 0.0;
    reflect = vote < 0.0;
    stride = nvc > CA_MAX_SAMPLES ? (nvc + CA_MAX_SAMPLES - 1) / CA_MAX_SAMPLES : 1;
    /* contacts */
    for (v = v0; v < v1; v += stride) {
        const float *q = &pool->verts[v * 3];
        float d2 = 0.0f;
        size_t nn = (size_t)pmap[KDTree_nearest(kd, q, &d2)];     /* the KD point order -> grid node */
        double d = 0.0, fu = 0.0, fv = 0.0, dot = 0.0;
        long fn = 0;
        int kk = 0;
        if (d2 > 36.0f) continue;
        fn = ca_project(g, nn, q, &d, &fu, &fv);
        if (fn < 0 || fabs(d) > o->on_vox) continue;
        if (g->state[fn] == CLOSEUP_NODE_FILL && g->dist[fn] > (float)o->foot_fill_px) continue;
        if (g->state[fn] == CLOSEUP_NODE_NONE) continue;
        for (kk = 0; kk < 3; kk++) dot += (double)pool->nrm[v * 3 + (size_t)kk] * g->nrm[(size_t)fn * 3 + (size_t)kk];
        if (fabs(dot) < cosn) continue;
        A[n * 2] = reflect ? -pool->uv[v * 2] : pool->uv[v * 2];
        A[n * 2 + 1] = pool->uv[v * 2 + 1];
        B[n * 2] = fu; B[n * 2 + 1] = fv;
        n++;
    }
    ev->p.corr = n;
    ev->p.reflect = reflect;
    if ((long)n < o->min_corr) return;
    ev->why = 1;
    for (k = 0; k < n; k++) { mean[0] += B[k * 2]; mean[1] += B[k * 2 + 1]; }
    mean[0] /= (double)n; mean[1] /= (double)n;
    for (k = 0; k < n; k++) {
        double x = B[k * 2] - mean[0], y = B[k * 2 + 1] - mean[1];
        sxx += x * x; syy += y * y; sxy += x * y;
    }
    sxx /= (double)n; syy /= (double)n; sxy /= (double)n;
    lam = 0.5 * (sxx + syy) + sqrt(0.25 * (sxx - syy) * (sxx - syy) + sxy * sxy);
    ev->p.extent = 4.0 * sqrt(lam);          /* ~ the full length of a uniform contact */
    if (ev->p.extent < o->min_extent_px) return;
    ev->why = 2;
    rms = ca_procrustes(A, B, NULL, n, &ev->p.c, &ev->p.s, &ev->p.tx, &ev->p.ty);
    {   /* one robust refit without the contacts beyond 3 x rms (at least 1 px) */
        double lim = fmax(3.0 * rms, 1.0);
        size_t kept = 0;
        for (k = 0; k < n; k++) {
            double px = ev->p.c * A[k * 2] - ev->p.s * A[k * 2 + 1] + ev->p.tx - B[k * 2];
            double py = ev->p.s * A[k * 2] + ev->p.c * A[k * 2 + 1] + ev->p.ty - B[k * 2 + 1];
            keep[k] = (uint8_t)(sqrt(px * px + py * py) <= lim);
            kept += keep[k];
        }
        if ((long)kept >= o->min_corr) rms = ca_procrustes(A, B, keep, n, &ev->p.c, &ev->p.s, &ev->p.tx, &ev->p.ty);
    }
    ev->p.rms = rms;
    if (rms > o->rms_px) return;
    ev->why = 3;
    /* validation: every sampled vertex where it lands */
    for (v = v0; v < v1; v += stride) {
        double u = 0.0, w = 0.0, p[3], dd = 0.0;
        long nd = 0;
        int kk = 0;
        CloseupPlacement_apply(&ev->p, pool->uv[v * 2], pool->uv[v * 2 + 1], &u, &w);
        nd = ca_node_at(g, u, w);
        if (nd < 0 || g->state[nd] == CLOSEUP_NODE_NONE) {
            ev->p.fresh++;
            ca_continuity(g, near_sup, nd, u, w, &pool->verts[v * 3], o, ev);
            continue;
        }
        if (g->state[nd] == CLOSEUP_NODE_FILL) ev->p.fresh++;
        (void)CloseupGrid_taylor(g, (size_t)nd % g->ni, (size_t)nd / g->ni, u, w, p);
        for (kk = 0; kk < 3; kk++) dd += (p[kk] - pool->verts[v * 3 + (size_t)kk]) * (p[kk] - pool->verts[v * 3 + (size_t)kk]);
        /* The contact search already refuses distant extrapolation. Apply that same
         * support rule to conflicts; an inferred border must not overrule observations. */
        if (g->state[nd] == CLOSEUP_NODE_FILL && g->dist[nd] > (float)o->foot_fill_px) {
            ev->ignored_fill_landed++;
            ev->ignored_fill_agree += sqrt(dd) <= o->conflict_vox;
            ev->ignored_fill_max_error = fmax(ev->ignored_fill_max_error, sqrt(dd));
            ca_continuity(g, near_sup, nd, u, w, &pool->verts[v * 3], o, ev);   /* unsupported: continuity instead */
            continue;
        }
        ev->p.landed++;
        if (sqrt(dd) <= o->conflict_vox) ev->p.agree++;
        if (g->state[nd] == CLOSEUP_NODE_DATA) {
            ev->data_landed++;
            ev->data_agree += sqrt(dd) <= o->conflict_vox;
            ev->data_max_error = fmax(ev->data_max_error, sqrt(dd));
        } else {
            ev->fill_landed++;
            ev->fill_agree += sqrt(dd) <= o->conflict_vox;
            ev->fill_max_error = fmax(ev->fill_max_error, sqrt(dd));
        }
    }
    if (ev->p.landed == 0 || (double)ev->p.agree < o->agree_frac * (double)ev->p.landed) return;
    ev->why = 8;           /* where it adds papyrus it must continue the surface around it */
    if (ev->cont_checked > 0 && (double)ev->cont_bad > (1.0 - o->agree_frac) * (double)ev->cont_checked) return;
    ev->why = 4;
    if ((double)ev->p.fresh < o->new_frac * (double)((nvc + stride - 1) / stride)) return;
    ev->why = 7;
    if (ev->metric_bad_share > o->metric_frac) return;    /* its own flattening would smear it on the close-up */
    ev->why = 9;            /* it must not stack onto the surface: a second layer across its normal, far along it */
    if (o->layer_vox > 0.0 && kd_data != NULL) {
        const size_t ns = (nvc + stride - 1) / stride;
        const size_t lstep = stride * (ns > CA_LAYER_SAMPLES ? (ns + CA_LAYER_SAMPLES - 1) / CA_LAYER_SAMPLES : 1);
        for (v = v0; v < v1; v += lstep) {
            double u = 0.0, w = 0.0;
            CloseupPlacement_apply(&ev->p, pool->uv[v * 2], pool->uv[v * 2 + 1], &u, &w);
            ev->layer_checked++;
            ev->layer_bad += (size_t)ca_layered(g, kd_data, pmap_data, &pool->verts[v * 3], &pool->nrm[v * 3], u, w, o);
        }
        if ((double)ev->layer_bad > o->layer_frac * (double)ev->layer_checked) return;
    }
    ev->why = 5;
    ev->ok = 1;
}

/* The 3-D neighbourhood of the surface as a bitmap of CA_CELL_VOX cells: every valid node marks its cell and
 * the 26 around it (the contact band, on_vox, is well inside).  A pool chart is evaluated only when at least
 * min_corr of its vertices fall in marked cells -- the other wraps (12-14 vox away) and the rest of the box are
 * never projected. */
#define CA_CELL_VOX 2.0
typedef struct CaScreen {
    double lo[3];
    size_t n[3];
    uint8_t *bits;
} CaScreen;

static int ca_screen_build(const CloseupGrid *g, CaScreen *sc)
{
    size_t nn = g->ni * g->nj, k = 0, tot = 0;
    double hi[3] = { -1e300, -1e300, -1e300 };
    int q = 0;
    memset(sc, 0, sizeof *sc);
    sc->lo[0] = sc->lo[1] = sc->lo[2] = 1e300;
    for (k = 0; k < nn; k++) {
        if (g->state[k] == CLOSEUP_NODE_NONE) continue;
        for (q = 0; q < 3; q++) {
            double x = g->pos[k * 3 + (size_t)q];
            if (x < sc->lo[q]) sc->lo[q] = x;
            if (x > hi[q]) hi[q] = x;
        }
    }
    if (hi[0] < sc->lo[0]) return -1;
    for (q = 0; q < 3; q++) {
        sc->lo[q] -= 2.0 * CA_CELL_VOX;
        sc->n[q] = (size_t)ceil((hi[q] + 2.0 * CA_CELL_VOX - sc->lo[q]) / CA_CELL_VOX) + 1;
    }
    tot = sc->n[0] * sc->n[1] * sc->n[2];
    sc->bits = (uint8_t *)calloc(tot / 8 + 1, 1);
    if (sc->bits == NULL) return -1;
    for (k = 0; k < nn; k++) {
        long c[3], dz = 0, dy = 0, dx = 0;
        if (g->state[k] == CLOSEUP_NODE_NONE) continue;
        for (q = 0; q < 3; q++) c[q] = (long)floor((g->pos[k * 3 + (size_t)q] - sc->lo[q]) / CA_CELL_VOX);
        for (dz = -1; dz <= 1; dz++)
            for (dy = -1; dy <= 1; dy++)
                for (dx = -1; dx <= 1; dx++) {
                    long z = c[0] + dz, y = c[1] + dy, x = c[2] + dx;
                    size_t b = 0;
                    if (z < 0 || y < 0 || x < 0 || (size_t)z >= sc->n[0] || (size_t)y >= sc->n[1] || (size_t)x >= sc->n[2])
                        continue;
                    b = ((size_t)z * sc->n[1] + (size_t)y) * sc->n[2] + (size_t)x;
                    sc->bits[b >> 3] |= (uint8_t)(1u << (b & 7));
                }
    }
    return 0;
}

static int ca_screen_hit(const CaScreen *sc, const float *p)
{
    long c[3];
    size_t b = 0;
    int q = 0;
    for (q = 0; q < 3; q++) {
        c[q] = (long)floor(((double)p[q] - sc->lo[q]) / CA_CELL_VOX);
        if (c[q] < 0 || (size_t)c[q] >= sc->n[q]) return 0;
    }
    b = ((size_t)c[0] * sc->n[1] + (size_t)c[1]) * sc->n[2] + (size_t)c[2];
    return (sc->bits[b >> 3] >> (b & 7)) & 1;
}

static int ca_cmp_corr(const void *x, const void *y)
{
    const CaEval *a = (const CaEval *)x, *b = (const CaEval *)y;
    if (a->p.corr != b->p.corr) return a->p.corr > b->p.corr ? -1 : 1;
    return a->p.chart < b->p.chart ? -1 : a->p.chart > b->p.chart ? 1 : 0;
}

int CloseupAdopt_round_audited(const CloseupGrid *g, const CloseupPool *pool, const uint8_t *done,
                              const CloseupAdoptOpts *oin, CloseupPlacement **out, size_t *nout,
                              CloseupAdoptStats *st, CloseupAdoptObserver observer, void *context)
{
    CloseupAdoptOpts o;
    CloseupAdoptStats s0;
    Arena_T arena = NULL;
    KDTree_T kd = NULL, kdd = NULL;
    CaRound rset = { 0 };
    int have_rset = 0;
    size_t nn = g->ni * g->nj, k = 0, nvalid = 0, nacc = 0, maxv = 0, ndata = 0;
    float *pts = NULL, *claim = NULL, *dpts = NULL;
    int32_t *pmap = NULL, *near_sup = NULL, *dmap = NULL;
    uint16_t *ccount = NULL;
    CaEval *ev = NULL;
    CaScreen sc;
    long ch = 0;
    int rc = -1, bad = 0;
    memset(&sc, 0, sizeof sc);
    if (oin == NULL) CloseupAdoptOpts_default(&o); else o = *oin;
    if (st == NULL) st = &s0;
    memset(st, 0, sizeof *st);
    *out = NULL; *nout = 0;
    for (k = 0; k < nn; k++) nvalid += g->state[k] != CLOSEUP_NODE_NONE;
    if (nvalid == 0 || pool->nch == 0) return 0;
    pts = (float *)malloc(nvalid * 3 * sizeof *pts);
    pmap = (int32_t *)malloc(nvalid * sizeof *pmap);
    ev = (CaEval *)calloc(pool->nch, sizeof *ev);
    if (pts == NULL || pmap == NULL || ev == NULL) goto done;
    for (k = 0, nvalid = 0; k < nn; k++) {
        if (g->state[k] == CLOSEUP_NODE_NONE) continue;
        memcpy(&pts[nvalid * 3], &g->pos[k * 3], 3 * sizeof *pts);
        pmap[nvalid++] = (int32_t)k;
    }
    arena = Arena_new();
    kd = KDTree_new(arena, pts, nvalid);
    if (o.layer_vox > 0.0) {   /* the layer gate's tree: the DATA nodes only (fill is an extrapolation, not papyrus) */
        for (k = 0; k < nn; k++) ndata += g->state[k] == CLOSEUP_NODE_DATA && isfinite(g->pos[k * 3]);
        if (ndata > 0) {
            dpts = (float *)malloc(ndata * 3 * sizeof *dpts);
            dmap = (int32_t *)malloc(ndata * sizeof *dmap);
            if (dpts == NULL || dmap == NULL) goto done;
            for (k = 0, ndata = 0; k < nn; k++) {
                if (g->state[k] != CLOSEUP_NODE_DATA || !isfinite(g->pos[k * 3])) continue;
                memcpy(&dpts[ndata * 3], &g->pos[k * 3], 3 * sizeof *dpts);
                dmap[ndata++] = (int32_t)k;
            }
            kdd = KDTree_new(arena, dpts, ndata);
        }
    }
    if (ca_screen_build(g, &sc) != 0) goto done;
    near_sup = ca_near_supported(g, &o);            /* the continuity gate's reference; NULL: gate off (no memory) */
    for (k = 0; k < pool->nch; k++) {
        size_t c = pool->vofs[k + 1] - pool->vofs[k];
        if (c > maxv) maxv = c;
    }
    if (maxv > CA_MAX_SAMPLES) maxv = CA_MAX_SAMPLES + 1;
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        double *A = (double *)malloc((maxv + 1) * 2 * sizeof *A), *B = (double *)malloc((maxv + 1) * 2 * sizeof *B);
        uint8_t *keep = (uint8_t *)malloc(maxv + 1);
        if (A == NULL || B == NULL || keep == NULL) {
#ifdef _OPENMP
#pragma omp atomic
#endif
            bad++;
        } else {
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 4)
#endif
            for (ch = 0; ch < (long)pool->nch; ch++) {
                CaEval *e = &ev[ch];
                size_t v = 0, hits = 0;
                e->p.chart = (size_t)ch;
                if (done != NULL && done[ch]) continue;
                for (v = pool->vofs[ch]; v < pool->vofs[ch + 1] && (long)hits < o.min_corr; v++)
                    hits += (size_t)ca_screen_hit(&sc, &pool->verts[v * 3]);
                e->screen_hits = hits;
                if ((long)hits < o.min_corr) continue;          /* never near the sheet */
                ca_eval_chart(g, kd, pmap, near_sup, kdd, dmap, pool, (size_t)ch, &o, A, B, keep, e);
                e->screen_hits = hits;   /* ca_eval_chart initializes its result */
            }
        }
        free(A); free(B); free(keep);
    }
    if (bad) goto done;
    for (ch = 0; ch < (long)pool->nch; ch++) {
        CaEval *e = &ev[ch];
        if (done != NULL && done[ch]) continue;
        st->candidates++;
        if ((long)e->p.corr >= o.min_corr) st->touched++;
        if (e->why == 1) st->refused_extent++;
        else if (e->why == 2) st->refused_rms++;
        else if (e->why == 3) st->refused_conflict++;
        else if (e->why == 4) st->refused_nothing_new++;
        else if (e->why == 7) st->refused_metric++;
        else if (e->why == 8) st->refused_continuity++;
        else if (e->why == 9) st->refused_layer++;
    }
    /* claims: the accepted charts in order of support; a claim that disagrees with an earlier one waits */
    qsort(ev, pool->nch, sizeof *ev, ca_cmp_corr);
    claim = (float *)calloc(nn * 3, sizeof *claim);
    ccount = (uint16_t *)calloc(nn, sizeof *ccount);
    *out = (CloseupPlacement *)malloc((pool->nch ? pool->nch : 1) * sizeof **out);
    if (claim == NULL || ccount == NULL || *out == NULL) goto done;
    if (o.layer_vox > 0.0) {   /* the same-round half of the layer gate: the accepted charts' sampled vertices */
        size_t nok = 0;
        for (k = 0; k < pool->nch; k++) nok += (size_t)ev[k].ok;
        have_rset = 1;                                   /* set first: a failed init frees what it got */
        if (ca_round_init(&rset, nok * (CA_LAYER_SAMPLES + 1) + 1) != 0) goto done;
    }
    for (k = 0; k < pool->nch; k++) {
        CaEval *e = &ev[k];
        size_t v0 = 0, v1 = 0, v = 0, stride = 1, clash = 0, fresh_in = 0;
        if (!e->ok) continue;
        v0 = pool->vofs[e->p.chart]; v1 = pool->vofs[e->p.chart + 1];
        stride = (v1 - v0) > CA_MAX_SAMPLES ? (v1 - v0 + CA_MAX_SAMPLES - 1) / CA_MAX_SAMPLES : 1;
        for (v = v0; v < v1; v += stride) {
            double u = 0.0, w = 0.0, dd = 0.0;
            long nd = 0;
            int kk = 0;
            CloseupPlacement_apply(&e->p, pool->uv[v * 2], pool->uv[v * 2 + 1], &u, &w);
            nd = ca_node_at(g, u, w);
            if (nd < 0 || g->state[nd] == CLOSEUP_NODE_DATA) continue;
            fresh_in++;
            if (ccount[nd] == 0) {
                /* no claim here: the nearest claim of an earlier chart of this round must continue it (two charts
                 * laid side by side in the frame but far apart in 3-D are a smear waiting to happen) */
                double duv = 0.0, d3 = 0.0;
                long nb = ca_claim_near(g, ccount, nd, u, w, o.continuity_px, &duv);
                if (nb >= 0) {
                    for (kk = 0; kk < 3; kk++) {
                        const double m = (double)claim[(size_t)nb * 3 + (size_t)kk] / ccount[nb];
                        d3 += (m - pool->verts[v * 3 + (size_t)kk]) * (m - pool->verts[v * 3 + (size_t)kk]);
                    }
                    if (sqrt(d3) > o.continuity_slope * duv + o.conflict_vox) clash++;
                }
                continue;
            }
            for (kk = 0; kk < 3; kk++) {
                double m = (double)claim[(size_t)nd * 3 + (size_t)kk] / ccount[nd];
                dd += (m - pool->verts[v * 3 + (size_t)kk]) * (m - pool->verts[v * 3 + (size_t)kk]);
            }
            if (sqrt(dd) > o.conflict_vox) clash++;
        }
        if (fresh_in > 0 && (double)clash > 0.1 * (double)fresh_in) { st->deferred++; e->why = 6; continue; }
        if (have_rset && rset.n > 0) {   /* stacked onto a chart accepted earlier in this round: it waits for the next */
            const size_t ns = (v1 - v0 + stride - 1) / stride;
            const size_t lstep = stride * (ns > CA_LAYER_SAMPLES ? (ns + CA_LAYER_SAMPLES - 1) / CA_LAYER_SAMPLES : 1);
            size_t checked = 0, stacked = 0;
            for (v = v0; v < v1; v += lstep) {
                double u = 0.0, w = 0.0;
                CloseupPlacement_apply(&e->p, pool->uv[v * 2], pool->uv[v * 2 + 1], &u, &w);
                checked++;
                stacked += (size_t)ca_round_layered(&rset, &pool->verts[v * 3], &pool->nrm[v * 3], u, w, &o);
            }
            if ((double)stacked > o.layer_frac * (double)checked) {
                st->deferred++; st->deferred_layer++; e->why = 6;
                continue;
            }
        }
        for (v = v0; v < v1; v += stride) {
            double u = 0.0, w = 0.0;
            long nd = 0;
            int kk = 0;
            CloseupPlacement_apply(&e->p, pool->uv[v * 2], pool->uv[v * 2 + 1], &u, &w);
            nd = ca_node_at(g, u, w);
            if (nd < 0 || g->state[nd] == CLOSEUP_NODE_DATA || ccount[nd] == UINT16_MAX) continue;
            for (kk = 0; kk < 3; kk++) claim[(size_t)nd * 3 + (size_t)kk] += pool->verts[v * 3 + (size_t)kk];
            ccount[nd]++;
        }
        if (have_rset) {
            const size_t ns = (v1 - v0 + stride - 1) / stride;
            const size_t lstep = stride * (ns > CA_LAYER_SAMPLES ? (ns + CA_LAYER_SAMPLES - 1) / CA_LAYER_SAMPLES : 1);
            for (v = v0; v < v1; v += lstep) {
                double u = 0.0, w = 0.0;
                CloseupPlacement_apply(&e->p, pool->uv[v * 2], pool->uv[v * 2 + 1], &u, &w);
                ca_round_add(&rset, &pool->verts[v * 3], u, w);
            }
        }
        (*out)[nacc++] = e->p;
    }
    st->placed = nacc;
    *nout = nacc;
    if (observer != NULL) {
        for (k = 0; k < pool->nch; k++) {
            const CaEval *e = &ev[k];
            CloseupAdoptAudit record;
            if (done != NULL && done[e->p.chart]) continue;
            memset(&record, 0, sizeof record);
            record.p = e->p; record.reason = e->why; record.screen_hits = e->screen_hits;
            record.data_landed = e->data_landed; record.data_agree = e->data_agree;
            record.fill_landed = e->fill_landed; record.fill_agree = e->fill_agree;
            record.data_max_error = e->data_max_error; record.fill_max_error = e->fill_max_error;
            record.ignored_fill_landed = e->ignored_fill_landed;
            record.ignored_fill_agree = e->ignored_fill_agree;
            record.ignored_fill_max_error = e->ignored_fill_max_error;
            record.metric_bad_share = e->metric_bad_share;
            record.cont_checked = e->cont_checked; record.cont_bad = e->cont_bad;
            record.layer_checked = e->layer_checked; record.layer_bad = e->layer_bad;
            if (observer(&record, context) != 0) goto done;
        }
    }
    rc = 0;
done:
    if (rc != 0) { free(*out); *out = NULL; *nout = 0; }
    free(pts); free(pmap); free(near_sup); free(ev); free(claim); free(ccount); free(sc.bits); free(dpts); free(dmap);
    if (have_rset) ca_round_free(&rset);
    Arena_dispose(&arena);
    return rc;
}

int CloseupAdopt_round(const CloseupGrid *g, const CloseupPool *pool, const uint8_t *done,
                       const CloseupAdoptOpts *o, CloseupPlacement **out, size_t *nout, CloseupAdoptStats *st)
{
    return CloseupAdopt_round_audited(g, pool, done, o, out, nout, st, NULL, NULL);
}

/* ------------------------------------------------------------------ selftest */

/* A cylinder sheet (radius 150) observed by one chart on u in [0, 200); the pool holds (a) its continuation
 * u in [201, 330) in a rotated, shifted, mirror-imaged chart frame, (b) the next wrap 13 vox inward over the
 * same u range as the observed chart, touching nothing, (c) a patch of the next wrap placed so it touches
 * the sheet at one corner.  (a) must be adopted with its placement within 0.5 px; (b) and (c) refused. */
static void ca_cyl(double R, double u, double v, float p[3])
{
    p[0] = (float)v; p[1] = (float)(R * sin(u / R)); p[2] = (float)(R * cos(u / R));
}

typedef struct CaTestAudit {
    CloseupAdoptAudit records[3];
    size_t seen;
} CaTestAudit;

static int ca_test_observer(const CloseupAdoptAudit *record, void *context)
{
    CaTestAudit *test = (CaTestAudit *)context;
    if (record->p.chart >= 3) return -1;
    test->records[record->p.chart] = *record;
    test->seen++;
    return 0;
}

/* ca_layered on a hand-made grid: DATA nodes (i < 60, j < 20) on the plane x = 50 at (z, y) = (j, i), h = 1.  A chart
 * vertex 7 vox across that plane (x 57, normal +x) stacks onto it when it is placed far along the surface (frame u 150),
 * not when placed beside the data it faces (u 36: a tight bend), not when it lies beside the data's edge, not 14 vox
 * across, not with the reach 6. */
static int ca_layer_test(void)
{
    CloseupGrid g;
    CloseupAdoptOpts o;
    Arena_T arena = NULL;
    KDTree_T kd = NULL;
    float *pts = NULL;
    int32_t *map = NULL;
    size_t i = 0, j = 0, n = 0;
    int fails = 0, r[5];
    const float nx[3] = { 0.0f, 0.0f, 1.0f }, across[3] = { 10.0f, 30.0f, 57.0f }, beside[3] = { 10.0f, 61.0f, 51.0f },
                far3[3] = { 10.0f, 30.0f, 64.0f };
    memset(&g, 0, sizeof g);
    if (CloseupGrid_alloc(&g, 0.0, 199.0, 0.0, 19.0, 1.0) != 0) return 1;
    pts = (float *)malloc(60 * g.nj * 3 * sizeof *pts);
    map = (int32_t *)malloc(60 * g.nj * sizeof *map);
    if (pts == NULL || map == NULL) { free(pts); free(map); CloseupGrid_free(&g); return 1; }
    for (j = 0; j < g.nj; j++)
        for (i = 0; i < 60 && i < g.ni; i++) {
            const size_t k = j * g.ni + i;
            g.pos[k * 3] = (float)j; g.pos[k * 3 + 1] = (float)i; g.pos[k * 3 + 2] = 50.0f;
            g.state[k] = CLOSEUP_NODE_DATA;
            memcpy(&pts[n * 3], &g.pos[k * 3], 3 * sizeof *pts);
            map[n++] = (int32_t)k;
        }
    arena = Arena_new();
    kd = KDTree_new(arena, pts, n);
    CloseupAdoptOpts_default(&o);
    r[0] = ca_layered(&g, kd, map, across, nx, 150.0, 10.0, &o);   /* far along the surface: a second layer */
    r[1] = ca_layered(&g, kd, map, across, nx, 36.0, 10.0, &o);    /* near along it: a tight bend */
    r[2] = ca_layered(&g, kd, map, beside, nx, 150.0, 10.0, &o);   /* beside the data's edge, not across */
    r[3] = ca_layered(&g, kd, map, far3, nx, 150.0, 10.0, &o);     /* 14 vox across: beyond reach */
    o.layer_vox = 6.0;
    r[4] = ca_layered(&g, kd, map, across, nx, 150.0, 10.0, &o);   /* reach 6 does not see 7 vox */
    fails += r[0] != 1 || r[1] != 0 || r[2] != 0 || r[3] != 0 || r[4] != 0;
    fprintf(stderr, "  closeup adopt layer gate: %d %d %d %d %d (want 1 0 0 0 0) %s\n", r[0], r[1], r[2], r[3], r[4],
            fails ? "FAIL" : "PASS");
    {   /* the same tests against a round's earlier chart (the plane's points as its accepted vertices) */
        CaRound rs;
        int q[4], f2 = 0;
        o.layer_vox = CLOSEUP_ADOPT_LAYER_VOX;
        if (ca_round_init(&rs, n + 1) != 0) { ca_round_free(&rs); f2 = 1; }
        else {
            size_t e = 0;
            for (e = 0; e < n; e++) {
                const size_t kk = (size_t)map[e];
                ca_round_add(&rs, &g.pos[kk * 3], (double)(kk % g.ni), (double)(kk / g.ni));
            }
            q[0] = ca_round_layered(&rs, across, nx, 150.0, 10.0, &o);
            q[1] = ca_round_layered(&rs, across, nx, 36.0, 10.0, &o);
            q[2] = ca_round_layered(&rs, beside, nx, 150.0, 10.0, &o);
            q[3] = ca_round_layered(&rs, far3, nx, 150.0, 10.0, &o);
            f2 = q[0] != 1 || q[1] != 0 || q[2] != 0 || q[3] != 0;
            fprintf(stderr, "  closeup adopt layer gate, same round: %d %d %d %d (want 1 0 0 0) %s\n", q[0], q[1], q[2], q[3],
                    f2 ? "FAIL" : "PASS");
            ca_round_free(&rs);
        }
        fails += f2;
    }
    Arena_dispose(&arena);
    free(pts); free(map);
    CloseupGrid_free(&g);
    return fails;
}

int CloseupAdopt_selftest(void)
{
    const double R = 150.0;
    CloseupObs obs;
    CloseupGrid g;
    CloseupPool pool;
    CloseupPlacement *pl = NULL;
    CloseupAdoptStats st;
    float *pv = NULL, *pn = NULL, *pin = NULL, *inw = NULL;
    double *puv = NULL;
    int32_t *pf = NULL;
    size_t vofs[4], fofs[4], nv = 0, nf = 0, npl = 0, k = 0;
    int x = 0, y = 0, fails = 0, ch = 0;
    double worst = 0.0;
    memset(&obs, 0, sizeof obs);
    memset(&pool, 0, sizeof pool);
    obs.cap = 200 * 120;
    obs.uv = (float *)malloc(obs.cap * 2 * sizeof(float)); obs.p = (float *)malloc(obs.cap * 3 * sizeof(float));
    obs.w = (float *)malloc(obs.cap * sizeof(float)); obs.chart = (int32_t *)malloc(obs.cap * sizeof(int32_t));
    if (obs.uv == NULL || obs.p == NULL || obs.w == NULL || obs.chart == NULL) { CloseupObs_free(&obs); return 1; }
    for (y = 0; y < 120; y++)
        for (x = 0; x < 200; x++) {
            ca_cyl(R, x + 0.5, y + 0.5, &obs.p[obs.n * 3]);
            obs.uv[obs.n * 2] = (float)(x + 0.5); obs.uv[obs.n * 2 + 1] = (float)(y + 0.5);
            obs.w[obs.n] = 1.0f; obs.chart[obs.n] = 0; obs.n++;
        }
    if (CloseupGrid_alloc(&g, -20.0, 420.0, -20.0, 140.0, 2.0) != 0 || CloseupSurface_fit(&g, &obs, NULL, NULL) != 0) {
        CloseupObs_free(&obs); CloseupGrid_free(&g); return 1;
    }
    inw = (float *)calloc(g.ni * g.nj * 3, sizeof *inw);
    if (inw == NULL) goto out;
    for (k = 0; k < g.ni * g.nj; k++) {
        double l = 0.0;
        if (g.state[k] == CLOSEUP_NODE_NONE) continue;
        l = sqrt((double)g.pos[k * 3 + 1] * g.pos[k * 3 + 1] + (double)g.pos[k * 3 + 2] * g.pos[k * 3 + 2]);
        inw[k * 3 + 1] = (float)(-g.pos[k * 3 + 1] / l); inw[k * 3 + 2] = (float)(-g.pos[k * 3 + 2] / l);
    }
    if (CloseupSurface_orient(&g, inw, NULL) != 1) { fails++; goto out; }
    /* the pool: three grid charts at 2 px spacing */
    pv = (float *)malloc(3 * 80 * 80 * 3 * sizeof *pv);
    pn = (float *)malloc(3 * 80 * 80 * 3 * sizeof *pn);
    pin = (float *)malloc(3 * 80 * 80 * 3 * sizeof *pin);
    puv = (double *)malloc(3 * 80 * 80 * 2 * sizeof *puv);
    pf = (int32_t *)malloc(3 * 80 * 80 * 6 * sizeof *pf);
    if (pv == NULL || pn == NULL || pin == NULL || puv == NULL || pf == NULL) goto out;
    for (ch = 0; ch < 3; ch++) {
        /* chart frames: (a) u in [201, 331), v in [0, 120); (b) inward wrap u [20, 180); (c) inward wrap near u 200 */
        const double ua = ch == 0 ? 201.0 : ch == 1 ? 20.0 : 150.0, ub = ch == 0 ? 331.0 : ch == 1 ? 180.0 : 230.0;
        /* (c) leaves the sheet ON the data at u 150 and reaches the inward wrap at u 230: its contacts are real
         * (u 150..162) but most of it lands on the surface 3-13 vox off -- the conflict gate must refuse it */
        const double Rw = ch == 0 ? R : R - 13.0;
        const int nx = 66, ny = 61;
        const double th = ch == 0 ? 0.7 : 0.2, cth = cos(th), sth = sin(th);
        size_t base = nv;
        vofs[ch] = nv; fofs[ch] = nf;
        for (y = 0; y < ny; y++)
            for (x = 0; x < nx; x++) {
                double u = ua + (ub - ua) * x / (nx - 1), v = 120.0 * y / (ny - 1), lu = 0.0, lv = 0.0, l = 0.0;
                double cu = ch == 2 ? u : u * (Rw / R);     /* arc length on the wrap itself */
                float p[3];
                if (ch == 2) {                        /* the touching patch bends from the sheet to the inward wrap */
                    double t = (u - ua) / (ub - ua);
                    double rr = R - 13.0 * t;
                    p[0] = (float)v; p[1] = (float)(rr * sin(u / R)); p[2] = (float)(rr * cos(u / R));
                } else {
                    ca_cyl(Rw, cu, v, p);
                }
                memcpy(&pv[nv * 3], p, sizeof p);
                l = sqrt((double)p[1] * p[1] + (double)p[2] * p[2]);
                pn[nv * 3] = 0.0f; pn[nv * 3 + 1] = (float)(p[1] / l); pn[nv * 3 + 2] = (float)(p[2] / l);
                pin[nv * 3] = 0.0f; pin[nv * 3 + 1] = (float)(-p[1] / l); pin[nv * 3 + 2] = (float)(-p[2] / l);
                /* the chart's own flattening: arc length, rotated by th, shifted, and mirror-imaged for (a) */
                lu = (ch == 0 ? -1.0 : 1.0) * (cth * cu - sth * v) + 1000.0 * ch;
                lv = sth * cu + cth * v - 500.0;
                puv[nv * 2] = lu; puv[nv * 2 + 1] = lv;
                nv++;
            }
        for (y = 0; y + 1 < ny; y++)
            for (x = 0; x + 1 < nx; x++) {
                int32_t a = (int32_t)(base + (size_t)(y * nx + x)), b = a + 1, c = a + nx, d = c + 1;
                pf[nf * 3] = a; pf[nf * 3 + 1] = b; pf[nf * 3 + 2] = d; nf++;
                pf[nf * 3] = a; pf[nf * 3 + 1] = d; pf[nf * 3 + 2] = c; nf++;
            }
    }
    vofs[3] = nv; fofs[3] = nf;
    pool.nv = nv; pool.nf = nf; pool.nch = 3; pool.verts = pv; pool.nrm = pn; pool.inward = pin; pool.uv = puv;
    pool.faces = pf; pool.vofs = vofs; pool.fofs = fofs;
    if (CloseupAdopt_round(&g, &pool, NULL, NULL, &pl, &npl, &st) != 0) { fails++; goto out; }
    for (k = 0; k < npl; k++) {
        size_t v = 0;
        if (pl[k].chart != 0) { fails++; continue; }
        for (v = vofs[0]; v < vofs[1]; v++) {
            double fu = 0.0, fv = 0.0, cu = 0.0, du = 0.0, dv = 0.0;
            CloseupPlacement_apply(&pl[k], puv[v * 2], puv[v * 2 + 1], &fu, &fv);
            /* truth: frame u = arc length, v = z */
            cu = atan2((double)pv[v * 3 + 1], (double)pv[v * 3 + 2]) * R;
            du = fu - cu; dv = fv - pv[v * 3];
            worst = fmax(worst, sqrt(du * du + dv * dv));
        }
    }
    fprintf(stderr, "  closeup adopt: %zu candidates, %zu touched, %zu placed (refused: extent %zu, rms %zu, conflict %zu, "
            "nothing new %zu; deferred %zu); continuation placed within %.3f px\n", st.candidates, st.touched, st.placed,
            st.refused_extent, st.refused_rms, st.refused_conflict, st.refused_nothing_new, st.deferred, worst);
    if (npl != 1 || worst > 0.5 || st.refused_conflict < 1) fails++;
    {
        CloseupPlacement *audited = NULL;
        CloseupAdoptStats ast;
        CaTestAudit audit;
        size_t count = 0;
        memset(&audit, 0, sizeof audit);
        if (CloseupAdopt_round_audited(&g, &pool, NULL, NULL, &audited, &count, &ast,
                                       ca_test_observer, &audit) != 0 || count != npl || audit.seen != 3) {
            fails++;
        } else {
            for (k = 0; k < npl; k++) if (memcmp(&pl[k], &audited[k], sizeof *pl) != 0) fails++;
            if (memcmp(&st, &ast, sizeof st) != 0 || audit.records[0].reason != 5 || audit.records[2].reason != 3)
                fails++;
            for (k = 0; k < 3; k++) {
                const CloseupAdoptAudit *r = &audit.records[k];
                if (r->p.landed != r->data_landed + r->fill_landed || r->p.agree != r->data_agree + r->fill_agree)
                    fails++;
            }
            if (audit.records[2].data_agree >= audit.records[2].data_landed) fails++;
        }
        free(audited);
        fprintf(stderr, "  closeup adopt audit: %zu records; unchanged placements and data/fill counts %s\n",
                audit.seen, fails ? "FAIL" : "PASS");
    }
    {   /* Corrupt only extrapolation outside the existing trusted contact band.
         * The true continuation must still join; the touching next wrap must not. */
        CloseupPlacement *checked = NULL;
        CloseupAdoptStats ast;
        CaTestAudit audit;
        CloseupAdoptOpts opts;
        size_t count = 0, changed = 0;
        CloseupAdoptOpts_default(&opts);
        memset(&audit, 0, sizeof audit);
        for (k = 0; k < g.ni * g.nj; k++) {
            int q = 0;
            if (g.state[k] != CLOSEUP_NODE_FILL || g.dist[k] <= (float)opts.foot_fill_px) continue;
            for (q = 0; q < 3; q++) g.pos[k * 3 + (size_t)q] += 20.0f * g.nrm[k * 3 + (size_t)q];
            changed++;
        }
        if (CloseupAdopt_round_audited(&g, &pool, NULL, NULL, &checked, &count, &ast,
                                       ca_test_observer, &audit) != 0 || count != 1 || changed == 0) {
            fails++;
        } else if (checked[0].chart != 0 || audit.records[2].reason != 3 ||
                   audit.records[0].ignored_fill_landed == 0 || audit.records[0].ignored_fill_max_error < 10.0) {
            fails++;
        } else {
            const CloseupAdoptAudit *r = &audit.records[0];
            /* The previous all-fill conflict gate must fail this same placement. */
            if ((double)(r->p.agree + r->ignored_fill_agree) >=
                opts.agree_frac * (double)(r->p.landed + r->ignored_fill_landed)) fails++;
        }
        free(checked);
        fprintf(stderr, "  closeup adopt distant fill: %zu poisoned nodes, %zu accepted; true sheet only %s\n",
                changed, count, fails ? "FAIL" : "PASS");
    }
    {   /* The metric gate (2026-10-07): the continuation (a) alone -- isometric, adopted, bad share ~0 -- then the
         * same chart with the far 15 of its 65 columns squeezed 8x in its OWN flattening.  Its contacts, fit,
         * conflicts and fresh share stay those of (a) (the squeezed part lands on unsupported fill), so only the
         * metric gate can refuse it. */
        CloseupPool one = pool;
        const double ua = 201.0, ub = 331.0, cth = cos(0.7), sth = sin(0.7);
        const double c50 = ua + (ub - ua) * 50.0 / 65.0;
        double *suv = (double *)malloc(vofs[1] * 2 * sizeof *suv);
        size_t vo[2], fo[2];
        int pass = 0;
        if (suv == NULL) { fails++; goto out; }
        vo[0] = 0; vo[1] = vofs[1]; fo[0] = 0; fo[1] = fofs[1];
        one.nch = 1; one.nv = vofs[1]; one.nf = fofs[1]; one.vofs = vo; one.fofs = fo;
        for (pass = 0; pass < 2; pass++) {
            CloseupPlacement *m = NULL;
            CloseupAdoptStats mst;
            CaTestAudit audit;
            size_t v = 0, count = 0;
            for (v = 0; v < vofs[1]; v++) {      /* chart (a) as generated above: x = v % 66, y = v / 66 */
                const double u = ua + (ub - ua) * (double)(v % 66) / 65.0, vv = 120.0 * (double)(v / 66) / 60.0;
                double cu = u;
                if (pass == 1 && cu > c50) cu = c50 + (cu - c50) / 8.0;
                suv[v * 2] = -(cth * cu - sth * vv);
                suv[v * 2 + 1] = sth * cu + cth * vv - 500.0;
            }
            one.uv = suv;
            memset(&audit, 0, sizeof audit);
            if (CloseupAdopt_round_audited(&g, &one, NULL, NULL, &m, &count, &mst, ca_test_observer, &audit) != 0 ||
                audit.seen != 1) {
                fails++;
            } else if (pass == 0 && (count != 1 || audit.records[0].reason != 5 || mst.refused_metric != 0 ||
                                     audit.records[0].metric_bad_share > 1e-3)) {
                fails++;
            } else if (pass == 1 && (count != 0 || audit.records[0].reason != 7 || mst.refused_metric != 1 ||
                                     audit.records[0].metric_bad_share < 0.05)) {
                fails++;
            }
            fprintf(stderr, "  closeup adopt metric gate, %s chart: %zu placed, reason %d, bad share %.3f %s\n",
                    pass ? "squeezed" : "isometric", count, audit.records[0].reason, audit.records[0].metric_bad_share,
                    fails ? "FAIL" : "PASS");
            free(m);
        }
        free(suv);
    }
    {   /* The continuity gate (2026-10-07): a second observed patch B (u in [345, 390)) beyond the end of the
         * continuation (a) (which ends at u 331: inside the continuity reach of B's border, outside B's supported
         * fill).  (i) B continues the sheet -- (a) is adopted; (ii) B lies 60 vox along the axis from where the frame
         * puts it (another sheet laid out next to this one) -- (a)'s end would be bridged to it, so (a) is refused
         * for continuity (reason 8), its contacts, fit, conflicts, fresh share and metric being those of (i). */
        int pass = 0;
        for (pass = 0; pass < 2; pass++) {
            CloseupObs ob;
            CloseupGrid gg;
            CloseupPool one = pool;
            CloseupPlacement *m = NULL;
            CloseupAdoptStats mst;
            CaTestAudit audit;
            float *iw = NULL;
            size_t vo[2], fo[2], count = 0, q = 0;
            int xx = 0, yy = 0, ok = 1;
            memset(&ob, 0, sizeof ob);
            memset(&gg, 0, sizeof gg);
            memset(&audit, 0, sizeof audit);
            ob.cap = 245 * 120;
            ob.uv = (float *)malloc(ob.cap * 2 * sizeof(float)); ob.p = (float *)malloc(ob.cap * 3 * sizeof(float));
            ob.w = (float *)malloc(ob.cap * sizeof(float)); ob.chart = (int32_t *)malloc(ob.cap * sizeof(int32_t));
            if (ob.uv == NULL || ob.p == NULL || ob.w == NULL || ob.chart == NULL) { CloseupObs_free(&ob); fails++; break; }
            for (yy = 0; yy < 120; yy++)
                for (xx = 0; xx < 245; xx++) {         /* 200 columns of the sheet, then B's 45 */
                    const double uu = xx < 200 ? xx + 0.5 : 345.0 + (double)(xx - 200) + 0.5;
                    const double zoff = (pass == 1 && xx >= 200) ? 60.0 : 0.0;
                    ca_cyl(R, uu, yy + 0.5 + zoff, &ob.p[ob.n * 3]);
                    ob.uv[ob.n * 2] = (float)uu; ob.uv[ob.n * 2 + 1] = (float)(yy + 0.5);
                    ob.w[ob.n] = 1.0f; ob.chart[ob.n] = xx < 200 ? 0 : 1; ob.n++;
                }
            if (CloseupGrid_alloc(&gg, -20.0, 420.0, -20.0, 140.0, 2.0) != 0 || CloseupSurface_fit(&gg, &ob, NULL, NULL) != 0 ||
                (iw = (float *)calloc(gg.ni * gg.nj * 3, sizeof *iw)) == NULL) {
                CloseupObs_free(&ob); CloseupGrid_free(&gg); fails++; break;
            }
            for (q = 0; q < gg.ni * gg.nj; q++) {
                double l = 0.0;
                if (gg.state[q] == CLOSEUP_NODE_NONE) continue;
                l = sqrt((double)gg.pos[q * 3 + 1] * gg.pos[q * 3 + 1] + (double)gg.pos[q * 3 + 2] * gg.pos[q * 3 + 2]);
                iw[q * 3 + 1] = (float)(-gg.pos[q * 3 + 1] / l); iw[q * 3 + 2] = (float)(-gg.pos[q * 3 + 2] / l);
            }
            (void)CloseupSurface_orient(&gg, iw, NULL);
            vo[0] = 0; vo[1] = vofs[1]; fo[0] = 0; fo[1] = fofs[1];
            one.nch = 1; one.nv = vofs[1]; one.nf = fofs[1]; one.vofs = vo; one.fofs = fo; one.uv = puv;
            if (CloseupAdopt_round_audited(&gg, &one, NULL, NULL, &m, &count, &mst, ca_test_observer, &audit) != 0 ||
                audit.seen != 1)
                ok = 0;
            else if (pass == 0 && (count != 1 || audit.records[0].reason != 5 || audit.records[0].cont_bad != 0))
                ok = 0;
            else if (pass == 1 && (count != 0 || audit.records[0].reason != 8 || mst.refused_continuity != 1))
                ok = 0;
            if (audit.records[0].cont_checked == 0) ok = 0;    /* the gate must have looked at (a)'s ends */
            if (!ok) fails++;
            fprintf(stderr, "  closeup adopt continuity, B %s: %zu placed, reason %d, %zu of %zu checked vertices break "
                    "it %s\n", pass ? "60 vox off" : "continuing the sheet", count, audit.records[0].reason,
                    audit.records[0].cont_bad, audit.records[0].cont_checked, ok ? "PASS" : "FAIL");
            free(m); free(iw);
            CloseupObs_free(&ob);
            CloseupGrid_free(&gg);
        }
    }
    {   /* The same-round claims (2026-10-07): P continues the sheet past u 201 over rows v [0, 56]; Q lies on the
         * sheet over u [150, 208) and bends inward past it (radius R - (u - 208) / 2) over rows v [64, 120].  Each
         * passes every gate alone (Q's bend stays within the continuity slope of the sheet's supported edge); in
         * one round the second to claim lands beside the first's claims, 8-24 px away in the frame but up to 60 vox
         * in 3-D, so it waits (reason 6).  Control: continuity_px 0 (the gate off) places both side by side. */
        CloseupObs ob;
        CloseupGrid gg;
        float *sv = NULL, *sn = NULL, *si = NULL, *iw = NULL;
        double *su = NULL;
        int32_t *sf = NULL;
        size_t svo[3], sfo[3], snv = 0, snf = 0, q = 0;
        int pass = 0, c2 = 0, xx = 0, yy = 0;
        memset(&ob, 0, sizeof ob);
        memset(&gg, 0, sizeof gg);
        ob.cap = 200 * 120;
        ob.uv = (float *)malloc(ob.cap * 2 * sizeof(float)); ob.p = (float *)malloc(ob.cap * 3 * sizeof(float));
        ob.w = (float *)malloc(ob.cap * sizeof(float)); ob.chart = (int32_t *)malloc(ob.cap * sizeof(int32_t));
        sv = (float *)malloc(2 * 91 * 29 * 3 * sizeof *sv); sn = (float *)malloc(2 * 91 * 29 * 3 * sizeof *sn);
        si = (float *)malloc(2 * 91 * 29 * 3 * sizeof *si); su = (double *)malloc(2 * 91 * 29 * 2 * sizeof *su);
        sf = (int32_t *)malloc(2 * 90 * 28 * 6 * sizeof *sf);
        if (ob.uv == NULL || ob.p == NULL || ob.w == NULL || ob.chart == NULL || sv == NULL || sn == NULL ||
            si == NULL || su == NULL || sf == NULL) {
            CloseupObs_free(&ob); free(sv); free(sn); free(si); free(su); free(sf); fails++; goto out;
        }
        for (yy = 0; yy < 120; yy++)
            for (xx = 0; xx < 200; xx++) {
                ca_cyl(R, xx + 0.5, yy + 0.5, &ob.p[ob.n * 3]);
                ob.uv[ob.n * 2] = (float)(xx + 0.5); ob.uv[ob.n * 2 + 1] = (float)(yy + 0.5);
                ob.w[ob.n] = 1.0f; ob.chart[ob.n] = 0; ob.n++;
            }
        if (CloseupGrid_alloc(&gg, -20.0, 420.0, -20.0, 140.0, 2.0) != 0 || CloseupSurface_fit(&gg, &ob, NULL, NULL) != 0 ||
            (iw = (float *)calloc(gg.ni * gg.nj * 3, sizeof *iw)) == NULL) {
            CloseupObs_free(&ob); CloseupGrid_free(&gg); free(sv); free(sn); free(si); free(su); free(sf);
            fails++; goto out;
        }
        for (q = 0; q < gg.ni * gg.nj; q++) {
            double l = 0.0;
            if (gg.state[q] == CLOSEUP_NODE_NONE) continue;
            l = sqrt((double)gg.pos[q * 3 + 1] * gg.pos[q * 3 + 1] + (double)gg.pos[q * 3 + 2] * gg.pos[q * 3 + 2]);
            iw[q * 3 + 1] = (float)(-gg.pos[q * 3 + 1] / l); iw[q * 3 + 2] = (float)(-gg.pos[q * 3 + 2] / l);
        }
        (void)CloseupSurface_orient(&gg, iw, NULL);
        for (c2 = 0; c2 < 2; c2++) {           /* P (u 201..331 x v 0..56), Q (u 150..330 x v 64..120), 2 px */
            const double ua = c2 ? 150.0 : 201.0, va = c2 ? 64.0 : 0.0;
            const int nx = c2 ? 91 : 66, ny = 29;
            const size_t base = snv;
            svo[c2] = snv; sfo[c2] = snf;
            for (yy = 0; yy < ny; yy++)
                for (xx = 0; xx < nx; xx++) {
                    const double u = ua + 2.0 * xx, v = va + 2.0 * yy;
                    const double r = (c2 && u > 208.0) ? R - 0.5 * (u - 208.0) : R;
                    float *p = &sv[snv * 3];
                    double l = 0.0;
                    p[0] = (float)v; p[1] = (float)(r * sin(u / R)); p[2] = (float)(r * cos(u / R));
                    l = sqrt((double)p[1] * p[1] + (double)p[2] * p[2]);
                    sn[snv * 3] = 0.0f; sn[snv * 3 + 1] = (float)(p[1] / l); sn[snv * 3 + 2] = (float)(p[2] / l);
                    si[snv * 3] = 0.0f; si[snv * 3 + 1] = (float)(-p[1] / l); si[snv * 3 + 2] = (float)(-p[2] / l);
                    su[snv * 2] = u + 1000.0 * (c2 + 1); su[snv * 2 + 1] = v - 500.0;   /* own UV: shifted frame */
                    snv++;
                }
            for (yy = 0; yy + 1 < ny; yy++)
                for (xx = 0; xx + 1 < nx; xx++) {
                    int32_t a = (int32_t)(base + (size_t)(yy * nx + xx)), b = a + 1, c = a + nx, d = c + 1;
                    sf[snf * 3] = a; sf[snf * 3 + 1] = b; sf[snf * 3 + 2] = d; snf++;
                    sf[snf * 3] = a; sf[snf * 3 + 1] = d; sf[snf * 3 + 2] = c; snf++;
                }
        }
        svo[2] = snv; sfo[2] = snf;
        for (pass = 0; pass < 2; pass++) {
            CloseupPool two;
            CloseupAdoptOpts opts;
            CloseupPlacement *m = NULL;
            CloseupAdoptStats mst;
            CaTestAudit audit;
            size_t count = 0;
            int ok = 1, r0 = 0, r1 = 0;
            memset(&two, 0, sizeof two);
            two.nv = snv; two.nf = snf; two.nch = 2; two.verts = sv; two.nrm = sn; two.inward = si; two.uv = su;
            two.faces = sf; two.vofs = svo; two.fofs = sfo;
            CloseupAdoptOpts_default(&opts);
            if (pass == 1) opts.continuity_px = 0.0;
            memset(&audit, 0, sizeof audit);
            if (CloseupAdopt_round_audited(&gg, &two, NULL, &opts, &m, &count, &mst, ca_test_observer, &audit) != 0 ||
                audit.seen != 2) {
                ok = 0;
            } else {
                r0 = audit.records[0].reason; r1 = audit.records[1].reason;
                if (pass == 0 && (count != 1 || mst.deferred != 1 || !((r0 == 5 && r1 == 6) || (r0 == 6 && r1 == 5))))
                    ok = 0;
                if (pass == 1 && (count != 2 || r0 != 5 || r1 != 5)) ok = 0;
                if (audit.records[0].cont_bad != 0 || audit.records[1].cont_bad != 0) ok = 0;   /* each fine alone */
            }
            if (!ok) fails++;
            fprintf(stderr, "  closeup adopt same-round continuity (%s): %zu placed, reasons P %d Q %d, deferred %zu %s\n",
                    pass ? "gate off" : "gate on", count, r0, r1, mst.deferred, ok ? "PASS" : "FAIL");
            free(m);
        }
        CloseupObs_free(&ob);
        CloseupGrid_free(&gg);
        free(sv); free(sn); free(si); free(su); free(sf); free(iw);
    }
    {   /* ca_claim_near (the same-round claims): the nearest claim within the radius, none beyond it */
        CloseupGrid tg;
        uint16_t *cc = NULL;
        double duv = 0.0, duv_miss = -1.0;
        long hit = 0, miss = 0;
        memset(&tg, 0, sizeof tg);
        if (CloseupGrid_alloc(&tg, 0.0, 38.0, 0.0, 38.0, 2.0) != 0 ||
            (cc = (uint16_t *)calloc(tg.ni * tg.nj, sizeof *cc)) == NULL) {
            CloseupGrid_free(&tg); fails++; goto out;
        }
        cc[10 * tg.ni + 10] = 1;                    /* one claim at frame (20, 20) */
        hit = ca_claim_near(&tg, cc, (long)(10 * tg.ni + 5), 10.0, 20.0, 24.0, &duv);
        miss = ca_claim_near(&tg, cc, (long)(10 * tg.ni + 5), 10.0, 20.0, 8.0, &duv_miss);
        if (hit != (long)(10 * tg.ni + 10) || fabs(duv - 10.0) > 1e-9 || miss != -1) fails++;
        fprintf(stderr, "  closeup adopt claim search: nearest claim %ld at %.1f px (radius 24), %ld within 8 px %s\n", hit,
                duv, miss, fails ? "FAIL" : "PASS");
        free(cc);
        CloseupGrid_free(&tg);
    }
    fails += ca_layer_test();
out:
    free(pl); free(pv); free(pn); free(pin); free(puv); free(pf); free(inw);
    CloseupObs_free(&obs);
    CloseupGrid_free(&g);
    if (fails) fprintf(stderr, "FAIL closeup adopt\n");
    return fails ? 1 : 0;
}
