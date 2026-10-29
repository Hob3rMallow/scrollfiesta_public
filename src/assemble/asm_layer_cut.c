#include "asm_layer_cut.h"

#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "../common/pipeline_constants.h"
#include "../common/ves_platform.h"

#include "../common/ves_omp.h"

/* The lift keeps branches -ALC_LIFT..ALC_LIFT of the angle about the axis:
 * a witness returns to branch 0, and a path that has wound twice already
 * is far longer than any fusion. */
#define ALC_LIFT 2
#define ALC_BRANCHES (2 * ALC_LIFT + 1)
#define ALC_TWO_PI 6.28318530717958647692

typedef struct AlcGraph {
    size_t   nc, ne;
    int32_t *ea, *eb;       /* [ne] endpoints */
    int32_t *rel;           /* [ne] relation index */
    uint8_t *alive;         /* [ne] */
    int32_t *first;         /* [nc+1] CSR offsets */
    int32_t *adj;           /* [2*ne] edge ids */
    const double *phi;      /* [nc] angle about the axis */
} AlcGraph;

typedef struct AlcScratch {
    int32_t  *prev;         /* [nc*B] reaching edge; -2 unvisited, -1 the start */
    int32_t  *queue;        /* [nc*B] */
    uint16_t *depth;        /* [nc*B] */
} AlcScratch;

static double alc_wrap(double x)
{
    x = fmod(x + ALC_TWO_PI / 2.0, ALC_TWO_PI);
    if (x < 0.0) x += ALC_TWO_PI;
    return x - ALC_TWO_PI / 2.0;
}

/* Branch crossing of the step u -> v: phi[u] + wrap(phi[v] - phi[u]) equals
 * phi[v] + 2 pi m with m in {-1, 0, 1}. */
static int alc_step(const double *phi, int32_t u, int32_t v)
{
    double lifted = phi[u] + alc_wrap(phi[v] - phi[u]);
    return (int)lround((lifted - phi[v]) / ALC_TWO_PI);
}

/* Shortest join path from a to b that returns to the start's branch (net
 * winding 0), at most max_hops joins.  Writes the edge ids (b back to a) and
 * returns their count, 0 when there is none.  The scratch arrives with every
 * prev entry -2 and is returned that way. */
static int alc_search(const AlcGraph *g, int32_t a, int32_t b, int max_hops, AlcScratch *w, int32_t *path)
{
    size_t qh = 0, qt = 0;
    int32_t start = a * ALC_BRANCHES + ALC_LIFT, goal = b * ALC_BRANCHES + ALC_LIFT;
    int found = 0, n = 0;
    w->prev[start] = -1; w->depth[start] = 0; w->queue[qt++] = start;
    while (qh < qt) {
        int32_t s = w->queue[qh++];
        if (s == goal) { found = 1; break; }
        if ((int)w->depth[s] >= max_hops) continue;
        int32_t u = s / ALC_BRANCHES;
        int k = s % ALC_BRANCHES - ALC_LIFT;
        for (int32_t j = g->first[u]; j < g->first[u + 1]; j++) {
            int32_t e = g->adj[j];
            if (!g->alive[e]) continue;
            int32_t v = g->ea[e] == u ? g->eb[e] : g->ea[e];
            int k2 = k + alc_step(g->phi, u, v);
            if (k2 < -ALC_LIFT || k2 > ALC_LIFT) continue;
            int32_t t = v * ALC_BRANCHES + k2 + ALC_LIFT;
            if (w->prev[t] != -2) continue;
            w->prev[t] = e; w->depth[t] = (uint16_t)(w->depth[s] + 1); w->queue[qt++] = t;
        }
    }
    if (found) {
        int32_t t = goal;
        while (w->prev[t] != -1 && n < max_hops) {
            int32_t e = w->prev[t], v = t / ALC_BRANCHES;
            int32_t u = g->ea[e] == v ? g->eb[e] : g->ea[e];
            int k = t % ALC_BRANCHES - ALC_LIFT - alc_step(g->phi, u, v);
            path[n++] = e;
            t = u * ALC_BRANCHES + k + ALC_LIFT;
        }
    }
    for (size_t i = 0; i < qt; i++) w->prev[w->queue[i]] = -2;
    return n;
}

static int alc_cmp_i64(const void *x, const void *y)
{
    int64_t a = *(const int64_t *)x, b = *(const int64_t *)y;
    return (a > b) - (a < b);
}

static int32_t alc_root(int32_t *p, int32_t x)
{
    while (p[x] != x) { p[x] = p[p[x]]; x = p[x]; }
    return x;
}

void AsmLayerCut_defaults(AsmLayerCutOpts *o)
{
    memset(o, 0, sizeof *o);
    o->max_hops = ASM_LAYER_CUT_MAX_HOPS;
    o->min_radius = ASM_LAYER_CUT_MIN_RADIUS;
    o->min_hits = ASM_LAYER_CUT_MIN_HITS;
    o->max_passes = ASM_LAYER_CUT_MAX_PASSES;
}

int AsmLayerCut_run(AsmRun *run, const AsmLayerCutOpts *o, AsmLayerCutStats *st, FILE *ledger)
{
    double t0 = ves_clock_sec();
    memset(st, 0, sizeof *st);
    if (!run || !o || o->max_hops < 1 || o->max_hops > 1000 || o->max_passes < 0) return -1;
    if (!o->axis || !run->n_charts || !run->n_layers || !run->n_rels) { st->sec = ves_clock_sec() - t0; return 0; }
    size_t nc = run->n_charts;
    if (nc > (size_t)INT32_MAX / ALC_BRANCHES || run->n_rels > (size_t)INT32_MAX / 2) return -1;
    Arena_T arena = run->arena;
    Arena_Mark mark = Arena_save(arena);
    int rc = -1;

    /* angles about the axis; charts too near it have none and are not traversed */
    double *phi = ARENA_ALLOC(arena, nc * sizeof *phi), *rad = ARENA_ALLOC(arena, nc * sizeof *rad);
    uint8_t *ok = ARENA_CALLOC(arena, nc, 1);
    for (size_t c = 0; c < nc; c++) {
        const AsmChart *ch = run->charts + c;
        phi[c] = rad[c] = 0.0;
        if (!AsmChart_in_layout(ch)) continue;
        double f[3];
        AsmAxis_frame(o->axis, ch->centroid, f);
        double r = hypot(f[1], f[2]);
        if (!isfinite(r) || r < o->min_radius) continue;
        phi[c] = atan2(f[2], f[1]);
        rad[c] = r;
        ok[c] = 1;
    }

    /* the join graph among traversable charts */
    AlcGraph g; memset(&g, 0, sizeof g);
    g.nc = nc; g.phi = phi;
    for (size_t i = 0; i < run->n_rels; i++) {
        const AsmRelation *r = run->rels + i;
        if (r->a < 0 || r->b < 0 || (size_t)r->a >= nc || (size_t)r->b >= nc) goto done;
        if (r->a != r->b && AsmRel_is_join(r) && ok[r->a] && ok[r->b]) g.ne++;
    }
    st->joins_first = g.ne;
    g.ea = ARENA_ALLOC(arena, (g.ne ? g.ne : 1) * sizeof(int32_t));
    g.eb = ARENA_ALLOC(arena, (g.ne ? g.ne : 1) * sizeof(int32_t));
    g.rel = ARENA_ALLOC(arena, (g.ne ? g.ne : 1) * sizeof(int32_t));
    g.alive = ARENA_ALLOC(arena, (g.ne ? g.ne : 1) * sizeof(uint8_t));
    g.first = ARENA_CALLOC(arena, nc + 1, sizeof(int32_t));
    g.adj = ARENA_ALLOC(arena, (2 * g.ne ? 2 * g.ne : 1) * sizeof(int32_t));
    {
        size_t e = 0;
        for (size_t i = 0; i < run->n_rels; i++) {
            const AsmRelation *r = run->rels + i;
            if (r->a == r->b || !AsmRel_is_join(r) || !ok[r->a] || !ok[r->b]) continue;
            g.ea[e] = r->a; g.eb[e] = r->b; g.rel[e] = (int32_t)i; g.alive[e] = 1; e++;
            g.first[r->a + 1]++; g.first[r->b + 1]++;
        }
        for (size_t c = 0; c < nc; c++) g.first[c + 1] += g.first[c];
        int32_t *fill = ARENA_CALLOC(arena, nc, sizeof(int32_t));
        for (size_t k = 0; k < g.ne; k++) {
            g.adj[g.first[g.ea[k]] + fill[g.ea[k]]++] = (int32_t)k;
            g.adj[g.first[g.eb[k]] + fill[g.eb[k]]++] = (int32_t)k;
        }
    }

    /* evidence: layer pairs, unordered and deduplicated, inside one join component */
    int64_t *keys = ARENA_ALLOC(arena, run->n_layers * sizeof(int64_t));
    size_t nk = 0;
    for (size_t l = 0; l < run->n_layers; l++) {
        const AsmLayerPair *lp = run->layers + l;
        if (lp->a < 0 || lp->b < 0 || (size_t)lp->a >= nc || (size_t)lp->b >= nc || lp->a == lp->b) continue;
        if (lp->count < o->min_hits || !ok[lp->a] || !ok[lp->b]) continue;
        int32_t a = lp->a < lp->b ? lp->a : lp->b, b = lp->a < lp->b ? lp->b : lp->a;
        keys[nk++] = ((int64_t)a << 32) | (int64_t)(uint32_t)b;
    }
    if (nk > 1) qsort(keys, nk, sizeof(int64_t), alc_cmp_i64);
    {
        size_t u = 0;
        for (size_t k = 0; k < nk; k++) if (k == 0 || keys[k] != keys[k - 1]) keys[u++] = keys[k];
        nk = u;
    }
    st->pairs = nk;
    int32_t *comp = ARENA_ALLOC(arena, nc * sizeof(int32_t));
    for (size_t c = 0; c < nc; c++) comp[c] = (int32_t)c;
    for (size_t e = 0; e < g.ne; e++) {
        int32_t x = alc_root(comp, g.ea[e]), y = alc_root(comp, g.eb[e]);
        if (x != y) comp[x > y ? x : y] = x < y ? x : y;
    }
    size_t np = 0;
    for (size_t k = 0; k < nk; k++) {
        int32_t a = (int32_t)(keys[k] >> 32), b = (int32_t)(keys[k] & 0xffffffff);
        if (alc_root(comp, a) == alc_root(comp, b)) keys[np++] = keys[k];
    }
    st->pairs_joined = np;
    if (np == 0 || np > (size_t)INT_MAX) { rc = np == 0 ? 0 : -1; goto done; }

    /* per-pair state: 0 search, 1 no witness, 2 witness with its path */
    size_t hops = (size_t)o->max_hops;
    uint8_t *status = ARENA_CALLOC(arena, np, 1);
    int32_t *path = ARENA_ALLOC(arena, np * hops * sizeof(int32_t));
    int32_t *path_n = ARENA_CALLOC(arena, np, sizeof(int32_t));
    uint32_t *count = ARENA_CALLOC(arena, g.ne ? g.ne : 1, sizeof(uint32_t));
    uint8_t *nominated = ARENA_CALLOC(arena, g.ne ? g.ne : 1, 1);
    int threads = 1;
#ifdef _OPENMP
    threads = omp_get_max_threads();
    if (threads < 1) threads = 1;
#endif
    size_t states = nc * ALC_BRANCHES;
    AlcScratch *scratch = ARENA_ALLOC(arena, (size_t)threads * sizeof(AlcScratch));
    for (int t = 0; t < threads; t++) {
        scratch[t].prev = ARENA_ALLOC(arena, states * sizeof(int32_t));
        scratch[t].queue = ARENA_ALLOC(arena, states * sizeof(int32_t));
        scratch[t].depth = ARENA_ALLOC(arena, states * sizeof(uint16_t));
        for (size_t s = 0; s < states; s++) scratch[t].prev[s] = -2;
    }
    if (ledger && fputs("pass,relation,a,b,witnesses,radial_step,seam_len,rms,n_corr,normal_gap\n", ledger) < 0) goto done;
    int np_int = (int)np;
    for (int pass = 0; ; pass++) {
        /* search every pair without a current witness path, in parallel;
         * each search is deterministic, so the result does not depend on
         * the schedule */
        int p = 0;
#pragma omp parallel for schedule(dynamic, 64) num_threads(threads)
        for (p = 0; p < np_int; p++) {
            if (status[p] != 0) continue;
            int t = 0;
#ifdef _OPENMP
            t = omp_get_thread_num();
#endif
            int32_t a = (int32_t)(keys[p] >> 32), b = (int32_t)(keys[p] & 0xffffffff);
            int n = alc_search(&g, a, b, o->max_hops, scratch + t, path + (size_t)p * hops);
            path_n[p] = n;
            status[p] = (uint8_t)(n ? 3 : 1);   /* 3 = a new path, counted below */
        }
        size_t witnesses = 0;
        for (size_t q = 0; q < np; q++) {
            if (status[q] == 3) {
                status[q] = 2;
                for (int32_t k = 0; k < path_n[q]; k++) count[path[q * hops + (size_t)k]]++;
            }
            if (status[q] == 2) witnesses++;
        }
        if (pass == 0) st->witnesses_first = witnesses;
        st->witnesses_last = witnesses;
        st->passes = pass;
        if (!witnesses || pass >= o->max_passes) break;
        /* Every witness path nominates its most used join (ties: the larger
         * radial step between the two charts, then the lower relation).  A
         * bridge is the bottleneck of all the paths through it; the ordinary
         * joins beside it carry fewer, so they are never nominated by those
         * paths, and independent bridges are cut in the same pass. */
        for (size_t q = 0; q < np; q++) if (status[q] == 2) {
            int32_t best = -1;
            double step_best = 0.0;
            for (int32_t k = 0; k < path_n[q]; k++) {
                int32_t e = path[q * hops + (size_t)k];
                double step = fabs(rad[g.ea[e]] - rad[g.eb[e]]);
                if (best < 0 || count[e] > count[best] || (count[e] == count[best] && step > step_best)) { best = e; step_best = step; }
            }
            if (best >= 0) nominated[best] = 1;
        }
        for (size_t e = 0; e < g.ne; e++) {
            if (!nominated[e]) continue;
            nominated[e] = 0;
            g.alive[e] = 0;
            AsmRelation *r = run->rels + g.rel[e];
            r->flags |= ASM_REL_DROPPED;
            st->dropped++;
            if (ledger && fprintf(ledger, "%d,%d,%d,%d,%u,%.2f,%.1f,%.3f,%zu,%.3f\n", pass + 1, g.rel[e], r->a, r->b,
                                  count[e], fabs(rad[r->a] - rad[r->b]), r->seam_len, r->rms, r->n_corr, r->normal_gap) < 0) goto done;
        }
        /* a witness through a dropped join is searched again; its old path
         * no longer counts */
        for (size_t q = 0; q < np; q++) if (status[q] == 2) {
            int broken = 0;
            for (int32_t k = 0; k < path_n[q] && !broken; k++) broken = !g.alive[path[q * hops + (size_t)k]];
            if (!broken) continue;
            for (int32_t k = 0; k < path_n[q]; k++) count[path[q * hops + (size_t)k]]--;
            status[q] = 0;
        }
    }
    rc = 0;
done:
    Arena_restore(arena, mark);
    st->sec = ves_clock_sec() - t0;
    return rc;
}

/* ---- selftest ---------------------------------------------------------------- */

static float alc_dummy_uv[2];

/* charts at (radius, angle) about the z axis, all in layout */
static void alc_chart(AsmChart *c, int32_t id, double r, double ang)
{
    memset(c, 0, sizeof *c);
    c->id = id; c->uv = alc_dummy_uv; c->nv = 1; c->nf = 1; c->area3d = 1000.0;
    c->centroid[0] = 0.0; c->centroid[1] = r * sin(ang); c->centroid[2] = r * cos(ang);
}

static void alc_join(AsmRelation *r, int32_t a, int32_t b)
{
    memset(r, 0, sizeof *r);
    r->a = a; r->b = b; r->seam_len = 100.0; r->rms = 0.5; r->n_corr = 20; r->robust_w = 1.0;
}

static void alc_layer(AsmLayerPair *l, int32_t a, int32_t b)
{
    memset(l, 0, sizeof *l);
    l->a = a; l->b = b; l->count = 8; l->side = 1; l->d_median = 12.0; l->k = 1; l->nsign = 1;
}

static int alc_case(int kind, size_t *dropped, int32_t *dropped_rel)
{
    Arena_T arena = Arena_new();
    AsmRun run; memset(&run, 0, sizeof run); run.arena = arena;
    double pt[3] = { 0.0, 0.0, 0.0 }, dz[3] = { 1.0, 0.0, 0.0 };
    const AsmAxis *axis = AsmAxis_line(arena, pt, dz, -100.0, 100.0);
    size_t n = 32;
    run.charts = ARENA_CALLOC(arena, n, sizeof(AsmChart)); run.n_charts = n;
    run.rels = ARENA_CALLOC(arena, 64, sizeof(AsmRelation));
    run.layers = ARENA_CALLOC(arena, 64, sizeof(AsmLayerPair));
    if (kind == 0) {
        /* two closed rings one wrap apart, one fusion join 0-16 */
        for (size_t i = 0; i < 16; i++) {
            alc_chart(run.charts + i, (int32_t)i, 300.0, ALC_TWO_PI * (double)i / 16.0);
            alc_chart(run.charts + 16 + i, (int32_t)(16 + i), 312.0, ALC_TWO_PI * (double)i / 16.0);
        }
        for (size_t i = 0; i < 16; i++) {
            alc_join(run.rels + run.n_rels++, (int32_t)i, (int32_t)((i + 1) % 16));
            alc_join(run.rels + run.n_rels++, (int32_t)(16 + i), (int32_t)(16 + (i + 1) % 16));
            alc_layer(run.layers + run.n_layers++, (int32_t)i, (int32_t)(16 + i));
        }
        alc_join(run.rels + run.n_rels++, 0, 16);
    } else {
        /* one sheet spiralling two turns, 16 charts per turn; kind 2 adds a fusion 3-19 */
        for (size_t i = 0; i < n; i++) {
            double t = ALC_TWO_PI * (double)i / 16.0;
            alc_chart(run.charts + i, (int32_t)i, 300.0 + 12.0 * t / ALC_TWO_PI, t);
        }
        for (size_t i = 0; i + 1 < n; i++) alc_join(run.rels + run.n_rels++, (int32_t)i, (int32_t)(i + 1));
        for (size_t i = 0; i + 16 < n; i++) alc_layer(run.layers + run.n_layers++, (int32_t)i, (int32_t)(i + 16));
        if (kind == 2) alc_join(run.rels + run.n_rels++, 3, 19);
    }
    AsmLayerCutOpts o; AsmLayerCut_defaults(&o);
    o.axis = axis; o.max_hops = 24; o.min_radius = 64.0;
    AsmLayerCutStats st;
    int rc = AsmLayerCut_run(&run, &o, &st, NULL);
    *dropped = 0; *dropped_rel = -1;
    for (size_t i = 0; i < run.n_rels; i++) if (run.rels[i].flags & ASM_REL_DROPPED) { (*dropped)++; *dropped_rel = (int32_t)i; }
    if (rc == 0 && st.dropped != *dropped) rc = -1;
    Arena_dispose(&arena);
    return rc;
}

int AsmLayerCut_selftest(void)
{
    int fails = 0;
    size_t dropped = 0; int32_t rel = -1;
    /* two rings: only the fusion join (the last relation, index 32) goes */
    if (alc_case(0, &dropped, &rel) != 0 || dropped != 1 || rel != 32) {
        fprintf(stderr, "  asm_layer_cut selftest FAIL: rings dropped %zu (relation %d), want the fusion join 32\n", dropped, rel);
        fails++;
    }
    /* a true spiral: its layer pairs are one full turn apart, no witness */
    if (alc_case(1, &dropped, &rel) != 0 || dropped != 0) {
        fprintf(stderr, "  asm_layer_cut selftest FAIL: clean spiral dropped %zu joins\n", dropped);
        fails++;
    }
    /* the spiral with a fusion 3-19 (relation 31): only it goes */
    if (alc_case(2, &dropped, &rel) != 0 || dropped != 1 || rel != 31) {
        fprintf(stderr, "  asm_layer_cut selftest FAIL: fused spiral dropped %zu (relation %d), want the fusion join 31\n", dropped, rel);
        fails++;
    }
    /* no axis: the pass is off and touches nothing */
    {
        Arena_T arena = Arena_new();
        AsmRun run; memset(&run, 0, sizeof run); run.arena = arena;
        AsmLayerCutOpts o; AsmLayerCut_defaults(&o);
        AsmLayerCutStats st;
        if (AsmLayerCut_run(&run, &o, &st, NULL) != 0 || st.dropped != 0) { fprintf(stderr, "  asm_layer_cut selftest FAIL: empty run\n"); fails++; }
        Arena_dispose(&arena);
    }
    return fails;
}
