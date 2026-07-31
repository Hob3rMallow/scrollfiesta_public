#include "atlas_sheet_split.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef ASS_2PI
#define ASS_2PI 6.283185307179586476925286766559
#endif

void AtlasSheetSplitOptions_default(AtlasSheetSplitOptions *opts)
{
    if (opts == NULL) return;
    opts->min_edge_pairs = 2;
    opts->align_max_residual = 0.35;
}

/* ------------------------------------------------------------------ helpers */

typedef struct {
    int32_t lo;
    int32_t hi;
    int32_t delta;   /* level[hi] - level[lo] */
} SsTriplet;

typedef struct {
    int32_t lo;
    int32_t hi;
    int32_t delta;
    int32_t weight;  /* radial pairs voting for this component pair */
} SsEdge;

static int ss_cmp_triplet(const void *a, const void *b)
{
    const SsTriplet *x = (const SsTriplet *)a;
    const SsTriplet *y = (const SsTriplet *)b;
    if (x->lo != y->lo) return x->lo < y->lo ? -1 : 1;
    if (x->hi != y->hi) return x->hi < y->hi ? -1 : 1;
    if (x->delta != y->delta) return x->delta < y->delta ? -1 : 1;
    return 0;
}

/* Descending weight: Kruskal keeps the best-supported edges, so whatever
 * contradiction remains is carried by the weakest evidence. */
static int ss_cmp_edge_weight(const void *a, const void *b)
{
    const SsEdge *x = (const SsEdge *)a;
    const SsEdge *y = (const SsEdge *)b;
    if (x->weight != y->weight) return x->weight > y->weight ? -1 : 1;
    if (x->lo != y->lo) return x->lo < y->lo ? -1 : 1;
    if (x->hi != y->hi) return x->hi < y->hi ? -1 : 1;
    return 0;
}

static int ss_cmp_double(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    if (x < y) return -1;
    if (x > y) return 1;
    return 0;
}

static double ss_median(double *v, size_t n)
{
    if (n == 0) return 0.0;
    qsort(v, n, sizeof(double), ss_cmp_double);
    return n % 2 ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

static int32_t ss_find(int32_t *parent, int32_t x)
{
    while (parent[x] != x) {
        parent[x] = parent[parent[x]];
        x = parent[x];
    }
    return x;
}

/* ---------------------------------------------------------------- the solve */

int AtlasSheetSplit_solve(Arena_T arena,
                          const AtlasRadialOrderPair *pairs, size_t npairs,
                          const int32_t *sample_component,
                          const int32_t *sample_support,
                          const double *sample_radius,
                          const double *sample_phi,
                          size_t nsamples,
                          size_t nmesh_components,
                          int32_t target_support,
                          double pitch,
                          double sense,
                          const AtlasSheetSplitOptions *opts,
                          AtlasSheetSplitResult *out)
{
    if (arena == NULL || out == NULL || sample_component == NULL ||
        sample_radius == NULL || sample_phi == NULL)
        return -1;
    if (!(pitch > 0.0)) return -1;
    AtlasSheetSplitOptions o;
    if (opts != NULL) o = *opts; else AtlasSheetSplitOptions_default(&o);

    memset(out, 0, sizeof *out);
    size_t nmesh = nmesh_components;
    out->nmesh = nmesh;
    out->sheet_of = (int32_t *)ARENA_ALLOC(
        arena, (nmesh ? nmesh : 1) * sizeof(int32_t));
    for (size_t i = 0; i < nmesh; i++) out->sheet_of[i] = -1;

    /* --- which mesh components belong to the target support ------------- */
    int32_t *node_of = (int32_t *)ARENA_ALLOC(
        arena, (nmesh ? nmesh : 1) * sizeof(int32_t));
    for (size_t i = 0; i < nmesh; i++) node_of[i] = -1;

    size_t *nsamp = (size_t *)ARENA_CALLOC(arena, nmesh ? nmesh : 1,
                                           sizeof(size_t));
    double *rsum = (double *)ARENA_CALLOC(arena, nmesh ? nmesh : 1,
                                          sizeof(double));
    double *psum = (double *)ARENA_CALLOC(arena, nmesh ? nmesh : 1,
                                          sizeof(double));
    double *plo = (double *)ARENA_ALLOC(arena, (nmesh ? nmesh : 1) *
                                               sizeof(double));
    double *phi_hi = (double *)ARENA_ALLOC(arena, (nmesh ? nmesh : 1) *
                                                  sizeof(double));
    double *rlo = (double *)ARENA_ALLOC(arena, (nmesh ? nmesh : 1) *
                                               sizeof(double));
    double *rhi = (double *)ARENA_ALLOC(arena, (nmesh ? nmesh : 1) *
                                               sizeof(double));
    for (size_t i = 0; i < nsamples; i++) {
        int32_t mc = sample_component[i];
        if (mc < 0 || (size_t)mc >= nmesh) continue;
        if (target_support >= 0 && sample_support != NULL &&
            sample_support[i] != target_support)
            continue;
        if (nsamp[mc] == 0) {
            plo[mc] = phi_hi[mc] = sample_phi[i];
            rlo[mc] = rhi[mc] = sample_radius[i];
        } else {
            if (sample_phi[i] < plo[mc]) plo[mc] = sample_phi[i];
            if (sample_phi[i] > phi_hi[mc]) phi_hi[mc] = sample_phi[i];
            if (sample_radius[i] < rlo[mc]) rlo[mc] = sample_radius[i];
            if (sample_radius[i] > rhi[mc]) rhi[mc] = sample_radius[i];
        }
        nsamp[mc]++;
        rsum[mc] += sample_radius[i];
        psum[mc] += sample_phi[i];
    }
    size_t nnode = 0;
    for (size_t i = 0; i < nmesh; i++)
        if (nsamp[i] > 0) node_of[i] = (int32_t)nnode++;
    out->stats.components = nnode;
    if (nnode == 0) {
        out->components = NULL;
        out->islands = NULL;
        return 0;
    }

    AtlasSheetComponent *comp = (AtlasSheetComponent *)ARENA_CALLOC(
        arena, nnode, sizeof(AtlasSheetComponent));
    for (size_t i = 0; i < nmesh; i++) {
        if (node_of[i] < 0) continue;
        AtlasSheetComponent *c = &comp[node_of[i]];
        c->mesh_component = (int32_t)i;
        c->samples = nsamp[i];
        c->mean_radius = rsum[i] / (double)nsamp[i];
        c->mean_phi = psum[i] / (double)nsamp[i];
        c->phi_min = plo[i];
        c->phi_max = phi_hi[i];
        c->radius_min = rlo[i];
        c->radius_max = rhi[i];
        c->turn_span = (phi_hi[i] - plo[i]) / ASS_2PI;
        c->island = -1;
        c->local_level = 0;
        c->sheet = -1;
        if (c->turn_span >= 1.0) out->stats.spiral_components++;
        if (c->turn_span > out->stats.turn_span_max)
            out->stats.turn_span_max = c->turn_span;
        out->stats.turn_span_total += c->turn_span;
    }

    /* --- radial pairs -> ordered component-pair votes -------------------- */
    SsTriplet *tri = (SsTriplet *)ARENA_ALLOC(
        arena, (npairs ? npairs : 1) * sizeof(SsTriplet));
    size_t ntri = 0;
    for (size_t i = 0; i < npairs; i++) {
        const AtlasRadialOrderPair *p = &pairs[i];
        if (p->inner < 0 || (size_t)p->inner >= nsamples) continue;
        if (p->outer < 0 || (size_t)p->outer >= nsamples) continue;
        int32_t mi = sample_component[p->inner];
        int32_t mo = sample_component[p->outer];
        if (mi < 0 || (size_t)mi >= nmesh || mo < 0 || (size_t)mo >= nmesh)
            continue;
        int32_t a = node_of[mi], b = node_of[mo];
        if (a < 0 || b < 0) continue;   /* outside the target support */
        out->stats.pairs_used++;
        if (a == b) {
            /* One connected piece of surface reaching a whole turn around to
             * touch itself.  A single sheet index cannot describe it; report
             * rather than silently mislabel. */
            out->stats.pairs_intra_component++;
            comp[a].pairs++;
            comp[a].pairs_internal++;
            continue;
        }
        comp[a].pairs++;
        comp[b].pairs++;
        SsTriplet *t = &tri[ntri++];
        if (a < b) { t->lo = a; t->hi = b; t->delta = p->turns; }
        else       { t->lo = b; t->hi = a; t->delta = -p->turns; }
    }
    for (size_t i = 0; i < nnode; i++)
        if (comp[i].pairs > 0) out->stats.components_with_pairs++;

    /* --- collapse to one modal edge per component pair ------------------- */
    qsort(tri, ntri, sizeof(SsTriplet), ss_cmp_triplet);
    SsEdge *edge = (SsEdge *)ARENA_ALLOC(
        arena, (ntri ? ntri : 1) * sizeof(SsEdge));
    size_t nedge = 0;
    {
        size_t i = 0;
        while (i < ntri) {
            size_t j = i;
            while (j < ntri && tri[j].lo == tri[i].lo && tri[j].hi == tri[i].hi)
                j++;
            /* modal delta within the group; the run is already sorted by it */
            int32_t best_delta = tri[i].delta;
            size_t best_n = 0, total = j - i, k = i;
            while (k < j) {
                size_t m = k;
                while (m < j && tri[m].delta == tri[k].delta) m++;
                if (m - k > best_n) { best_n = m - k; best_delta = tri[k].delta; }
                k = m;
            }
            out->stats.edges++;
            if ((int)total < o.min_edge_pairs) {
                out->stats.edges_dropped_weak++;
            } else {
                SsEdge *e = &edge[nedge++];
                e->lo = tri[i].lo;
                e->hi = tri[i].hi;
                e->delta = best_delta;
                e->weight = (int32_t)total;
            }
            i = j;
        }
    }

    /* --- max-weight spanning forest -------------------------------------- */
    qsort(edge, nedge, sizeof(SsEdge), ss_cmp_edge_weight);
    int32_t *parent = (int32_t *)ARENA_ALLOC(arena, nnode * sizeof(int32_t));
    for (size_t i = 0; i < nnode; i++) parent[i] = (int32_t)i;
    uint8_t *in_forest = (uint8_t *)ARENA_CALLOC(arena, nedge ? nedge : 1,
                                                 sizeof(uint8_t));
    for (size_t i = 0; i < nedge; i++) {
        int32_t ra = ss_find(parent, edge[i].lo);
        int32_t rb = ss_find(parent, edge[i].hi);
        if (ra == rb) continue;
        parent[ra] = rb;
        in_forest[i] = 1;
        out->stats.edges_in_forest++;
    }

    /* --- levels by BFS over the forest ----------------------------------- */
    size_t *deg = (size_t *)ARENA_CALLOC(arena, nnode + 1, sizeof(size_t));
    for (size_t i = 0; i < nedge; i++) {
        if (!in_forest[i]) continue;
        deg[(size_t)edge[i].lo + 1]++;
        deg[(size_t)edge[i].hi + 1]++;
    }
    for (size_t i = 0; i < nnode; i++) deg[i + 1] += deg[i];
    size_t nadj = deg[nnode];
    int32_t *adj = (int32_t *)ARENA_ALLOC(arena, (nadj ? nadj : 1) *
                                                 sizeof(int32_t));
    int32_t *adj_delta = (int32_t *)ARENA_ALLOC(arena, (nadj ? nadj : 1) *
                                                       sizeof(int32_t));
    size_t *cursor = (size_t *)ARENA_ALLOC(arena, nnode * sizeof(size_t));
    memcpy(cursor, deg, nnode * sizeof(size_t));
    for (size_t i = 0; i < nedge; i++) {
        if (!in_forest[i]) continue;
        size_t sa = cursor[edge[i].lo]++;
        adj[sa] = edge[i].hi;  adj_delta[sa] = edge[i].delta;
        size_t sb = cursor[edge[i].hi]++;
        adj[sb] = edge[i].lo;  adj_delta[sb] = -edge[i].delta;
    }

    int32_t *island = (int32_t *)ARENA_ALLOC(arena, nnode * sizeof(int32_t));
    int32_t *level = (int32_t *)ARENA_CALLOC(arena, nnode, sizeof(int32_t));
    for (size_t i = 0; i < nnode; i++) island[i] = -1;
    int32_t *queue = (int32_t *)ARENA_ALLOC(arena, nnode * sizeof(int32_t));
    int32_t nisl = 0;
    for (size_t s = 0; s < nnode; s++) {
        if (island[s] >= 0) continue;
        size_t head = 0, tail = 0;
        island[s] = nisl;
        level[s] = 0;
        queue[tail++] = (int32_t)s;
        while (head < tail) {
            size_t v = (size_t)queue[head++];
            for (size_t k = deg[v]; k < deg[v + 1]; k++) {
                size_t w = (size_t)adj[k];
                if (island[w] >= 0) continue;
                island[w] = nisl;
                level[w] = level[v] + adj_delta[k];
                queue[tail++] = (int32_t)w;
            }
        }
        nisl++;
    }
    out->stats.islands = (size_t)nisl;

    for (size_t i = 0; i < nedge; i++) {
        if (in_forest[i]) continue;
        if (level[edge[i].hi] - level[edge[i].lo] == edge[i].delta)
            out->stats.edges_consistent++;
        else
            out->stats.edges_contradictory++;
    }

    /* --- align the islands against each other by radial intercept -------- */
    AtlasSheetIsland *isl = (AtlasSheetIsland *)ARENA_CALLOC(
        arena, (size_t)nisl, sizeof(AtlasSheetIsland));
    double *scratch = (double *)ARENA_ALLOC(arena, nnode * sizeof(double));
    for (int32_t g = 0; g < nisl; g++) {
        size_t n = 0;
        AtlasSheetIsland *is = &isl[g];
        is->island = g;
        is->level_min = 0x7fffffff;
        is->level_max = -0x7fffffff;
        double rsum_isl = 0.0;
        size_t wsum = 0;
        for (size_t i = 0; i < nnode; i++) {
            if (island[i] != g) continue;
            scratch[n++] = comp[i].mean_radius - pitch * (double)level[i];
            is->components++;
            is->samples += comp[i].samples;
            if (level[i] < is->level_min) is->level_min = level[i];
            if (level[i] > is->level_max) is->level_max = level[i];
            rsum_isl += comp[i].mean_radius * (double)comp[i].samples;
            wsum += comp[i].samples;
        }
        is->intercept = ss_median(scratch, n);
        is->mean_radius = wsum ? rsum_isl / (double)wsum : 0.0;
    }
    for (int32_t g = 0; g < nisl; g++)
        if (isl[g].components == 1) out->stats.singleton_islands++;

    int32_t ref = 0;
    for (int32_t g = 1; g < nisl; g++)
        if (isl[g].samples > isl[ref].samples) ref = g;
    double sum_sq = 0.0;
    for (int32_t g = 0; g < nisl; g++) {
        double shift = (isl[g].intercept - isl[ref].intercept) / pitch;
        double rounded = floor(shift + 0.5);
        isl[g].offset = (int32_t)rounded;
        isl[g].residual = shift - rounded;
        sum_sq += isl[g].residual * isl[g].residual;
        if (fabs(isl[g].residual) > fabs(out->stats.align_residual_max))
            out->stats.align_residual_max = isl[g].residual;
        if (fabs(isl[g].residual) > o.align_max_residual)
            out->stats.align_unreliable++;
    }
    out->stats.align_residual_rms = nisl ? sqrt(sum_sq / (double)nisl) : 0.0;

    int32_t gmin = 0x7fffffff, gmax = -0x7fffffff;
    for (size_t i = 0; i < nnode; i++) {
        int32_t g = level[i] + isl[island[i]].offset;
        comp[i].local_level = level[i];
        comp[i].island = island[i];
        comp[i].sheet = g;
        if (g < gmin) gmin = g;
        if (g > gmax) gmax = g;
    }
    for (size_t i = 0; i < nnode; i++) {
        comp[i].sheet -= gmin;
        out->sheet_of[comp[i].mesh_component] = comp[i].sheet;
    }
    out->stats.sheet_min = 0;
    out->stats.sheet_max = gmax - gmin;
    {
        size_t span = (size_t)(gmax - gmin) + 1;
        uint8_t *seen = (uint8_t *)ARENA_CALLOC(arena, span, sizeof(uint8_t));
        for (size_t i = 0; i < nnode; i++) seen[comp[i].sheet] = 1;
        for (size_t i = 0; i < span; i++) if (seen[i]) out->stats.sheets++;
    }

    /* --- independent cross-check against phi ----------------------------- */
    {
        size_t best = 0;
        for (size_t i = 1; i < nnode; i++)
            if (comp[i].samples > comp[best].samples) best = i;
        double phi_ref = comp[best].mean_phi;
        int32_t sheet_ref = comp[best].sheet;
        for (size_t i = 0; i < nnode; i++) {
            double expect = (comp[i].mean_phi - phi_ref) / ASS_2PI * sense;
            int32_t e = (int32_t)floor(expect + 0.5);
            if (e == comp[i].sheet - sheet_ref) out->stats.phi_agree++;
            else out->stats.phi_disagree++;
        }
    }

    out->components = comp;
    out->ncomponents = nnode;
    out->islands = isl;
    out->nislands = (size_t)nisl;
    return 0;
}

/* ------------------------------------------------------------------ selftest */

#define SS_CHECK(cond, msg)                                                   \
    do {                                                                      \
        if (!(cond)) {                                                        \
            fprintf(stderr, "[atlas_sheet_split] FAIL: %s\n", (msg));         \
            failures++;                                                       \
        }                                                                     \
    } while (0)

/*
 * Four wraps of a spiral, chopped into eight fragments.  Fragments 0..3 sit at
 * angles that give them shared rays (one island); 4..7 are the same four wraps
 * on the far side with no ray in common with the first group (a second island),
 * so the only thing that can line the two groups up is the radial intercept.
 * A ninth fragment is joined to the first group by a SINGLE spurious pair --
 * an accidental bridge -- which must be dropped, not believed.
 */
int AtlasSheetSplit_selftest(void)
{
    int failures = 0;
    Arena_T arena = Arena_new();

    enum { NCOMP = 9, NSAMP = NCOMP * 4 };
    const double pitch = 12.0;
    const double r0 = 200.0;

    int32_t sample_component[NSAMP];
    int32_t sample_support[NSAMP];
    double  sample_radius[NSAMP];
    double  sample_phi[NSAMP];
    /* wrap index per fragment: two groups of four wraps, plus the bridge */
    const int wrap_of[NCOMP] = { 0, 1, 2, 3, 0, 1, 2, 3, 2 };

    for (int c = 0; c < NCOMP; c++)
        for (int k = 0; k < 4; k++) {
            int i = c * 4 + k;
            sample_component[i] = c;
            sample_support[i] = 7;
            sample_radius[i] = r0 + pitch * wrap_of[c] + 0.05 * k;
            sample_phi[i] = ASS_2PI * wrap_of[c] + 0.01 * k;
        }

    /* Island A: fragments 0-1-2-3 chained by well-supported pairs.
     * Island B: fragments 4-5-6-7, likewise, but never joined to A.
     * Bridge:   one lone pair claiming 8 is adjacent to 0. */
    AtlasRadialOrderPair pairs[64];
    size_t np = 0;
    for (int c = 0; c < 3; c++)
        for (int rep = 0; rep < 5; rep++) {
            AtlasRadialOrderPair *p = &pairs[np++];
            memset(p, 0, sizeof *p);
            p->inner = c * 4 + rep % 4;
            p->outer = (c + 1) * 4 + rep % 4;
            p->turns = 1;
        }
    for (int c = 4; c < 7; c++)
        for (int rep = 0; rep < 5; rep++) {
            AtlasRadialOrderPair *p = &pairs[np++];
            memset(p, 0, sizeof *p);
            p->inner = c * 4 + rep % 4;
            p->outer = (c + 1) * 4 + rep % 4;
            p->turns = 1;
        }
    {   /* the accidental bridge: one pair, claiming fragment 8 is a wrap out */
        AtlasRadialOrderPair *p = &pairs[np++];
        memset(p, 0, sizeof *p);
        p->inner = 0;
        p->outer = 8 * 4;
        p->turns = 1;
    }

    AtlasSheetSplitOptions opts;
    AtlasSheetSplitOptions_default(&opts);
    AtlasSheetSplitResult res;
    int rc = AtlasSheetSplit_solve(arena, pairs, np, sample_component,
                                   sample_support, sample_radius, sample_phi,
                                   (size_t)NSAMP, (size_t)NCOMP, 7, pitch,
                                   1.0, &opts, &res);
    SS_CHECK(rc == 0, "solve returns 0");
    SS_CHECK(res.ncomponents == NCOMP, "all fragments in the target support");
    SS_CHECK(res.stats.edges_dropped_weak == 1,
             "the single-pair accidental bridge is dropped");
    SS_CHECK(res.stats.pairs_intra_component == 0,
             "no fragment spans a whole turn");
    SS_CHECK(res.stats.edges_contradictory == 0,
             "the surviving evidence is self-consistent");
    /* Two chained islands plus the isolated bridge fragment. */
    SS_CHECK(res.stats.islands == 3, "two chains and one orphan");
    SS_CHECK(res.stats.sheets == 4, "four sheets recovered");
    SS_CHECK(res.stats.sheet_max == 3, "sheets are 0..3");
    SS_CHECK(fabs(res.stats.align_residual_max) < 0.05,
             "islands align on the radial intercept without a coin flip");
    for (size_t i = 0; i < res.ncomponents; i++) {
        int32_t mc = res.components[i].mesh_component;
        SS_CHECK(res.components[i].sheet == wrap_of[mc],
                 "each fragment lands on its own wrap");
    }

    /* A support component that is not the target contributes nothing. */
    AtlasSheetSplitResult empty;
    rc = AtlasSheetSplit_solve(arena, pairs, np, sample_component,
                               sample_support, sample_radius, sample_phi,
                               (size_t)NSAMP, (size_t)NCOMP, 99, pitch, 1.0,
                               &opts, &empty);
    SS_CHECK(rc == 0 && empty.ncomponents == 0,
             "an absent support component yields nothing, not a crash");

    /* Zero input must not crash either. */
    AtlasSheetSplitResult none;
    rc = AtlasSheetSplit_solve(arena, NULL, 0, sample_component,
                               sample_support, sample_radius, sample_phi,
                               (size_t)NSAMP, (size_t)NCOMP, 7, pitch, 1.0,
                               &opts, &none);
    SS_CHECK(rc == 0 && none.stats.islands == (size_t)NCOMP,
             "no evidence leaves every fragment its own island");

    Arena_dispose(&arena);
    fprintf(stderr, "[atlas_sheet_split] selftest %s (%d failures)\n",
            failures == 0 ? "PASSED" : "FAILED", failures);
    return failures;
}
