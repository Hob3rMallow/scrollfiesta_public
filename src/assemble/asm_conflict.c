/* asm_conflict.c -- contradiction audit + cleaning rounds. */
#include "asm_conflict.h"

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/pipeline_constants.h"
#include "../common/union_find.h"
#include "../common/ves_platform.h"

#include "../common/ves_omp.h"

void AsmConflict_default_opts(AsmConflictOpts *o)
{
    o->cell = ASM_CELL_VOX;
    o->contra_lo = ASM_CONTRA_LO;
    o->contra_hi = ASM_CONTRA_HI;
    o->rounds = ASM_CLEAN_ROUNDS;
    o->max_path = 12;
    o->drop_frac = 0.05;
    o->extras_ratio = 0.5;
    o->diag_pairs_csv = NULL;
}

/* (min, max) chart pair key for the diagnostic lookups of the audit ledger */
typedef struct AcPairKey { int64_t key; double value; uint32_t flags; } AcPairKey;
static int64_t ac_pair_key(int32_t a, int32_t b) { return a < b ? ((int64_t)a << 32) | (uint32_t)b : ((int64_t)b << 32) | (uint32_t)a; }
static int ac_cmp_pair_key(const void *x, const void *y)
{ const AcPairKey *p = x, *q = y; return (p->key > q->key) - (p->key < q->key); }
static const AcPairKey *ac_find_pair_key(const AcPairKey *keys, size_t n, int64_t key)
{
    size_t lo = 0, hi = n;
    while (lo < hi) { size_t mid = lo + (hi - lo) / 2; if (keys[mid].key < key) lo = mid + 1; else hi = mid; }
    return lo < n && keys[lo].key == key ? &keys[lo] : NULL;
}

/* ---- raster ------------------------------------------------------------------ */

typedef struct AcClaim {
    int64_t key;       /* (component << 40) | (cy+2^19) << 20 | (cx+2^19) */
    int32_t chart;
    int32_t face;      /* the chart face the claim was interpolated on */
    float   p[3];
    float   n[3];
} AcClaim;

static int64_t ac_key(int32_t component, int32_t cx, int32_t cy)
{
    return ((int64_t)component << 40) | ((int64_t)(cy + (1 << 19)) << 20) | (int64_t)(cx + (1 << 19));
}

static int ac_cmp_claim(const void *x, const void *y)
{
    const AcClaim *a = x, *b = y;
    if (a->key != b->key) return a->key < b->key ? -1 : 1;
    return (a->chart > b->chart) - (a->chart < b->chart);
}

static void ac_pose(const AsmChart *c, double u, double v, double *gx, double *gy)
{
    if (c->flags & ASM_CHART_MIRROR) u = -u;
    double ct = cos(c->pose_theta), st = sin(c->pose_theta);
    *gx = ct * u - st * v + c->pose_x;
    *gy = st * u + ct * v + c->pose_y;
}

static int av_cmp_double_ac(const void *x, const void *y) { double a = *(const double *)x, b = *(const double *)y; return (a > b) - (a < b); }

typedef struct AcBuf { Arena_T arena; AcClaim *c; size_t n, cap; } AcBuf;

static void ac_push(AcBuf *b, const AcClaim *cl)
{
    if (b->n == b->cap) {
        size_t ncap = b->cap ? b->cap * 2 : (1u << 16);
        AcClaim *na = ARENA_ALLOC(b->arena, ncap * sizeof(AcClaim));
        if (b->n) memcpy(na, b->c, b->n * sizeof(AcClaim));
        b->c = na; b->cap = ncap;
    }
    b->c[b->n++] = *cl;
}

/* Every vertex's layout point with exactly AsmChart_point's arithmetic, once
 * per chart: evaluated per face corner it cost three cos/sin pairs per face,
 * about 2.2e9 pairs per audit on the 21x21x21. xy holds 2*nv doubles. */
static void ac_chart_points(const AsmChart *c, double *xy)
{
    if (c->placed_uv) { for (size_t v = 0; v < 2*c->nv; v++) xy[v] = c->placed_uv[v]; return; }
    double ct = cos(c->pose_theta), sn = sin(c->pose_theta);
    for (size_t v = 0; v < c->nv; v++) {
        double u = c->uv[2*v], w = c->uv[2*v+1];
        if (c->flags & ASM_CHART_MIRROR) u = -u;
        xy[2*v] = ct*u-sn*w+c->pose_x; xy[2*v+1] = sn*u+ct*w+c->pose_y;
    }
}

/* rasterize one chart: cells whose centres fall inside a face. xy is the
 * caller's scratch for 2*nv doubles (see ac_chart_points). */
static void ac_raster_chart(const AsmChart *c, double cell, AcBuf *buf, double *xy)
{
    ac_chart_points(c, xy);
    for (size_t f = 0; f < c->nf; f++) {
        const int32_t *fv = &c->faces[f*3];
        double gx[3], gy[3];
        for (int k = 0; k < 3; k++) {gx[k]=xy[2*(size_t)fv[k]];gy[k]=xy[2*(size_t)fv[k]+1];}
        double minx = fmin(gx[0], fmin(gx[1], gx[2])), maxx = fmax(gx[0], fmax(gx[1], gx[2]));
        double miny = fmin(gy[0], fmin(gy[1], gy[2])), maxy = fmax(gy[0], fmax(gy[1], gy[2]));
        int cx0 = (int)floor(minx / cell), cx1 = (int)floor(maxx / cell);
        int cy0 = (int)floor(miny / cell), cy1 = (int)floor(maxy / cell);
        double area = (gx[1]-gx[0])*(gy[2]-gy[0]) - (gy[1]-gy[0])*(gx[2]-gx[0]);
        if (fabs(area) < 1e-9) continue;
        if (!(minx > -1e7 && maxx < 1e7 && miny > -1e7 && maxy < 1e7) || (cx1 - cx0) > 4096 || (cy1 - cy0) > 4096)
            continue;   /* a face spanning thousands of cells is a broken pose, not material */
        for (int cy = cy0; cy <= cy1; cy++) {
            double py = (cy + 0.5) * cell;
            for (int cx = cx0; cx <= cx1; cx++) {
                double px = (cx + 0.5) * cell;
                /* barycentric */
                double l0 = ((gx[1]-px)*(gy[2]-py) - (gy[1]-py)*(gx[2]-px)) / area;
                double l1 = ((gx[2]-px)*(gy[0]-py) - (gy[2]-py)*(gx[0]-px)) / area;
                double l2 = 1.0 - l0 - l1;
                if (l0 < -1e-6 || l1 < -1e-6 || l2 < -1e-6) continue;
                AcClaim cl;
                cl.key = ac_key(c->component, cx, cy);
                cl.chart = c->id;
                cl.face = (int32_t)f;
                for (int d = 0; d < 3; d++) {
                    cl.p[d] = (float)(l0 * c->xyz[(size_t)fv[0]*3+(size_t)d] + l1 * c->xyz[(size_t)fv[1]*3+(size_t)d] + l2 * c->xyz[(size_t)fv[2]*3+(size_t)d]);
                    cl.n[d] = (float)(l0 * c->nrm[(size_t)fv[0]*3+(size_t)d] + l1 * c->nrm[(size_t)fv[1]*3+(size_t)d] + l2 * c->nrm[(size_t)fv[2]*3+(size_t)d]);
                }
                float nl = sqrtf(cl.n[0]*cl.n[0] + cl.n[1]*cl.n[1] + cl.n[2]*cl.n[2]);
                if (nl > 0) { cl.n[0] /= nl; cl.n[1] /= nl; cl.n[2] /= nl; }
                ac_push(buf, &cl);
            }
        }
    }
}

/* Rasterize every placed chart in parallel into per-thread buffers, then
 * concatenate the claims in chart order: exactly the buffer the serial loop
 * built, so the sort that follows sees the same input and ties keep their
 * order. out->arena receives the concatenated claims. */
typedef struct AcSpan { int32_t thread; size_t start, count; } AcSpan;
static int ac_raster_all(const AsmRun *run, double cell, AcBuf *out)
{
    size_t nc = run->n_charts, maxnv = 0;
    if (nc > (size_t)INT_MAX) return -1;
    for (size_t i = 0; i < nc; i++) if (run->charts[i].nv > maxnv) maxnv = run->charts[i].nv;
    int threads = 1;
#ifdef _OPENMP
    threads = omp_get_max_threads(); if (threads < 1) threads = 1;
#endif
    AcBuf *tb = calloc((size_t)threads, sizeof *tb); double **xy = calloc((size_t)threads, sizeof *xy);
    AcSpan *span = calloc(nc ? nc : 1, sizeof *span);
    if (!tb || !xy || !span) { free(tb); free(xy); free(span); return -1; }
    for (int t = 0; t < threads; t++) { tb[t].arena = Arena_new(); xy[t] = ARENA_ALLOC(tb[t].arena, (maxnv ? 2*maxnv : 1) * sizeof(double)); }
    int n = (int)nc, i = 0;
#pragma omp parallel for schedule(dynamic,16) num_threads(threads) if(n >= 256)
    for (i = 0; i < n; i++) {
        const AsmChart *c = &run->charts[i];
        if (!AsmChart_in_layout(c) || !c->placed || c->component < 0) continue;
        int t = 0;
#ifdef _OPENMP
        t = omp_get_thread_num();
#endif
        span[i].thread = t; span[i].start = tb[t].n;
        ac_raster_chart(c, cell, &tb[t], xy[t]);
        span[i].count = tb[t].n - span[i].start;
    }
    size_t total = 0;
    for (size_t k = 0; k < nc; k++) total += span[k].count;
    out->c = ARENA_ALLOC(out->arena, (total ? total : 1) * sizeof(AcClaim)); out->n = out->cap = total;
    size_t at = 0;
    for (size_t k = 0; k < nc; k++) if (span[k].count) {
        memcpy(out->c + at, tb[span[k].thread].c + span[k].start, span[k].count * sizeof(AcClaim)); at += span[k].count;
    }
    for (int t = 0; t < threads; t++) Arena_dispose(&tb[t].arena);
    free(tb); free(xy); free(span);
    return 0;
}

static double ac_smoothstep(double x)
{
    if (x <= 0.0) return 0.0;
    if (x >= 1.0) return 1.0;
    return x * x * (3.0 - 2.0 * x);
}

typedef struct AcPair { int32_t a, b; double mass; double dn_sum; size_t cells; } AcPair;

static int ac_cmp_pair(const void *x, const void *y)
{
    const AcPair *p = x, *q = y;
    if (p->a != q->a) return p->a < q->a ? -1 : 1;
    return (p->b > q->b) - (p->b < q->b);
}

/* BFS shortest path over active relations between chart a and b; returns the
 * number of relations on the path (0 = none within max_path) and writes the
 * relation indices into path[]. */
/* prev_rel[] must arrive all -2 (unvisited) and is returned that way: only the
 * charts this search enqueued are reset, so a call costs what it explores
 * rather than the chart count.  Resetting every chart per call made one audit
 * 1.78 M pairs x 97,648 charts of stores on the 21x21x21. */
static int ac_witness_path(const AsmRun *run, const int32_t *deg, const int32_t *adj,
                           int32_t a, int32_t b, int max_path,
                           int32_t *prev_rel, int32_t *queue, int32_t *path)
{
    size_t nc = run->n_charts;
    size_t qh = 0, qt = 0;
    prev_rel[a] = -1;
    queue[qt++] = a;
    int32_t *depth = queue + nc;                          /* second half of the scratch */
    depth[a] = 0;
    while (qh < qt) {
        int32_t u = queue[qh++];
        if (u == b) break;
        if (depth[u] >= max_path) continue;
        for (int32_t k = deg[u]; k < deg[u+1]; k++) {
            int32_t r = adj[k];
            const AsmRelation *R = &run->rels[r];
            int32_t v = (R->a == u) ? R->b : R->a;
            if (prev_rel[v] != -2) continue;
            prev_rel[v] = r;
            depth[v] = depth[u] + 1;
            queue[qt++] = v;
        }
    }
    int n = 0;
    if (prev_rel[b] != -2) {
        int32_t v = b;
        while (v != a) {
            int32_t r = prev_rel[v];
            path[n++] = r;
            const AsmRelation *R = &run->rels[r];
            v = (R->a == v) ? R->b : R->a;
        }
    }
    for (size_t k = 0; k < qt; k++) prev_rel[queue[k]] = -2;
    return n;
}

static int ac_rel_active(const AsmRelation *R)
{
    return (R->flags & (ASM_REL_DROPPED | ASM_REL_CONTACT | ASM_REL_SWITCHED | ASM_REL_PLACEMENT_ONLY)) == 0;
}

int AsmConflict_audit(AsmRun *run, const AsmConflictOpts *o,
                      double *out_chart_mass, AsmConflictReport *rep)
{
    memset(rep, 0, sizeof *rep);
    Arena_T arena = run->arena;
    Arena_Mark mark = Arena_save(arena);
    AcBuf buf; memset(&buf, 0, sizeof buf); buf.arena = arena;
    if (ac_raster_all(run, o->cell, &buf)) { Arena_restore(arena, mark); return -1; }
    qsort(buf.c, buf.n, sizeof(AcClaim), ac_cmp_claim);
    for (size_t r = 0; r < run->n_rels; r++) run->rels[r].flag_mass = 0.0;
    if (out_chart_mass) memset(out_chart_mass, 0, run->n_charts * sizeof(double));

    /* per-cell contradiction, per-pair aggregation */
    size_t cap_pairs = 1024, n_pairs = 0;
    AcPair *pairs = ARENA_ALLOC(arena, cap_pairs * sizeof(AcPair));
    size_t cap_cells = 1024, n_cells = 0;
    AsmConflictCell *ccells = ARENA_ALLOC(arena, cap_cells * sizeof(AsmConflictCell));
    double cell_area = o->cell * o->cell;
    for (size_t i = 0; i < buf.n;) {
        size_t j = i;
        while (j < buf.n && buf.c[j].key == buf.c[i].key) j++;
        rep->cells++;
        /* distinct charts in the cell */
        size_t nd = 0;
        for (size_t k = i; k < j; k++) if (k == i || buf.c[k].chart != buf.c[k-1].chart) nd++;
        if (nd >= 2) {
            rep->cells_multi++;
            double wmax = 0.0, dn_at_max = 0.0;
            for (size_t k = i; k < j; k++) {
                for (size_t l = k + 1; l < j; l++) {
                    if (buf.c[k].chart == buf.c[l].chart) continue;
                    const AcClaim *A = &buf.c[k], *B = &buf.c[l];
                    double nx = A->n[0] + B->n[0], ny = A->n[1] + B->n[1], nz = A->n[2] + B->n[2];
                    double nl = sqrt(nx*nx + ny*ny + nz*nz);
                    double dn;
                    if (nl > 1e-6) {
                        dn = fabs(((A->p[0]-B->p[0])*nx + (A->p[1]-B->p[1])*ny + (A->p[2]-B->p[2])*nz) / nl);
                    } else {
                        /* anti-parallel normals: use the distance itself */
                        dn = sqrt((A->p[0]-B->p[0])*(A->p[0]-B->p[0]) + (A->p[1]-B->p[1])*(A->p[1]-B->p[1]) + (A->p[2]-B->p[2])*(A->p[2]-B->p[2]));
                    }
                    double dl = run->chart_layer_d ? run->chart_layer_d[A->chart] : 0.0;
                    if (dl <= 0.0 && run->chart_layer_d) dl = run->chart_layer_d[B->chart];
                    if (dl <= 0.0) dl = run->layer_d_global > 0.0 ? run->layer_d_global : ASM_LAYER_PROBE_VOX * 0.4;
                    double x = (dn / dl - o->contra_lo) / (o->contra_hi - o->contra_lo);
                    double w = ac_smoothstep(x);
                    if (w > 0.0) {
                        int32_t a = A->chart < B->chart ? A->chart : B->chart;
                        int32_t b = A->chart < B->chart ? B->chart : A->chart;
                        /* find or append the pair (linear in the small local set is fine; global sort later) */
                        if (n_pairs == cap_pairs) {
                            AcPair *np = ARENA_ALLOC(arena, cap_pairs * 2 * sizeof(AcPair));
                            memcpy(np, pairs, n_pairs * sizeof(AcPair)); pairs = np; cap_pairs *= 2;
                        }
                        pairs[n_pairs].a = a; pairs[n_pairs].b = b; pairs[n_pairs].mass = w * cell_area;
                        pairs[n_pairs].dn_sum = dn; pairs[n_pairs].cells = 1; n_pairs++;
                    }
                    if (w > wmax) { wmax = w; dn_at_max = dn; }
                }
            }
            if (wmax > 0.0) {
                rep->contra_area += wmax * cell_area;
                if (wmax > 0.5) rep->cells_contra++;
                if (n_cells == cap_cells) {
                    AsmConflictCell *nc2 = ARENA_ALLOC(arena, cap_cells * 2 * sizeof(AsmConflictCell));
                    memcpy(nc2, ccells, n_cells * sizeof(AsmConflictCell)); ccells = nc2; cap_cells *= 2;
                }
                int64_t key = buf.c[i].key;
                ccells[n_cells].component = (int32_t)(key >> 40);
                ccells[n_cells].cy = (int32_t)((key >> 20) & 0xFFFFF) - (1 << 19);
                ccells[n_cells].cx = (int32_t)(key & 0xFFFFF) - (1 << 19);
                ccells[n_cells].w = (float)wmax; ccells[n_cells].dn = (float)dn_at_max;
                n_cells++;
            }
        }
        i = j;
    }
    rep->covered_area = (double)rep->cells * cell_area;
    rep->score = rep->covered_area - rep->contra_area;
    /* aggregate pairs */
    qsort(pairs, n_pairs, sizeof(AcPair), ac_cmp_pair);
    size_t na = 0;
    for (size_t i = 0; i < n_pairs;) {
        size_t j = i;
        AcPair acc = pairs[i]; acc.mass = 0.0; acc.dn_sum = 0.0; acc.cells = 0;
        while (j < n_pairs && pairs[j].a == pairs[i].a && pairs[j].b == pairs[i].b) { acc.mass += pairs[j].mass; acc.dn_sum += pairs[j].dn_sum; acc.cells += pairs[j].cells; j++; }
        pairs[na++] = acc;
        i = j;
    }
    n_pairs = na;
    rep->pairs = n_pairs;

    /* witness attribution over the active relation graph */
    size_t nc = run->n_charts;
    int32_t *deg = ARENA_CALLOC(arena, nc + 1, sizeof(int32_t));
    for (size_t r = 0; r < run->n_rels; r++) if (ac_rel_active(&run->rels[r])) { deg[run->rels[r].a + 1]++; deg[run->rels[r].b + 1]++; }
    for (size_t i = 0; i < nc; i++) deg[i+1] += deg[i];
    int32_t *adj = ARENA_ALLOC(arena, ((size_t)deg[nc] + 1) * sizeof(int32_t));
    int32_t *fill = ARENA_CALLOC(arena, nc, sizeof(int32_t));
    for (size_t r = 0; r < run->n_rels; r++) {
        if (!ac_rel_active(&run->rels[r])) continue;
        adj[deg[run->rels[r].a] + fill[run->rels[r].a]++] = (int32_t)r;
        adj[deg[run->rels[r].b] + fill[run->rels[r].b]++] = (int32_t)r;
    }
    /* Witness searches are independent: run them in parallel with private
     * scratch and keep each pair's path, then attribute mass serially in pair
     * order below, so every sum accumulates exactly as a serial audit would. */
    if (n_pairs > (size_t)INT_MAX || o->max_path < 0) { Arena_restore(arena, mark); return -1; }
    int threads = 1;
#ifdef _OPENMP
    threads = omp_get_max_threads(); if (threads < 1) threads = 1;
#endif
    size_t stride = (size_t)o->max_path + 1;
    int32_t *scratch_prev = ARENA_ALLOC(arena, (size_t)threads * (nc ? nc : 1) * sizeof(int32_t));
    int32_t *scratch_queue = ARENA_ALLOC(arena, (size_t)threads * 2 * (nc ? nc : 1) * sizeof(int32_t));
    for (size_t i = 0; i < (size_t)threads * nc; i++) scratch_prev[i] = -2;
    int32_t *paths = ARENA_ALLOC(arena, (n_pairs ? n_pairs : 1) * stride * sizeof(int32_t));
    int32_t *path_n = ARENA_ALLOC(arena, (n_pairs ? n_pairs : 1) * sizeof(int32_t));
    int n_pairs_int = (int)n_pairs, pi = 0;
#pragma omp parallel for schedule(dynamic,256) num_threads(threads) if(n_pairs_int >= 1024)
    for (pi = 0; pi < n_pairs_int; pi++) {
        int t = 0;
#ifdef _OPENMP
        t = omp_get_thread_num();
#endif
        path_n[pi] = ac_witness_path(run, deg, adj, pairs[pi].a, pairs[pi].b, o->max_path,
                                     scratch_prev + (size_t)t * nc, scratch_queue + (size_t)t * 2 * nc,
                                     paths + (size_t)pi * stride);
    }
    /* the diagnostic ledger: every contradicting pair with its witness path length and the seam /
     * layer evidence between the two charts (any relation, dropped or not; any layer pair) */
    FILE *diag = o->diag_pairs_csv ? fopen(o->diag_pairs_csv, "wb") : NULL;
    AcPairKey *rel_keys = NULL, *lay_keys = NULL; size_t n_rk = 0, n_lk = 0;
    if (diag) {
        fputs("chart_a,chart_b,component,mass,cells,dn_mean,dn_over_d,path_len,rel_flags,layer_d,d3d_centroid,area_a,area_b\n", diag);
        rel_keys = ARENA_ALLOC(arena, (run->n_rels ? run->n_rels : 1) * sizeof(AcPairKey));
        for (size_t r = 0; r < run->n_rels; r++) { rel_keys[n_rk].key = ac_pair_key(run->rels[r].a, run->rels[r].b); rel_keys[n_rk].flags = run->rels[r].flags; rel_keys[n_rk].value = run->rels[r].rms; n_rk++; }
        qsort(rel_keys, n_rk, sizeof(AcPairKey), ac_cmp_pair_key);
        lay_keys = ARENA_ALLOC(arena, (run->n_layers ? run->n_layers : 1) * sizeof(AcPairKey));
        for (size_t l = 0; l < run->n_layers; l++) { lay_keys[n_lk].key = ac_pair_key(run->layers[l].a, run->layers[l].b); lay_keys[n_lk].value = run->layers[l].d_median; lay_keys[n_lk].flags = 0; n_lk++; }
        qsort(lay_keys, n_lk, sizeof(AcPairKey), ac_cmp_pair_key);
    }
    for (size_t p = 0; p < n_pairs; p++) {
        if (out_chart_mass) { out_chart_mass[pairs[p].a] += pairs[p].mass; out_chart_mass[pairs[p].b] += pairs[p].mass; }
        int n = path_n[p]; const int32_t *path = paths + p * stride;
        for (int k = 0; k < n; k++) run->rels[path[k]].flag_mass += pairs[p].mass;
        if (diag) {
            const AsmChart *A = &run->charts[pairs[p].a], *B = &run->charts[pairs[p].b];
            double dl = run->chart_layer_d ? run->chart_layer_d[pairs[p].a] : 0.0;
            if (dl <= 0.0 && run->chart_layer_d) dl = run->chart_layer_d[pairs[p].b];
            if (dl <= 0.0) dl = run->layer_d_global > 0.0 ? run->layer_d_global : ASM_LAYER_PROBE_VOX * 0.4;
            double dn_mean = pairs[p].cells ? pairs[p].dn_sum / (double)pairs[p].cells : 0.0;
            const AcPairKey *rk = ac_find_pair_key(rel_keys, n_rk, ac_pair_key(pairs[p].a, pairs[p].b));
            const AcPairKey *lk = ac_find_pair_key(lay_keys, n_lk, ac_pair_key(pairs[p].a, pairs[p].b));
            double d3 = sqrt((A->centroid[0]-B->centroid[0])*(A->centroid[0]-B->centroid[0]) + (A->centroid[1]-B->centroid[1])*(A->centroid[1]-B->centroid[1]) + (A->centroid[2]-B->centroid[2])*(A->centroid[2]-B->centroid[2]));
            fprintf(diag, "%d,%d,%d,%.6e,%zu,%.3f,%.4f,%d,%ld,%.2f,%.1f,%.4e,%.4e\n", pairs[p].a, pairs[p].b, A->component, pairs[p].mass, pairs[p].cells,
                    dn_mean, dn_mean / dl, n, rk ? (long)rk->flags : -1L, lk ? lk->value : -1.0, d3, A->area3d, B->area3d);
        }
    }
    if (diag) fclose(diag);
    /* the report keeps its cell list in the run arena (below the mark would be freed):
     * copy it above the restore point by allocating before restoring */
    AsmConflictCell *keep = ARENA_ALLOC(arena, (n_cells ? n_cells : 1) * sizeof(AsmConflictCell));
    memcpy(keep, ccells, n_cells * sizeof(AsmConflictCell));
    Arena_restore(arena, mark);
    rep->contra_cells = ARENA_ALLOC(arena, (n_cells ? n_cells : 1) * sizeof(AsmConflictCell));
    memcpy(rep->contra_cells, keep, n_cells * sizeof(AsmConflictCell));
    rep->n_contra_cells = n_cells;
    return 0;
}

/* ---- layout-confirmed joins ------------------------------------------------------- */

typedef struct AcJoinEv { int32_t a, b; int8_t agree, has_off, a_is_A; double dist, dn, dot, gx, gy, bx, by, offx, offy; } AcJoinEv;

/* Jacobian of chart c's face f: d xyz / d (u,v) of the triangle's affine map (J[d*2+k]) */
static int ac_face_jacobian(const AsmChart *c, int32_t f, double J[6])
{
    const int32_t *fv = &c->faces[(size_t)f*3];
    double u0 = c->uv[(size_t)fv[0]*2], v0 = c->uv[(size_t)fv[0]*2+1];
    double e1u = c->uv[(size_t)fv[1]*2] - u0, e1v = c->uv[(size_t)fv[1]*2+1] - v0;
    double e2u = c->uv[(size_t)fv[2]*2] - u0, e2v = c->uv[(size_t)fv[2]*2+1] - v0;
    double det = e1u * e2v - e1v * e2u;
    if (fabs(det) < 1e-9) return -1;
    double i00 = e2v / det, i01 = -e2u / det, i10 = -e1v / det, i11 = e1u / det;   /* [e1 e2]^-1 */
    for (int d = 0; d < 3; d++) {
        double E1 = c->xyz[(size_t)fv[1]*3+(size_t)d] - c->xyz[(size_t)fv[0]*3+(size_t)d];
        double E2 = c->xyz[(size_t)fv[2]*3+(size_t)d] - c->xyz[(size_t)fv[0]*3+(size_t)d];
        J[d*2+0] = E1 * i00 + E2 * i10;
        J[d*2+1] = E1 * i01 + E2 * i11;
    }
    return 0;
}

/* least-squares (u,v) for J (u,v) = r */
static int ac_jacobian_solve(const double J[6], const double r[3], double d[2])
{
    double a = 0.0, b = 0.0, c2 = 0.0, g0 = 0.0, g1 = 0.0;
    for (int k = 0; k < 3; k++) { a += J[k*2]*J[k*2]; b += J[k*2]*J[k*2+1]; c2 += J[k*2+1]*J[k*2+1]; g0 += J[k*2]*r[k]; g1 += J[k*2+1]*r[k]; }
    double det = a * c2 - b * b;
    if (fabs(det) < 1e-12) return -1;
    d[0] = (c2 * g0 - b * g1) / det; d[1] = (a * g1 - b * g0) / det;
    return 0;
}

/* In-plane offset (global 2-D) of claim B's material relative to where chart A's own map
 * would draw it, given the global lattice step (gstep) between the two cells: +delta is the
 * shift that aligns B's drawing with A's continuation. */
static int ac_inplane_offset(const AsmRun *run, const AcClaim *A, const AcClaim *B, double gstep_x, double gstep_y, double out[2])
{
    const AsmChart *ca = &run->charts[(size_t)A->chart];
    if (ca->uv == NULL) return -1;
    double J[6];
    if (ac_face_jacobian(ca, A->face, J) != 0) return -1;
    double ct = cos(ca->pose_theta), st = sin(ca->pose_theta);
    double su = ct * gstep_x + st * gstep_y, sv = -st * gstep_x + ct * gstep_y;   /* R(-theta) */
    if (ca->flags & ASM_CHART_MIRROR) su = -su;
    double r[3], dc[2];
    for (int d = 0; d < 3; d++) r[d] = ((double)B->p[d] - A->p[d]) - (J[d*2] * su + J[d*2+1] * sv);
    if (ac_jacobian_solve(J, r, dc) != 0) return -1;
    if (ca->flags & ASM_CHART_MIRROR) dc[0] = -dc[0];
    out[0] = ct * dc[0] - st * dc[1]; out[1] = st * dc[0] + ct * dc[1];   /* R(theta) */
    return 0;
}

typedef struct AcJoinPair {
    int32_t a, b;
    size_t ev_first, ev_count;
    size_t agree, disagree, n_off;
    double dsum, dnsum, dotsum, sx, sy, sxx, syy, offx, offy;
    uint32_t vflags;
    int outcome, seam_off;
    double theta, tx, ty;
} AcJoinPair;

/* Where a seam relation (rigidly fitting, gate-rejected: WEAK) puts chart b relative to the
 * layout's chart a: the offset (global 2-D) that moves b onto the seam's transform. */
static int ac_seam_offset(const AsmRun *run, const AsmRelation *R, int32_t b_of_pair, double out[2])
{
    const AsmChart *A = &run->charts[(size_t)R->a], *B = &run->charts[(size_t)R->b];
    int ma = (A->flags & ASM_CHART_MIRROR) != 0;
    double theta = ma ? -R->theta : R->theta, tx = ma ? -R->tx : R->tx, ty = R->ty;   /* asm_pose's ap_effective */
    (void)theta;
    double c = cos(A->pose_theta), s = sin(A->pose_theta);
    double ex = A->pose_x + c * tx - s * ty, ey = A->pose_y + s * tx + c * ty;   /* where the seam puts B */
    double dx = ex - B->pose_x, dy = ey - B->pose_y;                            /* the shift B needs */
    if (b_of_pair == R->b) { out[0] = dx; out[1] = dy; }
    else { out[0] = -dx; out[1] = -dy; }
    return 0;
}

static int ac_cmp_join_ev(const void *x, const void *y)
{
    const AcJoinEv *p = x, *q = y;
    if (p->a != q->a) return p->a < q->a ? -1 : 1;
    return (p->b > q->b) - (p->b < q->b);
}

static double ac_wrap_pi(double a)
{
    while (a > 3.14159265358979323846) a -= 2.0 * 3.14159265358979323846;
    while (a < -3.14159265358979323846) a += 2.0 * 3.14159265358979323846;
    return a;
}

/* +1 the two claims agree (one sheet), -1 they contradict, 0 neither */
static int ac_join_evidence(const AsmRun *run, const AsmConflictOpts *o, const AcClaim *A, const AcClaim *B,
                            double expected, double *dist, double *dn, double *dot)
{
    double nx = A->n[0] + B->n[0], ny = A->n[1] + B->n[1], nz = A->n[2] + B->n[2];
    double nl = sqrt(nx*nx + ny*ny + nz*nz);
    double d[3] = { (double)A->p[0] - B->p[0], (double)A->p[1] - B->p[1], (double)A->p[2] - B->p[2] };
    *dist = sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
    *dot = (double)A->n[0]*B->n[0] + (double)A->n[1]*B->n[1] + (double)A->n[2]*B->n[2];
    *dn = nl > 1e-6 ? fabs((d[0]*nx + d[1]*ny + d[2]*nz) / nl) : *dist;
    double dl = run->chart_layer_d ? run->chart_layer_d[A->chart] : 0.0;
    if (dl <= 0.0 && run->chart_layer_d) dl = run->chart_layer_d[B->chart];
    if (dl <= 0.0) dl = run->layer_d_global > 0.0 ? run->layer_d_global : ASM_LAYER_PROBE_VOX * 0.4;
    if (fabs(*dot) >= ASM_JOIN_NORMAL_DOT && *dn <= o->contra_lo * dl && *dist <= expected + ASM_JOIN_AGREE_REACH * ASM_JOIN_STEP_SLACK) return 1;
    if (*dn >= o->contra_hi * dl || *dist > expected + 2.0 * ASM_JOIN_STEP_SLACK || fabs(*dot) < 0.5) return -1;
    return 0;
}

static size_t ac_lower_key(const AcClaim *c, size_t lo, size_t hi, int64_t key)
{
    while (lo < hi) { size_t mid = (lo + hi) / 2; if (c[mid].key < key) lo = mid + 1; else hi = mid; }
    return lo;
}

int AsmConflict_confirm_joins(AsmRun *run, const AsmConflictOpts *o, const char *diag_csv, AsmJoinStats *st)
{
    memset(st, 0, sizeof *st);
    double t0 = ves_clock_sec();
    /* Certified placement owns transactional correspondence fits. A raster
     * proposal cannot shift those poses after their audits have passed. */
    int align = ASM_JOIN_ALIGN && !run->continuity_ready;
    Arena_T arena = run->arena;
    Arena_Mark mark = Arena_save(arena);
    size_t nc = run->n_charts;
    /* lineages over the joins; pairs with a decision already taken against them */
    UnionFind uf = UF_new(arena, (int32_t)nc);
    size_t nveto = 0, nweak = 0;
    int64_t *veto = ARENA_ALLOC(arena, (run->n_rels ? run->n_rels : 1) * sizeof(int64_t));
    uint32_t *veto_flags = ARENA_ALLOC(arena, (run->n_rels ? run->n_rels : 1) * sizeof(uint32_t));
    int64_t *weak = ARENA_ALLOC(arena, (run->n_rels ? run->n_rels : 1) * sizeof(int64_t));
    for (size_t r = 0; r < run->n_rels; r++) {
        const AsmRelation *R = &run->rels[r];
        if (AsmRel_is_join(R)) uf_union(&uf, R->a, R->b);
        int32_t a = R->a < R->b ? R->a : R->b, b = R->a < R->b ? R->b : R->a;
        if (R->flags & ASM_JOIN_VETO_FLAGS) { veto[nveto] = ((int64_t)a << 32) | (int64_t)b; veto_flags[nveto] = R->flags; nveto++; }
        if ((R->flags & ASM_REL_PLACEMENT_ONLY) && !(R->flags & ASM_JOIN_VETO_FLAGS)) weak[nweak++] = ((int64_t)a << 32) | (int64_t)b;
    }
    for (size_t x = 1; x < nveto; x++) { int64_t v = veto[x]; uint32_t f = veto_flags[x]; size_t y = x; while (y > 0 && veto[y-1] > v) { veto[y] = veto[y-1]; veto_flags[y] = veto_flags[y-1]; y--; } veto[y] = v; veto_flags[y] = f; }
    for (size_t x = 1; x < nweak; x++) { int64_t v = weak[x]; size_t y = x; while (y > 0 && weak[y-1] > v) { weak[y] = weak[y-1]; y--; } weak[y] = v; }
    /* lineage area (the smaller side of a join is the one that moves) */
    double *lin_area = ARENA_CALLOC(arena, nc ? nc : 1, sizeof(double));
    {
        uint8_t *seen = ARENA_CALLOC(arena, nc ? nc : 1, 1);
        for (size_t i = 0; i < nc; i++) {
            const AsmChart *c = &run->charts[i];
            if (!AsmChart_in_layout(c) || !c->placed || c->component < 0) continue;
            int32_t rt = uf_find(&uf, (int32_t)i);
            lin_area[rt] += c->area3d;
            if (!seen[rt]) { seen[rt] = 1; st->lineages_before++; }
        }
    }
    /* raster, as the audit */
    AcBuf buf; memset(&buf, 0, sizeof buf); buf.arena = arena;
    if (ac_raster_all(run, o->cell, &buf)) { Arena_restore(arena, mark); return -1; }
    qsort(buf.c, buf.n, sizeof(AcClaim), ac_cmp_claim);
    /* evidence per contact of two charts of different lineages: co-located, right, down */
    size_t cap_ev = 1u << 12, n_ev = 0;
    AcJoinEv *ev = ARENA_ALLOC(arena, cap_ev * sizeof(AcJoinEv));
    for (size_t i = 0; i < buf.n;) {
        size_t j = i;
        while (j < buf.n && buf.c[j].key == buf.c[i].key) j++;
        /* sides: 0 co-located; 1..G right by g cells; G+1..2G down by g cells */
        const int G = ASM_JOIN_GAP_CELLS;
        for (int side = 0; side <= 2 * G; side++) {
            int g = side == 0 ? 0 : (side <= G ? side : side - G);
            int right = side >= 1 && side <= G;
            int64_t key = side == 0 ? buf.c[i].key : (right ? buf.c[i].key + g : buf.c[i].key + ((int64_t)g << 20));
            double expected = (double)g * o->cell;
            double gsx = right ? expected : 0.0, gsy = right ? 0.0 : expected;
            size_t lo, hi;
            if (side == 0) { lo = i; hi = j; }
            else { lo = ac_lower_key(buf.c, j, buf.n, key); hi = lo; while (hi < buf.n && buf.c[hi].key == key) hi++; }
            if (lo >= hi) continue;
            for (size_t k = i; k < j; k++) {
                if (g >= 2) {
                    /* across a gap only from A's LAST cell before it: the cell just before B holds no A */
                    int64_t kb = right ? buf.c[i].key + (g - 1) : buf.c[i].key + ((int64_t)(g - 1) << 20);
                    size_t b0 = ac_lower_key(buf.c, j, buf.n, kb), b1 = b0;
                    int has_a = 0;
                    while (b1 < buf.n && buf.c[b1].key == kb) { if (buf.c[b1].chart == buf.c[k].chart) has_a = 1; b1++; }
                    if (has_a) continue;
                }
                for (size_t l = (side == 0 ? k + 1 : lo); l < hi; l++) {
                    const AcClaim *A = &buf.c[k], *B = &buf.c[l];
                    if (A->chart == B->chart) continue;
                    if (uf_find(&uf, A->chart) == uf_find(&uf, B->chart)) continue;
                    double dist, dn, dot;
                    int e = ac_join_evidence(run, o, A, B, expected, &dist, &dn, &dot);
                    if (e == 0) continue;
                    if (n_ev == cap_ev) { AcJoinEv *ne = ARENA_ALLOC(arena, cap_ev * 2 * sizeof(AcJoinEv)); memcpy(ne, ev, n_ev * sizeof(AcJoinEv)); ev = ne; cap_ev *= 2; }
                    AcJoinEv *q = &ev[n_ev++];
                    q->a = A->chart < B->chart ? A->chart : B->chart;
                    q->b = A->chart < B->chart ? B->chart : A->chart;
                    q->agree = (int8_t)e; q->dist = dist; q->dn = dn; q->dot = dot;
                    q->gx = ((double)((int32_t)(A->key & 0xFFFFF) - (1 << 19)) + 0.5) * o->cell;
                    q->gy = ((double)((int32_t)((A->key >> 20) & 0xFFFFF) - (1 << 19)) + 0.5) * o->cell;
                    q->bx = ((double)((int32_t)(B->key & 0xFFFFF) - (1 << 19)) + 0.5) * o->cell;
                    q->by = ((double)((int32_t)((B->key >> 20) & 0xFFFFF) - (1 << 19)) + 0.5) * o->cell;
                    q->a_is_A = (int8_t)(A->chart < B->chart);   /* the pair's a is the claim A */
                    q->has_off = 0; q->offx = q->offy = 0.0;
                    if (e > 0) {
                        double off[2];
                        if (ac_inplane_offset(run, A, B, gsx, gsy, off) == 0) {
                            double sgn = A->chart < B->chart ? 1.0 : -1.0;   /* the shift of b (the larger id) */
                            q->offx = sgn * off[0]; q->offy = sgn * off[1]; q->has_off = 1;
                        }
                    }
                }
            }
        }
        i = j;
    }
    qsort(ev, n_ev, sizeof(AcJoinEv), ac_cmp_join_ev);
    /* pass 1: aggregate per pair and decide */
    size_t npair = 0;
    AcJoinPair *pairs = ARENA_ALLOC(arena, (n_ev ? n_ev : 1) * sizeof(AcJoinPair));
    for (size_t i = 0; i < n_ev;) {
        size_t j = i;
        AcJoinPair P; memset(&P, 0, sizeof P);
        P.a = ev[i].a; P.b = ev[i].b;
        P.ev_first = i;
        while (j < n_ev && ev[j].a == ev[i].a && ev[j].b == ev[i].b) {
            if (ev[j].agree > 0) {
                P.agree++; P.dsum += ev[j].dist; P.dnsum += ev[j].dn; P.dotsum += fabs(ev[j].dot);
                P.sx += ev[j].gx; P.sy += ev[j].gy; P.sxx += ev[j].gx * ev[j].gx; P.syy += ev[j].gy * ev[j].gy;
                if (ev[j].has_off) { P.offx += ev[j].offx; P.offy += ev[j].offy; P.n_off++; }
            } else P.disagree++;
            j++;
        }
        P.ev_count = j - i;
        st->pairs_tested++;
        st->cells_agree += P.agree; st->cells_disagree += P.disagree;
        int vetoed = 0;
        {
            int64_t key = ((int64_t)P.a << 32) | (int64_t)P.b;
            size_t lo = 0, hi = nveto;
            while (lo < hi) { size_t mid = (lo + hi) / 2; if (veto[mid] < key) lo = mid + 1; else hi = mid; }
            vetoed = lo < nveto && veto[lo] == key;
            if (vetoed) P.vflags = veto_flags[lo];
        }
        size_t min_cells = P.disagree == 0 ? (size_t)ASM_JOIN_MIN_CELLS_CLEAN : (size_t)ASM_JOIN_MIN_CELLS;
        if (P.disagree == 0) {
            int64_t key = ((int64_t)P.a << 32) | (int64_t)P.b;
            size_t lo = 0, hi = nweak;
            while (lo < hi) { size_t mid = (lo + hi) / 2; if (weak[mid] < key) lo = mid + 1; else hi = mid; }
            if (lo < nweak && weak[lo] == key) min_cells = (size_t)ASM_JOIN_MIN_CELLS_WEAK;
        }
        if (P.agree < min_cells) { st->pairs_refused_short++; P.outcome = 1; }
        else if ((double)P.disagree > ASM_JOIN_MAX_DISAGREE * (double)P.agree) { st->pairs_refused_disagree++; P.outcome = 2; }
        else if (vetoed) { st->pairs_refused_veto++; P.outcome = 3; }
        else if (uf_find(&uf, P.a) != uf_find(&uf, P.b)) P.outcome = 0;
        else P.outcome = 4;
        if (P.n_off) { P.offx /= (double)P.n_off; P.offy /= (double)P.n_off; }
        /* a rigidly fitting seam between the two charts knows the alignment better than cells do */
        if (align) {
            for (size_t r = 0; r < run->n_rels; r++) {
                const AsmRelation *R = &run->rels[r];
                if (!((R->a == P.a && R->b == P.b) || (R->a == P.b && R->b == P.a))) continue;
                if ((R->flags & ASM_JOIN_VETO_FLAGS) || R->corr_count <= 0) continue;
                double off[2];
                if (ac_seam_offset(run, R, P.b, off) == 0 && hypot(off[0], off[1]) <= ASM_JOIN_AGREE_REACH * ASM_JOIN_STEP_SLACK) {
                    P.offx = off[0]; P.offy = off[1]; if (P.n_off == 0) P.n_off = 1; P.seam_off = 1;
                }
                break;
            }
        }
        pairs[npair++] = P;
        i = j;
    }
    /* pass 2: the smaller lineage of every join moves by a RIGID 2-D fit (weighted Procrustes)
     * of its drawn contact points onto where the evidence puts them: every agreeing cell pair
     * gives (drawn position of the moving side's cell, that position + its offset), and a seam
     * relation gives one pair at the contact's centre with the seam's offset. */
    if (align) {
        /* per lineage root: weighted sums for the fit */
        double *sw = ARENA_CALLOC(arena, nc ? nc : 1, sizeof(double));
        double *ssx = ARENA_CALLOC(arena, nc ? nc : 1, sizeof(double)), *ssy = ARENA_CALLOC(arena, nc ? nc : 1, sizeof(double));
        double *stx = ARENA_CALLOC(arena, nc ? nc : 1, sizeof(double)), *sty = ARENA_CALLOC(arena, nc ? nc : 1, sizeof(double));
        double *sxx = ARENA_CALLOC(arena, nc ? nc : 1, sizeof(double)), *sxy = ARENA_CALLOC(arena, nc ? nc : 1, sizeof(double));
        double *syx = ARENA_CALLOC(arena, nc ? nc : 1, sizeof(double)), *syy = ARENA_CALLOC(arena, nc ? nc : 1, sizeof(double));
        double *sss = ARENA_CALLOC(arena, nc ? nc : 1, sizeof(double));   /* sum w |source|^2, for the spread */
        for (size_t p = 0; p < npair; p++) {
            const AcJoinPair *P = &pairs[p];
            if (P->outcome != 0 || P->n_off == 0) continue;
            int32_t ra = uf_find(&uf, P->a), rb = uf_find(&uf, P->b);
            int b_moves = lin_area[rb] <= lin_area[ra];
            int32_t rt = b_moves ? rb : ra;
            double sgn = b_moves ? 1.0 : -1.0;   /* the offset is the shift of b; a moves by its negative */
            if (P->seam_off) {
                /* one point pair at the contact's centre, weighted like the cells */
                double cx = P->agree ? P->sx / (double)P->agree : 0.0, cy = P->agree ? P->sy / (double)P->agree : 0.0;
                double w = (double)(P->agree ? P->agree : 1);
                double tx2 = cx + sgn * P->offx, ty2 = cy + sgn * P->offy;
                sw[rt] += w; ssx[rt] += w * cx; ssy[rt] += w * cy; stx[rt] += w * tx2; sty[rt] += w * ty2;
                sxx[rt] += w * cx * tx2; sxy[rt] += w * cx * ty2; syx[rt] += w * cy * tx2; syy[rt] += w * cy * ty2;
                sss[rt] += w * (cx * cx + cy * cy);
                continue;
            }
            for (size_t e = P->ev_first; e < P->ev_first + P->ev_count; e++) {
                const AcJoinEv *q = &ev[e];
                if (q->agree <= 0 || !q->has_off) continue;
                /* the moving side's own cell: b's cell is B's when a_is_A, else A's */
                int moving_is_B = b_moves ? q->a_is_A : !q->a_is_A;
                double px = moving_is_B ? q->bx : q->gx, py = moving_is_B ? q->by : q->gy;
                double tx2 = px + sgn * q->offx, ty2 = py + sgn * q->offy;
                sw[rt] += 1.0; ssx[rt] += px; ssy[rt] += py; stx[rt] += tx2; sty[rt] += ty2;
                sxx[rt] += px * tx2; sxy[rt] += px * ty2; syx[rt] += py * tx2; syy[rt] += py * ty2;
                sss[rt] += px * px + py * py;
            }
        }
        double *rot = ARENA_CALLOC(arena, nc ? nc : 1, sizeof(double)), *trx = ARENA_CALLOC(arena, nc ? nc : 1, sizeof(double)), *try_ = ARENA_CALLOC(arena, nc ? nc : 1, sizeof(double));
        uint8_t *has = ARENA_CALLOC(arena, nc ? nc : 1, 1);
        for (size_t r = 0; r < nc; r++) {
            if (sw[r] <= 0.0) continue;
            double w = sw[r], msx = ssx[r] / w, msy = ssy[r] / w, mtx = stx[r] / w, mty = sty[r] / w;
            /* centred cross-covariance: rotation from source to target */
            double cxx = sxx[r] - w * msx * mtx, cxy = sxy[r] - w * msx * mty, cyx = syx[r] - w * msy * mtx, cyy = syy[r] - w * msy * mty;
            /* rotation only from evidence spread wide enough to carry one */
            double spread = sqrt(fmax(sss[r] / w - (msx * msx + msy * msy), 0.0));
            double th = (w >= 3.0 && spread >= ASM_JOIN_ALIGN_ROT_SPREAD) ? atan2(cxy - cyx, cxx + cyy) : 0.0;
            if (th > ASM_JOIN_ALIGN_ROT_MAX) th = ASM_JOIN_ALIGN_ROT_MAX;
            if (th < -ASM_JOIN_ALIGN_ROT_MAX) th = -ASM_JOIN_ALIGN_ROT_MAX;
            /* the centroid moves by (mt - ms), capped; the rotation is about the source centroid */
            double dx = mtx - msx, dy = mty - msy;
            double m = hypot(dx, dy), reach = ASM_JOIN_AGREE_REACH * ASM_JOIN_STEP_SLACK;
            if (m > reach) { dx *= reach / m; dy *= reach / m; }
            double c = cos(th), s = sin(th);
            double tx2 = (msx + dx) - (c * msx - s * msy), ty2 = (msy + dy) - (s * msx + c * msy);   /* x' = R (x - ms) + ms + d */
            rot[r] = th; trx[r] = tx2; try_[r] = ty2; has[r] = 1;
            if (m > st->shift_max) st->shift_max = m;
            if (fabs(th) > st->rot_max) st->rot_max = fabs(th);
        }
        for (size_t i = 0; i < nc; i++) {
            AsmChart *c = &run->charts[i];
            if (!AsmChart_in_layout(c) || !c->placed || c->component < 0) continue;
            int32_t rt = uf_find(&uf, (int32_t)i);
            if (!has[rt]) continue;
            double cs = cos(rot[rt]), sn = sin(rot[rt]);
            double x = cs * c->pose_x - sn * c->pose_y + trx[rt], y = sn * c->pose_x + cs * c->pose_y + try_[rt];
            c->pose_x = x; c->pose_y = y; c->pose_theta = ac_wrap_pi(c->pose_theta + rot[rt]);
            st->charts_shifted++;
        }
    }
    /* Keep appended relations outside the scratch arena's lifetime. */
    AsmRelation *pending = malloc((npair ? npair : 1)*sizeof *pending);
    if (!pending) { Arena_restore(arena, mark); return -1; }
    size_t n_pending = 0;
    /* pass 3: relations from the (aligned) poses, and the ledger */
    FILE *csv = diag_csv ? fopen(diag_csv, "wb") : NULL;
    if (csv) fprintf(csv, "chart_a,chart_b,agree_cells,disagree_cells,dist_mean,dn_mean,dot_mean,off_x,off_y,seam_off,theta,tx,ty,parity,veto_flags,outcome\n");   /* outcome: 0 joined, 1 short, 2 contradicting, 3 vetoed, 4 already one lineage */
    for (size_t p = 0; p < npair; p++) {
        AcJoinPair *P = &pairs[p];
        int32_t a = P->a, b = P->b;
        const AsmChart *A = &run->charts[(size_t)a], *B = &run->charts[(size_t)b];
        int ma = (A->flags & ASM_CHART_MIRROR) != 0, mb = (B->flags & ASM_CHART_MIRROR) != 0;
        /* the placed relative pose in the mirrored frames: m_a = R(theta) m_b + t */
        P->theta = ac_wrap_pi(B->pose_theta - A->pose_theta);
        {
            double dx = B->pose_x - A->pose_x, dy = B->pose_y - A->pose_y;
            double ca = cos(-A->pose_theta), sa = sin(-A->pose_theta);
            P->tx = ca * dx - sa * dy; P->ty = sa * dx + ca * dy;
        }
        if (P->outcome == 0) {
            AsmRelation r;
            memset(&r, 0, sizeof r);
            r.a = a; r.b = b;
            /* stored in the unmirrored frames (asm_pose's ap_effective inverts this) */
            r.theta = ma ? -P->theta : P->theta; r.tx = ma ? -P->tx : P->tx; r.ty = P->ty;
            r.flags = ASM_REL_LAYOUT | (ma != mb ? ASM_REL_PARITY : 0u);
            double n_equiv = 4.0 * (double)P->agree;   /* one 8-vox cell ~ four 2-vox correspondences */
            double rms = P->dnsum / (double)P->agree;
            double var = rms * rms; if (var < 0.25) var = 0.25;
            double s2 = (P->sxx - P->sx * P->sx / (double)P->agree) + (P->syy - P->sy * P->sy / (double)P->agree);
            r.w_xy = n_equiv / var;
            r.w_theta = 4.0 * s2 / var;
            r.rms = rms; r.seam_len = (double)P->agree * o->cell;
            r.normal_gap = P->dnsum / (double)P->agree; r.mean_gap = P->dsum / (double)P->agree;
            r.n_corr = (size_t)n_equiv; r.corr_first = -1; r.corr_count = 0;
            r.robust_w = 1.0;
            pending[n_pending++] = r;
            uf_union(&uf, a, b);
            st->pairs_joined++;
        }
        if (csv) fprintf(csv, "%d,%d,%zu,%zu,%.3f,%.3f,%.3f,%.3f,%.3f,%d,%.6f,%.3f,%.3f,%d,%u,%d\n", a, b, P->agree, P->disagree,
                         P->agree ? P->dsum / (double)P->agree : 0.0, P->agree ? P->dnsum / (double)P->agree : 0.0, P->agree ? P->dotsum / (double)P->agree : 0.0,
                         P->offx, P->offy, P->seam_off, P->theta, P->tx, P->ty, ma != mb, P->vflags, P->outcome);
    }
    if (csv) fclose(csv);
    {
        uint8_t *seen = ARENA_CALLOC(arena, nc ? nc : 1, 1);
        for (size_t i = 0; i < nc; i++) {
            const AsmChart *c = &run->charts[i];
            if (!AsmChart_in_layout(c) || !c->placed || c->component < 0) continue;
            int32_t rt = uf_find(&uf, (int32_t)i);
            if (!seen[rt]) { seen[rt] = 1; st->lineages_after++; }
        }
    }
    Arena_restore(arena, mark);
    for (size_t i = 0; i < n_pending; i++) AsmRun_push_rel(run, &pending[i]);
    free(pending);
    st->sec = ves_clock_sec() - t0;
    return 0;
}

/* ---- incremental grid ------------------------------------------------------------ */

struct AsmConflictGrid {
    Arena_T arena;          /* permanent: the table, the claim store (never restored) */
    Arena_T scratch;        /* rasters and per-probe scratch: marked and restored per call */
    const AsmRun *run;
    AsmConflictOpts o;
    /* open-addressing table: cell key -> head claim index (-1 = empty) */
    int64_t *keys; int32_t *head; size_t cap, count;
    /* claims with next links */
    AcClaim *claims; int32_t *next; size_t n, ncap;
    double *mass_scratch;   /* [n_charts] partner accumulation */
};

static int64_t acg_cell_key(int32_t cx, int32_t cy) { return ((int64_t)(cy + (1 << 19)) << 20) | (int64_t)(cx + (1 << 19)); }

static size_t acg_hash(int64_t key, size_t cap) { return (size_t)(((uint64_t)key * 0x9E3779B97F4A7C15ULL) >> 20) & (cap - 1); }

static void acg_rehash(AsmConflictGrid *g, size_t ncap)
{
    int64_t *nk = ARENA_ALLOC(g->arena, ncap * sizeof(int64_t));
    int32_t *nh = ARENA_ALLOC(g->arena, ncap * sizeof(int32_t));
    for (size_t i = 0; i < ncap; i++) nh[i] = -1;
    for (size_t i = 0; i < g->cap; i++) {
        if (g->head[i] < 0) continue;
        size_t h = acg_hash(g->keys[i], ncap);
        while (nh[h] >= 0) h = (h + 1) & (ncap - 1);
        nk[h] = g->keys[i]; nh[h] = g->head[i];
    }
    g->keys = nk; g->head = nh; g->cap = ncap;
}

static int32_t *acg_slot(AsmConflictGrid *g, int64_t key, int create)
{
    if (create && (g->count + 1) * 2 > g->cap) {
        acg_rehash(g, g->cap * 2);
    }
    size_t h = acg_hash(key, g->cap), steps = 0;
    for (;;) {
        if (g->head[h] < 0) {
            if (!create) return NULL;
            g->keys[h] = key; g->count++;
            return &g->head[h];
        }
        if (g->keys[h] == key) return &g->head[h];
        h = (h + 1) & (g->cap - 1);
        if (++steps > g->cap) { fprintf(stderr, "    [conflict grid] slot search exhausted (cap %zu count %zu)\n", g->cap, g->count); return NULL; }
    }
}

AsmConflictGrid *AsmConflictGrid_new(Arena_T arena, const AsmConflictOpts *o, const AsmRun *run)
{
    AsmConflictGrid *g = ARENA_CALLOC(arena, 1, sizeof *g);
    g->arena = arena; g->scratch = Arena_new(); g->run = run; g->o = *o;
    g->cap = 1u << 16;
    g->keys = ARENA_ALLOC(arena, g->cap * sizeof(int64_t));
    g->head = ARENA_ALLOC(arena, g->cap * sizeof(int32_t));
    for (size_t i = 0; i < g->cap; i++) g->head[i] = -1;
    g->ncap = 1u << 16;
    g->claims = ARENA_ALLOC(arena, g->ncap * sizeof(AcClaim));
    g->next = ARENA_ALLOC(arena, g->ncap * sizeof(int32_t));
    g->mass_scratch = ARENA_CALLOC(arena, run->n_charts ? run->n_charts : 1, sizeof(double));
    return g;
}

size_t AsmConflictGrid_cells(const AsmConflictGrid *g) { return g->count; }

void AsmConflictGrid_insert(AsmConflictGrid *g, const AsmChart *c)
{
    Arena_Mark mark = Arena_save(g->scratch);
    AcBuf buf; memset(&buf, 0, sizeof buf); buf.arena = g->scratch;
    ac_raster_chart(c, g->o.cell, &buf, ARENA_ALLOC(g->scratch, (c->nv ? 2*c->nv : 1) * sizeof(double)));
    size_t n = buf.n;
    if (g->n + n > g->ncap) {
        size_t ncap = g->ncap;
        while (ncap < g->n + n) ncap *= 2;
        AcClaim *nc = ARENA_ALLOC(g->arena, ncap * sizeof(AcClaim));
        int32_t *nn = ARENA_ALLOC(g->arena, ncap * sizeof(int32_t));
        memcpy(nc, g->claims, g->n * sizeof(AcClaim)); memcpy(nn, g->next, g->n * sizeof(int32_t));
        g->claims = nc; g->next = nn; g->ncap = ncap;
    }
    for (size_t i = 0; i < n; i++) {
        int64_t key = buf.c[i].key & ((((int64_t)1) << 40) - 1);   /* drop the component: one placed frame */
        int32_t *slot = acg_slot(g, key, 1);
        if (slot == NULL) break;
        g->claims[g->n] = buf.c[i]; g->claims[g->n].key = key;
        g->next[g->n] = *slot; *slot = (int32_t)g->n; g->n++;
    }
    Arena_restore(g->scratch, mark);
}

void AsmConflictGrid_dispose(AsmConflictGrid *g)
{
    if (g != NULL && g->scratch != NULL) Arena_dispose(&g->scratch);
}

static double acg_weight(const AsmConflictGrid *g, const AcClaim *A, const AcClaim *B)
{
    const AsmRun *run = g->run;
    double nx = A->n[0] + B->n[0], ny = A->n[1] + B->n[1], nz = A->n[2] + B->n[2];
    double nl = sqrt(nx*nx + ny*ny + nz*nz), dn;
    if (nl > 1e-6) dn = fabs(((A->p[0]-B->p[0])*nx + (A->p[1]-B->p[1])*ny + (A->p[2]-B->p[2])*nz) / nl);
    else dn = sqrt((A->p[0]-B->p[0])*(A->p[0]-B->p[0]) + (A->p[1]-B->p[1])*(A->p[1]-B->p[1]) + (A->p[2]-B->p[2])*(A->p[2]-B->p[2]));
    double dl = run->chart_layer_d ? run->chart_layer_d[A->chart] : 0.0;
    if (dl <= 0.0 && run->chart_layer_d) dl = run->chart_layer_d[B->chart];
    if (dl <= 0.0) dl = run->layer_d_global > 0.0 ? run->layer_d_global : ASM_LAYER_PROBE_VOX * 0.4;
    return ac_smoothstep((dn / dl - g->o.contra_lo) / (g->o.contra_hi - g->o.contra_lo));
}

double AsmConflictGrid_probe_excluding(AsmConflictGrid *g, const AsmChart *const *charts, size_t n,
                                       const uint8_t *excluded, int32_t *partner, double *partner_mass)
{
    Arena_Mark mark = Arena_save(g->scratch);
    size_t dbg_claims = 0, dbg_walk = 0;
    double cell_area = g->o.cell * g->o.cell, mass = 0.0;
    size_t nch = g->run->n_charts;
    /* partner accumulation: touched list to zero afterwards */
    int32_t *touched = ARENA_ALLOC(g->scratch, (nch ? nch : 1) * sizeof(int32_t));
    size_t ntouched = 0;
    for (size_t k = 0; k < n; k++) {
        Arena_Mark m2 = Arena_save(g->scratch);
        AcBuf buf; memset(&buf, 0, sizeof buf); buf.arena = g->scratch;
        ac_raster_chart(charts[k], g->o.cell, &buf, ARENA_ALLOC(g->scratch, (charts[k]->nv ? 2*charts[k]->nv : 1) * sizeof(double)));
        for (size_t i = 0; i < buf.n; i++) {
            int64_t key = buf.c[i].key & ((((int64_t)1) << 40) - 1);
            int32_t *slot = acg_slot(g, key, 0);
            dbg_claims++;
            if (slot == NULL) continue;
            double wmax = 0.0; int32_t wchart = -1;
            for (int32_t q = *slot; q >= 0; q = g->next[q]) {
                if (++dbg_walk > 50000000) { fprintf(stderr, "[conflict grid] runaway list walk at claim %zu (cycle?)\n", dbg_claims); break; }
                const AcClaim *B = &g->claims[q];
                if (excluded && excluded[B->chart]) continue;
                if (B->chart == buf.c[i].chart) continue;
                double w = acg_weight(g, &buf.c[i], B);
                if (w > wmax) { wmax = w; wchart = B->chart; }
            }
            if (wmax > 0.0) {
                mass += wmax * cell_area;
                if (wchart >= 0 && (size_t)wchart < nch) {
                    if (g->mass_scratch[wchart] == 0.0) touched[ntouched++] = wchart;
                    g->mass_scratch[wchart] += wmax * cell_area;
                }
            }
        }
        Arena_restore(g->scratch, m2);
    }
    int32_t best = -1; double bm = 0.0;
    for (size_t t = 0; t < ntouched; t++) {
        if (g->mass_scratch[touched[t]] > bm) { bm = g->mass_scratch[touched[t]]; best = touched[t]; }
        g->mass_scratch[touched[t]] = 0.0;
    }
    if (partner) *partner = best;
    if (partner_mass) *partner_mass = bm;
    Arena_restore(g->scratch, mark);
    return mass;
}

double AsmConflictGrid_probe(AsmConflictGrid *g, const AsmChart *const *charts, size_t n, int32_t *partner, double *partner_mass)
{
    return AsmConflictGrid_probe_excluding(g,charts,n,NULL,partner,partner_mass);
}

size_t AsmConflictGrid_stacked(const AsmConflictGrid *g, size_t *cells)
{
    size_t stacked = 0;
    for (size_t h = 0; h < g->cap; h++) {
        if (g->head[h] < 0) continue;
        int hit = 0;
        for (int32_t q = g->head[h]; q >= 0 && !hit; q = g->next[q])
            for (int32_t r = g->next[q]; r >= 0; r = g->next[r]) {
                if (g->claims[r].chart == g->claims[q].chart) continue;
                if (acg_weight(g, &g->claims[q], &g->claims[r]) > 0.5) { hit = 1; break; }
            }
        stacked += (size_t)hit;
    }
    if (cells) *cells = g->count;
    return stacked;
}

/* ---- cleaning rounds ----------------------------------------------------------- */

typedef struct AcFlag { double mass; int32_t rel; } AcFlag;

static int ac_cmp_flag_desc(const void *x, const void *y)
{
    const AcFlag *a = x, *b = y;
    if (a->mass != b->mass) return a->mass > b->mass ? -1 : 1;
    return (a->rel > b->rel) - (a->rel < b->rel);
}

/* A layout snapshot: everything the pose solve and the extras pass write.  The robust pose
 * solve is re-initialized from the chordal solve on every call and is MULTI-MODAL on a graph
 * with wrong joins: removing one 17k-vox^2 chart or three relations moved the familiar's
 * score from 7.9e6 to 7.0e6 and back (2026-09-08).  So a rejected step is undone by restoring
 * the snapshot, never by re-solving, and the loop finishes on the best audited state. */
typedef struct AcSnap { double *pose; int32_t *pc; uint32_t *cflags; uint32_t *rflags; double *rw; double score; } AcSnap;

static void ac_snap_alloc(Arena_T arena, AsmRun *run, AcSnap *s)
{
    size_t nch = run->n_charts ? run->n_charts : 1, nr = run->n_rels ? run->n_rels : 1;
    s->pose = ARENA_ALLOC(arena, nch * 3 * sizeof(double));
    s->pc = ARENA_ALLOC(arena, nch * 2 * sizeof(int32_t));
    s->cflags = ARENA_ALLOC(arena, nch * sizeof(uint32_t));
    s->rflags = ARENA_ALLOC(arena, nr * sizeof(uint32_t));
    s->rw = ARENA_ALLOC(arena, nr * 2 * sizeof(double));
    s->score = -1e300;
}
static void ac_snap_take(AsmRun *run, AcSnap *s, double score)
{
    for (size_t i = 0; i < run->n_charts; i++) {
        const AsmChart *c = &run->charts[i];
        s->pose[i*3] = c->pose_x; s->pose[i*3+1] = c->pose_y; s->pose[i*3+2] = c->pose_theta;
        s->pc[i*2] = c->placed; s->pc[i*2+1] = c->component; s->cflags[i] = c->flags;
    }
    for (size_t r = 0; r < run->n_rels; r++) { s->rflags[r] = run->rels[r].flags; s->rw[r*2] = run->rels[r].robust_w; s->rw[r*2+1] = run->rels[r].residual; }
    s->score = score;
}
static void ac_snap_restore(AsmRun *run, const AcSnap *s)
{
    for (size_t i = 0; i < run->n_charts; i++) {
        AsmChart *c = &run->charts[i];
        c->pose_x = s->pose[i*3]; c->pose_y = s->pose[i*3+1]; c->pose_theta = s->pose[i*3+2];
        c->placed = s->pc[i*2]; c->component = s->pc[i*2+1]; c->flags = s->cflags[i];
    }
    for (size_t r = 0; r < run->n_rels; r++) { run->rels[r].flags = s->rflags[r]; run->rels[r].robust_w = s->rw[r*2]; run->rels[r].residual = s->rw[r*2+1]; }
}

/* ---- post-placement seam re-solve ------------------------------------------------- */

typedef struct AcRefCand { int32_t rel; int kind; double gap, rot; int readmit; int reason; } AcRefCand;   /* reason: 0 none, 1 gap, 2 rotation, 3 parity, 4 no correspondences */

static double ac_wrap_angle(double a)
{
    while (a > 3.14159265358979323846) a -= 2.0 * 3.14159265358979323846;
    while (a < -3.14159265358979323846) a += 2.0 * 3.14159265358979323846;
    return a;
}

/* the median layout distance between a relation's correspondences: the slit's width */
static double ac_corr_gap(const AsmRun *run, const AsmRelation *R, Arena_T arena)
{
    if (R->corr_count <= 0 || R->corr_first < 0) return -1.0;
    Arena_Mark mark = Arena_save(arena);
    double *d = ARENA_ALLOC(arena, (size_t)R->corr_count * sizeof(double));
    const AsmChart *A = &run->charts[(size_t)R->a], *B = &run->charts[(size_t)R->b];
    size_t n = 0;
    for (int32_t k = 0; k < R->corr_count; k++) {
        const AsmCorr *cr = &run->corr[(size_t)R->corr_first + (size_t)k];
        if (cr->va < 0 || (size_t)cr->va >= A->nv || cr->vb < 0 || (size_t)cr->vb >= B->nv) continue;
        double ax, ay, bx, by;
        ac_pose(A, A->uv[(size_t)cr->va*2], A->uv[(size_t)cr->va*2+1], &ax, &ay);
        ac_pose(B, B->uv[(size_t)cr->vb*2], B->uv[(size_t)cr->vb*2+1], &bx, &by);
        d[n++] = hypot(ax - bx, ay - by);
    }
    double med = -1.0;
    if (n > 0) { qsort(d, n, sizeof(double), av_cmp_double_ac); med = d[n / 2]; }
    Arena_restore(arena, mark);
    return med;
}

static size_t ac_count_lineages(AsmRun *run, const uint8_t *placed, Arena_T arena, int32_t *root_out, double *lin_area_out)
{
    size_t nc = run->n_charts;
    UnionFind uf = UF_new(arena, (int32_t)nc);
    for (size_t r = 0; r < run->n_rels; r++) {
        const AsmRelation *R = &run->rels[r];
        if (!AsmRel_is_join(R) || !placed[R->a] || !placed[R->b]) continue;
        uf_union(&uf, R->a, R->b);
    }
    size_t n = 0;
    for (size_t i = 0; i < nc; i++) {
        if (!placed[i]) { if (root_out) root_out[i] = -1; continue; }
        int32_t rt = uf_find(&uf, (int32_t)i);
        if (root_out) root_out[i] = rt;
        if (lin_area_out) { if (lin_area_out[rt] == 0.0) n++; lin_area_out[rt] += run->charts[i].area3d; }
    }
    if (!lin_area_out) {
        uint8_t *seen = ARENA_CALLOC(arena, nc ? nc : 1, 1);
        for (size_t i = 0; i < nc; i++) { if (!placed[i]) continue; int32_t rt = uf_find(&uf, (int32_t)i); if (!seen[rt]) { seen[rt] = 1; n++; } }
    }
    return n;
}

int AsmConflict_refine_seams(AsmRun *run, const AsmConflictOpts *o, const AsmPoseOpts *po, const char *diag_csv, AsmRefineStats *st)
{
    memset(st, 0, sizeof *st);
    double t0 = ves_clock_sec();
    Arena_T arena = run->arena;
    Arena_Mark mark = Arena_save(arena);
    size_t nc = run->n_charts, nr = run->n_rels;
    /* the placed frame: the component holding the largest placed area */
    int32_t comp_primary = -1;
    {
        int32_t cmax = -1;
        for (size_t i = 0; i < nc; i++) if (run->charts[i].component > cmax) cmax = run->charts[i].component;
        if (cmax < 0) { Arena_restore(arena, mark); return 0; }
        double *carea = ARENA_CALLOC(arena, (size_t)cmax + 1, sizeof(double));
        for (size_t i = 0; i < nc; i++) { const AsmChart *c = &run->charts[i]; if (AsmChart_in_layout(c) && c->placed && c->component >= 0) carea[c->component] += c->area3d; }
        for (int32_t c = 0; c <= cmax; c++) if (comp_primary < 0 || carea[c] > carea[comp_primary]) comp_primary = c;
    }
    uint8_t *placed = ARENA_CALLOC(arena, nc ? nc : 1, 1);
    for (size_t i = 0; i < nc; i++) { const AsmChart *c = &run->charts[i]; placed[i] = AsmChart_in_layout(c) && c->placed && c->component == comp_primary; if (placed[i]) st->charts_placed++; }
    /* lineages over the joins; the largest by area is the gauge */
    int32_t *root = ARENA_ALLOC(arena, (nc ? nc : 1) * sizeof(int32_t));
    double *lin_area = ARENA_CALLOC(arena, nc ? nc : 1, sizeof(double));
    st->lineages_before = ac_count_lineages(run, placed, arena, root, lin_area);
    int32_t fixed_root = -1;
    for (size_t i = 0; i < nc; i++) if (placed[i] && (fixed_root < 0 || lin_area[root[i]] > lin_area[fixed_root])) fixed_root = root[i];
    uint8_t *fixed = ARENA_CALLOC(arena, nc ? nc : 1, 1);
    for (size_t i = 0; i < nc; i++) if (placed[i] && root[i] == fixed_root) { fixed[i] = 1; st->charts_fixed++; }
    /* candidates */
    AcRefCand *cand = ARENA_ALLOC(arena, (nr ? nr : 1) * sizeof(AcRefCand));
    int32_t *rel_idx = ARENA_ALLOC(arena, (nr ? nr : 1) * sizeof(int32_t));
    size_t ncand = 0, nsel = 0;
    for (size_t r = 0; r < nr; r++) {
        const AsmRelation *R = &run->rels[r];
        if (R->flags & ASM_REL_CONTACT) continue;
        if (!placed[R->a] || !placed[R->b]) continue;
        if (AsmRel_is_join(R)) {
            if (fixed[R->a] && fixed[R->b]) { st->internal_fixed++; continue; }
            st->cand_join++; rel_idx[nsel++] = (int32_t)r;
            continue;
        }
        int kind = (R->flags & ASM_REL_SHORT) ? 4 : (R->flags & ASM_REL_WEAK) ? 3 : (R->flags & ASM_REL_SWITCHED) ? 2 : (R->flags & ASM_REL_DROPPED) ? 1 : 0;
        if (kind == 0) continue;
        if (kind == 1 && !ASM_REFINE_READMIT_DROPPED) continue;
        st->cand_by_kind[kind]++;
        AcRefCand *cd = &cand[ncand++];
        cd->rel = (int32_t)r; cd->kind = kind; cd->gap = -1.0; cd->rot = 0.0; cd->readmit = 0; cd->reason = 0;
        const AsmChart *A = &run->charts[(size_t)R->a], *B = &run->charts[(size_t)R->b];
        int ma = (A->flags & ASM_CHART_MIRROR) != 0, mb = (B->flags & ASM_CHART_MIRROR) != 0;
        if ((ma ^ mb) != ((R->flags & ASM_REL_PARITY) != 0)) { cd->reason = 3; st->refused_parity++; continue; }
        cd->gap = ac_corr_gap(run, R, arena);
        if (cd->gap < 0.0) { cd->reason = 4; st->refused_nocorr++; continue; }
        double eff_theta = ma ? -R->theta : R->theta;
        cd->rot = fabs(ac_wrap_angle(B->pose_theta - A->pose_theta - eff_theta));
        if (cd->gap > ASM_READMIT_VOX) { cd->reason = 1; st->refused_gap++; continue; }
        if (cd->rot > ASM_READMIT_RAD) { cd->reason = 2; st->refused_rot++; continue; }
        cd->readmit = 1; st->readmitted++; st->readmit_by_kind[kind]++;
        rel_idx[nsel++] = (int32_t)r;
    }
    AcSnap snap; ac_snap_alloc(arena, run, &snap); ac_snap_take(run, &snap, 0.0);
    double *mass0 = ARENA_CALLOC(arena, nc ? nc : 1, sizeof(double)), *mass1 = ARENA_CALLOC(arena, nc ? nc : 1, sizeof(double));
    AsmConflictReport rep0, rep1;
    memset(&rep1, 0, sizeof rep1);
    AsmConflict_audit(run, o, mass0, &rep0);
    st->contra_before = rep0.contra_area; st->covered_before = rep0.covered_area;
    if (st->readmitted == 0) {
        st->contra_after = rep0.contra_area; st->covered_after = rep0.covered_area; st->lineages_after = st->lineages_before;
        st->sec = ves_clock_sec() - t0;
        Arena_restore(arena, mark);
        return 0;
    }
    for (size_t k = 0; k < ncand; k++) if (cand[k].readmit) { AsmRelation *R = &run->rels[(size_t)cand[k].rel]; R->flags &= ~(uint32_t)(ASM_REL_DROPPED | ASM_REL_SWITCHED | ASM_REL_WEAK | ASM_REL_SHORT); R->flags |= ASM_REL_READMIT; }
    double limit = rep0.contra_area * (1.0 + ASM_REFINE_MAX_CONTRA_GROWTH) + ASM_REFINE_CONTRA_FLOOR;
    for (int round = 1; round <= ASM_REFINE_ROUNDS; round++) {
        st->rounds = round;
        AsmPoseRefineStats ps;
        AsmPose_refine(run, po, rel_idx, nsel, fixed, &ps);
        st->iters = ps.iters; st->switched = ps.switched; st->charts_free = ps.charts_free; st->charts_pinned = ps.charts_pinned;
        /* a readmitted seam the robust solve switched off goes back to what it was: not a join */
        for (size_t k = 0; k < ncand; k++) {
            if (!cand[k].readmit) continue;
            AsmRelation *R = &run->rels[(size_t)cand[k].rel];
            if (R->flags & ASM_REL_SWITCHED) { R->flags = snap.rflags[(size_t)cand[k].rel]; cand[k].readmit = 0; st->readmit_switched++; }
        }
        AsmConflict_audit(run, o, mass1, &rep1);
        if (rep1.contra_area <= limit) { st->kept = 1; break; }
        if (round == ASM_REFINE_ROUNDS) break;
        /* round 2: back to the snapshot's poses; the seams touching a chart whose contradiction grew are un-readmitted */
        for (size_t i = 0; i < nc; i++) { AsmChart *c = &run->charts[i]; c->pose_x = snap.pose[i*3]; c->pose_y = snap.pose[i*3+1]; c->pose_theta = snap.pose[i*3+2]; }
        size_t removed = 0;
        for (size_t k = 0; k < ncand; k++) {
            if (!cand[k].readmit) continue;
            AsmRelation *R = &run->rels[(size_t)cand[k].rel];
            double ga = mass1[R->a] - mass0[R->a], gb = mass1[R->b] - mass0[R->b];
            if (ga > 0.25 * run->charts[(size_t)R->a].area3d || gb > 0.25 * run->charts[(size_t)R->b].area3d) {
                R->flags = snap.rflags[(size_t)cand[k].rel]; cand[k].readmit = 0; removed++; st->unreadmitted++;
            }
        }
        if (removed == 0) break;
        nsel = 0;
        for (size_t r = 0; r < nr; r++) {
            const AsmRelation *R = &run->rels[r];
            if ((R->flags & ASM_REL_CONTACT) || !placed[R->a] || !placed[R->b] || !AsmRel_is_join(R)) continue;
            if (fixed[R->a] && fixed[R->b]) continue;
            rel_idx[nsel++] = (int32_t)r;
        }
    }
    if (!st->kept) {
        ac_snap_restore(run, &snap);
        st->reverted = 1;
        for (size_t k = 0; k < ncand; k++) cand[k].readmit = 0;
        AsmConflict_audit(run, o, NULL, &rep1);
    }
    st->contra_after = rep1.contra_area; st->covered_after = rep1.covered_area;
    /* what moved */
    {
        double *sh = ARENA_ALLOC(arena, (nc ? nc : 1) * sizeof(double));
        size_t n = 0;
        for (size_t i = 0; i < nc; i++) {
            if (!placed[i] || fixed[i]) continue;
            const AsmChart *c = &run->charts[i];
            double d = hypot(c->pose_x - snap.pose[i*3], c->pose_y - snap.pose[i*3+1]);
            double dr = fabs(ac_wrap_angle(c->pose_theta - snap.pose[i*3+2]));
            if (d > 2.0) { st->moved++; st->moved_area += c->area3d; }
            if (d > st->shift_max) st->shift_max = d;
            if (dr > st->rot_max) st->rot_max = dr;
            sh[n++] = d;
        }
        if (n > 0) { qsort(sh, n, sizeof(double), av_cmp_double_ac); st->shift_p50 = sh[n/2]; st->shift_p90 = sh[(size_t)((double)(n-1)*0.9)]; }
    }
    st->lineages_after = ac_count_lineages(run, placed, arena, NULL, NULL);
    if (diag_csv && diag_csv[0]) {
        FILE *fp = fopen(diag_csv, "wb");
        if (fp) {
            fputs("chart_a,chart_b,kind,n_corr,rms,gap,rot,readmitted,reason,robust_w,switched,flags\n", fp);
            for (size_t k = 0; k < ncand; k++) {
                const AsmRelation *R = &run->rels[(size_t)cand[k].rel];
                fprintf(fp, "%d,%d,%d,%zu,%.3f,%.2f,%.4f,%d,%d,%.3f,%d,%u\n", R->a, R->b, cand[k].kind, R->n_corr, R->rms, cand[k].gap, cand[k].rot,
                        cand[k].readmit, cand[k].reason, R->robust_w, (R->flags & ASM_REL_SWITCHED) != 0, R->flags);
            }
            fclose(fp);
        }
    }
    st->sec = ves_clock_sec() - t0;
    Arena_restore(arena, mark);
    return 0;
}

/* The area a layout separated from the snapshot's components: for every component of the
 * snapshot, the area of its charts now outside the largest component they landed in (a dropped
 * bridge, or a switch the drop triggered; a chart that left the layout counts as separated). */
typedef struct AcPart { int32_t old_c, new_c; double area; } AcPart;
static int ac_cmp_part(const void *x, const void *y)
{
    const AcPart *a = x, *b = y;
    if (a->old_c != b->old_c) return (a->old_c > b->old_c) - (a->old_c < b->old_c);
    return (a->new_c > b->new_c) - (a->new_c < b->new_c);
}
static double ac_split_area(AsmRun *run, const AcSnap *s)
{
    size_t n = run->n_charts;
    Arena_Mark mark = Arena_save(run->arena);
    AcPart *p = ARENA_ALLOC(run->arena, (n ? n : 1) * sizeof(AcPart));
    size_t m = 0;
    for (size_t i = 0; i < n; i++) {
        const AsmChart *c = &run->charts[i];
        int32_t oc = s->pc[i*2+1];
        if (!s->pc[i*2] || oc < 0 || (s->cflags[i] & ASM_CHART_EXCLUDED_MASK) != 0u) continue;   /* not in the layout before the batch */
        int32_t nc = (AsmChart_in_layout(c) && c->placed && c->component >= 0) ? c->component : -1;
        p[m].old_c = oc; p[m].new_c = nc; p[m].area = c->area3d; m++;
    }
    if (m > 1) qsort(p, m, sizeof(AcPart), ac_cmp_part);
    double separated = 0.0;
    size_t i = 0;
    while (i < m) {
        size_t j = i;
        double total = 0.0, largest = 0.0;
        while (j < m && p[j].old_c == p[i].old_c) {
            size_t k = j;
            double part = 0.0;
            while (k < m && p[k].old_c == p[i].old_c && p[k].new_c == p[j].new_c) { part += p[k].area; k++; }
            total += part;
            if (p[j].new_c >= 0 && part > largest) largest = part;
            j = k;
        }
        separated += total - largest;
        i = j;
    }
    Arena_restore(run->arena, mark);
    return separated;
}

/* clean2's round budget: the configured rounds, then more while the last
 * round was kept with a gain of at least ASM_CLEAN_EXTEND_GAIN of the score,
 * within ASM_CLEAN_MAX_ROUNDS and ASM_CLEAN_EXTEND_SECONDS. */
static int ac_clean_extend(double gain, double score)
{
    return ASM_CLEAN_EXTEND_GAIN > 0.0 && isfinite(gain) && gain >= ASM_CLEAN_EXTEND_GAIN * fabs(score);
}

static int ac_clean_continue(int round, int rounds, int extend, double elapsed)
{
    return round < rounds || (extend && round < ASM_CLEAN_MAX_ROUNDS && elapsed < ASM_CLEAN_EXTEND_SECONDS);
}

int AsmConflict_clean(AsmRun *run, const AsmConflictOpts *o, const AsmPoseOpts *pose_opts,
                      AsmConflictReport *final_report, AsmConflictStats *st)
{
    memset(st, 0, sizeof *st);
    double t0 = ves_clock_sec();
    AsmPoseStats ps;
    AsmConflictReport rep;
    double *chart_mass = ARENA_CALLOC(run->arena, run->n_charts ? run->n_charts : 1, sizeof(double));
    uint8_t *protect = ARENA_CALLOC(run->arena, run->n_rels ? run->n_rels : 1, sizeof(uint8_t));   /* bridges the split guard refused */
    if (AsmConflict_audit(run, o, chart_mass, &rep) != 0) return -1;
    st->score_first = rep.score; st->contra_first = rep.contra_area; st->covered_first = rep.covered_area; st->pairs_first = rep.pairs;
    double score = rep.score;
    AcSnap cur, best;
    ac_snap_alloc(run->arena, run, &cur);
    ac_snap_alloc(run->arena, run, &best);
    ac_snap_take(run, &best, score);
    size_t batch_cap = 0;
    int round;
    /* Past o->rounds, cleaning continues only while the last round was kept
     * and gained at least ASM_CLEAN_EXTEND_GAIN of the score, within the
     * extension's time and round ceilings: a fixed round count sized for the
     * 10x10x10 (whose 8th round gained 0.06%) stopped the 21x21x21 while its
     * 8th round still gained 1.5%. */
    int extend = 0;
    for (round = 0; ac_clean_continue(round, o->rounds, extend, ves_clock_sec() - t0); round++) {
        if (rep.contra_area <= 0.0) break;
        extend = 0;
        /* candidate relations by witness mass */
        size_t nflag = 0;
        for (size_t r = 0; r < run->n_rels; r++) if (ac_rel_active(&run->rels[r]) && run->rels[r].flag_mass > 0.0 && !protect[r]) nflag++;
        if (nflag == 0) break;
        Arena_Mark mark = Arena_save(run->arena);
        AcFlag *fl = ARENA_ALLOC(run->arena, nflag * sizeof(AcFlag));
        size_t q = 0;
        for (size_t r = 0; r < run->n_rels; r++) {
            if (!ac_rel_active(&run->rels[r]) || run->rels[r].flag_mass <= 0.0 || protect[r]) continue;
            /* weak relations (low robust weight) are implicated more strongly */
            fl[q].mass = run->rels[r].flag_mass * (2.0 - run->rels[r].robust_w);
            fl[q].rel = (int32_t)r; q++;
        }
        qsort(fl, nflag, sizeof(AcFlag), ac_cmp_flag_desc);
        size_t batch = (size_t)ceil(o->drop_frac * (double)nflag);
        if (batch < 1) batch = 1;
        if (batch_cap > 0 && batch > batch_cap) batch = batch_cap;
        /* only relations carrying at least a quarter of the top mass */
        size_t nb = 0;
        for (size_t k = 0; k < batch && k < nflag; k++) { if (fl[k].mass < 0.25 * fl[0].mass) break; nb++; }
        if (nb < 1) nb = 1;
        ac_snap_take(run, &cur, score);
        for (size_t k = 0; k < nb; k++) run->rels[fl[k].rel].flags |= ASM_REL_DROPPED;
        AsmPose_solve(run, pose_opts, &ps);
        AsmConflictReport rep2;
        AsmConflict_audit(run, o, chart_mass, &rep2);
        /* SPLIT GUARD: the area this batch separated from its components against the score it gained */
        double split_area = ac_split_area(run, &cur);
        double gain = rep2.score - score;
        int split_refused = ASM_CLEAN_SPLIT_MIN_GAIN > 0.0 && split_area > 0.0 && gain < ASM_CLEAN_SPLIT_MIN_GAIN * split_area;
        int score_fell = rep2.score + 1e-9 < score;
        fprintf(stderr, "  [clean2 round %d] dropped %zu of %zu flagged: score %.4e -> %.4e (covered %.4e -> %.4e, contra %.4e -> %.4e, pairs %zu -> %zu, components %zu, separated %.3e vox^2)%s\n",
                round + 1, nb, nflag, score, rep2.score, rep.covered_area, rep2.covered_area, rep.contra_area, rep2.contra_area, rep.pairs, rep2.pairs, ps.components, split_area,
                split_refused ? " -> REFUSED: the gain does not pay for the split" : (score_fell ? " -> reverted: the score fell" : ""));
        size_t nprot = 0;
        for (size_t k = 0; k < nb; k++) {
            const AsmRelation *rl = &run->rels[fl[k].rel];
            int apart = run->charts[rl->a].component != run->charts[rl->b].component;
            fprintf(stderr, "  [clean2 round %d]   rel %d charts %d-%d: seam %.0f vox, %zu corr, rms %.2f, robust w %.2f, witness mass %.3e%s\n",
                    round + 1, fl[k].rel, rl->a, rl->b, rl->seam_len, rl->n_corr, rl->rms, rl->robust_w, fl[k].mass,
                    apart ? (split_refused ? " -> a bridge: its charts parted; PROTECTED" : " -> a bridge: its charts parted") : "");
            if (split_refused && apart) { protect[fl[k].rel] = 1; nprot++; }
        }
        Arena_restore(run->arena, mark);
        if (split_refused || score_fell) {
            ac_snap_restore(run, &cur);   /* exact: the layout before the batch, not a re-solve */
            st->rels_reverted += nb;
            if (split_refused) { st->splits_refused++; st->split_area_refused += split_area; }
            batch_cap = nb / 2;
            AsmConflict_audit(run, o, chart_mass, &rep);   /* rebuilds the cell list of the same layout */
            if (batch_cap == 0) {
                if (!split_refused || nprot == 0) break;
                batch_cap = 1;   /* the refused bridge is protected: the next candidate gets its turn */
            }
            continue;
        }
        st->rels_dropped += nb;
        extend = ac_clean_extend(gain, score);
        score = rep2.score;
        batch_cap = 0;
        AsmConflict_audit(run, o, chart_mass, &rep);
        if (score > best.score) ac_snap_take(run, &best, score);
    }
    st->rounds = round;
    /* SINGLE drops: the batch rounds stop at the first batch that lowers the score, and on a
     * graph with a wrong bridge join the solve is multi-modal enough that a batch of three
     * fails where one of them alone would not.  Every still-flagged relation is tried alone,
     * heaviest first, with the exact revert; a drop is kept when the score rises.  Bounded by
     * the flagged count and by consecutive failures (each try is one pose solve + audit). */
#if ASM_CLEAN_SINGLE_DROPS
    if (rep.contra_area > 0.0) {
        size_t nflag = 0;
        for (size_t r = 0; r < run->n_rels; r++) if (ac_rel_active(&run->rels[r]) && run->rels[r].flag_mass > 0.0) nflag++;
        Arena_Mark mark = Arena_save(run->arena);
        AcFlag *fl = ARENA_ALLOC(run->arena, (nflag ? nflag : 1) * sizeof(AcFlag));
        size_t q = 0;
        for (size_t r = 0; r < run->n_rels; r++) {
            if (!ac_rel_active(&run->rels[r]) || run->rels[r].flag_mass <= 0.0) continue;
            fl[q].mass = run->rels[r].flag_mass * (2.0 - run->rels[r].robust_w);
            fl[q].rel = (int32_t)r; q++;
        }
        qsort(fl, nflag, sizeof(AcFlag), ac_cmp_flag_desc);
        size_t tried = 0, kept = 0, fails = 0;
        for (size_t k = 0; k < nflag && fails < 12 && tried < 64; k++) {
            if (!ac_rel_active(&run->rels[fl[k].rel])) continue;   /* dropped or switched by a previous keep */
            ac_snap_take(run, &cur, score);
            run->rels[fl[k].rel].flags |= ASM_REL_DROPPED;
            AsmPose_solve(run, pose_opts, &ps);
            AsmConflictReport rep2;
            AsmConflict_audit(run, o, chart_mass, &rep2);
            tried++;
            if (rep2.score > score + 1e-9) {
                kept++; fails = 0;
                st->rels_dropped++;
                score = rep2.score;
                if (score > best.score) ac_snap_take(run, &best, score);
                fprintf(stderr, "  [clean2 single] dropped relation %d (charts %d-%d, mass %.3e): score %.4e contra %.4e pairs %zu\n",
                        fl[k].rel, run->rels[fl[k].rel].a, run->rels[fl[k].rel].b, fl[k].mass, rep2.score, rep2.contra_area, rep2.pairs);
                if (rep2.contra_area <= 0.0) break;
            } else {
                fails++;
                ac_snap_restore(run, &cur);
            }
        }
        Arena_restore(run->arena, mark);
        AsmConflict_audit(run, o, chart_mass, &rep);
        fprintf(stderr, "  [clean2 single] %zu flagged, tried %zu, kept %zu: score %.4e contra %.4e pairs %zu\n", nflag, tried, kept, rep.score, rep.contra_area, rep.pairs);
    }
#endif
    /* persistently contradicting charts -> extras, each removal kept only if the score holds */
    for (int pass = 0; pass < 32; pass++) {
        int32_t worst = -1; double worst_ratio = o->extras_ratio;
        for (size_t i = 0; i < run->n_charts; i++) {
            const AsmChart *c = &run->charts[i];
            if (!AsmChart_in_layout(c) || !c->placed || c->area3d <= 0.0) continue;
            double ratio = chart_mass[i] / c->area3d;
            if (ratio > worst_ratio) { worst_ratio = ratio; worst = (int32_t)i; }
        }
        if (worst < 0) break;
        ac_snap_take(run, &cur, score);
        run->charts[worst].flags |= ASM_CHART_EXTRAS;
        AsmPose_solve(run, pose_opts, &ps);
        AsmConflict_audit(run, o, chart_mass, &rep);
        /* UNCONDITIONAL (measured 2026-09-08): the extras pass is what removes the fused slabs
         * that make the big component non-sheet-like; on the familiar it dips first (7.7043e6 ->
         * 7.6973e6 on the first chart) and then climbs to 8.07e6 with contra 1.9e5 -> 4.5e3, and
         * placement goes from a 0.76M primary / 49% placed to a 2.55M primary / 88.5%.  Gating each
         * removal on the score stopped it at the dip. */
        fprintf(stderr, "  [clean2 extras] chart %d (cube %d, area %.3e, ratio %.2f) -> extras: score %.4e contra %.4e pairs %zu\n",
                worst, run->charts[worst].cube, run->charts[worst].area3d, worst_ratio, rep.score, rep.contra_area, rep.pairs);
        st->charts_extras++;
        score = rep.score;
        if (score > best.score) ac_snap_take(run, &best, score);
    }
    if (ASM_CLEAN_KEEP_BEST && best.score > score + 1e-9) {
        fprintf(stderr, "  [clean2] finishing on the best audited state: score %.4e (the last state scored %.4e)\n", best.score, score);
        ac_snap_restore(run, &best);
        AsmConflict_audit(run, o, chart_mass, &rep);
        score = rep.score;
    }
    st->score_last = rep.score; st->contra_last = rep.contra_area; st->covered_last = rep.covered_area; st->pairs_last = rep.pairs;
    st->sec = ves_clock_sec() - t0;
    *final_report = rep;
    return 0;
}

/* ac_witness_path resets only what it visited.  With one shared scratch, every
 * path must equal a search from freshly reset scratch, including a pair past
 * the depth cap and a pair in another component, and the scratch must come
 * back all -2 after every call. */
static int ac_witness_reset_control(void)
{
    enum { NC = 40, NR = 44 };
    int fails = 0;
    AsmRun run; memset(&run, 0, sizeof run); run.arena = Arena_new();
    AsmChart blank; memset(&blank, 0, sizeof blank);
    for (int i = 0; i < NC; i++) { blank.id = i; AsmRun_push_chart(&run, &blank); }
    /* a 30-chart chain with four chords, and a separate 10-chart chain */
    for (int i = 0; i + 1 < 30; i++) { AsmRelation r; memset(&r, 0, sizeof r); r.a = i; r.b = i + 1; AsmRun_push_rel(&run, &r); }
    const int chord[4][2] = { {0, 7}, {3, 19}, {12, 25}, {5, 28} };
    for (int k = 0; k < 4; k++) { AsmRelation r; memset(&r, 0, sizeof r); r.a = chord[k][0]; r.b = chord[k][1]; AsmRun_push_rel(&run, &r); }
    for (int i = 30; i + 1 < NC; i++) { AsmRelation r; memset(&r, 0, sizeof r); r.a = i; r.b = i + 1; AsmRun_push_rel(&run, &r); }
    int32_t deg[NC + 1] = {0}, fill[NC] = {0}, adj[2 * NR];
    for (size_t r = 0; r < run.n_rels; r++) { deg[run.rels[r].a + 1]++; deg[run.rels[r].b + 1]++; }
    for (int i = 0; i < NC; i++) deg[i + 1] += deg[i];
    for (size_t r = 0; r < run.n_rels; r++) {
        adj[deg[run.rels[r].a] + fill[run.rels[r].a]++] = (int32_t)r;
        adj[deg[run.rels[r].b] + fill[run.rels[r].b]++] = (int32_t)r;
    }
    int32_t shared_prev[NC], shared_queue[2 * NC], fresh_prev[NC], fresh_queue[2 * NC], p1[16], p2[16];
    for (int i = 0; i < NC; i++) shared_prev[i] = -2;
    const int pairs[][2] = { {0, 29}, {2, 20}, {29, 0}, {8, 9}, {0, 39}, {31, 38}, {14, 27}, {1, 26}, {6, 6 + 13} };
    for (size_t q = 0; q < sizeof pairs / sizeof pairs[0]; q++) {
        for (int i = 0; i < NC; i++) fresh_prev[i] = -2;
        int n_fresh = ac_witness_path(&run, deg, adj, pairs[q][0], pairs[q][1], 12, fresh_prev, fresh_queue, p2);
        int n_shared = ac_witness_path(&run, deg, adj, pairs[q][0], pairs[q][1], 12, shared_prev, shared_queue, p1);
        int same = n_fresh == n_shared && !memcmp(p1, p2, (size_t)n_fresh * sizeof(int32_t));
        int reset = 1;
        for (int i = 0; i < NC; i++) reset &= shared_prev[i] == -2;
        if (!same || !reset) {
            fprintf(stderr, "  asm_conflict selftest FAIL: witness path %d-%d shared scratch (length %d vs %d, reset %d)\n",
                    pairs[q][0], pairs[q][1], n_shared, n_fresh, reset);
            fails++;
        }
    }
    /* other components and pairs past the depth cap have no path */
    for (int i = 0; i < NC; i++) shared_prev[i] = -2;
    if (ac_witness_path(&run, deg, adj, 0, 39, 12, shared_prev, shared_queue, p1) != 0) { fprintf(stderr, "  asm_conflict selftest FAIL: witness path across components\n"); fails++; }
    if (ac_witness_path(&run, deg, adj, 8, 22, 2, shared_prev, shared_queue, p1) != 0) { fprintf(stderr, "  asm_conflict selftest FAIL: witness path past the depth cap\n"); fails++; }
    Arena_dispose(&run.arena);
    fprintf(stderr, "  asm_conflict witness scratch reset: %s (%d failures)\n", fails ? "FAIL" : "ok", fails);
    return fails;
}

int AsmConflict_selftest(void)
{
    int fails = 0;
    /* the round budget on measured rounds: the 10x10x10's 8th round (8.3030e7 -> 8.3080e7)
     * ends cleaning at the configured 8; the layer-cut 21x21x21's (7.7675e8 -> 7.8835e8)
     * continues, but never past the round or time ceiling, and a refused round stops it */
    if (ac_clean_extend(8.3080e7 - 8.3030e7, 8.3030e7) || !ac_clean_extend(7.8835e8 - 7.7675e8, 7.7675e8) ||
        ac_clean_extend(-1.0, 1.0e6) || ac_clean_extend(NAN, 1.0e6) ||
        !ac_clean_continue(7, 8, 0, 0.0) || ac_clean_continue(8, 8, 0, 0.0) || !ac_clean_continue(8, 8, 1, 10.0) ||
        ac_clean_continue(ASM_CLEAN_MAX_ROUNDS, 8, 1, 10.0) || ac_clean_continue(9, 8, 1, ASM_CLEAN_EXTEND_SECONDS + 1.0)) {
        fprintf(stderr, "  asm_conflict selftest FAIL: clean2 round budget\n"); fails++;
    }
    /* two square charts placed on top of each other, 12 vox apart along z (normal):
     * the audit must report a contradiction; side by side: none */
    AsmRun run; memset(&run, 0, sizeof run); run.arena = Arena_new();
    float uv[8] = { 0,0, 40,0, 40,40, 0,40 };
    int32_t faces[6] = { 0,1,2, 0,2,3 };
    float nrm[12] = { 1,0,0, 1,0,0, 1,0,0, 1,0,0 };
    float xyz_a[12] = { 0,0,0, 0,0,40, 0,40,40, 0,40,0 };
    float xyz_b[12] = { 12,0,0, 12,0,40, 12,40,40, 12,40,0 };
    for (int i = 0; i < 2; i++) {
        AsmChart c; memset(&c, 0, sizeof c);
        c.id = i; c.nv = 4; c.nf = 2; c.uv = uv; c.faces = faces; c.nrm = nrm; c.xyz = i ? xyz_b : xyz_a;
        c.area3d = 1600.0; c.placed = 1; c.component = 0; c.pose_x = 0.0; c.pose_y = 0.0;
        AsmRun_push_chart(&run, &c);
    }
    run.layer_d_global = 12.0;
    AsmConflictOpts o; AsmConflict_default_opts(&o);
    AsmConflictReport rep;
    AsmConflict_audit(&run, &o, NULL, &rep);
    if (!(rep.contra_area > 1000.0 && rep.pairs == 1)) { fprintf(stderr, "  asm_conflict selftest FAIL: stacked charts contra %.1f pairs %zu\n", rep.contra_area, rep.pairs); fails++; }
    {
        Arena_T arena = Arena_new();
        AsmConflictGrid *grid = AsmConflictGrid_new(arena,&o,&run);
        AsmConflictGrid_insert(grid,&run.charts[0]); AsmConflictGrid_insert(grid,&run.charts[1]);
        const AsmChart *candidate = &run.charts[0]; int32_t partner; double share;
        double baseline = AsmConflictGrid_probe(grid,&candidate,1,&partner,&share);
        uint8_t excluded[2] = {0,1};
        double without_old = AsmConflictGrid_probe_excluding(grid,&candidate,1,excluded,&partner,&share);
        excluded[0] = 1; excluded[1] = 0;
        double with_other = AsmConflictGrid_probe_excluding(grid,&candidate,1,excluded,&partner,&share);
        if (!(baseline > 1000 && without_old == 0 && with_other == baseline && partner == 1 && share == baseline &&
              AsmConflictGrid_probe(grid,&candidate,1,NULL,NULL) == baseline)) {
            fprintf(stderr,"  asm_conflict selftest FAIL: relocation sample exclusion\n"); fails++;
        }
        AsmConflictGrid_dispose(grid); Arena_dispose(&arena);
    }
    run.charts[1].pose_x = 100.0;   /* side by side */
    AsmConflict_audit(&run, &o, NULL, &rep);
    if (rep.contra_area != 0.0 || rep.pairs != 0) { fprintf(stderr, "  asm_conflict selftest FAIL: separated charts contra %.1f\n", rep.contra_area); fails++; }
    if (fabs(rep.covered_area - 2.0 * 1600.0) > 400.0) { fprintf(stderr, "  asm_conflict selftest FAIL: covered %.1f\n", rep.covered_area); fails++; }
    /* same place, same layer (duplicate coverage): no contradiction */
    run.charts[1].pose_x = 0.0; run.charts[1].xyz = xyz_a;
    AsmConflict_audit(&run, &o, NULL, &rep);
    if (rep.contra_area != 0.0) { fprintf(stderr, "  asm_conflict selftest FAIL: duplicate coverage flagged %.1f\n", rep.contra_area); fails++; }
    /* layout-confirmed joins: chart 1 continues chart 0 in uv (pose_x 40) and in 3-D (x + 40,
     * same plane) with no relation between them -> joined by the layout; one layer above -> not */
    {
        float xyz_c[12] = { 0,0,40, 0,0,80, 0,40,80, 0,40,40 };
        float xyz_d[12] = { 12,0,40, 12,0,80, 12,40,80, 12,40,40 };
        run.charts[1].pose_x = 40.0; run.charts[1].xyz = xyz_c;
        run.n_rels = 0;
        AsmJoinStats js;
        AsmConflict_confirm_joins(&run, &o, NULL, &js);
        if (js.pairs_joined != 1 || run.n_rels != 1 || !(run.rels[0].flags & ASM_REL_LAYOUT) ||
            fabs(run.rels[0].theta) > 1e-9 || fabs(run.rels[0].tx - 40.0) > 1e-6 || fabs(run.rels[0].ty) > 1e-6 ||
            js.lineages_before != 2 || js.lineages_after != 1) {
            fprintf(stderr, "  asm_conflict selftest FAIL: layout join (joined %zu rels %zu theta %.4f t %.2f %.2f lineages %zu -> %zu)\n",
                    js.pairs_joined, run.n_rels, run.n_rels ? run.rels[0].theta : 0.0, run.n_rels ? run.rels[0].tx : 0.0,
                    run.n_rels ? run.rels[0].ty : 0.0, js.lineages_before, js.lineages_after);
            fails++;
        }
        run.n_rels = 0;
        run.charts[1].xyz = xyz_d;
        AsmConflict_confirm_joins(&run, &o, NULL, &js);
        if (js.pairs_joined != 0 || run.n_rels != 0) { fprintf(stderr, "  asm_conflict selftest FAIL: a layer above must not join (joined %zu)\n", js.pairs_joined); fails++; }
        /* alignment: chart 1's material continues chart 0 at x + 40 but is drawn at pose_x 43;
         * the join measures the 3-vox offset and moves it back to 40 */
        run.n_rels = 0;
        run.charts[1].xyz = xyz_c; run.charts[1].pose_x = 43.0;
        AsmConflict_confirm_joins(&run, &o, NULL, &js);
        if (ASM_JOIN_ALIGN && (js.pairs_joined != 1 || fabs(run.charts[1].pose_x - 40.0) > 0.6 || fabs(run.charts[1].pose_y) > 0.6 || fabs(run.rels[0].tx - run.charts[1].pose_x) > 1e-6)) {
            fprintf(stderr, "  asm_conflict selftest FAIL: alignment (joined %zu, pose_x %.3f pose_y %.3f, tx %.3f, shifted %zu max %.2f)\n",
                    js.pairs_joined, run.charts[1].pose_x, run.charts[1].pose_y, run.n_rels ? run.rels[0].tx : 0.0, js.charts_shifted, js.shift_max);
            fails++;
        }
        run.charts[1].pose_x = 40.0; run.charts[1].pose_y = 0.0;
        run.n_rels = 0;
        /* a mirrored chart a stores the inverse of the effective relation */
        run.charts[1].xyz = xyz_c;
        run.charts[0].flags |= ASM_CHART_MIRROR; run.charts[0].pose_x = 40.0; run.charts[1].pose_x = 0.0; run.charts[1].flags |= ASM_CHART_MIRROR;
        AsmConflict_confirm_joins(&run, &o, NULL, &js);
        if (js.pairs_joined != 1 || fabs(run.rels[0].tx - 40.0) > 1e-6 || (run.rels[0].flags & ASM_REL_PARITY)) {
            fprintf(stderr, "  asm_conflict selftest FAIL: mirrored layout join (joined %zu tx %.2f flags %u)\n", js.pairs_joined, run.n_rels ? run.rels[0].tx : 0.0, run.n_rels ? run.rels[0].flags : 0u);
            fails++;
        }
        run.charts[0].flags &= ~(uint32_t)ASM_CHART_MIRROR; run.charts[1].flags &= ~(uint32_t)ASM_CHART_MIRROR;
        run.n_rels = 0;
    }
    Arena_dispose(&run.arena);
    /* SEAM RE-SOLVE (2026-09-10): chart 1 continues chart 0 (x + 40 in 3-D) but is drawn 20 vox further
     * along u with a WEAK seam between them (8 correspondences at the shared edge): readmitted, the slit
     * closes, lineages 2 -> 1; drawn 200 vox off: refused by the gap; a seam that pulls chart 1 onto
     * a third chart one layer away: the contradiction grows and the step is REVERTED */
    {
        AsmRun run5; memset(&run5, 0, sizeof run5); run5.arena = Arena_new();
        float uv5[8] = { 0,0, 80,0, 80,80, 0,80 };
        float xyz5a[12] = { 0,0,0, 0,0,80, 0,80,80, 0,80,0 }, xyz5b[12] = { 0,0,80, 0,0,160, 0,80,160, 0,80,80 }, xyz5c[12] = { 12,0,80, 12,0,160, 12,80,160, 12,80,80 };
        for (int i = 0; i < 3; i++) {
            AsmChart ch; memset(&ch, 0, sizeof ch);
            ch.id = i; ch.nv = 4; ch.nf = 2; ch.uv = uv5; ch.faces = faces; ch.nrm = nrm; ch.xyz = i == 0 ? xyz5a : (i == 1 ? xyz5b : xyz5c);
            ch.area3d = 6400.0 - 100.0 * i; ch.placed = 1; ch.component = 0; ch.pose_x = i == 0 ? 0.0 : (i == 1 ? 100.0 : 260.0); ch.pose_y = 0.0;
            AsmRun_push_chart(&run5, &ch);
        }
        run5.layer_d_global = 12.0;
        /* the seam 0-1: chart 1's left edge (vertices 0, 3) meets chart 0's right edge (1, 2); the true transform is tx 80 */
        AsmRelation r; memset(&r, 0, sizeof r);
        r.a = 0; r.b = 1; r.theta = 0.0; r.tx = 80.0; r.ty = 0.0; r.w_xy = 32.0; r.w_theta = 32.0 * 3200.0; r.robust_w = 1.0; r.rms = 0.5; r.seam_len = 80.0; r.n_corr = 8; r.flags = ASM_REL_WEAK;
        r.corr_first = (int32_t)run5.n_corr; r.corr_count = 8;
        for (int k = 0; k < 8; k++) { AsmCorr cc = { (k & 1) ? 2 : 1, (k & 1) ? 3 : 0 }; AsmRun_push_corr(&run5, &cc); }
        AsmRun_push_rel(&run5, &r);
        AsmPoseOpts po; AsmPose_default_opts(&po);
        AsmRefineStats rs;
        AsmConflict_refine_seams(&run5, &o, &po, NULL, &rs);
        if (!(rs.readmitted == 1 && rs.kept && !rs.reverted && fabs(run5.charts[1].pose_x - 80.0) < 0.5 && rs.lineages_before == 3 && rs.lineages_after == 2 && (run5.rels[0].flags & ASM_REL_READMIT))) {
            fprintf(stderr, "  asm_conflict selftest FAIL: seam re-solve (readmitted %zu kept %d reverted %d pose_x %.2f lineages %zu -> %zu flags %u)\n",
                    rs.readmitted, rs.kept, rs.reverted, run5.charts[1].pose_x, rs.lineages_before, rs.lineages_after, run5.rels[0].flags);
            fails++;
        }
        /* far off: refused by the gap, nothing moves */
        run5.rels[0].flags = ASM_REL_WEAK; run5.charts[1].pose_x = 380.0;
        AsmConflict_refine_seams(&run5, &o, &po, NULL, &rs);
        if (!(rs.readmitted == 0 && rs.refused_gap == 1 && fabs(run5.charts[1].pose_x - 380.0) < 1e-9)) {
            fprintf(stderr, "  asm_conflict selftest FAIL: seam re-solve far pair (readmitted %zu refused gap %zu pose_x %.2f)\n", rs.readmitted, rs.refused_gap, run5.charts[1].pose_x); fails++;
        }
        /* chart 1 drawn 20 vox off continuity (readmitted) but its seam claims it belongs where chart 2 (a
         * layer away) is drawn: the re-solve stacks them, the contradiction grows, the step is reverted */
        run5.charts[1].pose_x = 100.0; run5.charts[2].pose_x = 260.0;
        run5.rels[0].tx = 260.0; run5.rels[0].flags = ASM_REL_WEAK;
        AsmConflict_refine_seams(&run5, &o, &po, NULL, &rs);
        if (!(rs.readmitted == 1 && (rs.reverted || rs.unreadmitted == 1) && rs.rounds == 2 && fabs(run5.charts[1].pose_x - 100.0) < 1e-9 && !(run5.rels[0].flags & ASM_REL_READMIT) && rs.contra_after < 1.0)) {
            fprintf(stderr, "  asm_conflict selftest FAIL: seam re-solve revert (readmitted %zu kept %d reverted %d un-readmitted %zu rounds %d pose_x %.2f flags %u contra %.0f -> %.0f)\n",
                    rs.readmitted, rs.kept, rs.reverted, rs.unreadmitted, rs.rounds, run5.charts[1].pose_x, run5.rels[0].flags, rs.contra_before, rs.contra_after); fails++;
        }
        Arena_dispose(&run5.arena);
    }
    fails += ac_witness_reset_control();
    if (fails == 0) fprintf(stderr, "  asm_conflict selftest: all passed\n");
    return fails;
}
