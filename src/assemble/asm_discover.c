/* asm_discover.c -- the second pass over a layout.  See asm_discover.h. */
#include "asm_discover.h"
#include "asm_continuity.h"
#include "asm_report.h"
#include "asm_store.h"
#include "asm_audit.h"
#include "asm_field.h"
#include "../common/pipeline_constants.h"
#include "../common/ves_platform.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define AD_PI 3.14159265358979323846
#define AD_MAX_FITS 32          /* partner poses considered per chart (the best-supported ones) */
#define AD_SAMPLES 256          /* vertices a pose is applied to when two poses are compared */

/* A rigid pose in the global frame: g = R(theta) * m(uv) + (x, y), where m
 * mirrors u when `mirror` -- AsmChart_point's convention. */
typedef struct AdPose { double theta, x, y; int mirror; } AdPose;

static void ad_point(const AdPose *p, double u, double v, double g[2])
{
    if (p->mirror) u = -u;
    double ct = cos(p->theta), sn = sin(p->theta);
    g[0] = ct * u - sn * v + p->x;
    g[1] = sn * u + ct * v + p->y;
}

/* The stored transform is in unmirrored frames, p_a = R(theta) p_b + t; with
 * chart a mirrored the effective transform is (-theta, -tx, ty) (asm_pose.c). */
static void ad_effective(const AsmRelation *r, int mirror_a, double *theta, double *tx, double *ty)
{
    if (!mirror_a) { *theta = r->theta; *tx = r->tx; *ty = r->ty; }
    else { *theta = -r->theta; *tx = -r->tx; *ty = r->ty; }
}

/* The pose of chart `chart` implied by relation r and the pose of its other
 * endpoint.  Returns that partner. */
static int32_t ad_implied(const AsmRun *run, const AsmRelation *r, int32_t chart, AdPose *out)
{
    int32_t partner = r->a == chart ? r->b : r->a;
    const AsmChart *pc = run->charts + partner;
    int mirror_p = (pc->flags & ASM_CHART_MIRROR) != 0;
    int parity = (r->flags & ASM_REL_PARITY) != 0;
    int mirror_c = mirror_p ^ parity;
    double th = 0.0, tx = 0.0, ty = 0.0;
    if (r->a == partner) {
        /* m_P(p) = R(th) m_C(q) + t  ->  g = R(theta_P + th) m_C(q) + R(theta_P) t + t_P */
        ad_effective(r, mirror_p, &th, &tx, &ty);
        double ct = cos(pc->pose_theta), sn = sin(pc->pose_theta);
        out->theta = pc->pose_theta + th;
        out->x = pc->pose_x + ct * tx - sn * ty;
        out->y = pc->pose_y + sn * tx + ct * ty;
    } else {
        /* m_C(q) = R(th) m_P(p) + t  ->  g = R(theta_P - th) m_C(q) - R(theta_P - th) t + t_P */
        ad_effective(r, mirror_c, &th, &tx, &ty);
        out->theta = pc->pose_theta - th;
        double ct = cos(out->theta), sn = sin(out->theta);
        out->x = pc->pose_x - (ct * tx - sn * ty);
        out->y = pc->pose_y - (sn * tx + ct * ty);
    }
    out->mirror = mirror_c;
    return partner;
}

/* Relations that can carry discovery evidence: not dropped by the cleaning
 * loop and not tangency-rejected.  "Soft" ones (switched off by the pose
 * solve, or placement evidence only) may support an agreeing set but never
 * veto it. */
static int ad_usable(const AsmRelation *r)
{
    if (ASM_WRAP_VETO_DISCOVER && (r->flags & ASM_REL_CROSSWRAP)) return 0;   /* a sheet switch inside the seam is not a hypothesis */
    return (r->flags & (ASM_REL_DROPPED | ASM_REL_CONTACT)) == 0u && isfinite(r->theta) && isfinite(r->tx) && isfinite(r->ty);
}
static int ad_soft(const AsmRelation *r)
{
    return (r->flags & (ASM_REL_SWITCHED | ASM_REL_PLACEMENT_ONLY)) != 0u;
}
/* The same predicate the repair's admission uses for its source obligations. */
static int ad_source(const AsmRelation *r)
{
    return !(r->flags & ASM_REL_CONTACT) && (r->continuity & (ASM_CONT_SOURCE | ASM_CONT_CUT));
}

/* ---- registered-vertex index: global uv cells of ASM_DISCOVER_OVERLAP_VOX ---- */

typedef struct AdEntry { double u, v; int32_t chart, vert; int64_t key; } AdEntry;
typedef struct AdIndex { AdEntry *e; size_t n; double cell; } AdIndex;

static int64_t ad_key(double u, double v, double cell)
{
    long kx = (long)floor(u / cell), ky = (long)floor(v / cell);
    return ((int64_t)(kx + (1L << 24)) << 26) | (int64_t)(ky + (1L << 24));
}
static int ad_entry_order(const void *a, const void *b)
{
    const AdEntry *x = a, *y = b;
    if (x->key != y->key) return x->key < y->key ? -1 : 1;
    if (x->chart != y->chart) return x->chart < y->chart ? -1 : 1;
    return x->vert < y->vert ? -1 : x->vert > y->vert;
}
static int ad_index_build(Arena_T arena, const AsmRun *run, AdIndex *ix)
{
    size_t total = 0;
    for (size_t c = 0; c < run->n_charts; c++) if (AsmChart_registered(run->charts + c)) total += run->charts[c].nv;
    ix->cell = ASM_DISCOVER_OVERLAP_VOX; ix->n = 0;
    ix->e = ARENA_ALLOC(arena, (total ? total : 1) * sizeof(AdEntry));
    if (!ix->e) return -1;
    for (size_t c = 0; c < run->n_charts; c++) {
        const AsmChart *ch = run->charts + c;
        if (!AsmChart_registered(ch)) continue;
        for (size_t k = 0; k < ch->nv; k++) {
            double p[2] = {0, 0};
            AsmChart_point(ch, k, p);
            if (!isfinite(p[0]) || !isfinite(p[1])) continue;
            AdEntry *e = ix->e + ix->n++;
            e->u = p[0]; e->v = p[1]; e->chart = (int32_t)c; e->vert = (int32_t)k; e->key = ad_key(p[0], p[1], ix->cell);
        }
    }
    qsort(ix->e, ix->n, sizeof(AdEntry), ad_entry_order);
    return 0;
}
static size_t ad_lower(const AdIndex *ix, int64_t key)
{
    size_t lo = 0, hi = ix->n;
    while (lo < hi) { size_t mid = lo + (hi - lo) / 2; if (ix->e[mid].key < key) lo = mid + 1; else hi = mid; }
    return lo;
}

/* Registered vertices within `radius` (uv) of (u, v): their count, the one
 * nearest in uv, and the smallest 3-D distance from `xyz` to any of them. */
typedef struct AdNear { size_t hits; double uv_min, d3_min; int32_t chart, vert; } AdNear;
static void ad_near(const AdIndex *ix, const AsmRun *run, double u, double v, const float *xyz, double radius, AdNear *out)
{
    out->hits = 0; out->uv_min = DBL_MAX; out->d3_min = DBL_MAX; out->chart = -1; out->vert = -1;
    long kx0 = (long)floor((u - radius) / ix->cell), kx1 = (long)floor((u + radius) / ix->cell);
    long ky0 = (long)floor((v - radius) / ix->cell), ky1 = (long)floor((v + radius) / ix->cell);
    double r2 = radius * radius;
    for (long kx = kx0; kx <= kx1; kx++) for (long ky = ky0; ky <= ky1; ky++) {
        int64_t key = ((int64_t)(kx + (1L << 24)) << 26) | (int64_t)(ky + (1L << 24));
        for (size_t i = ad_lower(ix, key); i < ix->n && ix->e[i].key == key; i++) {
            const AdEntry *e = ix->e + i;
            double du = e->u - u, dv = e->v - v, d2 = du * du + dv * dv;
            if (d2 > r2) continue;
            out->hits++;
            const float *q = run->charts[e->chart].xyz + 3 * (size_t)e->vert;
            double dz = q[0] - xyz[0], dy = q[1] - xyz[1], dx = q[2] - xyz[2], d3 = sqrt(dz * dz + dy * dy + dx * dx);
            if (d2 < out->uv_min) { out->uv_min = d2; out->chart = e->chart; out->vert = e->vert; }
            if (d3 < out->d3_min) out->d3_min = d3;
        }
    }
    if (out->hits) out->uv_min = sqrt(out->uv_min);
}

/* ---- per-chart evidence ------------------------------------------------------ */

typedef struct AdFit { AdPose pose; int32_t rel, partner; int soft, quad; double rms; size_t corr; } AdFit;

static size_t ad_samples(const AsmChart *ch, size_t *idx, size_t cap)
{
    size_t n = 0, stride = ch->nv / cap ? ch->nv / cap : 1;
    for (size_t k = 0; k < ch->nv && n < cap; k += stride) idx[n++] = k;
    return n;
}
/* the largest vertex displacement between two poses of the chart, over the samples */
static double ad_displacement(const AsmChart *ch, const AdPose *a, const AdPose *b, const size_t *idx, size_t n)
{
    double worst = 0.0;
    for (size_t i = 0; i < n; i++) {
        double u = ch->uv[2 * idx[i]], v = ch->uv[2 * idx[i] + 1], ga[2], gb[2];
        ad_point(a, u, v, ga); ad_point(b, u, v, gb);
        double d = hypot(ga[0] - gb[0], ga[1] - gb[1]);
        if (d > worst) worst = d;
    }
    return worst;
}
static void ad_centroid(const AsmChart *ch, double c[2])
{
    c[0] = c[1] = 0.0;
    for (size_t k = 0; k < ch->nv; k++) { c[0] += ch->uv[2 * k]; c[1] += ch->uv[2 * k + 1]; }
    if (ch->nv) { c[0] /= (double)ch->nv; c[1] /= (double)ch->nv; }
}
/* which quadrant of the chart's own frame the partner sits in, under `pose` */
static int ad_quadrant(const AsmRun *run, const double centroid[2], const AdPose *pose, int32_t partner)
{
    const AsmChart *pc = run->charts + partner;
    double pcen[2], gp[2], gc[2];
    ad_centroid(pc, pcen);
    AdPose pp = { pc->pose_theta, pc->pose_x, pc->pose_y, (pc->flags & ASM_CHART_MIRROR) != 0 };
    ad_point(&pp, pcen[0], pcen[1], gp);
    ad_point(pose, centroid[0], centroid[1], gc);
    /* back into the chart's frame: undo the rotation (the mirror only swaps left and right) */
    double dx = gp[0] - gc[0], dy = gp[1] - gc[1], ct = cos(-pose->theta), sn = sin(-pose->theta);
    double lx = ct * dx - sn * dy, ly = sn * dx + ct * dy;
    if (pose->mirror) lx = -lx;
    double ang = atan2(ly, lx);
    int q = (int)floor((ang + AD_PI) / (AD_PI / 2.0));
    return ((q % 4) + 4) % 4;
}
static int ad_fit_order(const void *a, const void *b)
{
    const AdFit *x = a, *y = b;
    if (x->soft != y->soft) return x->soft < y->soft ? -1 : 1;
    if (x->corr != y->corr) return x->corr > y->corr ? -1 : 1;
    return x->rel < y->rel ? -1 : x->rel > y->rel;
}
static int ad_double_order(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return x < y ? -1 : x > y;
}

typedef struct AdRow {
    int32_t chart; int examined, edge, admitted, round;
    size_t partners, agreeing; int sides, required;
    double disc, overlap, sep_p50, rim_ok;
    AdPose pose; const char *reason;
    size_t selected_count;
    int32_t selected_rel[AD_MAX_FITS];
} AdRow;

/* Fill fits[] with the poses the chart's registered partners imply. */
static size_t ad_collect(const AsmRun *run, int32_t chart, const size_t *adj_off, const int32_t *adj_rel, AdFit *fits)
{
    size_t n = 0;
    for (size_t i = adj_off[chart]; i < adj_off[chart + 1]; i++) {
        const AsmRelation *r = run->rels + adj_rel[i];
        int32_t other = r->a == chart ? r->b : r->a;
        if (other < 0 || (size_t)other >= run->n_charts || other == chart) continue;
        if (!ad_usable(r) || !AsmChart_registered(run->charts + other)) continue;
        AdFit f; memset(&f, 0, sizeof f);
        f.partner = ad_implied(run, r, chart, &f.pose); f.rel = adj_rel[i]; f.soft = ad_soft(r); f.rms = r->rms; f.corr = r->n_corr;
        if (!isfinite(f.pose.theta) || !isfinite(f.pose.x) || !isfinite(f.pose.y)) continue;
        fits[n++] = f;
        if (n == AD_MAX_FITS * 4) break;
    }
    qsort(fits, n, sizeof(AdFit), ad_fit_order);
    return n > AD_MAX_FITS ? AD_MAX_FITS : n;
}

static void ad_examine(const AsmRun *run, const AdIndex *ix, int32_t chart, const size_t *adj_off, const int32_t *adj_rel,
                       int required_interior, int soft_min_sides, int edge, AdRow *row)
{
    const AsmChart *ch = run->charts + chart;
    AdFit fits[AD_MAX_FITS * 4];
    size_t idx[AD_SAMPLES];
    size_t n = ad_collect(run, chart, adj_off, adj_rel, fits), ns = ad_samples(ch, idx, AD_SAMPLES);
    double centroid[2]; ad_centroid(ch, centroid);
    row->chart = chart; row->examined = 1; row->edge = edge; row->admitted = 0; row->partners = n; row->agreeing = 0;
    row->sides = 0; row->required = edge ? 1 : required_interior; row->disc = 0.0; row->overlap = 0.0; row->sep_p50 = 0.0; row->rim_ok = 0.0;
    row->selected_count = 0;
    memset(&row->pose, 0, sizeof row->pose);
    if (!n) { row->reason = "no_registered_partner"; return; }
    /* agreement: the largest set of poses that all agree pairwise */
    unsigned char agree[AD_MAX_FITS][AD_MAX_FITS];
    for (size_t i = 0; i < n; i++) for (size_t k = 0; k < n; k++)
        agree[i][k] = i == k || ad_displacement(ch, &fits[i].pose, &fits[k].pose, idx, ns) <= ASM_DISCOVER_AGREE_VOX;
    unsigned char best[AD_MAX_FITS], cur[AD_MAX_FITS]; size_t best_n = 0;
    memset(best, 0, sizeof best);
    for (size_t i = 0; i < n; i++) {
        size_t cur_n = 1; memset(cur, 0, sizeof cur); cur[i] = 1;
        for (size_t k = 0; k < n; k++) {
            if (k == i) continue;
            int ok = 1;
            for (size_t m = 0; m < n && ok; m++) if (cur[m] && !agree[k][m]) ok = 0;
            if (ok) { cur[k] = 1; cur_n++; }
        }
        if (cur_n > best_n) { best_n = cur_n; memcpy(best, cur, sizeof best); }
    }
    row->agreeing = best_n;
    for (size_t i = 0; i < n; i++) if (best[i]) row->selected_rel[row->selected_count++] = fits[i].rel;
    int hard = 0, dissent = 0, soft_dissent = 0;
    for (size_t i = 0; i < n; i++) {
        if (best[i] && !fits[i].soft) hard = 1;
        if (!best[i]) { if (fits[i].soft) soft_dissent = 1; else dissent = 1; }
    }
    /* the mean pose of the agreeing set: mean rotation, translation fixed by the mean image of the centroid */
    double cs = 0.0, sn = 0.0, gx = 0.0, gy = 0.0; size_t support = 0; int32_t parent = -1; size_t parent_corr = 0;
    for (size_t i = 0; i < n; i++) if (best[i]) {
        double g[2]; ad_point(&fits[i].pose, centroid[0], centroid[1], g);
        cs += cos(fits[i].pose.theta); sn += sin(fits[i].pose.theta); gx += g[0]; gy += g[1]; support += fits[i].corr;
        if (parent < 0 || fits[i].corr > parent_corr) { parent = fits[i].partner; parent_corr = fits[i].corr; }
    }
    AdPose pose; pose.theta = atan2(sn, cs); pose.mirror = -1;
    for (size_t i = 0; i < n; i++) if (best[i]) { pose.mirror = fits[i].pose.mirror; break; }
    {
        double u = pose.mirror ? -centroid[0] : centroid[0], v = centroid[1], ct = cos(pose.theta), st = sin(pose.theta);
        pose.x = gx / (double)best_n - (ct * u - st * v);
        pose.y = gy / (double)best_n - (st * u + ct * v);
    }
    row->pose = pose;
    double disc = 0.0; unsigned char quads[4] = {0, 0, 0, 0};
    for (size_t i = 0; i < n; i++) if (best[i]) {
        double d = ad_displacement(ch, &fits[i].pose, &pose, idx, ns); if (d > disc) disc = d;
        quads[ad_quadrant(run, centroid, &fits[i].pose, fits[i].partner)] = 1;
    }
    row->disc = disc; row->sides = quads[0] + quads[1] + quads[2] + quads[3];
    /* the placed material where the chart would land: how much, and is it the same surface */
    size_t hits = 0, nsep = 0;
    double *sep = malloc((ch->nv ? ch->nv : 1) * sizeof(double));
    if (!sep) { row->reason = "out_of_memory"; return; }
    for (size_t k = 0; k < ch->nv; k++) {
        double g[2]; ad_point(&pose, ch->uv[2 * k], ch->uv[2 * k + 1], g);
        AdNear nb; ad_near(ix, run, g[0], g[1], ch->xyz + 3 * k, ASM_DISCOVER_OVERLAP_VOX, &nb);
        if (nb.hits) { hits++; sep[nsep++] = nb.d3_min; }
    }
    row->overlap = ch->nv ? (double)hits / (double)ch->nv : 0.0;
    if (nsep) { qsort(sep, nsep, sizeof(double), ad_double_order); row->sep_p50 = sep[nsep / 2]; }
    free(sep);
    /* the rim: boundary vertices whose nearest placed vertex in uv is also nb in 3-D */
    size_t rim = 0, rim_ok = 0;
    for (size_t k = 0; k < ch->nv; k++) {
        if (!ch->boundary || !ch->boundary[k]) continue;
        double g[2]; ad_point(&pose, ch->uv[2 * k], ch->uv[2 * k + 1], g);
        AdNear nb; ad_near(ix, run, g[0], g[1], ch->xyz + 3 * k, ASM_DISCOVER_RIM_VOX, &nb);
        if (!nb.hits) continue;
        const float *q = run->charts[nb.chart].xyz + 3 * (size_t)nb.vert, *p = ch->xyz + 3 * k;
        double d3 = sqrt((q[0] - p[0]) * (q[0] - p[0]) + (q[1] - p[1]) * (q[1] - p[1]) + (q[2] - p[2]) * (q[2] - p[2]));
        rim++; rim_ok += d3 <= ASM_DISCOVER_RIM_VOX;
    }
    row->rim_ok = rim ? (double)rim_ok / (double)rim : 1.0;
    if (dissent) { row->reason = "sides_disagree"; return; }
    /* Evidence from gate-rejected seams alone (switched off, weak, short) counts
     * only when it is unanimous and comes from two sides: two independent
     * rejected seams agreeing on one pose is the evidence, one is not. */
    if (!hard && soft_dissent) { row->reason = "soft_evidence_disagrees"; return; }
    if (!hard && row->sides < soft_min_sides) { row->reason = "soft_evidence_only"; return; }
    if (row->sides < row->required) { row->reason = "too_few_sides"; return; }
    if (row->overlap > ASM_DISCOVER_MAX_OVERLAP && row->sep_p50 > ASM_DISCOVER_SAME_SURFACE_VOX) { row->reason = "another_wrap"; return; }
    if (row->rim_ok < ASM_DISCOVER_MIN_RIM) { row->reason = "rim_inconsistent"; return; }
    row->reason = "admitted"; row->admitted = 1;
    (void)support; (void)parent;
}

static int ad_selected(const AdRow *row, int32_t relation)
{
    for (size_t i = 0; i < row->selected_count; i++) if (row->selected_rel[i] == relation) return 1;
    return 0;
}

/* Place the chart at the row's pose: the same bookkeeping as a seam placement,
 * plus the source obligations the repair will close. */
static void ad_admit(AsmRun *run, const AdRow *row, const size_t *adj_off, const int32_t *adj_rel, int selected_obligations)
{
    AsmChart *ch = run->charts + row->chart;
    int32_t parent = -1; size_t parent_corr = 0, support = 0;
    for (size_t i = adj_off[row->chart]; i < adj_off[row->chart + 1]; i++) {
        if (selected_obligations && !ad_selected(row, adj_rel[i])) continue;
        const AsmRelation *r = run->rels + adj_rel[i];
        int32_t other = r->a == row->chart ? r->b : r->a;
        if (!ad_usable(r) || other < 0 || (size_t)other >= run->n_charts || !AsmChart_registered(run->charts + other)) continue;
        support += r->n_corr;
        if (parent < 0 || r->n_corr > parent_corr) { parent = other; parent_corr = r->n_corr; }
    }
    ch->pose_theta = row->pose.theta; ch->pose_x = row->pose.x; ch->pose_y = row->pose.y;
    ch->flags = (ch->flags & ~(uint32_t)ASM_CHART_MIRROR) | (row->pose.mirror ? ASM_CHART_MIRROR : 0u);
    ch->placed = 1; ch->placement_state = ASM_PLACE_SEAM; ch->placed_uv = NULL;
    ch->placement_parent = parent; ch->placement_support = support; ch->placement_rms = row->disc;
    if (parent >= 0) ch->component = run->charts[parent].component;
    for (size_t i = adj_off[row->chart]; i < adj_off[row->chart + 1]; i++) {
        if (selected_obligations && !ad_selected(row, adj_rel[i])) continue;
        AsmRelation *r = run->rels + adj_rel[i];
        int32_t other = r->a == row->chart ? r->b : r->a;
        if (other >= 0 && (size_t)other < run->n_charts && AsmChart_registered(run->charts + other) && ad_source(r)) r->continuity |= ASM_CONT_REQUIRED;
    }
}

static int ad_is_edge(const AsmRun *run, const AsmChart *ch, const double lo[3], const double hi[3], int have_extent)
{
    (void)run;
    if (!have_extent) return 0;
    for (int d = 0; d < 3; d++) {
        if (ch->bbox_lo[d] <= lo[d] + ASM_DISCOVER_EDGE_MARGIN_VOX) return 1;
        if (ch->bbox_hi[d] >= hi[d] - ASM_DISCOVER_EDGE_MARGIN_VOX) return 1;
    }
    return 0;
}

void AsmDiscover_default_opts(AsmDiscoverOpts *o)
{
    o->enabled = ASM_DISCOVER;
    o->min_sides = ASM_DISCOVER_MIN_SIDES;
    o->soft_min_sides = ASM_DISCOVER_SOFT_MIN_SIDES;
    o->rounds = ASM_DISCOVER_ROUNDS;
    o->selected_obligations = 0;
}

int AsmDiscover_run(AsmRun *run, const AsmDiscoverOpts *o, const char *out_dir, AsmDiscoverStats *st)
{
    memset(st, 0, sizeof *st);
    if (!run || !o || o->min_sides < 1 || o->min_sides > 4 || o->soft_min_sides < 1 || o->soft_min_sides > 4 || o->rounds < 1 || o->rounds > 64 ||
        o->selected_obligations < 0 || o->selected_obligations > 1) return -1;
    double start = ves_clock_sec();
    Arena_T arena = Arena_new();
    size_t nc = run->n_charts;
    /* Stage-entry requirements make newly selected evidence reviewable. */
    uint8_t *required_before = NULL;
    if (out_dir) {
        required_before = ARENA_ALLOC(arena, run->n_rels ? run->n_rels : 1);
        for (size_t i = 0; i < run->n_rels; i++) required_before[i] = (run->rels[i].continuity & ASM_CONT_REQUIRED) != 0;
    }
    /* relation adjacency per chart */
    size_t *adj_off = ARENA_CALLOC(arena, nc + 2, sizeof(size_t));
    for (size_t i = 0; i < run->n_rels; i++) {
        const AsmRelation *r = run->rels + i;
        if (r->a < 0 || r->b < 0 || (size_t)r->a >= nc || (size_t)r->b >= nc) continue;
        adj_off[r->a + 1]++; adj_off[r->b + 1]++;
    }
    for (size_t c = 1; c <= nc; c++) adj_off[c] += adj_off[c - 1];
    int32_t *adj_rel = ARENA_ALLOC(arena, (adj_off[nc] ? adj_off[nc] : 1) * sizeof(int32_t));
    size_t *fill = ARENA_ALLOC(arena, (nc ? nc : 1) * sizeof(size_t));
    memcpy(fill, adj_off, nc * sizeof(size_t));
    for (size_t i = 0; i < run->n_rels; i++) {
        const AsmRelation *r = run->rels + i;
        if (r->a < 0 || r->b < 0 || (size_t)r->a >= nc || (size_t)r->b >= nc) continue;
        adj_rel[fill[r->a]++] = (int32_t)i; adj_rel[fill[r->b]++] = (int32_t)i;
    }
    /* charts replaced by cut children never re-enter */
    uint8_t *replaced = ARENA_CALLOC(arena, nc ? nc : 1, 1);
    for (size_t c = 0; c < nc; c++) if (run->charts[c].parent >= 0 && (size_t)run->charts[c].parent < nc) replaced[run->charts[c].parent] = 1;
    /* the pile's extent, for the crop-face rule */
    double lo[3] = {DBL_MAX, DBL_MAX, DBL_MAX}, hi[3] = {-DBL_MAX, -DBL_MAX, -DBL_MAX};
    int have_extent = 0;
    if (run->pile && run->n_cubes) {
        double edge = run->cube_size > 0 ? run->cube_size : 128.0;
        for (size_t i = 0; i < run->n_cubes; i++) {
            const MeshPileEntry *e = run->pile + i;
            if (!e->has_id) continue;
            double org[3] = {(double)e->oz, (double)e->oy, (double)e->ox};
            for (int d = 0; d < 3; d++) { if (org[d] < lo[d]) lo[d] = org[d]; if (org[d] + edge > hi[d]) hi[d] = org[d] + edge; }
            have_extent = 1;
        }
    }
    AdIndex ix; memset(&ix, 0, sizeof ix);
    if (ad_index_build(arena, run, &ix)) { Arena_dispose(&arena); return -1; }
    /* control: the placed charts' partner-implied poses against their actual poses */
    {
        double *disp = ARENA_ALLOC(arena, (run->n_rels ? 2 * run->n_rels : 1) * sizeof(double));
        size_t nd = 0, beyond = 0, idx[AD_SAMPLES];
        for (size_t c = 0; c < nc; c++) {
            const AsmChart *ch = run->charts + c;
            if (!AsmChart_registered(ch) || !ch->uv || !ch->nv) continue;
            AdPose actual = { ch->pose_theta, ch->pose_x, ch->pose_y, (ch->flags & ASM_CHART_MIRROR) != 0 };
            size_t ns = ad_samples(ch, idx, AD_SAMPLES);
            for (size_t i = adj_off[c]; i < adj_off[c + 1]; i++) {
                const AsmRelation *r = run->rels + adj_rel[i];
                int32_t other = r->a == (int32_t)c ? r->b : r->a;
                if (!ad_usable(r) || ad_soft(r) || !AsmChart_registered(run->charts + other)) continue;
                if (!AsmRel_is_join(r) || run->charts[other].component != ch->component) continue;
                AdPose implied; ad_implied(run, r, (int32_t)c, &implied);
                if (implied.mirror != actual.mirror) { disp[nd++] = DBL_MAX; beyond++; continue; }
                double d = ad_displacement(ch, &implied, &actual, idx, ns);
                disp[nd++] = d; beyond += d > ASM_DISCOVER_AGREE_VOX;
            }
        }
        if (nd) {
            qsort(disp, nd, sizeof(double), ad_double_order);
            st->control_p50 = disp[nd / 2]; st->control_p90 = disp[(nd * 9) / 10]; st->control_max = disp[nd - 1];
        }
        st->control_poses = nd; st->control_beyond = beyond;
    }
    AdRow *rows = ARENA_CALLOC(arena, nc ? nc : 1, sizeof(AdRow));
    for (size_t c = 0; c < nc; c++) { rows[c].chart = (int32_t)c; rows[c].reason = "not_examined"; }
    int round = 0;
    for (round = 1; round <= o->rounds; round++) {
        size_t admitted_now = 0;
        for (size_t c = 0; c < nc; c++) {
            const AsmChart *ch = run->charts + c;
            if (!AsmChart_in_layout(ch) || AsmChart_registered(ch) || replaced[c] || !ch->xyz || !ch->nv) continue;
            AdRow row; memset(&row, 0, sizeof row);
            ad_examine(run, &ix, (int32_t)c, adj_off, adj_rel, o->min_sides, o->soft_min_sides, ad_is_edge(run, ch, lo, hi, have_extent), &row);
            row.round = round;
            if (!rows[c].examined) st->examined++;
            if (row.partners && (!rows[c].examined || !rows[c].partners)) st->with_partner++;
            rows[c] = row;
        }
        for (size_t c = 0; c < nc; c++) if (rows[c].examined && rows[c].admitted && rows[c].round == round && !AsmChart_registered(run->charts + c)) {
            if (!o->enabled) continue;
            ad_admit(run, rows + c, adj_off, adj_rel, o->selected_obligations);
            st->admitted++; st->admitted_area += run->charts[c].area3d; admitted_now++;
            if (rows[c].overlap > ASM_DISCOVER_MAX_OVERLAP) st->admitted_overlapping++;
        }
        st->rounds = round;
        if (!admitted_now) break;
        if (ad_index_build(arena, run, &ix)) { Arena_dispose(&arena); return -1; }
    }
    for (size_t c = 0; c < nc; c++) if (rows[c].examined && !rows[c].admitted) {
        const char *why = rows[c].reason;
        if (!strcmp(why, "sides_disagree") || !strcmp(why, "soft_evidence_disagrees")) st->refused_disagree++;
        else if (!strcmp(why, "too_few_sides") || !strcmp(why, "soft_evidence_only")) st->refused_sides++;
        else if (!strcmp(why, "rim_inconsistent")) st->refused_rim++;
        else if (!strcmp(why, "another_wrap")) st->refused_wrap++;
    }
    int rc = 0;
    if (out_dir) {
        char path[2048];
        snprintf(path, sizeof path, "%s/stage5_discovery.csv", out_dir);
        FILE *fp = fopen(path, "wb");
        if (!fp) rc = -1;
        else {
            int ok = fprintf(fp, "chart,area,edge,partners,agreeing,sides,required,discrepancy,overlap,separation_p50,rim_ok,u,v,theta,mirror,round,admitted,reason\n") >= 0;
            for (size_t c = 0; c < nc && ok; c++) {
                const AdRow *r = rows + c;
                if (!r->examined) continue;
                double cen[2] = {0, 0}, g[2] = {0, 0};
                ad_centroid(run->charts + c, cen); ad_point(&r->pose, cen[0], cen[1], g);
                ok = fprintf(fp, "%d,%.6e,%d,%zu,%zu,%d,%d,%.3f,%.4f,%.3f,%.4f,%.2f,%.2f,%.6f,%d,%d,%d,%s\n",
                             r->chart, run->charts[c].area3d, r->edge, r->partners, r->agreeing, r->sides, r->required, r->disc, r->overlap, r->sep_p50, r->rim_ok,
                             r->partners ? g[0] : 0.0, r->partners ? g[1] : 0.0, r->partners ? r->pose.theta : 0.0, r->partners ? r->pose.mirror : 0, r->round, r->admitted, r->reason) >= 0;
            }
            if (fclose(fp) || !ok) rc = -1;
        }
        if (!rc) {
            snprintf(path, sizeof path, "%s/stage5_discovery_evidence.csv", out_dir);
            fp = fopen(path, "wb");
            if (!fp) rc = -1;
            else {
                int ok = fprintf(fp, "chart,round,relation,partner,selected,usable,flags,source,required_before,required_after,selected_obligations\n") >= 0;
                for (size_t c = 0; c < nc && ok; c++) if (rows[c].examined && rows[c].admitted) {
                    for (size_t i = adj_off[c]; i < adj_off[c + 1] && ok; i++) {
                        int32_t ri = adj_rel[i]; const AsmRelation *r = run->rels + ri;
                        int32_t other = r->a == (int32_t)c ? r->b : r->a;
                        ok = fprintf(fp, "%zu,%d,%d,%d,%d,%d,%u,%d,%d,%d,%d\n", c, rows[c].round, ri, other,
                            ad_selected(rows + c, ri), ad_usable(r), r->flags, ad_source(r) != 0, required_before[ri],
                            (r->continuity & ASM_CONT_REQUIRED) != 0, o->selected_obligations) >= 0;
                    }
                }
                if (fclose(fp) || !ok) rc = -1;
            }
        }
        if (!rc) {
            snprintf(path, sizeof path, "%s/stage5_discovered.png", out_dir);
            AsmChart *display = ARENA_ALLOC(arena, (nc ? nc : 1) * sizeof(AsmChart));
            memcpy(display, run->charts, nc * sizeof(AsmChart));
            for (size_t c = 0; c < nc; c++) { display[c].placed = AsmChart_registered(run->charts + c); display[c].component = 0; }
            if (nc && AsmReport_layout_png(arena, path, display, nc, NULL, 0, 8.0, 8192, 1)) rc = -1;
        }
    }
    st->sec = ves_clock_sec() - start;
    Arena_dispose(&arena);
    return rc;
}

/* ---- selftest ----------------------------------------------------------------- */

enum { ADT_N = 9, ADT_NV = ADT_N * ADT_N, ADT_NF = 2 * (ADT_N - 1) * (ADT_N - 1) };
typedef struct AdtChart { float xyz[3 * ADT_NV], uv[2 * ADT_NV]; int32_t faces[3 * ADT_NF]; uint8_t boundary[ADT_NV]; } AdtChart;

/* a square chart of ADT_N x ADT_N vertices: uv (x, y); xyz (z0, y0 + y, x0 + x) */
static void adt_chart(AdtChart *g, AsmChart *c, int32_t id, double z0, double y0, double x0)
{
    memset(c, 0, sizeof *c);
    for (int y = 0; y < ADT_N; y++) for (int x = 0; x < ADT_N; x++) {
        int v = y * ADT_N + x;
        g->xyz[3 * v] = (float)z0; g->xyz[3 * v + 1] = (float)(y0 + y); g->xyz[3 * v + 2] = (float)(x0 + x);
        g->uv[2 * v] = (float)x; g->uv[2 * v + 1] = (float)y;
        g->boundary[v] = (uint8_t)(x == 0 || y == 0 || x == ADT_N - 1 || y == ADT_N - 1);
        if (x + 1 < ADT_N && y + 1 < ADT_N) {
            int32_t f[6] = {v, v + 1, v + ADT_N + 1, v, v + ADT_N + 1, v + ADT_N};
            memcpy(g->faces + 6 * (y * (ADT_N - 1) + x), f, sizeof f);
        }
    }
    c->id = id; c->cube = 0; c->comp = id; c->parent = -1; c->component = id; c->nv = ADT_NV; c->nf = ADT_NF;
    c->xyz = g->xyz; c->uv = g->uv; c->faces = g->faces; c->boundary = g->boundary;
    c->area3d = (double)((ADT_N - 1) * (ADT_N - 1)); c->area_uv = c->area3d;
    c->bbox_lo[0] = (float)z0; c->bbox_lo[1] = (float)y0; c->bbox_lo[2] = (float)x0;
    c->bbox_hi[0] = (float)z0; c->bbox_hi[1] = (float)(y0 + ADT_N - 1); c->bbox_hi[2] = (float)(x0 + ADT_N - 1);
    c->centroid[0] = z0; c->centroid[1] = y0 + (ADT_N - 1) / 2.0; c->centroid[2] = x0 + (ADT_N - 1) / 2.0;
}
static void adt_register(AsmChart *c, double theta, double x, double y, int mirror)
{
    c->placed = 1; c->placement_state = ASM_PLACE_ROOT; c->pose_theta = theta; c->pose_x = x; c->pose_y = y;
    c->flags = mirror ? ASM_CHART_MIRROR : 0u;
}
/* the relation a-b consistent with the two charts' TRUE global poses:
 * m_a(p_a) = R(theta_b - theta_a) m_b(p_b) + R(-theta_a)(t_b - t_a), stored in unmirrored frames */
static void adt_relation(AsmRun *run, int32_t a, const AdPose *pa, int32_t b, const AdPose *pb, unsigned flags, unsigned continuity)
{
    AsmRelation r; memset(&r, 0, sizeof r);
    double th = pb->theta - pa->theta, dx = pb->x - pa->x, dy = pb->y - pa->y, ct = cos(-pa->theta), sn = sin(-pa->theta);
    double tx = ct * dx - sn * dy, ty = sn * dx + ct * dy;
    r.a = a; r.b = b; r.flags = flags | ((pa->mirror ^ pb->mirror) ? ASM_REL_PARITY : 0u);
    if (pa->mirror) { th = -th; tx = -tx; }
    r.theta = th; r.tx = tx; r.ty = ty; r.rms = 0.5; r.n_corr = 9; r.seam_len = 8.0; r.continuity = continuity; r.robust_w = 1.0;
    AsmRun_push_rel(run, &r);
}

/* Two valid anchors beat a higher-count weak alternative. Rejected, vetoed,
 * and not-yet-present partners must not become obligations in selected mode.
 * Existing requirements are immutable, including a disagreeing weak one. */
static int adt_obligations(int selected, int prior_required, int enabled)
{
    static AdtChart geo[7]; AdtChart saved[7];
    AsmRun run = {0}; run.arena = Arena_new(); AsmChart ch;
    AdPose truth[7] = {{0,100,100,0},{0,113,113,0},{0,100,113,0},
        {0,1000,1000,0},{0,1500,1500,0},{0,2000,2000,0},{0,600,600,0}};
    AdPose alternative = {0,3500,3500,0};
    for (int32_t c = 0; c < 7; c++) {
        adt_chart(geo + c, &ch, c, 0, truth[c].y - 100, truth[c].x - 100);
        if (c != 2 && c != 6) adt_register(&ch, 0, truth[c].x, truth[c].y, 0);
        AsmRun_push_chart(&run, &ch);
    }
    adt_relation(&run, 0, truth + 0, 2, truth + 2, 0, ASM_CONT_SOURCE);
    adt_relation(&run, 2, truth + 2, 1, truth + 1, 0, ASM_CONT_SOURCE | ASM_CONT_CUT);
    adt_relation(&run, 3, truth + 3, 2, &alternative, ASM_REL_WEAK,
        ASM_CONT_SOURCE | (prior_required ? ASM_CONT_REQUIRED : 0u));
    run.rels[2].n_corr = 99;
    adt_relation(&run, 4, truth + 4, 2, &alternative, ASM_REL_DROPPED, ASM_CONT_SOURCE);
    adt_relation(&run, 5, truth + 5, 2, &alternative, ASM_REL_CROSSWRAP | ASM_REL_WEAK, ASM_CONT_SOURCE);
    adt_relation(&run, 3, truth + 3, 2, &alternative, ASM_REL_CONTACT, ASM_CONT_SOURCE);
    adt_relation(&run, 0, truth + 0, 6, truth + 6, 0, ASM_CONT_SOURCE);
    adt_relation(&run, 2, truth + 2, 6, truth + 6, 0, ASM_CONT_SOURCE);
    AsmRelation before[8]; memcpy(before, run.rels, sizeof before); memcpy(saved, geo, sizeof saved);
    AsmDiscoverOpts opts; AsmDiscover_default_opts(&opts);
    opts.selected_obligations = selected; opts.enabled = enabled; opts.min_sides = 1; opts.rounds = 1;
    AsmDiscoverStats stats; int failed = AsmDiscover_run(&run, &opts, NULL, &stats) != 0;
    failed |= stats.admitted != (enabled ? 2u : 0u);
    for (int32_t c = 2; c <= 6; c += 4) {
        const AsmChart *got = run.charts + c;
        failed |= AsmChart_registered(got) != enabled;
        if (enabled) failed |= fabs(got->pose_x - truth[c].x) > 1e-9 || fabs(got->pose_y - truth[c].y) > 1e-9 ||
            fabs(got->pose_theta) > 1e-9 || (got->flags & ASM_CHART_MIRROR) != 0;
    }
    if (enabled) {
        const AsmChart *got = run.charts + 2;
        failed |= got->placement_parent != (selected ? 0 : 3);
        failed |= selected ? got->placement_support != 18 : got->placement_support <= 18;
        failed |= selected && run.charts[6].placement_support != 9;
    }
    for (size_t i = 0; i < 8; i++) {
        int chosen = selected ? (i == 0 || i == 1 || i == 6) : i != 5;
        AsmRelation want; memcpy(&want, before + i, sizeof want);
        if (enabled && chosen) want.continuity |= ASM_CONT_REQUIRED;
        failed |= memcmp(&want, run.rels + i, sizeof want) != 0;
    }
    failed |= memcmp(saved, geo, sizeof saved) != 0;
    Arena_dispose(&run.arena);
    if (failed) fprintf(stderr, "  asm_discover selftest FAIL: selected=%d prior_required=%d enabled=%d\n", selected, prior_required, enabled);
    return failed;
}

int AsmDiscover_selftest(void)
{
    int fails = 0;
    static AdtChart geo[8];
    AsmRun run; memset(&run, 0, sizeof run); run.arena = Arena_new();
    AsmChart c; AdPose truth[8];
    /* A at the origin and B up-right of it, both placed; L above A and left of B across 5-vox trim gaps; M mirrored, far right */
    truth[0] = (AdPose){0.0, 100.0, 100.0, 0}; truth[1] = (AdPose){0.0, 113.0, 113.0, 0}; truth[2] = (AdPose){0.0, 100.0, 113.0, 0};
    truth[3] = (AdPose){0.3, 300.0, 50.0, 0}; truth[4] = (AdPose){0.1, 330.0, 60.0, 1};
    adt_chart(geo + 0, &c, 0, 0.0, 0.0, 0.0); adt_register(&c, truth[0].theta, truth[0].x, truth[0].y, 0); AsmRun_push_chart(&run, &c);
    adt_chart(geo + 1, &c, 1, 0.0, 13.0, 13.0); adt_register(&c, truth[1].theta, truth[1].x, truth[1].y, 0); AsmRun_push_chart(&run, &c);
    adt_chart(geo + 2, &c, 2, 0.0, 13.0, 0.0); AsmRun_push_chart(&run, &c);                       /* L: unplaced, same surface, 5-vox trim gaps */
    adt_chart(geo + 3, &c, 3, 0.0, 0.0, 200.0); adt_register(&c, truth[3].theta, truth[3].x, truth[3].y, 0); AsmRun_push_chart(&run, &c);
    adt_chart(geo + 4, &c, 4, 0.0, 10.0, 200.0); AsmRun_push_chart(&run, &c);                     /* M: unplaced, mirrored parity */
    adt_chart(geo + 5, &c, 5, 15.0, 0.0, 0.0); AsmRun_push_chart(&run, &c);                       /* W: another wrap over A (15 vox up) */
    adt_chart(geo + 6, &c, 6, 2.0, 0.0, 0.0); AsmRun_push_chart(&run, &c);                        /* D: a duplicate surface over A (2 vox up) */
    adt_chart(geo + 7, &c, 7, 0.0, -10.0, 0.0); AsmRun_push_chart(&run, &c);                      /* X: below A, partners disagree */
    for (size_t i = 0; i < run.n_charts; i++) run.charts[i].component = (int32_t)(i < 2 ? 0 : (i == 3 ? 3 : (int32_t)i));
    adt_relation(&run, 0, truth + 0, 1, truth + 1, 0u, ASM_CONT_SOURCE);                           /* A-B: the control */
    adt_relation(&run, 0, truth + 0, 2, truth + 2, 0u, ASM_CONT_SOURCE);                           /* A-L and L-B: two sides */
    adt_relation(&run, 2, truth + 2, 1, truth + 1, 0u, ASM_CONT_SOURCE);
    adt_relation(&run, 3, truth + 3, 4, truth + 4, 0u, ASM_CONT_SOURCE);                           /* one side, mirrored: M is an interior chart */
    truth[5] = (AdPose){0.0, 100.0, 100.0, 0}; truth[6] = truth[5];                                /* W and D would land exactly on A */
    adt_relation(&run, 1, truth + 1, 5, truth + 5, 0u, ASM_CONT_SOURCE);
    adt_relation(&run, 5, truth + 5, 3, truth + 3, 0u, ASM_CONT_SOURCE);
    adt_relation(&run, 1, truth + 1, 6, truth + 6, 0u, ASM_CONT_SOURCE);
    adt_relation(&run, 6, truth + 6, 3, truth + 3, 0u, ASM_CONT_SOURCE);
    truth[7] = (AdPose){0.0, 100.0, 90.0, 0};
    adt_relation(&run, 0, truth + 0, 7, truth + 7, 0u, ASM_CONT_SOURCE);                           /* X from A: below A ... */
    { AdPose other = {0.0, 130.0, 90.0, 0}; adt_relation(&run, 7, &other, 1, truth + 1, 0u, ASM_CONT_SOURCE); } /* ... and 30 vox away from B */
    AsmDiscoverOpts o; AsmDiscover_default_opts(&o); o.min_sides = 2;
    AsmDiscoverStats st;
    if (AsmDiscover_run(&run, &o, NULL, &st)) { fprintf(stderr, "  asm_discover selftest FAIL: run\n"); fails++; }
    if (!(st.control_poses == 2 && st.control_max < 1e-6)) { fprintf(stderr, "  asm_discover selftest FAIL: control %zu poses, max %.3g\n", st.control_poses, st.control_max); fails++; }
    const AsmChart *L = run.charts + 2;
    if (!AsmChart_registered(L) || fabs(L->pose_x - truth[2].x) > 1e-6 || fabs(L->pose_y - truth[2].y) > 1e-6 || fabs(L->pose_theta) > 1e-9 || L->component != 0)
        { fprintf(stderr, "  asm_discover selftest FAIL: two agreeing sides not placed (registered %d pose %.3f %.3f %.3f comp %d)\n", AsmChart_registered(L), L->pose_x, L->pose_y, L->pose_theta, L->component); fails++; }
    if (AsmChart_registered(run.charts + 4)) { fprintf(stderr, "  asm_discover selftest FAIL: one interior side placed\n"); fails++; }
    if (AsmChart_registered(run.charts + 5)) { fprintf(stderr, "  asm_discover selftest FAIL: another wrap placed over A\n"); fails++; }
    if (!AsmChart_registered(run.charts + 6)) { fprintf(stderr, "  asm_discover selftest FAIL: same-surface duplicate refused\n"); fails++; }
    if (AsmChart_registered(run.charts + 7)) { fprintf(stderr, "  asm_discover selftest FAIL: disagreeing sides placed\n"); fails++; }
    if (!(st.admitted == 2 && st.refused_wrap == 1 && st.refused_disagree == 1 && st.refused_sides == 1 && st.admitted_overlapping == 1))
        { fprintf(stderr, "  asm_discover selftest FAIL: counts admitted %zu wrap %zu disagree %zu sides %zu overlapping %zu\n", st.admitted, st.refused_wrap, st.refused_disagree, st.refused_sides, st.admitted_overlapping); fails++; }
    /* the admitted chart's source obligations are required */
    { int req = 0; for (size_t i = 0; i < run.n_rels; i++) if ((run.rels[i].a == 2 || run.rels[i].b == 2) && (run.rels[i].continuity & ASM_CONT_REQUIRED)) req++;
      if (req != 2) { fprintf(stderr, "  asm_discover selftest FAIL: %d required obligations on the admitted chart\n", req); fails++; } }
    /* mirrored parity: with M an edge chart (one side suffices) it lands on its true, mirrored pose */
    o.min_sides = 1;
    if (AsmDiscover_run(&run, &o, NULL, &st)) { fprintf(stderr, "  asm_discover selftest FAIL: second run\n"); fails++; }
    const AsmChart *M = run.charts + 4;
    if (!AsmChart_registered(M) || !(M->flags & ASM_CHART_MIRROR) || fabs(M->pose_x - truth[4].x) > 1e-6 || fabs(M->pose_y - truth[4].y) > 1e-6 || fabs(M->pose_theta - truth[4].theta) > 1e-9)
        { fprintf(stderr, "  asm_discover selftest FAIL: mirrored chart (registered %d mirror %d pose %.3f %.3f %.4f)\n", AsmChart_registered(M), (M->flags & ASM_CHART_MIRROR) != 0, M->pose_x, M->pose_y, M->pose_theta); fails++; }
    /* report only: nothing moves */
    { AsmRun run2; memset(&run2, 0, sizeof run2); run2.arena = Arena_new();
      adt_chart(geo + 0, &c, 0, 0.0, 0.0, 0.0); adt_register(&c, 0.0, 100.0, 100.0, 0); AsmRun_push_chart(&run2, &c);
      adt_chart(geo + 2, &c, 1, 0.0, 13.0, 0.0); AsmRun_push_chart(&run2, &c);
      adt_relation(&run2, 0, truth + 0, 1, truth + 2, 0u, ASM_CONT_SOURCE);
      AsmDiscoverOpts off = o; off.enabled = 0;
      AsmDiscover_run(&run2, &off, NULL, &st);
      if (AsmChart_registered(run2.charts + 1) || st.admitted) { fprintf(stderr, "  asm_discover selftest FAIL: report-only placed a chart\n"); fails++; }
      Arena_dispose(&run2.arena); }
    Arena_dispose(&run.arena);
    for (int selected = 0; selected <= 1; selected++) for (int prior = 0; prior <= 1; prior++) for (int enabled = 0; enabled <= 1; enabled++)
        fails += adt_obligations(selected, prior, enabled);
    fprintf(stderr, "  asm_discover selftest: %s\n", fails ? "FAIL" : "ok");
    return fails;
}
