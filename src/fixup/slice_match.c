/* slice_match.c — per-plane endpoint matching. See slice_match.h. */

#include "slice_match.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/pipeline_constants.h"
#include "join_paint.h"

static inline size_t sm_idx(int y, int x, int W)
{
    return (size_t)y * (size_t)W + (size_t)x;
}

/* ---------------------------------------------------------------- */
/* Local corridor flood: adjacent-plane connectivity evidence       */
/* ---------------------------------------------------------------- */

/* In `plane` (binary H*W), can a's fragment reach b's fragment through
 * foreground restricted to the join corridor (endpoint bbox dilated by
 * FIXUP_ADJ_CORRIDOR)? Seeds are foreground within FIXUP_ADJ_SEED_R of a's
 * position; target likewise for b. 8-connected BFS on an arena-scratch
 * window. Returns 1 if connected. */
static int sm_corridor_connected(Arena_T arena, const uint8_t *plane,
                                 int H, int W,
                                 const SliceTraceEndpoint *a,
                                 const SliceTraceEndpoint *b)
{
    int pad = FIXUP_ADJ_CORRIDOR + FIXUP_ADJ_SEED_R;
    int y0 = (a->y < b->y ? a->y : b->y) - pad;
    int y1 = (a->y > b->y ? a->y : b->y) + pad;
    int x0 = (a->x < b->x ? a->x : b->x) - pad;
    int x1 = (a->x > b->x ? a->x : b->x) + pad;
    if (y0 < 0) y0 = 0;
    if (x0 < 0) x0 = 0;
    if (y1 >= H) y1 = H - 1;
    if (x1 >= W) x1 = W - 1;
    int wh = y1 - y0 + 1, ww = x1 - x0 + 1;
    if (wh <= 0 || ww <= 0) return 0;

    Arena_Mark mark = Arena_save(arena);
    size_t wsz = (size_t)wh * (size_t)ww;
    uint8_t *visited = ARENA_CALLOC(arena, wsz, 1);
    int32_t *queue = ARENA_ALLOC(arena, wsz * sizeof(int32_t));
    int32_t qhead = 0, qtail = 0;

    int r2 = FIXUP_ADJ_SEED_R * FIXUP_ADJ_SEED_R;

    /* Seed: foreground near a. */
    for (int dy = -FIXUP_ADJ_SEED_R; dy <= FIXUP_ADJ_SEED_R; dy++) {
        for (int dx = -FIXUP_ADJ_SEED_R; dx <= FIXUP_ADJ_SEED_R; dx++) {
            if (dy * dy + dx * dx > r2) continue;
            int py = a->y + dy, px = a->x + dx;
            if (py < y0 || py > y1 || px < x0 || px > x1) continue;
            if (!plane[sm_idx(py, px, W)]) continue;
            size_t widx = (size_t)(py - y0) * (size_t)ww + (size_t)(px - x0);
            if (visited[widx]) continue;
            visited[widx] = 1;
            queue[qtail++] = (int32_t)widx;
        }
    }

    static const int dy8[8] = {-1, -1, -1, 0, 0, 1, 1, 1};
    static const int dx8[8] = {-1,  0,  1, -1, 1, -1, 0, 1};

    int connected = 0;
    while (qhead < qtail && !connected) {
        int32_t widx = queue[qhead++];
        int wy = (int)(widx / ww), wx = (int)(widx % ww);
        int py = wy + y0, px = wx + x0;

        /* target reached? */
        int ty = py - b->y, tx = px - b->x;
        if (ty * ty + tx * tx <= r2) { connected = 1; break; }

        for (int d = 0; d < 8; d++) {
            int ny = wy + dy8[d], nx = wx + dx8[d];
            if (ny < 0 || ny >= wh || nx < 0 || nx >= ww) continue;
            size_t nidx = (size_t)ny * (size_t)ww + (size_t)nx;
            if (visited[nidx]) continue;
            if (!plane[sm_idx(ny + y0, nx + x0, W)]) continue;
            visited[nidx] = 1;
            queue[qtail++] = (int32_t)nidx;
        }
    }

    Arena_restore(arena, mark);
    return connected;
}

/* s_adj in [0,1]: weighted corridor connectivity over z-1,z+1 (weight 1)
 * and z-2,z+2 (weight 0.5), denominator fixed at 3 (missing planes at the
 * volume edge count as no evidence — conservative). */
static float sm_s_adj(Arena_T arena, const uint8_t *vol, int D, int H, int W,
                      int z, const SliceTraceEndpoint *a,
                      const SliceTraceEndpoint *b)
{
    static const int dz[4] = {-1, 1, -2, 2};
    static const float wz[4] = {1.0f, 1.0f, 0.5f, 0.5f};
    size_t plane_sz = (size_t)H * (size_t)W;
    float conn = 0.0f;
    for (int k = 0; k < 4; k++) {
        int zz = z + dz[k];
        if (zz < 0 || zz >= D) continue;
        if (sm_corridor_connected(arena, vol + (size_t)zz * plane_sz,
                                  H, W, a, b))
            conn += wz[k];
    }
    return conn / 3.0f;
}

/* ---------------------------------------------------------------- */
/* Candidate scoring + greedy stable matching                       */
/* ---------------------------------------------------------------- */

typedef struct {
    int32_t a, b;
    float score, dist, s_adj, s_track;
    uint8_t far_tier;
} sm_cand;

static int sm_cand_cmp(const void *pa, const void *pb)
{
    const sm_cand *ca = (const sm_cand *)pa;
    const sm_cand *cb = (const sm_cand *)pb;
    if (ca->score > cb->score) return -1;
    if (ca->score < cb->score) return 1;
    if (ca->a != cb->a) return (ca->a < cb->a) ? -1 : 1;
    if (ca->b != cb->b) return (ca->b < cb->b) ? -1 : 1;
    return 0;
}

enum { SM_MAX_REJECTS = 8192 };

static void sm_add_reject(SliceMatchReject *rej, int32_t *n_rej,
                          SliceMatchStats *stats,
                          int32_t a, int32_t b, uint8_t reason)
{
    if (stats) stats->rej_count[reason]++;
    if (rej && *n_rej < SM_MAX_REJECTS) {
        rej[*n_rej].a = a;
        rej[*n_rej].b = b;
        rej[*n_rej].reason = reason;
        (*n_rej)++;
    }
}

int SliceMatch_run(Arena_T arena,
                   const SliceTraceEndpoint *eps, int32_t n_eps,
                   const uint8_t *vol, const int32_t *all_labels,
                   int D, int H, int W, int z,
                   const int32_t *mask_labels,
                   const float *s_track,
                   const SliceMatchParams *params,
                   int32_t *match_out,
                   SliceMatchJoin **out_joins, int32_t *out_n_joins,
                   SliceMatchReject **out_rejects, int32_t *out_n_rejects,
                   SliceMatchStats *stats)
{
    assert(vol && all_labels && mask_labels && params && match_out);
    assert(out_joins && out_n_joins && out_rejects && out_n_rejects);

    *out_joins = NULL;
    *out_n_joins = 0;
    *out_rejects = NULL;
    *out_n_rejects = 0;
    if (stats) memset(stats, 0, sizeof(*stats));
    for (int32_t i = 0; i < n_eps; i++) match_out[i] = -1;
    if (n_eps < 2) return 0;

    SliceMatchReject *rej =
        ARENA_ALLOC(arena, (size_t)SM_MAX_REJECTS * sizeof(*rej));
    int32_t n_rej = 0;

    /* Pass 1: count pairs within reach (candidate array bound). */
    float reach2 = params->reach_max * params->reach_max;
    int32_t n_in_reach = 0;
    for (int32_t i = 0; i < n_eps; i++) {
        if (eps[i].excluded) continue;
        for (int32_t j = i + 1; j < n_eps; j++) {
            if (eps[j].excluded) continue;
            float dy = (float)(eps[j].y - eps[i].y);
            float dx = (float)(eps[j].x - eps[i].x);
            if (dy * dy + dx * dx <= reach2) n_in_reach++;
        }
    }
    if (stats) stats->n_pairs_in_reach = n_in_reach;
    if (n_in_reach == 0) {
        *out_rejects = rej;
        *out_n_rejects = 0;
        return 0;
    }

    sm_cand *cands = ARENA_ALLOC(arena, (size_t)n_in_reach * sizeof(*cands));
    int32_t n_cands = 0;

    /* Pass 2: gates + scoring. */
    for (int32_t i = 0; i < n_eps; i++) {
        if (eps[i].excluded) continue;
        for (int32_t j = i + 1; j < n_eps; j++) {
            if (eps[j].excluded) continue;
            float dy = (float)(eps[j].y - eps[i].y);
            float dx = (float)(eps[j].x - eps[i].x);
            float d2 = dy * dy + dx * dx;
            if (d2 > reach2) continue;
            float d = sqrtf(d2);
            if (d < 0.5f) continue;         /* coincident tips: nothing to join */

            if (eps[i].skel_cc == eps[j].skel_cc ||
                eps[i].mask_cc == eps[j].mask_cc) {
                sm_add_reject(rej, &n_rej, stats, i, j, SLICEMATCH_REJ_SAME_CC);
                continue;
            }

            /* facing + opposition gates, tiered: short hooks get 45 deg
             * (walk tangents are noisy at 4px gaps and the short chord +
             * radial + merger + support gates carry safety); far joins
             * must have clean geometry. */
            int is_far = (d > params->reach_safe);
            float need_cos = is_far ? FIXUP_FACING_MIN_COS
                                    : FIXUP_FACING_MIN_COS_NEAR;
            float need_opp = is_far ? FIXUP_OPPOSE_MAX_DOT
                                    : FIXUP_OPPOSE_MAX_DOT_NEAR;
            float dir_dy = dy / d, dir_dx = dx / d;
            float cos_a = eps[i].tan_dy * dir_dy + eps[i].tan_dx * dir_dx;
            float cos_b = -(eps[j].tan_dy * dir_dy + eps[j].tan_dx * dir_dx);
            float oppose = eps[i].tan_dy * eps[j].tan_dy +
                           eps[i].tan_dx * eps[j].tan_dx;
            if (cos_a < need_cos || cos_b < need_cos || oppose > need_opp) {
                sm_add_reject(rej, &n_rej, stats, i, j, SLICEMATCH_REJ_TANGENT);
                continue;
            }

            /* radial gate: |r_a - r_b| about the umbilicus */
            if (params->radial_armed) {
                float ray = (float)eps[i].y - params->umb_py;
                float rax = (float)eps[i].x - params->umb_px;
                float rby = (float)eps[j].y - params->umb_py;
                float rbx = (float)eps[j].x - params->umb_px;
                float ra = sqrtf(ray * ray + rax * rax);
                float rb = sqrtf(rby * rby + rbx * rbx);
                if (fabsf(ra - rb) > params->radial_dr_max) {
                    sm_add_reject(rej, &n_rej, stats, i, j,
                                  SLICEMATCH_REJ_RADIAL);
                    continue;
                }
            }

            float s_adj = sm_s_adj(arena, vol, D, H, W, z, &eps[i], &eps[j]);
            float st = s_track ? s_track[(size_t)i * (size_t)n_eps + (size_t)j]
                               : 0.0f;
            uint8_t far_tier = (uint8_t)is_far;
            /* Far-tier wrap-safety certificate: when the radial gate is
             * ARMED it is the certificate (a cross-wrap join must move ~one
             * pitch radially, forbidden at any reach) and evidence/track
             * stay advisory score terms. Unarmed runs (no umbilicus) keep
             * the strict corridor-evidence OR track-support requirement. */
            if (far_tier && !params->provisional && !params->radial_armed &&
                s_adj < FIXUP_ADJ_MIN_EVIDENCE &&
                st < FIXUP_TRACK_BOOST_MOD) {
                sm_add_reject(rej, &n_rej, stats, i, j,
                              SLICEMATCH_REJ_EVIDENCE);
                continue;
            }
            float s_dist = expf(-(d * d) /
                                (2.0f * FIXUP_DIST_SIGMA * FIXUP_DIST_SIGMA));
            float wa = (cos_a - need_cos) / (1.0f - need_cos);
            float wb = (cos_b - need_cos) / (1.0f - need_cos);
            if (wa > 1.0f) wa = 1.0f;
            if (wb > 1.0f) wb = 1.0f;
            float s_dir = 0.5f * (wa + wb);
            float score = FIXUP_W_DIST * s_dist + FIXUP_W_DIR * s_dir +
                          FIXUP_W_ADJ * s_adj + FIXUP_W_TRACK * st;

            /* provisional (round-1) floor is lower: long gaps must be able
             * to seed tracks or round 2 can never rescue them */
            float floor_score = params->min_score;
            if (params->provisional && FIXUP_MIN_SCORE_PROV < floor_score)
                floor_score = FIXUP_MIN_SCORE_PROV;
            if (score < floor_score) {
                sm_add_reject(rej, &n_rej, stats, i, j, SLICEMATCH_REJ_SCORE);
                continue;
            }
            if (!JoinPaint_arc_ratio_ok(&eps[i], &eps[j],
                                        FIXUP_MAX_ARC_RATIO)) {
                sm_add_reject(rej, &n_rej, stats, i, j, SLICEMATCH_REJ_ARC);
                continue;
            }
            if (!JoinPaint_merger_safe(mask_labels, H, W, &eps[i], &eps[j],
                                       FIXUP_MERGER_MARGIN)) {
                sm_add_reject(rej, &n_rej, stats, i, j, SLICEMATCH_REJ_MERGER);
                continue;
            }
            if (JoinPaint_crosses_adjacent_sheet(&eps[i], &eps[j], all_labels,
                                                 D, H, W, z,
                                                 FIXUP_CROSS_CHECK_RANGE)) {
                sm_add_reject(rej, &n_rej, stats, i, j, SLICEMATCH_REJ_XSHEET);
                continue;
            }

            sm_cand *c = &cands[n_cands++];
            c->a = i;
            c->b = j;
            c->score = score;
            c->dist = d;
            c->s_adj = s_adj;
            c->s_track = st;
            c->far_tier = far_tier;
        }
    }
    if (stats) stats->n_candidates = n_cands;

    /* Greedy stable matching on the symmetric score, with a live
     * no-crossing check against already-accepted joins. */
    qsort(cands, (size_t)n_cands, sizeof(*cands), sm_cand_cmp);

    SliceMatchJoin *joins =
        ARENA_ALLOC(arena, ((size_t)n_cands > 0 ? (size_t)n_cands : 1) *
                    sizeof(*joins));
    int32_t n_joins = 0;

    for (int32_t c = 0; c < n_cands; c++) {
        int32_t a = cands[c].a, b = cands[c].b;
        if (match_out[a] >= 0 || match_out[b] >= 0) {
            sm_add_reject(rej, &n_rej, stats, a, b, SLICEMATCH_REJ_OCCUPIED);
            continue;
        }
        int crosses = 0;
        float touch_dist = 2.0f * (float)FIXUP_PAINT_RADIUS + 1.0f;
        for (int32_t k = 0; k < n_joins && !crosses; k++) {
            crosses = JoinPaint_beziers_cross(&eps[a], &eps[b],
                                              &eps[joins[k].a],
                                              &eps[joins[k].b]) ||
                      JoinPaint_beziers_too_close(&eps[a], &eps[b],
                                                  &eps[joins[k].a],
                                                  &eps[joins[k].b],
                                                  touch_dist);
        }
        if (crosses) {
            sm_add_reject(rej, &n_rej, stats, a, b, SLICEMATCH_REJ_CROSSING);
            continue;
        }
        match_out[a] = b;
        match_out[b] = a;
        joins[n_joins].a = a;
        joins[n_joins].b = b;
        joins[n_joins].score = cands[c].score;
        joins[n_joins].dist = cands[c].dist;
        joins[n_joins].s_adj = cands[c].s_adj;
        joins[n_joins].s_track = cands[c].s_track;
        joins[n_joins].far_tier = cands[c].far_tier;
        n_joins++;
    }
    if (stats) stats->n_accepted = n_joins;

    *out_joins = joins;
    *out_n_joins = n_joins;
    *out_rejects = rej;
    *out_n_rejects = n_rej;
    return 0;
}

/* ---------------------------------------------------------------- */
/* Selftest                                                         */
/* ---------------------------------------------------------------- */

/* Paint a 1px-thick arc; th0..th1 radians, center (cy,cx), radius r. */
static void sm_test_arc(uint8_t *mask, int H, int W, float cy, float cx,
                        float r, float th0, float th1)
{
    int steps = (int)(r * 8.0f) + 16;
    for (int s = 0; s <= steps; s++) {
        float t = th0 + (th1 - th0) * (float)s / (float)steps;
        int iy = (int)(cy + r * sinf(t) + 0.5f);
        int ix = (int)(cx + r * cosf(t) + 0.5f);
        if (iy >= 0 && iy < H && ix >= 0 && ix < W)
            mask[sm_idx(iy, ix, W)] = 1;
    }
}

/* Full trace + match of plane z of vol. Returns joins via out params. */
static int sm_test_run(Arena_T arena, uint8_t *vol, int D, int H, int W, int z,
                       const SliceMatchParams *params,
                       SliceMatchJoin **out_joins, int32_t *out_n_joins,
                       SliceTraceEndpoint **out_eps, int32_t *out_n_eps,
                       SliceMatchStats *stats)
{
    size_t plane_sz = (size_t)H * (size_t)W;
    int32_t *all_labels = ARENA_ALLOC(arena, (size_t)D * plane_sz *
                                      sizeof(int32_t));
    for (int zz = 0; zz < D; zz++)
        SliceTrace_label_cc(arena, vol + (size_t)zz * plane_sz, H, W,
                            all_labels + (size_t)zz * plane_sz);
    SliceTraceEndpoint *eps = NULL;
    int32_t n_eps = 0;
    int rc = SliceTrace_run(arena, vol + (size_t)z * plane_sz, H, W, NULL,
                            all_labels + (size_t)z * plane_sz,
                            NULL, &eps, &n_eps);
    if (rc != 0) return rc;
    int32_t *match = ARENA_ALLOC(arena, ((size_t)n_eps > 0 ? (size_t)n_eps : 1)
                                 * sizeof(int32_t));
    SliceMatchReject *rejects = NULL;
    int32_t n_rejects = 0;
    rc = SliceMatch_run(arena, eps, n_eps, vol, all_labels, D, H, W, z,
                        all_labels + (size_t)z * plane_sz, NULL, params,
                        match, out_joins, out_n_joins, &rejects, &n_rejects,
                        stats);
    if (out_eps) *out_eps = eps;
    if (out_n_eps) *out_n_eps = n_eps;
    return rc;
}

int SliceMatch_selftest(void)
{
    int fails = 0;
    Arena_T arena = Arena_new();
    enum { H = 128, W = 128 };
    size_t plane_sz = (size_t)H * (size_t)W;

    SliceMatchParams params;
    memset(&params, 0, sizeof(params));
    params.reach_safe = FIXUP_REACH_SAFE_PX;
    params.reach_max = FIXUP_REACH_MAX_PX;
    params.min_score = FIXUP_MIN_SCORE;
    params.radial_armed = 0;
    params.radial_dr_max = FIXUP_RADIAL_DR_MAX;

    /* Case 1: broken C-arc, 3 planes (z0/z2 intact for evidence) -> 1 join.
     * The gap must split the curve into TWO components — a closed curve
     * with one gap stays one CC and the same-CC gate (correctly) refuses
     * to close the loop. A spiral cross-section is an open curve, so a
     * real gap always makes two pieces. */
    {
        Arena_Mark mark = Arena_save(arena);
        enum { D = 3 };
        uint8_t *vol = ARENA_CALLOC(arena, (size_t)D * plane_sz, 1);
        for (int zz = 0; zz < D; zz++) {
            uint8_t *pl = vol + (size_t)zz * plane_sz;
            if (zz == 1) {   /* two pieces: gap ~4.8px at angle pi. The C
                              * opening is kept > FIXUP_REACH_MAX_PX wide so
                              * only the pi gap is joinable. */
                sm_test_arc(pl, H, W, 64.0f, 64.0f, 40.0f, 0.50f, 3.08f);
                sm_test_arc(pl, H, W, 64.0f, 64.0f, 40.0f, 3.20f, 5.78f);
            } else {         /* intact C-arc */
                sm_test_arc(pl, H, W, 64.0f, 64.0f, 40.0f, 0.50f, 5.78f);
            }
        }
        SliceMatchJoin *joins = NULL;
        int32_t n_joins = 0;
        SliceMatchStats stats;
        int rc = sm_test_run(arena, vol, D, H, W, 1, &params, &joins, &n_joins,
                             NULL, NULL, &stats);
        int ok = (rc == 0 && n_joins == 1);
        fprintf(stderr, "[selftest] slice_match broken-arc (joins=%d) -> %s\n",
                (int)n_joins, ok ? "ok" : "FAIL");
        if (!ok) fails++;
        Arena_restore(arena, mark);
    }

    /* Case 2: two concentric arcs pitch 9.5 apart, aligned gaps, radial gate
     * armed -> 2 intra-arc joins, 0 cross-arc joins. */
    {
        Arena_Mark mark = Arena_save(arena);
        enum { D = 3 };
        uint8_t *vol = ARENA_CALLOC(arena, (size_t)D * plane_sz, 1);
        for (int zz = 0; zz < D; zz++) {
            uint8_t *pl = vol + (size_t)zz * plane_sz;
            if (zz == 1) {   /* aligned gaps at angle pi on both radii; the
                              * C openings stay wider than the reach cap */
                sm_test_arc(pl, H, W, 64.0f, 64.0f, 35.0f, 0.50f, 3.075f);
                sm_test_arc(pl, H, W, 64.0f, 64.0f, 35.0f, 3.205f, 5.78f);
                sm_test_arc(pl, H, W, 64.0f, 64.0f, 44.5f, 0.50f, 3.085f);
                sm_test_arc(pl, H, W, 64.0f, 64.0f, 44.5f, 3.195f, 5.78f);
            } else {
                sm_test_arc(pl, H, W, 64.0f, 64.0f, 35.0f, 0.50f, 5.78f);
                sm_test_arc(pl, H, W, 64.0f, 64.0f, 44.5f, 0.50f, 5.78f);
            }
        }
        SliceMatchParams p2 = params;
        p2.radial_armed = 1;
        p2.umb_py = 64.0f;
        p2.umb_px = 64.0f;
        SliceMatchJoin *joins = NULL;
        int32_t n_joins = 0;
        SliceTraceEndpoint *eps = NULL;
        int32_t n_eps = 0;
        int rc = sm_test_run(arena, vol, D, H, W, 1, &p2, &joins, &n_joins,
                             &eps, &n_eps, NULL);
        int n_cross_wrap = 0;
        for (int32_t k = 0; k < n_joins; k++) {
            float ra = sqrtf(((float)eps[joins[k].a].y - 64.0f) *
                             ((float)eps[joins[k].a].y - 64.0f) +
                             ((float)eps[joins[k].a].x - 64.0f) *
                             ((float)eps[joins[k].a].x - 64.0f));
            float rb = sqrtf(((float)eps[joins[k].b].y - 64.0f) *
                             ((float)eps[joins[k].b].y - 64.0f) +
                             ((float)eps[joins[k].b].x - 64.0f) *
                             ((float)eps[joins[k].b].x - 64.0f));
            if (fabsf(ra - rb) > 4.5f) n_cross_wrap++;
        }
        int ok = (rc == 0 && n_joins == 2 && n_cross_wrap == 0);
        if (!ok && rc == 0) {
            for (int32_t k = 0; k < n_joins; k++)
                fprintf(stderr,
                        "  [dbg] join %d: a(%d,%d) b(%d,%d) d=%.1f sc=%.2f\n",
                        (int)k,
                        (int)eps[joins[k].a].y, (int)eps[joins[k].a].x,
                        (int)eps[joins[k].b].y, (int)eps[joins[k].b].x,
                        (double)joins[k].dist, (double)joins[k].score);
            fprintf(stderr, "  [dbg] n_eps=%d\n", (int)n_eps);
        }
        fprintf(stderr,
                "[selftest] slice_match concentric-arcs (joins=%d cross=%d) -> %s\n",
                (int)n_joins, n_cross_wrap, ok ? "ok" : "FAIL");
        if (!ok) fails++;
        Arena_restore(arena, mark);
    }

    /* Case 3: same-CC reject — C-shape whose tips face each other. */
    {
        Arena_Mark mark = Arena_save(arena);
        enum { D = 1 };
        uint8_t *vol = ARENA_CALLOC(arena, plane_sz, 1);
        /* almost-closed circle: one CC, two tips ~5px apart */
        sm_test_arc(vol, H, W, 64.0f, 64.0f, 40.0f, 0.07f, 6.21f);
        SliceMatchJoin *joins = NULL;
        int32_t n_joins = 0;
        SliceMatchStats stats;
        int rc = sm_test_run(arena, vol, D, H, W, 0, &params, &joins, &n_joins,
                             NULL, NULL, &stats);
        int ok = (rc == 0 && n_joins == 0 &&
                  stats.rej_count[SLICEMATCH_REJ_SAME_CC] > 0);
        fprintf(stderr,
                "[selftest] slice_match same-cc-reject (joins=%d samecc=%d) -> %s\n",
                (int)n_joins, (int)stats.rej_count[SLICEMATCH_REJ_SAME_CC],
                ok ? "ok" : "FAIL");
        if (!ok) fails++;
        Arena_restore(arena, mark);
    }

    /* Case 4: far tier requires adjacent-plane evidence. Gap of ~9px:
     * with intact z+-1 -> join; with empty z+-1 -> REJ_EVIDENCE. */
    {
        Arena_Mark mark = Arena_save(arena);
        enum { D = 3 };
        uint8_t *vol = ARENA_CALLOC(arena, (size_t)D * plane_sz, 1);
        /* z=1: two collinear horizontal segments with a 9px gap */
        for (int x = 20; x <= 59; x++) vol[plane_sz + sm_idx(64, x, W)] = 1;
        for (int x = 69; x <= 108; x++) vol[plane_sz + sm_idx(64, x, W)] = 1;
        /* z=0, z=2: intact line */
        for (int x = 20; x <= 108; x++) {
            vol[sm_idx(64, x, W)] = 1;
            vol[2 * plane_sz + sm_idx(64, x, W)] = 1;
        }
        SliceMatchJoin *joins = NULL;
        int32_t n_joins = 0;
        int rc = sm_test_run(arena, vol, D, H, W, 1, &params, &joins, &n_joins,
                             NULL, NULL, NULL);
        int with_evidence = (rc == 0 && n_joins == 1 && joins[0].far_tier);

        /* clear the adjacent planes -> evidence gone -> no join */
        memset(vol, 0, plane_sz);
        memset(vol + 2 * plane_sz, 0, plane_sz);
        SliceMatchStats stats;
        rc = sm_test_run(arena, vol, D, H, W, 1, &params, &joins, &n_joins,
                         NULL, NULL, &stats);
        int without_evidence = (rc == 0 && n_joins == 0 &&
                                stats.rej_count[SLICEMATCH_REJ_EVIDENCE] > 0);
        int ok = (with_evidence && without_evidence);
        fprintf(stderr,
                "[selftest] slice_match far-tier-evidence (with=%d without=%d) -> %s\n",
                with_evidence, without_evidence, ok ? "ok" : "FAIL");
        if (!ok) fails++;
        Arena_restore(arena, mark);
    }

    /* Case 5: crossing ban — X configuration accepts only non-crossing set. */
    {
        Arena_Mark mark = Arena_save(arena);
        enum { D = 3 };
        uint8_t *vol = ARENA_CALLOC(arena, (size_t)D * plane_sz, 1);
        /* Four stubs pointing at the center of a 5px-radius clearing; the
         * natural pairs are left-right and top-bottom, which CROSS. Left-
         * right sits closer (higher score) so it wins; top-bottom must be
         * crossing-banned. Adjacent planes carry a full cross for evidence. */
        for (int zz = 0; zz < D; zz++) {
            uint8_t *pl = vol + (size_t)zz * plane_sz;
            int gap = (zz == 1) ? 3 : 0;
            for (int x = 24; x <= 64 - gap; x++) pl[sm_idx(64, x, W)] = 1;
            for (int x = 64 + gap; x <= 104; x++) pl[sm_idx(64, x, W)] = 1;
            for (int y = 24; y <= 64 - (gap + 1); y++) pl[sm_idx(y, 64, W)] = 1;
            for (int y = 64 + gap + 1; y <= 104; y++) pl[sm_idx(y, 64, W)] = 1;
        }
        SliceMatchJoin *joins = NULL;
        int32_t n_joins = 0;
        SliceMatchStats stats;
        int rc = sm_test_run(arena, vol, D, H, W, 1, &params, &joins, &n_joins,
                             NULL, NULL, &stats);
        /* at most one of the two crossing joins may be accepted */
        int ok = (rc == 0 && n_joins >= 1 &&
                  stats.rej_count[SLICEMATCH_REJ_CROSSING] +
                  stats.rej_count[SLICEMATCH_REJ_MERGER] +
                  stats.rej_count[SLICEMATCH_REJ_XSHEET] >= 1);
        if (ok && n_joins == 2) {
            /* if both survived they must not cross */
            SliceTraceEndpoint *eps2 = NULL;
            int32_t n2 = 0;
            sm_test_run(arena, vol, D, H, W, 1, &params, &joins, &n_joins,
                        &eps2, &n2, NULL);
            ok = !JoinPaint_beziers_cross(&eps2[joins[0].a], &eps2[joins[0].b],
                                          &eps2[joins[1].a], &eps2[joins[1].b]);
        }
        fprintf(stderr,
                "[selftest] slice_match crossing-ban (joins=%d banned=%d) -> %s\n",
                (int)n_joins, (int)stats.rej_count[SLICEMATCH_REJ_CROSSING],
                ok ? "ok" : "FAIL");
        if (!ok) fails++;
        Arena_restore(arena, mark);
    }

    Arena_dispose(&arena);
    fprintf(stderr, "=== slice_match selftest %s (%d failure%s) ===\n",
            fails ? "FAILED" : "passed", fails, fails == 1 ? "" : "s");
    return fails ? 3 : 0;
}
