/* slice_trace.c — per-plane skeleton tracing for the 2D prediction gap fixup.
 * See slice_trace.h. Core algorithms ported from the recovered 2026-05
 * prototype (git 744d0ba^:scripts/step0.5-2d-trace/gap_fill_2d.c). */

#include "slice_trace.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "../common/pipeline_constants.h"
#include "../common/union_find.h"

/* 8-connectivity neighbor offsets; first 4 are the "prior" neighbors
 * (NW, N, NE, W) used by the raster-scan union pass. */
static const int st_dy8[8] = {-1, -1, -1,  0, 0,  1, 1, 1};
static const int st_dx8[8] = {-1,  0,  1, -1, 1, -1, 0, 1};

static inline size_t st_idx(int y, int x, int W)
{
    return (size_t)y * (size_t)W + (size_t)x;
}

/* ---------------------------------------------------------------- */
/* Zhang-Suen thinning                                              */
/* ---------------------------------------------------------------- */

/* Neighbor ordering P2..P9 clockwise from north. */
static const int st_zs_dy[8] = {-1, -1, 0, 1, 1,  1,  0, -1};
static const int st_zs_dx[8] = { 0,  1, 1, 1, 0, -1, -1, -1};

static inline int st_zs_p(const uint8_t *img, int H, int W, int y, int x, int i)
{
    int ny = y + st_zs_dy[i], nx = x + st_zs_dx[i];
    if (ny < 0 || ny >= H || nx < 0 || nx >= W) return 0;
    return img[st_idx(ny, nx, W)] ? 1 : 0;
}

static inline int st_zs_b(const uint8_t *img, int H, int W, int y, int x)
{
    int count = 0;
    for (int i = 0; i < 8; i++) count += st_zs_p(img, H, W, y, x, i);
    return count;
}

static inline int st_zs_a(const uint8_t *img, int H, int W, int y, int x)
{
    int count = 0;
    for (int i = 0; i < 8; i++) {
        int cur = st_zs_p(img, H, W, y, x, i);
        int nxt = st_zs_p(img, H, W, y, x, (i + 1) % 8);
        if (!cur && nxt) count++;
    }
    return count;
}

/* In-place Zhang-Suen thinning. mark is scratch [H*W]. Returns px removed. */
static int st_thin(uint8_t *skel, int H, int W, uint8_t *mark)
{
    size_t sz = (size_t)H * (size_t)W;
    int total_removed = 0;
    int changed = 1;

    while (changed) {
        changed = 0;
        for (int sub = 0; sub < 2; sub++) {
            memset(mark, 0, sz);
            for (int y = 1; y < H - 1; y++) {
                for (int x = 1; x < W - 1; x++) {
                    if (!skel[st_idx(y, x, W)]) continue;
                    int b = st_zs_b(skel, H, W, y, x);
                    if (b < 2 || b > 6) continue;
                    if (st_zs_a(skel, H, W, y, x) != 1) continue;
                    int p2 = st_zs_p(skel, H, W, y, x, 0);
                    int p4 = st_zs_p(skel, H, W, y, x, 2);
                    int p6 = st_zs_p(skel, H, W, y, x, 4);
                    int p8 = st_zs_p(skel, H, W, y, x, 6);
                    if (sub == 0) {
                        if (p2 * p4 * p6 != 0) continue;
                        if (p4 * p6 * p8 != 0) continue;
                    } else {
                        if (p2 * p4 * p8 != 0) continue;
                        if (p2 * p6 * p8 != 0) continue;
                    }
                    mark[st_idx(y, x, W)] = 1;
                }
            }
            for (int y = 1; y < H - 1; y++) {
                for (int x = 1; x < W - 1; x++) {
                    if (mark[st_idx(y, x, W)]) {
                        skel[st_idx(y, x, W)] = 0;
                        changed = 1;
                        total_removed++;
                    }
                }
            }
        }
    }
    return total_removed;
}

/* ---------------------------------------------------------------- */
/* Endpoint predicates                                              */
/* ---------------------------------------------------------------- */

static inline int st_count_n8(const uint8_t *img, int H, int W, int y, int x)
{
    int count = 0;
    for (int i = 0; i < 8; i++) {
        int ny = y + st_dy8[i], nx = x + st_dx8[i];
        if (ny >= 0 && ny < H && nx >= 0 && nx < W && img[st_idx(ny, nx, W)])
            count++;
    }
    return count;
}

/* Crossing number: 0->1 transitions in the circular 8-neighborhood. An
 * endpoint has crossing_number == 1 (all foreground neighbors form one
 * contiguous cluster) — this catches Zhang-Suen staircase tips where the tip
 * pixel has 2 adjacent neighbors. */
static inline int st_crossing_number(const uint8_t *img, int H, int W,
                                     int y, int x)
{
    int v[8];
    for (int i = 0; i < 8; i++) {
        int ny = y + st_zs_dy[i], nx = x + st_zs_dx[i];
        v[i] = (ny >= 0 && ny < H && nx >= 0 && nx < W &&
                img[st_idx(ny, nx, W)]) ? 1 : 0;
    }
    int cn = 0;
    for (int i = 0; i < 8; i++)
        if (v[i] == 0 && v[(i + 1) & 7] == 1) cn++;
    return cn;
}

static inline int st_is_endpoint(const uint8_t *img, int H, int W, int y, int x)
{
    return st_count_n8(img, H, W, y, x) >= 1 &&
           st_crossing_number(img, H, W, y, x) == 1;
}

/* ---------------------------------------------------------------- */
/* Spur pruning                                                     */
/* ---------------------------------------------------------------- */

/* Delete endpoint->junction branches shorter than min_len (junction pixels
 * themselves are preserved). Iterates until stable. */
static int st_prune(uint8_t *skel, int H, int W, int min_len, uint8_t *is_junc)
{
    size_t sz = (size_t)H * (size_t)W;
    int pruned = 0;

    memset(is_junc, 0, sz);
    for (int y = 0; y < H; y++)
        for (int x = 0; x < W; x++)
            if (skel[st_idx(y, x, W)] && st_count_n8(skel, H, W, y, x) >= 3)
                is_junc[st_idx(y, x, W)] = 1;

    int changed = 1;
    while (changed) {
        changed = 0;
        for (int y = 0; y < H; y++) {
            for (int x = 0; x < W; x++) {
                if (!skel[st_idx(y, x, W)]) continue;
                if (is_junc[st_idx(y, x, W)]) continue;
                if (!st_is_endpoint(skel, H, W, y, x)) continue;

                int path_y[32], path_x[32];
                int path_len = 0;
                int cy = y, cx = x;
                int prev_y = -1, prev_x = -1;
                int reached_junction = 0;

                while (path_len < min_len + 2 && path_len < 30) {
                    path_y[path_len] = cy;
                    path_x[path_len] = cx;
                    path_len++;

                    if (is_junc[st_idx(cy, cx, W)]) { reached_junction = 1; break; }

                    int nn = 0, next_y = -1, next_x = -1;
                    for (int i = 0; i < 8; i++) {
                        int ny = cy + st_dy8[i], nx = cx + st_dx8[i];
                        if (ny >= 0 && ny < H && nx >= 0 && nx < W &&
                            skel[st_idx(ny, nx, W)] &&
                            !(ny == prev_y && nx == prev_x)) {
                            nn++;
                            next_y = ny;
                            next_x = nx;
                        }
                    }
                    if (nn == 0) break;               /* isolated segment */
                    if (nn >= 2) { reached_junction = 1; break; }
                    prev_y = cy; prev_x = cx;
                    cy = next_y; cx = next_x;
                }

                if (reached_junction && path_len <= min_len) {
                    for (int i = 0; i < path_len; i++) {
                        if (!is_junc[st_idx(path_y[i], path_x[i], W)]) {
                            skel[st_idx(path_y[i], path_x[i], W)] = 0;
                            pruned++;
                        }
                    }
                    changed = 1;
                }
            }
        }
    }
    return pruned;
}

/* ---------------------------------------------------------------- */
/* 8-connected CC labeling                                          */
/* ---------------------------------------------------------------- */

int32_t SliceTrace_label_cc(Arena_T arena, const uint8_t *mask, int H, int W,
                            int32_t *labels)
{
    assert(mask && labels);
    if (H <= 0 || W <= 0) return 0;

    int32_t n = (int32_t)((size_t)H * (size_t)W);
    Arena_Mark mark = Arena_save(arena);
    UnionFind uf = UF_new(arena, n);

    for (int y = 0; y < H; y++) {
        for (int x = 0; x < W; x++) {
            size_t idx = st_idx(y, x, W);
            if (!mask[idx]) continue;
            /* prior neighbors: NW, N, NE, W */
            for (int d = 0; d < 4; d++) {
                int ny = y + st_dy8[d], nx = x + st_dx8[d];
                if (ny >= 0 && ny < H && nx >= 0 && nx < W &&
                    mask[st_idx(ny, nx, W)])
                    uf_union(&uf, (int32_t)idx, (int32_t)st_idx(ny, nx, W));
            }
        }
    }

    int32_t *remap = ARENA_CALLOC(arena, (size_t)n, sizeof(int32_t));
    int32_t cc_count = 0;
    for (int32_t i = 0; i < n; i++) {
        labels[i] = 0;
        if (!mask[i]) continue;
        int32_t root = uf_find(&uf, i);
        if (remap[root] == 0) remap[root] = ++cc_count;
        labels[i] = remap[root];
    }

    Arena_restore(arena, mark);
    return cc_count;
}

/* ---------------------------------------------------------------- */
/* Tangent + curvature walk                                         */
/* ---------------------------------------------------------------- */

/* Walk up to `steps` skeleton pixels from (sy,sx), preferring the neighbor
 * most aligned with the current direction at branchings. Writes the final
 * position and the position after steps/2 (for the curvature chord). */
static void st_walk(const uint8_t *skel, int H, int W, int sy, int sx,
                    int steps, int *out_end_y, int *out_end_x,
                    int *out_mid_y, int *out_mid_x)
{
    int prev_y = -1, prev_x = -1;
    int walk_y = sy, walk_x = sx;
    float dir_y = 0.0f, dir_x = 0.0f;
    int mid_y = sy, mid_x = sx;
    int n_walked = 0;

    for (int step = 0; step < steps; step++) {
        int cands_y[8], cands_x[8], ncands = 0;
        for (int d = 0; d < 8; d++) {
            int ny = walk_y + st_dy8[d], nx = walk_x + st_dx8[d];
            if (ny >= 0 && ny < H && nx >= 0 && nx < W &&
                skel[st_idx(ny, nx, W)] &&
                !(ny == prev_y && nx == prev_x)) {
                cands_y[ncands] = ny;
                cands_x[ncands] = nx;
                ncands++;
            }
        }
        if (ncands == 0) break;

        int best = 0;
        if (ncands > 1 && (dir_y != 0.0f || dir_x != 0.0f)) {
            float best_dot = -2.0f;
            for (int c = 0; c < ncands; c++) {
                float dy = (float)(cands_y[c] - walk_y);
                float dx = (float)(cands_x[c] - walk_x);
                float dot = dy * dir_y + dx * dir_x;
                if (dot > best_dot) { best_dot = dot; best = c; }
            }
        }

        prev_y = walk_y; prev_x = walk_x;
        walk_y = cands_y[best];
        walk_x = cands_x[best];
        dir_y = (float)(walk_y - prev_y);
        dir_x = (float)(walk_x - prev_x);
        n_walked++;
        if (n_walked == steps / 2) { mid_y = walk_y; mid_x = walk_x; }
    }

    *out_end_y = walk_y;
    *out_end_x = walk_x;
    if (out_mid_y) *out_mid_y = mid_y;
    if (out_mid_x) *out_mid_x = mid_x;
}

static void st_compute_tangent(const uint8_t *skel, int H, int W,
                               SliceTraceEndpoint *ep)
{
    /* Curvature: angle between endpoint->midpoint and midpoint->end chords
     * of a longer walk (multi-pixel chords avoid pixel aliasing). */
    int end_y = ep->y, end_x = ep->x, mid_y = ep->y, mid_x = ep->x;
    st_walk(skel, H, W, ep->y, ep->x, FIXUP_CURV_WALK,
            &end_y, &end_x, &mid_y, &mid_x);
    {
        float a_dy = (float)(mid_y - ep->y), a_dx = (float)(mid_x - ep->x);
        float b_dy = (float)(end_y - mid_y), b_dx = (float)(end_x - mid_x);
        float a_len = sqrtf(a_dy * a_dy + a_dx * a_dx);
        float b_len = sqrtf(b_dy * b_dy + b_dx * b_dx);
        if (a_len > 0.001f && b_len > 0.001f) {
            float dot = (a_dy * b_dy + a_dx * b_dx) / (a_len * b_len);
            if (dot < -1.0f) dot = -1.0f;
            if (dot >  1.0f) dot =  1.0f;
            ep->curv = acosf(dot);
        } else {
            ep->curv = 0.0f;
        }
    }

    /* Tangent: shorter walk; outward direction = last walked pixel back to
     * the endpoint (points into the gap). */
    st_walk(skel, H, W, ep->y, ep->x, FIXUP_TANGENT_WALK,
            &end_y, &end_x, NULL, NULL);
    {
        float ddy = (float)(ep->y - end_y);
        float ddx = (float)(ep->x - end_x);
        float len = sqrtf(ddy * ddy + ddx * ddx);
        if (len > 0.001f) {
            ep->tan_dy = ddy / len;
            ep->tan_dx = ddx / len;
        } else {
            ep->tan_dy = 0.0f;
            ep->tan_dx = 1.0f;
        }
    }
}

/* Walk the endpoint outward along its tangent through the FOREGROUND MASK
 * (8-connected, most-aligned neighbor) to the actual mask edge. Corrects the
 * Zhang-Suen tip retraction so gap distances are mask-edge to mask-edge. */
static void st_snap_to_edge(const uint8_t *mask, int H, int W,
                            SliceTraceEndpoint *ep)
{
    float tdy = ep->tan_dy, tdx = ep->tan_dx;
    int cy = ep->y, cx = ep->x;
    int prev_y = -1, prev_x = -1;

    for (int step = 0; step < FIXUP_TANGENT_WALK; step++) {
        int best_d = -1;
        float best_dot = 0.0f;      /* must be in the tangent half-plane */
        for (int d = 0; d < 8; d++) {
            int ny = cy + st_dy8[d], nx = cx + st_dx8[d];
            if (ny < 0 || ny >= H || nx < 0 || nx >= W) continue;
            if (ny == prev_y && nx == prev_x) continue;
            if (!mask[st_idx(ny, nx, W)]) continue;
            float dot = (float)st_dy8[d] * tdy + (float)st_dx8[d] * tdx;
            if (dot > best_dot) { best_dot = dot; best_d = d; }
        }
        if (best_d < 0) break;
        prev_y = cy; prev_x = cx;
        cy += st_dy8[best_d];
        cx += st_dx8[best_d];
    }

    ep->y = cy;
    ep->x = cx;
}

/* ---------------------------------------------------------------- */
/* SliceTrace_run                                                   */
/* ---------------------------------------------------------------- */

int SliceTrace_run(Arena_T arena,
                   const uint8_t *mask, int H, int W,
                   const uint8_t *absent,
                   const int32_t *mask_labels,
                   uint8_t **out_skel,
                   SliceTraceEndpoint **out_eps, int32_t *out_n_eps)
{
    assert(mask && mask_labels && out_eps && out_n_eps);
    *out_eps = NULL;
    *out_n_eps = 0;
    if (out_skel) *out_skel = NULL;
    if (H < 3 || W < 3) return 0;

    size_t sz = (size_t)H * (size_t)W;
    uint8_t *skel = ARENA_ALLOC(arena, sz);
    memcpy(skel, mask, sz);

    Arena_Mark scratch_mark = Arena_save(arena);
    uint8_t *scratch = ARENA_ALLOC(arena, sz);

    st_thin(skel, H, W, scratch);
    st_prune(skel, H, W, FIXUP_PRUNE_LEN, scratch);
    st_thin(skel, H, W, scratch);       /* re-thin junction leftovers */
    Arena_restore(arena, scratch_mark);

    /* Label skeleton CCs + component sizes (small CCs yield no endpoints). */
    int32_t *skel_labels = ARENA_ALLOC(arena, sz * sizeof(int32_t));
    int32_t skel_cc = SliceTrace_label_cc(arena, skel, H, W, skel_labels);

    /* eps is a RESULT — allocate it before the scratch mark so the restore
     * below cannot recycle it out from under the caller. */
    SliceTraceEndpoint *eps =
        ARENA_ALLOC(arena, (size_t)FIXUP_MAX_EPS_PER_PLANE * sizeof(*eps));
    int32_t n_eps = 0;

    Arena_Mark size_mark = Arena_save(arena);
    int32_t *cc_size = ARENA_CALLOC(arena, (size_t)skel_cc + 1, sizeof(int32_t));
    for (size_t i = 0; i < sz; i++)
        if (skel_labels[i] > 0) cc_size[skel_labels[i]]++;

    /* Detect endpoints (crossing-number == 1), skipping small components. */

    for (int y = 0; y < H && n_eps < FIXUP_MAX_EPS_PER_PLANE; y++) {
        for (int x = 0; x < W && n_eps < FIXUP_MAX_EPS_PER_PLANE; x++) {
            size_t idx = st_idx(y, x, W);
            if (!skel[idx]) continue;
            if (!st_is_endpoint(skel, H, W, y, x)) continue;
            int32_t cc = skel_labels[idx];
            if (cc <= 0 || cc_size[cc] < FIXUP_MIN_CURVE_PX) continue;
            SliceTraceEndpoint *ep = &eps[n_eps++];
            ep->y = y;
            ep->x = x;
            ep->skel_cc = cc;
            ep->mask_cc = 0;
            ep->tan_dy = 0.0f;
            ep->tan_dx = 0.0f;
            ep->curv = 0.0f;
            ep->excluded = 0;
        }
    }
    Arena_restore(arena, size_mark);

    /* Tangents, then snap to the mask edge, then classify exclusions. */
    for (int32_t i = 0; i < n_eps; i++) {
        st_compute_tangent(skel, H, W, &eps[i]);
        st_snap_to_edge(mask, H, W, &eps[i]);
        eps[i].mask_cc = mask_labels[st_idx(eps[i].y, eps[i].x, W)];

        int b = FIXUP_BORDER_EXCLUDE;
        if (eps[i].y < b || eps[i].y >= H - b ||
            eps[i].x < b || eps[i].x >= W - b) {
            eps[i].excluded = 1;
        } else if (absent) {
            for (int dy = -b; dy <= b && !eps[i].excluded; dy++)
                for (int dx = -b; dx <= b && !eps[i].excluded; dx++)
                    if (absent[st_idx(eps[i].y + dy, eps[i].x + dx, W)])
                        eps[i].excluded = 1;
        }
    }

    if (out_skel) *out_skel = skel;
    *out_eps = eps;
    *out_n_eps = n_eps;
    return 0;
}

/* ---------------------------------------------------------------- */
/* Selftest                                                         */
/* ---------------------------------------------------------------- */

/* Paint a thick arc of a circle into a mask. */
static void st_test_arc(uint8_t *mask, int H, int W, float cy, float cx,
                        float r, float a0, float a1, int thick)
{
    int steps = (int)(r * 8.0f) + 16;
    for (int s = 0; s <= steps; s++) {
        float t = a0 + (a1 - a0) * (float)s / (float)steps;
        float py = cy + r * sinf(t);
        float px = cx + r * cosf(t);
        int iy = (int)(py + 0.5f), ix = (int)(px + 0.5f);
        for (int dy = -thick; dy <= thick; dy++)
            for (int dx = -thick; dx <= thick; dx++) {
                int ny = iy + dy, nx = ix + dx;
                if (ny >= 0 && ny < H && nx >= 0 && nx < W &&
                    dy * dy + dx * dx <= thick * thick)
                    mask[st_idx(ny, nx, W)] = 1;
            }
    }
}

int SliceTrace_selftest(void)
{
    int fails = 0;
    Arena_T arena = Arena_new();
    enum { H = 96, W = 96 };
    size_t sz = (size_t)H * (size_t)W;

    /* Case 1: circle with one gap -> exactly 2 endpoints, tangents facing. */
    {
        Arena_Mark mark = Arena_save(arena);
        uint8_t *mask = ARENA_CALLOC(arena, sz, 1);
        /* circle minus a ~10px gap around angle 0 (within join reach; a
         * 20px gap on r=30 tilts the walk-estimated tangents too far) */
        st_test_arc(mask, H, W, 48.0f, 48.0f, 30.0f, 0.17f, 6.11f, 1);
        int32_t *labels = ARENA_ALLOC(arena, sz * sizeof(int32_t));
        SliceTrace_label_cc(arena, mask, H, W, labels);
        SliceTraceEndpoint *eps = NULL;
        int32_t n_eps = 0;
        int rc = SliceTrace_run(arena, mask, H, W, NULL, labels,
                                NULL, &eps, &n_eps);
        int ok = (rc == 0 && n_eps == 2);
        if (ok) {
            /* tangents should point roughly at each other across the gap */
            float dy = (float)(eps[1].y - eps[0].y);
            float dx = (float)(eps[1].x - eps[0].x);
            float len = sqrtf(dy * dy + dx * dx);
            if (len < 1.0f) ok = 0;
            else {
                float c0 = (eps[0].tan_dy * dy + eps[0].tan_dx * dx) / len;
                float c1 = -(eps[1].tan_dy * dy + eps[1].tan_dx * dx) / len;
                if (c0 < 0.7f || c1 < 0.7f) {
                    fprintf(stderr,
                            "  [dbg] ep0 (%d,%d) tan (%.2f,%.2f)  "
                            "ep1 (%d,%d) tan (%.2f,%.2f)  c0=%.2f c1=%.2f\n",
                            (int)eps[0].y, (int)eps[0].x,
                            (double)eps[0].tan_dy, (double)eps[0].tan_dx,
                            (int)eps[1].y, (int)eps[1].x,
                            (double)eps[1].tan_dy, (double)eps[1].tan_dx,
                            (double)c0, (double)c1);
                    ok = 0;
                }
            }
            if (eps[0].skel_cc != eps[1].skel_cc) {
                /* broken circle is one open curve = one skeleton CC */
                fprintf(stderr, "  [dbg] skel_cc %d vs %d\n",
                        (int)eps[0].skel_cc, (int)eps[1].skel_cc);
                ok = 0;
            }
        } else if (rc == 0) {
            for (int32_t i = 0; i < n_eps && i < 8; i++)
                fprintf(stderr, "  [dbg] ep%d (%d,%d) cc=%d excl=%d\n",
                        (int)i, (int)eps[i].y, (int)eps[i].x,
                        (int)eps[i].skel_cc, (int)eps[i].excluded);
        }
        fprintf(stderr, "[selftest] slice_trace broken-circle (n=%d) -> %s\n",
                (int)n_eps, ok ? "ok" : "FAIL");
        if (!ok) fails++;
        Arena_restore(arena, mark);
    }

    /* Case 2: empty plane -> 0 endpoints, no crash. */
    {
        Arena_Mark mark = Arena_save(arena);
        uint8_t *mask = ARENA_CALLOC(arena, sz, 1);
        int32_t *labels = ARENA_ALLOC(arena, sz * sizeof(int32_t));
        SliceTrace_label_cc(arena, mask, H, W, labels);
        SliceTraceEndpoint *eps = NULL;
        int32_t n_eps = -1;
        int rc = SliceTrace_run(arena, mask, H, W, NULL, labels,
                                NULL, &eps, &n_eps);
        int ok = (rc == 0 && n_eps == 0);
        fprintf(stderr, "[selftest] slice_trace empty-plane -> %s\n",
                ok ? "ok" : "FAIL");
        if (!ok) fails++;
        Arena_restore(arena, mark);
    }

    /* Case 3: line running off the border -> its border endpoint excluded. */
    {
        Arena_Mark mark = Arena_save(arena);
        uint8_t *mask = ARENA_CALLOC(arena, sz, 1);
        for (int x = 0; x < 40; x++) {          /* row from border to x=39 */
            mask[st_idx(48, x, W)] = 1;
            mask[st_idx(49, x, W)] = 1;
        }
        int32_t *labels = ARENA_ALLOC(arena, sz * sizeof(int32_t));
        SliceTrace_label_cc(arena, mask, H, W, labels);
        SliceTraceEndpoint *eps = NULL;
        int32_t n_eps = 0;
        SliceTrace_run(arena, mask, H, W, NULL, labels, NULL, &eps, &n_eps);
        int n_excluded = 0, n_free = 0;
        for (int32_t i = 0; i < n_eps; i++) {
            if (eps[i].excluded) n_excluded++;
            else n_free++;
        }
        int ok = (n_eps == 2 && n_excluded == 1 && n_free == 1);
        fprintf(stderr,
                "[selftest] slice_trace border-exclusion (n=%d excl=%d) -> %s\n",
                (int)n_eps, n_excluded, ok ? "ok" : "FAIL");
        if (!ok) fails++;
        Arena_restore(arena, mark);
    }

    /* Case 4: absent-region exclusion. */
    {
        Arena_Mark mark = Arena_save(arena);
        uint8_t *mask = ARENA_CALLOC(arena, sz, 1);
        uint8_t *absent = ARENA_CALLOC(arena, sz, 1);
        for (int x = 20; x < 60; x++) {
            mask[st_idx(48, x, W)] = 1;
            mask[st_idx(49, x, W)] = 1;
        }
        for (int y = 0; y < H; y++)              /* absent right half */
            for (int x = 60; x < W; x++)
                absent[st_idx(y, x, W)] = 1;
        int32_t *labels = ARENA_ALLOC(arena, sz * sizeof(int32_t));
        SliceTrace_label_cc(arena, mask, H, W, labels);
        SliceTraceEndpoint *eps = NULL;
        int32_t n_eps = 0;
        SliceTrace_run(arena, mask, H, W, absent, labels, NULL, &eps, &n_eps);
        int n_excluded = 0;
        for (int32_t i = 0; i < n_eps; i++)
            if (eps[i].excluded) n_excluded++;
        int ok = (n_eps == 2 && n_excluded == 1);
        fprintf(stderr,
                "[selftest] slice_trace absent-exclusion (n=%d excl=%d) -> %s\n",
                (int)n_eps, n_excluded, ok ? "ok" : "FAIL");
        if (!ok) fails++;
        Arena_restore(arena, mark);
    }

    Arena_dispose(&arena);
    fprintf(stderr, "=== slice_trace selftest %s (%d failure%s) ===\n",
            fails ? "FAILED" : "passed", fails, fails == 1 ? "" : "s");
    return fails ? 3 : 0;
}
