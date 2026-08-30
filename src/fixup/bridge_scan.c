/* bridge_scan.c — thin-neck prediction-weld detector. See bridge_scan.h. */

#include "bridge_scan.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "../common/pipeline_constants.h"
#include "slice_trace.h"

static const int bs_dy8[8] = {-1, -1, -1,  0, 0,  1, 1, 1};
static const int bs_dx8[8] = {-1,  0,  1, -1, 1, -1, 0, 1};

static inline size_t bs_idx(int y, int x, int W)
{
    return (size_t)y * (size_t)W + (size_t)x;
}

/* 3-4 chamfer distance to background, scaled /3 to half-width in px. */
static void bs_chamfer(const uint8_t *mask, int H, int W, float *out)
{
    const float BIG = 1e6f;
    for (int y = 0; y < H; y++)
        for (int x = 0; x < W; x++)
            out[bs_idx(y, x, W)] = mask[bs_idx(y, x, W)] ? BIG : 0.0f;

    for (int y = 0; y < H; y++) {
        for (int x = 0; x < W; x++) {
            float v = out[bs_idx(y, x, W)];
            if (v == 0.0f) continue;
            if (y > 0) {
                float c = out[bs_idx(y - 1, x, W)] + 3.0f;
                if (c < v) v = c;
                if (x > 0) {
                    c = out[bs_idx(y - 1, x - 1, W)] + 4.0f;
                    if (c < v) v = c;
                }
                if (x < W - 1) {
                    c = out[bs_idx(y - 1, x + 1, W)] + 4.0f;
                    if (c < v) v = c;
                }
            }
            if (x > 0) {
                float c = out[bs_idx(y, x - 1, W)] + 3.0f;
                if (c < v) v = c;
            }
            out[bs_idx(y, x, W)] = v;
        }
    }
    for (int y = H - 1; y >= 0; y--) {
        for (int x = W - 1; x >= 0; x--) {
            float v = out[bs_idx(y, x, W)];
            if (v == 0.0f) continue;
            if (y < H - 1) {
                float c = out[bs_idx(y + 1, x, W)] + 3.0f;
                if (c < v) v = c;
                if (x > 0) {
                    c = out[bs_idx(y + 1, x - 1, W)] + 4.0f;
                    if (c < v) v = c;
                }
                if (x < W - 1) {
                    c = out[bs_idx(y + 1, x + 1, W)] + 4.0f;
                    if (c < v) v = c;
                }
            }
            if (x < W - 1) {
                float c = out[bs_idx(y, x + 1, W)] + 3.0f;
                if (c < v) v = c;
            }
            out[bs_idx(y, x, W)] = v;
        }
    }
    size_t n = (size_t)H * (size_t)W;
    for (size_t i = 0; i < n; i++) out[i] /= 3.0f;
}

int BridgeScan_run(Arena_T arena,
                   const uint8_t *mask, const uint8_t *skel, int H, int W,
                   int radial_armed, float umb_py, float umb_px,
                   BridgeScanHit **out_hits, int32_t *out_n_hits)
{
    assert(mask && skel && out_hits && out_n_hits);
    *out_hits = NULL;
    *out_n_hits = 0;
    if (H < 3 || W < 3) return 0;
    size_t n = (size_t)H * (size_t)W;

    /* results OUTLIVE the scratch mark — allocate them first */
    enum { MAX_HITS = 512 };
    BridgeScanHit *hits_tmp = ARENA_ALLOC(arena,
                                          (size_t)MAX_HITS * sizeof(*hits_tmp));
    int32_t n_hits = 0;

    Arena_Mark mark = Arena_save(arena);
    float *cham = ARENA_ALLOC(arena, n * sizeof(float));
    bs_chamfer(mask, H, W, cham);

    /* In-plane the predicted band is 1-3 px wide EVERYWHERE, so width can
     * never separate a weld from a sheet (measured: run2's thin-component
     * detector returned 3,391 junction crumbs, 0 real welds). A weld's true
     * signature is the RUNG-BETWEEN-RAILS pattern: a short skeleton segment
     * whose both ends are JUNCTIONS, whose own direction is radial, and
     * whose other branches at both junctions continue tangentially. */
    uint8_t *junc = ARENA_CALLOC(arena, n, 1);
    for (int y = 0; y < H; y++)
        for (int x = 0; x < W; x++) {
            size_t i = bs_idx(y, x, W);
            if (!skel[i]) continue;
            int deg = 0;
            for (int d = 0; d < 8; d++) {
                int ny = y + bs_dy8[d], nx = x + bs_dx8[d];
                if (ny >= 0 && ny < H && nx >= 0 && nx < W &&
                    skel[bs_idx(ny, nx, W)])
                    deg++;
            }
            junc[i] = (uint8_t)(deg >= 3);
        }

    uint8_t *seg_used = ARENA_CALLOC(arena, n, 1);

    for (int jy = 0; jy < H && n_hits < MAX_HITS; jy++) {
    for (int jx = 0; jx < W && n_hits < MAX_HITS; jx++) {
        if (!junc[bs_idx(jy, jx, W)]) continue;

        for (int d0 = 0; d0 < 8; d0++) {
            int sy = jy + bs_dy8[d0], sx = jx + bs_dx8[d0];
            if (sy < 0 || sy >= H || sx < 0 || sx >= W) continue;
            size_t si = bs_idx(sy, sx, W);
            if (!skel[si] || junc[si] || seg_used[si]) continue;

            /* walk the segment away from junction (jy,jx) */
            int path_y[FIXUP_BRIDGE_MAX_LEN + 2];
            int path_x[FIXUP_BRIDGE_MAX_LEN + 2];
            int len = 0;
            int cy = sy, cx = sx, py = jy, px = jx;
            int ky = -1, kx = -1;   /* terminating junction */
            while (len < FIXUP_BRIDGE_MAX_LEN) {
                path_y[len] = cy;
                path_x[len] = cx;
                len++;
                int ny_next = -1, nx_next = -1, n_next = 0;
                for (int d = 0; d < 8; d++) {
                    int ny = cy + bs_dy8[d], nx = cx + bs_dx8[d];
                    if (ny < 0 || ny >= H || nx < 0 || nx >= W) continue;
                    if (!skel[bs_idx(ny, nx, W)]) continue;
                    if (ny == py && nx == px) continue;
                    if (ny == jy && nx == jx) continue;
                    if (junc[bs_idx(ny, nx, W)]) { ky = ny; kx = nx; }
                    else { ny_next = ny; nx_next = nx; }
                    n_next++;
                }
                if (ky >= 0 || n_next != 1) break;
                py = cy; px = cx;
                cy = ny_next; cx = nx_next;
            }
            if (ky < 0) continue;               /* no junction at the far end */

            /* dedupe: claim the segment's interior pixels */
            for (int k = 0; k < len; k++)
                seg_used[bs_idx(path_y[k], path_x[k], W)] = 1;

            /* rung must be thin (hairline weld, not a confident sheet) */
            float wmax = 0.0f;
            for (int k = 0; k < len; k++)
                if (cham[bs_idx(path_y[k], path_x[k], W)] > wmax)
                    wmax = cham[bs_idx(path_y[k], path_x[k], W)];
            if (wmax > FIXUP_BRIDGE_MAX_WIDTH) continue;

            BridgeScanHit *h = &hits_tmp[n_hits];
            h->ay = jy; h->ax = jx;
            h->by = ky; h->bx = kx;
            h->my = (jy + ky) / 2;
            h->mx = (jx + kx) / 2;
            h->len = len;
            h->width = wmax;
            h->dr = -1.0f;
            h->radial_dot = -1.0f;
            h->certified = 0;

            if (radial_armed) {
                float ry = (float)h->my - umb_py;
                float rx = (float)h->mx - umb_px;
                float rl = hypotf(ry, rx);
                float ddy = (float)(ky - jy), ddx = (float)(kx - jx);
                float dl = hypotf(ddy, ddx);
                float ra = hypotf((float)jy - umb_py, (float)jx - umb_px);
                float rb = hypotf((float)ky - umb_py, (float)kx - umb_px);
                h->dr = fabsf(ra - rb);
                if (dl > 0.001f && rl > 0.001f) {
                    h->radial_dot = fabsf((ddy * ry + ddx * rx) / (dl * rl));

                    /* rails: at BOTH junctions, at least two OTHER branches
                     * must run tangentially */
                    int rails_ok = 1;
                    int ends_y[2] = { jy, ky }, ends_x[2] = { jx, kx };
                    for (int e = 0; e < 2 && rails_ok; e++) {
                        int ey = ends_y[e], ex = ends_x[e];
                        int n_tan = 0;
                        for (int d = 0; d < 8; d++) {
                            int ny = ey + bs_dy8[d], nx = ex + bs_dx8[d];
                            if (ny < 0 || ny >= H || nx < 0 || nx >= W)
                                continue;
                            if (!skel[bs_idx(ny, nx, W)]) continue;
                            /* skip the rung itself */
                            if ((ny == path_y[0] && nx == path_x[0]) ||
                                (ny == path_y[len - 1] &&
                                 nx == path_x[len - 1]))
                                continue;
                            if (ny == jy && nx == jx) continue;
                            if (ny == ky && nx == kx) continue;
                            /* walk 4 px along this branch for a direction */
                            int wy = ny, wx = nx, wpy = ey, wpx = ex;
                            for (int step = 0; step < 3; step++) {
                                int by2 = -1, bx2 = -1;
                                for (int dd = 0; dd < 8; dd++) {
                                    int zy = wy + bs_dy8[dd];
                                    int zx = wx + bs_dx8[dd];
                                    if (zy < 0 || zy >= H || zx < 0 ||
                                        zx >= W)
                                        continue;
                                    if (!skel[bs_idx(zy, zx, W)]) continue;
                                    if (zy == wpy && zx == wpx) continue;
                                    if (zy == ey && zx == ex) continue;
                                    by2 = zy; bx2 = zx;
                                    break;
                                }
                                if (by2 < 0) break;
                                wpy = wy; wpx = wx;
                                wy = by2; wx = bx2;
                            }
                            float bdy = (float)(wy - ey);
                            float bdx = (float)(wx - ex);
                            float bl = hypotf(bdy, bdx);
                            float ery = (float)ey - umb_py;
                            float erx = (float)ex - umb_px;
                            float erl = hypotf(ery, erx);
                            if (bl > 0.001f && erl > 0.001f) {
                                float bd = fabsf((bdy * ery + bdx * erx) /
                                                 (bl * erl));
                                if (bd <= 0.5f) n_tan++;
                            }
                        }
                        if (n_tan < 2) rails_ok = 0;
                    }

                    h->certified = (uint8_t)(rails_ok &&
                                             h->dr >= FIXUP_BRIDGE_MIN_DR &&
                                             h->radial_dot >=
                                                 FIXUP_BRIDGE_RADIAL_DOT);
                }
            }
            n_hits++;
            if (n_hits >= MAX_HITS) break;
        }
    }
    }

    Arena_restore(arena, mark);
    *out_hits = hits_tmp;
    *out_n_hits = n_hits;
    return 0;
}

int32_t BridgeScan_cut(uint8_t *mask, uint8_t *cut_mask, int H, int W,
                       const uint8_t *skel, const BridgeScanHit *hit)
{
    assert(mask && hit);
    (void)skel;
    float r = hit->width + 0.8f;
    float ay = (float)hit->ay, ax = (float)hit->ax;
    float by = (float)hit->by, bx = (float)hit->bx;
    float len = hypotf(by - ay, bx - ax);
    int steps = (int)(len * 2.0f) + 1;
    int32_t cleared = 0;
    int ri = (int)ceilf(r);
    for (int s = 0; s <= steps; s++) {
        float t = (float)s / (float)steps;
        float py = ay + t * (by - ay), px = ax + t * (bx - ax);
        int iy = (int)(py + 0.5f), ix = (int)(px + 0.5f);
        for (int dy = -ri; dy <= ri; dy++) {
            for (int dx = -ri; dx <= ri; dx++) {
                if ((float)(dy * dy + dx * dx) > r * r) continue;
                int ny = iy + dy, nx = ix + dx;
                if (ny < 0 || ny >= H || nx < 0 || nx >= W) continue;
                size_t idx = bs_idx(ny, nx, W);
                if (!mask[idx]) continue;
                mask[idx] = 0;
                if (cut_mask) cut_mask[idx] = 1;
                cleared++;
            }
        }
    }
    return cleared;
}

/* ---------------------------------------------------------------- */
/* Selftest                                                         */
/* ---------------------------------------------------------------- */

/* Trace helper: skeleton + labels for a synthetic mask. */
static int bst_scan(Arena_T arena, uint8_t *mask, int H, int W,
                    int armed, float uy, float ux,
                    BridgeScanHit **hits, int32_t *n_hits,
                    uint8_t **out_skel)
{
    int32_t *labels = ARENA_ALLOC(arena, (size_t)H * (size_t)W *
                                  sizeof(int32_t));
    SliceTrace_label_cc(arena, mask, H, W, labels);
    SliceTraceEndpoint *eps = NULL;
    int32_t n_eps = 0;
    uint8_t *skel = NULL;
    int rc = SliceTrace_run(arena, mask, H, W, NULL, labels, &skel,
                            &eps, &n_eps);
    if (rc != 0) return rc;
    if (out_skel) *out_skel = skel;
    return BridgeScan_run(arena, mask, skel, H, W, armed, uy, ux,
                          hits, n_hits);
}

static void bst_bar(uint8_t *mask, int W, int y0, int y1, int x0, int x1)
{
    for (int y = y0; y <= y1; y++)
        for (int x = x0; x <= x1; x++)
            mask[(size_t)y * (size_t)W + (size_t)x] = 1;
}

int BridgeScan_selftest(void)
{
    int fails = 0;
    Arena_T arena = Arena_new();
    enum { H = 96, W = 96 };
    size_t sz = (size_t)H * (size_t)W;

    /* Case 1: two thick bars + 1px radial bridge -> one CERTIFIED hit. */
    {
        Arena_Mark mark = Arena_save(arena);
        uint8_t *mask = ARENA_CALLOC(arena, sz, 1);
        bst_bar(mask, W, 32, 36, 10, 80);
        bst_bar(mask, W, 42, 46, 10, 80);
        bst_bar(mask, W, 37, 41, 40, 40);      /* the weld */
        BridgeScanHit *hits = NULL;
        int32_t n_hits = 0;
        int rc = bst_scan(arena, mask, H, W, 1, 200.0f, 40.0f,
                          &hits, &n_hits, NULL);
        int ok = (rc == 0 && n_hits == 1 && hits[0].certified &&
                  hits[0].dr > FIXUP_BRIDGE_MIN_DR &&
                  hits[0].radial_dot > FIXUP_BRIDGE_RADIAL_DOT);
        fprintf(stderr,
                "[selftest] bridge_scan radial-weld (n=%d cert=%d dr=%.1f "
                "dot=%.2f) -> %s\n",
                (int)n_hits, n_hits ? hits[0].certified : -1,
                n_hits ? (double)hits[0].dr : -1.0,
                n_hits ? (double)hits[0].radial_dot : -1.0,
                ok ? "ok" : "FAIL");
        if (!ok) fails++;
        Arena_restore(arena, mark);
    }

    /* Case 2: thick connector -> no thin neck, 0 hits. */
    {
        Arena_Mark mark = Arena_save(arena);
        uint8_t *mask = ARENA_CALLOC(arena, sz, 1);
        bst_bar(mask, W, 32, 36, 10, 80);
        bst_bar(mask, W, 42, 46, 10, 80);
        bst_bar(mask, W, 37, 41, 36, 44);      /* 9px-wide connector */
        BridgeScanHit *hits = NULL;
        int32_t n_hits = 0;
        int rc = bst_scan(arena, mask, H, W, 1, 200.0f, 40.0f,
                          &hits, &n_hits, NULL);
        int ok = (rc == 0 && n_hits == 0);
        fprintf(stderr, "[selftest] bridge_scan thick-connector (n=%d) -> %s\n",
                (int)n_hits, ok ? "ok" : "FAIL");
        if (!ok) fails++;
        Arena_restore(arena, mark);
    }

    /* Case 3: along-sheet thin section -> NOT a rung (no junctions bound
     * it), so zero hits. */
    {
        Arena_Mark mark = Arena_save(arena);
        uint8_t *mask = ARENA_CALLOC(arena, sz, 1);
        bst_bar(mask, W, 60, 64, 10, 80);
        for (int y = 60; y <= 64; y++)          /* pinch to 1px at x 38..44 */
            if (y != 62)
                for (int x = 38; x <= 44; x++)
                    mask[(size_t)y * (size_t)W + (size_t)x] = 0;
        BridgeScanHit *hits = NULL;
        int32_t n_hits = 0;
        int rc = bst_scan(arena, mask, H, W, 1, 200.0f, 40.0f,
                          &hits, &n_hits, NULL);
        int ok = (rc == 0 && n_hits == 0);
        (void)hits;
        fprintf(stderr,
                "[selftest] bridge_scan tangential-pinch (n=%d) -> %s\n",
                (int)n_hits, ok ? "ok" : "FAIL");
        if (!ok) fails++;
        Arena_restore(arena, mark);
    }

    /* Case 4: cutting the radial weld separates the bars. */
    {
        Arena_Mark mark = Arena_save(arena);
        uint8_t *mask = ARENA_CALLOC(arena, sz, 1);
        bst_bar(mask, W, 32, 36, 10, 80);
        bst_bar(mask, W, 42, 46, 10, 80);
        bst_bar(mask, W, 37, 41, 40, 40);
        int32_t *labels = ARENA_ALLOC(arena, sz * sizeof(int32_t));
        int32_t before = SliceTrace_label_cc(arena, mask, H, W, labels);
        BridgeScanHit *hits = NULL;
        int32_t n_hits = 0;
        uint8_t *skel = NULL;
        bst_scan(arena, mask, H, W, 1, 200.0f, 40.0f, &hits, &n_hits, &skel);
        int32_t cleared = 0;
        if (n_hits == 1)
            cleared = BridgeScan_cut(mask, NULL, H, W, skel, &hits[0]);
        int32_t after = SliceTrace_label_cc(arena, mask, H, W, labels);
        int ok = (before == 1 && after == 2 && cleared > 0 && cleared < 80);
        fprintf(stderr,
                "[selftest] bridge_scan cut (cc %d -> %d, cleared %d) -> %s\n",
                (int)before, (int)after, (int)cleared, ok ? "ok" : "FAIL");
        if (!ok) fails++;
        Arena_restore(arena, mark);
    }

    Arena_dispose(&arena);
    fprintf(stderr, "=== bridge_scan selftest %s (%d failure%s) ===\n",
            fails ? "FAILED" : "passed", fails, fails == 1 ? "" : "s");
    return fails ? 3 : 0;
}
