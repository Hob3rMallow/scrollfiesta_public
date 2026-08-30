/* fixup_viz.c — PNG viewers for the 2D prediction gap fixup.
 * See fixup_viz.h. Draw helpers ported from the recovered prototype. */

#include "fixup_viz.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "../common/pipeline_constants.h"
#include "../common/ves_png.h"
#include "join_paint.h"

static inline void fv_px(uint8_t *canvas, int cw, int ch, int x, int y,
                         uint8_t r, uint8_t g, uint8_t b)
{
    if (x < 0 || x >= cw || y < 0 || y >= ch) return;
    size_t idx = ((size_t)y * (size_t)cw + (size_t)x) * 3;
    canvas[idx] = r;
    canvas[idx + 1] = g;
    canvas[idx + 2] = b;
}

static void fv_dot(uint8_t *canvas, int cw, int ch, int x, int y, int radius,
                   uint8_t r, uint8_t g, uint8_t b)
{
    for (int dy = -radius; dy <= radius; dy++)
        for (int dx = -radius; dx <= radius; dx++)
            if (dy * dy + dx * dx <= radius * radius)
                fv_px(canvas, cw, ch, x + dx, y + dy, r, g, b);
}

static void fv_line(uint8_t *canvas, int cw, int ch,
                    int x0, int y0, int x1, int y1,
                    uint8_t r, uint8_t g, uint8_t b)
{
    int dx = x1 - x0;
    if (dx < 0) dx = -dx;
    int dy = y1 - y0;
    if (dy > 0) dy = -dy;
    int sx = x0 < x1 ? 1 : -1, sy = y0 < y1 ? 1 : -1;
    int err = dx + dy;
    for (;;) {
        fv_px(canvas, cw, ch, x0, y0, r, g, b);
        if (x0 == x1 && y0 == y1) break;
        int e2 = 2 * err;
        if (e2 >= dy) { err += dy; x0 += sx; }
        if (e2 <= dx) { err += dx; y0 += sy; }
    }
}

/* Draw a join's Bezier as a polyline. */
static void fv_join_curve(uint8_t *canvas, int cw, int ch, int scale,
                          const SliceTraceEndpoint *a,
                          const SliceTraceEndpoint *b,
                          uint8_t r, uint8_t g, uint8_t bl)
{
    float ys[JOINPAINT_CROSS_SAMPLES + 1], xs[JOINPAINT_CROSS_SAMPLES + 1];
    int n = JoinPaint_sample(a, b, ys, xs);
    for (int s = 0; s < n - 1; s++) {
        int x0 = (int)(xs[s] * (float)scale + 0.5f);
        int y0 = (int)(ys[s] * (float)scale + 0.5f);
        int x1 = (int)(xs[s + 1] * (float)scale + 0.5f);
        int y1 = (int)(ys[s + 1] * (float)scale + 0.5f);
        fv_line(canvas, cw, ch, x0, y0, x1, y1, r, g, bl);
    }
}

/* Render the overlay into an RGB canvas (arena). Exposed to the selftest. */
static uint8_t *fv_render_overlay(Arena_T arena,
                                  const uint8_t *plane_after,
                                  const uint8_t *painted,
                                  const uint8_t *skel,
                                  int H, int W,
                                  const SliceTraceEndpoint *eps, int32_t n_eps,
                                  const SliceMatchJoin *joins, int32_t n_joins,
                                  const int32_t *join_support,
                                  int32_t min_support,
                                  const SliceMatchReject *rejects,
                                  int32_t n_rejects,
                                  int scale, int *out_cw, int *out_ch)
{
    if (scale < 1) scale = 1;
    int cw = W * scale, ch = H * scale;
    uint8_t *canvas = ARENA_CALLOC(arena, (size_t)cw * (size_t)ch * 3, 1);

    /* base layers */
    for (int y = 0; y < H; y++) {
        for (int x = 0; x < W; x++) {
            size_t idx = (size_t)y * (size_t)W + (size_t)x;
            uint8_t r = 0, g = 0, b = 0;
            int is_painted = painted && painted[idx];
            if (is_painted) { r = 0; g = 255; b = 90; }
            else if (plane_after[idx]) { r = 150; g = 150; b = 150; }
            if (skel && skel[idx] && !is_painted) { r = 90; g = 110; b = 170; }
            if (r || g || b) {
                for (int sy = 0; sy < scale; sy++)
                    for (int sx = 0; sx < scale; sx++)
                        fv_px(canvas, cw, ch, x * scale + sx, y * scale + sy,
                              r, g, b);
            }
        }
    }

    /* reject lines under the accepted joins */
    for (int32_t k = 0; k < n_rejects; k++) {
        const SliceMatchReject *rj = &rejects[k];
        uint8_t r = 255, g = 150, b = 40;                 /* geometry: orange */
        switch (rj->reason) {
        case SLICEMATCH_REJ_MERGER:
        case SLICEMATCH_REJ_XSHEET:
        case SLICEMATCH_REJ_RADIAL:
            r = 255; g = 60; b = 60; break;               /* safety: red */
        case SLICEMATCH_REJ_CROSSING:
        case SLICEMATCH_REJ_OCCUPIED:
            r = 0; g = 180; b = 180; break;               /* competition: teal */
        default:
            break;
        }
        fv_line(canvas, cw, ch,
                eps[rj->a].x * scale, eps[rj->a].y * scale,
                eps[rj->b].x * scale, eps[rj->b].y * scale, r, g, b);
    }

    /* accepted joins (violet if dropped by the support filter) */
    for (int32_t k = 0; k < n_joins; k++) {
        int dropped = join_support && join_support[k] < min_support;
        if (dropped)
            fv_join_curve(canvas, cw, ch, scale,
                          &eps[joins[k].a], &eps[joins[k].b], 180, 80, 255);
        else
            fv_join_curve(canvas, cw, ch, scale,
                          &eps[joins[k].a], &eps[joins[k].b], 0, 255, 90);
    }

    /* endpoints on top */
    for (int32_t i = 0; i < n_eps; i++) {
        int ex = eps[i].x * scale, ey = eps[i].y * scale;
        int tx = ex + (int)(eps[i].tan_dx * 5.0f * (float)scale);
        int ty = ey + (int)(eps[i].tan_dy * 5.0f * (float)scale);
        fv_line(canvas, cw, ch, ex, ey, tx, ty, 0, 200, 255);
        if (eps[i].excluded)
            fv_dot(canvas, cw, ch, ex, ey, scale, 200, 110, 0);
        else
            fv_dot(canvas, cw, ch, ex, ey, scale, 255, 230, 40);
    }

    *out_cw = cw;
    *out_ch = ch;
    return canvas;
}

int FixupViz_plane_overlay(const char *path, Arena_T arena,
                           const uint8_t *plane_after, const uint8_t *painted,
                           const uint8_t *skel, int H, int W,
                           const SliceTraceEndpoint *eps, int32_t n_eps,
                           const SliceMatchJoin *joins, int32_t n_joins,
                           const int32_t *join_support, int32_t min_support,
                           const SliceMatchReject *rejects, int32_t n_rejects,
                           int scale)
{
    assert(path && plane_after);
    Arena_Mark mark = Arena_save(arena);
    int cw = 0, ch = 0;
    uint8_t *canvas = fv_render_overlay(arena, plane_after, painted, skel,
                                        H, W, eps, n_eps, joins, n_joins,
                                        join_support, min_support,
                                        rejects, n_rejects, scale, &cw, &ch);
    int rc = VesPng_write_rgb(path, canvas, cw, ch);
    Arena_restore(arena, mark);
    return rc;
}

int FixupViz_before_after(const char *path_before, const char *path_after,
                          Arena_T arena,
                          const uint8_t *plane_after, const uint8_t *painted,
                          int H, int W)
{
    assert(path_before && path_after && plane_after);
    Arena_Mark mark = Arena_save(arena);
    size_t sz = (size_t)H * (size_t)W;
    uint8_t *img = ARENA_ALLOC(arena, sz);
    for (size_t i = 0; i < sz; i++)
        img[i] = (uint8_t)((plane_after[i] && !(painted && painted[i]))
                           ? 255 : 0);
    int rc = VesPng_write_gray(path_before, img, W, H);
    for (size_t i = 0; i < sz; i++)
        img[i] = (uint8_t)(plane_after[i] ? 255 : 0);
    rc |= VesPng_write_gray(path_after, img, W, H);
    Arena_restore(arena, mark);
    return rc;
}

int FixupViz_pair_crop(const char *path, Arena_T arena,
                       const uint8_t *plane_after, const uint8_t *painted,
                       const uint8_t *skel, int H, int W,
                       const SliceTraceEndpoint *eps, int32_t n_eps,
                       int32_t a, int32_t b, int accepted,
                       uint8_t line_r, uint8_t line_g, uint8_t line_b,
                       int scale, int half)
{
    assert(path && plane_after && eps);
    if (scale < 1) scale = 1;
    if (half < 8) half = 8;

    int cy = (eps[a].y + eps[b].y) / 2;
    int cx = (eps[a].x + eps[b].x) / 2;
    int y0 = cy - half, y1 = cy + half;
    int x0 = cx - half, x1 = cx + half;
    if (y0 < 0) { y1 -= y0; y0 = 0; }
    if (x0 < 0) { x1 -= x0; x0 = 0; }
    if (y1 > H) { y0 -= (y1 - H); y1 = H; }
    if (x1 > W) { x0 -= (x1 - W); x1 = W; }
    if (y0 < 0) y0 = 0;
    if (x0 < 0) x0 = 0;
    int wh = y1 - y0, ww = x1 - x0;
    if (wh <= 0 || ww <= 0) return -1;

    enum { DIV = 4 };
    int pw = ww * scale, ph = wh * scale;
    int cw = pw * 2 + DIV, chh = ph;
    Arena_Mark mark = Arena_save(arena);
    uint8_t *canvas = ARENA_CALLOC(arena, (size_t)cw * (size_t)chh * 3, 1);

    /* divider */
    for (int y = 0; y < chh; y++)
        for (int d = 0; d < DIV; d++)
            fv_px(canvas, cw, chh, pw + d, y, 60, 60, 60);

    /* base layers: left = before, right = after */
    for (int y = y0; y < y1; y++) {
        for (int x = x0; x < x1; x++) {
            size_t idx = (size_t)y * (size_t)W + (size_t)x;
            int is_painted = painted && painted[idx];
            uint8_t lr2 = 0, lg2 = 0, lb2 = 0, rr = 0, rg = 0, rb = 0;
            if (plane_after[idx] && !is_painted) {
                lr2 = lg2 = lb2 = 150;
                rr = rg = rb = 150;
            }
            if (skel && skel[idx]) { lr2 = 90; lg2 = 110; lb2 = 170; }
            if (is_painted) { rr = 0; rg = 255; rb = 90; }
            for (int sy = 0; sy < scale; sy++) {
                for (int sx = 0; sx < scale; sx++) {
                    int py = (y - y0) * scale + sy;
                    int px = (x - x0) * scale + sx;
                    if (lr2 || lg2 || lb2)
                        fv_px(canvas, cw, chh, px, py, lr2, lg2, lb2);
                    if (rr || rg || rb)
                        fv_px(canvas, cw, chh, px + pw + DIV, py, rr, rg, rb);
                }
            }
        }
    }

    /* right panel: the connection */
    if (accepted) {
        float ys[JOINPAINT_CROSS_SAMPLES + 1], xs[JOINPAINT_CROSS_SAMPLES + 1];
        int n = JoinPaint_sample(&eps[a], &eps[b], ys, xs);
        for (int s = 0; s < n - 1; s++) {
            fv_line(canvas, cw, chh,
                    (int)((xs[s] - (float)x0) * (float)scale + 0.5f) + pw + DIV,
                    (int)((ys[s] - (float)y0) * (float)scale + 0.5f),
                    (int)((xs[s + 1] - (float)x0) * (float)scale + 0.5f) + pw + DIV,
                    (int)((ys[s + 1] - (float)y0) * (float)scale + 0.5f),
                    0, 255, 90);
        }
    } else {
        fv_line(canvas, cw, chh,
                (eps[a].x - x0) * scale + pw + DIV, (eps[a].y - y0) * scale,
                (eps[b].x - x0) * scale + pw + DIV, (eps[b].y - y0) * scale,
                line_r, line_g, line_b);
    }

    /* left panel: every endpoint in the window, ticks + dots */
    for (int32_t i = 0; i < n_eps; i++) {
        if (eps[i].y < y0 || eps[i].y >= y1 ||
            eps[i].x < x0 || eps[i].x >= x1) continue;
        int ex = (eps[i].x - x0) * scale, ey = (eps[i].y - y0) * scale;
        int tx = ex + (int)(eps[i].tan_dx * 5.0f * (float)scale);
        int ty = ey + (int)(eps[i].tan_dy * 5.0f * (float)scale);
        fv_line(canvas, cw, chh, ex, ey, tx, ty, 0, 200, 255);
        int hot = (i == a || i == b);
        if (eps[i].excluded)
            fv_dot(canvas, cw, chh, ex, ey, scale, 200, 110, 0);
        else
            fv_dot(canvas, cw, chh, ex, ey, hot ? scale + 1 : scale,
                   hot ? 255 : 255, hot ? 120 : 230, hot ? 120 : 40);
    }

    int rc = VesPng_write_rgb(path, canvas, cw, chh);
    Arena_restore(arena, mark);
    return rc;
}

int FixupViz_grid_summary(const char *path, Arena_T arena,
                          const uint8_t *vol, const uint8_t *painted,
                          int D, int H, int W)
{
    assert(path && vol);
    Arena_Mark mark = Arena_save(arena);
    size_t plane_sz = (size_t)H * (size_t)W;
    uint8_t *canvas = ARENA_CALLOC(arena, plane_sz * 3, 1);

    for (int z = 0; z < D; z++) {
        const uint8_t *pl = vol + (size_t)z * plane_sz;
        const uint8_t *pp = painted ? painted + (size_t)z * plane_sz : NULL;
        for (size_t i = 0; i < plane_sz; i++) {
            size_t c = i * 3;
            if (pp && pp[i]) {
                canvas[c] = 0;
                canvas[c + 1] = 255;
                canvas[c + 2] = 90;
            } else if (pl[i] && canvas[c + 1] != 255) {
                canvas[c] = 120;
                canvas[c + 1] = 120;
                canvas[c + 2] = 120;
            }
        }
    }

    int rc = VesPng_write_rgb(path, canvas, W, H);
    Arena_restore(arena, mark);
    return rc;
}

/* ---------------------------------------------------------------- */
/* Selftest (renders into memory; no file IO)                       */
/* ---------------------------------------------------------------- */

int FixupViz_selftest(void)
{
    int fails = 0;
    Arena_T arena = Arena_new();
    enum { H = 64, W = 64 };
    size_t sz = (size_t)H * (size_t)W;

    uint8_t *plane = ARENA_CALLOC(arena, sz, 1);
    uint8_t *painted = ARENA_CALLOC(arena, sz, 1);
    for (int x = 10; x <= 28; x++) plane[(size_t)32 * W + (size_t)x] = 1;
    for (int x = 36; x <= 54; x++) plane[(size_t)32 * W + (size_t)x] = 1;
    for (int x = 29; x <= 35; x++) {
        plane[(size_t)32 * W + (size_t)x] = 1;
        painted[(size_t)32 * W + (size_t)x] = 1;
    }
    /* one painted pixel off the tangent-tick row, for the color probe */
    plane[(size_t)31 * W + 32] = 1;
    painted[(size_t)31 * W + 32] = 1;

    SliceTraceEndpoint eps[2];
    memset(eps, 0, sizeof(eps));
    eps[0].y = 32; eps[0].x = 28; eps[0].tan_dx = 1.0f;
    eps[1].y = 32; eps[1].x = 36; eps[1].tan_dx = -1.0f;

    SliceMatchJoin join;
    memset(&join, 0, sizeof(join));
    join.a = 0;
    join.b = 1;

    int cw = 0, ch = 0;
    uint8_t *canvas = fv_render_overlay(arena, plane, painted, NULL, H, W,
                                        eps, 2, &join, 1, NULL, 0, NULL, 0,
                                        1, &cw, &ch);
    /* painted pixel green (probe off the tick row), mask gray, endpoint
     * yellow */
    size_t p_paint = (((size_t)31 * (size_t)W) + 32) * 3;
    size_t p_mask = (((size_t)32 * (size_t)W) + 15) * 3;
    int ok = (cw == W && ch == H &&
              canvas[p_paint] == 0 && canvas[p_paint + 1] == 255 &&
              canvas[p_mask] == 150);
    /* endpoint dot */
    size_t p_ep = (((size_t)32 * (size_t)W) + 28) * 3;
    if (!(canvas[p_ep] == 255 && canvas[p_ep + 1] == 230)) ok = 0;
    fprintf(stderr, "[selftest] fixup_viz overlay-colors -> %s\n",
            ok ? "ok" : "FAIL");
    if (!ok) fails++;

    /* scale=2 canvas doubles */
    canvas = fv_render_overlay(arena, plane, painted, NULL, H, W,
                               eps, 2, &join, 1, NULL, 0, NULL, 0,
                               2, &cw, &ch);
    ok = (cw == 2 * W && ch == 2 * H);
    fprintf(stderr, "[selftest] fixup_viz scale -> %s\n", ok ? "ok" : "FAIL");
    if (!ok) fails++;

    Arena_dispose(&arena);
    fprintf(stderr, "=== fixup_viz selftest %s (%d failure%s) ===\n",
            fails ? "FAILED" : "passed", fails, fails == 1 ? "" : "s");
    return fails ? 3 : 0;
}
