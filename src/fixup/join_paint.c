/* join_paint.c — cubic-Bezier join geometry, safety checks, painting.
 * See join_paint.h. Ported from the recovered 2026-05 prototype. */

#include "join_paint.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "../common/pipeline_constants.h"

static inline size_t jp_idx(int y, int x, int W)
{
    return (size_t)y * (size_t)W + (size_t)x;
}

/* Control points of the join Bezier. Returns chord length. */
static float jp_controls(const SliceTraceEndpoint *a,
                         const SliceTraceEndpoint *b,
                         float p[8] /* p0y,p0x,p1y,p1x,p2y,p2x,p3y,p3x */)
{
    float ddy = (float)(b->y - a->y);
    float ddx = (float)(b->x - a->x);
    float dist = sqrtf(ddy * ddy + ddx * ddx);
    p[0] = (float)a->y;
    p[1] = (float)a->x;
    p[2] = p[0] + (dist / 3.0f) * a->tan_dy;
    p[3] = p[1] + (dist / 3.0f) * a->tan_dx;
    p[6] = (float)b->y;
    p[7] = (float)b->x;
    p[4] = p[6] + (dist / 3.0f) * b->tan_dy;
    p[5] = p[7] + (dist / 3.0f) * b->tan_dx;
    return dist;
}

static inline void jp_eval(const float p[8], float t, float *py, float *px)
{
    float u = 1.0f - t;
    float uu = u * u, tt = t * t;
    *py = u * uu * p[0] + 3.0f * uu * t * p[2] + 3.0f * u * tt * p[4] + t * tt * p[6];
    *px = u * uu * p[1] + 3.0f * uu * t * p[3] + 3.0f * u * tt * p[5] + t * tt * p[7];
}

int JoinPaint_sample(const SliceTraceEndpoint *a, const SliceTraceEndpoint *b,
                     float *out_y, float *out_x)
{
    assert(a && b && out_y && out_x);
    float p[8];
    float dist = jp_controls(a, b, p);
    if (dist < 0.001f) {
        out_y[0] = (float)a->y; out_x[0] = (float)a->x;
        out_y[1] = (float)b->y; out_x[1] = (float)b->x;
        return 2;
    }
    int n = JOINPAINT_CROSS_SAMPLES + 1;
    for (int s = 0; s < n; s++) {
        float t = (float)s / (float)JOINPAINT_CROSS_SAMPLES;
        jp_eval(p, t, &out_y[s], &out_x[s]);
    }
    return n;
}

/* Proper segment crossing (parametric, shared endpoints excluded). */
static int jp_segments_cross(float ax, float ay, float bx, float by,
                             float cx, float cy, float dx, float dy)
{
    float d1x = bx - ax, d1y = by - ay;
    float d2x = dx - cx, d2y = dy - cy;
    float denom = d1x * d2y - d1y * d2x;
    if (fabsf(denom) < 1e-10f) return 0;    /* parallel */
    float t = ((cx - ax) * d2y - (cy - ay) * d2x) / denom;
    float u = ((cx - ax) * d1y - (cy - ay) * d1x) / denom;
    return (t > 0.001f && t < 0.999f && u > 0.001f && u < 0.999f);
}

int JoinPaint_beziers_cross(const SliceTraceEndpoint *a1,
                            const SliceTraceEndpoint *b1,
                            const SliceTraceEndpoint *a2,
                            const SliceTraceEndpoint *b2)
{
    float y1[JOINPAINT_CROSS_SAMPLES + 1], x1[JOINPAINT_CROSS_SAMPLES + 1];
    float y2[JOINPAINT_CROSS_SAMPLES + 1], x2[JOINPAINT_CROSS_SAMPLES + 1];
    int n1 = JoinPaint_sample(a1, b1, y1, x1);
    int n2 = JoinPaint_sample(a2, b2, y2, x2);

    for (int i = 0; i < n1 - 1; i++)
        for (int j = 0; j < n2 - 1; j++)
            if (jp_segments_cross(x1[i], y1[i], x1[i + 1], y1[i + 1],
                                  x2[j], y2[j], x2[j + 1], y2[j + 1]))
                return 1;
    return 0;
}

int JoinPaint_beziers_too_close(const SliceTraceEndpoint *a1,
                                const SliceTraceEndpoint *b1,
                                const SliceTraceEndpoint *a2,
                                const SliceTraceEndpoint *b2,
                                float min_dist)
{
    float y1[JOINPAINT_CROSS_SAMPLES + 1], x1[JOINPAINT_CROSS_SAMPLES + 1];
    float y2[JOINPAINT_CROSS_SAMPLES + 1], x2[JOINPAINT_CROSS_SAMPLES + 1];
    int n1 = JoinPaint_sample(a1, b1, y1, x1);
    int n2 = JoinPaint_sample(a2, b2, y2, x2);
    float d2 = min_dist * min_dist;

    for (int i = 0; i < n1; i++)
        for (int j = 0; j < n2; j++) {
            float dy = y1[i] - y2[j], dx = x1[i] - x2[j];
            if (dy * dy + dx * dx < d2) return 1;
        }
    return 0;
}

int JoinPaint_merger_safe(const int32_t *mask_labels, int H, int W,
                          const SliceTraceEndpoint *a,
                          const SliceTraceEndpoint *b, int margin)
{
    assert(mask_labels && a && b);
    float p[8];
    float dist = jp_controls(a, b, p);
    if (dist < 0.001f) return 1;

    int32_t cc_a = mask_labels[jp_idx(a->y, a->x, W)];
    int32_t cc_b = mask_labels[jp_idx(b->y, b->x, W)];

    int steps = (int)(dist * 2.0f) + 1;
    for (int s = 0; s <= steps; s++) {
        float t = (float)s / (float)steps;
        float py = 0.0f, px = 0.0f;
        jp_eval(p, t, &py, &px);
        int iy = (int)(py + 0.5f), ix = (int)(px + 0.5f);
        for (int my = -margin; my <= margin; my++) {
            for (int mx = -margin; mx <= margin; mx++) {
                int ny = iy + my, nx = ix + mx;
                if (ny < 0 || ny >= H || nx < 0 || nx >= W) continue;
                int32_t lbl = mask_labels[jp_idx(ny, nx, W)];
                if (lbl > 0 && lbl != cc_a && lbl != cc_b) return 0;
            }
        }
    }
    return 1;
}

int JoinPaint_arc_ratio_ok(const SliceTraceEndpoint *a,
                           const SliceTraceEndpoint *b, float max_ratio)
{
    float p[8];
    float chord = jp_controls(a, b, p);
    if (chord < 0.001f) return 1;

    enum { ARC_SAMPLES = 20 };
    float arc = 0.0f;
    float prev_y = p[0], prev_x = p[1];
    for (int s = 1; s <= ARC_SAMPLES; s++) {
        float t = (float)s / (float)ARC_SAMPLES;
        float py = 0.0f, px = 0.0f;
        jp_eval(p, t, &py, &px);
        float sy = py - prev_y, sx = px - prev_x;
        arc += sqrtf(sy * sy + sx * sx);
        prev_y = py;
        prev_x = px;
    }
    return arc / chord < max_ratio;
}

int JoinPaint_crosses_adjacent_sheet(const SliceTraceEndpoint *a,
                                     const SliceTraceEndpoint *b,
                                     const int32_t *all_labels,
                                     int D, int H, int W, int z, int range)
{
    assert(all_labels);
    float p[8];
    float dist = jp_controls(a, b, p);
    if (dist < 0.001f) return 0;

    size_t plane_sz = (size_t)H * (size_t)W;

    for (int zz = z - range; zz <= z + range; zz++) {
        if (zz < 0 || zz >= D || zz == z) continue;
        const int32_t *labels_zz = all_labels + (size_t)zz * plane_sz;

        int32_t label_a = labels_zz[jp_idx(a->y, a->x, W)];
        int32_t label_b = labels_zz[jp_idx(b->y, b->x, W)];

        /* Sample the middle 80% of the path (skip near the endpoints). */
        enum { NS = 20 };
        for (int s = 2; s <= NS - 2; s++) {
            float t = (float)s / (float)NS;
            float py = 0.0f, px = 0.0f;
            jp_eval(p, t, &py, &px);
            int iy = (int)(py + 0.5f), ix = (int)(px + 0.5f);
            if (iy < 0 || iy >= H || ix < 0 || ix >= W) continue;
            int32_t label_s = labels_zz[jp_idx(iy, ix, W)];
            if (label_s > 0 && label_s != label_a && label_s != label_b)
                return 1;
        }
    }
    return 0;
}

int32_t JoinPaint_draw(uint8_t *plane, uint8_t *painted, int H, int W,
                       const SliceTraceEndpoint *a,
                       const SliceTraceEndpoint *b, int radius)
{
    assert(plane && a && b);
    float p[8];
    float dist = jp_controls(a, b, p);
    if (dist < 0.001f) return 0;

    int steps = (int)(dist * 2.0f) + 1;
    int32_t added = 0;
    int r2 = radius * radius;

    for (int s = 0; s <= steps; s++) {
        float t = (float)s / (float)steps;
        float py = 0.0f, px = 0.0f;
        jp_eval(p, t, &py, &px);
        int iy = (int)(py + 0.5f), ix = (int)(px + 0.5f);
        for (int ry = -radius; ry <= radius; ry++) {
            for (int rx = -radius; rx <= radius; rx++) {
                if (ry * ry + rx * rx > r2) continue;
                int ny = iy + ry, nx = ix + rx;
                if (ny < 0 || ny >= H || nx < 0 || nx >= W) continue;
                size_t pos = jp_idx(ny, nx, W);
                if (!plane[pos]) {
                    plane[pos] = 1;
                    if (painted) painted[pos] = 1;
                    added++;
                }
            }
        }
    }
    return added;
}

/* ---------------------------------------------------------------- */
/* Selftest                                                         */
/* ---------------------------------------------------------------- */

static SliceTraceEndpoint jp_test_ep(int y, int x, float tdy, float tdx)
{
    SliceTraceEndpoint ep;
    memset(&ep, 0, sizeof(ep));
    ep.y = y;
    ep.x = x;
    float len = sqrtf(tdy * tdy + tdx * tdx);
    ep.tan_dy = (len > 0.0f) ? tdy / len : 0.0f;
    ep.tan_dx = (len > 0.0f) ? tdx / len : 1.0f;
    return ep;
}

int JoinPaint_selftest(void)
{
    int fails = 0;
    Arena_T arena = Arena_new();
    enum { H = 64, W = 64 };
    size_t sz = (size_t)H * (size_t)W;

    /* Case 1: crossing X vs parallel joins. (The X is deliberately
     * asymmetric so the intersection falls inside a sample segment, not
     * exactly on a shared sample point — the parametric test excludes
     * segment-endpoint contact, and the matcher layers the too-close
     * check on top for touching curves.) */
    {
        SliceTraceEndpoint a1 = jp_test_ep(20, 20, 1.0f, 1.0f);
        SliceTraceEndpoint b1 = jp_test_ep(40, 40, -1.0f, -1.0f);
        SliceTraceEndpoint a2 = jp_test_ep(44, 15, -1.0f, 1.0f);
        SliceTraceEndpoint b2 = jp_test_ep(17, 42, 1.0f, -1.0f);
        int cross = JoinPaint_beziers_cross(&a1, &b1, &a2, &b2);

        SliceTraceEndpoint a3 = jp_test_ep(20, 20, 0.0f, 1.0f);
        SliceTraceEndpoint b3 = jp_test_ep(20, 40, 0.0f, -1.0f);
        SliceTraceEndpoint a4 = jp_test_ep(40, 20, 0.0f, 1.0f);
        SliceTraceEndpoint b4 = jp_test_ep(40, 40, 0.0f, -1.0f);
        int par = JoinPaint_beziers_cross(&a3, &b3, &a4, &b4);

        int ok = (cross == 1 && par == 0);
        fprintf(stderr, "[selftest] join_paint beziers-cross (X=%d par=%d) -> %s\n",
                cross, par, ok ? "ok" : "FAIL");
        if (!ok) fails++;
    }

    /* Case 2: merger check — a third component between the endpoints. */
    {
        Arena_Mark mark = Arena_save(arena);
        int32_t *labels = ARENA_CALLOC(arena, sz, sizeof(int32_t));
        /* endpoint CCs at the two ends, an intruder CC in the middle */
        labels[jp_idx(32, 10, W)] = 1;
        labels[jp_idx(32, 50, W)] = 2;
        for (int y = 28; y <= 36; y++) labels[jp_idx(y, 30, W)] = 3;

        SliceTraceEndpoint a = jp_test_ep(32, 10, 0.0f, 1.0f);
        SliceTraceEndpoint b = jp_test_ep(32, 50, 0.0f, -1.0f);
        int unsafe = !JoinPaint_merger_safe(labels, H, W, &a, &b,
                                            FIXUP_MERGER_MARGIN);

        /* remove the intruder -> safe */
        for (int y = 28; y <= 36; y++) labels[jp_idx(y, 30, W)] = 0;
        int safe = JoinPaint_merger_safe(labels, H, W, &a, &b,
                                         FIXUP_MERGER_MARGIN);
        int ok = (unsafe && safe);
        fprintf(stderr, "[selftest] join_paint merger-check -> %s\n",
                ok ? "ok" : "FAIL");
        if (!ok) fails++;
        Arena_restore(arena, mark);
    }

    /* Case 3: arc ratio — facing tangents straight, hooked tangents bow. */
    {
        SliceTraceEndpoint a = jp_test_ep(32, 10, 0.0f, 1.0f);
        SliceTraceEndpoint b = jp_test_ep(32, 22, 0.0f, -1.0f);
        int straight = JoinPaint_arc_ratio_ok(&a, &b, FIXUP_MAX_ARC_RATIO);
        SliceTraceEndpoint c = jp_test_ep(32, 10, 1.0f, 0.3f);
        SliceTraceEndpoint d = jp_test_ep(32, 22, 1.0f, -0.3f);
        int hooked = JoinPaint_arc_ratio_ok(&c, &d, FIXUP_MAX_ARC_RATIO);
        int ok = (straight == 1 && hooked == 0);
        fprintf(stderr, "[selftest] join_paint arc-ratio (s=%d h=%d) -> %s\n",
                straight, hooked, ok ? "ok" : "FAIL");
        if (!ok) fails++;
    }

    /* Case 4: painting connects the two CCs and is additive-only. */
    {
        Arena_Mark mark = Arena_save(arena);
        uint8_t *plane = ARENA_CALLOC(arena, sz, 1);
        uint8_t *painted = ARENA_CALLOC(arena, sz, 1);
        for (int x = 4; x <= 20; x++) plane[jp_idx(32, x, W)] = 1;
        for (int x = 30; x <= 46; x++) plane[jp_idx(32, x, W)] = 1;
        size_t fg_before = 0;
        for (size_t i = 0; i < sz; i++) fg_before += plane[i];

        SliceTraceEndpoint a = jp_test_ep(32, 20, 0.0f, 1.0f);
        SliceTraceEndpoint b = jp_test_ep(32, 30, 0.0f, -1.0f);
        int32_t added = JoinPaint_draw(plane, painted, H, W, &a, &b,
                                       FIXUP_PAINT_RADIUS);

        size_t fg_after = 0, n_painted = 0;
        for (size_t i = 0; i < sz; i++) { fg_after += plane[i]; n_painted += painted[i]; }
        int ok = (added > 0 && fg_after == fg_before + (size_t)added &&
                  n_painted == (size_t)added);
        if (ok) {
            /* the two fragments must now be one 8-connected component */
            int32_t *labels = ARENA_ALLOC(arena, sz * sizeof(int32_t));
            int32_t ncc = SliceTrace_label_cc(arena, plane, H, W, labels);
            if (ncc != 1) ok = 0;
        }
        fprintf(stderr, "[selftest] join_paint draw-connects (added=%d) -> %s\n",
                (int)added, ok ? "ok" : "FAIL");
        if (!ok) fails++;
        Arena_restore(arena, mark);
    }

    /* Case 5: cross-sheet check — intruder in the adjacent plane corridor. */
    {
        Arena_Mark mark = Arena_save(arena);
        enum { D = 3 };
        int32_t *all = ARENA_CALLOC(arena, (size_t)D * sz, sizeof(int32_t));
        SliceTraceEndpoint a = jp_test_ep(32, 10, 0.0f, 1.0f);
        SliceTraceEndpoint b = jp_test_ep(32, 50, 0.0f, -1.0f);
        /* plane 0 (z-1 of z=1): third component in the middle of the path */
        for (int y = 30; y <= 34; y++)
            all[0 * (int)sz + (int)jp_idx(y, 30, W)] = 7;
        int hit = JoinPaint_crosses_adjacent_sheet(&a, &b, all, D, H, W, 1,
                                                   FIXUP_CROSS_CHECK_RANGE);
        /* clear -> no hit */
        for (int y = 30; y <= 34; y++)
            all[0 * (int)sz + (int)jp_idx(y, 30, W)] = 0;
        int clean = JoinPaint_crosses_adjacent_sheet(&a, &b, all, D, H, W, 1,
                                                     FIXUP_CROSS_CHECK_RANGE);
        int ok = (hit == 1 && clean == 0);
        fprintf(stderr, "[selftest] join_paint cross-sheet (hit=%d clean=%d) -> %s\n",
                hit, clean, ok ? "ok" : "FAIL");
        if (!ok) fails++;
        Arena_restore(arena, mark);
    }

    Arena_dispose(&arena);
    fprintf(stderr, "=== join_paint selftest %s (%d failure%s) ===\n",
            fails ? "FAILED" : "passed", fails, fails == 1 ? "" : "s");
    return fails ? 3 : 0;
}
