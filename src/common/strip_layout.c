/* strip_layout.c -- see strip_layout.h. */
#include "strip_layout.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "pipeline_constants.h"

#define SL_COLUMN_BLOCK 4096            /* columns per OpenMP work item in the profile */
#define SL_MIN_FLOW_WIDTH 64            /* narrowest line the aspect search tries, px */
#define SL_ASPECT_ITERATIONS 64

void StripLayout_params_default(StripParams *p, double du, double dv)
{
    assert(p != NULL);
    memset(p, 0, sizeof *p);
    p->du = du > 0.0 ? du : 1.0;
    p->dv = dv > 0.0 ? dv : 1.0;
    p->merge_gap_vox = 0.5 * ASM_READING_RIBBON_GAP;
    p->kibble_area_vox2 = STRIP_KIBBLE_AREA_VOX2;
    p->aspect = STRIP_PAGE_ASPECT;
    p->line_gap_vox = STRIP_LINE_GAP_VOX;
    p->piece_gap_vox = STRIP_PIECE_GAP_VOX;
    p->margin_vox = STRIP_MARGIN_VOX;
    p->tray_break_lines = STRIP_TRAY_BREAK_LINES;
    p->wrap_whole_fraction = STRIP_WRAP_WHOLE_FRACTION;
    p->soft_break_fraction = STRIP_SOFT_BREAK_FRACTION;
    p->min_room_fraction = STRIP_MIN_ROOM_FRACTION;
}

int StripLayout_columns(Arena_T arena, const uint8_t *cov, int32_t W, int32_t H,
                        StripColumns *out)
{
    if (arena == NULL || cov == NULL || out == NULL || W <= 0 || H <= 0) return -1;
    memset(out, 0, sizeof *out);
    int32_t *top = ARENA_ALLOC(arena, (size_t)W * sizeof *top);
    int32_t *bottom = ARENA_ALLOC(arena, (size_t)W * sizeof *bottom);
    int32_t *count = ARENA_ALLOC(arena, (size_t)W * sizeof *count);
    int nblocks = (int)(((int64_t)W + SL_COLUMN_BLOCK - 1) / SL_COLUMN_BLOCK);
    int b = 0;
#pragma omp parallel for schedule(dynamic, 1)
    for (b = 0; b < nblocks; b++) {
        int32_t x0 = (int32_t)b * SL_COLUMN_BLOCK;
        int32_t x1 = x0 + SL_COLUMN_BLOCK < W ? x0 + SL_COLUMN_BLOCK : W;
        for (int32_t x = x0; x < x1; x++) { top[x] = H; bottom[x] = 0; count[x] = 0; }
        for (int32_t y = 0; y < H; y++) {
            const uint8_t *row = cov + (size_t)y * (size_t)W;
            for (int32_t x = x0; x < x1; x++) {
                if (row[x] == 0) continue;
                if (top[x] == H) top[x] = y;
                bottom[x] = y + 1;
                count[x]++;
            }
        }
    }
    out->W = W; out->H = H; out->top = top; out->bottom = bottom; out->count = count;
    return 0;
}

/* Runs of covered columns; empty runs narrower than merge_px stay inside a
 * piece. out may be NULL (count only). */
static size_t sl_segments(const StripColumns *c, int32_t merge_px, double px_area, StripPiece *out)
{
    size_t n = 0;
    int open = 0;
    StripPiece cur;
    memset(&cur, 0, sizeof cur);
    for (int32_t x = 0; x < c->W; x++) {
        if (c->count[x] == 0) continue;
        if (open && x - cur.x1 < merge_px) {
            cur.x1 = x + 1;
            if (c->top[x] < cur.y0) cur.y0 = c->top[x];
            if (c->bottom[x] > cur.y1) cur.y1 = c->bottom[x];
            cur.covered += c->count[x];
            continue;
        }
        if (open) {
            cur.area_vox2 = (double)cur.covered * px_area;
            if (out != NULL) out[n] = cur;
            n++;
        }
        memset(&cur, 0, sizeof cur);
        cur.x0 = x; cur.x1 = x + 1; cur.y0 = c->top[x]; cur.y1 = c->bottom[x];
        cur.covered = c->count[x];
        open = 1;
    }
    if (open) {
        cur.area_vox2 = (double)cur.covered * px_area;
        if (out != NULL) out[n] = cur;
        n++;
    }
    return n;
}

typedef struct SlTrayKey { double area; int32_t x0; size_t index; } SlTrayKey;

static int sl_tray_order(const void *a, const void *b)
{
    const SlTrayKey *x = a, *y = b;
    if (x->area != y->area) return x->area > y->area ? -1 : 1;
    return (x->x0 > y->x0) - (x->x0 < y->x0);
}

/* One pass of the line flow at line width Lw. frags/lines NULL = dry run. */
typedef struct SlFlow {
    const StripColumns *c;
    const StripPiece *pieces;
    const size_t *order;        /* main pieces (ribbon order), then the tray (largest first) */
    size_t n, n_main;
    int32_t Lw, gap, min_room, wrap_whole, soft;
    StripFragment *frags;
    StripLine *lines;
    size_t n_frags, n_lines, n_split, tray_first_line;
    int64_t sum_heights;
    int32_t max_extent;
    int32_t last_main_extent;   /* extent of the last spiral (non-tray) line */
    /* current line */
    int32_t x, lh, lext;
    int line_open, tray;
} SlFlow;

static void sl_close_line(SlFlow *f)
{
    if (!f->line_open) return;
    if (f->lines != NULL) {
        StripLine *l = f->lines + f->n_lines;
        l->dy = 0; l->h = f->lh; l->tray = f->tray;
    }
    f->sum_heights += f->lh;
    if (f->lext > f->max_extent) f->max_extent = f->lext;
    if (!f->tray) f->last_main_extent = f->lext;
    f->n_lines++;
    f->x = 0; f->lh = 0; f->lext = 0; f->line_open = 0;
}

static void sl_add(SlFlow *f, size_t piece, int32_t sx0, int32_t sx1)
{
    assert(sx1 > sx0);
    int32_t y0 = INT32_MAX, y1 = 0;
    for (int32_t j = sx0; j < sx1; j++) {
        if (f->c->count[j] == 0) continue;
        if (f->c->top[j] < y0) y0 = f->c->top[j];
        if (f->c->bottom[j] > y1) y1 = f->c->bottom[j];
    }
    if (y1 <= y0) { y0 = 0; y1 = 0; }                        /* only empty columns */
    if (f->frags != NULL) {
        StripFragment *g = f->frags + f->n_frags;
        g->piece = (int32_t)piece; g->line = (int32_t)f->n_lines;
        g->sx0 = sx0; g->sx1 = sx1; g->sy0 = y0; g->sy1 = y1; g->dx = f->x;
    }
    f->n_frags++;
    if (y1 - y0 > f->lh) f->lh = y1 - y0;
    if (f->x + (sx1 - sx0) > f->lext) f->lext = f->x + (sx1 - sx0);
    f->line_open = 1;
}

static void sl_flow(SlFlow *f)
{
    f->n_frags = 0; f->n_lines = 0; f->n_split = 0; f->tray_first_line = SIZE_MAX;
    f->sum_heights = 0; f->max_extent = 0; f->last_main_extent = 0;
    f->x = 0; f->lh = 0; f->lext = 0; f->line_open = 0; f->tray = 0;
    for (size_t k = 0; k < f->n; k++) {
        size_t index = f->order[k];
        const StripPiece *pc = f->pieces + index;
        if (k == f->n_main) {                   /* the tray starts on a line of its own */
            sl_close_line(f);
            f->tray = 1;
            f->tray_first_line = f->n_lines;
        }
        int32_t sx = pc->x0;
        size_t parts = 0;
        while (sx < pc->x1) {
            int32_t w = pc->x1 - sx, room = f->Lw - f->x;
            if (w <= room) {
                sl_add(f, index, sx, pc->x1);
                parts++;
                f->x += w + f->gap;
                break;
            }
            int whole = f->tray ? w <= f->Lw : w <= f->wrap_whole;
            if (f->line_open && (room < f->min_room || whole)) {
                sl_close_line(f);
                continue;
            }
            /* Split: the break falls on the thinnest column near the end of
             * the room (rightmost of equals), so a cut crosses little material. */
            int32_t end = sx + room;                /* < pc->x1 because w > room */
            int32_t lo = end - f->soft;
            if (lo < sx + 1) lo = sx + 1;
            int32_t best = end, best_count = f->c->count[end];
            for (int32_t j = end - 1; j >= lo; j--) {
                if (f->c->count[j] < best_count) { best_count = f->c->count[j]; best = j; }
            }
            sl_add(f, index, sx, best);
            parts++;
            sx = best;
            sl_close_line(f);
        }
        if (parts > 1 && !f->tray) f->n_split++;
    }
    sl_close_line(f);
}

typedef struct SlGeometry { int32_t margin, line_gap, tray_gap; } SlGeometry;

static void sl_page_size(const SlFlow *f, const SlGeometry *g, int64_t *W, int64_t *H)
{
    size_t main_lines = f->tray_first_line == SIZE_MAX ? f->n_lines : f->tray_first_line;
    int64_t h = 2 * (int64_t)g->margin + f->sum_heights;
    if (f->n_lines > 1) h += (int64_t)g->line_gap * (int64_t)(f->n_lines - 1);
    if (main_lines > 0 && main_lines < f->n_lines) h += g->tray_gap;
    *W = 2 * (int64_t)g->margin + f->max_extent;
    *H = h;
}

static void sl_flow_setup(SlFlow *f, int32_t Lw, const StripParams *p, int32_t gap)
{
    f->Lw = Lw;
    f->gap = gap;
    f->min_room = (int32_t)lround(p->min_room_fraction * Lw);
    if (f->min_room < 1) f->min_room = 1;
    f->wrap_whole = (int32_t)lround(p->wrap_whole_fraction * Lw);
    f->soft = (int32_t)lround(p->soft_break_fraction * Lw);
    if (f->soft < 1) f->soft = 1;
}

int StripLayout_plan(Arena_T arena, const StripColumns *cols, const StripParams *p,
                     StripPage *out)
{
    if (arena == NULL || cols == NULL || p == NULL || out == NULL || cols->W <= 0 || cols->H <= 0 ||
        !(p->du > 0.0) || !(p->dv > 0.0) || !(p->aspect > 0.0)) return -1;
    memset(out, 0, sizeof *out);
    int32_t merge_px = (int32_t)lround(p->merge_gap_vox / p->du);
    if (merge_px < 1) merge_px = 1;
    double px_area = p->du * p->dv;
    size_t n = sl_segments(cols, merge_px, px_area, NULL);
    if (n == 0) return -1;
    StripPiece *pieces = ARENA_ALLOC(arena, n * sizeof *pieces);
    sl_segments(cols, merge_px, px_area, pieces);

    /* Main pieces in ribbon order, then the tray, largest first. */
    size_t *order = ARENA_ALLOC(arena, n * sizeof *order);
    SlTrayKey *keys = ARENA_ALLOC(arena, n * sizeof *keys);
    size_t n_main = 0, n_tray = 0;
    int64_t total_width = 0;
    for (size_t i = 0; i < n; i++) {
        StripPiece *pc = pieces + i;
        pc->tray = pc->area_vox2 < p->kibble_area_vox2;
        total_width += (int64_t)(pc->x1 - pc->x0);
        if (pc->tray) {
            keys[n_tray].area = pc->area_vox2; keys[n_tray].x0 = pc->x0; keys[n_tray].index = i;
            n_tray++;
            out->tray_area_vox2 += pc->area_vox2;
        } else {
            order[n_main++] = i;
            out->main_area_vox2 += pc->area_vox2;
        }
    }
    qsort(keys, n_tray, sizeof *keys, sl_tray_order);
    for (size_t i = 0; i < n_tray; i++) order[n_main + i] = keys[i].index;

    SlGeometry g;
    g.margin = (int32_t)lround(p->margin_vox / p->du);
    g.line_gap = (int32_t)lround(p->line_gap_vox / p->dv);
    g.tray_gap = (int32_t)lround(p->tray_break_lines * p->line_gap_vox / p->dv);
    int32_t gap = (int32_t)lround(p->piece_gap_vox / p->du);
    if (g.margin < 0 || g.line_gap < 0 || g.tray_gap < 0 || gap < 0) return -1;

    SlFlow f;
    memset(&f, 0, sizeof f);
    f.c = cols; f.pieces = pieces; f.order = order; f.n = n; f.n_main = n_main;

    /* Line width: bisect (log scale) towards the target aspect; keep the best. */
    int64_t hi64 = total_width + (int64_t)gap * (int64_t)n;
    if (hi64 > INT32_MAX / 2) hi64 = INT32_MAX / 2;
    int32_t lo = SL_MIN_FLOW_WIDTH, hi = (int32_t)hi64;
    if (hi < lo) hi = lo;
    int32_t best = hi;
    double best_err = HUGE_VAL;
    for (int it = 0; it < SL_ASPECT_ITERATIONS; it++) {
        int32_t mid = (int32_t)llround(sqrt((double)lo * (double)hi));
        if (mid < lo) mid = lo;
        if (mid > hi) mid = hi;
        sl_flow_setup(&f, mid, p, gap);
        sl_flow(&f);
        int64_t Wp = 0, Hp = 0;
        sl_page_size(&f, &g, &Wp, &Hp);
        double ratio = (double)Wp / (double)(Hp > 0 ? Hp : 1);
        double err = fabs(log(ratio / p->aspect));
        if (err < best_err) { best_err = err; best = mid; }
        if (hi - lo <= 1) break;
        if (ratio < p->aspect) lo = mid; else hi = mid;
    }
    /* Balance: near the best fit, the best aspect whose last spiral line is
     * at least STRIP_BALANCE_GOOD_FILL full, else the fullest last line. */
    if (n_main > 0) {
        int32_t chosen = best;
        double chosen_err = best_err, chosen_fill = -1.0;
        for (int i = 0; i <= STRIP_BALANCE_STEPS; i++) {
            double scale = 1.0 - STRIP_BALANCE_SPAN + 2.0 * STRIP_BALANCE_SPAN * (double)i / STRIP_BALANCE_STEPS;
            int32_t Lw = (int32_t)llround((double)best * scale);
            if (Lw < SL_MIN_FLOW_WIDTH) continue;
            sl_flow_setup(&f, Lw, p, gap);
            sl_flow(&f);
            int64_t Wp = 0, Hp = 0;
            sl_page_size(&f, &g, &Wp, &Hp);
            double err = fabs(log((double)Wp / (double)(Hp > 0 ? Hp : 1) / p->aspect));
            if (err > best_err + STRIP_BALANCE_ASPECT_SLACK) continue;
            double fill = fmin((double)f.last_main_extent / (double)Lw, STRIP_BALANCE_GOOD_FILL);
            if (fill > chosen_fill + 1e-9 || (fabs(fill - chosen_fill) <= 1e-9 && err < chosen_err)) {
                chosen = Lw; chosen_err = err; chosen_fill = fill;
            }
        }
        best = chosen;
    }

    /* Final flow: count, allocate, fill. */
    sl_flow_setup(&f, best, p, gap);
    sl_flow(&f);
    int64_t Wp = 0, Hp = 0;
    sl_page_size(&f, &g, &Wp, &Hp);
    if (Wp <= 0 || Hp <= 0 || Wp > INT32_MAX || Hp > INT32_MAX) return -1;
    f.frags = ARENA_ALLOC(arena, f.n_frags * sizeof *f.frags);
    f.lines = ARENA_ALLOC(arena, f.n_lines * sizeof *f.lines);
    sl_flow(&f);
    int32_t y = g.margin;
    for (size_t l = 0; l < f.n_lines; l++) {
        if (l == f.tray_first_line && l > 0) y += g.tray_gap;
        f.lines[l].dy = y;
        y += f.lines[l].h + g.line_gap;
    }
    out->pieces = pieces; out->n_pieces = n;
    out->frags = f.frags; out->n_frags = f.n_frags;
    out->lines = f.lines; out->n_lines = f.n_lines;
    out->W = (int32_t)Wp; out->H = (int32_t)Hp;
    out->flow_width = best; out->margin = g.margin;
    out->n_tray = n_tray; out->n_split = f.n_split;
    return 0;
}

int StripLayout_render(const StripPage *pg, const uint8_t *tex, const uint8_t *cov,
                       int32_t W, int32_t H, uint8_t *page, int64_t *out_copied)
{
    if (pg == NULL || tex == NULL || cov == NULL || page == NULL || W <= 0 || H <= 0) return -1;
    int64_t copied = 0;
    int bad = 0;
    int nf = (int)pg->n_frags, i = 0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+:copied, bad)
    for (i = 0; i < nf; i++) {
        const StripFragment *g = pg->frags + i;
        const StripLine *l = pg->lines + g->line;
        int32_t w = g->sx1 - g->sx0;
        int64_t px = (int64_t)pg->margin + g->dx;
        if (g->sx0 < 0 || g->sx1 > W || g->sy0 < 0 || g->sy1 > H || px + w > pg->W ||
            g->sy1 - g->sy0 > l->h || l->dy + l->h > pg->H) { bad++; continue; }
        for (int32_t sy = g->sy0; sy < g->sy1; sy++) {
            const uint8_t *st = tex + (size_t)sy * (size_t)W + (size_t)g->sx0;
            const uint8_t *sc = cov + (size_t)sy * (size_t)W + (size_t)g->sx0;
            uint8_t *d = page + (size_t)(l->dy + (sy - g->sy0)) * (size_t)pg->W + (size_t)px;
            for (int32_t x = 0; x < w; x++) {
                if (sc[x] == 0) continue;
                d[x] = st[x];
                copied++;
            }
        }
    }
    if (out_copied != NULL) *out_copied = copied;
    return bad ? -1 : 0;
}

int StripLayout_downsample(const uint8_t *img, int32_t W, int32_t H, int32_t f,
                           uint8_t *out, int32_t *out_w, int32_t *out_h)
{
    if (img == NULL || out == NULL || W <= 0 || H <= 0 || f < 1) return -1;
    int32_t ow = (W + f - 1) / f, oh = (H + f - 1) / f;
    int oy = 0;
#pragma omp parallel for schedule(dynamic, 16)
    for (oy = 0; oy < (int)oh; oy++) {
        int32_t y0 = (int32_t)oy * f, y1 = y0 + f < H ? y0 + f : H;
        for (int32_t ox = 0; ox < ow; ox++) {
            int32_t x0 = ox * f, x1 = x0 + f < W ? x0 + f : W;
            uint32_t sum = 0, cnt = 0;
            for (int32_t y = y0; y < y1; y++) {
                const uint8_t *row = img + (size_t)y * (size_t)W;
                for (int32_t x = x0; x < x1; x++) {
                    if (row[x] == 0) continue;
                    sum += row[x];
                    cnt++;
                }
            }
            out[(size_t)oy * (size_t)ow + (size_t)ox] = cnt ? (uint8_t)((sum + cnt / 2) / cnt) : 0;
        }
    }
    if (out_w != NULL) *out_w = ow;
    if (out_h != NULL) *out_h = oh;
    return 0;
}

int StripLayout_write_ledger(FILE *f, const StripPage *pg, const StripParams *p)
{
    if (f == NULL || pg == NULL || p == NULL) return -1;
    int ok = 1;
    ok &= fprintf(f, "\"params\":{\"du\":%.17g,\"dv\":%.17g,\"merge_gap_vox\":%.17g,\"kibble_area_vox2\":%.17g,"
                     "\"aspect\":%.17g,\"line_gap_vox\":%.17g,\"piece_gap_vox\":%.17g,\"margin_vox\":%.17g,"
                     "\"tray_break_lines\":%.17g,\"wrap_whole_fraction\":%.17g,\"soft_break_fraction\":%.17g,"
                     "\"min_room_fraction\":%.17g},\n",
                  p->du, p->dv, p->merge_gap_vox, p->kibble_area_vox2, p->aspect, p->line_gap_vox, p->piece_gap_vox,
                  p->margin_vox, p->tray_break_lines, p->wrap_whole_fraction, p->soft_break_fraction,
                  p->min_room_fraction) > 0;
    ok &= fprintf(f, "\"page_w\":%d,\"page_h\":%d,\"flow_width\":%d,\"margin\":%d,\"pieces\":%zu,\"main_pieces\":%zu,"
                     "\"tray_pieces\":%zu,\"split_main_pieces\":%zu,\"lines\":%zu,\"fragments\":%zu,"
                     "\"main_area_vox2\":%.17g,\"tray_area_vox2\":%.17g,\n",
                  (int)pg->W, (int)pg->H, (int)pg->flow_width, (int)pg->margin, pg->n_pieces,
                  pg->n_pieces - pg->n_tray, pg->n_tray, pg->n_split, pg->n_lines, pg->n_frags,
                  pg->main_area_vox2, pg->tray_area_vox2) > 0;
    ok &= fputs("\"piece_fields\":[\"x0\",\"x1\",\"y0\",\"y1\",\"covered_px\",\"area_vox2\",\"tray\"],\n\"piece_list\":[", f) >= 0;
    for (size_t i = 0; i < pg->n_pieces && ok; i++) {
        const StripPiece *pc = pg->pieces + i;
        ok &= fprintf(f, "%s[%d,%d,%d,%d,%lld,%.1f,%d]", i ? "," : "", (int)pc->x0, (int)pc->x1, (int)pc->y0,
                      (int)pc->y1, (long long)pc->covered, pc->area_vox2, (int)pc->tray) > 0;
    }
    ok &= fputs("],\n\"line_fields\":[\"dy\",\"h\",\"tray\"],\n\"line_list\":[", f) >= 0;
    for (size_t i = 0; i < pg->n_lines && ok; i++) {
        const StripLine *l = pg->lines + i;
        ok &= fprintf(f, "%s[%d,%d,%d]", i ? "," : "", (int)l->dy, (int)l->h, (int)l->tray) > 0;
    }
    ok &= fputs("],\n\"fragment_fields\":[\"piece\",\"line\",\"sx0\",\"sx1\",\"sy0\",\"sy1\",\"dx\"],\n"
                "\"fragment_list\":[", f) >= 0;
    for (size_t i = 0; i < pg->n_frags && ok; i++) {
        const StripFragment *g = pg->frags + i;
        ok &= fprintf(f, "%s[%d,%d,%d,%d,%d,%d,%d]", i ? "," : "", (int)g->piece, (int)g->line, (int)g->sx0,
                      (int)g->sx1, (int)g->sy0, (int)g->sy1, (int)g->dx) > 0;
    }
    ok &= fputs("]", f) >= 0;
    return ok ? 0 : -1;
}

/* ------------------------------------------------------------------ selftest */

static void sl_fill(uint8_t *cov, uint8_t *tex, int32_t W, int32_t x0, int32_t x1, int32_t y0, int32_t y1)
{
    for (int32_t y = y0; y < y1; y++) {
        for (int32_t x = x0; x < x1; x++) {
            cov[(size_t)y * (size_t)W + (size_t)x] = 255;
            tex[(size_t)y * (size_t)W + (size_t)x] = (uint8_t)(1 + (x * 7 + y * 13) % 250);
        }
    }
}

static int sl_check(int cond, const char *what, int *fails)
{
    if (!cond) {
        fprintf(stderr, "[strip_layout selftest]   FAIL: %s\n", what);
        (*fails)++;
    }
    return cond;
}

int StripLayout_selftest(void)
{
    int fails = 0;
    Arena_T arena = Arena_new();

    /* 1. Segmentation, tray order, conservation, line crop. Pieces (du = 1):
     *   A cols 10-59 rows 5-40 (1,800 px), B cols 100-139 rows 10-20 with an
     *   internal 5-column gap (one piece), C cols 170-179 rows 0-49 (500 px),
     *   D cols 220-229 rows 1-3 (30 px, kibble), E cols 260-269 rows 0-9
     *   (100 px, kibble): tray order E then D. */
    {
        const int32_t W = 400, H = 50;
        uint8_t *cov = ARENA_CALLOC(arena, (size_t)W * H, 1), *tex = ARENA_CALLOC(arena, (size_t)W * H, 1);
        sl_fill(cov, tex, W, 10, 60, 5, 41);
        sl_fill(cov, tex, W, 100, 120, 10, 21);
        sl_fill(cov, tex, W, 125, 140, 10, 21);
        sl_fill(cov, tex, W, 170, 180, 0, 50);
        sl_fill(cov, tex, W, 220, 230, 1, 4);
        sl_fill(cov, tex, W, 260, 270, 0, 10);
        StripColumns cols = {0};
        StripParams p;
        StripPage pg = {0};
        StripLayout_params_default(&p, 1.0, 1.0);
        p.kibble_area_vox2 = 200.0;
        p.line_gap_vox = 4.0; p.piece_gap_vox = 3.0; p.margin_vox = 2.0; p.tray_break_lines = 2.0;
        int rc = StripLayout_columns(arena, cov, W, H, &cols);
        rc = rc || StripLayout_plan(arena, &cols, &p, &pg);
        if (sl_check(rc == 0, "plan of five pieces", &fails)) {
            sl_check(pg.n_pieces == 5, "internal 5-column gap stays inside piece B (5 pieces)", &fails);
            sl_check(pg.n_tray == 2, "two kibble pieces", &fails);
            sl_check(pg.pieces[1].x0 == 100 && pg.pieces[1].x1 == 140, "piece B spans its gap", &fails);
            /* main fragments keep ribbon order; tray fragments follow, E before D */
            int32_t last_piece = -1, order_ok = 1, tray_seen = 0, e_line = -1, d_line = -1, e_dx = -1, d_dx = -1;
            for (size_t i = 0; i < pg.n_frags; i++) {
                const StripFragment *g = pg.frags + i;
                int tray = pg.pieces[g->piece].tray;
                if (tray) tray_seen = 1;
                else if (tray_seen || g->piece < last_piece) order_ok = 0;
                if (!tray) last_piece = g->piece;
                if (g->piece == 4) { e_line = g->line; e_dx = g->dx; }
                if (g->piece == 3) { d_line = g->line; d_dx = g->dx; }
            }
            sl_check(order_ok, "main pieces in ribbon order before the tray", &fails);
            sl_check(e_line >= 0 && d_line >= 0 && (e_line < d_line || (e_line == d_line && e_dx < d_dx)),
                     "tray is largest first (E before D)", &fails);
            sl_check(pg.lines[d_line].tray == 1, "kibble lines are tray lines", &fails);
            uint8_t *page = ARENA_CALLOC(arena, (size_t)pg.W * (size_t)pg.H, 1);
            int64_t copied = 0, covered = 0;
            int64_t sum_src = 0, sum_page = 0;
            for (size_t i = 0; i < (size_t)W * H; i++) { covered += cov[i] != 0; sum_src += tex[i]; }
            rc = StripLayout_render(&pg, tex, cov, W, H, page, &copied);
            for (size_t i = 0; i < (size_t)pg.W * pg.H; i++) sum_page += page[i];
            sl_check(rc == 0 && copied == covered, "every covered pixel copied once", &fails);
            sl_check(sum_page == sum_src, "pixel values unchanged", &fails);
            /* a line holding only E (rows 0-9) and D (rows 1-3) is 10 rows tall and
             * D hangs from its top by its own first row */
            if (e_line == d_line) {
                const StripLine *l = pg.lines + e_line;
                sl_check(l->h == 10, "tray line as tall as its tallest fragment", &fails);
                for (size_t i = 0; i < pg.n_frags; i++)
                    if (pg.frags[i].piece == 3)
                        sl_check(pg.frags[i].sy0 == 1 && pg.frags[i].sy1 == 4, "fragment rows are its covered rows", &fails);
            }
        }
    }

    /* 2. Flow rules at a fixed width (Lw 300): soft break on the thinnest
     *    column, whole moves for narrow main pieces, unsplit tray pieces. */
    {
        const int32_t W = 2000, H = 20;
        StripColumns c = {0};
        c.W = W; c.H = H;
        c.top = ARENA_ALLOC(arena, (size_t)W * sizeof *c.top);
        c.bottom = ARENA_ALLOC(arena, (size_t)W * sizeof *c.bottom);
        c.count = ARENA_ALLOC(arena, (size_t)W * sizeof *c.count);
        for (int32_t x = 0; x < W; x++) { c.top[x] = 0; c.bottom[x] = H; c.count[x] = 10; }
        c.count[290] = 1;                               /* thin column inside the window [282, 300] */
        /* The long piece breaks at 290, 590, 890; its last 275 columns end the
         * third line at x 280, so the narrow piece (30 < 12% of 300) finds 20
         * columns of room (>= the 15-column minimum) and moves whole. */
        StripPiece pc[3];
        memset(pc, 0, sizeof pc);
        pc[0].x0 = 0;    pc[0].x1 = 1165;               /* long main piece */
        pc[1].x0 = 1200; pc[1].x1 = 1230;               /* narrow main piece */
        pc[2].x0 = 1300; pc[2].x1 = 1500; pc[2].tray = 1; /* tray piece, 200 wide */
        size_t order[3] = {0, 1, 2};
        StripParams p;
        StripLayout_params_default(&p, 1.0, 1.0);
        p.soft_break_fraction = 0.06; p.wrap_whole_fraction = 0.12; p.min_room_fraction = 0.05;
        SlFlow f;
        memset(&f, 0, sizeof f);
        f.c = &c; f.pieces = pc; f.order = order; f.n = 3; f.n_main = 2;
        sl_flow_setup(&f, 300, &p, 5);
        sl_flow(&f);
        f.frags = ARENA_ALLOC(arena, f.n_frags * sizeof *f.frags);
        f.lines = ARENA_ALLOC(arena, f.n_lines * sizeof *f.lines);
        sl_flow(&f);
        sl_check(f.n_frags >= 2 && f.frags[0].sx0 == 0 && f.frags[0].sx1 == 290,
                 "first break falls on the thin column 290", &fails);
        sl_check(f.n_split == 1, "the long piece is the only split", &fails);
        int32_t covered_cols = 0;
        for (size_t i = 0; i < f.n_frags; i++) if (f.frags[i].piece == 0) covered_cols += f.frags[i].sx1 - f.frags[i].sx0;
        sl_check(covered_cols == 1165, "fragments of the long piece partition its columns", &fails);
        int narrow_whole = 0, tray_whole = 0, tray_new_line = 0;
        for (size_t i = 0; i < f.n_frags; i++) {
            if (f.frags[i].piece == 1 && f.frags[i].sx0 == 1200 && f.frags[i].sx1 == 1230 && f.frags[i].dx == 0 &&
                f.frags[i].line == 4) narrow_whole = 1;
            if (f.frags[i].piece == 2 && f.frags[i].sx0 == 1300 && f.frags[i].sx1 == 1500) {
                tray_whole = 1;
                tray_new_line = (size_t)f.frags[i].line == f.tray_first_line && f.frags[i].dx == 0;
            }
        }
        sl_check(narrow_whole, "narrow main piece moves whole to a new line", &fails);
        sl_check(tray_whole && tray_new_line, "tray piece unsplit, on a new line", &fails);
        for (size_t i = 0; i < f.n_frags; i++)
            sl_check(f.frags[i].dx + (f.frags[i].sx1 - f.frags[i].sx0) <= 300, "fragment within the line", &fails);
    }

    /* 3. Aspect fit on a long ribbon of equal pieces (60 x 100 px, 40 apart). */
    {
        const int32_t W = 40 * 100, H = 60;
        uint8_t *cov = ARENA_CALLOC(arena, (size_t)W * H, 1), *tex = ARENA_CALLOC(arena, (size_t)W * H, 1);
        for (int32_t k = 0; k < 40; k++) sl_fill(cov, tex, W, k * 100, k * 100 + 60, 0, 60);
        StripColumns cols = {0};
        StripParams p;
        StripPage pg = {0};
        StripLayout_params_default(&p, 1.0, 1.0);
        p.kibble_area_vox2 = 0.0; p.line_gap_vox = 10.0; p.piece_gap_vox = 10.0; p.margin_vox = 5.0;
        int rc = StripLayout_columns(arena, cov, W, H, &cols) || StripLayout_plan(arena, &cols, &p, &pg);
        if (sl_check(rc == 0, "plan of a long ribbon", &fails)) {
            double ratio = (double)pg.W / (double)pg.H;
            sl_check(pg.n_pieces == 40 && pg.n_tray == 0, "40 main pieces", &fails);
            sl_check(ratio > p.aspect / 1.35 && ratio < p.aspect * 1.35, "page near the target aspect", &fails);
            sl_check(pg.n_lines > 1, "the ribbon wraps", &fails);
        }
    }

    /* 3b. Balance: nine 100 x 50 pieces (gap 10, line gap 10, margin 5).  4
     *    per line (4+4+1) gives 440 x 180 (2.44), 3 per line 330 x 180 (1.83);
     *    at their geometric-mean aspect both fit equally and the balance must
     *    not leave the ninth piece alone on the last line. */
    {
        const int32_t W = 9 * 150, H = 50;
        uint8_t *cov = ARENA_CALLOC(arena, (size_t)W * H, 1), *tex = ARENA_CALLOC(arena, (size_t)W * H, 1);
        for (int32_t k = 0; k < 9; k++) sl_fill(cov, tex, W, k * 150, k * 150 + 100, 0, 50);
        StripColumns cols = {0};
        StripParams p;
        StripPage pg = {0};
        StripLayout_params_default(&p, 1.0, 1.0);
        p.kibble_area_vox2 = 0.0; p.line_gap_vox = 10.0; p.piece_gap_vox = 10.0; p.margin_vox = 5.0;
        p.wrap_whole_fraction = 0.3;                        /* pieces move whole */
        p.aspect = sqrt((440.0 / 180.0) * (330.0 / 180.0));
        int rc = StripLayout_columns(arena, cov, W, H, &cols) || StripLayout_plan(arena, &cols, &p, &pg);
        if (sl_check(rc == 0 && pg.n_pieces == 9, "plan of nine pieces", &fails)) {
            size_t last = pg.n_lines - 1, on_last = 0;
            for (size_t i = 0; i < pg.n_frags; i++) on_last += (size_t)pg.frags[i].line == last;
            sl_check(on_last >= 2, "no lone piece on the last line", &fails);
            sl_check(pg.n_split == 0, "equal pieces are never split", &fails);
        }
    }

    /* 3c. A diagonal band 600 columns wide (row = column / 10, 5 rows thick)
     *     split at a 400-column line: the continuation (rows ~38-65) hangs
     *     from its line's top, and the second line is only as tall as it. */
    {
        const int32_t W = 600, H = 70;
        uint8_t *cov = ARENA_CALLOC(arena, (size_t)W * H, 1), *tex = ARENA_CALLOC(arena, (size_t)W * H, 1);
        for (int32_t x = 0; x < W; x++) sl_fill(cov, tex, W, x, x + 1, x / 10, x / 10 + 5);
        StripColumns cols = {0};
        StripParams p;
        StripLayout_params_default(&p, 1.0, 1.0);
        int rc = StripLayout_columns(arena, cov, W, H, &cols);
        StripPiece band;
        memset(&band, 0, sizeof band);
        band.x0 = 0; band.x1 = W; band.y0 = 0; band.y1 = H;
        size_t order[1] = {0};
        SlFlow f;
        memset(&f, 0, sizeof f);
        f.c = &cols; f.pieces = &band; f.order = order; f.n = 1; f.n_main = 1;
        sl_flow_setup(&f, 400, &p, 5);
        sl_flow(&f);
        f.frags = ARENA_ALLOC(arena, f.n_frags * sizeof *f.frags);
        f.lines = ARENA_ALLOC(arena, f.n_lines * sizeof *f.lines);
        sl_flow(&f);
        if (sl_check(rc == 0 && f.n_frags == 2 && f.n_lines == 2, "diagonal band split once", &fails)) {
            const StripFragment *second = f.frags + 1;
            sl_check(second->sy0 == second->sx0 / 10, "continuation rows start at its own first covered row", &fails);
            sl_check(f.lines[1].h == second->sy1 - second->sy0 && f.lines[1].h < 30,
                     "second line only as tall as the lifted continuation", &fails);
        }
    }

    /* 4. Degenerate inputs, and a single covered pixel (tray only). */
    {
        uint8_t z[16] = {0};
        StripColumns cols = {0};
        StripParams p;
        StripPage pg = {0};
        StripLayout_params_default(&p, 2.0, 2.0);
        int rc = StripLayout_columns(arena, z, 4, 4, &cols);
        sl_check(rc == 0 && StripLayout_plan(arena, &cols, &p, &pg) == -1, "empty coverage is refused", &fails);
        sl_check(StripLayout_columns(arena, NULL, 4, 4, &cols) == -1, "NULL coverage is refused", &fails);
        sl_check(StripLayout_plan(NULL, &cols, &p, &pg) == -1, "NULL arena is refused", &fails);
        uint8_t one[16] = {0}, t1[16] = {0};
        one[5] = 255; t1[5] = 77;
        rc = StripLayout_columns(arena, one, 4, 4, &cols) || StripLayout_plan(arena, &cols, &p, &pg);
        if (sl_check(rc == 0 && pg.n_pieces == 1 && pg.n_tray == 1, "single pixel is one tray piece", &fails)) {
            uint8_t *page = ARENA_CALLOC(arena, (size_t)pg.W * (size_t)pg.H, 1);
            int64_t copied = 0;
            rc = StripLayout_render(&pg, t1, one, 4, 4, page, &copied);
            sl_check(rc == 0 && copied == 1, "single pixel copied", &fails);
        }
    }

    /* 5. Downsample: mean of nonzero pixels per block, 0 for empty blocks. */
    {
        const uint8_t img[4 * 3] = {10, 20, 0, 0,
                                    0, 30, 0, 0,
                                    5, 0, 0, 9};
        uint8_t out[4] = {0};
        int32_t ow = 0, oh = 0;
        int rc = StripLayout_downsample(img, 4, 3, 2, out, &ow, &oh);
        sl_check(rc == 0 && ow == 2 && oh == 2, "downsample size", &fails);
        sl_check(out[0] == 20 && out[1] == 0 && out[2] == 5 && out[3] == 9, "downsample means", &fails);
    }

    Arena_dispose(&arena);
    fprintf(stderr, "[strip_layout selftest] %s (%d failures)\n", fails == 0 ? "ok" : "FAILED", fails);
    return fails;
}
