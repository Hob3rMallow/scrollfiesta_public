#ifndef STRIP_LAYOUT_INCLUDED
#define STRIP_LAYOUT_INCLUDED

#include <stdint.h>
#include <stdio.h>
#include "arena.h"

/* ============================================================================
 * strip_layout.h -- piece-aware "strip" pages of a winding-ribbon bake.
 *
 * The winding ribbon (sheet_assemble --winding-order) lays every intact
 * piece in its own column range, in approximate winding (spiral) order and
 * ASM_READING_RIBBON_GAP apart; the bake (obj_bake_raw) rasters it into one
 * very wide sheet_rawtex.tif plus a coverage mask.  A strip page flows that
 * ribbon into lines, like text, so the whole sheet reads on one screen:
 *   - pieces are recovered from the coverage mask as runs of covered columns
 *     (empty runs narrower than half the ribbon gap stay inside a piece), so
 *     the page needs no geometry and is checked against ribbon_layout.json;
 *   - main pieces keep ribbon (spiral) order and flow across line breaks; a
 *     break falls on the piece's thinnest column near the end of the line;
 *   - kibble (covered area under STRIP_KIBBLE_AREA_VOX2) goes to a tray after
 *     the spiral, largest first, never split unless wider than a line;
 *   - every fragment hangs from the top of its line by its own first covered
 *     row (whole pieces touch the ribbon's top row anyway; the continuation
 *     of a split diagonal band is lifted instead of dangling below the line);
 *     a line is as tall as its tallest fragment, and the line width is chosen
 *     so that the page approaches STRIP_PAGE_ASPECT.
 * Pure re-arrangement: each covered source pixel lands on the page exactly
 * once with its value unchanged (StripLayout_render counts them).
 * ==========================================================================*/

typedef struct StripParams {
    double du, dv;              /* vox per source pixel (the bake's step_uv) */
    double merge_gap_vox;       /* empty column runs narrower than this stay inside a piece */
    double kibble_area_vox2;    /* pieces with less covered area go to the tray */
    double aspect;              /* target page width / height */
    double line_gap_vox;        /* blank rows between lines */
    double piece_gap_vox;       /* blank columns between pieces on a line */
    double margin_vox;          /* blank border around the page */
    double tray_break_lines;    /* extra line gaps between the spiral and the tray */
    double wrap_whole_fraction; /* a main piece narrower than this share of a line moves whole to the next line */
    double soft_break_fraction; /* a split falls on the thinnest column within this last share of the line */
    double min_room_fraction;   /* with less room than this share of a line left, start a new line */
} StripParams;

/* Per-column profile of a coverage mask (nonzero = covered). */
typedef struct StripColumns {
    int32_t  W, H;
    int32_t *top;               /* [W] first covered row, H when the column is empty */
    int32_t *bottom;            /* [W] last covered row + 1, 0 when empty */
    int32_t *count;             /* [W] covered pixels */
} StripColumns;

typedef struct StripPiece {
    int32_t x0, x1;             /* source columns [x0, x1): ribbon order is x order */
    int32_t y0, y1;             /* covered source rows [y0, y1) */
    int64_t covered;            /* covered source pixels */
    double  area_vox2;          /* covered * du * dv */
    int32_t tray;               /* 1 = kibble, laid out in the tray */
} StripPiece;

typedef struct StripFragment {
    int32_t piece;              /* index into StripPage.pieces */
    int32_t line;               /* index into StripPage.lines */
    int32_t sx0, sx1;           /* source columns [sx0, sx1) */
    int32_t sy0, sy1;           /* covered source rows of those columns [sy0, sy1) */
    int32_t dx;                 /* page column of sx0, from the left margin */
} StripFragment;

typedef struct StripLine {
    int32_t dy;                 /* page row of the line's top: every fragment's sy0 */
    int32_t h;                  /* rows: the tallest fragment */
    int32_t tray;               /* 1 = a tray line */
} StripLine;

typedef struct StripPage {
    StripPiece    *pieces;      size_t n_pieces;
    StripFragment *frags;       size_t n_frags;
    StripLine     *lines;       size_t n_lines;
    int32_t W, H;               /* page size, px */
    int32_t flow_width;         /* line width, px (between the margins) */
    int32_t margin;             /* px */
    size_t  n_tray;             /* kibble pieces */
    size_t  n_split;            /* main pieces split over more than one line */
    double  main_area_vox2, tray_area_vox2;
} StripPage;

/* Defaults from pipeline_constants.h (STRIP_*) for a bake of step du x dv. */
void StripLayout_params_default(StripParams *p, double du, double dv);

/* Column profile of the W x H row-major coverage mask. 0 = ok. */
int StripLayout_columns(Arena_T arena, const uint8_t *cov, int32_t W, int32_t H,
                        StripColumns *out);

/* Pieces, order, line flow and page geometry. 0 = ok; -1 = no covered pixel
 * or invalid input. */
int StripLayout_plan(Arena_T arena, const StripColumns *cols, const StripParams *p,
                     StripPage *out);

/* Copy every covered source pixel to the zeroed page[pg->W * pg->H]
 * (tex/cov are the W x H source). *out_copied = covered pixels written. */
int StripLayout_render(const StripPage *pg, const uint8_t *tex, const uint8_t *cov,
                       int32_t W, int32_t H, uint8_t *page, int64_t *out_copied);

/* Box-downsample by f (mean of the nonzero pixels of each f x f block, 0 when
 * none) into out[ceil(W/f) * ceil(H/f)]. 0 = ok. */
int StripLayout_downsample(const uint8_t *img, int32_t W, int32_t H, int32_t f,
                           uint8_t *out, int32_t *out_w, int32_t *out_h);

/* Pieces, lines and fragments as JSON members (no enclosing braces): page
 * pixel (x, y) of fragment f on line l shows source pixel
 * (f.sx0 + x - margin - f.dx, f.sy0 + y - l.dy). 0 = ok. */
int StripLayout_write_ledger(FILE *f, const StripPage *pg, const StripParams *p);

/* Synthetic ribbons: segmentation, tray order, split and wrap rules, pixel
 * conservation, aspect fit, degenerate inputs. 0 = pass. */
int StripLayout_selftest(void);

#endif
