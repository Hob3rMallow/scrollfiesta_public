#ifndef STRIP_PREVIEW_INCLUDED
#define STRIP_PREVIEW_INCLUDED

#include <stdint.h>
#include <stddef.h>
#include "../common/arena.h"

/* ============================================================================
 * strip_preview.h -- human-viewable renditions of the whole-grid strip TIFs
 * (1.19M x 511 px at 4x21x21 -- no image viewer opens that, and the .png
 * sibling is 114 MB). Everything streams rows via TiffIO_rows_* so peak
 * memory is one input row + the output image, never the 610 MB strip.
 *
 *   overview -- integer box-downsample (u by ds_u, v by ds_v), stacked into
 *               K rows so the PNG stays under max_w px wide.
 *   crop     -- full-res window [x0,x1) x [y0,y1).
 *   windows  -- K full-res +/-win column windows stacked vertically
 *               (the seam-gutter review instrument).
 *
 * Reduction semantics: gray = box MEAN with zeros included (empty stays
 * navy); class = MAX-PRIORITY reduce so one defect pixel survives any
 * downsample factor. Palettes mirror scripts/flatten/diag_render.c.
 * ==========================================================================*/

typedef struct StripPreviewOpts {
    int ds_u, ds_v;        /* integer box-filter factors, >=1 (def 16, 1) */
    int rows;              /* stacked rows; 0 = auto so width <= max_w */
    int max_w;             /* output width cap, px (def 20000) */
    int gap;               /* blank rows between stacks (def 8) */
    const char *mode;      /* "gray" | "class" | "coverage" (def "gray") */
    const char *overlay;   /* optional diagclass tif: tint multi-cover red
                              over the gray texture (gray mode only) */
    int verbose;           /* log geometry to stderr (def 1) */
} StripPreviewOpts;

void StripPreviewOpts_default(StripPreviewOpts *o);

/* Downsampled whole-strip overview -> RGB PNG. */
int StripPreview_overview(Arena_T arena, const char *tif, const char *png,
                          const StripPreviewOpts *o);

/* Full-res window [x0,x1) x [y0,y1) -> RGB PNG (y1<=0 means full height).
 * Width capped at 30000 px after clamping. */
int StripPreview_crop(Arena_T arena, const char *tif, const char *png,
                      long x0, long x1, int y0, int y1,
                      const StripPreviewOpts *o);

/* ncols full-res windows of +/-win columns centred on cols[], stacked
 * vertically with o->gap blank rows -> RGB PNG (seam zooms). */
int StripPreview_windows(Arena_T arena, const char *tif, const char *png,
                         const long *cols, size_t ncols, int win,
                         const StripPreviewOpts *o);

/* In-process unit tests (synthetic TIFs under output/_selftest_strip_preview/;
 * also exercises TiffIO_rows_*). Returns 0 if all pass. */
int StripPreview_selftest(void);

#endif
