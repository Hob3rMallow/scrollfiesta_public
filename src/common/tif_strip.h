#ifndef TIF_STRIP_INCLUDED
#define TIF_STRIP_INCLUDED

#include <stdint.h>

/* ============================================================================
 * tif_strip.h -- multi-line "strip" view of a very wide unrolled raster.
 *
 * A scroll unroll is thousands of columns wide but only a few hundred rows
 * tall, so the wide _rawtex.tif is unreadable at a glance. TifStrip_write_gray
 * slices it into ceil(W / target_strip_w) equal-width vertical strips and
 * stacks them top-to-bottom (with `gap` blank rows between) into a readable
 * near-square image -- the same layout scripts/flatten/diag_render.c produced
 * as a manual post-step, now emitted automatically beside every bake so every
 * pipeline stage delivers both <stage>_rawtex.tif and <stage>_rawtex_strip.tif.
 * ==========================================================================*/

/* Write the stacked strip image of the W x H row-major grayscale `img` to
 * `path` (a grayscale TIFF) plus a .png sibling at the same base name.
 * target_strip_w <= 0 -> default (~3000 px); gap < 0 -> default. Returns 0 on
 * success, -1 on failure / unreasonable size (caller may ignore the return:
 * the strip is a debug convenience, never load-bearing). */
int TifStrip_write_gray(const char *path, const uint8_t *img, int W, int H,
                        int target_strip_w, int gap);

/* In-process unit test of the stacking geometry + content. 0 = pass. */
int TifStrip_selftest(void);

#endif
