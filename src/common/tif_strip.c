/* tif_strip.c -- see tif_strip.h. Stacks a wide unrolled raster into a
 * readable multi-line strip image (ported from scripts/flatten/diag_render.c's
 * strip layout, factored into a reusable, unit-tested core). */
#include "tif_strip.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tiff_io.h"
#include "ves_png.h"

#define TIFSTRIP_DEFAULT_STRIP_W 3000
#define TIFSTRIP_DEFAULT_GAP     6

/* Pure stacker: slice W x H `img` into ceil(W/target_strip_w) equal-width
 * vertical strips stacked top-to-bottom with `gap` blank rows between. Returns
 * a freshly malloc'd row-major grayscale buffer (caller frees) and its dims
 * via *out_w/*out_h, or NULL on bad args / unreasonable size. */
static uint8_t *tifstrip_stack(const uint8_t *img, int W, int H,
                               int target_strip_w, int gap,
                               int *out_w, int *out_h)
{
    if (img == NULL || W <= 0 || H <= 0) return NULL;
    if (target_strip_w <= 0) target_strip_w = TIFSTRIP_DEFAULT_STRIP_W;
    if (gap < 0) gap = TIFSTRIP_DEFAULT_GAP;

    int nstrips = (W + target_strip_w - 1) / target_strip_w;
    if (nstrips < 1) nstrips = 1;
    int sw = (W + nstrips - 1) / nstrips;                 /* strip width */
    int64_t oh64 = (int64_t)nstrips * ((int64_t)H + gap) - gap;
    if (oh64 < H) oh64 = H;
    int64_t npix = (int64_t)sw * oh64;
    if (oh64 > 100000000 || npix > ((int64_t)1 << 30)) return NULL;  /* runaway */
    int oh = (int)oh64;

    uint8_t *out = (uint8_t *)calloc((size_t)sw * (size_t)oh, 1);
    if (out == NULL) return NULL;

    for (int s = 0; s < nstrips; s++) {
        int x0 = s * sw, oy0 = s * (H + gap);
        for (int y = 0; y < H; y++) {
            for (int x = 0; x < sw; x++) {
                int sx = x0 + x;
                if (sx >= W) continue;                    /* ragged last strip */
                out[(size_t)(oy0 + y) * (size_t)sw + (size_t)x]
                    = img[(size_t)y * (size_t)W + (size_t)sx];
            }
        }
    }
    *out_w = sw; *out_h = oh;
    return out;
}

int TifStrip_write_gray(const char *path, const uint8_t *img, int W, int H,
                        int target_strip_w, int gap)
{
    int sw = 0, oh = 0, rc = -1;
    uint8_t *out = NULL;
    if (path == NULL) return -1;

    out = tifstrip_stack(img, W, H, target_strip_w, gap, &sw, &oh);
    if (out == NULL) return -1;

    rc = TiffIO_save(path, out, 1, oh, sw);

    /* .png sibling at the same base (<path> with .tif -> .png), matching the
     * wide bake's grayscale PNG so the strip is viewable without a TIFF reader. */
    {
        char pngpath[2600];
        size_t plen = strlen(path);
        snprintf(pngpath, sizeof pngpath, "%s", path);
        if (plen > 4 && strcmp(pngpath + plen - 4, ".tif") == 0)
            snprintf(pngpath + plen - 4, sizeof pngpath - (plen - 4), ".png");
        else
            snprintf(pngpath + plen, sizeof pngpath - plen, ".png");
        VesPng_write_gray(pngpath, out, sw, oh);
    }

    free(out);
    return rc;
}

int TifStrip_selftest(void)
{
    int fails = 0;

    /* 10 x 3 image, pixel value = x. target_strip_w = 4, gap = 2:
     *   nstrips = ceil(10/4) = 3, sw = ceil(10/3) = 4, oh = 3*(3+2) - 2 = 13. */
    {
        uint8_t img[30];
        for (int y = 0; y < 3; y++)
            for (int x = 0; x < 10; x++)
                img[y * 10 + x] = (uint8_t)x;
        int w = 0, h = 0;
        uint8_t *out = tifstrip_stack(img, 10, 3, 4, 2, &w, &h);
        if (out == NULL) {
            fprintf(stderr, "[tif_strip selftest]   FAIL: stack returned NULL\n");
            fails++;
        } else {
            if (w != 4 || h != 13) {
                fprintf(stderr, "[tif_strip selftest]   FAIL: dims %dx%d (exp 4x13)\n", w, h);
                fails++;
            }
            /* strip 0 row 0 = cols 0..3 */
            if (!(out[0]==0 && out[1]==1 && out[2]==2 && out[3]==3)) {
                fprintf(stderr, "[tif_strip selftest]   FAIL: strip0 content\n"); fails++;
            }
            /* strip 1 starts at oy = 1*(3+2) = 5; row 0 = cols 4..7 */
            if (!(out[5*4+0]==4 && out[5*4+1]==5 && out[5*4+2]==6 && out[5*4+3]==7)) {
                fprintf(stderr, "[tif_strip selftest]   FAIL: strip1 content\n"); fails++;
            }
            /* strip 2 starts at oy = 10; row 0 = cols 8,9 then ragged (0,0) */
            if (!(out[10*4+0]==8 && out[10*4+1]==9 && out[10*4+2]==0 && out[10*4+3]==0)) {
                fprintf(stderr, "[tif_strip selftest]   FAIL: strip2 ragged tail\n"); fails++;
            }
            /* gap rows (oy 3,4) between strip 0 and strip 1 are blank */
            if (!(out[3*4+0]==0 && out[4*4+3]==0)) {
                fprintf(stderr, "[tif_strip selftest]   FAIL: gap not blank\n"); fails++;
            }
            free(out);
        }
    }

    /* single strip when the image is narrower than the target width */
    {
        uint8_t img[12];
        for (int i = 0; i < 12; i++) img[i] = (uint8_t)i;
        int w = 0, h = 0;
        uint8_t *out = tifstrip_stack(img, 4, 3, 3000, 6, &w, &h);
        if (out == NULL || w != 4 || h != 3) {
            fprintf(stderr, "[tif_strip selftest]   FAIL: single-strip dims %dx%d (exp 4x3)\n",
                    out ? w : -1, out ? h : -1);
            fails++;
        }
        if (out) free(out);
    }

    /* bad args -> NULL */
    {
        int w = 0, h = 0;
        if (tifstrip_stack(NULL, 4, 3, 0, 0, &w, &h) != NULL) {
            fprintf(stderr, "[tif_strip selftest]   FAIL: NULL img not rejected\n"); fails++;
        }
    }

    fprintf(stderr, "[tif_strip selftest] %s\n", fails == 0 ? "ok" : "FAILED");
    return fails;
}
