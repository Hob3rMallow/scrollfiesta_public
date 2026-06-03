#ifndef TIFF_IO_INCLUDED
#define TIFF_IO_INCLUDED

#include <stdint.h>
#include "arena.h"

/* Load multi-page TIFF as flat uint8 volume.
 * Arena-allocates *out_vol of size D*H*W.
 * Returns 0 on success, -1 on failure.
 * Page 0 = z=0, page 1 = z=1, etc. */
int TiffIO_load(Arena_T arena, const char *path,
                uint8_t **out_vol, int *out_D, int *out_H, int *out_W);

/* Save flat uint8 volume as multi-page TIFF.
 * Returns 0 on success, -1 on failure. */
int TiffIO_save(const char *path,
                const uint8_t *vol, int D, int H, int W);

/* Save a single 2D float image as an uncompressed 32-bit IEEE-float TIFF
 * (1 sample/pixel, MINISBLACK, top-left). Row-major img[H*W].
 * This is the on-disk shape of each tifxyz coordinate map (x.tif/y.tif/z.tif).
 * Returns 0 on success, -1 on failure. */
int TiffIO_save_float2d(const char *path, const float *img, int W, int H);

/* Load a single 2D float TIFF (as written by TiffIO_save_float2d).
 * Arena-allocates *out_img of size W*H. Returns 0 on success, -1 on failure. */
int TiffIO_load_float2d(Arena_T arena, const char *path,
                        float **out_img, int *out_W, int *out_H);

#endif
