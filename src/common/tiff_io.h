#ifndef TIFF_IO_INCLUDED
#define TIFF_IO_INCLUDED

#include <stdint.h>
#include <stdio.h>
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

/* ----------------------------------------------------------------------------
 * Streaming row access to a single-page UNCOMPRESSED 8-bit 1-sample gray TIFF
 * (the shape TiffIO_save writes: the whole-grid strip rasters are ~610 MB at
 * 4x21x21 scale, far too big to TiffIO_load just to look at a few rows).
 * Open copies the strip layout via libtiff, then reads rows by direct file
 * offset: peak memory = the caller's row buffer. */
typedef struct TiffRowReader {
    FILE     *f;             /* private stdio handle */
    int       W, H;
    uint32_t  rows_per_strip;
    uint64_t *strip_off;     /* [n_strips] byte offset of each strip (arena) */
    uint32_t  n_strips;
} TiffRowReader;

/* Open the FIRST directory of `path` for random row access. Fails (-1) unless
 * it is uncompressed, 8-bit, 1 sample/pixel -- caller falls back to
 * TiffIO_load. On success the reader owns a FILE handle until _close. */
int  TiffIO_rows_open(Arena_T arena, const char *path, TiffRowReader *rd);

/* Read `w` pixels of row `y` starting at column `x0` into dst[w].
 * Returns 0 on success, -1 on bounds/IO error. */
int  TiffIO_rows_read(TiffRowReader *rd, int y, int x0, int w, uint8_t *dst);

void TiffIO_rows_close(TiffRowReader *rd);

#endif
