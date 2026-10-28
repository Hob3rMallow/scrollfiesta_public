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

/* Save flat uint8 volume as multi-page TIFF; automatically uses BigTIFF when
 * classic 32-bit offsets cannot safely hold the raster and directory tables.
 * Returns 0 on success, -1 on failure. */
int TiffIO_save(const char *path,
                const uint8_t *vol, int D, int H, int W);

/* Save a single 2D float image as an uncompressed 32-bit IEEE-float TIFF
 * (1 sample/pixel, MINISBLACK, top-left). Row-major img[H*W].
 * This is the on-disk shape of each tifxyz coordinate map (x.tif/y.tif/z.tif).
 * Returns 0 on success, -1 on failure. */
int TiffIO_save_float2d(const char *path, const float *img, int W, int H);

/* TiffIO_save / TiffIO_save_float2d with every strip deflated (zlib level 1, no predictor; same strips, rows and
 * values).  2026-10-21: a block export's sheet canvas is ~90% filler (-1 / 0), and its 4.6 GB of plain grids held the
 * mechanical D: for ~78 s a block.  Every strip reader (TiffIO_load, TiffIO_load_float2d, TIFFReadScanline /
 * TIFFReadEncodedStrip, OpenCV, tifffile) reads the result unchanged; TiffIO_rows_open declines it (callers fall back
 * to TiffIO_load) and it cannot be memory-mapped. Builds with VES_TIFF_PARALLEL_DEFLATE
 * compress independent strips on their OpenMP workers, with bounded buffers;
 * libtiff calls and output order remain serialized. */
int TiffIO_save_deflate(const char *path, const uint8_t *vol, int D, int H, int W);
int TiffIO_save_float2d_deflate(const char *path, const float *img, int W, int H);
/* Lossless TIFF floating-point predictor (3), followed by deflate. Byte
 * shuffling/differencing preserves every float bit, including NaN payloads.
 * Explicit candidate API: existing deflate callers keep their current format. */
int TiffIO_save_float2d_predicted(const char *path, const float *img, int W, int H);

/* Fresh ABBA storage control for one native u8 or float32 TIFF. Reads pixels
 * once, writes four fresh files, verifies every decoded byte, and reports
 * write/verify times and file sizes. Float: deflate vs predictor+deflate;
 * uint8: plain vs deflate. Never changes the input or replaces an output. */
int TiffIO_storage_probe(const char *input, const char *fresh_directory);

/* Load a single 2D float TIFF (as written by TiffIO_save_float2d).
 * Arena-allocates *out_img of size W*H. Returns 0 on success, -1 on failure. */
int TiffIO_load_float2d(Arena_T arena, const char *path,
                        float **out_img, int *out_W, int *out_H);

/* Bounded random windows from a single float32 TIFF, stripped or tiled,
 * compressed or plain. Each decoded strip/tile must fit cache_bytes; up to
 * eight units are retained, within that total budget. No whole-plane fallback.
 * The opaque reader owns its TIFF handle and cache until close. Reads copy
 * native float bits into contiguous dst[h][w]; negative/out-of-bounds windows
 * and malformed or truncated units fail. Callers serialize use of a reader. */
typedef struct TiffFloatWindow TiffFloatWindow;
TiffFloatWindow *TiffIO_float_window_open(const char *path, size_t cache_bytes, int *W, int *H);
int TiffIO_float_window_read(TiffFloatWindow *reader, int y, int x, int h, int w, float *dst);
void TiffIO_float_window_close(TiffFloatWindow *reader);

/* Load a single 2D signed 32-bit integer TIFF.  This is used by tifxyz turn /
 * topology labels, whose -1 background cannot be represented by TiffIO_load. */
int TiffIO_load_int32_2d(Arena_T arena, const char *path,
                         int32_t **out_img, int *out_W, int *out_H);

/* ----------------------------------------------------------------------------
 * Streaming row access to a single-page UNCOMPRESSED 8-bit 1-sample gray TIFF
 * (the shape TiffIO_save writes: the whole-grid strip rasters are ~610 MB at
 * 4x21x21 scale, far too big to TiffIO_load just to look at a few rows).
 * Open copies the strip layout via libtiff, then reads rows by direct file
 * offset: peak memory = the caller's row buffer. */
typedef struct TiffRowReader {
    FILE     *f;             /* private stdio handle */
    int       W, H;
    int       px_bytes;      /* 1: 8-bit gray; 4: 32-bit IEEE float (a tifxyz coordinate plane) */
    uint32_t  rows_per_strip;
    uint64_t *strip_off;     /* [n_strips] byte offset of each strip (arena) */
    uint32_t  n_strips;
} TiffRowReader;

/* Open the FIRST directory of `path` for random row access. Fails (-1) unless
 * it is uncompressed, 1 sample/pixel and either 8-bit (px_bytes 1) or 32-bit
 * IEEE float (px_bytes 4, the tifxyz planes TiffIO_save_float2d writes) --
 * caller falls back to TiffIO_load. On success the reader owns a FILE handle
 * until _close. */
int  TiffIO_rows_open(Arena_T arena, const char *path, TiffRowReader *rd);

/* Read `w` pixels of row `y` starting at column `x0` into dst[w] (8-bit
 * readers only). Returns 0 on success, -1 on bounds/IO error. */
int  TiffIO_rows_read(TiffRowReader *rd, int y, int x0, int w, uint8_t *dst);

/* The same for a 32-bit float reader: dst[w] floats. */
int  TiffIO_rows_read_f32(TiffRowReader *rd, int y, int x0, int w, float *dst);

void TiffIO_rows_close(TiffRowReader *rd);

/* Native I/O regression fixtures, written beneath an existing directory. */
int TiffIO_selftest(const char *dir);

#endif
