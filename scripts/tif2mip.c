/*
 * tif2mip.c — Maximum Intensity Projection visualization tool
 *
 * Loads 3 multi-page TIFFs (raw CT, GT labels, predictions), computes
 * MIP along Z/Y/X axes, composites a 3×3 grid (rows: CT/GT/Pred,
 * cols: Z-proj/Y-proj/X-proj), and writes PNG.
 *
 * For binary masks (GT/pred), MIP = "any foreground along this ray".
 * For raw CT, MIP = brightest voxel along each ray.
 *
 * Usage: tif2mip <raw.tif> <gt.tif> <pred.tif> <output.png> [cube_id] [--diff] [--contour]
 *
 * Options:
 *   --diff      Add row: GT vs Pred MIP diff (green=TP, red=FP, blue=FN)
 *   --contour   Add row: Pred MIP boundary contour on CT MIP
 *
 * Build: gcc -O2 -o tif2mip scripts/tif2mip.c -ltiff -lpng16 -lm
 *
 * No dependencies on pipeline modules. Only needs libtiff + libpng + libm.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <tiffio.h>
#include <png.h>

/* ================================================================
 * Embedded 8x8 bitmap font (ASCII 32-127, 96 glyphs)
 * Each glyph is 8 bytes, one per row, MSB-first.
 * ================================================================ */

static const uint8_t font8x8[96][8] = {
    {0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00},
    {0x18,0x3C,0x3C,0x18,0x18,0x00,0x18,0x00},
    {0x6C,0x6C,0x24,0x00,0x00,0x00,0x00,0x00},
    {0x6C,0xFE,0x6C,0x6C,0xFE,0x6C,0x00,0x00},
    {0x18,0x7E,0xC0,0x7C,0x06,0xFC,0x18,0x00},
    {0x00,0xC6,0xCC,0x18,0x30,0x66,0xC6,0x00},
    {0x38,0x6C,0x38,0x76,0xDC,0xCC,0x76,0x00},
    {0x18,0x18,0x30,0x00,0x00,0x00,0x00,0x00},
    {0x0C,0x18,0x30,0x30,0x30,0x18,0x0C,0x00},
    {0x30,0x18,0x0C,0x0C,0x0C,0x18,0x30,0x00},
    {0x00,0x66,0x3C,0xFF,0x3C,0x66,0x00,0x00},
    {0x00,0x18,0x18,0x7E,0x18,0x18,0x00,0x00},
    {0x00,0x00,0x00,0x00,0x00,0x18,0x18,0x30},
    {0x00,0x00,0x00,0x7E,0x00,0x00,0x00,0x00},
    {0x00,0x00,0x00,0x00,0x00,0x18,0x18,0x00},
    {0x06,0x0C,0x18,0x30,0x60,0xC0,0x80,0x00},
    {0x7C,0xC6,0xCE,0xDE,0xF6,0xE6,0x7C,0x00},
    {0x18,0x38,0x78,0x18,0x18,0x18,0x7E,0x00},
    {0x7C,0xC6,0x06,0x1C,0x30,0x66,0xFE,0x00},
    {0x7C,0xC6,0x06,0x3C,0x06,0xC6,0x7C,0x00},
    {0x1C,0x3C,0x6C,0xCC,0xFE,0x0C,0x1E,0x00},
    {0xFE,0xC0,0xFC,0x06,0x06,0xC6,0x7C,0x00},
    {0x38,0x60,0xC0,0xFC,0xC6,0xC6,0x7C,0x00},
    {0xFE,0xC6,0x0C,0x18,0x30,0x30,0x30,0x00},
    {0x7C,0xC6,0xC6,0x7C,0xC6,0xC6,0x7C,0x00},
    {0x7C,0xC6,0xC6,0x7E,0x06,0x0C,0x78,0x00},
    {0x00,0x18,0x18,0x00,0x00,0x18,0x18,0x00},
    {0x00,0x18,0x18,0x00,0x00,0x18,0x18,0x30},
    {0x06,0x0C,0x18,0x30,0x18,0x0C,0x06,0x00},
    {0x00,0x00,0x7E,0x00,0x7E,0x00,0x00,0x00},
    {0x60,0x30,0x18,0x0C,0x18,0x30,0x60,0x00},
    {0x7C,0xC6,0x0C,0x18,0x18,0x00,0x18,0x00},
    {0x7C,0xC6,0xDE,0xDE,0xDE,0xC0,0x78,0x00},
    {0x38,0x6C,0xC6,0xC6,0xFE,0xC6,0xC6,0x00},
    {0xFC,0x66,0x66,0x7C,0x66,0x66,0xFC,0x00},
    {0x3C,0x66,0xC0,0xC0,0xC0,0x66,0x3C,0x00},
    {0xF8,0x6C,0x66,0x66,0x66,0x6C,0xF8,0x00},
    {0xFE,0x62,0x68,0x78,0x68,0x62,0xFE,0x00},
    {0xFE,0x62,0x68,0x78,0x68,0x60,0xF0,0x00},
    {0x3C,0x66,0xC0,0xC0,0xCE,0x66,0x3E,0x00},
    {0xC6,0xC6,0xC6,0xFE,0xC6,0xC6,0xC6,0x00},
    {0x3C,0x18,0x18,0x18,0x18,0x18,0x3C,0x00},
    {0x1E,0x0C,0x0C,0x0C,0xCC,0xCC,0x78,0x00},
    {0xE6,0x66,0x6C,0x78,0x6C,0x66,0xE6,0x00},
    {0xF0,0x60,0x60,0x60,0x62,0x66,0xFE,0x00},
    {0xC6,0xEE,0xFE,0xFE,0xD6,0xC6,0xC6,0x00},
    {0xC6,0xE6,0xF6,0xDE,0xCE,0xC6,0xC6,0x00},
    {0x7C,0xC6,0xC6,0xC6,0xC6,0xC6,0x7C,0x00},
    {0xFC,0x66,0x66,0x7C,0x60,0x60,0xF0,0x00},
    {0x7C,0xC6,0xC6,0xC6,0xD6,0xDE,0x7C,0x06},
    {0xFC,0x66,0x66,0x7C,0x6C,0x66,0xE6,0x00},
    {0x7C,0xC6,0xE0,0x7C,0x0E,0xC6,0x7C,0x00},
    {0x7E,0x5A,0x18,0x18,0x18,0x18,0x3C,0x00},
    {0xC6,0xC6,0xC6,0xC6,0xC6,0xC6,0x7C,0x00},
    {0xC6,0xC6,0xC6,0xC6,0x6C,0x38,0x10,0x00},
    {0xC6,0xC6,0xD6,0xFE,0xFE,0xEE,0xC6,0x00},
    {0xC6,0x6C,0x38,0x38,0x6C,0xC6,0xC6,0x00},
    {0x66,0x66,0x66,0x3C,0x18,0x18,0x3C,0x00},
    {0xFE,0xC6,0x8C,0x18,0x32,0x66,0xFE,0x00},
    {0x3C,0x30,0x30,0x30,0x30,0x30,0x3C,0x00},
    {0xC0,0x60,0x30,0x18,0x0C,0x06,0x02,0x00},
    {0x3C,0x0C,0x0C,0x0C,0x0C,0x0C,0x3C,0x00},
    {0x10,0x38,0x6C,0xC6,0x00,0x00,0x00,0x00},
    {0x00,0x00,0x00,0x00,0x00,0x00,0x00,0xFF},
    {0x30,0x18,0x0C,0x00,0x00,0x00,0x00,0x00},
    {0x00,0x00,0x78,0x0C,0x7C,0xCC,0x76,0x00},
    {0xE0,0x60,0x60,0x7C,0x66,0x66,0xDC,0x00},
    {0x00,0x00,0x7C,0xC6,0xC0,0xC6,0x7C,0x00},
    {0x1C,0x0C,0x0C,0x7C,0xCC,0xCC,0x76,0x00},
    {0x00,0x00,0x7C,0xC6,0xFE,0xC0,0x7C,0x00},
    {0x38,0x6C,0x60,0xF0,0x60,0x60,0xF0,0x00},
    {0x00,0x00,0x76,0xCC,0xCC,0x7C,0x0C,0xF8},
    {0xE0,0x60,0x6C,0x76,0x66,0x66,0xE6,0x00},
    {0x18,0x00,0x38,0x18,0x18,0x18,0x3C,0x00},
    {0x06,0x00,0x06,0x06,0x06,0x66,0x66,0x3C},
    {0xE0,0x60,0x66,0x6C,0x78,0x6C,0xE6,0x00},
    {0x38,0x18,0x18,0x18,0x18,0x18,0x3C,0x00},
    {0x00,0x00,0xEC,0xFE,0xD6,0xD6,0xD6,0x00},
    {0x00,0x00,0xDC,0x66,0x66,0x66,0x66,0x00},
    {0x00,0x00,0x7C,0xC6,0xC6,0xC6,0x7C,0x00},
    {0x00,0x00,0xDC,0x66,0x66,0x7C,0x60,0xF0},
    {0x00,0x00,0x76,0xCC,0xCC,0x7C,0x0C,0x1E},
    {0x00,0x00,0xDC,0x76,0x66,0x60,0xF0,0x00},
    {0x00,0x00,0x7C,0xC0,0x7C,0x06,0xFC,0x00},
    {0x10,0x30,0x7C,0x30,0x30,0x34,0x18,0x00},
    {0x00,0x00,0xCC,0xCC,0xCC,0xCC,0x76,0x00},
    {0x00,0x00,0xC6,0xC6,0x6C,0x38,0x10,0x00},
    {0x00,0x00,0xC6,0xD6,0xD6,0xFE,0x6C,0x00},
    {0x00,0x00,0xC6,0x6C,0x38,0x6C,0xC6,0x00},
    {0x00,0x00,0xC6,0xC6,0x6C,0x38,0x30,0x60},
    {0x00,0x00,0xFE,0x8C,0x18,0x32,0xFE,0x00},
    {0x0E,0x18,0x18,0x70,0x18,0x18,0x0E,0x00},
    {0x18,0x18,0x18,0x00,0x18,0x18,0x18,0x00},
    {0x70,0x18,0x18,0x0E,0x18,0x18,0x70,0x00},
    {0x76,0xDC,0x00,0x00,0x00,0x00,0x00,0x00},
    {0x00,0x10,0x38,0x6C,0xC6,0xFE,0x00,0x00},
};

#define FONT_SCALE 2
#define GLYPH_W (8 * FONT_SCALE)
#define GLYPH_H (8 * FONT_SCALE)

static void draw_char(uint8_t *canvas, int cw, int ch,
                      int px, int py, char c,
                      uint8_t r, uint8_t g, uint8_t b)
{
    if (c < 32 || c > 127) c = '?';
    const uint8_t *glyph = font8x8[c - 32];
    for (int gy = 0; gy < 8; gy++) {
        uint8_t row = glyph[gy];
        for (int gx = 0; gx < 8; gx++) {
            if (row & (0x80 >> gx)) {
                for (int sy = 0; sy < FONT_SCALE; sy++) {
                    for (int sx = 0; sx < FONT_SCALE; sx++) {
                        int cx_pos = px + gx * FONT_SCALE + sx;
                        int cy_pos = py + gy * FONT_SCALE + sy;
                        if (cx_pos >= 0 && cx_pos < cw && cy_pos >= 0 && cy_pos < ch) {
                            size_t idx = ((size_t)cy_pos * (size_t)cw + (size_t)cx_pos) * 3;
                            canvas[idx + 0] = r;
                            canvas[idx + 1] = g;
                            canvas[idx + 2] = b;
                        }
                    }
                }
            }
        }
    }
}

static void draw_string(uint8_t *canvas, int cw, int ch,
                        int px, int py, const char *str,
                        uint8_t r, uint8_t g, uint8_t b)
{
    for (int i = 0; str[i]; i++)
        draw_char(canvas, cw, ch, px + i * GLYPH_W, py, str[i], r, g, b);
}

static void draw_string_centered(uint8_t *canvas, int cw, int ch,
                                 int py, const char *str,
                                 uint8_t r, uint8_t g, uint8_t b)
{
    int len = (int)strlen(str);
    int px = (cw - len * GLYPH_W) / 2;
    draw_string(canvas, cw, ch, px, py, str, r, g, b);
}

/* ================================================================
 * TIFF Loading
 * ================================================================ */

static uint8_t *load_tiff(const char *path, int *out_D, int *out_H, int *out_W)
{
    TIFF *tif = TIFFOpen(path, "r");
    if (!tif) {
        fprintf(stderr, "ERROR: cannot open TIFF: %s\n", path);
        return NULL;
    }

    uint32_t w = 0, h = 0;
    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);

    int depth = 0;
    do { depth++; } while (TIFFReadDirectory(tif));

    size_t vol_size = (size_t)depth * (size_t)h * (size_t)w;
    uint8_t *vol = (uint8_t *)calloc(vol_size, 1);
    if (!vol) {
        fprintf(stderr, "ERROR: cannot allocate %zu bytes\n", vol_size);
        TIFFClose(tif);
        return NULL;
    }

    TIFFSetDirectory(tif, 0);
    uint8_t *row_buf = (uint8_t *)malloc((size_t)TIFFScanlineSize(tif));
    if (!row_buf) { free(vol); TIFFClose(tif); return NULL; }

    for (int z = 0; z < depth; z++) {
        for (uint32_t y = 0; y < h; y++) {
            TIFFReadScanline(tif, row_buf, y, 0);
            memcpy(vol + (size_t)z * h * w + (size_t)y * w, row_buf, w);
        }
        if (z < depth - 1) TIFFReadDirectory(tif);
    }

    free(row_buf);
    TIFFClose(tif);
    *out_D = depth;
    *out_H = (int)h;
    *out_W = (int)w;
    return vol;
}

/* ================================================================
 * PNG Writing
 * ================================================================ */

static int write_png(const char *path, const uint8_t *canvas, int w, int h)
{
    FILE *fp = fopen(path, "wb");
    if (!fp) { fprintf(stderr, "ERROR: cannot open: %s\n", path); return -1; }

    png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!png) { fclose(fp); return -1; }

    png_infop info = png_create_info_struct(png);
    if (!info) { png_destroy_write_struct(&png, NULL); fclose(fp); return -1; }

    if (setjmp(png_jmpbuf(png))) {
        png_destroy_write_struct(&png, &info); fclose(fp); return -1;
    }

    png_init_io(png, fp);
    png_set_compression_level(png, 6);
    png_set_IHDR(png, info, (png_uint_32)w, (png_uint_32)h,
                 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
    png_write_info(png, info);

    for (int y = 0; y < h; y++)
        png_write_row(png, (const png_bytep)(canvas + (size_t)y * (size_t)w * 3));

    png_write_end(png, NULL);
    png_destroy_write_struct(&png, &info);
    fclose(fp);
    return 0;
}

/* ================================================================
 * MIP computation
 * ================================================================ */

/* Z-projection: for each (y,x), max over all z. Output: H x W */
static void mip_z(const uint8_t *vol, int D, int H, int W, uint8_t *out)
{
    memset(out, 0, (size_t)H * (size_t)W);
    for (int z = 0; z < D; z++) {
        const uint8_t *slice = vol + (size_t)z * (size_t)H * (size_t)W;
        for (int i = 0; i < H * W; i++) {
            if (slice[i] > out[i]) out[i] = slice[i];
        }
    }
}

/* Y-projection: for each (z,x), max over all y. Output: D x W */
static void mip_y(const uint8_t *vol, int D, int H, int W, uint8_t *out)
{
    memset(out, 0, (size_t)D * (size_t)W);
    for (int z = 0; z < D; z++) {
        for (int y = 0; y < H; y++) {
            const uint8_t *row = vol + (size_t)z * (size_t)H * (size_t)W
                                     + (size_t)y * (size_t)W;
            uint8_t *dst = out + (size_t)z * (size_t)W;
            for (int x = 0; x < W; x++) {
                if (row[x] > dst[x]) dst[x] = row[x];
            }
        }
    }
}

/* X-projection: for each (z,y), max over all x. Output: D x H */
static void mip_x(const uint8_t *vol, int D, int H, int W, uint8_t *out)
{
    memset(out, 0, (size_t)D * (size_t)H);
    for (int z = 0; z < D; z++) {
        for (int y = 0; y < H; y++) {
            const uint8_t *row = vol + (size_t)z * (size_t)H * (size_t)W
                                     + (size_t)y * (size_t)W;
            uint8_t mx = 0;
            for (int x = 0; x < W; x++) {
                if (row[x] > mx) mx = row[x];
            }
            out[(size_t)z * (size_t)H + (size_t)y] = mx;
        }
    }
}

/* ================================================================
 * Helpers
 * ================================================================ */

static inline uint8_t clamp_u8(float v)
{
    if (v < 0.0f) return 0;
    if (v > 255.0f) return 255;
    return (uint8_t)(v + 0.5f);
}

static inline int is_boundary(const uint8_t *img, int rows, int cols, int r, int c)
{
    if (r == 0 || r == rows - 1 || c == 0 || c == cols - 1) return 1;
    if (!img[(r - 1) * cols + c]) return 1;
    if (!img[(r + 1) * cols + c]) return 1;
    if (!img[r * cols + (c - 1)]) return 1;
    if (!img[r * cols + (c + 1)]) return 1;
    return 0;
}

/* Row types */
enum {
    ROW_CT = 0,
    ROW_GT,
    ROW_PRED,
    ROW_DIFF,
    ROW_CONTOUR,
    ROW_MAX
};

/* ================================================================
 * Main
 * ================================================================ */

int main(int argc, char **argv)
{
    if (argc < 5) {
        fprintf(stderr,
            "Usage: %s <raw.tif> <gt.tif> <pred.tif> <output.png> [cube_id] [--diff] [--contour]\n",
            argv[0]);
        return 1;
    }

    const char *raw_path  = argv[1];
    const char *gt_path   = argv[2];
    const char *pred_path = argv[3];
    const char *out_path  = argv[4];
    const char *cube_id   = "unknown";
    int opt_diff = 0;
    int opt_contour = 0;

    for (int i = 5; i < argc; i++) {
        if (strcmp(argv[i], "--diff") == 0)         opt_diff = 1;
        else if (strcmp(argv[i], "--contour") == 0)  opt_contour = 1;
        else if (argv[i][0] != '-')                  cube_id = argv[i];
        else { fprintf(stderr, "Unknown option: %s\n", argv[i]); return 1; }
    }

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    /* Load volumes */
    int rD, rH, rW, gD, gH, gW, pD, pH, pW;
    uint8_t *raw  = load_tiff(raw_path,  &rD, &rH, &rW);  if (!raw)  return 1;
    uint8_t *gt   = load_tiff(gt_path,   &gD, &gH, &gW);  if (!gt)   { free(raw); return 1; }
    uint8_t *pred = load_tiff(pred_path, &pD, &pH, &pW);   if (!pred) { free(raw); free(gt); return 1; }

    if (rD != gD || rD != pD || rH != gH || rH != pH || rW != gW || rW != pW) {
        fprintf(stderr, "ERROR: dimension mismatch\n");
        free(raw); free(gt); free(pred); return 1;
    }

    int D = rD, H = rH, W = rW;
    fprintf(stderr, "  loaded 3 volumes: %d x %d x %d\n", D, H, W);

    struct timespec t_load;
    clock_gettime(CLOCK_MONOTONIC, &t_load);
    double load_sec = (t_load.tv_sec - t0.tv_sec) + (t_load.tv_nsec - t0.tv_nsec) * 1e-9;
    fprintf(stderr, "  load time: %.2fs\n", load_sec);

    /* Compute 9 MIPs: 3 volumes × 3 axes */
    /* Z-proj: H×W, Y-proj: D×W, X-proj: D×H */
    uint8_t *raw_mip_z  = (uint8_t *)malloc((size_t)H * (size_t)W);
    uint8_t *raw_mip_y  = (uint8_t *)malloc((size_t)D * (size_t)W);
    uint8_t *raw_mip_x  = (uint8_t *)malloc((size_t)D * (size_t)H);
    uint8_t *gt_mip_z   = (uint8_t *)malloc((size_t)H * (size_t)W);
    uint8_t *gt_mip_y   = (uint8_t *)malloc((size_t)D * (size_t)W);
    uint8_t *gt_mip_x   = (uint8_t *)malloc((size_t)D * (size_t)H);
    uint8_t *pred_mip_z = (uint8_t *)malloc((size_t)H * (size_t)W);
    uint8_t *pred_mip_y = (uint8_t *)malloc((size_t)D * (size_t)W);
    uint8_t *pred_mip_x = (uint8_t *)malloc((size_t)D * (size_t)H);

    mip_z(raw,  D, H, W, raw_mip_z);
    mip_y(raw,  D, H, W, raw_mip_y);
    mip_x(raw,  D, H, W, raw_mip_x);
    mip_z(gt,   D, H, W, gt_mip_z);
    mip_y(gt,   D, H, W, gt_mip_y);
    mip_x(gt,   D, H, W, gt_mip_x);
    mip_z(pred, D, H, W, pred_mip_z);
    mip_y(pred, D, H, W, pred_mip_y);
    mip_x(pred, D, H, W, pred_mip_x);

    /* Free full volumes — only need MIPs from here */
    free(raw); free(gt); free(pred);

    struct timespec t_mip;
    clock_gettime(CLOCK_MONOTONIC, &t_mip);
    double mip_sec = (t_mip.tv_sec - t_load.tv_sec) + (t_mip.tv_nsec - t_load.tv_nsec) * 1e-9;
    fprintf(stderr, "  MIP computation: %.2fs\n", mip_sec);

    /* MIP image arrays per projection axis:
     * col 0 (Z-proj): H rows × W cols
     * col 1 (Y-proj): D rows × W cols
     * col 2 (X-proj): D rows × H cols
     */
    int col_w[3]    = { W, W, H };
    int col_h[3]    = { H, D, D };
    uint8_t *ct_mips[3]   = { raw_mip_z,  raw_mip_y,  raw_mip_x  };
    uint8_t *gt_mips[3]   = { gt_mip_z,   gt_mip_y,   gt_mip_x   };
    uint8_t *pred_mips[3] = { pred_mip_z, pred_mip_y, pred_mip_x };

    static const char *col_names[3] = { "Z-proj", "Y-proj", "X-proj" };

    /* Build row list */
    int row_types[ROW_MAX];
    int n_rows = 0;
    row_types[n_rows++] = ROW_CT;
    row_types[n_rows++] = ROW_GT;
    row_types[n_rows++] = ROW_PRED;
    if (opt_diff)    row_types[n_rows++] = ROW_DIFF;
    if (opt_contour) row_types[n_rows++] = ROW_CONTOUR;

    static const char *row_name_table[] = {
        [ROW_CT]      = "CT",
        [ROW_GT]      = "GT",
        [ROW_PRED]    = "Pred",
        [ROW_DIFF]    = "Diff",
        [ROW_CONTOUR] = "Contour",
    };

    /* Canvas: find max cell dimensions per column */
    int pad = 6;
    int title_h = GLYPH_H + 4;
    int label_h = GLYPH_H + 2;

    /* Each column has its own width; all rows in a column share it.
     * Row height varies per column too, but we use max height per row. */
    int max_cell_h = 0;
    int total_cell_w = 0;
    for (int c = 0; c < 3; c++) {
        if (col_h[c] > max_cell_h) max_cell_h = col_h[c];
        total_cell_w += col_w[c];
    }

    int canvas_w = pad + 3 * pad + total_cell_w;
    int canvas_h = title_h + n_rows * (label_h + max_cell_h + pad) + pad;

    size_t canvas_size = (size_t)canvas_w * (size_t)canvas_h * 3;
    uint8_t *canvas = (uint8_t *)malloc(canvas_size);
    if (!canvas) {
        fprintf(stderr, "ERROR: cannot allocate canvas\n");
        return 1;
    }
    memset(canvas, 0x30, canvas_size);

    /* Composite */
    for (int c = 0; c < 3; c++) {
        /* Compute cell_x: pad + sum of previous column widths + previous pads */
        int cell_x = pad;
        for (int cc = 0; cc < c; cc++) cell_x += col_w[cc] + pad;

        int img_h = col_h[c];
        int img_w = col_w[c];
        const uint8_t *ct_img   = ct_mips[c];
        const uint8_t *gt_img   = gt_mips[c];
        const uint8_t *pred_img = pred_mips[c];

        for (int row = 0; row < n_rows; row++) {
            int rtype = row_types[row];
            int cell_y = title_h + row * (label_h + max_cell_h + pad) + label_h;
            /* Center vertically if this column is shorter than max */
            int y_off = (max_cell_h - img_h) / 2;

            /* Label */
            char label[64];
            snprintf(label, sizeof(label), "%s %s", row_name_table[rtype], col_names[c]);
            int label_y = title_h + row * (label_h + max_cell_h + pad);
            draw_string(canvas, canvas_w, canvas_h,
                       cell_x + 2, label_y + 1, label, 0xCC, 0xCC, 0xCC);

            /* Pixels */
            for (int r = 0; r < img_h; r++) {
                for (int x = 0; x < img_w; x++) {
                    size_t mi = (size_t)r * (size_t)img_w + (size_t)x;
                    float ct = (float)ct_img[mi];
                    uint8_t rv, gv, bv;

                    switch (rtype) {
                    case ROW_CT:
                        rv = gv = bv = (uint8_t)ct;
                        break;

                    case ROW_GT:
                        if (gt_img[mi] >= 1) {
                            rv = clamp_u8(ct * 0.3f);
                            gv = clamp_u8(ct * 0.4f + 255.0f * 0.6f);
                            bv = clamp_u8(ct * 0.3f);
                        } else {
                            rv = gv = bv = (uint8_t)ct;
                        }
                        break;

                    case ROW_PRED:
                        if (pred_img[mi] >= 1) {
                            rv = clamp_u8(ct * 0.4f + 255.0f * 0.6f);
                            gv = clamp_u8(ct * 0.3f);
                            bv = clamp_u8(ct * 0.3f);
                        } else {
                            rv = gv = bv = (uint8_t)ct;
                        }
                        break;

                    case ROW_DIFF: {
                        int gf = (gt_img[mi] >= 1);
                        int pf = (pred_img[mi] >= 1);
                        if (gf && pf) {
                            rv = clamp_u8(ct * 0.3f);
                            gv = clamp_u8(ct * 0.4f + 255.0f * 0.6f);
                            bv = clamp_u8(ct * 0.3f);
                        } else if (pf && !gf) {
                            rv = clamp_u8(ct * 0.4f + 255.0f * 0.6f);
                            gv = clamp_u8(ct * 0.3f);
                            bv = clamp_u8(ct * 0.3f);
                        } else if (gf && !pf) {
                            rv = clamp_u8(ct * 0.3f);
                            gv = clamp_u8(ct * 0.3f);
                            bv = clamp_u8(ct * 0.4f + 255.0f * 0.6f);
                        } else {
                            rv = gv = bv = (uint8_t)ct;
                        }
                        break;
                    }

                    case ROW_CONTOUR:
                        if (pred_img[mi] >= 1 && is_boundary(pred_img, img_h, img_w, r, x)) {
                            rv = clamp_u8(ct * 0.3f + 255.0f * 0.7f);
                            gv = clamp_u8(ct * 0.2f);
                            bv = clamp_u8(ct * 0.2f);
                        } else {
                            rv = gv = bv = (uint8_t)ct;
                        }
                        break;

                    default:
                        rv = gv = bv = (uint8_t)ct;
                        break;
                    }

                    int cx = cell_x + x;
                    int cy = cell_y + y_off + r;
                    size_t idx = ((size_t)cy * (size_t)canvas_w + (size_t)cx) * 3;
                    canvas[idx + 0] = rv;
                    canvas[idx + 1] = gv;
                    canvas[idx + 2] = bv;
                }
            }
        }
    }

    /* Title */
    {
        char title[128];
        snprintf(title, sizeof(title), "MIP: %s  (%dx%dx%d)", cube_id, D, H, W);
        draw_string_centered(canvas, canvas_w, canvas_h, 2, title, 0xFF, 0xFF, 0xFF);
    }

    struct timespec t_comp;
    clock_gettime(CLOCK_MONOTONIC, &t_comp);
    double comp_sec = (t_comp.tv_sec - t_mip.tv_sec) + (t_comp.tv_nsec - t_mip.tv_nsec) * 1e-9;
    fprintf(stderr, "  composited %d x %d canvas in %.2fs\n", canvas_w, canvas_h, comp_sec);

    /* Write */
    int rc = write_png(out_path, canvas, canvas_w, canvas_h);
    free(canvas);
    free(raw_mip_z); free(raw_mip_y); free(raw_mip_x);
    free(gt_mip_z);  free(gt_mip_y);  free(gt_mip_x);
    free(pred_mip_z); free(pred_mip_y); free(pred_mip_x);

    if (rc != 0) return 1;

    clock_gettime(CLOCK_MONOTONIC, &t1);
    double total_sec = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) * 1e-9;
    fprintf(stderr, "  wrote %s (%.2fs PNG)\n", out_path,
            (t1.tv_sec - t_comp.tv_sec) + (t1.tv_nsec - t_comp.tv_nsec) * 1e-9);
    fprintf(stderr, "  total: %.2fs\n", total_sec);

    return 0;
}
