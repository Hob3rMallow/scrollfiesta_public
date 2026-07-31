/*
 * tif2png.c — Standalone TIFF slice visualization tool
 *
 * Loads 3 multi-page TIFFs (raw CT, GT labels, predictions), composites
 * 10 equidistant Z-slices across 3 rows (CT / GT overlay / Pred overlay),
 * renders text labels with an embedded 8x8 bitmap font, and writes PNG.
 *
 * Usage: tif2png <raw.tif> <gt.tif> <pred.tif> <output.png> [cube_id] [--diff] [--contour]
 *
 * Options:
 *   --diff      Add row: GT vs Pred diff (green=TP, red=FP, blue=FN, gray=TN)
 *   --contour   Add row: Pred boundary contour on CT (red outline only)
 *
 * Build: gcc -O2 -o tif2png scripts/tif2png.c -ltiff -lpng16 -lm
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
    /* 32 ' ' */ {0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00},
    /* 33 '!' */ {0x18,0x3C,0x3C,0x18,0x18,0x00,0x18,0x00},
    /* 34 '"' */ {0x6C,0x6C,0x24,0x00,0x00,0x00,0x00,0x00},
    /* 35 '#' */ {0x6C,0xFE,0x6C,0x6C,0xFE,0x6C,0x00,0x00},
    /* 36 '$' */ {0x18,0x7E,0xC0,0x7C,0x06,0xFC,0x18,0x00},
    /* 37 '%' */ {0x00,0xC6,0xCC,0x18,0x30,0x66,0xC6,0x00},
    /* 38 '&' */ {0x38,0x6C,0x38,0x76,0xDC,0xCC,0x76,0x00},
    /* 39 ''' */ {0x18,0x18,0x30,0x00,0x00,0x00,0x00,0x00},
    /* 40 '(' */ {0x0C,0x18,0x30,0x30,0x30,0x18,0x0C,0x00},
    /* 41 ')' */ {0x30,0x18,0x0C,0x0C,0x0C,0x18,0x30,0x00},
    /* 42 '*' */ {0x00,0x66,0x3C,0xFF,0x3C,0x66,0x00,0x00},
    /* 43 '+' */ {0x00,0x18,0x18,0x7E,0x18,0x18,0x00,0x00},
    /* 44 ',' */ {0x00,0x00,0x00,0x00,0x00,0x18,0x18,0x30},
    /* 45 '-' */ {0x00,0x00,0x00,0x7E,0x00,0x00,0x00,0x00},
    /* 46 '.' */ {0x00,0x00,0x00,0x00,0x00,0x18,0x18,0x00},
    /* 47 '/' */ {0x06,0x0C,0x18,0x30,0x60,0xC0,0x80,0x00},
    /* 48 '0' */ {0x7C,0xC6,0xCE,0xDE,0xF6,0xE6,0x7C,0x00},
    /* 49 '1' */ {0x18,0x38,0x78,0x18,0x18,0x18,0x7E,0x00},
    /* 50 '2' */ {0x7C,0xC6,0x06,0x1C,0x30,0x66,0xFE,0x00},
    /* 51 '3' */ {0x7C,0xC6,0x06,0x3C,0x06,0xC6,0x7C,0x00},
    /* 52 '4' */ {0x1C,0x3C,0x6C,0xCC,0xFE,0x0C,0x1E,0x00},
    /* 53 '5' */ {0xFE,0xC0,0xFC,0x06,0x06,0xC6,0x7C,0x00},
    /* 54 '6' */ {0x38,0x60,0xC0,0xFC,0xC6,0xC6,0x7C,0x00},
    /* 55 '7' */ {0xFE,0xC6,0x0C,0x18,0x30,0x30,0x30,0x00},
    /* 56 '8' */ {0x7C,0xC6,0xC6,0x7C,0xC6,0xC6,0x7C,0x00},
    /* 57 '9' */ {0x7C,0xC6,0xC6,0x7E,0x06,0x0C,0x78,0x00},
    /* 58 ':' */ {0x00,0x18,0x18,0x00,0x00,0x18,0x18,0x00},
    /* 59 ';' */ {0x00,0x18,0x18,0x00,0x00,0x18,0x18,0x30},
    /* 60 '<' */ {0x06,0x0C,0x18,0x30,0x18,0x0C,0x06,0x00},
    /* 61 '=' */ {0x00,0x00,0x7E,0x00,0x7E,0x00,0x00,0x00},
    /* 62 '>' */ {0x60,0x30,0x18,0x0C,0x18,0x30,0x60,0x00},
    /* 63 '?' */ {0x7C,0xC6,0x0C,0x18,0x18,0x00,0x18,0x00},
    /* 64 '@' */ {0x7C,0xC6,0xDE,0xDE,0xDE,0xC0,0x78,0x00},
    /* 65 'A' */ {0x38,0x6C,0xC6,0xC6,0xFE,0xC6,0xC6,0x00},
    /* 66 'B' */ {0xFC,0x66,0x66,0x7C,0x66,0x66,0xFC,0x00},
    /* 67 'C' */ {0x3C,0x66,0xC0,0xC0,0xC0,0x66,0x3C,0x00},
    /* 68 'D' */ {0xF8,0x6C,0x66,0x66,0x66,0x6C,0xF8,0x00},
    /* 69 'E' */ {0xFE,0x62,0x68,0x78,0x68,0x62,0xFE,0x00},
    /* 70 'F' */ {0xFE,0x62,0x68,0x78,0x68,0x60,0xF0,0x00},
    /* 71 'G' */ {0x3C,0x66,0xC0,0xC0,0xCE,0x66,0x3E,0x00},
    /* 72 'H' */ {0xC6,0xC6,0xC6,0xFE,0xC6,0xC6,0xC6,0x00},
    /* 73 'I' */ {0x3C,0x18,0x18,0x18,0x18,0x18,0x3C,0x00},
    /* 74 'J' */ {0x1E,0x0C,0x0C,0x0C,0xCC,0xCC,0x78,0x00},
    /* 75 'K' */ {0xE6,0x66,0x6C,0x78,0x6C,0x66,0xE6,0x00},
    /* 76 'L' */ {0xF0,0x60,0x60,0x60,0x62,0x66,0xFE,0x00},
    /* 77 'M' */ {0xC6,0xEE,0xFE,0xFE,0xD6,0xC6,0xC6,0x00},
    /* 78 'N' */ {0xC6,0xE6,0xF6,0xDE,0xCE,0xC6,0xC6,0x00},
    /* 79 'O' */ {0x7C,0xC6,0xC6,0xC6,0xC6,0xC6,0x7C,0x00},
    /* 80 'P' */ {0xFC,0x66,0x66,0x7C,0x60,0x60,0xF0,0x00},
    /* 81 'Q' */ {0x7C,0xC6,0xC6,0xC6,0xD6,0xDE,0x7C,0x06},
    /* 82 'R' */ {0xFC,0x66,0x66,0x7C,0x6C,0x66,0xE6,0x00},
    /* 83 'S' */ {0x7C,0xC6,0xE0,0x7C,0x0E,0xC6,0x7C,0x00},
    /* 84 'T' */ {0x7E,0x5A,0x18,0x18,0x18,0x18,0x3C,0x00},
    /* 85 'U' */ {0xC6,0xC6,0xC6,0xC6,0xC6,0xC6,0x7C,0x00},
    /* 86 'V' */ {0xC6,0xC6,0xC6,0xC6,0x6C,0x38,0x10,0x00},
    /* 87 'W' */ {0xC6,0xC6,0xD6,0xFE,0xFE,0xEE,0xC6,0x00},
    /* 88 'X' */ {0xC6,0x6C,0x38,0x38,0x6C,0xC6,0xC6,0x00},
    /* 89 'Y' */ {0x66,0x66,0x66,0x3C,0x18,0x18,0x3C,0x00},
    /* 90 'Z' */ {0xFE,0xC6,0x8C,0x18,0x32,0x66,0xFE,0x00},
    /* 91 '[' */ {0x3C,0x30,0x30,0x30,0x30,0x30,0x3C,0x00},
    /* 92 '\' */ {0xC0,0x60,0x30,0x18,0x0C,0x06,0x02,0x00},
    /* 93 ']' */ {0x3C,0x0C,0x0C,0x0C,0x0C,0x0C,0x3C,0x00},
    /* 94 '^' */ {0x10,0x38,0x6C,0xC6,0x00,0x00,0x00,0x00},
    /* 95 '_' */ {0x00,0x00,0x00,0x00,0x00,0x00,0x00,0xFF},
    /* 96 '`' */ {0x30,0x18,0x0C,0x00,0x00,0x00,0x00,0x00},
    /* 97 'a' */ {0x00,0x00,0x78,0x0C,0x7C,0xCC,0x76,0x00},
    /* 98 'b' */ {0xE0,0x60,0x60,0x7C,0x66,0x66,0xDC,0x00},
    /* 99 'c' */ {0x00,0x00,0x7C,0xC6,0xC0,0xC6,0x7C,0x00},
    /*100 'd' */ {0x1C,0x0C,0x0C,0x7C,0xCC,0xCC,0x76,0x00},
    /*101 'e' */ {0x00,0x00,0x7C,0xC6,0xFE,0xC0,0x7C,0x00},
    /*102 'f' */ {0x38,0x6C,0x60,0xF0,0x60,0x60,0xF0,0x00},
    /*103 'g' */ {0x00,0x00,0x76,0xCC,0xCC,0x7C,0x0C,0xF8},
    /*104 'h' */ {0xE0,0x60,0x6C,0x76,0x66,0x66,0xE6,0x00},
    /*105 'i' */ {0x18,0x00,0x38,0x18,0x18,0x18,0x3C,0x00},
    /*106 'j' */ {0x06,0x00,0x06,0x06,0x06,0x66,0x66,0x3C},
    /*107 'k' */ {0xE0,0x60,0x66,0x6C,0x78,0x6C,0xE6,0x00},
    /*108 'l' */ {0x38,0x18,0x18,0x18,0x18,0x18,0x3C,0x00},
    /*109 'm' */ {0x00,0x00,0xEC,0xFE,0xD6,0xD6,0xD6,0x00},
    /*110 'n' */ {0x00,0x00,0xDC,0x66,0x66,0x66,0x66,0x00},
    /*111 'o' */ {0x00,0x00,0x7C,0xC6,0xC6,0xC6,0x7C,0x00},
    /*112 'p' */ {0x00,0x00,0xDC,0x66,0x66,0x7C,0x60,0xF0},
    /*113 'q' */ {0x00,0x00,0x76,0xCC,0xCC,0x7C,0x0C,0x1E},
    /*114 'r' */ {0x00,0x00,0xDC,0x76,0x66,0x60,0xF0,0x00},
    /*115 's' */ {0x00,0x00,0x7C,0xC0,0x7C,0x06,0xFC,0x00},
    /*116 't' */ {0x10,0x30,0x7C,0x30,0x30,0x34,0x18,0x00},
    /*117 'u' */ {0x00,0x00,0xCC,0xCC,0xCC,0xCC,0x76,0x00},
    /*118 'v' */ {0x00,0x00,0xC6,0xC6,0x6C,0x38,0x10,0x00},
    /*119 'w' */ {0x00,0x00,0xC6,0xD6,0xD6,0xFE,0x6C,0x00},
    /*120 'x' */ {0x00,0x00,0xC6,0x6C,0x38,0x6C,0xC6,0x00},
    /*121 'y' */ {0x00,0x00,0xC6,0xC6,0x6C,0x38,0x30,0x60},
    /*122 'z' */ {0x00,0x00,0xFE,0x8C,0x18,0x32,0xFE,0x00},
    /*123 '{' */ {0x0E,0x18,0x18,0x70,0x18,0x18,0x0E,0x00},
    /*124 '|' */ {0x18,0x18,0x18,0x00,0x18,0x18,0x18,0x00},
    /*125 '}' */ {0x70,0x18,0x18,0x0E,0x18,0x18,0x70,0x00},
    /*126 '~' */ {0x76,0xDC,0x00,0x00,0x00,0x00,0x00,0x00},
    /*127 DEL */ {0x00,0x10,0x38,0x6C,0xC6,0xFE,0x00,0x00},
};

#define FONT_SCALE 2   /* 2x scale: 16x16 per glyph */
#define GLYPH_W (8 * FONT_SCALE)
#define GLYPH_H (8 * FONT_SCALE)

/* Draw a single character at (px, py) on RGB canvas */
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
                /* Fill FONT_SCALE x FONT_SCALE block */
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

/* Draw a string at (px, py) */
static void draw_string(uint8_t *canvas, int cw, int ch,
                        int px, int py, const char *str,
                        uint8_t r, uint8_t g, uint8_t b)
{
    for (int i = 0; str[i]; i++) {
        draw_char(canvas, cw, ch, px + i * GLYPH_W, py, str[i], r, g, b);
    }
}

/* Draw a string centered horizontally at y */
static void draw_string_centered(uint8_t *canvas, int cw, int ch,
                                 int py, const char *str,
                                 uint8_t r, uint8_t g, uint8_t b)
{
    int len = (int)strlen(str);
    int tw = len * GLYPH_W;
    int px = (cw - tw) / 2;
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

    /* Count pages */
    int depth = 0;
    do { depth++; } while (TIFFReadDirectory(tif));

    size_t vol_size = (size_t)depth * (size_t)h * (size_t)w;
    uint8_t *vol = (uint8_t *)calloc(vol_size, 1);
    if (!vol) {
        fprintf(stderr, "ERROR: cannot allocate %zu bytes for volume\n", vol_size);
        TIFFClose(tif);
        return NULL;
    }

    /* Re-read from beginning */
    TIFFSetDirectory(tif, 0);
    uint8_t *row_buf = (uint8_t *)malloc((size_t)TIFFScanlineSize(tif));
    if (!row_buf) {
        fprintf(stderr, "ERROR: cannot allocate scanline buffer\n");
        free(vol);
        TIFFClose(tif);
        return NULL;
    }

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
    if (!fp) {
        fprintf(stderr, "ERROR: cannot open output: %s\n", path);
        return -1;
    }

    png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING,
                                              NULL, NULL, NULL);
    if (!png) {
        fclose(fp);
        return -1;
    }

    png_infop info = png_create_info_struct(png);
    if (!info) {
        png_destroy_write_struct(&png, NULL);
        fclose(fp);
        return -1;
    }

    if (setjmp(png_jmpbuf(png))) {
        png_destroy_write_struct(&png, &info);
        fclose(fp);
        return -1;
    }

    png_init_io(png, fp);
    png_set_compression_level(png, 6);

    png_set_IHDR(png, info, (png_uint_32)w, (png_uint_32)h,
                 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
    png_write_info(png, info);

    for (int y = 0; y < h; y++) {
        png_write_row(png, (const png_bytep)(canvas + (size_t)y * (size_t)w * 3));
    }

    png_write_end(png, NULL);
    png_destroy_write_struct(&png, &info);
    fclose(fp);
    return 0;
}

/* ================================================================
 * Clamp helper
 * ================================================================ */

static inline uint8_t clamp_u8(float v)
{
    if (v < 0.0f) return 0;
    if (v > 255.0f) return 255;
    return (uint8_t)(v + 0.5f);
}

/* Check if pixel (y,x) in a 2D slice is on the boundary of the mask.
 * Boundary = mask is 1 here and at least one 4-neighbor is 0 or out of bounds. */
static inline int is_boundary(const uint8_t *slice, int H, int W, int y, int x)
{
    if (y == 0 || y == H - 1 || x == 0 || x == W - 1) return 1;
    if (!slice[(y - 1) * W + x]) return 1;
    if (!slice[(y + 1) * W + x]) return 1;
    if (!slice[y * W + (x - 1)]) return 1;
    if (!slice[y * W + (x + 1)]) return 1;
    return 0;
}

/* Row type enum */
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

    /* Parse optional positional + flags */
    for (int i = 5; i < argc; i++) {
        if (strcmp(argv[i], "--diff") == 0) {
            opt_diff = 1;
        } else if (strcmp(argv[i], "--contour") == 0) {
            opt_contour = 1;
        } else if (argv[i][0] != '-') {
            cube_id = argv[i];
        } else {
            fprintf(stderr, "Unknown option: %s\n", argv[i]);
            return 1;
        }
    }

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    /* Load 3 volumes */
    int rD = 0, rH = 0, rW = 0;
    uint8_t *raw = load_tiff(raw_path, &rD, &rH, &rW);
    if (!raw) return 1;
    fprintf(stderr, "  raw:  %d x %d x %d\n", rD, rH, rW);

    int gD = 0, gH = 0, gW = 0;
    uint8_t *gt = load_tiff(gt_path, &gD, &gH, &gW);
    if (!gt) { free(raw); return 1; }
    fprintf(stderr, "  gt:   %d x %d x %d\n", gD, gH, gW);

    int pD = 0, pH = 0, pW = 0;
    uint8_t *pred = load_tiff(pred_path, &pD, &pH, &pW);
    if (!pred) { free(raw); free(gt); return 1; }
    fprintf(stderr, "  pred: %d x %d x %d\n", pD, pH, pW);

    /* Validate dimensions match */
    if (rD != gD || rD != pD || rH != gH || rH != pH || rW != gW || rW != pW) {
        fprintf(stderr, "ERROR: dimension mismatch: raw=%dx%dx%d gt=%dx%dx%d pred=%dx%dx%d\n",
                rD, rH, rW, gD, gH, gW, pD, pH, pW);
        free(raw); free(gt); free(pred);
        return 1;
    }

    int D = rD, H = rH, W = rW;

    struct timespec t_load;
    clock_gettime(CLOCK_MONOTONIC, &t_load);
    double load_sec = (t_load.tv_sec - t0.tv_sec) + (t_load.tv_nsec - t0.tv_nsec) * 1e-9;
    fprintf(stderr, "  loaded 3 volumes in %.2fs\n", load_sec);

    /* Compute 10 equidistant Z-slice positions, keep middle 8 (drop endpoints) */
    int n_slices = 8;
    int slices[8];
    for (int i = 0; i < n_slices; i++) {
        slices[i] = (D - 1) * (i + 1) / 9;
    }

    /* Build row list: always CT, GT, Pred; optionally Diff and Contour */
    int row_types[ROW_MAX];
    int n_rows = 0;
    row_types[n_rows++] = ROW_CT;
    row_types[n_rows++] = ROW_GT;
    row_types[n_rows++] = ROW_PRED;
    if (opt_diff)    row_types[n_rows++] = ROW_DIFF;
    if (opt_contour) row_types[n_rows++] = ROW_CONTOUR;

    /* Canvas layout */
    int pad = 4;           /* padding between cells */
    int title_h = GLYPH_H + 4;     /* title bar height */
    int label_h = GLYPH_H + 2;     /* per-cell label height */
    int n_cols = n_slices;

    int canvas_w = pad + n_cols * (W + pad);
    int canvas_h = title_h + n_rows * (label_h + H + pad) + pad;

    size_t canvas_size = (size_t)canvas_w * (size_t)canvas_h * 3;
    uint8_t *canvas = (uint8_t *)malloc(canvas_size);
    if (!canvas) {
        fprintf(stderr, "ERROR: cannot allocate canvas (%d x %d = %zu bytes)\n",
                canvas_w, canvas_h, canvas_size);
        free(raw); free(gt); free(pred);
        return 1;
    }

    /* Fill with dark gray background (#303030) */
    memset(canvas, 0x30, canvas_size);

    /* Composite cells */
    static const char *row_name_table[] = {
        [ROW_CT]      = "CT",
        [ROW_GT]      = "GT",
        [ROW_PRED]    = "Pred",
        [ROW_DIFF]    = "Diff",
        [ROW_CONTOUR] = "Contour",
    };

    for (int col = 0; col < n_cols; col++) {
        int z = slices[col];
        int cell_x = pad + col * (W + pad);
        const uint8_t *pred_slice = pred + (size_t)z * (size_t)H * (size_t)W;

        for (int row = 0; row < n_rows; row++) {
            int rtype = row_types[row];
            int cell_y = title_h + row * (label_h + H + pad) + label_h;

            /* Draw label */
            char label[64];
            snprintf(label, sizeof(label), "%s z=%d", row_name_table[rtype], z);
            int label_y = title_h + row * (label_h + H + pad);
            draw_string(canvas, canvas_w, canvas_h,
                       cell_x + 2, label_y + 1, label,
                       0xCC, 0xCC, 0xCC);

            /* Composite pixels */
            for (int y = 0; y < H; y++) {
                for (int x = 0; x < W; x++) {
                    size_t vol_idx = (size_t)z * (size_t)H * (size_t)W
                                   + (size_t)y * (size_t)W + (size_t)x;
                    float ct = (float)raw[vol_idx];
                    uint8_t r, g, b;

                    switch (rtype) {
                    case ROW_CT:
                        /* Pure grayscale CT */
                        r = g = b = (uint8_t)ct;
                        break;

                    case ROW_GT: {
                        /* GT overlay: green=fg, yellow=unlabeled */
                        uint8_t gv = gt[vol_idx];
                        if (gv == 1) {
                            r = clamp_u8(ct * 0.3f);
                            g = clamp_u8(ct * 0.4f + 255.0f * 0.6f);
                            b = clamp_u8(ct * 0.3f);
                        } else if (gv == 2) {
                            r = clamp_u8(ct * 0.6f + 180.0f * 0.4f);
                            g = clamp_u8(ct * 0.6f + 160.0f * 0.4f);
                            b = clamp_u8(ct * 0.4f);
                        } else {
                            r = g = b = (uint8_t)ct;
                        }
                        break;
                    }

                    case ROW_PRED: {
                        /* Pred overlay: red fill */
                        uint8_t pv = pred[vol_idx];
                        if (pv == 1) {
                            r = clamp_u8(ct * 0.4f + 255.0f * 0.6f);
                            g = clamp_u8(ct * 0.3f);
                            b = clamp_u8(ct * 0.3f);
                        } else {
                            r = g = b = (uint8_t)ct;
                        }
                        break;
                    }

                    case ROW_DIFF: {
                        /* Diff: compare GT(==1) vs Pred(==1)
                         * TP (both 1): green
                         * FP (pred=1, gt!=1): red
                         * FN (gt=1, pred!=1): blue
                         * TN: plain grayscale */
                        uint8_t gv = gt[vol_idx];
                        uint8_t pv = pred[vol_idx];
                        int gt_fg = (gv == 1);
                        int pred_fg = (pv == 1);
                        if (gt_fg && pred_fg) {
                            /* TP: green */
                            r = clamp_u8(ct * 0.3f);
                            g = clamp_u8(ct * 0.4f + 255.0f * 0.6f);
                            b = clamp_u8(ct * 0.3f);
                        } else if (pred_fg && !gt_fg) {
                            /* FP: red */
                            r = clamp_u8(ct * 0.4f + 255.0f * 0.6f);
                            g = clamp_u8(ct * 0.3f);
                            b = clamp_u8(ct * 0.3f);
                        } else if (gt_fg && !pred_fg) {
                            /* FN: blue */
                            r = clamp_u8(ct * 0.3f);
                            g = clamp_u8(ct * 0.3f);
                            b = clamp_u8(ct * 0.4f + 255.0f * 0.6f);
                        } else {
                            r = g = b = (uint8_t)ct;
                        }
                        break;
                    }

                    case ROW_CONTOUR: {
                        /* Pred boundary contour: red outline where pred=1 meets pred=0 */
                        uint8_t pv = pred[vol_idx];
                        if (pv == 1 && is_boundary(pred_slice, H, W, y, x)) {
                            r = clamp_u8(ct * 0.3f + 255.0f * 0.7f);
                            g = clamp_u8(ct * 0.2f);
                            b = clamp_u8(ct * 0.2f);
                        } else {
                            r = g = b = (uint8_t)ct;
                        }
                        break;
                    }

                    default:
                        r = g = b = (uint8_t)ct;
                        break;
                    }

                    int cx_pos = cell_x + x;
                    int cy_pos = cell_y + y;
                    size_t idx = ((size_t)cy_pos * (size_t)canvas_w + (size_t)cx_pos) * 3;
                    canvas[idx + 0] = r;
                    canvas[idx + 1] = g;
                    canvas[idx + 2] = b;
                }
            }
        }
    }

    /* Draw title */
    {
        char title[128];
        snprintf(title, sizeof(title), "Cube %s  (%dx%dx%d)", cube_id, D, H, W);
        draw_string_centered(canvas, canvas_w, canvas_h, 2, title,
                           0xFF, 0xFF, 0xFF);
    }

    struct timespec t_comp;
    clock_gettime(CLOCK_MONOTONIC, &t_comp);
    double comp_sec = (t_comp.tv_sec - t_load.tv_sec) + (t_comp.tv_nsec - t_load.tv_nsec) * 1e-9;
    fprintf(stderr, "  composited %d x %d canvas in %.2fs\n", canvas_w, canvas_h, comp_sec);

    /* Write PNG */
    int rc = write_png(out_path, canvas, canvas_w, canvas_h);

    free(canvas);
    free(raw);
    free(gt);
    free(pred);

    if (rc != 0) return 1;

    clock_gettime(CLOCK_MONOTONIC, &t1);
    double total_sec = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) * 1e-9;
    double png_sec = (t1.tv_sec - t_comp.tv_sec) + (t1.tv_nsec - t_comp.tv_nsec) * 1e-9;
    fprintf(stderr, "  wrote %s (%.2fs)\n", out_path, png_sec);
    fprintf(stderr, "  total: %.2fs\n", total_sec);

    return 0;
}
