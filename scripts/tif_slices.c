/*
 * tif_slices.c — Extract z-slices from a 3D binary TIFF as PNGs
 *
 * Usage: tif_slices <input.tif> <output_dir>
 *        tif_slices <input.tif> <input2.tif> <output_dir>   (side-by-side)
 *        tif_slices --overlay <raw.tif> <mask_a.tif> <mask_b.tif> <output_dir>
 *
 * When two TIFFs are given, produces side-by-side comparison PNGs:
 *   left = first TIFF (white), right = second TIFF (white/green for diff)
 *
 * --overlay mode: raw CT background (grayscale) with two mask overlays
 *   side-by-side. mask_a shown in cyan, mask_b shown in cyan.
 *   For GT labels, only value==1 is treated as foreground (value 2 = ignore).
 *
 * Open output_dir in eog and use arrow keys to scrub through slices.
 *
 * Build: cc -O2 -o tif_slices tif_slices.c -ltiff -lpng -lm
 */

#define _POSIX_C_SOURCE 200809L

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <sys/stat.h>
#include <tiffio.h>
#include <png.h>

/* ------------------------------------------------------------------ */
/* TIFF I/O                                                           */
/* ------------------------------------------------------------------ */

static uint8_t *load_tiff(const char *path, int *D, int *H, int *W)
{
    TIFFSetWarningHandler(NULL);
    TIFF *tif = TIFFOpen(path, "r");
    if (!tif) {
        fprintf(stderr, "ERROR: cannot open %s\n", path);
        return NULL;
    }

    uint32_t w = 0, h = 0;
    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);

    int d = 0;
    do { d++; } while (TIFFReadDirectory(tif));

    *D = d; *H = (int)h; *W = (int)w;

    size_t slice_sz = (size_t)h * (size_t)w;
    uint8_t *vol = (uint8_t *)malloc((size_t)d * slice_sz);
    if (!vol) {
        fprintf(stderr, "ERROR: cannot allocate %zu MB\n",
                (size_t)d * slice_sz / (1024 * 1024));
        TIFFClose(tif);
        return NULL;
    }

    TIFFSetDirectory(tif, 0);
    for (int z = 0; z < d; z++) {
        for (int y = 0; y < (int)h; y++)
            TIFFReadScanline(tif, vol + (size_t)z * slice_sz + (size_t)y * w, y, 0);
        TIFFReadDirectory(tif);
    }

    TIFFClose(tif);
    return vol;
}

static uint16_t *load_tiff_u16(const char *path, int *D, int *H, int *W)
{
    TIFFSetWarningHandler(NULL);
    TIFF *tif = TIFFOpen(path, "r");
    if (!tif) {
        fprintf(stderr, "ERROR: cannot open %s\n", path);
        return NULL;
    }

    uint32_t w = 0, h = 0;
    uint16_t bps = 0;
    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
    TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bps);

    int d = 0;
    do { d++; } while (TIFFReadDirectory(tif));

    *D = d; *H = (int)h; *W = (int)w;

    size_t slice_sz = (size_t)h * (size_t)w;
    uint16_t *vol = (uint16_t *)malloc((size_t)d * slice_sz * sizeof(uint16_t));
    if (!vol) {
        fprintf(stderr, "ERROR: cannot allocate %zu MB\n",
                (size_t)d * slice_sz * 2 / (1024 * 1024));
        TIFFClose(tif);
        return NULL;
    }

    TIFFSetDirectory(tif, 0);
    if (bps == 16) {
        for (int z = 0; z < d; z++) {
            for (int y = 0; y < (int)h; y++)
                TIFFReadScanline(tif, vol + (size_t)z * slice_sz + (size_t)y * w, y, 0);
            TIFFReadDirectory(tif);
        }
    } else {
        /* 8-bit: promote to 16-bit */
        uint8_t *row = (uint8_t *)malloc(w);
        for (int z = 0; z < d; z++) {
            for (int y = 0; y < (int)h; y++) {
                TIFFReadScanline(tif, row, y, 0);
                for (uint32_t x = 0; x < w; x++)
                    vol[(size_t)z * slice_sz + (size_t)y * w + x] = (uint16_t)row[x] * 257;
            }
            TIFFReadDirectory(tif);
        }
        free(row);
    }

    TIFFClose(tif);
    return vol;
}

/* ------------------------------------------------------------------ */
/* PNG writer                                                         */
/* ------------------------------------------------------------------ */

static int write_png_rgb(const char *path, const uint8_t *rgb, int w, int h)
{
    FILE *fp = fopen(path, "wb");
    if (!fp) return -1;

    png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING,
                                              NULL, NULL, NULL);
    png_infop info = png_create_info_struct(png);
    if (setjmp(png_jmpbuf(png))) {
        fclose(fp);
        png_destroy_write_struct(&png, &info);
        return -1;
    }

    png_init_io(png, fp);
    png_set_IHDR(png, info, (uint32_t)w, (uint32_t)h, 8,
                 PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
    png_set_compression_level(png, 1); /* fast */
    png_write_info(png, info);

    for (int y = 0; y < h; y++)
        png_write_row(png, rgb + (size_t)y * (size_t)w * 3);

    png_write_end(png, NULL);
    png_destroy_write_struct(&png, &info);
    fclose(fp);
    return 0;
}

/* ------------------------------------------------------------------ */
/* Overlay mode                                                       */
/* ------------------------------------------------------------------ */

static int run_overlay(const char *raw_path, const char *mask_a_path,
                       const char *mask_b_path, const char *out_dir)
{
    /* Load raw CT as uint16 */
    int Dr = 0, Hr = 0, Wr = 0;
    uint16_t *raw = load_tiff_u16(raw_path, &Dr, &Hr, &Wr);
    if (!raw) return 1;
    fprintf(stderr, "Loaded raw %s: %d x %d x %d\n", raw_path, Dr, Hr, Wr);

    /* Find min/max for contrast stretching */
    size_t vol_sz = (size_t)Dr * (size_t)Hr * (size_t)Wr;
    uint16_t rmin = 65535, rmax = 0;
    for (size_t i = 0; i < vol_sz; i++) {
        if (raw[i] < rmin) rmin = raw[i];
        if (raw[i] > rmax) rmax = raw[i];
    }
    float rscale = (rmax > rmin) ? 255.0f / (float)(rmax - rmin) : 1.0f;
    fprintf(stderr, "  raw range: [%u, %u]\n", rmin, rmax);

    /* Load masks as uint8 */
    int Da = 0, Ha = 0, Wa = 0;
    uint8_t *mask_a = load_tiff(mask_a_path, &Da, &Ha, &Wa);
    if (!mask_a) { free(raw); return 1; }
    fprintf(stderr, "Loaded mask_a %s: %d x %d x %d\n", mask_a_path, Da, Ha, Wa);

    int Db = 0, Hb = 0, Wb = 0;
    uint8_t *mask_b = load_tiff(mask_b_path, &Db, &Hb, &Wb);
    if (!mask_b) { free(raw); free(mask_a); return 1; }
    fprintf(stderr, "Loaded mask_b %s: %d x %d x %d\n", mask_b_path, Db, Hb, Wb);

    if (Dr != Da || Hr != Ha || Wr != Wa || Dr != Db || Hr != Hb || Wr != Wb) {
        fprintf(stderr, "ERROR: dimensions don't match\n");
        free(raw); free(mask_a); free(mask_b);
        return 1;
    }

    mkdir(out_dir, 0755);

    int D = Dr, H = Hr, W = Wr;
    size_t slice_sz = (size_t)H * (size_t)W;
    int gap = 4;
    int out_w = W * 2 + gap;
    int out_h = H;

    uint8_t *rgb = (uint8_t *)malloc((size_t)out_w * (size_t)out_h * 3);
    if (!rgb) {
        fprintf(stderr, "ERROR: cannot allocate RGB buffer\n");
        free(raw); free(mask_a); free(mask_b);
        return 1;
    }

    fprintf(stderr, "Writing %d slices to %s/\n", D, out_dir);

    for (int z = 0; z < D; z++) {
        const uint16_t *sr = raw + (size_t)z * slice_sz;
        const uint8_t *sa = mask_a + (size_t)z * slice_sz;
        const uint8_t *sb = mask_b + (size_t)z * slice_sz;

        for (int y = 0; y < H; y++) {
            for (int x = 0; x < W; x++) {
                /* Grayscale background from raw CT */
                uint8_t g = (uint8_t)((float)(sr[y * W + x] - rmin) * rscale);

                /* Left panel: raw + mask_a overlay */
                size_t px = ((size_t)y * (size_t)out_w + (size_t)x) * 3;
                int a_fg = (sa[y * W + x] == 1); /* only value 1 = foreground */
                if (a_fg) {
                    /* Cyan overlay blended with background */
                    rgb[px]     = (uint8_t)(g / 3);
                    rgb[px + 1] = (uint8_t)(g / 3 + 170);
                    rgb[px + 2] = (uint8_t)(g / 3 + 170);
                } else {
                    rgb[px] = rgb[px + 1] = rgb[px + 2] = g;
                }

                /* Right panel: raw + mask_b overlay */
                size_t px2 = ((size_t)y * (size_t)out_w + (size_t)(W + gap + x)) * 3;
                int b_fg = (sb[y * W + x] != 0);
                if (b_fg) {
                    rgb[px2]     = (uint8_t)(g / 3);
                    rgb[px2 + 1] = (uint8_t)(g / 3 + 170);
                    rgb[px2 + 2] = (uint8_t)(g / 3 + 170);
                } else {
                    rgb[px2] = rgb[px2 + 1] = rgb[px2 + 2] = g;
                }
            }

            /* Gap stripe */
            for (int gi = 0; gi < gap; gi++) {
                size_t px = ((size_t)y * (size_t)out_w + (size_t)(W + gi)) * 3;
                rgb[px] = rgb[px + 1] = rgb[px + 2] = 40;
            }
        }

        char out_path[1024];
        snprintf(out_path, sizeof(out_path), "%s/z_%03d.png", out_dir, z);
        if (write_png_rgb(out_path, rgb, out_w, out_h) != 0) {
            fprintf(stderr, "ERROR: failed to write %s\n", out_path);
        }

        if ((z + 1) % 50 == 0 || z == D - 1)
            fprintf(stderr, "  %d/%d slices\n", z + 1, D);
    }

    fprintf(stderr, "Done. Open with: eog %s/z_000.png\n", out_dir);

    free(rgb);
    free(raw);
    free(mask_a);
    free(mask_b);
    return 0;
}

/* ------------------------------------------------------------------ */
/* Main                                                               */
/* ------------------------------------------------------------------ */

int main(int argc, char **argv)
{
    /* --overlay mode */
    if (argc >= 2 && strcmp(argv[1], "--overlay") == 0) {
        if (argc != 6) {
            fprintf(stderr, "Usage: %s --overlay <raw.tif> <mask_a.tif> <mask_b.tif> <output_dir>\n", argv[0]);
            fprintf(stderr, "  Left panel:  raw CT + mask_a (cyan, value==1 only)\n");
            fprintf(stderr, "  Right panel: raw CT + mask_b (cyan, any nonzero)\n");
            return 1;
        }
        return run_overlay(argv[2], argv[3], argv[4], argv[5]);
    }

    if (argc < 3 || argc > 4) {
        fprintf(stderr, "Usage: %s <input.tif> [input2.tif] <output_dir>\n", argv[0]);
        fprintf(stderr, "       %s --overlay <raw.tif> <mask_a.tif> <mask_b.tif> <output_dir>\n", argv[0]);
        fprintf(stderr, "  Single TIFF:  white = foreground, black = background\n");
        fprintf(stderr, "  Two TIFFs:    side-by-side, green = added in second\n");
        fprintf(stderr, "  --overlay:    raw CT background + two mask overlays side-by-side\n");
        return 1;
    }

    int two_inputs = (argc == 4);
    const char *path_a = argv[1];
    const char *path_b = two_inputs ? argv[2] : NULL;
    const char *out_dir = argv[argc - 1];

    /* Load first TIFF */
    int Da = 0, Ha = 0, Wa = 0;
    uint8_t *vol_a = load_tiff(path_a, &Da, &Ha, &Wa);
    if (!vol_a) return 1;
    fprintf(stderr, "Loaded %s: %d x %d x %d\n", path_a, Da, Ha, Wa);

    /* Load optional second TIFF */
    int Db = 0, Hb = 0, Wb = 0;
    uint8_t *vol_b = NULL;
    if (two_inputs) {
        vol_b = load_tiff(path_b, &Db, &Hb, &Wb);
        if (!vol_b) { free(vol_a); return 1; }
        fprintf(stderr, "Loaded %s: %d x %d x %d\n", path_b, Db, Hb, Wb);
        if (Da != Db || Ha != Hb || Wa != Wb) {
            fprintf(stderr, "ERROR: dimensions don't match\n");
            free(vol_a); free(vol_b); return 1;
        }
    }

    /* Create output directory */
    mkdir(out_dir, 0755);

    int D = Da, H = Ha, W = Wa;
    size_t slice_sz = (size_t)H * (size_t)W;

    /* Output image dimensions */
    int gap = two_inputs ? 4 : 0;
    int out_w = two_inputs ? W * 2 + gap : W;
    int out_h = H;

    uint8_t *rgb = (uint8_t *)malloc((size_t)out_w * (size_t)out_h * 3);
    if (!rgb) {
        fprintf(stderr, "ERROR: cannot allocate RGB buffer\n");
        free(vol_a); free(vol_b); return 1;
    }

    fprintf(stderr, "Writing %d slices to %s/\n", D, out_dir);

    for (int z = 0; z < D; z++) {
        const uint8_t *sa = vol_a + (size_t)z * slice_sz;
        const uint8_t *sb = vol_b ? vol_b + (size_t)z * slice_sz : NULL;

        memset(rgb, 0, (size_t)out_w * (size_t)out_h * 3);

        for (int y = 0; y < H; y++) {
            for (int x = 0; x < W; x++) {
                /* Left panel (or only panel): first TIFF in white */
                size_t px = ((size_t)y * (size_t)out_w + (size_t)x) * 3;
                if (sa[y * W + x]) {
                    rgb[px] = rgb[px + 1] = rgb[px + 2] = 255;
                }

                /* Right panel: second TIFF with diff coloring */
                if (sb) {
                    size_t px2 = ((size_t)y * (size_t)out_w + (size_t)(W + gap + x)) * 3;
                    int a_val = sa[y * W + x];
                    int b_val = sb[y * W + x];
                    if (b_val && !a_val) {
                        /* Added in second: green */
                        rgb[px2]     = 0;
                        rgb[px2 + 1] = 220;
                        rgb[px2 + 2] = 80;
                    } else if (a_val && !b_val) {
                        /* Removed in second: red */
                        rgb[px2]     = 220;
                        rgb[px2 + 1] = 60;
                        rgb[px2 + 2] = 60;
                    } else if (b_val) {
                        /* Unchanged: white */
                        rgb[px2] = rgb[px2 + 1] = rgb[px2 + 2] = 255;
                    }
                }
            }

            /* Gap stripe between panels: dark gray */
            if (two_inputs) {
                for (int g = 0; g < gap; g++) {
                    size_t px = ((size_t)y * (size_t)out_w + (size_t)(W + g)) * 3;
                    rgb[px] = rgb[px + 1] = rgb[px + 2] = 40;
                }
            }
        }

        char out_path[1024];
        snprintf(out_path, sizeof(out_path), "%s/z_%03d.png", out_dir, z);
        if (write_png_rgb(out_path, rgb, out_w, out_h) != 0) {
            fprintf(stderr, "ERROR: failed to write %s\n", out_path);
        }

        if ((z + 1) % 50 == 0 || z == D - 1)
            fprintf(stderr, "  %d/%d slices\n", z + 1, D);
    }

    fprintf(stderr, "Done. Open with: eog %s/z_000.png\n", out_dir);

    free(rgb);
    free(vol_a);
    free(vol_b);
    return 0;
}
