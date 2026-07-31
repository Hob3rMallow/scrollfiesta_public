/* strip_preview -- human-viewable renditions of whole-grid strip TIFs.
 * Entry point only: parse args, call StripPreview_*, print timing.
 *
 *   strip_preview <in.tif> <out.png> [--mode gray|class|coverage]
 *                 [--ds N] [--ds-v M]
 *                 [--rows K] [--max-w W] [--gap G] [--overlay class.tif]
 *                 [--crop x0 x1 [y0 y1]]      full-res window (strip px)
 *                 [--cols file.txt --win N]   stacked seam windows
 *                 [--quiet]
 *   strip_preview --selftest
 *
 * The default (no --crop/--cols) is the downsampled whole-strip overview:
 * at 4x21x21 (--ds 16) that is a ~18,650 x ~2,068 PNG from the 610 MB TIF,
 * streamed row-by-row (peak memory = output image, never the input). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../common/arena.h"
#include "../common/ves_platform.h"
#include "../unroll/strip_preview.h"

static int read_cols(const char *path, long **out_cols, size_t *out_n)
{
    FILE *f = fopen(path, "r");
    if (f == NULL) {
        fprintf(stderr, "strip_preview: cannot open cols file %s\n", path);
        return -1;
    }
    size_t cap = 64, n = 0;
    long *cols = (long *)malloc(cap * sizeof(long));
    char line[256];
    while (cols != NULL && fgets(line, sizeof(line), f) != NULL) {
        if (line[0] == '#' || line[0] == '\n' || line[0] == '\r')
            continue;
        long v = strtol(line, NULL, 10);
        if (n == cap) {
            cap *= 2;
            long *nc = (long *)realloc(cols, cap * sizeof(long));
            if (nc == NULL) { free(cols); cols = NULL; break; }
            cols = nc;
        }
        cols[n++] = v;
    }
    fclose(f);
    if (cols == NULL)
        return -1;
    *out_cols = cols;
    *out_n = n;
    return 0;
}

static void usage(const char *argv0)
{
    fprintf(stderr,
        "usage: %s <in.tif> <out.png> [options]\n"
        "       %s --selftest\n"
        "  --mode gray|class|coverage\n"
        "                      palette (coverage = ribbon omission reasons)\n"
        "  --ds N              u box-downsample factor (default 16)\n"
        "  --ds-v M            v box-downsample factor (default 1)\n"
        "  --rows K            stacked rows (default: auto to fit --max-w)\n"
        "  --max-w W           output width cap for auto rows (default 20000)\n"
        "  --gap G             blank rows between stacks (default 8)\n"
        "  --overlay f.tif     gray mode: tint multi-cover (class 4) red\n"
        "  --crop x0 x1 [y0 y1]  full-res window instead of overview\n"
        "  --cols f.txt --win N  stacked full-res +/-N windows at the listed\n"
        "                        columns (seam zooms; one strip-px col/line)\n"
        "  --quiet             suppress geometry logging\n",
        argv0, argv0);
}

int main(int argc, char **argv)
{
    if (argc >= 2 && strcmp(argv[1], "--selftest") == 0)
        return StripPreview_selftest() == 0 ? 0 : 1;
    if (argc < 3) {
        usage(argv[0]);
        return 1;
    }

    const char *in_tif = argv[1];
    const char *out_png = argv[2];
    StripPreviewOpts o;
    StripPreviewOpts_default(&o);

    int do_crop = 0;
    long crop_x0 = 0, crop_x1 = 0;
    int crop_y0 = 0, crop_y1 = 0;
    const char *cols_file = NULL;
    int win = 256;

    for (int a = 3; a < argc; a++) {
        if (strcmp(argv[a], "--mode") == 0 && a + 1 < argc) {
            o.mode = argv[++a];
        } else if (strcmp(argv[a], "--ds") == 0 && a + 1 < argc) {
            o.ds_u = atoi(argv[++a]);
        } else if (strcmp(argv[a], "--ds-v") == 0 && a + 1 < argc) {
            o.ds_v = atoi(argv[++a]);
        } else if (strcmp(argv[a], "--rows") == 0 && a + 1 < argc) {
            o.rows = atoi(argv[++a]);
        } else if (strcmp(argv[a], "--max-w") == 0 && a + 1 < argc) {
            o.max_w = atoi(argv[++a]);
        } else if (strcmp(argv[a], "--gap") == 0 && a + 1 < argc) {
            o.gap = atoi(argv[++a]);
        } else if (strcmp(argv[a], "--overlay") == 0 && a + 1 < argc) {
            o.overlay = argv[++a];
        } else if (strcmp(argv[a], "--crop") == 0 && a + 2 < argc) {
            do_crop = 1;
            crop_x0 = strtol(argv[++a], NULL, 10);
            crop_x1 = strtol(argv[++a], NULL, 10);
            if (a + 2 < argc && argv[a + 1][0] != '-') {
                crop_y0 = atoi(argv[++a]);
                crop_y1 = atoi(argv[++a]);
            }
        } else if (strcmp(argv[a], "--cols") == 0 && a + 1 < argc) {
            cols_file = argv[++a];
        } else if (strcmp(argv[a], "--win") == 0 && a + 1 < argc) {
            win = atoi(argv[++a]);
        } else if (strcmp(argv[a], "--quiet") == 0) {
            o.verbose = 0;
        } else {
            fprintf(stderr, "strip_preview: unknown option %s\n", argv[a]);
            usage(argv[0]);
            return 1;
        }
    }
    if (o.ds_u < 1) o.ds_u = 1;
    if (o.ds_v < 1) o.ds_v = 1;
    if (o.max_w < 256) o.max_w = 256;
    if (o.gap < 0) o.gap = 0;
    if (strcmp(o.mode, "gray") != 0 && strcmp(o.mode, "class") != 0 &&
        strcmp(o.mode, "coverage") != 0) {
        fprintf(stderr,
                "strip_preview: mode must be gray, class, or coverage\n");
        return 1;
    }

    Arena_T arena = Arena_new();
    double t0 = ves_clock_sec();
    int rc = 0;

    if (cols_file != NULL) {
        long *cols = NULL;
        size_t ncols = 0;
        rc = read_cols(cols_file, &cols, &ncols);
        if (rc == 0)
            rc = StripPreview_windows(arena, in_tif, out_png,
                                      cols, ncols, win, &o);
        free(cols);
    } else if (do_crop) {
        rc = StripPreview_crop(arena, in_tif, out_png,
                               crop_x0, crop_x1, crop_y0, crop_y1, &o);
    } else {
        rc = StripPreview_overview(arena, in_tif, out_png, &o);
    }

    if (rc == 0)
        fprintf(stderr, "strip_preview: OK (%.1fs)\n", ves_clock_sec() - t0);
    else
        fprintf(stderr, "strip_preview: FAILED\n");
    Arena_dispose(&arena);
    return rc == 0 ? 0 : 1;
}
