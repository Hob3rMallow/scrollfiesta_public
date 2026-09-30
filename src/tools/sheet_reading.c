/* sheet_reading -- a reading view of a sheet image: mirror-imaged pieces flipped top to bottom.
 *
 *   sheet_reading <in.tif> <handedness.json> <out_stem> [--strip <ledger.json>]
 *   sheet_reading --selftest
 *
 * <in.tif> is an 8-bit gray image (stripped or tiled) on the big sheet's grid (the review bake / sheet_layers
 * grid: a CT sheet, an ink map), or with --strip a strip page of that sheet (sheet_strip's <stem>.tif, CT or
 * --texture) together with its ledger.  <handedness.json> is sheet_layers' <id>_handedness.json for the same
 * sheet: per piece (a run of painted columns) the pixels showing the recto in the handedness of the team's
 * reading canvases and those showing it mirror-imaged.  The review layouts place some pieces with the scroll
 * axis running up the rows and others down, which is exactly the mirror: every piece with more mirror-imaged
 * than team-handed pixels is flipped top to bottom -- on the sheet grid within its own painted rows
 * [row0, row1), on a strip page within each of its fragments' page rectangles (page pixel (x, y) of fragment f
 * on line l shows source (f.sx0 + x - margin - f.dx, f.sy0 + y - l.dy)).  Pixels only move.  Which way is up
 * (0 or 180 degrees) is not decided here.
 * Writes <out_stem>.tif, <out_stem>.png (downsampled to fit the PNG writer) and <out_stem>.json. */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/arena.h"
#include "../common/json_read.h"
#include "../common/strip_layout.h"
#include "../common/tiff_io.h"
#include "../common/ves_platform.h"
#include "../common/ves_png.h"

#include "tiffio.h"

typedef struct { long col0, col1, row0, row1, team, mirror; } Piece;

/* Single-page 8-bit gray TIFF, stripped (TiffIO_load) or tiled (libtiff tiles), as sheet_strip reads textures. */
static int load_gray_any(Arena_T arena, const char *path, uint8_t **img, int *H, int *W)
{
    TIFF *tif = TIFFOpen(path, "r");
    uint32_t w = 0, h = 0, tw = 0, th = 0, x = 0, y = 0, r = 0;
    uint16_t bps = 0, spp = 0;
    uint8_t *out = NULL, *tile = NULL;
    if (tif == NULL) { fprintf(stderr, "sheet_reading: cannot open %s\n", path); return -1; }
    if (!TIFFIsTiled(tif)) {
        int D = 0;
        TIFFClose(tif);
        if (TiffIO_load(arena, path, img, &D, H, W) != 0 || D != 1) {
            fprintf(stderr, "sheet_reading: cannot load a single-page 8-bit TIFF %s\n", path);
            return -1;
        }
        return 0;
    }
    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
    TIFFGetFieldDefaulted(tif, TIFFTAG_BITSPERSAMPLE, &bps);
    TIFFGetFieldDefaulted(tif, TIFFTAG_SAMPLESPERPIXEL, &spp);
    TIFFGetField(tif, TIFFTAG_TILEWIDTH, &tw);
    TIFFGetField(tif, TIFFTAG_TILELENGTH, &th);
    if (w == 0 || h == 0 || w > INT32_MAX || h > INT32_MAX || bps != 8 || spp != 1 || tw == 0 || th == 0 ||
        TIFFNumberOfDirectories(tif) != 1) {
        fprintf(stderr, "sheet_reading: %s is not a single-page 8-bit gray TIFF\n", path);
        TIFFClose(tif);
        return -1;
    }
    out = (uint8_t *)ARENA_CALLOC(arena, (size_t)w * (size_t)h, 1);
    tile = (uint8_t *)ARENA_ALLOC(arena, (size_t)TIFFTileSize(tif));
    for (y = 0; y < h; y += th)
        for (x = 0; x < w; x += tw) {
            uint32_t cw = w - x < tw ? w - x : tw, ch = h - y < th ? h - y : th;
            if (TIFFReadTile(tif, tile, x, y, 0, 0) < 0) {
                fprintf(stderr, "sheet_reading: tile read failed in %s at (%u, %u)\n", path, x, y);
                TIFFClose(tif);
                return -1;
            }
            for (r = 0; r < ch; r++) memcpy(out + (size_t)(y + r) * w + x, tile + (size_t)r * tw, cw);
        }
    TIFFClose(tif);
    *img = out; *H = (int)h; *W = (int)w;
    return 0;
}

/* Reverse the rows [y0, y1) of columns [x0, x1) of a W-wide image in place. */
static void flip_rows(uint8_t *img, size_t W, long x0, long x1, long y0, long y1, uint8_t *tmp)
{
    long a = y0, b = y1 - 1;
    size_t n = (size_t)(x1 - x0);
    for (; a < b; a++, b--) {
        uint8_t *ra = img + (size_t)a * W + (size_t)x0, *rb = img + (size_t)b * W + (size_t)x0;
        memcpy(tmp, ra, n); memcpy(ra, rb, n); memcpy(rb, tmp, n);
    }
}

static int read_pieces(Arena_T arena, const char *path, long *gh, long *gw, Piece **out, size_t *n)
{
    const char *err = NULL;
    const JsonValue *j = Json_parse_file(arena, path, &err), *pl = NULL, *g = NULL;
    size_t i = 0, k = 0;
    Piece *p = NULL;
    if (j == NULL) { fprintf(stderr, "sheet_reading: %s: %s\n", path, err ? err : "parse error"); return -1; }
    g = Json_object_get(j, "grid_hw");
    pl = Json_object_get(j, "pieces");
    if (g == NULL || Json_array_len(g) != 2 || pl == NULL || Json_type(pl) != JSON_ARRAY) {
        fprintf(stderr, "sheet_reading: %s is not a sheet-handedness file\n", path);
        return -1;
    }
    *gh = Json_as_long(Json_array_get(g, 0), -1);
    *gw = Json_as_long(Json_array_get(g, 1), -1);
    *n = Json_array_len(pl);
    p = (Piece *)ARENA_CALLOC(arena, *n ? *n : 1, sizeof *p);
    for (i = 0; i < *n; i++) {
        const JsonValue *e = Json_array_get(pl, i);
        long v[6] = { 0, 0, 0, 0, 0, 0 };
        if (Json_array_len(e) != 6) { fprintf(stderr, "sheet_reading: bad piece %zu in %s\n", i, path); return -1; }
        for (k = 0; k < 6; k++) v[k] = Json_as_long(Json_array_get(e, k), -1);
        p[i].col0 = v[0]; p[i].col1 = v[1]; p[i].row0 = v[2]; p[i].row1 = v[3]; p[i].team = v[4]; p[i].mirror = v[5];
    }
    *out = p;
    return 0;
}

/* The mirror-imaged piece whose columns meet [c0, c1), or -1. */
static long mirrored_piece(const Piece *p, size_t n, long c0, long c1)
{
    size_t i = 0;
    for (i = 0; i < n; i++)
        if (p[i].col0 < c1 && p[i].col1 > c0 && p[i].mirror > p[i].team) return (long)i;
    return -1;
}

/* Flip the mirror-imaged pieces of img (H x W) in place.  ledger == NULL: img is on the sheet grid (gh x gw);
 * else img is that sheet's strip page.  *nflip = pieces (sheet) or fragments (strip) flipped.  0 = ok. */
static int reading_view(Arena_T arena, uint8_t *img, int H, int W, const Piece *p, size_t n, long gh, long gw,
                        const char *ledger, size_t *nflip, size_t *nfrag)
{
    uint8_t *tmp = (uint8_t *)ARENA_ALLOC(arena, (size_t)W);
    size_t i = 0;
    *nflip = 0; *nfrag = 0;
    if (ledger == NULL) {
        if (H != gh || W != gw) {
            fprintf(stderr, "sheet_reading: image %dx%d is not the sheet grid %ldx%ld\n", W, H, gw, gh);
            return -1;
        }
        for (i = 0; i < n; i++) {
            if (!(p[i].mirror > p[i].team)) continue;
            if (p[i].col0 < 0 || p[i].col1 > W || p[i].row0 < 0 || p[i].row1 > H || p[i].col0 >= p[i].col1) return -1;
            flip_rows(img, (size_t)W, p[i].col0, p[i].col1, p[i].row0, p[i].row1, tmp);
            (*nflip)++;
        }
        return 0;
    } else {
        const char *err = NULL;
        const JsonValue *j = Json_parse_file(arena, ledger, &err), *fl = NULL, *ll = NULL;
        long margin = 0, sw = 0, sh = 0;
        if (j == NULL) { fprintf(stderr, "sheet_reading: %s: %s\n", ledger, err ? err : "parse error"); return -1; }
        margin = Json_member_long(j, "margin", -1);
        sw = Json_member_long(j, "source_w", -1);
        sh = Json_member_long(j, "source_h", -1);
        fl = Json_object_get(j, "fragment_list");
        ll = Json_object_get(j, "line_list");
        if (margin < 0 || fl == NULL || ll == NULL || Json_member_long(j, "page_w", -1) != W ||
            Json_member_long(j, "page_h", -1) != H || sw != gw || sh != gh) {
            fprintf(stderr, "sheet_reading: ledger %s does not describe this page of this sheet\n", ledger);
            return -1;
        }
        *nfrag = Json_array_len(fl);
        for (i = 0; i < *nfrag; i++) {
            const JsonValue *f = Json_array_get(fl, i);
            long line = Json_as_long(Json_array_get(f, 1), -1), sx0 = Json_as_long(Json_array_get(f, 2), -1);
            long sx1 = Json_as_long(Json_array_get(f, 3), -1), sy0 = Json_as_long(Json_array_get(f, 4), -1);
            long sy1 = Json_as_long(Json_array_get(f, 5), -1), dx = Json_as_long(Json_array_get(f, 6), -1);
            long dy = 0, x0 = 0, y0 = 0;
            if (line < 0 || (size_t)line >= Json_array_len(ll) || sx0 < 0 || sx1 <= sx0 || sy0 < 0 || sy1 < sy0) return -1;
            if (mirrored_piece(p, n, sx0, sx1) < 0) continue;
            dy = Json_as_long(Json_array_get(Json_array_get(ll, (size_t)line), 0), -1);
            x0 = margin + dx; y0 = dy;
            if (dy < 0 || x0 < 0 || x0 + (sx1 - sx0) > W || y0 + (sy1 - sy0) > H) return -1;
            flip_rows(img, (size_t)W, x0, x0 + (sx1 - sx0), y0, y0 + (sy1 - sy0), tmp);
            (*nflip)++;
        }
        return 0;
    }
}

static int selftest(void)
{
    int fails = 0;
    Arena_T arena = Arena_new();
    const int W = 12, H = 8;
    uint8_t img[12 * 8], orig[12 * 8];
    Piece p[2];
    size_t nf = 0, nfr = 0;
    int x = 0, y = 0, bad = 0;
    long sum0 = 0, sum1 = 0;
    const char *ledger = "sheet_reading_selftest_ledger.json";
    FILE *fp = NULL;
    /* sheet grid: piece 0 cols [0,4) rows [1,7) team-handed; piece 1 cols [6,10) rows [0,8) mirror-imaged */
    for (y = 0; y < H; y++)
        for (x = 0; x < W; x++) { img[y * W + x] = (uint8_t)(10 * y + x); sum0 += img[y * W + x]; }
    memcpy(orig, img, sizeof img);
    p[0].col0 = 0; p[0].col1 = 4; p[0].row0 = 1; p[0].row1 = 7; p[0].team = 20; p[0].mirror = 1;
    p[1].col0 = 6; p[1].col1 = 10; p[1].row0 = 0; p[1].row1 = 8; p[1].team = 2; p[1].mirror = 30;
    if (reading_view(arena, img, H, W, p, 2, H, W, NULL, &nf, &nfr) != 0 || nf != 1) fails++;
    for (y = 0; y < H; y++)
        for (x = 0; x < W; x++) {
            int want = (x >= 6 && x < 10) ? orig[(H - 1 - y) * W + x] : orig[y * W + x];
            if (img[y * W + x] != want) bad++;
            sum1 += img[y * W + x];
        }
    if (bad || sum1 != sum0) { fprintf(stderr, "  FAIL: sheet-grid flip (%d wrong px)\n", bad); fails++; }
    /* strip page 20 x 10, margin 1, one line at dy 1: fragment of piece 0 (source cols 0-3, rows 1-6) at dx 0,
     * fragment of piece 1 (source cols 6-9, rows 0-7) at dx 6 -> page cols 7-10, rows 1-8 flip */
    fp = fopen(ledger, "wb");
    if (fp == NULL) fails++;
    else {
        uint8_t page[20 * 10], pg0[20 * 10];
        fprintf(fp, "{\"schema\":\"vesuvius-sheet-strip-v1\",\"source_w\":12,\"source_h\":8,\"page_w\":20,"
                    "\"page_h\":10,\"margin\":1,\"line_list\":[[1,8,0]],"
                    "\"fragment_list\":[[0,0,0,4,1,7,0],[1,0,6,10,0,8,6]]}\n");
        fclose(fp);
        for (y = 0; y < 10; y++) for (x = 0; x < 20; x++) page[y * 20 + x] = (uint8_t)(y * 20 + x);
        memcpy(pg0, page, sizeof page);
        if (reading_view(arena, page, 10, 20, p, 2, H, W, ledger, &nf, &nfr) != 0 || nf != 1 || nfr != 2) fails++;
        bad = 0;
        for (y = 0; y < 10; y++)
            for (x = 0; x < 20; x++) {
                int in = x >= 7 && x < 11 && y >= 1 && y < 9;
                int want = in ? pg0[(1 + 8 - 1 - (y - 1)) * 20 + x] : pg0[y * 20 + x];
                if (page[y * 20 + x] != want) bad++;
            }
        if (bad) { fprintf(stderr, "  FAIL: strip-page flip (%d wrong px)\n", bad); fails++; }
        /* a ledger of another page size is refused */
        if (reading_view(arena, page, 9, 20, p, 2, H, W, ledger, &nf, &nfr) == 0) {
            fprintf(stderr, "  FAIL: mismatched ledger accepted\n");
            fails++;
        }
        remove(ledger);
    }
    /* the sheet grid must match the handedness grid */
    if (reading_view(arena, img, H, W - 1, p, 2, H, W, NULL, &nf, &nfr) == 0) {
        fprintf(stderr, "  FAIL: mismatched grid accepted\n");
        fails++;
    }
    Arena_dispose(&arena);
    fprintf(stderr, "[selftest] sheet_reading %s (%d failures)\n", fails ? "FAIL" : "PASS", fails);
    return fails ? 1 : 0;
}

int main(int argc, char **argv)
{
    const char *in = NULL, *hand = NULL, *stem = NULL, *ledger = NULL;
    char path[2400];
    Arena_T arena = NULL;
    uint8_t *img = NULL;
    int H = 0, W = 0, i = 0, rc = 1, f = 0, png_f = 0;
    long gh = 0, gw = 0;
    Piece *p = NULL;
    size_t n = 0, nflip = 0, nfrag = 0, k = 0, mp = 0;
    double t0 = ves_clock_sec();
    FILE *fp = NULL;

    if (argc >= 2 && strcmp(argv[1], "--selftest") == 0) return selftest();
    if (argc < 4) {
        fprintf(stderr, "usage: sheet_reading <in.tif> <handedness.json> <out_stem> [--strip <ledger.json>]\n"
                        "       sheet_reading --selftest\n");
        return 1;
    }
    in = argv[1]; hand = argv[2]; stem = argv[3];
    for (i = 4; i < argc; i++) {
        if (strcmp(argv[i], "--strip") == 0 && i + 1 < argc) ledger = argv[++i];
        else { fprintf(stderr, "unknown option %s\n", argv[i]); return 1; }
    }
    arena = Arena_new();
    if (read_pieces(arena, hand, &gh, &gw, &p, &n) != 0 || load_gray_any(arena, in, &img, &H, &W) != 0) goto out;
    for (k = 0; k < n; k++) mp += p[k].mirror > p[k].team;
    if (reading_view(arena, img, H, W, p, n, gh, gw, ledger, &nflip, &nfrag) != 0) {
        fprintf(stderr, "ERROR: reading view failed (grid, ledger or piece bounds)\n");
        goto out;
    }
    snprintf(path, sizeof path, "%s.tif", stem);
    if (ves_ensure_parent_dir(path) != 0 || TiffIO_save(path, img, 1, H, W) != 0) {
        fprintf(stderr, "ERROR: cannot write %s\n", path);
        goto out;
    }
    for (f = 1; f <= 16 && png_f == 0; f *= 2) {       /* PNG at full size, or the smallest power-of-2 downsample */
        Arena_Mark mark = Arena_save(arena);
        int32_t ow = 0, oh = 0;
        uint8_t *small = (uint8_t *)ARENA_ALLOC(arena, (size_t)((W + f - 1) / f) * (size_t)((H + f - 1) / f));
        if (StripLayout_downsample(img, W, H, f, small, &ow, &oh) == 0) {
            snprintf(path, sizeof path, "%s.png", stem);
            if (VesPng_write_gray(path, small, ow, oh) == 0) png_f = f;
        }
        Arena_restore(arena, mark);
    }
    snprintf(path, sizeof path, "%s.json", stem);
    fp = fopen(path, "wb");
    if (fp == NULL) { fprintf(stderr, "ERROR: cannot write %s\n", path); goto out; }
    fprintf(fp, "{\"schema\":\"sheet-reading-v1\",\"input\":\"");
    for (k = 0; in[k]; k++) fprintf(fp, in[k] == '\\' ? "/" : "%c", in[k]);
    fprintf(fp, "\",\"mode\":\"%s\",\"pieces\":%zu,\"mirrored_pieces\":%zu,\"flipped\":%zu,\"fragments\":%zu,"
                "\"png_downsample\":%d,\"meaning\":\"mirror-imaged pieces flipped top to bottom; which way is up "
                "is not decided\"}\n", ledger ? "strip" : "sheet", n, mp, nflip, nfrag, png_f);
    fclose(fp);
    fprintf(stderr, "[sheet_reading] %s: %zu of %zu pieces mirror-imaged; flipped %zu %s; %dx%d, png 1/%d, %.1f s\n",
            stem, mp, n, nflip, ledger ? "fragments" : "pieces", W, H, png_f, ves_clock_sec() - t0);
    rc = 0;
out:
    if (arena != NULL) Arena_dispose(&arena);
    return rc;
}
