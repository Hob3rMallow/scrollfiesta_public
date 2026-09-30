/* sheet_strip -- strip page of a winding-ribbon bake: the ribbon flowed into
 * lines in spiral order, the kibble in a tray at the end (strip_layout.h).
 * Entry point only: read the bake, call StripLayout_*, write the page.
 *
 *   sheet_strip <layout_dir> <out_stem> [--texture <same_grid.tif>]
 *   sheet_strip --selftest
 *
 * --texture moves the pixels of another single-page 8-bit raster on the bake's grid (e.g. a reader-v2 ink
 * map, stripped or tiled) through the identical plan, which always comes from the bake's coverage; the ledger
 * then names the texture and otherwise equals the bake's page ledger.
 *
 * <layout_dir> is a review layout (qc_review) holding ribbon_layout.json and
 * bake/{sheet_rawtex.tif, sheet_rawtex_coverage.tif, sheet_rawtex_grid.json}.
 * Writes <out_stem>.tif (full resolution), <out_stem>.png (full resolution,
 * or the smallest integer downsample a PNG can hold), <out_stem>_preview.png
 * (longest side <= STRIP_PREVIEW_MAX_PX) and <out_stem>.json (inputs, counts,
 * the pixel-conservation check and the fragment ledger). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/arena.h"
#include "../common/json_read.h"
#include "../common/pipeline_constants.h"
#include "../common/strip_layout.h"
#include "../common/tiff_io.h"
#include "../common/ves_platform.h"
#include "../common/ves_png.h"

#include "tiffio.h"

static int load_gray(Arena_T arena, const char *path, uint8_t **img, int *H, int *W)
{
    int D = 0;
    if (TiffIO_load(arena, path, img, &D, H, W) != 0 || D != 1) {
        fprintf(stderr, "sheet_strip: cannot load a single-page 8-bit TIFF %s\n", path);
        return -1;
    }
    return 0;
}

/* --texture rasters (e.g. an ink map inferred on the bake's grid) may be TILED, which TiffIO_load's strip
 * reader refuses: single-page, 8-bit, 1 sample, any libtiff codec, stripped or tiled. */
static int load_gray_any(Arena_T arena, const char *path, uint8_t **img, int *H, int *W)
{
    TIFF *tif = TIFFOpen(path, "r");
    uint32_t w = 0, h = 0, tw = 0, th = 0;
    uint16_t bps = 0, spp = 0;
    if (tif == NULL) {
        fprintf(stderr, "sheet_strip: cannot open %s\n", path);
        return -1;
    }
    if (!TIFFIsTiled(tif)) {
        TIFFClose(tif);
        return load_gray(arena, path, img, H, W);
    }
    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
    TIFFGetFieldDefaulted(tif, TIFFTAG_BITSPERSAMPLE, &bps);
    TIFFGetFieldDefaulted(tif, TIFFTAG_SAMPLESPERPIXEL, &spp);
    TIFFGetField(tif, TIFFTAG_TILEWIDTH, &tw);
    TIFFGetField(tif, TIFFTAG_TILELENGTH, &th);
    if (w == 0 || h == 0 || w > INT32_MAX || h > INT32_MAX || bps != 8 || spp != 1 || tw == 0 || th == 0 ||
        TIFFNumberOfDirectories(tif) != 1) {
        fprintf(stderr, "sheet_strip: %s is not a single-page 8-bit gray TIFF\n", path);
        TIFFClose(tif);
        return -1;
    }
    uint8_t *out = ARENA_CALLOC(arena, (size_t)w * (size_t)h, 1);
    uint8_t *tile = ARENA_ALLOC(arena, (size_t)TIFFTileSize(tif));
    for (uint32_t y = 0; y < h; y += th)
        for (uint32_t x = 0; x < w; x += tw) {
            if (TIFFReadTile(tif, tile, x, y, 0, 0) < 0) {
                fprintf(stderr, "sheet_strip: tile read failed in %s at (%u, %u)\n", path, x, y);
                TIFFClose(tif);
                return -1;
            }
            uint32_t cw = w - x < tw ? w - x : tw, ch = h - y < th ? h - y : th;
            for (uint32_t r = 0; r < ch; r++)
                memcpy(out + (size_t)(y + r) * w + x, tile + (size_t)r * tw, cw);
        }
    TIFFClose(tif);
    *img = out;
    *H = (int)h;
    *W = (int)w;
    return 0;
}

/* Tiled round trip for load_gray_any: a 40 x 21 ramp in 16 x 16 tiles (partial edge tiles). 0 = pass. */
static int texture_selftest(void)
{
    const char *path = "sheet_strip_texture_selftest.tif";
    const uint32_t w = 40, h = 21, tsz = 16;
    int fails = 0;
    TIFF *tif = TIFFOpen(path, "w");
    if (tif == NULL) return 1;
    TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, w);
    TIFFSetField(tif, TIFFTAG_IMAGELENGTH, h);
    TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 8);
    TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);
    TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
    TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_LZW);
    TIFFSetField(tif, TIFFTAG_TILEWIDTH, tsz);
    TIFFSetField(tif, TIFFTAG_TILELENGTH, tsz);
    uint8_t tile[16 * 16];
    for (uint32_t y = 0; y < h; y += tsz)
        for (uint32_t x = 0; x < w; x += tsz) {
            for (uint32_t r = 0; r < tsz; r++)
                for (uint32_t c = 0; c < tsz; c++)
                    tile[r * tsz + c] = (uint8_t)((y + r) < h && (x + c) < w ? ((y + r) * 7 + (x + c) * 3) & 255 : 0);
            if (TIFFWriteTile(tif, tile, x, y, 0, 0) < 0) fails++;
        }
    TIFFClose(tif);
    Arena_T arena = Arena_new();
    uint8_t *img = NULL;
    int H = 0, W = 0;
    if (load_gray_any(arena, path, &img, &H, &W) != 0 || H != (int)h || W != (int)w) fails++;
    else
        for (uint32_t y = 0; y < h; y++)
            for (uint32_t x = 0; x < w; x++)
                if (img[(size_t)y * w + x] != (uint8_t)((y * 7 + x * 3) & 255)) fails++;
    Arena_dispose(&arena);
    remove(path);
    fprintf(stderr, "[selftest] sheet_strip tiled texture load %s (%d failures)\n", fails ? "FAIL" : "PASS", fails);
    return fails ? 1 : 0;
}

/* texture == NULL: the page of the bake itself.  Otherwise the pixels of `texture` (a raster on the bake's
 * grid, e.g. an ink map) move through the identical plan -- the plan always comes from the bake's coverage
 * -- and the ledger records the texture. */
static int run(const char *layout, const char *stem, const char *texture)
{
    double t0 = ves_clock_sec();
    Arena_T arena = Arena_new();
    char path[4096];
    int rc = 1;
    const char *err = NULL;

    snprintf(path, sizeof path, "%s/ribbon_layout.json", layout);
    const JsonValue *ribbon = Json_parse_file(arena, path, &err);
    if (ribbon == NULL || !Json_as_bool(Json_object_get(ribbon, "kept"), 0)) {
        fprintf(stderr, "sheet_strip: %s is not a kept winding ribbon (%s)\n", path, err ? err : "kept false");
        goto done;
    }
    long regions = Json_member_long(ribbon, "regions", -1);
    double ribbon_gap = Json_member_double(ribbon, "gap", ASM_READING_RIBBON_GAP);
    snprintf(path, sizeof path, "%s/bake/sheet_rawtex_grid.json", layout);
    const JsonValue *grid = Json_parse_file(arena, path, &err);
    const JsonValue *step = grid ? Json_object_get(grid, "step_uv") : NULL;
    if (step == NULL || Json_array_len(step) != 2) {
        fprintf(stderr, "sheet_strip: %s has no step_uv\n", path);
        goto done;
    }
    double du = Json_as_double(Json_array_get(step, 0), 0.0), dv = Json_as_double(Json_array_get(step, 1), 0.0);

    uint8_t *tex = NULL, *cov = NULL;
    int H = 0, W = 0, Hc = 0, Wc = 0;
    if (texture != NULL) {
        if (load_gray_any(arena, texture, &tex, &H, &W)) goto done;
    } else {
        snprintf(path, sizeof path, "%s/bake/sheet_rawtex.tif", layout);
        if (load_gray(arena, path, &tex, &H, &W)) goto done;
    }
    snprintf(path, sizeof path, "%s/bake/sheet_rawtex_coverage.tif", layout);
    if (load_gray(arena, path, &cov, &Hc, &Wc)) goto done;
    if (H != Hc || W != Wc) {
        fprintf(stderr, "sheet_strip: texture %dx%d and coverage %dx%d differ\n", W, H, Wc, Hc);
        goto done;
    }
    double t_load = ves_clock_sec();

    StripParams p;
    StripLayout_params_default(&p, du, dv);
    if (ribbon_gap > 0.0) p.merge_gap_vox = 0.5 * ribbon_gap;
    StripColumns cols = {0};
    StripPage pg = {0};
    if (StripLayout_columns(arena, cov, W, H, &cols) || StripLayout_plan(arena, &cols, &p, &pg)) {
        fprintf(stderr, "sheet_strip: no layout (empty coverage?)\n");
        goto done;
    }
    int64_t covered = 0;
    for (int32_t x = 0; x < W; x++) covered += cols.count[x];
    double t_plan = ves_clock_sec();

    uint8_t *page = ARENA_CALLOC(arena, (size_t)pg.W * (size_t)pg.H, 1);
    int64_t copied = 0;
    if (StripLayout_render(&pg, tex, cov, W, H, page, &copied) || copied != covered) {
        fprintf(stderr, "sheet_strip: pixel conservation failed: %lld covered, %lld copied\n",
                (long long)covered, (long long)copied);
        goto done;
    }
    double t_render = ves_clock_sec();

    snprintf(path, sizeof path, "%s.tif", stem);
    if (TiffIO_save(path, page, 1, pg.H, pg.W)) {
        fprintf(stderr, "sheet_strip: cannot write %s\n", path);
        goto done;
    }
    /* PNG at full resolution, or the smallest downsample the PNG writer holds. */
    int png_factor = 0;
    for (int32_t f = 1; f <= 8 && png_factor == 0; f++) {
        const uint8_t *img = page;
        int32_t ow = pg.W, oh = pg.H;
        Arena_Mark mark = Arena_save(arena);
        if (f > 1) {
            uint8_t *small = ARENA_ALLOC(arena, (size_t)((pg.W + f - 1) / f) * (size_t)((pg.H + f - 1) / f));
            StripLayout_downsample(page, pg.W, pg.H, f, small, &ow, &oh);
            img = small;
        }
        snprintf(path, sizeof path, "%s.png", stem);
        if (VesPng_write_gray(path, img, ow, oh) == 0) png_factor = (int)f;
        Arena_restore(arena, mark);
    }
    int32_t longest = pg.W > pg.H ? pg.W : pg.H;
    int32_t pf = (longest + STRIP_PREVIEW_MAX_PX - 1) / STRIP_PREVIEW_MAX_PX, pw = 0, ph = 0;
    if (pf < 1) pf = 1;
    uint8_t *preview = ARENA_ALLOC(arena, (size_t)((pg.W + pf - 1) / pf) * (size_t)((pg.H + pf - 1) / pf));
    StripLayout_downsample(page, pg.W, pg.H, pf, preview, &pw, &ph);
    snprintf(path, sizeof path, "%s_preview.png", stem);
    int preview_ok = VesPng_write_gray(path, preview, pw, ph) == 0;
    double t_write = ves_clock_sec();

    snprintf(path, sizeof path, "%s.json", stem);
    FILE *f = fopen(path, "wb");
    if (f == NULL) {
        fprintf(stderr, "sheet_strip: cannot write %s\n", path);
        goto done;
    }
    int ok = fprintf(f, "{\"schema\":\"vesuvius-sheet-strip-v1\",\"layout\":\"") > 0;
    for (const char *s = layout; *s; s++) {                 /* JSON string: '/' separators, quotes escaped */
        if (*s == '"') ok &= fputs("\\\"", f) >= 0;
        else ok &= fputc(*s == '\\' ? '/' : *s, f) != EOF;
    }
    if (texture != NULL) {                                   /* default ledgers stay byte-identical */
        ok &= fputs("\",\"texture\":\"", f) >= 0;
        for (const char *s = texture; *s; s++) {
            if (*s == '"') ok &= fputs("\\\"", f) >= 0;
            else ok &= fputc(*s == '\\' ? '/' : *s, f) != EOF;
        }
    }
    ok &= fprintf(f, "\",\"source_w\":%d,\"source_h\":%d,\"ribbon_regions\":%ld,\"ribbon_gap_vox\":%.17g,"
                     "\"segments_match_ribbon_regions\":%s,\"covered_px\":%lld,\"copied_px\":%lld,"
                     "\"pixels_conserved\":%s,\"png_downsample\":%d,\"preview_downsample\":%d,\"preview_written\":%s,"
                     "\"seconds\":{\"load\":%.1f,\"plan\":%.1f,\"render\":%.1f,\"write\":%.1f},\n",
                  W, H, regions, ribbon_gap, regions == (long)pg.n_pieces ? "true" : "false",
                  (long long)covered, (long long)copied, copied == covered ? "true" : "false", png_factor, (int)pf,
                  preview_ok ? "true" : "false", t_load - t0, t_plan - t_load, t_render - t_plan,
                  t_write - t_render) > 0;
    ok &= StripLayout_write_ledger(f, &pg, &p) == 0;
    ok &= fputs("}\n", f) >= 0;
    ok &= fclose(f) == 0;
    if (!ok) {
        fprintf(stderr, "sheet_strip: ledger write failed\n");
        goto done;
    }
    fprintf(stderr, "[sheet_strip] %s: %zu pieces (ribbon %ld), %zu main + %zu tray (%.1f%% of the area); "
                    "%zu lines, %zu splits; page %dx%d px (line %d); png 1/%d; %.1f s\n",
            stem, pg.n_pieces, regions, pg.n_pieces - pg.n_tray, pg.n_tray,
            100.0 * pg.tray_area_vox2 / (pg.main_area_vox2 + pg.tray_area_vox2), pg.n_lines, pg.n_split,
            (int)pg.W, (int)pg.H, (int)pg.flow_width, png_factor, ves_clock_sec() - t0);
    if (regions != (long)pg.n_pieces)
        fprintf(stderr, "[sheet_strip] WARNING: %zu column runs but %ld ribbon regions\n", pg.n_pieces, regions);
    rc = png_factor ? 0 : 1;
done:
    Arena_dispose(&arena);
    return rc;
}

int main(int argc, char **argv)
{
    if (argc == 2 && strcmp(argv[1], "--selftest") == 0)
        return StripLayout_selftest() != 0 || texture_selftest() != 0;
    if (argc == 5 && strcmp(argv[3], "--texture") == 0) return run(argv[1], argv[2], argv[4]);
    if (argc != 3) {
        fprintf(stderr, "usage: sheet_strip <layout_dir> <out_stem> [--texture <same_grid.tif>]\n"
                        "       sheet_strip --selftest\n");
        return 2;
    }
    return run(argv[1], argv[2], NULL);
}
