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
#include "../common/block_list.h"
#include "../common/closeup_window.h"
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

/* --plan (2026-10-13, user: "get rid of the big CT TIFF and just do the annotator strips first"): the page PLAN alone
 * of a sheet exported in --window blocks -- each block's sheet-space mask (sheet_layers --tifxyz *_sheet_mask.tif) at
 * its first column of the full H x W sheet grid -- with no texture, no whole-sheet raster and no page.  The column
 * profile is per column, so the blocks' profiles stitched side by side ARE the profile of the composed coverage the
 * --texture run reads (region_ink_all strip_layout: the composed sheet mask on a 1 px/vox grid), and the plan, a pure
 * function of that profile and the ribbon, is the one --texture draws.  The ledger carries the page ledger's members
 * with "plan_only": true and no page (copied_px 0, pixels_conserved false: page publishers refuse it). */

/* A block's column profile placed at col0 of the full profile. 0 = ok; -1 = a height mismatch, an overhang or a column
 * another block already holds (blocks are cut in empty columns and never overlap). */
static int stitch_columns(StripColumns *full, const StripColumns *blk, int32_t col0, uint8_t *used)
{
    if (blk->H > full->H || col0 < 0 || (int64_t)col0 + blk->W > full->W) return -1;
    for (int32_t x = 0; x < blk->W; x++)
        if (used[col0 + x]) return -1;
    for (int32_t x = 0; x < blk->W; x++) {
        used[col0 + x] = 1;
        full->top[col0 + x] = blk->count[x] ? blk->top[x] : full->H;
        full->bottom[col0 + x] = blk->bottom[x];
        full->count[col0 + x] = blk->count[x];
    }
    return 0;
}

/* An empty H x W profile: every column uncovered (top H, bottom 0), as a zero-filled composed raster reads. */
static int empty_columns(Arena_T arena, int32_t W, int32_t H, StripColumns *out)
{
    memset(out, 0, sizeof *out);
    out->top = ARENA_ALLOC(arena, (size_t)W * sizeof *out->top);
    out->bottom = ARENA_ALLOC(arena, (size_t)W * sizeof *out->bottom);
    out->count = ARENA_ALLOC(arena, (size_t)W * sizeof *out->count);
    if (out->top == NULL || out->bottom == NULL || out->count == NULL) return -1;
    for (int32_t x = 0; x < W; x++) { out->top[x] = H; out->bottom[x] = 0; out->count[x] = 0; }
    out->W = W; out->H = H;
    return 0;
}

/* --repack (2026-10-24, the user: "specify that a given chart in a sheet is the wrong orientation and must be rotated
 * 90 degrees (and then all charts in the sheet should shuffle to make room for it)"): --plan with some pieces turned.
 * turns.json = {"turns": [[piece, rot, flip], ...]} (piece = the plan's piece index, as the lines' fragments name
 * it; rot 0/90/180/270 clockwise after a left-right flip). The same blocks give the same pieces; a turned piece's
 * page profile comes from its sub-mask read out of the blocks (StripLayout_columns_turned), and the ledger carries
 * piece_turns / fragment_turns (readerv2_annotate.fragments reads them). With no turns file this is --plan. */
static int read_turns(Arena_T arena, const char *path, size_t n, StripTurn *turns, int *n_turned)
{
    const char *err = NULL;
    const JsonValue *j = Json_parse_file(arena, path, &err);
    const JsonValue *arr = j ? Json_object_get(j, "turns") : NULL;
    if (arr == NULL || Json_type(arr) != JSON_ARRAY) {
        fprintf(stderr, "sheet_strip: %s: no \"turns\" list (%s)\n", path, err ? err : "not an object");
        return -1;
    }
    *n_turned = 0;
    for (size_t i = 0; i < Json_array_len(arr); i++) {
        const JsonValue *e = Json_array_get(arr, i);
        if (Json_type(e) != JSON_ARRAY || Json_array_len(e) < 2) { fprintf(stderr, "sheet_strip: %s: turn %zu is not [piece, rot, flip]\n", path, i); return -1; }
        long piece = Json_as_long(Json_array_get(e, 0), -1), rot = Json_as_long(Json_array_get(e, 1), 0);
        long flip = Json_array_len(e) > 2 ? (Json_type(Json_array_get(e, 2)) == JSON_BOOL ? Json_as_bool(Json_array_get(e, 2), 0)
                                                                                           : Json_as_long(Json_array_get(e, 2), 0)) : 0;
        rot = ((rot % 360) + 360) % 360;
        if (piece < 0 || (size_t)piece >= n || rot % 90 || (flip != 0 && flip != 1)) {
            fprintf(stderr, "sheet_strip: %s: turn %zu names piece %ld of %zu, rot %ld, flip %ld: not a turn\n", path, i, piece, n, rot, flip);
            return -1;
        }
        turns[piece].rot = (int32_t)rot;
        turns[piece].flip = (int32_t)flip;
        *n_turned += !StripTurn_identity(turns[piece]);
    }
    return 0;
}

static int plan_run(const char *layout, const char *stem, int32_t H, int32_t W, int nb, char **blocks, const char *turns_path)
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
    StripColumns cols;
    uint8_t *used = ARENA_CALLOC(arena, (size_t)W, 1);
    if (W <= 0 || H <= 0 || used == NULL || empty_columns(arena, W, H, &cols)) {
        fprintf(stderr, "sheet_strip: bad sheet grid %d x %d\n", (int)W, (int)H);
        goto done;
    }
    for (int i = 0; i < nb; i++) {
        Arena_Mark mark = Arena_save(arena);
        char *end = NULL;
        long col0 = strtol(blocks[2 * i], &end, 10);
        uint8_t *m = NULL;
        int bh = 0, bw = 0;
        StripColumns bc;
        if (end == blocks[2 * i] || *end != '\0' || load_gray(arena, blocks[2 * i + 1], &m, &bh, &bw) ||
            StripLayout_columns(arena, m, bw, bh, &bc) || stitch_columns(&cols, &bc, (int32_t)col0, used)) {
            fprintf(stderr, "sheet_strip: block %s at column %s does not fit the %d x %d sheet grid "
                            "(height, overhang or overlap)\n", blocks[2 * i + 1], blocks[2 * i], (int)W, (int)H);
            goto done;
        }
        Arena_restore(arena, mark);                     /* the block raster and its profile: copied into cols */
    }
    double t_load = ves_clock_sec();
    StripParams p;
    StripLayout_params_default(&p, 1.0, 1.0);           /* block masks sit on the 1 px/vox sheet grid */
    const char *pixel_budget = getenv("VESUVIUS_STRIP_LINE_PIXELS");
    if (pixel_budget && *pixel_budget) {
        char *end = NULL;
        p.max_line_pixels = strtoll(pixel_budget, &end, 10);
        if (end == pixel_budget || *end || p.max_line_pixels <= 0) {
            fprintf(stderr, "sheet_strip: VESUVIUS_STRIP_LINE_PIXELS must be a positive integer\n");
            goto done;
        }
    }
    if (ribbon_gap > 0.0) p.merge_gap_vox = 0.5 * ribbon_gap;
    StripPage pg = {0};
    StripPiece *pieces = NULL;
    size_t n_pieces = 0;
    if (StripLayout_segments(arena, &cols, &p, &pieces, &n_pieces)) {
        fprintf(stderr, "sheet_strip: no layout (empty coverage?)\n");
        goto done;
    }
    StripTurn *turns = NULL;
    const StripColumns **prof = NULL;
    int n_turned = 0;
    if (turns_path != NULL) {
        turns = ARENA_CALLOC(arena, n_pieces, sizeof *turns);
        prof = ARENA_CALLOC(arena, n_pieces, sizeof *prof);
        if (turns == NULL || prof == NULL || read_turns(arena, turns_path, n_pieces, turns, &n_turned)) goto done;
        /* the turned pieces' sub-masks, gathered block by block (each block read once more), then their page profiles */
        uint8_t **sub = ARENA_CALLOC(arena, n_pieces, sizeof *sub);
        if (sub == NULL) goto done;
        for (size_t i = 0; i < n_pieces; i++) {
            if (StripTurn_identity(turns[i])) continue;
            const int64_t area = (int64_t)(pieces[i].x1 - pieces[i].x0) * (int64_t)(pieces[i].y1 - pieces[i].y0);
            if (area <= 0 || (sub[i] = ARENA_CALLOC(arena, (size_t)area, 1)) == NULL) {
                fprintf(stderr, "sheet_strip: piece %zu: cannot hold its %lld px sub-mask\n", i, (long long)area);
                goto done;
            }
        }
        for (int b = 0; b < nb; b++) {
            Arena_Mark mark = Arena_save(arena);
            const long col0 = strtol(blocks[2 * b], NULL, 10);
            uint8_t *m = NULL;
            int bh = 0, bw = 0;
            if (n_turned > 0 && load_gray(arena, blocks[2 * b + 1], &m, &bh, &bw)) goto done;
            for (size_t i = 0; m != NULL && i < n_pieces; i++) {
                if (sub[i] == NULL) continue;
                const StripPiece *pc = pieces + i;
                const int32_t x0 = pc->x0 > (int32_t)col0 ? pc->x0 : (int32_t)col0;
                const int32_t x1 = pc->x1 < (int32_t)col0 + bw ? pc->x1 : (int32_t)col0 + bw;
                if (x0 >= x1) continue;
                for (int32_t y = pc->y0; y < pc->y1 && y < bh; y++)
                    memcpy(sub[i] + (size_t)(y - pc->y0) * (size_t)(pc->x1 - pc->x0) + (size_t)(x0 - pc->x0),
                           m + (size_t)y * (size_t)bw + (size_t)(x0 - col0), (size_t)(x1 - x0));
            }
            Arena_restore(arena, mark);
        }
        for (size_t i = 0; i < n_pieces; i++) {
            if (sub[i] == NULL) continue;
            StripColumns *c = ARENA_CALLOC(arena, 1, sizeof *c);
            if (c == NULL || StripLayout_columns_turned(arena, sub[i], pieces[i].x1 - pieces[i].x0, pieces[i].y1 - pieces[i].y0,
                                                        turns[i], c)) {
                fprintf(stderr, "sheet_strip: piece %zu: no turned profile\n", i);
                goto done;
            }
            prof[i] = c;
        }
    }
    if (StripLayout_plan_pieces(arena, &cols, pieces, n_pieces, turns, prof, &p, &pg)) {
        fprintf(stderr, "sheet_strip: no layout (%s)\n", turns_path ? "a turned piece's profile" : "empty coverage?");
        goto done;
    }
    int64_t covered = 0;
    for (int32_t x = 0; x < W; x++) covered += cols.count[x];
    double t_plan = ves_clock_sec();
    snprintf(path, sizeof path, "%s.json", stem);
    FILE *f = fopen(path, "wb");
    if (f == NULL) {
        fprintf(stderr, "sheet_strip: cannot write %s\n", path);
        goto done;
    }
    int ok = fprintf(f, "{\"schema\":\"vesuvius-sheet-strip-v1\",\"plan_only\":true,\"layout\":\"") > 0;
    for (const char *s = layout; *s; s++) {
        if (*s == '"') ok &= fputs("\\\"", f) >= 0;
        else ok &= fputc(*s == '\\' ? '/' : *s, f) != EOF;
    }
    ok &= fprintf(f, "\",\"blocks\":%d,\"source_w\":%d,\"source_h\":%d,\"ribbon_regions\":%ld,\"ribbon_gap_vox\":%.17g,"
                     "\"segments_match_ribbon_regions\":%s,\"covered_px\":%lld,\"copied_px\":0,"
                     "\"pixels_conserved\":false,\"png_downsample\":0,\"preview_downsample\":0,\"preview_written\":false,"
                     "\"seconds\":{\"load\":%.1f,\"plan\":%.1f,\"render\":0.0,\"write\":0.0},\n",
                  nb, (int)W, (int)H, regions, ribbon_gap, regions == (long)pg.n_pieces ? "true" : "false",
                  (long long)covered, t_load - t0, t_plan - t_load) > 0;
    ok &= StripLayout_write_ledger(f, &pg, &p) == 0;
    ok &= fputs("}\n", f) >= 0;
    ok &= fclose(f) == 0;
    if (!ok) {
        fprintf(stderr, "sheet_strip: ledger write failed\n");
        goto done;
    }
    fprintf(stderr, "[sheet_strip] %s: PLAN ONLY from %d blocks: %zu pieces (ribbon %ld), %zu main + %zu tray; %zu lines, "
                    "%zu splits; page %dx%d px; %d turned; %.1f s\n", stem, nb, pg.n_pieces, regions, pg.n_pieces - pg.n_tray,
            pg.n_tray, pg.n_lines, pg.n_split, (int)pg.W, (int)pg.H, n_turned, ves_clock_sec() - t0);
    rc = 0;
done:
    Arena_dispose(&arena);
    return rc;
}

/* --plan's premise, in memory: a coverage cut into column blocks gives, stitched, the whole raster's column profile
 * and so the same plan; overlaps, overhangs and height mismatches are refused. 0 = pass. */
static int plan_selftest(void)
{
    const int32_t W = 90, H = 12, cut = 37;
    int fails = 0;
    Arena_T arena = Arena_new();
    uint8_t *cov = ARENA_CALLOC(arena, (size_t)W * (size_t)H, 1);
    for (int32_t x = 0; x < W; x++)                        /* three pieces of varying height, gaps between them */
        for (int32_t y = 0; y < H; y++)
            if ((x >= 3 && x < 30 && y >= 2 && y < 9 - (x % 3)) || (x >= 41 && x < 66 && y >= 1 + (x % 4) && y < 11) ||
                (x >= 70 && x < 74 && y >= 4 && y < 6))
                cov[(size_t)y * (size_t)W + (size_t)x] = 255;
    StripColumns whole, a, b, stitched;
    uint8_t *left = ARENA_ALLOC(arena, (size_t)cut * (size_t)H), *right = ARENA_ALLOC(arena, (size_t)(W - cut) * (size_t)H);
    for (int32_t y = 0; y < H; y++) {
        memcpy(left + (size_t)y * (size_t)cut, cov + (size_t)y * (size_t)W, (size_t)cut);
        memcpy(right + (size_t)y * (size_t)(W - cut), cov + (size_t)y * (size_t)W + (size_t)cut, (size_t)(W - cut));
    }
    uint8_t *used = ARENA_CALLOC(arena, (size_t)W, 1);
    if (StripLayout_columns(arena, cov, W, H, &whole) || StripLayout_columns(arena, left, cut, H, &a) ||
        StripLayout_columns(arena, right, W - cut, H, &b) || empty_columns(arena, W, H, &stitched) ||
        stitch_columns(&stitched, &b, cut, used) || stitch_columns(&stitched, &a, 0, used))   /* any block order */
        fails++;
    else {
        for (int32_t x = 0; x < W; x++)
            fails += whole.top[x] != stitched.top[x] || whole.bottom[x] != stitched.bottom[x] ||
                     whole.count[x] != stitched.count[x];
        StripParams p;
        StripPage p1 = {0}, p2 = {0};
        StripLayout_params_default(&p, 1.0, 1.0);
        if (StripLayout_plan(arena, &whole, &p, &p1) || StripLayout_plan(arena, &stitched, &p, &p2)) fails++;
        else {
            fails += p1.n_pieces != p2.n_pieces || p1.n_frags != p2.n_frags || p1.n_lines != p2.n_lines ||
                     p1.W != p2.W || p1.H != p2.H || p1.n_tray != p2.n_tray;
            for (size_t i = 0; i < p1.n_frags && i < p2.n_frags; i++)
                fails += memcmp(&p1.frags[i], &p2.frags[i], sizeof p1.frags[i]) != 0;
        }
        fails += stitch_columns(&stitched, &a, 0, used) != -1;                 /* the same columns twice */
        StripColumns tall = a;
        tall.H = H + 1;
        fails += stitch_columns(&stitched, &tall, 0, used) != -1;               /* another height */
        uint8_t *fresh = ARENA_CALLOC(arena, (size_t)W, 1);
        fails += stitch_columns(&stitched, &b, cut + 1, fresh) != -1;          /* past the sheet's last column */
        /* Trailing zero rows need no storage, including empty columns whose
         * sentinel top must still equal the logical sheet height. */
        StripColumns short_a;
        memset(fresh, 0, (size_t)W);
        if (StripLayout_columns(arena, left, cut, 9, &short_a) ||
            stitch_columns(&stitched, &short_a, 0, fresh)) fails++;
        else for (int32_t x = 0; x < cut; x++)
            fails += whole.top[x] != stitched.top[x] || whole.bottom[x] != stitched.bottom[x] ||
                     whole.count[x] != stitched.count[x];
        /* A bounded line layout preserves every covered source pixel once,
         * even when the budget requires lines narrower than the old minimum. */
        for (int budget = 30; budget <= 120; budget *= 2) {
            StripPage bounded = {0}; p.max_line_pixels = budget;
            p.kibble_area_vox2 = 0;
            uint8_t *seen = ARENA_CALLOC(arena, (size_t)W*H, 1);
            if (StripLayout_plan(arena, &whole, &p, &bounded)) { fails++; continue; }
            for (size_t l = 0; l < bounded.n_lines; l++)
                fails += (int64_t)bounded.flow_width * bounded.lines[l].h > budget;
            for (size_t i = 0; i < bounded.n_frags; i++) {
                const StripFragment *g = bounded.frags+i;
                for (int y = g->sy0; y < g->sy1; y++)
                    for (int x = g->sx0; x < g->sx1; x++)
                        if (cov[y*W+x]) seen[y*W+x]++;
            }
            for (int i = 0; i < W*H; i++) fails += seen[i] != (cov[i] != 0);
        }
    }
    Arena_dispose(&arena);
    fprintf(stderr, "[selftest] sheet_strip --plan from blocks %s (%d failures)\n", fails ? "FAIL" : "PASS", fails);
    return fails ? 1 : 0;
}

/* --repack end to end (2026-10-24): the plan_selftest coverage as two block masks on disk, piece 1 turned 90 degrees
 * through a turns.json, the ledger read back: piece_turns / fragment_turns present, the turned piece's fragments
 * its height wide and its rows its columns, the other pieces as --plan lays them. 0 = pass. */
static int repack_selftest(void)
{
    const int32_t W = 90, H = 12, cut = 37;
    int fails = 0;
    char dir[1024], left[1200], right[1200], turns[1200], ribbon[1200], stem[1200], plan[1200];
    if (ves_selftest_dir(dir, sizeof dir, "sheet_strip_repack") != 0) { fprintf(stderr, "[selftest] sheet_strip --repack: no scratch dir\n"); return 1; }
    Arena_T arena = Arena_new();
    uint8_t *cov = ARENA_CALLOC(arena, (size_t)W * (size_t)H, 1);
    for (int32_t x = 0; x < W; x++)
        for (int32_t y = 0; y < H; y++)
            if ((x >= 3 && x < 30 && y >= 2 && y < 9 - (x % 3)) || (x >= 41 && x < 66 && y >= 1 + (x % 4) && y < 11) ||
                (x >= 70 && x < 74 && y >= 4 && y < 6))
                cov[(size_t)y * (size_t)W + (size_t)x] = 255;
    uint8_t *l = ARENA_ALLOC(arena, (size_t)cut * (size_t)H), *r = ARENA_ALLOC(arena, (size_t)(W - cut) * (size_t)H);
    for (int32_t y = 0; y < H; y++) {
        memcpy(l + (size_t)y * (size_t)cut, cov + (size_t)y * (size_t)W, (size_t)cut);
        memcpy(r + (size_t)y * (size_t)(W - cut), cov + (size_t)y * (size_t)W + (size_t)cut, (size_t)(W - cut));
    }
    snprintf(left, sizeof left, "%s/left.tif", dir);
    snprintf(right, sizeof right, "%s/right.tif", dir);
    snprintf(turns, sizeof turns, "%s/turns.json", dir);
    snprintf(ribbon, sizeof ribbon, "%s/ribbon_layout.json", dir);
    snprintf(stem, sizeof stem, "%s/plan", dir);
    snprintf(plan, sizeof plan, "%s/plan.json", dir);
    FILE *f = fopen(turns, "wb");
    int ok = f != NULL && fputs("{\"turns\": [[1, 90, 0]]}\n", f) >= 0 && fclose(f) == 0;
    f = fopen(ribbon, "wb");
    ok = ok && f != NULL && fputs("{\"kept\": true, \"gap\": 8, \"regions\": 3}\n", f) >= 0 && fclose(f) == 0;
    ok = ok && TiffIO_save(left, l, 1, H, cut) == 0 && TiffIO_save(right, r, 1, H, W - cut) == 0;
    char col_l[8] = "0", col_r[8];
    snprintf(col_r, sizeof col_r, "%d", (int)cut);
    char *blocks[4] = {col_l, left, col_r, right};
    if (!ok || plan_run(dir, stem, H, W, 2, blocks, turns) != 0) fails++;
    else {
        const char *err = NULL;
        const JsonValue *j = Json_parse_file(arena, plan, &err);
        const JsonValue *pt = j ? Json_object_get(j, "piece_turns") : NULL, *ft = j ? Json_object_get(j, "fragment_turns") : NULL;
        const JsonValue *fl = j ? Json_object_get(j, "fragment_list") : NULL, *pl = j ? Json_object_get(j, "piece_list") : NULL;
        if (pt == NULL || ft == NULL || fl == NULL || pl == NULL || Json_array_len(pt) != 3 || Json_array_len(ft) != Json_array_len(fl)) fails++;
        else {
            fails += Json_as_long(Json_array_get(Json_array_get(pt, 1), 0), -1) != 90 || Json_as_long(Json_array_get(Json_array_get(pt, 0), 0), -1) != 0;
            const JsonValue *B = Json_array_get(pl, 1);                  /* [x0, x1, y0, y1, covered, area, tray] */
            long bx0 = Json_as_long(Json_array_get(B, 0), -1), bx1 = Json_as_long(Json_array_get(B, 1), -1);
            long by0 = Json_as_long(Json_array_get(B, 2), -1), by1 = Json_as_long(Json_array_get(B, 3), -1);
            int nb = 0;
            for (size_t i = 0; i < Json_array_len(fl); i++) {
                const JsonValue *g = Json_array_get(fl, i), *t = Json_array_get(ft, i);
                long piece = Json_as_long(Json_array_get(g, 0), -1), sx0 = Json_as_long(Json_array_get(g, 2), -1);
                long sx1 = Json_as_long(Json_array_get(g, 3), -1), sy0 = Json_as_long(Json_array_get(g, 4), -1), sy1 = Json_as_long(Json_array_get(g, 5), -1);
                long rot = Json_as_long(Json_array_get(t, 0), -1);
                if (piece == 1) {                                       /* the whole piece (11 columns wide on the page) on one line */
                    nb++;
                    fails += rot != 90 || sy0 != by0 || sy1 != by1 || sx0 < bx0 || sx1 > bx1;
                } else fails += rot != 0;
            }
            fails += nb != 1;
        }
    }
    remove(left); remove(right); remove(turns); remove(ribbon); remove(plan);
    ves_selftest_clean(dir);
    Arena_dispose(&arena);
    fprintf(stderr, "[selftest] sheet_strip --repack %s (%d failures)\n", fails ? "FAIL" : "PASS", fails);
    return fails ? 1 : 0;
}

/* --plan / --repack with the blocks in a file: the same plan_run on the parsed list. */
static int plan_run_list(const char *layout, const char *stem, int32_t H, int32_t W, const char *list, const char *turns)
{
    int n = 0, rc = 0;
    char **blocks = NULL;
    if (BlockList_read(list, &n, &blocks) != 0) return 2;
    fprintf(stderr, "[sheet_strip] %d blocks from %s\n", n, list);
    rc = plan_run(layout, stem, H, W, (size_t)n, blocks, turns);
    free(blocks);
    return rc;
}

/* --plan through a block list (block_list.h): the repack fixture's two masks named in a file, any order, give the
 * fixture's three pieces. 0 = pass. */
static int block_list_selftest(void)
{
    int fails = 0;
    char dir[1024], list[1200], left[1200], right[1200], ribbon[1200], stem[1200], plan[1200];
    FILE *f = NULL;
    if (ves_selftest_dir(dir, sizeof dir, "sheet_strip_blocklist") != 0) { fprintf(stderr, "[selftest] sheet_strip block list: no scratch dir\n"); return 1; }
    snprintf(list, sizeof list, "%s/blocks.txt", dir);
    {                                                      /* end to end: the repack fixture's masks through a list */
        const int32_t W = 90, H = 12, cut = 37;
        Arena_T arena = Arena_new();
        uint8_t *cov = ARENA_CALLOC(arena, (size_t)W * (size_t)H, 1);
        uint8_t *l = ARENA_ALLOC(arena, (size_t)cut * (size_t)H), *r = ARENA_ALLOC(arena, (size_t)(W - cut) * (size_t)H);
        int ok = 1;
        for (int32_t x = 0; x < W; x++)
            for (int32_t y = 0; y < H; y++)
                if ((x >= 3 && x < 30 && y >= 2 && y < 9 - (x % 3)) || (x >= 41 && x < 66 && y >= 1 + (x % 4) && y < 11) ||
                    (x >= 70 && x < 74 && y >= 4 && y < 6))
                    cov[(size_t)y * (size_t)W + (size_t)x] = 255;
        for (int32_t y = 0; y < H; y++) {
            memcpy(l + (size_t)y * (size_t)cut, cov + (size_t)y * (size_t)W, (size_t)cut);
            memcpy(r + (size_t)y * (size_t)(W - cut), cov + (size_t)y * (size_t)W + (size_t)cut, (size_t)(W - cut));
        }
        snprintf(left, sizeof left, "%s/left.tif", dir);
        snprintf(right, sizeof right, "%s/right.tif", dir);
        snprintf(ribbon, sizeof ribbon, "%s/ribbon_layout.json", dir);
        snprintf(stem, sizeof stem, "%s/plan", dir);
        snprintf(plan, sizeof plan, "%s/plan.json", dir);
        f = fopen(ribbon, "wb");
        ok = f != NULL && fputs("{\"kept\": true, \"gap\": 8, \"regions\": 3}\n", f) >= 0 && fclose(f) == 0;
        ok = ok && TiffIO_save(left, l, 1, H, cut) == 0 && TiffIO_save(right, r, 1, H, W - cut) == 0;
        f = fopen(list, "wb");
        ok = ok && f != NULL && fprintf(f, "%d %s\n0 %s\n", (int)cut, right, left) > 0 && fclose(f) == 0;   /* any block order */
        if (!ok || plan_run_list(dir, stem, H, W, list, NULL) != 0) { fprintf(stderr, "  FAIL: --plan through the list\n"); fails++; }
        else {
            const char *err = NULL;
            const JsonValue *j = Json_parse_file(arena, plan, &err);
            const JsonValue *pl = j ? Json_object_get(j, "piece_list") : NULL;
            if (pl == NULL || Json_array_len(pl) != 3) { fprintf(stderr, "  FAIL: the list plan has %zu pieces, not 3\n", pl ? Json_array_len(pl) : (size_t)0); fails++; }
        }
        remove(left); remove(right); remove(ribbon); remove(plan);
        Arena_dispose(&arena);
    }
    remove(list);
    if (fails == 0) ves_selftest_clean(dir);
    fprintf(stderr, "[selftest] sheet_strip block list %s (%d failures)\n", fails ? "FAIL" : "PASS", fails);
    return fails ? 1 : 0;
}

int main(int argc, char **argv)
{
    if (argc == 2 && strcmp(argv[1], "--selftest") == 0)
        return StripLayout_selftest() != 0 || texture_selftest() != 0 || plan_selftest() != 0 || repack_selftest() != 0 ||
               BlockList_selftest() != 0 || block_list_selftest() != 0 || CloseupWindow_selftest() != 0;
    if (argc == 3 && strcmp(argv[1], "--closeup-window") == 0) return CloseupWindow_run(argv[2]);
    if (argc == 7 && strcmp(argv[1], "--plan") == 0 && argv[6][0] == '@')
        return plan_run_list(argv[2], argv[3], (int32_t)atol(argv[4]), (int32_t)atol(argv[5]), argv[6] + 1, NULL);
    if (argc == 8 && strcmp(argv[1], "--repack") == 0 && argv[7][0] == '@')
        return plan_run_list(argv[2], argv[3], (int32_t)atol(argv[4]), (int32_t)atol(argv[5]), argv[7] + 1, argv[6]);
    if (argc >= 8 && strcmp(argv[1], "--plan") == 0 && (argc - 6) % 2 == 0)
        return plan_run(argv[2], argv[3], (int32_t)atol(argv[4]), (int32_t)atol(argv[5]), (argc - 6) / 2, argv + 6, NULL);
    if (argc >= 9 && strcmp(argv[1], "--repack") == 0 && (argc - 7) % 2 == 0)
        return plan_run(argv[2], argv[3], (int32_t)atol(argv[4]), (int32_t)atol(argv[5]), (argc - 7) / 2, argv + 7, argv[6]);
    if (argc == 5 && strcmp(argv[3], "--texture") == 0) return run(argv[1], argv[2], argv[4]);
    if (argc != 3) {
        fprintf(stderr, "usage: sheet_strip <layout_dir> <out_stem> [--texture <same_grid.tif>]\n"
                        "       sheet_strip --plan <layout_dir> <out_stem> <H> <W> <col0> <block_mask.tif> [...]\n"
                        "       sheet_strip --plan <layout_dir> <out_stem> <H> <W> @<blocks.txt>   (lines: <col0> <block_mask.tif>)\n"
                        "       sheet_strip --repack <layout_dir> <out_stem> <H> <W> <turns.json> <col0> <block_mask.tif> [...]\n"
                        "       sheet_strip --repack <layout_dir> <out_stem> <H> <W> <turns.json> @<blocks.txt>\n"
                        "       sheet_strip --closeup-window <placement.json>\n"
                        "       sheet_strip --selftest\n");
        return 2;
    }
    return run(argv[1], argv[2], NULL);
}
