/* ink_closeup -- a continuous, isometric close-up of one ink region of an assembled sheet.
 *
 *   ink_closeup <config.json> [--stop-after surface|render]
 *   ink_closeup --selftest
 *
 * The config (configs/ink/closeup_*.json) names the region and its sources; tuning lives in the CLOSEUP_*
 * constants (pipeline_constants.h), never on the command line:
 *   "id"        experiment id: every output file is <out_dir>/<id>_*
 *   "out_dir"   output directory (created)
 *   "layout"    a review layout dir: stage6_repaired.vmesh (+ _chart.i32) and sheet_index.json
 *   "raw"       the RAW source (uncompressed uint8 Zarr v2, 128^3 chunks, or cubes_RAW)
 *   "axis"      the z,y,x scroll-axis table
 *   "voxel_um"  scan pitch (reports and scale bars)
 *   "seed_box"  [x0, y0, x1, y1] on the layout's sheet grid (= its UV, 1 px per vox): the region
 *   "grow_px"   the close-up takes the seed pieces' faces within this margin of the box
 *   "complete_gaps" optional wider segmentation-gap fill, with observed/fitted provenance retained
 *   "adopt_audit" optional per-round CSV diagnostics; does not change admission or fitting
 *   "surface_source" optional existing close-up stem (before _surface.vmesh): refit its observed
 *               faces and complete gaps in its existing UV frame; layout and seed are then unused
 *   "grow"      optional, with surface_source, "growth_uv" (the target box in its frame) and "pool": frontier growth
 *               (ic_stage_grow) -- the pools read again and again near the data boundary, chart by chart, until the
 *               sheet, the box or the adoption gates stop it; {margin_vox, rounds, max_reads, reflatten}
 *   "z_top_to_bottom"  optional, default true: the scan's catalog flag.  The frame puts the scroll's TOP at the top
 *               of the rows: increasing z down the rows when true (PHerc0139/1447/0211), decreasing z down the rows
 *               when false (PHerc0343 20250521140437 runs z from the scroll's bottom to its top -- with the default
 *               frame its letters read upside down).  A 180-degree turn: the reading handedness is unchanged.
 *
 * Stages (each resumable from its files):
 *   surface  the pieces owning the seed box; their faces within the grown box become observations in a
 *            FRAME (axis down the rows, team reading handedness); the fused surface (closeup_surface: robust
 *            local quadrics, holes and borders filled up to CLOSEUP_FILL_PX); its isometric flattening
 *            (coarse SLIM with a scaffold), rotated so the axial coordinate runs down the rows ->
 *              <id>_surface.vmesh (+ _nrm.f32 oriented unit normals, _attr.u8 1 data / 2 fill,
 *              _chart.i32 source chart per vertex), <id>_surface_report.json, <id>_frame_state.png
 *   render   the reader-v2 surface volume of the close-up at 1 px per vox (rawtex_layers: 21 layers inward,
 *            the context apron, the geometry volume), its CT bake and a provenance map ->
 *              <id>_layers.zarr, <id>_layers_mask.tif, <id>_geometry.zarr, <id>_ct.tif/.png,
 *              <id>_provenance.png, <id>_render_report.json
 * Threads: OMP_NUM_THREADS. */
#include <errno.h>
#include "../assemble/asm_axis.h"
#include "../assemble/asm_contacts.h"
#include "../common/arena.h"
#include "../common/json_read.h"
#include "../common/mesh_bin.h"
#include "../common/mesh_chart_index.h"
#include "../common/pipeline_constants.h"
#include "../common/raw_sample.h"
#include "../common/tiff_io.h"
#include "../common/uv_guard.h"
#include "../common/ves_platform.h"
#include "../common/ves_png.h"
#include "../common/mesh_normals.h"
#include "../common/kdtree.h"
#include "../flatten/chart_pool.h"
#include "../flatten/closeup_adopt.h"
#include "../flatten/closeup_cttrack.h"
#include "../flatten/closeup_f0bridge.h"
#include "../flatten/closeup_fiber.h"
#include "../flatten/closeup_fiber_volume.h"
#include "../flatten/closeup_surface.h"
#include "../flatten/closeup_import.h"
#include "../flatten/fiber_map.h"
#include "../flatten/ink_letters.h"
#include "../flatten/rawtex_bake.h"
#include "../flatten/rawtex_layers.h"
#include "../flatten/rawtex_tifxyz.h"
#include "../flatten/closeup_seed.h"
#include "../flatten/sparse_solve.h"
#include "../flatten/slim_refine.h"
#include "../flatten/seam_cut.h"

#include <math.h>
#include <float.h>
#include <stdarg.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/ves_omp.h"

typedef struct IcSource { char path[2048]; int surface; } IcSource;

typedef struct IcConfig {
    char id[256], out_dir[2048], layout[2048], raw[2048], axis[2048];
    char surface_source[2048];
    char fiber_volume[2048];       /* optional provenance-pinned physical 3-D field */
    IcSource *pool;
    int npool;
    int candidate_manifest;
    int z_top_to_bottom;           /* the scan's catalog flag: the frame puts the scroll's top at the top of the rows */
    int complete_gaps;             /* optional wider fill, still distinct from observed geometry */
    int seed_all_pieces;           /* whole-region requests retain small pieces beside large ones */
    int adopt_audit;               /* read-only per-chart admission diagnostics */
    int fiber_constraints;         /* measured physical CT arms in recovery SLIM */
    CloseupAdoptOpts adopt;        /* admission gates (CloseupAdoptOpts_default + config "adopt" overrides) */
    int adopt_rounds;              /* adoption rounds per pool (CLOSEUP_ADOPT_ROUNDS unless "adopt" overrides) */
    double voxel_um, seed[4], grow_px;
    double growth_uv[4], search_margin;
    int has_growth;
    double ct_track_box[6], ct_track_seed_z;
    int ct_track_configured;
    /* gated bridging (2026-10-13): "wrap_charts" = wrap_index.exe's wrap_charts.csv (chart -> wrap) and "bridge" =
     * true let observed blobs on one wrap within 2 x CLOSEUP_BRIDGE_PX of each other share one surface component */
    char wrap_charts[2048];
    int bridge;
    int32_t *chart_wrap;
    int32_t nchart_wrap;
    /* frontier growth (2026-10-14, the PHerc0211 candidates grown 3x; ic_stage_grow): "grow" {margin_vox, rounds,
     * max_reads, reflatten} with surface_source + growth_uv + pool */
    int has_grow, grow_rounds, grow_max_reads, grow_reflatten;
    double grow_margin_vox;
    double grow_min_gain;          /* two reads in a row growing the data by less than this share end it (0: off) */
    double grow_tile_px;           /* > 0: the pools' charts judged as tiles of this many px of their own UV */
    /* the F0 bridge of a growth (closeup_f0bridge.c): "grow" {"f0": [cubes_PRED dirs], "f0_depth_px", "f0_rings",
     * "f0_passes"}; no dirs = off */
    int nf0, f0_rings, f0_passes;
    char f0_dir[16][1024];
    double f0_depth_px;
} IcConfig;

/* wrap_charts.csv (chart,...,wrap,...) -> c->chart_wrap[chart] (INT32_MIN where unknown); 0 on success */
static int ic_load_wrap_charts(IcConfig *c)
{
    FILE *fp = fopen(c->wrap_charts, "rb");
    char line[4096], *f[64];
    int nf = 0, cc = -1, cw = -1, k = 0;
    long maxc = -1;
    if (fp == NULL) return -1;
    if (fgets(line, sizeof line, fp) == NULL) { fclose(fp); return -1; }
    {   char *p = line;
        for (nf = 0; nf < 64 && p != NULL; nf++) { char *q = strchr(p, ','); f[nf] = p; if (q) *q = 0; p = q ? q + 1 : NULL; }
        {   char *e = f[nf - 1] + strlen(f[nf - 1]);
            while (e > f[nf - 1] && (e[-1] == '\n' || e[-1] == '\r')) *--e = 0; }
        for (k = 0; k < nf; k++) { if (!strcmp(f[k], "chart")) cc = k; if (!strcmp(f[k], "wrap")) cw = k; }
    }
    if (cc < 0 || cw < 0) { fclose(fp); return -1; }
    while (fgets(line, sizeof line, fp) != NULL) {
        char *p = line; int n = 0; long ch = 0;
        for (n = 0; n < 64 && p != NULL; n++) { char *q = strchr(p, ','); f[n] = p; if (q) *q = 0; p = q ? q + 1 : NULL; }
        if (cc >= n || cw >= n) continue;
        ch = strtol(f[cc], NULL, 10);
        if (ch > maxc) maxc = ch;
    }
    if (maxc < 0 || maxc >= 100000000L) { fclose(fp); return -1; }
    c->nchart_wrap = (int32_t)(maxc + 1);
    c->chart_wrap = (int32_t *)malloc((size_t)c->nchart_wrap * sizeof *c->chart_wrap);
    if (c->chart_wrap == NULL) { fclose(fp); return -1; }
    for (k = 0; k < c->nchart_wrap; k++) c->chart_wrap[k] = INT32_MIN;
    rewind(fp);
    if (fgets(line, sizeof line, fp) == NULL) { fclose(fp); return -1; }
    while (fgets(line, sizeof line, fp) != NULL) {
        char *p = line; int n = 0; long ch = 0;
        for (n = 0; n < 64 && p != NULL; n++) { char *q = strchr(p, ','); f[n] = p; if (q) *q = 0; p = q ? q + 1 : NULL; }
        if (cc >= n || cw >= n || !f[cw][0]) continue;
        ch = strtol(f[cc], NULL, 10);
        if (ch < 0 || ch >= c->nchart_wrap) continue;
        c->chart_wrap[ch] = (int32_t)strtol(f[cw], NULL, 10);
    }
    fclose(fp);
    return 0;
}

/* A JSON number with one decimal into buf, or "null" when v is not finite (JSON has no nan / inf). */
static const char *ic_json_num(char *buf, size_t cap, double v)
{
    if (!isfinite(v)) return "null";
    snprintf(buf, cap, "%.1f", v);
    return buf;
}

/* JSON string literal (Windows paths carry backslashes) into buf. */
static void ic_json_str(char *buf, size_t cap, const char *s)
{
    size_t o = 0;
    if (cap < 3) { if (cap) buf[0] = 0; return; }
    buf[o++] = '"';
    for (; *s && o + 8 < cap; s++) {
        if (*s == '"' || *s == '\\') { buf[o++] = '\\'; buf[o++] = *s; }
        else if ((unsigned char)*s < 0x20) o += (size_t)snprintf(buf + o, cap - o, "\\u%04x", (unsigned)(unsigned char)*s);
        else buf[o++] = *s;
    }
    buf[o++] = '"';
    buf[o] = 0;
}

static int ic_copy_str(char *dst, size_t cap, const JsonValue *obj, const char *key)
{
    const char *s = Json_as_string(Json_object_get(obj, key));
    if (s == NULL || strlen(s) + 1 > cap) { fprintf(stderr, "ERROR: config needs a string \"%s\"\n", key); return -1; }
    memcpy(dst, s, strlen(s) + 1);
    return 0;
}

static int ic_config(const char *path, IcConfig *c)
{
    Arena_T a = Arena_new();
    const char *err = NULL;
    const JsonValue *j = Json_parse_file(a, path, &err), *sb = NULL;
    int rc = -1, k = 0;
    memset(c, 0, sizeof *c);
    if (j == NULL) { fprintf(stderr, "ERROR: cannot parse %s (%s)\n", path, err ? err : "?"); goto out; }
    if (ic_copy_str(c->id, sizeof c->id, j, "id") || ic_copy_str(c->out_dir, sizeof c->out_dir, j, "out_dir") ||
        ic_copy_str(c->layout, sizeof c->layout, j, "layout") || ic_copy_str(c->raw, sizeof c->raw, j, "raw") ||
        ic_copy_str(c->axis, sizeof c->axis, j, "axis"))
        goto out;
    c->voxel_um = Json_member_double(j, "voxel_um", 0.0);
    c->grow_px = Json_member_double(j, "grow_px", 0.0);
    c->search_margin = Json_member_double(j, "search_margin_vox", CLOSEUP_ADOPT_MARGIN_PX);
    if (!(c->search_margin > 0) || !isfinite(c->search_margin)) goto out;
    sb = Json_object_get(j, "growth_uv");
    if (sb != NULL) {
        if (Json_array_len(sb) != 4) goto out;
        for (k = 0; k < 4; k++) {
            c->growth_uv[k] = Json_as_double(Json_array_get(sb, (size_t)k), NAN);
            if (!isfinite(c->growth_uv[k])) goto out;
        }
        if (c->growth_uv[2] <= c->growth_uv[0] || c->growth_uv[3] <= c->growth_uv[1]) goto out;
        c->has_growth = 1;
    }
    c->z_top_to_bottom = Json_as_bool(Json_object_get(j, "z_top_to_bottom"), 1);
    c->complete_gaps = Json_as_bool(Json_object_get(j, "complete_gaps"), 0);
    c->seed_all_pieces = Json_as_bool(Json_object_get(j, "seed_all_pieces"), 0);
    c->adopt_audit = Json_as_bool(Json_object_get(j, "adopt_audit"), 0);
    c->fiber_constraints = Json_as_bool(Json_object_get(j, "fiber_constraints"), 0);
    CloseupAdoptOpts_default(&c->adopt);
    c->adopt_rounds = CLOSEUP_ADOPT_ROUNDS;
    {   /* "adopt": optional overrides of the admission gates for one run (e.g. a relaxed completion pass over a
         * region's holes); keys = CloseupAdoptOpts fields + "rounds"; absent keys keep the defaults */
        const JsonValue *ad = Json_object_get(j, "adopt");
        if (ad != NULL) {
            c->adopt.rms_px = Json_member_double(ad, "rms_px", c->adopt.rms_px);
            c->adopt.normal_deg = Json_member_double(ad, "normal_deg", c->adopt.normal_deg);
            c->adopt.conflict_vox = Json_member_double(ad, "conflict_vox", c->adopt.conflict_vox);
            c->adopt.agree_frac = Json_member_double(ad, "agree_frac", c->adopt.agree_frac);
            c->adopt.on_vox = Json_member_double(ad, "on_vox", c->adopt.on_vox);
            c->adopt.min_extent_px = Json_member_double(ad, "min_extent_px", c->adopt.min_extent_px);
            c->adopt.new_frac = Json_member_double(ad, "new_frac", c->adopt.new_frac);
            c->adopt.min_corr = (int)Json_member_double(ad, "min_corr", (double)c->adopt.min_corr);
            c->adopt.metric_ratio = Json_member_double(ad, "metric_ratio", c->adopt.metric_ratio);
            c->adopt.metric_frac = Json_member_double(ad, "metric_frac", c->adopt.metric_frac);
            c->adopt.continuity_px = Json_member_double(ad, "continuity_px", c->adopt.continuity_px);
            c->adopt.continuity_slope = Json_member_double(ad, "continuity_slope", c->adopt.continuity_slope);
            c->adopt.layer_vox = Json_member_double(ad, "layer_vox", c->adopt.layer_vox);
            c->adopt.layer_frac = Json_member_double(ad, "layer_frac", c->adopt.layer_frac);
            c->adopt.foot_fill_px = Json_member_double(ad, "foot_fill_px", c->adopt.foot_fill_px);
            c->adopt_rounds = (int)Json_member_double(ad, "rounds", (double)c->adopt_rounds);
            if (!(c->adopt.rms_px > 0.0) || !(c->adopt.normal_deg > 0.0) || !(c->adopt.conflict_vox > 0.0) ||
                !(c->adopt.metric_ratio > 1.0) || !(c->adopt.metric_frac >= 0.0 && c->adopt.metric_frac <= 1.0) ||
                !(c->adopt.continuity_px >= 0.0) || !(c->adopt.continuity_slope > 0.0) ||
                !(c->adopt.layer_vox >= 0.0 && c->adopt.layer_vox <= 64.0) ||
                !(c->adopt.layer_frac >= 0.0 && c->adopt.layer_frac <= 1.0) ||
                !(c->adopt.foot_fill_px >= 0.0 && c->adopt.foot_fill_px <= 256.0) ||
                c->adopt_rounds < 1 || c->adopt_rounds > 64) {
                fprintf(stderr, "ERROR: config \"adopt\" out of range\n");
                goto out;
            }
            fprintf(stderr, "[ink_closeup] adopt overrides: rms %.2f px, normal %.1f deg, conflict %.2f vox, "
                    "agree %.2f, foot fill %.1f px, rounds %d\n", c->adopt.rms_px, c->adopt.normal_deg,
                    c->adopt.conflict_vox, c->adopt.agree_frac, c->adopt.foot_fill_px, c->adopt_rounds);
        }
    }
    {
        const JsonValue *box = Json_object_get(j, "ct_track_box");
        if (box != NULL) {
            if (Json_array_len(box) != 6 || Json_object_get(j, "ct_track_seed_z") == NULL) goto out;
            for (k = 0; k < 6; k++) c->ct_track_box[k] = Json_as_double(Json_array_get(box, (size_t)k), NAN);
            c->ct_track_seed_z = Json_member_double(j, "ct_track_seed_z", NAN);
            c->ct_track_configured = 1;
        }
    }
    if (Json_object_get(j, "surface_source") != NULL &&
        ic_copy_str(c->surface_source, sizeof c->surface_source, j, "surface_source") != 0) goto out;
    if (Json_object_get(j, "fiber_volume") != NULL &&
        ic_copy_str(c->fiber_volume, sizeof c->fiber_volume, j, "fiber_volume") != 0) goto out;
    {   /* "pool": layouts whose charts may join the close-up (the reduced placement); optional */
        const JsonValue *pl = Json_object_get(j, "pool");
        const JsonValue *sources = NULL;
        const char *manifest = Json_as_string(Json_object_get(j, "candidate_sources"));
        size_t q = 0, n;
        if (manifest) {
            c->candidate_manifest = 1;
            const JsonValue *mj = Json_parse_file(a, manifest, &err);
            if (!mj) { fprintf(stderr, "ERROR: candidate sources %s (%s)\n", manifest, err ? err : "?"); goto out; }
            sources = Json_object_get(mj, "sources");
        }
        n = Json_array_len(pl) + Json_array_len(sources);
        if (n > INT_MAX || n > SIZE_MAX / sizeof *c->pool) goto out;
        c->pool = (IcSource *)calloc(n ? n : 1, sizeof *c->pool);
        if (!c->pool) goto out;
        for (q = 0; q < n; q++) {
            const JsonValue *item = q < Json_array_len(pl) ? Json_array_get(pl, q) : Json_array_get(sources, q - Json_array_len(pl));
            const char *ps = Json_as_string(item);
            int surface = 0;
            if (!ps) {
                ps = Json_as_string(Json_object_get(item, "surface")); surface = ps != NULL;
                if (!ps) ps = Json_as_string(Json_object_get(item, "layout"));
            }
            if (!ps || strlen(ps) + 1 > sizeof c->pool[0].path) {
                fprintf(stderr, "ERROR: invalid candidate source %zu\n", q); goto out;
            }
            memcpy(c->pool[c->npool].path, ps, strlen(ps) + 1);
            c->pool[c->npool++].surface = surface;
        }
    }
    sb = Json_object_get(j, "seed_box");
    if (Json_array_len(sb) != 4 || !(c->voxel_um > 0.0)) {
        fprintf(stderr, "ERROR: config needs \"seed_box\": [x0, y0, x1, y1] and \"voxel_um\"\n");
        goto out;
    }
    for (k = 0; k < 4; k++) c->seed[k] = Json_as_double(Json_array_get(sb, (size_t)k), 0.0);
    if (!(c->seed[2] > c->seed[0]) || !(c->seed[3] > c->seed[1])) { fprintf(stderr, "ERROR: empty seed_box\n"); goto out; }
    {   /* "grow": frontier growth over growth_uv from surface_source (ic_stage_grow); absent keys keep the defaults */
        const JsonValue *gr = Json_object_get(j, "grow");
        if (gr != NULL) {
            c->has_grow = 1;
            c->grow_margin_vox = Json_member_double(gr, "margin_vox", CLOSEUP_GROW_MARGIN_VOX);
            c->grow_rounds = (int)Json_member_double(gr, "rounds", (double)CLOSEUP_GROW_ROUNDS);
            c->grow_max_reads = (int)Json_member_double(gr, "max_reads", (double)CLOSEUP_GROW_MAX_READS);
            c->grow_reflatten = Json_as_bool(Json_object_get(gr, "reflatten"), 1);
            c->grow_min_gain = Json_member_double(gr, "min_gain", CLOSEUP_GROW_MIN_GAIN);
            c->grow_tile_px = Json_member_double(gr, "tile_px", 0.0);
            c->f0_depth_px = Json_member_double(gr, "f0_depth_px", CLOSEUP_F0_DEPTH_PX);
            c->f0_rings = (int)Json_member_double(gr, "f0_rings", (double)CLOSEUP_F0_RINGS);
            c->f0_passes = (int)Json_member_double(gr, "f0_passes", (double)CLOSEUP_F0_PASSES);
            {
                const JsonValue *fd = Json_object_get(gr, "f0");
                size_t q = 0;
                c->nf0 = 0;
                if (Json_array_len(fd) > 16) { fprintf(stderr, "ERROR: \"grow\".\"f0\" holds more than 16 dirs\n"); goto out; }
                for (q = 0; q < Json_array_len(fd); q++) {
                    const char *d = Json_as_string(Json_array_get(fd, q));
                    if (d == NULL || strlen(d) + 1 > sizeof c->f0_dir[0]) { fprintf(stderr, "ERROR: bad \"grow\".\"f0\" dir\n"); goto out; }
                    memcpy(c->f0_dir[c->nf0++], d, strlen(d) + 1);
                }
            }
            if (!(c->grow_margin_vox > 0.0) || !isfinite(c->grow_margin_vox) || c->grow_rounds < 1 || c->grow_rounds > 64 ||
                c->grow_max_reads < 1 || c->grow_max_reads > 100000 || !(c->grow_min_gain >= 0.0 && c->grow_min_gain < 1.0) ||
                !(c->grow_tile_px == 0.0 || (c->grow_tile_px >= 16.0 && c->grow_tile_px <= 4096.0)) ||
                !(c->f0_depth_px >= 0.0 && c->f0_depth_px <= 100000.0) || c->f0_rings < 1 || c->f0_rings > 64 ||
                c->f0_passes < 1 || c->f0_passes > 1000) {
                fprintf(stderr, "ERROR: config \"grow\" out of range\n");
                goto out;
            }
        }
    }
    c->bridge = Json_as_bool(Json_object_get(j, "bridge"), 0);
    if (Json_object_get(j, "wrap_charts") != NULL) {
        if (ic_copy_str(c->wrap_charts, sizeof c->wrap_charts, j, "wrap_charts")) goto out;
        if (ic_load_wrap_charts(c) != 0) { fprintf(stderr, "ERROR: cannot read wrap_charts %s\n", c->wrap_charts); goto out; }
    }
    rc = 0;
out:
    if (rc) { free(c->pool); c->pool = NULL; c->npool = 0; free(c->chart_wrap); c->chart_wrap = NULL; c->nchart_wrap = 0; }
    Arena_dispose(&a);
    return rc;
}

static int ic_seed_selected(size_t count, size_t total, int all)
{
    return all ? count > 0 : count >= 1000 && (double)count >= 0.02 * (double)total;
}

/* A whole binary file (malloc'd); NULL when absent, unreadable, or not `expect` bytes (0 = any size). */
static void *ic_read_file(const char *path, size_t expect, size_t *got)
{
    FILE *f = fopen(path, "rb");
    void *buf = NULL;
    long long n = 0;
    if (f == NULL) return NULL;
#ifdef _WIN32
    if (_fseeki64(f, 0, SEEK_END) != 0 || (n = _ftelli64(f)) < 0 || _fseeki64(f, 0, SEEK_SET) != 0) { fclose(f); return NULL; }
#else
    if (fseeko(f, 0, SEEK_END) != 0 || (n = (long long)ftello(f)) < 0 || fseeko(f, 0, SEEK_SET) != 0) { fclose(f); return NULL; }
#endif
    if (expect && (size_t)n != expect) { fclose(f); return NULL; }
    *got = (size_t)n;
    buf = malloc(*got ? *got : 1);
    if (buf == NULL || fread(buf, 1, *got, f) != *got) { free(buf); fclose(f); return NULL; }
    fclose(f);
    return buf;
}

static int ic_write_file(const char *path, const void *buf, size_t n)
{
    FILE *f = fopen(path, "wb");
    if (f == NULL) return -1;
    if (n && fwrite(buf, 1, n, f) != n) { fclose(f); return -1; }
    return fclose(f) == 0 ? 0 : -1;
}

/* The frames run their axial field down the rows: negate it on a scan whose z runs from the scroll's bottom to its
 * top, so the scroll's top lands at the top of the close-up (a 180-degree turn of the frame; handedness unchanged). */
static void ic_axial_sense(const IcConfig *c, double *axial, size_t n)
{
    size_t i = 0;
    if (c->z_top_to_bottom) return;
    for (i = 0; i < n; i++) axial[i] = -axial[i];
}

#define IC_AXIS_INDEX_MIN 4096   /* points at which the nearest-foot index pays for its build */

/* Unit direction to the axis foot and the axial arc length of every point.  Without the arc length only the foot is
 * needed, and the axis's nearest-foot index returns AsmAxis_project's own foot (the asm_axis selftest compares them
 * bit for bit) without scanning the ~500 segments of each point's z window: 44 s of a completion pass's 120 M pool
 * vertices on 16 threads (2026-10-30). */
/* One point of ic_axis_fields: the foot by the index when there is one and no arc length is wanted. */
static void ic_axis_point(const AsmAxis *axis, const AsmAxisProjector *index, const float *p, double *axial,
                          float *inward)
{
    double q[3], s = 0.0, foot[3], d[3], l = 0.0;
    int k = 0;
    for (k = 0; k < 3; k++) q[k] = (double)p[k];
    if (!isfinite(q[0])) { if (axial) *axial = 0.0; if (inward) memset(inward, 0, 3 * sizeof *inward); return; }
    if (index != NULL && axial == NULL) AsmAxis_project_foot(index, q, foot);
    else AsmAxis_project(axis, q, &s, NULL, NULL, foot);
    for (k = 0; k < 3; k++) { d[k] = foot[k] - q[k]; l += d[k] * d[k]; }
    l = sqrt(l);
    if (axial) *axial = s;
    if (inward)
        for (k = 0; k < 3; k++) inward[k] = l > 1e-9 ? (float)(d[k] / l) : 0.0f;
}

static void ic_axis_fields(const AsmAxis *axis, const float *p, size_t n, double *axial, float *inward)
{
    Arena_T arena = NULL;
    const AsmAxisProjector *index = NULL;
    long i = 0;
    if (axial == NULL && n >= IC_AXIS_INDEX_MIN) {
        arena = Arena_new_sized((size_t)1 << 20);
        index = AsmAxis_projector(arena, axis);
    }
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (i = 0; i < (long)n; i++)
        ic_axis_point(axis, index, &p[(size_t)i * 3], axial ? &axial[i] : NULL, inward ? &inward[(size_t)i * 3] : NULL);
    if (arena != NULL) Arena_dispose(&arena);
}

/* ic_axis_fields' inward for the points whose position changed since prev[] (bitwise; prev is updated): the same
 * function of the same bits, so the field is the one a full pass computes.  A completion's refits move a few percent
 * of the grid's nodes. */
static void ic_axis_inward_moved(const AsmAxis *axis, const AsmAxisProjector *index, const float *p, size_t n,
                                 float *prev, float *inward)
{
    long i = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096)
#endif
    for (i = 0; i < (long)n; i++) {
        if (memcmp(&p[(size_t)i * 3], &prev[(size_t)i * 3], 3 * sizeof *p) == 0) continue;
        ic_axis_point(axis, index, &p[(size_t)i * 3], NULL, &inward[(size_t)i * 3]);
        memcpy(&prev[(size_t)i * 3], &p[(size_t)i * 3], 3 * sizeof *p);
    }
}


static uint32_t ic_hash(uint32_t x)
{
    x ^= x >> 16; x *= 0x7feb352dU; x ^= x >> 15; x *= 0x846ca68bU; x ^= x >> 16;
    return x;
}

/* ------------------------------------------------------------------ pool */

/* A pool's arrays, by kind: what a chain keeps from one pass's pool for the next (ic_pool_keep) instead of handing
 * ~9 GB back to the system and faulting it in again, zeroed, a few seconds later. */
enum { IB_VOFS, IB_FOFS, IB_CHART_ID, IB_FACES, IB_FACE_CHART, IB_VERTS, IB_UV, IB_INWARD, IB_NRM, IB_READY, IB_COUNT };
typedef struct IcPoolBufs {
    void *ptr[IB_COUNT];
    size_t bytes[IB_COUNT];
} IcPoolBufs;

typedef struct IcPool {
    CloseupPool p;
    float *verts, *nrm, *inward;
    double *uv;
    int32_t *faces, *face_chart, *chart_id;
    size_t *vofs, *fofs;
    Arena_T arena;                  /* the lazy fields' axis index (ic_pool_load lazy) */
    const AsmAxisProjector *index;
    uint8_t *chart_ready;           /* [nch] */
    size_t bytes[IB_COUNT];         /* each array's allocation */
} IcPool;

static void ic_pool_free(IcPool *q)
{
    free(q->verts); free(q->nrm); free(q->inward); free(q->uv); free(q->faces); free(q->face_chart);
    free(q->chart_id); free(q->vofs); free(q->fofs); free(q->chart_ready);
    if (q->arena != NULL) Arena_dispose(&q->arena);
    memset(q, 0, sizeof *q);
}

static void ic_bufs_free(IcPoolBufs *b)
{
    int k = 0;
    for (k = 0; k < IB_COUNT; k++) free(b->ptr[k]);
    memset(b, 0, sizeof *b);
}

/* A buffer of at least `bytes` for array kind k: the spare one when it is large enough, else a new one.  The
 * allocation's size lands in *got.  Contents are undefined. */
static void *ic_buf(IcPoolBufs *spare, int k, size_t bytes, size_t *got)
{
    void *p = NULL;
    if (bytes == 0) bytes = 1;
    if (spare != NULL && spare->ptr[k] != NULL && spare->bytes[k] >= bytes) {
        p = spare->ptr[k]; *got = spare->bytes[k];
        spare->ptr[k] = NULL; spare->bytes[k] = 0;
        return p;
    }
    if (spare != NULL) { free(spare->ptr[k]); spare->ptr[k] = NULL; spare->bytes[k] = 0; }
    p = malloc(bytes);
    *got = p != NULL ? bytes : 0;
    return p;
}

/* The pool's arrays into spare for the next pass's pool; the rest freed. */
static void ic_pool_keep(IcPool *q, IcPoolBufs *spare)
{
    void *arrays[IB_COUNT];
    int k = 0;
    arrays[IB_VOFS] = q->vofs; arrays[IB_FOFS] = q->fofs; arrays[IB_CHART_ID] = q->chart_id; arrays[IB_FACES] = q->faces;
    arrays[IB_FACE_CHART] = q->face_chart; arrays[IB_VERTS] = q->verts; arrays[IB_UV] = q->uv;
    arrays[IB_INWARD] = q->inward; arrays[IB_NRM] = q->nrm; arrays[IB_READY] = q->chart_ready;
    for (k = 0; k < IB_COUNT; k++) {
        free(spare->ptr[k]);
        spare->ptr[k] = arrays[k]; spare->bytes[k] = arrays[k] != NULL ? q->bytes[k] : 0;
    }
    if (q->arena != NULL) Arena_dispose(&q->arena);
    memset(q, 0, sizeof *q);
}

/* CloseupPool.chart_fill of a pool whose normals and inward field are computed per chart on first use (a completion
 * round evaluates ~3,000 of its ~640,000 charts): the chart's run of MeshNormals_compute_runs and its points of
 * ic_axis_fields, the values a full pass gives them.  ctx is the IcPool. */
static void ic_pool_chart_fill(const void *ctx, size_t chart)
{
    const IcPool *q = (const IcPool *)ctx;
    const size_t v0 = q->vofs[chart], v1 = q->vofs[chart + 1];
    size_t v = 0;
    MeshNormals_run(q->verts, q->faces, q->fofs[chart], q->fofs[chart + 1], v0, v1, q->nrm);
    for (v = v0; v < v1; v++) ic_axis_point(NULL, q->index, &q->verts[v * 3], NULL, &q->inward[v * 3]);
}

/* The pool's near region (CLOSEUP_POOL_NEAR_SURFACE): cells holding a supported node of the fused surface, dilated
 * (box metric) by the node grid's reach x CLOSEUP_POOL_NEAR_SLACK. */
typedef struct IcNear {
    double lo[3], cell;
    size_t n[3];
    uint8_t *on;                    /* [n0 * n1 * n2]: within the margin of a supported node */
} IcNear;

static void ic_near_free(IcNear *nr)
{
    free(nr->on);
    memset(nr, 0, sizeof *nr);
}

/* 1-D dilation of a line of `len` cells `stride` apart by `r` cells (two distance sweeps) */
static void ic_near_dilate_line(uint8_t *on, size_t len, size_t stride, size_t r, size_t *dist)
{
    size_t i = 0, d = r + 1;
    for (i = 0; i < len; i++) {                      /* distance to the nearest marked cell on the left */
        d = on[i * stride] == 1 ? 0 : (d <= r ? d + 1 : r + 1);
        dist[i] = d;
    }
    for (i = len, d = r + 1; i-- > 0;) {             /* ... and on the right */
        d = on[i * stride] == 1 ? 0 : (d <= r ? d + 1 : r + 1);
        if (d < dist[i]) dist[i] = d;
    }
    for (i = 0; i < len; i++) if (dist[i] <= r && on[i * stride] == 0) on[i * stride] = 2;
}

/* The region of the supported nodes (state != NONE) dilated by `reach` vox, over the search box [lo, hi]: 0 built, 1
 * nothing to build (no supported node: the search box alone decides, as it always did), -1 no memory.  Cells are
 * CLOSEUP_POOL_CELL_VOX, coarser when the box would need more than 2^27 of them. */
static int ic_near_build(const float *pos, const uint8_t *state, size_t nn, const double lo[3], const double hi[3],
                         double reach, IcNear *nr)
{
    size_t k = 0, total = 1, v = 0, r = 0, a = 0, b = 0, maxlen = 0, *dist = NULL;
    double vol = 1.0;
    memset(nr, 0, sizeof *nr);
    for (k = 0; k < 3; k++) {
        double ext = hi[k] - lo[k];
        if (!(ext >= 0.0) || !isfinite(ext)) return 1;
        vol *= ext + CLOSEUP_POOL_CELL_VOX;
    }
    if (!(reach >= 0.0) || !isfinite(reach)) return 1;
    nr->cell = fmax(CLOSEUP_POOL_CELL_VOX, cbrt(vol / (double)((size_t)1 << 27)));
    for (k = 0; k < 3; k++) {
        nr->lo[k] = lo[k];
        nr->n[k] = (size_t)floor((hi[k] - lo[k]) / nr->cell) + 1;
        if (nr->n[k] > maxlen) maxlen = nr->n[k];
        total *= nr->n[k];
    }
    nr->on = (uint8_t *)calloc(total, 1);
    dist = (size_t *)malloc(maxlen * sizeof *dist);
    if (nr->on == NULL || dist == NULL) { free(dist); ic_near_free(nr); return -1; }
    for (v = 0; v < nn; v++) {
        size_t c[3];
        int inside = 1;
        if (state[v] == CLOSEUP_NODE_NONE) continue;
        for (k = 0; k < 3; k++) {
            double t = ((double)pos[v * 3 + k] - nr->lo[k]) / nr->cell;
            if (!(t >= 0.0) || t >= (double)nr->n[k]) { inside = 0; break; }
            c[k] = (size_t)t;
        }
        if (inside) nr->on[(c[0] * nr->n[1] + c[1]) * nr->n[2] + c[2]] = 1;
    }
    r = (size_t)ceil(reach / nr->cell);
    /* separable box dilation: along x, then y, then z (a marked cell is 1, a dilated one 2; each pass reads both) */
    for (a = 0; a < nr->n[0]; a++)
        for (b = 0; b < nr->n[1]; b++) ic_near_dilate_line(nr->on + (a * nr->n[1] + b) * nr->n[2], nr->n[2], 1, r, dist);
    for (k = 0; k < total; k++) if (nr->on[k]) nr->on[k] = 1;
    for (a = 0; a < nr->n[0]; a++)
        for (b = 0; b < nr->n[2]; b++) ic_near_dilate_line(nr->on + a * nr->n[1] * nr->n[2] + b, nr->n[1], nr->n[2], r, dist);
    for (k = 0; k < total; k++) if (nr->on[k]) nr->on[k] = 1;
    for (a = 0; a < nr->n[1]; a++)
        for (b = 0; b < nr->n[2]; b++) ic_near_dilate_line(nr->on + a * nr->n[2] + b, nr->n[0], nr->n[1] * nr->n[2], r, dist);
    for (k = 0; k < total; k++) if (nr->on[k]) nr->on[k] = 1;
    free(dist);
    return 0;
}

/* The farthest node of the grid from a supported node, in frame px: adopted papyrus lands on the grid's nodes, and an
 * isometric chart lies no farther from the sheet in 3-D than in the frame.  A two-pass chamfer distance (steps h and
 * h sqrt 2) is a path length in the node graph, never below the Euclidean distance.  0 without a supported node, -1 on
 * no memory. */
static double ic_grid_reach(const CloseupGrid *g)
{
    size_t ni = g->ni, nj = g->nj, i = 0, j = 0, k = 0;
    const double h = g->h, hd = g->h * sqrt(2.0);
    double r = 0.0, *d = (double *)malloc((ni * nj ? ni * nj : 1) * sizeof *d);
    if (d == NULL) return -1.0;
    for (k = 0; k < ni * nj; k++) d[k] = g->state[k] != CLOSEUP_NODE_NONE ? 0.0 : INFINITY;
    for (j = 0; j < nj; j++)
        for (i = 0; i < ni; i++) {
            double *p = &d[j * ni + i];
            if (i > 0) *p = fmin(*p, d[j * ni + i - 1] + h);
            if (j > 0) {
                *p = fmin(*p, d[(j - 1) * ni + i] + h);
                if (i > 0) *p = fmin(*p, d[(j - 1) * ni + i - 1] + hd);
                if (i + 1 < ni) *p = fmin(*p, d[(j - 1) * ni + i + 1] + hd);
            }
        }
    for (j = nj; j-- > 0;)
        for (i = ni; i-- > 0;) {
            double *p = &d[j * ni + i];
            if (i + 1 < ni) *p = fmin(*p, d[j * ni + i + 1] + h);
            if (j + 1 < nj) {
                *p = fmin(*p, d[(j + 1) * ni + i] + h);
                if (i + 1 < ni) *p = fmin(*p, d[(j + 1) * ni + i + 1] + hd);
                if (i > 0) *p = fmin(*p, d[(j + 1) * ni + i - 1] + hd);
            }
        }
    for (k = 0; k < ni * nj; k++) if (isfinite(d[k]) && d[k] > r) r = d[k];
    free(d);
    return r;
}

/* Build the near region for a pool read (CLOSEUP_POOL_NEAR_SURFACE): *use is NULL when there is none to apply. */
static int ic_near_for_pool(const CloseupGrid *g, const double lo[3], const double hi[3], IcNear *nr, const IcNear **use)
{
    double reach = 0.0;
    int rc = 0;
    memset(nr, 0, sizeof *nr);
    *use = NULL;
    if (!CLOSEUP_POOL_NEAR_SURFACE) return 0;
    reach = ic_grid_reach(g);
    if (reach < 0.0) return -1;
    reach *= CLOSEUP_POOL_NEAR_SLACK;
    rc = ic_near_build(g->pos, g->state, g->ni * g->nj, lo, hi, reach, nr);
    if (rc < 0) return -1;
    if (rc == 0) {
        *use = nr;
        fprintf(stderr, "  pool near region: within %.0f vox of the fused surface (%.0f px grid reach x %.2f), cells of "
                "%.0f vox\n", reach, reach / CLOSEUP_POOL_NEAR_SLACK, CLOSEUP_POOL_NEAR_SLACK, nr->cell);
    }
    return 0;
}

/* a box [lo, hi] meets the region: some cell it overlaps is in it */
static int ic_near_meets(const IcNear *nr, const double lo[3], const double hi[3])
{
    size_t c0[3], c1[3], a = 0, b = 0, d = 0, k = 0;
    for (k = 0; k < 3; k++) {
        double t0 = floor((lo[k] - nr->lo[k]) / nr->cell), t1 = floor((hi[k] - nr->lo[k]) / nr->cell);
        if (!(t1 >= 0.0) || t0 >= (double)nr->n[k] || t1 < t0) return 0;
        c0[k] = t0 < 0.0 ? 0 : (size_t)t0;
        c1[k] = t1 >= (double)nr->n[k] ? nr->n[k] - 1 : (size_t)t1;
    }
    for (a = c0[0]; a <= c1[0]; a++)
        for (b = c0[1]; b <= c1[1]; b++)
            for (d = c0[2]; d <= c1[2]; d++)
                if (nr->on[(a * nr->n[1] + b) * nr->n[2] + d]) return 1;
    return 0;
}

/* A pool chart is a candidate when the 3-D box of all its vertices (blo, bhi) meets the search box [lo, hi] and, with
 * a near region, the region.  One test for the chart index's records and for a whole read's charts, so a paged and a
 * whole read select the same charts; a candidate's faces are then clipped to the search box as they always were. */
static int ic_pool_candidate(const double blo[3], const double bhi[3], const double lo[3], const double hi[3],
                             const IcNear *nreg)
{
    const double eps = 1e-3;            /* face centroids are thirds summed in double: far inside this */
    double a[3], b[3];
    int k = 0;
    for (k = 0; k < 3; k++) {
        a[k] = blo[k] - eps;
        b[k] = bhi[k] + eps;
        if (!(a[k] <= b[k]) || b[k] < lo[k] || a[k] > hi[k]) return 0;
    }
    return nreg == NULL || ic_near_meets(nreg, a, b);
}

/* After a whole read of a mesh that had no current chart index: write it from the data in hand (the raw per-vertex
 * charts), so the next close-up on this mesh pages instead of reading it all.  A failure only costs that. */
static void ic_index_after_read(const MeshChartStamp *st, const char *vmesh, const char *chart, const MeshBinData *m,
                                const int32_t *vchart)
{
    char ipath[2700];
    double t0 = ves_clock_sec();
    int rc = 0;
    if (MeshChartIndex_path(vmesh, ipath, sizeof ipath) != 0) return;
    rc = MeshChartIndex_build_mem(st, vmesh, chart, m, vchart, ipath);
    if (rc == 0) fprintf(stderr, "  wrote the chart index %s for the next close-up (%.1fs)\n", ipath, ves_clock_sec() - t0);
    else if (rc == 1) fprintf(stderr, "  %s is not chart-contiguous: no chart index\n", vmesh);
    else fprintf(stderr, "  no chart index written for %s (it changed during the read, or no memory)\n", vmesh);
}

/* The layout's charts that can hold a face of the grown seed box, when the layout has a current chart index: a face
 * whose UV centroid lies in the box belongs to a chart whose UV box meets it, and the paged mesh keeps the file's
 * order, so the surface stage's own face tests see exactly the faces of a whole-mesh read.  1 = no current index
 * (read it all), 0 = paged, -1 = error. */
static int ic_layout_paged(const IcConfig *c, const char *vmesh, MeshBinData *m, int32_t **vchart)
{
    char ipath[2700];
    MeshChartIndex ix;
    size_t *sel = NULL, ns = 0, k = 0;
    const double eps = 1e-6;
    double gb[4], t0 = ves_clock_sec();
    int rc = 1;
    gb[0] = c->seed[0] - c->grow_px; gb[1] = c->seed[1] - c->grow_px;
    gb[2] = c->seed[2] + c->grow_px; gb[3] = c->seed[3] + c->grow_px;
    if (MeshChartIndex_path(vmesh, ipath, sizeof ipath) != 0 || MeshChartIndex_load(ipath, vmesh, &ix) != 0) return 1;
    if (!(ix.mesh.flags & MESH_BIN_HAS_UV)) { MeshChartIndex_free(&ix); return 1; }
    sel = (size_t *)malloc((ix.n ? ix.n : 1) * sizeof *sel);
    if (sel == NULL) { MeshChartIndex_free(&ix); return -1; }
    for (k = 0; k < ix.n; k++) {
        const MeshChartRec *r = &ix.rec[k];
        if (r->f1 == r->f0) continue;
        if (r->uvhi[0] < gb[0] - eps || r->uvlo[0] >= gb[2] + eps || r->uvhi[1] < gb[1] - eps || r->uvlo[1] >= gb[3] + eps)
            continue;
        sel[ns++] = k;
    }
    rc = MeshChartIndex_read(&ix, vmesh, sel, ns, m, vchart) == 0 ? 0 : -1;
    if (rc == 0)
        fprintf(stderr, "  layout paged (%s): %zu of %zu chart runs, %zu of %llu verts, %zu of %llu faces (%.1fs)\n",
                ipath, ns, ix.n, m->nv, (unsigned long long)ix.mesh.nv, m->nf, (unsigned long long)ix.mesh.nf,
                ves_clock_sec() - t0);
    free(sel);
    MeshChartIndex_free(&ix);
    return rc;
}

typedef struct IcAdoptAudit {
    const IcPool *pool;
    char *buf;                     /* the CSV, built in memory (fprintf's bytes) and written in one go */
    size_t n, cap;
    int failed;
} IcAdoptAudit;

/* A file's bytes, written later (a chain's pass files, ic_write_job_run). */
typedef struct IcFileBytes {
    char path[2700];
    char *data;
    size_t n;
} IcFileBytes;

static void ic_audit_printf(IcAdoptAudit *a, const char *fmt, ...)
{
    va_list ap;
    int need = 0;
    if (a->failed) return;
    for (;;) {
        va_start(ap, fmt);
        need = vsnprintf(a->buf + a->n, a->cap - a->n, fmt, ap);
        va_end(ap);
        if (need < 0) { a->failed = 1; return; }
        if ((size_t)need < a->cap - a->n) { a->n += (size_t)need; return; }
        {
            size_t cap = a->cap * 2 + (size_t)need + 1;
            char *grown = (char *)realloc(a->buf, cap);
            if (grown == NULL) { a->failed = 1; return; }
            a->buf = grown; a->cap = cap;
        }
    }
}

static int ic_adopt_observer(const CloseupAdoptAudit *record, void *context)
{
    IcAdoptAudit *audit = (IcAdoptAudit *)context;
    const IcPool *pool = audit->pool;
    const CloseupPlacement *p = &record->p;
    double lo[3] = {1e300, 1e300, 1e300}, hi[3] = {-1e300, -1e300, -1e300};
    double ulo = 1e300, uhi = -1e300, vlo = 1e300, vhi = -1e300;
    size_t v = 0;
    int q = 0;
    if (record->screen_hits == 0 && p->corr == 0) return 0;
    for (v = pool->vofs[p->chart]; v < pool->vofs[p->chart + 1]; v++) {
        for (q = 0; q < 3; q++) {
            double x = pool->verts[v * 3 + (size_t)q];
            lo[q] = fmin(lo[q], x); hi[q] = fmax(hi[q], x);
        }
        if (record->reason >= 2) {
            double u = 0.0, w = 0.0;
            CloseupPlacement_apply(p, pool->uv[v * 2], pool->uv[v * 2 + 1], &u, &w);
            ulo = fmin(ulo, u); uhi = fmax(uhi, u); vlo = fmin(vlo, w); vhi = fmax(vhi, w);
        }
    }
    if (record->reason < 2) ulo = uhi = vlo = vhi = NAN;
    ic_audit_printf(audit, "%zu,%d,%d,%zu,%zu,%.9g,%.9g,%zu,%zu,%zu,%zu,%.9g,%.9g,%zu,"
                   "%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%zu,%zu,%.9g\n",
                   p->chart, pool->chart_id[p->chart], record->reason, record->screen_hits, p->corr, p->rms, p->extent,
                   record->data_landed, record->data_agree, record->fill_landed, record->fill_agree,
                   record->data_max_error, record->fill_max_error, p->fresh,
                   lo[0], hi[0], lo[1], hi[1], lo[2], hi[2], ulo, uhi, vlo, vhi,
                   record->ignored_fill_landed, record->ignored_fill_agree, record->ignored_fill_max_error);
    return audit->failed ? -1 : 0;
}

/* later: NULL writes the audit CSV now; else it receives the CSV's bytes for the caller to write (a chain's writer). */
static int ic_adopt_round(const IcConfig *c, const CloseupGrid *g, const IcPool *pool, const uint8_t *done,
                          int pool_index, int round, CloseupPlacement **out, size_t *nout, CloseupAdoptStats *st,
                          IcFileBytes *later)
{
    char path[2600];
    IcAdoptAudit audit;
    int rc = 0;
    if (!c->adopt_audit) return CloseupAdopt_round(g, &pool->p, done, &c->adopt, out, nout, st);
    snprintf(path, sizeof path, "%s/%s_adopt_pool%d_round%d.csv", c->out_dir, c->id, pool_index, round);
    memset(&audit, 0, sizeof audit);
    audit.pool = pool;
    audit.cap = (size_t)1 << 20;
    audit.buf = (char *)malloc(audit.cap);
    if (audit.buf == NULL) return -1;
    ic_audit_printf(&audit, "pool_chart,source_chart,reason,screen_hits,corr,rms,extent,data_landed,data_agree,"
                    "fill_landed,fill_agree,data_max_error,fill_max_error,fresh,z_min,z_max,y_min,y_max,x_min,x_max,"
                    "u_min,u_max,v_min,v_max,ignored_fill_landed,ignored_fill_agree,ignored_fill_max_error\n");
    if (audit.failed) rc = -1;
    if (rc == 0) rc = CloseupAdopt_round_audited(g, &pool->p, done, &c->adopt, out, nout, st, ic_adopt_observer, &audit);
    if (rc == 0 && later != NULL) {
        snprintf(later->path, sizeof later->path, "%s", path);
        later->data = audit.buf; later->n = audit.n;
        audit.buf = NULL;
    } else if (rc == 0) {
        FILE *f = NULL;
        ves_ensure_parent_dir(path);
        f = fopen(path, "wb");
        if (f == NULL || (audit.n && fwrite(audit.buf, 1, audit.n, f) != audit.n)) rc = -1;
        if (f != NULL && fclose(f) != 0) rc = -1;
    }
    free(audit.buf);
    if (rc != 0) { free(*out); *out = NULL; *nout = 0; }
    return rc;
}

/* The pool's candidate charts (ic_pool_candidate), when the pool mesh has a current chart index.  A chart id >= 0 owns
 * one record, whose box is the box of all its vertices, and the paged mesh keeps the file's order, so the selection in
 * ic_pool_load sees exactly the faces of a whole-mesh read.  1 = no current index, 0 = paged (m, vchart and, for a
 * surface pool, attr filled), -1 = error. */
static int ic_pool_paged(const char *vmesh, const char *attr_path, const double lo[3], const double hi[3],
                         const IcNear *nreg, const uint8_t *exclude, size_t nexcl, MeshBinData *m, int32_t **vchart,
                         uint8_t **attr)
{
    char ipath[2700];
    MeshChartIndex ix;
    size_t *sel = NULL, ns = 0, k = 0;
    uint64_t box_faces = 0;
    double t0 = ves_clock_sec();
    int rc = 1;
    if (MeshChartIndex_path(vmesh, ipath, sizeof ipath) != 0 || MeshChartIndex_load(ipath, vmesh, &ix) != 0) return 1;
    if (!(ix.mesh.flags & MESH_BIN_HAS_UV)) { MeshChartIndex_free(&ix); return 1; }
    sel = (size_t *)malloc((ix.n ? ix.n : 1) * sizeof *sel);
    if (sel == NULL) { MeshChartIndex_free(&ix); return -1; }
    for (k = 0; k < ix.n; k++) {
        const MeshChartRec *r = &ix.rec[k];
        if (r->f1 == r->f0 || r->chart < 0 || ((size_t)r->chart < nexcl && exclude[r->chart])) continue;
        if (!ic_pool_candidate(r->lo, r->hi, lo, hi, NULL)) continue;
        box_faces += r->f1 - r->f0;                  /* the search box alone (diagnostic) */
        if (nreg != NULL && !ic_pool_candidate(r->lo, r->hi, lo, hi, nreg)) continue;
        sel[ns++] = k;
    }
    if (nreg != NULL)
        fprintf(stderr, "  pool near region: charts meeting the search box hold %llu faces\n", (unsigned long long)box_faces);
    rc = MeshChartIndex_read(&ix, vmesh, sel, ns, m, vchart) == 0 ? 0 : -1;
    if (rc == 0 && attr_path != NULL && MeshChartIndex_read_side(&ix, attr_path, 1, sel, ns, (void **)attr) != 0) rc = -1;
    if (rc == 0)
        fprintf(stderr, "  pool paged (%s): %zu of %zu chart runs, %zu of %llu verts, %zu of %llu faces (%.1fs)\n",
                ipath, ns, ix.n, m->nv, (unsigned long long)ix.mesh.nv, m->nf, (unsigned long long)ix.mesh.nf,
                ves_clock_sec() - t0);
    else if (rc < 0) { MeshBin_dispose(m); free(*vchart); *vchart = NULL; }
    free(sel);
    MeshChartIndex_free(&ix);
    return rc;
}

/* A pool layout as read over a search box: its paged (or whole) mesh, per-vertex charts and, for a close-up surface
 * pool, its attributes. */
typedef struct IcPoolMesh {
    MeshBinData m;
    int32_t *vchart;
    uint8_t *attr;
} IcPoolMesh;

static void ic_pool_mesh_free(IcPoolMesh *pm)
{
    MeshBin_dispose(&pm->m); free(pm->vchart); free(pm->attr);
    memset(pm, 0, sizeof *pm);
}

/* The pool's mesh over [lo, hi]: the chart index's runs meeting it (ic_pool_paged), else the whole mesh, then its
 * index for the next close-up.  0 read, -1 error. */
static int ic_pool_read(const char *dir, int surface, const double lo[3], const double hi[3], const IcNear *nreg,
                        const uint8_t *exclude, size_t nexcl, IcPoolMesh *pm)
{
    char path[2600], apath[2600];
    size_t got = 0;
    int pg = 1;
    double tpl = ves_clock_sec();
    memset(pm, 0, sizeof *pm);
    snprintf(path, sizeof path, surface ? "%s_surface.vmesh" : "%s/stage6_repaired.vmesh", dir);
    snprintf(apath, sizeof apath, "%s_surface_attr.u8", dir);
    pg = ic_pool_paged(path, surface ? apath : NULL, lo, hi, nreg, exclude, nexcl, &pm->m, &pm->vchart, &pm->attr);
    if (pg < 0) { fprintf(stderr, "ERROR: cannot page the pool layout %s\n", path); return -1; }
    fprintf(stderr, "  [complete-timing] pool paged read %.2fs\n", ves_clock_sec() - tpl);
    if (pg > 0) {                                    /* no current chart index: the whole mesh, then its index */
        char cpath[2600];
        MeshChartStamp stamp;
        int stamped = 0;
        snprintf(cpath, sizeof cpath, surface ? "%s_surface_chart.i32" : "%s/stage6_repaired_chart.i32", dir);
        stamped = !surface && MeshChartIndex_stamp(path, cpath, &stamp) == 0;    /* a close-up surface is small */
        if (MeshBin_read_precise_malloc(path, &pm->m) != 0 || (pm->m.uv == NULL && pm->m.uv64 == NULL)) {
            fprintf(stderr, "ERROR: cannot read the pool layout %s\n", path);
            ic_pool_mesh_free(pm);
            return -1;
        }
        pm->vchart = (int32_t *)ic_read_file(cpath, pm->m.nv * sizeof(int32_t), &got);
        if (pm->vchart == NULL) { fprintf(stderr, "ERROR: %s missing or short\n", cpath); ic_pool_mesh_free(pm); return -1; }
        if (stamped) ic_index_after_read(&stamp, path, cpath, &pm->m, pm->vchart);  /* the raw charts, before any mask */
        if (surface) {
            pm->attr = (uint8_t *)ic_read_file(apath, pm->m.nv, &got);
            if (!pm->attr) { ic_pool_mesh_free(pm); return -1; }
        }
    }
    return 0;
}

/* The pool over [lo, hi] from its mesh: charts with a face whose centroid lies in the box, in chart-id order, their
 * faces in file order and vertices numbered as those faces first name them.  A mesh read over a larger box selects
 * the same pool: a run that misses [lo, hi] holds no such face (its centroids lie in its vertices' box).
 * lazy_inward: the inward field is filled per chart when an adopt round first evaluates it (CloseupPool.inward_fill),
 * instead of for every pool vertex here; the values are the same either way. */
static int ic_pool_select(const IcPoolMesh *pm, int surface, const AsmAxis *axis, const double lo[3], const double hi[3],
                          const IcNear *nreg, const uint8_t *exclude, size_t nexcl, int lazy_inward, IcPoolBufs *spare,
                          IcPool *q)
{
    const MeshBinData m = pm->m;
    const uint8_t *attr = pm->attr;
    int32_t *vchart = pm->vchart, *masked = NULL, *vmap = NULL;
    uint8_t *elig = NULL, *infc = NULL;
    size_t f = 0, v = 0, maxid = 0, *cnt = NULL, *pos = NULL, nsel = 0, *order = NULL, nch = 0, k = 0, nelig = 0;
    int rc = -1, own = 0;
    double tpl = ves_clock_sec();
    memset(q, 0, sizeof *q);
    if (nreg != NULL) {                /* the candidate charts, from the boxes of all their vertices (before the mask) */
        double *cbox = NULL;
        for (v = 0; v < m.nv; v++) if (vchart[v] >= 0 && (size_t)vchart[v] + 1 > nelig) nelig = (size_t)vchart[v] + 1;
        cbox = (double *)malloc((nelig ? nelig : 1) * 6 * sizeof *cbox);
        elig = (uint8_t *)calloc(nelig ? nelig : 1, 1);
        if (cbox == NULL || elig == NULL) { free(cbox); goto out; }
        for (k = 0; k < nelig; k++) {
            cbox[k * 6] = cbox[k * 6 + 1] = cbox[k * 6 + 2] = INFINITY;
            cbox[k * 6 + 3] = cbox[k * 6 + 4] = cbox[k * 6 + 5] = -INFINITY;
        }
        for (v = 0; v < m.nv; v++) {
            size_t b = 0, kk = 0;
            if (vchart[v] < 0) continue;
            b = (size_t)vchart[v] * 6;
            for (kk = 0; kk < 3; kk++) {
                double x = (double)m.verts[v * 3 + kk];
                if (x < cbox[b + kk]) cbox[b + kk] = x;
                if (x > cbox[b + 3 + kk]) cbox[b + 3 + kk] = x;
            }
        }
        for (k = 0; k < nelig; k++) elig[k] = (uint8_t)ic_pool_candidate(&cbox[k * 6], &cbox[k * 6 + 3], lo, hi, nreg);
        free(cbox);
    }
    if (surface) {                     /* Filled/extrapolated geometry is never evidence for a connection. */
        masked = (int32_t *)malloc((m.nv ? m.nv : 1) * sizeof *masked);
        if (masked == NULL) goto out;
        for (v = 0; v < m.nv; v++) masked[v] = attr[v] != 1 && attr[v] != 3 ? -1 : vchart[v];
        vchart = masked;
    }
    for (v = 0; v < m.nv; v++) if (vchart[v] >= 0 && (size_t)vchart[v] + 1 > maxid) maxid = (size_t)vchart[v] + 1;
    cnt = (size_t *)calloc(maxid + 1, sizeof *cnt);
    pos = (size_t *)calloc(maxid + 1, sizeof *pos);
    order = (size_t *)malloc((m.nf ? m.nf : 1) * sizeof *order);
    vmap = (int32_t *)malloc(m.nv * sizeof *vmap);
    infc = (uint8_t *)malloc(m.nf ? m.nf : 1);
    if (cnt == NULL || pos == NULL || order == NULL || vmap == NULL || infc == NULL) goto out;
    {   /* each face's verdict in parallel (a face's centroid test reads only its own corners); the counts and the
         * chart-ordered face lists stay one pass in file order */
        long lf = 0;
        int mixed = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 65536) reduction(| : mixed)
#endif
        for (lf = 0; lf < (long)m.nf; lf++) {
            const size_t ff = (size_t)lf;
            int32_t ch = vchart[m.faces[ff * 3]];
            double c3[3] = { 0.0, 0.0, 0.0 };
            int kk = 0, in = 1;
            infc[ff] = 0;
            if (ch < 0 || vchart[m.faces[ff * 3 + 1]] < 0 || vchart[m.faces[ff * 3 + 2]] < 0 ||
                ((size_t)ch < nexcl && exclude[ch]) || (elig != NULL && !elig[ch])) continue;
            mixed |= vchart[m.faces[ff * 3 + 1]] != ch || vchart[m.faces[ff * 3 + 2]] != ch;
            for (kk = 0; kk < 3; kk++) {
                size_t w = (size_t)m.faces[ff * 3 + (size_t)kk];
                c3[0] += m.verts[w * 3] / 3.0; c3[1] += m.verts[w * 3 + 1] / 3.0; c3[2] += m.verts[w * 3 + 2] / 3.0;
            }
            for (kk = 0; kk < 3; kk++) if (c3[kk] < lo[kk] || c3[kk] > hi[kk]) in = 0;
            infc[ff] = (uint8_t)in;
        }
        own = !mixed;      /* every candidate face on its own chart's vertices: charts number their vertices alone */
    }
    for (f = 0; f < m.nf; f++) if (infc[f]) cnt[vchart[m.faces[f * 3]]]++;
    for (k = 0; k < maxid; k++) { pos[k] = nsel; nsel += cnt[k]; if (cnt[k]) nch++; }
    for (f = 0; f < m.nf; f++) if (infc[f]) { int32_t ch = vchart[m.faces[f * 3]]; order[pos[ch]++] = f; }
    fprintf(stderr, "  [complete-timing] pool select %.2fs\n", ves_clock_sec() - tpl);
    /* charts in id order, vertices numbered in chart order (charts share no vertices) */
    q->vofs = (size_t *)ic_buf(spare, IB_VOFS, (nch + 1) * sizeof *q->vofs, &q->bytes[IB_VOFS]);
    q->fofs = (size_t *)ic_buf(spare, IB_FOFS, (nch + 1) * sizeof *q->fofs, &q->bytes[IB_FOFS]);
    q->chart_id = (int32_t *)ic_buf(spare, IB_CHART_ID, (nch ? nch : 1) * sizeof *q->chart_id, &q->bytes[IB_CHART_ID]);
    q->faces = (int32_t *)ic_buf(spare, IB_FACES, (nsel ? nsel : 1) * 3 * sizeof *q->faces, &q->bytes[IB_FACES]);
    q->face_chart = (int32_t *)ic_buf(spare, IB_FACE_CHART, (nsel ? nsel : 1) * sizeof *q->face_chart,
                                      &q->bytes[IB_FACE_CHART]);
    if (q->vofs == NULL || q->fofs == NULL || q->chart_id == NULL || q->faces == NULL || q->face_chart == NULL) goto out;
    {
        long lv = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 65536)
#endif
        for (lv = 0; lv < (long)m.nv; lv++) vmap[lv] = -1;
    }
    if (own) {   /* the serial numbering below, chart by chart in parallel: a chart's ids are its first-seen order */
        size_t nv = 0, c = 0, fi = 0;
        size_t *nloc = (size_t *)malloc((nch + 1) * sizeof *nloc);
        long lc = 0;
        if (nloc == NULL) goto out;
        for (k = 0; k < maxid; k++) {
            if (cnt[k] == 0) continue;
            q->fofs[c] = fi; q->chart_id[c] = (int32_t)k;
            fi += cnt[k];
            c++;
        }
        q->fofs[c] = fi;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 64)
#endif
        for (lc = 0; lc < (long)nch; lc++) {
            const size_t ck = (size_t)q->chart_id[lc];
            size_t a = 0, local = 0;
            for (a = pos[ck] - cnt[ck]; a < pos[ck]; a++) {
                int kk = 0;
                for (kk = 0; kk < 3; kk++) {
                    const size_t w = (size_t)m.faces[order[a] * 3 + (size_t)kk];
                    if (vmap[w] < 0) vmap[w] = (int32_t)local++;
                }
            }
            nloc[lc] = local;
        }
        for (c = 0; c < nch; c++) { q->vofs[c] = nv; nv += nloc[c]; }
        q->vofs[nch] = nv;
        free(nloc);
        q->verts = (float *)ic_buf(spare, IB_VERTS, (nv ? nv : 1) * 3 * sizeof *q->verts, &q->bytes[IB_VERTS]);
        q->uv = (double *)ic_buf(spare, IB_UV, (nv ? nv : 1) * 2 * sizeof *q->uv, &q->bytes[IB_UV]);
        q->inward = (float *)ic_buf(spare, IB_INWARD, (nv ? nv : 1) * 3 * sizeof *q->inward, &q->bytes[IB_INWARD]);
        if (q->verts == NULL || q->uv == NULL || q->inward == NULL) goto out;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 64)
#endif
        for (lc = 0; lc < (long)nch; lc++) {   /* global ids; each vertex copied as its first face names it */
            const size_t ck = (size_t)q->chart_id[lc], base = q->vofs[lc];
            size_t a = 0, fo = q->fofs[lc];
            for (a = pos[ck] - cnt[ck]; a < pos[ck]; a++, fo++) {
                int kk = 0;
                for (kk = 0; kk < 3; kk++) {
                    const size_t w = (size_t)m.faces[order[a] * 3 + (size_t)kk];
                    int32_t id = vmap[w];
                    if (id >= 0) {
                        const size_t gi = base + (size_t)id;
                        memcpy(&q->verts[gi * 3], &m.verts[w * 3], 3 * sizeof *q->verts);
                        q->uv[gi * 2] = MeshBin_uv(&m, w * 2); q->uv[gi * 2 + 1] = MeshBin_uv(&m, w * 2 + 1);
                        vmap[w] = -(int32_t)gi - 2;
                        id = (int32_t)gi;
                    } else {
                        id = -id - 2;
                    }
                    q->faces[fo * 3 + (size_t)kk] = id;
                }
                q->face_chart[fo] = (int32_t)ck;
            }
        }
        goto numbered;
    }
    {
        size_t nv = 0, fi = 0, c = 0;
        for (k = 0; k < maxid; k++) {
            size_t a = 0;
            if (cnt[k] == 0) continue;
            q->vofs[c] = nv; q->fofs[c] = fi; q->chart_id[c] = (int32_t)k;
            for (a = pos[k] - cnt[k]; a < pos[k]; a++) {
                size_t ff = order[a];
                int kk = 0;
                for (kk = 0; kk < 3; kk++) {
                    size_t w = (size_t)m.faces[ff * 3 + (size_t)kk];
                    if (vmap[w] < 0) vmap[w] = (int32_t)nv++;
                    q->faces[fi * 3 + (size_t)kk] = vmap[w];
                }
                q->face_chart[fi] = (int32_t)k;
                fi++;
            }
            c++;
        }
        q->vofs[c] = nv; q->fofs[c] = fi;
        q->verts = (float *)ic_buf(spare, IB_VERTS, (nv ? nv : 1) * 3 * sizeof *q->verts, &q->bytes[IB_VERTS]);
        q->uv = (double *)ic_buf(spare, IB_UV, (nv ? nv : 1) * 2 * sizeof *q->uv, &q->bytes[IB_UV]);
        q->inward = (float *)ic_buf(spare, IB_INWARD, (nv ? nv : 1) * 3 * sizeof *q->inward, &q->bytes[IB_INWARD]);
        if (q->verts == NULL || q->uv == NULL || q->inward == NULL) goto out;
        {
            long lv = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 65536)
#endif
            for (lv = 0; lv < (long)m.nv; lv++) {      /* vmap is one-to-one: every target written once */
                size_t w = 0, vv = (size_t)lv;
                if (vmap[vv] < 0) continue;
                w = (size_t)vmap[vv];
                memcpy(&q->verts[w * 3], &m.verts[vv * 3], 3 * sizeof *q->verts);
                q->uv[w * 2] = MeshBin_uv(&m, vv * 2); q->uv[w * 2 + 1] = MeshBin_uv(&m, vv * 2 + 1);
            }
        }
    }
numbered:
    {
        const size_t nv = q->vofs[nch], fi = q->fofs[nch];
        {
            double tp = ves_clock_sec();
            if (lazy_inward) {           /* both fields per chart when a round first evaluates it */
                q->arena = Arena_new_sized((size_t)1 << 20);
                q->index = AsmAxis_projector(q->arena, axis);
                q->nrm = (float *)ic_buf(spare, IB_NRM, (nv ? nv : 1) * 3 * sizeof *q->nrm, &q->bytes[IB_NRM]);
                q->chart_ready = (uint8_t *)ic_buf(spare, IB_READY, nch ? nch : 1, &q->bytes[IB_READY]);
                if (q->nrm == NULL || q->chart_ready == NULL || q->index == NULL) goto out;
                memset(q->chart_ready, 0, nch ? nch : 1);
                q->p.chart_fill = ic_pool_chart_fill; q->p.chart_ctx = q; q->p.chart_ready = q->chart_ready;
            } else {
                q->nrm = MeshNormals_compute_runs(q->verts, nv, q->faces, q->vofs, q->fofs, nch);
                if (q->nrm == NULL) goto out;
                ic_axis_fields(axis, q->verts, nv, NULL, q->inward);
            }
            fprintf(stderr, "  [complete-timing] pool normals + axis fields %.2fs (%zu verts %zu faces%s)\n",
                    ves_clock_sec() - tp, nv, fi, lazy_inward ? ", per chart on first use" : "");
        }
        q->p.nv = nv; q->p.nf = fi; q->p.nch = nch;
    }
    q->p.verts = q->verts; q->p.nrm = q->nrm; q->p.inward = q->inward; q->p.uv = q->uv; q->p.faces = q->faces;
    q->p.vofs = q->vofs; q->p.fofs = q->fofs;
    rc = 0;
out:
    if (rc != 0) ic_pool_free(q);
    free(masked); free(elig); free(vmap); free(cnt); free(pos); free(order); free(infc);
    return rc;
}

static int ic_pool_load(const char *dir, int surface, const AsmAxis *axis, const double lo[3], const double hi[3],
                        const IcNear *nreg, const uint8_t *exclude, size_t nexcl, int lazy_inward, IcPool *q)
{
    IcPoolMesh pm;
    int rc = ic_pool_read(dir, surface, lo, hi, nreg, exclude, nexcl, &pm);
    memset(q, 0, sizeof *q);
    if (rc == 0) rc = ic_pool_select(&pm, surface, axis, lo, hi, nreg, exclude, nexcl, lazy_inward, NULL, q);
    ic_pool_mesh_free(&pm);
    return rc;
}

/* The completion chain's pool meshes (ic_stage_complete_chain): each read once over its first pass's search box plus
 * IC_POOL_CACHE_SLACK_VOX, then every later pass selects from it (ic_pool_select), reading again only when a pass's
 * box leaves it.  The HDD-resident pools of a 32^3 take-in cost minutes per cold read (2026-10-30). */
#define IC_POOL_CACHE_SLOTS 4
#define IC_POOL_CACHE_SLACK_VOX 512.0
typedef struct IcPoolCache {
    IcPoolBufs spare;              /* the last pass's pool arrays, for the next pass's (ic_pool_keep) */
    IcPoolMesh pm[IC_POOL_CACHE_SLOTS];
    char dir[IC_POOL_CACHE_SLOTS][2600];
    int surface[IC_POOL_CACHE_SLOTS], loaded[IC_POOL_CACHE_SLOTS];
    double lo[IC_POOL_CACHE_SLOTS][3], hi[IC_POOL_CACHE_SLOTS][3];
} IcPoolCache;

static void ic_pool_cache_free(IcPoolCache *pc)
{
    int k = 0;
    for (k = 0; k < IC_POOL_CACHE_SLOTS; k++) ic_pool_mesh_free(&pc->pm[k]);
    ic_bufs_free(&pc->spare);
    memset(pc, 0, sizeof *pc);
}

static int ic_pool_load_cached(IcPoolCache *pc, const char *dir, int surface, const AsmAxis *axis, const double lo[3],
                               const double hi[3], int lazy_inward, IcPool *q)
{
    int slot = -1, k = 0, d = 0;
    for (k = 0; k < IC_POOL_CACHE_SLOTS && slot < 0; k++)
        if (pc->loaded[k] && pc->surface[k] == surface && strcmp(pc->dir[k], dir) == 0) slot = k;
    if (slot >= 0) {
        int inside = 1;
        for (d = 0; d < 3; d++) inside = inside && lo[d] >= pc->lo[slot][d] && hi[d] <= pc->hi[slot][d];
        if (inside) {
            fprintf(stderr, "  pool %s from the chain's copy (read over %.0f..%.0f %.0f..%.0f %.0f..%.0f)\n", dir,
                    pc->lo[slot][0], pc->hi[slot][0], pc->lo[slot][1], pc->hi[slot][1], pc->lo[slot][2], pc->hi[slot][2]);
            return ic_pool_select(&pc->pm[slot], surface, axis, lo, hi, NULL, NULL, 0, lazy_inward, &pc->spare, q);
        }
    } else {
        for (k = 0; k < IC_POOL_CACHE_SLOTS && slot < 0; k++) if (!pc->loaded[k]) slot = k;
        if (slot < 0) slot = 0;
    }
    ic_pool_mesh_free(&pc->pm[slot]);
    pc->loaded[slot] = 0;
    for (d = 0; d < 3; d++) { pc->lo[slot][d] = lo[d] - IC_POOL_CACHE_SLACK_VOX; pc->hi[slot][d] = hi[d] + IC_POOL_CACHE_SLACK_VOX; }
    if (ic_pool_read(dir, surface, pc->lo[slot], pc->hi[slot], NULL, NULL, 0, &pc->pm[slot]) != 0) return -1;
    snprintf(pc->dir[slot], sizeof pc->dir[slot], "%s", dir);
    pc->surface[slot] = surface; pc->loaded[slot] = 1;
    fprintf(stderr, "  pool %s read once for the chain over the search box + %.0f vox\n", dir, IC_POOL_CACHE_SLACK_VOX);
    return ic_pool_select(&pc->pm[slot], surface, axis, lo, hi, NULL, NULL, 0, lazy_inward, &pc->spare, q);
}

/* ------------------------------------------------------------------ surface */

static int ic_stage_surface(const IcConfig *c)
{
    char path[2600], js[4800];
    MeshBinData m;
    Arena_T arena = NULL;
    const JsonValue *idx = NULL, *charts = NULL;
    const char *err = NULL;
    int32_t *vchart = NULL, *region = NULL, *face_chart = NULL, *vmap = NULL, *sfaces = NULL, *node_of = NULL;
    int32_t *mfaces = NULL, *vtx_chart = NULL;
    uint8_t *seed_region = NULL, *keep = NULL, *vattr = NULL;
    float *sverts = NULL, *sinward = NULL, *ninward = NULL, *mverts = NULL, *mnrm = NULL;
    double *suv = NULL, *fuv = NULL, *saxial = NULL, *uvn = NULL, *muv = NULL, *maxial = NULL;
    size_t got = 0, nchart_ids = 0, f = 0, v = 0, snv = 0, snf = 0, in_box = 0, mnv = 0, mnf = 0;
    size_t *reg_count = NULL;
    long max_region = -1;
    AsmAxis *axis = NULL;
    CloseupFrame fr, fr2;
    CloseupObs obs;
    CloseupObsGate gate;
    CloseupGrid g;
    CloseupFitStats fs;
    CloseupFitOpts fit;
    CloseupFlatStats flat;
    double t0 = ves_clock_sec(), gb[4], uvmin[2], uvmax[2], agree = 0.0;
    size_t used = 0, gated = 0, nseed = 0, nadopt_total = 0, adopt_obs = 0, nadopted_ids = 0;
    uint8_t *seedchart = NULL, *adopted = NULL;
    char adopt_json[8192], frame_json[4096];
    size_t adopt_len = 0, frame_len = 0;
    int orient = 0, rc = 1;
    FILE *fp = NULL;

    memset(&m, 0, sizeof m);
    memset(&obs, 0, sizeof obs);
    memset(&g, 0, sizeof g);
    if (c->complete_gaps) CloseupFitOpts_complete(&fit); else CloseupFitOpts_default(&fit);
    fit.bridge = c->bridge && c->chart_wrap != NULL; fit.chart_wrap = c->chart_wrap; fit.nchart_wrap = c->nchart_wrap;
    fit.bridge_px = fit.fill_px;    /* the bridge reaches as far as the fill in effect (24 px, or the completed-gap 160 px) */
    adopt_json[0] = 0;
    frame_json[0] = 0;
    arena = Arena_new();
    axis = AsmAxis_load(arena, c->axis);
    if (axis == NULL) { fprintf(stderr, "ERROR: cannot load axis table %s\n", c->axis); goto out; }
    snprintf(path, sizeof path, "%s/stage6_repaired.vmesh", c->layout);
    {
        int pg = ic_layout_paged(c, path, &m, &vchart);
        if (pg < 0) { fprintf(stderr, "ERROR: cannot page the layout %s\n", path); goto out; }
        if (pg > 0) {                                /* no current chart index: the whole layout, then its index */
            char cpath[2600];
            MeshChartStamp stamp;
            int stamped = 0;
            snprintf(cpath, sizeof cpath, "%s/stage6_repaired_chart.i32", c->layout);
            stamped = MeshChartIndex_stamp(path, cpath, &stamp) == 0;
            if (MeshBin_read_precise_malloc(path, &m) != 0 || (m.uv == NULL && m.uv64 == NULL)) {
                fprintf(stderr, "ERROR: cannot read a UV-carrying layout %s\n", path);
                goto out;
            }
            vchart = (int32_t *)ic_read_file(cpath, m.nv * sizeof(int32_t), &got);
            if (vchart == NULL) { fprintf(stderr, "ERROR: %s missing or not %zu vertices\n", cpath, m.nv); goto out; }
            if (stamped) ic_index_after_read(&stamp, path, cpath, &m, vchart);
        }
    }
    snprintf(path, sizeof path, "%s/sheet_index.json", c->layout);
    idx = Json_parse_file(arena, path, &err);
    charts = idx != NULL ? Json_object_get(idx, "charts") : NULL;
    if (charts == NULL || Json_type(charts) != JSON_ARRAY) {
        fprintf(stderr, "ERROR: cannot read the chart index %s (%s)\n", path, err ? err : "no charts");
        goto out;
    }
    for (f = 0; f < Json_array_len(charts); f++) {
        long id = Json_member_long(Json_array_get(charts, f), "id", -1);
        if (id >= 0 && (size_t)id + 1 > nchart_ids) nchart_ids = (size_t)id + 1;
    }
    region = (int32_t *)malloc((nchart_ids ? nchart_ids : 1) * sizeof *region);
    if (region == NULL) goto out;
    for (f = 0; f < nchart_ids; f++) region[f] = -1;
    for (f = 0; f < Json_array_len(charts); f++) {
        const JsonValue *ch = Json_array_get(charts, f);
        long id = Json_member_long(ch, "id", -1), rg = Json_member_long(ch, "region", -1);
        if (id < 0 || !Json_as_bool(Json_object_get(ch, "placed"), 0)) continue;
        region[id] = (int32_t)rg;
        if (rg > max_region) max_region = rg;
    }
    fprintf(stderr, "[ink_closeup] %s: layout %zu verts %zu faces, %zu charts indexed (%.1fs)\n", c->id, m.nv, m.nf,
            nchart_ids, ves_clock_sec() - t0);

    /* the seed pieces: regions owning at least 2% (and 1,000) of the faces centred in the seed box */
    reg_count = (size_t *)calloc((size_t)(max_region + 2), sizeof *reg_count);
    seed_region = (uint8_t *)calloc((size_t)(max_region + 2), 1);
    keep = (uint8_t *)calloc(m.nf, 1);
    if (reg_count == NULL || seed_region == NULL || keep == NULL) goto out;
    for (f = 0; f < m.nf; f++) {
        double fu = 0.0, fv = 0.0;
        int32_t ch = vchart[m.faces[f * 3]];
        int k = 0;
        for (k = 0; k < 3; k++) {
            fu += MeshBin_uv(&m, (size_t)m.faces[f * 3 + (size_t)k] * 2) / 3.0;
            fv += MeshBin_uv(&m, (size_t)m.faces[f * 3 + (size_t)k] * 2 + 1) / 3.0;
        }
        if (fu < c->seed[0] || fu >= c->seed[2] || fv < c->seed[1] || fv >= c->seed[3]) continue;
        in_box++;
        if (ch >= 0 && (size_t)ch < nchart_ids && region[ch] >= 0) reg_count[region[ch]]++;
    }
    for (f = 0; f <= (size_t)max_region; f++)
        if (ic_seed_selected(reg_count[f], in_box, c->seed_all_pieces)) {
            seed_region[f] = 1;
            nseed++;
            fprintf(stderr, "  seed piece: region %zu (%zu of %zu faces in the box)\n", f, reg_count[f], in_box);
        }
    if (nseed == 0) { fprintf(stderr, "ERROR: no piece owns the seed box\n"); goto out; }
    gb[0] = c->seed[0] - c->grow_px; gb[1] = c->seed[1] - c->grow_px;
    gb[2] = c->seed[2] + c->grow_px; gb[3] = c->seed[3] + c->grow_px;
    vmap = (int32_t *)malloc(m.nv * sizeof *vmap);
    if (vmap == NULL) goto out;
    for (v = 0; v < m.nv; v++) vmap[v] = -1;
    for (f = 0; f < m.nf; f++) {
        double fu = 0.0, fv = 0.0;
        int32_t ch = vchart[m.faces[f * 3]];
        int k = 0;
        if (ch < 0 || (size_t)ch >= nchart_ids || region[ch] < 0 || !seed_region[region[ch]]) continue;
        for (k = 0; k < 3; k++) {
            fu += MeshBin_uv(&m, (size_t)m.faces[f * 3 + (size_t)k] * 2) / 3.0;
            fv += MeshBin_uv(&m, (size_t)m.faces[f * 3 + (size_t)k] * 2 + 1) / 3.0;
        }
        if (fu < gb[0] || fu >= gb[2] || fv < gb[1] || fv >= gb[3]) continue;
        keep[f] = 1;
        snf++;
        for (k = 0; k < 3; k++) {
            size_t w = (size_t)m.faces[f * 3 + (size_t)k];
            if (vmap[w] < 0) vmap[w] = (int32_t)snv++;
        }
    }
    sverts = (float *)malloc(snv * 3 * sizeof *sverts);
    suv = (double *)malloc(snv * 2 * sizeof *suv);
    fuv = (double *)malloc(snv * 2 * sizeof *fuv);
    saxial = (double *)malloc(snv * sizeof *saxial);
    sinward = (float *)malloc(snv * 3 * sizeof *sinward);
    sfaces = (int32_t *)malloc(snf * 3 * sizeof *sfaces);
    face_chart = (int32_t *)malloc(snf * sizeof *face_chart);
    if (sverts == NULL || suv == NULL || fuv == NULL || saxial == NULL || sinward == NULL || sfaces == NULL ||
        face_chart == NULL)
        goto out;
    for (v = 0; v < m.nv; v++) {
        size_t w = 0;
        if (vmap[v] < 0) continue;
        w = (size_t)vmap[v];
        memcpy(&sverts[w * 3], &m.verts[v * 3], 3 * sizeof *sverts);
        suv[w * 2] = MeshBin_uv(&m, v * 2); suv[w * 2 + 1] = MeshBin_uv(&m, v * 2 + 1);
    }
    for (f = 0, snf = 0; f < m.nf; f++) {
        int k = 0;
        if (!keep[f]) continue;
        for (k = 0; k < 3; k++) sfaces[snf * 3 + (size_t)k] = vmap[m.faces[f * 3 + (size_t)k]];
        face_chart[snf] = vchart[m.faces[f * 3]];
        snf++;
    }
    fprintf(stderr, "  close-up source: %zu faces, %zu verts within %.0f px of the seed box (%.1fs)\n", snf, snv,
            c->grow_px, ves_clock_sec() - t0);
    MeshBin_dispose(&m);
    free(vchart); vchart = NULL;

    /* the frame, PER SEED PIECE: the ribbon places pieces mirror-imaged or rotated independently, so each piece
     * gets its own reflection (to the team's reading handedness) and rotation (scroll axis down the rows); the
     * pieces are then laid side by side, in their layout order, CLOSEUP_PIECE_GAP_PX apart */
    ic_axis_fields(axis, sverts, snv, saxial, sinward);
    ic_axial_sense(c, saxial, snv);
    if (!c->z_top_to_bottom)
        fprintf(stderr, "  scan z runs from the scroll's bottom to its top: frames put high z at the top of the rows\n");
    uvmin[0] = uvmin[1] = 1e300; uvmax[0] = uvmax[1] = -1e300;
    {
        uint8_t *fuse = (uint8_t *)calloc(snf + 1, 1), *vdone = (uint8_t *)calloc(snv + 1, 1);
        size_t *order = (size_t *)malloc(((size_t)max_region + 2) * sizeof *order), no = 0, a = 0, b = 0;
        double *mu = (double *)calloc((size_t)max_region + 2, sizeof *mu), cursor = 0.0;
        if (fuse == NULL || vdone == NULL || order == NULL || mu == NULL) {
            free(fuse); free(vdone); free(order); free(mu);
            goto out;
        }
        for (f = 0; f <= (size_t)max_region; f++) {
            size_t cnt = 0;
            if (!seed_region[f]) continue;
            for (a = 0; a < snf; a++)
                if (face_chart[a] >= 0 && region[face_chart[a]] == (int32_t)f) {
                    mu[f] += suv[(size_t)sfaces[a * 3] * 2];
                    cnt++;
                }
            mu[f] = cnt ? mu[f] / (double)cnt : 0.0;
            order[no++] = f;
        }
        for (a = 0; a < no; a++)             /* layout order: left to right */
            for (b = a + 1; b < no; b++)
                if (mu[order[b]] < mu[order[a]]) { size_t t = order[a]; order[a] = order[b]; order[b] = t; }
        for (a = 0; a < no; a++) {
            size_t rg = order[a];
            double lo2[2] = { 1e300, 1e300 }, hi2[2] = { -1e300, -1e300 }, bu = 0.0, bv = 0.0, shift_u = 0.0, shift_v = 0.0;
            for (f = 0; f < snf; f++) fuse[f] = (uint8_t)(face_chart[f] >= 0 && region[face_chart[f]] == (int32_t)rg);
            for (f = 0; f < snf; f++) {
                int q = 0;
                if (!fuse[f]) continue;
                for (q = 0; q < 3; q++) {
                    size_t w = (size_t)sfaces[f * 3 + (size_t)q];
                    if (suv[w * 2] < lo2[0]) lo2[0] = suv[w * 2];
                    if (suv[w * 2] > hi2[0]) hi2[0] = suv[w * 2];
                    if (suv[w * 2 + 1] < lo2[1]) lo2[1] = suv[w * 2 + 1];
                    if (suv[w * 2 + 1] > hi2[1]) hi2[1] = suv[w * 2 + 1];
                }
            }
            bu = 0.5 * (lo2[0] + hi2[0]); bv = 0.5 * (lo2[1] + hi2[1]);
            if (CloseupSurface_frame(sverts, suv, sfaces, snf, fuse, saxial, sinward, bu, bv, &fr) != 0) {
                fprintf(stderr, "ERROR: seed piece %zu gives no frame\n", rg);
                free(fuse); free(vdone); free(order); free(mu);
                goto out;
            }
            lo2[0] = lo2[1] = 1e300; hi2[0] = hi2[1] = -1e300;
            for (f = 0; f < snf; f++) {
                int q = 0;
                if (!fuse[f]) continue;
                for (q = 0; q < 3; q++) {
                    size_t w = (size_t)sfaces[f * 3 + (size_t)q];
                    if (!vdone[w]) {
                        CloseupFrame_apply(&fr, suv[w * 2], suv[w * 2 + 1], &fuv[w * 2], &fuv[w * 2 + 1]);
                        vdone[w] = 1;
                    }
                    if (fuv[w * 2] < lo2[0]) lo2[0] = fuv[w * 2];
                    if (fuv[w * 2] > hi2[0]) hi2[0] = fuv[w * 2];
                    if (fuv[w * 2 + 1] < lo2[1]) lo2[1] = fuv[w * 2 + 1];
                    if (fuv[w * 2 + 1] > hi2[1]) hi2[1] = fuv[w * 2 + 1];
                }
            }
            shift_u = cursor - lo2[0];
            shift_v = -lo2[1];
            for (f = 0; f < snf; f++) {
                int q = 0;
                if (!fuse[f]) continue;
                for (q = 0; q < 3; q++) {
                    size_t w = (size_t)sfaces[f * 3 + (size_t)q];
                    if (vdone[w] != 1) continue;
                    fuv[w * 2] += shift_u; fuv[w * 2 + 1] += shift_v;
                    vdone[w] = 2;
                }
            }
            cursor += hi2[0] - lo2[0] + CLOSEUP_PIECE_GAP_PX;
            uvmin[0] = fmin(uvmin[0], lo2[0] + shift_u); uvmax[0] = fmax(uvmax[0], hi2[0] + shift_u);
            uvmin[1] = fmin(uvmin[1], 0.0); uvmax[1] = fmax(uvmax[1], hi2[1] - lo2[1]);
            fprintf(stderr, "  frame of piece %zu: layout handedness %s (agree %.3f), reflect %d, rotation %.2f deg, "
                    "at u %.0f..%.0f\n", rg, fr.hand_vote > 0.0 ? "team" : "MIRROR", fr.hand_agree, fr.reflect,
                    fr.theta * 180.0 / 3.14159265358979, lo2[0] + shift_u, hi2[0] + shift_u);
            frame_len += (size_t)snprintf(frame_json + frame_len, sizeof frame_json - frame_len,
                "%s{\"region\": %zu, \"reflect\": %d, \"theta_deg\": %.4f, \"hand_agree\": %.5f, \"shift\": [%.3f, %.3f]}",
                a ? ", " : "", rg, fr.reflect, fr.theta * 180.0 / 3.14159265358979, fr.hand_agree, shift_u, shift_v);
        }
        free(fuse); free(vdone); free(order); free(mu);
    }

    /* observations and the fused surface */
    CloseupObsGate_default(&gate);
    if (CloseupObs_add_faces(&obs, sverts, fuv, sfaces, snf, NULL, face_chart, &gate, &used, &gated) != 0) {
        fprintf(stderr, "ERROR: out of memory for observations\n");
        goto out;
    }
    fprintf(stderr, "  observations: %zu from %zu faces (%zu gated) (%.1fs)\n", obs.n, used, gated, ves_clock_sec() - t0);
    {   /* room for adopted papyrus around the seed pieces */
        double mg = fit.fill_px + 4.0 + (c->npool > 0 ? CLOSEUP_ADOPT_MARGIN_PX : 0.0);
        if (CloseupGrid_alloc(&g, uvmin[0] - mg, uvmax[0] + mg, uvmin[1] - mg, uvmax[1] + mg, CLOSEUP_GRID_PX) != 0) {
            fprintf(stderr, "ERROR: cannot allocate the node grid\n");
            goto out;
        }
    }
    if (CloseupSurface_fit(&g, &obs, &fit, &fs) != 0) { fprintf(stderr, "ERROR: the fused fit failed\n"); goto out; }
    fprintf(stderr, "  fused: %zu x %zu nodes: %zu data (%zu failed), %zu fill (%zu rejected); observations %zu "
            "(%zu downweighted), |obs - surface| p50 %.3f p90 %.3f p99 %.3f vox; torn %zu stretched + %zu floating; "
            "blobs %zu (%zu with a wrap), bridges %zu (%.1fs)\n", g.ni, g.nj,
            fs.data_nodes, fs.data_failed, fs.fill_nodes, fs.fill_rejected, fs.obs, fs.obs_downweighted,
            fs.resid_p50, fs.resid_p90, fs.resid_p99, fs.torn_stretch, fs.torn_float, fs.blobs, fs.blobs_with_wrap, fs.bridges, fs.seconds);
    ninward = (float *)malloc(g.ni * g.nj * 3 * sizeof *ninward);
    if (ninward == NULL) goto out;
    ic_axis_fields(axis, g.pos, g.ni * g.nj, NULL, ninward);
    orient = CloseupSurface_orient(&g, ninward, &agree);
    fprintf(stderr, "  normals: inward sign %d (agree %.4f)%s\n", orient, agree,
            orient < 0 ? " -- the fused frame is MIRROR-imaged" : "");

    /* the reduced placement: pool charts that continue the fused sheet join it by their 3-D contact */
    seedchart = (uint8_t *)calloc(nchart_ids + 1, 1);
    adopted = (uint8_t *)calloc(nchart_ids + 1, 1);
    if (seedchart == NULL || adopted == NULL) goto out;
    for (f = 0; f < snf; f++) if (face_chart[f] >= 0 && (size_t)face_chart[f] < nchart_ids) seedchart[face_chart[f]] = 1;
    adopt_len += (size_t)snprintf(adopt_json + adopt_len, sizeof adopt_json - adopt_len, "[");
    if (adopt_len >= sizeof adopt_json) adopt_len = sizeof adopt_json - 1;   /* never past the buffer */
    {
        int pi = 0;
        for (pi = 0; pi < c->npool; pi++) {
            IcPool pool;
            uint8_t *done = NULL;
            double lo3[3] = { 1e300, 1e300, 1e300 }, hi3[3] = { -1e300, -1e300, -1e300 }, *pfuv = NULL;
            int round = 0, k3 = 0;
            for (f = 0; f < g.ni * g.nj; f++) {
                if (g.state[f] == CLOSEUP_NODE_NONE) continue;
                for (k3 = 0; k3 < 3; k3++) {
                    if (g.pos[f * 3 + (size_t)k3] < lo3[k3]) lo3[k3] = g.pos[f * 3 + (size_t)k3];
                    if (g.pos[f * 3 + (size_t)k3] > hi3[k3]) hi3[k3] = g.pos[f * 3 + (size_t)k3];
                }
            }
            for (k3 = 0; k3 < 3; k3++) { lo3[k3] -= c->search_margin; hi3[k3] += c->search_margin; }
            {
                IcNear nreg;
                const IcNear *nuse = NULL;
                int lr = 0;
                if (ic_near_for_pool(&g, lo3, hi3, &nreg, &nuse) != 0) goto out;
                lr = ic_pool_load(c->pool[pi].path, c->pool[pi].surface, axis, lo3, hi3, nuse,
                                  c->pool[pi].surface || c->candidate_manifest ? NULL : seedchart,
                                  c->pool[pi].surface || c->candidate_manifest ? 0 : nchart_ids, 0, &pool);
                ic_near_free(&nreg);
                if (lr != 0) goto out;
            }
            fprintf(stderr, "  pool %s: %zu charts, %zu faces, %zu verts within %.0f vox of the surface (%.1fs)\n",
                    c->pool[pi].path, pool.p.nch, pool.p.nf, pool.p.nv, c->search_margin, ves_clock_sec() - t0);
            done = (uint8_t *)calloc(pool.p.nch + 1, 1);
            pfuv = (double *)malloc((pool.p.nv + 1) * 2 * sizeof *pfuv);
            if (done == NULL || pfuv == NULL) { free(done); free(pfuv); ic_pool_free(&pool); goto out; }
            for (round = 0; round < c->adopt_rounds; round++) {
                CloseupPlacement *pl = NULL;
                CloseupAdoptStats ast;
                size_t npl = 0, k = 0, obs0 = obs.n, area_faces = 0;
                if (ic_adopt_round(c, &g, &pool, done, pi, round, &pl, &npl, &ast, NULL) != 0) {
                    free(done); free(pfuv); ic_pool_free(&pool); goto out;
                }
                for (k = 0; k < npl; k++) {
                    size_t ch = pl[k].chart, a = 0, u2 = 0, g2 = 0;
                    done[ch] = 1;
                    for (a = pool.vofs[ch]; a < pool.vofs[ch + 1]; a++)
                        CloseupPlacement_apply(&pl[k], pool.uv[a * 2], pool.uv[a * 2 + 1], &pfuv[a * 2], &pfuv[a * 2 + 1]);
                    if (CloseupObs_add_faces(&obs, pool.verts, pfuv, pool.faces + pool.fofs[ch] * 3,
                                             pool.fofs[ch + 1] - pool.fofs[ch], NULL, pool.face_chart + pool.fofs[ch],
                                             &gate, &u2, &g2) != 0) {
                        free(pl); free(done); free(pfuv); ic_pool_free(&pool); goto out;
                    }
                    area_faces += pool.fofs[ch + 1] - pool.fofs[ch];
                    if (pool.chart_id[ch] >= 0 && (size_t)pool.chart_id[ch] < nchart_ids) adopted[pool.chart_id[ch]] = 1;
                }
                fprintf(stderr, "  adopt round %d: %zu candidates, %zu touch the sheet, %zu placed (+%zu observations from "
                        "%zu faces); refused: extent %zu, rms %zu, conflict %zu, nothing new %zu, metric %zu, "
                        "continuity %zu, layered %zu; deferred %zu (layered %zu)\n",
                        round, ast.candidates, ast.touched, npl, obs.n - obs0, area_faces, ast.refused_extent,
                        ast.refused_rms, ast.refused_conflict, ast.refused_nothing_new, ast.refused_metric,
                        ast.refused_continuity, ast.refused_layer, ast.deferred, ast.deferred_layer);
                adopt_len += (size_t)snprintf(adopt_json + adopt_len, sizeof adopt_json - adopt_len,
                    "%s\n    {\"pool\": %d, \"round\": %d, \"candidates\": %zu, \"touched\": %zu, \"placed\": %zu, "
                    "\"observations\": %zu, \"refused_extent\": %zu, \"refused_rms\": %zu, \"refused_conflict\": %zu, "
                    "\"refused_nothing_new\": %zu, \"refused_metric\": %zu, \"refused_continuity\": %zu, "
                    "\"refused_layer\": %zu, \"deferred\": %zu, \"deferred_layer\": %zu}",
                    (nadopt_total || round || pi) ? "," : "", pi, round,
                    ast.candidates, ast.touched, npl, obs.n - obs0, ast.refused_extent, ast.refused_rms,
                    ast.refused_conflict, ast.refused_nothing_new, ast.refused_metric, ast.refused_continuity,
                    ast.refused_layer, ast.deferred, ast.deferred_layer);
                if (adopt_len >= sizeof adopt_json) adopt_len = sizeof adopt_json - 1;   /* never past the buffer */
                nadopt_total += npl;
                adopt_obs += obs.n - obs0;
                free(pl);
                if (npl == 0) break;
                if (CloseupSurface_fit(&g, &obs, &fit, &fs) != 0) {
                    free(done); free(pfuv); ic_pool_free(&pool);
                    fprintf(stderr, "ERROR: refit failed\n");
                    goto out;
                }
                ic_axis_fields(axis, g.pos, g.ni * g.nj, NULL, ninward);
                orient = CloseupSurface_orient(&g, ninward, &agree);
                fprintf(stderr, "    refit: %zu data, %zu fill; |obs - surface| p50 %.3f p90 %.3f p99 %.3f vox; torn %zu "
                        "stretched + %zu floating; inward %d (agree %.4f) (%.1fs)\n", fs.data_nodes, fs.fill_nodes,
                        fs.resid_p50, fs.resid_p90, fs.resid_p99, fs.torn_stretch, fs.torn_float, orient, agree,
                        ves_clock_sec() - t0);
            }
            free(done); free(pfuv);
            ic_pool_free(&pool);
            fprintf(stderr, "  pool released (%.1fs)\n", ves_clock_sec() - t0);
        }
        for (f = 0; f < nchart_ids; f++) nadopted_ids += adopted[f];
    }
    adopt_len += (size_t)snprintf(adopt_json + adopt_len, sizeof adopt_json - adopt_len, "\n  ]");
    if (adopt_len >= sizeof adopt_json) adopt_len = sizeof adopt_json - 1;   /* never past the buffer */
    {   /* the frame's node-state map: data grey by chart, fill blue, none black */
        uint8_t *rgb = (uint8_t *)calloc(g.ni * g.nj * 3, 1);
        if (rgb != NULL) {
            for (f = 0; f < g.ni * g.nj; f++) {
                if (g.state[f] == CLOSEUP_NODE_DATA) {
                    uint32_t hsh = ic_hash((uint32_t)g.chart[f] + 1u);
                    int adp = g.chart[f] >= 0 && (size_t)g.chart[f] < nchart_ids && adopted[g.chart[f]];
                    rgb[f * 3] = (uint8_t)(adp ? 200 + (hsh & 55) : 90 + (hsh & 127));
                    rgb[f * 3 + 1] = (uint8_t)(adp ? 90 + ((hsh >> 8) & 60) : 90 + ((hsh >> 8) & 127));
                    rgb[f * 3 + 2] = (uint8_t)(adp ? 20 : 90 + ((hsh >> 16) & 127));
                } else if (g.state[f] == CLOSEUP_NODE_FILL) {
                    rgb[f * 3] = 30; rgb[f * 3 + 1] = 60; rgb[f * 3 + 2] = 200;
                }
            }
            snprintf(path, sizeof path, "%s/%s_frame_state.png", c->out_dir, c->id);
            ves_ensure_parent_dir(path);
            VesPng_write_rgb(path, rgb, (int)g.ni, (int)g.nj);
            fprintf(stderr, "  frame state map written (%.1fs)\n", ves_clock_sec() - t0);
            free(rgb);
        }
    }

    /* the isometric flattening, rotated so the axial coordinate runs down the rows */
    uvn = (double *)malloc(g.ni * g.nj * 2 * sizeof *uvn);
    if (uvn == NULL) goto out;
    fprintf(stderr, "  flatten start (%.1fs)\n", ves_clock_sec() - t0);
    if (CloseupSurface_flatten(&g, CLOSEUP_FLAT_PX, CLOSEUP_FLAT_ITERS, CLOSEUP_FLAT_TOL, uvn, &flat) != 0) {
        fprintf(stderr, "ERROR: flattening failed\n");
        goto out;
    }
    fprintf(stderr, "  flattened: %zu coarse nodes, %zu tris, %zu components, %zu iterations, energy %.4g -> %.4g; "
            "|log sigma| p50 %.4f p99 %.4f -> p50 %.4f p99 %.4f; area ratio %.4f (%.1fs)\n", flat.coarse_nodes,
            flat.coarse_tris, flat.components, flat.iterations, flat.energy_first, flat.energy_last,
            flat.logs_p50_before, flat.logs_p99_before, flat.logs_p50_after, flat.logs_p99_after, flat.area_ratio,
            flat.seconds);
    if (CloseupSurface_mesh(&g, uvn, &mverts, &muv, &mnrm, &mfaces, &node_of, &mnv, &mnf) != 0 || mnf == 0) {
        fprintf(stderr, "ERROR: the close-up surface has no faces\n");
        goto out;
    }
    maxial = (double *)malloc(mnv * sizeof *maxial);
    {
        float *minward = (float *)malloc(mnv * 3 * sizeof *minward);
        if (maxial == NULL || minward == NULL) { free(minward); goto out; }
        ic_axis_fields(axis, mverts, mnv, maxial, minward);
        ic_axial_sense(c, maxial, mnv);
        if (CloseupSurface_frame(mverts, muv, mfaces, mnf, NULL, maxial, minward, 0.0, 0.0, &fr2) != 0) {
            free(minward);
            fprintf(stderr, "ERROR: no frame on the flattened surface\n");
            goto out;
        }
        free(minward);
    }
    if (fr2.reflect) fprintf(stderr, "WARNING: the flattened surface is mirror-imaged (vote %.3g)\n", fr2.hand_vote);
    uvmin[0] = uvmin[1] = 1e300;
    for (v = 0; v < mnv; v++) {
        double a = 0.0, b = 0.0;
        CloseupFrame_apply(&fr2, muv[v * 2], muv[v * 2 + 1], &a, &b);
        muv[v * 2] = a; muv[v * 2 + 1] = b;
        if (a < uvmin[0]) uvmin[0] = a;
        if (b < uvmin[1]) uvmin[1] = b;
    }
    for (v = 0; v < mnv; v++) {
        muv[v * 2] += (double)CLOSEUP_RASTER_MARGIN_PX - floor(uvmin[0]);
        muv[v * 2 + 1] += (double)CLOSEUP_RASTER_MARGIN_PX - floor(uvmin[1]);
    }
    vattr = (uint8_t *)malloc(mnv);
    vtx_chart = (int32_t *)malloc(mnv * sizeof *vtx_chart);
    if (vattr == NULL || vtx_chart == NULL) goto out;
    for (v = 0; v < mnv; v++) {
        size_t n = (size_t)node_of[v];
        vattr[v] = g.state[n];
        vtx_chart[v] = g.chart[n];
        if (g.state[n] == CLOSEUP_NODE_DATA && g.chart[n] >= 0 && (size_t)g.chart[n] < nchart_ids && adopted[g.chart[n]])
            vattr[v] = 3;
    }
    snprintf(path, sizeof path, "%s/%s_surface.vmesh", c->out_dir, c->id);
    fprintf(stderr, "  surface mesh ready (%.1fs)\n", ves_clock_sec() - t0);
    if (MeshBin_write_uv64(path, mverts, mnv, mfaces, mnf, muv) != 0) { fprintf(stderr, "ERROR: cannot write %s\n", path); goto out; }
    snprintf(path, sizeof path, "%s/%s_surface_nrm.f32", c->out_dir, c->id);
    if (ic_write_file(path, mnrm, mnv * 3 * sizeof *mnrm) != 0) goto out;
    snprintf(path, sizeof path, "%s/%s_surface_attr.u8", c->out_dir, c->id);
    if (ic_write_file(path, vattr, mnv) != 0) goto out;
    snprintf(path, sizeof path, "%s/%s_surface_chart.i32", c->out_dir, c->id);
    if (ic_write_file(path, vtx_chart, mnv * sizeof *vtx_chart) != 0) goto out;
    snprintf(path, sizeof path, "%s/%s_surface_report.json", c->out_dir, c->id);
    if ((fp = fopen(path, "wb")) == NULL) goto out;
    ic_json_str(js, sizeof js, c->layout);
    fprintf(fp, "{\n  \"schema\": \"ink-closeup-surface-v1\",\n  \"id\": \"%s\",\n  \"layout\": %s,\n", c->id, js);
    fprintf(fp, "  \"seed_box\": [%.1f, %.1f, %.1f, %.1f],\n  \"grow_px\": %.1f,\n  \"seed_regions\": [", c->seed[0],
            c->seed[1], c->seed[2], c->seed[3], c->grow_px);
    {
        int first = 1;
        for (f = 0; f <= (size_t)max_region; f++)
            if (seed_region[f]) { fprintf(fp, "%s%zu", first ? "" : ", ", f); first = 0; }
    }
    fprintf(fp, "],\n  \"source\": {\"faces\": %zu, \"verts\": %zu, \"observations\": %zu, \"faces_gated\": %zu},\n",
            snf, snv, obs.n, gated);
    fprintf(fp, "  \"frames\": [%s],\n  \"z_top_to_bottom\": %s,\n", frame_json, c->z_top_to_bottom ? "true" : "false");
    fprintf(fp, "  \"complete_gaps\": %s,\n  \"fill_policy\": {\"reach_px\": %.1f, \"fit_max_px\": %.1f, "
            "\"cartesian_thin_plate\": %s, \"holes\": %zu, \"filled\": %zu, \"refused\": %zu, \"nodes\": %zu},\n",
            c->complete_gaps ? "true" : "false", fit.fill_px, fit.fill_fit_max_px,
            fit.variational_holes ? "true" : "false", fs.variational_holes, fs.variational_filled,
            fs.variational_refused, fs.variational_nodes);
    fprintf(fp, "  \"fit\": {\"grid_px\": %.1f, \"nodes\": [%zu, %zu], \"data\": %zu, \"data_failed\": %zu, "
            "\"fill\": %zu, \"fill_rejected\": %zu, \"obs_downweighted\": %zu, \"resid_p50\": %.4f, "
            "\"resid_p90\": %.4f, \"resid_p99\": %.4f, \"torn_stretch\": %zu, \"torn_float\": %zu, "
            "\"blobs\": %zu, \"blobs_with_wrap\": %zu, \"bridges\": %zu, \"bridge\": %s, "
            "\"orient\": %d, \"orient_agree\": %.5f, \"seconds\": %.2f},\n",
            CLOSEUP_GRID_PX, g.ni, g.nj, fs.data_nodes, fs.data_failed, fs.fill_nodes, fs.fill_rejected,
            fs.obs_downweighted, fs.resid_p50, fs.resid_p90, fs.resid_p99, fs.torn_stretch, fs.torn_float,
            fs.blobs, fs.blobs_with_wrap, fs.bridges, fit.bridge ? "true" : "false", orient,
            agree, fs.seconds);
    fprintf(fp, "  \"flatten\": {\"method\": \"slim_scaf_fringe_v2\", \"coarse_px\": %.1f, \"coarse_nodes\": %zu, \"tris\": %zu, \"components\": %zu, "
            "\"iterations\": %zu, \"energy_first\": %.6g, \"energy_last\": %.6g, \"logs_p50_before\": %.5f, "
            "\"logs_p99_before\": %.5f, \"logs_p50_after\": %.5f, \"logs_p99_after\": %.5f, \"area_ratio\": %.5f, "
            "\"final_rotation_deg\": %.4f, \"seconds\": %.2f},\n", CLOSEUP_FLAT_PX, flat.coarse_nodes, flat.coarse_tris,
            flat.components, flat.iterations, flat.energy_first, flat.energy_last, flat.logs_p50_before,
            flat.logs_p99_before, flat.logs_p50_after, flat.logs_p99_after, flat.area_ratio,
            fr2.theta * 180.0 / 3.14159265358979, flat.seconds);
    {
        size_t q = 0;
        fprintf(fp, "  \"flatten_energy\": [");
        for (q = 0; q < flat.iterations && q < 64; q++) fprintf(fp, "%s%.6g", q ? ", " : "", flat.energy_hist[q]);
        fprintf(fp, "],\n");
    }
    fprintf(fp, "  \"adopt\": {\"charts\": %zu, \"chart_ids\": %zu, \"observations\": %zu, \"rounds\": %s},\n",
            nadopt_total, nadopted_ids, adopt_obs, adopt_json);
    fprintf(fp, "  \"surface\": {\"verts\": %zu, \"faces\": %zu},\n  \"seconds\": %.2f\n}\n", mnv, mnf,
            ves_clock_sec() - t0);
    fclose(fp); fp = NULL;
    fprintf(stderr, "[ink_closeup] SUMMARY surface id=%s verts=%zu faces=%zu data=%zu fill=%zu logs_p99=%.4f->%.4f (%.1fs)\n",
            c->id, mnv, mnf, fs.data_nodes, fs.fill_nodes, flat.logs_p99_before, flat.logs_p99_after,
            ves_clock_sec() - t0);
    rc = 0;
out:
    if (fp != NULL) fclose(fp);
    MeshBin_dispose(&m);
    CloseupObs_free(&obs);
    CloseupGrid_free(&g);
    free(vchart); free(region); free(face_chart); free(vmap); free(sfaces); free(node_of); free(mfaces);
    free(vtx_chart); free(seed_region); free(keep); free(vattr); free(sverts); free(sinward); free(ninward);
    free(mverts); free(mnrm); free(suv); free(fuv); free(saxial); free(uvn); free(muv); free(maxial); free(reg_count);
    free(seedchart); free(adopted);
    Arena_dispose(&arena);
    return rc;
}

/* PROFILE (closeup_complete_speed_20261030): one stderr line per phase with its wall time. */
static void ic_lap(const char *what, double *t)
{
    double now = ves_clock_sec();
    fprintf(stderr, "  [complete-timing] %-28s %8.2fs\n", what, now - *t);
    *t = now;
}

/* A completion pass's surface as its files hold it: the next pass of a chain reads it here instead of from disk. */
typedef struct IcSurfaceMem {
    MeshBinData m;                 /* verts, uv64, faces */
    uint8_t *attr;
    int32_t *chart;
} IcSurfaceMem;

static void ic_surface_mem_free(IcSurfaceMem *s)
{
    MeshBin_dispose(&s->m); free(s->attr); free(s->chart);
    memset(s, 0, sizeof *s);
}

/* One pass's files, written off the computing thread: the mesh and its per-vertex files, then the report (last, so a
 * report on disk means a complete pass).  The job owns its arrays. */
#define IC_WRITE_EXTRA 64
typedef struct IcWriteJob {
    IcFileBytes extra[IC_WRITE_EXTRA];   /* the pass's audit CSVs */
    int nextra;
    char stem[2700];               /* <out_dir>/<id> */
    float *verts, *nrm;
    double *muv;
    int32_t *faces, *oc;
    uint8_t *oa;
    size_t nv, nf;
    char *report;                  /* the report JSON up to its "seconds" */
    double t0;                     /* the pass's start: "seconds" run to the report's write */
    int rc;
    int no_mesh;                   /* a chain's intermediate pass: its CSVs and report, not its mesh */
} IcWriteJob;

typedef struct IcWriter {
    /* the chain's plan for the pass about to queue its files: keep every mesh, or only the final pass's -- the last
     * pass, or one gaining less than min_gain of the previous vertex count (the chain's own stop rule) */
    int keep_all, last;
    double previous, min_gain;
    int queued_no_mesh;            /* the last pass queued without its mesh */
#ifdef _WIN32
    uintptr_t thread;
#else
    pthread_t thread;
#endif
    int running, failed;
    IcWriteJob *job;
} IcWriter;

static void ic_write_job_free(IcWriteJob *j)
{
    int k = 0;
    if (j == NULL) return;
    for (k = 0; k < j->nextra; k++) free(j->extra[k].data);
    free(j->verts); free(j->nrm); free(j->muv); free(j->faces); free(j->oc); free(j->oa); free(j->report);
    free(j);
}

static int ic_write_job_run(IcWriteJob *j)
{
    char path[2800];
    FILE *fp = NULL;
    int k = 0;
    for (k = 0; k < j->nextra; k++) {
        ves_ensure_parent_dir(j->extra[k].path);
        if (ic_write_file(j->extra[k].path, j->extra[k].data, j->extra[k].n) != 0) return -1;
    }
    if (!j->no_mesh) {
        snprintf(path, sizeof path, "%s_surface.vmesh", j->stem);
        ves_ensure_parent_dir(path);
        if (MeshBin_write_uv64(path, j->verts, j->nv, j->faces, j->nf, j->muv) != 0) return -1;
        snprintf(path, sizeof path, "%s_surface_nrm.f32", j->stem);
        if (ic_write_file(path, j->nrm, j->nv * 3 * sizeof *j->nrm) != 0) return -1;
        snprintf(path, sizeof path, "%s_surface_attr.u8", j->stem);
        if (ic_write_file(path, j->oa, j->nv) != 0) return -1;
        snprintf(path, sizeof path, "%s_surface_chart.i32", j->stem);
        if (ic_write_file(path, j->oc, j->nv * sizeof *j->oc) != 0) return -1;
    }
    snprintf(path, sizeof path, "%s_surface_report.json", j->stem);
    fp = fopen(path, "wb");
    if (fp == NULL) return -1;
    if (fputs(j->report, fp) < 0 || fprintf(fp, "  \"seconds\": %.2f\n}\n", ves_clock_sec() - j->t0) < 0) {
        fclose(fp);
        return -1;
    }
    return fclose(fp) == 0 ? 0 : -1;
}

#ifdef _WIN32
static unsigned __stdcall ic_writer_main(void *arg)
#else
static void *ic_writer_main(void *arg)
#endif
{
    IcWriteJob *j = (IcWriteJob *)arg;
    j->rc = ic_write_job_run(j);
    if (j->rc != 0) fprintf(stderr, "ERROR: writing %s_surface.* failed\n", j->stem);
#ifdef _WIN32
    return 0;
#else
    return NULL;
#endif
}

/* The queued pass's files are on disk (or failed: -1, and every later wait fails too). */
static int ic_writer_wait(IcWriter *w)
{
    double t0 = ves_clock_sec();
    if (w->running) {
#ifdef _WIN32
        WaitForSingleObject((HANDLE)w->thread, INFINITE);
        CloseHandle((HANDLE)w->thread);
#else
        pthread_join(w->thread, NULL);
#endif
        w->running = 0;
        if (ves_clock_sec() - t0 > 0.05)
            fprintf(stderr, "  [complete-timing] waited %.2fs for the previous pass's files\n", ves_clock_sec() - t0);
    }
    if (w->job != NULL) {
        if (w->job->rc != 0) w->failed = 1;
        ic_write_job_free(w->job);
        w->job = NULL;
    }
    return w->failed ? -1 : 0;
}

/* Queue a pass's files (the writer takes the job); written inline when no thread starts. */
static int ic_writer_start(IcWriter *w, IcWriteJob *j)
{
    if (ic_writer_wait(w) != 0) { ic_write_job_free(j); return -1; }
    w->job = j;
#ifdef _WIN32
    w->thread = _beginthreadex(NULL, 0, ic_writer_main, j, 0, NULL);
    w->running = w->thread != 0;
#else
    w->running = pthread_create(&w->thread, NULL, ic_writer_main, j) == 0;
#endif
    if (!w->running) (void)ic_writer_main(j);
    return 0;
}

/* Complete an already flattened close-up when the original assembly is unavailable. Only faces
 * whose three corners were observed support the fit; old FILL never becomes a new observation.
 * Keep the existing UV frame plus an integer translation, so pixels compare without a glyph warp.
 * In a chain (ic_stage_complete_chain): src (taken) is the previous pass's surface in memory, out receives this one,
 * pc keeps the pool meshes between passes and wr writes the files behind the next pass; NULL each: the files. */
static int ic_complete_pass(const IcConfig *c, IcSurfaceMem *src, IcSurfaceMem *out, IcPoolCache *pc, IcWriter *wr)
{
    char path[2600], js[4800];
    MeshBinData m;
    CloseupObs obs;
    CloseupObsGate gate;
    CloseupGrid g;
    CloseupFitOpts fit;
    CloseupFitStats fs;
    Arena_T arena = NULL;
    AsmAxis *axis = NULL;
    CloseupFitCache *fitcache = NULL;   /* the adopt rounds' refits redo only what their charts reach */
    uint8_t *attr = NULL, *use = NULL, *oa = NULL;
    int32_t *chart = NULL, *fc = NULL, *faces = NULL, *nodes = NULL, *oc = NULL;
    float *inward = NULL, *verts = NULL, *nrm = NULL;
    double *uv = NULL, *uvn = NULL, *muv = NULL;
    double lo[2] = {1e300, 1e300}, hi[2] = {-1e300, -1e300}, shift = 0.0, agree = 0.0;
    double t0 = ves_clock_sec();
    size_t got = 0, f = 0, v = 0, used = 0, gated = 0, nv = 0, nf = 0, nn = 0;
    size_t nadopt = 0, obs_source = 0;
    int32_t adopt_first = 0, chart_offset = 0;
    int k = 0, orient = 0, rc = 1, ncsv = 0, q0 = 0;
    FILE *fp = NULL;
    char *report = NULL;
    IcFileBytes csv[IC_WRITE_EXTRA];
    AsmAxisProjector *axis_index = NULL;
    float *inward_at = NULL;           /* the positions inward was last measured at (ic_axis_inward_moved) */
    double tl = t0;
    memset(&m, 0, sizeof m); memset(&obs, 0, sizeof obs); memset(&g, 0, sizeof g);
    for (q0 = 0; q0 < IC_WRITE_EXTRA; q0++) { csv[q0].data = NULL; csv[q0].n = 0; csv[q0].path[0] = 0; }
    if (out != NULL) memset(out, 0, sizeof *out);
    if (c->complete_gaps) CloseupFitOpts_complete(&fit); else CloseupFitOpts_default(&fit);
    fit.bridge = c->bridge && c->chart_wrap != NULL; fit.chart_wrap = c->chart_wrap; fit.nchart_wrap = c->nchart_wrap;
    fit.bridge_px = fit.fill_px;    /* the bridge reaches as far as the fill in effect (24 px, or the completed-gap 160 px) */
    fprintf(stderr, "  [complete-timing] bridge=%d threads=%d\n", fit.bridge,
#ifdef _OPENMP
            omp_get_max_threads()
#else
            1
#endif
            );
    if (src != NULL) {                 /* the previous pass's files, as they would read back */
        m = src->m; attr = src->attr; chart = src->chart;
        memset(src, 0, sizeof *src);
        fprintf(stderr, "[ink_closeup] completion: the source %s from the chain's memory\n", c->surface_source);
    } else {
        snprintf(path, sizeof path, "%s_surface.vmesh", c->surface_source);
        if (MeshBin_read_precise_malloc(path, &m) != 0 || (m.uv == NULL && m.uv64 == NULL)) goto out;
        snprintf(path, sizeof path, "%s_surface_attr.u8", c->surface_source);
        attr = (uint8_t *)ic_read_file(path, m.nv, &got);
        snprintf(path, sizeof path, "%s_surface_chart.i32", c->surface_source);
        chart = (int32_t *)ic_read_file(path, m.nv * sizeof *chart, &got);
    }
    ic_lap("read source", &tl);
    uv = (double *)malloc(m.nv * 2 * sizeof *uv);
    use = (uint8_t *)calloc(m.nf, 1);
    fc = (int32_t *)malloc(m.nf * sizeof *fc);
    if (attr == NULL || chart == NULL || uv == NULL || use == NULL || fc == NULL) goto out;
    for (v = 0; v < m.nv; v++)
    {
        if (chart[v] >= adopt_first) {
            if (chart[v] == INT32_MAX) goto out;
            adopt_first = chart[v] + 1;
        }
        for (k = 0; k < 2; k++) {
            double a = MeshBin_uv(&m, v * 2 + (size_t)k);
            if (!isfinite(a)) goto out;
            uv[v * 2 + (size_t)k] = a;
            if (attr[v] == 1 || attr[v] == 3) { lo[k] = fmin(lo[k], a); hi[k] = fmax(hi[k], a); }
        }
    }
    for (f = 0; f < m.nf; f++) {
        use[f] = 1;
        fc[f] = chart[m.faces[f * 3]];
        for (k = 0; k < 3; k++) {
            uint8_t a = attr[m.faces[f * 3 + (size_t)k]];
            if (a != 1 && a != 3) use[f] = 0;
        }
    }
    CloseupObsGate_default(&gate);
    if (CloseupObs_add_faces(&obs, m.verts, uv, m.faces, m.nf, use, fc, &gate, &used, &gated) != 0 || obs.n == 0) goto out;
    obs_source = obs.n;
    fprintf(stderr, "[ink_closeup] completion: %zu observations from %zu observed faces (%zu gated)\n", obs.n, used, gated);
    ic_lap("source observations", &tl);
    {
        double margin = fit.fill_px + 4.0 + (c->npool ? CLOSEUP_ADOPT_MARGIN_PX : 0.0);
        if (c->has_growth) {
            lo[0] = fmin(lo[0], c->growth_uv[0]); lo[1] = fmin(lo[1], c->growth_uv[1]);
            hi[0] = fmax(hi[0], c->growth_uv[2]); hi[1] = fmax(hi[1], c->growth_uv[3]);
            margin = 0;
        }
        if (CloseupGrid_alloc(&g, lo[0] - margin, hi[0] + margin, lo[1] - margin, hi[1] + margin, CLOSEUP_GRID_PX) != 0) goto out;
    }
    fprintf(stderr, "  [complete-timing] grid %zu x %zu = %zu nodes\n", g.ni, g.nj, g.ni * g.nj);
    ic_lap("grid alloc", &tl);
    fitcache = CloseupFitCache_new();
    if (fitcache == NULL || CloseupSurface_fit_cached(&g, &obs, &fit, &fs, fitcache) != 0) goto out;
    ic_lap("fit 0", &tl);
    nn = g.ni * g.nj;
    inward = (float *)malloc(nn * 3 * sizeof *inward);
    uvn = (double *)malloc(nn * 2 * sizeof *uvn);
    if (inward == NULL || uvn == NULL) goto out;
    arena = Arena_new(); axis = AsmAxis_load(arena, c->axis);
    if (axis == NULL) goto out;
    axis_index = AsmAxis_projector(arena, axis);
    inward_at = (float *)malloc(nn * 3 * sizeof *inward_at);
    if (inward_at == NULL) goto out;
    memset(inward_at, 0xff, nn * 3 * sizeof *inward_at);     /* a NaN no position carries: every node measured */
    ic_axis_inward_moved(axis, axis_index, g.pos, nn, inward_at, inward);
    orient = CloseupSurface_orient(&g, inward, &agree);
    ic_lap("axis + orient 0", &tl);
    chart_offset = adopt_first;
    for (k = 0; k < c->npool; k++) {
        IcPool pool;
        double lo3[3] = {1e300, 1e300, 1e300}, hi3[3] = {-1e300, -1e300, -1e300};
        double *puv = NULL;
        uint8_t *done = NULL;
        int q = 0, round = 0;
        for (v = 0; v < nn; v++) {
            if (g.state[v] == CLOSEUP_NODE_NONE) continue;
            for (q = 0; q < 3; q++) { lo3[q] = fmin(lo3[q], g.pos[v * 3 + (size_t)q]); hi3[q] = fmax(hi3[q], g.pos[v * 3 + (size_t)q]); }
        }
        for (q = 0; q < 3; q++) { lo3[q] -= c->search_margin; hi3[q] += c->search_margin; }
        fprintf(stderr, "  [complete-timing] search box %.0f..%.0f %.0f..%.0f %.0f..%.0f\n",
                lo3[0], hi3[0], lo3[1], hi3[1], lo3[2], hi3[2]);
        {
            IcNear nreg;
            const IcNear *nuse = NULL;
            int lr = 0;
            if (ic_near_for_pool(&g, lo3, hi3, &nreg, &nuse) != 0) goto out;
            ic_lap("near region", &tl);
            lr = pc != NULL && nuse == NULL
                     ? ic_pool_load_cached(pc, c->pool[k].path, c->pool[k].surface, axis, lo3, hi3, 1, &pool)
                     : ic_pool_load(c->pool[k].path, c->pool[k].surface, axis, lo3, hi3, nuse, NULL, 0, 1, &pool);
            ic_near_free(&nreg);
            if (lr != 0) goto out;
            ic_lap("pool load", &tl);
        }
        puv = (double *)malloc((pool.p.nv + 1) * 2 * sizeof *puv);
        done = (uint8_t *)calloc(pool.p.nch + 1, 1);
        if (puv == NULL || done == NULL) { free(puv); free(done); ic_pool_free(&pool); goto out; }
        /* Each source gets its own chart namespace; unrelated cached chart IDs may collide. */
        for (v = 0; v < pool.p.nf; v++) {
            if (pool.face_chart[v] > INT32_MAX - chart_offset) { free(puv); free(done); ic_pool_free(&pool); goto out; }
            pool.face_chart[v] += chart_offset;
        }
        for (v = 0; v < pool.p.nf; v++) if (pool.face_chart[v] >= chart_offset) {
            if (pool.face_chart[v] == INT32_MAX) { free(puv); free(done); ic_pool_free(&pool); goto out; }
            chart_offset = pool.face_chart[v] + 1;
        }
        fprintf(stderr, "  completion pool %d: %zu charts, %zu faces\n", k, pool.p.nch, pool.p.nf);
        for (round = 0; round < c->adopt_rounds; round++) {
            CloseupPlacement *pl = NULL;
            CloseupAdoptStats ast;
            size_t npl = 0, p = 0, before = obs.n;
            if (ic_adopt_round(c, &g, &pool, done, k, round, &pl, &npl, &ast,
                               wr != NULL && ncsv < IC_WRITE_EXTRA ? &csv[ncsv] : NULL) != 0) {
                free(puv); free(done); ic_pool_free(&pool); goto out;
            }
            if (wr != NULL && ncsv < IC_WRITE_EXTRA && csv[ncsv].data != NULL) ncsv++;
            ic_lap("adopt round", &tl);
            for (p = 0; p < npl; p++) {
                size_t ch = pl[p].chart, a = 0, nu = 0, ng = 0;
                done[ch] = 1;
                for (a = pool.vofs[ch]; a < pool.vofs[ch + 1]; a++)
                    CloseupPlacement_apply(&pl[p], pool.uv[a * 2], pool.uv[a * 2 + 1], &puv[a * 2], &puv[a * 2 + 1]);
                if (CloseupObs_add_faces(&obs, pool.verts, puv, pool.faces + pool.fofs[ch] * 3,
                                         pool.fofs[ch + 1] - pool.fofs[ch], NULL, pool.face_chart + pool.fofs[ch],
                                         &gate, &nu, &ng) != 0) {
                    free(pl); free(puv); free(done); ic_pool_free(&pool); goto out;
                }
            }
            fprintf(stderr, "  completion adopt %d/%d: %zu placed, +%zu observations; refused rms %zu, conflict %zu; deferred %zu\n",
                    k, round, npl, obs.n - before, ast.refused_rms, ast.refused_conflict, ast.deferred);
            ic_lap("adopted observations", &tl);
            free(pl); nadopt += npl;
            if (npl == 0) break;
            if (CloseupSurface_fit_cached(&g, &obs, &fit, &fs, fitcache) != 0) { free(puv); free(done); ic_pool_free(&pool); goto out; }
            ic_lap("fit", &tl);
            ic_axis_inward_moved(axis, axis_index, g.pos, nn, inward_at, inward);
            orient = CloseupSurface_orient(&g, inward, &agree);
            ic_lap("axis + orient", &tl);
        }
        free(puv); free(done);
        if (pc != NULL) ic_pool_keep(&pool, &pc->spare); else ic_pool_free(&pool);
        ic_lap("pool free", &tl);
    }
    if (CLOSEUP_APRON_TRIM) {   /* the outward apron: FILL only between observed papyrus (2026-10-15) */
        CloseupApronStats ap;
        memset(&ap, 0, sizeof ap);
        if (CloseupSurface_trim_apron(&g, fit.fill_px, &ap) != 0) goto out;
        fprintf(stderr, "[ink_closeup] completion: outward apron: %zu of %zu FILL nodes had papyrus on one side only (reach "
                "%.0f px): dropped (%.1fs)\n", ap.dropped, ap.fill_before, fit.fill_px, ap.seconds);
    }
    ic_lap("apron", &tl);
    {   /* duplicated papyrus: FILL extrapolated over papyrus the surface holds elsewhere goes (2026-10-15) */
        CloseupDedupStats dd;
        memset(&dd, 0, sizeof dd);
        if (CloseupSurface_dedup_fill(&g, CLOSEUP_DEDUP_VOX, CLOSEUP_DEDUP_FRAME_PX, CLOSEUP_DEDUP_COS, &dd) != 0) goto out;
        fprintf(stderr, "[ink_closeup] completion: duplicated papyrus: %zu node pairs (%zu DATA-DATA, kept); %zu FILL nodes "
                "dropped\n", dd.pairs, dd.data_pairs, dd.dropped_fill);
    }
    ic_lap("dedup", &tl);
    shift = ceil((CLOSEUP_RASTER_MARGIN_PX - fmin(g.u0, g.v0)) / g.h) * g.h;
    if (shift < 0.0) shift = 0.0;
    for (v = 0; v < nn; v++) {
        uvn[v * 2] = g.state[v] ? g.u0 + g.h * (double)(v % g.ni) + shift : NAN;
        uvn[v * 2 + 1] = g.state[v] ? g.v0 + g.h * (double)(v / g.ni) + shift : NAN;
    }
    if (CloseupSurface_mesh(&g, uvn, &verts, &muv, &nrm, &faces, &nodes, &nv, &nf) != 0 || nf == 0) goto out;
    oa = (uint8_t *)malloc(nv); oc = (int32_t *)malloc(nv * sizeof *oc);
    if (oa == NULL || oc == NULL) goto out;
    for (v = 0; v < nv; v++) {
        oa[v] = g.state[nodes[v]]; oc[v] = g.chart[nodes[v]];
        if (oa[v] == CLOSEUP_NODE_DATA && oc[v] >= adopt_first) oa[v] = 3;
    }
    ic_lap("mesh", &tl);
    {   /* the report, all but its "seconds" (ic_write_job_run adds them as it writes it, after the other files) */
        const size_t cap = 16384;
        int at = 0;
        report = (char *)malloc(cap);
        if (report == NULL) goto out;
        ic_json_str(js, sizeof js, c->surface_source);
        at += snprintf(report + at, cap - (size_t)at, "{\n  \"schema\": \"ink-closeup-completion-v1\",\n  \"surface_source\": %s,\n", js);
        at += snprintf(report + at, cap - (size_t)at, "  \"uv_translation\": [%.1f, %.1f],\n  \"frame\": \"source UV retained; integer translation only\",\n", shift, shift);
        at += snprintf(report + at, cap - (size_t)at, "  \"complete_gaps\": %s,\n  \"fill_policy\": {\"reach_px\": %.1f, \"fit_max_px\": %.1f, "
                "\"cartesian_thin_plate\": %s, \"holes\": %zu, \"filled\": %zu, \"refused\": %zu, \"nodes\": %zu},\n",
                c->complete_gaps ? "true" : "false", fit.fill_px, fit.fill_fit_max_px,
                fit.variational_holes ? "true" : "false", fs.variational_holes, fs.variational_filled,
                fs.variational_refused, fs.variational_nodes);
        at += snprintf(report + at, cap - (size_t)at, "  \"adopt\": {\"charts\": %zu, \"observations\": %zu, \"first_chart_id\": %d},\n",
                nadopt, obs.n - obs_source, adopt_first);
        at += snprintf(report + at, cap - (size_t)at, "  \"fit\": {\"observations\": %zu, \"observed_faces\": %zu, \"faces_gated\": %zu, \"data\": %zu, "
                "\"fill\": %zu, \"fill_rejected\": %zu, \"resid_p50\": %.6f, \"resid_p99\": %.6f, "
                "\"torn_stretch\": %zu, \"torn_float\": %zu, "
                "\"orient\": %d, \"orient_agree\": %.6f},\n", obs.n, used, gated, fs.data_nodes, fs.fill_nodes,
                fs.fill_rejected, fs.resid_p50, fs.resid_p99, fs.torn_stretch, fs.torn_float, orient, agree);
        at += snprintf(report + at, cap - (size_t)at, "  \"surface\": {\"verts\": %zu, \"faces\": %zu},\n", nv, nf);
        if (at <= 0 || (size_t)at >= cap) goto out;
    }
    {   /* the files: behind the next pass (a chain) or now */
        IcWriteJob *job = (IcWriteJob *)calloc(1, sizeof *job);
        if (job == NULL) goto out;
        snprintf(job->stem, sizeof job->stem, "%s/%s", c->out_dir, c->id);
        job->nv = nv; job->nf = nf; job->t0 = t0;
        for (q0 = 0; q0 < ncsv; q0++) { job->extra[q0] = csv[q0]; csv[q0].data = NULL; }
        job->nextra = ncsv; ncsv = 0;
        job->nrm = nrm; nrm = NULL;
        job->report = report; report = NULL;
        if (wr != NULL && out != NULL && !wr->keep_all && !wr->last &&
            !(wr->previous > 0.0 && (double)nv / wr->previous - 1.0 < wr->min_gain)) {
            job->no_mesh = 1;          /* an intermediate pass: the next pass takes the arrays, no file needs them */
            free(job->nrm); job->nrm = NULL;
        }
        if (wr != NULL) wr->queued_no_mesh = job->no_mesh;
        if (out != NULL && job->no_mesh) {
            out->m.verts = verts; out->m.uv64 = muv; out->m.faces = faces; out->m.nv = nv; out->m.nf = nf;
            out->attr = oa; out->chart = oc;
            verts = NULL; muv = NULL; faces = NULL; oa = NULL; oc = NULL;
        } else if (out != NULL) {      /* the next pass takes the arrays; the writer, copies */
            job->verts = (float *)malloc(nv * 3 * sizeof *verts); job->muv = (double *)malloc(nv * 2 * sizeof *muv);
            job->faces = (int32_t *)malloc(nf * 3 * sizeof *faces); job->oa = (uint8_t *)malloc(nv);
            job->oc = (int32_t *)malloc(nv * sizeof *oc);
            if (!job->verts || !job->muv || !job->faces || !job->oa || !job->oc) { ic_write_job_free(job); goto out; }
            memcpy(job->verts, verts, nv * 3 * sizeof *verts); memcpy(job->muv, muv, nv * 2 * sizeof *muv);
            memcpy(job->faces, faces, nf * 3 * sizeof *faces); memcpy(job->oa, oa, nv); memcpy(job->oc, oc, nv * sizeof *oc);
            out->m.verts = verts; out->m.uv64 = muv; out->m.faces = faces; out->m.nv = nv; out->m.nf = nf;
            out->attr = oa; out->chart = oc;
            verts = NULL; muv = NULL; faces = NULL; oa = NULL; oc = NULL;
        } else {
            job->verts = verts; job->muv = muv; job->faces = faces; job->oa = oa; job->oc = oc;
            verts = NULL; muv = NULL; faces = NULL; oa = NULL; oc = NULL;
        }
        if (wr != NULL) {
            if (ic_writer_start(wr, job) != 0) goto out;
        } else {
            int wrc = ic_write_job_run(job);
            ic_write_job_free(job);
            if (wrc != 0) goto out;
        }
    }
    ic_lap(wr != NULL ? "queue the files" : "write", &tl);
    fprintf(stderr, "[ink_closeup] completion: %zu data / %zu fill (%zu rejected); residual p50 %.4f p99 %.4f vox; "
            "torn %zu stretched + %zu floating; UV shift %.0f; normals %d agree %.4f; %zu verts %zu faces (%.1fs)\n",
            fs.data_nodes, fs.fill_nodes, fs.fill_rejected, fs.resid_p50, fs.resid_p99, fs.torn_stretch,
            fs.torn_float, shift, orient, agree, nv, nf, ves_clock_sec() - t0);
    rc = 0;
out:
    if (fp != NULL) fclose(fp);
    if (rc) fprintf(stderr, "ERROR: surface completion failed for %s\n", c->surface_source);
    if (rc && out != NULL) ic_surface_mem_free(out);
    for (q0 = 0; q0 < IC_WRITE_EXTRA; q0++) free(csv[q0].data);
    free(report); free(inward_at);
    MeshBin_dispose(&m); CloseupObs_free(&obs); CloseupGrid_free(&g); Arena_dispose(&arena);
    CloseupFitCache_free(&fitcache);
    free(attr); free(use); free(oa); free(chart); free(fc); free(faces); free(nodes); free(oc);
    free(inward); free(verts); free(nrm); free(uv); free(uvn); free(muv);
    return rc;
}

static int ic_stage_complete_surface(const IcConfig *c)
{
    return ic_complete_pass(c, NULL, NULL, NULL, NULL);
}

/* "1,234,567" (a pass's progress line, as closeup_fibers.py printed it) */
static void ic_thousands(char *buf, size_t cap, size_t value)
{
    char digits[32];
    size_t n = (size_t)snprintf(digits, sizeof digits, "%zu", value), i = 0, o = 0;
    for (i = 0; i < n && o + 2 < cap; i++) {
        if (i > 0 && (n - i) % 3 == 0) buf[o++] = ',';
        buf[o++] = digits[i];
    }
    buf[o] = 0;
}

/* Whether two stems name one file set (separators and, on Windows, case aside). */
static int ic_same_stem(const char *a, const char *b)
{
    for (; *a && *b; a++, b++) {
        int x = (unsigned char)*a, y = (unsigned char)*b;
        if (x == '\\') x = '/';
        if (y == '\\') y = '/';
#ifdef _WIN32
        if (x >= 'A' && x <= 'Z') x += 'a' - 'A';
        if (y >= 'A' && y <= 'Z') y += 'a' - 'A';
#endif
        if (x != y) return 0;
    }
    return *a == *b;
}

/* ---- the completion stage's own files (closeup_fibers.py wrote these before 2026-10-30): pass configs, summary */

#ifdef _WIN32
#define IC_SEP '\\'
#else
#define IC_SEP '/'
#endif

/* A JSON number the way Python's json writes one: the shortest text that reads back as the same double, with ".0" on
 * an integral value written as a float (as_float); null when not finite. */
static void ic_json_shortest(char *buf, size_t cap, double v, int as_float)
{
    int p = 0;
    if (!isfinite(v)) { snprintf(buf, cap, "null"); return; }
    for (p = 1; p <= 17; p++) {
        double back = 0.0;
        snprintf(buf, cap, "%.*g", p, v);
        back = strtod(buf, NULL);
        if (memcmp(&back, &v, sizeof v) == 0) break;
    }
    if (as_float && strpbrk(buf, ".eEn") == NULL && strlen(buf) + 3 <= cap) strcat(buf, ".0");
}

static void ic_json_put_str(FILE *f, const char *s)
{
    fputc('"', f);
    for (; *s; s++) {
        const unsigned ch = (unsigned char)*s;
        if (ch == '"' || ch == '\\') { fputc('\\', f); fputc((int)ch, f); }
        else if (ch == '\n') fputs("\\n", f);
        else if (ch == '\r') fputs("\\r", f);
        else if (ch == '\t') fputs("\\t", f);
        else if (ch < 0x20) fprintf(f, "\\u%04x", ch);
        else fputc((int)ch, f);
    }
    fputc('"', f);
}

/* v as Python's json.dump(..., indent=1) lays it out, on a line at depth d. */
static void ic_json_put(FILE *f, const JsonValue *v, int d)
{
    char num[48];
    size_t k = 0, n = 0;
    switch (Json_type(v)) {
    case JSON_BOOL: fputs(Json_as_bool(v, 0) ? "true" : "false", f); break;
    case JSON_NUMBER: ic_json_shortest(num, sizeof num, Json_as_double(v, NAN), 0); fputs(num, f); break;
    case JSON_STRING: ic_json_put_str(f, Json_as_string(v)); break;
    case JSON_ARRAY:
        n = Json_array_len(v);
        if (n == 0) { fputs("[]", f); break; }
        fputs("[\n", f);
        for (k = 0; k < n; k++) {
            fprintf(f, "%*s", d + 1, "");
            ic_json_put(f, Json_array_get(v, k), d + 1);
            fputs(k + 1 < n ? ",\n" : "\n", f);
        }
        fprintf(f, "%*s]", d, "");
        break;
    case JSON_OBJECT:
        n = Json_object_len(v);
        if (n == 0) { fputs("{}", f); break; }
        fputs("{\n", f);
        for (k = 0; k < n; k++) {
            fprintf(f, "%*s", d + 1, "");
            ic_json_put_str(f, Json_object_key(v, k));
            fputs(": ", f);
            ic_json_put(f, Json_object_value_at(v, k), d + 1);
            fputs(k + 1 < n ? ",\n" : "\n", f);
        }
        fprintf(f, "%*s}", d, "");
        break;
    case JSON_NULL:
    default: fputs("null", f); break;
    }
}

/* The same double, bit for bit (a value written and read back must be this exactly). */
static int ic_same_bits(double a, double b)
{
    return memcmp(&a, &b, sizeof a) == 0;
}

/* Two parsed values the same: types, numbers bit for bit, strings, members in order. */
static int ic_json_same(const JsonValue *a, const JsonValue *b)
{
    size_t k = 0;
    if (Json_type(a) != Json_type(b)) return 0;
    switch (Json_type(a)) {
    case JSON_BOOL: return Json_as_bool(a, 0) == Json_as_bool(b, 0);
    case JSON_NUMBER: return ic_same_bits(Json_as_double(a, 0.0), Json_as_double(b, 0.0));
    case JSON_STRING: return strcmp(Json_as_string(a), Json_as_string(b)) == 0;
    case JSON_ARRAY:
        if (Json_array_len(a) != Json_array_len(b)) return 0;
        for (k = 0; k < Json_array_len(a); k++) if (!ic_json_same(Json_array_get(a, k), Json_array_get(b, k))) return 0;
        return 1;
    case JSON_OBJECT:
        if (Json_object_len(a) != Json_object_len(b)) return 0;
        for (k = 0; k < Json_object_len(a); k++)
            if (strcmp(Json_object_key(a, k), Json_object_key(b, k)) != 0 ||
                !ic_json_same(Json_object_value_at(a, k), Json_object_value_at(b, k))) return 0;
        return 1;
    default: return 1;
    }
}

/* tmp takes path's place in one step: a reader sees the old file or the new one, never half of one. */
static int ic_replace_file(const char *tmp, const char *path)
{
#ifdef _WIN32
    return MoveFileExA(tmp, path, MOVEFILE_REPLACE_EXISTING | MOVEFILE_WRITE_THROUGH) ? 0 : -1;
#else
    return rename(tmp, path) == 0 ? 0 : -1;
#endif
}

static FILE *ic_open_tmp(const char *path, char *tmp, size_t cap)
{
    if (snprintf(tmp, cap, "%s.tmp", path) >= (int)cap || ves_ensure_parent_dir(path) != 0) return NULL;
    return fopen(tmp, "wb");
}

static int ic_close_tmp(FILE *f, const char *tmp, const char *path)
{
    int ok = !ferror(f);
    if (fclose(f) != 0) ok = 0;
    if (!ok || ic_replace_file(tmp, path) != 0) {
        remove(tmp);
        fprintf(stderr, "ERROR: cannot write %s\n", path);
        return -1;
    }
    return 0;
}

/* One pass's config as the per-pass stage wrote it: its id and folder, the request's pass settings, its source. */
static int ic_chain_write_config(const char *path, const char *id, const char *out_dir, const JsonValue *settings,
                                 const char *source)
{
    char tmp[2900];
    FILE *f = ic_open_tmp(path, tmp, sizeof tmp);
    size_t k = 0;
    if (f == NULL) { fprintf(stderr, "ERROR: cannot write %s\n", path); return -1; }
    fputs("{\n \"id\": ", f); ic_json_put_str(f, id);
    fputs(",\n \"out_dir\": ", f); ic_json_put_str(f, out_dir);
    for (k = 0; k < Json_object_len(settings); k++) {
        fputs(",\n ", f);
        ic_json_put_str(f, Json_object_key(settings, k));
        fputs(": ", f);
        ic_json_put(f, Json_object_value_at(settings, k), 1);
    }
    fputs(",\n \"surface_source\": ", f); ic_json_put_str(f, source);
    fputs("\n}", f);
    return ic_close_tmp(f, tmp, path);
}

/* A surface's vertex count as its report records it ("vertices", else surface.verts); 0 without one. */
static double ic_report_vertices(const char *stem)
{
    char path[2800];
    const char *err = NULL;
    Arena_T a = NULL;
    const JsonValue *r = NULL;
    FILE *probe = NULL;
    double v = 0.0;
    if (snprintf(path, sizeof path, "%s_surface_report.json", stem) >= (int)sizeof path) return 0.0;
    probe = fopen(path, "rb");
    if (probe == NULL) return 0.0;
    fclose(probe);
    a = Arena_new();
    r = Json_parse_file(a, path, &err);
    if (r != NULL) {
        v = Json_member_double(r, "vertices", 0.0);
        if (!(v > 0.0)) v = Json_member_double(Json_object_get(r, "surface"), "verts", 0.0);
    }
    Arena_dispose(&a);
    return isfinite(v) && v > 0.0 ? v : 0.0;
}

/* What the summary records of one pass: vertices and gain when it ends; adopted charts and seconds from its report. */
typedef struct IcChainPass {
    double verts, gain;            /* gain NAN: no vertex count before it */
    double adopted, seconds;       /* NAN: not in the report */
    int facts;                     /* adopted / seconds read */
} IcChainPass;

static int ic_chain_pass_facts(IcChainPass *p, const char *stem)
{
    char path[2800];
    const char *err = NULL;
    Arena_T a = NULL;
    const JsonValue *r = NULL, *ad = NULL;
    if (p->facts) return 0;
    if (snprintf(path, sizeof path, "%s_surface_report.json", stem) >= (int)sizeof path) return -1;
    a = Arena_new();
    r = Json_parse_file(a, path, &err);
    if (r == NULL) { fprintf(stderr, "ERROR: cannot read %s (%s)\n", path, err ? err : "?"); Arena_dispose(&a); return -1; }
    ad = Json_object_get(Json_object_get(r, "adopt"), "charts");
    p->adopted = Json_type(ad) == JSON_NUMBER ? Json_as_double(ad, NAN) : NAN;
    p->seconds = Json_type(Json_object_get(r, "seconds")) == JSON_NUMBER ? Json_member_double(r, "seconds", NAN) : NAN;
    p->facts = 1;
    Arena_dispose(&a);
    return 0;
}

/* The close-up summary closeup_fibers.py's finish() wrote: the final pass, every pass so far, why they ended. */
static int ic_chain_write_summary(const JsonValue *sj, const char *final_dir, const char *final_stem,
                                  const IcChainPass *p, size_t n, const char *stop)
{
    const char *path = Json_as_string(Json_object_get(sj, "path"));
    char tmp[2900], num[48];
    FILE *f = path != NULL ? ic_open_tmp(path, tmp, sizeof tmp) : NULL;
    size_t k = 0;
    if (f == NULL) { fprintf(stderr, "ERROR: cannot write the close-up summary %s\n", path ? path : "(none)"); return -1; }
    fputs("{\n \"area\": ", f); ic_json_put_str(f, Json_as_string(Json_object_get(sj, "area")));
    fputs(",\n \"operation\": ", f); ic_json_put_str(f, Json_as_string(Json_object_get(sj, "operation")));
    fputs(",\n \"final_dir\": ", f); ic_json_put_str(f, final_dir);
    fputs(",\n \"final_stem\": ", f); ic_json_put_str(f, final_stem);
    fputs(",\n \"seed_box\": ", f); ic_json_put(f, Json_object_get(sj, "seed_box"), 1);
    fputs(",\n \"passes\": ", f);
    if (n == 0) fputs("[]", f);
    else {
        fputs("[\n", f);
        for (k = 0; k < n; k++) {
            fprintf(f, "  {\n   \"k\": %zu,\n   \"verts\": %.0f,\n   \"gain\": ", k + 1, p[k].verts);
            ic_json_shortest(num, sizeof num, p[k].gain, 1); fputs(num, f);
            fputs(",\n   \"adopted\": ", f);
            if (isfinite(p[k].adopted)) fprintf(f, "%.0f", p[k].adopted); else fputs("null", f);
            fputs(",\n   \"seconds\": ", f);
            ic_json_shortest(num, sizeof num, p[k].seconds, 1); fputs(num, f);
            fputs(k + 1 < n ? "\n  },\n" : "\n  }\n", f);
        }
        fputs(" ]", f);
    }
    fputs(",\n \"stop_reason\": ", f); ic_json_put_str(f, stop);
    fputs("\n}", f);
    return ic_close_tmp(f, tmp, path);
}

/* The stage job's shape; NULL when it is usable, else what is wrong with it. */
static const char *ic_chain_job_problem(const JsonValue *j)
{
    const JsonValue *settings = Json_object_get(j, "pass_config"), *sj = Json_object_get(j, "summary");
    const char *schema = Json_as_string(Json_object_get(j, "schema"));
    const double np = Json_member_double(j, "passes", NAN), gain = Json_member_double(j, "min_gain", NAN);
    static const char *const req[] = { "base", "stem", "source" };
    static const char *const sreq[] = { "path", "area", "operation" };
    static const char *const own[] = { "id", "out_dir", "surface_source" };
    size_t k = 0;
    if (schema == NULL || strcmp(schema, "ink-closeup-complete-chain-v2") != 0) return "schema is not ink-closeup-complete-chain-v2";
    for (k = 0; k < 3; k++) {
        const char *v = Json_as_string(Json_object_get(j, req[k]));
        if (v == NULL || v[0] == 0) return "base, stem and source must be strings";
        v = Json_as_string(Json_object_get(sj, sreq[k]));
        if (v == NULL || v[0] == 0) return "summary needs path, area and operation";
        if (Json_object_get(settings, own[k]) != NULL) return "pass_config may not name id, out_dir or surface_source";
    }
    if (!(np >= 1.0 && np <= 64.0) || floor(np) != np) return "passes must be a whole number from 1 to 64";
    if (!isfinite(gain)) return "min_gain must be a number";
    if (Json_type(settings) != JSON_OBJECT) return "pass_config must be an object";
    if (Json_object_get(sj, "seed_box") == NULL) return "summary needs seed_box";
    if (Json_object_get(j, "keep_pass_meshes") != NULL && Json_type(Json_object_get(j, "keep_pass_meshes")) != JSON_BOOL)
        return "keep_pass_meshes must be true or false";
    return NULL;
}

static int ic_file_exists(const char *stem, const char *suffix)
{
    char path[2900];
    FILE *f = NULL;
    if (snprintf(path, sizeof path, "%s%s", stem, suffix) >= (int)sizeof path) return 0;
    f = fopen(path, "rb");
    if (f == NULL) return 0;
    fclose(f);
    return 1;
}

/* Where a retry resumes: the passes recorded in order (their reports), up to the last of them whose mesh is on disk --
 * the next pass's source.  A pass recorded after it has no mesh (an intermediate pass) and runs again. */
static size_t ic_chain_resume_point(char (*stems)[2048], size_t np)
{
    size_t k = 0, upto = 0;
    for (k = 0; k < np; k++) {
        if (!ic_file_exists(stems[k], "_surface_report.json")) break;
        if (ic_file_exists(stems[k], "_surface.vmesh")) upto = k + 1;
    }
    return upto;
}

/* closeup_fibers.py's fiber_complete stage in one process (2026-10-30; the user: the stage "took 30.2 minutes. Get that
 * to 2 minutes", and "make sure that everything is C and not python").  job.json {"schema":
 * "ink-closeup-complete-chain-v2", "base", "stem", "source", "passes" (1..64), "min_gain", "pass_config": the request's
 * pass settings, "summary": {"path", "area", "operation", "seed_box"}}.  Pass k runs in <base>/grow<k> as
 * <stem>_grow<k> on the pass before it (pass 1 on "source"), from a config written there as the per-pass stage wrote
 * it.  A pass whose report exists is done (a retry resumes after it); otherwise it runs on the previous pass's surface
 * in memory, the pool meshes read once, and its files are written behind the next pass.  The passes end where a pass
 * gains less than min_gain of the vertices before it (pass 1: the source's report), else at the pass limit.  Progress:
 * the stage's two status lines per pass.  Only the final pass's mesh is written unless "keep_pass_meshes" (the user,
 * 2026-10-30: write only the final completion pass's mesh -- "yes"); every pass keeps its config, report and adoption
 * CSVs.  A retry resumes after the last recorded pass whose mesh is on disk.  The close-up summary names a pass only
 * once its files are complete and its mesh is on disk (at the latest, the final pass at the end), with every pass's
 * vertices, gain, adopted charts and seconds, and why the passes ended. */
static int ic_stage_complete_chain(const char *job_path)
{
    Arena_T arena = Arena_new();
    const char *err = NULL, *problem = NULL, *base = NULL, *stem = NULL, *source0 = NULL;
    const JsonValue *j = Json_parse_file(arena, job_path, &err), *settings = NULL, *sj = NULL;
    IcPoolCache *pc = (IcPoolCache *)calloc(1, sizeof *pc);
    IcConfig *c = (IcConfig *)calloc(1, sizeof *c);
    IcChainPass *pass = (IcChainPass *)calloc(64, sizeof *pass);
    IcSurfaceMem src;
    IcWriter wr;
    char prev_stem[2048] = "", report[2800], count[48], folder[2048], name[256], cfg[2400];
    char (*stems)[2048] = (char (*)[2048])calloc(64, sizeof *stems);
    double previous = 0.0, min_gain = 0.0, t0 = ves_clock_sec();
    size_t np = 0, k = 0, ran = 0, landed = 0, reached = 0, resume = 0;
    int have_src = 0, rc = 1;
    const char *stop = "pass_limit";
    memset(&src, 0, sizeof src); memset(&wr, 0, sizeof wr);
    if (j == NULL || pc == NULL || c == NULL || pass == NULL || stems == NULL) {
        fprintf(stderr, "ERROR: cannot read the completion job %s (%s)\n", job_path, err ? err : "no memory");
        goto out;
    }
    if ((problem = ic_chain_job_problem(j)) != NULL) { fprintf(stderr, "ERROR: %s: %s\n", job_path, problem); goto out; }
    base = Json_as_string(Json_object_get(j, "base"));
    stem = Json_as_string(Json_object_get(j, "stem"));
    source0 = Json_as_string(Json_object_get(j, "source"));
    settings = Json_object_get(j, "pass_config");
    sj = Json_object_get(j, "summary");
    np = (size_t)Json_member_double(j, "passes", 1.0);
    min_gain = Json_member_double(j, "min_gain", 0.01);
    wr.keep_all = Json_as_bool(Json_object_get(j, "keep_pass_meshes"), 0);
    wr.min_gain = min_gain;
    previous = ic_report_vertices(source0);
    for (k = 0; k < np; k++)
        if (snprintf(folder, sizeof folder, "%s%cgrow%zu", base, IC_SEP, k + 1) >= (int)sizeof folder ||
            snprintf(name, sizeof name, "%s_grow%zu", stem, k + 1) >= (int)sizeof name ||
            snprintf(stems[k], sizeof stems[k], "%s%c%s", folder, IC_SEP, name) >= (int)sizeof stems[k]) {
            fprintf(stderr, "ERROR: pass %zu: path too long under %s\n", k + 1, base);
            goto out;
        }
    resume = ic_chain_resume_point(stems, np);
    for (k = 0; k < np; k++) {
        double verts = 0.0;
        int ran_here = 0;
        snprintf(folder, sizeof folder, "%s%cgrow%zu", base, IC_SEP, k + 1);
        snprintf(name, sizeof name, "%s_grow%zu", stem, k + 1);
        if (snprintf(cfg, sizeof cfg, "%s%c%s.json", folder, IC_SEP, name) >= (int)sizeof cfg) {
            fprintf(stderr, "ERROR: pass %zu: path too long under %s\n", k + 1, base);
            goto out;
        }
        if (ic_chain_write_config(cfg, name, folder, settings, k == 0 ? source0 : prev_stem) != 0) goto out;
        free(c->chart_wrap);
        if (ic_config(cfg, c) != 0) { fprintf(stderr, "ERROR: pass %zu: cannot read %s\n", k + 1, cfg); goto out; }
        printf("[orch] status Close-up completion pass %zu of %zu\n", k + 1, np);
        fflush(stdout);
        snprintf(report, sizeof report, "%s_surface_report.json", stems[k]);
        if (k < resume) {              /* done before (a retry): its files, not memory, are the next source */
            verts = ic_report_vertices(stems[k]);
            if (!(verts > 0.0)) { fprintf(stderr, "ERROR: %s holds no surface\n", report); goto out; }
            if (have_src) { ic_surface_mem_free(&src); have_src = 0; }
            if (ic_writer_wait(&wr) != 0) goto out;
            fprintf(stderr, "[ink_closeup] chain: pass %zu was complete (%s)\n", k + 1, report);
        } else {
            IcSurfaceMem next;
            const int mem = have_src && ic_same_stem(c->surface_source, prev_stem);
            if (have_src && !mem) { ic_surface_mem_free(&src); have_src = 0; }
            if (!mem && ic_writer_wait(&wr) != 0) goto out;      /* a source on disk must have landed */
            wr.last = k + 1 == np;
            wr.previous = previous;
            if (ic_complete_pass(c, mem ? &src : NULL, &next, pc, &wr) != 0) goto out;
            src = next; have_src = 1; ran++; ran_here = 1;
            verts = (double)next.m.nv;
        }
        reached = k + 1;
        memcpy(prev_stem, stems[k], strlen(stems[k]) + 1);
        pass[k].verts = verts;
        pass[k].gain = previous > 0.0 ? verts / previous - 1.0 : NAN;
        pass[k].adopted = pass[k].seconds = NAN;
        ic_thousands(count, sizeof count, (size_t)verts);
        if (previous > 0.0) printf("[orch] status Close-up pass %zu: %s vertices; %+.1f%% change\n", k + 1, count,
                                   pass[k].gain * 100.0);
        else printf("[orch] status Close-up pass %zu: %s vertices; initial completion\n", k + 1, count);
        fflush(stdout);
        {   /* the passes whose files are complete: before this one when it was queued here, else this one too */
            const size_t now = ran_here ? k : k + 1;
            if (now > landed && ic_file_exists(stems[now - 1], "_surface.vmesh")) {
                size_t q = 0;
                for (q = 0; q < now; q++) if (ic_chain_pass_facts(&pass[q], stems[q]) != 0) goto out;
                landed = now;
                snprintf(folder, sizeof folder, "%s%cgrow%zu", base, IC_SEP, landed);
                snprintf(name, sizeof name, "%s_grow%zu", stem, landed);
                if (ic_chain_write_summary(sj, folder, name, pass, landed, "continuing") != 0) goto out;
            }
        }
        if (previous > 0.0 && pass[k].gain < min_gain) {
            stop = "min_gain";
            if (ran_here && wr.queued_no_mesh) { fprintf(stderr, "ERROR: pass %zu ended the passes without its mesh\n", k + 1); goto out; }
            break;
        }
        previous = verts;
    }
    if (ic_writer_wait(&wr) != 0) goto out;
    if (reached == 0 || !ic_file_exists(stems[reached - 1], "_surface.vmesh")) {
        fprintf(stderr, "ERROR: the final pass %zu has no mesh on disk\n", reached);
        goto out;
    }
    for (k = 0; k < reached; k++) if (ic_chain_pass_facts(&pass[k], stems[k]) != 0) goto out;
    snprintf(folder, sizeof folder, "%s%cgrow%zu", base, IC_SEP, reached);
    snprintf(name, sizeof name, "%s_grow%zu", stem, reached);
    if (ic_chain_write_summary(sj, folder, name, pass, reached, stop) != 0) goto out;
    fprintf(stderr, "[ink_closeup] chain: %zu of %zu passes run here, ended by %s (%.1fs)\n", ran, np, stop,
            ves_clock_sec() - t0);
    rc = 0;
out:
    if (ic_writer_wait(&wr) != 0) rc = 1;
    if (have_src) ic_surface_mem_free(&src);
    if (pc != NULL) ic_pool_cache_free(pc);
    if (c != NULL) free(c->chart_wrap);
    free(pc); free(c); free(pass); free(stems);
    Arena_dispose(&arena);
    return rc;
}

/* The stage's files without a close-up: pass configs, report vertex counts, the summary, numbers, the job check. */
static int ic_chain_selftest(void)
{
    char dir[600], path[800], stem[800], num[48];
    const char *err = NULL;
    int fails = 0;
    size_t k = 0;
    Arena_T a = Arena_new();
    const char *settings_text =
        "{\"layout\": \"D:\\\\vesuvius\\\\takein\\\\x\\\\closeup_context\", \"voxel_um\": 8.64, \"z_top_to_bottom\": true,"
        " \"seed_box\": [15636359.760440797, 14.0, 1e-7, -3], \"grow_px\": 0, \"complete_gaps\": true,"
        " \"adopt\": {\"rounds\": 3}, \"empty\": {}, \"none\": [], \"q\": \"a \\\"quoted\\\" \\\\ name\\n\", \"n\": null}";
    const JsonValue *settings = Json_parse(a, settings_text, &err), *back = NULL;
    static const double nums[] = { 0.1, 1.0 / 3.0, 1e-7, 262.57, 14.0, 1e300, 6800516.0, -0.0, 0.14372101917319613,
                                   -0.02, 1e16, 123456.0, 5e-324 };
    if (ves_selftest_dir(dir, sizeof dir, "ink_closeup_chain") != 0 || settings == NULL) { Arena_dispose(&a); return 1; }
    for (k = 0; k < sizeof nums / sizeof nums[0]; k++) {   /* shortest text, read back bit for bit */
        double v = 0.0;
        ic_json_shortest(num, sizeof num, nums[k], 1);
        v = strtod(num, NULL);
        if (memcmp(&v, &nums[k], sizeof v) != 0) fails++;
    }
    ic_json_shortest(num, sizeof num, 14.0, 1); if (strcmp(num, "14.0") != 0) fails++;
    ic_json_shortest(num, sizeof num, 14.0, 0); if (strcmp(num, "14") != 0) fails++;
    ic_json_shortest(num, sizeof num, 1e-5, 1); if (strcmp(num, "1e-05") != 0) fails++;
    ic_json_shortest(num, sizeof num, 0.1, 1); if (strcmp(num, "0.1") != 0) fails++;
    ic_json_shortest(num, sizeof num, NAN, 1); if (strcmp(num, "null") != 0) fails++;
    snprintf(path, sizeof path, "%s/grow1/x_grow1.json", dir);
    if (ic_chain_write_config(path, "x_grow1", "D:\\out\\grow1", settings, "D:\\src \"stem\"") != 0) fails++;
    else if ((back = Json_parse_file(a, path, &err)) == NULL) fails++;
    else {
        if (strcmp(Json_as_string(Json_object_get(back, "id")), "x_grow1") != 0 ||
            strcmp(Json_as_string(Json_object_get(back, "out_dir")), "D:\\out\\grow1") != 0 ||
            strcmp(Json_as_string(Json_object_get(back, "surface_source")), "D:\\src \"stem\"") != 0) fails++;
        if (Json_object_len(back) != Json_object_len(settings) + 3) fails++;
        for (k = 0; k < Json_object_len(settings); k++)
            if (strcmp(Json_object_key(back, k + 2), Json_object_key(settings, k)) != 0 ||
                !ic_json_same(Json_object_value_at(back, k + 2), Json_object_value_at(settings, k))) fails++;
    }
    {   /* report vertex counts: "vertices" first, then surface.verts, else none */
        static const char *const reports[] = { "{\"vertices\": 5945957, \"surface\": {\"verts\": 1}}",
                                               "{\"surface\": {\"verts\": 77}}", "{\"vertices\": 0, \"surface\": {\"verts\": 5}}",
                                               "{}" };
        static const double want[] = { 5945957.0, 77.0, 5.0, 0.0 };
        for (k = 0; k < 4; k++) {
            FILE *f = NULL;
            snprintf(stem, sizeof stem, "%s/r%zu", dir, k);
            snprintf(path, sizeof path, "%s_surface_report.json", stem);
            f = fopen(path, "wb");
            if (f == NULL) { fails++; continue; }
            fputs(reports[k], f);
            fclose(f);
            if (!ic_same_bits(ic_report_vertices(stem), want[k])) fails++;
        }
        snprintf(stem, sizeof stem, "%s/absent", dir);
        if (!ic_same_bits(ic_report_vertices(stem), 0.0)) fails++;
    }
    {   /* the summary: two passes, a null gain and a null adopted count, read back */
        IcChainPass p[2];
        char sjt[1200], bs[900], ps[900];
        const JsonValue *sj = NULL, *sb = NULL, *pp = NULL;
        FILE *f = NULL;
        memset(p, 0, sizeof p);
        p[0].verts = 6800516; p[0].gain = NAN; p[0].adopted = 300; p[0].seconds = 262.57;
        p[1].verts = 7564718; p[1].gain = 0.11237411984619983; p[1].adopted = NAN; p[1].seconds = 262.0;
        snprintf(path, sizeof path, "%s/cu_x_closeup_summary.json", dir);
        ic_json_str(ps, sizeof ps, path);
        snprintf(sjt, sizeof sjt, "{\"path\": %s, \"area\": \"cu_x\", \"operation\": \"align-closeup\", "
                 "\"seed_box\": [15636359.760440797, 14.0, 15645124.760440797, 3594.0]}", ps);
        sj = Json_parse(a, sjt, &err);
        if (sj == NULL || ic_chain_write_summary(sj, "D:\\out\\grow2", "x_grow2", p, 2, "min_gain") != 0) fails++;
        else if ((sb = Json_parse_file(a, path, &err)) == NULL) fails++;
        else {
            pp = Json_object_get(sb, "passes");
            if (strcmp(Json_as_string(Json_object_get(sb, "area")), "cu_x") != 0 ||
                strcmp(Json_as_string(Json_object_get(sb, "operation")), "align-closeup") != 0 ||
                strcmp(Json_as_string(Json_object_get(sb, "final_dir")), "D:\\out\\grow2") != 0 ||
                strcmp(Json_as_string(Json_object_get(sb, "final_stem")), "x_grow2") != 0 ||
                strcmp(Json_as_string(Json_object_get(sb, "stop_reason")), "min_gain") != 0 ||
                !ic_json_same(Json_object_get(sb, "seed_box"), Json_object_get(sj, "seed_box")) ||
                Json_array_len(pp) != 2) fails++;
            else {
                const JsonValue *p0 = Json_array_get(pp, 0), *p1 = Json_array_get(pp, 1);
                if (!ic_same_bits(Json_member_double(p0, "k", 0), 1.0) ||
                    !ic_same_bits(Json_member_double(p0, "verts", 0), 6800516.0) ||
                    Json_type(Json_object_get(p0, "gain")) != JSON_NULL ||
                    !ic_same_bits(Json_member_double(p0, "adopted", 0), 300.0) ||
                    !ic_same_bits(Json_member_double(p0, "seconds", 0), 262.57)) fails++;
                if (!ic_same_bits(Json_member_double(p1, "k", 0), 2.0) ||
                    !ic_same_bits(Json_member_double(p1, "gain", 0), 0.11237411984619983) ||
                    Json_type(Json_object_get(p1, "adopted")) != JSON_NULL ||
                    !ic_same_bits(Json_member_double(p1, "seconds", 0), 262.0)) fails++;
            }
        }
        f = fopen(path, "rb");     /* no leftover temporary beside it */
        if (f != NULL) fclose(f);
        snprintf(bs, sizeof bs, "%s.tmp", path);
        f = fopen(bs, "rb");
        if (f != NULL) { fclose(f); fails++; }
    }
    {   /* the resume point: reports in order, up to the last pass with its mesh */
        static char st[5][2048];
        static const int has_report[5] = { 1, 1, 1, 1, 0 }, has_mesh[5] = { 0, 1, 0, 0, 1 };
        size_t q = 0;
        for (q = 0; q < 5; q++) {
            FILE *f = NULL;
            snprintf(st[q], sizeof st[q], "%s/rp%zu", dir, q);
            snprintf(path, sizeof path, "%s_surface_report.json", st[q]);
            if (has_report[q] && (f = fopen(path, "wb")) != NULL) { fputs("{}", f); fclose(f); }
            snprintf(path, sizeof path, "%s_surface.vmesh", st[q]);
            if (has_mesh[q] && (f = fopen(path, "wb")) != NULL) fclose(f);
        }
        if (ic_chain_resume_point(st, 5) != 2) fails++;           /* pass 2 has the last mesh; 3, 4 run again */
        if (ic_chain_resume_point(st, 1) != 0) fails++;           /* pass 1: recorded, no mesh */
        snprintf(path, sizeof path, "%s_surface.vmesh", st[3]);
        { FILE *f = fopen(path, "wb"); if (f) fclose(f); }
        if (ic_chain_resume_point(st, 5) != 4) fails++;           /* pass 5 has a mesh but no report: not recorded */
    }
    {   /* the job check */
        static const char *const bad[] = {
            "{\"schema\": \"ink-closeup-complete-chain-v1\"}",
            "{\"schema\": \"ink-closeup-complete-chain-v2\", \"base\": \"b\", \"stem\": \"s\", \"source\": \"x\", \"passes\": 0,"
            " \"min_gain\": 0.01, \"pass_config\": {}, \"summary\": {\"path\": \"p\", \"area\": \"a\", \"operation\": \"o\", \"seed_box\": []}}",
            "{\"schema\": \"ink-closeup-complete-chain-v2\", \"base\": \"b\", \"stem\": \"s\", \"source\": \"x\", \"passes\": 2.5,"
            " \"min_gain\": 0.01, \"pass_config\": {}, \"summary\": {\"path\": \"p\", \"area\": \"a\", \"operation\": \"o\", \"seed_box\": []}}",
            "{\"schema\": \"ink-closeup-complete-chain-v2\", \"base\": \"b\", \"stem\": \"s\", \"source\": \"x\", \"passes\": 65,"
            " \"min_gain\": 0.01, \"pass_config\": {}, \"summary\": {\"path\": \"p\", \"area\": \"a\", \"operation\": \"o\", \"seed_box\": []}}",
            "{\"schema\": \"ink-closeup-complete-chain-v2\", \"base\": \"b\", \"stem\": \"s\", \"source\": \"x\", \"passes\": 6,"
            " \"min_gain\": 0.01, \"pass_config\": {\"id\": \"no\"}, \"summary\": {\"path\": \"p\", \"area\": \"a\", \"operation\": \"o\", \"seed_box\": []}}",
            "{\"schema\": \"ink-closeup-complete-chain-v2\", \"base\": \"b\", \"stem\": \"s\", \"source\": \"x\", \"passes\": 6,"
            " \"min_gain\": 0.01, \"pass_config\": {}, \"summary\": {\"path\": \"p\", \"area\": \"a\", \"operation\": \"o\"}}",
            "{\"schema\": \"ink-closeup-complete-chain-v2\", \"base\": \"b\", \"stem\": \"s\", \"source\": \"x\", \"passes\": 6,"
            " \"min_gain\": 0.01, \"keep_pass_meshes\": 1, \"pass_config\": {}, \"summary\": {\"path\": \"p\", \"area\": \"a\","
            " \"operation\": \"o\", \"seed_box\": []}}" };
        const char *good = "{\"schema\": \"ink-closeup-complete-chain-v2\", \"base\": \"b\", \"stem\": \"s\", \"source\": \"x\","
                           " \"passes\": 6, \"min_gain\": 0.01, \"pass_config\": {\"layout\": \"l\"}, \"summary\": {\"path\": \"p\","
                           " \"area\": \"a\", \"operation\": \"o\", \"seed_box\": [1, 2, 3, 4]}}";
        for (k = 0; k < sizeof bad / sizeof bad[0]; k++) {
            const JsonValue *jb = Json_parse(a, bad[k], &err);
            if (jb == NULL || ic_chain_job_problem(jb) == NULL) fails++;
        }
        if (ic_chain_job_problem(Json_parse(a, good, &err)) != NULL) fails++;
    }
    Arena_dispose(&a);
    fprintf(stderr, "  completion stage files (pass configs, report vertices, summary, numbers, resume point, job check): %s\n",
            fails ? "FAIL" : "ok");
    if (!fails) (void)ves_selftest_clean(dir);
    return fails;
}

/* ------------------------------------------------------------------ grow */

/* The data boundary: DATA nodes with a 4-neighbour inside the grid that is not DATA (where the sheet can still grow;
 * the grid's own edge is the target box's and does not count).  mark[k] = CLOSEUP_NODE_DATA there, else NONE; the
 * count. */
static size_t ic_data_boundary(const CloseupGrid *g, uint8_t *mark)
{
    size_t i = 0, j = 0, n = 0, ni = g->ni, nj = g->nj;
    memset(mark, CLOSEUP_NODE_NONE, ni * nj);
    for (j = 0; j < nj; j++)
        for (i = 0; i < ni; i++) {
            size_t k = j * ni + i;
            int edge = 0;
            if (g->state[k] != CLOSEUP_NODE_DATA) continue;
            edge = (i > 0 && g->state[k - 1] != CLOSEUP_NODE_DATA) || (i + 1 < ni && g->state[k + 1] != CLOSEUP_NODE_DATA) ||
                   (j > 0 && g->state[k - ni] != CLOSEUP_NODE_DATA) || (j + 1 < nj && g->state[k + ni] != CLOSEUP_NODE_DATA);
            if (edge) { mark[k] = CLOSEUP_NODE_DATA; n++; }
        }
    return n;
}

/* Charts of one pool already adopted, by pool chart id (ic_pool_load's exclusion list). */
typedef struct IcAdopted {
    uint8_t *on;
    size_t n;
} IcAdopted;

static int ic_adopted_add(IcAdopted *a, int32_t id)
{
    if (id < 0) return 0;
    if ((size_t)id >= a->n) {
        size_t n = a->n ? a->n : 4096;
        uint8_t *grown = NULL;
        while (n <= (size_t)id) n *= 2;
        grown = (uint8_t *)realloc(a->on, n);
        if (grown == NULL) return -1;
        memset(grown + a->n, 0, n - a->n);
        a->on = grown;
        a->n = n;
    }
    a->on[id] = 1;
    return 0;
}

/* Tiles already adopted, by key (parent chart id << 24 | tile cell): an open-addressing set of nonnegative int64. */
typedef struct IcKeySet {
    int64_t *slot;           /* -1 = empty */
    size_t cap, n;
} IcKeySet;

static size_t ic_key_hash(int64_t k, size_t cap)
{
    uint64_t x = (uint64_t)k * 0x9E3779B97F4A7C15ull;
    return (size_t)(x >> 17) & (cap - 1);
}

static int ic_keyset_has(const IcKeySet *s, int64_t k)
{
    size_t i = 0;
    if (s->cap == 0) return 0;
    for (i = ic_key_hash(k, s->cap); s->slot[i] != -1; i = (i + 1) & (s->cap - 1))
        if (s->slot[i] == k) return 1;
    return 0;
}

static size_t ic_keysets_total(const IcKeySet *s, int n)
{
    size_t t = 0;
    int k = 0;
    if (s == NULL) return 0;
    for (k = 0; k < n; k++) t += s[k].n;
    return t;
}

static int ic_keyset_add(IcKeySet *s, int64_t k)
{
    size_t i = 0;
    if (k < 0) return -1;
    if (2 * (s->n + 1) > s->cap) {                    /* keep the load under 1/2: grow and rehash */
        size_t cap = s->cap ? s->cap * 2 : 4096, j = 0;
        int64_t *slot = (int64_t *)malloc(cap * sizeof *slot);
        if (slot == NULL) return -1;
        for (j = 0; j < cap; j++) slot[j] = -1;
        for (j = 0; j < s->cap; j++) {
            size_t t = 0;
            if (s->slot[j] == -1) continue;
            for (t = ic_key_hash(s->slot[j], cap); slot[t] != -1; t = (t + 1) & (cap - 1)) {}
            slot[t] = s->slot[j];
        }
        free(s->slot);
        s->slot = slot;
        s->cap = cap;
    }
    for (i = ic_key_hash(k, s->cap); s->slot[i] != -1; i = (i + 1) & (s->cap - 1))
        if (s->slot[i] == k) return 0;
    s->slot[i] = k;
    s->n++;
    return 0;
}

/* Every chart of a loaded pool split into TILES (tile_px x tile_px cells of its own UV, faces by their UV centroid;
 * vertex-disjoint: a vertex on a cell border is copied into each tile that uses it).  A tile keeps its parent's
 * flattening and id (chart_id, face_chart) and gets keys[t] = parent id << 24 | its cell, the same on every read of
 * that chart.  closeup_adopt then judges each tile on its own contacts: a chart that runs onto the neighbouring wrap in
 * a dense bundle, or whose flattening drifts over its length, still gives the tiles that continue the sheet (the
 * cleaning doctrine's "cut the chart at the neck", done at the tile scale).  0, -1 on no memory or a pool chart
 * wider than 4096 tiles (left as it was). */
static int ic_pool_tile(IcPool *q, double tile_px, int64_t **keys_out)
{
    size_t c = 0, nt = 0, cap_t = 0, nv2 = 0, nf2 = 0, cap_v = 0, maxv = 0, maxf = 0;
    float *verts = NULL, *nrm = NULL, *inward = NULL;
    double *uv = NULL;
    int32_t *faces = NULL, *fchart = NULL, *chart_id = NULL, *vmap = NULL;
    size_t *vofs = NULL, *fofs = NULL, *order = NULL;
    uint32_t *cell = NULL;
    int64_t *keys = NULL;
    int rc = -1;
    *keys_out = NULL;
    if (!(tile_px > 0.0) || q->p.nch == 0) return 0;
    for (c = 0; c < q->p.nch; c++) {
        if (q->vofs[c + 1] - q->vofs[c] > maxv) maxv = q->vofs[c + 1] - q->vofs[c];
        if (q->fofs[c + 1] - q->fofs[c] > maxf) maxf = q->fofs[c + 1] - q->fofs[c];
    }
    cap_v = q->p.nv + q->p.nv / 4 + 16;
    cap_t = q->p.nch * 4 + 16;
    verts = (float *)malloc(cap_v * 3 * sizeof *verts); nrm = (float *)malloc(cap_v * 3 * sizeof *nrm);
    inward = (float *)malloc(cap_v * 3 * sizeof *inward); uv = (double *)malloc(cap_v * 2 * sizeof *uv);
    faces = (int32_t *)malloc((q->p.nf ? q->p.nf : 1) * 3 * sizeof *faces);
    fchart = (int32_t *)malloc((q->p.nf ? q->p.nf : 1) * sizeof *fchart);
    vofs = (size_t *)malloc((cap_t + 1) * sizeof *vofs); fofs = (size_t *)malloc((cap_t + 1) * sizeof *fofs);
    chart_id = (int32_t *)malloc(cap_t * sizeof *chart_id); keys = (int64_t *)malloc(cap_t * sizeof *keys);
    vmap = (int32_t *)malloc((maxv ? maxv : 1) * sizeof *vmap);
    order = (size_t *)malloc((maxf ? maxf : 1) * sizeof *order);
    cell = (uint32_t *)malloc((maxf ? maxf : 1) * sizeof *cell);
    if (!verts || !nrm || !inward || !uv || !faces || !fchart || !vofs || !fofs || !chart_id || !keys || !vmap || !order ||
        !cell)
        goto out;
    for (c = 0; c < q->p.nch; c++) {
        const size_t v0 = q->vofs[c], v1 = q->vofs[c + 1], f0 = q->fofs[c], f1 = q->fofs[c + 1], nfc = f1 - f0;
        double umin = 1e300, wmin = 1e300;
        size_t k = 0, a = 0, v = 0;
        for (v = v0; v < v1; v++) { umin = fmin(umin, q->uv[v * 2]); wmin = fmin(wmin, q->uv[v * 2 + 1]); }
        for (k = 0; k < nfc; k++) {
            const int32_t *fv = &q->faces[(f0 + k) * 3];
            const double cu = (q->uv[(size_t)fv[0] * 2] + q->uv[(size_t)fv[1] * 2] + q->uv[(size_t)fv[2] * 2]) / 3.0;
            const double cw = (q->uv[(size_t)fv[0] * 2 + 1] + q->uv[(size_t)fv[1] * 2 + 1] + q->uv[(size_t)fv[2] * 2 + 1]) / 3.0;
            const double ci = floor((cu - umin) / tile_px), cj = floor((cw - wmin) / tile_px);
            if (!(ci >= 0.0 && ci < 4096.0 && cj >= 0.0 && cj < 4096.0)) goto out;
            cell[k] = (uint32_t)ci | ((uint32_t)cj << 12);
            order[k] = k;
        }
        /* the chart's faces grouped by cell, in file order within a cell (a shell sort of the index by (cell, k):
         * charts hold 10^2..10^5 faces) */
        {
            size_t gap = nfc / 2;
            for (; gap > 0; gap /= 2)
                for (k = gap; k < nfc; k++) {
                    size_t t = order[k], j = k;
                    while (j >= gap && (cell[order[j - gap]] > cell[t] || (cell[order[j - gap]] == cell[t] && order[j - gap] > t))) {
                        order[j] = order[j - gap];
                        j -= gap;
                    }
                    order[j] = t;
                }
        }
        for (v = 0; v < v1 - v0; v++) vmap[v] = -1;
        for (a = 0; a < nfc;) {
            size_t b = a, base_v = nv2, z = 0;
            while (b < nfc && cell[order[b]] == cell[order[a]]) b++;
            if (nt + 1 > cap_t) {
                size_t nc = cap_t * 2;
                size_t *vo = (size_t *)realloc(vofs, (nc + 1) * sizeof *vofs), *fo = NULL;
                int32_t *ci2 = NULL;
                int64_t *ke = NULL;
                if (vo == NULL) goto out;
                vofs = vo;
                fo = (size_t *)realloc(fofs, (nc + 1) * sizeof *fofs); if (fo == NULL) goto out; fofs = fo;
                ci2 = (int32_t *)realloc(chart_id, nc * sizeof *chart_id); if (ci2 == NULL) goto out; chart_id = ci2;
                ke = (int64_t *)realloc(keys, nc * sizeof *keys); if (ke == NULL) goto out; keys = ke;
                cap_t = nc;
            }
            vofs[nt] = nv2; fofs[nt] = nf2;
            chart_id[nt] = q->chart_id[c];
            keys[nt] = ((int64_t)(q->chart_id[c] < 0 ? 0 : q->chart_id[c]) << 24) | (int64_t)cell[order[a]];
            for (z = a; z < b; z++) {
                const size_t f = f0 + order[z];
                int kk = 0;
                for (kk = 0; kk < 3; kk++) {
                    const size_t ov = (size_t)q->faces[f * 3 + (size_t)kk], lv = ov - v0;
                    if (vmap[lv] < 0 || (size_t)vmap[lv] < base_v) {      /* not yet in THIS tile */
                        if (nv2 + 1 > cap_v) {
                            size_t ncv = cap_v + cap_v / 2 + 1024;
                            float *t1 = (float *)realloc(verts, ncv * 3 * sizeof *verts), *t2 = NULL, *t3 = NULL;
                            double *t4 = NULL;
                            if (t1 == NULL) goto out;
                            verts = t1;
                            t2 = (float *)realloc(nrm, ncv * 3 * sizeof *nrm); if (t2 == NULL) goto out; nrm = t2;
                            t3 = (float *)realloc(inward, ncv * 3 * sizeof *inward); if (t3 == NULL) goto out; inward = t3;
                            t4 = (double *)realloc(uv, ncv * 2 * sizeof *uv); if (t4 == NULL) goto out; uv = t4;
                            cap_v = ncv;
                        }
                        memcpy(&verts[nv2 * 3], &q->verts[ov * 3], 3 * sizeof *verts);
                        memcpy(&nrm[nv2 * 3], &q->nrm[ov * 3], 3 * sizeof *nrm);
                        memcpy(&inward[nv2 * 3], &q->inward[ov * 3], 3 * sizeof *inward);
                        uv[nv2 * 2] = q->uv[ov * 2]; uv[nv2 * 2 + 1] = q->uv[ov * 2 + 1];
                        vmap[lv] = (int32_t)nv2++;
                    }
                    faces[nf2 * 3 + (size_t)kk] = vmap[lv];
                }
                fchart[nf2] = q->face_chart[f];
                nf2++;
            }
            nt++;
            a = b;
        }
    }
    vofs[nt] = nv2; fofs[nt] = nf2;
    free(q->verts); free(q->nrm); free(q->inward); free(q->uv); free(q->faces); free(q->face_chart);
    free(q->chart_id); free(q->vofs); free(q->fofs);
    q->verts = verts; q->nrm = nrm; q->inward = inward; q->uv = uv; q->faces = faces; q->face_chart = fchart;
    q->chart_id = chart_id; q->vofs = vofs; q->fofs = fofs;
    q->p.nv = nv2; q->p.nf = nf2; q->p.nch = nt;
    q->p.verts = q->verts; q->p.nrm = q->nrm; q->p.inward = q->inward; q->p.uv = q->uv; q->p.faces = q->faces;
    q->p.vofs = q->vofs; q->p.fofs = q->fofs;
    verts = nrm = inward = NULL; uv = NULL; faces = fchart = chart_id = NULL; vofs = fofs = NULL;
    *keys_out = keys;
    keys = NULL;
    rc = 0;
out:
    free(verts); free(nrm); free(inward); free(uv); free(faces); free(fchart); free(chart_id); free(vofs); free(fofs);
    free(keys); free(vmap); free(order); free(cell);
    return rc;
}

/* The method name the surface report gives the growth's final flattening (CloseupSurface_flatten_slim). */
#define IC_GROW_FLATTEN_METHOD "slim_scaf_fringe_v2"

/* The grown surface's last steps, shared by the growth and --reflatten (2026-10-15, so a later flattening can be tried on
 * a finished growth): the final flattening (CloseupSurface_flatten_slim; the seed's frame when reflatten is off or the
 * flattening fails), the mesh, rows along the axis, and the surface files <out_dir>/<id>_surface{.vmesh, _nrm.f32,
 * _attr.u8, _chart.i32, _f0.u8 when have_f0}.  *flat, *reflattened, *nv_out, *nf_out report it; 0 on success. */
static int ic_grow_surface_files(const IcConfig *c, CloseupGrid *g, const AsmAxis *axis, int32_t adopt_first,
                                 int32_t f0_chart, int have_f0, double t0, CloseupFlatStats *flat, int *reflattened,
                                 size_t *nv_out, size_t *nf_out, CloseupDedupStats *dd)
{
    char path[2600];
    const size_t nn = g->ni * g->nj;
    uint8_t *oa = NULL, *of0 = NULL;
    int32_t *faces = NULL, *nodes = NULL, *oc = NULL;
    float *verts = NULL, *nrm = NULL;
    double *uvn = NULL, *muv = NULL, shift = 0.0;
    size_t v = 0, nv = 0, nf = 0;
    int rc = 1;
    memset(flat, 0, sizeof *flat);
    memset(dd, 0, sizeof *dd);
    *reflattened = 0;
    if (CLOSEUP_APRON_TRIM) {   /* the outward apron: FILL only between observed papyrus (2026-10-15) */
        CloseupFitOpts fo;
        CloseupApronStats ap;
        memset(&ap, 0, sizeof ap);
        if (c->complete_gaps) CloseupFitOpts_complete(&fo); else CloseupFitOpts_default(&fo);
        if (CloseupSurface_trim_apron(g, fo.fill_px, &ap) != 0) goto out;
        fprintf(stderr, "  grow: outward apron: %zu of %zu FILL nodes had papyrus on one side only (reach %.0f px): dropped "
                "(%.1fs)\n", ap.dropped, ap.fill_before, fo.fill_px, ap.seconds);
    }
    /* duplicated papyrus: FILL extrapolated over papyrus the surface holds elsewhere goes before the flattening */
    if (CloseupSurface_dedup_fill(g, CLOSEUP_DEDUP_VOX, CLOSEUP_DEDUP_FRAME_PX, CLOSEUP_DEDUP_COS, dd) != 0) goto out;
    fprintf(stderr, "  grow: duplicated papyrus: %zu node pairs within %.1f vox but > %.0f px apart in the frame (%zu DATA-DATA, "
            "kept); %zu FILL nodes dropped (%.1fs)\n", dd->pairs, CLOSEUP_DEDUP_VOX, CLOSEUP_DEDUP_FRAME_PX, dd->data_pairs,
            dd->dropped_fill, dd->seconds);
    uvn = (double *)malloc(nn * 2 * sizeof *uvn);
    if (uvn == NULL) goto out;
    if (c->grow_reflatten) {
        fprintf(stderr, "  grow: flatten start (%.1fs)\n", ves_clock_sec() - t0);
        if (CloseupSurface_flatten_slim(g, CLOSEUP_FLAT_PX, CLOSEUP_GROW_SLIM_ITERS, uvn, flat) == 0) {
            *reflattened = 1;
            fprintf(stderr, "  grow: flattened %zu coarse nodes (+%zu fringe), %zu components, %zu iterations; |log sigma| p99 "
                    "%.4f -> %.4f; area ratio %.4f (%.1fs)\n", flat->coarse_nodes, flat->coarse_fringe, flat->components,
                    flat->iterations, flat->logs_p99_before, flat->logs_p99_after, flat->area_ratio, flat->seconds);
        } else {
            fprintf(stderr, "ERROR: SLIM flattening failed; no replacement written\n");
            goto out;
        }
    }
    if (!*reflattened)
        for (v = 0; v < nn; v++) {
            uvn[v * 2] = g->state[v] ? g->u0 + g->h * (double)(v % g->ni) : NAN;
            uvn[v * 2 + 1] = g->state[v] ? g->v0 + g->h * (double)(v / g->ni) : NAN;
        }
    if (CloseupSurface_mesh(g, uvn, &verts, &muv, &nrm, &faces, &nodes, &nv, &nf) != 0 || nf == 0) goto out;
    if (*reflattened) {                           /* rows along the axis, as the surface stage leaves them */
        CloseupFrame fr2;
        double *maxial = (double *)malloc(nv * sizeof *maxial), umin[2] = {1e300, 1e300};
        float *minward = (float *)malloc(nv * 3 * sizeof *minward);
        if (maxial == NULL || minward == NULL) { free(maxial); free(minward); goto out; }
        ic_axis_fields(axis, verts, nv, maxial, minward);
        ic_axial_sense(c, maxial, nv);
        if (CloseupSurface_frame(verts, muv, faces, nf, NULL, maxial, minward, 0.0, 0.0, &fr2) != 0) {
            free(maxial); free(minward);
            fprintf(stderr, "ERROR: no frame on the grown surface\n");
            goto out;
        }
        free(maxial); free(minward);
        if (fr2.reflect) fprintf(stderr, "WARNING: the grown surface is mirror-imaged (vote %.3g)\n", fr2.hand_vote);
        for (v = 0; v < nv; v++) {
            double a = 0.0, b = 0.0;
            CloseupFrame_apply(&fr2, muv[v * 2], muv[v * 2 + 1], &a, &b);
            muv[v * 2] = a; muv[v * 2 + 1] = b;
            umin[0] = fmin(umin[0], a); umin[1] = fmin(umin[1], b);
        }
        for (v = 0; v < nv; v++) {
            muv[v * 2] += (double)CLOSEUP_RASTER_MARGIN_PX - floor(umin[0]);
            muv[v * 2 + 1] += (double)CLOSEUP_RASTER_MARGIN_PX - floor(umin[1]);
        }
    } else {                                     /* the seed's frame, an integer translation off the raster margin */
        double umin[2] = {1e300, 1e300};
        for (v = 0; v < nv; v++) { umin[0] = fmin(umin[0], muv[v * 2]); umin[1] = fmin(umin[1], muv[v * 2 + 1]); }
        shift = ceil((CLOSEUP_RASTER_MARGIN_PX - fmin(umin[0], umin[1])) / g->h) * g->h;
        if (shift < 0.0) shift = 0.0;
        for (v = 0; v < nv * 2; v++) muv[v] += shift;
    }
    oa = (uint8_t *)malloc(nv); oc = (int32_t *)malloc(nv * sizeof *oc);
    if (oa == NULL || oc == NULL) goto out;
    for (v = 0; v < nv; v++) {
        oa[v] = g->state[nodes[v]]; oc[v] = g->chart[nodes[v]];
        if (oa[v] == CLOSEUP_NODE_DATA && oc[v] >= adopt_first) oa[v] = 3;
    }
    snprintf(path, sizeof path, "%s/%s_surface.vmesh", c->out_dir, c->id);
    ves_ensure_parent_dir(path);
    if (MeshBin_write_uv64(path, verts, nv, faces, nf, muv) != 0) goto out;
    snprintf(path, sizeof path, "%s/%s_surface_nrm.f32", c->out_dir, c->id);
    if (ic_write_file(path, nrm, nv * 3 * sizeof *nrm) != 0) goto out;
    snprintf(path, sizeof path, "%s/%s_surface_attr.u8", c->out_dir, c->id);
    if (ic_write_file(path, oa, nv) != 0) goto out;
    snprintf(path, sizeof path, "%s/%s_surface_chart.i32", c->out_dir, c->id);
    if (ic_write_file(path, oc, nv * sizeof *oc) != 0) goto out;
    if (have_f0) {   /* provenance: 1 where the node's nearest observation is the F0 bridge's */
        of0 = (uint8_t *)malloc(nv ? nv : 1);
        if (of0 == NULL) goto out;
        for (v = 0; v < nv; v++) of0[v] = (uint8_t)(oa[v] == 3 && oc[v] == f0_chart);
        snprintf(path, sizeof path, "%s/%s_surface_f0.u8", c->out_dir, c->id);
        if (ic_write_file(path, of0, nv) != 0) goto out;
    }
    rc = 0;
out:
    free(uvn); free(verts); free(muv); free(nrm); free(faces); free(nodes); free(oa); free(oc); free(of0);
    *nv_out = nv; *nf_out = nf;
    return rc;
}

/* Grow an existing close-up (surface_source) over the target box growth_uv in its own UV frame: its observed faces
 * seed the fit; then, read after read, the pools' charts near the DATA boundary (ic_data_boundary, dilated by
 * margin_vox: CLOSEUP_GROW_*) join by closeup_adopt's 3-D contact, round by round, the surface refitted after every
 * round that placed a chart.  Nothing is joined that the adoption gates refuse (contact, rigid fit, other-wrap
 * conflict, continuity, metric): the growth stops where the sheet does, at the box, or when a read places nothing.
 * reflatten: the grown surface flattened injectively (CloseupSurface_flatten_slim: SLIM with a scaffold, its rows
 * along the axis); else the source frame kept, as a completion pass does.  The grid as the flattening sees it goes
 * to <id>_grid.cgrid first (--reflatten).  Writes the surface files and <id>_surface_report.json
 * (schema ink-closeup-grow-v1, every read's numbers). */
static int ic_stage_grow(const IcConfig *c)
{
    char path[2600], js[4800];
    MeshBinData m;
    CloseupObs obs;
    CloseupObsGate gate;
    CloseupGrid g;
    CloseupFitOpts fit;
    CloseupFitStats fs;
    CloseupFlatStats flat;
    CloseupDedupStats dedup;
    Arena_T arena = NULL;
    AsmAxis *axis = NULL;
    IcAdopted *adopted = NULL;
    IcKeySet *tiles_done = NULL;               /* grow_tile_px > 0: the adopted tiles of each pool, by key */
    CloseupF0 f0 = NULL;                       /* the F0 bridge's cube cache (grow "f0") */
    CloseupF0BridgeOpts bo;
    CloseupF0BridgeStats f0tot;
    int32_t f0_chart = -1;                     /* the bridge's observations: one id past the last pool namespace */
    size_t f0_snapped = 0, f0_passes_run = 0;
    int64_t *tkeys = NULL;                     /* the loaded pool's tile keys (ic_pool_tile) */
    uint8_t *attr = NULL, *use = NULL, *mark = NULL;
    int32_t *chart = NULL, *fc = NULL;
    float *inward = NULL;
    double *uv = NULL;
    double lo[2] = {1e300, 1e300}, hi[2] = {-1e300, -1e300}, agree = 0.0, t0 = ves_clock_sec();
    size_t got = 0, f = 0, v = 0, used = 0, gated = 0, nv = 0, nf = 0, nn = 0, obs_source = 0, nadopt = 0;
    size_t target_nodes = 0, reads_done = 0, data0 = 0, data_prev = 0;
    int prev_low = 0;
    const char *stop_reason = "max_reads";
    size_t reads_since_pool = 0, last_placed = 1;   /* adaptive pool reads (CLOSEUP_GROW_POOL_EVERY) */
    int force_pool = 1;
    char *reads_json = NULL;
    size_t reads_len = 0, reads_cap = 0;
    int32_t adopt_first = 0;
    int k = 0, orient = 0, rc = 1, reflattened = 0;
    FILE *fp = NULL;
    memset(&m, 0, sizeof m); memset(&obs, 0, sizeof obs); memset(&g, 0, sizeof g); memset(&flat, 0, sizeof flat);
    memset(&dedup, 0, sizeof dedup);
    if (!c->has_growth || c->npool < 1) {
        fprintf(stderr, "ERROR: \"grow\" needs \"growth_uv\" (the target box) and at least one \"pool\"\n");
        return 1;
    }
    if (c->complete_gaps) CloseupFitOpts_complete(&fit); else CloseupFitOpts_default(&fit);
    fit.bridge = c->bridge && c->chart_wrap != NULL; fit.chart_wrap = c->chart_wrap; fit.nchart_wrap = c->nchart_wrap;
    fit.bridge_px = fit.fill_px;
    /* the seed: the source close-up's observed faces in its own frame (old FILL never observes) */
    snprintf(path, sizeof path, "%s_surface.vmesh", c->surface_source);
    if (MeshBin_read_precise_malloc(path, &m) != 0 || (m.uv == NULL && m.uv64 == NULL)) {
        fprintf(stderr, "ERROR: cannot read the source surface %s\n", path);
        goto out;
    }
    snprintf(path, sizeof path, "%s_surface_attr.u8", c->surface_source);
    attr = (uint8_t *)ic_read_file(path, m.nv, &got);
    snprintf(path, sizeof path, "%s_surface_chart.i32", c->surface_source);
    chart = (int32_t *)ic_read_file(path, m.nv * sizeof *chart, &got);
    uv = (double *)malloc(m.nv * 2 * sizeof *uv);
    use = (uint8_t *)calloc(m.nf ? m.nf : 1, 1);
    fc = (int32_t *)malloc((m.nf ? m.nf : 1) * sizeof *fc);
    adopted = (IcAdopted *)calloc((size_t)c->npool, sizeof *adopted);
    tiles_done = (IcKeySet *)calloc((size_t)c->npool, sizeof *tiles_done);
    if (attr == NULL || chart == NULL || uv == NULL || use == NULL || fc == NULL || adopted == NULL || tiles_done == NULL)
        goto out;
    for (v = 0; v < m.nv; v++) {
        if (chart[v] >= adopt_first) {
            if (chart[v] == INT32_MAX) goto out;
            adopt_first = chart[v] + 1;
        }
        for (k = 0; k < 2; k++) {
            double a = MeshBin_uv(&m, v * 2 + (size_t)k);
            if (!isfinite(a)) goto out;
            uv[v * 2 + (size_t)k] = a;
            if (attr[v] == 1 || attr[v] == 3) { lo[k] = fmin(lo[k], a); hi[k] = fmax(hi[k], a); }
        }
    }
    if ((double)adopt_first + (double)c->npool * (double)CLOSEUP_GROW_POOL_ID_STRIDE + 1.0 > (double)INT32_MAX) goto out;
    f0_chart = adopt_first + (int32_t)c->npool * CLOSEUP_GROW_POOL_ID_STRIDE;
    memset(&f0tot, 0, sizeof f0tot);
    CloseupF0BridgeOpts_default(&bo);
    bo.depth_px = c->f0_depth_px; bo.rings = c->f0_rings;
    for (k = 0; k < 4; k++) bo.box[k] = c->growth_uv[k];
    if (c->nf0 > 0) {
        const char *dirs[16];
        for (k = 0; k < c->nf0; k++) dirs[k] = c->f0_dir[k];
        f0 = CloseupF0_open(dirs, (size_t)c->nf0, CLOSEUP_F0_MAX_CUBES);
        if (f0 == NULL) goto out;
    }
    for (f = 0; f < m.nf; f++) {
        use[f] = 1;
        fc[f] = chart[m.faces[f * 3]];
        for (k = 0; k < 3; k++) {
            uint8_t a = attr[m.faces[f * 3 + (size_t)k]];
            if (a != 1 && a != 3) use[f] = 0;
        }
    }
    CloseupObsGate_default(&gate);
    if (CloseupObs_add_faces(&obs, m.verts, uv, m.faces, m.nf, use, fc, &gate, &used, &gated) != 0 || obs.n == 0) goto out;
    obs_source = obs.n;
    MeshBin_dispose(&m);
    free(attr); attr = NULL; free(chart); chart = NULL; free(uv); uv = NULL; free(use); use = NULL; free(fc); fc = NULL;
    fprintf(stderr, "[ink_closeup] grow: %zu observations from %zu observed faces of %s (%zu gated); target box "
            "[%.0f, %.0f, %.0f, %.0f], the seed's [%.0f, %.0f, %.0f, %.0f]\n", obs.n, used, c->surface_source, gated,
            c->growth_uv[0], c->growth_uv[1], c->growth_uv[2], c->growth_uv[3], lo[0], lo[1], hi[0], hi[1]);
    lo[0] = fmin(lo[0], c->growth_uv[0]); lo[1] = fmin(lo[1], c->growth_uv[1]);
    hi[0] = fmax(hi[0], c->growth_uv[2]); hi[1] = fmax(hi[1], c->growth_uv[3]);
    if (CloseupGrid_alloc(&g, lo[0], hi[0], lo[1], hi[1], CLOSEUP_GRID_PX) != 0) goto out;
    nn = g.ni * g.nj;
    for (v = 0; v < nn; v++) {
        double u = g.u0 + g.h * (double)(v % g.ni), w = g.v0 + g.h * (double)(v / g.ni);
        target_nodes += u >= c->growth_uv[0] && u <= c->growth_uv[2] && w >= c->growth_uv[1] && w <= c->growth_uv[3];
    }
    inward = (float *)malloc(nn * 3 * sizeof *inward);
    mark = (uint8_t *)malloc(nn);
    if (inward == NULL || mark == NULL) goto out;
    if (CloseupSurface_fit(&g, &obs, &fit, &fs) != 0) goto out;
    data0 = fs.data_nodes;
    arena = Arena_new(); axis = AsmAxis_load(arena, c->axis);
    if (axis == NULL) { fprintf(stderr, "ERROR: cannot load axis table %s\n", c->axis); goto out; }
    ic_axis_fields(axis, g.pos, nn, NULL, inward);
    orient = CloseupSurface_orient(&g, inward, &agree);
    fprintf(stderr, "  grow: grid %zu x %zu nodes (%zu in the target box); seed fit %zu data / %zu fill; inward %d "
            "(agree %.4f) (%.1fs)\n", g.ni, g.nj, target_nodes, fs.data_nodes, fs.fill_nodes, orient, agree,
            ves_clock_sec() - t0);
    reads_cap = 4096;
    reads_json = (char *)malloc(reads_cap);
    if (reads_json == NULL) goto out;
    reads_json[0] = 0;
    for (reads_done = 0; reads_done < (size_t)c->grow_max_reads; reads_done++) {
        double lo3[3] = {1e300, 1e300, 1e300}, hi3[3] = {-1e300, -1e300, -1e300}, tr = ves_clock_sec();
        size_t nb = ic_data_boundary(&g, mark), placed = 0, loaded = 0, data_now = 0, q = 0, placed_or_snapped = 0;
        IcNear nreg;
        const IcNear *nuse = NULL;
        int nr = 0, rounds = 0, read_pools = 1;
        char row[512];
        if (reads_done == 0) data_prev = data0;
        if (nb == 0) { stop_reason = "no_boundary"; break; }
        {   /* a graceful stop: <out_dir>/STOP ends the growth here; the surface grown so far is written as usual */
            char stop[2600];
            FILE *sf = NULL;
            snprintf(stop, sizeof stop, "%s/STOP", c->out_dir);
            sf = fopen(stop, "rb");
            if (sf != NULL) {
                fclose(sf);
                fprintf(stderr, "  grow: %s found: the growth ends after %zu reads\n", stop, reads_done);
                stop_reason = "stop_file";
                break;
            }
        }
        {   /* adaptive: after a read that placed no chart, the next reads are F0 bridge passes only; the pools are read
             * again every CLOSEUP_GROW_POOL_EVERY reads, or once more before the growth stops */
            if (f0 != NULL && !force_pool && last_placed == 0 && reads_since_pool + 1 < (size_t)CLOSEUP_GROW_POOL_EVERY)
                read_pools = 0;
            force_pool = 0;
        }
        for (v = 0; v < nn; v++) {
            if (mark[v] == CLOSEUP_NODE_NONE) continue;
            for (q = 0; q < 3; q++) {
                lo3[q] = fmin(lo3[q], (double)g.pos[v * 3 + q]);
                hi3[q] = fmax(hi3[q], (double)g.pos[v * 3 + q]);
            }
        }
        for (q = 0; q < 3; q++) { lo3[q] -= c->grow_margin_vox; hi3[q] += c->grow_margin_vox; }
        memset(&nreg, 0, sizeof nreg);
        if (read_pools) {
            nr = ic_near_build(g.pos, mark, nn, lo3, hi3, c->grow_margin_vox, &nreg);
            if (nr < 0) goto out;
            if (nr == 0) nuse = &nreg;
        }
        for (k = 0; k < c->npool; k++) {
            /* the near region picks the charts; each is read WHOLE (no clip to a search box: an adopted chart is never
             * read again, so a clipped one would lose its far part for good) */
            const double all_lo[3] = {-INFINITY, -INFINITY, -INFINITY}, all_hi[3] = {INFINITY, INFINITY, INFINITY};
            IcPool pool;
            double *puv = NULL;
            uint8_t *done = NULL;
            int32_t base = adopt_first + (int32_t)k * CLOSEUP_GROW_POOL_ID_STRIDE;
            int round = 0;
            if (nuse == NULL) break;                         /* nothing to read near: never the whole pool */
            /* tiled: a chart with adopted tiles is read again for the rest, its adopted tiles marked done */
            if (ic_pool_load(c->pool[k].path, c->pool[k].surface, axis, all_lo, all_hi, nuse,
                             c->grow_tile_px > 0.0 ? NULL : adopted[k].on, c->grow_tile_px > 0.0 ? 0 : adopted[k].n,
                             0, &pool) != 0) {
                ic_near_free(&nreg);
                goto out;
            }
            free(tkeys); tkeys = NULL;
            if (c->grow_tile_px > 0.0 && ic_pool_tile(&pool, c->grow_tile_px, &tkeys) != 0) {
                fprintf(stderr, "ERROR: cannot tile pool %s at %.1f px\n", c->pool[k].path, c->grow_tile_px);
                ic_pool_free(&pool); ic_near_free(&nreg);
                goto out;
            }
            loaded += pool.p.nch;
            puv = (double *)malloc((pool.p.nv + 1) * 2 * sizeof *puv);
            done = (uint8_t *)calloc(pool.p.nch + 1, 1);
            if (puv == NULL || done == NULL) { free(puv); free(done); ic_pool_free(&pool); ic_near_free(&nreg); goto out; }
            if (tkeys != NULL)
                for (f = 0; f < pool.p.nch; f++) done[f] = (uint8_t)ic_keyset_has(&tiles_done[k], tkeys[f]);
            for (f = 0; f < pool.p.nf; f++) {
                if (pool.face_chart[f] < 0 || pool.face_chart[f] >= CLOSEUP_GROW_POOL_ID_STRIDE) {
                    fprintf(stderr, "ERROR: pool %s holds chart id %d (>= %d)\n", c->pool[k].path, pool.face_chart[f],
                            CLOSEUP_GROW_POOL_ID_STRIDE);
                    free(puv); free(done); ic_pool_free(&pool); ic_near_free(&nreg);
                    goto out;
                }
                pool.face_chart[f] += base;                  /* each pool its own chart namespace */
            }
            for (round = 0; round < c->grow_rounds && pool.p.nch > 0; round++) {
                CloseupPlacement *pl = NULL;
                CloseupAdoptStats ast;
                size_t npl = 0, p = 0;
                if (ic_adopt_round(c, &g, &pool, done, k, (int)reads_done * 64 + round, &pl, &npl, &ast, NULL) != 0) {
                    free(puv); free(done); ic_pool_free(&pool); ic_near_free(&nreg);
                    goto out;
                }
                for (p = 0; p < npl; p++) {
                    size_t ch = pl[p].chart, a = 0, nu = 0, ng = 0;
                    done[ch] = 1;
                    for (a = pool.vofs[ch]; a < pool.vofs[ch + 1]; a++)
                        CloseupPlacement_apply(&pl[p], pool.uv[a * 2], pool.uv[a * 2 + 1], &puv[a * 2], &puv[a * 2 + 1]);
                    if (CloseupObs_add_faces(&obs, pool.verts, puv, pool.faces + pool.fofs[ch] * 3,
                                             pool.fofs[ch + 1] - pool.fofs[ch], NULL, pool.face_chart + pool.fofs[ch],
                                             &gate, &nu, &ng) != 0 ||
                        (tkeys != NULL ? ic_keyset_add(&tiles_done[k], tkeys[ch])
                                       : ic_adopted_add(&adopted[k], pool.chart_id[ch])) != 0) {
                        free(pl); free(puv); free(done); ic_pool_free(&pool); ic_near_free(&nreg);
                        goto out;
                    }
                }
                free(pl);
                rounds++;
                placed += npl;
                fprintf(stderr, "    read %zu pool %d round %d: %zu of %zu candidates placed (%zu touch); refused rms %zu, "
                        "conflict %zu, continuity %zu, metric %zu, layered %zu; deferred %zu (layered %zu)\n", reads_done, k,
                        round, npl, ast.candidates, ast.touched, ast.refused_rms, ast.refused_conflict, ast.refused_continuity,
                        ast.refused_metric, ast.refused_layer, ast.deferred, ast.deferred_layer);
                if (npl == 0) break;
                if (CloseupSurface_fit(&g, &obs, &fit, &fs) != 0) {
                    free(puv); free(done); ic_pool_free(&pool); ic_near_free(&nreg);
                    goto out;
                }
                ic_axis_fields(axis, g.pos, nn, NULL, inward);
                orient = CloseupSurface_orient(&g, inward, &agree);
            }
            free(puv); free(done); ic_pool_free(&pool);
            free(tkeys); tkeys = NULL;
        }
        ic_near_free(&nreg);
        nadopt += placed;
        if (read_pools) {
            reads_since_pool = 0;
            last_placed = placed;
        } else {
            reads_since_pool++;
        }
        {   /* the F0 bridge: where no chart reaches across a gap, the prediction carries the data on (bounded depth) */
            size_t pass = 0, snapped_read = 0;
            for (pass = 0; f0 != NULL && pass < (size_t)c->f0_passes; pass++) {
                CloseupF0BridgeStats bs;
                double tb = ves_clock_sec();
                if (CloseupF0Bridge_pass(&g, f0, &bo, f0_chart, &obs, &bs) != 0) goto out;
                f0_passes_run++;
                f0tot.candidates += bs.candidates; f0tot.snapped += bs.snapped; f0tot.no_cube += bs.no_cube;
                f0tot.no_run += bs.no_run; f0tot.shifted += bs.shifted; f0tot.thick += bs.thick; f0tot.broken += bs.broken;
                f0tot.deep += bs.deep; f0tot.faces += bs.faces; f0tot.observations += bs.observations;
                f0tot.dup += bs.dup; f0tot.ambiguous += bs.ambiguous; f0tot.layered += bs.layered;
                fprintf(stderr, "    read %zu f0 bridge pass %zu: %zu of %zu ring nodes snapped (no cube %zu, no run %zu, shifted "
                        "%zu, thick %zu, broken %zu, duplicate %zu, ambiguous %zu, layered %zu; %zu too deep); %zu faces, %zu "
                        "observations (%.1fs)\n", reads_done, pass, bs.snapped, bs.candidates, bs.no_cube, bs.no_run, bs.shifted,
                        bs.thick, bs.broken, bs.dup, bs.ambiguous, bs.layered, bs.deep, bs.faces, bs.observations,
                        ves_clock_sec() - tb);
                if (bs.snapped == 0 || bs.observations == 0) break;
                snapped_read += bs.snapped;
                if (CloseupSurface_fit(&g, &obs, &fit, &fs) != 0) goto out;
                ic_axis_fields(axis, g.pos, nn, NULL, inward);
                orient = CloseupSurface_orient(&g, inward, &agree);
            }
            f0_snapped += snapped_read;
            placed_or_snapped = placed + snapped_read;
        }
        for (v = 0; v < nn; v++) data_now += g.state[v] == CLOSEUP_NODE_DATA;
        fprintf(stderr, "  grow read %zu: data boundary %zu nodes, %zu charts read%s, %zu placed in %d rounds; data %zu nodes "
                "(%.1f%% of the target box; seed %zu) (%.1fs, total %.1fs)\n", reads_done, nb, loaded,
                read_pools ? "" : " (pools skipped: bridge only)", placed, rounds,
                data_now, target_nodes ? 100.0 * (double)data_now / (double)target_nodes : 0.0, data0,
                ves_clock_sec() - tr, ves_clock_sec() - t0);
        snprintf(row, sizeof row, "%s\n    {\"read\": %zu, \"boundary_nodes\": %zu, \"pools_read\": %d, \"charts_read\": %zu, "
                 "\"placed\": %zu, \"rounds\": %d, \"f0_snapped\": %zu, \"data_nodes\": %zu, \"seconds\": %.1f}",
                 reads_done ? "," : "", reads_done, nb, read_pools, loaded, placed, rounds, placed_or_snapped - placed,
                 data_now, ves_clock_sec() - tr);
        if (reads_len + strlen(row) + 1 > reads_cap) {
            char *grown = (char *)realloc(reads_json, reads_cap * 2 + strlen(row));
            if (grown == NULL) goto out;
            reads_json = grown;
            reads_cap = reads_cap * 2 + strlen(row);
        }
        memcpy(reads_json + reads_len, row, strlen(row) + 1);
        reads_len += strlen(row);
        if (placed_or_snapped == 0) {
            if (!read_pools) { force_pool = 1; continue; }     /* the pools once more before the growth stops */
            reads_done++;
            stop_reason = "nothing_placed";
            break;
        }
        {   /* measured progress: two reads in a row that each grew the data by less than min_gain of it end the growth,
             * after one more read of the pools (CLOSEUP_GROW_MIN_GAIN) */
            const int low = c->grow_min_gain > 0.0 && data_prev > 0 &&
                            (double)data_now < (1.0 + c->grow_min_gain) * (double)data_prev;
            data_prev = data_now;
            if (low && prev_low) {
                if (!read_pools) { force_pool = 1; continue; }
                fprintf(stderr, "  grow: two reads in a row grew the data by less than %.1f%% each: the growth ends after %zu "
                        "reads\n", 100.0 * c->grow_min_gain, reads_done + 1);
                reads_done++;
                stop_reason = "min_gain";
                break;
            }
            prev_low = low;
        }
    }
    /* the grid as the final flattening sees it: --reflatten redoes the steps after this on it (2026-10-15) */
    snprintf(path, sizeof path, "%s/%s_grid.cgrid", c->out_dir, c->id);
    ves_ensure_parent_dir(path);
    {
        int32_t meta[3];
        meta[0] = adopt_first; meta[1] = f0_chart; meta[2] = f0 != NULL;
        if (CloseupGrid_save(path, &g, meta, 3) != 0) fprintf(stderr, "WARNING: cannot write %s\n", path);
    }
    if (ic_grow_surface_files(c, &g, axis, adopt_first, f0_chart, f0 != NULL, t0, &flat, &reflattened, &nv, &nf, &dedup) != 0)
        goto out;
    {
        size_t data_end = 0;
        for (v = 0; v < nn; v++) data_end += g.state[v] == CLOSEUP_NODE_DATA;
        snprintf(path, sizeof path, "%s/%s_surface_report.json", c->out_dir, c->id);
        fp = fopen(path, "wb"); if (fp == NULL) goto out;
        ic_json_str(js, sizeof js, c->surface_source);
        fprintf(fp, "{\n  \"schema\": \"ink-closeup-grow-v1\",\n  \"id\": \"%s\",\n  \"surface_source\": %s,\n", c->id, js);
        fprintf(fp, "  \"target_uv\": [%.1f, %.1f, %.1f, %.1f],\n  \"target_nodes\": %zu,\n", c->growth_uv[0],
                c->growth_uv[1], c->growth_uv[2], c->growth_uv[3], target_nodes);
        fprintf(fp, "  \"grow\": {\"margin_vox\": %.1f, \"rounds\": %d, \"max_reads\": %d, \"reflatten\": %s, \"reflattened\": %s, "
                "\"tile_px\": %.1f, \"tiles_adopted\": %zu, \"min_gain\": %.4f, \"stop_reason\": \"%s\"},\n",
                c->grow_margin_vox, c->grow_rounds, c->grow_max_reads, c->grow_reflatten ? "true" : "false",
                reflattened ? "true" : "false", c->grow_tile_px, ic_keysets_total(tiles_done, c->npool), c->grow_min_gain,
                stop_reason);
        fprintf(fp, "  \"gates\": {\"on_vox\": %.2f, \"foot_fill_px\": %.1f, \"min_extent_px\": %.1f, \"rms_px\": %.2f, "
                "\"conflict_vox\": %.2f, \"agree_frac\": %.3f, \"normal_deg\": %.1f, \"min_corr\": %d, "
                "\"continuity_px\": %.1f},\n", c->adopt.on_vox, c->adopt.foot_fill_px, c->adopt.min_extent_px,
                c->adopt.rms_px, c->adopt.conflict_vox, c->adopt.agree_frac, c->adopt.normal_deg, c->adopt.min_corr,
                c->adopt.continuity_px);
        fprintf(fp, "  \"pools\": [");
        for (k = 0; k < c->npool; k++) {
            ic_json_str(js, sizeof js, c->pool[k].path);
            fprintf(fp, "%s%s", k ? ", " : "", js);
        }
        fprintf(fp, "],\n  \"complete_gaps\": %s,\n", c->complete_gaps ? "true" : "false");
        fprintf(fp, "  \"seed\": {\"observations\": %zu, \"data_nodes\": %zu},\n", obs_source, data0);
        fprintf(fp, "  \"adopt\": {\"charts\": %zu, \"observations\": %zu, \"first_chart_id\": %d, \"pool_id_stride\": %d},\n",
                nadopt, obs.n - obs_source, adopt_first, CLOSEUP_GROW_POOL_ID_STRIDE);
        fprintf(fp, "  \"reads\": [%s\n  ],\n", reads_json);
        fprintf(fp, "  \"f0\": {\"dirs\": [");
        for (k = 0; k < c->nf0; k++) {
            ic_json_str(js, sizeof js, c->f0_dir[k]);
            fprintf(fp, "%s%s", k ? ", " : "", js);
        }
        fprintf(fp, "], \"chart\": %d, \"depth_px\": %.1f, \"rings\": %d, \"passes_max\": %d, \"passes_run\": %zu, "
                "\"window_vox\": %.2f, \"shift_vox\": %.2f, \"max_run_vox\": %.2f, \"candidates\": %zu, \"snapped\": %zu, "
                "\"no_cube\": %zu, \"no_run\": %zu, \"shifted\": %zu, \"thick\": %zu, \"broken\": %zu, \"duplicate\": %zu, "
                "\"ambiguous\": %zu, \"layered\": %zu, \"deep\": %zu, \"dup_vox\": %.2f, \"ambig_vox\": %.2f, \"layer_vox\": %.2f, "
                "\"faces\": %zu, \"observations\": %zu},\n",
                f0_chart, bo.depth_px, bo.rings, c->f0_passes, f0_passes_run, bo.window_vox, bo.shift_vox, bo.max_run_vox,
                f0tot.candidates, f0tot.snapped, f0tot.no_cube, f0tot.no_run, f0tot.shifted, f0tot.thick, f0tot.broken,
                f0tot.dup, f0tot.ambiguous, f0tot.layered, f0tot.deep, bo.dup_vox, bo.ambig_vox, bo.layer_vox, f0tot.faces,
                f0tot.observations);
        fprintf(fp, "  \"fit\": {\"data\": %zu, \"fill\": %zu, \"fill_rejected\": %zu, \"resid_p50\": %.6f, \"resid_p99\": %.6f, "
                "\"torn_stretch\": %zu, \"torn_float\": %zu, \"orient\": %d, \"orient_agree\": %.6f, "
                "\"target_coverage\": %.5f},\n", data_end, fs.fill_nodes, fs.fill_rejected, fs.resid_p50, fs.resid_p99,
                fs.torn_stretch, fs.torn_float, orient, agree,
                target_nodes ? (double)data_end / (double)target_nodes : 0.0);
        fprintf(fp, "  \"flatten\": {\"method\": \"%s\", \"coarse_nodes\": %zu, \"fringe_nodes\": %zu, \"components\": %zu, "
                "\"iterations\": %zu, \"logs_p99_before\": %.5f, "
                "\"logs_p99_after\": %.5f, \"area_ratio\": %.5f, \"seconds\": %.2f},\n",
                reflattened ? IC_GROW_FLATTEN_METHOD : "frame", flat.coarse_nodes, flat.coarse_fringe, flat.components,
                flat.iterations, flat.logs_p99_before, flat.logs_p99_after, flat.area_ratio, flat.seconds);
        fprintf(fp, "  \"dedup\": {\"vox\": %.2f, \"frame_px\": %.1f, \"pairs\": %zu, \"data_pairs\": %zu, \"dropped_fill\": %zu},\n",
                CLOSEUP_DEDUP_VOX, CLOSEUP_DEDUP_FRAME_PX, dedup.pairs, dedup.data_pairs, dedup.dropped_fill);
        fprintf(fp, "  \"surface\": {\"verts\": %zu, \"faces\": %zu},\n  \"seconds\": %.2f\n}\n", nv, nf, ves_clock_sec() - t0);
        fclose(fp); fp = NULL;
        fprintf(stderr, "[ink_closeup] SUMMARY grow id=%s reads=%zu adopted=%zu f0_snapped=%zu data %zu -> %zu nodes (%.1f%% of "
                "the target box) verts=%zu faces=%zu (%.1fs)\n", c->id, reads_done, nadopt, f0_snapped, data0, data_end,
                target_nodes ? 100.0 * (double)data_end / (double)target_nodes : 0.0, nv, nf, ves_clock_sec() - t0);
    }
    rc = 0;
out:
    if (fp != NULL) fclose(fp);
    if (rc) fprintf(stderr, "ERROR: growth failed for %s\n", c->surface_source);
    MeshBin_dispose(&m); CloseupObs_free(&obs); CloseupGrid_free(&g);
    if (arena != NULL) Arena_dispose(&arena);
    if (adopted != NULL) for (k = 0; k < c->npool; k++) free(adopted[k].on);
    if (tiles_done != NULL) for (k = 0; k < c->npool; k++) free(tiles_done[k].slot);
    free(tiles_done); free(tkeys);
    CloseupF0_close(&f0);
    free(adopted); free(attr); free(use); free(chart); free(fc);
    free(inward); free(uv); free(mark); free(reads_json);
    return rc;
}

/* --reflatten (2026-10-15): a finished growth's last steps again, from its grid (<source id>_grid.cgrid, written before
 * the final flattening), with this build's flattening, into the config's out_dir and id: the surface files, and the
 * source's surface report with this flattening's "flatten", "surface" and "seconds", "id" and a "reflatten_of" line. */
static int ic_stage_reflatten(const IcConfig *c0, const char *grid_path)
{
    char path[2600], src[2600], js[4800];
    const char *suffix = "_grid.cgrid";
    IcConfig c = *c0;
    CloseupGrid g;
    CloseupFlatStats flat;
    CloseupDedupStats dedup;
    Arena_T arena = NULL;
    AsmAxis *axis = NULL;
    int32_t meta[3] = { 0, -1, 0 };
    size_t nm = 0, nv = 0, nf = 0, len = 0, got = 0;
    char *rep = NULL, *p = NULL;
    int reflattened = 0, rc = 1;
    double t0 = ves_clock_sec();
    FILE *fp = NULL;
    memset(&g, 0, sizeof g);
    memset(&dedup, 0, sizeof dedup);
    len = strlen(grid_path);
    if (len <= strlen(suffix) || strcmp(grid_path + len - strlen(suffix), suffix) != 0 || len - strlen(suffix) + 32 > sizeof src) {
        fprintf(stderr, "ERROR: --reflatten wants <source id>_grid.cgrid, not %s\n", grid_path);
        return 1;
    }
    snprintf(src, sizeof src, "%.*s_surface_report.json", (int)(len - strlen(suffix)), grid_path);
    rep = (char *)ic_read_file(src, 0, &got);
    if (rep == NULL || got == 0) { fprintf(stderr, "ERROR: cannot read the source report %s\n", src); goto out; }
    {
        char *t = (char *)realloc(rep, got + 1);          /* NUL-terminated: read line by line */
        if (t == NULL) goto out;
        rep = t;
        rep[got] = '\0';
    }
    if (CloseupGrid_load(grid_path, &g, meta, 3, &nm) != 0 || nm != 3) {
        fprintf(stderr, "ERROR: cannot read the grid %s\n", grid_path);
        goto out;
    }
    arena = Arena_new();
    axis = arena != NULL ? AsmAxis_load(arena, c.axis) : NULL;
    if (axis == NULL) { fprintf(stderr, "ERROR: cannot read the axis %s\n", c.axis); goto out; }
    c.grow_reflatten = 1;
    fprintf(stderr, "[ink_closeup] reflatten %s: grid %zu x %zu -> %s/%s\n", grid_path, g.ni, g.nj, c.out_dir, c.id);
    if (ic_grow_surface_files(&c, &g, axis, meta[0], meta[1], meta[2] != 0, t0, &flat, &reflattened, &nv, &nf, &dedup) != 0)
        goto out;
    snprintf(path, sizeof path, "%s/%s_surface_report.json", c.out_dir, c.id);
    fp = fopen(path, "wb");
    if (fp == NULL) goto out;
    for (p = rep; p < rep + got && *p;) {
        char *e = strchr(p, '\n');
        size_t n = e != NULL ? (size_t)(e - p) + 1 : strlen(p);
        if (strncmp(p, "  \"id\": ", 8) == 0) fprintf(fp, "  \"id\": \"%s\",\n", c.id);
        else if (strncmp(p, "  \"flatten\": ", 13) == 0) {
            fprintf(fp, "  \"flatten\": {\"method\": \"%s\", \"coarse_nodes\": %zu, \"fringe_nodes\": %zu, \"components\": %zu, "
                    "\"iterations\": %zu, \"logs_p99_before\": %.5f, \"logs_p99_after\": %.5f, \"area_ratio\": %.5f, "
                    "\"seconds\": %.2f},\n",
                    reflattened ? IC_GROW_FLATTEN_METHOD : "frame", flat.coarse_nodes, flat.coarse_fringe, flat.components,
                    flat.iterations,
                    flat.logs_p99_before, flat.logs_p99_after, flat.area_ratio, flat.seconds);
            fprintf(fp, "  \"dedup\": {\"vox\": %.2f, \"frame_px\": %.1f, \"pairs\": %zu, \"data_pairs\": %zu, "
                    "\"dropped_fill\": %zu},\n", CLOSEUP_DEDUP_VOX, CLOSEUP_DEDUP_FRAME_PX, dedup.pairs, dedup.data_pairs,
                    dedup.dropped_fill);
        } else if (strncmp(p, "  \"dedup\": ", 11) == 0) {
            /* the source growth's: replaced by this run's, written with the flatten line */
        } else if (strncmp(p, "  \"surface\": ", 13) == 0) fprintf(fp, "  \"surface\": {\"verts\": %zu, \"faces\": %zu},\n", nv, nf);
        else if (strncmp(p, "  \"seconds\": ", 13) == 0) fprintf(fp, "  \"seconds\": %.2f\n", ves_clock_sec() - t0);
        else {
            const char *rf = reflattened ? "\"reflattened\": true" : "\"reflattened\": false";
            char *q = strstr(p, "\"reflattened\": ");
            if (q != NULL && q < p + n && strncmp(p, "  \"grow\": ", 10) == 0) {
                const char *after = strchr(q, ',');
                fprintf(fp, "%.*s%s%.*s", (int)(q - p), p, rf, after != NULL ? (int)(p + n - after) : 0, after != NULL ? after : "");
            } else {
                fwrite(p, 1, n, fp);
            }
            if (strncmp(p, "  \"schema\": ", 12) == 0) {
                ic_json_str(js, sizeof js, grid_path);
                fprintf(fp, "  \"reflatten_of\": %s,\n", js);
            }
        }
        p += n;
    }
    if (fclose(fp) != 0) { fp = NULL; goto out; }
    fp = NULL;
    fprintf(stderr, "[ink_closeup] SUMMARY reflatten id=%s %s: %zu coarse nodes, %zu iterations; |log sigma| p99 %.4f -> %.4f; "
            "verts=%zu faces=%zu (%.1fs)\n", c.id, reflattened ? IC_GROW_FLATTEN_METHOD : "frame (the flattening failed)",
            flat.coarse_nodes, flat.iterations, flat.logs_p99_before, flat.logs_p99_after, nv, nf, ves_clock_sec() - t0);
    rc = reflattened ? 0 : 1;
out:
    if (fp != NULL) fclose(fp);
    if (arena != NULL) Arena_dispose(&arena);
    CloseupGrid_free(&g);
    free(rep);
    return rc;
}

/* --reuv (2026-10-15; the user, on the grown sheets' bridged lobes: "align themselves correctly with fibers after
 * landing so that they're *straight*, using the fiber path; and then ... expand their content once they're aligned
 * correctly everywhere with the rest of the sheet using the same strategy as we use for close-up"): a grown surface
 * (surface_source) carried into the fibre-aligned frame of a fiber_flatten run on its export.  Every vertex's canvas
 * position (its UV less the export's uv_origin, times its scale) goes through the forward map (fiber_map_u/v.tif at the
 * node grid of fiber_flatten.json): bilinear over the solved nodes of its cell, else the displacement of the nearest
 * solved node within IC_REUV_NEAR_CELLS, else the map's global affine; new UV = (u' - u_min, v' - v_min), the new page's
 * pixels.  3-D, faces and labels unchanged: <out_dir>/<id>_surface{.vmesh, _attr.u8, _chart.i32, _nrm.f32, _f0.u8}
 * + <id>_surface_report.json -- the seed of a completion pass (the close-up strategy) in that frame. */
#define IC_REUV_NEAR_CELLS 16

typedef struct IcReuvGuardStats {
    size_t unsafe_paths, limited_vertices;
    int passes;
    double mean_step, max_correction;
} IcReuvGuardStats;

/* Keep the full fibre proposal wherever possible. As in SLIM's local step,
 * damp only vertices incident on an unsafe triangle, then certify the whole
 * nonuniform path again. The affine fit supplies an orientation-preserving
 * starting frame; the original material vertices and faces never change. */
static int ic_reuv_guard(Arena_T arena, const double *safe, double *target,
                         size_t nv, const int32_t *faces, size_t nf,
                         IcReuvGuardStats *st)
{
    double *alpha = (double *)ARENA_ALLOC(arena, nv * sizeof *alpha);
    double *trial = (double *)ARENA_ALLOC(arena, nv * 2 * sizeof *trial);
    uint8_t *bad_face = (uint8_t *)ARENA_ALLOC(arena, nf);
    uint8_t *swap_face = (uint8_t *)ARENA_ALLOC(arena, nf);
    uint8_t *damp = (uint8_t *)ARENA_ALLOC(arena, nv);
    memset(st, 0, sizeof *st);
    for (size_t v = 0; v < nv; v++) alpha[v] = 1.0;
    for (size_t f = 0; f < nf; f++) {
        int32_t rev[3] = { faces[3*f], faces[3*f+2], faces[3*f+1] };
        swap_face[f] = (uint8_t)!UvGuard_positive(safe, faces+3*f);
        if (swap_face[f] && !UvGuard_positive(safe, rev)) return -1;
    }
    for (int pass = 0; pass < CLOSEUP_REUV_GUARD_PASSES; pass++) {
        size_t bad = 0;
        int64_t fi = 0;
        for (size_t v = 0; v < nv; v++)
            for (size_t k = 0; k < 2; k++)
                trial[2*v+k] = safe[2*v+k] + alpha[v] * (target[2*v+k]-safe[2*v+k]);
        #pragma omp parallel for schedule(static) reduction(+:bad) if(nf > 4096)
        for (fi = 0; fi < (int64_t)nf; fi++) {
            size_t f = (size_t)fi;
            int32_t idx[3] = { faces[3*f], faces[3*f+1], faces[3*f+2] };
            if (swap_face[f]) { int32_t t = idx[1]; idx[1] = idx[2]; idx[2] = t; }
            bad_face[f] = (uint8_t)!UvGuard_interval(safe, trial, idx, 1.0);
            bad += bad_face[f];
        }
        if (pass == 0) st->unsafe_paths = bad;
        st->passes = pass+1;
        if (bad == 0) {
            for (size_t v = 0; v < nv; v++) {
                double dx = trial[2*v]-target[2*v], dy = trial[2*v+1]-target[2*v+1];
                double movement = hypot(dx,dy);
                st->limited_vertices += alpha[v] < 1.0;
                st->mean_step += alpha[v];
                if (movement > st->max_correction) st->max_correction = movement;
            }
            if (nv) st->mean_step /= (double)nv;
            memcpy(target,trial,nv*2*sizeof *target);
            return 0;
        }
        if (pass == 0 || (pass+1)%10 == 0) {
            fprintf(stderr,"[orch] status Fibre mapping guard pass %d: %zu unsafe triangle paths\n",pass+1,bad);
            fflush(stderr);
        }
        memset(damp,0,nv);
        for (size_t f = 0; f < nf; f++) if (bad_face[f])
            for (size_t k = 0; k < 3; k++) damp[(size_t)faces[3*f+k]] = 1;
        for (size_t v = 0; v < nv; v++) if (damp[v]) alpha[v] *= 0.5;
    }
    return -1;
}

static int ic_copy_file(const char *from, const char *to)
{
    size_t got = 0;
    void *buf = ic_read_file(from, 0, &got);
    int rc = -1;
    if (buf == NULL) return -1;
    rc = ic_write_file(to, buf, got);
    free(buf);
    return rc;
}

static int ic_stage_reuv(const IcConfig *c, const char *fdir)
{
    char path[2600], js[4800];
    Arena_T a = Arena_new();
    const char *err = NULL, *seg = NULL;
    const JsonValue *fj = NULL, *gr = NULL, *np = NULL, *mj = NULL, *uo = NULL, *sc = NULL, *mp = NULL, *af = NULL;
    MeshBinData m;
    float *MU = NULL, *MV = NULL;
    double *uv = NULL, *canvas = NULL, aff[6] = { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0 }, u0 = 0.0, v0 = 0.0, su = 1.0, sv = 1.0;
    double row0 = 0.0, col0 = 0.0, cell = 8.0, umin = 0.0, vmin = 0.0, t0 = ves_clock_sec();
    int gw = 0, gh = 0, mw = 0, mh = 0, k = 0, rc = 1;
    size_t v = 0, n_bil = 0, n_near = 0, n_aff = 0, f = 0, pos = 0, neg = 0;
    FILE *fp = NULL;
    FiberMapStats guard;
    FiberMapResult guarded_map;
    memset(&m, 0, sizeof m);
    if (a == NULL) return 1;
    snprintf(path, sizeof path, "%s/fiber_flatten.json", fdir);
    fj = Json_parse_file(a, path, &err);
    if (fj == NULL) { fprintf(stderr, "ERROR: cannot parse %s (%s)\n", path, err ? err : "?"); goto out; }
    seg = Json_as_string(Json_object_get(fj, "seg"));
    gr = Json_object_get(fj, "grid");
    np = Json_object_get(fj, "new_page");
    mp = Json_object_get(fj, "map");
    if (seg == NULL || gr == NULL || np == NULL) { fprintf(stderr, "ERROR: %s lacks seg/grid/new_page\n", path); goto out; }
    gh = (int)Json_member_long(gr, "rows", 0); gw = (int)Json_member_long(gr, "cols", 0);
    cell = Json_member_double(gr, "cell", 8.0); row0 = Json_member_double(gr, "row0", 0.0); col0 = Json_member_double(gr, "col0", 0.0);
    umin = Json_member_double(np, "u_min", 0.0); vmin = Json_member_double(np, "v_min", 0.0);
    af = mp != NULL ? Json_object_get(mp, "affine") : NULL;
    if (af != NULL && Json_array_len(af) == 2)
        for (k = 0; k < 3; k++) {
            aff[k] = Json_as_double(Json_array_get(Json_array_get(af, 0), (size_t)k), aff[k]);
            aff[3 + k] = Json_as_double(Json_array_get(Json_array_get(af, 1), (size_t)k), aff[3 + k]);
        }
    snprintf(path, sizeof path, "%s/meta.json", seg);
    mj = Json_parse_file(a, path, &err);
    uo = mj != NULL ? Json_object_get(mj, "uv_origin") : NULL;
    sc = mj != NULL ? Json_object_get(mj, "scale") : NULL;
    if (uo == NULL || Json_array_len(uo) != 2) { fprintf(stderr, "ERROR: %s has no uv_origin\n", path); goto out; }
    u0 = Json_as_double(Json_array_get(uo, 0), 0.0); v0 = Json_as_double(Json_array_get(uo, 1), 0.0);
    if (sc != NULL && Json_array_len(sc) == 2) { su = Json_as_double(Json_array_get(sc, 0), 1.0); sv = Json_as_double(Json_array_get(sc, 1), 1.0); }
    snprintf(path, sizeof path, "%s/fiber_map_u.tif", fdir);
    if (TiffIO_load_float2d(a, path, &MU, &mw, &mh) != 0 || mw != gw || mh != gh) { fprintf(stderr, "ERROR: %s\n", path); goto out; }
    snprintf(path, sizeof path, "%s/fiber_map_v.tif", fdir);
    if (TiffIO_load_float2d(a, path, &MV, &mw, &mh) != 0 || mw != gw || mh != gh) { fprintf(stderr, "ERROR: %s\n", path); goto out; }
    snprintf(path, sizeof path, "%s_surface.vmesh", c->surface_source);
    if (MeshBin_read_precise_malloc(path, &m) != 0 || (m.uv == NULL && m.uv64 == NULL)) {
        fprintf(stderr, "ERROR: cannot read %s\n", path);
        goto out;
    }
    uv = (double *)malloc(m.nv * 2 * sizeof *uv);
    if (uv == NULL) goto out;
    canvas = (double *)ARENA_ALLOC(a,m.nv*2*sizeof *canvas);
    for (v = 0; v < m.nv; v++) {
        canvas[2*v] = (MeshBin_uv(&m,2*v)-u0)*su;
        canvas[2*v+1] = (MeshBin_uv(&m,2*v+1)-v0)*sv;
    }
    {
        FiberMapGrid grid = {0};
        FiberMapDirections directions = {0};
        grid.rows=gh; grid.cols=gw; grid.row0=row0; grid.col0=col0; grid.cell=cell;
        grid.u=MU; grid.v=MV; memcpy(grid.affine,aff,sizeof aff);
        grid.source_scale[0]=su; grid.source_scale[1]=sv;
        snprintf(path,sizeof path,"%s/fiber_directions.csv",fdir);
        FILE *probe=fopen(path,"rb");
        if(probe){
            fclose(probe);
            if(FiberMap_read_directions(a,path,&directions)){
                fprintf(stderr,"ERROR: invalid measured fibre directions: %s\n",path);goto out;
            }
            grid.directions=&directions;
        }else if(errno!=ENOENT){fprintf(stderr,"ERROR: cannot read fibre directions: %s\n",path);goto out;}
        if (FiberMap_fit(a,&grid,m.verts,canvas,m.nv,m.faces,m.nf,uv,&guarded_map,&guard) != 0) {
            fprintf(stderr,"ERROR: continuous fibre map could not preserve triangle orientation and physical metric; no replacement written\n");
            goto out;
        }
    }
    n_bil=m.nv;
    for (v = 0; v < m.nv; v++) { uv[2*v]-=umin; uv[2*v+1]-=vmin; }
    for (f = 0; f < m.nf; f++) {
        const double *A = &uv[(size_t)m.faces[f * 3] * 2], *B = &uv[(size_t)m.faces[f * 3 + 1] * 2];
        const double *Q = &uv[(size_t)m.faces[f * 3 + 2] * 2];
        const double ar = (B[0] - A[0]) * (Q[1] - A[1]) - (B[1] - A[1]) * (Q[0] - A[0]);
        pos += ar > 0.0; neg += ar < 0.0;
    }
    snprintf(path, sizeof path, "%s/%s_surface.vmesh", c->out_dir, c->id);
    ves_ensure_parent_dir(path);
    if (MeshBin_write_uv64(path, m.verts, m.nv, m.faces, m.nf, uv) != 0) { fprintf(stderr, "ERROR: cannot write %s\n", path); goto out; }
    {
        static const char *const side[4] = { "_surface_attr.u8", "_surface_chart.i32", "_surface_nrm.f32", "_surface_f0.u8" };
        char from[2600];
        for (k = 0; k < 4; k++) {
            FILE *probe = NULL;
            snprintf(from, sizeof from, "%s%s", c->surface_source, side[k]);
            probe = fopen(from, "rb");
            if (probe == NULL) {
                if (k < 2) { fprintf(stderr, "ERROR: no %s\n", from); goto out; }
                continue;
            }
            fclose(probe);
            snprintf(path, sizeof path, "%s/%s%s", c->out_dir, c->id, side[k]);
            if (ic_copy_file(from, path) != 0) { fprintf(stderr, "ERROR: cannot copy %s\n", from); goto out; }
        }
    }
    snprintf(path, sizeof path, "%s/%s_fiber_map.f64", c->out_dir, c->id);
    if (ic_write_file(path,guarded_map.uv,(size_t)guarded_map.rows*(size_t)guarded_map.cols*2*sizeof(double))) goto out;
    snprintf(path, sizeof path, "%s/%s_surface_report.json", c->out_dir, c->id);
    fp = fopen(path, "wb");
    if (fp == NULL) goto out;
    ic_json_str(js, sizeof js, c->surface_source);
    fprintf(fp, "{\n  \"schema\": \"ink-closeup-reuv-v1\",\n  \"id\": \"%s\",\n  \"surface_source\": %s,\n", c->id, js);
    ic_json_str(js, sizeof js, fdir);
    fprintf(fp, "  \"fiber_flatten\": %s,\n", js);
    fprintf(fp, "  \"vertices\": %zu,\n  \"mapped\": {\"bilinear\": %zu, \"nearest_node\": %zu, \"affine\": %zu},\n", m.nv, n_bil,
            n_near, n_aff);
    fprintf(fp,"  \"orientation_guard\": {\"method\": \"%s\", \"unsafe_paths_before\": %zu, "
               "\"unsafe_cells_before\": %zu, \"source_nodes\": %zu, \"extended_nodes\": %zu, \"relaxed_nodes\": %zu, "
               "\"baseline_pinned_nodes\": %zu, \"passes\": %d, \"node_correction_rms_px\": %.9f, \"node_correction_max_px\": %.9f, \"whole_path_certified\": true},\n",
            guard.baseline_pinned_nodes?"continuous_pinned_harmonic_map":guard.metric_iterations?"continuous_metric_constrained_map":"continuous_harmonic_map",
            guard.unsafe_faces_before,guard.unsafe_cells_before,guard.source_nodes,guard.extended_nodes,guard.relaxed_nodes,
            guard.baseline_pinned_nodes,guard.passes,guard.node_correction_rms,guard.node_correction_max);
    fprintf(fp,"  \"metric_guard\": {\"reference\": \"original_world_triangles\", \"source_bad_faces\": %zu, "
               "\"affine_bad_faces\": %zu, \"rigid_baseline\": %s, \"failed_faces_before\": %zu, "
               "\"failed_faces_after\": %zu, \"endpoint_certified\": true},\n",
            guard.source_metric_bad,guard.affine_metric_bad,guard.rigid_baseline?"true":"false",
            guard.metric_faces_before,guard.metric_faces_after);
    fprintf(fp,"  \"metric_fit\": {\"iterations\": %d, \"target_error_before\": %.17g, \"target_error_after\": %.17g, \"translation_removed\": true},\n",
            guard.metric_iterations,guard.metric_fit_before,guard.metric_fit_after);
    fprintf(fp,"  \"direction_fit\": {\"enabled\": %s, \"original_samples\": %zu, \"mean_angle_squared_before_rad2\": %.17g, "
               "\"mean_angle_squared_after_rad2\": %.17g, \"measure\": \"physical_area_evidence_weighted_both_arms_all_incident_cells\", \"assignment\": \"unordered_pair_quarter_turn\"},\n",
            guard.direction_fit?"true":"false",guard.direction_samples,guard.direction_angle_before,guard.direction_angle_after);
    fprintf(fp,"  \"guarded_map\": {\"rows\": %d, \"cols\": %d, \"row0\": %.17g, \"col0\": %.17g, "
               "\"cell\": %.17g, \"output_uv_offset\": [%.17g, %.17g], \"storage\": \"interleaved_float64_uv\"},\n",
            guarded_map.rows,guarded_map.cols,guarded_map.row0,guarded_map.col0,guarded_map.cell,umin,vmin);
    fprintf(fp, "  \"faces\": %zu,\n  \"faces_flipped\": %zu,\n  \"seconds\": %.2f\n}\n", m.nf, pos < neg ? pos : neg,
            ves_clock_sec() - t0);
    if (fclose(fp) != 0) { fp = NULL; goto out; }
    fp = NULL;
    fprintf(stderr, "[ink_closeup] SUMMARY reuv id=%s: %zu vertices (%zu bilinear, %zu nearest node, %zu affine); %zu of %zu faces "
            "flipped (%.1fs)\n", c->id, m.nv, n_bil, n_near, n_aff, pos < neg ? pos : neg, m.nf, ves_clock_sec() - t0);
    rc = 0;
out:
    if (fp != NULL) fclose(fp);
    MeshBin_dispose(&m);
    free(uv);
    Arena_dispose(&a);
    return rc;
}

/* ------------------------------------------------------------------ render */

/* The close-up grid is already in physical pixels. Use its checked area as the bake budget:
 * the generic review-bake default (2^28 pixels) rejects a whole strip such as PHerc1203 line 7.
 * Never shrink the image to fit that default; layers, geometry and CT must have the same scale. */
static int ic_render_plan(RawtexUv uv, size_t nv, RawtexWindow *win, size_t *W, size_t *H)
{
    double umax = 0.0, vmax = 0.0;
    size_t v = 0;
    RawtexPlan plan;
    for (v = 0; v < nv; v++) {
        double u = Rawtex_uv(uv, v * 2), w = Rawtex_uv(uv, v * 2 + 1);
        if (!isfinite(u) || !isfinite(w)) return 1;
        if (u > umax) umax = u;
        if (w > vmax) vmax = w;
    }
    if (umax > INT_MAX - CLOSEUP_RASTER_MARGIN_PX || vmax > INT_MAX - CLOSEUP_RASTER_MARGIN_PX) return 1;
    *W = (size_t)ceil(umax) + (size_t)CLOSEUP_RASTER_MARGIN_PX;
    *H = (size_t)ceil(vmax) + (size_t)CLOSEUP_RASTER_MARGIN_PX;
    if (*W == 0 || *H == 0 || *W > SIZE_MAX / *H / sizeof(double)) return 1;
    win->umin = 0.0; win->vmin = 0.0; win->umax = (double)*W; win->vmax = (double)*H;
    return Rawtex_plan_field(uv, nv, 1.0, 1.0, *W * *H, win, &plan) != 0 || !plan.ok;
}

/* Render the surface <surf>_surface.vmesh as <stem>_layers.zarr (depth layers, 1 vox apart, inward) + mask;
 * full: also the geometry volume, the CT bake, the provenance map and the report. */
static int ic_stage_render(const IcConfig *c, const char *surf, const char *stem, int depth, int full)
{
    char path[2600], zarr[2600], mpath[2600], gpath[2600], attrs[4096];
    MeshBinData m;
    Arena_T arena = NULL;
    CubeTable ct;
    RawtexLayersOpts opt;
    RawtexLayersStats ls;
    RawtexHandedness hand;
    RawtexWindow win;
    float *nrm = NULL;
    uint8_t *attr = NULL, *prov = NULL, *chan = NULL, *vcol = NULL, *has = NULL;
    int32_t *vchart = NULL;
    double *val = NULL, lo = 0.0, hi = 255.0, reach = 0.0, t0 = ves_clock_sec(), bfill = 0.0;
    size_t bmulti = 0;
    size_t got = 0, v = 0, W = 0, H = 0, i = 0, team = 0, mirror = 0, bw = 0, bh = 0, px = 0, pfill = 0, pdata = 0;
    size_t padopt = 0;
    int rc = 1;
    FILE *fp = NULL;

    memset(&m, 0, sizeof m);
    memset(&hand, 0, sizeof hand);
    snprintf(path, sizeof path, "%s/%s_surface.vmesh", c->out_dir, surf);
    if (MeshBin_read_precise_malloc(path, &m) != 0 || m.uv64 == NULL) {
        fprintf(stderr, "ERROR: run the surface stage first (%s)\n", path);
        return 1;
    }
    snprintf(path, sizeof path, "%s/%s_surface_nrm.f32", c->out_dir, surf);
    nrm = (float *)ic_read_file(path, m.nv * 3 * sizeof(float), &got);
    snprintf(path, sizeof path, "%s/%s_surface_attr.u8", c->out_dir, surf);
    attr = (uint8_t *)ic_read_file(path, m.nv, &got);
    snprintf(path, sizeof path, "%s/%s_surface_chart.i32", c->out_dir, surf);
    vchart = (int32_t *)ic_read_file(path, m.nv * sizeof(int32_t), &got);
    if (nrm == NULL || attr == NULL || vchart == NULL) { fprintf(stderr, "ERROR: surface companions missing\n"); goto out; }
    if (ic_render_plan((RawtexUv){ NULL, m.uv64 }, m.nv, &win, &W, &H) != 0) {
        fprintf(stderr, "ERROR: close-up raster dimensions cannot be represented at one pixel per voxel\n");
        goto out;
    }
    arena = Arena_new();
    reach = 0.5 * (double)(depth - 1) * SHEET_LAYERS_STEP_VOX + 2.0 + 1.6 * SHEET_LAYERS_APRON_PX + 2.0;
    if (cubetable_init(&ct, arena, c->raw, SHEET_LAYERS_CHUNK, m.verts, m.nv, reach) != 0) {
        fprintf(stderr, "ERROR: cube table init failed for %s\n", c->raw);
        goto out;
    }
    fprintf(stderr,"[orch] stage raw: Loading CT chunks\n");fflush(stderr);
    cubetable_prewarm_all(&ct);
    fprintf(stderr, "[ink_closeup] %s render: %zu x %zu px, %zu verts %zu faces; RAW %d chunks loaded, %d missing (%.1fs)\n",
            stem, W, H, m.nv, m.nf, ct.n_loaded, ct.n_missing, ves_clock_sec() - t0);

    memset(&opt, 0, sizeof opt);
    opt.depth = depth; opt.step = SHEET_LAYERS_STEP_VOX; opt.offset = SHEET_LAYERS_OFFSET_VOX;
    opt.tile_cols = SHEET_LAYERS_TILE_COLS; opt.chunk = SHEET_LAYERS_CHUNK;
    opt.stretch_ratio = SHEET_LAYERS_STRETCH_RATIO; opt.stretch_floor = SHEET_LAYERS_STRETCH_FLOOR;
    opt.apron_px = SHEET_LAYERS_APRON_PX; opt.apron_fit_px = SHEET_LAYERS_APRON_FIT_PX;
    opt.apron_cell_px = SHEET_LAYERS_APRON_CELL_PX; opt.apron_step_px = SHEET_LAYERS_APRON_STEP_PX;
    opt.window = &win; opt.hand = &hand;
    snprintf(zarr, sizeof zarr, "%s/%s_layers.zarr", c->out_dir, stem);
    snprintf(mpath, sizeof mpath, "%s/%s_layers_mask.tif", c->out_dir, stem);
    snprintf(gpath, sizeof gpath, "%s/%s_geometry.zarr", c->out_dir, stem);
    opt.geom_dir = full ? gpath : NULL;
    snprintf(attrs, sizeof attrs, "{\"source\": \"ink_closeup render\", \"id\": \"%s\", \"depth\": %d, \"step_vox\": %.3f, "
             "\"layer_index\": \"inward\", \"voxel_um\": %.4f}", stem, depth, SHEET_LAYERS_STEP_VOX,
             c->voxel_um);
    fprintf(stderr,"[orch] stage layers: Rendering depth layers\n");fflush(stderr);
    if (RawtexLayers_write(zarr, mpath, &ct, m.verts, (RawtexUv){ NULL, m.uv64 }, m.nv, m.faces, m.nf, nrm, NULL, &opt,
                           attrs, &ls) != 0) {
        fprintf(stderr, "ERROR: layer render failed (an existing %s is refused)\n", zarr);
        goto out;
    }
    for (i = 0; i < hand.npieces; i++) { team += hand.pieces[i].px_team; mirror += hand.pieces[i].px_mirror; }
    fprintf(stderr, "  layers: %zu painted px (%zu complete), %zu apron px (%zu rejected); handedness %zu team / %zu "
            "mirror px (%.1fs)\n", ls.painted_px, ls.complete_px, ls.apron_px, ls.apron_rejected, team, mirror, ls.seconds);

    if (!full) {
        fprintf(stderr, "[ink_closeup] SUMMARY render id=%s %zux%zu depth=%d painted=%zu (%.1fs)\n", stem, W, H, depth,
                ls.painted_px, ves_clock_sec() - t0);
        rc = 0;
        goto out;
    }
    /* the CT bake: the review bake's normal-max composite and contrast window */
    val = (double *)malloc(m.nv * sizeof *val);
    has = (uint8_t *)calloc(m.nv, 1);
    if (val == NULL || has == NULL) goto out;
    {
        long k = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (k = 0; k < (long)m.nv; k++) {
            double s = sample_vertex(&ct, &m.verts[(size_t)k * 3], &nrm[(size_t)k * 3], 4.0, 5);
            val[k] = s >= 0.0 ? s : 0.0;
            has[k] = s >= 0.0;
        }
    }
    Rawtex_stretch_window(val, has, m.nv, 1.0, 99.0, &lo, &hi);
    snprintf(path, sizeof path, "%s/%s_ct.tif", c->out_dir, stem);
    if (Rawtex_write_ct_tiled(path, &ct, m.verts, (RawtexUv){ NULL, m.uv64 }, m.nv, m.faces, m.nf, nrm,
                             lo, hi, &win, W*H, 0, &bw, &bh, &bfill, &bmulti) != 0 || bw != W || bh != H) {
        fprintf(stderr, "ERROR: CT bake failed (%zu x %zu vs %zu x %zu)\n", bw, bh, W, H);
        goto out;
    }
    fprintf(stderr, "  CT bake: window [%.0f, %.0f] -> %s\n", lo, hi, path);

    /* provenance: data pixels coloured by source chart, fill blue */
    fprintf(stderr,"[orch] stage provenance: Preserving source provenance\n");fflush(stderr);
    prov = (uint8_t *)calloc(W * H * 3, 1);
    chan = (uint8_t *)calloc(W * H, 1);
    vcol = (uint8_t *)malloc(m.nv);
    if (prov == NULL || chan == NULL || vcol == NULL) goto out;
    for (i = 0; i < 3; i++) {
        memset(chan, 0, W * H);
        for (v = 0; v < m.nv; v++) {
            uint32_t hsh = ic_hash((uint32_t)vchart[v] + 1u);
            if (attr[v] == CLOSEUP_NODE_FILL) vcol[v] = (uint8_t)(i == 0 ? 30 : i == 1 ? 60 : 200);
            else if (attr[v] == 3) vcol[v] = (uint8_t)(i == 0 ? 200 + (hsh & 55) : i == 1 ? 90 + ((hsh >> 8) & 60) : 20);
            else vcol[v] = (uint8_t)(90 + ((hsh >> (8 * i)) & 127));
        }
        CloseupSurface_raster_attr(m.uv64, m.faces, m.nf, vcol, 0.0, 0.0, W, H, chan);
        for (px = 0; px < W * H; px++) prov[px * 3 + i] = chan[px];
    }
    memset(chan, 0, W * H);
    CloseupSurface_raster_attr(m.uv64, m.faces, m.nf, attr, 0.0, 0.0, W, H, chan);
    for (px = 0; px < W * H; px++) {
        pfill += chan[px] == CLOSEUP_NODE_FILL;
        pdata += chan[px] == CLOSEUP_NODE_DATA;
        padopt += chan[px] == 3;
    }
    snprintf(path, sizeof path, "%s/%s_provenance.png", c->out_dir, stem);
    VesPng_write_rgb(path, prov, (int)W, (int)H);

    snprintf(path, sizeof path, "%s/%s_render_report.json", c->out_dir, stem);
    if ((fp = fopen(path, "wb")) == NULL) goto out;
    fprintf(fp, "{\n  \"schema\": \"ink-closeup-render-v1\",\n  \"id\": \"%s\",\n  \"grid\": {\"origin_uv\": [0, 0], "
            "\"shape_hw\": [%zu, %zu], \"pixel_center_rule\": \"origin + (index + 0.5)\", \"px_per_vox\": 1},\n", stem, H, W);
    fprintf(fp, "  \"voxel_um\": %.4f,\n  \"layers\": {\"depth\": %d, \"step_vox\": %.3f, \"painted_px\": %zu, "
            "\"complete_px\": %zu, \"apron_px\": %zu, \"apron_rejected\": %zu},\n", c->voxel_um, SHEET_LAYERS_DEPTH,
            SHEET_LAYERS_STEP_VOX, ls.painted_px, ls.complete_px, ls.apron_px, ls.apron_rejected);
    fprintf(fp, "  \"provenance_px\": {\"data\": %zu, \"adopted\": %zu, \"fill\": %zu},\n  \"handedness_px\": "
            "{\"team\": %zu, \"mirror\": %zu},\n", pdata, padopt, pfill, team, mirror);
    fprintf(fp, "  \"ct_window\": [%.1f, %.1f],\n  \"seconds\": %.2f\n}\n", lo, hi, ves_clock_sec() - t0);
    fclose(fp); fp = NULL;
    fprintf(stderr, "[ink_closeup] SUMMARY render id=%s %zux%zu painted=%zu data=%zu adopted=%zu fill=%zu team=%zu "
            "mirror=%zu (%.1fs)\n", stem, W, H, ls.painted_px, pdata, padopt, pfill, team, mirror, ves_clock_sec() - t0);
    rc = 0;
out:
    if (fp != NULL) fclose(fp);
    RawtexHandedness_free(&hand);
    MeshBin_dispose(&m);
    free(nrm); free(attr); free(vchart); free(val); free(has); free(prov); free(chan); free(vcol);
    Arena_dispose(&arena);
    return rc;
}

/* ------------------------------------------------------------------ letters */

/* Bilinear sample of an 8-bit raster (0 outside). */
static double ic_bilinear(const uint8_t *img, size_t W, size_t H, double x, double y)
{
    double fx = x - 0.5, fy = y - 0.5, ax = 0.0, ay = 0.0, v = 0.0;
    long i = (long)floor(fx), j = (long)floor(fy);
    int q = 0;
    ax = fx - (double)i; ay = fy - (double)j;
    for (q = 0; q < 4; q++) {
        long ii = i + (q & 1), jj = j + (q >> 1);
        double w = ((q & 1) ? ax : 1.0 - ax) * ((q >> 1) ? ay : 1.0 - ay);
        if (ii < 0 || jj < 0 || (size_t)ii >= W || (size_t)jj >= H) continue;
        v += w * img[(size_t)jj * W + (size_t)ii];
    }
    return v;
}

static double ic_row_y(const InkRow *r, double x) { return r->a + r->b * x + r->c * x * x; }

static const uint8_t ic_row_rgb[6][3] = { {255, 220, 0}, {0, 230, 255}, {120, 255, 80}, {255, 120, 255}, {255, 160, 60},
                                          {160, 160, 255} };

/* ink_closeup <cfg> --letters <ink.tif> base|rect: letter boxes and rows of an ink map inferred on <stem>'s grid;
 * with base, also the rectification warp and the rectified surface <id>_rect_surface.vmesh. */
static int ic_stage_letters(const IcConfig *c, const char *which, const char *ink_path)
{
    char stem[512], path[2600];
    Arena_T arena = Arena_new();
    uint8_t *ink = NULL, *mask = NULL, *ct = NULL, *rgb = NULL, *gal = NULL, *attr = NULL;
    float *pos = NULL, *nrm = NULL;
    int32_t *vchart = NULL;
    size_t W = 0, H = 0, n = 0, nr = 0, k = 0, got = 0;
    int D = 0, Hm = 0, Wm = 0, rc = 1, make_warp = strcmp(which, "base") == 0;
    InkLettersOpts lo;
    InkLetter *L = NULL;
    InkRow *R = NULL;
    InkWarp wp;
    double hm = 0.0, t0 = ves_clock_sec(), px_mm = 1000.0 / c->voxel_um;
    MeshBinData m;
    AsmAxis *axis = NULL;
    FILE *fp = NULL;
    memset(&wp, 0, sizeof wp);
    memset(&m, 0, sizeof m);
    snprintf(stem, sizeof stem, make_warp ? "%s" : "%s_rect", c->id);
    if (RawtexTifxyz_load_gray(ink_path, &ink, &W, &H) != 0) { fprintf(stderr, "ERROR: cannot read %s\n", ink_path); goto out; }
    snprintf(path, sizeof path, "%s/%s_layers_mask.tif", c->out_dir, stem);
    if (TiffIO_load(arena, path, &mask, &D, &Hm, &Wm) != 0 || (size_t)Wm != W || (size_t)Hm != H) {
        fprintf(stderr, "ERROR: %s missing or not %zu x %zu\n", path, W, H);
        goto out;
    }
    snprintf(path, sizeof path, "%s/%s_ct.tif", c->out_dir, stem);
    if (TiffIO_load(arena, path, &ct, &D, &Hm, &Wm) != 0 || (size_t)Wm != W || (size_t)Hm != H) {
        fprintf(stderr, "ERROR: %s missing or not %zu x %zu\n", path, W, H);
        goto out;
    }
    InkLettersOpts_default(&lo);
    if (InkLetters_detect(ink, mask, W, H, &lo, &L, &n) != 0 || InkLetters_rows(L, n, &lo, &R, &nr, &hm) != 0) {
        fprintf(stderr, "ERROR: letter detection failed\n");
        goto out;
    }
    fprintf(stderr, "[ink_closeup] %s letters: %zu letters, %zu rows, median letter height %.0f px = %.2f mm (%.1fs)\n",
            stem, n, nr, hm, hm / px_mm, ves_clock_sec() - t0);
    for (k = 0; k < nr; k++)
        fprintf(stderr, "  row %zu: %zu letters, x %.0f..%.0f, height %.0f px, centre rms %.1f px about y = %.1f + %.3g x + %.3g x^2\n",
                k, R[k].n, R[k].x0, R[k].x1, R[k].height, R[k].rms, R[k].a, R[k].b, R[k].c);

    /* 3-D position of every pixel of the close-up (geometry only), for the letters' scroll coordinates */
    snprintf(path, sizeof path, "%s/%s_surface.vmesh", c->out_dir, stem);
    if (MeshBin_read_precise_malloc(path, &m) != 0 || m.uv64 == NULL) { fprintf(stderr, "ERROR: %s\n", path); goto out; }
    snprintf(path, sizeof path, "%s/%s_surface_nrm.f32", c->out_dir, stem);
    nrm = (float *)ic_read_file(path, m.nv * 3 * sizeof(float), &got);
    pos = (float *)malloc(3 * W * H * sizeof *pos);
    axis = AsmAxis_load(arena, c->axis);
    if (nrm == NULL || pos == NULL || axis == NULL) goto out;
    for (k = 0; k < 3 * W * H; k++) pos[k] = NAN;
    {
        RawtexLayersOpts opt;
        RawtexLayersStats ls;
        RawtexWindow win = { 0.0, (double)W, 0.0, (double)H };
        memset(&opt, 0, sizeof opt);
        opt.depth = SHEET_LAYERS_DEPTH; opt.step = SHEET_LAYERS_STEP_VOX; opt.tile_cols = SHEET_LAYERS_TILE_COLS;
        opt.chunk = SHEET_LAYERS_CHUNK; opt.stretch_ratio = SHEET_LAYERS_STRETCH_RATIO;
        opt.stretch_floor = SHEET_LAYERS_STRETCH_FLOOR; opt.window = &win;
        opt.pos_sink = pos; opt.sink_w = W; opt.sink_h = H;
        if (RawtexLayers_write(NULL, NULL, NULL, m.verts, (RawtexUv){ NULL, m.uv64 }, m.nv, m.faces, m.nf, nrm, NULL,
                               &opt, NULL, &ls) != 0) {
            fprintf(stderr, "ERROR: geometry raster failed\n");
            goto out;
        }
    }

    /* letters.json */
    snprintf(path, sizeof path, "%s/%s_letters.json", c->out_dir, stem);
    if ((fp = fopen(path, "wb")) == NULL) goto out;
    fprintf(fp, "{\n  \"schema\": \"ink-closeup-letters-v1\",\n  \"id\": \"%s\",\n  \"ink\": ", stem);
    {
        char js[4800];
        ic_json_str(js, sizeof js, ink_path);
        fprintf(fp, "%s,\n", js);
    }
    fprintf(fp, "  \"voxel_um\": %.4f,\n  \"px_per_mm\": %.3f,\n  \"median_letter_height_px\": %.1f,\n"
            "  \"median_letter_height_mm\": %.3f,\n  \"detector\": {\"smooth_px\": %.1f, \"t_hi\": %.2f, \"t_lo\": %.2f, "
            "\"close_px\": %d, \"min_area\": %.0f, \"max_area\": %.0f},\n  \"rows\": [", c->voxel_um, px_mm, hm, hm / px_mm,
            lo.smooth_px, lo.t_hi, lo.t_lo, lo.close_px, lo.min_area, lo.max_area);
    for (k = 0; k < nr; k++)
        fprintf(fp, "%s\n    {\"row\": %zu, \"letters\": %zu, \"x0\": %.1f, \"x1\": %.1f, \"height_px\": %.1f, "
                "\"centre_rms_px\": %.2f, \"curve\": [%.6g, %.6g, %.6g]}", k ? "," : "", k, R[k].n, R[k].x0, R[k].x1,
                R[k].height, R[k].rms, R[k].a, R[k].b, R[k].c);
    fprintf(fp, "\n  ],\n  \"letters\": [");
    for (k = 0; k < n; k++) {
        size_t xi = (size_t)(L[k].cx < 0 ? 0 : L[k].cx), yi = (size_t)(L[k].cy < 0 ? 0 : L[k].cy);
        double q[3] = { NAN, NAN, NAN }, sa = NAN, rr = NAN, fr[3] = { NAN, NAN, NAN };
        if (xi < W && yi < H) {
            size_t pl = W * H, px = yi * W + xi;
            q[0] = pos[px]; q[1] = pos[pl + px]; q[2] = pos[2 * pl + px];
            if (isfinite(q[0])) {
                AsmAxis_project(axis, q, &sa, &rr, NULL, NULL);
                AsmAxis_frame(axis, q, fr);
            }
        }
        {   /* a letter whose centre has no geometry has no 3-D position: JSON null, never a bare nan */
            char b0[32], b1[32], b2[32], b3[32], b4[32], b5[32];
            double th = isfinite(fr[1]) ? atan2(fr[2], fr[1]) * 180.0 / 3.14159265358979 : NAN;
            fprintf(fp, "%s\n    {\"id\": %zu, \"row\": %d, \"cx\": %.1f, \"cy\": %.1f, \"angle_deg\": %.1f, \"half\": [%.1f, %.1f], "
                    "\"box\": [%.0f, %.0f, %.0f, %.0f], \"area\": %.0f, \"mean_p\": %.3f, \"max_p\": %.3f, "
                    "\"zyx\": [%s, %s, %s], \"axial\": %s, \"radius\": %s, \"theta_deg\": %s}", k ? "," : "", k,
                    L[k].row, L[k].cx, L[k].cy, L[k].angle * 180.0 / 3.14159265358979, L[k].ha, L[k].hb, L[k].bx0, L[k].by0,
                    L[k].bx1, L[k].by1, L[k].area, L[k].mean_p, L[k].max_p, ic_json_num(b0, sizeof b0, q[0]),
                    ic_json_num(b1, sizeof b1, q[1]), ic_json_num(b2, sizeof b2, q[2]), ic_json_num(b3, sizeof b3, sa),
                    ic_json_num(b4, sizeof b4, rr), ic_json_num(b5, sizeof b5, th));
        }
    }
    fprintf(fp, "\n  ]\n}\n");
    fclose(fp); fp = NULL;

    /* the overlay: CT grey, ink red, letter boxes by row, row curves, 1 cm scale bar */
    rgb = (uint8_t *)malloc(W * H * 3);
    if (rgb == NULL) goto out;
    for (k = 0; k < W * H; k++) {
        double pv = (double)ink[k] / 255.0, a = (pv - 0.35) / 0.4, g = (double)ct[k];
        if (a < 0.0 || !mask[k]) a = 0.0;
        if (a > 1.0) a = 1.0;
        rgb[k * 3] = (uint8_t)(g * (1.0 - a) + 255.0 * a);
        rgb[k * 3 + 1] = (uint8_t)(g * (1.0 - 0.8 * a));
        rgb[k * 3 + 2] = (uint8_t)(g * (1.0 - 0.8 * a));
    }
    for (k = 0; k < nr; k++) {
        const uint8_t *col = ic_row_rgb[k % 6];
        double x = 0.0, x0 = R[k].x0 - hm, x1 = R[k].x1 + hm;
        for (x = x0; x < x1; x += 8.0)
            InkDraw_line(rgb, W, H, x, ic_row_y(&R[k], x), x + 8.0, ic_row_y(&R[k], x + 8.0), 3, col[0], col[1], col[2]);
    }
    for (k = 0; k < n; k++) {
        const uint8_t *col = L[k].row >= 0 ? ic_row_rgb[L[k].row % 6] : (const uint8_t *)"\xb0\xb0\xb0";
        char lab[32];
        InkDraw_obb(rgb, W, H, &L[k], 3, col[0], col[1], col[2]);
        snprintf(lab, sizeof lab, "%zu", k);
        InkDraw_text(rgb, W, H, (int)(L[k].bx0), (int)(L[k].by0) - 30, 4, lab, col[0], col[1], col[2]);
    }
    {
        int x0 = 60, y0 = (int)H - 60, len = (int)floor(10.0 * px_mm + 0.5);
        char lab[160];
        InkDraw_line(rgb, W, H, x0, y0, x0 + len, y0, 14, 255, 255, 255);
        InkDraw_text(rgb, W, H, x0, y0 - 50, 5, "1 cm", 255, 255, 255);
        snprintf(lab, sizeof lab, "%s  reader-v2  %.2f um/vox  %zu letters  %zu rows  letter height %.2f mm", stem,
                 c->voxel_um, n, nr, hm / px_mm);
        InkDraw_text(rgb, W, H, 30, 20, 4, lab, 255, 255, 255);
    }
    snprintf(path, sizeof path, "%s/%s_letters_overlay.png", c->out_dir, stem);
    VesPng_write_rgb(path, rgb, (int)W, (int)H);

    /* the gallery: each letter upright along its row (row tangent), CT above ink, in reading order */
    {
        const int T = 176, lab_h = 30, per_line = 12;
        size_t *ord = (size_t *)malloc((n + 1) * sizeof *ord), no = 0, r = 0, GW = 0, GH = 0, lines = 0;
        double S = 1.9 * hm;
        if (ord == NULL) goto out;
        for (r = 0; r < nr; r++) {           /* rows top to bottom, letters left to right */
            size_t a = 0, b = 0;
            size_t first = no;
            for (a = 0; a < n; a++) if (L[a].row == (int)r) ord[no++] = a;
            for (a = first; a < no; a++)
                for (b = a + 1; b < no; b++)
                    if (L[ord[b]].cx < L[ord[a]].cx) { size_t t = ord[a]; ord[a] = ord[b]; ord[b] = t; }
        }
        for (k = 0; k < n; k++) if (L[k].row < 0) ord[no++] = k;
        lines = (no + (size_t)per_line - 1) / (size_t)per_line;
        GW = (size_t)per_line * (size_t)(T + 8) + 8;
        GH = lines * (size_t)(2 * T + lab_h + 16) + 8;
        gal = (uint8_t *)calloc(GW * GH * 3, 1);
        if (gal == NULL || no == 0) { free(ord); goto gallery_done; }
        for (k = 0; k < no; k++) {
            const InkLetter *e = &L[ord[k]];
            double th = e->row >= 0 ? atan(R[e->row].b + 2.0 * R[e->row].c * e->cx) : 0.0, ca = cos(th), sa = sin(th);
            size_t gx = 8 + (k % (size_t)per_line) * (size_t)(T + 8), gy = 8 + (k / (size_t)per_line) * (size_t)(2 * T + lab_h + 16);
            int tx = 0, ty = 0, pass = 0;
            char lab[32];
            for (pass = 0; pass < 2; pass++)
                for (ty = 0; ty < T; ty++)
                    for (tx = 0; tx < T; tx++) {
                        double lx = ((tx + 0.5) / T - 0.5) * S, ly = ((ty + 0.5) / T - 0.5) * S;
                        double x = e->cx + lx * ca - ly * sa, y = e->cy + lx * sa + ly * ca;
                        double v = ic_bilinear(pass == 0 ? ct : ink, W, H, x, y);
                        size_t o = ((gy + (size_t)lab_h + (size_t)(pass * T) + (size_t)ty) * GW + gx + (size_t)tx) * 3;
                        if (pass == 0) { gal[o] = gal[o + 1] = gal[o + 2] = (uint8_t)v; }
                        else {
                            double a = (v / 255.0 - 0.3) / 0.45;
                            if (a < 0.0) a = 0.0;
                            if (a > 1.0) a = 1.0;
                            gal[o] = (uint8_t)(255.0 * a); gal[o + 1] = (uint8_t)(255.0 * a * a); gal[o + 2] = (uint8_t)(90.0 * a * a * a);
                        }
                    }
            snprintf(lab, sizeof lab, "%zu r%d", ord[k], e->row);
            InkDraw_text(gal, GW, GH, (int)gx + 2, (int)gy + 2, 3, lab, 255, 255, 255);
        }
        snprintf(path, sizeof path, "%s/%s_letters_gallery.png", c->out_dir, stem);
        VesPng_write_rgb(path, gal, (int)GW, (int)GH);
        free(ord);
    }
gallery_done:

    /* base: the rectification warp, applied to the surface UV (then re-rendered from 3-D) */
    if (make_warp) {
        double umin = 1e300, vmin = 1e300;
        size_t v = 0;
        if (nr == 0) { fprintf(stderr, "  no row: no rectification\n"); rc = 0; goto out; }
        if (InkLetters_warp(L, n, R, nr, W, H, &lo, &wp) != 0) { fprintf(stderr, "ERROR: warp failed\n"); goto out; }
        fprintf(stderr, "  rectification: letter centres about their row levels rms %.1f -> %.2f px; local scale %.3f..%.3f, "
                "rotation <= %.1f deg\n", wp.rms_before, wp.rms_after, wp.scale_min, wp.scale_max, wp.rot_max_deg);
        for (v = 0; v < m.nv; v++) {
            double dx = 0.0, dy = 0.0;
            InkWarp_eval(&wp, m.uv64[v * 2], m.uv64[v * 2 + 1], &dx, &dy);
            m.uv64[v * 2] += dx; m.uv64[v * 2 + 1] += dy;
            if (m.uv64[v * 2] < umin) umin = m.uv64[v * 2];
            if (m.uv64[v * 2 + 1] < vmin) vmin = m.uv64[v * 2 + 1];
        }
        for (v = 0; v < m.nv; v++) {
            m.uv64[v * 2] += (double)CLOSEUP_RASTER_MARGIN_PX - floor(umin);
            m.uv64[v * 2 + 1] += (double)CLOSEUP_RASTER_MARGIN_PX - floor(vmin);
        }
        snprintf(path, sizeof path, "%s/%s_rect_surface.vmesh", c->out_dir, c->id);
        if (MeshBin_write_uv64(path, m.verts, m.nv, m.faces, m.nf, m.uv64) != 0) goto out;
        snprintf(path, sizeof path, "%s/%s_rect_surface_nrm.f32", c->out_dir, c->id);
        if (ic_write_file(path, nrm, m.nv * 3 * sizeof *nrm) != 0) goto out;
        snprintf(path, sizeof path, "%s/%s_surface_attr.u8", c->out_dir, c->id);
        attr = (uint8_t *)ic_read_file(path, m.nv, &got);
        snprintf(path, sizeof path, "%s/%s_surface_chart.i32", c->out_dir, c->id);
        vchart = (int32_t *)ic_read_file(path, m.nv * sizeof(int32_t), &got);
        if (attr == NULL || vchart == NULL) goto out;
        snprintf(path, sizeof path, "%s/%s_rect_surface_attr.u8", c->out_dir, c->id);
        if (ic_write_file(path, attr, m.nv) != 0) goto out;
        snprintf(path, sizeof path, "%s/%s_rect_surface_chart.i32", c->out_dir, c->id);
        if (ic_write_file(path, vchart, m.nv * sizeof *vchart) != 0) goto out;
        snprintf(path, sizeof path, "%s/%s_rect_warp.json", c->out_dir, c->id);
        if ((fp = fopen(path, "wb")) == NULL) goto out;
        fprintf(fp, "{\n  \"schema\": \"ink-closeup-warp-v1\",\n  \"grid_px\": %.1f,\n  \"rows\": %zu,\n  \"letters\": %zu,\n"
                "  \"rms_before_px\": %.3f,\n  \"rms_after_px\": %.3f,\n  \"scale_min\": %.4f,\n  \"scale_max\": %.4f,\n"
                "  \"rotation_max_deg\": %.2f,\n  \"weights\": {\"row\": %.3g, \"conformal\": %.3g, \"smooth\": %.3g, "
                "\"anchor\": %.3g}\n}\n", wp.step, nr, n, wp.rms_before, wp.rms_after, wp.scale_min, wp.scale_max,
                wp.rot_max_deg, lo.w_row, lo.w_conf, lo.w_smooth, lo.w_anchor);
        fclose(fp); fp = NULL;
    }
    fprintf(stderr, "[ink_closeup] SUMMARY letters id=%s letters=%zu rows=%zu h=%.0fpx (%.1fs)\n", stem, n, nr, hm,
            ves_clock_sec() - t0);
    rc = 0;
out:
    if (fp != NULL) fclose(fp);
    InkWarp_free(&wp);
    MeshBin_dispose(&m);
    free(ink); free(rgb); free(gal); free(pos); free(nrm); free(attr); free(vchart); free(L); free(R);
    Arena_dispose(&arena);
    return rc;
}

/* ------------------------------------------------------------------ main */

static int ic_adopt_opts_equal(const CloseupAdoptOpts *x, const CloseupAdoptOpts *y)   /* field by field: no padding */
{
    return x->on_vox == y->on_vox && x->foot_fill_px == y->foot_fill_px && x->min_extent_px == y->min_extent_px &&
           x->rms_px == y->rms_px && x->conflict_vox == y->conflict_vox && x->agree_frac == y->agree_frac &&
           x->normal_deg == y->normal_deg && x->new_frac == y->new_frac && x->min_corr == y->min_corr &&
           x->metric_ratio == y->metric_ratio && x->metric_frac == y->metric_frac &&
           x->continuity_px == y->continuity_px && x->continuity_slope == y->continuity_slope &&
           x->layer_vox == y->layer_vox && x->layer_frac == y->layer_frac && x->layer_margin_px == y->layer_margin_px;
}

/* The "adopt" overrides: absent = the compiled defaults (byte-identical runs), present = applied, out of range =
 * refused.  Configs are written under VES_SELFTEST_DIR (the repo's rule for selftest scratch, 2026-10-12). */
static int ic_selftest_adopt_config(void)
{
    static const char *base = "{\"id\": \"t\", \"out_dir\": \"o\", \"layout\": \"l\", \"raw\": \"r\", \"axis\": \"a\", "
                              "\"voxel_um\": 8.64, \"seed_box\": [0, 0, 10, 10]";
    const char *path = VES_SELFTEST_DIR "/ink_closeup_selftest_adopt.json";
    const char *tails[3] = {"}", ", \"adopt\": {\"rms_px\": 4.0, \"normal_deg\": 30, \"rounds\": 8}}",
                            ", \"adopt\": {\"rms_px\": -1}}"};
    CloseupAdoptOpts def;
    IcConfig c;
    int k = 0, rc[3] = {0, 0, 0}, bad = 0;
    FILE *fp = NULL;
    CloseupAdoptOpts_default(&def);
    ves_ensure_parent_dir(path);
    for (k = 0; k < 3; k++) {
        fp = fopen(path, "wb");
        if (fp == NULL) return 1;
        fprintf(fp, "%s%s\n", base, tails[k]);
        fclose(fp);
        rc[k] = ic_config(path, &c);
        if (k == 0 && (rc[k] != 0 || !ic_adopt_opts_equal(&c.adopt, &def) || c.adopt_rounds != CLOSEUP_ADOPT_ROUNDS))
            bad = 1;
        if (k == 1 && (rc[k] != 0 || fabs(c.adopt.rms_px - 4.0) > 1e-12 || fabs(c.adopt.normal_deg - 30.0) > 1e-12 ||
                       c.adopt_rounds != 8 || fabs(c.adopt.conflict_vox - def.conflict_vox) > 1e-12))
            bad = 1;
        if (k == 2 && rc[k] == 0) bad = 1;
        free(c.pool);
    }
    remove(path);
    return bad;
}

/* The near region: a node's own cell and every cell within the reach (box metric) of it, nothing beyond; a chart
 * whose box reaches into it is a candidate, one beside it, or outside the search box, is not. */
static int ic_near_selftest(void)
{
    float pos[6] = { 100.0f, 100.0f, 100.0f,  NAN, NAN, NAN };
    uint8_t state[2] = { CLOSEUP_NODE_DATA, CLOSEUP_NODE_NONE };
    const double reach = 3.0 * CLOSEUP_POOL_CELL_VOX;
    double lo[3] = { 0.0, 0.0, 0.0 }, hi[3] = { 400.0, 400.0, 400.0 };
    double node[3] = { 100.0, 100.0, 100.0 }, p_in[3] = { 100.0 + reach, 100.0 - reach, 100.0 };
    double p_out[3] = { 100.0, 100.0, 100.0 + 2.0 * reach + 2.0 * CLOSEUP_POOL_CELL_VOX };
    double blo[3] = { 300.0, 300.0, 300.0 }, bhi[3] = { 310.0, 310.0, 310.0 };
    double llo[3] = { 0.0, 0.0, 90.0 }, lhi[3] = { 390.0, 390.0, 100.0 };
    double olo[3] = { 500.0, 100.0, 100.0 }, ohi[3] = { 600.0, 100.0, 100.0 };
    IcNear nr;
    int bad = 0;
    if (ic_near_build(pos, state, 2, lo, hi, reach, &nr) != 0) return 1;
    if (!ic_pool_candidate(node, node, lo, hi, &nr) || !ic_pool_candidate(p_in, p_in, lo, hi, &nr) ||
        ic_pool_candidate(p_out, p_out, lo, hi, &nr)) bad = 1;
    if (ic_pool_candidate(blo, bhi, lo, hi, &nr) || !ic_pool_candidate(llo, lhi, lo, hi, &nr)) bad = 1;
    if (ic_pool_candidate(olo, ohi, lo, hi, NULL) || !ic_pool_candidate(blo, bhi, lo, hi, NULL)) bad = 1;
    ic_near_free(&nr);
    {   /* the grid's reach: corner to far corner of 5 x 3 nodes 2 px apart = 2 diagonal + 2 straight steps */
        uint8_t st[15];
        CloseupGrid gg;
        memset(st, 0, sizeof st);
        memset(&gg, 0, sizeof gg);
        gg.ni = 5; gg.nj = 3; gg.h = 2.0; gg.state = st;
        if (ic_grid_reach(&gg) != 0.0) bad = 1;                     /* no supported node */
        st[0] = CLOSEUP_NODE_DATA;
        if (fabs(ic_grid_reach(&gg) - (4.0 + 4.0 * sqrt(2.0))) > 1e-9) bad = 1;
        st[14] = CLOSEUP_NODE_FILL;
        if (fabs(ic_grid_reach(&gg) - (2.0 + 2.0 * sqrt(2.0))) > 1e-9) bad = 1;
    }
    {   /* no supported node: no region (the search box alone decides) */
        uint8_t none[2] = { CLOSEUP_NODE_NONE, CLOSEUP_NODE_NONE };
        double ilo[3] = { 1e300, 1e300, 1e300 }, ihi[3] = { -1e300, -1e300, -1e300 };
        if (ic_near_build(pos, none, 2, ilo, ihi, reach, &nr) != 1) bad = 1;
    }
    return bad;
}

/* A sheet of the cylinder of radius R about the z axis: (u, v) = (arc length, z) -> world (z, y, x). */
static void ic_test_cyl(double R, double u, double v, float p[3])
{
    p[0] = (float)v; p[1] = (float)(R * sin(u / R)); p[2] = (float)(R * cos(u / R));
}

/* A regular patch [u0, u0 + side] x [v0, v0 + side] of that cylinder, 2 px apart, appended to verts/uv/faces (vertex
 * k's UV = its (u, v) rotated by theta and shifted by (tu, tv): any rigid frame of the chart's own flattening). */
static size_t ic_test_patch(double R, double u0, double v0, double side, double theta, double tu, double tv,
                            float *verts, double *uv, int32_t *faces, size_t *nv, size_t *nf)
{
    size_t n = (size_t)(side / 2.0) + 1, i = 0, j = 0, base = *nv;
    double cs = cos(theta), sn = sin(theta);
    for (j = 0; j < n; j++)
        for (i = 0; i < n; i++) {
            double u = u0 + 2.0 * (double)i, v = v0 + 2.0 * (double)j;
            ic_test_cyl(R, u, v, &verts[*nv * 3]);
            uv[*nv * 2] = cs * u - sn * v + tu;
            uv[*nv * 2 + 1] = sn * u + cs * v + tv;
            (*nv)++;
        }
    for (j = 0; j + 1 < n; j++)
        for (i = 0; i + 1 < n; i++) {
            int32_t a = (int32_t)(base + j * n + i), b = a + 1, cc = a + (int32_t)n, d = cc + 1;
            faces[*nf * 3] = a; faces[*nf * 3 + 1] = b; faces[*nf * 3 + 2] = d; (*nf)++;
            faces[*nf * 3] = a; faces[*nf * 3 + 1] = d; faces[*nf * 3 + 2] = cc; (*nf)++;
        }
    return n * n;
}

/* The growth loop end to end on a cylinder sheet (R 300): a 160 px seed close-up, a pool of 25 tiles of 96 px (each in
 * its own rotated frame) over the 480 px target box around it; the pool is read whole once (and indexed), then paged.
 * The grown surface must cover most of the target box, keep the inward orientation, and its data must lie on the
 * cylinder.  Scratch under selftest/ (removed on a pass). */
static int ic_reuv_guard_selftest(void)
{
    Arena_T arena = Arena_new();
    double safe[50], target[50], expected[50];
    int32_t faces[96];
    size_t nf = 0;
    int fails = 0;
    IcReuvGuardStats st;
    for (int y = 0; y < 5; y++) for (int x = 0; x < 5; x++) {
        size_t v = (size_t)(5*y+x);
        safe[2*v] = x; safe[2*v+1] = y;
        expected[2*v] = x + .05*y; expected[2*v+1] = y;
    }
    for (int y = 0; y < 4; y++) for (int x = 0; x < 4; x++) {
        int32_t a = 5*y+x, b = a+1, c = a+5, d = c+1;
        faces[3*nf] = a; faces[3*nf+1] = b; faces[3*nf+2] = d; nf++;
        faces[3*nf] = a; faces[3*nf+1] = d; faces[3*nf+2] = c; nf++;
    }
    memcpy(target,expected,sizeof target);
    if (ic_reuv_guard(arena,safe,target,25,faces,nf,&st) != 0 || st.limited_vertices != 0 ||
        memcmp(target,expected,sizeof target) != 0) fails++;
    target[24] = 4.5; target[25] = 1.5;       /* local inversion in an otherwise valid guide */
    if (ic_reuv_guard(arena,safe,target,25,faces,nf,&st) != 0 || st.unsafe_paths == 0 ||
        st.limited_vertices == 0 || st.limited_vertices >= 25 || !UvGuard_mesh(target,faces,nf) ||
        target[0] != expected[0] || target[1] != expected[1]) fails++;
    for (size_t v = 0; v < 25; v++) { target[2*v] = -safe[2*v]; target[2*v+1] = -safe[2*v+1]; }
    if (ic_reuv_guard(arena,safe,target,25,faces,nf,&st) != 0 || st.unsafe_paths == 0) fails++;
    for (size_t f = 0; f < nf; f++) if (!UvGuard_interval(safe,target,faces+3*f,1.0)) fails++;
    Arena_dispose(&arena);
    fprintf(stderr,"[ink_closeup] reuv orientation guard selftest: %s (%d failures)\n",fails ? "FAIL" : "ok",fails);
    return fails;
}

static int ic_reuv_selftest(void)
{
    char dir[512], path[700], stem[700], cfg[700], fdir[700], seg[700];
    float *verts = NULL, *MU = NULL, *MV = NULL;
    double *uv = NULL;
    int32_t *faces = NULL, *chart = NULL;
    uint8_t *attr = NULL;
    size_t nv = 0, nf = 0, i = 0, j = 0, bad = 0;
    const int gw = 21, gh = 16;
    const double th = 10.0 * 3.14159265358979323846 / 180.0, cs = cos(th), sn = sin(th);
    int fails = ic_reuv_guard_selftest() + FiberMap_selftest();
    IcConfig c;
    MeshBinData m;
    FILE *fp = NULL;
    memset(&m, 0, sizeof m);
    memset(&c, 0, sizeof c);
    if (ves_selftest_dir(dir, sizeof dir, "ink_closeup_reuv") != 0) return 1;
    verts = (float *)malloc(40 * 30 * 3 * sizeof *verts); uv = (double *)malloc(40 * 30 * 2 * sizeof *uv);
    faces = (int32_t *)malloc(39 * 29 * 6 * sizeof *faces); chart = (int32_t *)malloc(40 * 30 * sizeof *chart);
    attr = (uint8_t *)malloc(40 * 30); MU = (float *)malloc((size_t)gw * gh * sizeof *MU); MV = (float *)malloc((size_t)gw * gh * sizeof *MV);
    if (!verts || !uv || !faces || !chart || !attr || !MU || !MV) { fails++; goto done; }
    for (j = 0; j < 30; j++)                      /* UV = (100 + 4 i, 50 + 4 j): canvas = UV - (100, 50) */
        for (i = 0; i < 40; i++) {
            verts[nv * 3] = 0.0f; verts[nv * 3 + 1] = (float)(4 * j); verts[nv * 3 + 2] = (float)(4 * i);
            uv[nv * 2] = 100.0 + 4.0 * (double)i; uv[nv * 2 + 1] = 50.0 + 4.0 * (double)j;
            chart[nv] = 7; attr[nv] = 1; nv++;
        }
    for (j = 0; j + 1 < 30; j++)
        for (i = 0; i + 1 < 40; i++) {
            const int32_t a = (int32_t)(j * 40 + i), b = a + 1, cc = a + 40, d = cc + 1;
            faces[nf * 3] = a; faces[nf * 3 + 1] = b; faces[nf * 3 + 2] = d; nf++;
            faces[nf * 3] = a; faces[nf * 3 + 1] = d; faces[nf * 3 + 2] = cc; nf++;
        }
    snprintf(stem, sizeof stem, "%s/src", dir);
    snprintf(path, sizeof path, "%s_surface.vmesh", stem); fails += MeshBin_write_uv64(path, verts, nv, faces, nf, uv) != 0;
    snprintf(path, sizeof path, "%s_surface_attr.u8", stem); fails += ic_write_file(path, attr, nv) != 0;
    snprintf(path, sizeof path, "%s_surface_chart.i32", stem); fails += ic_write_file(path, chart, nv * sizeof *chart) != 0;
    snprintf(seg, sizeof seg, "%s/seg", dir); (void)ves_mkdir(seg);
    snprintf(path, sizeof path, "%s/meta.json", seg);
    fp = fopen(path, "wb"); if (fp == NULL) { fails++; goto done; }
    fprintf(fp, "{\"type\": \"seg\", \"scale\": [1.0, 1.0], \"uv_origin\": [100.0, 50.0]}\n"); fclose(fp); fp = NULL;
    snprintf(fdir, sizeof fdir, "%s/fiber", dir); (void)ves_mkdir(fdir);
    for (j = 0; j < (size_t)gh; j++)               /* node (i, j) at canvas (row 8 j, col 8 i): u' = R (col, row) + (500, 300) */
        for (i = 0; i < (size_t)gw; i++) {
            const double col = 8.0 * (double)i, row = 8.0 * (double)j;
            MU[j * (size_t)gw + i] = (float)(cs * col - sn * row + 500.0);
            MV[j * (size_t)gw + i] = (float)(sn * col + cs * row + 300.0);
        }
    MU[0] = MV[0] = NAN;                          /* one unsolved node: its cell falls back to the nearest node */
    snprintf(path, sizeof path, "%s/fiber_map_u.tif", fdir); fails += TiffIO_save_float2d(path, MU, gw, gh) != 0;
    snprintf(path, sizeof path, "%s/fiber_map_v.tif", fdir); fails += TiffIO_save_float2d(path, MV, gw, gh) != 0;
    snprintf(path, sizeof path, "%s/fiber_flatten.json", fdir);
    fp = fopen(path, "wb"); if (fp == NULL) { fails++; goto done; }
    fprintf(fp, "{\"seg\": \"%s\", \"grid\": {\"rows\": %d, \"cols\": %d, \"cell\": 8, \"row0\": 0, \"col0\": 0}, "
            "\"map\": {\"affine\": [[%.9f, %.9f, 500.0], [%.9f, %.9f, 300.0]]}, \"new_page\": {\"u_min\": 450.0, \"v_min\": 250.0}}\n",
            seg, gh, gw, cs, -sn, sn, cs);
    fclose(fp); fp = NULL;
    snprintf(cfg, sizeof cfg, "%s/reuv.json", dir);
    fp = fopen(cfg, "wb"); if (fp == NULL) { fails++; goto done; }
    fprintf(fp, "{\"id\": \"aligned\", \"out_dir\": \"%s/out\", \"layout\": \"%s/none\", \"raw\": \"%s/none\", \"axis\": \"%s/none\", "
            "\"voxel_um\": 9.362, \"seed_box\": [0, 0, 1, 1], \"surface_source\": \"%s\"}\n", dir, dir, dir, dir, stem);
    fclose(fp); fp = NULL;
    if (ic_config(cfg, &c) != 0 || ic_stage_reuv(&c, fdir) != 0) { fprintf(stderr, "  reuv selftest FAIL: the stage\n"); fails++; goto done; }
    snprintf(path, sizeof path, "%s/out/aligned_surface.vmesh", dir);
    if (MeshBin_read_precise_malloc(path, &m) != 0 || m.uv64 == NULL || m.nv != nv || m.nf != nf) { fails++; goto done; }
    for (i = 0; i < nv; i++) {                    /* expected: R (canvas) + (500, 300) - (450, 250) */
        const double col = uv[i * 2] - 100.0, row = uv[i * 2 + 1] - 50.0;
        const double eu = cs * col - sn * row + 50.0, ev = sn * col + cs * row + 50.0;
        if (fabs(m.uv64[i * 2] - eu) > 1e-3 || fabs(m.uv64[i * 2 + 1] - ev) > 1e-3) bad++;
    }
    snprintf(path, sizeof path, "%s/out/aligned_surface_chart.i32", dir);
    {
        size_t got = 0;
        int32_t *ch2 = (int32_t *)ic_read_file(path, nv * sizeof *ch2, &got);
        if (ch2 == NULL || ch2[5] != 7) fails++;
        free(ch2);
    }
    if (bad != 0) { fprintf(stderr, "  reuv selftest FAIL: %zu of %zu vertices off the rotated frame\n", bad, nv); fails++; }
    if (c.pool) free(c.pool);
done:
    MeshBin_dispose(&m);
    free(verts); free(uv); free(faces); free(chart); free(attr); free(MU); free(MV);
    fprintf(stderr, "[ink_closeup] reuv selftest: %s\n", fails ? "FAIL" : "ok");
    if (fails == 0) (void)ves_selftest_clean(dir);
    return fails ? 1 : 0;
}

static int ic_grow_selftest(double tile_px, int stop)
{
    char dir[512], path[700], stem[700], cfg[700], pool[700], tiles[64];
    const double R = 300.0;
    float *verts = NULL;
    double *uv = NULL;
    int32_t *faces = NULL, *chart = NULL;
    uint8_t *attr = NULL;
    size_t nv = 0, nf = 0, cap_v = 0, cap_f = 0, v = 0, t = 0, data = 0;
    int fails = 0;
    IcConfig c;
    FILE *fp = NULL;
    MeshBinData m;
    memset(&m, 0, sizeof m);
    if (ves_selftest_dir(dir, sizeof dir, stop ? "ink_closeup_grow_stop" : tile_px > 0.0 ? "ink_closeup_grow_tiled"
                                                                                     : "ink_closeup_grow") != 0) return 1;
    tiles[0] = 0;
    if (tile_px > 0.0) snprintf(tiles, sizeof tiles, ", \"tile_px\": %.1f", tile_px);
    cap_v = 26u * 81u * 81u; cap_f = 26u * 2u * 80u * 80u;
    verts = (float *)malloc(cap_v * 3 * sizeof *verts);
    uv = (double *)malloc(cap_v * 2 * sizeof *uv);
    faces = (int32_t *)malloc(cap_f * 3 * sizeof *faces);
    chart = (int32_t *)malloc(cap_v * sizeof *chart);
    attr = (uint8_t *)malloc(cap_v);
    if (verts == NULL || uv == NULL || faces == NULL || chart == NULL || attr == NULL) { fails++; goto done; }
    /* the seed close-up: [0, 160]^2 in its frame = (u, v), every vertex observed, chart 0 */
    ic_test_patch(R, 0.0, 0.0, 160.0, 0.0, 0.0, 0.0, verts, uv, faces, &nv, &nf);
    for (v = 0; v < nv; v++) { attr[v] = 1; chart[v] = 0; }
    snprintf(stem, sizeof stem, "%s/seed", dir);
    snprintf(path, sizeof path, "%s_surface.vmesh", stem);
    fails += MeshBin_write_uv64(path, verts, nv, faces, nf, uv) != 0;
    snprintf(path, sizeof path, "%s_surface_attr.u8", stem);
    fails += ic_write_file(path, attr, nv) != 0;
    snprintf(path, sizeof path, "%s_surface_chart.i32", stem);
    fails += ic_write_file(path, chart, nv * sizeof *chart) != 0;
    /* the pool: 5 x 5 tiles of 96 px from -160, each its own chart and frame */
    nv = nf = 0;
    for (t = 0; t < 25; t++) {
        size_t n0 = nv, k = 0;
        ic_test_patch(R, -160.0 + 96.0 * (double)(t % 5), -160.0 + 96.0 * (double)(t / 5), 96.0, 0.37 * (double)t,
                      31.0 * (double)t, -17.0 * (double)t, verts, uv, faces, &nv, &nf);
        for (k = n0; k < nv; k++) chart[k] = (int32_t)t;
    }
    snprintf(pool, sizeof pool, "%s/pool", dir);
    (void)ves_mkdir(pool);
    snprintf(path, sizeof path, "%s/stage6_repaired.vmesh", pool);
    fails += MeshBin_write_uv64(path, verts, nv, faces, nf, uv) != 0;
    snprintf(path, sizeof path, "%s/stage6_repaired_chart.i32", pool);
    fails += ic_write_file(path, chart, nv * sizeof *chart) != 0;
    snprintf(path, sizeof path, "%s/axis.csv", dir);
    fp = fopen(path, "wb");
    if (fp == NULL) { fails++; goto done; }
    fprintf(fp, "# test axis\nZ_l0,y_l0,x_l0\n-1024,0,0\n-512,0,0\n0,0,0\n512,0,0\n1024,0,0\n");
    fclose(fp);
    snprintf(cfg, sizeof cfg, "%s/grow.json", dir);
    fp = fopen(cfg, "wb");
    if (fp == NULL) { fails++; goto done; }
    fprintf(fp, "{\"id\": \"grown\", \"out_dir\": \"%s/out\", \"layout\": \"%s/none\", \"raw\": \"%s/none\", "
            "\"axis\": \"%s/axis.csv\", \"voxel_um\": 9.362, \"seed_box\": [0, 0, 1, 1], \"surface_source\": \"%s\", "
            "\"growth_uv\": [-160, -160, 320, 320], \"pool\": [\"%s\"], "
            "\"grow\": {\"margin_vox\": 128, \"rounds\": 6, \"max_reads\": 20, \"reflatten\": false%s}}\n",
            dir, dir, dir, dir, stem, pool, tiles);
    fclose(fp);
    if (ic_config(cfg, &c) != 0 || !c.has_grow || fabs(c.grow_tile_px - tile_px) > 1e-9) { fails++; goto done; }
    if (stop) {   /* the graceful stop: <out_dir>/STOP before the first read -> the seed's surface, nothing grown */
        snprintf(path, sizeof path, "%s/out", dir);
        (void)ves_mkdir(path);
        snprintf(path, sizeof path, "%s/out/STOP", dir);
        fp = fopen(path, "wb");
        if (fp == NULL) { fails++; goto done; }
        fclose(fp);
    }
    if (ic_stage_grow(&c) != 0) { fprintf(stderr, "  grow selftest FAIL: the stage failed\n"); fails++; }
    snprintf(path, sizeof path, "%s/out/grown_surface.vmesh", dir);
    if (MeshBin_read_precise_malloc(path, &m) != 0 || m.uv64 == NULL) {
        fprintf(stderr, "  grow selftest FAIL: no grown surface\n");
        fails++;
    } else {
        size_t got = 0, off = 0;
        uint8_t *oa = NULL;
        double worst = 0.0;
        snprintf(path, sizeof path, "%s/out/grown_surface_attr.u8", dir);
        oa = (uint8_t *)ic_read_file(path, m.nv, &got);
        if (oa == NULL) fails++;
        for (v = 0; oa != NULL && v < m.nv; v++) {
            double r = hypot((double)m.verts[v * 3 + 1], (double)m.verts[v * 3 + 2]);
            if (oa[v] != 1 && oa[v] != 3) continue;
            data++;
            worst = fmax(worst, fabs(r - R));
            off += fabs(r - R) > 0.5;
        }
        /* 241 x 241 target nodes; the seed alone holds 81 x 81 */
        if (stop ? (data < 6000 || data > 7200) : (data < (size_t)(0.85 * 241.0 * 241.0) || off > data / 100)) {
            fprintf(stderr, "  grow selftest FAIL: %zu data nodes (seed 6561, target 58081), %zu off the sheet (worst %.2f vox)\n",
                    data, off, worst);
            fails++;
        }
        free(oa);
    }
    if (!stop && !(tile_px > 0.0) && fails == 0) {   /* --reflatten: the grid dump flattened again, by SLIM, into a new id */
        IcConfig c2;
        MeshBinData m2;
        char cfg2[700], *rep = NULL;
        size_t got = 0, pos = 0, neg = 0;
        memset(&m2, 0, sizeof m2);
        memset(&c2, 0, sizeof c2);
        snprintf(cfg2, sizeof cfg2, "%s/regrow.json", dir);
        fp = fopen(cfg2, "wb");
        if (fp == NULL) { fails++; goto done; }
        fprintf(fp, "{\"id\": \"regrown\", \"out_dir\": \"%s/out2\", \"layout\": \"%s/none\", \"raw\": \"%s/none\", "
                "\"axis\": \"%s/axis.csv\", \"voxel_um\": 9.362, \"seed_box\": [0, 0, 1, 1], \"surface_source\": \"%s\", "
                "\"growth_uv\": [-160, -160, 320, 320], \"pool\": [\"%s\"], "
                "\"grow\": {\"margin_vox\": 128, \"rounds\": 6, \"max_reads\": 20, \"reflatten\": true}}\n",
                dir, dir, dir, dir, stem, pool);
        fclose(fp);
        fp = NULL;
        snprintf(path, sizeof path, "%s/out/grown_grid.cgrid", dir);
        if (ic_config(cfg2, &c2) != 0 || ic_stage_reflatten(&c2, path) != 0) {
            fprintf(stderr, "  grow selftest FAIL: --reflatten failed\n");
            fails++;
        } else {
            snprintf(path, sizeof path, "%s/out2/regrown_surface.vmesh", dir);
            if (MeshBin_read_precise_malloc(path, &m2) != 0 || m2.uv64 == NULL || m2.nv != m.nv || m2.nf != m.nf ||
                memcmp(m2.verts, m.verts, m.nv * 3 * sizeof *m.verts) != 0) {
                fprintf(stderr, "  grow selftest FAIL: the reflattened surface is not the grown one (%zu/%zu verts)\n", m2.nv, m.nv);
                fails++;
            } else {
                for (t = 0; t < m2.nf; t++) {   /* one orientation throughout: no face flipped */
                    const double *A = &m2.uv64[(size_t)m2.faces[t * 3] * 2], *B = &m2.uv64[(size_t)m2.faces[t * 3 + 1] * 2];
                    const double *Q = &m2.uv64[(size_t)m2.faces[t * 3 + 2] * 2];
                    const double ar = (B[0] - A[0]) * (Q[1] - A[1]) - (B[1] - A[1]) * (Q[0] - A[0]);
                    pos += ar > 0.0; neg += ar < 0.0;
                }
                if (pos != 0 && neg != 0) {
                    fprintf(stderr, "  grow selftest FAIL: the reflattened UV flips %zu of %zu faces\n", pos < neg ? pos : neg, m2.nf);
                    fails++;
                }
            }
            snprintf(path, sizeof path, "%s/out2/regrown_surface_report.json", dir);
            rep = (char *)ic_read_file(path, 0, &got);
            if (rep == NULL) fails++;
            else {
                char *z = (char *)realloc(rep, got + 1);
                if (z == NULL) fails++;
                else {
                    rep = z;
                    rep[got] = '\0';
                    if (strstr(rep, "\"reflatten_of\": ") == NULL || strstr(rep, "\"id\": \"regrown\"") == NULL ||
                        strstr(rep, "\"method\": \"" IC_GROW_FLATTEN_METHOD "\"") == NULL ||
                        strstr(rep, "\"reflattened\": true") == NULL || strstr(rep, "\"reads\": [") == NULL ||
                        strstr(rep, "\"dedup\": {\"vox\": ") == NULL) {
                        fprintf(stderr, "  grow selftest FAIL: the reflattened report lacks its lines\n");
                        fails++;
                    }
                }
            }
        }
        free(rep);
        MeshBin_dispose(&m2);
        if (c2.pool) free(c2.pool);
    }
    if (!stop) {   /* the pool got its chart index on the first (whole) read (a stopped run reads nothing) */
        char idx[800];
        MeshChartIndex ix;
        snprintf(path, sizeof path, "%s/stage6_repaired.vmesh", pool);
        if (MeshChartIndex_path(path, idx, sizeof idx) != 0 || MeshChartIndex_load(idx, path, &ix) != 0) {
            fprintf(stderr, "  grow selftest FAIL: the pool was not indexed\n");
            fails++;
        } else {
            MeshChartIndex_free(&ix);
        }
    }
    {   /* the data boundary of a 5 x 4 grid with a 3 x 2 data block: its 6 nodes all border non-data */
        CloseupGrid g;
        uint8_t mk[20];
        size_t i = 0, nb = 0;
        memset(&g, 0, sizeof g);
        if (CloseupGrid_alloc(&g, 0.0, 8.0, 0.0, 6.0, 2.0) != 0 || g.ni != 5 || g.nj != 4) fails++;
        else {
            memset(g.state, CLOSEUP_NODE_NONE, 20);
            for (i = 0; i < 20; i++) if (i % 5 >= 1 && i % 5 <= 3 && i / 5 >= 1 && i / 5 <= 2) g.state[i] = CLOSEUP_NODE_DATA;
            nb = ic_data_boundary(&g, mk);
            if (nb != 6 || mk[6] != CLOSEUP_NODE_DATA || mk[0] != CLOSEUP_NODE_NONE) fails++;
            g.state[2] = g.state[16] = CLOSEUP_NODE_DATA;          /* the grid's edge never counts as non-data */
            for (i = 0; i < 20; i++) if (i / 5 == 0 || i / 5 == 3 || i % 5 == 0 || i % 5 == 4) g.state[i] = CLOSEUP_NODE_DATA;
            if (ic_data_boundary(&g, mk) != 0) fails++;
        }
        CloseupGrid_free(&g);
    }
    {   /* the adopted set grows on demand */
        IcAdopted a;
        memset(&a, 0, sizeof a);
        if (ic_adopted_add(&a, 5) != 0 || ic_adopted_add(&a, 100000) != 0 || a.n <= 100000 || !a.on[5] || !a.on[100000] ||
            a.on[6] || ic_adopted_add(&a, -1) != 0) fails++;
        free(a.on);
    }
    if (c.pool) { free(c.pool); c.pool = NULL; }
done:
    MeshBin_dispose(&m);
    free(verts); free(uv); free(faces); free(chart); free(attr);
    fprintf(stderr, "[ink_closeup] grow selftest%s: %s (%zu data nodes)\n",
            stop ? " (STOP before the first read)" : tile_px > 0.0 ? " (tiled)" : "", fails ? "FAIL" : "ok", data);
    if (fails == 0) (void)ves_selftest_clean(dir);
    return fails ? 1 : 0;
}

/* The growth with the F0 bridge end to end: a cylinder sheet (radius 60, axis (y, x) = (128, 128), 3 vox thick) in four
 * 128^3 F0 cubes (z 0..127), a 40 x 80 px seed close-up on it and a pool holding one far chart only (nothing to adopt):
 * the bridge alone must carry the data over most of the target box's rows z 32..96 (280 px wide), on the sheet; reads
 * after one that placed no chart skip the pools (CLOSEUP_GROW_POOL_EVERY); the box reaches z 150, past the cubes, so the
 * bridge stalls there and the last read reads the pools once more before the growth stops. */
static int ic_grow_f0_selftest(void)
{
    char dir[512], path[700], stem[700], cfg[700], pool[700], cubes[700];
    const double R = 60.0;
    const size_t nu = 21, nvv = 41;
    float *verts = NULL;
    double *uv = NULL;
    int32_t *faces = NULL, *chart = NULL;
    uint8_t *attr = NULL, *vol = NULL;
    size_t nv = 0, nf = 0, i = 0, j = 0, data = 0, off = 0, skipped = 0, nreads = 0, v = 0;
    int fails = 0, cy = 0, cx = 0, last_pools = -1;
    IcConfig c;
    FILE *fp = NULL;
    MeshBinData m;
    Arena_T arena = NULL;
    memset(&m, 0, sizeof m); memset(&c, 0, sizeof c);
    if (ves_selftest_dir(dir, sizeof dir, "ink_closeup_grow_f0") != 0) return 1;
    snprintf(cubes, sizeof cubes, "%s/cubes", dir);
    (void)ves_mkdir(cubes);
    vol = (uint8_t *)malloc((size_t)128 * 128 * 128);
    verts = (float *)malloc(nu * nvv * 3 * sizeof *verts);
    uv = (double *)malloc(nu * nvv * 2 * sizeof *uv);
    faces = (int32_t *)malloc((nu - 1) * (nvv - 1) * 6 * sizeof *faces);
    chart = (int32_t *)calloc(nu * nvv, sizeof *chart);
    attr = (uint8_t *)malloc(nu * nvv);
    if (vol == NULL || verts == NULL || uv == NULL || faces == NULL || chart == NULL || attr == NULL) { fails++; goto done; }
    for (cy = 0; cy < 2; cy++)
        for (cx = 0; cx < 2; cx++) {
            size_t z = 0, y = 0, x = 0;
            for (z = 0; z < 128; z++)
                for (y = 0; y < 128; y++)
                    for (x = 0; x < 128; x++) {
                        const double Y = (double)(cy * 128) + (double)y - 128.0, X = (double)(cx * 128) + (double)x - 128.0;
                        vol[(z * 128 + y) * 128 + x] = fabs(hypot(Y, X) - R) <= 1.5 ? 255 : 0;
                    }
            snprintf(path, sizeof path, "%s/z00000_y%05d_x%05d.tif", cubes, cy * 128, cx * 128);
            fails += TiffIO_save(path, vol, 128, 128, 128) != 0;
        }
    /* the seed close-up: u in [0, 40] (theta = u / R), v = z in [24, 104], every 2 px, all observed, chart 0 */
    for (j = 0; j < nvv; j++)
        for (i = 0; i < nu; i++) {
            const double u = 2.0 * (double)i, z = 24.0 + 2.0 * (double)j;
            verts[nv * 3] = (float)z; verts[nv * 3 + 1] = (float)(128.0 + R * sin(u / R));
            verts[nv * 3 + 2] = (float)(128.0 + R * cos(u / R));
            uv[nv * 2] = u; uv[nv * 2 + 1] = z;
            attr[nv] = 1;
            nv++;
        }
    for (j = 0; j + 1 < nvv; j++)
        for (i = 0; i + 1 < nu; i++) {
            const int32_t a = (int32_t)(j * nu + i);
            faces[nf * 3] = a; faces[nf * 3 + 1] = a + 1; faces[nf * 3 + 2] = a + (int32_t)nu + 1; nf++;
            faces[nf * 3] = a; faces[nf * 3 + 1] = a + (int32_t)nu + 1; faces[nf * 3 + 2] = a + (int32_t)nu; nf++;
        }
    snprintf(stem, sizeof stem, "%s/seed", dir);
    snprintf(path, sizeof path, "%s_surface.vmesh", stem);
    fails += MeshBin_write_uv64(path, verts, nv, faces, nf, uv) != 0;
    snprintf(path, sizeof path, "%s_surface_attr.u8", stem);
    fails += ic_write_file(path, attr, nv) != 0;
    snprintf(path, sizeof path, "%s_surface_chart.i32", stem);
    fails += ic_write_file(path, chart, nv * sizeof *chart) != 0;
    {   /* the pool: one triangle far from everything (z 900) */
        float pv[9] = { 900.0f, 10.0f, 10.0f, 900.0f, 10.0f, 12.0f, 902.0f, 10.0f, 10.0f };
        double pu[6] = { 0.0, 0.0, 2.0, 0.0, 0.0, 2.0 };
        int32_t pf[3] = { 0, 1, 2 }, pc[3] = { 0, 0, 0 };
        snprintf(pool, sizeof pool, "%s/pool", dir);
        (void)ves_mkdir(pool);
        snprintf(path, sizeof path, "%s/stage6_repaired.vmesh", pool);
        fails += MeshBin_write_uv64(path, pv, 3, pf, 1, pu) != 0;
        snprintf(path, sizeof path, "%s/stage6_repaired_chart.i32", pool);
        fails += ic_write_file(path, pc, sizeof pc) != 0;
    }
    snprintf(path, sizeof path, "%s/axis.csv", dir);
    fp = fopen(path, "wb");
    if (fp == NULL) { fails++; goto done; }
    fprintf(fp, "# test axis\nZ_l0,y_l0,x_l0\n-512,128,128\n0,128,128\n512,128,128\n1024,128,128\n");
    fclose(fp);
    snprintf(cfg, sizeof cfg, "%s/grow.json", dir);
    fp = fopen(cfg, "wb");
    if (fp == NULL) { fails++; goto done; }
    fprintf(fp, "{\"id\": \"grownf0\", \"out_dir\": \"%s/out\", \"layout\": \"%s/none\", \"raw\": \"%s/none\", "
            "\"axis\": \"%s/axis.csv\", \"voxel_um\": 9.362, \"seed_box\": [0, 0, 1, 1], \"surface_source\": \"%s\", "
            "\"growth_uv\": [-120, 24, 160, 150], \"pool\": [\"%s\"], \"complete_gaps\": true, "
            "\"grow\": {\"margin_vox\": 64, \"rounds\": 2, \"max_reads\": 40, \"reflatten\": false, \"f0\": [\"%s\"], "
            "\"f0_depth_px\": 100000, \"f0_rings\": 8, \"f0_passes\": 2}}\n", dir, dir, dir, dir, stem, pool, cubes);
    fclose(fp);
    if (ic_config(cfg, &c) != 0 || !c.has_grow || c.nf0 != 1) { fails++; goto done; }
    if (ic_stage_grow(&c) != 0) { fprintf(stderr, "  grow f0 selftest FAIL: the stage failed\n"); fails++; goto done; }
    snprintf(path, sizeof path, "%s/out/grownf0_surface.vmesh", dir);
    if (MeshBin_read_precise_malloc(path, &m) != 0) { fails++; goto done; }
    {
        size_t got = 0;
        uint8_t *oa = NULL;
        snprintf(path, sizeof path, "%s/out/grownf0_surface_attr.u8", dir);
        oa = (uint8_t *)ic_read_file(path, m.nv, &got);
        if (oa == NULL) { fails++; goto done; }
        for (v = 0; v < m.nv; v++) {
            const double z = m.verts[v * 3], r = hypot((double)m.verts[v * 3 + 1] - 128.0, (double)m.verts[v * 3 + 2] - 128.0);
            if ((oa[v] != 1 && oa[v] != 3) || z < 32.0 || z > 96.0) continue;
            data++;
            off += fabs(r - R) > 1.2;
        }
        free(oa);
    }
    {   /* the reads: some bridge-only, the last one with the pools */
        const char *err = NULL;
        const JsonValue *rep = NULL, *reads = NULL;
        arena = Arena_new();
        snprintf(path, sizeof path, "%s/out/grownf0_surface_report.json", dir);
        rep = arena ? Json_parse_file(arena, path, &err) : NULL;
        reads = rep ? Json_object_get(rep, "reads") : NULL;
        nreads = Json_array_len(reads);
        for (i = 0; i < nreads; i++) {
            const int pr = (int)Json_member_double(Json_array_get(reads, i), "pools_read", -1.0);
            skipped += pr == 0;
            last_pools = pr;
        }
    }
    /* the box rows 32..96: 33 rows x 141 columns = 4,653 nodes */
    if (data < (size_t)(0.8 * 4653.0) || off > data / 100 || skipped == 0 || last_pools != 1) {
        fprintf(stderr, "  grow f0 selftest FAIL: %zu data nodes in rows 32..96 (want >= 3722), %zu off the sheet; %zu of %zu "
                "reads bridge-only, the last read pools %d\n", data, off, skipped, nreads, last_pools);
        fails++;
    }
    {   /* measured progress: the same growth with min_gain 0.3 ends earlier, by "min_gain", its last read with the pools */
        const char *err = NULL, *why = NULL;
        const JsonValue *rep = NULL, *reads = NULL;
        size_t nreads2 = 0;
        int last2 = -1;
        IcConfig c2;
        memset(&c2, 0, sizeof c2);
        snprintf(cfg, sizeof cfg, "%s/grow2.json", dir);
        fp = fopen(cfg, "wb");
        if (fp == NULL) { fails++; goto done; }
        fprintf(fp, "{\"id\": \"grownf0b\", \"out_dir\": \"%s/out2\", \"layout\": \"%s/none\", \"raw\": \"%s/none\", "
                "\"axis\": \"%s/axis.csv\", \"voxel_um\": 9.362, \"seed_box\": [0, 0, 1, 1], \"surface_source\": \"%s\", "
                "\"growth_uv\": [-120, 24, 160, 150], \"pool\": [\"%s\"], \"complete_gaps\": true, "
                "\"grow\": {\"margin_vox\": 64, \"rounds\": 2, \"max_reads\": 40, \"reflatten\": false, \"f0\": [\"%s\"], "
                "\"f0_depth_px\": 100000, \"f0_rings\": 8, \"f0_passes\": 2, \"min_gain\": 0.3}}\n", dir, dir, dir, dir, stem,
                pool, cubes);
        fclose(fp);
        if (ic_config(cfg, &c2) != 0 || ic_stage_grow(&c2) != 0) { fprintf(stderr, "  grow f0 selftest FAIL: min_gain run\n"); fails++; }
        if (c2.pool) free(c2.pool);
        snprintf(path, sizeof path, "%s/out2/grownf0b_surface_report.json", dir);
        rep = arena ? Json_parse_file(arena, path, &err) : NULL;
        reads = rep ? Json_object_get(rep, "reads") : NULL;
        nreads2 = Json_array_len(reads);
        if (nreads2 > 0) last2 = (int)Json_member_double(Json_array_get(reads, nreads2 - 1), "pools_read", -1.0);
        why = rep ? Json_as_string(Json_object_get(Json_object_get(rep, "grow"), "stop_reason")) : NULL;
        if (nreads2 == 0 || nreads2 >= nreads || last2 != 1 || why == NULL || strcmp(why, "min_gain") != 0) {
            fprintf(stderr, "  grow f0 selftest FAIL: min_gain 0.3: %zu reads (want fewer than %zu), the last pools %d, "
                    "stop %s (want min_gain)\n", nreads2, nreads, last2, why ? why : "?");
            fails++;
        }
    }
done:
    if (c.pool) free(c.pool);
    MeshBin_dispose(&m);
    if (arena != NULL) Arena_dispose(&arena);
    free(vol); free(verts); free(uv); free(faces); free(chart); free(attr);
    fprintf(stderr, "[ink_closeup] grow f0 selftest: %s (%zu data nodes, %zu of %zu reads bridge-only)\n",
            fails ? "FAIL" : "ok", data, skipped, nreads);
    if (fails == 0) (void)ves_selftest_clean(dir);
    return fails ? 1 : 0;
}

/* ic_pool_tile on a two-chart pool: a 10 x 10 vertex grid at 10 px (chart 7) and a lone triangle (chart 3), tiled at
 * 32 px: 9 + 1 tiles, every face kept once on its own tile's vertices, positions and UV carried together, the parents'
 * ids kept, the keys distinct and the same on a second tiling; then the tile key set. */
static int ic_tile_selftest(void)
{
    IcPool q;
    int64_t *keys[2] = { NULL, NULL };
    const size_t nv = 103, nf = 2 * 81 + 1;
    size_t i = 0, j = 0, t = 0, f = 0, tiles7 = 0;
    int fails = 0, pass = 0, kk = 0;
    for (pass = 0; pass < 2; pass++) {
        memset(&q, 0, sizeof q);
        q.verts = (float *)malloc(nv * 3 * sizeof *q.verts); q.nrm = (float *)malloc(nv * 3 * sizeof *q.nrm);
        q.inward = (float *)malloc(nv * 3 * sizeof *q.inward); q.uv = (double *)malloc(nv * 2 * sizeof *q.uv);
        q.faces = (int32_t *)malloc(nf * 3 * sizeof *q.faces); q.face_chart = (int32_t *)malloc(nf * sizeof *q.face_chart);
        q.chart_id = (int32_t *)malloc(2 * sizeof *q.chart_id);
        q.vofs = (size_t *)malloc(3 * sizeof *q.vofs); q.fofs = (size_t *)malloc(3 * sizeof *q.fofs);
        if (!q.verts || !q.nrm || !q.inward || !q.uv || !q.faces || !q.face_chart || !q.chart_id || !q.vofs || !q.fofs) {
            ic_pool_free(&q);
            fails++;
            break;
        }
        for (j = 0; j < 10; j++)
            for (i = 0; i < 10; i++) {
                const size_t v = j * 10 + i;
                q.verts[v * 3] = 0.0f; q.verts[v * 3 + 1] = (float)(10 * j); q.verts[v * 3 + 2] = (float)(10 * i);
                q.nrm[v * 3] = 1.0f; q.nrm[v * 3 + 1] = 0.0f; q.nrm[v * 3 + 2] = 0.0f;
                q.inward[v * 3] = 1.0f; q.inward[v * 3 + 1] = 0.0f; q.inward[v * 3 + 2] = 0.0f;
                q.uv[v * 2] = 10.0 * (double)i + 5.0; q.uv[v * 2 + 1] = 10.0 * (double)j - 3.0;
            }
        for (j = 0, f = 0; j < 9; j++)
            for (i = 0; i < 9; i++, f += 2) {
                const int32_t a = (int32_t)(j * 10 + i);
                q.faces[f * 3] = a; q.faces[f * 3 + 1] = a + 1; q.faces[f * 3 + 2] = a + 11;
                q.faces[f * 3 + 3] = a; q.faces[f * 3 + 4] = a + 11; q.faces[f * 3 + 5] = a + 10;
                q.face_chart[f] = q.face_chart[f + 1] = 7;
            }
        for (i = 100; i < 103; i++) {
            q.verts[i * 3] = 5.0f; q.verts[i * 3 + 1] = 500.0f + (i == 102 ? 10.0f : 0.0f);
            q.verts[i * 3 + 2] = 500.0f + (i == 101 ? 10.0f : 0.0f);
            for (kk = 0; kk < 3; kk++) q.nrm[i * 3 + (size_t)kk] = q.inward[i * 3 + (size_t)kk] = kk == 0 ? 1.0f : 0.0f;
            q.uv[i * 2] = i == 101 ? 10.0 : 0.0; q.uv[i * 2 + 1] = i == 102 ? 10.0 : 0.0;
        }
        q.faces[162 * 3] = 100; q.faces[162 * 3 + 1] = 101; q.faces[162 * 3 + 2] = 102; q.face_chart[162] = 3;
        q.chart_id[0] = 7; q.chart_id[1] = 3;
        q.vofs[0] = 0; q.vofs[1] = 100; q.vofs[2] = 103; q.fofs[0] = 0; q.fofs[1] = 162; q.fofs[2] = 163;
        q.p.nv = nv; q.p.nf = nf; q.p.nch = 2;
        q.p.verts = q.verts; q.p.nrm = q.nrm; q.p.inward = q.inward; q.p.uv = q.uv; q.p.faces = q.faces;
        q.p.vofs = q.vofs; q.p.fofs = q.fofs;
        if (ic_pool_tile(&q, 32.0, &keys[pass]) != 0 || keys[pass] == NULL) { fails++; ic_pool_free(&q); break; }
        if (q.p.nch != 10 || q.p.nf != nf || q.p.nv <= nv || q.vofs[10] != q.p.nv || q.fofs[10] != nf) fails++;
        for (t = 0, tiles7 = 0; t < q.p.nch && t < 10; t++) {
            tiles7 += q.chart_id[t] == 7;
            if (q.chart_id[t] != 7 && q.chart_id[t] != 3) fails++;
            if (q.fofs[t + 1] <= q.fofs[t]) fails++;                                 /* no empty tile */
            for (f = q.fofs[t]; f < q.fofs[t + 1]; f++) {
                if (q.face_chart[f] != q.chart_id[t]) fails++;
                for (kk = 0; kk < 3; kk++) {
                    const size_t v = (size_t)q.faces[f * 3 + (size_t)kk];
                    if (v < q.vofs[t] || v >= q.vofs[t + 1]) { fails++; continue; }
                    if (q.chart_id[t] == 7 && (fabs((double)q.verts[v * 3 + 2] - (q.uv[v * 2] - 5.0)) > 1e-6 ||
                                               fabs((double)q.verts[v * 3 + 1] - (q.uv[v * 2 + 1] + 3.0)) > 1e-6)) fails++;
                }
            }
            for (j = 0; j < t; j++) if (keys[pass][j] == keys[pass][t]) fails++;
        }
        if (tiles7 != 9) fails++;
        ic_pool_free(&q);
    }
    if (fails == 0 && keys[0] != NULL && keys[1] != NULL)
        for (t = 0; t < 10; t++) if (keys[0][t] != keys[1][t]) fails++;
    free(keys[0]); free(keys[1]);
    {   /* the key set: 20,000 keys through several rehashes, members found, non-members not, duplicates not counted */
        IcKeySet ks;
        int64_t k = 0;
        memset(&ks, 0, sizeof ks);
        for (k = 0; k < 20000; k++) if (ic_keyset_add(&ks, k * 7919 + ((int64_t)k << 30)) != 0) fails++;
        if (ic_keyset_add(&ks, 7919 + ((int64_t)1 << 30)) != 0 || ks.n != 20000 || ic_keyset_add(&ks, -1) == 0) fails++;
        for (k = 0; k < 20000; k++) if (!ic_keyset_has(&ks, k * 7919 + ((int64_t)k << 30))) { fails++; break; }
        if (ic_keyset_has(&ks, 3) || ic_keysets_total(&ks, 1) != 20000) fails++;
        free(ks.slot);
    }
    fprintf(stderr, "[ink_closeup] pool tiling selftest: %s\n", fails ? "FAIL" : "ok");
    return fails ? 1 : 0;
}

static int ic_selftest_rebuild_material(void);
static int ic_material_selftest(void);
static int ic_projection_selftest(void);
static int ic_align_global_selftest(void);
static int ic_selftest(void)
{
    if (ic_align_global_selftest() || ic_projection_selftest() || CloseupFiber_selftest() || CloseupFiberVolume_selftest() || ic_selftest_rebuild_material() ||
        ic_material_selftest() || AsmContacts_selftest()) return 1;
    int a = CloseupSurface_selftest(), b = Sparse_selftest(), d = Json_read_selftest(), e = AsmAxis_selftest();
    int h = CloseupAdopt_selftest(), l = InkLetters_selftest();
    int ct = CloseupCtTrack_selftest() || CloseupImport_selftest() || CloseupSeed_selftest();
    int ac = ic_selftest_adopt_config();
    int tiles = Rawtex_tiled_selftest();
    int speed = KDTree_selftest() != 0 || MeshNormals_selftest() != 0;   /* the completion speed paths (2026-10-30) */
    /* A full-line selection includes small islands, regardless of a neighbour's area. */
    if (ic_seed_selected(500000,32000000,0) || !ic_seed_selected(500000,32000000,1) ||
        ic_seed_selected(0,32000000,1) || !ic_seed_selected(1,32000000,1) ||
        !ic_seed_selected(1000,2000,0)) ac = 1;
    int rp = 0;
    {   /* PHerc1203ztfix2 line 7: above the default cap, with no lost physical pixels. */
        double uv[4] = {0, 0, 41624 - CLOSEUP_RASTER_MARGIN_PX, 7306 - CLOSEUP_RASTER_MARGIN_PX};
        RawtexWindow win;
        RawtexPlan plan;
        size_t test_w = 0, test_h = 0;
        rp = ic_render_plan((RawtexUv){ NULL, uv }, 2, &win, &test_w, &test_h) != 0 || test_w != 41624 || test_h != 7306;
        rp |= Rawtex_plan_field((RawtexUv){ NULL, uv }, 2, 1, 1, 0, &win, &plan) != 0 || plan.ok;
        uv[2] = INFINITY;
        rp |= ic_render_plan((RawtexUv){ NULL, uv }, 2, &win, &test_w, &test_h) == 0;
        uv[2] = (double)INT_MAX;
        rp |= ic_render_plan((RawtexUv){ NULL, uv }, 2, &win, &test_w, &test_h) == 0;
    }
    {   /* letters.json numbers: finite values printed, nan / inf as null */
        char nb[32];
        if (strcmp(ic_json_num(nb, sizeof nb, NAN), "null") != 0 || strcmp(ic_json_num(nb, sizeof nb, INFINITY), "null") != 0 ||
            strcmp(ic_json_num(nb, sizeof nb, 12.04), "12.0") != 0) ac = 1;
    }
    fprintf(stderr, "[ink_closeup] selftest: closeup_surface %s, closeup_adopt %s, ink_letters %s, sparse %s, json %s, "
            "axis %s\n", a ? "FAIL" : "ok", h ? "FAIL" : "ok", l ? "FAIL" : "ok", b ? "FAIL" : "ok", d ? "FAIL" : "ok",
            e ? "FAIL" : "ok");
    fprintf(stderr, "[ink_closeup] ct tracking selftest: %s; adopt config overrides: %s\n", ct ? "FAIL" : "ok",
            ac ? "FAIL" : "ok");
    fprintf(stderr, "[ink_closeup] full-line render plan: %s; parallel kd-tree and pool normals: %s\n", rp ? "FAIL" : "ok",
            speed ? "FAIL" : "ok");
    if (speed) rp = 1;
    {
        int mci = MeshChartIndex_selftest(), nrt = ic_near_selftest(), chn = ic_chain_selftest();
        fprintf(stderr, "[ink_closeup] chart index: %s; pool near region: %s; completion stage files: %s\n",
                mci ? "FAIL" : "ok", nrt ? "FAIL" : "ok", chn ? "FAIL" : "ok");
        if (mci || nrt || chn) rp = 1;
    }
    {
        int cp = ChartPool_selftest(), gr = ic_grow_selftest(0.0, 0), tl = ic_tile_selftest(), gt = ic_grow_selftest(48.0, 0);
        int fb = CloseupF0Bridge_selftest(), gs = ic_grow_selftest(0.0, 1), gf = ic_grow_f0_selftest();
        int ru = ic_reuv_selftest();
        fprintf(stderr, "[ink_closeup] chart pool: %s; frontier growth: %s; pool tiling: %s; tiled growth: %s; F0 bridge: %s; "
                "graceful stop: %s; bridged growth: %s; fibre-frame reuv: %s\n", cp ? "FAIL" : "ok", gr ? "FAIL" : "ok",
                tl ? "FAIL" : "ok", gt ? "FAIL" : "ok", fb ? "FAIL" : "ok", gs ? "FAIL" : "ok", gf ? "FAIL" : "ok", ru ? "FAIL" : "ok");
        if (cp || gr || tl || gt || fb || gs || gf || ru) rp = 1;
    }
    return a || b || d || e || h || l || ct || ac || rp || tiles ? 1 : 0;
}

/* --build-index: write <vmesh>.chartidx, so a close-up reads only the charts near its box.  Exit 0 written, 2 the mesh
 * is not chart-contiguous (its readers keep the whole-mesh read), 1 error. */
static int ic_build_index(const char *vmesh, const char *chart)
{
    char ipath[2700];
    double t0 = ves_clock_sec();
    MeshChartIndex ix;
    int rc = 0;
    if (MeshChartIndex_path(vmesh, ipath, sizeof ipath) != 0) return 1;
    rc = MeshChartIndex_build(vmesh, chart, ipath);
    if (rc == 1) { fprintf(stderr, "[ink_closeup] %s is not chart-contiguous: no index\n", vmesh); return 2; }
    if (rc != 0 || MeshChartIndex_load(ipath, vmesh, &ix) != 0) {
        fprintf(stderr, "ERROR: cannot index %s with %s\n", vmesh, chart);
        return 1;
    }
    fprintf(stderr, "[ink_closeup] wrote %s: %zu chart runs over %llu verts, %llu faces (%.1fs)\n", ipath, ix.n,
            (unsigned long long)ix.mesh.nv, (unsigned long long)ix.mesh.nf, ves_clock_sec() - t0);
    MeshChartIndex_free(&ix);
    return 0;
}

#include "ink_closeup_reparameterize.inc"
#include "ink_closeup_complete_fixed.inc"
#include "ink_closeup_rebuild.inc"

int main(int argc, char **argv)
{
    IcConfig c;
    const char *stop = "render";
    char probe[2600];
    int i = 0;
    ves_prof_start_from_env();         /* VES_PROF_OUT=<file>: the sampling profiler (diagnostics only) */
    if (argc >= 2 && strcmp(argv[1], "--selftest") == 0) return ic_selftest();
    if (argc == 2 && strcmp(argv[1], "--reparameterize-selftest") == 0) return ic_reparameterize_selftest();
    if (argc == 2 && strcmp(argv[1], "--slim-selftest") == 0) return SlimRefine_selftest();
#ifdef VES_CLOSEUP_LOCAL_GLOBAL
    if (argc == 2 && strcmp(argv[1], "--local-global-selftest") == 0) return CloseupLocalGlobal_selftest();
#endif
    if (argc == 2 && strcmp(argv[1], "--physical-align-selftest") == 0) return ic_physical_alignment_selftest();
    if (argc == 2 && strcmp(argv[1], "--align-module-selftest") == 0) return CloseupAlign_selftest();
    if (argc == 2 && strcmp(argv[1], "--projection-selftest") == 0) return ic_projection_selftest();
    if (argc == 2 && strcmp(argv[1], "--global-contact-selftest") == 0) return ic_align_global_selftest();
    if (argc == 2 && strcmp(argv[1], "--fiber-measure-selftest") == 0) return CloseupFiber_selftest();
    if (argc == 2 && strcmp(argv[1], "--fiber-volume-selftest") == 0) return CloseupFiberVolume_selftest();
    if (argc == 2 && strcmp(argv[1], "--rebuild-material-selftest") == 0) return ic_selftest_rebuild_material();
    if (argc == 2 && strcmp(argv[1], "--rebuild-core-selftest") == 0) return CloseupSurface_rebuild_selftest();
    if (argc == 2 && strcmp(argv[1], "--material-contact-selftest") == 0) return ic_material_selftest() || AsmContacts_selftest();
    if ((argc == 5 || argc == 6) && strcmp(argv[1], "--audit-material") == 0)
        return ic_material_audit(argv[2],argv[3],argv[4],argc == 6 ? argv[5] : NULL);
    if (argc == 2 && strcmp(argv[1], "--seed-selftest") == 0) return CloseupSeed_selftest();
    if (argc == 3 && strcmp(argv[1], "--seed-rendered") == 0) return CloseupSeed_write(argv[2]) != 0;
    if (argc == 2 && strcmp(argv[1], "--import-selftest") == 0) return CloseupImport_selftest();
    if ((argc == 5 || argc == 6) && strcmp(argv[1], "--import-surface") == 0)
        return CloseupImport_write(argv[2], argv[3], argv[4], argc == 6 ? argv[5] : NULL) != 0;
    if (argc == 2 && strcmp(argv[1], "--ct-track-selftest") == 0) return CloseupCtTrack_selftest();
    if (argc == 4 && strcmp(argv[1], "--build-index") == 0) return ic_build_index(argv[2], argv[3]);
    if ((argc == 4 || argc == 10) && strcmp(argv[1], "--pool-from-store") == 0) {
        long box[6];
        ChartPoolStats ps;
        for (i = 0; argc == 10 && i < 6; i++) box[i] = strtol(argv[4 + i], NULL, 10);
        i = ChartPool_from_store(argv[2], argv[3], argc == 10 ? box : NULL, &ps);
        if (i == 1) fprintf(stderr, "[ink_closeup] %s: no chart qualifies (%zu cubes)\n", argv[2], ps.cubes);
        return i == 0 ? 0 : 1;
    }
    if (argc < 2) {
        fprintf(stderr, "usage: ink_closeup <config.json> [--stop-after surface|render]\n"
                        "       ink_closeup <config.json> --letters <ink.tif> base|rect\n"
                        "       ink_closeup <config.json> --render-rect\n"
                        "       ink_closeup <config.json> --render-sweep     (CLOSEUP_SWEEP_DEPTH layers: the depth sweep)\n"
                        "       ink_closeup <config.json> --ct-track-audit   (CT-only traces; no surface modification)\n"
                        "       ink_closeup <config.json> --reflatten <source id>_grid.cgrid   (a growth's last steps again)\n"
                        "       ink_closeup <config.json> --reuv <fiber_flatten out dir>   (the surface into the fibre frame)\n"
                        "       ink_closeup <job.json> --complete-chain     (the fiber_complete stage: every pass in one process)\n"
                        "       ink_closeup --build-index <mesh.vmesh> <chart.i32>   (the <mesh>.chartidx pager)\n"
                        "       ink_closeup --import-surface <reading.vmesh> <chart.i32> <new output stem> [charts.csv]\n"
                        "       ink_closeup --pool-from-store <chart store dir> <out dir> [z0 z1 y0 y1 x0 x1]   (an adoption pool)\n"
                        "       ink_closeup --selftest\n");
        return 1;
    }
    if (argc == 3 && strcmp(argv[2], "--complete-chain") == 0) return ic_stage_complete_chain(argv[1]) != 0;
    if (argc == 3 && strcmp(argv[2], "--ct-track-audit") == 0) {
        char stem[2600];
        if (ic_config(argv[1], &c) != 0 || !c.ct_track_configured || !c.surface_source[0]) return 1;
        snprintf(stem, sizeof stem, "%s/%s", c.out_dir, c.id);
        return CloseupCtTrack_audit(c.surface_source, c.raw, stem, c.ct_track_box, c.ct_track_seed_z) == 0 ? 0 : 1;
    }
    if (argc == 4 && strcmp(argv[2], "--reuv") == 0) {
        if (ic_config(argv[1], &c) != 0 || !c.surface_source[0]) return 1;
        return ic_stage_reuv(&c, argv[3]);
    }
    if (argc == 3 && strcmp(argv[2], "--measure-fibers") == 0) {
        if (ic_config(argv[1], &c) != 0 || !c.surface_source[0]) return 1;
        return ic_stage_measure_fibers(&c);
    }
    if (argc == 3 && (!strcmp(argv[2], "--reparameterize") || !strcmp(argv[2], "--reparameterize-conformal") ||
                     !strcmp(argv[2], "--reparameterize-cuts"))) {
        if (ic_config(argv[1], &c) != 0 || !c.surface_source[0]) return 1;
        ic_reparameterize_mode = !strcmp(argv[2], "--reparameterize-cuts") ? 2 : !strcmp(argv[2], "--reparameterize-conformal");
        return ic_stage_reparameterize(&c);
    }
    if ((argc == 3 || argc == 4) && (!strcmp(argv[2], "--repair-completed") ||
                                  !strcmp(argv[2], "--repair-completed-cuts") ||
                                  !strcmp(argv[2], "--repair-completed-provisional") ||
                                  !strcmp(argv[2], "--repair-completed-provisional-cuts"))) {
        if (ic_config(argv[1], &c) != 0 || !c.surface_source[0]) return 1;
        ic_reparameterize_mode = strstr(argv[2], "-cuts") ? 2 : 0;
        return ic_stage_repair_completed(&c,argc == 4 ? argv[3] : NULL,strstr(argv[2], "-provisional") != NULL);
    }
    if (argc == 3 && (!strcmp(argv[2], "--rebuild-completed") || !strcmp(argv[2], "--rebuild-completed-cuts"))) {
        if (ic_config(argv[1], &c) != 0 || !c.surface_source[0]) return 1;
        ic_reparameterize_mode = strstr(argv[2], "-cuts") ? 2 : 0;
        return ic_stage_rebuild_completed(&c);
    }
    if (argc == 4 && strcmp(argv[2], "--reflatten") == 0) {
        if (ic_config(argv[1], &c) != 0) return 1;
        return ic_stage_reflatten(&c, argv[3]);
    }
    if (argc == 5 && strcmp(argv[2], "--letters") == 0) {
        if (ic_config(argv[1], &c) != 0) return 1;
        if (strcmp(argv[4], "base") != 0 && strcmp(argv[4], "rect") != 0) { fprintf(stderr, "base or rect\n"); return 1; }
        return ic_stage_letters(&c, argv[4], argv[3]);
    }
    if (argc == 3 && strcmp(argv[2], "--render-rect") == 0) {
        char stem[512];
        if (ic_config(argv[1], &c) != 0) return 1;
        snprintf(stem, sizeof stem, "%s_rect", c.id);
        return ic_stage_render(&c, stem, stem, SHEET_LAYERS_DEPTH, 1);
    }
    if (argc == 3 && strcmp(argv[2], "--render-sweep") == 0) {
        char stem[512];
        if (ic_config(argv[1], &c) != 0) return 1;
        snprintf(stem, sizeof stem, "%s_sweep", c.id);
        return ic_stage_render(&c, c.id, stem, CLOSEUP_SWEEP_DEPTH, 0);
    }
    for (i = 2; i < argc; i++) {
        if (strcmp(argv[i], "--stop-after") == 0 && i + 1 < argc) stop = argv[++i];
        else { fprintf(stderr, "unknown option %s\n", argv[i]); return 1; }
    }
    if (strcmp(stop, "surface") != 0 && strcmp(stop, "render") != 0) { fprintf(stderr, "bad --stop-after %s\n", stop); return 1; }
    if (ic_config(argv[1], &c) != 0) return 1;
    snprintf(probe, sizeof probe, "%s/%s_surface.vmesh", c.out_dir, c.id);
    if (MeshBin_looks_complete(probe)) fprintf(stderr, "[ink_closeup] surface stage: resuming from %s\n", probe);
    else if (c.surface_source[0] && c.has_grow) { if (ic_stage_grow(&c) != 0) return 1; }
    else if (c.surface_source[0]) { if (ic_stage_complete_surface(&c) != 0) return 1; }
    else if (ic_stage_surface(&c) != 0) return 1;
    if (strcmp(stop, "surface") == 0) return 0;
    snprintf(probe, sizeof probe, "%s/%s_render_report.json", c.out_dir, c.id);
    {
        FILE *f = fopen(probe, "rb");
        if (f != NULL) { fclose(f); fprintf(stderr, "[ink_closeup] render stage: done (%s)\n", probe); return 0; }
    }
    return ic_stage_render(&c, c.id, c.id, SHEET_LAYERS_DEPTH, 1);
}
