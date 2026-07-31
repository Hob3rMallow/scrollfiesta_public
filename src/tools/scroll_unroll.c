/*
 * scroll_unroll.c -- whole-grid back half: concatenate a scroll_whole placed
 * dir into one flat piece set, then chain the v2 stages over it, baking a
 * strip TIF after every stage (the user reviews TIFs, never meshes):
 *
 *   step1_ribbon   base bake (v1)
 *   step2_join     UV-weld (z- AND u-seams) + banded relax        [v4]
 *   step3_overlap  SeamOwn ownership: multicut/quilting de-blur
 *   step4_snap     dark-vertex snap (GLOBAL mode after the weld)
 *   step5_relax    final light banded relax
 *
 * v4 rationale: the trim-inset-0 dump kills the manufactured 2.75-vox z-seam
 * gap, so the standard weld gates now close the "cut into 4 strips" v-seams
 * -- weld+flatten runs FIRST (step2_join), before any ownership/snap work,
 * exactly as on the u-atlas gutters. Post-weld invariants: SeamOwn exempts
 * mesh-adjacent cross-cube pairs; snap runs whole-mesh (global_mode).
 *
 * Per stage: <out>/<id>_<step>_{rawtex,diagclass}.tif (+ _strip.tif, .png
 * when narrow), <id>_<step>_rawtex_preview.png / _diagclass_preview.png /
 * _seamzoom.png (streamed previews), metrics row in the end summary, and a
 * block in <id>_pipeline_stats.json. All stages share one contrast window
 * (from step1) so stage TIFs difference cleanly, and one RAW CubeTable.
 *
 * Usage:
 *   scroll_unroll <placed_dir> <out_dir> --raw <cubes_RAW_dir>
 *       [--steps S] [--id X] [--du F] [--dv F] [--band-cols N]
 *       [--raw-chunk N] [--range F] [--nsteps N] [--threads N]
 *       [--stretch-ratio F] [--stretch-floor F] [--max-edge F]
 *       [--strip-w N] [--dark-thresh N] [--seam-top N]
 *       [--no-diag] [--no-preview]
 *   scroll_unroll --selftest
 *
 *   --steps S: stage subset as a digit string (default "12345"; stages not
 *   built into this binary are skipped with a warning). NOTE: the v1 flag
 *   "--steps N" (samples along the normal) is now --nsteps.
 */
#include "../common/ves_platform.h"

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/arena.h"
#include "../common/raw_sample.h"
#include "../common/tiff_io.h"
#include "../unroll/export_atlas.h"
#include "../unroll/piece_set.h"
#include "../unroll/relax_grid.h"
#include "../unroll/scaffold.h"
#include "../unroll/scroll_raster.h"
#include "../unroll/seam_own.h"
#include "../unroll/snap_grid.h"
#include "../unroll/strip_metrics.h"
#include "../unroll/strip_preview.h"
#include "../unroll/tifxyz_export.h"
#include "../unroll/uv_weld.h"

#define HAVE_STAGE2 1
#define HAVE_STAGE3 1
#define HAVE_STAGE4 1
#define HAVE_STAGE5 1

static FILE *g_log = NULL;

static void write_json_string(FILE *fp, const char *value)
{
    fputc('"', fp);
    for (const unsigned char *p = (const unsigned char *)value; *p; p++) {
        switch (*p) {
        case '"': fputs("\\\"", fp); break;
        case '\\': fputs("\\\\", fp); break;
        case '\b': fputs("\\b", fp); break;
        case '\f': fputs("\\f", fp); break;
        case '\n': fputs("\\n", fp); break;
        case '\r': fputs("\\r", fp); break;
        case '\t': fputs("\\t", fp); break;
        default:
            if (*p < 0x20) fprintf(fp, "\\u%04x", (unsigned)*p);
            else fputc(*p, fp);
        }
    }
    fputc('"', fp);
}

static int write_json_string_selftest(void)
{
    static const char input[] = "C:\\raw\\x\"\n";
    static const char expected[] = "\"C:\\\\raw\\\\x\\\"\\n\"";
    char got[sizeof expected + 8];
    FILE *fp = tmpfile();
    if (fp == NULL) return 1;
    write_json_string(fp, input);
    fflush(fp);
    rewind(fp);
    size_t n = fread(got, 1, sizeof got - 1, fp);
    got[n] = '\0';
    fclose(fp);
    int fail = strcmp(got, expected) != 0;
    fprintf(stderr, "write_json_string_selftest: %s\n",
            fail ? "FAIL" : "PASS");
    return fail;
}

static void logf_both(const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    if (g_log != NULL) {
        va_start(ap, fmt);
        vfprintf(g_log, fmt, ap);
        va_end(ap);
        fflush(g_log);
    }
}

static double xyz_color_component(double value, double lo, double hi)
{
    if (!(hi > lo)) return 0.5;
    double q = (value - lo) / (hi - lo);
    if (q < 0.0) q = 0.0;
    if (q > 1.0) q = 1.0;
    return q;
}

static int write_xyz_face_category_obj(const char *path,
                                       const char *group_name,
                                       const PieceSet *ps,
                                       const uint8_t *face_keep,
                                       int want_kept,
                                       int flattened_uv,
                                       const double xyz_lo[3],
                                       const double xyz_hi[3],
                                       size_t *out_vertices,
                                       size_t *out_faces)
{
    if (path == NULL || group_name == NULL || ps == NULL ||
        face_keep == NULL || xyz_lo == NULL || xyz_hi == NULL)
        return -1;
    int32_t *remap = (int32_t *)malloc(ps->nv * sizeof(int32_t));
    if (remap == NULL) return -1;
    for (size_t i = 0; i < ps->nv; i++) remap[i] = -1;

    size_t nfaces = 0;
    for (size_t f = 0; f < ps->nf; f++) {
        if ((face_keep[f] != 0) != (want_kept != 0)) continue;
        nfaces++;
        for (int e = 0; e < 3; e++) {
            int32_t v = ps->faces[f * 3 + (size_t)e];
            if (v < 0 || (size_t)v >= ps->nv) {
                free(remap);
                return -1;
            }
            remap[v] = 0;
        }
    }
    size_t nvertices = 0;
    for (size_t i = 0; i < ps->nv; i++)
        if (remap[i] == 0) remap[i] = (int32_t)nvertices++;

    FILE *fp = fopen(path, "wb");
    if (fp == NULL) {
        free(remap);
        return -1;
    }
    fprintf(fp,
        "# %s: SeamOwn %s faces before PieceSet compaction\n"
        "# positions are %s; RGB is normalized world (x,y,z)\n"
        "# Both position and face selection are from the same post-weld, "
        "post-relax PieceSet passed to SeamOwn.\n"
        "# vertices=%zu faces=%zu\n"
        "g %s\n",
        group_name, want_kept ? "kept/main-spiral" : "dropped/overlap",
        flattened_uv ? "classifier-space flattening (u,v,0)"
                     : "original retained 3-D coordinates (z,y,x)",
        nvertices, nfaces, group_name);
    for (size_t i = 0; i < ps->nv; i++) {
        if (remap[i] < 0) continue;
        const float *p = &ps->verts[i * 3];
        double r = xyz_color_component((double)p[2], xyz_lo[2], xyz_hi[2]);
        double g = xyz_color_component((double)p[1], xyz_lo[1], xyz_hi[1]);
        double b = xyz_color_component((double)p[0], xyz_lo[0], xyz_hi[0]);
        if (flattened_uv)
            fprintf(fp, "v %.9g %.9g 0 %.6f %.6f %.6f\n",
                    (double)ps->uv[i * 2 + 0],
                    (double)ps->uv[i * 2 + 1], r, g, b);
        else
            fprintf(fp, "v %.9g %.9g %.9g %.6f %.6f %.6f\n",
                    (double)p[0], (double)p[1], (double)p[2], r, g, b);
    }
    for (size_t f = 0; f < ps->nf; f++) {
        if ((face_keep[f] != 0) != (want_kept != 0)) continue;
        int32_t a = remap[ps->faces[f * 3 + 0]];
        int32_t b = remap[ps->faces[f * 3 + 1]];
        int32_t c = remap[ps->faces[f * 3 + 2]];
        if (a < 0 || b < 0 || c < 0) {
            fclose(fp);
            free(remap);
            return -1;
        }
        fprintf(fp, "f %d %d %d\n", a + 1, b + 1, c + 1);
    }
    if (fclose(fp) != 0) {
        free(remap);
        return -1;
    }
    free(remap);
    if (out_vertices != NULL) *out_vertices = nvertices;
    if (out_faces != NULL) *out_faces = nfaces;
    return 0;
}

static int write_ownership_split_objs(const char *out_dir, const char *id,
                                      const PieceSet *ps,
                                      const SeamOwnResult *ownership)
{
    if (out_dir == NULL || id == NULL || ps == NULL || ownership == NULL ||
        ownership->face_keep == NULL)
        return -1;
    double lo[3] = {HUGE_VAL, HUGE_VAL, HUGE_VAL};
    double hi[3] = {-HUGE_VAL, -HUGE_VAL, -HUGE_VAL};
    double uv_lo[2] = {HUGE_VAL, HUGE_VAL};
    double uv_hi[2] = {-HUGE_VAL, -HUGE_VAL};
    for (size_t f = 0; f < ps->nf; f++) {
        for (int e = 0; e < 3; e++) {
            int32_t v = ps->faces[f * 3 + (size_t)e];
            if (v < 0 || (size_t)v >= ps->nv) return -1;
            const float *p = &ps->verts[(size_t)v * 3];
            for (int d = 0; d < 3; d++) {
                if ((double)p[d] < lo[d]) lo[d] = (double)p[d];
                if ((double)p[d] > hi[d]) hi[d] = (double)p[d];
            }
            for (int d = 0; d < 2; d++) {
                double q = (double)ps->uv[(size_t)v * 2 + (size_t)d];
                if (q < uv_lo[d]) uv_lo[d] = q;
                if (q > uv_hi[d]) uv_hi[d] = q;
            }
        }
    }
    char scroll_xyz_path[1024], leftover_xyz_path[1024];
    char scroll_uv_path[1024], leftover_uv_path[1024], stats_path[1024];
    snprintf(scroll_xyz_path, sizeof(scroll_xyz_path),
             "%s/%s_scroll_xyz.obj", out_dir, id);
    snprintf(leftover_xyz_path, sizeof(leftover_xyz_path),
             "%s/%s_leftover_xyz.obj", out_dir, id);
    snprintf(scroll_uv_path, sizeof(scroll_uv_path),
             "%s/%s_scroll_uv.obj", out_dir, id);
    snprintf(leftover_uv_path, sizeof(leftover_uv_path),
             "%s/%s_leftover_uv.obj", out_dir, id);
    size_t main_vertices = 0, main_faces = 0;
    size_t overlap_vertices = 0, overlap_faces = 0;
    if (write_xyz_face_category_obj(
            scroll_xyz_path, "scroll", ps, ownership->face_keep, 1, 0,
            lo, hi, &main_vertices, &main_faces) != 0 ||
        write_xyz_face_category_obj(
            leftover_xyz_path, "leftover", ps, ownership->face_keep, 0, 0,
            lo, hi, &overlap_vertices, &overlap_faces) != 0 ||
        write_xyz_face_category_obj(
            scroll_uv_path, "scroll", ps, ownership->face_keep, 1, 1,
            lo, hi, NULL, NULL) != 0 ||
        write_xyz_face_category_obj(
            leftover_uv_path, "leftover", ps, ownership->face_keep, 0, 1,
            lo, hi, NULL, NULL) != 0)
        return -1;

    size_t decision[7] = {0, 0, 0, 0, 0, 0, 0};
    if (ownership->face_dec != NULL)
        for (size_t f = 0; f < ps->nf; f++)
            if (ownership->face_dec[f] < 7)
                decision[ownership->face_dec[f]]++;
    snprintf(stats_path, sizeof(stats_path), "%s/%s_ownership_split.json",
             out_dir, id);
    FILE *fp = fopen(stats_path, "wb");
    if (fp != NULL) {
        fprintf(fp,
            "{\n"
            "  \"classification\": \"SeamOwn pre-compaction face_keep\",\n"
            "  \"classification_applied_to_source\": false,\n"
            "  \"xyz_positions\": \"original retained world coordinates\",\n"
            "  \"xyz_position_order\": \"z,y,x\",\n"
            "  \"uv_positions\": \"post-weld/post-relax classifier uv\",\n"
            "  \"uv_position_order\": \"u,v,0\",\n"
            "  \"color_source\": \"world position\",\n"
            "  \"color_order\": \"x,y,z\",\n"
            "  \"uv_bbox\": [[%.9g, %.9g], [%.9g, %.9g]],\n"
            "  \"world_bbox_zyx\": [[%.9g, %.9g, %.9g], "
            "[%.9g, %.9g, %.9g]],\n"
            "  \"main_spiral\": {\"vertices\": %zu, \"faces\": %zu},\n"
            "  \"overlap\": {\"vertices\": %zu, \"faces\": %zu},\n"
            "  \"face_decisions\": {\"none\": %zu, \"overlap_kept\": %zu, "
            "\"overlap_dropped\": %zu, \"seam_kept\": %zu, "
            "\"seam_dropped\": %zu, \"rehomed\": %zu}\n"
            "}\n",
            uv_lo[0], uv_lo[1], uv_hi[0], uv_hi[1],
            lo[0], lo[1], lo[2], hi[0], hi[1], hi[2],
            main_vertices, main_faces, overlap_vertices, overlap_faces,
            decision[0], decision[1], decision[3], decision[4], decision[5],
            decision[6]);
        fclose(fp);
    }
    logf_both("[step3] diagnostic OBJ split: scroll=%zu faces "
              "leftover=%zu faces; XYZ and UV pairs written "
              "(%s, %s, %s, %s)\n",
              main_faces, overlap_faces, scroll_xyz_path,
              leftover_xyz_path, scroll_uv_path, leftover_uv_path);
    return 0;
}

static void usage(void)
{
    fprintf(stderr,
        "usage: scroll_unroll <placed_dir> <out_dir> --raw <cubes_RAW_dir>\n"
        "       [--steps S] [--id X] [--du F] [--dv F] [--band-cols N]\n"
        "       [--raw-chunk N] [--range F] [--nsteps N] [--threads N]\n"
        "       [--snap-tensor-weight F] [--snap-tensor-radius F] (legacy)\n"
        "       [--snap-recto-iters N] [--snap-recto-range F]\n"
        "       [--stretch-ratio F] [--stretch-floor F] [--max-edge F]\n"
        "       [--synth-max-edge F]\n"
        "       [--strip-w N] [--dark-thresh N] [--seam-top N]\n"
        "       [--apply-ownership]  (unsafe legacy face deletion; default is\n"
        "                              diagnostic-only ownership)\n"
        "       [--no-diag] [--no-preview] [--no-xyzmap]\n"
        "       [--export-tifxyz DIR] [--tifxyz-du F] [--tifxyz-dv F]\n"
        "       [--tifxyz-flip-u] [--tifxyz-flip-v] [--tifxyz-winding]\n"
        "       [--export-atlas DIR] [--atlas-wraps N] [--atlas-slab F]\n"
        "       [--z-range Z0 Z1]  (only cubes with z-origin in [Z0,Z1);\n"
        "                           slab/strip export from a shared solve)\n"
        "       scroll_unroll --selftest\n"
        "  --export-atlas DIR   L3 atlas: split the strip into (wrap x z-slab)\n"
        "              tifxyz segments under DIR + atlas.json (uuid=DIR/<seg>);\n"
        "              --atlas-wraps wraps per piece (1), --atlas-slab z-slab\n"
        "              height vox (4096). Works with --steps 0 (raw placement)\n"
        "  --no-xyzmap turn off the per-stage world-pos->RGB xyzmap PNGs\n"
        "  --max-edge F  optional 3D-edge cap for real/source faces (default\n"
        "              disabled; their valid scale follows the remesher)\n"
        "  --synth-max-edge F  synthetic-fill membrane cap (default 6 vox)\n"
        "  --steps S   stage digits, default 12345 (1 bake, 2 join =\n"
        "              weld+relax, 3 overlap, 4 snap, 5 final relax);\n"
        "              unbuilt stages are skipped; --steps 0 runs NO stages\n"
        "              (valid only with --export-tifxyz: export the placed\n"
        "              dir's registered UVs as-is, no bake)\n"
        "  --export-tifxyz DIR  after the last stage, write the strip as a\n"
        "              VC3D tifxyz segment (x/y/z.tif + mask.tif + meta.json)\n"
        "              into DIR (meta uuid = DIR basename); --tifxyz-du/-dv\n"
        "              set the segment grid step (default --du/--dv),\n"
        "              --tifxyz-flip-u/-v mirror the grid, --tifxyz-winding\n"
        "              adds winding.tif = phi/2pi\n");
}

/* one row of the end summary + one JSON block */
typedef struct StageRow {
    char name[32];
    RasterStats rs;
    StripMetrics sm;
    double secs;
    int have_sm;
} StageRow;

typedef struct RunCtx {
    const char *out_dir;
    const char *id;
    Arena_T arena;
    const int32_t *seam_cols;
    size_t n_seam_cols;
    int dark_thresh;
    int seam_top;
    int do_preview;
    double v0;           /* step1 canvas origin (v == world z) */
    long vchunk;         /* cube-grid pitch for v-seam rows (128; 0 = off) */
    StageRow rows[8];
    int n_rows;
} RunCtx;

/* previews + metrics + summary row for the stage just baked */
static void stage_report(RunCtx *rc, const char *step, const RasterStats *rs,
                         double secs)
{
    StageRow *row = &rc->rows[rc->n_rows];
    memset(row, 0, sizeof(*row));
    snprintf(row->name, sizeof(row->name), "%s", step);
    row->rs = *rs;
    row->secs = secs;

    char tif[1024], cls[1024], png[1024];
    snprintf(tif, sizeof(tif), "%s/%s_%s_rawtex.tif", rc->out_dir, rc->id, step);
    snprintf(cls, sizeof(cls), "%s/%s_%s_diagclass.tif", rc->out_dir, rc->id, step);

    double *col_dcol = NULL;
    if (StripMetrics_compute(rc->arena, tif, cls, rc->seam_cols,
                             rc->n_seam_cols, rc->dark_thresh, &col_dcol,
                             &row->sm) == 0) {
        row->have_sm = 1;
        logf_both("[%s] metrics: fill=%.1f%% multi=%.2f%% dark=%.2f%% "
                  "seam_dcol=%.2f base_dcol=%.2f seam_ratio=%.2f "
                  "(%zu seam cols)\n", step, 100.0 * row->sm.fill,
                  100.0 * row->sm.multi_frac, 100.0 * row->sm.dark_frac,
                  row->sm.seam_dcol_mean, row->sm.base_dcol_mean,
                  row->sm.seam_ratio, row->sm.n_seam_cols);
        if (rc->vchunk > 0
            && StripMetrics_vseams(rc->arena, tif, cls, rc->v0, rc->vchunk,
                                   &row->sm) == 0
            && row->sm.n_vseam_rows > 0)
            logf_both("[%s] vseams: %zu rows gap_fill=%.2f ratio=%.2f\n",
                      step, row->sm.n_vseam_rows, row->sm.vseam_gap_fill,
                      row->sm.vseam_ratio);
    } else {
        logf_both("[%s] WARN: metrics pass failed (missing diagclass?)\n",
                  step);
    }

    if (rc->do_preview) {
        StripPreviewOpts po;
        StripPreviewOpts_default(&po);
        po.verbose = 0;
        po.overlay = cls;
        snprintf(png, sizeof(png), "%s/%s_%s_rawtex_preview.png",
                 rc->out_dir, rc->id, step);
        if (StripPreview_overview(rc->arena, tif, png, &po) != 0)
            logf_both("[%s] WARN: rawtex preview failed\n", step);
        po.overlay = NULL;
        po.mode = "class";
        snprintf(png, sizeof(png), "%s/%s_%s_diagclass_preview.png",
                 rc->out_dir, rc->id, step);
        if (StripPreview_overview(rc->arena, cls, png, &po) != 0)
            logf_both("[%s] WARN: diagclass preview failed\n", step);

        /* seam zoom: top-K worst seam columns by per-col |dcol| */
        if (col_dcol != NULL && rc->n_seam_cols > 0 && rc->seam_top > 0) {
            size_t k = (size_t)rc->seam_top;
            if (k > rc->n_seam_cols) k = rc->n_seam_cols;
            long *top = (long *)ARENA_ALLOC(rc->arena, k * sizeof(long));
            uint8_t *used = (uint8_t *)ARENA_CALLOC(rc->arena,
                                                    rc->n_seam_cols, 1);
            size_t got = 0;
            for (size_t r = 0; r < k; r++) {
                double best = -1.0;
                size_t bi = (size_t)-1;
                for (size_t i = 0; i < rc->n_seam_cols; i++) {
                    if (used[i] || col_dcol[i] < 0.0) continue;
                    if (col_dcol[i] > best) { best = col_dcol[i]; bi = i; }
                }
                if (bi == (size_t)-1) break;
                used[bi] = 1;
                top[got++] = (long)rc->seam_cols[bi];
            }
            if (got > 0) {
                StripPreviewOpts zo;
                StripPreviewOpts_default(&zo);
                zo.verbose = 0;
                zo.overlay = cls;
                snprintf(png, sizeof(png), "%s/%s_%s_seamzoom.png",
                         rc->out_dir, rc->id, step);
                if (StripPreview_windows(rc->arena, tif, png, top, got, 256,
                                         &zo) != 0)
                    logf_both("[%s] WARN: seam zoom failed\n", step);
            }
        }
    }
    rc->n_rows++;
}

static void write_pipeline_json(const RunCtx *rc, const char *placed_dir,
                                const char *raw_dir, const PieceSet *ps,
                                double total_secs)
{
    char jp[1024];
    snprintf(jp, sizeof(jp), "%s/%s_pipeline_stats.json", rc->out_dir, rc->id);
    FILE *jf = fopen(jp, "w");
    if (jf == NULL)
        return;
    fputs("{\n"
          "  \"tool\": \"scroll_unroll\",\n"
          "  \"placed_dir\": ", jf);
    write_json_string(jf, placed_dir);
    fputs(", \"raw_dir\": ", jf);
    write_json_string(jf, raw_dir);
    fprintf(jf, ",\n"
            "  \"n_cubes\": %zu, \"nv\": %zu, \"nf\": %zu,\n"
            "  \"n_seam_cols\": %zu,\n"
            "  \"stages\": [\n",
            ps->n_cubes, ps->nv, ps->nf,
            rc->n_seam_cols);
    for (int i = 0; i < rc->n_rows; i++) {
        const StageRow *r = &rc->rows[i];
        fprintf(jf, "    {\n"
                "      \"name\": \"%s\",\n"
                "      \"W\": %zu, \"H\": %zu, \"bands\": %zu,\n"
                "      \"window\": [%.1f, %.1f],\n"
                "      \"fill\": %.4f, \"filled_px\": %zu, \"multi_px\": %zu,\n"
                "      \"skip_uv_faces\": %zu, \"skip_3d_faces\": %zu, "
                "\"skip_own_faces\": %zu,\n"
                "      \"multi_frac\": %.4f, \"dark_frac\": %.4f,\n"
                "      \"seam_dcol_mean\": %.3f, \"base_dcol_mean\": %.3f, "
                "\"seam_ratio\": %.3f,\n"
                "      \"vseam_gap_fill\": %.3f, \"vseam_ratio\": %.3f, "
                "\"n_vseam_rows\": %zu,\n"
                "      \"synth_px\": %zu, \"synth_frac\": %.4f, "
                "\"fill_total\": %.4f,\n"
                "      \"seconds\": %.2f\n"
                "    }%s\n",
                r->name, r->rs.W, r->rs.H, r->rs.n_bands, r->rs.lo, r->rs.hi,
                r->rs.fill, r->rs.filled, r->rs.multi, r->rs.skip_uv,
                r->rs.skip_3d, r->rs.skip_own,
                r->have_sm ? r->sm.multi_frac : 0.0,
                r->have_sm ? r->sm.dark_frac : 0.0,
                r->have_sm ? r->sm.seam_dcol_mean : 0.0,
                r->have_sm ? r->sm.base_dcol_mean : 0.0,
                r->have_sm ? r->sm.seam_ratio : 0.0,
                r->have_sm ? r->sm.vseam_gap_fill : 0.0,
                r->have_sm ? r->sm.vseam_ratio : 0.0,
                r->have_sm ? r->sm.n_vseam_rows : (size_t)0,
                r->have_sm ? r->sm.synth_px : (size_t)0,
                r->have_sm ? r->sm.synth_frac : 0.0,
                r->have_sm ? r->sm.fill_total : 0.0,
                r->secs, i + 1 < rc->n_rows ? "," : "");
    }
    fprintf(jf, "  ],\n  \"total_seconds\": %.2f\n}\n", total_secs);
    fclose(jf);
}

int main(int argc, char **argv)
{
    if (argc >= 2 && strcmp(argv[1], "--selftest") == 0) {
        int f = write_json_string_selftest();
        f += PieceSet_selftest();
        f += ScrollRaster_selftest();
        f += StripMetrics_selftest();
        f += StripPreview_selftest();
        f += SeamOwn_selftest();
        f += SnapGrid_selftest();
        f += UVWeld_selftest();
        f += RelaxGrid_selftest();
        f += TifxyzExport_selftest();
        f += Scaffold_selftest();
        f += ExportAtlas_selftest();
        fprintf(stderr, "scroll_unroll --selftest: %s (%d failures)\n",
                f == 0 ? "PASS" : "FAIL", f);
        return f == 0 ? 0 : 1;
    }
    if (argc < 4) {
        usage();
        return 2;
    }
    const char *placed_dir = argv[1];
    const char *out_dir = argv[2];
    const char *raw_dir = NULL;
    char id[256] = "";
    char steps[16] = "12345";
    long raw_chunk = 128;
    int threads = 32;
    int dark_thresh = 40;
    int seam_top = 8;
    int do_preview = 1;
    double snap_tensor_weight = 0.0, snap_tensor_radius = 2.0;
    int snap_recto_iters = 4;
    double snap_recto_range = 3.0;
    double own_radius_gate = 0.0, own_cell = 0.0;
    long own_region_cap = 0;
    int own_grid_cap = 0, own_seam_cut = 1, own_halfband = 0;
    int own_rehome = 1;
    int own_apply = 0;
    const char *tifxyz_dir = NULL;
    double tifxyz_du = 0.0, tifxyz_dv = 0.0;   /* 0 = inherit --du/--dv */
    int tifxyz_flip_u = 0, tifxyz_flip_v = 0, tifxyz_winding = 0;
    const char *atlas_dir = NULL;              /* L3 export atlas root */
    int atlas_wraps = 1;
    double atlas_slab = 4096.0;
    int atlas_relax = 0;
    long z_lo = LONG_MIN, z_hi = LONG_MAX;     /* z-slab / strip filter (vox) */
    RasterOpts ro;
    RasterOpts_default(&ro);

    for (int i = 3; i < argc; i++) {
        if (strcmp(argv[i], "--raw") == 0 && i + 1 < argc)
            raw_dir = argv[++i];
        else if (strcmp(argv[i], "--id") == 0 && i + 1 < argc)
            snprintf(id, sizeof(id), "%s", argv[++i]);
        else if (strcmp(argv[i], "--steps") == 0 && i + 1 < argc)
            snprintf(steps, sizeof(steps), "%s", argv[++i]);
        else if (strcmp(argv[i], "--du") == 0 && i + 1 < argc)
            ro.du = atof(argv[++i]);
        else if (strcmp(argv[i], "--dv") == 0 && i + 1 < argc)
            ro.dv = atof(argv[++i]);
        else if (strcmp(argv[i], "--band-cols") == 0 && i + 1 < argc)
            ro.band_cols = (size_t)atoll(argv[++i]);
        else if (strcmp(argv[i], "--raw-chunk") == 0 && i + 1 < argc)
            raw_chunk = atol(argv[++i]);
        else if (strcmp(argv[i], "--range") == 0 && i + 1 < argc)
            ro.range = atof(argv[++i]);
        else if (strcmp(argv[i], "--nsteps") == 0 && i + 1 < argc)
            ro.nsteps = atoi(argv[++i]);
        else if (strcmp(argv[i], "--threads") == 0 && i + 1 < argc)
            threads = atoi(argv[++i]);
        else if (strcmp(argv[i], "--snap-tensor-weight") == 0 && i + 1 < argc)
            snap_tensor_weight = atof(argv[++i]);
        else if (strcmp(argv[i], "--snap-tensor-radius") == 0 && i + 1 < argc)
            snap_tensor_radius = atof(argv[++i]);
        else if (strcmp(argv[i], "--snap-recto-iters") == 0 && i + 1 < argc)
            snap_recto_iters = atoi(argv[++i]);
        else if (strcmp(argv[i], "--snap-recto-range") == 0 && i + 1 < argc)
            snap_recto_range = atof(argv[++i]);
        else if (strcmp(argv[i], "--stretch-ratio") == 0 && i + 1 < argc)
            ro.stretch_ratio = atof(argv[++i]);
        else if (strcmp(argv[i], "--stretch-floor") == 0 && i + 1 < argc)
            ro.stretch_floor = atof(argv[++i]);
        else if (strcmp(argv[i], "--max-edge") == 0 && i + 1 < argc)
            ro.max_edge3d = atof(argv[++i]);
        else if (strcmp(argv[i], "--synth-max-edge") == 0 && i + 1 < argc)
            ro.synth_max_edge3d = atof(argv[++i]);
        else if (strcmp(argv[i], "--strip-w") == 0 && i + 1 < argc)
            ro.strip_w = atoi(argv[++i]);
        else if (strcmp(argv[i], "--dark-thresh") == 0 && i + 1 < argc)
            dark_thresh = atoi(argv[++i]);
        else if (strcmp(argv[i], "--seam-top") == 0 && i + 1 < argc)
            seam_top = atoi(argv[++i]);
        else if (strcmp(argv[i], "--no-diag") == 0)
            ro.write_diag = 0;
        else if (strcmp(argv[i], "--no-preview") == 0)
            do_preview = 0;
        else if (strcmp(argv[i], "--radius-gate") == 0 && i + 1 < argc)
            own_radius_gate = atof(argv[++i]);
        else if (strcmp(argv[i], "--own-cell") == 0 && i + 1 < argc)
            own_cell = atof(argv[++i]);
        else if (strcmp(argv[i], "--region-cap") == 0 && i + 1 < argc)
            own_region_cap = atol(argv[++i]);
        else if (strcmp(argv[i], "--grid-cap") == 0 && i + 1 < argc)
            own_grid_cap = atoi(argv[++i]);
        else if (strcmp(argv[i], "--seam-halfband") == 0 && i + 1 < argc)
            own_halfband = atoi(argv[++i]);
        else if (strcmp(argv[i], "--no-seam-cut") == 0)
            own_seam_cut = 0;
        else if (strcmp(argv[i], "--no-rehome") == 0)
            own_rehome = 0;
        else if (strcmp(argv[i], "--apply-ownership") == 0)
            own_apply = 1;
        else if (strcmp(argv[i], "--no-xyzmap") == 0)
            ro.write_xyzmap = 0;
        else if (strcmp(argv[i], "--export-tifxyz") == 0 && i + 1 < argc)
            tifxyz_dir = argv[++i];
        else if (strcmp(argv[i], "--tifxyz-du") == 0 && i + 1 < argc)
            tifxyz_du = atof(argv[++i]);
        else if (strcmp(argv[i], "--tifxyz-dv") == 0 && i + 1 < argc)
            tifxyz_dv = atof(argv[++i]);
        else if (strcmp(argv[i], "--tifxyz-flip-u") == 0)
            tifxyz_flip_u = 1;
        else if (strcmp(argv[i], "--tifxyz-flip-v") == 0)
            tifxyz_flip_v = 1;
        else if (strcmp(argv[i], "--tifxyz-winding") == 0)
            tifxyz_winding = 1;
        else if (strcmp(argv[i], "--export-atlas") == 0 && i + 1 < argc)
            atlas_dir = argv[++i];
        else if (strcmp(argv[i], "--atlas-wraps") == 0 && i + 1 < argc)
            atlas_wraps = atoi(argv[++i]);
        else if (strcmp(argv[i], "--atlas-slab") == 0 && i + 1 < argc)
            atlas_slab = atof(argv[++i]);
        else if (strcmp(argv[i], "--atlas-relax") == 0)
            atlas_relax = 1;
        else if (strcmp(argv[i], "--z-range") == 0 && i + 2 < argc) {
            z_lo = strtol(argv[++i], NULL, 10);
            z_hi = strtol(argv[++i], NULL, 10);
        }
        else
            fprintf(stderr, "WARN: ignoring unknown arg %s\n", argv[i]);
    }
    if (raw_dir == NULL) {
        usage();
        return 2;
    }
    if (id[0] == '\0') {
        const char *base = ves_path_basename(placed_dir);
        snprintf(id, sizeof(id), "%s", base[0] ? base : "whole");
    }
    (void)threads;   /* consumed by stages 2/4/5 (relax + snap) */

    int want[6] = { 0, 0, 0, 0, 0, 0 };
    for (const char *s = steps; *s; s++)
        if (*s >= '1' && *s <= '5') want[*s - '0'] = 1;
    int have[6] = { 0, 1, HAVE_STAGE2, HAVE_STAGE3, HAVE_STAGE4, HAVE_STAGE5 };
    for (int s = 2; s <= 5; s++)
        if (want[s] && !have[s]) {
            fprintf(stderr, "WARN: stage %d not built into this binary -- "
                    "skipped\n", s);
            want[s] = 0;
        }
    int export_only = 0;
    if (!want[1] && !want[2] && !want[3] && !want[4] && !want[5]) {
        if ((tifxyz_dir != NULL || atlas_dir != NULL) && strcmp(steps, "0") == 0) {
            export_only = 1;   /* --steps 0: export the placed UVs as-is */
        } else {
            fprintf(stderr, "ERROR: no runnable stages in --steps %s\n", steps);
            return 2;
        }
    }

    {
        char probe[1024];
        snprintf(probe, sizeof(probe), "%s/x", out_dir);
        ves_ensure_parent_dir(probe);
        char lp[1024];
        snprintf(lp, sizeof(lp), "%s/scroll_unroll.log", out_dir);
        g_log = fopen(lp, "w");
    }
    double t_all = ves_clock_sec();
    logf_both("scroll_unroll: placed=%s out=%s raw=%s id=%s steps=%d%d%d%d%d "
              "du=%.2f dv=%.2f band=%zu\n", placed_dir, out_dir, raw_dir, id,
              want[1], want[2] ? 2 : 0, want[3] ? 3 : 0, want[4] ? 4 : 0,
              want[5] ? 5 : 0, ro.du, ro.dv, ro.band_cols);

    Arena_T arena = Arena_new();
    PieceSet ps;
    double t0 = ves_clock_sec();
    if (PieceSet_build_z(arena, placed_dir, z_lo, z_hi, &ps) != 0) {
        logf_both("ERROR: no complete placed cubes under %s%s\n", placed_dir,
                  (z_lo != LONG_MIN || z_hi != LONG_MAX) ? " (in z-range)" : "");
        return 1;
    }
    if (z_lo != LONG_MIN || z_hi != LONG_MAX)
        logf_both("[z-range] slab z in [%ld,%ld): %zu cubes\n",
                  z_lo, z_hi, ps.n_cubes);
    logf_both("[pieces] %zu cubes -> %zu verts, %zu kept faces; "
              "u=[%.0f,%.0f] (%.0f vox) v=[%.0f,%.0f] (%.2fs)\n",
              ps.n_cubes, ps.nv, ps.nf, ps.u_min, ps.u_max,
              ps.u_max - ps.u_min, ps.v_min, ps.v_max, ves_clock_sec() - t0);

    /* ONE RAW table for every stage. Lives in the main arena; never spanned
     * by an Arena_save/restore (scroll_raster.c aliasing rule). Pre-warmed
     * when a parallel-sampling stage will run. */
    CubeTable ct;
    if (cubetable_init(&ct, arena, raw_dir, raw_chunk, ps.verts, ps.nv,
                       ro.range + 2.0) != 0) {
        logf_both("ERROR: cubetable_init failed (%s)\n", raw_dir);
        return 1;
    }
    if (want[2] || want[3] || want[4] || want[5]) {
        t0 = ves_clock_sec();
        int n = cubetable_prewarm_all(&ct);
        logf_both("[raw] pre-warmed %d cubes (%d missing) (%.1fs)\n",
                  n, ct.n_missing, ves_clock_sec() - t0);
    }
    ro.ct = &ct;

    RunCtx rc;
    memset(&rc, 0, sizeof(rc));
    rc.out_dir = out_dir;
    rc.id = id;
    rc.arena = arena;
    rc.dark_thresh = dark_thresh;
    rc.seam_top = seam_top;
    rc.do_preview = do_preview;

    char prefix[1024];
    RasterStats st1;
    memset(&st1, 0, sizeof(st1));
    double t_raster = 0.0;

    /* ---- step1: base bake (skipped by --steps 0 export-only) -------------- */
    if (export_only) {
        logf_both("[step1] skipped (--steps 0 export-only)\n");
    } else {
        snprintf(prefix, sizeof(prefix), "%s/%s_step1_ribbon", out_dir, id);
        t0 = ves_clock_sec();
        if (ScrollRaster_run(arena, &ps, raw_dir, raw_chunk, prefix, &ro,
                             &st1) != 0) {
            logf_both("ERROR: step1 raster failed\n");
            return 1;
        }
        t_raster = ves_clock_sec() - t0;
        logf_both("[step1] %zux%zu px (origin u=%.0f v=%.0f), %zu band(s): "
                  "%.1f%% filled, multi=%zu (%.2f%% of filled), "
                  "skip_uv=%zu skip_3d=%zu faces, window=[%.0f,%.0f], "
                  "raw cubes %d loaded / %d missing (%.1fs)\n",
                  st1.W, st1.H, st1.u0, st1.v0, st1.n_bands, 100.0 * st1.fill,
                  st1.multi,
                  st1.filled ? 100.0 * (double)st1.multi / (double)st1.filled
                             : 0.0,
                  st1.skip_uv, st1.skip_3d, st1.lo, st1.hi,
                  st1.raw_loaded, st1.raw_missing, t_raster);

        /* later stages re-bake into the SAME window so stage TIFs difference */
        ro.lo_fixed = st1.lo;
        ro.hi_fixed = st1.hi;

        /* seam columns (needs the canvas mapping from the bake) */
        {
            int32_t *cols = NULL;
            size_t n = 0;
            StripMetrics_seam_cols(arena, &ps, st1.u0, ro.du, st1.W, &cols,
                                   &n);
            rc.seam_cols = cols;
            rc.n_seam_cols = n;
            logf_both("[seams] %zu seam columns from %zu cubes\n", n,
                      ps.n_cubes);
            char sp[1024];
            snprintf(sp, sizeof(sp), "%s/%s_seamcols.txt", out_dir, id);
            FILE *sf = fopen(sp, "w");
            if (sf != NULL) {
                for (size_t i = 0; i < n; i++)
                    fprintf(sf, "%d\n", cols[i]);
                fclose(sf);
            }
        }
        /* v-seam rows: world-z cube boundaries mapped into the canvas (v == z;
         * the prediction grid pitch is 128, independent of --raw-chunk) */
        rc.v0 = st1.v0;
        rc.vchunk = 128;
        stage_report(&rc, "step1_ribbon", &st1, t_raster);
    }

    /* v1-compat stats JSON for the step1 bake */
    if (!export_only) {
        char jp[1024];
        snprintf(jp, sizeof(jp), "%s/%s_strip_stats.json", out_dir, id);
        FILE *jf = fopen(jp, "w");
        if (jf != NULL) {
            fputs("{\n"
                  "  \"tool\": \"scroll_unroll\",\n"
                  "  \"placed_dir\": ", jf);
            write_json_string(jf, placed_dir);
            fputs(", \"raw_dir\": ", jf);
            write_json_string(jf, raw_dir);
            fprintf(jf, ",\n"
                    "  \"n_cubes\": %zu, \"nv\": %zu, \"nf\": %zu,\n"
                    "  \"u_range\": [%.2f, %.2f], \"v_range\": [%.2f, %.2f],\n"
                    "  \"du\": %.3f, \"dv\": %.3f,\n"
                    "  \"W\": %zu, \"H\": %zu, \"origin\": [%.2f, %.2f],\n"
                    "  \"bands\": %zu, \"band_cols\": %zu,\n"
                    "  \"fill\": %.4f, \"filled_px\": %zu, \"multi_px\": %zu,\n"
                    "  \"skip_uv_faces\": %zu, \"skip_3d_faces\": %zu,\n"
                    "  \"window\": [%.1f, %.1f],\n"
                    "  \"raw_loaded\": %d, \"raw_missing\": %d,\n"
                    "  \"seconds\": %.2f\n}\n",
                    ps.n_cubes, ps.nv, ps.nf,
                    ps.u_min, ps.u_max, ps.v_min, ps.v_max, ro.du, ro.dv,
                    st1.W, st1.H, st1.u0, st1.v0, st1.n_bands,
                    ro.band_cols, st1.fill, st1.filled, st1.multi,
                    st1.skip_uv, st1.skip_3d, st1.lo, st1.hi,
                    st1.raw_loaded, st1.raw_missing, t_raster);
            fclose(jf);
        }
    }

    /* ---- step2: JOIN -- weld z- and u-seams, then flatten ------------------ */
    /* The dv gap the old chain could not close is gone at the source
     * (trim-inset-0 dump), so the standard weld gates stitch BOTH seam
     * families; the banded relax then diffuses the weld strain ("merging at
     * the geometry level") before any ownership/snap work sees the mesh. */
    int ran_join = 0;
    if (want[2]) {
        double t2 = ves_clock_sec();
        UVWeldOpts wo;
        UVWeldOpts_default(&wo);
        wo.chunk = 128;
        wo.verbose = 1;
        UVWeldStats wst;
        if (UVWeld_run(arena, &ps, &wo, &wst) != 0) {
            logf_both("[step2] WARN: weld failed -- relaxing unwelded\n");
            memset(&wst, 0, sizeof(wst));
        } else {
            ran_join = 1;   /* cross-cube edges exist from here on */
            logf_both("[step2] weld: shell=%zu pairs=%zu clusters=%zu "
                      "merged-refs=%zu degen=%zu rej(cube/uv)=%zu/%zu "
                      "|du|med=%.2f\n", wst.n_boundary, wst.n_pairs,
                      wst.n_clusters, wst.n_merged_verts, wst.n_degen_faces,
                      wst.rej_cube, wst.rej_uv, wst.du_absmed);
            PieceSet_refresh_normals(&ps);   /* weld moved verts to means */
        }

        char fib_tif[1024];
        snprintf(fib_tif, sizeof(fib_tif), "%s/%s_step1_ribbon_rawtex.tif",
                 out_dir, id);
        RelaxGridOpts rgo;
        RelaxGridOpts_default(&rgo);
        rgo.threads = threads < 16 ? threads : 16;
        rgo.verbose = 1;
        RelaxGridStats rst;
        if (RelaxGrid_run(arena, &ps, fib_tif, st1.u0, st1.v0, ro.du, ro.dv,
                          &rgo, &rst) != 0) {
            logf_both("[step2] WARN: relax failed -- uv unchanged\n");
            memset(&rst, 0, sizeof(rst));
        } else {
            logf_both("[step2] relax: bands=%zu run=%zu reverted=%zu "
                      "no_fiber=%zu moved=%zu stretch %.3f->%.3f "
                      "disp mean=%.2f max=%.2f (%.1fs)\n",
                      rst.n_bands, rst.bands_run, rst.bands_reverted,
                      rst.bands_no_fiber, rst.n_moved, rst.stretch_before,
                      rst.stretch_after, rst.mean_disp, rst.max_disp,
                      rst.seconds);
        }

        snprintf(prefix, sizeof(prefix), "%s/%s_step2_join", out_dir, id);
        RasterStats st2;
        t0 = ves_clock_sec();
        if (ScrollRaster_run(arena, &ps, raw_dir, raw_chunk, prefix, &ro,
                             &st2) != 0) {
            logf_both("ERROR: step2 raster failed\n");
            return 1;
        }
        logf_both("[step2] bake: %.1f%% filled, multi=%zu (%.1fs)\n",
                  100.0 * st2.fill, st2.multi, ves_clock_sec() - t0);
        stage_report(&rc, "step2_join", &st2, ves_clock_sec() - t2);

        char jp[1024];
        snprintf(jp, sizeof(jp), "%s/%s_step2_join_stats.json", out_dir, id);
        FILE *jf = fopen(jp, "w");
        if (jf != NULL) {
            fprintf(jf, "{\n"
                "  \"tool\": \"scroll_unroll step2 uv_weld+relax_grid\",\n"
                "  \"weld\": { \"shell\": %zu, \"pairs\": %zu, "
                "\"clusters\": %zu, \"merged_refs\": %zu, \"degen\": %zu,\n"
                "    \"rej_cube\": %zu, \"rej_uv\": %zu, "
                "\"du_absmed\": %.3f },\n"
                "  \"relax\": { \"bands\": %zu, \"run\": %zu, "
                "\"reverted\": %zu, \"no_fiber\": %zu, \"moved\": %zu,\n"
                "    \"stretch_before\": %.4f, \"stretch_after\": %.4f,\n"
                "    \"mean_disp\": %.3f, \"max_disp\": %.3f },\n"
                "  \"seconds\": %.2f\n}\n",
                wst.n_boundary, wst.n_pairs, wst.n_clusters,
                wst.n_merged_verts, wst.n_degen_faces, wst.rej_cube,
                wst.rej_uv, wst.du_absmed, rst.n_bands, rst.bands_run,
                rst.bands_reverted, rst.bands_no_fiber, rst.n_moved,
                rst.stretch_before, rst.stretch_after, rst.mean_disp,
                rst.max_disp, ves_clock_sec() - t2);
            fclose(jf);
        }
    }

    /* ---- step3: SeamOwn ownership + re-bake ------------------------------- */
    if (want[3]) {
        double t3 = ves_clock_sec();
        SeamOwnOpts so;
        SeamOwnOpts_default(&so);
        SeamOwn_read_index(placed_dir, &so);
        so.raw_dir = raw_dir;
        so.raw_chunk = raw_chunk;
        so.ct = &ct;
        so.verbose = 1;
        if (own_radius_gate > 0.0) so.radius_gate = own_radius_gate;
        if (own_cell > 0.0) so.cell = own_cell;
        if (own_region_cap > 0) so.region_cap = (size_t)own_region_cap;
        if (own_grid_cap > 0) so.grid_cap = own_grid_cap;
        if (own_halfband > 0) so.seam_halfband = own_halfband;
        so.seam_cut = own_seam_cut;
        /* Diagnostic ownership must be read-only. Rehome changes ps->uv, so
         * it is permitted only with the explicit legacy apply switch. */
        so.rehome = own_apply ? own_rehome : 0;

        SeamOwnResult sr;
        if (SeamOwn_run(arena, &ps, &so, &sr) != 0) {
            logf_both("[step3] WARN: SeamOwn failed -- continuing on the "
                      "unfiltered faces\n");
        } else {
            logf_both("[step3] seam_own: pairs true=%zu seam=%zu mystery=%zu; "
                      "%zu regions (%zu layers, max %zu f) picks "
                      "vote/energy/largest=%zu/%zu/%zu multicut=%zu "
                      "skip cap/grid=%zu/%zu; dropped=%zu of %zu; seams=%zu "
                      "(abut %zu, skip %zu) drop=%zu restore=%zu; "
                      "cells %zu->%zu (%.1f%%) E=%.2f (%.1fs)\n",
                      sr.n_true_pairs, sr.n_seam_pairs, sr.n_mystery_pairs,
                      sr.n_regions, sr.n_layers_total, sr.max_region_faces,
                      sr.n_vote_picks, sr.n_energy_picks, sr.n_largest_picks,
                      sr.n_multicut_fallbacks, sr.n_regions_skipped_cap,
                      sr.n_regions_skipped_grid, sr.n_dropped + sr.n_seam_dropped,
                      ps.nf, sr.n_seam_regions, sr.n_seams_abut_only,
                      sr.n_seams_skipped, sr.n_seam_dropped, sr.n_seam_restored,
                      sr.multi_cells_before, sr.multi_cells_after,
                      sr.multi_cells_before
                          ? 100.0 * (double)sr.multi_cells_after
                                / (double)sr.multi_cells_before : 0.0,
                      sr.energy_mean, ves_clock_sec() - t3);
            if (write_ownership_split_objs(out_dir, id, &ps, &sr) != 0)
                logf_both("[step3] WARN: XYZ ownership OBJ split failed\n");
            if (sr.n_rehome_layers + sr.n_rehome_dup + sr.n_rehome_blocked
                + sr.n_rehome_adjacent + sr.n_rehome_incoherent > 0)
                logf_both("[step3] rehome: %zu layers (%zu faces) found "
                          "homes; dup=%zu blocked=%zu adjacent=%zu "
                          "small=%zu incoherent=%zu\n",
                          sr.n_rehome_layers, sr.n_rehomed_faces,
                          sr.n_rehome_dup, sr.n_rehome_blocked,
                          sr.n_rehome_adjacent, sr.n_rehome_small,
                          sr.n_rehome_incoherent);

            RasterOpts ro2 = ro;
            ro2.face_keep = own_apply ? sr.face_keep : NULL;
            snprintf(prefix, sizeof(prefix), "%s/%s_step3_overlap", out_dir,
                     id);
            RasterStats st3;
            t0 = ves_clock_sec();
            if (ScrollRaster_run(arena, &ps, raw_dir, raw_chunk, prefix, &ro2,
                                 &st3) != 0) {
                logf_both("ERROR: step3 raster failed\n");
                return 1;
            }
            double t3r = ves_clock_sec() - t0;
            logf_both("[step3] bake%s: %.1f%% filled, multi=%zu (%.2f%% of "
                      "filled, was %.2f%%), skip_own=%zu faces; "
                      "diagnostic_would_drop=%zu (%.1fs)\n",
                      own_apply ? " (ownership applied)"
                                : " (diagnostic-only; source preserved)",
                      100.0 * st3.fill, st3.multi,
                      st3.filled ? 100.0 * (double)st3.multi
                                       / (double)st3.filled : 0.0,
                      st1.filled ? 100.0 * (double)st1.multi
                                       / (double)st1.filled : 0.0,
                      st3.skip_own, sr.n_dropped + sr.n_seam_dropped, t3r);
            stage_report(&rc, "step3_overlap", &st3,
                         ves_clock_sec() - t3);

            /* step3 stats JSON (SeamOwn counters + bake) */
            {
                char jp[1024];
                snprintf(jp, sizeof(jp), "%s/%s_step3_overlap_stats.json",
                         out_dir, id);
                FILE *jf = fopen(jp, "w");
                if (jf != NULL) {
                    fprintf(jf, "{\n"
                        "  \"tool\": \"scroll_unroll step3 seam_own\",\n"
                        "  \"ownership_applied\": %s,\n"
                        "  \"pairs\": { \"true\": %zu, \"seam\": %zu, "
                        "\"mystery\": %zu },\n"
                        "  \"regions\": %zu, \"layers\": %zu, "
                        "\"max_region_faces\": %zu,\n"
                        "  \"picks\": { \"vote\": %zu, \"energy\": %zu, "
                        "\"largest\": %zu, \"multicut\": %zu },\n"
                        "  \"skipped\": { \"cap\": %zu, \"grid\": %zu },\n"
                        "  \"dropped\": %zu, \"seam_regions\": %zu,\n"
                        "  \"seam\": { \"abut_only\": %zu, \"skipped\": %zu, "
                        "\"dropped\": %zu, \"restored\": %zu },\n"
                        "  \"multi_cells\": [%zu, %zu],\n"
                        "  \"rehome\": { \"layers\": %zu, \"faces\": %zu, "
                        "\"dup\": %zu, \"blocked\": %zu, \"adjacent\": %zu, "
                        "\"small\": %zu, \"incoherent\": %zu },\n"
                        "  \"energy_mean\": %.3f, \"have_texture\": %d,\n"
                        "  \"bake\": { \"filled_px\": %zu, \"multi_px\": %zu, "
                        "\"skip_own_faces\": %zu },\n"
                        "  \"step1_multi_px\": %zu,\n"
                        "  \"seconds\": %.2f\n}\n",
                        own_apply ? "true" : "false",
                        sr.n_true_pairs, sr.n_seam_pairs, sr.n_mystery_pairs,
                        sr.n_regions, sr.n_layers_total, sr.max_region_faces,
                        sr.n_vote_picks, sr.n_energy_picks,
                        sr.n_largest_picks, sr.n_multicut_fallbacks,
                        sr.n_regions_skipped_cap, sr.n_regions_skipped_grid,
                        sr.n_dropped, sr.n_seam_regions,
                        sr.n_seams_abut_only, sr.n_seams_skipped,
                        sr.n_seam_dropped, sr.n_seam_restored,
                        sr.multi_cells_before, sr.multi_cells_after,
                        sr.n_rehome_layers, sr.n_rehomed_faces,
                        sr.n_rehome_dup, sr.n_rehome_blocked,
                        sr.n_rehome_adjacent, sr.n_rehome_small,
                        sr.n_rehome_incoherent,
                        sr.energy_mean, sr.have_texture,
                        st3.filled, st3.multi, st3.skip_own, st1.multi,
                        ves_clock_sec() - t3);
                    fclose(jf);
                }
            }

            if (own_apply) {
                /* Legacy opt-in only: downstream stages see kept-only faces. */
                size_t nf_before = ps.nf;
                PieceSet_apply_facekeep(&ps, sr.face_keep);
                PieceSet_refresh_normals(&ps);
                logf_both("[step3] ownership explicitly applied: piece set "
                          "compacted %zu -> %zu faces\n", nf_before, ps.nf);
            } else {
                logf_both("[step3] ownership quarantined: preserved all %zu "
                          "faces and immutable UV for downstream stages\n",
                          ps.nf);
            }
        }
    }

    /* ---- step4: dark-vertex snap + re-bake -------------------------------- */
    if (want[4]) {
        double t4 = ves_clock_sec();
        SnapGridOpts sgo;
        SnapGridOpts_default(&sgo);
        sgo.tensor_weight = snap_tensor_weight;
        sgo.tensor_radius = snap_tensor_radius;
        sgo.recto_iters = snap_recto_iters;
        sgo.recto_range = snap_recto_range;
        {   /* axis from the placed index (same source as SeamOwn) */
            SeamOwnOpts axo;
            SeamOwnOpts_default(&axo);
            if (SeamOwn_read_index(placed_dir, &axo) == 0) {
                memcpy(sgo.axis_point, axo.axis_point, sizeof(sgo.axis_point));
                memcpy(sgo.axis_dir, axo.axis_dir, sizeof(sgo.axis_dir));
            }
        }
        sgo.threads = threads;
        /* after step2_join the mesh has cross-cube edges: per-cube slicing
         * is no longer exact -- run the whole mesh as one unit */
        sgo.global_mode = ran_join;
        sgo.verbose = 1;
        SnapGridStats sst;
        if (SnapGrid_run(arena, &ps, &ct, &sgo, &sst) != 0) {
            logf_both("[step4] WARN: snap failed -- continuing unmoved\n");
        } else {
            logf_both("[step4] two-pass snap%s: dark=%zu fixable=%zu crack=%zu "
                      "rejected=%zu; regions=%zu (fix=%zu anch=%zu) max=%zu; "
                      "repair moved=%zu reverted=%zu mean=%.2f max=%.2f "
                      "quilt=%.1f target_dist=%.2f qfb=%zu; "
                      "recto support=%zu moved=%zu reverted=%zu slope=%zu "
                      "mean=%.2f max=%.2f iter=%d (%.1fs)\n",
                      sgo.global_mode ? " (global)" : "",
                      sst.n_dark, sst.n_fixable, sst.n_crack, sst.n_rejected,
                      sst.nreg, sst.n_reg_fixable, sst.n_reg_anchorless,
                      sst.max_region_size, sst.n_moved, sst.n_reverted,
                      sst.mean_disp, sst.max_disp, sst.mean_quilt_cost,
                      sst.mean_target_dist, sst.n_quilt_fallback,
                      sst.n_recto_supported, sst.n_recto_moved,
                      sst.n_recto_reverted, sst.n_recto_slope_limited,
                      sst.recto_mean_disp, sst.recto_max_disp,
                      sst.recto_iterations,
                      ves_clock_sec() - t4);
            PieceSet_refresh_normals(&ps);

            snprintf(prefix, sizeof(prefix), "%s/%s_step4_snap", out_dir, id);
            RasterStats st4;
            t0 = ves_clock_sec();
            if (ScrollRaster_run(arena, &ps, raw_dir, raw_chunk, prefix, &ro,
                                 &st4) != 0) {
                logf_both("ERROR: step4 raster failed\n");
                return 1;
            }
            logf_both("[step4] bake: %.1f%% filled, multi=%zu (%.1fs)\n",
                      100.0 * st4.fill, st4.multi, ves_clock_sec() - t0);
            stage_report(&rc, "step4_snap", &st4, ves_clock_sec() - t4);

            char jp[1024];
            snprintf(jp, sizeof(jp), "%s/%s_step4_snap_stats.json", out_dir,
                     id);
            FILE *jf = fopen(jp, "w");
            if (jf != NULL) {
                fprintf(jf, "{\n"
                    "  \"tool\": \"scroll_unroll step4 snap_grid\",\n"
                    "  \"global_mode\": %d,\n"
                    "  \"sampled\": %zu, \"dark\": %zu, \"fixable\": %zu, "
                    "\"crack\": %zu, \"rejected\": %zu,\n"
                    "  \"regions\": { \"n\": %zu, \"fixable\": %zu, "
                    "\"crack\": %zu, \"anchorless\": %zu, \"mixed\": %zu, "
                    "\"reject\": %zu, \"max_size\": %zu },\n"
                    "  \"repair\": { \"moved\": %zu, \"reverted\": %zu, "
                    "\"mean_disp\": %.3f, \"max_disp\": %.3f, "
                    "\"quilt_cost_mean\": %.4f, \"target_dist_mean\": %.4f, "
                    "\"quilt_fallback\": %zu },\n"
                    "  \"recto\": { \"iterations\": %d, \"supported\": %zu, "
                    "\"moved\": %zu, \"reverted\": %zu, \"slope_limited\": %zu, "
                    "\"mean_disp\": %.4f, \"max_disp\": %.4f, "
                    "\"mean_gradient\": %.4f },\n"
                    "  \"window\": [%.1f, %.1f],\n"
                    "  \"units\": %zu, \"gco_fallback\": %zu,\n"
                    "  \"seconds\": %.2f\n}\n",
                    sgo.global_mode,
                    sst.n_sampled, sst.n_dark, sst.n_fixable, sst.n_crack,
                    sst.n_rejected, sst.nreg, sst.n_reg_fixable,
                    sst.n_reg_crack, sst.n_reg_anchorless, sst.n_reg_mixed,
                    sst.n_reg_reject, sst.max_region_size, sst.n_moved,
                    sst.n_reverted, sst.mean_disp, sst.max_disp,
                    sst.mean_quilt_cost, sst.mean_target_dist,
                    sst.n_quilt_fallback, sst.recto_iterations,
                    sst.n_recto_supported, sst.n_recto_moved,
                    sst.n_recto_reverted, sst.n_recto_slope_limited,
                    sst.recto_mean_disp, sst.recto_max_disp,
                    sst.recto_mean_gradient,
                    sst.win_lo, sst.win_hi, sst.n_cubes_run,
                    sst.n_cubes_gco_fallback, ves_clock_sec() - t4);
                fclose(jf);
            }
        }
    }

    /* ---- step5: final light banded relax + re-bake ------------------------- */
    if (want[5]) {
        double t5 = ves_clock_sec();
        /* fiber source = the newest baked rawtex */
        char fib_tif[1024];
        if (want[4])
            snprintf(fib_tif, sizeof(fib_tif), "%s/%s_step4_snap_rawtex.tif",
                     out_dir, id);
        else if (want[3])
            snprintf(fib_tif, sizeof(fib_tif),
                     "%s/%s_step3_overlap_rawtex.tif", out_dir, id);
        else if (want[2])
            snprintf(fib_tif, sizeof(fib_tif),
                     "%s/%s_step2_join_rawtex.tif", out_dir, id);
        else
            snprintf(fib_tif, sizeof(fib_tif),
                     "%s/%s_step1_ribbon_rawtex.tif", out_dir, id);

        RelaxGridOpts rgo;
        RelaxGridOpts_default(&rgo);
        rgo.base.sweeps = 8;   /* light polish: step2_join did the heavy relax */
        rgo.threads = threads < 16 ? threads : 16;
        rgo.verbose = 1;
        RelaxGridStats rst;
        if (RelaxGrid_run(arena, &ps, fib_tif, st1.u0, st1.v0, ro.du, ro.dv,
                          &rgo, &rst) != 0) {
            logf_both("[step5] WARN: relax failed -- uv unchanged\n");
            memset(&rst, 0, sizeof(rst));
        } else {
            logf_both("[step5] relax: bands=%zu run=%zu reverted=%zu "
                      "no_fiber=%zu moved=%zu stretch %.3f->%.3f "
                      "disp mean=%.2f max=%.2f (%.1fs)\n",
                      rst.n_bands, rst.bands_run, rst.bands_reverted,
                      rst.bands_no_fiber, rst.n_moved, rst.stretch_before,
                      rst.stretch_after, rst.mean_disp, rst.max_disp,
                      rst.seconds);
        }

        snprintf(prefix, sizeof(prefix), "%s/%s_step5_relax", out_dir, id);
        RasterStats st5;
        t0 = ves_clock_sec();
        if (ScrollRaster_run(arena, &ps, raw_dir, raw_chunk, prefix, &ro,
                             &st5) != 0) {
            logf_both("ERROR: step5 raster failed\n");
            return 1;
        }
        logf_both("[step5] bake: %.1f%% filled, multi=%zu (%.1fs)\n",
                  100.0 * st5.fill, st5.multi, ves_clock_sec() - t0);
        stage_report(&rc, "step5_relax", &st5, ves_clock_sec() - t5);

        char jp[1024];
        snprintf(jp, sizeof(jp), "%s/%s_step5_relax_stats.json", out_dir, id);
        FILE *jf = fopen(jp, "w");
        if (jf != NULL) {
            fprintf(jf, "{\n"
                "  \"tool\": \"scroll_unroll step5 relax_grid\",\n"
                "  \"relax\": { \"bands\": %zu, \"run\": %zu, "
                "\"reverted\": %zu, \"no_fiber\": %zu, \"moved\": %zu,\n"
                "    \"stretch_before\": %.4f, \"stretch_after\": %.4f,\n"
                "    \"mean_disp\": %.3f, \"max_disp\": %.3f },\n"
                "  \"seconds\": %.2f\n}\n",
                rst.n_bands, rst.bands_run,
                rst.bands_reverted, rst.bands_no_fiber, rst.n_moved,
                rst.stretch_before, rst.stretch_after, rst.mean_disp,
                rst.max_disp, ves_clock_sec() - t5);
            fclose(jf);
        }
    }

    /* ---- optional tifxyz export (VC3D segment handoff) --------------------- */
    int export_failed = 0;
    if (tifxyz_dir != NULL) {
        double te = ves_clock_sec();
        TifxyzOpts to;
        TifxyzOpts_default(&to);
        to.du = tifxyz_du > 0.0 ? tifxyz_du : ro.du;
        to.dv = tifxyz_dv > 0.0 ? tifxyz_dv : ro.dv;
        to.stretch_ratio = ro.stretch_ratio;   /* gate parity with the bakes */
        to.stretch_floor = ro.stretch_floor;
        to.max_edge3d = ro.max_edge3d;
        to.flip_u = tifxyz_flip_u;
        to.flip_v = tifxyz_flip_v;
        to.write_winding = tifxyz_winding;
        const char *uuid = ves_path_basename(tifxyz_dir);
        if (uuid[0] == '\0') uuid = id;
        TifxyzStats ts;
        if (TifxyzExport_run(arena, &ps, tifxyz_dir, uuid, &to, &ts) != 0) {
            logf_both("[tifxyz] ERROR: export to %s failed\n", tifxyz_dir);
            export_failed = 1;
        } else {
            logf_both("[tifxyz] wrote %s: %zux%zu px scale=[%.4f,%.4f] "
                      "valid=%zu (%.1f%%) multi=%zu conflicts=%zu "
                      "gated uv/3d/own=%zu/%zu/%zu area=%.0f vx2 "
                      "bbox x[%.0f,%.0f] y[%.0f,%.0f] z[%.0f,%.0f] (%.1fs)\n",
                      tifxyz_dir, ts.W, ts.H, 1.0 / to.du, 1.0 / to.dv,
                      ts.valid,
                      100.0 * (double)ts.valid / (double)(ts.W * ts.H),
                      ts.multi, ts.conflicts,
                      ts.skip_uv, ts.skip_3d, ts.skip_own, ts.area_vx2,
                      ts.bbox_lo[0], ts.bbox_hi[0], ts.bbox_lo[1],
                      ts.bbox_hi[1], ts.bbox_lo[2], ts.bbox_hi[2],
                      ves_clock_sec() - te);
        }
    }

    /* ---- optional export atlas (L3: per-(wrap x z-slab) tifxyz segments) --- */
    if (atlas_dir != NULL) {
        double ta = ves_clock_sec();
        ScaffoldCalib cal;
        if (Scaffold_read_calib(placed_dir, &cal) != 0)
            logf_both("[atlas] no placed_index.json calibration; using defaults "
                      "(sense=%d)\n", cal.sense);
        AtlasOpts ao; AtlasOpts_default(&ao);
        ao.piece_wraps = atlas_wraps > 0 ? atlas_wraps : 1;
        ao.slab_v = atlas_slab;
        ao.du = tifxyz_du > 0.0 ? tifxyz_du : ro.du;
        ao.dv = tifxyz_dv > 0.0 ? tifxyz_dv : ro.dv;
        ao.stretch_ratio = ro.stretch_ratio;
        ao.stretch_floor = ro.stretch_floor;
        ao.max_edge3d = ro.max_edge3d;
        ao.write_winding = 1;
        ao.relax = atlas_relax;
        ao.verbose = 0;
        AtlasStats as;
        if (ExportAtlas_run(arena, &ps, &cal, atlas_dir, id, &ao, &as) != 0) {
            logf_both("[atlas] ERROR: export to %s failed\n", atlas_dir);
            export_failed = 1;
        } else {
            logf_both("[atlas] %s: %zu pieces written (%zu empty, %zu over-cap, "
                      "%zu quarantine faces); valid=%zu conflict=%zu "
                      "(%.3f%%) in %.1fs\n",
                      atlas_dir, as.n_written, as.n_empty, as.n_over_cap,
                      as.n_quarantine_faces, as.total_valid_px,
                      as.total_conflict_px, 100.0 * as.conflict_frac,
                      ves_clock_sec() - ta);
        }
    }

    /* ---- end summary ------------------------------------------------------ */
    write_pipeline_json(&rc, placed_dir, raw_dir, &ps,
                        ves_clock_sec() - t_all);
    logf_both("---- summary ----------------------------------------------\n");
    logf_both("%-14s %10s %8s %8s %8s %7s %7s %7s %10s %8s\n",
              "stage", "fill%", "multi%", "dark%", "seamR", "vfill", "vseamR",
              "synth%", "drop_own", "secs");
    for (int i = 0; i < rc.n_rows; i++) {
        const StageRow *r = &rc.rows[i];
        logf_both("%-14s %9.1f%% %7.2f%% %7.2f%% %8.2f %7.2f %7.2f %6.2f%% "
                  "%10zu %8.1f\n",
                  r->name, 100.0 * (r->have_sm ? r->sm.fill : r->rs.fill),
                  100.0 * (r->have_sm ? r->sm.multi_frac : 0.0),
                  100.0 * (r->have_sm ? r->sm.dark_frac : 0.0),
                  r->have_sm ? r->sm.seam_ratio : 0.0,
                  r->have_sm ? r->sm.vseam_gap_fill : 0.0,
                  r->have_sm ? r->sm.vseam_ratio : 0.0,
                  100.0 * (r->have_sm ? r->sm.synth_frac : 0.0),
                  r->rs.skip_own, r->secs);
    }
    logf_both("scroll_unroll: %s (total %.1fs)\n",
              export_failed ? "DONE with tifxyz export FAILURE" : "OK",
              ves_clock_sec() - t_all);
    if (g_log) fclose(g_log);
    Arena_dispose(&arena);
    return export_failed ? 1 : 0;
}
