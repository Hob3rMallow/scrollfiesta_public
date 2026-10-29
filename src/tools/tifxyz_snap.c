/* ============================================================================
 * tifxyz_snap -- snap a tifxyz segment onto the papyrus RECTO edge of the RAW CT
 * (src/flatten/tifxyz_recto.c; 2026-10-09, PHerc0211 close-ups; 2026-10-10
 * whole-sheet blocks).
 *
 *   tifxyz_snap <in_seg_dir> <raw_root> <out_seg_dir> [--midline | --inset] [--smooth]
 *               [--no-prior | --prior <prior.json>] [--prior-only] [--window C0 C1]
 *   tifxyz_snap --selftest
 *
 * in_seg_dir  x.tif y.tif z.tif [mask.tif] meta.json: a VC3D tifxyz segment at
 *             1 node per voxel (vmesh_tifxyz, sheet_layers --tifxyz) -- a
 *             close-up or a whole-sheet block (4295 x 70k nodes).
 * raw_root    the uint8 Zarr v2 root (<root>/0/.zarray, uncompressed) or a
 *             directory of per-cube RAW TIFFs (raw_sample.h).
 * out_seg_dir created, refused if it already holds a meta.json: the moved grid,
 *             meta.json with the provenance.  Beside it:
 *             <out_seg_dir>_prior.json (the rays' reading of the wraps and the
 *             parameters derived from it), <out_seg_dir>_snap_report.json and
 *             the diagnostic PNGs <out_seg_dir>_recto_before.png (the recto
 *             offset of every node before the move), _disp.png (the move),
 *             _recto_after.png (the offset measured again afterwards), and the
 *             same through the median window, downsampled TXR_DIAG_SCALE x:
 *             blue = outward, red = inward, saturating at +-TXR_DIAG_FULL_VOX,
 *             dark grey = no edge found, black = invalid.  A run in one window
 *             (a close-up) also writes the rasters in full as float TIFFs.
 *
 * Big segments run in windows (TifxyzRecto_run_windows): cores of
 * TXR_WINDOW_COLS columns, each measured and solved with TXR_WINDOW_HALO
 * columns of the INPUT either side, sharing one prior (the segment's) and one
 * CT table loaded only where the nodes reach (cubetable_prewarm_points).
 *   --window C0 C1  only the columns [C0, C1): read by rows with the halo,
 *                   written alone (meta.json "window" says where they sit).
 *   --prior FILE    the parameters of that record (another segment's prior --
 *                   one per scroll), instead of shooting rays here.
 * Threads: OMP_NUM_THREADS.  No tuning flags: the policy is the TXR_* block of
 * pipeline_constants.h.
 * ==========================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/arena.h"
#include "../common/json_read.h"
#include "../common/raw_sample.h"
#include "../common/tiff_io.h"
#include "../common/ves_platform.h"
#include "../common/ves_png.h"
#include "../flatten/rawtex_tifxyz.h"
#include "../flatten/tifxyz_recto.h"

static int file_exists(const char *path)
{
    FILE *fp = fopen(path, "rb");
    if (fp == NULL) return 0;
    fclose(fp);
    return 1;
}

/* A path with forward slashes only, so it can sit inside a JSON string. */
static const char *json_path(const char *path, char *buf, size_t len)
{
    size_t i = 0;
    for (i = 0; i + 1 < len && path[i] != '\0'; i++) buf[i] = path[i] == '\\' ? '/' : path[i];
    buf[i] = '\0';
    return buf;
}

/* The chunk edge of a Zarr root's level 0 (its .zarray "chunks"), else 128. */
static long raw_chunk(Arena_T arena, const char *raw)
{
    char path[2600];
    const char *err = NULL;
    const JsonValue *j = NULL, *ch = NULL;
    snprintf(path, sizeof path, "%s/0/.zarray", raw);
    if (!file_exists(path)) return 128;
    j = Json_parse_file(arena, path, &err);
    if (j == NULL) return 128;
    ch = Json_object_get(j, "chunks");
    if (ch == NULL || Json_array_len(ch) < 1) return 128;
    return Json_as_long(Json_array_get(ch, 0), 128);
}

/* How far from an input node any CT sample can land: the prior's rays (+-TXR_PRIOR_RANGE_VOX and their smoothing),
 * or a node moved by the largest cap the prior may derive and searching the largest range from there. */
static double snap_reach(void)
{
    const double rays = TXR_PRIOR_RANGE_VOX + 3.0 * TXR_PRIOR_SIGMA;
    const double snap = (TXR_SEARCH_MAX_VOX - 1.0) + TXR_SEARCH_MAX_VOX + 3.0 * TXR_SIGMA_MAX;
    const double prof = (TXR_SEARCH_MAX_VOX - 1.0) + (double)TXR_PROFILE_HALF;
    double r = rays > snap ? rays : snap;
    if (prof > r) r = prof;
    return r + 2.0;
}

/* One plane's columns [a, b) of every row into dst ((b - a) x H): by rows when the TIFF is uncompressed 32-bit float
 * strips, else whole and cropped.  *W / *H: the plane's full size (checked against nonzero inputs). */
static int read_plane_cols(Arena_T io, const char *dir, const char *name, long a, long b, float **dst, long *W, long *H)
{
    char path[2600];
    TiffRowReader rd;
    long y = 0;
    snprintf(path, sizeof path, "%s/%s", dir, name);
    if (TiffIO_rows_open(io, path, &rd) == 0 && rd.px_bytes == 4) {
        if ((*W != 0 && rd.W != *W) || (*H != 0 && rd.H != *H)) {
            fprintf(stderr, "ERROR: %s is %d x %d, not %ld x %ld\n", path, rd.W, rd.H, *W, *H);
            TiffIO_rows_close(&rd);
            return -1;
        }
        *W = rd.W;
        *H = rd.H;
        if (b < 0) b = rd.W;
        if (a < 0 || b > rd.W || b <= a) {
            TiffIO_rows_close(&rd);
            return -1;
        }
        *dst = (float *)ARENA_ALLOC(io, (size_t)(b - a) * (size_t)rd.H * sizeof(float));
        for (y = 0; y < rd.H; y++)
            if (TiffIO_rows_read_f32(&rd, (int)y, (int)a, (int)(b - a), *dst + (size_t)y * (size_t)(b - a)) != 0) {
                fprintf(stderr, "ERROR: cannot read row %ld of %s\n", y, path);
                TiffIO_rows_close(&rd);
                return -1;
            }
        TiffIO_rows_close(&rd);
        return 0;
    }
    if (rd.f != NULL) TiffIO_rows_close(&rd);
    {
        float *img = NULL;
        int w = 0, h = 0;
        if (TiffIO_load_float2d(io, path, &img, &w, &h) != 0) {
            fprintf(stderr, "ERROR: cannot read %s\n", path);
            return -1;
        }
        if ((*W != 0 && w != *W) || (*H != 0 && h != *H)) {
            fprintf(stderr, "ERROR: %s is %d x %d, not %ld x %ld\n", path, w, h, *W, *H);
            return -1;
        }
        *W = w;
        *H = h;
        if (b < 0) b = w;
        if (a < 0 || b > w || b <= a) return -1;
        if (a == 0 && b == w) {
            *dst = img;
            return 0;
        }
        *dst = (float *)ARENA_ALLOC(io, (size_t)(b - a) * (size_t)h * sizeof(float));
        for (y = 0; y < h; y++)
            memcpy(*dst + (size_t)y * (size_t)(b - a), img + (size_t)y * (size_t)w + (size_t)a, (size_t)(b - a) * sizeof(float));
        return 0;
    }
}

/* The mask's columns [a, b) (NULL when the segment has none). */
static int read_mask_cols(Arena_T io, const char *dir, long a, long b, long W, long H, uint8_t **dst)
{
    char path[2600];
    TiffRowReader rd;
    long y = 0;
    *dst = NULL;
    snprintf(path, sizeof path, "%s/mask.tif", dir);
    if (!file_exists(path)) return 0;
    *dst = (uint8_t *)ARENA_ALLOC(io, (size_t)(b - a) * (size_t)H);
    if (TiffIO_rows_open(io, path, &rd) == 0 && rd.px_bytes == 1 && rd.W == W && rd.H == H) {
        for (y = 0; y < H; y++)
            if (TiffIO_rows_read(&rd, (int)y, (int)a, (int)(b - a), *dst + (size_t)y * (size_t)(b - a)) != 0) {
                TiffIO_rows_close(&rd);
                return -1;
            }
        TiffIO_rows_close(&rd);
        return 0;
    }
    if (rd.f != NULL) TiffIO_rows_close(&rd);
    {
        uint8_t *img = NULL;
        size_t mW = 0, mH = 0;
        if (RawtexTifxyz_load_gray(path, &img, &mW, &mH) != 0 || (long)mW != W || (long)mH != H) {
            fprintf(stderr, "ERROR: cannot read %s at %ld x %ld\n", path, W, H);
            free(img);
            return -1;
        }
        for (y = 0; y < H; y++)
            memcpy(*dst + (size_t)y * (size_t)(b - a), img + (size_t)y * mW + (size_t)a, (size_t)(b - a));
        free(img);
        return 0;
    }
}

/* A prior record (tifxyz-snap-prior-v1) back into the parameters the snap reads. */
static int read_prior(Arena_T arena, const char *path, TifxyzRectoPrior *pr)
{
    const char *err = NULL;
    const JsonValue *j = Json_parse_file(arena, path, &err), *d = NULL, *lv = NULL, *rays = NULL;
    if (j == NULL) return -1;
    d = Json_object_get(j, "derived");
    lv = Json_object_get(j, "levels");
    rays = Json_object_get(j, "rays");
    if (d == NULL || lv == NULL || rays == NULL) return -1;
    memset(pr, 0, sizeof *pr);
    pr->rays = (size_t)Json_member_long(rays, "shot", 0);
    pr->rays_with_band = (size_t)Json_member_long(rays, "with_band", 0);
    pr->rays_with_recto = (size_t)Json_member_long(rays, "with_recto", 0);
    pr->rays_with_next = (size_t)Json_member_long(rays, "with_next_wrap", 0);
    pr->rays_with_seam = (size_t)Json_member_long(rays, "with_seam", 0);
    pr->search_vox = Json_member_double(d, "search_vox", TXR_SEARCH_VOX);
    pr->max_total_vox = Json_member_double(d, "max_total_vox", TXR_MAX_TOTAL_VOX);
    pr->min_contrast = Json_member_double(d, "min_contrast", TXR_MIN_CONTRAST);
    pr->gate_error = Json_member_double(d, "gate_error", 0.0);
    pr->min_edge = Json_member_double(d, "min_edge", TXR_MIN_EDGE);
    pr->edge_full = Json_member_double(d, "edge_full", TXR_EDGE_FULL);
    pr->sigma = Json_member_double(d, "sigma", TXR_SMOOTH_SIGMA);
    pr->thick_min = Json_member_double(d, "thick_min", TXR_MIDLINE_MIN_THICK);
    pr->thick_max = Json_member_double(d, "thick_max", TXR_MIDLINE_MAX_THICK);
    pr->crowded = Json_as_bool(Json_object_get(d, "crowded"), 0);
    pr->air_level = Json_member_double(lv, "air", 0.0);
    pr->band_level = Json_member_double(lv, "band", 0.0);
    pr->band_min_level = Json_member_double(lv, "band_min", 0.0);
    pr->level_separability = Json_member_double(lv, "separability", 0.0);
    /* a record says whether the rays were good enough to be used; one that was not runs on the constants */
    if (!Json_as_bool(Json_object_get(d, "used"), 1)) pr->rays_with_recto = 0;
    return 0;
}

static void hist_json(FILE *fp, const char *key, const TifxyzRectoHist *h, const char *tail)
{
    fprintf(fp, "  \"%s\": {\"n\": %zu, \"p5\": %.3f, \"p25\": %.3f, \"p50\": %.3f, \"p75\": %.3f, \"p95\": %.3f, "
                "\"mean_abs\": %.3f, \"max_abs\": %.3f, \"share_over_1p5\": %.4f, \"share_over_3\": %.4f}%s\n",
            key, h->n, h->p5, h->p25, h->p50, h->p75, h->p95, h->mean_abs, h->max_abs, h->share_over_1p5,
            h->share_over_3, tail);
}

static void profile_json(FILE *fp, const char *key, const double *prof, const char *tail)
{
    int k = 0;
    fprintf(fp, "  \"%s\": [", key);
    for (k = 0; k < TXR_PROFILE_N; k++) {
        if (isfinite(prof[k])) fprintf(fp, "%s%.2f", k ? ", " : "", prof[k]);
        else fprintf(fp, "%snull", k ? ", " : "");
    }
    fprintf(fp, "]%s\n", tail);
}

/* ------------------------------------------------------------ windows' diagnostics */

typedef struct SnapWin {
    size_t k, lc0, ww, c0, c1;      /* loaded-grid columns: the window [lc0, lc0 + ww), its core [lc0 + c0, lc0 + c1) */
    size_t n_valid, n_supported, n_beyond, n_capped;
    double offset_med_first_p50, offset_med_last_p50, disp_p50, disp_p95, seconds;
    int iters;
} SnapWin;

/* The output frame (OW columns from loaded column col0) in cells of TXR_DIAG_SCALE: per raster the sum and count
 * of its finite values, and the valid nodes. */
typedef struct SnapDiag {
    size_t OW, H, cw, ch, col0;
    float *sum[5];
    uint8_t *cnt[5], *nvalid;
    SnapWin *win;
    size_t n_win, cap_win;
    const char *out;
    double t_last;
} SnapDiag;

enum { DG_FIRST, DG_LAST, DG_FIRST_MED, DG_LAST_MED, DG_DISP };

static void snap_window(void *ctx, size_t k, size_t n_windows, size_t lc0, size_t ww, size_t H, size_t c0, size_t c1,
                        const TifxyzRectoStats *st, const float *tf, const float *tl, const float *tfm, const float *tlm,
                        const float *dp)
{
    SnapDiag *d = (SnapDiag *)ctx;
    const float *ras[5] = { tf, tl, tfm, tlm, dp };
    const size_t s = TXR_DIAG_SCALE;
    const double now = ves_clock_sec();
    size_t r = 0, c = 0;
    int q = 0;
    for (r = 0; r < H; r++)
        for (c = c0; c < c1; c++) {
            const size_t i = r * ww + c, oc = lc0 + c - d->col0, cell = (r / s) * d->cw + oc / s;
            if (isfinite(dp[i])) d->nvalid[cell]++;
            for (q = 0; q < 5; q++)
                if (isfinite(ras[q][i])) {
                    d->sum[q][cell] += ras[q][i];
                    d->cnt[q][cell]++;
                }
        }
    if (d->n_win < d->cap_win) {
        SnapWin *w = &d->win[d->n_win++];
        w->k = k;
        w->lc0 = lc0;
        w->ww = ww;
        w->c0 = c0;
        w->c1 = c1;
        w->n_valid = st->n_valid;
        w->n_supported = st->n_supported_first;
        w->n_beyond = st->n_beyond_first;
        w->n_capped = st->n_capped_total;
        w->offset_med_first_p50 = st->target_first_med.p50;
        w->offset_med_last_p50 = st->target_last_med.p50;
        w->disp_p50 = st->disp.p50;
        w->disp_p95 = st->disp.p95;
        w->seconds = now - d->t_last;
        w->iters = st->outer_iters;
    }
    fprintf(stderr, "[tifxyz_snap] window %zu/%zu: core columns [%zu, %zu) of the output (loaded [%zu, %zu)), %zu valid, "
                    "supported %.1f%% (beyond reach %.1f%%), median offset %+.2f -> %+.2f, moved p50 %+.2f p95 %+.2f, "
                    "capped %zu, %d iters (%.1fs)\n",
            k + 1, n_windows, lc0 + c0 - d->col0, lc0 + c1 - d->col0, lc0, lc0 + ww, st->n_valid,
            st->n_valid ? 100.0 * (double)st->n_supported_first / (double)st->n_valid : 0.0,
            st->n_valid ? 100.0 * (double)st->n_beyond_first / (double)st->n_valid : 0.0, st->target_first_med.p50,
            st->target_last_med.p50, st->disp.p50, st->disp.p95, st->n_capped_total, st->outer_iters, now - d->t_last);
    d->t_last = now;
    /* one window covers the output: the rasters in full, as float TIFFs (NaN = undefined), for analysis */
    if (n_windows == 1) {
        const char *names[3] = { "_recto_before.tif", "_recto_after.tif", "_disp.tif" };
        const float *src[3] = { tf, tl, dp };
        float *buf = (float *)malloc(d->OW * H * sizeof(float));
        char path[2600];
        if (buf == NULL) return;
        for (q = 0; q < 3; q++) {
            for (r = 0; r < H; r++)
                memcpy(buf + r * d->OW, src[q] + r * ww + (d->col0 - lc0), d->OW * sizeof(float));
            snprintf(path, sizeof path, "%s%s", d->out, names[q]);
            if (TiffIO_save_float2d(path, buf, (int)d->OW, (int)H) != 0) fprintf(stderr, "WARNING: cannot write %s\n", path);
        }
        free(buf);
    }
}

/* A raster's cells as a diverging PNG. */
static int write_cells_png(const char *path, const SnapDiag *d, int q)
{
    uint8_t *img = (uint8_t *)calloc(d->cw * d->ch * 3, 1);
    size_t cell = 0;
    int rc = -1;
    if (img == NULL) return -1;
    for (cell = 0; cell < d->cw * d->ch; cell++) {
        uint8_t *px = img + cell * 3;
        double f = 0.0;
        if (d->nvalid[cell] == 0) continue;                  /* black: invalid */
        if (d->cnt[q][cell] == 0) {
            px[0] = px[1] = px[2] = 48;                      /* dark grey: nothing measured */
            continue;
        }
        f = ((double)d->sum[q][cell] / (double)d->cnt[q][cell]) / TXR_DIAG_FULL_VOX;
        if (f > 1.0) f = 1.0;
        if (f < -1.0) f = -1.0;
        if (f >= 0.0) {                                      /* light grey -> red */
            px[0] = (uint8_t)(235.0 + 20.0 * f);
            px[1] = (uint8_t)(235.0 - 160.0 * f);
            px[2] = (uint8_t)(235.0 - 185.0 * f);
        } else {                                             /* light grey -> blue */
            px[0] = (uint8_t)(235.0 - 185.0 * -f);
            px[1] = (uint8_t)(235.0 - 130.0 * -f);
            px[2] = (uint8_t)(235.0 + 20.0 * -f);
        }
    }
    rc = VesPng_write_rgb(path, img, (int)d->cw, (int)d->ch);
    free(img);
    return rc;
}

/* ------------------------------------------------------------ the prior's record */

static int write_prior_json(const char *path, const char *in, const char *raw, const TifxyzRectoPrior *prior, int no_prior,
                            double seconds)
{
    FILE *fp = fopen(path, "wb");
    int k = 0;
    if (fp == NULL) return -1;
    fprintf(fp, "{\n  \"schema\": \"tifxyz-snap-prior-v1\",\n  \"in\": \"%s\",\n  \"raw\": \"%s\",\n"
                "  \"rays\": {\"stride\": %d, \"range_vox\": %.1f, \"smoothing_sigma\": %.2f, \"floor_grey\": %.1f, \"air_rel\": %.2f, "
                "\"shot\": %zu, \"with_band\": %zu, \"with_recto\": %zu, \"with_next_wrap\": %zu, \"with_seam\": %zu},\n",
            in, raw, TXR_PRIOR_STRIDE, TXR_PRIOR_RANGE_VOX, TXR_PRIOR_SIGMA, TXR_PRIOR_FLOOR, TXR_PRIOR_AIR_REL, prior->rays,
            prior->rays_with_band, prior->rays_with_recto, prior->rays_with_next, prior->rays_with_seam);
    {
        const char *names[8] = { "recto_offset_vox", "wrap_spacing_vox", "band_thickness_vox", "gap_width_vox", "gap_depth_grey",
                                 "seam_depth_grey", "recto_slope_grey_per_vox", "edge_width_vox" };
        const TifxyzRectoQ *qs[8] = { &prior->recto_offset, &prior->spacing, &prior->thickness, &prior->gap_width,
                                      &prior->gap_depth, &prior->seam_depth, &prior->slope, &prior->edge_width };
        for (k = 0; k < 8; k++)
            fprintf(fp, "  \"%s\": {\"n\": %zu, \"p01\": %.3f, \"p05\": %.3f, \"p25\": %.3f, \"p50\": %.3f, \"p75\": %.3f, "
                        "\"p95\": %.3f, \"p99\": %.3f},\n",
                    names[k], qs[k]->n, qs[k]->p01, qs[k]->p05, qs[k]->p25, qs[k]->p50, qs[k]->p75, qs[k]->p95, qs[k]->p99);
    }
    fprintf(fp, "  \"levels\": {\"air\": %.2f, \"band\": %.2f, \"band_min\": %.2f, \"separability\": %.3f, \"trough_p05_p50_p95\": [%.1f, %.1f, %.1f], "
                "\"crest_p05_p50_p95\": [%.1f, %.1f, %.1f]},\n",
            prior->air_level, prior->band_level, prior->band_min_level, prior->level_separability, prior->trough_level.p05,
            prior->trough_level.p50, prior->trough_level.p95, prior->crest_level.p05, prior->crest_level.p50, prior->crest_level.p95);
    fprintf(fp, "  \"trough_widths\": {\"gap_p05_p50_p95\": [%.2f, %.2f, %.2f], \"seam_p05_p50_p95\": [%.2f, %.2f, %.2f], "
                "\"split\": %.2f, \"separability\": %.3f, \"error_vs_levels\": %.3f, \"n_gap\": %zu, \"n_seam\": %zu},\n",
            prior->gap_trough_width.p05, prior->gap_trough_width.p50, prior->gap_trough_width.p95, prior->seam_trough_width.p05,
            prior->seam_trough_width.p50, prior->seam_trough_width.p95, prior->trough_width_split, prior->trough_width_separability,
            prior->trough_width_error, prior->gap_trough_width.n, prior->seam_trough_width.n);
    fprintf(fp, "  \"derived\": {\"search_vox\": %.3f, \"max_total_vox\": %.3f, \"min_contrast\": %.2f, \"gate_error\": %.4f, "
                "\"min_edge\": %.2f, \"edge_full\": %.2f, \"sigma\": %.3f, \"thick_min\": %.2f, \"thick_max\": %.2f, \"crowded\": %s, "
                "\"used\": %s},\n",
            prior->search_vox, prior->max_total_vox, prior->min_contrast, prior->gate_error, prior->min_edge, prior->edge_full,
            prior->sigma, prior->thick_min, prior->thick_max, prior->crowded ? "true" : "false",
            !no_prior && prior->rays_with_recto >= 100 ? "true" : "false");
    fprintf(fp, "  \"rules\": {\"search_vox\": \"max(p95 |recto offset| + 1, 4), at most 0.9 x p05 wrap spacing and %.0f\", "
                "\"max_total_vox\": \"search - 1\", \"min_contrast\": \"the depth that best separates seam dips (inside a band) "
                "from gaps (air, within %.2f of the ray's contrast above its darkest sample); 0.3 x p05 gap depth when no seams; "
                "6..80\", \"min_edge\": \"p05 recto slope, 2..30\", \"edge_full\": \"p50 recto slope\", \"sigma\": \"sqrt((p50 edge "
                "width / 2.56)^2 - %.2f^2), 0.4..%.1f\", \"thick_min\": \"p01 thickness - 0.5, 1..6\", \"thick_max\": \"p99 thickness "
                "+ 1, 8..24\"},\n  \"seconds\": %.1f\n}\n",
            TXR_SEARCH_MAX_VOX, TXR_PRIOR_AIR_REL, TXR_PRIOR_SIGMA, TXR_SIGMA_MAX, seconds);
    fclose(fp);
    return 0;
}

static void log_prior(const TifxyzRectoPrior *prior, int no_prior, const char *from, double seconds)
{
    fprintf(stderr, "[tifxyz_snap] PRIOR %s%zu rays (stride %d, +-%.0f vox): levels air %.0f, band %.0f (separability %.2f; troughs "
                    "p05/p50/p95 %.0f/%.0f/%.0f, crests %.0f/%.0f/%.0f); band on %.1f%%, recto %.1f%%, next wrap %.1f%%, seam "
                    "%.1f%%; recto offset p50 %+.2f (p5 %+.2f, p95 %+.2f), spacing p05/p50 %.1f/%.1f, thickness p50 %.1f "
                    "(p01 %.1f, p99 %.1f), gap width p50 %.1f, gap depth p05/p50 %.0f/%.0f, seam depth p50/p95 %.0f/%.0f, "
                    "slope p05/p50 %.0f/%.0f, edge width p50 %.2f -> search %.2f, cap %.2f, gate %.0f (error %.3f), edge "
                    "%.0f/%.0f, sigma %.2f, thick %.1f..%.1f%s%s (%.1fs)\n",
            from, prior->rays, TXR_PRIOR_STRIDE, TXR_PRIOR_RANGE_VOX, prior->air_level, prior->band_level,
            prior->level_separability, prior->trough_level.p05, prior->trough_level.p50, prior->trough_level.p95,
            prior->crest_level.p05, prior->crest_level.p50, prior->crest_level.p95,
            prior->rays ? 100.0 * (double)prior->rays_with_band / (double)prior->rays : 0.0,
            prior->rays ? 100.0 * (double)prior->rays_with_recto / (double)prior->rays : 0.0,
            prior->rays ? 100.0 * (double)prior->rays_with_next / (double)prior->rays : 0.0,
            prior->rays ? 100.0 * (double)prior->rays_with_seam / (double)prior->rays : 0.0,
            prior->recto_offset.p50, prior->recto_offset.p05, prior->recto_offset.p95, prior->spacing.p05, prior->spacing.p50,
            prior->thickness.p50, prior->thickness.p01, prior->thickness.p99, prior->gap_width.p50, prior->gap_depth.p05,
            prior->gap_depth.p50, prior->seam_depth.p50, prior->seam_depth.p95, prior->slope.p05, prior->slope.p50,
            prior->edge_width.p50, prior->search_vox, prior->max_total_vox, prior->min_contrast, prior->gate_error,
            prior->min_edge, prior->edge_full, prior->sigma, prior->thick_min, prior->thick_max,
            prior->crowded ? " CROWDED (the recto's reach exceeds the wrap spacing)" : "",
            no_prior ? " -- NOT USED (--no-prior)" : prior->rays_with_recto < 100 ? " -- too few rays: the constants run" : "",
            seconds);
}

/* ------------------------------------------------------------ main */

int main(int argc, char **argv)
{
    Arena_T arena = NULL, io = NULL;
    const char *in = NULL, *raw = NULL, *out = NULL, *err = NULL, *prior_path = NULL;
    char path[2600], uuid[600], extra[4096], jin[2048], jraw[2048], jout[2048], jprior[2048];
    float *pos = NULL, *pos_out = NULL, *pl[3] = { NULL, NULL, NULL };
    uint8_t *mask = NULL;
    size_t W = 0, H = 0, n = 0, i = 0, n_valid = 0, OW = 0, col0 = 0, col1 = 0;
    long full_w = 0, full_h = 0, win_c0 = 0, win_c1 = 0, load_c0 = 0, load_c1 = 0;
    int a = 0, k = 0, windowed = 0, no_prior = 0, prior_only = 0, use_prior = 0;
    float corners[6];
    double lo[3] = { 1e30, 1e30, 1e30 }, hi[3] = { -1e30, -1e30, -1e30 }, t0 = 0.0, t_load = 0.0, t_raw = 0.0;
    double t_prior = 0.0, reach = snap_reach();
    CubeTable ct;
    TifxyzRectoOpts o;
    TifxyzRectoStats st;
    TifxyzRectoPrior prior;
    SnapDiag diag;
    const JsonValue *meta = NULL;
    FILE *fp = NULL;

    if (argc >= 2 && strcmp(argv[1], "--selftest") == 0) {
        char tiff_dir[256];
        int fails = TifxyzRecto_selftest();
        fails += RawtexTifxyz_selftest() != 0;
        ves_selftest_dir(tiff_dir, sizeof tiff_dir, "tifxyz_snap_tiff");
        if (TiffIO_selftest(tiff_dir) != 0) fails++;
        else ves_selftest_clean(tiff_dir);
        fprintf(stderr, "[selftest] tifxyz_snap %s (%d failure%s)\n", fails == 0 ? "PASS" : "FAIL", fails,
                fails == 1 ? "" : "s");
        return fails == 0 ? 0 : 1;
    }
    if (argc < 4) {
        fprintf(stderr, "usage: tifxyz_snap <in_seg_dir> <raw_root> <out_seg_dir> [--midline | --inset] [--smooth]\n"
                        "                   [--no-prior | --prior <prior.json>] [--prior-only] [--window C0 C1]\n"
                        "       tifxyz_snap --selftest\n"
                        "  --inset:      the recto edge minus TXR_RECTO_INSET_VOX (the surface stops inside the papyrus)\n"
                        "  --smooth:     the TXR_MEDIAN_HALF_SMOOTH window (the slow drift only, no fibre weave)\n"
                        "  --no-prior:   run on the TXR_* constants instead of the parameters the rays measure\n"
                        "  --prior:      run on the parameters of that record (another segment's: one per scroll)\n"
                        "  --prior-only: shoot the rays, write <out_seg_dir>_prior.json and stop (no segment)\n"
                        "  --window:     only the columns [C0, C1), read by rows with TXR_WINDOW_HALO columns either side\n");
        return 1;
    }
    in = json_path(argv[1], jin, sizeof jin);
    raw = json_path(argv[2], jraw, sizeof jraw);
    out = json_path(argv[3], jout, sizeof jout);
    memset(&o, 0, sizeof o);
    o.verbose = 1;
    for (a = 4; a < argc; a++) {
        if (strcmp(argv[a], "--midline") == 0) o.midline = 1;
        else if (strcmp(argv[a], "--inset") == 0) o.inset = 1;
        else if (strcmp(argv[a], "--smooth") == 0) o.smooth = 1;
        else if (strcmp(argv[a], "--no-prior") == 0) no_prior = 1;
        else if (strcmp(argv[a], "--prior-only") == 0) prior_only = 1;
        else if (strcmp(argv[a], "--prior") == 0 && a + 1 < argc) prior_path = json_path(argv[++a], jprior, sizeof jprior);
        else if (strcmp(argv[a], "--window") == 0 && a + 2 < argc) {
            win_c0 = atol(argv[++a]);
            win_c1 = atol(argv[++a]);
            windowed = 1;
        } else {
            fprintf(stderr, "unknown argument %s\n", argv[a]);
            return 1;
        }
    }
    if (no_prior && prior_path != NULL) {
        fprintf(stderr, "ERROR: --no-prior and --prior exclude each other\n");
        return 1;
    }
    if (windowed && (win_c0 < 0 || win_c1 <= win_c0)) {
        fprintf(stderr, "ERROR: --window %ld %ld is empty\n", win_c0, win_c1);
        return 1;
    }
    snprintf(path, sizeof path, "%s/meta.json", out);
    if (!prior_only && file_exists(path)) {
        fprintf(stderr, "ERROR: %s already holds a segment; refusing to replace it\n", out);
        return 1;
    }
    t0 = ves_clock_sec();
    arena = Arena_new();
    /* the records beside the segment (<out>_prior.json, _snap_report.json, the PNGs) need its parent folder */
    snprintf(path, sizeof path, "%s_prior.json", out);
    ves_ensure_parent_dir(path);

    /* the grid: the window's columns with the halo, or all of them; plane-major (z, y, x), -1 = invalid */
    io = Arena_new();
    {
        const char *names[3] = { "z.tif", "y.tif", "x.tif" };
        long a0 = -1, a1 = -1;
        if (windowed) {
            /* the plane's width first (the row reader opens without reading) */
            TiffRowReader rd;
            snprintf(path, sizeof path, "%s/x.tif", in);
            if (TiffIO_rows_open(io, path, &rd) == 0) {
                full_w = rd.W;
                full_h = rd.H;
                TiffIO_rows_close(&rd);
            } else {
                float *img = NULL;
                int w = 0, h = 0;
                if (TiffIO_load_float2d(io, path, &img, &w, &h) != 0) {
                    fprintf(stderr, "ERROR: cannot read %s\n", path);
                    return 1;
                }
                full_w = w;
                full_h = h;
            }
            if (win_c1 > full_w) {
                fprintf(stderr, "ERROR: --window %ld %ld beyond the %ld columns of %s\n", win_c0, win_c1, full_w, in);
                return 1;
            }
            a0 = win_c0 - TXR_WINDOW_HALO > 0 ? win_c0 - TXR_WINDOW_HALO : 0;
            a1 = win_c1 + TXR_WINDOW_HALO < full_w ? win_c1 + TXR_WINDOW_HALO : full_w;
        }
        for (k = 0; k < 3; k++)
            if (read_plane_cols(io, in, names[k], windowed ? a0 : 0, windowed ? a1 : -1, &pl[k], &full_w, &full_h) != 0) {
                fprintf(stderr, "ERROR: cannot read %s/%s\n", in, names[k]);
                return 1;
            }
        load_c0 = windowed ? a0 : 0;
        load_c1 = windowed ? a1 : full_w;
        if (read_mask_cols(io, in, load_c0, load_c1, full_w, full_h, &mask) != 0) return 1;
    }
    W = (size_t)(load_c1 - load_c0);
    H = (size_t)full_h;
    n = W * H;
    col0 = windowed ? (size_t)(win_c0 - load_c0) : 0;
    col1 = windowed ? (size_t)(win_c1 - load_c0) : W;
    OW = col1 - col0;
    pos = (float *)ARENA_ALLOC(arena, 3 * n * sizeof(float));
    for (i = 0; i < n; i++) {
        const int ok = pl[0][i] > -0.5f && pl[1][i] > -0.5f && pl[2][i] > -0.5f && (mask == NULL || mask[i] != 0);
        for (k = 0; k < 3; k++) pos[(size_t)k * n + i] = ok ? pl[k][i] : -1.0f;
        if (!ok) continue;
        n_valid++;
        for (k = 0; k < 3; k++) {
            const double v = (double)pos[(size_t)k * n + i];
            if (v < lo[k]) lo[k] = v;
            if (v > hi[k]) hi[k] = v;
        }
    }
    Arena_dispose(&io);
    mask = NULL;
    snprintf(path, sizeof path, "%s/meta.json", in);
    meta = Json_parse_file(arena, path, &err);
    snprintf(uuid, sizeof uuid, "%s_recto", meta != NULL && Json_as_string(Json_object_get(meta, "uuid")) != NULL
                                                 ? Json_as_string(Json_object_get(meta, "uuid")) : ves_path_basename(in));
    t_load = ves_clock_sec() - t0;
    fprintf(stderr, "[tifxyz_snap] %s: %ld x %ld nodes; columns [%ld, %ld) loaded, [%ld, %ld) snapped; %zu valid, bbox z %.0f..%.0f "
                    "y %.0f..%.0f x %.0f..%.0f (%.1fs)\n",
            in, full_w, full_h, load_c0, load_c1, load_c0 + (long)col0, load_c0 + (long)col1, n_valid, lo[0], hi[0], lo[1],
            hi[1], lo[2], hi[2], t_load);
    if (n_valid == 0) {
        fprintf(stderr, "ERROR: no valid node in %s\n", in);
        return 1;
    }

    /* the CT: the bbox's table, loaded only where some node reaches */
    for (k = 0; k < 3; k++) {
        corners[k] = (float)lo[k];
        corners[3 + k] = (float)hi[k];
    }
    if (cubetable_init(&ct, arena, raw, raw_chunk(arena, raw), corners, 2, reach) != 0) {
        fprintf(stderr, "ERROR: cube table init failed for %s\n", raw);
        return 1;
    }
    t_raw = ves_clock_sec();
    if (cubetable_prewarm_points(&ct, pos, pos + n, pos + 2 * n, n, reach) < 0) {
        fprintf(stderr, "ERROR: the CT prewarm failed\n");
        return 1;
    }
    t_raw = ves_clock_sec() - t_raw;
    fprintf(stderr, "[tifxyz_snap] RAW %s: %d chunks loaded (%.1f GB), %d missing, %zu out of reach (reach %.1f vox; %.1fs)\n",
            raw, ct.n_loaded, (double)ct.n_loaded * (double)(ct.chunk * ct.chunk * ct.chunk) / 1e9, ct.n_missing,
            ct.n_untouched, reach, t_raw);

    /* the prior: the rays read the wrap structure and set the parameters -- or a record of that */
    memset(&prior, 0, sizeof prior);
    t_prior = ves_clock_sec();
    if (prior_path != NULL) {
        if (read_prior(arena, prior_path, &prior) != 0) {
            fprintf(stderr, "ERROR: cannot read the prior %s\n", prior_path);
            return 1;
        }
    } else if (TifxyzRecto_prior(arena, &ct, pos, W, H, TXR_PRIOR_STRIDE, &prior) != 0) {
        fprintf(stderr, "ERROR: the prior failed\n");
        return 1;
    }
    t_prior = ves_clock_sec() - t_prior;
    log_prior(&prior, no_prior, prior_path != NULL ? "(from the record) " : "", t_prior);
    if (prior_path == NULL) {
        fprintf(stderr, "[tifxyz_snap] PRIOR trough widths at half depth: gaps (fallen to the air level) p05/p50/p95 %.2f/%.2f/%.2f "
                        "(%zu), seams p05/p50/p95 %.2f/%.2f/%.2f (%zu); the narrow/wide split %.2f vox (separability %.2f, "
                        "disagreement with the levels %.3f)\n",
                prior.gap_trough_width.p05, prior.gap_trough_width.p50, prior.gap_trough_width.p95, prior.gap_trough_width.n,
                prior.seam_trough_width.p05, prior.seam_trough_width.p50, prior.seam_trough_width.p95, prior.seam_trough_width.n,
                prior.trough_width_split, prior.trough_width_separability, prior.trough_width_error);
        snprintf(path, sizeof path, "%s_prior.json", out);
        if (write_prior_json(path, in, raw, &prior, no_prior, t_prior) != 0) {
            fprintf(stderr, "ERROR: cannot write %s\n", path);
            return 1;
        }
        if (prior_only) {
            fprintf(stderr, "[tifxyz_snap] prior written: %s\n", path);
            Arena_dispose(&arena);
            return 0;
        }
    } else if (prior_only) {
        fprintf(stderr, "ERROR: --prior-only with --prior has nothing to do\n");
        return 1;
    }
    use_prior = !no_prior && prior.rays_with_recto >= 100;

    /* the snap, in windows of the output's columns */
    memset(&diag, 0, sizeof diag);
    diag.OW = OW;
    diag.H = H;
    diag.cw = (OW + TXR_DIAG_SCALE - 1) / TXR_DIAG_SCALE;
    diag.ch = (H + TXR_DIAG_SCALE - 1) / TXR_DIAG_SCALE;
    diag.col0 = col0;
    diag.out = out;
    diag.t_last = ves_clock_sec();
    for (k = 0; k < 5; k++) {
        diag.sum[k] = (float *)ARENA_CALLOC(arena, diag.cw * diag.ch, sizeof(float));
        diag.cnt[k] = (uint8_t *)ARENA_CALLOC(arena, diag.cw * diag.ch, 1);
    }
    diag.nvalid = (uint8_t *)ARENA_CALLOC(arena, diag.cw * diag.ch, 1);
    diag.cap_win = OW / TXR_WINDOW_COLS + 2;
    diag.win = (SnapWin *)ARENA_CALLOC(arena, diag.cap_win, sizeof(SnapWin));
    pos_out = (float *)ARENA_ALLOC(arena, 3 * n * sizeof(float));
    if (TifxyzRecto_run_windows(arena, &ct, pos, pos_out, W, H, col0, col1, TXR_WINDOW_COLS, TXR_WINDOW_HALO, &o,
                                use_prior ? &prior : NULL, &st, snap_window, &diag) != 0) {
        fprintf(stderr, "ERROR: the snap failed\n");
        return 1;
    }

    /* the segment: the output's columns of the moved grid; invalid nodes stay invalid */
    snprintf(extra, sizeof extra,
             "\"source_segment\": \"%s\", \"raw\": \"%s\", \"snap\": {\"tool\": \"tifxyz_snap\", \"target\": \"%s\", "
             "\"outer_iters\": %d, \"supported_share\": %.4f, \"disp_p50\": %.3f, \"disp_p95_abs\": %.3f, "
             "\"recto_offset_p50_before\": %.3f, \"recto_offset_p50_after\": %.3f, \"search_vox\": %.2f, "
             "\"smooth_mu\": %.2f, \"max_total_vox\": %.2f, \"parameters\": \"%s\", \"windows\": %d}",
             in, raw, o.midline ? "midline" : o.inset ? "recto-inset" : "recto", st.outer_iters,
             st.n_valid > 0 ? (double)st.n_supported_first / (double)st.n_valid : 0.0, st.disp.p50,
             st.disp.max_abs, st.target_first.p50, st.target_last.p50, st.prm_search, TXR_SMOOTH_MU,
             st.prm_max_total, use_prior ? (prior_path != NULL ? prior_path : "prior") : "constants", st.windows);
    if (windowed)
        snprintf(extra + strlen(extra), sizeof extra - strlen(extra),
                 ", \"window\": {\"col0\": %ld, \"col1\": %ld, \"source_width\": %ld, \"halo\": %d, "
                 "\"rule\": \"column c of this segment = column col0 + c of the source segment\"}",
                 win_c0, win_c1, full_w, TXR_WINDOW_HALO);
    {
        float *wr = pos_out;
        if (col0 != 0 || col1 != W) {
            size_t r = 0;
            wr = (float *)ARENA_ALLOC(arena, 3 * OW * H * sizeof(float));
            for (k = 0; k < 3; k++)
                for (r = 0; r < H; r++)
                    memcpy(wr + (size_t)k * OW * H + r * OW, pos_out + (size_t)k * n + r * W + col0, OW * sizeof(float));
        }
        for (i = 0; i < 3 * OW * H; i++)
            if (wr[i] < -0.5f) wr[i] = NAN;      /* the writer's invalid marker */
        if (RawtexTifxyz_write(out, uuid, wr, OW, H, extra) != 0) return 1;
    }

    /* the report */
    snprintf(path, sizeof path, "%s_snap_report.json", out);
    fp = fopen(path, "wb");
    if (fp == NULL) {
        fprintf(stderr, "ERROR: cannot write %s\n", path);
        return 1;
    }
    fprintf(fp, "{\n  \"schema\": \"tifxyz-snap-report-v1\",\n  \"tool\": \"tifxyz_snap\",\n  \"target\": \"%s\",\n"
                "  \"in\": \"%s\",\n  \"raw\": \"%s\",\n  \"out\": \"%s\",\n  \"W\": %zu,\n  \"H\": %zu,\n"
                "  \"n_valid\": %zu,\n  \"n_dir\": %zu,\n",
            o.midline ? "midline" : o.inset ? "recto-inset" : "recto", in, raw, out, OW, H, st.n_valid, st.n_dir);
    fprintf(fp, "  \"frame\": {\"source_width\": %ld, \"loaded\": [%ld, %ld], \"columns\": [%ld, %ld], \"windowed\": %s},\n",
            full_w, load_c0, load_c1, load_c0 + (long)col0, load_c0 + (long)col1, windowed ? "true" : "false");
    fprintf(fp, "  \"inset_vox\": %.3f,\n  \"smooth\": %s,\n", o.inset ? TXR_RECTO_INSET_VOX : 0.0, o.smooth ? "true" : "false");
    fprintf(fp, "  \"first\": {\"supported\": %zu, \"edge_inward\": %zu, \"edge_outward\": %zu, \"band_beyond_reach\": %zu},\n"
                "  \"last\": {\"supported\": %zu, \"edge_inward\": %zu, \"edge_outward\": %zu},\n",
            st.n_supported_first, st.n_in_band_first, st.n_in_gap_first, st.n_beyond_first, st.n_supported_last,
            st.n_in_band_last, st.n_in_gap_last);
    fprintf(fp, "  \"outer_iters\": %d,\n  \"mean_abs_step\": [", st.outer_iters);
    for (k = 0; k < st.outer_iters; k++) fprintf(fp, "%s%.4f", k ? ", " : "", st.mean_abs_step[k]);
    fprintf(fp, "],\n  \"pcg_iters\": [");
    for (k = 0; k < st.outer_iters; k++) fprintf(fp, "%s%d", k ? ", " : "", st.pcg_iters[k]);
    fprintf(fp, "],\n  \"pcg_resid\": [");
    for (k = 0; k < st.outer_iters; k++) fprintf(fp, "%s%.2e", k ? ", " : "", st.pcg_resid[k]);
    fprintf(fp, "],\n  \"n_capped_total\": %zu,\n", st.n_capped_total);
    hist_json(fp, "recto_offset_before", &st.target_first, ",");
    hist_json(fp, "recto_offset_after", &st.target_last, ",");
    fprintf(fp, "  \"median_window\": {\"half\": %d, \"filtered_first\": %zu, \"filtered_last\": %zu},\n",
            o.smooth ? TXR_MEDIAN_HALF_SMOOTH : TXR_MEDIAN_HALF, st.n_filtered_first, st.n_filtered_last);
    hist_json(fp, "recto_offset_before_median", &st.target_first_med, ",");
    hist_json(fp, "recto_offset_after_median", &st.target_last_med, ",");
    hist_json(fp, "displacement", &st.disp, ",");
    fprintf(fp, "  \"slope\": {\"p99\": %.4f, \"max\": %.4f},\n", st.slope_p99, st.slope_max);
    fprintf(fp, "  \"profile_vox\": {\"half\": %d, \"step\": 1, \"stride\": %d, \"nodes\": %zu, \"sign\": \"+ inward\"},\n",
            TXR_PROFILE_HALF, TXR_PROFILE_STRIDE, st.profile_n);
    profile_json(fp, "profile_before", st.profile_first, ",");
    profile_json(fp, "profile_after", st.profile_last, ",");
    fprintf(fp, "  \"parameters_used\": {\"source\": \"%s\", \"prior_record\": %s%s%s, \"search_vox\": %.3f, \"max_total_vox\": %.3f, "
                "\"min_contrast\": %.2f, \"min_edge\": %.2f, \"edge_full\": %.2f, \"sigma\": %.3f, \"thick_min\": %.2f, "
                "\"thick_max\": %.2f, \"band_rule\": \"%s\", \"air_level\": %.2f, \"band_min_level\": %.2f},\n",
            use_prior ? "prior" : "constants", prior_path != NULL ? "\"" : "", prior_path != NULL ? prior_path : "null",
            prior_path != NULL ? "\"" : "", st.prm_search, st.prm_max_total, st.prm_min_contrast, st.prm_min_edge,
            st.prm_edge_full, st.prm_sigma, st.prm_thick_min, st.prm_thick_max,
            st.prm_levels ? "a band ends where a trough falls to the air level" : "a band ends at a trough the gate below its crest",
            st.prm_air_level, st.prm_band_min_level);
    fprintf(fp, "  \"constants\": {\"search_vox\": %.2f, \"sample_step\": %.2f, \"smooth_sigma\": %.2f, "
                "\"min_contrast\": %.1f, \"min_edge\": %.1f, \"edge_full\": %.1f, \"data_w\": %.2f, "
                "\"anchor_w\": %.3f, \"smooth_mu\": %.2f, \"irls_rounds\": %d, \"robust_scale_vox\": %.2f, "
                "\"median_half\": %d, \"median_min_support\": %.2f, \"max_step_vox\": %.2f, "
                "\"max_total_vox\": %.2f, \"outer_iters\": %d, \"step_tol\": %.3f, \"pcg_tol\": %.1e, "
                "\"window_cols\": %d, \"window_halo\": %d},\n",
            TXR_SEARCH_VOX, TXR_SAMPLE_STEP, TXR_SMOOTH_SIGMA, TXR_MIN_CONTRAST, TXR_MIN_EDGE, TXR_EDGE_FULL,
            TXR_DATA_W, TXR_ANCHOR_W, TXR_SMOOTH_MU, TXR_IRLS_ROUNDS, TXR_ROBUST_SCALE_VOX, TXR_MEDIAN_HALF,
            TXR_MEDIAN_MIN_SUPPORT, TXR_MAX_STEP_VOX, TXR_MAX_TOTAL_VOX, TXR_OUTER_ITERS, TXR_STEP_TOL, TXR_PCG_TOL,
            TXR_WINDOW_COLS, TXR_WINDOW_HALO);
    fprintf(fp, "  \"windows\": {\"run\": %d, \"list\": [", st.windows);
    for (i = 0; i < diag.n_win; i++) {
        const SnapWin *w = &diag.win[i];
        fprintf(fp, "%s\n    {\"k\": %zu, \"loaded\": [%zu, %zu], \"core\": [%zu, %zu], \"valid\": %zu, \"supported\": %zu, "
                    "\"beyond_reach\": %zu, \"capped\": %zu, \"offset_median_p50_before\": %.3f, \"offset_median_p50_after\": %.3f, "
                    "\"disp_p50\": %.3f, \"disp_p95\": %.3f, \"iters\": %d, \"seconds\": %.1f}",
                i ? "," : "", w->k, (size_t)load_c0 + w->lc0, (size_t)load_c0 + w->lc0 + w->ww, (size_t)load_c0 + w->lc0 + w->c0,
                (size_t)load_c0 + w->lc0 + w->c1, w->n_valid, w->n_supported, w->n_beyond, w->n_capped,
                w->offset_med_first_p50, w->offset_med_last_p50, w->disp_p50, w->disp_p95, w->iters, w->seconds);
    }
    fprintf(fp, "\n  ]},\n");
    fprintf(fp, "  \"raw_chunks\": {\"loaded\": %d, \"missing\": %d, \"out_of_reach\": %zu, \"reach_vox\": %.1f},\n",
            ct.n_loaded, ct.n_missing, ct.n_untouched, reach);
    fprintf(fp, "  \"seconds\": {\"load\": %.1f, \"raw\": %.1f, \"prior\": %.1f, \"measure\": %.1f, \"solve\": %.1f, "
                "\"snap_total\": %.1f, \"total\": %.1f}\n}\n",
            t_load, t_raw, t_prior, st.sec_measure, st.sec_solve, st.sec_total, ves_clock_sec() - t0);
    fclose(fp);

    {
        const char *names[5] = { "_recto_before.png", "_recto_after.png", "_recto_before_median.png", "_recto_after_median.png",
                                 "_disp.png" };
        for (k = 0; k < 5; k++) {
            snprintf(path, sizeof path, "%s%s", out, names[k]);
            if (write_cells_png(path, &diag, k) != 0) fprintf(stderr, "WARNING: cannot write %s\n", path);
        }
    }

    fprintf(stderr, "[tifxyz_snap] SUMMARY %s: %zu valid nodes in %d window%s (search %.2f, cap %.2f, gate %.0f, sigma %.2f from the %s), "
                    "supported %.1f%% (edge inward %.1f%%, outward %.1f%%; abstained, band beyond reach %.1f%%); recto offset per "
                    "node before p50 %+.2f (p5 %+.2f, p95 %+.2f; >1.5 vox %.1f%%) -> after p50 %+.2f (p5 %+.2f, p95 "
                    "%+.2f; >1.5 vox %.1f%%); %d-node median before p50 %+.2f (p5 %+.2f, p95 %+.2f; >1.5 vox %.1f%%) "
                    "-> after p50 %+.2f (p5 %+.2f, p95 %+.2f; >1.5 vox %.1f%%); moved p50 %+.2f, |max| %.2f, capped "
                    "%zu; slope p99 %.3f; %d iters (pcg %d its, resid %.1e); profile before %.0f/%.0f/%.0f after "
                    "%.0f/%.0f/%.0f at -4/0/+4 (%.1fs)\n",
            out, st.n_valid, st.windows, st.windows == 1 ? "" : "s", st.prm_search, st.prm_max_total, st.prm_min_contrast,
            st.prm_sigma, use_prior ? (prior_path != NULL ? "prior record" : "prior") : "constants",
            100.0 * (double)st.n_supported_first / (double)st.n_valid,
            100.0 * (double)st.n_in_band_first / (double)st.n_valid, 100.0 * (double)st.n_in_gap_first / (double)st.n_valid,
            100.0 * (double)st.n_beyond_first / (double)st.n_valid,
            st.target_first.p50, st.target_first.p5, st.target_first.p95, 100.0 * st.target_first.share_over_1p5,
            st.target_last.p50, st.target_last.p5, st.target_last.p95, 100.0 * st.target_last.share_over_1p5,
            o.smooth ? 2 * TXR_MEDIAN_HALF_SMOOTH + 1 : 2 * TXR_MEDIAN_HALF + 1, st.target_first_med.p50,
            st.target_first_med.p5, st.target_first_med.p95,
            100.0 * st.target_first_med.share_over_1p5, st.target_last_med.p50, st.target_last_med.p5,
            st.target_last_med.p95, 100.0 * st.target_last_med.share_over_1p5,
            st.disp.p50, st.disp.max_abs, st.n_capped_total, st.slope_p99, st.outer_iters,
            st.outer_iters > 0 ? st.pcg_iters[st.outer_iters - 1] : 0,
            st.outer_iters > 0 ? st.pcg_resid[st.outer_iters - 1] : 0.0,
            st.profile_first[TXR_PROFILE_HALF - 4], st.profile_first[TXR_PROFILE_HALF], st.profile_first[TXR_PROFILE_HALF + 4],
            st.profile_last[TXR_PROFILE_HALF - 4], st.profile_last[TXR_PROFILE_HALF], st.profile_last[TXR_PROFILE_HALF + 4],
            ves_clock_sec() - t0);
    Arena_dispose(&arena);
    return 0;
}
