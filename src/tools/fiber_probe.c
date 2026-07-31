/* ============================================================================
 * fiber_probe.c -- Phase-A verification spike for the ribbon-relax stage.
 *
 * Reads a baked rawtex texture (the (u,v) raster of RAW CT intensity over the
 * unrolled sheet) and visualizes the papyrus fiber orientation field recovered
 * by the structure tensor (src/flatten/fiber_field.c), so we can confirm --
 * before building the optimizer -- that classical (no-ML) fiber detection
 * actually captures the cross-hatch, how far off-axis fibers currently sit, and
 * where coherence collapses (cracks/holes) so alignment must be gated off.
 *
 * Default mode is 4-RoSy (--rosy 4): the two orthogonal cross-hatch families
 * reinforce into one grid-orientation field (0 = axes aligned), which is what
 * the aligner consumes. --rosy 2 gives the plain line-orientation field.
 *
 * Outputs (given -o PREFIX), full-width and stacked-strip PNGs:
 *   PREFIX_fiber_orient.png[/_strip]  HSV overlay: hue=orientation, sat=coherence,
 *                                     val=texture intensity (gray where coh<coh_min).
 *   PREFIX_fiber_coh.png / _strip     grayscale coherence (bright=strongly oriented).
 *   PREFIX_fiber_quiver.png[/_strip]  strokes along the fiber direction(s) over a
 *                                     dim texture, green=on-axis / red=off-axis.
 *   PREFIX_fiber_stats.json           histogram + peaks + coverage/axis stats.
 *
 * --passthrough just re-emits a uint8 TIF as PNG + strip (to view the existing
 * _diagstretch/_diagsquash distortion diagnostics, which are written TIF-only).
 * ==========================================================================*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/arena.h"
#include "../common/tiff_io.h"
#include "../common/ves_png.h"
#include "../common/tif_strip.h"
#include "../flatten/fiber_field.h"

#define PROBE_PI 3.14159265358979323846

/* ---- small raster helpers ------------------------------------------------ */

static void hsv2rgb(double h, double s, double v, uint8_t *r, uint8_t *g, uint8_t *b)
{
    double c = v * s, hp = h / 60.0, x = c * (1.0 - fabs(fmod(hp, 2.0) - 1.0));
    double m = v - c, rr = 0, gg = 0, bb = 0;
    if      (hp < 1) { rr = c; gg = x; }
    else if (hp < 2) { rr = x; gg = c; }
    else if (hp < 3) { gg = c; bb = x; }
    else if (hp < 4) { gg = x; bb = c; }
    else if (hp < 5) { rr = x; bb = c; }
    else             { rr = c; bb = x; }
    *r = (uint8_t)((rr + m) * 255.0 + 0.5);
    *g = (uint8_t)((gg + m) * 255.0 + 0.5);
    *b = (uint8_t)((bb + m) * 255.0 + 0.5);
}

static void draw_line_rgb(uint8_t *rgb, int W, int H, int x0, int y0, int x1, int y1,
                          uint8_t r, uint8_t g, uint8_t b)
{
    int dx = abs(x1 - x0), sx = x0 < x1 ? 1 : -1;
    int dy = -abs(y1 - y0), sy = y0 < y1 ? 1 : -1;
    int err = dx + dy, e2;
    for (;;) {
        if (x0 >= 0 && x0 < W && y0 >= 0 && y0 < H) {
            size_t o = ((size_t)y0 * W + x0) * 3;
            rgb[o] = r; rgb[o + 1] = g; rgb[o + 2] = b;
        }
        if (x0 == x1 && y0 == y1) break;
        e2 = 2 * err;
        if (e2 >= dy) { err += dy; x0 += sx; }
        if (e2 <= dx) { err += dx; y0 += sy; }
    }
}

/* Stack a wide RGB image into strips (mirror of tif_strip's grayscale stacker)
 * and write it as a PNG. Returns 0 on success. */
static int rgb_strip_png(const char *path, const uint8_t *rgb, int W, int H,
                         int target_w, int gap)
{
    int nstrips, sw, oh, s, x, y, rc;
    int64_t oh64, npix;
    uint8_t *out = NULL;
    if (target_w <= 0) target_w = 3000;
    if (gap < 0) gap = 6;
    nstrips = (W + target_w - 1) / target_w; if (nstrips < 1) nstrips = 1;
    sw = (W + nstrips - 1) / nstrips;
    oh64 = (int64_t)nstrips * ((int64_t)H + gap) - gap; if (oh64 < H) oh64 = H;
    npix = (int64_t)sw * oh64 * 3;
    if (oh64 > 100000000 || npix > ((int64_t)1 << 31)) return -1;
    oh = (int)oh64;
    out = (uint8_t *)calloc((size_t)sw * oh * 3, 1);
    if (out == NULL) return -1;
    for (s = 0; s < nstrips; s++) {
        int x0 = s * sw, oy0 = s * (H + gap);
        for (y = 0; y < H; y++)
            for (x = 0; x < sw; x++) {
                int sx = x0 + x;
                if (sx >= W) continue;
                memcpy(out + ((size_t)(oy0 + y) * sw + x) * 3,
                       rgb + ((size_t)y * W + sx) * 3, 3);
            }
    }
    rc = VesPng_write_rgb(path, out, sw, oh);
    free(out);
    return rc;
}

/* ---- output builders ----------------------------------------------------- */

static void write_orient(const char *prefix, const uint8_t *img, const FiberField *f,
                         double coh_min, int strip_w)
{
    int W = f->W, H = f->H;
    size_t N = (size_t)W * H, p;
    double range = f->range_deg * PROBE_PI / 180.0;   /* theta domain in radians */
    uint8_t *rgb = (uint8_t *)malloc(N * 3);
    char path[2600];
    if (rgb == NULL) return;
    for (p = 0; p < N; p++) {
        uint8_t r, g, b;
        if (f->valid[p] && f->coh[p] >= coh_min) {
            double hue = (double)f->theta[p] / range * 360.0;         /* [0,360) */
            hsv2rgb(hue, fmin(1.0, (double)f->coh[p] * 1.5), (double)img[p] / 255.0, &r, &g, &b);
        } else { r = g = b = (uint8_t)((double)img[p] * 0.5); }        /* dim gray */
        rgb[p * 3] = r; rgb[p * 3 + 1] = g; rgb[p * 3 + 2] = b;
    }
    snprintf(path, sizeof path, "%s_fiber_orient.png", prefix);       VesPng_write_rgb(path, rgb, W, H);
    snprintf(path, sizeof path, "%s_fiber_orient_strip.png", prefix); rgb_strip_png(path, rgb, W, H, strip_w, -1);
    free(rgb);
}

static void write_coh(const char *prefix, const FiberField *f, int strip_w)
{
    int W = f->W, H = f->H;
    size_t N = (size_t)W * H, p;
    uint8_t *g = (uint8_t *)malloc(N);
    char path[2600];
    if (g == NULL) return;
    for (p = 0; p < N; p++) g[p] = (uint8_t)((double)f->coh[p] * 255.0 + 0.5);
    snprintf(path, sizeof path, "%s_fiber_coh.png", prefix);        VesPng_write_gray(path, g, W, H);
    snprintf(path, sizeof path, "%s_fiber_coh_strip.tif", prefix);  TifStrip_write_gray(path, g, W, H, strip_w, -1);
    free(g);
}

static void stroke(uint8_t *rgb, int W, int H, int cx, int cy, double ang, double hl,
                   double coh, double axis_tol)
{
    double deg = ang * 180.0 / PROBE_PI, d90 = fmod(fmod(deg, 90.0) + 90.0, 90.0);
    double err = fmin(d90, 90.0 - d90), t = err / 45.0;
    double br = 0.4 + 0.6 * coh;
    double R = (err <= axis_tol) ? 40.0 : (60.0 + 195.0 * (t > 1 ? 1 : t));
    double G = (err <= axis_tol) ? 235.0 : (200.0 * (1.0 - (t > 1 ? 1 : t)));
    uint8_t r = (uint8_t)(R * br), g = (uint8_t)(G * br), b = (uint8_t)(40.0 * br);
    draw_line_rgb(rgb, W, H, (int)(cx - hl * cos(ang) + 0.5), (int)(cy - hl * sin(ang) + 0.5),
                  (int)(cx + hl * cos(ang) + 0.5), (int)(cy + hl * sin(ang) + 0.5), r, g, b);
}

static void write_quiver(const char *prefix, const uint8_t *img, const FiberField *f,
                         int step, double axis_tol, double coh_min, int strip_w)
{
    int W = f->W, H = f->H, cx, cy;
    size_t N = (size_t)W * H, p;
    uint8_t *rgb = (uint8_t *)malloc(N * 3);
    char path[2600];
    double hl = step * 0.45;
    if (rgb == NULL) return;
    for (p = 0; p < N; p++) { uint8_t v = (uint8_t)((double)img[p] * 0.35); rgb[p*3]=v; rgb[p*3+1]=v; rgb[p*3+2]=v; }
    for (cy = step / 2; cy < H; cy += step)
        for (cx = step / 2; cx < W; cx += step) {
            size_t c = (size_t)cy * W + cx;
            double th, coh;
            if (!f->valid[c] || f->coh[c] < coh_min) continue;
            th = (double)f->theta[c]; coh = (double)f->coh[c];
            stroke(rgb, W, H, cx, cy, th, hl, coh, axis_tol);
            if (f->rosy == 4) stroke(rgb, W, H, cx, cy, th + PROBE_PI / 2.0, hl, coh, axis_tol);  /* other family */
        }
    snprintf(path, sizeof path, "%s_fiber_quiver.png", prefix);       VesPng_write_rgb(path, rgb, W, H);
    snprintf(path, sizeof path, "%s_fiber_quiver_strip.png", prefix); rgb_strip_png(path, rgb, W, H, strip_w, -1);
    free(rgb);
}

/* Aggregate the field over cell x cell blocks (coherence-weighted doubled-angle
 * mean -- the same reduction the relax stage does per mesh face over the texels
 * a face covers) and render the coarse field: blocky hue overlay + one stroke
 * per cell. This previews what the optimizer actually consumes. */
static void write_cells(const char *prefix, const uint8_t *img, const FiberField *f,
                        int cell, double axis_tol, int strip_w)
{
    int W = f->W, H = f->H, ncx = (W + cell - 1) / cell, ncy = (H + cell - 1) / cell;
    int ci, cj, x, y;
    double range = f->range_deg * PROBE_PI / 180.0, k = (double)f->rosy, gate = 0.20;
    double aerr_sum = 0.0, aw = 0.0;
    size_t nc = (size_t)ncx * ncy, ncov = 0;
    double *ctheta = (double *)malloc(nc * sizeof *ctheta);
    double *ccoh   = (double *)malloc(nc * sizeof *ccoh);
    uint8_t *cvalid = (uint8_t *)calloc(nc, 1);
    uint8_t *rgb = (uint8_t *)malloc((size_t)W * H * 3);
    char path[2600];
    if (ctheta == NULL || ccoh == NULL || cvalid == NULL || rgb == NULL) {
        free(ctheta); free(ccoh); free(cvalid); free(rgb); return;
    }
    for (cj = 0; cj < ncy; cj++)
        for (ci = 0; ci < ncx; ci++) {
            double sx = 0.0, sy = 0.0, wc = 0.0;
            int nv = 0, ntot = 0, x0 = ci * cell, y0 = cj * cell;
            for (y = y0; y < y0 + cell && y < H; y++)
                for (x = x0; x < x0 + cell && x < W; x++) {
                    size_t p = (size_t)y * W + x;
                    ntot++;
                    if (!f->valid[p]) continue;
                    { double c = f->coh[p], th = f->theta[p];
                      sx += c * cos(k * th); sy += c * sin(k * th); wc += c; nv++; }
                }
            if (nv > 0 && nv >= ntot / 4 && wc > 1e-9) {
                double ang = atan2(sy, sx) / k;
                ang = fmod(ang, range); if (ang < 0) ang += range;
                ctheta[(size_t)cj * ncx + ci] = ang;
                ccoh[(size_t)cj * ncx + ci] = sqrt(sx * sx + sy * sy) / wc;   /* within-cell consistency */
                cvalid[(size_t)cj * ncx + ci] = 1;
            }
        }
    /* blocky hue overlay */
    for (y = 0; y < H; y++)
        for (x = 0; x < W; x++) {
            size_t p = (size_t)y * W + x, cidx = (size_t)(y / cell) * ncx + (x / cell);
            uint8_t r, g, b;
            if (cvalid[cidx] && ccoh[cidx] >= gate) {
                double hue = ctheta[cidx] / range * 360.0;
                hsv2rgb(hue, fmin(1.0, ccoh[cidx] * 1.3), (double)img[p] / 255.0, &r, &g, &b);
            } else { r = g = b = (uint8_t)((double)img[p] * 0.5); }
            rgb[p * 3] = r; rgb[p * 3 + 1] = g; rgb[p * 3 + 2] = b;
        }
    snprintf(path, sizeof path, "%s_fiber_cellorient.png", prefix);       VesPng_write_rgb(path, rgb, W, H);
    snprintf(path, sizeof path, "%s_fiber_cellorient_strip.png", prefix); rgb_strip_png(path, rgb, W, H, strip_w, -1);
    /* one stroke per cell over dim texture */
    for (y = 0; y < H; y++)
        for (x = 0; x < W; x++) { size_t p=(size_t)y*W+x; uint8_t v=(uint8_t)((double)img[p]*0.35); rgb[p*3]=v; rgb[p*3+1]=v; rgb[p*3+2]=v; }
    for (cj = 0; cj < ncy; cj++)
        for (ci = 0; ci < ncx; ci++) {
            size_t cidx = (size_t)cj * ncx + ci;
            double th, co, d90, err;
            int cx, cy;
            if (!cvalid[cidx] || ccoh[cidx] < gate) continue;
            th = ctheta[cidx]; co = ccoh[cidx];
            cx = ci * cell + cell / 2; cy = cj * cell + cell / 2;
            stroke(rgb, W, H, cx, cy, th, cell * 0.45, co, axis_tol);
            if (f->rosy == 4) stroke(rgb, W, H, cx, cy, th + PROBE_PI / 2.0, cell * 0.45, co, axis_tol);
            d90 = fmod(th * 180.0 / PROBE_PI, 90.0); err = fmin(d90, 90.0 - d90);
            aerr_sum += co * err; aw += co; ncov++;
        }
    snprintf(path, sizeof path, "%s_fiber_cellquiver.png", prefix);       VesPng_write_rgb(path, rgb, W, H);
    snprintf(path, sizeof path, "%s_fiber_cellquiver_strip.png", prefix); rgb_strip_png(path, rgb, W, H, strip_w, -1);
    fprintf(stderr, "  [fiber] cells=%dpx: %zu/%zu covered, mean_cell_axis_err=%.2fdeg (coh-wtd)\n",
            cell, ncov, nc, aw > 0 ? aerr_sum / aw : 0.0);
    free(ctheta); free(ccoh); free(cvalid); free(rgb);
}

static void report_stats(const char *prefix, const FiberField *f)
{
    char path[2600];
    double binw = f->range_deg / FIBER_NBINS;
    int b, peak1 = 0, peak2 = -1, sep = (int)(30.0 / binw + 0.5);
    FILE *fp;
    for (b = 1; b < FIBER_NBINS; b++) if (f->hist[b] > f->hist[peak1]) peak1 = b;
    for (b = 0; b < FIBER_NBINS; b++) {
        if (abs(b - peak1) < sep) continue;
        if (peak2 < 0 || f->hist[b] > f->hist[peak2]) peak2 = b;
    }
    fprintf(stderr, "  [fiber] rosy=%d range=%.0fdeg  valid=%.1f%%  mean_coh=%.3f  peak=%.1fdeg  "
            "axis_frac=%.1f%%  mean_axis_err=%.2fdeg\n",
            f->rosy, f->range_deg, 100.0 * f->frac_valid, f->mean_coh, f->peak_deg,
            100.0 * f->axis_frac, f->mean_axis_err_deg);
    fprintf(stderr, "  [fiber] orientation peaks: %.1fdeg (%.3f)  %.1fdeg (%.3f)\n",
            (peak1 + 0.5) * binw, f->hist[peak1],
            peak2 >= 0 ? (peak2 + 0.5) * binw : -1.0, peak2 >= 0 ? f->hist[peak2] : 0.0);
    fprintf(stderr, "  [fiber] histogram (deg: weight), 0=aligned:\n    ");
    for (b = 0; b < FIBER_NBINS; b++) {
        int nb = (int)(f->hist[b] * 200.0 + 0.5), k;
        fprintf(stderr, "%3d:", (int)(b * binw));
        for (k = 0; k < nb && k < 20; k++) fputc('#', stderr);
        fputc(b % 6 == 5 ? '\n' : ' ', stderr);
        if (b % 6 == 5) fprintf(stderr, "    ");
    }
    fprintf(stderr, "\n");

    snprintf(path, sizeof path, "%s_fiber_stats.json", prefix);
    fp = fopen(path, "wb");
    if (fp != NULL) {
        fprintf(fp, "{\n  \"W\": %d, \"H\": %d, \"rosy\": %d, \"range_deg\": %.1f,\n", f->W, f->H, f->rosy, f->range_deg);
        fprintf(fp, "  \"frac_valid\": %.4f, \"mean_coh\": %.4f,\n", f->frac_valid, f->mean_coh);
        fprintf(fp, "  \"peak_deg\": %.2f, \"axis_frac\": %.4f, \"mean_axis_err_deg\": %.3f,\n",
                f->peak_deg, f->axis_frac, f->mean_axis_err_deg);
        fprintf(fp, "  \"hist\": [");
        for (b = 0; b < FIBER_NBINS; b++) fprintf(fp, "%s%.5f", b ? ", " : "", f->hist[b]);
        fprintf(fp, "]\n}\n");
        fclose(fp);
        fprintf(stderr, "  [fiber] wrote %s\n", path);
    }
}

/* ---- passthrough: uint8 TIF -> PNG + strip ------------------------------- */

static int do_passthrough(Arena_T arena, const char *in, const char *prefix)
{
    uint8_t *vol = NULL; int D = 0, H = 0, W = 0;
    char path[2600];
    if (TiffIO_load(arena, in, &vol, &D, &H, &W) != 0) { fprintf(stderr, "load fail: %s\n", in); return -1; }
    snprintf(path, sizeof path, "%s.png", prefix);       VesPng_write_gray(path, vol, W, H);
    snprintf(path, sizeof path, "%s_strip.tif", prefix); TifStrip_write_gray(path, vol, W, H, 0, -1);
    fprintf(stderr, "  [passthrough] %s (%dx%d) -> %s.png (+_strip)\n", in, W, H, prefix);
    return 0;
}

static void usage(void)
{
    fprintf(stderr,
        "usage: fiber_probe --selftest\n"
        "       fiber_probe <in_rawtex.tif> -o <out_prefix> [options]\n"
        "       fiber_probe --passthrough <in.tif> -o <out_prefix>\n"
        "options: --rosy 2|4 (4)  --grad-sigma f (1.5)  --tensor-sigma f (6.0)\n"
        "         --axis-tol f (15)  --coh-min f (0.05)  --quiver-step n (16)\n"
        "         --cell n (0=off; per-cell aggregate preview)  --strip-w n (3000)\n");
}

int main(int argc, char **argv)
{
    const char *in = NULL, *prefix = NULL;
    double grad_sigma = 1.5, tensor_sigma = 6.0, axis_tol = 15.0, coh_min = 0.05;
    int rosy = 4, quiver_step = 16, cell = 0, strip_w = 3000, passthrough = 0, i;
    Arena_T arena = NULL;
    uint8_t *vol = NULL; int D = 0, H = 0, W = 0;
    FiberField f;

    if (argc >= 2 && strcmp(argv[1], "--selftest") == 0) {
        int fails = FiberField_selftest() + TifStrip_selftest();
        fprintf(stderr, "[fiber_probe selftest] %s\n", fails == 0 ? "ALL OK" : "FAILED");
        return fails == 0 ? 0 : 1;
    }

    for (i = 1; i < argc; i++) {
        if      (strcmp(argv[i], "-o") == 0 && i + 1 < argc)             prefix = argv[++i];
        else if (strcmp(argv[i], "--passthrough") == 0)                  passthrough = 1;
        else if (strcmp(argv[i], "--rosy") == 0 && i + 1 < argc)         rosy = atoi(argv[++i]);
        else if (strcmp(argv[i], "--grad-sigma") == 0 && i + 1 < argc)   grad_sigma = atof(argv[++i]);
        else if (strcmp(argv[i], "--tensor-sigma") == 0 && i + 1 < argc) tensor_sigma = atof(argv[++i]);
        else if (strcmp(argv[i], "--axis-tol") == 0 && i + 1 < argc)     axis_tol = atof(argv[++i]);
        else if (strcmp(argv[i], "--coh-min") == 0 && i + 1 < argc)      coh_min = atof(argv[++i]);
        else if (strcmp(argv[i], "--quiver-step") == 0 && i + 1 < argc)  quiver_step = atoi(argv[++i]);
        else if (strcmp(argv[i], "--cell") == 0 && i + 1 < argc)         cell = atoi(argv[++i]);
        else if (strcmp(argv[i], "--strip-w") == 0 && i + 1 < argc)      strip_w = atoi(argv[++i]);
        else if (argv[i][0] != '-')                                      in = argv[i];
        else { fprintf(stderr, "unknown arg: %s\n", argv[i]); usage(); return 2; }
    }
    if (in == NULL || prefix == NULL) { usage(); return 2; }

    arena = Arena_new();
    if (passthrough) { int rc = do_passthrough(arena, in, prefix); Arena_dispose(&arena); return rc == 0 ? 0 : 1; }

    if (TiffIO_load(arena, in, &vol, &D, &H, &W) != 0) {
        fprintf(stderr, "load fail: %s\n", in); Arena_dispose(&arena); return 1;
    }
    if (D != 1) fprintf(stderr, "  [fiber] note: %d pages, using page 0\n", D);
    fprintf(stderr, "[fiber_probe] %s (%dx%d)  rosy=%d grad_sigma=%.2f tensor_sigma=%.2f coh_min=%.2f\n",
            in, W, H, rosy, grad_sigma, tensor_sigma, coh_min);

    if (FiberField_compute(arena, vol, W, H, grad_sigma, tensor_sigma, axis_tol, rosy, &f) != 0) {
        fprintf(stderr, "fiber compute fail\n"); Arena_dispose(&arena); return 1;
    }
    write_orient(prefix, vol, &f, coh_min, strip_w);
    write_coh(prefix, &f, strip_w);
    write_quiver(prefix, vol, &f, quiver_step, axis_tol, coh_min, strip_w);
    if (cell > 1) write_cells(prefix, vol, &f, cell, axis_tol, strip_w);
    report_stats(prefix, &f);

    Arena_dispose(&arena);
    return 0;
}
