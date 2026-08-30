/* marble_diag.c -- see marble_diag.h. */

#include "marble_diag.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "rawtex_bake.h"
#include "../common/arena.h"
#include "../common/tiff_io.h"
#include "../common/ves_png.h"

typedef struct MarbleSample {
    float p[3];
    float tu[3];
    float tv[3];
    float n[3];
    uint8_t owner;             /* 0 absent, 1 boundary, 2 strict interior */
} MarbleSample;

typedef struct MarbleProbe {
    float offset;
    float center_support;
    float best_support;
    float axis_bad;
    uint8_t valid;
} MarbleProbe;

static void *md_calloc(size_t count, size_t size)
{
    if (count > 0 && size > SIZE_MAX / count) return NULL;
    return calloc(count > 0 ? count : 1, size);
}

/* PGM is intentionally emitted beside the human-facing PNG: it is a tiny,
 * dependency-free contract between diagnosis and the geometry repair stage.
 * The repairer must never infer its weights from a coloured visualization. */
static int write_pgm(const char *path, const uint8_t *pixels, size_t W, size_t H)
{
    FILE *f;
    if (path == NULL || pixels == NULL || W == 0 || H == 0) return -1;
    f = fopen(path, "wb");
    if (f == NULL) return -1;
    if (fprintf(f, "P5\n%zu %zu\n255\n", W, H) < 0 ||
        fwrite(pixels, 1, W * H, f) != W * H || fclose(f) != 0)
        return -1;
    return 0;
}

static double clamp01(double x)
{
    if (x < 0.0) return 0.0;
    if (x > 1.0) return 1.0;
    return x;
}

void MarbleDiag_defaults(MarbleDiagOpts *opts)
{
    if (opts == NULL) return;
    memset(opts, 0, sizeof *opts);
    opts->stride_pixels = 16;
    opts->depth_range = 4.0;
    opts->depth_step = 1.0;
    opts->tensor_radius = 2.0;
    opts->score_threshold = 0.40;
    opts->minimum_region_cells = 12;
}

/* Affine UV->XYZ frame for one triangle. UV has already been converted to
 * baked-raster coordinates, which changes only the lengths, not directions. */
static int face_frame(const float *verts, size_t a, size_t b, size_t c,
                      double ua, double va, double ub, double vb,
                      double uc, double vc, const float *vertex_normals,
                      double l0, double l1, double l2,
                      MarbleSample *sample)
{
    double A2 = (ub - ua) * (vc - va) - (vb - va) * (uc - ua);
    double nn = 0.0, align = 0.0;
    double vn[3] = {0.0, 0.0, 0.0};
    int k;
    if (fabs(A2) < 1e-12) return -1;
    for (k = 0; k < 3; k++) {
        double pa = (double)verts[a * 3 + (size_t)k];
        double pb = (double)verts[b * 3 + (size_t)k];
        double pc = (double)verts[c * 3 + (size_t)k];
        sample->p[k] = (float)(l0 * pa + l1 * pb + l2 * pc);
        sample->tu[k] = (float)((pa * (vb - vc) + pb * (vc - va)
                                 + pc * (va - vb)) / A2);
        sample->tv[k] = (float)((pa * (uc - ub) + pb * (ua - uc)
                                 + pc * (ub - ua)) / A2);
        if (vertex_normals != NULL)
            vn[k] = l0 * (double)vertex_normals[a * 3 + (size_t)k]
                  + l1 * (double)vertex_normals[b * 3 + (size_t)k]
                  + l2 * (double)vertex_normals[c * 3 + (size_t)k];
    }
    sample->n[0] = sample->tu[1] * sample->tv[2]
                 - sample->tu[2] * sample->tv[1];
    sample->n[1] = sample->tu[2] * sample->tv[0]
                 - sample->tu[0] * sample->tv[2];
    sample->n[2] = sample->tu[0] * sample->tv[1]
                 - sample->tu[1] * sample->tv[0];
    for (k = 0; k < 3; k++) nn += (double)sample->n[k] * sample->n[k];
    nn = sqrt(nn);
    if (nn < 1e-9) return -1;
    for (k = 0; k < 3; k++) sample->n[k] = (float)((double)sample->n[k] / nn);
    if (vertex_normals != NULL) {
        for (k = 0; k < 3; k++) align += (double)sample->n[k] * vn[k];
        if (align < 0.0)
            for (k = 0; k < 3; k++) sample->n[k] = -sample->n[k];
    }
    return 0;
}

static int probe_alignment(CubeTable *ct, const MarbleSample *sample,
                           const MarbleDiagOpts *opts, MarbleProbe *probe,
                           float *profile, int half)
{
    int di, have_center = 0;
    RawTangentTensor center;
    memset(probe, 0, sizeof *probe);
    memset(&center, 0, sizeof center);
    for (di = -half; di <= half; di++) {
        double t = (double)di * opts->depth_step;
        double q[3], tu[3], tv[3];
        RawTangentTensor tensor;
        int k;
        if (fabs(t) > opts->depth_range + 1e-9) continue;
        for (k = 0; k < 3; k++) {
            q[k] = (double)sample->p[k] + t * (double)sample->n[k];
            tu[k] = sample->tu[k];
            tv[k] = sample->tv[k];
        }
        if (sample_tangent_tensor(ct, q, tu, tv, opts->tensor_radius,
                                  &tensor) != 0)
            continue;
        {
            double strength = tensor.rms_gradient /
                              (tensor.rms_gradient + 8.0);
            /* Locate the RAW sheet independently of fiber direction.  Axis
             * agreement is meaningful only after a sheet displacement has
             * been established; using it here made arbitrary diagonal fibers
             * win the depth search. */
            profile[di + half] = (float)clamp01(
                tensor.ridge_center * strength *
                (0.25 + 0.75 * tensor.coherence));
        }
        if (di == 0) {
            center = tensor;
            have_center = 1;
        }
    }
    if (!have_center) return -1;
    {
        double strength = center.rms_gradient / (center.rms_gradient + 8.0);
        probe->axis_bad = (float)(center.ridge_center * strength
                                  * center.coherence
                                  * (1.0 - center.axis_alignment));
        probe->valid = 1;
    }
    return 0;
}

static double percentile_from_hist(const size_t hist[256], size_t total,
                                   double percentile)
{
    size_t target, sum = 0;
    int i;
    if (total == 0) return 0.0;
    target = (size_t)ceil(percentile * 0.01 * (double)total);
    if (target < 1) target = 1;
    for (i = 0; i < 256; i++) {
        sum += hist[i];
        if (sum >= target) return (double)i / 255.0;
    }
    return 1.0;
}

static void write_gray_map(const char *prefix, const char *suffix,
                           const float *value, const uint8_t *valid,
                           size_t n, int W, int H)
{
    char path[2600];
    uint8_t *img = (uint8_t *)md_calloc(n, 1);
    size_t i;
    if (img == NULL) return;
    for (i = 0; i < n; i++)
        if (valid[i]) img[i] = (uint8_t)(255.0 * clamp01(value[i]) + 0.5);
    snprintf(path, sizeof path, "%s_%s.png", prefix, suffix);
    VesPng_write_gray(path, img, W, H);
    free(img);
}

static void write_offset_map(const char *prefix, const MarbleProbe *probe,
                             size_t n, int W, int H, double depth_range)
{
    char path[2600];
    uint8_t *rgb = (uint8_t *)md_calloc(n * 3, 1);
    size_t i;
    if (rgb == NULL) return;
    for (i = 0; i < n; i++) {
        double t, conf, mag, base;
        if (!probe[i].valid) continue;
        t = depth_range > 0.0 ? (double)probe[i].offset / depth_range : 0.0;
        if (t < -1.0) t = -1.0; else if (t > 1.0) t = 1.0;
        conf = clamp01((double)probe[i].best_support / 0.20);
        mag = fabs(t);
        base = 32.0 + 96.0 * conf * (1.0 - mag);
        rgb[i * 3 + 0] = (uint8_t)(255.0 * clamp01(
            (base + (t > 0.0 ? 127.0 * conf * mag : 0.0)) / 255.0) + 0.5);
        rgb[i * 3 + 1] = (uint8_t)base;
        rgb[i * 3 + 2] = (uint8_t)(255.0 * clamp01(
            (base + (t < 0.0 ? 127.0 * conf * mag : 0.0)) / 255.0) + 0.5);
    }
    snprintf(path, sizeof path, "%s_diagmarble_offset.png", prefix);
    VesPng_write_rgb(path, rgb, W, H);
    free(rgb);
}

static void remove_small_regions(uint8_t *mask, size_t W, size_t H,
                                 int minimum_cells, size_t *out_regions,
                                 size_t *out_cells)
{
    size_t n = W * H, i, regions = 0, kept = 0;
    uint8_t *seen = (uint8_t *)md_calloc(n, 1);
    size_t *queue = (size_t *)md_calloc(n, sizeof *queue);
    if (seen == NULL || queue == NULL) {
        free(seen); free(queue); return;
    }
    for (i = 0; i < n; i++) {
        size_t head = 0, tail = 0;
        if (!mask[i] || seen[i]) continue;
        seen[i] = 1; queue[tail++] = i;
        while (head < tail) {
            size_t p = queue[head++], x = p % W, y = p / W;
#define MD_PUSH(Q) do { size_t _q=(Q); if(mask[_q]&&!seen[_q]){seen[_q]=1;queue[tail++]=_q;} } while(0)
            if (x > 0) MD_PUSH(p - 1);
            if (x + 1 < W) MD_PUSH(p + 1);
            if (y > 0) MD_PUSH(p - W);
            if (y + 1 < H) MD_PUSH(p + W);
#undef MD_PUSH
        }
        if (tail < (size_t)minimum_cells) {
            size_t q;
            for (q = 0; q < tail; q++) mask[queue[q]] = 0;
        } else {
            regions++; kept += tail;
        }
    }
    free(seen); free(queue);
    if (out_regions != NULL) *out_regions = regions;
    if (out_cells != NULL) *out_cells = kept;
}

int MarbleDiag_write(const char *prefix, const char *rawtex_path,
                     CubeTable *ct,
                     const float *verts, const float *uv, size_t nv,
                     const int32_t *faces, size_t nf,
                     const float *vertex_normals,
                     double raster_du, double raster_dv, size_t raster_max_px,
                     const MarbleDiagOpts *input_opts, MarbleDiagStats *stats)
{
    MarbleDiagOpts opts;
    RawtexPlan plan;
    MarbleSample *sample = NULL;
    MarbleProbe *probe = NULL;
    float *profile = NULL, *smooth_profile = NULL;
    float *support_center = NULL, *support_best = NULL;
    float *gain = NULL, *axis = NULL, *axis_broad = NULL;
    float *disorder = NULL, *miss = NULL;
    float *raw_score = NULL, *score = NULL;
    uint8_t *valid = NULL, *mask = NULL, *full_mask = NULL, *rgb = NULL;
    size_t W, H, MW, MH, mn, f, i, sampled = 0, highlighted = 0;
    size_t nd;
    size_t hist[256];
    Arena_T image_arena = NULL;
    uint8_t *texture = NULL;
    int depth = 0, image_h = 0, image_w = 0;
    int half;
    int rc = -1;
    char path[2600];
    MarbleDiagStats result;

    if (stats != NULL) memset(stats, 0, sizeof *stats);
    memset(&result, 0, sizeof result);
    if (input_opts != NULL) opts = *input_opts; else MarbleDiag_defaults(&opts);
    if (prefix == NULL || rawtex_path == NULL || ct == NULL || verts == NULL ||
        uv == NULL || faces == NULL || nv == 0 || nf == 0 ||
        opts.stride_pixels < 2 || opts.depth_range <= 0.0 ||
        opts.depth_step <= 0.0 || opts.tensor_radius <= 0.0 ||
        opts.score_threshold < 0.0 || opts.score_threshold > 1.0 ||
        opts.minimum_region_cells < 1)
        return -1;
    if (Rawtex_plan(uv, nv, raster_du, raster_dv, raster_max_px, &plan) != 0 ||
        !plan.ok)
        return -1;
    W = plan.W; H = plan.H;
    MW = (W + (size_t)opts.stride_pixels - 1) / (size_t)opts.stride_pixels;
    MH = (H + (size_t)opts.stride_pixels - 1) / (size_t)opts.stride_pixels;
    mn = MW * MH;
    half = (int)floor(opts.depth_range / opts.depth_step + 1e-9);
    if (half < 1 || half > (INT32_MAX - 1) / 2) return -1;
    nd = (size_t)(2 * half + 1);
    if (mn == 0 || mn > (size_t)INT32_MAX || mn > SIZE_MAX / nd) return -1;
    sample = (MarbleSample *)md_calloc(mn, sizeof *sample);
    probe = (MarbleProbe *)md_calloc(mn, sizeof *probe);
    profile = (float *)md_calloc(mn * nd, sizeof *profile);
    smooth_profile = (float *)md_calloc(mn * nd, sizeof *smooth_profile);
    support_center = (float *)md_calloc(mn, sizeof *support_center);
    support_best = (float *)md_calloc(mn, sizeof *support_best);
    gain = (float *)md_calloc(mn, sizeof *gain);
    axis = (float *)md_calloc(mn, sizeof *axis);
    axis_broad = (float *)md_calloc(mn, sizeof *axis_broad);
    disorder = (float *)md_calloc(mn, sizeof *disorder);
    miss = (float *)md_calloc(mn, sizeof *miss);
    raw_score = (float *)md_calloc(mn, sizeof *raw_score);
    score = (float *)md_calloc(mn, sizeof *score);
    valid = (uint8_t *)md_calloc(mn, 1);
    mask = (uint8_t *)md_calloc(mn, 1);
    if (!sample || !probe || !profile || !smooth_profile || !support_center ||
        !support_best || !gain || !axis || !axis_broad || !disorder || !miss ||
        !raw_score || !score || !valid || !mask)
        goto cleanup;
    for (i = 0; i < mn * nd; i++) {
        profile[i] = -1.0f;
        smooth_profile[i] = -1.0f;
    }

    /* Sparse UV rasterization: one barycentric source point and affine frame
     * per audit cell. A strict triangle interior supersedes an edge owner. */
    for (f = 0; f < nf; f++) {
        size_t a = (size_t)faces[f * 3], b = (size_t)faces[f * 3 + 1];
        size_t c = (size_t)faces[f * 3 + 2];
        double ua, va, ub, vb, uc, vc, A2, lox, hix, loy, hiy;
        long mx0, mx1, my0, my1, mx, my;
        if (a >= nv || b >= nv || c >= nv) continue;
        ua = ((double)uv[a * 2] - plan.umin) / raster_du;
        va = ((double)uv[a * 2 + 1] - plan.vmin) / raster_dv;
        ub = ((double)uv[b * 2] - plan.umin) / raster_du;
        vb = ((double)uv[b * 2 + 1] - plan.vmin) / raster_dv;
        uc = ((double)uv[c * 2] - plan.umin) / raster_du;
        vc = ((double)uv[c * 2 + 1] - plan.vmin) / raster_dv;
        A2 = (ub - ua) * (vc - va) - (vb - va) * (uc - ua);
        if (fabs(A2) < 1e-12) continue;
        lox = fmin(ua, fmin(ub, uc)); hix = fmax(ua, fmax(ub, uc));
        loy = fmin(va, fmin(vb, vc)); hiy = fmax(va, fmax(vb, vc));
        mx0 = (long)floor(lox / opts.stride_pixels);
        mx1 = (long)floor(hix / opts.stride_pixels);
        my0 = (long)floor(loy / opts.stride_pixels);
        my1 = (long)floor(hiy / opts.stride_pixels);
        if (mx0 < 0) mx0 = 0; if (my0 < 0) my0 = 0;
        if (mx1 >= (long)MW) mx1 = (long)MW - 1;
        if (my1 >= (long)MH) my1 = (long)MH - 1;
        for (my = my0; my <= my1; my++) for (mx = mx0; mx <= mx1; mx++) {
            size_t xstart = (size_t)mx * (size_t)opts.stride_pixels;
            size_t ystart = (size_t)my * (size_t)opts.stride_pixels;
            size_t xcount = xstart + (size_t)opts.stride_pixels <= W
                          ? (size_t)opts.stride_pixels : W - xstart;
            size_t ycount = ystart + (size_t)opts.stride_pixels <= H
                          ? (size_t)opts.stride_pixels : H - ystart;
            size_t xpixel = xstart + (xcount - 1) / 2;
            size_t ypixel = ystart + (ycount - 1) / 2;
            double pu = (double)xpixel + 0.5, pv = (double)ypixel + 0.5;
            double l1 = ((pu - ua) * (vc - va) - (pv - va) * (uc - ua)) / A2;
            double l2 = ((ub - ua) * (pv - va) - (vb - va) * (pu - ua)) / A2;
            double l0 = 1.0 - l1 - l2;
            int boundary;
            size_t mi;
            MarbleSample candidate;
            if (l0 < -1e-9 || l1 < -1e-9 || l2 < -1e-9) continue;
            boundary = l0 <= 1e-9 || l1 <= 1e-9 || l2 <= 1e-9;
            mi = (size_t)my * MW + (size_t)mx;
            if (sample[mi].owner == 2 || (boundary && sample[mi].owner != 0))
                continue;
            memset(&candidate, 0, sizeof candidate);
            if (face_frame(verts, a, b, c, ua, va, ub, vb, uc, vc,
                           vertex_normals, l0, l1, l2, &candidate) != 0)
                continue;
            candidate.owner = (uint8_t)(boundary ? 1 : 2);
            sample[mi] = candidate;
        }
    }

    {
        int kk;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 32)
#endif
        for (kk = 0; kk < (int)mn; kk++) {
            if (sample[kk].owner != 0)
                (void)probe_alignment(ct, &sample[kk], &opts, &probe[kk],
                                      profile + (size_t)kk * nd, half);
        }
    }
    for (i = 0; i < mn; i++) if (probe[i].valid) valid[i] = 1;

    /* Average each candidate-depth plane in UV before choosing a depth.  The
     * old order (argmax first, average offsets second) turned harmless local
     * CT texture into random +/- depth choices.  A real sheet miss persists
     * across neighboring UV cells and therefore survives this operation. */
    {
        int kk;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (kk = 0; kk < (int)mn; kk++) {
            size_t x = (size_t)kk % MW, y = (size_t)kk / MW;
            int di, dy, dx;
            if (!valid[kk]) continue;
            for (di = 0; di < (int)nd; di++) {
                double sum = 0.0, wsum = 0.0;
                for (dy = -2; dy <= 2; dy++) for (dx = -2; dx <= 2; dx++) {
                    long xx = (long)x + dx, yy = (long)y + dy;
                    size_t q;
                    double w;
                    float v;
                    if (xx < 0 || yy < 0 || xx >= (long)MW || yy >= (long)MH)
                        continue;
                    q = (size_t)yy * MW + (size_t)xx;
                    if (!valid[q]) continue;
                    v = profile[q * nd + (size_t)di];
                    if (v < 0.0f) continue;
                    w = (double)(3 - abs(dx)) * (double)(3 - abs(dy));
                    sum += w * (double)v;
                    wsum += w;
                }
                if (wsum > 0.0)
                    smooth_profile[(size_t)kk * nd + (size_t)di] =
                        (float)(sum / wsum);
            }
        }
    }
    for (i = 0; i < mn; i++) {
        int di, best_di = half;
        float center, best;
        if (!valid[i]) continue;
        center = smooth_profile[i * nd + (size_t)half];
        if (center < 0.0f) { valid[i] = 0; continue; }
        best = center;
        for (di = 0; di < (int)nd; di++) {
            float candidate = smooth_profile[i * nd + (size_t)di];
            if (candidate < 0.0f) continue;
            /* Equal evidence prefers the mesh surface, then the nearer depth. */
            if (candidate > best + 1e-6f ||
                (fabs((double)candidate - (double)best) <= 1e-6 &&
                 abs(di - half) < abs(best_di - half))) {
                best = candidate;
                best_di = di;
            }
        }
        probe[i].offset = (float)((best_di - half) * opts.depth_step);
        probe[i].center_support = center;
        probe[i].best_support = best;
        support_center[i] = center;
        support_best[i] = best;
        sampled++;
    }

    /* Convert the coherent depth profiles to independent, inspectable scores.
     * - gain: a nearby normal offset has more RAW sheet support than the mesh;
     * - axis: centered coherent RAW fibers disagree with the intended U/V axes;
     * - disorder: supported neighboring cells prefer changing normal offsets.
     * Axis and disorder only corroborate an actual displaced peak.  Darkness
     * and image busyness alone contribute nothing. */
    for (i = 0; i < mn; i++) {
        double best_conf, delta, relative_gain, gain_abs, gain_rel;
        double gain_score, offset_score, axis_score;
        if (!valid[i]) continue;
        best_conf = clamp01((double)probe[i].best_support / 0.20);
        delta = fmax(0.0, (double)probe[i].best_support
                          - (double)probe[i].center_support);
        relative_gain = delta / (0.05 + (double)probe[i].center_support);
        gain_abs = clamp01((delta - 0.015) / 0.085);
        gain_rel = clamp01((relative_gain - 0.10) / 0.50);
        gain_score = sqrt(gain_abs * gain_rel);
        gain[i] = (float)gain_score;
        offset_score = clamp01((fabs((double)probe[i].offset) - 0.5)
                                / fmax(0.5, opts.depth_range - 0.5));
        axis_score = clamp01((double)probe[i].axis_bad / 0.12);
        axis[i] = (float)axis_score;
        raw_score[i] = (float)(best_conf *
            pow(gain_score * offset_score, 1.0 / 3.0) *
            (0.85 + 0.15 * axis_score));
    }

    /* Require neighboring UV cells to support the displaced peak.  This is a
     * second coherence gate, independent of the profile smoothing above: an
     * isolated high score cannot light a region, while a changing but locally
     * supported offset field still exposes the edge of a marbled patch. */
    for (i = 0; i < mn; i++) {
        size_t x = i % MW, y = i / MW;
        double valid_w = 0.0, occupied_w = 0.0, signal_w = 0.0;
        double agree_w = 0.0, mean_num = 0.0, variance_num = 0.0;
        double coverage, agreement, local_signal, mean_offset;
        double coherent_score, disorder_score;
        int dy, dx;
        if (!valid[i]) continue;
        for (dy = -2; dy <= 2; dy++) for (dx = -2; dx <= 2; dx++) {
            long xx = (long)x + dx, yy = (long)y + dy;
            size_t q;
            double w, s, od;
            if (xx < 0 || yy < 0 || xx >= (long)MW || yy >= (long)MH) continue;
            q = (size_t)yy * MW + (size_t)xx;
            if (!valid[q]) continue;
            w = (double)(3 - abs(dx)) * (double)(3 - abs(dy));
            s = clamp01((double)raw_score[q]);
            valid_w += w;
            occupied_w += w * clamp01(s / 0.20);
            signal_w += w * s;
            od = fabs((double)probe[q].offset - (double)probe[i].offset);
            agree_w += w * s * clamp01(1.0 - od /
                                       fmax(2.0 * opts.depth_step, 1.0));
            mean_num += w * s * (double)probe[q].offset;
        }
        coverage = valid_w > 0.0 ? occupied_w / valid_w : 0.0;
        local_signal = valid_w > 0.0 ? signal_w / valid_w : 0.0;
        agreement = signal_w > 0.0 ? agree_w / signal_w : 0.0;
        mean_offset = signal_w > 0.0 ? mean_num / signal_w : 0.0;
        if (signal_w > 0.0) {
            for (dy = -2; dy <= 2; dy++) for (dx = -2; dx <= 2; dx++) {
                long xx = (long)x + dx, yy = (long)y + dy;
                size_t q;
                double w, s, d;
                if (xx < 0 || yy < 0 || xx >= (long)MW || yy >= (long)MH)
                    continue;
                q = (size_t)yy * MW + (size_t)xx;
                if (!valid[q]) continue;
                w = (double)(3 - abs(dx)) * (double)(3 - abs(dy));
                s = clamp01((double)raw_score[q]);
                d = (double)probe[q].offset - mean_offset;
                variance_num += w * s * d * d;
            }
        }
        coherent_score = (double)raw_score[i] * sqrt(coverage) *
                         (0.50 + 0.50 * agreement);
        disorder_score = signal_w > 0.0 ?
            clamp01(sqrt(variance_num / signal_w) /
                    fmax(1.5 * opts.depth_step, 1.0)) * local_signal : 0.0;
        disorder[i] = (float)disorder_score;
        miss[i] = (float)fmax(coherent_score, 0.70 * disorder_score);
    }

    /* Marble is not the same class as a surface miss.  A miss has weak RAW
     * support at the mesh and a better nearby normal depth (the usual black
     * crack).  Marble has RAW surface support at the mesh, but its tangent
     * structure is broadly inconsistent with the mesh's U/V frame.  Average
     * that disagreement over a 17x17 audit neighborhood so a crack edge or a
     * single diagonal fiber cannot become a marbled region. */
    for (i = 0; i < mn; i++) {
        size_t x = i % MW, y = i / MW;
        double sum = 0.0, wsum = 0.0;
        int dy, dx;
        if (!valid[i]) continue;
        for (dy = -8; dy <= 8; dy++) for (dx = -8; dx <= 8; dx++) {
            long xx = (long)x + dx, yy = (long)y + dy;
            size_t q;
            if (xx < 0 || yy < 0 || xx >= (long)MW || yy >= (long)MH)
                continue;
            q = (size_t)yy * MW + (size_t)xx;
            if (!valid[q]) continue;
            sum += (double)axis[q];
            wsum += 1.0;
        }
        axis_broad[i] = (float)(wsum > 0.0 ? sum / wsum : 0.0);
        {
            double frame_miss = clamp01(((double)axis_broad[i] - 0.32) / 0.25);
            double on_surface = clamp01((double)support_center[i] / 0.20);
            raw_score[i] = (float)(frame_miss * on_surface);
        }
    }
    for (i = 0; i < mn; i++) {
        size_t x = i % MW, y = i / MW;
        double sum = 0.0; int count = 0, dy, dx;
        if (!valid[i]) continue;
        for (dy = -1; dy <= 1; dy++) for (dx = -1; dx <= 1; dx++) {
            long xx = (long)x + dx, yy = (long)y + dy;
            size_t q;
            if (xx < 0 || yy < 0 || xx >= (long)MW || yy >= (long)MH) continue;
            q = (size_t)yy * MW + (size_t)xx;
            if (!valid[q]) continue;
            sum += raw_score[q]; count++;
        }
        score[i] = (float)clamp01(0.70 * raw_score[i]
                                  + 0.30 * (count ? sum / count : 0.0));
        if (score[i] >= opts.score_threshold) mask[i] = 1;
    }
    remove_small_regions(mask, MW, MH, opts.minimum_region_cells,
                         &result.highlighted_regions, &highlighted);

    memset(hist, 0, sizeof hist);
    for (i = 0; i < mn; i++) if (valid[i]) {
        int b = (int)(255.0 * clamp01(score[i]) + 0.5);
        hist[b]++;
    }
    result.raster_width = W; result.raster_height = H;
    result.audit_width = MW; result.audit_height = MH;
    result.sampled_cells = sampled;
    result.highlighted_cells = highlighted;
    result.score_p50 = percentile_from_hist(hist, sampled, 50.0);
    result.score_p90 = percentile_from_hist(hist, sampled, 90.0);
    result.score_p95 = percentile_from_hist(hist, sampled, 95.0);
    result.score_p99 = percentile_from_hist(hist, sampled, 99.0);

    write_offset_map(prefix, probe, mn, (int)MW, (int)MH, opts.depth_range);
    write_gray_map(prefix, "diagmarble_support_center", support_center, valid,
                   mn, (int)MW, (int)MH);
    write_gray_map(prefix, "diagmarble_support_best", support_best, valid,
                   mn, (int)MW, (int)MH);
    write_gray_map(prefix, "diagmarble_gain", gain, valid, mn, (int)MW, (int)MH);
    write_gray_map(prefix, "diagmarble_axis", axis, valid, mn, (int)MW, (int)MH);
    write_gray_map(prefix, "diagmarble_axis_broad", axis_broad, valid,
                   mn, (int)MW, (int)MH);
    write_gray_map(prefix, "diagmarble_disorder", disorder, valid, mn, (int)MW, (int)MH);
    write_gray_map(prefix, "diagmarble_surface_miss", miss, valid,
                   mn, (int)MW, (int)MH);
    write_gray_map(prefix, "diagmarble_score", score, valid, mn, (int)MW, (int)MH);
    {
        uint8_t *coarse_mask = (uint8_t *)md_calloc(mn, 1);
        if (coarse_mask == NULL) goto cleanup;
        for (i = 0; i < mn; i++) coarse_mask[i] = mask[i] ? 255 : 0;
        snprintf(path, sizeof path, "%s_diagmarble_mask.png", prefix);
        VesPng_write_gray(path, coarse_mask, (int)MW, (int)MH);
        free(coarse_mask);
    }

    image_arena = Arena_new();
    if (TiffIO_load(image_arena, rawtex_path, &texture, &depth,
                    &image_h, &image_w) != 0 || depth != 1 ||
        image_w != (int)W || image_h != (int)H) {
        fprintf(stderr, "  [marble] cannot load matching rawtex %s\n", rawtex_path);
        goto cleanup;
    }
    full_mask = (uint8_t *)md_calloc(W * H, 1);
    rgb = (uint8_t *)md_calloc(W * H * 3, 1);
    if (full_mask == NULL || rgb == NULL) goto cleanup;
    for (i = 0; i < W * H; i++) {
        size_t x = i % W, y = i / W;
        size_t mi = (y / (size_t)opts.stride_pixels) * MW
                  + x / (size_t)opts.stride_pixels;
        double a = mask[mi] ? 0.25 + 0.35 * clamp01(score[mi]) : 0.0;
        double g = texture[i];
        full_mask[i] = mask[mi] ? 255 : 0;
        rgb[i * 3] = (uint8_t)((1.0 - a) * g + a * 16.0 + 0.5);
        rgb[i * 3 + 1] = (uint8_t)((1.0 - a) * g + a * 255.0 + 0.5);
        rgb[i * 3 + 2] = (uint8_t)((1.0 - a) * g + a * 48.0 + 0.5);
        if (mask[mi]) result.highlighted_pixels++;
    }
    snprintf(path, sizeof path, "%s_marble_mask.png", prefix);
    VesPng_write_gray(path, full_mask, (int)W, (int)H);
    snprintf(path, sizeof path, "%s_marble_mask.pgm", prefix);
    if (write_pgm(path, full_mask, W, H) != 0) goto cleanup;
    snprintf(path, sizeof path, "%s_marble_overlay.png", prefix);
    VesPng_write_rgb(path, rgb, (int)W, (int)H);
    snprintf(path, sizeof path, "%s_marble_report.json", prefix);
    {
        FILE *report = fopen(path, "wb");
        if (report == NULL) goto cleanup;
        fprintf(report,
            "{\n \"schema\":\"raw-mesh-uv-marble-audit-v2\",\n"
            " \"definition\":\"RAW/mesh UV misalignment; darkness is not evidence\",\n"
            " \"signals\":{\"centered_surface_support\":\"RAW tangent structure confirms that the mesh is on a surface\","
            "\"broad_uv_frame_disagreement\":\"primary marble signal: on-surface RAW tangent structure disagrees with mesh U/V throughout a broad UV neighborhood\","
            "\"surface_miss\":\"separate non-marble diagnostic: a coherent nearby normal depth has more RAW support than mesh depth\"},\n"
            " \"classification\":\"marble requires on-mesh RAW support plus broad UV-frame disagreement; normal-offset gain is reported as a crack/surface miss and never makes marble by itself\",\n"
            " \"ordering\":\"smooth every candidate-depth plane in UV, then select depth; never select noisy per-cell depths first\",\n"
            " \"policy\":{\"stride_pixels\":%d,\"depth_range\":%.9g,"
            "\"depth_step\":%.9g,\"tensor_radius\":%.9g,"
            "\"score_threshold\":%.9g,\"minimum_region_cells\":%d},\n"
            " \"raster\":{\"width\":%zu,\"height\":%zu},\n"
            " \"audit_grid\":{\"width\":%zu,\"height\":%zu,"
            "\"sampled_cells\":%zu,\"highlighted_cells\":%zu,"
            "\"highlighted_regions\":%zu},\n"
            " \"highlighted_pixels\":%zu,\n"
            " \"score_percentiles\":{\"p50\":%.9g,\"p90\":%.9g,"
            "\"p95\":%.9g,\"p99\":%.9g}\n}\n",
            opts.stride_pixels, opts.depth_range, opts.depth_step,
            opts.tensor_radius, opts.score_threshold, opts.minimum_region_cells,
            W, H, MW, MH, sampled, highlighted,
            result.highlighted_regions, result.highlighted_pixels,
            result.score_p50, result.score_p90, result.score_p95,
            result.score_p99);
        fclose(report);
    }
    fprintf(stderr,
        "  [marble] sampled %zu/%zu audit cells; highlighted %zu cells in "
        "%zu regions (%.2f%% of sheet pixels); score p50/p90/p95/p99 "
        "%.3f/%.3f/%.3f/%.3f\n",
        sampled, mn, highlighted, result.highlighted_regions,
        W * H ? 100.0 * (double)result.highlighted_pixels / (double)(W * H) : 0.0,
        result.score_p50, result.score_p90, result.score_p95, result.score_p99);
    if (stats != NULL) *stats = result;
    rc = 0;

cleanup:
    if (image_arena != NULL) Arena_dispose(&image_arena);
    free(sample); free(probe); free(profile); free(smooth_profile);
    free(support_center); free(support_best);
    free(gain); free(axis); free(axis_broad); free(disorder); free(miss);
    free(raw_score); free(score); free(valid); free(mask);
    free(full_mask); free(rgb);
    return rc;
}
