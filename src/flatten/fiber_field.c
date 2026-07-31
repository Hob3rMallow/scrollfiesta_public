/* fiber_field.c -- see fiber_field.h. Structure-tensor fiber orientation +
 * coherence over a baked rawtex texture, with a mod-180 (line) and a mod-90
 * (4-RoSy grid) mode. Cold path (analysis/diagnostics): scratch is malloc/
 * free'd; only the output arrays come from the arena. */
#include "fiber_field.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define FIBER_PI 3.14159265358979323846

/* ---- separable Gaussian blur (reflect-101 border) ------------------------ */

static float *gauss_kernel(double sigma, int *out_r)
{
    int r = (int)ceil(3.0 * sigma), i;
    float *k = NULL;
    double s = 0.0;
    if (r < 1) r = 1;
    k = (float *)malloc((size_t)(2 * r + 1) * sizeof *k);
    if (k == NULL) { *out_r = 0; return NULL; }
    for (i = -r; i <= r; i++) {
        double v = exp(-(double)(i * i) / (2.0 * sigma * sigma));
        k[i + r] = (float)v;
        s += v;
    }
    for (i = 0; i < 2 * r + 1; i++) k[i] = (float)(k[i] / s);
    *out_r = r;
    return k;
}

static int reflect101(int i, int n)
{
    if (n == 1) return 0;
    while (i < 0 || i >= n) {
        if (i < 0) i = -i;
        if (i >= n) i = 2 * n - 2 - i;
    }
    return i;
}

/* dst != src; tmp is a caller-provided W*H scratch buffer. */
static void gauss_blur(const float *src, float *dst, int W, int H,
                       double sigma, float *tmp)
{
    int r = 0, x, y, t;
    float *k = gauss_kernel(sigma, &r);
    if (k == NULL) { memcpy(dst, src, (size_t)W * (size_t)H * sizeof *dst); return; }
    for (y = 0; y < H; y++) {                          /* horizontal: src -> tmp */
        const float *row = src + (size_t)y * W;
        float *orow = tmp + (size_t)y * W;
        for (x = 0; x < W; x++) {
            double a = 0.0;
            for (t = -r; t <= r; t++)
                a += (double)k[t + r] * row[reflect101(x + t, W)];
            orow[x] = (float)a;
        }
    }
    for (y = 0; y < H; y++) {                           /* vertical: tmp -> dst */
        for (x = 0; x < W; x++) {
            double a = 0.0;
            for (t = -r; t <= r; t++)
                a += (double)k[t + r] * tmp[(size_t)reflect101(y + t, H) * W + x];
            dst[(size_t)y * W + x] = (float)a;
        }
    }
    free(k);
}

/* ---- core ---------------------------------------------------------------- */

int FiberField_compute(Arena_T arena, const uint8_t *img, int W, int H,
                       double grad_sigma, double tensor_sigma, double axis_tol_deg,
                       int rosy, FiberField *out)
{
    size_t N = (size_t)W * (size_t)H, p;
    float *fimg = NULL, *Is = NULL, *gx = NULL, *gy = NULL;
    float *Ac = NULL, *As = NULL, *En = NULL, *cov = NULL, *covs = NULL, *tmp = NULL;
    int x, y, b;
    double wsum = 0.0, cohsum = 0.0, axiserr_sum = 0.0, axisfrac_sum = 0.0;
    double range, binw;
    size_t nvalid = 0;

    if (arena == NULL || img == NULL || out == NULL || W <= 0 || H <= 0) return -1;
    if (grad_sigma   <= 0.0) grad_sigma   = 1.0;
    if (tensor_sigma <= 0.0) tensor_sigma = 4.0;
    if (axis_tol_deg <= 0.0) axis_tol_deg = 15.0;
    if (rosy != 4) rosy = 2;
    range = (rosy == 4) ? (FIBER_PI / 2.0) : FIBER_PI;    /* orientation domain */
    binw = (range * 180.0 / FIBER_PI) / FIBER_NBINS;

    memset(out, 0, sizeof *out);
    out->W = W; out->H = H; out->rosy = rosy; out->range_deg = range * 180.0 / FIBER_PI;
    out->theta = (float   *)Arena_calloc(arena, N, sizeof(float),   __FILE__, __LINE__);
    out->coh   = (float   *)Arena_calloc(arena, N, sizeof(float),   __FILE__, __LINE__);
    out->valid = (uint8_t *)Arena_calloc(arena, N, sizeof(uint8_t), __FILE__, __LINE__);

    tmp  = (float *)malloc(N * sizeof *tmp);
    fimg = (float *)malloc(N * sizeof *fimg);
    if (tmp == NULL || fimg == NULL) { free(tmp); free(fimg); return -1; }
    for (p = 0; p < N; p++) fimg[p] = (float)img[p];

    /* 1. denoise, gradient */
    Is = (float *)malloc(N * sizeof *Is);
    if (Is == NULL) { free(tmp); free(fimg); return -1; }
    gauss_blur(fimg, Is, W, H, grad_sigma, tmp);
    free(fimg); fimg = NULL;

    gx = (float *)malloc(N * sizeof *gx);
    gy = (float *)malloc(N * sizeof *gy);
    if (gx == NULL || gy == NULL) { free(tmp); free(Is); free(gx); free(gy); return -1; }
    for (y = 0; y < H; y++) {
        int yu = (y > 0) ? y - 1 : 0, yd = (y < H - 1) ? y + 1 : H - 1;
        for (x = 0; x < W; x++) {
            int xl = (x > 0) ? x - 1 : 0, xr = (x < W - 1) ? x + 1 : W - 1;
            gx[(size_t)y * W + x] = 0.5f * (Is[(size_t)y * W + xr] - Is[(size_t)y * W + xl]);
            gy[(size_t)y * W + x] = 0.5f * (Is[(size_t)yd * W + x] - Is[(size_t)yu * W + x]);
        }
    }
    free(Is); Is = NULL;

    /* 2. structure tensor J = [[gx^2, gx gy],[gx gy, gy^2]], stored as the
     *    doubled-angle vector (Ac,As) = (Sxx-Syy, 2 Sxy) plus energy En = Sxx+Syy.
     *    (The grid/rosy4 folding of the recovered orientation happens in step 4;
     *    folding the *tensor* itself would blend the two families to the diagonal.) */
    Ac = (float *)malloc(N * sizeof *Ac);
    As = (float *)malloc(N * sizeof *As);
    En = (float *)malloc(N * sizeof *En);
    if (Ac == NULL || As == NULL || En == NULL) {
        free(tmp); free(gx); free(gy); free(Ac); free(As); free(En); return -1;
    }
    for (p = 0; p < N; p++) {
        double a = gx[p], g = gy[p];
        En[p] = (float)(a * a + g * g);      /* gradient energy E = |grad I|^2 */
        Ac[p] = (float)(a * a - g * g);      /* E cos2phi */
        As[p] = (float)(2.0 * a * g);        /* E sin2phi */
    }
    free(gx); gx = NULL; free(gy); gy = NULL;

    /* 3. smooth the accumulators at the integration scale (in place via tmp) */
    {
        float *sc = (float *)malloc(N * sizeof *sc);
        if (sc == NULL) { free(tmp); free(Ac); free(As); free(En); return -1; }
        gauss_blur(Ac, sc, W, H, tensor_sigma, tmp); memcpy(Ac, sc, N * sizeof *Ac);
        gauss_blur(As, sc, W, H, tensor_sigma, tmp); memcpy(As, sc, N * sizeof *As);
        gauss_blur(En, sc, W, H, tensor_sigma, tmp); memcpy(En, sc, N * sizeof *En);
        free(sc);
    }

    /* coverage mask smoothed at the same scale -> drop under-supported texels */
    cov  = (float *)malloc(N * sizeof *cov);
    covs = (float *)malloc(N * sizeof *covs);
    if (cov == NULL || covs == NULL) {
        free(tmp); free(Ac); free(As); free(En); free(cov); free(covs); return -1;
    }
    for (p = 0; p < N; p++) cov[p] = (img[p] > 0) ? 1.0f : 0.0f;
    gauss_blur(cov, covs, W, H, tensor_sigma, tmp);
    free(cov); cov = NULL;

    /* 4. extract orientation + coherence + stats */
    for (p = 0; p < N; p++) {
        double A = Ac[p], B = As[p], E = En[p];
        double mag = sqrt(A * A + B * B);
        double coh = (E > 1e-12) ? mag / E : 0.0;
        double phi, theta, deg, d90, axis_err;
        int valid = (img[p] > 0) && (covs[p] > 0.9f);

        if (coh < 0.0) coh = 0.0; else if (coh > 1.0) coh = 1.0;
        phi = 0.5 * atan2(B, A);                    /* gradient orientation (mod pi) */
        theta = phi + FIBER_PI / 2.0;               /* fiber runs perpendicular */
        theta = fmod(theta, FIBER_PI); if (theta < 0.0) theta += FIBER_PI;
        if (rosy == 4) theta = fmod(theta, FIBER_PI / 2.0);   /* fold to grid orientation (mod 90) */

        if (!valid) { out->theta[p] = 0.0f; out->coh[p] = 0.0f; out->valid[p] = 0; continue; }
        out->theta[p] = (float)theta;
        out->coh[p]   = (float)coh;
        out->valid[p] = 1;

        deg = theta * 180.0 / FIBER_PI;                 /* [0, range_deg) */
        b = (int)(deg / binw); if (b < 0) b = 0; else if (b >= FIBER_NBINS) b = FIBER_NBINS - 1;
        out->hist[b] += coh;
        d90 = fmod(deg, 90.0); axis_err = fmin(d90, 90.0 - d90);   /* distance to nearest axis */
        wsum += coh; cohsum += coh; nvalid++;
        axiserr_sum += coh * axis_err;
        if (axis_err <= axis_tol_deg) axisfrac_sum += coh;
    }

    if (wsum > 0.0) {
        int best = 0;
        for (b = 0; b < FIBER_NBINS; b++) out->hist[b] /= wsum;
        for (b = 1; b < FIBER_NBINS; b++) if (out->hist[b] > out->hist[best]) best = b;
        out->peak_deg = (best + 0.5) * binw;
        out->axis_frac = axisfrac_sum / wsum;
        out->mean_axis_err_deg = axiserr_sum / wsum;
    }
    out->mean_coh = (nvalid > 0) ? cohsum / (double)nvalid : 0.0;
    out->frac_valid = (double)nvalid / (double)N;

    free(tmp); free(Ac); free(As); free(En); free(covs);
    return 0;
}

/* ---- selftest ------------------------------------------------------------ */

/* coherence-weighted mean orientation (deg) via order-`rosy` doubled angle */
static double ff_mean_orient_deg(const FiberField *f)
{
    double sx = 0.0, sy = 0.0, range = f->range_deg * FIBER_PI / 180.0;
    size_t p, N = (size_t)f->W * (size_t)f->H;
    double k = (double)f->rosy;
    for (p = 0; p < N; p++) {
        if (!f->valid[p]) continue;
        sx += (double)f->coh[p] * cos(k * (double)f->theta[p]);
        sy += (double)f->coh[p] * sin(k * (double)f->theta[p]);
    }
    {
        double m = atan2(sy, sx) / k * 180.0 / FIBER_PI;
        double rd = range * 180.0 / FIBER_PI;
        while (m < 0.0) m += rd;
        while (m >= rd) m -= rd;
        return m;
    }
}

static int ff_check(const char *name, double got, double want, double tol, double wrap)
{
    double d = fabs(got - want);
    if (wrap > 0.0) { double d2 = wrap - d; if (d2 < d) d = d2; }
    if (d > tol) {
        fprintf(stderr, "[fiber_field selftest]   FAIL: %s got %.2f want %.2f (tol %.2f)\n",
                name, got, want, tol);
        return 1;
    }
    return 0;
}

static void ff_grating(uint8_t *img, int W, int H, double dirx, double diry, double lambda)
{
    int x, y;
    for (y = 0; y < H; y++)
        for (x = 0; x < W; x++) {
            double phase = 2.0 * FIBER_PI * ((double)x * dirx + (double)y * diry) / lambda;
            img[(size_t)y * W + x] = (uint8_t)(128.0 + 100.0 * sin(phase));
        }
}

/* additive orthogonal gratings -> a cross-hatch grid rotated by `deg` */
static void ff_crosshatch(uint8_t *img, int W, int H, double deg, double lambda)
{
    double a = deg * FIBER_PI / 180.0, cx = cos(a), sx = sin(a);
    int x, y;
    for (y = 0; y < H; y++)
        for (x = 0; x < W; x++) {
            double u = (double)x * cx + (double)y * sx;      /* along grid */
            double v = -(double)x * sx + (double)y * cx;     /* across grid */
            double val = 128.0 + 65.0 * sin(2.0 * FIBER_PI * u / lambda)   /* dominant family */
                               + 30.0 * sin(2.0 * FIBER_PI * v / lambda);  /* weaker family */
            if (val < 0) val = 0; if (val > 255) val = 255;
            img[(size_t)y * W + x] = (uint8_t)val;
        }
}

int FiberField_selftest(void)
{
    int fails = 0, W = 160, H = 120;
    uint8_t *img = (uint8_t *)malloc((size_t)W * H);
    Arena_T arena = Arena_new();
    FiberField f;

    if (img == NULL) { fprintf(stderr, "[fiber_field selftest] FAIL alloc\n"); return 1; }

    /* rosy2 line orientation: horizontal fibers -> ~0 deg */
    ff_grating(img, W, H, 0.0, 1.0, 12.0);
    if (FiberField_compute(arena, img, W, H, 1.0, 4.0, 15.0, 2, &f) != 0) { fprintf(stderr, "compute fail\n"); fails++; }
    else {
        fails += ff_check("r2 horiz", ff_mean_orient_deg(&f), 0.0, 6.0, 180.0);
        if (f.mean_coh < 0.5) { fprintf(stderr, "[fiber_field selftest]   FAIL: r2 horiz coh %.2f\n", f.mean_coh); fails++; }
    }
    /* rosy2 vertical -> ~90 */
    ff_grating(img, W, H, 1.0, 0.0, 12.0);
    if (FiberField_compute(arena, img, W, H, 1.0, 4.0, 15.0, 2, &f) == 0)
        fails += ff_check("r2 vert", ff_mean_orient_deg(&f), 90.0, 6.0, 180.0);
    /* rosy2 diagonal -> ~135 */
    ff_grating(img, W, H, 0.70710678, 0.70710678, 12.0);
    if (FiberField_compute(arena, img, W, H, 1.0, 4.0, 15.0, 2, &f) == 0)
        fails += ff_check("r2 diag", ff_mean_orient_deg(&f), 135.0, 8.0, 180.0);

    /* rosy4 GRID: axis-aligned cross-hatch -> grid orientation ~0 (mod 90) */
    ff_crosshatch(img, W, H, 0.0, 12.0);
    if (FiberField_compute(arena, img, W, H, 1.0, 5.0, 15.0, 4, &f) == 0) {
        fails += ff_check("r4 aligned", ff_mean_orient_deg(&f), 0.0, 6.0, 90.0);
        if (f.mean_coh < 0.3) { fprintf(stderr, "[fiber_field selftest]   FAIL: r4 aligned coh %.2f\n", f.mean_coh); fails++; }
    }
    /* rosy4 GRID rotated 20 deg -> grid orientation ~20 */
    ff_crosshatch(img, W, H, 20.0, 12.0);
    if (FiberField_compute(arena, img, W, H, 1.0, 5.0, 15.0, 4, &f) == 0)
        fails += ff_check("r4 rot20", ff_mean_orient_deg(&f), 20.0, 7.0, 90.0);

    /* flat: no orientation -> coherence ~0 */
    memset(img, 128, (size_t)W * H);
    if (FiberField_compute(arena, img, W, H, 1.0, 4.0, 15.0, 4, &f) == 0)
        if (f.mean_coh > 0.1) { fprintf(stderr, "[fiber_field selftest]   FAIL: flat coh %.2f\n", f.mean_coh); fails++; }

    if (FiberField_compute(arena, NULL, W, H, 1.0, 4.0, 15.0, 2, &f) != -1) { fprintf(stderr, "null not rejected\n"); fails++; }

    free(img);
    Arena_dispose(&arena);
    fprintf(stderr, "[fiber_field selftest] %s\n", fails == 0 ? "ok" : "FAILED");
    return fails;
}
