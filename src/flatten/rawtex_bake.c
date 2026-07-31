/* rawtex_bake.c -- see rawtex_bake.h. Extracted verbatim from obj_bake_raw.c
 * (the rasterizer + its private helpers), plus an optional face_skip mask. */
#include "rawtex_bake.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#include "../common/tiff_io.h"
#include "../common/ves_png.h"
#include "../common/tif_strip.h"
#include "../common/raw_sample.h"

/* ------------------------------------------------------------------ util
 * Trivial malloc wrappers + the rasterizer's private math helpers. These are
 * intentionally file-local duplicates of the same statics in the standalone diagnostic baker
 * (5-line wrappers not worth a shared alloc header). */
static void *xmalloc(size_t nbytes)
{
    void *p = malloc(nbytes > 0 ? nbytes : 1);
    if (p == NULL) {
        fprintf(stderr, "ERROR: out of memory (%zu bytes)\n", nbytes);
        exit(1);
    }
    return p;
}

static void *xcalloc(size_t count, size_t size)
{
    void *p = calloc(count > 0 ? count : 1, size);
    if (p == NULL) {
        fprintf(stderr, "ERROR: out of memory (%zu x %zu)\n", count, size);
        exit(1);
    }
    return p;
}

static double gray_of(double v, double lo, double hi)
{
    double g = (v - lo) / (hi - lo);
    if (g < 0.0) g = 0.0;
    if (g > 1.0) g = 1.0;
    return g;
}

static double face_sigma(const float *verts, const float *uv,
                         size_t a, size_t b, size_t c, double *out_smin)
{
    /* Sander et al. 2001 singular values of the UV->3D map */
    double s1 = (double)uv[a * 2 + 0], t1 = (double)uv[a * 2 + 1];
    double s2 = (double)uv[b * 2 + 0], t2 = (double)uv[b * 2 + 1];
    double s3 = (double)uv[c * 2 + 0], t3 = (double)uv[c * 2 + 1];
    double A2 = (s2 - s1) * (t3 - t1) - (s3 - s1) * (t2 - t1);
    double Ss[3], St[3];
    double aa = 0.0, bb = 0.0, cc = 0.0, disc = 0.0;
    int k = 0;
    *out_smin = 0.0;
    if (fabs(A2) < 1e-12) return 1e6;   /* UV-degenerate: infinite compression */
    for (k = 0; k < 3; k++) {
        double q1 = (double)verts[a * 3 + k], q2 = (double)verts[b * 3 + k];
        double q3 = (double)verts[c * 3 + k];
        Ss[k] = (q1 * (t2 - t3) + q2 * (t3 - t1) + q3 * (t1 - t2)) / A2;
        St[k] = (q1 * (s3 - s2) + q2 * (s1 - s3) + q3 * (s2 - s1)) / A2;
    }
    aa = Ss[0] * Ss[0] + Ss[1] * Ss[1] + Ss[2] * Ss[2];
    bb = Ss[0] * St[0] + Ss[1] * St[1] + Ss[2] * St[2];
    cc = St[0] * St[0] + St[1] * St[1] + St[2] * St[2];
    disc = sqrt((aa - cc) * (aa - cc) + 4.0 * bb * bb);
    *out_smin = sqrt(((aa + cc - disc) > 0.0 ? (aa + cc - disc) : 0.0) / 2.0);
    return sqrt((aa + cc + disc) / 2.0);
}

static uint8_t log2_to_u8(double sigma)
{
    double t = 0.0;
    if (sigma < 1e-9) return 0;
    t = 128.0 + 42.5 * (log(sigma) / 0.6931471805599453);
    if (t < 0.0) t = 0.0;
    if (t > 255.0) t = 255.0;
    return (uint8_t)(t + 0.5);
}

static int cmp_double(const void *x, const void *y)
{
    double a = *(const double *)x, b = *(const double *)y;
    return (a < b) ? -1 : (a > b) ? 1 : 0;
}

void Rawtex_stretch_window(const double *val, const uint8_t *has, size_t nv,
                           double pct_lo, double pct_hi,
                           double *out_lo, double *out_hi)
{
    size_t hist[256];
    size_t i = 0, n = 0, cum = 0, tlo = 0, thi = 0;
    int b = 0, lo = 0, hi = 255;

    memset(hist, 0, sizeof hist);
    for (i = 0; i < nv; i++) {
        int bin = 0;
        if (!has[i]) continue;
        bin = (int)(val[i] + 0.5);
        if (bin < 0) bin = 0;
        if (bin > 255) bin = 255;
        hist[bin]++;
        n++;
    }
    *out_lo = 0.0; *out_hi = 255.0;
    if (n == 0) return;
    tlo = (size_t)(pct_lo / 100.0 * (double)n);
    thi = (size_t)(pct_hi / 100.0 * (double)n);
    cum = 0;
    for (b = 0; b < 256; b++) {
        cum += hist[b];
        if (cum > tlo) { lo = b; break; }
    }
    cum = 0;
    for (b = 0; b < 256; b++) {
        cum += hist[b];
        if (cum >= thi) { hi = b; break; }
    }
    if (hi <= lo) { lo = 0; hi = 255; }
    *out_lo = (double)lo;
    *out_hi = (double)hi;
}

int Rawtex_write_tif(const char *path, CubeTable *ct,
                     const float *verts, const float *uv, size_t nv,
                     const int32_t *faces, size_t nf,
                     const float *normals,
                     const uint8_t *face_skip,
                     double range, int nsteps,
                     double du, double dv, double lo, double hi,
                     double stretch_ratio, double stretch_floor,
                     double max_edge3d,
                     const DiagOpts *diag,
                     size_t *out_W, size_t *out_H,
                     double *out_fill, size_t *out_multi,
                     size_t *out_skip_uv, size_t *out_skip_3d)
{
    size_t skip_uv = 0, skip_3d = 0;
    double umax = 0.0, vmax = 0.0;
    size_t W = 0, H = 0, i = 0, f = 0, px = 0, filled = 0, multi = 0;
    double *sum = NULL;
    uint32_t *cnt = NULL;
    uint8_t *img = NULL;
    /* diag buffers */
    uint8_t *dsmax = NULL, *dsmin = NULL, *dverr = NULL, *dcls = NULL;
    uint8_t *cover_uv = NULL, *cover_3d = NULL;
    uint8_t *is_smear = NULL;   /* per-face skip-uv flag for the 3D dump */
    double zoff = 0.0;
    int rc = 0;

    *out_W = 0; *out_H = 0; *out_fill = 0.0; *out_multi = 0;
    if (out_skip_uv != NULL) *out_skip_uv = 0;
    if (out_skip_3d != NULL) *out_skip_3d = 0;
    if (nv == 0 || nf == 0 || uv == NULL || du <= 0.0 || dv <= 0.0) return -1;
    for (i = 0; i < nv; i++) {
        if ((double)uv[i * 2 + 0] > umax) umax = (double)uv[i * 2 + 0];
        if ((double)uv[i * 2 + 1] > vmax) vmax = (double)uv[i * 2 + 1];
    }
    W = (size_t)floor(umax / du + 0.5) + 1;
    H = (size_t)floor(vmax / dv + 0.5) + 1;
    if (W > (size_t)1 << 20 || H > (size_t)1 << 20 || W * H > (size_t)1 << 26) {
        fprintf(stderr, "ERROR: raster %zux%zu unreasonable (du/dv too small?)\n",
                W, H);
        return -1;
    }
    sum = (double *)xcalloc(W * H, sizeof(double));
    cnt = (uint32_t *)xcalloc(W * H, sizeof(uint32_t));
    img = (uint8_t *)xcalloc(W * H, sizeof(uint8_t));
    if (diag != NULL) {
        size_t n_probe = 0, stride = 0;
        double *zv = NULL;
        dsmax = (uint8_t *)xcalloc(W * H, 1);
        dsmin = (uint8_t *)xcalloc(W * H, 1);
        dverr = (uint8_t *)xcalloc(W * H, 1);
        dcls = (uint8_t *)xcalloc(W * H, 1);
        cover_uv = (uint8_t *)xcalloc(W * H, 1);
        cover_3d = (uint8_t *)xcalloc(W * H, 1);
        /* zmin estimate: median of (z - vt.v) over a vertex sample -- exact
         * on correctly-mapped verts, robust to the broken ones */
        stride = (nv > 200000) ? nv / 200000 : 1;
        zv = (double *)xmalloc((nv / stride + 1) * sizeof(double));
        for (i = 0; i < nv; i += stride)
            zv[n_probe++] = (double)verts[i * 3 + 0] - (double)uv[i * 2 + 1];
        qsort(zv, n_probe, sizeof(double), cmp_double);
        zoff = zv[n_probe / 2];
        free(zv);
        if (diag->smear_obj != NULL)
            is_smear = (uint8_t *)xcalloc(nf, 1);
    }

    for (f = 0; f < nf; f++) {
        size_t a = (size_t)faces[f * 3 + 0];
        size_t b = (size_t)faces[f * 3 + 1];
        size_t c = (size_t)faces[f * 3 + 2];
        double ua = (double)uv[a * 2 + 0] / du, va = (double)uv[a * 2 + 1] / dv;
        double ub = (double)uv[b * 2 + 0] / du, vb = (double)uv[b * 2 + 1] / dv;
        double uc = (double)uv[c * 2 + 0] / du, vc = (double)uv[c * 2 + 1] / dv;
        double A2 = (ub - ua) * (vc - va) - (vb - va) * (uc - ua);
        double lox = ua < ub ? (ua < uc ? ua : uc) : (ub < uc ? ub : uc);
        double hix = ua > ub ? (ua > uc ? ua : uc) : (ub > uc ? ub : uc);
        double loy = va < vb ? (va < vc ? va : vc) : (vb < vc ? vb : vc);
        double hiy = va > vb ? (va > vc ? va : vc) : (vb > vc ? vb : vc);
        long x0 = 0, x1 = 0, y0 = 0, y1 = 0, xx = 0, yy = 0;
        int bad_uv = 0, bad_3d = 0;
        uint8_t f_smax = 0, f_smin = 0, f_verr = 0;
        if (face_skip != NULL && face_skip[f]) continue;   /* omitted entirely */
        if (fabs(A2) < 1e-12) continue;   /* degenerate in UV */
        {   /* face gates: UV stretch + max 3D edge */
            size_t vi[3];
            int e = 0;
            vi[0] = a; vi[1] = b; vi[2] = c;
            for (e = 0; e < 3; e++) {
                size_t p = vi[e], q = vi[(e + 1) % 3];
                double d3 = 0.0, dueg = 0.0, dveg = 0.0, e2d = 0.0;
                int k = 0;
                for (k = 0; k < 3; k++) {
                    double dd = (double)verts[q * 3 + k]
                                - (double)verts[p * 3 + k];
                    d3 += dd * dd;
                }
                d3 = sqrt(d3);
                dueg = (double)uv[q * 2 + 0] - (double)uv[p * 2 + 0];
                dveg = (double)uv[q * 2 + 1] - (double)uv[p * 2 + 1];
                e2d = sqrt(dueg * dueg + dveg * dveg);
                if (max_edge3d > 0.0 && d3 > max_edge3d) bad_3d = 1;
                if (stretch_ratio > 0.0
                    && e2d > (stretch_ratio * d3 > stretch_floor
                              ? stretch_ratio * d3 : stretch_floor)) bad_uv = 1;
            }
            if (bad_3d) skip_3d++;
            else if (bad_uv) skip_uv++;
            if (is_smear != NULL && bad_uv && !bad_3d) is_smear[f] = 1;
            if ((bad_3d || bad_uv) && diag == NULL) continue;
        }
        if (diag != NULL) {
            double smin = 0.0, smax = face_sigma(verts, uv, a, b, c, &smin);
            double e0 = fabs(((double)verts[a * 3 + 0] - zoff)
                             - (double)uv[a * 2 + 1]);
            double e1 = fabs(((double)verts[b * 3 + 0] - zoff)
                             - (double)uv[b * 2 + 1]);
            double e2 = fabs(((double)verts[c * 3 + 0] - zoff)
                             - (double)uv[c * 2 + 1]);
            double emax = e0 > e1 ? (e0 > e2 ? e0 : e2) : (e1 > e2 ? e1 : e2);
            f_smax = log2_to_u8(smax);
            f_smin = log2_to_u8(smin);
            emax *= 32.0;
            f_verr = (uint8_t)(emax > 255.0 ? 255.0 : emax);
        }
        x0 = (long)ceil(lox - 0.001); x1 = (long)floor(hix + 0.001);
        y0 = (long)ceil(loy - 0.001); y1 = (long)floor(hiy + 0.001);
        if (x0 < 0) x0 = 0;
        if (y0 < 0) y0 = 0;
        if (x1 >= (long)W) x1 = (long)W - 1;
        if (y1 >= (long)H) y1 = (long)H - 1;
        for (yy = y0; yy <= y1; yy++) {
            for (xx = x0; xx <= x1; xx++) {
                double pu = (double)xx, pv = (double)yy;
                double l1 = ((pu - ua) * (vc - va) - (pv - va) * (uc - ua)) / A2;
                double l2 = ((ub - ua) * (pv - va) - (vb - va) * (pu - ua)) / A2;
                double l0 = 1.0 - l1 - l2;
                size_t pi = 0;
                float p3[3], n3[3];
                double nn = 0.0, s = 0.0;
                int k = 0;
                if (l0 < -1e-6 || l1 < -1e-6 || l2 < -1e-6) continue;
                pi = (size_t)yy * W + (size_t)xx;
                if (diag != NULL) {
                    if (f_smax > dsmax[pi]) dsmax[pi] = f_smax;
                    if (dsmin[pi] == 0 || (f_smin != 0 && f_smin < dsmin[pi]))
                        dsmin[pi] = f_smin;
                    if (f_verr > dverr[pi]) dverr[pi] = f_verr;
                    if (bad_uv) cover_uv[pi] = 1;
                    else if (bad_3d) cover_3d[pi] = 1;
                }
                if (bad_uv || bad_3d) continue;   /* diagnosed, not painted */
                for (k = 0; k < 3; k++) {
                    p3[k] = (float)(l0 * (double)verts[a * 3 + k]
                                    + l1 * (double)verts[b * 3 + k]
                                    + l2 * (double)verts[c * 3 + k]);
                    n3[k] = (normals != NULL)
                            ? (float)(l0 * (double)normals[a * 3 + k]
                                      + l1 * (double)normals[b * 3 + k]
                                      + l2 * (double)normals[c * 3 + k])
                            : 0.0f;
                }
                nn = sqrt((double)n3[0] * n3[0] + (double)n3[1] * n3[1]
                          + (double)n3[2] * n3[2]);
                if (nn > 1e-6) {
                    n3[0] = (float)((double)n3[0] / nn);
                    n3[1] = (float)((double)n3[1] / nn);
                    n3[2] = (float)((double)n3[2] / nn);
                }
                s = sample_vertex(ct, p3, nn > 0.5 ? n3 : NULL, range, nsteps);
                if (s < 0.0) continue;
                sum[pi] += s;
                cnt[pi]++;
            }
        }
    }
    for (px = 0; px < W * H; px++) {
        if (cnt[px] > 0) {
            double g = gray_of(sum[px] / (double)cnt[px], lo, hi);
            img[px] = (uint8_t)(g * 255.0 + 0.5);
            filled++;
            if (cnt[px] > 1) multi++;   /* overlapping wraps (rarely: on-edge) */
        }
    }
    rc = TiffIO_save(path, img, 1, (int)H, (int)W);
    {   /* also a full-res grayscale PNG sibling (<path>.png), written from C
         * via stb -- no GDI+ round-trip (which truncates >~3k-col images). */
        char pngpath[2600];
        size_t plen = strlen(path);
        snprintf(pngpath, sizeof pngpath, "%s", path);
        if (plen > 4 && strcmp(pngpath + plen - 4, ".tif") == 0)
            snprintf(pngpath + plen - 4, sizeof pngpath - (plen - 4), ".png");
        else
            snprintf(pngpath + plen, sizeof pngpath - plen, ".png");
        VesPng_write_gray(pngpath, img, (int)W, (int)H);
    }
    {   /* multi-line strip sibling (<base>_strip.tif + .png): a readable
         * stacked view of the very wide unroll (see tif_strip.h). Every bake
         * stage thus delivers both <stage>_rawtex.tif and _rawtex_strip.tif. */
        char sp[2600];
        size_t plen = strlen(path);
        snprintf(sp, sizeof sp, "%s", path);
        if (plen > 4 && strcmp(sp + plen - 4, ".tif") == 0)
            snprintf(sp + plen - 4, sizeof sp - (plen - 4), "_strip.tif");
        else
            snprintf(sp + plen, sizeof sp - plen, "_strip.tif");
        TifStrip_write_gray(sp, img, (int)W, (int)H, 0, -1);
    }
    *out_W = W; *out_H = H;
    *out_fill = (double)filled / (double)(W * H);
    *out_multi = multi;
    if (out_skip_uv != NULL) *out_skip_uv = skip_uv;
    if (out_skip_3d != NULL) *out_skip_3d = skip_3d;

    if (diag != NULL && rc == 0) {
        char dpath[2600];
        size_t n_bg = 0, n_ok = 0, n_suv = 0, n_s3d = 0, n_multi = 0, n_dark = 0;
        for (px = 0; px < W * H; px++) {
            uint8_t cl = 0;
            if (cover_uv[px]) cl = 2;
            else if (cover_3d[px]) cl = 3;
            else if (cnt[px] > 1) cl = 4;
            else if (cnt[px] == 1)
                cl = (img[px] < diag->dark_thresh) ? 5 : 1;
            dcls[px] = cl;
            switch (cl) {
            case 0: n_bg++; break;
            case 1: n_ok++; break;
            case 2: n_suv++; break;
            case 3: n_s3d++; break;
            case 4: n_multi++; break;
            default: n_dark++; break;
            }
        }
        fprintf(stderr, "  [diag] class px: ok=%zu (%.2f%%) skip-uv=%zu "
                "(%.2f%%) skip-3d=%zu (%.2f%%) dark=%zu (%.2f%%) multi=%zu "
                "bg=%zu  (%% of non-bg non-multi)\n",
                n_ok, 100.0 * (double)n_ok
                      / (double)(n_ok + n_suv + n_s3d + n_dark + 1),
                n_suv, 100.0 * (double)n_suv
                       / (double)(n_ok + n_suv + n_s3d + n_dark + 1),
                n_s3d, 100.0 * (double)n_s3d
                       / (double)(n_ok + n_suv + n_s3d + n_dark + 1),
                n_dark, 100.0 * (double)n_dark
                        / (double)(n_ok + n_suv + n_s3d + n_dark + 1),
                n_multi, n_bg);
        {   /* worst-tile table: defect fraction per 256x64-px tile */
            size_t tw = 256, th = 64;
            size_t ntx = (W + tw - 1) / tw, nty = (H + th - 1) / th;
            typedef struct { double bad; size_t tx, ty, nbad, ncov; } Tile;
            Tile *tiles = (Tile *)xcalloc(ntx * nty, sizeof(Tile));
            size_t ti = 0, shown = 0;
            for (px = 0; px < W * H; px++) {
                size_t x = px % W, y = px / W;
                uint8_t cl = dcls[px];
                Tile *T = &tiles[(y / th) * ntx + (x / tw)];
                T->tx = x / tw; T->ty = y / th;
                if (cl == 1) T->ncov++;
                else if (cl == 2 || cl == 3 || cl == 5) { T->ncov++; T->nbad++; }
            }
            for (ti = 0; ti < ntx * nty; ti++)
                tiles[ti].bad = (tiles[ti].ncov > 300)
                    ? (double)tiles[ti].nbad / (double)tiles[ti].ncov : 0.0;
            /* selection sort of top 12 (tiny n) */
            for (shown = 0; shown < 12 && shown < ntx * nty; shown++) {
                size_t best = shown;
                Tile tmp;
                for (ti = shown; ti < ntx * nty; ti++)
                    if (tiles[ti].bad > tiles[best].bad) best = ti;
                tmp = tiles[shown]; tiles[shown] = tiles[best];
                tiles[best] = tmp;
                if (tiles[shown].bad <= 0.0) break;
                fprintf(stderr, "  [diag] hot tile %2zu: px(%zu..%zu, %zu..%zu)"
                        " u=[%.0f,%.0f] v=[%.0f,%.0f]  bad=%.1f%% (%zu/%zu)\n",
                        shown,
                        tiles[shown].tx * tw, (tiles[shown].tx + 1) * tw,
                        tiles[shown].ty * th, (tiles[shown].ty + 1) * th,
                        (double)(tiles[shown].tx * tw) * du,
                        (double)((tiles[shown].tx + 1) * tw) * du,
                        (double)(tiles[shown].ty * th) * dv,
                        (double)((tiles[shown].ty + 1) * th) * dv,
                        100.0 * tiles[shown].bad,
                        tiles[shown].nbad, tiles[shown].ncov);
            }
            free(tiles);
        }
        snprintf(dpath, sizeof dpath, "%s_diagclass.tif", diag->prefix);
        TiffIO_save(dpath, dcls, 1, (int)H, (int)W);
        snprintf(dpath, sizeof dpath, "%s_diagstretch.tif", diag->prefix);
        TiffIO_save(dpath, dsmax, 1, (int)H, (int)W);
        snprintf(dpath, sizeof dpath, "%s_diagsquash.tif", diag->prefix);
        TiffIO_save(dpath, dsmin, 1, (int)H, (int)W);
        snprintf(dpath, sizeof dpath, "%s_diagverr.tif", diag->prefix);
        TiffIO_save(dpath, dverr, 1, (int)H, (int)W);
        fprintf(stderr, "  [diag] wrote %s_diag{class,stretch,squash,verr}.tif"
                " (zmin est %.2f)\n", diag->prefix, zoff);
    }
    if (is_smear != NULL) {
        /* 3D dump of the skip-uv (smear) faces + a z-slice histogram, to
         * localize where the u-stretch defects sit on the scroll. */
        FILE *fp = fopen(diag->smear_obj, "wb");
        size_t nsm = 0;
        long zhist[16];
        double zmin = 1e300, zmax = -1e300;
        int hb = 0;
        memset(zhist, 0, sizeof zhist);
        for (i = 0; i < nv; i++) {
            double z = (double)verts[i * 3 + 0];
            if (z < zmin) zmin = z;
            if (z > zmax) zmax = z;
        }
        /* winding-jump test: for each smear face, is max|du| across its edges
         * close to an integer multiple of 2*pi*r (a winding-index disagreement)
         * or arbitrary (a QP/registration offset)? r = mean radial distance of
         * the face from the axis (axis_point z,y,x; scroll axis = Z). */
        long ratio_hist[8];   /* du/(2*pi*r) rounded: 0,1,2,3,4,5,6,>=7 */
        memset(ratio_hist, 0, sizeof ratio_hist);
        if (fp != NULL) {
            for (i = 0; i < nv; i++)
                fprintf(fp, "v %.4f %.4f %.4f\n", (double)verts[i * 3 + 0],
                        (double)verts[i * 3 + 1], (double)verts[i * 3 + 2]);
            for (f = 0; f < nf; f++) {
                double zc = 0.0, umn = 1e300, umx = -1e300, rsum = 0.0, twopir = 0.0;
                int bidx = 0, m = 0, rr = 0;
                if (!is_smear[f]) continue;
                fprintf(fp, "f %d %d %d\n", faces[f * 3 + 0] + 1,
                        faces[f * 3 + 1] + 1, faces[f * 3 + 2] + 1);
                nsm++;
                zc = ((double)verts[(size_t)faces[f*3+0]*3]
                      + (double)verts[(size_t)faces[f*3+1]*3]
                      + (double)verts[(size_t)faces[f*3+2]*3]) / 3.0;
                bidx = (zmax > zmin)
                     ? (int)(15.0 * (zc - zmin) / (zmax - zmin)) : 0;
                if (bidx < 0) bidx = 0;
                if (bidx > 15) bidx = 15;
                zhist[bidx]++;
                for (m = 0; m < 3; m++) {
                    size_t vi = (size_t)faces[f*3+m];
                    double uu = (double)uv[vi*2+0];
                    double dy = (double)verts[vi*3+1] - 3405.0;
                    double dx = (double)verts[vi*3+2] - 2878.0;
                    if (uu < umn) umn = uu;
                    if (uu > umx) umx = uu;
                    rsum += sqrt(dy*dy + dx*dx);
                }
                twopir = 6.283185307 * (rsum / 3.0);
                rr = (twopir > 1.0) ? (int)((umx - umn) / twopir + 0.5) : 7;
                if (rr < 0) rr = 0;
                if (rr > 7) rr = 7;
                ratio_hist[rr]++;
            }
            fclose(fp);
            fprintf(stderr, "  [diag] wrote %s (%zu smear faces, 3D)\n",
                    diag->smear_obj, nsm);
            fprintf(stderr, "  [diag] smear z-histogram over [%.0f,%.0f] "
                    "(16 bands):\n    ", zmin, zmax);
            for (hb = 0; hb < 16; hb++)
                fprintf(stderr, "%ld ", zhist[hb]);
            fprintf(stderr, "\n  [diag] smear du/(2pi r) rounded 0..6,>=7: ");
            for (hb = 0; hb < 8; hb++)
                fprintf(stderr, "%ld ", ratio_hist[hb]);
            fprintf(stderr, "  (>=1 = winding-index jump; 0 = sub-wrap "
                    "registration/QP)\n");
        }
        free(is_smear);
    }
    free(sum); free(cnt); free(img);
    free(dsmax); free(dsmin); free(dverr); free(dcls);
    free(cover_uv); free(cover_3d);
    return rc;
}
