/* ink_letters.c -- see ink_letters.h. */
#include "ink_letters.h"
#include "sparse_solve.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/ves_omp.h"

void InkLettersOpts_default(InkLettersOpts *o)
{
    o->smooth_px = 2.0;
    o->t_hi = 0.60;
    o->t_lo = 0.45;
    o->close_px = 6;
    o->min_area = 600.0;
    o->max_area = 120000.0;
    o->row_link_h = 8.0;
    o->row_tol_h = 0.5;
    o->row_min = 3;
    o->warp_px = 32.0;
    o->w_row = 1.0;
    o->w_conf = 1.0;
    o->w_smooth = 0.25;
    o->w_anchor = 1e-6;
}

/* ------------------------------------------------------------------ detection */

static void il_blur(const float *in, float *out, size_t W, size_t H, double sigma)
{
    int R = (int)ceil(3.0 * sigma), k = 0;
    double *g = NULL, gs = 0.0;
    float *tmp = NULL;
    long y = 0;
    if (!(sigma > 0.0)) { memcpy(out, in, W * H * sizeof *out); return; }
    g = (double *)malloc((size_t)(2 * R + 1) * sizeof *g);
    tmp = (float *)malloc(W * H * sizeof *tmp);
    if (g == NULL || tmp == NULL) { free(g); free(tmp); memcpy(out, in, W * H * sizeof *out); return; }
    for (k = -R; k <= R; k++) { g[k + R] = exp(-0.5 * (double)k * k / (sigma * sigma)); gs += g[k + R]; }
    for (k = 0; k <= 2 * R; k++) g[k] /= gs;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (y = 0; y < (long)H; y++) {
        size_t x = 0;
        for (x = 0; x < W; x++) {
            double acc = 0.0, ws = 0.0;
            int q = 0;
            for (q = -R; q <= R; q++) {
                long xx = (long)x + q;
                if (xx < 0 || xx >= (long)W) continue;
                acc += g[q + R] * in[(size_t)y * W + (size_t)xx]; ws += g[q + R];
            }
            tmp[(size_t)y * W + x] = (float)(acc / ws);
        }
    }
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (y = 0; y < (long)H; y++) {
        size_t x = 0;
        for (x = 0; x < W; x++) {
            double acc = 0.0, ws = 0.0;
            int q = 0;
            for (q = -R; q <= R; q++) {
                long yy = y + q;
                if (yy < 0 || yy >= (long)H) continue;
                acc += g[q + R] * tmp[(size_t)yy * W + x]; ws += g[q + R];
            }
            out[(size_t)y * W + x] = (float)(acc / ws);
        }
    }
    free(g); free(tmp);
}

/* 8-connected labelling of on[] by union-find over the raster; labels 1.. into lab (0 = off). Returns count. */
static int32_t il_find(int32_t *par, int32_t a)
{
    while (par[a] != a) { par[a] = par[par[a]]; a = par[a]; }
    return a;
}

static size_t il_label(const uint8_t *on, size_t W, size_t H, int32_t *lab)
{
    size_t n = W * H, k = 0, next = 0;
    int32_t *par = (int32_t *)malloc(n * sizeof *par);
    int32_t *map = NULL;
    if (par == NULL) return 0;
    for (k = 0; k < n; k++) par[k] = (int32_t)k;
    for (k = 0; k < n; k++) {
        size_t x = k % W, y = k / W;
        if (!on[k]) continue;
        if (x > 0 && on[k - 1]) { int32_t a = il_find(par, (int32_t)k), b = il_find(par, (int32_t)(k - 1)); if (a != b) par[a] = b; }
        if (y > 0) {
            size_t u = k - W;
            if (on[u]) { int32_t a = il_find(par, (int32_t)k), b = il_find(par, (int32_t)u); if (a != b) par[a] = b; }
            if (x > 0 && on[u - 1]) { int32_t a = il_find(par, (int32_t)k), b = il_find(par, (int32_t)(u - 1)); if (a != b) par[a] = b; }
            if (x + 1 < W && on[u + 1]) { int32_t a = il_find(par, (int32_t)k), b = il_find(par, (int32_t)(u + 1)); if (a != b) par[a] = b; }
        }
    }
    map = (int32_t *)calloc(n, sizeof *map);
    if (map == NULL) { free(par); return 0; }
    for (k = 0; k < n; k++) {
        int32_t r = 0;
        if (!on[k]) { lab[k] = 0; continue; }
        r = il_find(par, (int32_t)k);
        if (map[r] == 0) map[r] = (int32_t)++next;
        lab[k] = map[r];
    }
    free(par); free(map);
    return next;
}

int InkLetters_detect(const uint8_t *ink, const uint8_t *mask, size_t W, size_t H, const InkLettersOpts *oin,
                      InkLetter **letters, size_t *nout)
{
    InkLettersOpts o;
    size_t n = W * H, k = 0, ncomp = 0, nmerge = 0, nl = 0;
    float *p = NULL, *s = NULL;
    uint8_t *on = NULL, *dil = NULL, *seeded = NULL;
    int32_t *lab = NULL, *mlab = NULL;
    double *acc = NULL;          /* per merged label: w, wx, wy, wxx, wyy, wxy, cnt, maxp, bx0, by0, bx1, by1 */
    InkLetter *L = NULL;
    int rc = -1;
    if (oin == NULL) InkLettersOpts_default(&o); else o = *oin;
    *letters = NULL; *nout = 0;
    p = (float *)malloc(n * sizeof *p);
    s = (float *)malloc(n * sizeof *s);
    on = (uint8_t *)calloc(n, 1);
    dil = (uint8_t *)calloc(n, 1);
    lab = (int32_t *)malloc(n * sizeof *lab);
    mlab = (int32_t *)malloc(n * sizeof *mlab);
    if (p == NULL || s == NULL || on == NULL || dil == NULL || lab == NULL || mlab == NULL) goto done;
    for (k = 0; k < n; k++) p[k] = (float)ink[k] / 255.0f;
    il_blur(p, s, W, H, o.smooth_px);
    /* hysteresis: components of s >= t_lo (inside the mask) holding a seed s >= t_hi */
    for (k = 0; k < n; k++) on[k] = (uint8_t)(s[k] >= (float)o.t_lo && (mask == NULL || mask[k]));
    ncomp = il_label(on, W, H, lab);
    seeded = (uint8_t *)calloc(ncomp + 1, 1);
    if (seeded == NULL) goto done;
    for (k = 0; k < n; k++) if (lab[k] && s[k] >= (float)o.t_hi) seeded[lab[k]] = 1;
    for (k = 0; k < n; k++) on[k] = (uint8_t)(lab[k] && seeded[lab[k]]);
    /* strokes of one letter: merge components within close_px (a square dilation of the kept pixels) */
    {
        long y = 0;
        int r = o.close_px > 0 ? o.close_px : 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (y = 0; y < (long)H; y++) {
            size_t x = 0;
            for (x = 0; x < W; x++) {
                long yy = 0, xx = 0;
                int hit = 0;
                for (yy = y - r; yy <= y + r && !hit; yy++) {
                    if (yy < 0 || yy >= (long)H) continue;
                    for (xx = (long)x - r; xx <= (long)x + r; xx++) {
                        if (xx < 0 || xx >= (long)W) continue;
                        if (on[(size_t)yy * W + (size_t)xx]) { hit = 1; break; }
                    }
                }
                dil[(size_t)y * W + x] = (uint8_t)hit;
            }
        }
    }
    nmerge = il_label(dil, W, H, mlab);
    acc = (double *)calloc((nmerge + 1) * 12, sizeof *acc);
    if (acc == NULL) goto done;
    for (k = 0; k <= nmerge; k++) { acc[k * 12 + 8] = 1e300; acc[k * 12 + 9] = 1e300; acc[k * 12 + 10] = -1e300; acc[k * 12 + 11] = -1e300; }
    for (k = 0; k < n; k++) {
        double x = (double)(k % W) + 0.5, y = (double)(k / W) + 0.5, w = 0.0;
        double *a = NULL;
        if (!on[k]) continue;
        a = &acc[(size_t)mlab[k] * 12];
        w = s[k];
        a[0] += w; a[1] += w * x; a[2] += w * y; a[3] += w * x * x; a[4] += w * y * y; a[5] += w * x * y;
        a[6] += 1.0;
        if (s[k] > a[7]) a[7] = s[k];
        if (x - 0.5 < a[8]) a[8] = x - 0.5;
        if (y - 0.5 < a[9]) a[9] = y - 0.5;
        if (x + 0.5 > a[10]) a[10] = x + 0.5;
        if (y + 0.5 > a[11]) a[11] = y + 0.5;
    }
    L = (InkLetter *)calloc(nmerge + 1, sizeof *L);
    if (L == NULL) goto done;
    {
        int32_t *slot = (int32_t *)malloc((nmerge + 1) * sizeof *slot);
        double *ext = NULL;
        if (slot == NULL) goto done;
        for (k = 1; k <= nmerge; k++) {
            double *a = &acc[k * 12], mx = 0.0, my = 0.0, cxx = 0.0, cyy = 0.0, cxy = 0.0;
            InkLetter *e = &L[nl];
            slot[k] = -1;
            if (a[6] < o.min_area || a[6] > o.max_area || !(a[0] > 0.0)) continue;
            mx = a[1] / a[0]; my = a[2] / a[0];
            cxx = a[3] / a[0] - mx * mx; cyy = a[4] / a[0] - my * my; cxy = a[5] / a[0] - mx * my;
            e->cx = mx; e->cy = my;
            e->angle = 0.5 * atan2(2.0 * cxy, cxx - cyy);
            e->bx0 = a[8]; e->by0 = a[9]; e->bx1 = a[10]; e->by1 = a[11];
            e->area = a[6]; e->mean_p = a[0] / a[6]; e->max_p = a[7];
            e->row = -1;
            e->ha = e->hb = 0.0;
            slot[k] = (int32_t)nl++;
        }
        /* exact OBB extents: min/max projections of the letter's pixels on its axes */
        ext = (double *)malloc((nl + 1) * 4 * sizeof *ext);
        if (ext == NULL) { free(slot); goto done; }
        for (k = 0; k < nl; k++) { ext[k * 4] = ext[k * 4 + 2] = 1e300; ext[k * 4 + 1] = ext[k * 4 + 3] = -1e300; }
        for (k = 0; k < n; k++) {
            int32_t sl = 0;
            double x = 0.0, y = 0.0, ca = 0.0, sa = 0.0, pa = 0.0, pb = 0.0;
            if (!on[k] || (sl = slot[mlab[k]]) < 0) continue;
            x = (double)(k % W) + 0.5 - L[sl].cx; y = (double)(k / W) + 0.5 - L[sl].cy;
            ca = cos(L[sl].angle); sa = sin(L[sl].angle);
            pa = x * ca + y * sa; pb = -x * sa + y * ca;
            if (pa < ext[sl * 4]) ext[sl * 4] = pa;
            if (pa > ext[sl * 4 + 1]) ext[sl * 4 + 1] = pa;
            if (pb < ext[sl * 4 + 2]) ext[sl * 4 + 2] = pb;
            if (pb > ext[sl * 4 + 3]) ext[sl * 4 + 3] = pb;
        }
        for (k = 0; k < nl; k++) {
            double ca = cos(L[k].angle), sa = sin(L[k].angle);
            double ma = 0.5 * (ext[k * 4] + ext[k * 4 + 1]), mb = 0.5 * (ext[k * 4 + 2] + ext[k * 4 + 3]);
            L[k].ha = 0.5 * (ext[k * 4 + 1] - ext[k * 4]) + 0.5;
            L[k].hb = 0.5 * (ext[k * 4 + 3] - ext[k * 4 + 2]) + 0.5;
            /* the OBB is centred on its extent, not the weighted centre */
            L[k].cx += ma * ca - mb * sa;
            L[k].cy += ma * sa + mb * ca;
        }
        free(ext);
        free(slot);
    }
    *letters = L; *nout = nl; L = NULL;
    rc = 0;
done:
    free(p); free(s); free(on); free(dil); free(lab); free(mlab); free(acc); free(seeded); free(L);
    return rc;
}

/* ------------------------------------------------------------------ rows */

static int il_dcmp(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return x < y ? -1 : x > y ? 1 : 0;
}

/* Weighted robust fit y = a + b x + c x^2 (c = 0 with fewer than 5 points); x centred on xm. */
static void il_rowfit(const InkLetter *L, size_t n, int row, double scale, InkRow *R)
{
    double coef[3] = { 0.0, 0.0, 0.0 }, xm = 0.0, wsum = 0.0;
    size_t k = 0, m = 0;
    int it = 0, deg = 2;
    for (k = 0; k < n; k++) if (L[k].row == row) { xm += L[k].cx; m++; }
    if (m == 0) return;
    xm /= (double)m;
    if (m < 8) deg = 1;
    {
        double lo = 1e300, hi = -1e300, hh = 0.0;
        for (k = 0; k < n; k++) if (L[k].row == row) {
            if (L[k].cx < lo) lo = L[k].cx;
            if (L[k].cx > hi) hi = L[k].cx;
            hh += (L[k].by1 - L[k].by0) / (double)m;
        }
        if (hi - lo < 10.0 * hh) deg = 1;
    }
    for (it = 0; it < 4; it++) {
        double A[9], B[3];
        int i = 0, j = 0;
        memset(A, 0, sizeof A); memset(B, 0, sizeof B);
        for (k = 0; k < n; k++) {
            double x = 0.0, phi[3], w = 0.0, r = 0.0;
            if (L[k].row != row) continue;
            x = (L[k].cx - xm) / 1000.0;
            phi[0] = 1.0; phi[1] = x; phi[2] = x * x;
            r = L[k].cy - (coef[0] + coef[1] * x + coef[2] * x * x);
            w = sqrt(L[k].area);
            if (it > 0 && fabs(r) > scale) w *= scale / fabs(r);        /* Huber */
            for (i = 0; i <= deg; i++) {
                for (j = 0; j <= deg; j++) A[i * 3 + j] += w * phi[i] * phi[j];
                B[i] += w * phi[i] * L[k].cy;
            }
        }
        /* Cramer / Gauss on the (deg+1) system */
        {
            double M[3][4];
            int d = deg + 1, c = 0, r2 = 0;
            for (i = 0; i < d; i++) { for (j = 0; j < d; j++) M[i][j] = A[i * 3 + j]; M[i][d] = B[i]; }
            for (c = 0; c < d; c++) {
                int piv = c;
                for (r2 = c + 1; r2 < d; r2++) if (fabs(M[r2][c]) > fabs(M[piv][c])) piv = r2;
                if (fabs(M[piv][c]) < 1e-12) break;
                if (piv != c) for (j = 0; j <= d; j++) { double t = M[c][j]; M[c][j] = M[piv][j]; M[piv][j] = t; }
                for (r2 = 0; r2 < d; r2++) {
                    double f = 0.0;
                    if (r2 == c) continue;
                    f = M[r2][c] / M[c][c];
                    for (j = c; j <= d; j++) M[r2][j] -= f * M[c][j];
                }
            }
            for (i = 0; i < 3; i++) coef[i] = i < d && fabs(M[i][i]) > 1e-12 ? M[i][d] / M[i][i] : 0.0;
        }
    }
    /* back to raw x: y = c0 + c1 (x - xm)/1000 + c2 ((x - xm)/1000)^2 */
    {
        double s1 = coef[1] / 1000.0, s2 = coef[2] / 1e6;
        R->a = coef[0] - s1 * xm + s2 * xm * xm;
        R->b = s1 - 2.0 * s2 * xm;
        R->c = s2;
    }
    R->x0 = 1e300; R->x1 = -1e300; R->n = m;
    {
        double *h = (double *)malloc(m * sizeof *h), e = 0.0;
        size_t q = 0;
        for (k = 0; k < n; k++) {
            double r = 0.0;
            if (L[k].row != row) continue;
            if (L[k].cx < R->x0) R->x0 = L[k].cx;
            if (L[k].cx > R->x1) R->x1 = L[k].cx;
            r = L[k].cy - (R->a + R->b * L[k].cx + R->c * L[k].cx * L[k].cx);
            e += r * r;
            if (h != NULL) h[q++] = L[k].by1 - L[k].by0;
            wsum += 1.0;
        }
        R->rms = sqrt(e / wsum);
        if (h != NULL) { qsort(h, q, sizeof *h, il_dcmp); R->height = h[q / 2]; free(h); }
    }
}

int InkLetters_rows(InkLetter *L, size_t n, const InkLettersOpts *oin, InkRow **rows, size_t *nrows, double *h_med)
{
    InkLettersOpts o;
    double *h = NULL, hm = 0.0;
    int32_t *par = NULL, *rid = NULL;
    size_t *cnt = NULL, k = 0, j = 0, nr = 0;
    InkRow *R = NULL;
    int pass = 0;
    if (oin == NULL) InkLettersOpts_default(&o); else o = *oin;
    *rows = NULL; *nrows = 0; if (h_med) *h_med = 0.0;
    if (n == 0) return 0;
    h = (double *)malloc(n * sizeof *h);
    par = (int32_t *)malloc(n * sizeof *par);
    rid = (int32_t *)malloc(n * sizeof *rid);
    cnt = (size_t *)calloc(n, sizeof *cnt);
    R = (InkRow *)calloc(n, sizeof *R);
    if (h == NULL || par == NULL || rid == NULL || cnt == NULL || R == NULL) {
        free(h); free(par); free(rid); free(cnt); free(R);
        return -1;
    }
    for (k = 0; k < n; k++) h[k] = L[k].by1 - L[k].by0;
    qsort(h, n, sizeof *h, il_dcmp);
    hm = h[n / 2];
    for (k = 0; k < n; k++) par[k] = (int32_t)k;
    for (k = 0; k < n; k++)
        for (j = k + 1; j < n; j++) {
            double gx = fabs(L[k].cx - L[j].cx);
            if (gx > o.row_link_h * hm || fabs(L[k].cy - L[j].cy) > o.row_tol_h * hm + 0.05 * gx) continue;
            {
                int32_t a = (int32_t)k, b = (int32_t)j;
                while (par[a] != a) a = par[a] = par[par[a]];
                while (par[b] != b) b = par[b] = par[par[b]];
                if (a != b) par[a] = b;
            }
        }
    for (k = 0; k < n; k++) {
        int32_t a = (int32_t)k;
        while (par[a] != a) a = par[a];
        cnt[a]++;
        rid[k] = a;
    }
    for (k = 0; k < n; k++) L[k].row = -1;
    {   /* rows numbered by first appearance */
        int32_t *num = (int32_t *)malloc(n * sizeof *num);
        if (num == NULL) { free(h); free(par); free(rid); free(cnt); free(R); return -1; }
        for (k = 0; k < n; k++) num[k] = -1;
        for (k = 0; k < n; k++) {
            if (cnt[rid[k]] < (size_t)o.row_min) continue;
            if (num[rid[k]] < 0) num[rid[k]] = (int32_t)nr++;
            L[k].row = num[rid[k]];
        }
        free(num);
    }
    /* fit; merge fragments of one line (each one's curve predicts the other's letters within row_tol_h on
     * average); reassign every letter to the nearest row curve within row_tol_h; refit */
    for (pass = 0; pass < 2; pass++) {
        size_t r = 0;
        for (r = 0; r < nr; r++) il_rowfit(L, n, (int)r, 0.25 * hm, &R[r]);
        if (pass == 0) {
            int merged = 1;
            while (merged && nr > 1) {
                size_t r1 = 0, r2 = 0;
                merged = 0;
                for (r1 = 0; r1 < nr && !merged; r1++)
                    for (r2 = r1 + 1; r2 < nr && !merged; r2++) {
                        double e = 0.0, cnt = 0.0;
                        for (k = 0; k < n; k++) {
                            const InkRow *other = NULL;
                            if (L[k].row == (int)r1) other = &R[r2];
                            else if (L[k].row == (int)r2) other = &R[r1];
                            else continue;
                            e += fabs(L[k].cy - (other->a + other->b * L[k].cx + other->c * L[k].cx * L[k].cx));
                            cnt += 1.0;
                        }
                        if (cnt > 0.0 && e / cnt <= o.row_tol_h * hm) {
                            for (k = 0; k < n; k++) {
                                if (L[k].row == (int)r2) L[k].row = (int)r1;
                                else if (L[k].row > (int)r2) L[k].row--;
                            }
                            nr--;
                            for (r = 0; r < nr; r++) il_rowfit(L, n, (int)r, 0.25 * hm, &R[r]);
                            merged = 1;
                        }
                    }
            }
        }
        if (pass == 1) break;
        for (k = 0; k < n; k++) {
            double best = o.row_tol_h * hm;
            int br = -1;
            for (r = 0; r < nr; r++) {
                double y = R[r].a + R[r].b * L[k].cx + R[r].c * L[k].cx * L[k].cx;
                if (L[k].cx < R[r].x0 - o.row_link_h * hm || L[k].cx > R[r].x1 + o.row_link_h * hm) continue;
                if (fabs(L[k].cy - y) < best) { best = fabs(L[k].cy - y); br = (int)r; }
            }
            L[k].row = br;
        }
    }
    free(h); free(par); free(rid); free(cnt);
    *rows = R; *nrows = nr;
    if (h_med) *h_med = hm;
    return 0;
}

/* ------------------------------------------------------------------ warp */

void InkWarp_free(InkWarp *w)
{
    if (w == NULL) return;
    free(w->dx); free(w->dy); free(w->level);
    memset(w, 0, sizeof *w);
}

void InkWarp_eval(const InkWarp *w, double x, double y, double *dx, double *dy)
{
    double fx = x / w->step, fy = y / w->step, ax = 0.0, ay = 0.0;
    long i = 0, j = 0;
    size_t a = 0;
    if (fx < 0.0) fx = 0.0;
    if (fy < 0.0) fy = 0.0;
    if (fx > (double)(w->nx - 1) - 1e-9) fx = (double)(w->nx - 1) - 1e-9;
    if (fy > (double)(w->ny - 1) - 1e-9) fy = (double)(w->ny - 1) - 1e-9;
    i = (long)floor(fx); j = (long)floor(fy);
    ax = fx - (double)i; ay = fy - (double)j;
    a = (size_t)j * w->nx + (size_t)i;
    *dx = (1 - ax) * (1 - ay) * w->dx[a] + ax * (1 - ay) * w->dx[a + 1] + (1 - ax) * ay * w->dx[a + w->nx] +
          ax * ay * w->dx[a + w->nx + 1];
    *dy = (1 - ax) * (1 - ay) * w->dy[a] + ax * (1 - ay) * w->dy[a + 1] + (1 - ax) * ay * w->dy[a + w->nx] +
          ax * ay * w->dy[a + w->nx + 1];
}

/* A sparse least-squares accumulator: normal equations as lower-triangle COO. */
typedef struct IlLsq {
    size_t n, cap, nt;
    int *r, *c;
    double *v, *rhs;
} IlLsq;

static int il_lsq_add(IlLsq *q, const int *idx, const double *coef, int m, double target, double w)
{
    int a = 0, b = 0;
    if (q->nt + (size_t)(m * m) > q->cap) {
        size_t nc = (q->cap ? q->cap * 2 : 1 << 16) + (size_t)(m * m);
        int *r = (int *)realloc(q->r, nc * sizeof *r), *c = NULL;
        double *v = NULL;
        if (r == NULL) return -1;
        q->r = r;
        c = (int *)realloc(q->c, nc * sizeof *c);
        if (c == NULL) return -1;
        q->c = c;
        v = (double *)realloc(q->v, nc * sizeof *v);
        if (v == NULL) return -1;
        q->v = v;
        q->cap = nc;
    }
    for (a = 0; a < m; a++) {
        q->rhs[idx[a]] += w * coef[a] * target;
        for (b = 0; b < m; b++) {
            if (idx[b] > idx[a]) continue;
            q->r[q->nt] = idx[a]; q->c[q->nt] = idx[b]; q->v[q->nt] = w * coef[a] * coef[b];
            q->nt++;
        }
    }
    return 0;
}

int InkLetters_warp(const InkLetter *L, size_t n, const InkRow *rows, size_t nrows, size_t W, size_t H,
                    const InkLettersOpts *oin, InkWarp *w)
{
    InkLettersOpts o;
    IlLsq q;
    size_t nx = 0, ny = 0, nn = 0, k = 0, i = 0, j = 0, nun = 0;
    double *x = NULL, st = 0.0, cell_w = 0.0;
    int rc = -1;
    if (oin == NULL) InkLettersOpts_default(&o); else o = *oin;
    memset(w, 0, sizeof *w);
    memset(&q, 0, sizeof q);
    st = o.warp_px;
    nx = (size_t)ceil((double)W / st) + 1;
    ny = (size_t)ceil((double)H / st) + 1;
    nn = nx * ny;
    nun = 2 * nn + nrows;
    w->nx = nx; w->ny = ny; w->step = st; w->nrows = nrows;
    w->dx = (double *)calloc(nn, sizeof *w->dx);
    w->dy = (double *)calloc(nn, sizeof *w->dy);
    w->level = (double *)calloc(nrows + 1, sizeof *w->level);
    q.n = nun;
    q.rhs = (double *)calloc(nun, sizeof *q.rhs);
    x = (double *)calloc(nun, sizeof *x);
    if (w->dx == NULL || w->dy == NULL || w->level == NULL || q.rhs == NULL || x == NULL) goto done;
#define IL_DX(ii, jj) ((int)((jj) * nx + (ii)))
#define IL_DY(ii, jj) ((int)(nn + (jj) * nx + (ii)))
    /* row constraints: every row's fitted curve, sampled each step over its letters' extent (+- half a row
     * height), maps to one height: f(x) + Dy(x, f(x)) - level_r = 0.  The curve, not the letter centres: a
     * letter's centroid depends on its shape (omicron vs lambda), the line it was written on does not. */
    {
        size_t r = 0;
        (void)L; (void)n;
        for (r = 0; r < nrows; r++) {
            double x = 0.0, xa = rows[r].x0 - 0.5 * rows[r].height, xb = rows[r].x1 + 0.5 * rows[r].height;
            for (x = xa; x <= xb; x += 0.5 * st) {
                double y = rows[r].a + rows[r].b * x + rows[r].c * x * x, fx = x / st, fy = y / st, ax = 0.0, ay = 0.0, cf[5];
                int id[5];
                long ii = (long)floor(fx), jj = (long)floor(fy);
                if (ii < 0 || jj < 0 || (size_t)ii + 1 >= nx || (size_t)jj + 1 >= ny) continue;
                ax = fx - (double)ii; ay = fy - (double)jj;
                id[0] = IL_DY(ii, jj); cf[0] = (1 - ax) * (1 - ay);
                id[1] = IL_DY(ii + 1, jj); cf[1] = ax * (1 - ay);
                id[2] = IL_DY(ii, jj + 1); cf[2] = (1 - ax) * ay;
                id[3] = IL_DY(ii + 1, jj + 1); cf[3] = ax * ay;
                id[4] = (int)(2 * nn + r); cf[4] = -1.0;
                if (il_lsq_add(&q, id, cf, 5, -y, o.w_row * st) != 0) goto done;
            }
        }
    }
    /* per cell: Cauchy-Riemann (conformal) residuals of grad D, weighted by the cell area */
    cell_w = o.w_conf * st * st;
    for (j = 0; j + 1 < ny; j++)
        for (i = 0; i + 1 < nx; i++) {
            int id[8];
            double cf[8], h2 = 0.5 / st;
            /* dxDx = (Dx(i+1,j) + Dx(i+1,j+1) - Dx(i,j) - Dx(i,j+1)) / 2st ; dyDx likewise; same for Dy */
            /* r1 = dxDx - dyDy */
            id[0] = IL_DX(i + 1, j); cf[0] = h2; id[1] = IL_DX(i + 1, j + 1); cf[1] = h2;
            id[2] = IL_DX(i, j); cf[2] = -h2; id[3] = IL_DX(i, j + 1); cf[3] = -h2;
            id[4] = IL_DY(i, j + 1); cf[4] = -h2; id[5] = IL_DY(i + 1, j + 1); cf[5] = -h2;
            id[6] = IL_DY(i, j); cf[6] = h2; id[7] = IL_DY(i + 1, j); cf[7] = h2;
            if (il_lsq_add(&q, id, cf, 8, 0.0, cell_w) != 0) goto done;
            /* r2 = dyDx + dxDy */
            id[0] = IL_DX(i, j + 1); cf[0] = h2; id[1] = IL_DX(i + 1, j + 1); cf[1] = h2;
            id[2] = IL_DX(i, j); cf[2] = -h2; id[3] = IL_DX(i + 1, j); cf[3] = -h2;
            id[4] = IL_DY(i + 1, j); cf[4] = h2; id[5] = IL_DY(i + 1, j + 1); cf[5] = h2;
            id[6] = IL_DY(i, j); cf[6] = -h2; id[7] = IL_DY(i, j + 1); cf[7] = -h2;
            if (il_lsq_add(&q, id, cf, 8, 0.0, cell_w) != 0) goto done;
        }
    /* smoothness: second differences of both components, and the anchor */
    for (j = 0; j < ny; j++)
        for (i = 0; i < nx; i++) {
            int id[3], comp = 0;
            double cf[3] = { 1.0 / st, -2.0 / st, 1.0 / st };
            for (comp = 0; comp < 2; comp++) {
                int base = comp ? (int)nn : 0;
                if (i > 0 && i + 1 < nx) {
                    id[0] = base + IL_DX(i - 1, j); id[1] = base + IL_DX(i, j); id[2] = base + IL_DX(i + 1, j);
                    if (il_lsq_add(&q, id, cf, 3, 0.0, o.w_smooth * st * st) != 0) goto done;
                }
                if (j > 0 && j + 1 < ny) {
                    id[0] = base + IL_DX(i, j - 1); id[1] = base + IL_DX(i, j); id[2] = base + IL_DX(i, j + 1);
                    if (il_lsq_add(&q, id, cf, 3, 0.0, o.w_smooth * st * st) != 0) goto done;
                }
                {
                    int a1 = base + IL_DX(i, j);
                    double one = 1.0;
                    if (il_lsq_add(&q, &a1, &one, 1, 0.0, o.w_anchor) != 0) goto done;
                }
            }
        }
    {   /* the row levels: a weak anchor at their fitted mean height */
        size_t r = 0;
        for (r = 0; r < nrows; r++) {
            int a1 = (int)(2 * nn + r);
            double one = 1.0, mid = 0.5 * (rows[r].x0 + rows[r].x1);
            if (il_lsq_add(&q, &a1, &one, 1, rows[r].a + rows[r].b * mid + rows[r].c * mid * mid, o.w_anchor) != 0)
                goto done;
        }
    }
#undef IL_DX
#undef IL_DY
    if (Sparse_solve_sym((int)nun, (int)q.nt, q.r, q.c, q.v, q.rhs, x, SPARSE_SPD) != 0) {
        fprintf(stderr, "[ink_letters] warp solve failed (%zu unknowns)\n", nun);
        goto done;
    }
    memcpy(w->dx, x, nn * sizeof *x);
    memcpy(w->dy, x + nn, nn * sizeof *x);
    memcpy(w->level, x + 2 * nn, nrows * sizeof *x);
    {   /* diagnostics: straightness before (about each row's own mean height) and after; local scale/rotation */
        double e0 = 0.0, e1 = 0.0, cnt = 0.0, *mean = (double *)calloc(nrows + 1, sizeof *mean);
        double *mc = (double *)calloc(nrows + 1, sizeof *mc);
        w->scale_min = 1e300; w->scale_max = 0.0; w->rot_max_deg = 0.0;
        if (mean != NULL && mc != NULL) {
            for (k = 0; k < n; k++) if (L[k].row >= 0 && (size_t)L[k].row < nrows) { mean[L[k].row] += L[k].cy; mc[L[k].row] += 1.0; }
            for (k = 0; k < n; k++) {
                double ddx = 0.0, ddy = 0.0, e = 0.0, jx[4], sc = 0.0, rot = 0.0;
                int r = L[k].row;
                if (r < 0 || (size_t)r >= nrows) continue;
                InkWarp_eval(w, L[k].cx, L[k].cy, &ddx, &ddy);
                e = L[k].cy - mean[r] / mc[r];
                e0 += e * e;
                e = L[k].cy + ddy - w->level[r];
                e1 += e * e;
                cnt += 1.0;
                {   /* J = I + grad D by central differences */
                    double ax = 0.0, ay = 0.0, bx = 0.0, by = 0.0, hx = 0.5 * st;
                    InkWarp_eval(w, L[k].cx + hx, L[k].cy, &ax, &ay);
                    InkWarp_eval(w, L[k].cx - hx, L[k].cy, &bx, &by);
                    jx[0] = 1.0 + (ax - bx) / (2 * hx); jx[2] = (ay - by) / (2 * hx);
                    InkWarp_eval(w, L[k].cx, L[k].cy + hx, &ax, &ay);
                    InkWarp_eval(w, L[k].cx, L[k].cy - hx, &bx, &by);
                    jx[1] = (ax - bx) / (2 * hx); jx[3] = 1.0 + (ay - by) / (2 * hx);
                    sc = sqrt(fabs(jx[0] * jx[3] - jx[1] * jx[2]));
                    rot = atan2(jx[2] - jx[1], jx[0] + jx[3]) * 180.0 / 3.14159265358979;
                    if (sc < w->scale_min) w->scale_min = sc;
                    if (sc > w->scale_max) w->scale_max = sc;
                    if (fabs(rot) > w->rot_max_deg) w->rot_max_deg = fabs(rot);
                }
            }
        }
        /* rms_before: letters about their own row's mean height (0 for a level row); rms_after: about the level */
        w->rms_before = cnt > 0.0 ? sqrt(e0 / cnt) : 0.0;
        w->rms_after = cnt > 0.0 ? sqrt(e1 / cnt) : 0.0;
        if (w->scale_min > w->scale_max) w->scale_min = w->scale_max = 1.0;
        free(mean); free(mc);
    }
    rc = 0;
done:
    free(q.r); free(q.c); free(q.v); free(q.rhs); free(x);
    if (rc != 0) InkWarp_free(w);
    return rc;
}

/* ------------------------------------------------------------------ drawing */

static void il_put(uint8_t *rgb, size_t W, size_t H, long x, long y, uint8_t r, uint8_t g, uint8_t b)
{
    if (x < 0 || y < 0 || (size_t)x >= W || (size_t)y >= H) return;
    rgb[((size_t)y * W + (size_t)x) * 3] = r;
    rgb[((size_t)y * W + (size_t)x) * 3 + 1] = g;
    rgb[((size_t)y * W + (size_t)x) * 3 + 2] = b;
}

void InkDraw_line(uint8_t *rgb, size_t W, size_t H, double x0, double y0, double x1, double y1, int thick,
                  uint8_t r, uint8_t g, uint8_t b)
{
    double len = sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0));
    long steps = (long)ceil(len * 2.0) + 1, s = 0;
    int h = thick / 2;
    for (s = 0; s <= steps; s++) {
        double t = (double)s / (double)steps;
        long cx = (long)floor(x0 + t * (x1 - x0)), cy = (long)floor(y0 + t * (y1 - y0)), dx = 0, dy = 0;
        for (dy = -h; dy <= h; dy++)
            for (dx = -h; dx <= h; dx++) il_put(rgb, W, H, cx + dx, cy + dy, r, g, b);
    }
}

void InkDraw_obb(uint8_t *rgb, size_t W, size_t H, const InkLetter *L, int thick, uint8_t r, uint8_t g, uint8_t b)
{
    double ca = cos(L->angle), sa = sin(L->angle), px[4], py[4];
    int k = 0;
    for (k = 0; k < 4; k++) {
        double a = (k == 0 || k == 3) ? -L->ha : L->ha, c = (k < 2) ? -L->hb : L->hb;
        px[k] = L->cx + a * ca - c * sa;
        py[k] = L->cy + a * sa + c * ca;
    }
    for (k = 0; k < 4; k++) InkDraw_line(rgb, W, H, px[k], py[k], px[(k + 1) % 4], py[(k + 1) % 4], thick, r, g, b);
}

/* 5x7 glyphs: 7 rows, bit 4 = leftmost column */
static const char il_glyph_chars[] = "0123456789abcdefghijklmnopqrstuvwxyz .-=:/%()+,";
static const uint8_t il_glyphs[][7] = {
    {0x0E,0x11,0x13,0x15,0x19,0x11,0x0E}, {0x04,0x0C,0x04,0x04,0x04,0x04,0x0E}, {0x0E,0x11,0x01,0x02,0x04,0x08,0x1F},
    {0x1F,0x02,0x04,0x02,0x01,0x11,0x0E}, {0x02,0x06,0x0A,0x12,0x1F,0x02,0x02}, {0x1F,0x10,0x1E,0x01,0x01,0x11,0x0E},
    {0x06,0x08,0x10,0x1E,0x11,0x11,0x0E}, {0x1F,0x01,0x02,0x04,0x08,0x08,0x08}, {0x0E,0x11,0x11,0x0E,0x11,0x11,0x0E},
    {0x0E,0x11,0x11,0x0F,0x01,0x02,0x0C},
    {0x00,0x00,0x0E,0x01,0x0F,0x11,0x0F}, {0x10,0x10,0x16,0x19,0x11,0x11,0x1E}, {0x00,0x00,0x0E,0x10,0x10,0x11,0x0E},
    {0x01,0x01,0x0D,0x13,0x11,0x11,0x0F}, {0x00,0x00,0x0E,0x11,0x1F,0x10,0x0E}, {0x06,0x09,0x08,0x1C,0x08,0x08,0x08},
    {0x00,0x0F,0x11,0x11,0x0F,0x01,0x0E}, {0x10,0x10,0x16,0x19,0x11,0x11,0x11}, {0x04,0x00,0x0C,0x04,0x04,0x04,0x0E},
    {0x02,0x00,0x06,0x02,0x02,0x12,0x0C}, {0x10,0x10,0x12,0x14,0x18,0x14,0x12}, {0x0C,0x04,0x04,0x04,0x04,0x04,0x0E},
    {0x00,0x00,0x1A,0x15,0x15,0x11,0x11}, {0x00,0x00,0x16,0x19,0x11,0x11,0x11}, {0x00,0x00,0x0E,0x11,0x11,0x11,0x0E},
    {0x00,0x00,0x1E,0x11,0x1E,0x10,0x10}, {0x00,0x00,0x0D,0x13,0x0F,0x01,0x01}, {0x00,0x00,0x16,0x19,0x10,0x10,0x10},
    {0x00,0x00,0x0E,0x10,0x0E,0x01,0x1E}, {0x08,0x08,0x1C,0x08,0x08,0x09,0x06}, {0x00,0x00,0x11,0x11,0x11,0x13,0x0D},
    {0x00,0x00,0x11,0x11,0x11,0x0A,0x04}, {0x00,0x00,0x11,0x11,0x15,0x15,0x0A}, {0x00,0x00,0x11,0x0A,0x04,0x0A,0x11},
    {0x00,0x00,0x11,0x11,0x0F,0x01,0x0E}, {0x00,0x00,0x1F,0x02,0x04,0x08,0x1F},
    {0x00,0x00,0x00,0x00,0x00,0x00,0x00}, {0x00,0x00,0x00,0x00,0x00,0x0C,0x0C}, {0x00,0x00,0x00,0x1F,0x00,0x00,0x00},
    {0x00,0x00,0x1F,0x00,0x1F,0x00,0x00}, {0x00,0x0C,0x0C,0x00,0x0C,0x0C,0x00}, {0x00,0x01,0x02,0x04,0x08,0x10,0x00},
    {0x18,0x19,0x02,0x04,0x08,0x13,0x03}, {0x02,0x04,0x08,0x08,0x08,0x04,0x02}, {0x08,0x04,0x02,0x02,0x02,0x04,0x08},
    {0x00,0x04,0x04,0x1F,0x04,0x04,0x00}, {0x00,0x00,0x00,0x00,0x0C,0x04,0x08}
};

void InkDraw_text(uint8_t *rgb, size_t W, size_t H, int x, int y, int scale, const char *s,
                  uint8_t r, uint8_t g, uint8_t b)
{
    int cx = x;
    for (; *s; s++) {
        char ch = (*s >= 'A' && *s <= 'Z') ? (char)(*s - 'A' + 'a') : *s;
        const char *p = strchr(il_glyph_chars, ch);
        int gi = p != NULL && ch != 0 ? (int)(p - il_glyph_chars) : -1, row = 0, col = 0, dy = 0, dx = 0;
        if (gi >= 0)
            for (row = 0; row < 7; row++)
                for (col = 0; col < 5; col++) {
                    if (!((il_glyphs[gi][row] >> (4 - col)) & 1)) continue;
                    for (dy = 0; dy < scale; dy++)
                        for (dx = 0; dx < scale; dx++)
                            il_put(rgb, W, H, cx + col * scale + dx, y + row * scale + dy, r, g, b);
                }
        cx += 6 * scale;
    }
}

/* ------------------------------------------------------------------ selftest */

int InkLetters_selftest(void)
{
    /* A curved row of 12 synthetic letters (discs and bars, 60-110 px) on a 2400 x 900 map, the row
     * y = 330 + 6e-5 (x - 1200)^2 (66 px of sag), and a second row 280 px below bending the same way.  Detection must find
     * the 24 letters, rows must chain 2 x 12, and the warp must bring the letter centres onto their levels
     * within 2 px while staying near-conformal (local scale within 5%): one smooth bend of the sheet is
     * straightened by a near-rotation, whereas rows of different curvature cannot all be (conformally). */
    const size_t W = 2400, H = 900;
    uint8_t *ink = (uint8_t *)malloc(W * H);
    InkLetter *L = NULL;
    InkRow *R = NULL;
    InkWarp w;
    size_t n = 0, nr = 0, k = 0, y = 0, x = 0;
    double hm = 0.0;
    int fails = 0, i = 0;
    memset(&w, 0, sizeof w);
    if (ink == NULL) return 1;
    for (k = 0; k < W * H; k++) ink[k] = 70;
    for (i = 0; i < 24; i++) {
        int row = i / 12, c = i % 12;
        /* both rows bend the same way (a smooth distortion of the sheet, as a flattening error is) */
        double cx = 150.0 + 180.0 * c, cy = (row == 0 ? 330.0 : 610.0) + 6e-5 * (cx - 1200.0) * (cx - 1200.0);
        int shape = c % 3;
        for (y = (size_t)(cy - 70); y < (size_t)(cy + 70); y++)
            for (x = (size_t)(cx - 70); x < (size_t)(cx + 70); x++) {
                double dx = (double)x + 0.5 - cx, dy = (double)y + 0.5 - cy, rr = sqrt(dx * dx + dy * dy);
                int on = shape == 0 ? (rr < 45.0 && rr > 25.0) : shape == 1 ? (fabs(dx) < 12.0 && fabs(dy) < 50.0)
                                                                            : (rr < 35.0);
                if (on) ink[y * W + x] = 220;
            }
    }
    if (InkLetters_detect(ink, NULL, W, H, NULL, &L, &n) != 0) { free(ink); return 1; }
    if (n != 24) { fprintf(stderr, "FAIL ink_letters: %zu letters detected (want 24)\n", n); fails++; }
    if (InkLetters_rows(L, n, NULL, &R, &nr, &hm) != 0 || nr != 2) {
        fprintf(stderr, "FAIL ink_letters: %zu rows (want 2)\n", nr);
        fails++;
    }
    for (k = 0; k < nr && !fails; k++) if (R[k].n != 12) { fprintf(stderr, "FAIL ink_letters: row %zu has %zu\n", k, R[k].n); fails++; }
    if (!fails && InkLetters_warp(L, n, R, nr, W, H, NULL, &w) != 0) fails++;
    if (!fails) {
        fprintf(stderr, "  ink_letters: %zu letters, %zu rows (h %.0f px); straightness rms %.1f -> %.2f px; scale %.3f..%.3f, "
                "rotation <= %.1f deg\n", n, nr, hm, w.rms_before, w.rms_after, w.scale_min, w.scale_max, w.rot_max_deg);
        if (w.rms_after > 2.0 || w.rms_before < 10.0 || w.scale_min < 0.95 || w.scale_max > 1.05) fails++;
    }
    {   /* the font and the drawing stay inside the raster */
        uint8_t *rgb = (uint8_t *)calloc(200 * 40 * 3, 1);
        if (rgb != NULL) {
            InkDraw_text(rgb, 200, 40, 2, 2, 2, "row 1: 10 mm", 255, 255, 255);
            InkDraw_line(rgb, 200, 40, -10, -10, 300, 60, 3, 255, 0, 0);
            free(rgb);
        }
    }
    free(ink); free(L); free(R);
    InkWarp_free(&w);
    if (fails) fprintf(stderr, "FAIL ink_letters\n");
    return fails ? 1 : 0;
}
