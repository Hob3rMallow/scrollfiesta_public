/* fiber_frame.c -- see fiber_frame.h.  A C port of the probe that measured it (PHerc0139, 2026-10-14:
 * line 11's 2,089 px fibre climb -> 30 px). */
#include "fiber_frame.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/pipeline_constants.h"
#include "../common/ves_platform.h"
#include "../common/zarr_u8.h"
#include "sparse_solve.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

struct FiberVolume {
    ZarrU8 pres, nx, ny;
};

int FiberVolume_open(const char *root, FiberVolume_T *out)
{
    char path[1100];
    FiberVolume_T fv = NULL;
    if (root == NULL || out == NULL) return -1;
    *out = NULL;
    fv = (FiberVolume_T)calloc(1, sizeof *fv);
    if (fv == NULL) return -1;
    snprintf(path, sizeof path, "%s/presence", root);
    if (ZarrU8_open(&fv->pres, path) != 0) { free(fv); return -1; }
    snprintf(path, sizeof path, "%s/nx", root);
    if (ZarrU8_open(&fv->nx, path) != 0) { ZarrU8_close(&fv->pres); free(fv); return -1; }
    snprintf(path, sizeof path, "%s/ny", root);
    if (ZarrU8_open(&fv->ny, path) != 0) { ZarrU8_close(&fv->pres); ZarrU8_close(&fv->nx); free(fv); return -1; }
    *out = fv;
    return 0;
}

void FiberVolume_close(FiberVolume_T *fv)
{
    if (fv == NULL || *fv == NULL) return;
    fprintf(stderr, "[fiber_volume] cache %d chunks/plane: presence %zu hits/%zu loads (%zu missing), "
                    "nx %zu/%zu, ny %zu/%zu\n", ZARR_U8_CACHE_SLOTS,
            (*fv)->pres.hits, (*fv)->pres.misses, (*fv)->pres.missing,
            (*fv)->nx.hits, (*fv)->nx.misses, (*fv)->ny.hits, (*fv)->ny.misses);
    ZarrU8_close(&(*fv)->pres);
    ZarrU8_close(&(*fv)->nx);
    ZarrU8_close(&(*fv)->ny);
    free(*fv);
    *fv = NULL;
}

double FiberVolume_sample(FiberVolume_T fv, double x, double y, double z, double d[3])
{
    long iz = (long)floor(z + 0.5), iy = (long)floor(y + 0.5), ix = (long)floor(x + 0.5);
    uint8_t p = ZarrU8_get(&fv->pres, iz, iy, ix);
    double dx = 0.0, dy = 0.0, s = 0.0;
    d[0] = d[1] = d[2] = 0.0;
    if (p == 0) return 0.0;
    dx = ((double)ZarrU8_get(&fv->nx, iz, iy, ix) - 128.0) / 127.0;
    dy = ((double)ZarrU8_get(&fv->ny, iz, iy, ix) - 128.0) / 127.0;
    s = 1.0 - dx * dx - dy * dy;
    d[0] = dx;
    d[1] = dy;
    d[2] = s > 0.0 ? sqrt(s) : 0.0;
    return (double)p / 255.0;
}

/* ------------------------------------------------------------------------------------------------ small helpers */

static int ff_cmp_double(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return x < y ? -1 : x > y ? 1 : 0;
}

/* numpy-style linear-interpolated percentile of v[0..n) (sorts v in place); NaN when n == 0 */
static double ff_pct_sorted(const double *v, size_t n, double q)
{
    double pos = 0.0, f = 0.0;
    size_t i = 0;
    if (n == 0) return NAN;
    pos = q / 100.0 * (double)(n - 1);
    i = (size_t)floor(pos);
    if (i + 1 >= n) return v[n - 1];
    f = pos - (double)i;
    return v[i] * (1.0 - f) + v[i + 1] * f;
}

static FiberQuant ff_quant(double *v, size_t n)
{
    FiberQuant q;
    qsort(v, n, sizeof *v, ff_cmp_double);
    q.p50 = ff_pct_sorted(v, n, 50.0);
    q.p90 = ff_pct_sorted(v, n, 90.0);
    q.p99 = ff_pct_sorted(v, n, 99.0);
    return q;
}

static int ff_tangents(const float *X, const uint8_t *valid, int H, int W, int r, int c, double pu[3], double pv[3])
{
    size_t a = 0, b = 0, e = 0, f = 0;
    if (r < 1 || c < 1 || r > H - 2 || c > W - 2) return -1;
    a = (size_t)r * (size_t)W + (size_t)(c + 1);
    b = (size_t)r * (size_t)W + (size_t)(c - 1);
    e = (size_t)(r + 1) * (size_t)W + (size_t)c;
    f = (size_t)(r - 1) * (size_t)W + (size_t)c;
    if (!valid[a] || !valid[b] || !valid[e] || !valid[f]) return -1;
    for (int k = 0; k < 3; k++) {
        pu[k] = ((double)X[a * 3 + (size_t)k] - (double)X[b * 3 + (size_t)k]) / 2.0;
        pv[k] = ((double)X[e * 3 + (size_t)k] - (double)X[f * 3 + (size_t)k]) / 2.0;
    }
    return 0;
}

/* Orthonormal PHYSICAL tangent axes. In these axes the canvas differential is
 * B=[length,along; 0,height]. A sheared canvas is not an orthogonal frame. */
static int ff_physical_frame(const double pu[3], const double pv[3], double e1[3], double e2[3], double metric[3])
{
    double cross[3],length2=0,along=0;
    for(size_t k=0;k<3;k++){
        cross[k]=pu[(k+1)%3]*pv[(k+2)%3]-pu[(k+2)%3]*pv[(k+1)%3];
        length2+=pu[k]*pu[k];along+=pu[k]*pv[k];
    }
    double area=hypot(hypot(cross[0],cross[1]),cross[2]),length=sqrt(length2);
    if(!(length>0) || !(area*area>1e-9) || !isfinite(length) || !isfinite(area) || !isfinite(along))return -1;
    metric[0]=length;metric[1]=along/length;metric[2]=area/length;
    for(size_t k=0;k<3;k++)e1[k]=pu[k]/length;
    for(size_t k=0;k<3;k++)e2[k]=(cross[(k+1)%3]*e1[(k+2)%3]-cross[(k+2)%3]*e1[(k+1)%3])/area;
    return 0;
}

static double ff_physical_angle(const double direction[3], const double e1[3], const double e2[3])
{
    double a=0,b=0;
    for(size_t k=0;k<3;k++){a+=direction[k]*e1[k];b+=direction[k]*e2[k];}
    return atan2(b,a);
}

/* Preserve theta's public meaning as a direction in canvas coordinates.
 * Fourfold averaging happens before this pullback, in the physical frame. */
static double ff_canvas_angle(double physical_angle, const double metric[3])
{
    double v=sin(physical_angle)/metric[2];
    double u=(cos(physical_angle)-metric[1]*v)/metric[0];
    return atan2(v,u);
}

/* separable Gaussian, scipy.ndimage mode 'reflect' (d c b a | a b c d), truncate 4 sigma */
static void ff_gauss(double *img, int h, int w, double sigma, double *tmp)
{
    int R = (int)(4.0 * sigma + 0.5);
    double k[64], ks = 0.0;
    if (R > 31) R = 31;
    for (int i = -R; i <= R; i++) {
        k[i + R] = exp(-0.5 * (double)(i * i) / (sigma * sigma));
        ks += k[i + R];
    }
    for (int i = 0; i <= 2 * R; i++) k[i] /= ks;
    for (int y = 0; y < h; y++)
        for (int x = 0; x < w; x++) {
            double s = 0.0;
            for (int i = -R; i <= R; i++) {
                int xx = x + i;
                while (xx < 0 || xx >= w) xx = xx < 0 ? -xx - 1 : 2 * w - xx - 1;
                s += k[i + R] * img[(size_t)y * (size_t)w + (size_t)xx];
            }
            tmp[(size_t)y * (size_t)w + (size_t)x] = s;
        }
    for (int y = 0; y < h; y++)
        for (int x = 0; x < w; x++) {
            double s = 0.0;
            for (int i = -R; i <= R; i++) {
                int yy = y + i;
                while (yy < 0 || yy >= h) yy = yy < 0 ? -yy - 1 : 2 * h - yy - 1;
                s += k[i + R] * tmp[(size_t)yy * (size_t)w + (size_t)x];
            }
            img[(size_t)y * (size_t)w + (size_t)x] = s;
        }
}

/* binary dilation / erosion with a disk of radius r; outside the array counts as 0 (scipy's border_value 0) */
static void ff_morph(const uint8_t *in, uint8_t *out, int h, int w, int r, int dilate)
{
    for (int y = 0; y < h; y++)
        for (int x = 0; x < w; x++) {
            int hit = dilate ? 0 : 1;
            for (int dy = -r; dy <= r && (dilate ? !hit : hit); dy++)
                for (int dx = -r; dx <= r; dx++) {
                    int yy = y + dy, xx = x + dx, v = 0;
                    if (dx * dx + dy * dy > r * r) continue;
                    v = (yy >= 0 && yy < h && xx >= 0 && xx < w) ? in[(size_t)yy * (size_t)w + (size_t)xx] : 0;
                    if (dilate && v) { hit = 1; break; }
                    if (!dilate && !v) { hit = 0; break; }
                }
            out[(size_t)y * (size_t)w + (size_t)x] = (uint8_t)hit;
        }
}

/* scipy binary_fill_holes: background 4-connected to the border stays background, the rest becomes foreground */
static int ff_fill_holes(uint8_t *m, int h, int w)
{
    size_t n = (size_t)h * (size_t)w, head = 0, tail = 0;
    size_t *q = (size_t *)malloc(n * sizeof *q);
    uint8_t *out = (uint8_t *)calloc(n, 1);
    if (q == NULL || out == NULL) { free(q); free(out); return -1; }
    for (int y = 0; y < h; y++)
        for (int x = 0; x < w; x++) {
            size_t i = (size_t)y * (size_t)w + (size_t)x;
            if ((y == 0 || x == 0 || y == h - 1 || x == w - 1) && !m[i] && !out[i]) { out[i] = 1; q[tail++] = i; }
        }
    while (head < tail) {
        size_t i = q[head++];
        int y = (int)(i / (size_t)w), x = (int)(i % (size_t)w);
        const int ny[4] = { y - 1, y + 1, y, y }, nx[4] = { x, x, x - 1, x + 1 };
        for (int k = 0; k < 4; k++) {
            size_t j = 0;
            if (ny[k] < 0 || ny[k] >= h || nx[k] < 0 || nx[k] >= w) continue;
            j = (size_t)ny[k] * (size_t)w + (size_t)nx[k];
            if (!m[j] && !out[j]) { out[j] = 1; q[tail++] = j; }
        }
    }
    for (size_t i = 0; i < n; i++) m[i] = (uint8_t)(m[i] || !out[i]);
    free(q);
    free(out);
    return 0;
}

/* axial angles in degrees (period 180): the doubled-angle mean direction, then the median of the residuals */
static double ff_circ_median_deg(double *a, size_t n)
{
    double sx = 0.0, sy = 0.0, m = 0.0;
    for (size_t i = 0; i < n; i++) {
        sx += cos(2.0 * a[i] * M_PI / 180.0);
        sy += sin(2.0 * a[i] * M_PI / 180.0);
    }
    m = atan2(sy, sx) * 180.0 / M_PI / 2.0;
    for (size_t i = 0; i < n; i++) {
        double r = fmod(a[i] - m + 90.0, 180.0);
        if (r < 0.0) r += 180.0;
        a[i] = r - 90.0;
    }
    qsort(a, n, sizeof *a, ff_cmp_double);
    return m + ff_pct_sorted(a, n, 50.0);
}

/* np.gradient along one axis of a [h*w] grid (central inside, one-sided at the ends), divided by `step` */
static void ff_gradient(const double *g, double *out, int h, int w, int axis, double step)
{
    for (int y = 0; y < h; y++)
        for (int x = 0; x < w; x++) {
            size_t i = (size_t)y * (size_t)w + (size_t)x;
            int n = axis ? h : w, k = axis ? y : x;
            size_t s = axis ? (size_t)w : 1;
            double v = NAN;
            if (n >= 2) {
                if (k == 0) v = g[i + s] - g[i];
                else if (k == n - 1) v = g[i] - g[i - s];
                else v = (g[i + s] - g[i - s]) / 2.0;
            }
            out[i] = v / step;
        }
}

/* ------------------------------------------------------------------------------------------------------ field */

int FiberFrame_field(Arena_T arena, const float *X, const uint8_t *valid, int H, int W,
                     int x0, int y0, int x1, int y1, FiberVolume_T fv, FiberFrame *ff)
{
    const int S = FIBER_FLATTEN_SAMPLE_PX, D = FIBER_FLATTEN_HALF_DEPTH, cell = FIBER_FLATTEN_CELL_PX;
    size_t nn = 0;
    double *zr = NULL, *zi = NULL, *w = NULL, *tmp = NULL;
    if (X == NULL || valid == NULL || fv == NULL || ff == NULL || H < 3 || W < 3) return -1;
    if (x0 < 0) x0 = 0;
    if (y0 < 0) y0 = 0;
    if (x1 > W) x1 = W;
    if (y1 > H) y1 = H;
    if (x1 - x0 < 2 * cell || y1 - y0 < 2 * cell) return -1;
    memset(ff, 0, sizeof *ff);
    ff->W = W; ff->H = H; ff->x0 = x0; ff->y0 = y0; ff->x1 = x1; ff->y1 = y1;
    ff->cell = cell; ff->gy0 = y0; ff->gx0 = x0;
    ff->gh = (y1 - y0) / cell + 1;
    ff->gw = (x1 - x0) / cell + 1;
    nn = (size_t)ff->gh * (size_t)ff->gw;
    zr = (double *)calloc(nn, sizeof *zr);
    zi = (double *)calloc(nn, sizeof *zi);
    w = (double *)calloc(nn, sizeof *w);
    tmp = (double *)calloc(nn, sizeof *tmp);
    ff->theta = (float *)ARENA_ALLOC(arena, nn * sizeof *ff->theta);
    ff->weight = (float *)ARENA_ALLOC(arena, nn * sizeof *ff->weight);
    ff->coherence = (float *)ARENA_ALLOC(arena, nn * sizeof *ff->coherence);
    ff->real = (uint8_t *)ARENA_ALLOC(arena, nn);
    ff->closed = (uint8_t *)ARENA_ALLOC(arena, nn);
    ff->pu = (double *)ARENA_ALLOC(arena, nn * 3 * sizeof *ff->pu);
    ff->pv = (double *)ARENA_ALLOC(arena, nn * 3 * sizeof *ff->pv);
    ff->U = (double *)ARENA_ALLOC(arena, nn * sizeof *ff->U);
    ff->V = (double *)ARENA_ALLOC(arena, nn * sizeof *ff->V);
    if (zr == NULL || zi == NULL || w == NULL || tmp == NULL || ff->theta == NULL || ff->U == NULL) {
        free(zr); free(zi); free(w); free(tmp);
        return -1;
    }
    for (int r = y0; r < y1; r += S)
        for (int c = x0; c < x1; c += S) {
            size_t i = (size_t)r * (size_t)W + (size_t)c;
            double pu[3], pv[3], N[3], e1[3], e2[3], metric[3], ar = 0.0, ai = 0.0, aw = 0.0;
            long gr = 0, gc = 0;
            int any = 0;
            if (!valid[i] || ff_tangents(X, valid, H, W, r, c, pu, pv) != 0) continue;
            if(ff_physical_frame(pu,pv,e1,e2,metric))continue;
            for(size_t k=0;k<3;k++)N[k]=-(e1[(k+1)%3]*e2[(k+2)%3]-e1[(k+2)%3]*e2[(k+1)%3]);
            ff->samples++;
            for (int o = -D; o <= D; o++) {
                double d[3], p = 0.0, psi = 0.0;
                p = FiberVolume_sample(fv, (double)X[i * 3] + o * N[0], (double)X[i * 3 + 1] + o * N[1],
                                       (double)X[i * 3 + 2] + o * N[2], d);
                if (p <= 0.0) continue;
                psi = ff_physical_angle(d,e1,e2);
                ar += p * cos(4.0 * psi);
                ai += p * sin(4.0 * psi);
                aw += p;
                any = 1;
            }
            if (any) ff->samples_with_fibres++;
            gr = (long)rint((double)(r - y0) / cell);       /* numpy rint: half to even */
            gc = (long)rint((double)(c - x0) / cell);
            if (gr < 0 || gc < 0 || gr >= ff->gh || gc >= ff->gw) continue;
            zr[(size_t)gr * (size_t)ff->gw + (size_t)gc] += ar;
            zi[(size_t)gr * (size_t)ff->gw + (size_t)gc] += ai;
            w[(size_t)gr * (size_t)ff->gw + (size_t)gc] += aw;
        }
    ff_gauss(zr, ff->gh, ff->gw, FIBER_FLATTEN_SIGMA_PX / cell, tmp);
    ff_gauss(zi, ff->gh, ff->gw, FIBER_FLATTEN_SIGMA_PX / cell, tmp);
    ff_gauss(w, ff->gh, ff->gw, FIBER_FLATTEN_SIGMA_PX / cell, tmp);
    for (int gi = 0; gi < ff->gh; gi++)
        for (int gj = 0; gj < ff->gw; gj++) {
            size_t n = (size_t)gi * (size_t)ff->gw + (size_t)gj;
            int r = y0 + gi * cell, c = x0 + gj * cell;
            int inside = r < y1 && c < x1;
            int rc = r < H - 1 ? r : H - 1, cc = c < W - 1 ? c : W - 1;
            double pu[3] = { 0, 0, 0 }, pv[3] = { 0, 0, 0 }, e1[3],e2[3],metric[3];
            double physical_angle=atan2(zi[n],zr[n])/4.0;
            ff->weight[n] = (float)w[n];
            ff->coherence[n] = (float)(w[n] > 0 ? fmin(1.0,hypot(zr[n],zi[n])/w[n]) : 0.0);
            ff->real[n] = (uint8_t)(inside && valid[(size_t)rc * (size_t)W + (size_t)cc]
                                    && ff_tangents(X, valid, H, W, rc, cc, pu, pv) == 0
                                    && ff_physical_frame(pu,pv,e1,e2,metric)==0);
            ff->theta[n]=(float)(ff->real[n]?ff_canvas_angle(physical_angle,metric):0);
            for (int k = 0; k < 3; k++) {
                ff->pu[n * 3 + (size_t)k] = ff->real[n] ? pu[k] : 0.0;
                ff->pv[n * 3 + (size_t)k] = ff->real[n] ? pv[k] : 0.0;
            }
            ff->U[n] = ff->V[n] = NAN;
        }
    free(zr); free(zi); free(w); free(tmp);
    return 0;
}

/* ------------------------------------------------------------------------------------------------------ solve */

int FiberFrame_solve(Arena_T arena, FiberFrame *ff, FiberFrameStats *st)
{
    const int gh = ff->gh, gw = ff->gw, cell = ff->cell;
    const size_t nn = (size_t)gh * (size_t)gw;
    uint8_t *dil = NULL, *clo = NULL;
    int *ids = NULL, *rows = NULL, *cols = NULL;
    double *vals = NULL, *rhs = NULL, *sol = NULL, *gu = NULL, *gv = NULL, *tu = NULL, *tv = NULL, *we = NULL;
    int *ea = NULL, *eb = NULL;
    size_t n = 0, ne = 0, nt = 0, nreal_e = 0;
    SparseFactor_T fac = NULL;
    double su = 0.0, sv = 0.0;
    double *lh = NULL, *lv = NULL;
    size_t nl = 0;
    int rc = -1;
    (void)arena;
    memset(st, 0, sizeof *st);
    dil = (uint8_t *)malloc(nn);
    clo = (uint8_t *)malloc(nn);
    ids = (int *)malloc(nn * sizeof *ids);
    gu = (double *)malloc(nn * 2 * sizeof *gu);
    gv = (double *)malloc(nn * 2 * sizeof *gv);
    lh = (double *)malloc(nn * sizeof *lh);
    lv = (double *)malloc(nn * sizeof *lv);
    if (dil == NULL || clo == NULL || ids == NULL || gu == NULL || gv == NULL || lh == NULL || lv == NULL) goto done;
    ff_morph(ff->real, dil, gh, gw, FIBER_FLATTEN_GAP_CELLS, 1);
    ff_morph(dil, clo, gh, gw, FIBER_FLATTEN_GAP_CELLS, 0);
    for (size_t i = 0; i < nn; i++) clo[i] = (uint8_t)(clo[i] || ff->real[i]);
    if (ff_fill_holes(clo, gh, gw) != 0) goto done;
    memcpy(ff->closed, clo, nn);
    for (size_t i = 0; i < nn; i++) {
        lh[i]=lv[i]=gu[2*i]=gu[2*i+1]=gv[2*i]=gv[2*i+1]=0;
        if(ff->real[i]){
            double e1[3],e2[3],metric[3],t=(double)ff->theta[i];
            if(ff_physical_frame(ff->pu+3*i,ff->pv+3*i,e1,e2,metric))goto done;
            double hx=metric[0]*cos(t)+metric[1]*sin(t),hy=metric[2]*sin(t),length=hypot(hx,hy);
            if(!(length>0) || !isfinite(length))goto done;
            hx/=length;hy/=length;
            /* du'=H.dP and dv'=V.dP. These are covectors, not the canvas
             * tangent directions scaled by their lengths. V is perpendicular
             * to H in the physical plane, including the off-diagonal metric. */
            gu[2*i]=hx*metric[0]*cell;gu[2*i+1]=(hx*metric[1]+hy*metric[2])*cell;
            gv[2*i]=-hy*metric[0]*cell;gv[2*i+1]=(-hy*metric[1]+hx*metric[2])*cell;
            lh[i]=length;
            double vy=hx/metric[2],vx=(-hy-metric[1]*vy)/metric[0];
            lv[i]=1/hypot(vx,vy);
        }
        ids[i] = clo[i] ? (int)n++ : -1;
        if (ff->real[i]) { lh[nl] = lh[i]; lv[nl] = lv[i]; nl++; }   /* compacted copies for the medians (nl <= i) */
        if (clo[i] && !ff->real[i]) st->nodes_gap++;
        if (ff->real[i]) st->nodes_real++;
    }
    /* the edges */
    ea = (int *)malloc(nn * 2 * sizeof *ea);
    eb = (int *)malloc(nn * 2 * sizeof *eb);
    tu = (double *)malloc(nn * 2 * sizeof *tu);
    tv = (double *)malloc(nn * 2 * sizeof *tv);
    we = (double *)malloc(nn * 2 * sizeof *we);
    if (ea == NULL || eb == NULL || tu == NULL || tv == NULL || we == NULL) goto done;
    for (int dir = 0; dir < 2; dir++) {
        int dy = dir == 0 ? 0 : 1, dx = dir == 0 ? 1 : 0, comp = dx ? 0 : 1;
        for (int y = 0; y + dy < gh; y++)
            for (int x = 0; x + dx < gw; x++) {
                size_t a = (size_t)y * (size_t)gw + (size_t)x, b = (size_t)(y + dy) * (size_t)gw + (size_t)(x + dx);
                int real = 0;
                if (ids[a] < 0 || ids[b] < 0) continue;
                real = ff->real[a] && ff->real[b];
                ea[ne] = ids[a];
                eb[ne] = ids[b];
                tu[ne] = real ? (gu[a * 2 + (size_t)comp] + gu[b * 2 + (size_t)comp]) / 2.0 : (dx ? (double)cell : 0.0);
                tv[ne] = real ? (gv[a * 2 + (size_t)comp] + gv[b * 2 + (size_t)comp]) / 2.0 : (dx ? 0.0 : (double)cell);
                we[ne] = real ? 1.0 : FIBER_FLATTEN_GAP_WEIGHT;
                ne++;
            }
    }
    st->edges = ne;
    nt = 3 * ne + n;
    rows = (int *)malloc(nt * sizeof *rows);
    cols = (int *)malloc(nt * sizeof *cols);
    vals = (double *)malloc(nt * sizeof *vals);
    rhs = (double *)calloc(n * 2, sizeof *rhs);
    sol = (double *)calloc(n * 2, sizeof *sol);
    if (rows == NULL || cols == NULL || vals == NULL || rhs == NULL || sol == NULL) goto done;
    nt = 0;
    for (size_t e = 0; e < ne; e++) {
        int a = ea[e], b = eb[e];
        rows[nt] = a; cols[nt] = a; vals[nt++] = we[e];
        rows[nt] = b; cols[nt] = b; vals[nt++] = we[e];
        rows[nt] = a > b ? a : b; cols[nt] = a > b ? b : a; vals[nt++] = -we[e];
        rhs[a] -= we[e] * tu[e];       rhs[b] += we[e] * tu[e];
        rhs[n + (size_t)a] -= we[e] * tv[e];  rhs[n + (size_t)b] += we[e] * tv[e];
    }
    for (int gi = 0; gi < gh; gi++)
        for (int gj = 0; gj < gw; gj++) {
            size_t i = (size_t)gi * (size_t)gw + (size_t)gj;
            if (ids[i] < 0) continue;
            rows[nt] = ids[i]; cols[nt] = ids[i]; vals[nt++] = FIBER_FLATTEN_EPS;
            rhs[ids[i]] += FIBER_FLATTEN_EPS * (double)(ff->gx0 + gj * cell);
            rhs[n + (size_t)ids[i]] += FIBER_FLATTEN_EPS * (double)(ff->gy0 + gi * cell);
        }
    if (n == 0 || Sparse_factor_spd((int)n, (int)nt, rows, cols, vals, &fac) != 0) goto done;
    if (Sparse_factor_solve_multi(fac, rhs, sol, 2) != 0) goto done;
    for (size_t i = 0; i < nn; i++) {
        ff->U[i] = ids[i] >= 0 ? sol[ids[i]] : NAN;
        ff->V[i] = ids[i] >= 0 ? sol[n + (size_t)ids[i]] : NAN;
    }
    for (size_t e = 0; e < ne; e++) {
        double ru = 0.0, rv = 0.0;
        if (we[e] != 1.0) continue;
        ru = sol[eb[e]] - sol[ea[e]] - tu[e];
        rv = sol[n + (size_t)eb[e]] - sol[n + (size_t)ea[e]] - tv[e];
        su += ru * ru;
        sv += rv * rv;
        nreal_e++;
    }
    st->rms_residual_u = nreal_e ? sqrt(su / (double)nreal_e) : 0.0;
    st->rms_residual_v = nreal_e ? sqrt(sv / (double)nreal_e) : 0.0;
    qsort(lh, nl, sizeof *lh, ff_cmp_double);
    qsort(lv, nl, sizeof *lv, ff_cmp_double);
    st->L_h_p50 = ff_pct_sorted(lh, nl, 50.0);
    st->L_v_p50 = ff_pct_sorted(lv, nl, 50.0);
    rc = 0;
done:
    Sparse_factor_free(&fac);
    free(dil); free(clo); free(ids); free(gu); free(gv); free(lh); free(lv);
    free(ea); free(eb); free(tu); free(tv); free(we);
    free(rows); free(cols); free(vals); free(rhs); free(sol);
    return rc;
}

/* ------------------------------------------------------------------------------------------------------ audit */

/* the 2-D rotation (rad) taking the centred a's to the centred b's, and the residuals' magnitudes into res[] */
static double ff_rigid(const double *ax, const double *ay, const double *bx, const double *by, size_t n, double *res)
{
    double cax = 0.0, cay = 0.0, cbx = 0.0, cby = 0.0, sxx = 0.0, sxy = 0.0, th = 0.0, c = 0.0, s = 0.0;
    for (size_t i = 0; i < n; i++) { cax += ax[i]; cay += ay[i]; cbx += bx[i]; cby += by[i]; }
    cax /= (double)n; cay /= (double)n; cbx /= (double)n; cby /= (double)n;
    for (size_t i = 0; i < n; i++) {
        double px = ax[i] - cax, py = ay[i] - cay, qx = bx[i] - cbx, qy = by[i] - cby;
        sxx += px * qx + py * qy;
        sxy += px * qy - py * qx;
    }
    th = atan2(sxy, sxx);
    c = cos(th);
    s = sin(th);
    if (res != NULL)
        for (size_t i = 0; i < n; i++) {
            double px = ax[i] - cax, py = ay[i] - cay;
            double rx = (bx[i] - cbx) - (c * px - s * py), ry = (by[i] - cby) - (s * px + c * py);
            res[i] = sqrt(rx * rx + ry * ry);
        }
    return th;
}

typedef struct { double x, a; } FFPair;

static int ff_cmp_pair(const void *p, const void *q)
{
    double x = ((const FFPair *)p)->x, y = ((const FFPair *)q)->x;
    return x < y ? -1 : x > y ? 1 : 0;
}

/* The fibre-climb profile: bins of `step` px along x; per bin with >= MIN nodes the circular median angle; empty
 * bins interpolated; height = cumsum(tan(angle) * step).  Fills range / bow / end / width and the p10/p90 angles. */
static int ff_climb(FFPair *pp, size_t n, double step, double *range, double *bow, double *end, double *width,
                    double *p10, double *p90, FILE *csv, const char *kind)
{
    double lo = 0.0, hi = 0.0, *med = NULL, *h = NULL, *xb = NULL, *buf = NULL, *ok = NULL;
    size_t nb = 0, j = 0, k = 0, nok = 0;
    double sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0, slope = 0.0, icpt = 0.0, rmin = 0.0, rmax = 0.0, hmin = 0.0, hmax = 0.0;
    if (n == 0) return -1;
    qsort(pp, n, sizeof *pp, ff_cmp_pair);
    lo = floor(pp[0].x / step);
    hi = ceil(pp[n - 1].x / step);
    nb = (size_t)(hi - lo) + 1;
    med = (double *)malloc(nb * sizeof *med);
    h = (double *)malloc(nb * sizeof *h);
    xb = (double *)malloc(nb * sizeof *xb);
    buf = (double *)malloc(n * sizeof *buf);
    ok = (double *)malloc(nb * sizeof *ok);
    if (med == NULL || h == NULL || xb == NULL || buf == NULL || ok == NULL) {
        free(med); free(h); free(xb); free(buf); free(ok);
        return -1;
    }
    for (size_t b = 0; b < nb; b++) { med[b] = NAN; xb[b] = (lo + (double)b) * step; }
    k = 0;
    for (size_t b = 0; b < nb; b++) {
        size_t m = 0;
        double right = (lo + (double)b + 1.0) * step;
        while (k < n && (b + 1 == nb || pp[k].x < right)) buf[m++] = pp[k++].a;
        if (m >= FIBER_FLATTEN_PROFILE_MIN) med[b] = ff_circ_median_deg(buf, m);
    }
    for (size_t b = 0; b < nb; b++) if (!isnan(med[b])) ok[nok++] = med[b];
    if (nok == 0) { free(med); free(h); free(xb); free(buf); free(ok); return -1; }
    qsort(ok, nok, sizeof *ok, ff_cmp_double);
    *p10 = ff_pct_sorted(ok, nok, 10.0);
    *p90 = ff_pct_sorted(ok, nok, 90.0);
    /* interpolate the empty bins (np.interp: clamp at the ends) */
    for (size_t b = 0; b < nb; b++) {
        size_t l = b, r = b;
        if (!isnan(med[b])) continue;
        while (l > 0 && isnan(med[l])) l--;
        while (r + 1 < nb && isnan(med[r])) r++;
        if (isnan(med[l]) && isnan(med[r])) continue;
        if (isnan(med[l])) ok[0] = med[r];
        else if (isnan(med[r])) ok[0] = med[l];
        else ok[0] = med[l] + (med[r] - med[l]) * (double)(b - l) / (double)(r - l);
        h[b] = ok[0];                         /* staged: a filled bin must not feed its neighbours' search */
    }
    for (size_t b = 0; b < nb; b++) if (isnan(med[b])) med[b] = h[b];
    for (j = 0; j < nb; j++) {
        double t = isnan(med[j]) ? 0.0 : tan(med[j] * M_PI / 180.0) * step;
        h[j] = (j ? h[j - 1] : 0.0) + t;
    }
    for (j = 0; j < nb; j++) { sx += xb[j]; sy += h[j]; sxx += xb[j] * xb[j]; sxy += xb[j] * h[j]; }
    slope = (nb > 1) ? ((double)nb * sxy - sx * sy) / ((double)nb * sxx - sx * sx) : 0.0;
    icpt = (sy - slope * sx) / (double)nb;
    rmin = rmax = h[0] - (slope * xb[0] + icpt);
    hmin = hmax = h[0];
    for (j = 0; j < nb; j++) {
        double r = h[j] - (slope * xb[j] + icpt);
        if (r < rmin) rmin = r;
        if (r > rmax) rmax = r;
        if (h[j] < hmin) hmin = h[j];
        if (h[j] > hmax) hmax = h[j];
        if (csv != NULL) fprintf(csv, "%s,%.1f,%.3f,%.4f\n", kind, xb[j], h[j], med[j]);
    }
    *range = hmax - hmin;
    *bow = rmax - rmin;
    *end = h[nb - 1] - h[0];
    *width = xb[nb - 1] - xb[0];
    free(med); free(h); free(xb); free(buf); free(ok);
    return 0;
}

int FiberFrame_audit(Arena_T arena, const FiberFrame *ff, FiberFrameStats *st, const char *profile_csv)
{
    const int gh = ff->gh, gw = ff->gw, cell = ff->cell;
    const size_t nn = (size_t)gh * (size_t)gw;
    double *Uc = NULL, *Ur = NULL, *Vc = NULL, *Vr = NULL, *tmp = NULL;
    double *ax = NULL, *ay = NULL, *bx = NULL, *by = NULL, *res = NULL, *tres = NULL, *trot = NULL;
    double *ho = NULL, *hn = NULL, *area = NULL, wmed = 0.0, sho = 0.0, shn = 0.0;
    FFPair *po = NULL, *pn = NULL;
    size_t nok = 0, nho = 0, nhn = 0, narea = 0, ntres = 0, ntrot = 0, npo = 0, npn = 0;
    FILE *csv = NULL;
    int rc = -1;
    (void)arena;
    Uc = (double *)malloc(nn * sizeof(double)); Ur = (double *)malloc(nn * sizeof(double));
    Vc = (double *)malloc(nn * sizeof(double)); Vr = (double *)malloc(nn * sizeof(double));
    tmp = (double *)malloc(nn * sizeof(double));
    ax = (double *)malloc(nn * sizeof(double)); ay = (double *)malloc(nn * sizeof(double));
    bx = (double *)malloc(nn * sizeof(double)); by = (double *)malloc(nn * sizeof(double));
    res = (double *)malloc(nn * sizeof(double)); tres = (double *)malloc(nn * sizeof(double));
    trot = (double *)malloc(nn * sizeof(double));
    ho = (double *)malloc(nn * sizeof(double)); hn = (double *)malloc(nn * sizeof(double));
    area = (double *)malloc(nn * sizeof(double));
    po = (FFPair *)malloc(nn * sizeof *po); pn = (FFPair *)malloc(nn * sizeof *pn);
    if (Uc == NULL || Ur == NULL || Vc == NULL || Vr == NULL || tmp == NULL || ax == NULL || ay == NULL || bx == NULL ||
        by == NULL || res == NULL || tres == NULL || trot == NULL || ho == NULL || hn == NULL || area == NULL ||
        po == NULL || pn == NULL) goto done;
    ff_gradient(ff->U, Uc, gh, gw, 0, cell);
    ff_gradient(ff->U, Ur, gh, gw, 1, cell);
    ff_gradient(ff->V, Vc, gh, gw, 0, cell);
    ff_gradient(ff->V, Vr, gh, gw, 1, cell);
    /* the evidence threshold: FIBER_FLATTEN_EVIDENCE_REL x the median smoothed weight on real nodes */
    for (size_t i = 0; i < nn; i++) if (ff->real[i]) tmp[nok++] = (double)ff->weight[i];
    qsort(tmp, nok, sizeof *tmp, ff_cmp_double);
    wmed = ff_pct_sorted(tmp, nok, 50.0);
    nok = 0;
    for (int gi = 0; gi < gh; gi++)
        for (int gj = 0; gj < gw; gj++) {
            size_t i = (size_t)gi * (size_t)gw + (size_t)gj;
            int ev = 0;
            double det = 0.0, c = 0.0, s = 0.0, an = 0.0, tdeg = 0.0;
            if (!ff->real[i] || isnan(ff->U[i]) || isnan(ff->V[i])) continue;
            ax[nok] = (double)(ff->gx0 + gj * cell);
            ay[nok] = (double)(ff->gy0 + gi * cell);
            bx[nok] = ff->U[i];
            by[nok] = ff->V[i];
            nok++;
            det = Uc[i] * Vr[i] - Ur[i] * Vc[i];
            if (!isnan(det)) {
                area[narea++] = det;
                if (det <= 0.0) st->flipped++;
            }
            ev = (double)ff->weight[i] > FIBER_FLATTEN_EVIDENCE_REL * wmed;
            if (!ev) continue;
            tdeg = (double)ff->theta[i] * 180.0 / M_PI;
            ho[nho++] = fabs(tdeg);
            sho += tdeg;
            po[npo].x = (double)(ff->gx0 + gj * cell);
            po[npo++].a = tdeg;
            c = cos((double)ff->theta[i]);
            s = sin((double)ff->theta[i]);
            an = atan2(c * Vc[i] + s * Vr[i], c * Uc[i] + s * Ur[i]) * 180.0 / M_PI;
            if (isnan(an)) continue;
            hn[nhn++] = fabs(an);
            shn += an;
            pn[npn].x = ff->U[i];
            pn[npn++].a = an;
        }
    if (nok < 3) goto done;
    st->evidence_share = (double)nho / (double)nok;
    st->H_old_signed_mean = nho ? sho / (double)nho : NAN;
    st->H_new_signed_mean = nhn ? shn / (double)nhn : NAN;
    st->H_old_abs = ff_quant(ho, nho);
    st->H_new_abs = ff_quant(hn, nhn);
    qsort(area, narea, sizeof *area, ff_cmp_double);
    st->area_p1 = ff_pct_sorted(area, narea, 1.0);
    st->area_p10 = ff_pct_sorted(area, narea, 10.0);
    st->area_p50 = ff_pct_sorted(area, narea, 50.0);
    st->area_p90 = ff_pct_sorted(area, narea, 90.0);
    st->area_p99 = ff_pct_sorted(area, narea, 99.0);
    /* the rigid part and what is left */
    {
        double th = ff_rigid(ax, ay, bx, by, nok, res);
        st->rotation_deg = th * 180.0 / M_PI;
        memcpy(tmp, res, nok * sizeof *tmp);
        st->nonrigid_px = ff_quant(tmp, nok);
        st->nonrigid_max = tmp[nok - 1];
        /* FIBER_FLATTEN_TILE_PX patches: nodes sorted by tile key, a rigid fit per tile with >= 300 nodes */
        {
            FFPair *tk = (FFPair *)malloc(nok * sizeof *tk);
            double *qx = NULL, *qy = NULL, *rx = NULL, *ry = NULL;
            size_t a = 0;
            if (tk == NULL) goto done;
            qx = ax; qy = ay; rx = bx; ry = by;
            for (size_t i = 0; i < nok; i++) {
                tk[i].x = floor(qy[i] / FIBER_FLATTEN_TILE_PX) * 1e6 + floor(qx[i] / FIBER_FLATTEN_TILE_PX);
                tk[i].a = (double)i;
            }
            qsort(tk, nok, sizeof *tk, ff_cmp_pair);
            while (a < nok) {
                size_t b = a, m = 0;
                while (b < nok && tk[b].x == tk[a].x) b++;
                m = b - a;
                if (m >= 300) {
                    double *t1 = (double *)malloc(m * 5 * sizeof(double));
                    if (t1 != NULL) {
                        double tth = 0.0, rel = 0.0;
                        for (size_t k = 0; k < m; k++) {
                            size_t id = (size_t)tk[a + k].a;
                            t1[k] = qx[id]; t1[m + k] = qy[id]; t1[2 * m + k] = rx[id]; t1[3 * m + k] = ry[id];
                        }
                        tth = ff_rigid(t1, t1 + m, t1 + 2 * m, t1 + 3 * m, m, t1 + 4 * m);
                        for (size_t k = 0; k < m; k++) tres[ntres++] = t1[4 * m + k];
                        rel = fmod(tth * 180.0 / M_PI - st->rotation_deg + 90.0, 180.0);
                        if (rel < 0.0) rel += 180.0;
                        trot[ntrot++] = fabs(rel - 90.0);
                        free(t1);
                    }
                }
                a = b;
            }
            free(tk);
            st->tiles = ntrot;
            st->tile_nonrigid_px = ff_quant(tres, ntres);
            st->tile_rotation_deg = ff_quant(trot, ntrot);
        }
    }
    /* the affine fit (col, row, 1) -> (u', v'): normal equations, 3x3 */
    {
        double M[9] = { 0 }, Ru[3] = { 0 }, Rv[3] = { 0 }, inv[9], dt = 0.0;
        for (size_t i = 0; i < nok; i++) {
            double r3[3] = { ax[i], ay[i], 1.0 };
            for (int p = 0; p < 3; p++) {
                for (int q = 0; q < 3; q++) M[p * 3 + q] += r3[p] * r3[q];
                Ru[p] += r3[p] * bx[i];
                Rv[p] += r3[p] * by[i];
            }
        }
        inv[0] = M[4] * M[8] - M[5] * M[7];
        inv[1] = M[2] * M[7] - M[1] * M[8];
        inv[2] = M[1] * M[5] - M[2] * M[4];
        inv[3] = M[5] * M[6] - M[3] * M[8];
        inv[4] = M[0] * M[8] - M[2] * M[6];
        inv[5] = M[2] * M[3] - M[0] * M[5];
        inv[6] = M[3] * M[7] - M[4] * M[6];
        inv[7] = M[1] * M[6] - M[0] * M[7];
        inv[8] = M[0] * M[4] - M[1] * M[3];
        dt = M[0] * inv[0] + M[1] * inv[3] + M[2] * inv[6];
        if (fabs(dt) < 1e-300) goto done;
        for (int p = 0; p < 3; p++) {
            double cu = 0.0, cv = 0.0;
            for (int q = 0; q < 3; q++) { cu += inv[p * 3 + q] * Ru[q]; cv += inv[p * 3 + q] * Rv[q]; }
            cu /= dt; cv /= dt;
            if (p < 2) { st->A[p] = cu; st->A[2 + p] = cv; } else { st->b[0] = cu; st->b[1] = cv; }
        }
    }
    if (profile_csv != NULL) {
        csv = fopen(profile_csv, "w");
        if (csv != NULL) fprintf(csv, "kind,x_px,height_px,angle_deg\n");
    }
    if (npo) ff_climb(po, npo, (double)cell, &st->climb_old_range, &st->climb_old_bow, &st->climb_old_end,
                      &st->climb_old_width, &st->col_angle_old_p10, &st->col_angle_old_p90, csv, "old");
    if (npn) ff_climb(pn, npn, (double)cell, &st->climb_new_range, &st->climb_new_bow, &st->climb_new_end,
                      &st->climb_new_width, &st->col_angle_new_p10, &st->col_angle_new_p90, csv, "new");
    if (csv != NULL) fclose(csv);
    rc = 0;
done:
    free(Uc); free(Ur); free(Vc); free(Vr); free(tmp); free(ax); free(ay); free(bx); free(by); free(res);
    free(tres); free(trot); free(ho); free(hn); free(area); free(po); free(pn);
    return rc;
}

/* ---------------------------------------------------------------------------------------------------- inverse */

static int ff_grid_bilinear(const double *G, int gh, int gw, double gy, double gx, double *out)
{
    int y0 = 0, x0 = 0;
    double fy = 0.0, fx = 0.0;
    if (!(gy >= 0.0) || !(gx >= 0.0) || gy > (double)(gh - 1) || gx > (double)(gw - 1)) return -1;
    if (gy > (double)gh - 1.000001) gy = (double)gh - 1.000001;
    if (gx > (double)gw - 1.000001) gx = (double)gw - 1.000001;
    y0 = (int)floor(gy);
    x0 = (int)floor(gx);
    fy = gy - y0;
    fx = gx - x0;
    *out = G[(size_t)y0 * (size_t)gw + (size_t)x0] * (1 - fy) * (1 - fx)
         + G[(size_t)(y0 + 1) * (size_t)gw + (size_t)x0] * fy * (1 - fx)
         + G[(size_t)y0 * (size_t)gw + (size_t)(x0 + 1)] * (1 - fy) * fx
         + G[(size_t)(y0 + 1) * (size_t)gw + (size_t)(x0 + 1)] * fy * fx;
    return isnan(*out) ? -1 : 0;
}

int FiberFrame_inverse(Arena_T arena, const FiberFrame *ff, const FiberFrameStats *st, int threads,
                       int *Hn, int *Wn, double *u_min, double *v_min, float **old_rc)
{
    const int gh = ff->gh, gw = ff->gw, cell = ff->cell;
    const size_t nn = (size_t)gh * (size_t)gw;
    double *Uf = NULL, *Vf = NULL, umn = INFINITY, umx = -INFINITY, vmn = INFINITY, vmx = -INFINITY, Li[4], dt = 0.0;
    long *src = NULL;
    size_t *q = NULL, head = 0, tail = 0;
    int hn = 0, wn = 0, pr_row = 0;
    float *orc = NULL;
    if (threads <= 0) threads = 1;
    Uf = (double *)malloc(nn * sizeof *Uf);
    Vf = (double *)malloc(nn * sizeof *Vf);
    src = (long *)malloc(nn * sizeof *src);
    q = (size_t *)malloc(nn * sizeof *q);
    if (Uf == NULL || Vf == NULL || src == NULL || q == NULL) { free(Uf); free(Vf); free(src); free(q); return -1; }
    /* past the solved nodes: the nearest solved node's displacement (multi-source BFS, 8-connected) */
    for (size_t i = 0; i < nn; i++) {
        src[i] = -1;
        if (ff->closed[i] && !isnan(ff->U[i]) && !isnan(ff->V[i])) { src[i] = (long)i; q[tail++] = i; }
    }
    while (head < tail) {
        size_t i = q[head++];
        int y = (int)(i / (size_t)gw), x = (int)(i % (size_t)gw);
        for (int dy = -1; dy <= 1; dy++)
            for (int dx = -1; dx <= 1; dx++) {
                size_t j = 0;
                if ((dy == 0 && dx == 0) || y + dy < 0 || y + dy >= gh || x + dx < 0 || x + dx >= gw) continue;
                j = (size_t)(y + dy) * (size_t)gw + (size_t)(x + dx);
                if (src[j] < 0) { src[j] = src[i]; q[tail++] = j; }
            }
    }
    for (int gi = 0; gi < gh; gi++)
        for (int gj = 0; gj < gw; gj++) {
            size_t i = (size_t)gi * (size_t)gw + (size_t)gj;
            double gx = (double)(ff->gx0 + gj * cell), gy = (double)(ff->gy0 + gi * cell);
            if (src[i] < 0) { Uf[i] = Vf[i] = NAN; continue; }
            {
                size_t s = (size_t)src[i];
                double sgx = (double)(ff->gx0 + (int)(s % (size_t)gw) * cell);
                double sgy = (double)(ff->gy0 + (int)(s / (size_t)gw) * cell);
                Uf[i] = gx + (ff->U[s] - sgx);
                Vf[i] = gy + (ff->V[s] - sgy);
            }
            if (ff->closed[i] && !isnan(ff->U[i])) {
                if (ff->U[i] < umn) umn = ff->U[i];
                if (ff->U[i] > umx) umx = ff->U[i];
                if (ff->V[i] < vmn) vmn = ff->V[i];
                if (ff->V[i] > vmx) vmx = ff->V[i];
            }
        }
    free(src);
    free(q);
    if (!isfinite(umn) || !isfinite(vmn)) { free(Uf); free(Vf); return -1; }
    umn = floor(umn) - FIBER_FLATTEN_PAD_CELLS * cell;
    vmn = floor(vmn) - FIBER_FLATTEN_PAD_CELLS * cell;
    umx = ceil(umx) + FIBER_FLATTEN_PAD_CELLS * cell;
    vmx = ceil(vmx) + FIBER_FLATTEN_PAD_CELLS * cell;
    hn = (int)(vmx - vmn) + 1;
    wn = (int)(umx - umn) + 1;
    dt = st->A[0] * st->A[3] - st->A[1] * st->A[2];
    if (fabs(dt) < 1e-12 || hn <= 0 || wn <= 0) { free(Uf); free(Vf); return -1; }
    Li[0] = st->A[3] / dt;  Li[1] = -st->A[1] / dt;
    Li[2] = -st->A[2] / dt; Li[3] = st->A[0] / dt;
    orc = (float *)ARENA_ALLOC(arena, (size_t)hn * (size_t)wn * 2 * sizeof *orc);
    if (orc == NULL) { free(Uf); free(Vf); return -1; }
#pragma omp parallel for schedule(dynamic, 8) num_threads(threads)
    for (pr_row = 0; pr_row < hn; pr_row++) {
        const int r = pr_row;
        for (int c = 0; c < wn; c++) {
            double tu = umn + c, tv = vmn + r;
            double pc = Li[0] * (tu - st->b[0]) + Li[1] * (tv - st->b[1]);
            double pr = Li[2] * (tu - st->b[0]) + Li[3] * (tv - st->b[1]);
            double res = INFINITY;
            int it = 0, out = 0;
            for (it = 0; it < FIBER_FLATTEN_INV_ITERS; it++) {
                double fu = 0.0, fv = 0.0, eu = 0.0, ev = 0.0, sc = 0.0, sr = 0.0;
                double gy = (pr - ff->gy0) / cell, gx = (pc - ff->gx0) / cell;
                if (ff_grid_bilinear(Uf, gh, gw, gy, gx, &fu) != 0 || ff_grid_bilinear(Vf, gh, gw, gy, gx, &fv) != 0) {
                    out = 1;
                    break;
                }
                eu = fu - tu;
                ev = fv - tv;
                res = sqrt(eu * eu + ev * ev);
                sc = Li[0] * eu + Li[1] * ev;
                sr = Li[2] * eu + Li[3] * ev;
                pc -= sc;
                pr -= sr;
                if (fabs(sc) < FIBER_FLATTEN_INV_STEP && fabs(sr) < FIBER_FLATTEN_INV_STEP) break;
            }
            {
                size_t o = ((size_t)r * (size_t)wn + (size_t)c) * 2;
                if (out || !(res < FIBER_FLATTEN_INV_TOL)) { orc[o] = NAN; orc[o + 1] = NAN; }
                else { orc[o] = (float)pr; orc[o + 1] = (float)pc; }
            }
        }
    }
    free(Uf);
    free(Vf);
    *Hn = hn; *Wn = wn; *u_min = umn; *v_min = vmn; *old_rc = orc;
    return 0;
}

int FiberFrame_resample(Arena_T arena, const float *X, const uint8_t *valid, int H, int W, int x0, int y0, int x1, int y1,
                        const float *old_rc, int Hn, int Wn, float **Xn, uint8_t **vn, size_t *n_valid)
{
    size_t np = (size_t)Hn * (size_t)Wn, cnt = 0;
    int rs_row = 0;
    float *xo = (float *)ARENA_ALLOC(arena, np * 3 * sizeof *xo);
    uint8_t *v0 = (uint8_t *)ARENA_ALLOC(arena, np), *v1 = (uint8_t *)ARENA_ALLOC(arena, np);
    if (xo == NULL || v0 == NULL || v1 == NULL) return -1;
    if (x1 > W) x1 = W;
    if (y1 > H) y1 = H;
#pragma omp parallel for schedule(static)
    for (rs_row = 0; rs_row < Hn; rs_row++)
        for (int c = 0; c < Wn; c++) {
            size_t o = (size_t)rs_row * (size_t)Wn + (size_t)c;
            double pr = (double)old_rc[o * 2], pc = (double)old_rc[o * 2 + 1];
            int ok = 0;
            v0[o] = 0;
            xo[o * 3] = xo[o * 3 + 1] = xo[o * 3 + 2] = -1.0f;
            if (!isnan(pr) && !isnan(pc) && pr >= y0 && pc >= x0 && pr <= y1 - 1 && pc <= x1 - 1) {
                int r0 = (int)floor(pr), c0 = (int)floor(pc), r1 = 0, c1 = 0;
                double fy = 0.0, fx = 0.0;
                if (r0 > y1 - 2) r0 = y1 - 2;
                if (c0 > x1 - 2) c0 = x1 - 2;
                r1 = r0 + 1; c1 = c0 + 1;
                fy = pr - r0; fx = pc - c0;
                ok = r0 >= y0 && c0 >= x0 &&
                     valid[(size_t)r0 * (size_t)W + (size_t)c0] && valid[(size_t)r1 * (size_t)W + (size_t)c0] &&
                     valid[(size_t)r0 * (size_t)W + (size_t)c1] && valid[(size_t)r1 * (size_t)W + (size_t)c1];
                if (ok) {
                    for (int k = 0; k < 3; k++) {
                        double a = X[((size_t)r0 * (size_t)W + (size_t)c0) * 3 + (size_t)k];
                        double b = X[((size_t)r1 * (size_t)W + (size_t)c0) * 3 + (size_t)k];
                        double cc = X[((size_t)r0 * (size_t)W + (size_t)c1) * 3 + (size_t)k];
                        double d = X[((size_t)r1 * (size_t)W + (size_t)c1) * 3 + (size_t)k];
                        xo[o * 3 + (size_t)k] = (float)(a * (1 - fy) * (1 - fx) + b * fy * (1 - fx) + cc * (1 - fy) * fx + d * fy * fx);
                    }
                    v0[o] = 1;
                }
            }
        }
    /* erode the valid mask by one pixel (3 x 3, outside = invalid) */
    for (int r = 0; r < Hn; r++)
        for (int c = 0; c < Wn; c++) {
            size_t o = (size_t)r * (size_t)Wn + (size_t)c;
            int keep = v0[o];
            for (int dy = -1; dy <= 1 && keep; dy++)
                for (int dx = -1; dx <= 1; dx++) {
                    int rr = r + dy, cc = c + dx;
                    if (rr < 0 || rr >= Hn || cc < 0 || cc >= Wn || !v0[(size_t)rr * (size_t)Wn + (size_t)cc]) { keep = 0; break; }
                }
            v1[o] = (uint8_t)keep;
            if (!keep) xo[o * 3] = xo[o * 3 + 1] = xo[o * 3 + 2] = -1.0f;
            else cnt++;
        }
    *Xn = xo;
    *vn = v1;
    *n_valid = cnt;
    return 0;
}

/* --------------------------------------------------------------------------------------------------- selftest */

#define FFT_CHECK(cond, msg) do { if (!(cond)) { fprintf(stderr, "[selftest] fiber_frame FAIL: %s\n", msg); fails++; } } while (0)

static int ff_write_zarray(const char *dir, const long shape[3], int chunk)
{
    char path[1200];
    FILE *f = NULL;
    snprintf(path, sizeof path, "%s/0/.zarray", dir);
    f = fopen(path, "wb");
    if (f == NULL) return -1;
    fprintf(f, "{\"chunks\": [%d, %d, %d], \"compressor\": null, \"dimension_separator\": \"/\", \"dtype\": \"|u1\", "
               "\"fill_value\": 0, \"filters\": null, \"order\": \"C\", \"shape\": [%ld, %ld, %ld], \"zarr_format\": 2}",
            chunk, chunk, chunk, shape[0], shape[1], shape[2]);
    fclose(f);
    return 0;
}

#include "fiber_frame_map_audit.inc"
#include "fiber_frame_metric_test.inc"
#include "fiber_frame_branch_audit.inc"

int FiberFrame_selftest(void)
{
    int fails = ff_metric_selftest() + ff_map_audit_selftest() + ff_branch_selftest();
    char root[512], path[1200];
    /* helpers */
    {
        uint8_t m[7 * 7], d[7 * 7], e[7 * 7];
        double a[4] = { 89.0, -89.0, 88.0, -88.0 }, g[5 * 5], t[5 * 5];
        memset(m, 0, sizeof m);
        for (int y = 1; y < 6; y++) for (int x = 1; x < 6; x++) m[y * 7 + x] = 1;
        m[3 * 7 + 3] = 0;                                  /* an enclosed hole */
        ff_fill_holes(m, 7, 7);
        FFT_CHECK(m[3 * 7 + 3] == 1 && m[0] == 0, "fill_holes fills the enclosed hole only");
        memset(m, 0, sizeof m);
        for (int y = 0; y < 7; y++) for (int x = 0; x < 7; x++) m[y * 7 + x] = (uint8_t)(x != 3);   /* a 1-cell crack */
        ff_morph(m, d, 7, 7, 1, 1);
        ff_morph(d, e, 7, 7, 1, 0);
        FFT_CHECK(e[3 * 7 + 3] == 1, "closing fills a one-cell crack");
        FFT_CHECK(fabs(fabs(ff_circ_median_deg(a, 4)) - 90.0) < 2.0, "the axial median wraps at +-90 deg");
        for (int i = 0; i < 25; i++) g[i] = 3.0;
        ff_gauss(g, 5, 5, 1.5, t);
        FFT_CHECK(fabs(g[12] - 3.0) < 1e-9 && fabs(g[0] - 3.0) < 1e-9, "the Gaussian keeps a constant");
    }
    /* a planted bowed fibre field on a tilted plane */
    if (ves_selftest_dir(root, sizeof root, "fiber_frame") != 0) return 1;
    {
        const int H = 96, W = 320, CH = 32;
        const double kappa = 1.5625e-3;                  /* rad per px: a ~20 px climb over 320 px */
        double e1[3] = { 1.0, 0.1, 0.0 }, e2[3] = { 0.0, 0.05, 1.0 }, n3[3], O[3] = { 40.0, 60.0, 20.0 }, l = 0.0, dd = 0.0;
        long shape[3] = { 160, 128, 384 };
        const char *names[3] = { "presence", "nx", "ny" };
        Arena_T ar = Arena_new();
        float *X = (float *)calloc((size_t)H * W * 3, sizeof(float));
        uint8_t *valid = (uint8_t *)calloc((size_t)H * W, 1);
        uint8_t *chunk[3];
        FiberVolume_T fv = NULL;
        FiberFrame ff;
        FiberFrameStats st;
        l = sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);
        for (int k = 0; k < 3; k++) e1[k] /= l;
        dd = e2[0] * e1[0] + e2[1] * e1[1] + e2[2] * e1[2];
        for (int k = 0; k < 3; k++) e2[k] -= dd * e1[k];
        l = sqrt(e2[0] * e2[0] + e2[1] * e2[1] + e2[2] * e2[2]);
        for (int k = 0; k < 3; k++) e2[k] /= l;
        n3[0] = e1[1] * e2[2] - e1[2] * e2[1];
        n3[1] = e1[2] * e2[0] - e1[0] * e2[2];
        n3[2] = e1[0] * e2[1] - e1[1] * e2[0];
        for (int r = 0; r < H; r++)
            for (int c = 0; c < W; c++)
                for (int k = 0; k < 3; k++) {
                    X[((size_t)r * W + (size_t)c) * 3 + (size_t)k] = (float)(O[k] + c * e1[k] + r * e2[k]);
                    valid[(size_t)r * W + (size_t)c] = 1;
                }
        for (int a = 0; a < 3; a++) {
            snprintf(path, sizeof path, "%s/f", root); ves_mkdir(path);
            snprintf(path, sizeof path, "%s/f/%s", root, names[a]); ves_mkdir(path);
            snprintf(path, sizeof path, "%s/f/%s/0", root, names[a]); ves_mkdir(path);
            snprintf(path, sizeof path, "%s/f/%s", root, names[a]);
            ff_write_zarray(path, shape, CH);
            chunk[a] = (uint8_t *)malloc((size_t)CH * CH * CH);
        }
        for (long cz = 0; cz < shape[0] / CH; cz++)
            for (long cy = 0; cy < shape[1] / CH; cy++)
                for (long cx = 0; cx < shape[2] / CH; cx++) {
                    int any = 0;
                    for (int a = 0; a < 3; a++) memset(chunk[a], 0, (size_t)CH * CH * CH);
                    for (int z = 0; z < CH; z++)
                        for (int y = 0; y < CH; y++)
                            for (int x = 0; x < CH; x++) {
                                double q[3] = { (double)(cx * CH + x), (double)(cy * CH + y), (double)(cz * CH + z) };
                                double rel[3] = { q[0] - O[0], q[1] - O[1], q[2] - O[2] };
                                double dist = rel[0] * n3[0] + rel[1] * n3[1] + rel[2] * n3[2];
                                double u = rel[0] * e1[0] + rel[1] * e1[1] + rel[2] * e1[2];
                                double v = rel[0] * e2[0] + rel[1] * e2[1] + rel[2] * e2[2];
                                double al = 0.0, d3[3];
                                size_t o = ((size_t)z * CH + (size_t)y) * CH + (size_t)x;
                                if (fabs(dist) > 4.0 || u < -4.0 || u > W + 4.0 || v < -4.0 || v > H + 4.0) continue;
                                al = kappa * (u - W / 2.0);
                                for (int k = 0; k < 3; k++) d3[k] = cos(al) * e1[k] + sin(al) * e2[k];
                                if (d3[2] < 0.0) for (int k = 0; k < 3; k++) d3[k] = -d3[k];
                                chunk[0][o] = 255;
                                chunk[1][o] = (uint8_t)lround(127.0 * d3[0] + 128.0);
                                chunk[2][o] = (uint8_t)lround(127.0 * d3[1] + 128.0);
                                any = 1;
                            }
                    if (!any) continue;
                    for (int a = 0; a < 3; a++) {
                        FILE *f = NULL;
                        snprintf(path, sizeof path, "%s/f/%s/0/%ld", root, names[a], cz); ves_mkdir(path);
                        snprintf(path, sizeof path, "%s/f/%s/0/%ld/%ld", root, names[a], cz, cy); ves_mkdir(path);
                        snprintf(path, sizeof path, "%s/f/%s/0/%ld/%ld/%ld", root, names[a], cz, cy, cx);
                        f = fopen(path, "wb");
                        if (f != NULL) { fwrite(chunk[a], 1, (size_t)CH * CH * CH, f); fclose(f); }
                    }
                }
        snprintf(path, sizeof path, "%s/f", root);
        FFT_CHECK(FiberVolume_open(path, &fv) == 0, "open the planted field");
        if (fv != NULL) {
            int hn = 0, wn = 0;
            double umin = 0.0, vmin = 0.0;
            float *orc = NULL, *Xn = NULL;
            uint8_t *vn = NULL;
            size_t nv = 0, checked = 0, bad = 0;
            FFT_CHECK(FiberFrame_field(ar, X, valid, H, W, 0, 0, W, H, fv, &ff) == 0, "field");
            FFT_CHECK(ff.samples_with_fibres > ff.samples / 2, "most samples see the planted fibres");
            FFT_CHECK(FiberFrame_solve(ar, &ff, &st) == 0, "solve (TAUCS)");
            FFT_CHECK(FiberFrame_audit(ar, &ff, &st, NULL) == 0, "audit");
            fprintf(stderr, "[selftest] fiber_frame planted bow: climb %.1f -> %.1f px, H |p50| %.2f -> %.2f deg, area %.3f..%.3f, "
                            "flipped %zu\n", st.climb_old_range, st.climb_new_range, st.H_old_abs.p50, st.H_new_abs.p50,
                    st.area_p1, st.area_p99, st.flipped);
            FFT_CHECK(st.climb_old_range > 12.0 && st.climb_old_range < 30.0, "the planted bow reads ~20 px");
            FFT_CHECK(st.climb_new_range < 3.0, "the re-flattened page is straight");
            FFT_CHECK(st.H_new_abs.p50 < 1.0, "H lies along the new rows");
            FFT_CHECK(st.flipped == 0, "no folded nodes");
            FFT_CHECK(st.area_p1 > 0.85 && st.area_p99 < 1.15, "area within 15%");
            FFT_CHECK(FiberFrame_inverse(ar, &ff, &st, 2, &hn, &wn, &umin, &vmin, &orc) == 0, "inverse");
            if (orc != NULL) {
                for (int r = 0; r < hn; r += 7)
                    for (int c = 0; c < wn; c += 7) {
                        size_t o = ((size_t)r * (size_t)wn + (size_t)c) * 2;
                        double fu = 0.0, fv2 = 0.0;
                        if (isnan(orc[o])) continue;
                        if (ff_grid_bilinear(ff.U, ff.gh, ff.gw, ((double)orc[o] - ff.gy0) / ff.cell,
                                             ((double)orc[o + 1] - ff.gx0) / ff.cell, &fu) != 0 ||
                            ff_grid_bilinear(ff.V, ff.gh, ff.gw, ((double)orc[o] - ff.gy0) / ff.cell,
                                             ((double)orc[o + 1] - ff.gx0) / ff.cell, &fv2) != 0) continue;
                        checked++;
                        if (fabs(fu - (umin + c)) > 0.06 || fabs(fv2 - (vmin + r)) > 0.06) bad++;
                    }
                FFT_CHECK(checked > 100 && bad == 0, "the inverse maps back onto the forward map");
                FFT_CHECK(FiberFrame_resample(ar, X, valid, H, W, 0, 0, W, H, orc, hn, wn, &Xn, &vn, &nv) == 0, "resample");
                if (Xn != NULL) {
                    size_t nfar = 0, seen = 0;
                    for (size_t o = 0; o < (size_t)hn * (size_t)wn; o += 13) {
                        double pr = (double)orc[o * 2], pc = (double)orc[o * 2 + 1];
                        if (!vn[o]) continue;
                        seen++;
                        for (int k = 0; k < 3; k++)
                            if (fabs((double)Xn[o * 3 + (size_t)k] - (O[k] + pc * e1[k] + pr * e2[k])) > 1e-3) { nfar++; break; }
                    }
                    FFT_CHECK(seen > 100 && nfar == 0, "the new tifxyz samples the plane where the inverse says");
                    FFT_CHECK((double)nv > 0.8 * H * W, "the new page keeps most of the surface");
                }
            }
            FiberVolume_close(&fv);
        }
        for (int a = 0; a < 3; a++) free(chunk[a]);
        free(X);
        free(valid);
        Arena_dispose(&ar);
    }
    if (fails == 0) ves_selftest_clean(root);
    fprintf(stderr, "[selftest] fiber_frame %s (%d failures)\n", fails ? "FAIL" : "PASS", fails);
    return fails;
}
