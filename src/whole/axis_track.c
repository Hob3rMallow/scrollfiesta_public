#include "axis_track.h"
#include "axis_warp.h"
#include "../common/pipeline_constants.h"
#include "../common/png_write.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---- reservoirs -----------------------------------------------------------*/

static uint32_t at_rng(uint32_t *s)
{
    uint32_t x = *s;
    x ^= x << 13; x ^= x >> 17; x ^= x << 5;
    *s = x;
    return x;
}

int AxisSlab_init(AxisSlab *s, double z0, double z1,
                  size_t face_cap, size_t vert_cap, uint32_t seed)
{
    if (s == NULL || face_cap == 0 || vert_cap == 0) return -1;
    memset(s, 0, sizeof *s);
    s->z0 = z0; s->z1 = z1;
    s->fy  = (double *)malloc(face_cap * sizeof(double));
    s->fx  = (double *)malloc(face_cap * sizeof(double));
    s->fny = (double *)malloc(face_cap * sizeof(double));
    s->fnx = (double *)malloc(face_cap * sizeof(double));
    s->fa  = (double *)malloc(face_cap * sizeof(double));
    s->vy  = (double *)malloc(vert_cap * sizeof(double));
    s->vx  = (double *)malloc(vert_cap * sizeof(double));
    if (s->fy == NULL || s->fx == NULL || s->fny == NULL || s->fnx == NULL ||
        s->fa == NULL || s->vy == NULL || s->vx == NULL) {
        AxisSlab_dispose(s);
        return -1;
    }
    s->nf_cap = face_cap;
    s->nv_cap = vert_cap;
    s->rng = seed ? seed : 0x9E3779B9u;
    return 0;
}

void AxisSlab_dispose(AxisSlab *s)
{
    if (s == NULL) return;
    free(s->fy); free(s->fx); free(s->fny); free(s->fnx); free(s->fa);
    free(s->vy); free(s->vx);
    memset(s, 0, sizeof *s);
}

void AxisSlab_add_face(AxisSlab *s, const float a[3], const float b[3],
                       const float c[3])
{
    /* face normal (z,y,x) = (b-a) x (c-a); in-plane part = (ny, nx) */
    double e1[3] = { (double)b[0] - a[0], (double)b[1] - a[1], (double)b[2] - a[2] };
    double e2[3] = { (double)c[0] - a[0], (double)c[1] - a[1], (double)c[2] - a[2] };
    double nz = e1[1] * e2[2] - e1[2] * e2[1];
    double ny = e1[2] * e2[0] - e1[0] * e2[2];
    double nx = e1[0] * e2[1] - e1[1] * e2[0];
    double len = sqrt(nz * nz + ny * ny + nx * nx);
    double lp = sqrt(ny * ny + nx * nx);
    size_t slot = 0;
    if (s == NULL || !(len > 1e-12)) return;
    if (lp < 0.3 * len) return;       /* near-horizontal face: no centre information */
    s->nf_seen++;
    if (s->nf < s->nf_cap) slot = s->nf++;
    else {
        slot = (size_t)(at_rng(&s->rng) % (uint32_t)s->nf_seen);
        if (slot >= s->nf_cap) return;
    }
    s->fy[slot]  = ((double)a[1] + b[1] + c[1]) / 3.0;
    s->fx[slot]  = ((double)a[2] + b[2] + c[2]) / 3.0;
    s->fny[slot] = ny / lp;
    s->fnx[slot] = nx / lp;
    s->fa[slot]  = 0.5 * len;
}

void AxisSlab_add_vertex(AxisSlab *s, const float p[3])
{
    size_t slot = 0;
    if (s == NULL) return;
    s->nv_seen++;
    if (s->nv < s->nv_cap) slot = s->nv++;
    else {
        slot = (size_t)(at_rng(&s->rng) % (uint32_t)s->nv_seen);
        if (slot >= s->nv_cap) return;
    }
    s->vy[slot] = (double)p[1];
    s->vx[slot] = (double)p[2];
}

/* ---- primary: normal-line least squares -----------------------------------*/

static int at_cmp_double(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return x < y ? -1 : (x > y);
}

static double at_quantile(double *v, size_t n, double q)
{
    size_t k = 0;
    if (n == 0) return 0.0;
    qsort(v, n, sizeof *v, at_cmp_double);
    k = (size_t)(q * (double)(n - 1) + 0.5);
    if (k >= n) k = n - 1;
    return v[k];
}

int AxisTrack_primary(const AxisSlab *s, double prior_y, double prior_x,
                      AxisRow *row)
{
    size_t *idx = NULL;
    double *res = NULL, *tmp = NULL;
    unsigned char *keep = NULL;
    size_t n_in = 0, i = 0;
    double cy = prior_y, cx = prior_x;
    double lmin = 0.0, lmax = 0.0;
    int ok = 0;
    if (s == NULL || row == NULL) return -1;
    row->primary_ok = 0; row->n_in = 0; row->cond = 0.0; row->resid_p50 = 0.0;
    row->sectors = 0; row->void_ratio = 1.0;
    row->y1 = prior_y; row->x1 = prior_x; row->jump = 0.0;
    idx  = (size_t *)malloc((s->nf + 1) * sizeof *idx);
    res  = (double *)malloc((s->nf + 1) * sizeof *res);
    tmp  = (double *)malloc((s->nf + 1) * sizeof *tmp);
    keep = (unsigned char *)malloc(s->nf + 1);
    if (idx == NULL || res == NULL || tmp == NULL || keep == NULL) goto done;
    for (i = 0; i < s->nf; i++) {
        double dy = s->fy[i] - prior_y, dx = s->fx[i] - prior_x;
        if (dy * dy + dx * dx <= AXIS_TRACE_LOCAL_R * AXIS_TRACE_LOCAL_R)
            idx[n_in++] = i;
    }
    row->n_in = n_in;
    if (n_in < AXIS_TRACE_MIN_SAMPLES) goto done;
    {   /* conditioning of the whole in-window normal scatter */
        double a11 = 0.0, a12 = 0.0, a22 = 0.0, tr = 0.0, disc = 0.0, l1 = 0.0, l2 = 0.0;
        for (i = 0; i < n_in; i++) {
            double my = s->fny[idx[i]], mx = s->fnx[idx[i]];
            a11 += 1.0 - my * my;
            a12 += -my * mx;
            a22 += 1.0 - mx * mx;
        }
        tr = a11 + a22;
        disc = sqrt(fmax(0.0, 0.25 * (a11 - a22) * (a11 - a22) + a12 * a12));
        l1 = 0.5 * tr + disc; l2 = 0.5 * tr - disc;
        row->cond = l1 > 0.0 ? l2 / l1 : 0.0;
    }
    memset(keep, 1, n_in);
    for (int round = 0; round <= AXIS_TRACE_IRLS_ROUNDS; round++) {
        /* normal equations of sum_i w_i (c-p_i)^T (I - m m^T) (c-p_i) */
        double a11 = 0.0, a12 = 0.0, a22 = 0.0, b1 = 0.0, b2 = 0.0;
        double tr = 0.0, det = 0.0, disc = 0.0, thr = 0.0;
        size_t nk = 0;
        for (i = 0; i < n_in; i++) {
            size_t j = idx[i];
            double my = s->fny[j], mx = s->fnx[j];
            double p11 = 1.0 - my * my, p12 = -my * mx, p22 = 1.0 - mx * mx;
            if (!keep[i]) continue;
            a11 += p11; a12 += p12; a22 += p22;
            b1 += p11 * s->fy[j] + p12 * s->fx[j];
            b2 += p12 * s->fy[j] + p22 * s->fx[j];
            nk++;
        }
        det = a11 * a22 - a12 * a12;
        if (!(det > 1e-9) || nk < 3) goto done;
        cy = ( a22 * b1 - a12 * b2) / det;
        cx = (-a12 * b1 + a11 * b2) / det;
        tr = a11 + a22;
        disc = sqrt(fmax(0.0, 0.25 * (a11 - a22) * (a11 - a22) + a12 * a12));
        lmax = 0.5 * tr + disc; lmin = 0.5 * tr - disc;
        /* residual = distance from c to every in-window normal line */
        for (i = 0; i < n_in; i++) {
            size_t j = idx[i];
            double dy = cy - s->fy[j], dx = cx - s->fx[j];
            double along = dy * s->fny[j] + dx * s->fnx[j];
            double ry = dy - along * s->fny[j], rx = dx - along * s->fnx[j];
            res[i] = sqrt(ry * ry + rx * rx);
        }
        if (round == AXIS_TRACE_IRLS_ROUNDS) break;
        memcpy(tmp, res, n_in * sizeof(double));
        thr = at_quantile(tmp, n_in, AXIS_TRACE_TRIM_FRAC);
        for (i = 0; i < n_in; i++) keep[i] = res[i] <= thr;
    }
    memcpy(tmp, res, n_in * sizeof(double));
    row->resid_p50 = at_quantile(tmp, n_in, 0.5);
    (void)lmin; (void)lmax;
    row->y1 = cy; row->x1 = cx;
    row->jump = hypot(cy - prior_y, cx - prior_x);
    {   /* angular coverage (the wraps must surround the centre) and the
         * void test (wrap-area density inside the disc vs the annulus) */
        int hit[12] = {0};
        int ns = 0;
        double a_in = 0.0, a_ann = 0.0;
        const double r0 = AXIS_TRACE_VOID_R, r1 = 3.0 * AXIS_TRACE_VOID_R;
        for (i = 0; i < n_in; i++) {
            size_t j = idx[i];
            double dy = s->fy[j] - cy, dx = s->fx[j] - cx, rr = sqrt(dy * dy + dx * dx);
            int sct = 0;
            if (rr < r0) a_in += s->fa[j];
            else if (rr < r1) a_ann += s->fa[j];
            if (!keep[i]) continue;
            sct = (int)((atan2(dx, dy) + M_PI) * (6.0 / M_PI));
            if (sct < 0) sct = 0;
            if (sct > 11) sct = 11;
            hit[sct] = 1;
        }
        for (int k = 0; k < 12; k++) ns += hit[k];
        row->sectors = ns;
        {
            double d_in = a_in / (M_PI * r0 * r0);
            double d_ann = a_ann / (M_PI * (r1 * r1 - r0 * r0));
            row->void_ratio = d_ann > 1e-9 ? d_in / d_ann : 1.0;
        }
    }
    ok = row->cond >= AXIS_TRACE_MIN_COND &&
         row->sectors >= AXIS_TRACE_MIN_SECTORS &&
         row->resid_p50 <= AXIS_TRACE_MAX_RESID_P50 &&
         row->jump <= AXIS_TRACE_MAX_JUMP;
    row->primary_ok = ok;
done:
    free(idx); free(res); free(tmp); free(keep);
    return ok ? 0 : -1;
}

/* ---- secondary: polar first-harmonic recentring ---------------------------*/

#define AT_NTH   ((int)(360.0 / AXIS_TRACE_POLAR_DTH))
#define AT_NR    ((int)(AXIS_TRACE_POLAR_RMAX / AXIS_TRACE_POLAR_DR))
#define AT_PEAKS_PER_BIN 64
#define AT_MIN_SUPPORT 24.0      /* vox^2 of wrap area in the smoothed (1,2,1) cell sum */
#define AT_VERTEX_AREA 4.0       /* vox^2 credited to a bare vertex sample */
#define AT_LINK_TOL 5.0          /* vox: peak-to-peak link tolerance across theta bins */
#define AT_LINK_MISSES 2         /* consecutive empty bins a chain may cross (cube seams) */
#define AT_MERGE_TOL 3.0         /* vox: thick-shell front/back peaks merge */
#define AT_MAX_RIDGES 512
#define AT_MAX_ITERS 24
#define AT_CONVERGED 1.0         /* vox: mean first-harmonic shift below this = converged */
#define AT_SETTLED 3.0           /* vox: at the iteration cap a shift below this still counts */
int AxisTrack_verbose = 0;

typedef struct { float r; unsigned char claimed; } AtPeak;

typedef struct {
    double R, k, a, b;   /* r(theta) = R + k theta + a cos theta + b sin theta */
    double c2, s2;       /* + c2 cos 2theta + s2 sin 2theta (ellipticity) */
    double rmean;
    int span;            /* theta bins */
} AtRidge;

/* Solve an n x n (n <= 6) system by Gaussian elimination with pivoting. */
#define AT_NPARAM 6
static int at_solve(int n, double A[AT_NPARAM][AT_NPARAM], double b[AT_NPARAM],
                    double x[AT_NPARAM])
{
    int p[AT_NPARAM];
    for (int i = 0; i < n; i++) p[i] = i;
    for (int c = 0; c < n; c++) {
        int best = c;
        for (int r = c + 1; r < n; r++)
            if (fabs(A[p[r]][c]) > fabs(A[p[best]][c])) best = r;
        { int t = p[c]; p[c] = p[best]; p[best] = t; }
        if (fabs(A[p[c]][c]) < 1e-12) return -1;
        for (int r = c + 1; r < n; r++) {
            double f = A[p[r]][c] / A[p[c]][c];
            for (int k = c; k < n; k++) A[p[r]][k] -= f * A[p[c]][k];
            b[p[r]] -= f * b[p[c]];
        }
    }
    for (int c = n - 1; c >= 0; c--) {
        double acc = b[p[c]];
        for (int k = c + 1; k < n; k++) acc -= A[p[c]][k] * x[k];
        x[c] = acc / A[p[c]][c];
    }
    return 0;
}

/* r(theta) = R + k theta + a cos + b sin + c cos 2theta + d sin 2theta: the
 * second harmonic absorbs the wraps' ellipticity, which over a partial arc
 * would otherwise leak into the first harmonic (the centre) and the slope
 * (the pitch/sense). */
static int at_fit_ridge(const double *th, const double *r, int n, AtRidge *out)
{
    double A[AT_NPARAM][AT_NPARAM] = {{0}}, b[AT_NPARAM] = {0}, x[AT_NPARAM] = {0};
    for (int i = 0; i < n; i++) {
        double f[AT_NPARAM] = { 1.0, th[i], cos(th[i]), sin(th[i]),
                                cos(2.0 * th[i]), sin(2.0 * th[i]) };
        for (int p = 0; p < AT_NPARAM; p++) {
            for (int q = 0; q < AT_NPARAM; q++) A[p][q] += f[p] * f[q];
            b[p] += f[p] * r[i];
        }
    }
    if (at_solve(AT_NPARAM, A, b, x) != 0) return -1;
    out->R = x[0]; out->k = x[1]; out->a = x[2]; out->b = x[3];
    out->c2 = x[4]; out->s2 = x[5];
    out->span = n;
    out->rmean = 0.0;
    for (int i = 0; i < n; i++) out->rmean += r[i];
    out->rmean /= (double)n;
    return 0;
}

/* One recentring pass about (cy, cx): rasterize, extract ridges, fit.
 * Returns the number of ridges spanning the minimum arc, fills shift. */
static int at_polar_pass(const AxisSlab *s, double cy, double cx, int stride,
                         float *acc, AtPeak *peaks, int *npeaks,
                         AtRidge *ridges, int max_ridges,
                         double *shift_y, double *shift_x)
{
    const int nth = AT_NTH, nr = AT_NR;
    const double dth = AXIS_TRACE_POLAR_DTH * M_PI / 180.0;
    const int min_span = (int)(AXIS_TRACE_RIDGE_MIN_ARC / AXIS_TRACE_POLAR_DTH);
    int nridge = 0;
    double *th_buf = NULL, *r_buf = NULL;
    AtPeak **pk_buf = NULL;
    double sum_a = 0.0, sum_b = 0.0, sum_w = 0.0;
    int total_peaks = 0, longest = 0;
    memset(acc, 0, (size_t)nth * (size_t)nr * sizeof *acc);
    /* area-weighted: a wrap crossing a cell credits its surface area, so
     * the raster is independent of how densely the mesh is triangulated */
    for (size_t i = 0; i < s->nf; i += (size_t)stride) {
        double dy = s->fy[i] - cy, dx = s->fx[i] - cx;
        double r = sqrt(dy * dy + dx * dx);
        int ir = (int)(r / AXIS_TRACE_POLAR_DR);
        int it = (int)((atan2(dx, dy) + M_PI) / dth);
        if (ir < 0 || ir >= nr) continue;
        if (it < 0) it = 0;
        if (it >= nth) it = nth - 1;
        acc[(size_t)it * (size_t)nr + (size_t)ir] += (float)(s->fa[i] * (double)stride);
    }
    if (s->nf == 0) {
        for (size_t i = 0; i < s->nv; i += (size_t)stride) {
            double dy = s->vy[i] - cy, dx = s->vx[i] - cx;
            double r = sqrt(dy * dy + dx * dx);
            int ir = (int)(r / AXIS_TRACE_POLAR_DR);
            int it = (int)((atan2(dx, dy) + M_PI) / dth);
            if (ir < 0 || ir >= nr) continue;
            if (it < 0) it = 0;
            if (it >= nth) it = nth - 1;
            acc[(size_t)it * (size_t)nr + (size_t)ir] += (float)(AT_VERTEX_AREA * stride);
        }
    }
    /* peaks along r per theta bin (3-tap smoothed, parabolic refinement) */
    for (int it = 0; it < nth; it++) {
        const float *col = acc + (size_t)it * (size_t)nr;
        int np = 0;
        double last_r = -1e9;
        for (int ir = 1; ir < nr - 1; ir++) {
            double sm = (double)col[ir - 1] + 2.0 * col[ir] + col[ir + 1];
            double sl = (ir >= 2 ? (double)col[ir - 2] : 0.0) + 2.0 * col[ir - 1] + col[ir];
            double sr = (double)col[ir] + 2.0 * col[ir + 1] + (ir + 2 < nr ? (double)col[ir + 2] : 0.0);
            double denom = 0.0, rr = 0.0;
            if (sm < AT_MIN_SUPPORT || sm <= sl || sm < sr) continue;
            denom = sl - 2.0 * sm + sr;
            rr = (double)ir + (fabs(denom) > 1e-9 ? 0.5 * (sl - sr) / denom : 0.0);
            rr *= AXIS_TRACE_POLAR_DR;
            if (rr - last_r < AT_MERGE_TOL) { last_r = rr; continue; }  /* thick shell */
            if (np < AT_PEAKS_PER_BIN) {
                peaks[(size_t)it * AT_PEAKS_PER_BIN + (size_t)np].r = (float)rr;
                peaks[(size_t)it * AT_PEAKS_PER_BIN + (size_t)np].claimed = 0;
                np++;
            }
            last_r = rr;
        }
        npeaks[it] = np;
        total_peaks += np;
    }
    th_buf = (double *)malloc((size_t)(2 * nth) * sizeof(double));
    r_buf  = (double *)malloc((size_t)(2 * nth) * sizeof(double));
    pk_buf = (AtPeak **)malloc((size_t)(2 * nth) * sizeof(AtPeak *));
    if (th_buf == NULL || r_buf == NULL || pk_buf == NULL) {
        free(th_buf); free(r_buf); free(pk_buf);
        return 0;
    }
    /* chain peaks across theta (two laps, starts in the first lap only) */
    for (int it0 = 0; it0 < nth; it0++) {
        for (int p0 = 0; p0 < npeaks[it0]; p0++) {
            AtPeak *start = &peaks[(size_t)it0 * AT_PEAKS_PER_BIN + (size_t)p0];
            int len = 0, misses = 0, it = it0;
            double r_prev = 0.0;
            if (start->claimed) continue;
            r_prev = (double)start->r;
            start->claimed = 1;
            th_buf[len] = (double)it0 * dth - M_PI;
            r_buf[len] = r_prev;
            pk_buf[len] = start;
            len++;
            for (it = it0 + 1; it < it0 + nth; it++) {
                int itm = it % nth, best = -1;
                double best_d = AT_LINK_TOL;
                for (int p = 0; p < npeaks[itm]; p++) {
                    AtPeak *pk = &peaks[(size_t)itm * AT_PEAKS_PER_BIN + (size_t)p];
                    double d = fabs((double)pk->r - r_prev);
                    if (!pk->claimed && d < best_d) { best_d = d; best = p; }
                }
                if (best < 0) {
                    if (++misses > AT_LINK_MISSES) break;
                    continue;
                }
                misses = 0;
                {
                    AtPeak *pk = &peaks[(size_t)itm * AT_PEAKS_PER_BIN + (size_t)best];
                    pk->claimed = 1;
                    r_prev = (double)pk->r;
                    th_buf[len] = (double)it * dth - M_PI;
                    r_buf[len] = r_prev;
                    pk_buf[len] = pk;
                    len++;
                }
            }
            if (len > longest) longest = len;
            {
                AtRidge rg;
                int kept = 0;
                if (len >= min_span && nridge < max_ridges &&
                    at_fit_ridge(th_buf, r_buf, len, &rg) == 0) {
                    ridges[nridge++] = rg;
                    sum_a += rg.a * (double)len;
                    sum_b += rg.b * (double)len;
                    sum_w += (double)len;
                    kept = 1;
                }
                /* a rejected (short) chain must not steal peaks from the
                 * longer ridge that starts after the coverage gap */
                if (!kept)
                    for (int j = 0; j < len; j++) pk_buf[j]->claimed = 0;
            }
        }
    }
    free(th_buf); free(r_buf); free(pk_buf);
    if (AxisTrack_verbose > 1)
        fprintf(stderr, "      pass about (%.1f, %.1f): %d peaks, longest chain %d bins, %d ridges >= %d bins\n",
                cy, cx, total_peaks, longest, nridge, min_span);
    *shift_y = sum_w > 0.0 ? sum_a / sum_w : 0.0;
    *shift_x = sum_w > 0.0 ? sum_b / sum_w : 0.0;
    return nridge;
}

static int at_secondary_impl(const AxisSlab *s, double sy, double sx, int stride,
                             int max_iters, AxisRow *row)
{
    float *acc = NULL;
    AtPeak *peaks = NULL;
    int *npeaks = NULL;
    AtRidge *ridges = NULL;
    double cy = sy, cx = sx, shy = 0.0, shx = 0.0;
    int nridge = 0, iter = 0, converged = 0;
    row->secondary_ok = 0; row->ridges = 0; row->d1 = 0.0; row->pitch = 0.0;
    row->r_wall = 0.0; row->y2 = sy; row->x2 = sx;
    if (s == NULL || s->nv + s->nf < 100) return -1;
    acc    = (float *)malloc((size_t)AT_NTH * (size_t)AT_NR * sizeof *acc);
    peaks  = (AtPeak *)malloc((size_t)AT_NTH * AT_PEAKS_PER_BIN * sizeof *peaks);
    npeaks = (int *)calloc((size_t)AT_NTH, sizeof *npeaks);
    ridges = (AtRidge *)malloc(AT_MAX_RIDGES * sizeof *ridges);
    if (acc == NULL || peaks == NULL || npeaks == NULL || ridges == NULL) {
        free(acc); free(peaks); free(npeaks); free(ridges);
        return -1;
    }
    for (iter = 0; iter < max_iters; iter++) {
        nridge = at_polar_pass(s, cy, cx, stride, acc, peaks, npeaks,
                               ridges, AT_MAX_RIDGES, &shy, &shx);
        if (AxisTrack_verbose > 1)
            fprintf(stderr, "      iter %d: centre (%.1f, %.1f) ridges=%d shift (%.2f, %.2f)\n",
                    iter, cy, cx, nridge, shy, shx);
        if (nridge == 0) break;
        cy += shy; cx += shx;
        if (hypot(shy, shx) < AT_CONVERGED) { converged = 1; break; }
        if (iter + 1 == max_iters && hypot(shy, shx) < AT_SETTLED) converged = 1;
    }
    if (AxisTrack_verbose > 0)
        fprintf(stderr, "    secondary from (%.1f, %.1f): %d iter -> (%.1f, %.1f) ridges=%d last shift %.2f %s\n",
                sy, sx, iter, cy, cx, nridge, hypot(shy, shx), converged ? "converged" : "NOT converged");
    if (nridge > 0 && converged) {
        double *d = NULL;
        double ksum = 0.0, kw = 0.0, rmin = 1e30;
        /* One final pass at the converged centre supplies the row's numbers.
         * Its ridge count can be larger than the pass that declared
         * convergence, so size diagnostics from the final count, not the
         * stale pre-update count. */
        nridge = at_polar_pass(s, cy, cx, stride, acc, peaks, npeaks,
                               ridges, AT_MAX_RIDGES, &shy, &shx);
        if (nridge > 0)
            d = (double *)malloc((size_t)nridge * sizeof(double));
        for (int i = 0; i < nridge; i++) {
            if (d != NULL) d[i] = hypot(ridges[i].a, ridges[i].b);
            ksum += ridges[i].k * (double)ridges[i].span;
            kw += (double)ridges[i].span;
            if (ridges[i].rmean < rmin) rmin = ridges[i].rmean;
        }
        row->d1 = d != NULL && nridge > 0 ? at_quantile(d, (size_t)nridge, 0.5) : 0.0;
        row->pitch = kw > 0.0 ? 2.0 * M_PI * ksum / kw : 0.0;
        row->r_wall = nridge > 0 ? rmin : 0.0;
        free(d);
    }
    row->y2 = cy; row->x2 = cx;
    row->ridges = (size_t)(nridge > 0 ? nridge : 0);
    row->secondary_ok = converged && nridge >= 3;
    free(acc); free(peaks); free(npeaks); free(ridges);
    return row->secondary_ok ? 0 : -1;
}

int AxisTrack_polar_debug(const AxisSlab *s, double cy, double cx,
                          int r_max, const char *png_path)
{
    const int nth = AT_NTH, nr = AT_NR;
    float *acc = NULL;
    AtPeak *peaks = NULL;
    int *npeaks = NULL;
    AtRidge *ridges = NULL;
    uint8_t *rgb = NULL;
    double shy = 0.0, shx = 0.0;
    int nridge = 0, W = nth, H = r_max < nr ? r_max : nr, rc = -1;
    if (s == NULL || png_path == NULL || H <= 0) return -1;
    acc    = (float *)malloc((size_t)nth * (size_t)nr * sizeof *acc);
    peaks  = (AtPeak *)malloc((size_t)nth * AT_PEAKS_PER_BIN * sizeof *peaks);
    npeaks = (int *)calloc((size_t)nth, sizeof *npeaks);
    ridges = (AtRidge *)malloc(AT_MAX_RIDGES * sizeof *ridges);
    rgb    = (uint8_t *)calloc((size_t)W * (size_t)H * 3, 1);
    if (acc == NULL || peaks == NULL || npeaks == NULL || ridges == NULL || rgb == NULL)
        goto done;
    nridge = at_polar_pass(s, cy, cx, 1, acc, peaks, npeaks, ridges,
                           AT_MAX_RIDGES, &shy, &shx);
    for (int ir = 0; ir < H; ir++)
        for (int it = 0; it < W; it++) {
            double v = acc[(size_t)it * (size_t)nr + (size_t)ir];
            uint8_t g = (uint8_t)(v > 0.0 ? fmin(255.0, 40.0 + 45.0 * log(1.0 + v)) : 0);
            uint8_t *px = rgb + ((size_t)ir * (size_t)W + (size_t)it) * 3;
            px[0] = g; px[1] = g; px[2] = g;
        }
    for (int it = 0; it < W; it++)
        for (int p = 0; p < npeaks[it]; p++) {
            const AtPeak *pk = &peaks[(size_t)it * AT_PEAKS_PER_BIN + (size_t)p];
            int ir = (int)((double)pk->r / AXIS_TRACE_POLAR_DR + 0.5);
            uint8_t *px = NULL;
            if (ir < 0 || ir >= H) continue;
            px = rgb + ((size_t)ir * (size_t)W + (size_t)it) * 3;
            if (pk->claimed) { px[0] = 40; px[1] = 255; px[2] = 60; }
            else { px[0] = 255; px[1] = 50; px[2] = 50; }
        }
    rc = PngWrite_rgb(png_path, W, H, rgb);
    fprintf(stderr, "  polar debug about (%.1f, %.1f): %d ridges, mean shift (%.2f, %.2f) -> %s\n",
            cy, cx, nridge, shy, shx, png_path);
done:
    free(acc); free(peaks); free(npeaks); free(ridges); free(rgb);
    return rc;
}

int AxisTrack_secondary(const AxisSlab *s, double start_y, double start_x,
                        AxisRow *row)
{
    if (row == NULL) return -1;
    return at_secondary_impl(s, start_y, start_x, 1, AT_MAX_ITERS, row);
}

/* ---- seed ------------------------------------------------------------------*/

int AxisTrack_seed(const AxisSlab *s, double *seed_y, double *seed_x,
                   double *cl_y, double *cl_x, size_t *cl_members,
                   double *cl_d1, size_t *cl_ridges, size_t max_clusters,
                   size_t *n_clusters)
{
    return AxisTrack_seed_in(s, 0.0, 0.0, 0.0, seed_y, seed_x, cl_y, cl_x,
                             cl_members, cl_d1, cl_ridges, max_clusters,
                             n_clusters);
}

int AxisTrack_seed_in(const AxisSlab *s, double win_y, double win_x, double win_r,
                      double *seed_y, double *seed_x,
                      double *cl_y, double *cl_x, size_t *cl_members,
                      double *cl_d1, size_t *cl_ridges, size_t max_clusters,
                      size_t *n_clusters)
{
    double ylo = 1e30, yhi = -1e30, xlo = 1e30, xhi = -1e30;
    size_t nc = 0, best = 0;
    int found = 0;
    if (s == NULL || seed_y == NULL || seed_x == NULL || n_clusters == NULL)
        return -1;
    *n_clusters = 0;
    for (size_t i = 0; i < s->nv; i++) {
        if (s->vy[i] < ylo) ylo = s->vy[i];
        if (s->vy[i] > yhi) yhi = s->vy[i];
        if (s->vx[i] < xlo) xlo = s->vx[i];
        if (s->vx[i] > xhi) xhi = s->vx[i];
    }
    if (win_r > 0.0) {
        ylo = fmax(ylo, win_y - win_r); yhi = fmin(yhi, win_y + win_r);
        xlo = fmax(xlo, win_x - win_r); xhi = fmin(xhi, win_x + win_r);
    }
    if (!(yhi > ylo) || !(xhi > xlo)) return -1;
    for (double y = ylo + 0.5 * AXIS_TRACE_SEED_GRID; y < yhi; y += AXIS_TRACE_SEED_GRID) {
        for (double x = xlo + 0.5 * AXIS_TRACE_SEED_GRID; x < xhi; x += AXIS_TRACE_SEED_GRID) {
            AxisRow r;
            double py = y, px = x;
            size_t c = 0;
            int ok = 0;
            memset(&r, 0, sizeof r);
            /* the local window follows the estimate: iterate to a fixed point */
            for (int it = 0; it < 6; it++) {
                ok = AxisTrack_primary(s, py, px, &r) == 0;
                if (r.n_in < AXIS_TRACE_MIN_SAMPLES) { ok = 0; break; }
                if (hypot(r.y1 - py, r.x1 - px) < 1.0) break;
                py = r.y1; px = r.x1;
            }
            if (AxisTrack_verbose > 0)
                fprintf(stderr, "    seed start (%.0f, %.0f) -> (%.1f, %.1f) n=%zu cond=%.2f sectors=%d void=%.2f resid=%.1f %s\n",
                        y, x, r.y1, r.x1, r.n_in, r.cond, r.sectors, r.void_ratio, r.resid_p50,
                        ok ? "lock" : "abstain");
            if (!ok) continue;
            for (c = 0; c < nc; c++)
                if (hypot(cl_y[c] - r.y1, cl_x[c] - r.x1) <= 30.0) break;
            if (c == nc) {
                if (nc >= max_clusters) continue;
                cl_y[nc] = r.y1; cl_x[nc] = r.x1; cl_members[nc] = 0;
                cl_d1[nc] = r.resid_p50; cl_ridges[nc] = r.n_in;
                nc++;
            }
            cl_members[c]++;
            /* the cluster keeps its lowest-residual member; cl_ridges carries
             * that member's void ratio x 1000 (0 = emptiest disc) */
            if (r.resid_p50 < cl_d1[c] || cl_members[c] == 1) {
                cl_ridges[c] = (size_t)(fmin(r.void_ratio, 9.999) * 1000.0);
                cl_d1[c] = r.resid_p50;
                cl_y[c] = r.y1; cl_x[c] = r.x1;
            }
        }
    }
    *n_clusters = nc;
    /* the umbilicus sits in a VOID: candidates whose inner disc is emptier
     * than AXIS_TRACE_VOID_MAX of the annulus qualify; the lowest residual
     * among them wins.  Without any qualifying cluster nothing is seeded. */
    for (size_t c = 0; c < nc; c++) {
        if ((double)cl_ridges[c] / 1000.0 >= AXIS_TRACE_VOID_MAX) continue;
        if (!found || cl_d1[c] < cl_d1[best]) { best = c; found = 1; }
    }
    if (!found) return -1;
    *seed_y = cl_y[best]; *seed_x = cl_x[best];
    return 0;
}

/* ---- finish: fuse, median, slope limit, gaps, extrapolation ----------------*/

size_t AxisTrack_finish(AxisRow *rows, size_t n)
{
    size_t nlock = 0, last_lock = 0;
    int have_lock = 0;
    if (rows == NULL || n == 0) return 0;
    /* 1. fuse: the primary is the measured position authority.  The polar
     * harmonic estimator is known to lock sparsely on PHerc0139's folded,
     * non-circular wraps, so its disagreement is reported but cannot veto an
     * otherwise valid primary row. */
    for (size_t i = 0; i < n; i++) {
        AxisRow *r = &rows[i];
        r->agree = r->secondary_ok ? hypot(r->y1 - r->y2, r->x1 - r->x2) : -1.0;
        if (r->primary_ok) {
            r->y = r->y1; r->x = r->x1; r->status = AXIS_ROW_LOCK;
        } else {
            r->y = 0.0; r->x = 0.0; r->status = AXIS_ROW_ABSTAIN;
        }
    }
    /* 2. slope limit on the raw fused rows: abstain the later row of a
     * too-steep segment (before the median, which would hide the jump) */
    for (size_t i = 0; i < n; i++) {
        if (rows[i].status != AXIS_ROW_LOCK) continue;
        if (have_lock) {
            double dz = rows[i].z - rows[last_lock].z;
            double slope = dz > 0.0 ? hypot(rows[i].y - rows[last_lock].y,
                                            rows[i].x - rows[last_lock].x) / dz : 0.0;
            if (slope > AXIS_TABLE_MAX_SLOPE) { rows[i].status = AXIS_ROW_ABSTAIN; continue; }
        }
        last_lock = i; have_lock = 1;
    }
    /* 3. (no smoothing) The 5-sample median measured 2026-09-02 on the
     * 21x3x3 column turned a monotone drift with sparse locks into a
     * staircase (equal values on neighbouring rows); the raw locked
     * estimates stay, guarded by the slope limit and the forward/backward
     * agreement. */
    /* 4. interpolate short gaps between locks */
    have_lock = 0;
    for (size_t i = 0; i < n; i++) {
        if (rows[i].status != AXIS_ROW_LOCK) continue;
        if (have_lock && i - last_lock > 1 && i - last_lock - 1 <= AXIS_TABLE_MAX_GAP_ROWS) {
            for (size_t j = last_lock + 1; j < i; j++) {
                double t = (rows[j].z - rows[last_lock].z) / (rows[i].z - rows[last_lock].z);
                rows[j].y = rows[last_lock].y + t * (rows[i].y - rows[last_lock].y);
                rows[j].x = rows[last_lock].x + t * (rows[i].x - rows[last_lock].x);
                rows[j].status = AXIS_ROW_INTERP;
            }
        }
        last_lock = i; have_lock = 1;
    }
    /* 5. explicit extrapolation past the first/last lock (trend of 3 locks) */
    if (have_lock) {
        size_t first = 0, cnt = 0;
        for (first = 0; first < n && rows[first].status != AXIS_ROW_LOCK; first++) ;
        /* trailing */
        {
            size_t locks[3] = {0, 0, 0};
            for (size_t i = n; i-- > 0 && cnt < 3;)
                if (rows[i].status == AXIS_ROW_LOCK) locks[cnt++] = i;
            if (cnt >= 2) {
                double dz = rows[locks[0]].z - rows[locks[cnt - 1]].z;
                double sy = dz > 0.0 ? (rows[locks[0]].y - rows[locks[cnt - 1]].y) / dz : 0.0;
                double sx = dz > 0.0 ? (rows[locks[0]].x - rows[locks[cnt - 1]].x) / dz : 0.0;
                for (size_t i = locks[0] + 1; i < n; i++) {
                    double d = rows[i].z - rows[locks[0]].z;
                    if (d > AXIS_TABLE_EXTRAP_VOX || rows[i].status != AXIS_ROW_ABSTAIN) break;
                    rows[i].y = rows[locks[0]].y + sy * d;
                    rows[i].x = rows[locks[0]].x + sx * d;
                    rows[i].status = AXIS_ROW_EXTRAP;
                }
            }
        }
        /* leading */
        {
            size_t locks[3] = {0, 0, 0};
            cnt = 0;
            for (size_t i = 0; i < n && cnt < 3; i++)
                if (rows[i].status == AXIS_ROW_LOCK) locks[cnt++] = i;
            if (cnt >= 2) {
                double dz = rows[locks[cnt - 1]].z - rows[locks[0]].z;
                double sy = dz > 0.0 ? (rows[locks[cnt - 1]].y - rows[locks[0]].y) / dz : 0.0;
                double sx = dz > 0.0 ? (rows[locks[cnt - 1]].x - rows[locks[0]].x) / dz : 0.0;
                for (size_t i = locks[0]; i-- > 0;) {
                    double d = rows[i].z - rows[locks[0]].z;   /* negative */
                    if (-d > AXIS_TABLE_EXTRAP_VOX || rows[i].status != AXIS_ROW_ABSTAIN) break;
                    rows[i].y = rows[locks[0]].y + sy * d;
                    rows[i].x = rows[locks[0]].x + sx * d;
                    rows[i].status = AXIS_ROW_EXTRAP;
                }
            }
        }
        (void)first;
    }
    for (size_t i = 0; i < n; i++) if (rows[i].status == AXIS_ROW_LOCK) nlock++;
    return nlock;
}

static const char *at_status_name(int st)
{
    switch (st) {
    case AXIS_ROW_LOCK:   return "lock";
    case AXIS_ROW_INTERP: return "interp";
    case AXIS_ROW_EXTRAP: return "extrap";
    default:              return "abstain";
    }
}

int AxisTrack_write_csv(const char *path, const AxisRow *rows, size_t n,
                        const char *header_lines)
{
    FILE *f = NULL;
    if (path == NULL || rows == NULL) return -1;
    f = fopen(path, "wb");
    if (f == NULL) return -1;
    fprintf(f, "# scroll axis table: z,y,x = the umbilicus centre per z-slab (voxels, (z,y,x) order)\n");
    if (header_lines != NULL) {
        const char *p = header_lines;
        while (*p != '\0') {
            const char *e = strchr(p, '\n');
            size_t len = e != NULL ? (size_t)(e - p) : strlen(p);
            fprintf(f, "# %.*s\n", (int)len, p);
            if (e == NULL) break;
            p = e + 1;
        }
    }
    fprintf(f, "# abstained rows (not written):");
    for (size_t i = 0; i < n; i++)
        if (rows[i].status == AXIS_ROW_ABSTAIN) fprintf(f, " %.0f", rows[i].z);
    fprintf(f, "\n# columns: z,y,x,status,n,cond,sectors,void_ratio,resid_p50,d1,agree,y1,x1,y2,x2,pitch,r_wall,ridges\n");
    for (size_t i = 0; i < n; i++) {
        const AxisRow *r = &rows[i];
        if (r->status == AXIS_ROW_ABSTAIN) continue;
        fprintf(f, "%.1f,%.2f,%.2f,%s,%zu,%.3f,%d,%.3f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.3f,%.1f,%zu\n",
                r->z, r->y, r->x, at_status_name(r->status), r->n_in, r->cond,
                r->sectors, r->void_ratio, r->resid_p50, r->d1, r->agree, r->y1, r->x1,
                r->y2, r->x2, r->pitch, r->r_wall, r->ridges);
    }
    return fclose(f) == 0 ? 0 : -1;
}

int AxisTrack_sense(const AxisRow *rows, size_t n, double *agree)
{
    long pos = 0, neg = 0;
    if (agree != NULL) *agree = 0.0;
    if (rows == NULL) return 0;
    for (size_t i = 0; i < n; i++) {
        if (rows[i].status != AXIS_ROW_LOCK || fabs(rows[i].pitch) < 1.0) continue;
        if (rows[i].pitch > 0.0) pos++; else neg++;
    }
    if (pos + neg == 0) return 0;
    if (agree != NULL)
        *agree = (double)(pos > neg ? pos : neg) / (double)(pos + neg);
    return pos > neg ? 1 : (neg > pos ? -1 : 0);
}

/* ---- selftest -------------------------------------------------------------*/

/* Elliptical Archimedean spiral about (cy, cx): r(theta) = r0 + p*theta/2pi,
 * ellipticity e along y; theta covers [0, cover_deg) per turn (a gap of
 * 360 - cover_deg).  Faces are quads between successive arc samples and
 * z rows, split into triangles. */
static void at_synth_spiral(AxisSlab *s, double cy, double cx, double r0,
                            double pitch, int turns, double ellip,
                            double cover_deg, double z0, double z1,
                            double jitter, uint32_t *rng)
{
    const double dz = 4.0, ds = 3.0;
    double total = 2.0 * M_PI * (double)turns;
    double th = 0.0;
    float prev[2][3];
    int have_prev = 0;
    while (th < total) {
        double r = r0 + pitch * th / (2.0 * M_PI);
        double lap = fmod(th, 2.0 * M_PI) * 180.0 / M_PI;
        double dth = ds / fmax(r, 1.0);
        if (lap < cover_deg) {
            for (double z = z0; z < z1; z += dz) {
                float p[2][3];
                double jy = jitter * ((double)(at_rng(rng) % 2001) / 1000.0 - 1.0);
                double jx = jitter * ((double)(at_rng(rng) % 2001) / 1000.0 - 1.0);
                p[0][0] = (float)z;
                p[0][1] = (float)(cy + ellip * r * cos(th) + jy);
                p[0][2] = (float)(cx + r * sin(th) + jx);
                p[1][0] = (float)(z + dz);
                p[1][1] = p[0][1];
                p[1][2] = p[0][2];
                AxisSlab_add_vertex(s, p[0]);
                if (have_prev && z > z0) {
                    AxisSlab_add_face(s, prev[0], p[0], p[1]);
                    AxisSlab_add_face(s, prev[0], p[1], prev[1]);
                }
                (void)p;
                if (z + dz >= z1) {
                    memcpy(prev, p, sizeof prev);
                }
            }
            have_prev = 1;
        } else {
            have_prev = 0;
        }
        th += dth;
    }
}

int AxisTrack_selftest(void)
{
    int fails = 0;
    uint32_t rng = 12345u;
    const double cy = 3400.0, cx = 2900.0;
    AxisSlab s;
    AxisRow row;
    memset(&row, 0, sizeof row);
    if (AxisSlab_init(&s, 0.0, 64.0, 200000, 500000, 7u) != 0) return 1;
    at_synth_spiral(&s, cy, cx, 60.0, 9.5, 30, 0.85, 300.0, 0.0, 64.0, 1.0, &rng);
    fprintf(stderr, "[selftest] axis_track synthetic slab: %zu faces, %zu verts\n",
            s.nf, s.nv);
    /* primary from a 50-vox-off prior */
    if (AxisTrack_primary(&s, cy + 40.0, cx - 30.0, &row) != 0 ||
        hypot(row.y1 - cy, row.x1 - cx) > 10.0) {
        fprintf(stderr, "  primary: (%.1f, %.1f) want (%.0f, %.0f) n=%zu cond=%.2f resid=%.1f -> FAIL\n",
                row.y1, row.x1, cy, cx, row.n_in, row.cond, row.resid_p50);
        fails++;
    } else
        fprintf(stderr, "  primary: (%.1f, %.1f) n=%zu cond=%.2f resid=%.1f -> ok\n",
                row.y1, row.x1, row.n_in, row.cond, row.resid_p50);
    /* secondary from a 64-vox-off start */
    if (AxisTrack_secondary(&s, cy + 50.0, cx + 40.0, &row) != 0 ||
        hypot(row.y2 - cy, row.x2 - cx) > 5.0 || row.pitch <= 0.0) {
        fprintf(stderr, "  secondary: (%.1f, %.1f) ridges=%zu d1=%.1f pitch=%.2f -> FAIL\n",
                row.y2, row.x2, row.ridges, row.d1, row.pitch);
        fails++;
    } else
        fprintf(stderr, "  secondary: (%.1f, %.1f) ridges=%zu d1=%.1f pitch=%.2f -> ok\n",
                row.y2, row.x2, row.ridges, row.d1, row.pitch);
    /* seed from the grid */
    {
        double sy = 0.0, sx = 0.0, cl_y[16], cl_x[16], cl_d1[16];
        size_t cl_m[16], cl_r[16], nc = 0;
        if (AxisTrack_seed(&s, &sy, &sx, cl_y, cl_x, cl_m, cl_d1, cl_r, 16, &nc) != 0 ||
            hypot(sy - cy, sx - cx) > 15.0) {
            fprintf(stderr, "  seed: (%.1f, %.1f) clusters=%zu -> FAIL\n", sy, sx, nc);
            fails++;
        } else
            fprintf(stderr, "  seed: (%.1f, %.1f) clusters=%zu -> ok\n", sy, sx, nc);
    }
    AxisSlab_dispose(&s);
    /* one-sided slab: both estimators abstain */
    if (AxisSlab_init(&s, 0.0, 64.0, 200000, 500000, 9u) != 0) return fails + 1;
    at_synth_spiral(&s, cy, cx, 60.0, 9.5, 30, 0.85, 110.0, 0.0, 64.0, 1.0, &rng);
    if (AxisTrack_secondary(&s, cy + 5.0, cx - 5.0, &row) == 0) {
        fprintf(stderr, "  one-sided slab: secondary locked (ridges=%zu) -> FAIL\n", row.ridges);
        fails++;
    } else
        fprintf(stderr, "  one-sided slab: secondary abstains -> ok\n");
    AxisSlab_dispose(&s);
    /* finish: slope limit and gap interpolation */
    {
        AxisRow rows[8];
        size_t nlock = 0;
        memset(rows, 0, sizeof rows);
        for (int i = 0; i < 8; i++) {
            rows[i].z = 4384.0 + 64.0 * i;
            rows[i].y1 = rows[i].y2 = 3400.0 + 10.0 * i;
            rows[i].x1 = rows[i].x2 = 2900.0 + 5.0 * i;
            rows[i].primary_ok = rows[i].secondary_ok = 1;
            rows[i].sectors = 12;
            rows[i].pitch = 9.5;
        }
        rows[3].primary_ok = 0;                   /* one abstain: interpolated */
        rows[6].y1 = rows[6].y2 = 3400.0 + 60.0 + 91.0;  /* 91-vox jump: slope 1.4 */
        nlock = AxisTrack_finish(rows, 8);
        if (rows[3].status != AXIS_ROW_INTERP || fabs(rows[3].y - 3430.0) > 1e-6 ||
            rows[6].status != AXIS_ROW_INTERP || fabs(rows[6].y - 3460.0) > 1e-6 ||
            nlock != 6) {
            fprintf(stderr, "  finish: row3=%d y=%.1f row6=%d nlock=%zu -> FAIL\n",
                    rows[3].status, rows[3].y, rows[6].status, nlock);
            fails++;
        } else
            fprintf(stderr, "  finish: gap interpolated, jump abstained then interpolated, %zu locks -> ok\n", nlock);
        if (AxisTrack_sense(rows, 8, NULL) != 1) { fprintf(stderr, "  sense -> FAIL\n"); fails++; }
        /* A sparse polar lock is diagnostic only: it must not veto the
         * normal-line position authority. */
        {
            AxisRow authority;
            memset(&authority, 0, sizeof authority);
            authority.z = 4384.0;
            authority.y1 = 3400.0; authority.x1 = 2900.0;
            authority.y2 = 3500.0; authority.x2 = 3000.0;
            authority.primary_ok = authority.secondary_ok = 1;
            if (AxisTrack_finish(&authority, 1) != 1 ||
                authority.status != AXIS_ROW_LOCK ||
                authority.y != authority.y1 || authority.x != authority.x1) {
                fprintf(stderr, "  finish primary authority -> FAIL\n");
                fails++;
            } else
                fprintf(stderr, "  finish primary authority despite polar disagreement -> ok\n");
        }
        /* csv round trip through AxisWarp_load_csv */
        {
            const char *tmp = "axis_track_selftest.axis.csv";
            AxisWarp w;
            double y = 0.0, x = 0.0;
            AxisWarp_init(&w);
            if (AxisTrack_write_csv(tmp, rows, 8, "selftest\nsecond line") != 0 ||
                AxisWarp_load_csv(&w, tmp) != 0 || w.n < 6) {
                fprintf(stderr, "  csv round trip -> FAIL (n=%zu)\n", w.n);
                fails++;
            } else {
                AxisWarp_eval(&w, 4384.0 + 64.0, &y, &x);
                if (fabs(y - 3410.0) > 1e-6 || fabs(x - 2905.0) > 1e-6) {
                    fprintf(stderr, "  csv eval (%.2f, %.2f) -> FAIL\n", y, x);
                    fails++;
                } else
                    fprintf(stderr, "  csv round trip: %zu rows -> ok\n", w.n);
            }
            AxisWarp_dispose(&w);
            remove(tmp);
        }
    }
    fprintf(stderr, "[selftest] axis_track %s (%d failures)\n",
            fails == 0 ? "PASS" : "FAIL", fails);
    return fails;
}
