#include "tifxyz_recto.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/ves_omp.h"

#include "../common/tiff_io.h"
#include "../common/ves_platform.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define TXR_NS_MAX   ((int)(2.0 * TXR_SEARCH_MAX_VOX / TXR_SAMPLE_STEP + 0.5) + 1)   /* snap profile samples, bound */
#define TXR_KH_MAX   ((int)(2.0 * TXR_SIGMA_MAX / TXR_SAMPLE_STEP + 0.5))            /* smoothing half-width, bound */
#define TXR_PNS      ((int)(2.0 * TXR_PRIOR_RANGE_VOX / TXR_SAMPLE_STEP + 0.5) + 1)  /* prior ray samples */
#define TXR_PKH      ((int)(2.0 * TXR_PRIOR_SIGMA / TXR_SAMPLE_STEP + 0.5))           /* prior smoothing half-width */
#define TXR_RAW_MAX  (TXR_PNS + 2 * TXR_PKH > TXR_NS_MAX + 2 * TXR_KH_MAX ? TXR_PNS + 2 * TXR_PKH : TXR_NS_MAX + 2 * TXR_KH_MAX)
#define TXR_HIST_BINS 6400

/* ------------------------------------------------------------ histograms */

typedef struct TxrHist {
    double lo, hi;
    size_t bin[TXR_HIST_BINS];
    size_t n, over_1p5, over_3;
    double sum_abs, max_abs;
} TxrHist;

static void txr_hist_init(TxrHist *h, double lo, double hi)
{
    memset(h, 0, sizeof *h);
    h->lo = lo;
    h->hi = hi;
}

static void txr_hist_add(TxrHist *h, double v)
{
    double f = (v - h->lo) / (h->hi - h->lo) * (double)TXR_HIST_BINS;
    long b = (long)floor(f);
    if (b < 0) b = 0;
    if (b >= TXR_HIST_BINS) b = TXR_HIST_BINS - 1;
    h->bin[b]++;
    h->n++;
    h->sum_abs += fabs(v);
    if (fabs(v) > h->max_abs) h->max_abs = fabs(v);
    if (fabs(v) > 1.5) h->over_1p5++;
    if (fabs(v) > 3.0) h->over_3++;
}

static void txr_hist_merge(TxrHist *into, const TxrHist *from)
{
    size_t b = 0;
    for (b = 0; b < TXR_HIST_BINS; b++) into->bin[b] += from->bin[b];
    into->n += from->n;
    into->over_1p5 += from->over_1p5;
    into->over_3 += from->over_3;
    into->sum_abs += from->sum_abs;
    if (from->max_abs > into->max_abs) into->max_abs = from->max_abs;
}

static double txr_hist_pct(const TxrHist *h, double p)
{
    size_t want = 0, acc = 0;
    long b = 0;
    if (h->n == 0) return 0.0;
    want = (size_t)floor(p / 100.0 * (double)(h->n - 1));
    for (b = 0; b < TXR_HIST_BINS; b++) {
        acc += h->bin[b];
        if (acc > want) break;
    }
    if (b >= TXR_HIST_BINS) b = TXR_HIST_BINS - 1;
    return h->lo + ((double)b + 0.5) * (h->hi - h->lo) / (double)TXR_HIST_BINS;
}

/* The share of the histogram's mass at or above v. */
static double txr_hist_share_ge(const TxrHist *h, double v)
{
    double f = (v - h->lo) / (h->hi - h->lo) * (double)TXR_HIST_BINS;
    long b0 = (long)floor(f), b = 0;
    size_t acc = 0;
    if (h->n == 0) return 0.0;
    if (b0 < 0) b0 = 0;
    for (b = b0; b < TXR_HIST_BINS; b++) acc += h->bin[b];
    return (double)acc / (double)h->n;
}

static void txr_hist_out(const TxrHist *h, TifxyzRectoHist *o)
{
    memset(o, 0, sizeof *o);
    o->n = h->n;
    if (h->n == 0) return;
    o->p5 = txr_hist_pct(h, 5.0);
    o->p25 = txr_hist_pct(h, 25.0);
    o->p50 = txr_hist_pct(h, 50.0);
    o->p75 = txr_hist_pct(h, 75.0);
    o->p95 = txr_hist_pct(h, 95.0);
    o->mean_abs = h->sum_abs / (double)h->n;
    o->max_abs = h->max_abs;
    o->share_over_1p5 = (double)h->over_1p5 / (double)h->n;
    o->share_over_3 = (double)h->over_3 / (double)h->n;
}

static void txr_hist_q(const TxrHist *h, TifxyzRectoQ *o)
{
    memset(o, 0, sizeof *o);
    o->n = h->n;
    if (h->n == 0) return;
    o->p01 = txr_hist_pct(h, 1.0);
    o->p05 = txr_hist_pct(h, 5.0);
    o->p25 = txr_hist_pct(h, 25.0);
    o->p50 = txr_hist_pct(h, 50.0);
    o->p75 = txr_hist_pct(h, 75.0);
    o->p95 = txr_hist_pct(h, 95.0);
    o->p99 = txr_hist_pct(h, 99.0);
}

static double txr_clamp(double v, double lo, double hi)
{
    return v < lo ? lo : v > hi ? hi : v;
}

/* ------------------------------------------------------------ the grid */

static int txr_valid(const float *pos, size_t plane, size_t i)
{
    return pos[i] > -0.5f && pos[plane + i] > -0.5f && pos[2 * plane + i] > -0.5f;
}

/* The inward unit direction of node (r, c): normalize(dP/dcol x dP/drow) in
 * (z, y, x) order = -N of vc_render (rawtex_tifxyz.h), central differences,
 * one-sided beside an invalid node; 0 when an axis has no valid neighbour. */
static int txr_inward(const float *pos, size_t W, size_t H, size_t r, size_t c, double n[3])
{
    const size_t plane = W * H, i = r * W + c;
    double du[3], dv[3], len = 0.0;
    size_t a = 0, b = 0;
    int k = 0;
    n[0] = n[1] = n[2] = 0.0;
    if (!txr_valid(pos, plane, i)) return 0;
    a = (c > 0 && txr_valid(pos, plane, i - 1)) ? i - 1 : i;
    b = (c + 1 < W && txr_valid(pos, plane, i + 1)) ? i + 1 : i;
    if (a == b) return 0;
    for (k = 0; k < 3; k++) du[k] = (double)pos[(size_t)k * plane + b] - (double)pos[(size_t)k * plane + a];
    a = (r > 0 && txr_valid(pos, plane, i - W)) ? i - W : i;
    b = (r + 1 < H && txr_valid(pos, plane, i + W)) ? i + W : i;
    if (a == b) return 0;
    for (k = 0; k < 3; k++) dv[k] = (double)pos[(size_t)k * plane + b] - (double)pos[(size_t)k * plane + a];
    n[0] = du[1] * dv[2] - du[2] * dv[1];
    n[1] = du[2] * dv[0] - du[0] * dv[2];
    n[2] = du[0] * dv[1] - du[1] * dv[0];
    len = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    if (len < 1e-9) {
        n[0] = n[1] = n[2] = 0.0;
        return 0;
    }
    for (k = 0; k < 3; k++) n[k] /= len;
    return 1;
}

/* The directions of every valid node (3 floats each, 0 when none). */
static size_t txr_directions(const float *pos, size_t W, size_t H, float *dir)
{
    long long n_dir = 0;
    long rr = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:n_dir)
#endif
    for (rr = 0; rr < (long)H; rr++) {
        const size_t r = (size_t)rr;
        size_t c = 0;
        for (c = 0; c < W; c++) {
            const size_t i = r * W + c;
            double n[3];
            if (txr_inward(pos, W, H, r, c, n)) {
                dir[i * 3 + 0] = (float)n[0];
                dir[i * 3 + 1] = (float)n[1];
                dir[i * 3 + 2] = (float)n[2];
                n_dir++;
            } else {
                dir[i * 3 + 0] = dir[i * 3 + 1] = dir[i * 3 + 2] = 0.0f;
            }
        }
    }
    return (size_t)n_dir;
}

/* ------------------------------------------------------------ parameters and profiles */

/* The parameters a run reads its profiles with: the prior's values, or the constants. */
typedef struct TxrPrm {
    double range, step, sigma, min_contrast, min_edge, edge_full, thick_min, thick_max, max_step, max_total;
    int levels;                                 /* 1: a band ends where a trough falls to air_level (+ the floor) and a
                                                 * crest is a band from band_min_level; 0: the prominence rule */
    double air_level, band_min_level;
    int ns, kh;
    double kernel[2 * TXR_KH_MAX + 1];
} TxrPrm;

static void txr_kernel(double sigma, double step, int kh, double *kernel)
{
    const double s = sigma / step;
    double sum = 0.0;
    int k = 0;
    for (k = -kh; k <= kh; k++) {
        kernel[k + kh] = exp(-0.5 * (double)k * (double)k / (s * s));
        sum += kernel[k + kh];
    }
    for (k = 0; k < 2 * kh + 1; k++) kernel[k] /= sum;
}

static void txr_prm_init(TxrPrm *p, const TifxyzRectoPrior *prior)
{
    memset(p, 0, sizeof *p);
    p->range = TXR_SEARCH_VOX;
    p->step = TXR_SAMPLE_STEP;
    p->sigma = TXR_SMOOTH_SIGMA;
    p->min_contrast = TXR_MIN_CONTRAST;
    p->min_edge = TXR_MIN_EDGE;
    p->edge_full = TXR_EDGE_FULL;
    p->thick_min = TXR_MIDLINE_MIN_THICK;
    p->thick_max = TXR_MIDLINE_MAX_THICK;
    p->max_step = TXR_MAX_STEP_VOX;
    p->max_total = TXR_MAX_TOTAL_VOX;
    if (prior != NULL && prior->rays_with_recto >= 100) {
        p->range = prior->search_vox;
        p->max_total = prior->max_total_vox;
        p->min_contrast = prior->min_contrast;
        p->min_edge = prior->min_edge;
        p->edge_full = prior->edge_full;
        p->sigma = prior->sigma;
        p->thick_min = prior->thick_min;
        p->thick_max = prior->thick_max;
        p->levels = 1;
        p->air_level = prior->air_level;
        p->band_min_level = prior->band_min_level;
    }
    p->range = txr_clamp(p->range, 2.0, TXR_SEARCH_MAX_VOX);
    p->sigma = txr_clamp(p->sigma, 0.25, TXR_SIGMA_MAX);
    p->ns = 2 * (int)floor(p->range / p->step + 0.5) + 1;             /* odd: a sample at t = 0 */
    if (p->ns > TXR_NS_MAX) p->ns = TXR_NS_MAX - ((TXR_NS_MAX + 1) % 2);
    if (p->ns < 5) p->ns = 5;
    p->range = 0.5 * (double)(p->ns - 1) * p->step;
    p->kh = (int)floor(2.0 * p->sigma / p->step + 0.5);
    if (p->kh > TXR_KH_MAX) p->kh = TXR_KH_MAX;
    if (p->kh < 1) p->kh = 1;
    txr_kernel(p->sigma, p->step, p->kh, p->kernel);
}

static double txr_t(int s, double range, double step)
{
    return -range + (double)s * step;
}

/* The smoothed CT profile q[0..ns) at t = -range + s * step along n from p;
 * 0 when a sample has no RAW datum. */
static int txr_profile(CubeTable *ct, const double p[3], const double n[3], double range, double step, int ns, int kh,
                       const double *kernel, double *q)
{
    double raw[TXR_RAW_MAX];
    int s = 0, k = 0;
    for (s = 0; s < ns + 2 * kh; s++) {
        const double t = -range + (double)(s - kh) * step;
        const double v = sample_trilinear(ct, p[0] + t * n[0], p[1] + t * n[1], p[2] + t * n[2]);
        if (v < 0.0) return 0;
        raw[s] = v;
    }
    for (s = 0; s < ns; s++) {
        double acc = 0.0;
        for (k = -kh; k <= kh; k++) acc += kernel[k + kh] * raw[s + kh + k];
        q[s] = acc;
    }
    return 1;
}

/* The steepest descent (sign -1) or ascent (+1) of d[] over samples lo..hi, at
 * a sub-sample position from the parabola through the three samples around
 * it; *strength = |d| there. */
static int txr_steepest(const double *d, int ns, double range, double step, int lo, int hi, int sign, double *t_out,
                        double *strength)
{
    int s = 0, sm = -1;
    double best = 0.0, off = 0.0;
    if (lo < 1) lo = 1;
    if (hi > ns - 2) hi = ns - 2;
    for (s = lo; s <= hi; s++) {
        const double v = (double)sign * d[s];
        if (sm < 0 || v > best) {
            best = v;
            sm = s;
        }
    }
    if (sm < 0) return 0;
    if (sm - 1 >= 1 && sm + 1 <= ns - 2) {
        const double a = (double)sign * d[sm - 1], b = (double)sign * d[sm], c = (double)sign * d[sm + 1];
        const double denom = a - 2.0 * b + c;
        if (denom < -1e-9) {
            off = 0.5 * (a - c) / denom;
            if (off > 0.5) off = 0.5;
            if (off < -0.5) off = -0.5;
        }
    }
    *t_out = txr_t(sm, range, step) + off * step;
    *strength = best;
    return 1;
}

/* Local extrema with the plateau rule: one crest per flat top (its last
 * sample), one trough per flat bottom; the first sample counts only when it
 * differs from the next, the last one as the top or the bottom the profile
 * rose or fell into (a flat run ending there is that run's last sample). */
#define TXR_EQ_GREY 1e-6     /* smoothed samples of a flat run differ by rounding only: equal within this */

static int txr_cmp(double a, double b)
{
    return a > b + TXR_EQ_GREY ? 1 : a < b - TXR_EQ_GREY ? -1 : 0;
}

static void txr_extrema(const double *q, int ns, int *is_max, int *is_min)
{
    int s = 0;
    for (s = 0; s < ns; s++) is_max[s] = is_min[s] = 0;
    for (s = 1; s + 1 < ns; s++) {
        is_max[s] = txr_cmp(q[s], q[s - 1]) >= 0 && txr_cmp(q[s], q[s + 1]) > 0;
        is_min[s] = txr_cmp(q[s], q[s - 1]) <= 0 && txr_cmp(q[s], q[s + 1]) < 0;
    }
    if (ns >= 2) {
        is_max[0] = txr_cmp(q[0], q[1]) > 0;
        is_min[0] = txr_cmp(q[0], q[1]) < 0;
        s = ns - 2;
        while (s > 0 && txr_cmp(q[s], q[ns - 1]) == 0) s--;
        is_max[ns - 1] = txr_cmp(q[ns - 1], q[s]) > 0;
        is_min[ns - 1] = txr_cmp(q[ns - 1], q[s]) < 0;
    }
}

/* ------------------------------------------------------------ the snap's target */

/* The target offset along inward of one node (t > 0 inward), 1 when found, 0
 * when the node abstains.
 *
 * A crest counts as a BAND when it stands the gate (min_contrast) above the
 * profile's darkest sample (a flat plateau's rounding or CT noise does not).
 * The node belongs to the band whose crest is NEAREST to it -- the band's body,
 * not an edge: a node in the dark gap between two wraps is assigned to the
 * wrap it is closer to.  Only that band offers a recto candidate: the first
 * trough inward of its crest that stands the gate below it (shallower dips --
 * the seam between the sheet's two fibre layers, noise -- are walked over; a
 * rise above the crest ends the walk), and the steepest descent between the
 * two, if at least min_edge.  Otherwise the node ABSTAINS and follows its
 * neighbours (*beyond = 1 when the band is there but its recto is out of
 * reach).  It never takes another band's recto: choosing the recto edge
 * nearest the node pulled 2.5% of cand5's nodes 3-6 vox OUTWARD onto the
 * neighbouring outer wrap's recto when their own band lay inward of them with
 * its recto out of reach, and snapped others to the fibre-layer seam
 * (2026-10-10).  kind 1: the edge is inward of the node (the node sits in or
 * outward of its band); 2: outward (the node sits in the gap inside its band).
 * midline: the midpoint of the recto and the band's verso edge (the steepest
 * ascent from the trough outward of the crest), for a band thick_min..max. */
static int txr_node_target(CubeTable *ct, const double p[3], const double n[3], int midline, double inset,
                           const TxrPrm *prm, double *out_t, double *out_conf, int *out_kind, int *beyond)
{
    double q[TXR_NS_MAX], d[TXR_NS_MAX];
    int is_max[TXR_NS_MAX], is_min[TXR_NS_MAX];
    const int ns = prm->ns, s0 = (ns - 1) / 2;
    double qmin = 0.0, tr = 0.0, g = 0.0, t_cand = 0.0;
    int s = 0, sc = -1, st = -1, sv = -1;

    *beyond = 0;
    if (!txr_profile(ct, p, n, prm->range, prm->step, ns, prm->kh, prm->kernel, q)) return 0;
    for (s = 1; s + 1 < ns; s++) d[s] = (q[s + 1] - q[s - 1]) / (2.0 * prm->step);
    d[0] = d[1];
    d[ns - 1] = d[ns - 2];
    qmin = q[0];
    for (s = 1; s < ns; s++)
        if (q[s] < qmin) qmin = q[s];
    txr_extrema(q, ns, is_max, is_min);
    /* a band: a crest above the scroll's band_min_level (the prior), or standing the gate above the darkest
     * sample (the constants); its end: a trough fallen to the air level (+ the noise floor), or standing the gate
     * below the crest */
    for (s = 0; s < ns; s++) {
        const int band = prm->levels ? q[s] >= prm->band_min_level : q[s] - qmin >= prm->min_contrast;
        if (is_max[s] && band && (sc < 0 || abs(s - s0) < abs(sc - s0))) sc = s;
    }
    if (sc < 0) return 0;
    for (s = sc + 1; s < ns; s++) {
        const int gap = prm->levels ? q[s] <= prm->air_level + TXR_PRIOR_FLOOR && q[sc] - q[s] >= TXR_PRIOR_FLOOR
                                    : q[sc] - q[s] >= prm->min_contrast;
        if (q[s] > q[sc]) break;
        if (is_min[s] && gap) { st = s; break; }
    }
    if (st < 0) {
        *beyond = 1;
        return 0;
    }
    if (!txr_steepest(d, ns, prm->range, prm->step, sc, st, -1, &tr, &g) || g < prm->min_edge) return 0;
    t_cand = tr - inset;
    if (midline) {
        double tv = 0.0, gv = 0.0, thick = 0.0;
        for (s = sc - 1; s >= 0; s--) {
            const int gap = prm->levels ? q[s] <= prm->air_level + TXR_PRIOR_FLOOR && q[sc] - q[s] >= TXR_PRIOR_FLOOR
                                        : q[sc] - q[s] >= prm->min_contrast;
            if (q[s] > q[sc]) break;
            if (is_min[s] && gap) { sv = s; break; }
        }
        if (sv < 0) return 0;
        if (!txr_steepest(d, ns, prm->range, prm->step, sv, sc, 1, &tv, &gv) || gv < prm->min_edge) return 0;
        thick = tr - tv;
        if (thick < prm->thick_min || thick > prm->thick_max) return 0;
        t_cand = 0.5 * (tr + tv);
    }
    *out_conf = txr_clamp((g - prm->min_edge) / (prm->edge_full - prm->min_edge), 0.0, 1.0);
    *out_kind = t_cand >= 0.0 ? 1 : 2;
    *out_t = t_cand;
    return 1;
}

/* Targets and weights of every node; counts[0] supported, [1] the edge inward of the node, [2] outward, [3]
 * abstained because the node's band has its recto beyond the range. */
/* Is node i in the counted columns [cc0, cc1)? */
static int txr_in_core(size_t i, size_t W, size_t cc0, size_t cc1)
{
    const size_t c = i % W;
    return c >= cc0 && c < cc1;
}

/* Nodes with a nonzero flag / an inward direction within the counted columns. */
static size_t txr_count_flags(const uint8_t *f, size_t W, size_t H, size_t cc0, size_t cc1)
{
    size_t r = 0, c = 0, cnt = 0;
    for (r = 0; r < H; r++)
        for (c = cc0; c < cc1; c++)
            if (f[r * W + c]) cnt++;
    return cnt;
}

static size_t txr_count_dir(const float *dir, size_t W, size_t H, size_t cc0, size_t cc1)
{
    size_t r = 0, c = 0, cnt = 0;
    for (r = 0; r < H; r++)
        for (c = cc0; c < cc1; c++) {
            const size_t i = r * W + c;
            if (dir[i * 3 + 0] != 0.0f || dir[i * 3 + 1] != 0.0f || dir[i * 3 + 2] != 0.0f) cnt++;
        }
    return cnt;
}

static void txr_measure(CubeTable *ct, const float *pos, const float *dir, size_t W, size_t H, int midline, double inset,
                        const TxrPrm *prm, float *target, float *weight, uint8_t *kind, size_t cc0, size_t cc1,
                        size_t counts[4])
{
    const size_t plane = W * H;
    long long n_sup = 0, n_in = 0, n_out = 0, n_beyond = 0;
    long rr = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 8) reduction(+:n_sup,n_in,n_out,n_beyond)
#endif
    for (rr = 0; rr < (long)H; rr++) {
        const size_t r = (size_t)rr;
        size_t c = 0;
        for (c = 0; c < W; c++) {
            const size_t i = r * W + c;
            double p[3], n[3], t = 0.0, conf = 0.0;
            int k = 0, beyond = 0;
            target[i] = 0.0f;
            kind[i] = 0;
            weight[i] = txr_valid(pos, plane, i) ? (float)TXR_ANCHOR_W : 0.0f;
            if (dir[i * 3 + 0] == 0.0f && dir[i * 3 + 1] == 0.0f && dir[i * 3 + 2] == 0.0f) continue;
            p[0] = pos[i];
            p[1] = pos[plane + i];
            p[2] = pos[2 * plane + i];
            n[0] = dir[i * 3 + 0];
            n[1] = dir[i * 3 + 1];
            n[2] = dir[i * 3 + 2];
            if (!txr_node_target(ct, p, n, midline, inset, prm, &t, &conf, &k, &beyond)) {
                if (beyond && c >= cc0 && c < cc1) n_beyond++;
                continue;
            }
            target[i] = (float)t;
            weight[i] = (float)(TXR_DATA_W * (0.35 + 0.65 * conf));
            kind[i] = (uint8_t)k;
            if (c < cc0 || c >= cc1) continue;
            n_sup++;
            if (k == 1) n_in++; else n_out++;
        }
    }
    counts[0] = (size_t)n_sup;
    counts[1] = (size_t)n_in;
    counts[2] = (size_t)n_out;
    counts[3] = (size_t)n_beyond;
}

/* ------------------------------------------------------------ the median filter */

/* The k-th smallest of a[0..n) (0-based); a is permuted. */
static float txr_select(float *a, int n, int k)
{
    int lo = 0, hi = n - 1;
    while (lo < hi) {
        const float pivot = a[(lo + hi) / 2];
        int i = lo, j = hi;
        while (i <= j) {
            while (a[i] < pivot) i++;
            while (a[j] > pivot) j--;
            if (i <= j) {
                const float t = a[i];
                a[i] = a[j];
                a[j] = t;
                i++;
                j--;
            }
        }
        if (k <= j) hi = j;
        else if (k >= i) lo = i;
        else return a[k];
    }
    return a[lo];
}

#define TXR_MEDIAN_BUF ((2 * TXR_MEDIAN_HALF + 1) * (2 * TXR_MEDIAN_HALF + 1) > \
                        (2 * (TXR_MEDIAN_HALF_SMOOTH / TXR_MEDIAN_STRIDE_SMOOTH) + 1) * (2 * (TXR_MEDIAN_HALF_SMOOTH / TXR_MEDIAN_STRIDE_SMOOTH) + 1) \
                        ? (2 * TXR_MEDIAN_HALF + 1) * (2 * TXR_MEDIAN_HALF + 1) \
                        : (2 * (TXR_MEDIAN_HALF_SMOOTH / TXR_MEDIAN_STRIDE_SMOOTH) + 1) * (2 * (TXR_MEDIAN_HALF_SMOOTH / TXR_MEDIAN_STRIDE_SMOOTH) + 1))

/* The targets through a median over the supported nodes of each valid node's
 * window.  The window is POINT-SYMMETRIC: an offset counts only together with
 * its mirror (both valid, both supported), so the median is unbiased for a
 * linear field at the grid's border or beside a hole too (a one-sided window
 * leans toward the interior, and the outer loop's fixed point then never
 * reaches the border node).  The weight is the window's mean data weight
 * times its supported share; fkind 0 none, else the node's own kind or 3
 * when only the window supports it.  Counts the filtered nodes. */
static size_t txr_median_filter(const float *target, const float *weight, const uint8_t *kind, const uint8_t *valid,
                                size_t W, size_t H, int smooth, float *ftarget, float *fweight, uint8_t *fkind)
{
    const long half = smooth ? TXR_MEDIAN_HALF_SMOOTH : TXR_MEDIAN_HALF;
    const long stride = smooth ? TXR_MEDIAN_STRIDE_SMOOTH : 1;
    long long n_filtered = 0;
    long rr = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 4) reduction(+:n_filtered)
#endif
    for (rr = 0; rr < (long)H; rr++) {
        float buf[TXR_MEDIAN_BUF + 2];
        const size_t r = (size_t)rr;
        size_t c = 0;
        for (c = 0; c < W; c++) {
            const size_t i = r * W + c;
            long dr = 0, dc = 0;
            int n = 0, nwin = 0;
            double wsum = 0.0;
            fkind[i] = 0;
            ftarget[i] = 0.0f;
            fweight[i] = valid[i] ? (float)TXR_ANCHOR_W : 0.0f;
            if (!valid[i]) continue;
            nwin = 1;
            if (kind[i] != 0) {
                buf[n++] = target[i];
                wsum += weight[i];
            }
            for (dr = 0; dr <= half; dr += stride) {
                const long r1 = (long)r + dr, r2 = (long)r - dr;
                if (r1 >= (long)H || r2 < 0) continue;
                for (dc = (dr == 0 ? stride : -half); dc <= half; dc += stride) {
                    const long c1 = (long)c + dc, c2 = (long)c - dc;
                    size_t j1 = 0, j2 = 0;
                    if (c1 < 0 || c1 >= (long)W || c2 < 0 || c2 >= (long)W) continue;
                    j1 = (size_t)r1 * W + (size_t)c1;
                    j2 = (size_t)r2 * W + (size_t)c2;
                    if (!valid[j1] || !valid[j2]) continue;
                    nwin += 2;
                    if (kind[j1] == 0 || kind[j2] == 0) continue;
                    buf[n++] = target[j1];
                    buf[n++] = target[j2];
                    wsum += weight[j1] + weight[j2];
                }
            }
            if (n == 0 || (double)n < TXR_MEDIAN_MIN_SUPPORT * (double)nwin) continue;
            ftarget[i] = txr_select(buf, n, (n - 1) / 2);
            fweight[i] = (float)((wsum / (double)n) * ((double)n / (double)nwin));
            fkind[i] = kind[i] != 0 ? kind[i] : (uint8_t)3;
            n_filtered++;
        }
    }
    return (size_t)n_filtered;
}

/* The mean CT profile at -TXR_PROFILE_HALF..+HALF vox along inward over every
 * TXR_PROFILE_STRIDE-th node with a direction and a complete profile. */
static size_t txr_mean_profile(CubeTable *ct, const float *pos, const float *dir, size_t W, size_t H,
                               size_t cc0, size_t cc1, double out[TXR_PROFILE_N])
{
    const size_t plane = W * H;
    double sum[TXR_PROFILE_N];
    size_t r = 0, c = 0, n = 0;
    int k = 0;
    for (k = 0; k < TXR_PROFILE_N; k++) sum[k] = 0.0;
    for (r = 0; r < H; r += TXR_PROFILE_STRIDE) {
        for (c = 0; c < W; c += TXR_PROFILE_STRIDE) {
            const size_t i = r * W + c;
            double v[TXR_PROFILE_N];
            int ok = 1;
            if (c < cc0 || c >= cc1) continue;
            if (dir[i * 3 + 0] == 0.0f && dir[i * 3 + 1] == 0.0f && dir[i * 3 + 2] == 0.0f) continue;
            for (k = 0; k < TXR_PROFILE_N && ok; k++) {
                const double t = (double)(k - TXR_PROFILE_HALF);
                v[k] = sample_trilinear(ct, pos[i] + t * dir[i * 3 + 0], pos[plane + i] + t * dir[i * 3 + 1],
                                        pos[2 * plane + i] + t * dir[i * 3 + 2]);
                if (v[k] < 0.0) ok = 0;
            }
            if (!ok) continue;
            for (k = 0; k < TXR_PROFILE_N; k++) sum[k] += v[k];
            n++;
        }
    }
    for (k = 0; k < TXR_PROFILE_N; k++) out[k] = n > 0 ? sum[k] / (double)n : NAN;
    return n;
}

/* ------------------------------------------------------------ the solve */

/* y = (W + mu L) x over the compact valid nodes. */
static void txr_matvec(const int32_t *nodes, const int32_t *cidx, size_t nv, size_t W, size_t H, const float *w,
                       const double *x, double *y)
{
    long long k = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (k = 0; k < (long long)nv; k++) {
        const size_t i = (size_t)nodes[k], r = i / W, c = i % W;
        double acc = (double)w[i] * x[k];
        int32_t j = -1;
        if (c > 0 && (j = cidx[i - 1]) >= 0) acc += TXR_SMOOTH_MU * (x[k] - x[j]);
        if (c + 1 < W && (j = cidx[i + 1]) >= 0) acc += TXR_SMOOTH_MU * (x[k] - x[j]);
        if (r > 0 && (j = cidx[i - W]) >= 0) acc += TXR_SMOOTH_MU * (x[k] - x[j]);
        if (r + 1 < H && (j = cidx[i + W]) >= 0) acc += TXR_SMOOTH_MU * (x[k] - x[j]);
        y[k] = acc;
    }
}

static double txr_diag(const int32_t *nodes, const int32_t *cidx, size_t k, size_t W, size_t H, const float *w)
{
    const size_t i = (size_t)nodes[k], r = i / W, c = i % W;
    double deg = 0.0, diag = 0.0;
    if (c > 0 && cidx[i - 1] >= 0) deg += 1.0;
    if (c + 1 < W && cidx[i + 1] >= 0) deg += 1.0;
    if (r > 0 && cidx[i - W] >= 0) deg += 1.0;
    if (r + 1 < H && cidx[i + W] >= 0) deg += 1.0;
    diag = (double)w[i] + TXR_SMOOTH_MU * deg;
    return diag < 1e-12 ? 1e-12 : diag;
}

static double txr_dot(const double *a, const double *b, size_t n)
{
    double s = 0.0;
    long long i = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:s)
#endif
    for (i = 0; i < (long long)n; i++) s += a[i] * b[i];
    return s;
}

/* Jacobi-preconditioned CG for (W + mu L) d = W target over the compact
 * nodes, from the d given (zero, or a warm start).  Returns the iterations
 * taken; *resid the relative residual reached. */
static int txr_pcg(const int32_t *nodes, const int32_t *cidx, size_t nv, size_t W, size_t H, const float *w,
                   const float *target, double *d, double *r, double *z, double *p, double *Ap, double *resid)
{
    double bnorm = 0.0, rz = 0.0, rz_new = 0.0, alpha = 0.0, beta = 0.0, pAp = 0.0, rnorm = 0.0;
    long long k = 0;
    int it = 0;

    if (nv == 0) {
        *resid = 0.0;
        return 0;
    }
    txr_matvec(nodes, cidx, nv, W, H, w, d, Ap);               /* r = b - A d */
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:bnorm)
#endif
    for (k = 0; k < (long long)nv; k++) {
        const size_t i = (size_t)nodes[k];
        const double b = (double)w[i] * (double)target[i];
        r[k] = b - Ap[k];
        bnorm += b * b;
    }
    bnorm = sqrt(bnorm);
    if (bnorm <= 0.0) {
        *resid = 0.0;
        return 0;
    }
    rnorm = sqrt(txr_dot(r, r, nv));
    if (rnorm <= TXR_PCG_TOL * bnorm) {
        *resid = rnorm / bnorm;
        return 0;
    }
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (k = 0; k < (long long)nv; k++) {
        z[k] = r[k] / txr_diag(nodes, cidx, (size_t)k, W, H, w);
        p[k] = z[k];
    }
    rz = txr_dot(r, z, nv);
    for (it = 0; it < TXR_PCG_MAX_ITERS; it++) {
        txr_matvec(nodes, cidx, nv, W, H, w, p, Ap);
        pAp = txr_dot(p, Ap, nv);
        if (pAp <= 0.0) break;
        alpha = rz / pAp;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (k = 0; k < (long long)nv; k++) {
            d[k] += alpha * p[k];
            r[k] -= alpha * Ap[k];
        }
        rnorm = sqrt(txr_dot(r, r, nv));
        if (rnorm <= TXR_PCG_TOL * bnorm) {
            it++;
            break;
        }
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (k = 0; k < (long long)nv; k++) z[k] = r[k] / txr_diag(nodes, cidx, (size_t)k, W, H, w);
        rz_new = txr_dot(r, z, nv);
        beta = rz_new / rz;
        rz = rz_new;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (k = 0; k < (long long)nv; k++) p[k] = z[k] + beta * p[k];
    }
    *resid = rnorm / bnorm;
    return it;
}

/* ------------------------------------------------------------ the prior */

/* One ray's reading of the wrap structure around its node. */
#define TXR_RAY_TROUGHS 24
typedef struct TxrRay {
    int has_band, has_recto, has_verso, has_next_in, has_next_out, has_seam;
    double t_crest, t_recto, t_verso, spacing_in, spacing_out, gap_width_in, gap_depth_in, gap_depth_out,
           seam_depth_max, slope, width;
    int n_troughs;                       /* every pruned trough on the ray: its half-depth width and whether it reads as a gap */
    double trough_width[TXR_RAY_TROUGHS];
    int trough_is_gap[TXR_RAY_TROUGHS];
} TxrRay;

/* A trough's width at half its depth below the lower flanking crest. */
static double txr_trough_width(const double *q, int ns, double step, int trough, double depth)
{
    const double level = q[trough] + 0.5 * depth;
    int a = trough, b = trough;
    double tl = 0.0, tr = 0.0;
    while (a > 0 && q[a] < level) a--;
    while (b + 1 < ns && q[b] < level) b++;
    tl = a < trough && q[a + 1] != q[a] ? (double)a + (level - q[a]) / (q[a + 1] - q[a]) : (double)a;
    tr = b > trough && q[b - 1] != q[b] ? (double)b - (level - q[b]) / (q[b - 1] - q[b]) : (double)b;
    return (tr - tl) * step;
}

/* The extrema of a ray as an alternating sequence of crests and troughs,
 * pruned of pairs closer than TXR_PRIOR_FLOOR (noise); returns their count. */
static int txr_ray_extrema(const double *q, int ns, int *ext, int *is_crest)
{
    int is_max[TXR_PNS], is_min[TXR_PNS];
    int ne = 0, s = 0;
    txr_extrema(q, ns, is_max, is_min);
    for (s = 0; s < ns; s++) {
        int mx = is_max[s], mn = is_min[s];
        if (!mx && !mn) continue;
        if (mx && mn) continue;                     /* a lone flat sample between equals: neither */
        if (ne > 0 && is_crest[ne - 1] == mx) {     /* two of a kind in a row: keep the better one */
            if ((mx && q[s] > q[ext[ne - 1]]) || (mn && q[s] < q[ext[ne - 1]])) ext[ne - 1] = s;
            continue;
        }
        ext[ne] = s;
        is_crest[ne] = mx;
        ne++;
    }
    for (;;) {
        int i = 0, best = -1;
        double bd = 1e30;
        if (ne < 3) break;
        for (i = 0; i + 1 < ne; i++) {
            const double dd = fabs(q[ext[i]] - q[ext[i + 1]]);
            if (dd < bd) { bd = dd; best = i; }
        }
        if (best < 0 || bd >= TXR_PRIOR_FLOOR) break;
        memmove(ext + best, ext + best + 2, (size_t)(ne - best - 2) * sizeof(int));
        memmove(is_crest + best, is_crest + best + 2, (size_t)(ne - best - 2) * sizeof(int));
        ne -= 2;
    }
    return ne;
}

/* The 10%..90% distance of the descent from level a (at sample sa) to level b (at sample sb > sa). */
static double txr_edge_width(const double *q, int sa, int sb, double step)
{
    const double a = q[sa], b = q[sb], l10 = a - 0.1 * (a - b), l90 = a - 0.9 * (a - b);
    int s = 0, s10 = -1, s90 = -1;
    for (s = sa; s <= sb; s++) {
        if (s10 < 0 && q[s] <= l10) s10 = s;
        if (s90 < 0 && q[s] <= l90) { s90 = s; break; }
    }
    if (s10 < 0 || s90 < 0 || s90 <= s10) return step;
    return (double)(s90 - s10) * step;
}

/* The prior's ray: raw samples where RAW exists, smoothed where the kernel
 * fits, over the contiguous span around the node's own sample; *lo..*hi are
 * the smoothed samples' indices (q[s] at t = -range + s * step).  0 when the
 * node itself has no datum.  A ray leaving the volume or the masked scroll is
 * truncated, not lost. */
static int txr_ray_profile(CubeTable *ct, const double p[3], const double n[3], double range, double step, int ns,
                           int kh, const double *kernel, double *q, int *lo, int *hi)
{
    double raw[TXR_RAW_MAX];
    int have[TXR_RAW_MAX];
    const int c = kh + (ns - 1) / 2, nraw = ns + 2 * kh;
    int s = 0, a = 0, b = 0, k = 0;
    for (s = 0; s < nraw; s++) {
        const double t = -range + (double)(s - kh) * step;
        const double v = sample_trilinear(ct, p[0] + t * n[0], p[1] + t * n[1], p[2] + t * n[2]);
        have[s] = v >= 0.0;
        raw[s] = have[s] ? v : 0.0;
    }
    if (!have[c]) return 0;
    a = c;
    while (a > 0 && have[a - 1]) a--;
    b = c;
    while (b + 1 < nraw && have[b + 1]) b++;
    *lo = a;                          /* q[s] reads raw[s .. s + 2 kh] */
    *hi = b - 2 * kh;
    if (*hi < *lo) return 0;
    for (s = *lo; s <= *hi; s++) {
        double acc = 0.0;
        for (k = -kh; k <= kh; k++) acc += kernel[k + kh] * raw[s + kh + k];
        q[s] = acc;
    }
    return 1;
}

/* q[0..ns): the ray at t = (s - s0) * step = -range + s * step. */
/* The levels pass: every pruned crest's and trough's level into the histograms (the scroll's air and band levels
 * come from all rays together, not from each ray's own darkest sample: a ray that crosses a void makes every gap
 * on it look shallow). */
static void txr_ray_levels(const double *q, int ns, TxrHist *crests, TxrHist *troughs)
{
    int ext[TXR_PNS], is_crest[TXR_PNS];
    int ne = 0, j = 0;
    double qmin = q[0], qmax = q[0];
    int s = 0;
    for (s = 1; s < ns; s++) {
        if (q[s] < qmin) qmin = q[s];
        if (q[s] > qmax) qmax = q[s];
    }
    if (qmax - qmin < TXR_PRIOR_FLOOR) return;
    ne = txr_ray_extrema(q, ns, ext, is_crest);
    for (j = 0; j < ne; j++) txr_hist_add(is_crest[j] ? crests : troughs, q[ext[j]]);
}

/* Otsu's split of a level histogram: the threshold that maximizes the between-class variance; *separability =
 * its share of the total variance (1 = two clean modes, ~0 = one mode). */
static double txr_otsu(const TxrHist *h, double *separability)
{
    double total = 0.0, sum = 0.0, sumb = 0.0, wb = 0.0, best = -1.0, thr = h->lo, mean = 0.0, var = 0.0;
    size_t b = 0;
    if (h->n == 0) { *separability = 0.0; return h->lo; }
    for (b = 0; b < TXR_HIST_BINS; b++) {
        const double v = h->lo + ((double)b + 0.5) * (h->hi - h->lo) / (double)TXR_HIST_BINS;
        total += (double)h->bin[b];
        sum += v * (double)h->bin[b];
    }
    mean = sum / total;
    for (b = 0; b < TXR_HIST_BINS; b++) {
        const double v = h->lo + ((double)b + 0.5) * (h->hi - h->lo) / (double)TXR_HIST_BINS;
        var += (v - mean) * (v - mean) * (double)h->bin[b];
    }
    var /= total;
    for (b = 0; b < TXR_HIST_BINS; b++) {
        const double v = h->lo + ((double)b + 0.5) * (h->hi - h->lo) / (double)TXR_HIST_BINS;
        double wf = 0.0, mb = 0.0, mf = 0.0, between = 0.0;
        wb += (double)h->bin[b];
        sumb += v * (double)h->bin[b];
        if (wb <= 0.0) continue;
        wf = total - wb;
        if (wf <= 0.0) break;
        mb = sumb / wb;
        mf = (sum - sumb) / wf;
        between = wb * wf * (mb - mf) * (mb - mf) / (total * total);
        if (between > best) {
            best = between;
            thr = v;
        }
    }
    *separability = var > 0.0 ? best / var : 0.0;
    return thr;
}

static void txr_ray_read(const double *q, int ns, int s0, double range, double step, double air_level, TxrRay *o)
{
    double d[TXR_PNS];
    int ext[TXR_PNS], is_crest[TXR_PNS];
    double qmin = q[0], qmax = q[0], contrast = 0.0;
    int s = 0, ne = 0, k = -1, j = 0;

    memset(o, 0, sizeof *o);
    for (s = 1; s < ns; s++) {
        if (q[s] < qmin) qmin = q[s];
        if (q[s] > qmax) qmax = q[s];
    }
    contrast = qmax - qmin;
    if (contrast < TXR_PRIOR_FLOOR) return;
    for (s = 1; s + 1 < ns; s++) d[s] = (q[s + 1] - q[s - 1]) / (2.0 * step);
    d[0] = d[1];
    d[ns - 1] = d[ns - 2];
    ne = txr_ray_extrema(q, ns, ext, is_crest);
    for (j = 0; j < ne; j++) {
        if (is_crest[j] || o->n_troughs >= TXR_RAY_TROUGHS) continue;
        {
            const double cl = j > 0 ? q[ext[j - 1]] : 1e30, cr = j + 1 < ne ? q[ext[j + 1]] : 1e30;
            const double lower = cl < cr ? cl : cr;
            const double depth = lower < 1e29 ? lower - q[ext[j]] : 0.0;
            if (depth < TXR_PRIOR_FLOOR) continue;
            o->trough_width[o->n_troughs] = txr_trough_width(q, ns, step, ext[j], depth);
            o->trough_is_gap[o->n_troughs] = q[ext[j]] <= air_level + TXR_PRIOR_FLOOR;
            o->n_troughs++;
        }
    }
    for (j = 0; j < ne; j++)
        if (is_crest[j] && (k < 0 || abs(ext[j] - s0) < abs(ext[k] - s0))) k = j;
    if (k < 0) return;
    o->has_band = 1;
    o->t_crest = txr_t(ext[k], range, step);
    /* inward over the sequence: troughs at k + 1, k + 3, ...; a trough fallen to the air level is a GAP and ends
     * the band; a higher one is a dip inside it (the fibre-layer seam) and the band goes on */
    for (j = k + 1; j < ne; j += 2) {
        const int trough = ext[j], crest_before = ext[j - 1];
        const double depth_before = q[crest_before] - q[trough];
        const double depth = j + 1 < ne ? (q[ext[j + 1]] < q[crest_before] ? q[ext[j + 1]] - q[trough] : depth_before)
                                        : depth_before;
        if (q[trough] <= air_level + TXR_PRIOR_FLOOR) {
            double tr = 0.0, g = 0.0;
            if (txr_steepest(d, ns, range, step, crest_before, trough, -1, &tr, &g)) {
                o->has_recto = 1;
                o->t_recto = tr;
                o->slope = g;
                o->width = txr_edge_width(q, crest_before, trough, step);
                o->gap_depth_in = depth;
                if (j + 1 < ne) {
                    double tv = 0.0, gv = 0.0;
                    o->has_next_in = 1;
                    o->spacing_in = txr_t(ext[j + 1], range, step) - o->t_crest;
                    if (txr_steepest(d, ns, range, step, trough, ext[j + 1], 1, &tv, &gv)) o->gap_width_in = tv - tr;
                }
            }
            break;
        }
        o->has_seam = 1;
        if (depth > o->seam_depth_max) o->seam_depth_max = depth;
    }
    for (j = k - 1; j >= 0; j -= 2) {
        const int trough = ext[j], crest_after = ext[j + 1];
        const double depth_after = q[crest_after] - q[trough];
        const double depth = j - 1 >= 0 ? (q[ext[j - 1]] < q[crest_after] ? q[ext[j - 1]] - q[trough] : depth_after)
                                        : depth_after;
        if (q[trough] <= air_level + TXR_PRIOR_FLOOR) {
            double tv = 0.0, gv = 0.0;
            if (txr_steepest(d, ns, range, step, trough, crest_after, 1, &tv, &gv)) {
                o->has_verso = 1;
                o->t_verso = tv;
                o->gap_depth_out = depth;
                if (j - 1 >= 0) {
                    o->has_next_out = 1;
                    o->spacing_out = o->t_crest - txr_t(ext[j - 1], range, step);
                }
            }
            break;
        }
        o->has_seam = 1;
        if (depth > o->seam_depth_max) o->seam_depth_max = depth;
    }
}

enum { PH_RECTO, PH_RECTO_ABS, PH_SPACING, PH_THICK, PH_GAPW, PH_GAPD, PH_SEAMD, PH_SLOPE, PH_WIDTH, PH_CREST, PH_TROUGH,
       PH_VERSO_ABS, PH_TWID_GAP, PH_TWID_SEAM, PH_TWID_ALL, PH_COUNT };

int TifxyzRecto_prior(Arena_T arena, CubeTable *ct, const float *pos, size_t W, size_t H, int stride,
                      TifxyzRectoPrior *out)
{
    const size_t plane = W * H;
    const double range = TXR_PRIOR_RANGE_VOX, step = TXR_SAMPLE_STEP;
    const int ns = TXR_PNS, kh = TXR_PKH;
    double kernel[2 * TXR_PKH + 1];
    int nthreads = 1, t = 0;
    TxrHist *hists = NULL;
    Arena_Mark mark;
    long long rays = 0, with_band = 0, with_recto = 0, with_next = 0, with_seam = 0;
    long rr = 0;
    const double qlo[PH_COUNT] = { -TXR_PRIOR_RANGE_VOX, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    const double qhi[PH_COUNT] = { TXR_PRIOR_RANGE_VOX, TXR_PRIOR_RANGE_VOX, 2.0 * TXR_PRIOR_RANGE_VOX, TXR_PRIOR_RANGE_VOX,
                                   TXR_PRIOR_RANGE_VOX, 256.0, 256.0, 400.0, 8.0, 256.0, 256.0, TXR_PRIOR_RANGE_VOX, 16.0, 16.0, 16.0 };
    double air_level = 0.0, band_level = 0.0, separability = 0.0;
    int pass = 0;

    assert(arena && ct && out);
    memset(out, 0, sizeof *out);
    out->range_vox = range;
    out->stride = stride > 0 ? stride : TXR_PRIOR_STRIDE;
    /* the defaults, until the rays say otherwise */
    out->search_vox = TXR_SEARCH_VOX;
    out->max_total_vox = TXR_MAX_TOTAL_VOX;
    out->min_contrast = TXR_MIN_CONTRAST;
    out->min_edge = TXR_MIN_EDGE;
    out->edge_full = TXR_EDGE_FULL;
    out->sigma = TXR_SMOOTH_SIGMA;
    out->thick_min = TXR_MIDLINE_MIN_THICK;
    out->thick_max = TXR_MIDLINE_MAX_THICK;
    if (W == 0 || H == 0 || pos == NULL) return 0;
    txr_kernel(TXR_PRIOR_SIGMA, step, kh, kernel);
#ifdef _OPENMP
    nthreads = omp_get_max_threads();
#endif
    if (nthreads < 1) nthreads = 1;
    mark = Arena_save(arena);
    hists = (TxrHist *)ARENA_CALLOC(arena, (size_t)nthreads * PH_COUNT, sizeof(TxrHist));
    for (t = 0; t < nthreads; t++) {
        int h = 0;
        for (h = 0; h < PH_COUNT; h++) txr_hist_init(&hists[(size_t)t * PH_COUNT + (size_t)h], qlo[h], qhi[h]);
    }
    /* pass 0: the levels of every crest and trough on every ray -> the scroll's air and band levels; pass 1: the
     * wraps' structure read against them */
    for (pass = 0; pass < 2; pass++) {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 4) reduction(+:rays,with_band,with_recto,with_next,with_seam)
#endif
    for (rr = 0; rr < (long)H; rr += (long)out->stride) {
        const size_t r = (size_t)rr;
        size_t c = 0;
        int tid = 0;
        TxrHist *hh = NULL;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        hh = hists + (size_t)tid * PH_COUNT;
        for (c = 0; c < W; c += (size_t)out->stride) {
            const size_t i = r * W + c;
            const int s0 = (ns - 1) / 2, min_side = (int)floor(8.0 / step + 0.5);   /* 8 vox each side at least */
            double p[3], n[3], q[TXR_PNS];
            int lo = 0, hi = 0;
            TxrRay ray;
            if (!txr_inward(pos, W, H, r, c, n)) continue;
            p[0] = pos[i];
            p[1] = pos[plane + i];
            p[2] = pos[2 * plane + i];
            if (!txr_ray_profile(ct, p, n, range, step, ns, kh, kernel, q, &lo, &hi)) continue;
            if (s0 - lo < min_side || hi - s0 < min_side) continue;
            if (pass == 0) {
                txr_ray_levels(q + lo, hi - lo + 1, &hh[PH_CREST], &hh[PH_TROUGH]);
                continue;
            }
            rays++;
            txr_ray_read(q + lo, hi - lo + 1, s0 - lo, (double)(s0 - lo) * step, step, air_level, &ray);
            {
                int tj = 0;
                for (tj = 0; tj < ray.n_troughs; tj++) {
                    txr_hist_add(&hh[PH_TWID_ALL], ray.trough_width[tj]);
                    txr_hist_add(ray.trough_is_gap[tj] ? &hh[PH_TWID_GAP] : &hh[PH_TWID_SEAM], ray.trough_width[tj]);
                }
            }
            if (!ray.has_band) continue;
            with_band++;
            if (ray.has_seam) {
                with_seam++;
                txr_hist_add(&hh[PH_SEAMD], ray.seam_depth_max);
            }
            if (ray.has_recto) {
                with_recto++;
                txr_hist_add(&hh[PH_RECTO], ray.t_recto);
                txr_hist_add(&hh[PH_RECTO_ABS], fabs(ray.t_recto));
                txr_hist_add(&hh[PH_GAPD], ray.gap_depth_in);
                txr_hist_add(&hh[PH_SLOPE], ray.slope);
                txr_hist_add(&hh[PH_WIDTH], ray.width);
                if (ray.has_verso) txr_hist_add(&hh[PH_THICK], ray.t_recto - ray.t_verso);
                if (ray.has_next_in) {
                    txr_hist_add(&hh[PH_SPACING], ray.spacing_in);
                    txr_hist_add(&hh[PH_GAPW], ray.gap_width_in);
                }
            }
            if (ray.has_verso) {
                txr_hist_add(&hh[PH_GAPD], ray.gap_depth_out);
                txr_hist_add(&hh[PH_VERSO_ABS], fabs(ray.t_verso));
            }
            if (ray.has_next_out) txr_hist_add(&hh[PH_SPACING], ray.spacing_out);
            if (ray.has_next_in || ray.has_next_out) with_next++;
        }
    }
    if (pass == 0) {
        for (t = 1; t < nthreads; t++) {
            txr_hist_merge(&hists[PH_CREST], &hists[(size_t)t * PH_COUNT + PH_CREST]);
            txr_hist_merge(&hists[PH_TROUGH], &hists[(size_t)t * PH_COUNT + PH_TROUGH]);
        }
        txr_hist_q(&hists[PH_TROUGH], &out->trough_level);
        txr_hist_q(&hists[PH_CREST], &out->crest_level);
        air_level = txr_otsu(&hists[PH_TROUGH], &separability);
        band_level = out->crest_level.p50;
        if (hists[PH_TROUGH].n < 100 || air_level >= band_level) {   /* no structure: the darkest troughs */
            air_level = out->trough_level.p05;
            separability = 0.0;
        }
        out->air_level = air_level;
        out->band_level = band_level;
        out->band_min_level = 0.5 * (air_level + band_level);
        out->level_separability = separability;
    }
    }
    for (t = 1; t < nthreads; t++) {
        int h = 0;
        for (h = 0; h < PH_COUNT; h++)
            if (h != PH_CREST && h != PH_TROUGH) txr_hist_merge(&hists[h], &hists[(size_t)t * PH_COUNT + (size_t)h]);
    }
    out->rays = (size_t)rays;
    out->rays_with_band = (size_t)with_band;
    out->rays_with_recto = (size_t)with_recto;
    out->rays_with_next = (size_t)with_next;
    out->rays_with_seam = (size_t)with_seam;
    txr_hist_q(&hists[PH_RECTO], &out->recto_offset);
    txr_hist_q(&hists[PH_SPACING], &out->spacing);
    txr_hist_q(&hists[PH_THICK], &out->thickness);
    txr_hist_q(&hists[PH_GAPW], &out->gap_width);
    txr_hist_q(&hists[PH_GAPD], &out->gap_depth);
    txr_hist_q(&hists[PH_SEAMD], &out->seam_depth);
    txr_hist_q(&hists[PH_SLOPE], &out->slope);
    txr_hist_q(&hists[PH_WIDTH], &out->edge_width);
    txr_hist_q(&hists[PH_TWID_GAP], &out->gap_trough_width);
    txr_hist_q(&hists[PH_TWID_SEAM], &out->seam_trough_width);
    if (hists[PH_TWID_ALL].n >= 100) {
        double sep = 0.0;
        out->trough_width_split = txr_otsu(&hists[PH_TWID_ALL], &sep);
        out->trough_width_separability = sep;
        /* how the level classes fall either side of the width split */
        out->trough_width_error = 0.5 * ((1.0 - txr_hist_share_ge(&hists[PH_TWID_GAP], out->trough_width_split))
                                         + txr_hist_share_ge(&hists[PH_TWID_SEAM], out->trough_width_split));
    }
    if (with_recto >= 100) {
        /* the band gate: the depth that best separates the dips inside bands from the gaps between them */
        double gate = 0.0, best_err = 1e30, g = 0.0;
        double reach = 0.0, reach_recto = 0.0, bound = TXR_SEARCH_MAX_VOX, w = 0.0, s2 = 0.0;
        TifxyzRectoQ abs_recto, abs_verso;
        txr_hist_q(&hists[PH_RECTO_ABS], &abs_recto);
        txr_hist_q(&hists[PH_VERSO_ABS], &abs_verso);
        if (hists[PH_SEAMD].n >= 50 && hists[PH_GAPD].n >= 50) {
            for (g = 2.0; g <= 200.0; g += 1.0) {
                const double err = txr_hist_share_ge(&hists[PH_SEAMD], g) + (1.0 - txr_hist_share_ge(&hists[PH_GAPD], g));
                if (err < best_err) { best_err = err; gate = g; }
            }
            out->gate_error = 0.5 * best_err;
        } else {
            gate = 0.3 * out->gap_depth.p05;
            if (gate < TXR_MIN_CONTRAST) gate = TXR_MIN_CONTRAST;
            out->gate_error = 0.0;
        }
        out->min_contrast = txr_clamp(gate, 6.0, 80.0);
        out->min_edge = txr_clamp(out->slope.p05, 2.0, 30.0);
        out->edge_full = txr_clamp(out->slope.p50, out->min_edge + 4.0, 120.0);
        w = out->edge_width.p50 / 2.56;                        /* a Gaussian-blurred step's 10..90% width is 2.56 sigma */
        s2 = w * w - TXR_PRIOR_SIGMA * TXR_PRIOR_SIGMA;        /* less the rays' own smoothing */
        out->sigma = txr_clamp(s2 > 0.0 ? sqrt(s2) : 0.4, 0.4, TXR_SIGMA_MAX);
        /* the range must reach the recto (p95 of its offset) with the gap trough beyond it -- the edge's width,
         * the smoothing's spread and a sample of slack -- and, for the midline, the verso too; the wrap spacing
         * bounds it, and the recto's reach is what matters when both cannot fit */
        {
            const double margin = out->edge_width.p50 + 2.0 * out->sigma + 1.0;
            reach_recto = abs_recto.p95 + margin;
            reach = abs_verso.n >= 100 && abs_verso.p95 + margin > reach_recto ? abs_verso.p95 + margin : reach_recto;
        }
        if (reach < 4.0) reach = 4.0;
        if (with_next >= 100) bound = 0.9 * out->spacing.p05;
        if (bound > TXR_SEARCH_MAX_VOX) bound = TXR_SEARCH_MAX_VOX;
        out->crowded = reach_recto > bound;
        out->search_vox = txr_clamp(reach < bound ? reach : bound, 3.0, TXR_SEARCH_MAX_VOX);
        out->max_total_vox = out->search_vox - 1.0;
        if (out->max_total_vox < 2.0) out->max_total_vox = 2.0;
        out->thick_min = txr_clamp(out->thickness.p01 - 0.5, 1.0, 6.0);
        out->thick_max = txr_clamp(out->thickness.p99 + 1.0, 8.0, 24.0);
    }
    Arena_restore(arena, mark);
    return 0;
}

/* ------------------------------------------------------------ the run */

int TifxyzRecto_run(Arena_T arena, CubeTable *ct, float *pos, size_t W, size_t H,
                    const TifxyzRectoOpts *opts, const TifxyzRectoPrior *prior, TifxyzRectoStats *st,
                    float *target_first, float *target_last, float *target_first_med,
                    float *target_last_med, float *disp)
{
    const size_t n = W * H, plane = W * H;
    Arena_Mark mark;
    uint8_t *valid = NULL, *kind = NULL, *fkind = NULL;
    int32_t *nodes = NULL, *cidx = NULL;
    float *dir = NULL, *target = NULL, *weight = NULL, *total = NULL, *ftarget = NULL, *fweight = NULL, *w_round = NULL;
    double *d = NULL, *r = NULL, *z = NULL, *p = NULL, *Ap = NULL;
    TxrHist *hist = NULL;
    TxrPrm prm;
    size_t i = 0, n_valid = 0, n_core = 0, cc0 = 0, cc1 = 0, counts[4] = { 0, 0, 0, 0 };
    int outer = 0, midline = 0, verbose = 0, smooth = 0;
    double t0 = 0.0, t1 = 0.0, inset = 0.0;

    assert(arena && ct && st);
    memset(st, 0, sizeof *st);
    if (W == 0 || H == 0 || pos == NULL) return 0;
    if (n > (size_t)0x7fffffff) return -1;
    midline = opts != NULL && opts->midline;
    verbose = opts != NULL && opts->verbose;
    smooth = opts != NULL && opts->smooth;
    inset = opts != NULL && opts->inset ? TXR_RECTO_INSET_VOX : 0.0;
    cc0 = opts != NULL && opts->count_c1 > opts->count_c0 ? opts->count_c0 : 0;
    cc1 = opts != NULL && opts->count_c1 > opts->count_c0 ? opts->count_c1 : W;
    if (cc1 > W) cc1 = W;
    if (cc0 > cc1) cc0 = cc1;
    txr_prm_init(&prm, prior);
    st->prm_search = prm.range;
    st->prm_max_total = prm.max_total;
    st->prm_min_contrast = prm.min_contrast;
    st->prm_min_edge = prm.min_edge;
    st->prm_edge_full = prm.edge_full;
    st->prm_sigma = prm.sigma;
    st->prm_thick_min = prm.thick_min;
    st->prm_thick_max = prm.thick_max;
    st->prm_levels = prm.levels;
    st->prm_air_level = prm.air_level;
    st->prm_band_min_level = prm.band_min_level;
    t0 = ves_clock_sec();
    st->n_nodes = n;
    for (i = 0; i < n; i++) {
        if (!txr_valid(pos, plane, i)) continue;
        n_valid++;
        if (txr_in_core(i, W, cc0, cc1)) n_core++;
    }
    st->n_valid = n_core;
    if (target_first != NULL) for (i = 0; i < n; i++) target_first[i] = NAN;
    if (target_last != NULL) for (i = 0; i < n; i++) target_last[i] = NAN;
    if (target_first_med != NULL) for (i = 0; i < n; i++) target_first_med[i] = NAN;
    if (target_last_med != NULL) for (i = 0; i < n; i++) target_last_med[i] = NAN;
    if (disp != NULL) for (i = 0; i < n; i++) disp[i] = NAN;
    if (n_valid == 0) return 0;

    mark = Arena_save(arena);
    valid = (uint8_t *)ARENA_CALLOC(arena, n, 1);
    kind = (uint8_t *)ARENA_CALLOC(arena, n, 1);
    fkind = (uint8_t *)ARENA_CALLOC(arena, n, 1);
    dir = (float *)ARENA_CALLOC(arena, n, 3 * sizeof(float));
    target = (float *)ARENA_CALLOC(arena, n, sizeof(float));
    weight = (float *)ARENA_CALLOC(arena, n, sizeof(float));
    ftarget = (float *)ARENA_CALLOC(arena, n, sizeof(float));
    fweight = (float *)ARENA_CALLOC(arena, n, sizeof(float));
    w_round = (float *)ARENA_CALLOC(arena, n, sizeof(float));
    total = (float *)ARENA_CALLOC(arena, n, sizeof(float));
    nodes = (int32_t *)ARENA_CALLOC(arena, n_valid, sizeof(int32_t));
    cidx = (int32_t *)ARENA_CALLOC(arena, n, sizeof(int32_t));
    d = (double *)ARENA_CALLOC(arena, n_valid, sizeof(double));
    r = (double *)ARENA_CALLOC(arena, n_valid, sizeof(double));
    z = (double *)ARENA_CALLOC(arena, n_valid, sizeof(double));
    p = (double *)ARENA_CALLOC(arena, n_valid, sizeof(double));
    Ap = (double *)ARENA_CALLOC(arena, n_valid, sizeof(double));
    hist = (TxrHist *)ARENA_CALLOC(arena, 1, sizeof(TxrHist));
    {
        size_t k = 0;
        for (i = 0; i < n; i++) {
            valid[i] = (uint8_t)txr_valid(pos, plane, i);
            cidx[i] = valid[i] ? (int32_t)k : -1;
            if (valid[i]) nodes[k++] = (int32_t)i;
        }
    }

    for (outer = 0; outer < TXR_OUTER_ITERS; outer++) {
        double ta = ves_clock_sec(), tb = 0.0, resid = 0.0, step_sum = 0.0;
        long long moved = 0, capped = 0;
        long long ii = 0;
        txr_directions(pos, W, H, dir);
        if (outer == 0) st->n_dir = txr_count_dir(dir, W, H, cc0, cc1);
        txr_measure(ct, pos, dir, W, H, midline, inset, &prm, target, weight, kind, cc0, cc1, counts);
        txr_median_filter(target, weight, kind, valid, W, H, smooth, ftarget, fweight, fkind);
        if (outer == 0) st->n_filtered_first = txr_count_flags(fkind, W, H, cc0, cc1);
        tb = ves_clock_sec();
        st->sec_measure += tb - ta;
        if (outer == 0) {
            st->n_supported_first = counts[0];
            st->n_in_band_first = counts[1];
            st->n_in_gap_first = counts[2];
            st->n_beyond_first = counts[3];
            txr_hist_init(hist, -16.0, 16.0);
            for (i = 0; i < n; i++) {
                if (kind[i] == 0) continue;
                if (txr_in_core(i, W, cc0, cc1)) txr_hist_add(hist, (double)target[i]);
                if (target_first != NULL) target_first[i] = target[i];
            }
            txr_hist_out(hist, &st->target_first);
            txr_hist_init(hist, -16.0, 16.0);
            for (i = 0; i < n; i++) {
                if (fkind[i] == 0) continue;
                if (txr_in_core(i, W, cc0, cc1)) txr_hist_add(hist, (double)ftarget[i]);
                if (target_first_med != NULL) target_first_med[i] = ftarget[i];
            }
            txr_hist_out(hist, &st->target_first_med);
            st->profile_n = txr_mean_profile(ct, pos, dir, W, H, cc0, cc1, st->profile_first);
        }
        /* least squares first, then Cauchy-reweighted by the residual to that smooth field (IRLS): a target far
         * from what its neighbourhood agrees on -- CT noise, another band -- fades continuously */
        {
            int round = 0;
            memcpy(w_round, fweight, n * sizeof(float));
            for (ii = 0; ii < (long long)n_valid; ii++) d[ii] = 0.0;
            for (round = 0; round < TXR_IRLS_ROUNDS; round++) {
                if (round > 0) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
                    for (ii = 0; ii < (long long)n_valid; ii++) {
                        const size_t k = (size_t)nodes[ii];
                        double res = 0.0;
                        if (fkind[k] == 0) continue;
                        res = (d[ii] - (double)ftarget[k]) / TXR_ROBUST_SCALE_VOX;
                        w_round[k] = (float)((double)fweight[k] / (1.0 + res * res));
                    }
                }
                st->pcg_iters[outer] += txr_pcg(nodes, cidx, n_valid, W, H, w_round, ftarget, d, r, z, p, Ap, &resid);
            }
        }
        st->pcg_resid[outer] = resid;
        st->sec_solve += ves_clock_sec() - tb;
        /* the move: capped per iteration and from the input grid */
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:moved,capped,step_sum)
#endif
        for (ii = 0; ii < (long long)n_valid; ii++) {
            const size_t k = (size_t)nodes[ii];
            double step = 0.0, newtotal = 0.0;
            if (dir[k * 3 + 0] == 0.0f && dir[k * 3 + 1] == 0.0f && dir[k * 3 + 2] == 0.0f) continue;
            step = d[ii];
            if (step > prm.max_step) step = prm.max_step;
            if (step < -prm.max_step) step = -prm.max_step;
            newtotal = (double)total[k] + step;
            if (newtotal > prm.max_total) { newtotal = prm.max_total; capped += txr_in_core(k, W, cc0, cc1); }
            if (newtotal < -prm.max_total) { newtotal = -prm.max_total; capped += txr_in_core(k, W, cc0, cc1); }
            step = newtotal - (double)total[k];
            if (fabs(step) < 1e-7) continue;
            pos[k] = (float)((double)pos[k] + step * (double)dir[k * 3 + 0]);
            pos[plane + k] = (float)((double)pos[plane + k] + step * (double)dir[k * 3 + 1]);
            pos[2 * plane + k] = (float)((double)pos[2 * plane + k] + step * (double)dir[k * 3 + 2]);
            total[k] = (float)newtotal;
            moved++;
            step_sum += fabs(step);
        }
        st->n_capped_total += (size_t)capped;
        st->mean_abs_step[outer] = n_valid > 0 ? step_sum / (double)n_valid : 0.0;
        st->outer_iters = outer + 1;
        if (verbose)
            fprintf(stderr, "[tifxyz_recto] iter %d: supported %zu (edge inward %zu, outward %zu; band beyond reach "
                            "%zu) of %zu valid, filtered %zu; pcg %d its resid %.1e; moved %lld, mean |step| %.3f "
                            "vox, capped %lld (%.1fs)\n",
                    outer + 1, counts[0], counts[1], counts[2], counts[3], n_valid,
                    outer == 0 ? st->n_filtered_first : (size_t)0, st->pcg_iters[outer], resid, moved,
                    st->mean_abs_step[outer], capped, ves_clock_sec() - ta);
        if (moved == 0 || st->mean_abs_step[outer] < TXR_STEP_TOL) break;
    }

    /* measure again on the moved grid: the residual */
    t1 = ves_clock_sec();
    txr_directions(pos, W, H, dir);
    txr_measure(ct, pos, dir, W, H, midline, inset, &prm, target, weight, kind, cc0, cc1, counts);
    txr_median_filter(target, weight, kind, valid, W, H, smooth, ftarget, fweight, fkind);
    st->n_filtered_last = txr_count_flags(fkind, W, H, cc0, cc1);
    st->sec_measure += ves_clock_sec() - t1;
    st->n_supported_last = counts[0];
    st->n_in_band_last = counts[1];
    st->n_in_gap_last = counts[2];
    txr_hist_init(hist, -16.0, 16.0);
    for (i = 0; i < n; i++) {
        if (kind[i] == 0) continue;
        if (txr_in_core(i, W, cc0, cc1)) txr_hist_add(hist, (double)target[i]);
        if (target_last != NULL) target_last[i] = target[i];
    }
    txr_hist_out(hist, &st->target_last);
    txr_hist_init(hist, -16.0, 16.0);
    for (i = 0; i < n; i++) {
        if (fkind[i] == 0) continue;
        if (txr_in_core(i, W, cc0, cc1)) txr_hist_add(hist, (double)ftarget[i]);
        if (target_last_med != NULL) target_last_med[i] = ftarget[i];
    }
    txr_hist_out(hist, &st->target_last_med);
    txr_mean_profile(ct, pos, dir, W, H, cc0, cc1, st->profile_last);
    txr_hist_init(hist, -16.0, 16.0);
    for (i = 0; i < n; i++) {
        if (!valid[i]) continue;
        if (txr_in_core(i, W, cc0, cc1)) txr_hist_add(hist, (double)total[i]);
        if (disp != NULL) disp[i] = total[i];
    }
    txr_hist_out(hist, &st->disp);
    txr_hist_init(hist, -16.0, 16.0);
    for (i = 0; i < n; i++) {
        const size_t rr = i / W, c = i % W;
        if (!valid[i] || c < cc0 || c >= cc1) continue;
        if (c + 1 < W && valid[i + 1]) txr_hist_add(hist, fabs((double)total[i] - (double)total[i + 1]));
        if (rr + 1 < H && valid[i + W]) txr_hist_add(hist, fabs((double)total[i] - (double)total[i + W]));
    }
    st->slope_p99 = hist->n > 0 ? txr_hist_pct(hist, 99.0) : 0.0;
    st->slope_max = hist->max_abs;
    st->sec_total = ves_clock_sec() - t0;
    Arena_restore(arena, mark);
    return 0;
}

/* ------------------------------------------------------------ windows */

/* The combined statistics of a run in windows: the distributions are rebuilt from each window's core rasters (the
 * summaries of TifxyzRectoHist do not merge); counts add; profiles are weighted by their node counts. */
enum { TXA_FIRST, TXA_LAST, TXA_FIRST_MED, TXA_LAST_MED, TXA_DISP, TXA_SLOPE, TXA_COUNT };

typedef struct TxrAcc {
    TxrHist h[TXA_COUNT];
    double prof_first[TXR_PROFILE_N], prof_last[TXR_PROFILE_N];
    size_t prof_n;
    double step_w[TXR_OUTER_ITERS];
    size_t step_n[TXR_OUTER_ITERS];
} TxrAcc;

static void txr_acc_add(TxrAcc *a, TifxyzRectoStats *tot, const TifxyzRectoStats *st, const float *const ras[5], size_t ww,
                        size_t H, size_t c0, size_t c1)
{
    size_t r = 0, c = 0;
    int k = 0;
    tot->n_nodes += H * (c1 - c0);
    tot->n_valid += st->n_valid;
    tot->n_dir += st->n_dir;
    tot->n_supported_first += st->n_supported_first;
    tot->n_in_band_first += st->n_in_band_first;
    tot->n_in_gap_first += st->n_in_gap_first;
    tot->n_beyond_first += st->n_beyond_first;
    tot->n_supported_last += st->n_supported_last;
    tot->n_in_band_last += st->n_in_band_last;
    tot->n_in_gap_last += st->n_in_gap_last;
    tot->n_capped_total += st->n_capped_total;
    tot->n_filtered_first += st->n_filtered_first;
    tot->n_filtered_last += st->n_filtered_last;
    if (st->outer_iters > tot->outer_iters) tot->outer_iters = st->outer_iters;
    for (k = 0; k < st->outer_iters && k < TXR_OUTER_ITERS; k++) {
        a->step_w[k] += st->mean_abs_step[k] * (double)st->n_valid;
        a->step_n[k] += st->n_valid;
        tot->pcg_iters[k] += st->pcg_iters[k];
        if (st->pcg_resid[k] > tot->pcg_resid[k]) tot->pcg_resid[k] = st->pcg_resid[k];
    }
    if (st->profile_n > 0) {
        for (k = 0; k < TXR_PROFILE_N; k++) {
            if (isfinite(st->profile_first[k])) a->prof_first[k] += st->profile_first[k] * (double)st->profile_n;
            if (isfinite(st->profile_last[k])) a->prof_last[k] += st->profile_last[k] * (double)st->profile_n;
        }
        a->prof_n += st->profile_n;
    }
    tot->sec_measure += st->sec_measure;
    tot->sec_solve += st->sec_solve;
    tot->sec_total += st->sec_total;
    for (r = 0; r < H; r++) {
        for (c = c0; c < c1; c++) {
            const size_t i = r * ww + c;
            for (k = 0; k < 5; k++)
                if (isfinite(ras[k][i])) txr_hist_add(&a->h[k], (double)ras[k][i]);
            if (!isfinite(ras[4][i])) continue;
            /* the slope over the grid edges leaving a core node (a right neighbour in the halo is the window's own) */
            if (c + 1 < ww && isfinite(ras[4][i + 1])) txr_hist_add(&a->h[TXA_SLOPE], fabs((double)ras[4][i] - (double)ras[4][i + 1]));
            if (r + 1 < H && isfinite(ras[4][i + ww])) txr_hist_add(&a->h[TXA_SLOPE], fabs((double)ras[4][i] - (double)ras[4][i + ww]));
        }
    }
}

static void txr_acc_out(const TxrAcc *a, TifxyzRectoStats *tot)
{
    int k = 0;
    txr_hist_out(&a->h[TXA_FIRST], &tot->target_first);
    txr_hist_out(&a->h[TXA_LAST], &tot->target_last);
    txr_hist_out(&a->h[TXA_FIRST_MED], &tot->target_first_med);
    txr_hist_out(&a->h[TXA_LAST_MED], &tot->target_last_med);
    txr_hist_out(&a->h[TXA_DISP], &tot->disp);
    tot->slope_p99 = a->h[TXA_SLOPE].n > 0 ? txr_hist_pct(&a->h[TXA_SLOPE], 99.0) : 0.0;
    tot->slope_max = a->h[TXA_SLOPE].max_abs;
    for (k = 0; k < TXR_OUTER_ITERS; k++) tot->mean_abs_step[k] = a->step_n[k] > 0 ? a->step_w[k] / (double)a->step_n[k] : 0.0;
    for (k = 0; k < TXR_PROFILE_N; k++) {
        tot->profile_first[k] = a->prof_n > 0 ? a->prof_first[k] / (double)a->prof_n : NAN;
        tot->profile_last[k] = a->prof_n > 0 ? a->prof_last[k] / (double)a->prof_n : NAN;
    }
    tot->profile_n = a->prof_n;
}

int TifxyzRecto_run_windows(Arena_T arena, CubeTable *ct, const float *pos_in, float *pos_out, size_t W, size_t H,
                            size_t col0, size_t col1, size_t window_cols, size_t halo,
                            const TifxyzRectoOpts *opts, const TifxyzRectoPrior *prior, TifxyzRectoStats *st,
                            TifxyzRectoWindowFn fn, void *fn_ctx)
{
    const size_t n = W * H;
    TxrAcc *acc = NULL;
    TxrPrm prm;
    size_t nwin = 0, k = 0;
    int h = 0;

    assert(arena && ct && st);
    memset(st, 0, sizeof *st);
    txr_prm_init(&prm, prior);
    st->prm_search = prm.range;
    st->prm_max_total = prm.max_total;
    st->prm_min_contrast = prm.min_contrast;
    st->prm_min_edge = prm.min_edge;
    st->prm_edge_full = prm.edge_full;
    st->prm_sigma = prm.sigma;
    st->prm_thick_min = prm.thick_min;
    st->prm_thick_max = prm.thick_max;
    st->prm_levels = prm.levels;
    st->prm_air_level = prm.air_level;
    st->prm_band_min_level = prm.band_min_level;
    for (k = 0; k < TXR_PROFILE_N; k++) st->profile_first[k] = st->profile_last[k] = NAN;
    if (W == 0 || H == 0 || pos_in == NULL || pos_out == NULL) return 0;
    if (n > (size_t)0x7fffffff / 3) return -1;
    memcpy(pos_out, pos_in, 3 * n * sizeof(float));
    if (col1 > W) col1 = W;
    if (col0 >= col1) return 0;
    if (window_cols == 0 || window_cols > col1 - col0) window_cols = col1 - col0;
    nwin = (col1 - col0 + window_cols - 1) / window_cols;
    acc = (TxrAcc *)ARENA_CALLOC(arena, 1, sizeof(TxrAcc));
    for (h = 0; h < TXA_COUNT; h++) txr_hist_init(&acc->h[h], -16.0, 16.0);
    for (k = 0; k < nwin; k++) {
        const size_t c0 = col0 + k * window_cols, c1 = c0 + window_cols < col1 ? c0 + window_cols : col1;
        const size_t lc0 = c0 > halo ? c0 - halo : 0, lc1 = c1 + halo < W ? c1 + halo : W, ww = lc1 - lc0, wn = ww * H;
        Arena_Mark mark = Arena_save(arena);
        TifxyzRectoOpts wo;
        TifxyzRectoStats ws;
        float *wpos = NULL, *ras[5] = { NULL, NULL, NULL, NULL, NULL };
        size_t r = 0, c = 0, any = 0;
        int p = 0;
        for (r = 0; r < H && !any; r++)
            for (c = c0; c < c1 && !any; c++)
                any = txr_valid(pos_in, n, r * W + c);
        if (!any) continue;
        if (opts != NULL) wo = *opts;
        else memset(&wo, 0, sizeof wo);
        wo.count_c0 = c0 - lc0;
        wo.count_c1 = c1 - lc0;
        wpos = (float *)ARENA_ALLOC(arena, 3 * wn * sizeof(float));
        for (p = 0; p < 5; p++) ras[p] = (float *)ARENA_ALLOC(arena, wn * sizeof(float));
        for (p = 0; p < 3; p++)
            for (r = 0; r < H; r++)
                memcpy(wpos + (size_t)p * wn + r * ww, pos_in + (size_t)p * n + r * W + lc0, ww * sizeof(float));
        if (TifxyzRecto_run(arena, ct, wpos, ww, H, &wo, prior, &ws, ras[0], ras[1], ras[2], ras[3], ras[4]) != 0) {
            Arena_restore(arena, mark);
            return -1;
        }
        for (p = 0; p < 3; p++)
            for (r = 0; r < H; r++)
                memcpy(pos_out + (size_t)p * n + r * W + c0, wpos + (size_t)p * wn + r * ww + (c0 - lc0),
                       (c1 - c0) * sizeof(float));
        txr_acc_add(acc, st, &ws, (const float *const *)ras, ww, H, c0 - lc0, c1 - lc0);
        st->windows++;
        if (fn != NULL) fn(fn_ctx, k, nwin, lc0, ww, H, c0 - lc0, c1 - lc0, &ws, ras[0], ras[1], ras[2], ras[3], ras[4]);
        Arena_restore(arena, mark);
    }
    txr_acc_out(acc, st);
    return 0;
}

/* ------------------------------------------------------------ self-test */

#define TXR_CHECK(cond, msg) \
    do { \
        if (!(cond)) { \
            fprintf(stderr, "  FAIL: %s\n", (msg)); \
            fails++; \
        } \
    } while (0)

/* A W x H grid with rows along +z and columns along +x at y = base + amp *
 * sin(2 pi r / H) cos(2 pi c / W); the node (10..12, 10..12) block is a hole. */
static void txr_test_grid(float *pos, size_t W, size_t H, double base, double amp)
{
    const size_t plane = W * H;
    size_t r = 0, c = 0;
    for (r = 0; r < H; r++) {
        for (c = 0; c < W; c++) {
            const size_t i = r * W + c;
            const double bump = amp * sin(2.0 * M_PI * (double)r / (double)H) * cos(2.0 * M_PI * (double)c / (double)W);
            if (r >= 10 && r < 13 && c >= 10 && c < 13) {
                pos[i] = pos[plane + i] = pos[2 * plane + i] = -1.0f;
                continue;
            }
            pos[i] = (float)(8 + r);                  /* z */
            pos[plane + i] = (float)(base + bump);    /* y: the across-sheet axis */
            pos[2 * plane + i] = (float)(8 + c);      /* x */
        }
    }
}

int TifxyzRecto_selftest(void)
{
    enum { CHUNK = 128, SIDE = 96 };    /* the bend's wavelength (SIDE) is 6 x the median window, as a real drift is */
    char root[256], dir[300];
    const size_t W = SIDE, H = SIDE, n = (size_t)SIDE * SIDE, plane = n;
    ves_selftest_dir(root, sizeof root, "tifxyz_recto");
    snprintf(dir, sizeof dir, "%s/raw", root);
    Arena_T arena = Arena_new();
    uint8_t *vol = (uint8_t *)malloc((size_t)CHUNK * CHUNK * CHUNK);
    float *pos = (float *)malloc(3 * n * sizeof(float));
    float *tf = (float *)malloc(n * sizeof(float));
    float *dp = (float *)malloc(n * sizeof(float));
    float *tl = (float *)malloc(n * sizeof(float));
    float corners[6] = { 8.0f, 10.0f, 8.0f, 8.0f + SIDE, 40.0f, 8.0f + SIDE };
    CubeTable ct;
    TifxyzRectoOpts o;
    TifxyzRectoPrior prior;
    TifxyzRectoStats st, st2, st3;
    char path[512];
    size_t i = 0;
    int fails = 0, z = 0, y = 0, x = 0;

    fprintf(stderr, "[selftest] tifxyz_recto\n");
    if (vol == NULL || pos == NULL || tf == NULL || dp == NULL || tl == NULL) {
        free(vol); free(pos); free(tf); free(dp); free(tl);
        Arena_dispose(&arena);
        return 1;
    }
    /* two flat bands (200) with one brighter centre voxel each (220): the own wrap at y = 19..25 (7 vox, the
     * cand5 median: its inward (+y) edge at 25.5, the verso at 18.5, the midline at 22) and the neighbouring OUTER
     * wrap at y = 9..15 (its recto at 15.5), a 3 vox gap between; the edges stay exactly at the voxel steps and
     * each band shows a crest at its centre (a flat plateau alone hides a band's crest: its profile only rises
     * into the search range's end) */
    for (z = 0; z < CHUNK; z++)
        for (y = 0; y < CHUNK; y++)
            for (x = 0; x < CHUNK; x++) {
                const int own = y >= 19 && y <= 25, outer = y >= 9 && y <= 15;
                vol[((size_t)z * CHUNK + (size_t)y) * CHUNK + (size_t)x] =
                    (uint8_t)((y == 22 || y == 12) ? 220 : (own || outer) ? 200 : 20);
            }
    snprintf(path, sizeof path, "%s/z00000_y00000_x00000.tif", dir);
    ves_ensure_parent_dir(path);
    TXR_CHECK(TiffIO_save(path, vol, CHUNK, CHUNK, CHUNK) == 0, "fixture write");
    free(vol);
    TXR_CHECK(cubetable_init(&ct, arena, dir, CHUNK, corners, 2, 30.0) == 0, "cubetable init");
    cubetable_prewarm_all(&ct);
    memset(&o, 0, sizeof o);

    /* 0. the prior reads the wraps from a grid inside the own band */
    txr_test_grid(pos, W, H, 22.0, 1.2);
    TXR_CHECK(TifxyzRecto_prior(arena, &ct, pos, W, H, 2, &prior) == 0, "prior");
    fprintf(stderr, "  prior: rays %zu, with band %zu, recto %zu, next wrap %zu, seam %zu; air %.0f band %.0f (separability %.2f), "
                    "thickness p50 %.2f, spacing p50 %.2f, gap depth p50 %.0f, slope p50 %.0f, edge width p50 %.2f; gate %.0f "
                    "(error %.3f), search %.2f, cap %.2f, sigma %.2f, thick %.1f..%.1f%s\n",
            prior.rays, prior.rays_with_band, prior.rays_with_recto, prior.rays_with_next, prior.rays_with_seam,
            prior.air_level, prior.band_level, prior.level_separability, prior.thickness.p50, prior.spacing.p50,
            prior.gap_depth.p50, prior.slope.p50, prior.edge_width.p50, prior.min_contrast, prior.gate_error,
            prior.search_vox, prior.max_total_vox, prior.sigma, prior.thick_min, prior.thick_max,
            prior.crowded ? " (crowded)" : "");
    TXR_CHECK(prior.air_level < 40.0 && prior.band_level > 190.0, "the air and band levels are read");
    TXR_CHECK(prior.rays_with_recto >= 100 && prior.rays_with_recto == prior.rays_with_band, "every ray finds its band's recto");
    TXR_CHECK(fabs(prior.thickness.p50 - 7.0) < 0.6, "the band is 7 vox thick");
    TXR_CHECK(fabs(prior.spacing.p50 - 10.0) < 0.6, "the wraps are 10 vox apart");
    TXR_CHECK(prior.gap_depth.p50 > 150.0 && prior.slope.p50 > 60.0, "the gap is deep and the edge steep");
    TXR_CHECK(prior.min_contrast > 10.0 && prior.min_contrast < 170.0, "the gate lies between the seam dips and the gaps");
    TXR_CHECK(prior.search_vox >= 4.0 && prior.search_vox <= 9.0 && !prior.crowded,
              "the search range covers the recto's reach and stays inside the wrap spacing");

    /* 1. the recto, on the prior's parameters: a grid bent from y 22.5 (inside the band) to 26.5 (in the gap
     * beyond its edge) */
    txr_test_grid(pos, W, H, 24.5, 2.0);
    TXR_CHECK(TifxyzRecto_run(arena, &ct, pos, W, H, &o, &prior, &st, tf, tl, NULL, NULL, dp) == 0, "run");
    {
        size_t bad = 0, nv = 0, jumps = 0, holes = 0;
        double sum = 0.0, maxj = 0.0;
        for (i = 0; i < n; i++) {
            const size_t r = i / W, c = i % W;
            if (pos[i] < -0.5f) { holes++; continue; }
            nv++;
            sum += fabs((double)pos[plane + i] - 25.5);
            if (fabs((double)pos[plane + i] - 25.5) > 0.2) bad++;
            if (c + 1 < W && pos[i + 1] > -0.5f) {
                const double j = fabs((double)pos[plane + i] - (double)pos[plane + i + 1]);
                if (j > maxj) maxj = j;
                if (j > 0.15) jumps++;
            }
            if (r + 1 < H && pos[i + W] > -0.5f) {
                const double j = fabs((double)pos[plane + i] - (double)pos[plane + i + W]);
                if (j > maxj) maxj = j;
                if (j > 0.15) jumps++;
            }
        }
        fprintf(stderr, "  recto: valid %zu, supported first %zu (edge inward %zu, outward %zu), iters %d, mean |y - 25.5| %.3f, "
                        "off by > 0.2: %zu, max neighbour jump %.3f, target first p50 %.2f p5 %.2f p95 %.2f, last p50 %.3f "
                        "mean|.| %.3f, disp p50 %.2f, profile(+8) %.1f profile(0) %.1f; ran on search %.2f gate %.0f sigma %.2f\n",
                nv, st.n_supported_first, st.n_in_band_first, st.n_in_gap_first, st.outer_iters, sum / (double)nv,
                bad, maxj, st.target_first.p50, st.target_first.p5, st.target_first.p95, st.target_last.p50,
                st.target_last.mean_abs, st.disp.p50, st.profile_first[TXR_PROFILE_HALF + 8],
                st.profile_first[TXR_PROFILE_HALF], st.prm_search, st.prm_min_contrast, st.prm_sigma);
        for (i = 0; i < n; i++)            /* the worst nodes, for the record */
            if (pos[i] > -0.5f && fabs((double)pos[plane + i] - 25.5) > 0.1)
                fprintf(stderr, "    node (%zu, %zu): y %.3f, target first %.3f, last %.3f, disp %.3f\n", i / W, i % W,
                        (double)pos[plane + i], (double)tf[i], (double)tl[i], (double)dp[i]);
        TXR_CHECK(holes == 9, "the hole stays invalid");
        TXR_CHECK(st.n_supported_first == nv, "every valid node finds the recto edge");
        TXR_CHECK(st.n_in_band_first > 0 && st.n_in_gap_first > 0, "nodes start both inside the band and in the gap");
        TXR_CHECK(sum / (double)nv < 0.08, "the grid lands on the inward edge of the band");
        TXR_CHECK(bad == 0, "no node stays more than 0.2 vox off the edge (a narrow residual peak at a border, "
                            "narrower than half the median window, is tolerated by design)");
        TXR_CHECK(maxj < 0.15 && jumps == 0, "the moved grid is flat");
        TXR_CHECK(st.target_last.mean_abs < 0.08, "the residual offset is small");
        TXR_CHECK(st.target_last_med.mean_abs < 0.05, "the median-filtered residual is smaller still");
        TXR_CHECK(st.profile_first[TXR_PROFILE_HALF + 8] < st.profile_first[TXR_PROFILE_HALF] - 50.0,
                  "the mean profile is dark 8 vox INWARD of the bright band");
        TXR_CHECK(st.profile_first[TXR_PROFILE_HALF - 7] < st.profile_first[TXR_PROFILE_HALF] - 30.0,
                  "... and darker 7 vox outward (the gap before the outer wrap)");
        TXR_CHECK(isfinite(tf[0]) && isfinite(dp[0]) && isfinite(tl[0]) && !isfinite(dp[10 * W + 10]),
                  "diagnostic rasters: defined on valid nodes, NaN in the hole");
        TXR_CHECK(fabs(st.prm_search - prior.search_vox) < 0.3 && fabs(st.prm_min_contrast - prior.min_contrast) < 0.01,
                  "the run used the prior's parameters");
    }
    /* 2. idempotent */
    TXR_CHECK(TifxyzRecto_run(arena, &ct, pos, W, H, &o, &prior, &st2, NULL, NULL, NULL, NULL, NULL) == 0, "second run");
    fprintf(stderr, "  again: disp mean |.| %.3f max %.3f\n", st2.disp.mean_abs, st2.disp.max_abs);
    TXR_CHECK(st2.disp.mean_abs < 0.05 && st2.disp.max_abs < 0.2, "a second pass barely moves");

    /* 2b. two wraps: a patch of nodes in the dark gap OUTSIDE the own band (y = 17.5: 1 vox outside its verso, the
     * outer wrap's recto 2 vox outward; its own band's crest 4.5 vox inward is nearer than the outer wrap's 5.5 vox
     * outward, and its own recto lies 8 vox away, beyond the range) must follow its neighbourhood inward, never
     * the outer wrap's edge; the rest of the grid (y = 22 +- 1.2, in the band) lands on the edge.  A node nearer to
     * the other wrap's body than to its own is, by the rule, that wrap's: the fixture does not go there.  The
     * patch's core (rows and columns 36..59, 24 wide: its middle is beyond a median half-window from the ramp) is
     * 4.5 vox outward, with a raised-cosine ramp over 6 nodes around it (a real grid has no cliffs: the as-built
     * close-ups' slope p99 is 0.3-0.9 vox per node).  On the constants, not the prior (the range is what matters). */
    txr_test_grid(pos, W, H, 22.0, 1.2);
    for (i = 0; i < n; i++) {
        const double r = (double)(i / W), c = (double)(i % W);
        const double dr = r < 36.0 ? 36.0 - r : r > 59.0 ? r - 59.0 : 0.0;
        const double dc = c < 36.0 ? 36.0 - c : c > 59.0 ? c - 59.0 : 0.0;
        const double dist = sqrt(dr * dr + dc * dc);
        const double w = dist >= 6.0 ? 0.0 : 0.5 * (1.0 + cos(M_PI * dist / 6.0));
        if (pos[i] > -0.5f) pos[plane + i] = (float)((double)pos[plane + i] - 4.5 * w);
    }
    {
        TifxyzRectoStats st7;
        size_t nv = 0, bad = 0, patch = 0, patch_in = 0, patch_out = 0, neg = 0, shown = 0;
        double sum = 0.0, patch_move = 0.0, patch_min = 1e9, patch_max = -1e9;
        TXR_CHECK(TifxyzRecto_run(arena, &ct, pos, W, H, &o, NULL, &st7, NULL, NULL, NULL, NULL, dp) == 0, "two-wrap run");
        for (i = 0; i < n; i++) {        /* every target in this fixture is inward: any outward move is a wrong assignment */
            if (pos[i] > -0.5f && dp[i] < -0.05f) {
                neg++;
                if (shown++ < 6)
                    fprintf(stderr, "    moved outward: node (%zu, %zu) y %.2f disp %.2f\n", i / W, i % W, (double)pos[plane + i],
                            (double)dp[i]);
            }
        }
        for (i = 0; i < n; i++) {
            const size_t r = i / W, c = i % W;
            const int in_patch = r >= 44 && r < 52 && c >= 44 && c < 52;        /* the core's middle, 8 x 8 */
            const int near_patch = r + 14 >= 36 && r < 74 && c + 14 >= 36 && c < 74;
            if (pos[i] < -0.5f) continue;
            if (in_patch) {
                const double mv = (double)pos[plane + i] - 17.5;
                patch++;
                patch_move += mv;
                if (mv < patch_min) patch_min = mv;
                if (mv > patch_max) patch_max = mv;
                if (mv >= 4.0) patch_in++;
                if (mv < 0.0) patch_out++;
                continue;
            }
            if (near_patch) continue;
            nv++;
            sum += fabs((double)pos[plane + i] - 25.5);
            if (fabs((double)pos[plane + i] - 25.5) > 0.2) bad++;
        }
        fprintf(stderr, "  two wraps: band beyond reach %zu; nodes moved outward %zu; patch middle %zu nodes moved %.2f vox "
                        "inward on average (min %.2f, max %.2f; >= 4.0: %zu); the far grid: mean |y - 25.5| %.3f, off > 0.2: %zu\n",
                st7.n_beyond_first, neg, patch, patch_move / (double)patch, patch_min, patch_max, patch_in, sum / (double)nv, bad);
        TXR_CHECK(st7.n_beyond_first > 0, "the patch's own band is in view with its recto out of reach");
        TXR_CHECK(neg == 0, "no node moves outward: nothing takes the outer wrap's edge");
        TXR_CHECK(patch_out == 0 && patch_in == patch, "the patch's middle follows its neighbourhood inward, to the cap");
        TXR_CHECK(sum / (double)nv < 0.08 && bad == 0, "the far grid lands on its own recto edge");
    }

    /* 3. the midline: the band centre at y = 22 */
    txr_test_grid(pos, W, H, 22.0, 1.2);
    o.midline = 1;
    TXR_CHECK(TifxyzRecto_run(arena, &ct, pos, W, H, &o, &prior, &st3, NULL, NULL, NULL, NULL, NULL) == 0, "midline run");
    {
        size_t nv = 0, bad = 0;
        double sum = 0.0;
        for (i = 0; i < n; i++) {
            if (pos[i] < -0.5f) continue;
            nv++;
            sum += fabs((double)pos[plane + i] - 22.0);
            if (fabs((double)pos[plane + i] - 22.0) > 0.2) bad++;
        }
        fprintf(stderr, "  midline: supported %zu of %zu, mean |y - 22| %.3f, off by > 0.2: %zu\n",
                st3.n_supported_first, nv, sum / (double)nv, bad);
        TXR_CHECK(st3.n_supported_first == nv, "every node finds both edges");
        TXR_CHECK(sum / (double)nv < 0.1 && bad == 0, "the grid lands on the band's midline");
    }
    /* 3b. the inset: the recto edge minus TXR_RECTO_INSET_VOX, inside the band */
    txr_test_grid(pos, W, H, 24.5, 2.0);
    o.midline = 0;
    o.inset = 1;
    {
        TifxyzRectoStats st5;
        size_t nv = 0, bad = 0;
        double sum = 0.0;
        const double want = 25.5 - TXR_RECTO_INSET_VOX;
        TXR_CHECK(TifxyzRecto_run(arena, &ct, pos, W, H, &o, &prior, &st5, NULL, NULL, NULL, NULL, NULL) == 0, "inset run");
        for (i = 0; i < n; i++) {
            if (pos[i] < -0.5f) continue;
            nv++;
            sum += fabs((double)pos[plane + i] - want);
            if (fabs((double)pos[plane + i] - want) > 0.2) bad++;
        }
        fprintf(stderr, "  inset %.2f: supported %zu of %zu, mean |y - %.2f| %.3f, off by > 0.2: %zu\n", TXR_RECTO_INSET_VOX,
                st5.n_supported_first, nv, want, sum / (double)nv, bad);
        TXR_CHECK(st5.n_supported_first == nv && sum / (double)nv < 0.1 && bad == 0,
                  "the grid stops TXR_RECTO_INSET_VOX outward of the edge");
    }
    /* 3c. the smooth window: the drift only; the fixture's bend is a drift, so the grid still heads for the edge */
    txr_test_grid(pos, W, H, 24.5, 2.0);
    o.inset = 0;
    o.smooth = 1;
    {
        TifxyzRectoStats st6;
        size_t nv = 0;
        double sum = 0.0;
        TXR_CHECK(TifxyzRecto_run(arena, &ct, pos, W, H, &o, &prior, &st6, NULL, NULL, NULL, NULL, NULL) == 0, "smooth run");
        for (i = 0; i < n; i++) {
            if (pos[i] < -0.5f) continue;
            nv++;
            sum += fabs((double)pos[plane + i] - 25.5);
        }
        fprintf(stderr, "  smooth (%d-node window, stride %d): filtered %zu of %zu, mean |y - 25.5| %.3f, iters %d\n",
                2 * TXR_MEDIAN_HALF_SMOOTH + 1, TXR_MEDIAN_STRIDE_SMOOTH, st6.n_filtered_first, nv, sum / (double)nv, st6.outer_iters);
        /* the fixture's bend (wavelength 96) is sharper than any real drift against a 63-node window: the median's
         * curvature bias makes the fixed point slow (0.21 vox left of ~1.0 after 6 iterations, 2026-10-10) */
        TXR_CHECK(st6.n_filtered_first == nv && sum / (double)nv < 0.3, "the wide window follows most of the bend to the edge");
    }
    o.smooth = 0;
    /* 5. windows: 32-column cores with a 32-column halo of the input give the whole run's grid in every core (the
     * halo carries the median window and the solve's reach), and the cores' counts add up to the whole's; one
     * window is the plain run, exactly; a column subset leaves every other column as it came */
    {
        float *pin = (float *)malloc(3 * n * sizeof(float)), *pw = (float *)malloc(3 * n * sizeof(float));
        float *pwin = (float *)malloc(3 * n * sizeof(float)), *pplain = (float *)malloc(3 * n * sizeof(float));
        TifxyzRectoStats sw1, sw2, sw3, sp;
        size_t nv = 0, outside_moved = 0;
        double maxd = 0.0, maxp = 0.0, maxs = 0.0;
        if (pin == NULL || pw == NULL || pwin == NULL || pplain == NULL) {
            fails++;
        } else {
            txr_test_grid(pin, W, H, 24.5, 2.0);
            memcpy(pplain, pin, 3 * n * sizeof(float));
            TXR_CHECK(TifxyzRecto_run(arena, &ct, pplain, W, H, &o, &prior, &sp, NULL, NULL, NULL, NULL, NULL) == 0, "plain run");
            TXR_CHECK(TifxyzRecto_run_windows(arena, &ct, pin, pw, W, H, 0, W, 0, 0, &o, &prior, &sw1, NULL, NULL) == 0
                      && sw1.windows == 1, "one window");
            TXR_CHECK(TifxyzRecto_run_windows(arena, &ct, pin, pwin, W, H, 0, W, 32, 32, &o, &prior, &sw2, NULL, NULL) == 0
                      && sw2.windows == 3, "three windows");
            for (i = 0; i < 3 * n; i++) {
                if (pin[i % n] < -0.5f) continue;
                if (i < n) nv++;
                if (fabs((double)pw[i] - (double)pplain[i]) > maxp) maxp = fabs((double)pw[i] - (double)pplain[i]);
                if (fabs((double)pwin[i] - (double)pw[i]) > maxd) maxd = fabs((double)pwin[i] - (double)pw[i]);
            }
            fprintf(stderr, "  windows: one window vs the plain run max |d| %.2e; 3 x 32 columns + 32 halo vs whole max |d| %.4f "
                            "vox; supported first %zu / %zu (whole %zu), valid %zu / %zu, capped %zu / %zu\n",
                    maxp, maxd, sw2.n_supported_first, sw1.n_supported_first, sp.n_supported_first, sw2.n_valid, sw1.n_valid,
                    sw2.n_capped_total, sw1.n_capped_total);
            TXR_CHECK(maxp < 1e-4, "one window is the plain run (up to the reductions' order)");
            TXR_CHECK(sw1.n_supported_first == sp.n_supported_first && sw1.n_valid == sp.n_valid
                      && fabs(sw1.target_first.p50 - sp.target_first.p50) < 1e-9 && fabs(sw1.slope_p99 - sp.slope_p99) < 1e-3,
                      "one window's statistics are the plain run's");
            TXR_CHECK(maxd < 0.05, "the windows' cores match the whole run");
            TXR_CHECK(sw2.n_supported_first == sw1.n_supported_first && sw2.n_valid == sw1.n_valid && sw2.n_valid == nv,
                      "the cores' counts add up to the whole's");
            /* the columns [40, 56) only: the rest untouched, the subset as the whole run left it (its halo reaches the
             * grid's edges, so the subset's window is the whole grid) */
            TXR_CHECK(TifxyzRecto_run_windows(arena, &ct, pin, pwin, W, H, 40, 56, 0, 64, &o, &prior, &sw3, NULL, NULL) == 0
                      && sw3.windows == 1, "a column subset");
            for (i = 0; i < 3 * n; i++) {
                const size_t c = (i % n) % W;
                if (pin[i % n] < -0.5f) continue;
                if (c < 40 || c >= 56) {
                    if (pwin[i] != pin[i]) outside_moved++;
                } else if (fabs((double)pwin[i] - (double)pw[i]) > maxs) {
                    maxs = fabs((double)pwin[i] - (double)pw[i]);
                }
            }
            fprintf(stderr, "  column subset [40, 56): %zu values moved outside it, inside vs whole max |d| %.2e, valid %zu\n",
                    outside_moved, maxs, sw3.n_valid);
            TXR_CHECK(outside_moved == 0 && maxs < 1e-4 && sw3.n_valid == 16 * H, "the subset moves only its columns");
        }
        free(pin); free(pw); free(pwin); free(pplain);
    }
    /* 6. the selective prewarm: a table over two chunks, points within reach of the first only -> the second is
     * marked absent unread, the first serves samples */
    {
        char dir2[300];
        snprintf(dir2, sizeof dir2, "%s/raw2", root);
        uint8_t *v2 = (uint8_t *)calloc((size_t)CHUNK * CHUNK * CHUNK, 1);
        float corners2[6] = { 8.0f, 8.0f, 8.0f, 100.0f, 100.0f, 200.0f };
        float pz[3] = { 50.0f, -1.0f, 60.0f }, py[3] = { 50.0f, 20.0f, 60.0f }, px[3] = { 50.0f, 300.0f, 90.0f };
        CubeTable ct2;
        int got = 0;
        if (v2 == NULL) {
            fails++;
        } else {
            for (i = 0; i < (size_t)CHUNK * CHUNK * CHUNK; i++) v2[i] = (uint8_t)(i % 251);
            snprintf(path, sizeof path, "%s/z00000_y00000_x00000.tif", dir2);
            ves_ensure_parent_dir(path);
            TXR_CHECK(TiffIO_save(path, v2, CHUNK, CHUNK, CHUNK) == 0, "prewarm fixture 0");
            snprintf(path, sizeof path, "%s/z00000_y00000_x00128.tif", dir2);
            TXR_CHECK(TiffIO_save(path, v2, CHUNK, CHUNK, CHUNK) == 0, "prewarm fixture 1");
            free(v2);
            TXR_CHECK(cubetable_init(&ct2, arena, dir2, CHUNK, corners2, 2, 2.0) == 0 && ct2.nx == 2, "two-chunk table");
            /* the second point is invalid (z -1) and would reach the second chunk */
            got = cubetable_prewarm_points(&ct2, pz, py, px, 3, 6.0);
            fprintf(stderr, "  selective prewarm: loaded %d, untouched %zu, missing %d; sample in the first %d, in the second %d\n",
                    got, ct2.n_untouched, ct2.n_missing, cube_fetch(&ct2, 10, 10, 10), cube_fetch(&ct2, 10, 10, 140));
            TXR_CHECK(got == 1 && ct2.n_loaded == 1 && ct2.n_untouched == 1 && ct2.n_missing == 0,
                      "only the reached chunk loads; the other is marked absent unread");
            TXR_CHECK(cube_fetch(&ct2, 10, 10, 10) == (int)((((size_t)10 * CHUNK + 10) * CHUNK + 10) % 251)
                      && cube_fetch(&ct2, 10, 10, 140) < 0 && ct2.n_loaded == 1,
                      "the loaded chunk serves samples; the untouched one reads no datum and stays unloaded");
        }
    }
    /* 4. empty input */
    {
        TifxyzRectoStats st4;
        TifxyzRectoPrior pr4;
        TXR_CHECK(TifxyzRecto_run(arena, &ct, pos, 0, 0, &o, NULL, &st4, NULL, NULL, NULL, NULL, NULL) == 0,
                  "empty input is a no-op");
        for (i = 0; i < 3 * n; i++) pos[i] = -1.0f;
        TXR_CHECK(TifxyzRecto_run(arena, &ct, pos, W, H, &o, NULL, &st4, NULL, NULL, NULL, NULL, NULL) == 0
                  && st4.n_valid == 0, "an all-invalid grid is a no-op");
        TXR_CHECK(TifxyzRecto_prior(arena, &ct, pos, W, H, 2, &pr4) == 0 && pr4.rays == 0
                  && fabs(pr4.search_vox - TXR_SEARCH_VOX) < 1e-9, "the prior of an all-invalid grid keeps the defaults");
    }
    free(pos); free(tf); free(dp); free(tl);
    Arena_dispose(&arena);
    if (fails == 0) ves_selftest_clean(root);
    fprintf(stderr, "[selftest] tifxyz_recto %s (%d failure%s)\n", fails == 0 ? "PASS" : "FAIL", fails,
            fails == 1 ? "" : "s");
    return fails;
}
