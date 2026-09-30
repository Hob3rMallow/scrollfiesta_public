/* ridge_track.c -- see ridge_track.h.
 *
 * Per cell we sample the CT along the cell's lattice normal at K discrete
 * offsets and build a data cost
 *      d(c, k) = (255 - CT(c, k)) / 255 + prior * |t_k|
 * then choose an offset field t(c) minimising
 *      sum_c d(c, t(c)) + smooth_u * |t(c) - t(c-1)|   (along a row)
 *                      + smooth_v * |t(c) - t(c-W)|    (across rows)
 * The u term is solved exactly per row by a dynamic program whose L1
 * transition is evaluated with a forward/backward distance transform
 * (Felzenszwalb & Huttenlocher), which keeps it O(K) per cell.  The v term is
 * handled by iterating: after the first round each cell's data cost gains a
 * pull toward the mean of its neighbouring rows' chosen offsets.
 *
 * Memory is per ROW, not per lattice: two W x K buffers, so a 3,700-column run
 * costs about 2 MB regardless of its height.
 */
#include "ridge_track.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/pipeline_constants.h"

#define RT_MAX_K 257

void RidgeTrack_defaults(RidgeTrackOpts *opts)
{
    opts->reach = QS_RIDGE_REACH;
    opts->step = QS_RIDGE_STEP;
    opts->prior = QS_RIDGE_PRIOR;
    opts->smooth_u = QS_RIDGE_SMOOTH_U;
    opts->smooth_v = QS_RIDGE_SMOOTH_V;
    opts->ridge_min = QS_RIDGE_MIN_U8;
    opts->iters = QS_RIDGE_ITERS;
    opts->max_step = QS_RIDGE_MAX_STEP;
    opts->edge_limit = QS_RIDGE_EDGE_LIMIT;
    opts->flatten = QS_RIDGE_FLATTEN;
    opts->flatten_win = QS_RIDGE_FLATTEN_WIN;
    opts->profile_win = QS_RIDGE_PROFILE_WIN;
    opts->fill_reach = QS_RIDGE_FILL_REACH;
    opts->fill_prior = QS_RIDGE_FILL_PRIOR;
}

/* median of a small array (copies, insertion sort: the window is <= 15) */
static double rt_median(double *a, int n)
{
    int i = 0, j = 0;
    if (n <= 0) return 0.0;
    for (i = 1; i < n; i++) {
        double v = a[i];
        for (j = i; j > 0 && a[j - 1] > v; j--) a[j] = a[j - 1];
        a[j] = v;
    }
    return a[n / 2];
}

/* Separable local median of `d` over valid cells: along u, then along v.  The
 * verdict uses a square window; separating it is close enough to steer the
 * offsets and costs 2*win instead of win^2 per cell. */
static void rt_local_median(const double *d, const uint8_t *has, int W, int H,
                            int win, double *tmp, double *out)
{
    int hw = win / 2, i = 0, j = 0, k = 0;
    double buf[17];
    for (j = 0; j < H; j++)
        for (i = 0; i < W; i++) {
            size_t c = (size_t)j * (size_t)W + (size_t)i;
            int n = 0;
            if (!has[c]) { tmp[c] = 0.0; continue; }
            for (k = -hw; k <= hw; k++) {
                int ii = i + k;
                size_t o = 0;
                if (ii < 0 || ii >= W) continue;
                o = (size_t)j * (size_t)W + (size_t)ii;
                if (!has[o]) continue;
                buf[n++] = d[o];
                if (n >= 17) break;
            }
            tmp[c] = n > 0 ? rt_median(buf, n) : d[c];
        }
    for (j = 0; j < H; j++)
        for (i = 0; i < W; i++) {
            size_t c = (size_t)j * (size_t)W + (size_t)i;
            int n = 0;
            if (!has[c]) { out[c] = 0.0; continue; }
            for (k = -hw; k <= hw; k++) {
                int jj = j + k;
                size_t o = 0;
                if (jj < 0 || jj >= H) continue;
                o = (size_t)jj * (size_t)W + (size_t)i;
                if (!has[o]) continue;
                buf[n++] = tmp[o];
                if (n >= 17) break;
            }
            out[c] = n > 0 ? rt_median(buf, n) : tmp[c];
        }
}

/* lattice normal at cell (i, j) from its neighbours; 0 when unavailable */
static int rt_normal(const float *pos, const uint8_t *valid, int W, int H,
                     int i, int j, double n[3])
{
    size_t c = (size_t)j * (size_t)W + (size_t)i;
    double du[3] = { 0.0, 0.0, 0.0 }, dv[3] = { 0.0, 0.0, 0.0 }, len = 0.0;
    int have_u = 0, have_v = 0;
    int k = 0;
    if (i + 1 < W && valid[c + 1] && i > 0 && valid[c - 1]) {
        for (k = 0; k < 3; k++)
            du[k] = (double)pos[(c + 1) * 3 + (size_t)k] - (double)pos[(c - 1) * 3 + (size_t)k];
        have_u = 1;
    } else if (i + 1 < W && valid[c + 1]) {
        for (k = 0; k < 3; k++)
            du[k] = (double)pos[(c + 1) * 3 + (size_t)k] - (double)pos[c * 3 + (size_t)k];
        have_u = 1;
    } else if (i > 0 && valid[c - 1]) {
        for (k = 0; k < 3; k++)
            du[k] = (double)pos[c * 3 + (size_t)k] - (double)pos[(c - 1) * 3 + (size_t)k];
        have_u = 1;
    }
    if (j + 1 < H && valid[c + (size_t)W] && j > 0 && valid[c - (size_t)W]) {
        for (k = 0; k < 3; k++)
            dv[k] = (double)pos[(c + (size_t)W) * 3 + (size_t)k] - (double)pos[(c - (size_t)W) * 3 + (size_t)k];
        have_v = 1;
    } else if (j + 1 < H && valid[c + (size_t)W]) {
        for (k = 0; k < 3; k++)
            dv[k] = (double)pos[(c + (size_t)W) * 3 + (size_t)k] - (double)pos[c * 3 + (size_t)k];
        have_v = 1;
    } else if (j > 0 && valid[c - (size_t)W]) {
        for (k = 0; k < 3; k++)
            dv[k] = (double)pos[c * 3 + (size_t)k] - (double)pos[(c - (size_t)W) * 3 + (size_t)k];
        have_v = 1;
    }
    if (!have_u || !have_v) return 0;
    n[0] = du[1] * dv[2] - du[2] * dv[1];
    n[1] = du[2] * dv[0] - du[0] * dv[2];
    n[2] = du[0] * dv[1] - du[1] * dv[0];
    len = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    if (!(len > 1e-9)) return 0;
    n[0] /= len; n[1] /= len; n[2] /= len;
    return 1;
}

/* in-place L1 distance transform: f[k] <- min_l (f[l] + lambda * |k - l|) */
static void rt_dt_l1(double *f, int K, double lambda)
{
    int k = 0;
    for (k = 1; k < K; k++)
        if (f[k] > f[k - 1] + lambda) f[k] = f[k - 1] + lambda;
    for (k = K - 2; k >= 0; k--)
        if (f[k] > f[k + 1] + lambda) f[k] = f[k + 1] + lambda;
}

/* BOUNDED transition: g[k] = min over |d| <= B of (f[k + d] + lambda * |d|).
 * The bound is what keeps an emitted edge from growing: two neighbouring cells
 * can never choose offsets more than B steps apart. */
static void rt_min_bounded(const double *f, double *g, int K, int B, double lambda)
{
    int k = 0, d = 0;
    for (k = 0; k < K; k++) {
        double best = f[k];
        for (d = 1; d <= B; d++) {
            if (k - d >= 0 && f[k - d] + lambda * (double)d < best) best = f[k - d] + lambda * (double)d;
            if (k + d < K && f[k + d] + lambda * (double)d < best) best = f[k + d] + lambda * (double)d;
        }
        g[k] = best;
    }
}

static int rt_cmp_double(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return x < y ? -1 : (x > y ? 1 : 0);
}

/* length of the edge (a, b), measured in the frame the verdict uses, after both
 * endpoints move along their normals */
static double rt_edge_len(const float *pos, const double *nrm, const double *cur,
                          size_t a, size_t b, RidgeTrackToWorld to_world, void *ctx)
{
    double pa[3], pb[3], wa[3], wb[3], d[3];
    int k = 0;
    for (k = 0; k < 3; k++) {
        pa[k] = (double)pos[a * 3 + (size_t)k] + cur[a] * nrm[a * 3 + (size_t)k];
        pb[k] = (double)pos[b * 3 + (size_t)k] + cur[b] * nrm[b * 3 + (size_t)k];
    }
    if (to_world != NULL) {
        to_world(ctx, pa, wa);
        to_world(ctx, pb, wb);
    } else {
        for (k = 0; k < 3; k++) { wa[k] = pa[k]; wb[k] = pb[k]; }
    }
    for (k = 0; k < 3; k++) d[k] = wb[k] - wa[k];
    return sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
}

/* the same edge before any displacement */
static double rt_edge_len0(const float *pos, size_t a, size_t b,
                           RidgeTrackToWorld to_world, void *ctx)
{
    double pa[3], pb[3], wa[3], wb[3], d[3];
    int k = 0;
    for (k = 0; k < 3; k++) {
        pa[k] = (double)pos[a * 3 + (size_t)k];
        pb[k] = (double)pos[b * 3 + (size_t)k];
    }
    if (to_world != NULL) {
        to_world(ctx, pa, wa);
        to_world(ctx, pb, wb);
    } else {
        for (k = 0; k < 3; k++) { wa[k] = pa[k]; wb[k] = pb[k]; }
    }
    for (k = 0; k < 3; k++) d[k] = wb[k] - wa[k];
    return sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
}

int RidgeTrack_solve(Arena_T arena, const RidgeTrackOpts *opts,
                     RidgeTrackSampler sample, RidgeTrackToWorld to_world, void *sample_ctx,
                     float *pos, const uint8_t *valid, const uint8_t *wide, int W, int H,
                     const int32_t *tri, size_t ntri, const double *depth,
                     RidgeTrackReport *rep)
{
    size_t HW = (size_t)W * (size_t)H;
    /* the offset axis spans the WIDEST class; a cell's own reach forbids the rest */
    double reach_max = (wide != NULL && opts->fill_reach > opts->reach) ? opts->fill_reach : opts->reach;
    int K = (int)(2.0 * reach_max / opts->step) + 1;
    double *data = NULL;      /* [W * K] this row's data cost */
    double *fwd = NULL;       /* [W * K] forward messages, for the backtrack */
    double *cur = NULL;       /* [HW] chosen offset */
    double *nrm = NULL;       /* [HW * 3] cell normals */
    double *scratch = NULL;   /* [HW] shifts / CT samples */
    double *dcur = NULL;      /* [HW] current depth (depth + offset) */
    double *dmed = NULL;      /* [HW] its local median */
    double *dtmp = NULL;      /* [HW] separable-median scratch */
    double *target = NULL;    /* [HW] offset that would level the cell */
    uint8_t *has = NULL;      /* [HW] cell has a normal */
    size_t nshift = 0;
    int it = 0, i = 0, j = 0, k = 0;
    Arena_Mark mark;

    if (arena == NULL || opts == NULL || sample == NULL || pos == NULL || valid == NULL || rep == NULL)
        return -1;
    memset(rep, 0, sizeof *rep);
    if (W < 3 || H < 3 || K < 3 || K > RT_MAX_K) return -1;
    mark = Arena_save(arena);
    data = (double *)ARENA_ALLOC(arena, (size_t)W * (size_t)K * sizeof *data);
    fwd = (double *)ARENA_ALLOC(arena, (size_t)W * (size_t)K * sizeof *fwd);
    cur = (double *)ARENA_CALLOC(arena, HW, sizeof *cur);
    nrm = (double *)ARENA_CALLOC(arena, HW * 3, sizeof *nrm);
    scratch = (double *)ARENA_ALLOC(arena, HW * sizeof *scratch);
    has = (uint8_t *)ARENA_CALLOC(arena, HW, 1);
    if (depth != NULL && opts->flatten > 0.0) {
        dcur = (double *)ARENA_ALLOC(arena, HW * sizeof *dcur);
        dmed = (double *)ARENA_ALLOC(arena, HW * sizeof *dmed);
        dtmp = (double *)ARENA_ALLOC(arena, HW * sizeof *dtmp);
        target = (double *)ARENA_CALLOC(arena, HW, sizeof *target);
        if (dcur == NULL || dmed == NULL || dtmp == NULL || target == NULL) {
            Arena_restore(arena, mark);
            return -1;
        }
    }
    if (data == NULL || fwd == NULL || cur == NULL || nrm == NULL || scratch == NULL || has == NULL) {
        Arena_restore(arena, mark);
        return -1;
    }
    /* normals and the "before" statistics */
    {
        size_t nb = 0;
        for (j = 0; j < H; j++)
            for (i = 0; i < W; i++) {
                size_t c = (size_t)j * (size_t)W + (size_t)i;
                double n[3];
                double v = 0.0;
                if (!valid[c]) continue;
                if (!rt_normal(pos, valid, W, H, i, j, n)) continue;
                has[c] = 1;
                nrm[c * 3] = n[0]; nrm[c * 3 + 1] = n[1]; nrm[c * 3 + 2] = n[2];
                v = sample(sample_ctx, (double)pos[c * 3], (double)pos[c * 3 + 1],
                           (double)pos[c * 3 + 2]);
                if (v >= opts->ridge_min) rep->on_ridge_before++;
                if (wide != NULL && wide[c]) { rep->fill_cells++; if (v >= opts->ridge_min) rep->fill_on_before++; }
                scratch[nb++] = v;
                rep->cells++;
            }
        if (nb > 0) {
            qsort(scratch, nb, sizeof *scratch, rt_cmp_double);
            rep->ct_before_p50 = scratch[nb / 2];
        }
    }
    if (rep->cells == 0) { Arena_restore(arena, mark); return 0; }

    for (it = 0; it < opts->iters; it++) {
        if (target != NULL) {
            /* the offset that would put each cell level with its neighbourhood */
            size_t c = 0, nh = 0;
            for (c = 0; c < HW; c++) dcur[c] = has[c] ? depth[c] + cur[c] : 0.0;
            rt_local_median(dcur, has, W, H, opts->flatten_win, dtmp, dmed);
            for (c = 0; c < HW; c++) {
                if (!has[c]) { target[c] = 0.0; continue; }
                target[c] = cur[c] - (dcur[c] - dmed[c]);
                if (it == 0) scratch[nh++] = fabs(dcur[c] - dmed[c]);
            }
            if (it == 0 && nh > 0) {
                qsort(scratch, nh, sizeof *scratch, rt_cmp_double);
                rep->hp_p90_before = scratch[(size_t)((double)nh * 0.9)];
            }
        }
        for (j = 0; j < H; j++) {
            int any = 0, last = -1;
            for (i = 0; i < W; i++) {
                size_t c = (size_t)j * (size_t)W + (size_t)i;
                if (!has[c]) continue;
                any = 1;
                last = i;
                for (k = 0; k < K; k++) {
                    double t = -reach_max + (double)k * opts->step;
                    double reach_c = (wide != NULL && wide[c]) ? opts->fill_reach : opts->reach;
                    double prior_c = (wide != NULL && wide[c]) ? opts->fill_prior : opts->prior;
                    double v = sample(sample_ctx,
                                      (double)pos[c * 3] + t * nrm[c * 3],
                                      (double)pos[c * 3 + 1] + t * nrm[c * 3 + 1],
                                      (double)pos[c * 3 + 2] + t * nrm[c * 3 + 2]);
                    double d = (255.0 - v) / 255.0 + prior_c * fabs(t);
                    if (fabs(t) > reach_c + 1e-9) d += 1e3;   /* outside this cell's own search */
                    if (target != NULL) d += opts->flatten * fabs(t - target[c]);
                    if (it > 0 && opts->smooth_v > 0.0) {
                        double acc = 0.0;
                        int nn = 0;
                        if (j > 0 && has[c - (size_t)W]) { acc += cur[c - (size_t)W]; nn++; }
                        if (j + 1 < H && has[c + (size_t)W]) { acc += cur[c + (size_t)W]; nn++; }
                        if (nn > 0) d += opts->smooth_v * fabs(t - acc / (double)nn);
                    }
                    data[(size_t)i * (size_t)K + (size_t)k] = d;
                }
            }
            if (!any) continue;
            /* A cell's CT profile is a per-cell maximum and flickers between the
             * two fibre layers of a wrap; averaging it over a few cells along u
             * turns the ridge into a surface fitted to a neighbourhood, which is
             * what the tracker actually wants to follow. */
            if (opts->profile_win > 1) {
                int hwp = opts->profile_win / 2;
                for (i = 0; i < W; i++) {
                    size_t c = (size_t)j * (size_t)W + (size_t)i;
                    if (!has[c]) continue;
                    for (k = 0; k < K; k++) {
                        double acc = 0.0;
                        int n = 0, q = 0;
                        for (q = -hwp; q <= hwp; q++) {
                            int ii = i + q;
                            size_t o = 0;
                            if (ii < 0 || ii >= W) continue;
                            o = (size_t)j * (size_t)W + (size_t)ii;
                            if (!has[o]) continue;
                            acc += data[(size_t)ii * (size_t)K + (size_t)k];
                            n++;
                        }
                        if (n > 0) fwd[(size_t)i * (size_t)K + (size_t)k] = acc / (double)n;
                    }
                }
                for (i = 0; i < W; i++) {
                    size_t c = (size_t)j * (size_t)W + (size_t)i;
                    if (!has[c]) continue;
                    for (k = 0; k < K; k++)
                        data[(size_t)i * (size_t)K + (size_t)k] = fwd[(size_t)i * (size_t)K + (size_t)k];
                }
            }
            /* forward messages along u; a gap in the row restarts the chain */
            {
                int prev = -1;
                for (i = 0; i < W; i++) {
                    size_t c = (size_t)j * (size_t)W + (size_t)i;
                    if (!has[c]) { prev = -1; continue; }
                    if (prev < 0) {
                        for (k = 0; k < K; k++)
                            fwd[(size_t)i * (size_t)K + (size_t)k] = data[(size_t)i * (size_t)K + (size_t)k];
                    } else {
                        double *m = fwd + (size_t)i * (size_t)K;
                        if (opts->max_step > 0)
                            rt_min_bounded(fwd + (size_t)prev * (size_t)K, m, K, opts->max_step,
                                           opts->smooth_u * opts->step);
                        else {
                            for (k = 0; k < K; k++) m[k] = fwd[(size_t)prev * (size_t)K + (size_t)k];
                            rt_dt_l1(m, K, opts->smooth_u * opts->step);
                        }
                        for (k = 0; k < K; k++) m[k] += data[(size_t)i * (size_t)K + (size_t)k];
                    }
                    prev = i;
                }
            }
            /* backtrack from the last valid cell of the row */
            {
                int best = 0, right = -1;
                for (i = last; i >= 0; i--) {
                    size_t c = (size_t)j * (size_t)W + (size_t)i;
                    if (!has[c]) { right = -1; continue; }
                    if (right < 0) {
                        best = 0;
                        for (k = 1; k < K; k++)
                            if (fwd[(size_t)i * (size_t)K + (size_t)k] <
                                fwd[(size_t)i * (size_t)K + (size_t)best]) best = k;
                    } else {
                        int bk = -1;
                        double bv = 0.0;
                        int klo = opts->max_step > 0 ? best - opts->max_step : 0;
                        int khi = opts->max_step > 0 ? best + opts->max_step : K - 1;
                        if (klo < 0) klo = 0;
                        if (khi > K - 1) khi = K - 1;
                        for (k = klo; k <= khi; k++) {
                            double v = fwd[(size_t)i * (size_t)K + (size_t)k] +
                                       opts->smooth_u * opts->step * fabs((double)k - (double)best);
                            if (bk < 0 || v < bv) { bv = v; bk = k; }
                        }
                        best = bk;
                    }
                    cur[c] = -reach_max + (double)best * opts->step;
                    right = i;
                }
            }
        }
    }
    /* LIPSCHITZ PROJECTION over the 4-neighbourhood.  The row DP bounds the
     * offset change along u exactly; across v it is only a soft pull, so two
     * vertically adjacent cells could still differ by the whole search width
     * and stretch an emitted edge.  Raster sweeps of
     *     t[c] <- min(t[c], t[n] + L),  t[c] <- max(t[c], t[n] - L)
     * converge to the nearest field obeying |t[c] - t[n]| <= L everywhere. */
    if (opts->max_step > 0) {
        double L = opts->smooth_u > 0.0 ? (double)opts->max_step * opts->step
                                        : (double)opts->max_step * opts->step;
        int sweep = 0;
        for (sweep = 0; sweep < 4; sweep++) {
            for (j = 0; j < H; j++)
                for (i = 0; i < W; i++) {
                    size_t c = (size_t)j * (size_t)W + (size_t)i;
                    if (!has[c]) continue;
                    if (i > 0 && has[c - 1]) {
                        if (cur[c] > cur[c - 1] + L) cur[c] = cur[c - 1] + L;
                        if (cur[c] < cur[c - 1] - L) cur[c] = cur[c - 1] - L;
                    }
                    if (j > 0 && has[c - (size_t)W]) {
                        if (cur[c] > cur[c - (size_t)W] + L) cur[c] = cur[c - (size_t)W] + L;
                        if (cur[c] < cur[c - (size_t)W] - L) cur[c] = cur[c - (size_t)W] - L;
                    }
                }
            for (j = H - 1; j >= 0; j--)
                for (i = W - 1; i >= 0; i--) {
                    size_t c = (size_t)j * (size_t)W + (size_t)i;
                    if (!has[c]) continue;
                    if (i + 1 < W && has[c + 1]) {
                        if (cur[c] > cur[c + 1] + L) cur[c] = cur[c + 1] + L;
                        if (cur[c] < cur[c + 1] - L) cur[c] = cur[c + 1] - L;
                    }
                    if (j + 1 < H && has[c + (size_t)W]) {
                        if (cur[c] > cur[c + (size_t)W] + L) cur[c] = cur[c + (size_t)W] + L;
                        if (cur[c] < cur[c + (size_t)W] - L) cur[c] = cur[c + (size_t)W] - L;
                    }
                }
        }
        for (j = 0; j < H; j++)
            for (i = 0; i < W; i++) {
                size_t c = (size_t)j * (size_t)W + (size_t)i;
                if (!has[c]) continue;
                if (i + 1 < W && has[c + 1]) {
                    double dd = fabs(cur[c] - cur[c + 1]);
                    if (dd > rep->neighbour_max) rep->neighbour_max = dd;
                }
                if (j + 1 < H && has[c + (size_t)W]) {
                    double dd = fabs(cur[c] - cur[c + (size_t)W]);
                    if (dd > rep->neighbour_max) rep->neighbour_max = dd;
                }
            }
    }
    /* EDGE CLAMP over the faces the caller will EMIT.  Two neighbouring cells
     * move along different normals, so even a bounded offset difference can
     * stretch the edge between them where the surface curves.  The verdict
     * tolerates no emitted edge over 6 voxels, so shrink the offsets on any
     * emitted edge that would exceed the limit until it fits.  Offsets only
     * shrink here, so a few sweeps converge.  An edge that was ALREADY over the
     * limit is left alone: it is the fit's, not ours. */
    if (opts->edge_limit > 0.0 && tri != NULL && ntri > 0) {
        int sweep = 0;
        for (sweep = 0; sweep < 8; sweep++) {
            size_t clamped = 0, t = 0;
            for (t = 0; t < ntri; t++) {
                int e = 0;
                for (e = 0; e < 3; e++) {
                    size_t a = (size_t)tri[t * 3 + (size_t)e];
                    size_t b = (size_t)tri[t * 3 + (size_t)((e + 1) % 3)];
                    double len0 = 0.0, len = 0.0;
                    if (a >= HW || b >= HW || (!has[a] && !has[b])) continue;
                    len0 = rt_edge_len0(pos, a, b, to_world, sample_ctx);
                    /* an edge already at or past the limit may not GROW; a shorter
                     * one may grow up to the limit.  Comparing against the limit
                     * alone let edges between the limit and the gate slip through
                     * and cross it (measured 2026-09-03: 159 such edges). */
                    len = rt_edge_len(pos, nrm, cur, a, b, to_world, sample_ctx);
                    /* 0.01 vox of margin: the clamp works in double on `pos`,
                     * the verdict measures the float32 file, and an edge sitting
                     * exactly on the 6.0 wrap gate rounds over it on the write. */
                    if (len0 > opts->edge_limit) {
                        /* an edge already at the fit's own wrap-gate maximum must not
                         * move at all: scaling leaves float-rounding residue that the
                         * verdict counts (measured: one edge at exactly 6.000) */
                        if (cur[a] != 0.0 || cur[b] != 0.0) {
                            cur[a] = 0.0;
                            cur[b] = 0.0;
                            clamped++;
                        }
                    } else if (len > opts->edge_limit) {
                        cur[a] *= 0.6;
                        cur[b] *= 0.6;
                        clamped++;
                    }
                }
            }
            rep->edge_clamped += clamped;
            if (clamped == 0) break;
        }
        {
            size_t t = 0;
            for (t = 0; t < ntri; t++) {
                int e = 0;
                for (e = 0; e < 3; e++) {
                    size_t a = (size_t)tri[t * 3 + (size_t)e];
                    size_t b = (size_t)tri[t * 3 + (size_t)((e + 1) % 3)];
                    double len = 0.0;
                    if (a >= HW || b >= HW || (!has[a] && !has[b])) continue;
                    len = rt_edge_len(pos, nrm, cur, a, b, to_world, sample_ctx);
                    if (len > rep->edge_max) rep->edge_max = len;
                }
            }
        }
    }
    if (target != NULL) {
        size_t c = 0, nh = 0;
        for (c = 0; c < HW; c++) dcur[c] = has[c] ? depth[c] + cur[c] : 0.0;
        rt_local_median(dcur, has, W, H, opts->flatten_win, dtmp, dmed);
        for (c = 0; c < HW; c++)
            if (has[c]) scratch[nh++] = fabs(dcur[c] - dmed[c]);
        if (nh > 0) {
            qsort(scratch, nh, sizeof *scratch, rt_cmp_double);
            rep->hp_p90_after = scratch[(size_t)((double)nh * 0.9)];
        }
    }
    /* apply */
    {
        size_t na = 0;
        double *shift = (double *)ARENA_ALLOC(arena, HW * sizeof *shift);
        if (shift == NULL) { Arena_restore(arena, mark); return -1; }
        for (size_t c = 0; c < HW; c++) {
            double t = cur[c], v = 0.0;
            if (!has[c]) continue;
            if (fabs(t) >= opts->step) rep->moved++;
            pos[c * 3] = (float)((double)pos[c * 3] + t * nrm[c * 3]);
            pos[c * 3 + 1] = (float)((double)pos[c * 3 + 1] + t * nrm[c * 3 + 1]);
            pos[c * 3 + 2] = (float)((double)pos[c * 3 + 2] + t * nrm[c * 3 + 2]);
            shift[nshift++] = fabs(t);
            v = sample(sample_ctx, (double)pos[c * 3], (double)pos[c * 3 + 1],
                       (double)pos[c * 3 + 2]);
            if (v >= opts->ridge_min) rep->on_ridge++;
            if (wide != NULL && wide[c] && v >= opts->ridge_min) rep->fill_on_after++;
            scratch[na++] = v;
        }
        if (na > 0) {
            qsort(scratch, na, sizeof *scratch, rt_cmp_double);
            rep->ct_after_p50 = scratch[na / 2];
        }
        if (nshift > 0) {
            qsort(shift, nshift, sizeof *shift, rt_cmp_double);
            rep->shift_p50 = shift[nshift / 2];
            rep->shift_p90 = shift[(size_t)((double)nshift * 0.9)];
            rep->shift_max = shift[nshift - 1];
        }
    }
    Arena_restore(arena, mark);
    return 0;
}

int RidgeTrack_selftest(void)
{
    int failures = 0;
    fprintf(stderr, "[selftest] ridge_track: L1 distance transform\n");
    {
        double f[5] = { 0.0, 10.0, 10.0, 10.0, 10.0 };
        int k = 0;
        rt_dt_l1(f, 5, 1.0);
        for (k = 0; k < 5; k++)
            if (fabs(f[k] - (double)k) > 1e-9) {
                fprintf(stderr, "  dt[%d] = %.3f, want %.3f\n", k, f[k], (double)k);
                failures++;
            }
    }
    {
        double f[4] = { 2.0, 2.0, 2.0, 2.0 };
        int k = 0;
        rt_dt_l1(f, 4, 1.0);
        for (k = 0; k < 4; k++)
            if (fabs(f[k] - 2.0) > 1e-9) { fprintf(stderr, "  flat dt broken\n"); failures++; break; }
    }
    fprintf(stderr, "[selftest] ridge_track: bounded transition\n");
    {
        double f[7] = { 9.0, 9.0, 9.0, 0.0, 9.0, 9.0, 9.0 };
        double g[7];
        int k = 0;
        rt_min_bounded(f, g, 7, 1, 1.0);
        /* only the cells within one step of the minimum may take it */
        if (fabs(g[3] - 0.0) > 1e-9 || fabs(g[2] - 1.0) > 1e-9 || fabs(g[4] - 1.0) > 1e-9 ||
            fabs(g[0] - 9.0) > 1e-9) {
            fprintf(stderr, "  bounded min wrong: ");
            for (k = 0; k < 7; k++) fprintf(stderr, "%.1f ", g[k]);
            fprintf(stderr, "\n");
            failures++;
        }
    }
    fprintf(stderr, "[selftest] ridge_track: normal of a planar patch\n");
    {
        /* a 3x3 patch in the plane z = 5, spaced 2 apart: normal must be +/-Z */
        float pos[27 * 3];
        uint8_t valid[9];
        double n[3];
        int i = 0, j = 0;
        for (j = 0; j < 3; j++)
            for (i = 0; i < 3; i++) {
                size_t c = (size_t)j * 3u + (size_t)i;
                pos[c * 3] = 5.0f;
                pos[c * 3 + 1] = (float)(2 * j);
                pos[c * 3 + 2] = (float)(2 * i);
                valid[c] = 1;
            }
        if (!rt_normal(pos, valid, 3, 3, 1, 1, n)) {
            fprintf(stderr, "  no normal on a planar patch\n");
            failures++;
        } else if (fabs(fabs(n[0]) - 1.0) > 1e-9) {
            fprintf(stderr, "  planar normal = (%.3f, %.3f, %.3f), want +/-Z\n", n[0], n[1], n[2]);
            failures++;
        }
    }
    fprintf(stderr, "[selftest] ridge_track %s (%d failure(s))\n",
            failures == 0 ? "PASS" : "FAIL", failures);
    return failures;
}
