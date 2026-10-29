/* rawtex_layers.c -- see rawtex_layers.h. */
#include "rawtex_layers.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/arena.h"
#include "../common/mesh_normals.h"
#include "../common/tiff_io.h"
#include "../common/ves_platform.h"

#include "../common/ves_omp.h"

/* The bake's total-pixel cap for 1 px/vox sheets (obj_bake_raw --raster-max-px
 * in the review stage); only the 2^20 per-axis cap can refuse a sheet here. */
#define RL_MAX_PX ((size_t)6000000000ULL)

/* ------------------------------------------------------------ components */

static int32_t rl_find(int32_t *parent, int32_t v)
{
    while (parent[v] != v) {
        parent[v] = parent[parent[v]];
        v = parent[v];
    }
    return v;
}

static void rl_union(int32_t *parent, int32_t a, int32_t b)
{
    int32_t ra = rl_find(parent, a), rb = rl_find(parent, b);
    if (ra == rb) return;
    if (ra < rb) parent[rb] = ra;
    else parent[ra] = rb;
}

static void rl_normalize(float *n)
{
    double l = sqrt((double)n[0] * n[0] + (double)n[1] * n[1] + (double)n[2] * n[2]);
    if (l > 1e-12) {
        n[0] = (float)((double)n[0] / l);
        n[1] = (float)((double)n[1] / l);
        n[2] = (float)((double)n[2] / l);
    }
}

int RawtexLayers_orient(float *normals, const float *verts, size_t nv,
                        const int32_t *faces, size_t nf, const float *inward,
                        int smooth_iters, RawtexOrientStats *st)
{
    int32_t *parent = NULL;
    uint8_t *used = NULL;
    double *sum = NULL, *mag = NULL, total = 0.0, agree = 0.0, flipped_mag = 0.0;
    size_t f = 0, v = 0;

    if (st != NULL) memset(st, 0, sizeof *st);
    if (nv == 0 || nf == 0) return 0;
    if (normals == NULL || verts == NULL || faces == NULL || inward == NULL ||
        nv > (size_t)INT32_MAX)
        return -1;
    parent = (int32_t *)malloc(nv * sizeof *parent);
    used = (uint8_t *)calloc(nv, 1);
    sum = (double *)calloc(nv, sizeof *sum);
    mag = (double *)calloc(nv, sizeof *mag);
    if (parent == NULL || used == NULL || sum == NULL || mag == NULL) {
        free(parent); free(used); free(sum); free(mag);
        return -1;
    }
    for (v = 0; v < nv; v++) parent[v] = (int32_t)v;
    for (f = 0; f < nf; f++) {
        int32_t a = faces[f * 3], b = faces[f * 3 + 1], c = faces[f * 3 + 2];
        if (a < 0 || b < 0 || c < 0 || (size_t)a >= nv || (size_t)b >= nv || (size_t)c >= nv) {
            free(parent); free(used); free(sum); free(mag);
            return -1;
        }
        used[a] = used[b] = used[c] = 1;
        rl_union(parent, a, b);
        rl_union(parent, a, c);
    }
    /* per-face vote: area x cosine between the winding normal and the mean
     * inward direction of its corners */
    for (f = 0; f < nf; f++) {
        size_t a = (size_t)faces[f * 3], b = (size_t)faces[f * 3 + 1], c = (size_t)faces[f * 3 + 2];
        double e1[3], e2[3], cr[3], in[3], il = 0.0, vote = 0.0;
        int k = 0;
        for (k = 0; k < 3; k++) {
            e1[k] = (double)verts[b * 3 + (size_t)k] - (double)verts[a * 3 + (size_t)k];
            e2[k] = (double)verts[c * 3 + (size_t)k] - (double)verts[a * 3 + (size_t)k];
            in[k] = (double)inward[a * 3 + (size_t)k] + (double)inward[b * 3 + (size_t)k]
                  + (double)inward[c * 3 + (size_t)k];
            il += in[k] * in[k];
        }
        cr[0] = e1[1] * e2[2] - e1[2] * e2[1];
        cr[1] = e1[2] * e2[0] - e1[0] * e2[2];
        cr[2] = e1[0] * e2[1] - e1[1] * e2[0];
        il = sqrt(il);
        if (il < 1e-12) continue;
        vote = 0.5 * (cr[0] * in[0] + cr[1] * in[1] + cr[2] * in[2]) / il;
        v = (size_t)rl_find(parent, (int32_t)a);
        sum[v] += vote;
        mag[v] += fabs(vote);
    }
    for (v = 0; v < nv; v++) {
        if (!used[v] || parent[v] != (int32_t)v) continue;
        if (st != NULL) {
            st->components++;
            if (sum[v] < 0.0) st->flipped++;
            if (fabs(sum[v]) < 0.25 * mag[v]) st->weak++;
        }
        total += mag[v];
        agree += 0.5 * (mag[v] + fabs(sum[v]));   /* |vote| on the majority side */
        if (sum[v] < 0.0) flipped_mag += mag[v];
    }
    for (v = 0; v < nv; v++) {
        if (!used[v]) continue;
        if (sum[rl_find(parent, (int32_t)v)] < 0.0) {
            normals[v * 3] = -normals[v * 3];
            normals[v * 3 + 1] = -normals[v * 3 + 1];
            normals[v * 3 + 2] = -normals[v * 3 + 2];
        }
    }
    if (st != NULL) {
        st->agree_area = total > 0.0 ? agree / total : 1.0;
        st->flipped_area = total > 0.0 ? flipped_mag / total : 0.0;
    }
    /* one-ring smoothing of the oriented normals; components share no vertex,
     * so a smoothing step never mixes two orientation decisions */
    if (smooth_iters > 0) {
        float *acc = (float *)malloc(nv * 3 * sizeof *acc);
        int it = 0;
        if (acc == NULL) { free(parent); free(used); free(sum); free(mag); return -1; }
        for (it = 0; it < smooth_iters; it++) {
            memset(acc, 0, nv * 3 * sizeof *acc);
            for (f = 0; f < nf; f++) {
                size_t a = (size_t)faces[f * 3], b = (size_t)faces[f * 3 + 1], c = (size_t)faces[f * 3 + 2];
                int k = 0;
                for (k = 0; k < 3; k++) {
                    float s = normals[a * 3 + (size_t)k] + normals[b * 3 + (size_t)k]
                            + normals[c * 3 + (size_t)k];
                    acc[a * 3 + (size_t)k] += s;
                    acc[b * 3 + (size_t)k] += s;
                    acc[c * 3 + (size_t)k] += s;
                }
            }
            for (v = 0; v < nv; v++) {
                if (!used[v]) continue;
                rl_normalize(&acc[v * 3]);
                memcpy(&normals[v * 3], &acc[v * 3], 3 * sizeof *acc);
            }
        }
        free(acc);
    }
    free(parent); free(used); free(sum); free(mag);
    return 0;
}

/* ---------------------------------------------------------------- render */

/* Trilinear sample; *full = 1 only when all eight corners hold RAW. */
static double rl_sample(CubeTable *ct, double z, double y, double x, int *full)
{
    double fz = floor(z), fy = floor(y), fx = floor(x);
    long iz = (long)fz, iy = (long)fy, ix = (long)fx;
    double dz = z - fz, dy = y - fy, dx = x - fx;
    double acc = 0.0, wsum = 0.0;
    int k = 0, all = 1;
    for (k = 0; k < 8; k++) {
        long oz = k & 1, oy = (k >> 1) & 1, ox = (k >> 2) & 1;
        int val = cube_fetch(ct, iz + oz, iy + oy, ix + ox);
        double w = 0.0;
        if (val < 0) { all = 0; continue; }
        w = (oz ? dz : 1.0 - dz) * (oy ? dy : 1.0 - dy) * (ox ? dx : 1.0 - dx);
        acc += w * (double)val;
        wsum += w;
    }
    *full = all;
    if (wsum < 1e-9) return -1.0;
    return acc / wsum;
}

typedef struct RlFace {
    double ua, va, ub, vb, uc, vc, A2;  /* pixel-space UV and doubled signed area */
    long x0, x1, y0, y1;                /* pixel-centre bbox, clamped to the raster */
} RlFace;

/* The bake's face prologue: pixel-space corners, UV-stretch gate, centre bbox.
 * Returns 1 when the face can paint, 0 when gated/degenerate. */
static int rl_face(const float *verts, RawtexUv uv, const int32_t *faces, size_t f,
                   double umin, double vmin, size_t W, size_t H,
                   double stretch_ratio, double stretch_floor, RlFace *o, int *bad_uv)
{
    size_t vi[3];
    int e = 0;
    double lox, hix, loy, hiy;
    vi[0] = (size_t)faces[f * 3]; vi[1] = (size_t)faces[f * 3 + 1]; vi[2] = (size_t)faces[f * 3 + 2];
    *bad_uv = 0;
    o->ua = Rawtex_uv(uv, vi[0] * 2) - umin; o->va = Rawtex_uv(uv, vi[0] * 2 + 1) - vmin;
    o->ub = Rawtex_uv(uv, vi[1] * 2) - umin; o->vb = Rawtex_uv(uv, vi[1] * 2 + 1) - vmin;
    o->uc = Rawtex_uv(uv, vi[2] * 2) - umin; o->vc = Rawtex_uv(uv, vi[2] * 2 + 1) - vmin;
    o->A2 = (o->ub - o->ua) * (o->vc - o->va) - (o->vb - o->va) * (o->uc - o->ua);
    if (fabs(o->A2) < 1e-12) return 0;
    for (e = 0; e < 3; e++) {
        size_t p = vi[e], q = vi[(e + 1) % 3];
        double d3 = 0.0, du = 0.0, dv = 0.0, e2d = 0.0;
        int k = 0;
        for (k = 0; k < 3; k++) {
            double dd = (double)verts[q * 3 + (size_t)k] - (double)verts[p * 3 + (size_t)k];
            d3 += dd * dd;
        }
        d3 = sqrt(d3);
        du = Rawtex_uv(uv, q * 2) - Rawtex_uv(uv, p * 2);
        dv = Rawtex_uv(uv, q * 2 + 1) - Rawtex_uv(uv, p * 2 + 1);
        e2d = sqrt(du * du + dv * dv);
        if (stretch_ratio > 0.0 &&
            e2d > (stretch_ratio * d3 > stretch_floor ? stretch_ratio * d3 : stretch_floor))
            *bad_uv = 1;
    }
    if (*bad_uv) return 0;
    lox = fmin(o->ua, fmin(o->ub, o->uc)); hix = fmax(o->ua, fmax(o->ub, o->uc));
    loy = fmin(o->va, fmin(o->vb, o->vc)); hiy = fmax(o->va, fmax(o->vb, o->vc));
    /* centres at index + 0.5: the first covered centre is ceil(lo - 0.5)
     * (rawtex_bake.c uses the same -0.501 / +0.001 bounds) */
    o->x0 = (long)ceil(lox - 0.501); o->x1 = (long)floor(hix + 0.001);
    o->y0 = (long)ceil(loy - 0.501); o->y1 = (long)floor(hiy + 0.001);
    if (o->x0 < 0) o->x0 = 0;
    if (o->y0 < 0) o->y0 = 0;
    if (o->x1 >= (long)W) o->x1 = (long)W - 1;
    if (o->y1 >= (long)H) o->y1 = (long)H - 1;
    return o->x0 <= o->x1 && o->y0 <= o->y1;
}

/* One seam-fill pass over a w x h raster (vmesh_tifxyz v2's rule, on positions
 * AND normals).  state: 0 empty, 1/2 painted (boundary / interior), 3 filled.
 * snap is scratch [w*h].  Returns the number of pixels filled. */
static size_t rl_seam_fill(uint8_t *state, float *pos, float *nrm, size_t w, size_t h,
                           int reach, double tol, uint8_t *snap)
{
    size_t n = w * h, filled = 0, r = 0, c = 0;
    memcpy(snap, state, n);
    for (r = 0; r < h; r++)
        for (c = 0; c < w; c++) {
            size_t i = r * w + c, ba = 0, bb = 0;
            int best_span = reach + 2, da_best = 0, db_best = 0, dir = 0, q = 0;
            if (snap[i]) continue;
            for (dir = 0; dir < 2; dir++) {
                long dc = dir == 0 ? 1 : 0, dr = dir == 0 ? 0 : 1;
                int da = 0, db = 0, k = 0;
                size_t ia = 0, ib = 0;
                double ex = 0.0, ey = 0.0, ez = 0.0;
                for (k = 1; k <= reach && !da; k++) {
                    long rr = (long)r - k * dr, cc = (long)c - k * dc;
                    size_t j = 0;
                    if (rr < 0 || cc < 0) break;
                    j = (size_t)rr * w + (size_t)cc;
                    if (snap[j]) { da = k; ia = j; }
                }
                for (k = 1; k <= reach && !db; k++) {
                    long rr = (long)r + k * dr, cc = (long)c + k * dc;
                    size_t j = 0;
                    if (rr >= (long)h || cc >= (long)w) break;
                    j = (size_t)rr * w + (size_t)cc;
                    if (snap[j]) { db = k; ib = j; }
                }
                if (!da || !db || da + db - 1 > reach || da + db >= best_span) continue;
                ex = (double)pos[ia * 3] - (double)pos[ib * 3];
                ey = (double)pos[ia * 3 + 1] - (double)pos[ib * 3 + 1];
                ez = (double)pos[ia * 3 + 2] - (double)pos[ib * 3 + 2];
                if (ex * ex + ey * ey + ez * ez > tol * tol) continue;
                best_span = da + db;
                ba = ia; bb = ib; da_best = da; db_best = db;
            }
            if (best_span > reach + 1) continue;
            {
                double wa = (double)db_best / (double)(da_best + db_best), wb = 1.0 - wa;
                for (q = 0; q < 3; q++) {
                    pos[i * 3 + (size_t)q] = (float)(wa * (double)pos[ba * 3 + (size_t)q] + wb * (double)pos[bb * 3 + (size_t)q]);
                    nrm[i * 3 + (size_t)q] = (float)(wa * (double)nrm[ba * 3 + (size_t)q] + wb * (double)nrm[bb * 3 + (size_t)q]);
                }
            }
            state[i] = 3;
            filled++;
        }
    return filled;
}

/* ------------------------------------------------------------ apron */

static inline void rl_try(int32_t *near_px, float *d2, size_t i, size_t j, long x, long y, long w)
{
    if (near_px[j] >= 0) {
        long fy = near_px[j] / w, fx = near_px[j] % w;
        float dd = (float)((x - fx) * (x - fx) + (y - fy) * (y - fy));
        if (dd < d2[i]) { d2[i] = dd; near_px[i] = near_px[j]; }
    }
}

/* Nearest painted pixel (state 1..3) of every pixel of a w x h raster: two raster passes propagating the
 * nearest feature through the 8-neighbourhood (8SSEDT; exact up to rare ties).  near_px = -1 where none. */
static void rl_nearest(const uint8_t *state, size_t w, size_t h, int32_t *near_px, float *d2)
{
    long x = 0, y = 0;
    const long W = (long)w, H = (long)h;
    size_t i = 0;
    for (i = 0; i < w * h; i++) {
        int painted = state[i] >= 1 && state[i] <= 3;
        near_px[i] = painted ? (int32_t)i : -1;
        d2[i] = painted ? 0.0f : 3.0e38f;
    }
    for (y = 0; y < H; y++) {
        for (x = 0; x < W; x++) {
            i = (size_t)(y * W + x);
            if (x > 0) rl_try(near_px, d2, i, i - 1, x, y, W);
            if (y > 0) {
                rl_try(near_px, d2, i, i - (size_t)W, x, y, W);
                if (x > 0) rl_try(near_px, d2, i, i - (size_t)W - 1, x, y, W);
                if (x < W - 1) rl_try(near_px, d2, i, i - (size_t)W + 1, x, y, W);
            }
        }
        for (x = W - 2; x >= 0; x--) {
            i = (size_t)(y * W + x);
            rl_try(near_px, d2, i, i + 1, x, y, W);
        }
    }
    for (y = H - 1; y >= 0; y--) {
        for (x = W - 1; x >= 0; x--) {
            i = (size_t)(y * W + x);
            if (x < W - 1) rl_try(near_px, d2, i, i + 1, x, y, W);
            if (y < H - 1) {
                rl_try(near_px, d2, i, i + (size_t)W, x, y, W);
                if (x < W - 1) rl_try(near_px, d2, i, i + (size_t)W + 1, x, y, W);
                if (x > 0) rl_try(near_px, d2, i, i + (size_t)W - 1, x, y, W);
            }
        }
        for (x = 1; x < W; x++) {
            i = (size_t)(y * W + x);
            rl_try(near_px, d2, i, i - 1, x, y, W);
        }
    }
}

/* 6x6 symmetric positive definite solve for three right-hand sides (Cholesky).  0 on success. */
static int rl_chol6(double A[36], double B[18])
{
    double L[36];
    int i = 0, j = 0, k = 0, r = 0;
    memset(L, 0, sizeof L);
    for (i = 0; i < 6; i++)
        for (j = 0; j <= i; j++) {
            double s = A[i * 6 + j];
            for (k = 0; k < j; k++) s -= L[i * 6 + k] * L[j * 6 + k];
            if (i == j) {
                if (!(s > 1e-12)) return -1;
                L[i * 6 + i] = sqrt(s);
            } else {
                L[i * 6 + j] = s / L[j * 6 + j];
            }
        }
    for (r = 0; r < 3; r++) {
        double yv[6], xv[6];
        for (i = 0; i < 6; i++) {
            double s = B[r * 6 + i];
            for (k = 0; k < i; k++) s -= L[i * 6 + k] * yv[k];
            yv[i] = s / L[i * 6 + i];
        }
        for (i = 5; i >= 0; i--) {
            double s = yv[i];
            for (k = i + 1; k < 6; k++) s -= L[k * 6 + i] * xv[k];
            xv[i] = s / L[i * 6 + i];
        }
        for (i = 0; i < 6; i++) B[r * 6 + i] = xv[i];
    }
    return 0;
}

int RawtexLayers_quadric_fit(const uint8_t *state, const float *pos, long W, long H, long cx, long cy,
                             int fit_px, int step_px, double coef[18])
{
    double A[36], B[18], wsum = 0.0;
    const double inv = 1.0 / (double)fit_px, s2 = 2.0 * (0.5 * fit_px) * (0.5 * fit_px);
    long x = 0, y = 0;
    size_t nsamp = 0;
    int j = 0, k = 0, q = 0;
    if (fit_px < 2 || step_px < 1) return -1;
    memset(A, 0, sizeof A);
    memset(B, 0, sizeof B);
    for (y = cy - fit_px; y <= cy + fit_px; y += step_px) {
        if (y < 0 || y >= H) continue;
        for (x = cx - fit_px; x <= cx + fit_px; x += step_px) {
            long dx = x - cx, dy = y - cy;
            size_t i = 0;
            double r2 = (double)(dx * dx + dy * dy), wgt = 0.0, phi[6];
            if (x < 0 || x >= W || r2 > (double)fit_px * fit_px) continue;
            i = (size_t)(y * W + x);
            if (state[i] < 1 || state[i] > 3) continue;
            wgt = exp(-r2 / s2);
            phi[0] = 1.0; phi[1] = dx * inv; phi[2] = dy * inv;
            phi[3] = phi[1] * phi[1]; phi[4] = phi[1] * phi[2]; phi[5] = phi[2] * phi[2];
            for (j = 0; j < 6; j++) {
                for (k = 0; k <= j; k++) A[j * 6 + k] += wgt * phi[j] * phi[k];
                for (q = 0; q < 3; q++) B[q * 6 + j] += wgt * phi[j] * (double)pos[i * 3 + (size_t)q];
            }
            wsum += wgt;
            nsamp++;
        }
    }
    if (nsamp < 12) return -1;
    for (j = 0; j < 6; j++)
        for (k = j + 1; k < 6; k++) A[j * 6 + k] = A[k * 6 + j];
    /* a RELATIVE ridge on the curvature terms keeps one-sided (border) fits of
     * thin strips well posed without biasing a sheet's curvature: an absolute
     * 1e-3 * wsum ridge doubled the 17 px error on a radius-100 turn, because
     * the Gaussian weights keep sum(w u^4) far below wsum */
    for (j = 3; j < 6; j++) A[j * 6 + j] += 1e-3 * A[j * 6 + j] + 1e-9 * wsum;
    if (rl_chol6(A, B) != 0) return -1;
    memcpy(coef, B, sizeof B);
    return 0;
}

/* Fill the context apron of one w x h tile raster (see RawtexLayersOpts): pixels within apron_px of the
 * painted sheet in columns [c0, c1) get state 4 and extrapolated position / normal.  Scratch is allocated
 * here.  Returns the number of apron pixels; *rejected / *fits count failed extrapolations and fits. */
/* Fit cells are GLOBAL (column gx0 + x): a cell's fit is centred on its painted pixel nearest the cell centre
 * (ties: lowest index), so the extrapolation does not depend on the tiling. */
static long rl_cell_centre(const uint8_t *state, long W, long H, long gx0, long gcx, long cyc, long cell)
{
    long best = -1, x = 0, y = 0;
    double bd = 1e300, mx = (gcx + 0.5) * cell - 0.5, my = (cyc + 0.5) * cell - 0.5;
    for (y = cyc * cell; y < (cyc + 1) * cell && y < H; y++)
        for (x = gcx * cell - gx0; x < (gcx + 1) * cell - gx0; x++) {
            size_t i = 0;
            double dd = 0.0;
            if (x < 0 || x >= W) continue;
            i = (size_t)(y * W + x);
            if (state[i] < 1 || state[i] > 3) continue;
            dd = ((double)(x + gx0) - mx) * ((double)(x + gx0) - mx) + ((double)y - my) * ((double)y - my);
            if (dd < bd) { bd = dd; best = (long)i; }
        }
    return best;
}

static size_t rl_apron(uint8_t *state, float *pos, float *nrm, size_t w, size_t h, size_t gx0,
                       size_t c0, size_t c1, const RawtexLayersOpts *o, size_t *rejected, size_t *fits)
{
    const long W = (long)w, H = (long)h, cell = o->apron_cell_px;
    const long gc0 = (long)gx0 / cell;
    const size_t cw = (size_t)((long)(gx0 + w - 1) / cell - gc0 + 1), ch = (h + (size_t)cell - 1) / (size_t)cell;
    const double reach2 = (double)o->apron_px * o->apron_px, inv = 1.0 / (double)o->apron_fit_px;
    int32_t *near_px = (int32_t *)malloc(w * h * sizeof *near_px);
    float *d2 = (float *)malloc(w * h * sizeof *d2);
    uint8_t *fstate = (uint8_t *)calloc(cw * ch, 1);
    double *fc = (double *)malloc(cw * ch * 18 * sizeof *fc);
    int32_t *fcentre = (int32_t *)malloc(cw * ch * sizeof *fcentre);
    size_t done = 0;
    long x = 0, y = 0;
    if (near_px == NULL || d2 == NULL || fstate == NULL || fc == NULL || fcentre == NULL) {
        free(near_px); free(d2); free(fstate); free(fc); free(fcentre);
        return 0;
    }
    rl_nearest(state, w, h, near_px, d2);
    for (y = 0; y < H; y++)
        for (x = (long)c0; x < (long)c1; x++) {
            size_t i = (size_t)(y * W + x), b = 0, cidx = 0;
            long bx = 0, by = 0, cxp = 0, cyp = 0;
            double u = 0.0, v = 0.0, phi[6], p[3], du[3], dv[3], nn[3], nl = 0.0, dist3 = 0.0, dpx = 0.0;
            double lu = 0.0, lv = 0.0, dotn = 0.0;
            const double *c = NULL;
            int k = 0, j = 0;
            if (state[i] != 0 || near_px[i] < 0 || (double)d2[i] > reach2) continue;
            b = (size_t)near_px[i];
            by = (long)(b / w); bx = (long)(b % w);
            cidx = (size_t)(by / cell) * cw + (size_t)(((long)gx0 + bx) / cell - gc0);
            if (fstate[cidx] == 0) {
                long cc = rl_cell_centre(state, W, H, (long)gx0, ((long)gx0 + bx) / cell, by / cell, cell);
                fstate[cidx] = 2;
                if (cc >= 0) {
                    fcentre[cidx] = (int32_t)cc;
                    if (RawtexLayers_quadric_fit(state, pos, W, H, cc % W, cc / W, o->apron_fit_px,
                                                 o->apron_step_px, &fc[cidx * 18]) == 0)
                        fstate[cidx] = 1;
                }
                (*fits)++;
            }
            if (fstate[cidx] != 1) { (*rejected)++; continue; }
            c = &fc[cidx * 18];
            cyp = (long)((size_t)fcentre[cidx] / w); cxp = (long)((size_t)fcentre[cidx] % w);
            u = (double)(x - cxp) * inv; v = (double)(y - cyp) * inv;
            phi[0] = 1.0; phi[1] = u; phi[2] = v; phi[3] = u * u; phi[4] = u * v; phi[5] = v * v;
            for (k = 0; k < 3; k++) {
                p[k] = 0.0;
                for (j = 0; j < 6; j++) p[k] += c[k * 6 + j] * phi[j];
                du[k] = (c[k * 6 + 1] + 2.0 * c[k * 6 + 3] * u + c[k * 6 + 4] * v) * inv;
                dv[k] = (c[k * 6 + 2] + c[k * 6 + 4] * u + 2.0 * c[k * 6 + 5] * v) * inv;
            }
            nn[0] = du[1] * dv[2] - du[2] * dv[1];
            nn[1] = du[2] * dv[0] - du[0] * dv[2];
            nn[2] = du[0] * dv[1] - du[1] * dv[0];
            nl = sqrt(nn[0] * nn[0] + nn[1] * nn[1] + nn[2] * nn[2]);
            lu = sqrt(du[0] * du[0] + du[1] * du[1] + du[2] * du[2]);
            lv = sqrt(dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2]);
            for (k = 0; k < 3; k++) {
                double e = p[k] - (double)pos[b * 3 + (size_t)k];
                dist3 += e * e;
            }
            dist3 = sqrt(dist3);
            dpx = sqrt((double)d2[i]);
            if (nl > 1e-9) {
                for (k = 0; k < 3; k++) nn[k] /= nl;
                dotn = nn[0] * nrm[b * 3] + nn[1] * nrm[b * 3 + 1] + nn[2] * nrm[b * 3 + 2];
                dotn /= sqrt((double)nrm[b * 3] * nrm[b * 3] + (double)nrm[b * 3 + 1] * nrm[b * 3 + 1] +
                             (double)nrm[b * 3 + 2] * nrm[b * 3 + 2]) + 1e-12;
                if (dotn < 0.0) { for (k = 0; k < 3; k++) nn[k] = -nn[k]; dotn = -dotn; }
            }
            if (!(nl > 1e-9) || dotn < 0.5 || lu < 0.5 || lu > 2.0 || lv < 0.5 || lv > 2.0 ||
                dist3 > 1.6 * dpx + 2.0 || dist3 < 0.5 * dpx - 2.0) {
                (*rejected)++;
                continue;
            }
            for (k = 0; k < 3; k++) {
                pos[i * 3 + (size_t)k] = (float)p[k];
                nrm[i * 3 + (size_t)k] = (float)nn[k];
            }
            state[i] = 4;
            done++;
        }
    free(near_px); free(d2); free(fstate); free(fc); free(fcentre);
    return done;
}

/* ------------------------------------------------------------ surface regularisation */

/* Kernel half-width of the regularisation: 3 sigma, at least one pixel. */
static int rl_smooth_radius(double sigma)
{
    int r = (int)ceil(3.0 * sigma);
    return r < 1 ? 1 : r;
}

/* One separable pass of an unnormalised Gaussian over a W x H float image,
 * zero outside: along rows (horizontal) or columns (vertical). */
static void rl_gauss_pass(const float *in, float *out, size_t W, size_t H, const double *g, int R, int vertical)
{
    size_t x = 0, y = 0;
    if (!vertical) {
        for (y = 0; y < H; y++) {
            const float *row = in + y * W;
            float *o = out + y * W;
            for (x = 0; x < W; x++) {
                long k0 = (long)x - R < 0 ? -(long)x : -R;
                long k1 = (long)x + R >= (long)W ? (long)W - 1 - (long)x : R;
                double s = 0.0;
                long k = 0;
                for (k = k0; k <= k1; k++) s += g[k + R] * (double)row[(long)x + k];
                o[x] = (float)s;
            }
        }
    } else {
        for (y = 0; y < H; y++) {
            long k0 = (long)y - R < 0 ? -(long)y : -R;
            long k1 = (long)y + R >= (long)H ? (long)H - 1 - (long)y : R;
            long k = 0;
            float *o = out + y * W;
            for (x = 0; x < W; x++) o[x] = 0.0f;
            for (k = k0; k <= k1; k++) {
                const float *row = in + (size_t)((long)y + k) * W;
                const float gk = (float)g[k + R];
                for (x = 0; x < W; x++) o[x] += gk * row[x];
            }
        }
    }
}

int RawtexLayers_smooth(const uint8_t *state, float *pos, float *nrm, size_t W, size_t H,
                        double sigma_px, double cap_vox, size_t c0, size_t c1,
                        size_t *npx, size_t *ncapped, double *shift2, double *turn)
{
    const size_t n = W * H;
    const int R = rl_smooth_radius(sigma_px);
    double *g = NULL;
    float *tmp = NULL, *acc = NULL;
    double ref[3] = { 0.0, 0.0, 0.0 };
    size_t i = 0;
    int ch = 0, have_ref = 0;
    if (!(sigma_px > 0.0) || !(cap_vox >= 0.0)) return -1;
    if (n == 0) return 0;
    g = (double *)malloc((size_t)(2 * R + 1) * sizeof *g);
    tmp = (float *)malloc(2 * n * sizeof *tmp);
    acc = (float *)malloc(7 * n * sizeof *acc);      /* mask, point (3, about ref), normal (3) */
    if (g == NULL || tmp == NULL || acc == NULL) { free(g); free(tmp); free(acc); return -1; }
    for (ch = -R; ch <= R; ch++) g[ch + R] = exp(-0.5 * (double)ch * (double)ch / (sigma_px * sigma_px));
    /* points about a tile reference keep float sums exact to ~1e-3 vox at scroll coordinates */
    for (i = 0; i < n && !have_ref; i++)
        if (state[i] >= 1 && state[i] <= 3) {
            ref[0] = pos[i * 3]; ref[1] = pos[i * 3 + 1]; ref[2] = pos[i * 3 + 2];
            have_ref = 1;
        }
    if (!have_ref) { free(g); free(tmp); free(acc); return 0; }
    for (ch = 0; ch < 7; ch++) {
        float *src = tmp, *mid = tmp + n;
        for (i = 0; i < n; i++) {
            int on = state[i] >= 1 && state[i] <= 3;
            float v = 0.0f;
            if (on) {
                if (ch == 0) v = 1.0f;
                else if (ch <= 3) v = (float)((double)pos[i * 3 + (size_t)(ch - 1)] - ref[ch - 1]);
                else {
                    const float *nn = &nrm[i * 3];
                    double nl = sqrt((double)nn[0] * nn[0] + (double)nn[1] * nn[1] + (double)nn[2] * nn[2]);
                    v = nl > 1e-12 ? (float)((double)nn[ch - 4] / nl) : 0.0f;
                }
            }
            src[i] = v;
        }
        rl_gauss_pass(src, mid, W, H, g, R, 0);
        rl_gauss_pass(mid, acc + (size_t)ch * n, W, H, g, R, 1);
    }
    for (i = 0; i < n; i++) {
        const float m = acc[i];
        double pb[3], nb[3], no[3], nl = 0.0, nol = 0.0, sh = 0.0, c = 0.0;
        int k = 0, capped = 0;
        if (state[i] < 1 || state[i] > 3 || !(m > 1e-6f)) continue;
        for (k = 0; k < 3; k++) {
            pb[k] = ref[k] + (double)acc[(size_t)(1 + k) * n + i] / (double)m;
            nb[k] = (double)acc[(size_t)(4 + k) * n + i] / (double)m;
            no[k] = (double)nrm[i * 3 + (size_t)k];
        }
        nl = sqrt(nb[0] * nb[0] + nb[1] * nb[1] + nb[2] * nb[2]);
        nol = sqrt(no[0] * no[0] + no[1] * no[1] + no[2] * no[2]);
        if (nl < 0.5 || nol < 1e-12) continue;          /* normals cancel (a fold): leave the pixel */
        for (k = 0; k < 3; k++) { nb[k] /= nl; no[k] /= nol; }
        c = nb[0] * no[0] + nb[1] * no[1] + nb[2] * no[2];
        if (c < 0.0) continue;                           /* never turn a pixel's depth axis over */
        for (k = 0; k < 3; k++) sh += (pb[k] - (double)pos[i * 3 + (size_t)k]) * nb[k];
        if (sh > cap_vox) { sh = cap_vox; capped = 1; }
        if (sh < -cap_vox) { sh = -cap_vox; capped = 1; }
        for (k = 0; k < 3; k++) {
            pos[i * 3 + (size_t)k] = (float)((double)pos[i * 3 + (size_t)k] + sh * nb[k]);
            nrm[i * 3 + (size_t)k] = (float)nb[k];
        }
        if ((i % W) >= c0 && (i % W) < c1) {
            if (npx != NULL) (*npx)++;
            if (ncapped != NULL && capped) (*ncapped)++;
            if (shift2 != NULL) *shift2 += sh * sh;
            if (turn != NULL) *turn += acos(c > 1.0 ? 1.0 : c) * 180.0 / 3.14159265358979323846;
        }
    }
    free(g); free(tmp); free(acc);
    return 0;
}

/* ------------------------------------------------------------ step blending */

static double rl_smoothstep01(double t)
{
    t = t < 0.0 ? 0.0 : (t > 1.0 ? 1.0 : t);
    return t * t * (3.0 - 2.0 * t);
}

int RawtexLayers_stepblend(const uint8_t *state, float *pos, float *nrm, size_t W, size_t H,
                           double sigma_px, double spread_px, size_t c0, size_t c1,
                           size_t *nmoved, double *move2)
{
    const size_t n = W * H;
    const int R = rl_smooth_radius(sigma_px), Rs = rl_smooth_radius(spread_px);
    double *g = NULL, *gs = NULL, ref[3] = { 0.0, 0.0, 0.0 };
    float *str = NULL, *wgt = NULL, *tmp = NULL;
    size_t i = 0;
    int ch = 0, any = 0, have_ref = 0;
    const double line_norm = 1.0 / (sqrt(2.0 * 3.14159265358979323846) * spread_px);
    if (!(sigma_px > 0.0) || !(spread_px > 0.0)) return -1;
    if (n == 0) return 0;
    str = (float *)calloc(n, sizeof *str);
    if (str == NULL) return -1;
    /* step strength: the largest smoothstep((|dP| - 1) / 2) to a painted 4-neighbour */
    for (i = 0; i < n; i++) {
        size_t x = i % W, y = i / W, nb[4];
        int k = 0, m = 0;
        double s = 0.0;
        if (state[i] < 1 || state[i] > 3) continue;
        if (x > 0) nb[m++] = i - 1;
        if (x + 1 < W) nb[m++] = i + 1;
        if (y > 0) nb[m++] = i - W;
        if (y + 1 < H) nb[m++] = i + W;
        for (k = 0; k < m; k++) {
            size_t j = nb[k];
            double d = 0.0, e = 0.0;
            int q = 0;
            if (state[j] < 1 || state[j] > 3) continue;
            for (q = 0; q < 3; q++) {
                e = (double)pos[j * 3 + (size_t)q] - (double)pos[i * 3 + (size_t)q];
                d += e * e;
            }
            d = rl_smoothstep01((sqrt(d) - 1.0) / 2.0);
            if (d > s) s = d;
        }
        str[i] = (float)s;
        if (s > 0.0) any = 1;
    }
    if (!any) { free(str); return 0; }
    /* the blend target is a Gaussian-weighted LOCAL LINEAR fit of the painted points around each pixel (its
     * value at the pixel), not their mean: a mean is biased toward the interior near any border, a linear
     * fit reproduces a flat parametrisation exactly; its slopes give the target normal */
    {
        /* kernels: G, G*k, G*k^2 (k = pixel offset) for the fit, Gs for spreading the step strength */
        double *gk = NULL, *gk2 = NULL;
        float *mom = NULL;                 /* 15 channels: m{1,x,y,xx,xy,yy}, P_c{1,x,y} for c = 0..2 */
        const double *kx[15], *ky[15];
        int kk = 0;
        g = (double *)malloc((size_t)(2 * R + 1) * sizeof *g);
        gk = (double *)malloc((size_t)(2 * R + 1) * sizeof *gk);
        gk2 = (double *)malloc((size_t)(2 * R + 1) * sizeof *gk2);
        gs = (double *)malloc((size_t)(2 * Rs + 1) * sizeof *gs);
        wgt = (float *)malloc(n * sizeof *wgt);
        tmp = (float *)malloc(2 * n * sizeof *tmp);
        mom = (float *)malloc(15 * n * sizeof *mom);
        if (g == NULL || gk == NULL || gk2 == NULL || gs == NULL || wgt == NULL || tmp == NULL || mom == NULL) {
            free(str); free(g); free(gk); free(gk2); free(gs); free(wgt); free(tmp); free(mom);
            return -1;
        }
        for (ch = -R; ch <= R; ch++) {
            double e = exp(-0.5 * (double)ch * (double)ch / (sigma_px * sigma_px));
            g[ch + R] = e; gk[ch + R] = e * (double)ch; gk2[ch + R] = e * (double)ch * (double)ch;
        }
        for (ch = -Rs; ch <= Rs; ch++) gs[ch + Rs] = exp(-0.5 * (double)ch * (double)ch / (spread_px * spread_px));
        /* blend weight: the strength spread by a Gaussian, a straight seam line weighing 1 */
        rl_gauss_pass(str, tmp, W, H, gs, Rs, 0);
        rl_gauss_pass(tmp, wgt, W, H, gs, Rs, 1);
        for (i = 0; i < n; i++) {
            double w = (double)wgt[i] * line_norm;
            wgt[i] = (float)(w > 1.0 ? 1.0 : w);
        }
        for (i = 0; i < n && !have_ref; i++)
            if (state[i] >= 1 && state[i] <= 3) {
                ref[0] = pos[i * 3]; ref[1] = pos[i * 3 + 1]; ref[2] = pos[i * 3 + 2];
                have_ref = 1;
            }
        /* channel kernels (horizontal x, vertical y): mask moments, then per point component 1, x, y */
        kx[0] = g;   ky[0] = g;   kx[1] = gk;  ky[1] = g;   kx[2] = g;   ky[2] = gk;
        kx[3] = gk2; ky[3] = g;   kx[4] = gk;  ky[4] = gk;  kx[5] = g;   ky[5] = gk2;
        for (kk = 0; kk < 3; kk++) {
            kx[6 + 3 * kk] = g;  ky[6 + 3 * kk] = g;
            kx[7 + 3 * kk] = gk; ky[7 + 3 * kk] = g;
            kx[8 + 3 * kk] = g;  ky[8 + 3 * kk] = gk;
        }
        for (ch = 0; ch < 15; ch++) {
            float *src = tmp, *mid = tmp + n;
            const int comp = ch < 6 ? -1 : (ch - 6) / 3;
            for (i = 0; i < n; i++) {
                int on = state[i] >= 1 && state[i] <= 3;
                src[i] = !on ? 0.0f : (comp < 0 ? 1.0f : (float)((double)pos[i * 3 + (size_t)comp] - ref[comp]));
            }
            rl_gauss_pass(src, mid, W, H, kx[ch], R, 0);
            rl_gauss_pass(mid, mom + (size_t)ch * n, W, H, ky[ch], R, 1);
        }
        for (i = 0; i < n; i++) {
            const double w = (double)wgt[i], m0 = (double)mom[i];
            double mx = 0.0, my = 0.0, cxx = 0.0, cxy = 0.0, cyy = 0.0, det = 0.0, a[3], bx[3], by[3], mv = 0.0;
            double no[3], nf[3], ol = 0.0, fl = 0.0;
            int k = 0;
            if (state[i] < 1 || state[i] > 3 || w < 1e-3 || !(m0 > 1e-6)) continue;
            mx = (double)mom[n + i] / m0; my = (double)mom[2 * n + i] / m0;
            cxx = (double)mom[3 * n + i] / m0 - mx * mx;
            cxy = (double)mom[4 * n + i] / m0 - mx * my;
            cyy = (double)mom[5 * n + i] / m0 - my * my;
            det = cxx * cyy - cxy * cxy;
            if (!(det > 1e-3 * (cxx + cyy) * (cxx + cyy) && cxx > 0.25 && cyy > 0.25)) continue;  /* no 2-D support */
            for (k = 0; k < 3; k++) {
                double s0 = (double)mom[(size_t)(6 + 3 * k) * n + i] / m0;
                double sx = (double)mom[(size_t)(7 + 3 * k) * n + i] / m0 - s0 * mx;
                double sy = (double)mom[(size_t)(8 + 3 * k) * n + i] / m0 - s0 * my;
                bx[k] = (cyy * sx - cxy * sy) / det;
                by[k] = (cxx * sy - cxy * sx) / det;
                a[k] = ref[k] + s0 - bx[k] * mx - by[k] * my;       /* the fit's value at the pixel */
            }
            for (k = 0; k < 3; k++) {
                double d = w * (a[k] - (double)pos[i * 3 + (size_t)k]);
                pos[i * 3 + (size_t)k] = (float)((double)pos[i * 3 + (size_t)k] + d);
                mv += d * d;
                no[k] = (double)nrm[i * 3 + (size_t)k];
            }
            /* the fit's normal, in the stored (z, y, x) order, oriented like the pixel's own */
            nf[0] = bx[1] * by[2] - bx[2] * by[1];
            nf[1] = bx[2] * by[0] - bx[0] * by[2];
            nf[2] = bx[0] * by[1] - bx[1] * by[0];
            ol = sqrt(no[0] * no[0] + no[1] * no[1] + no[2] * no[2]);
            fl = sqrt(nf[0] * nf[0] + nf[1] * nf[1] + nf[2] * nf[2]);
            if (ol > 1e-12 && fl > 1e-12) {
                double c = 0.0, mix[3], ml = 0.0;
                for (k = 0; k < 3; k++) { no[k] /= ol; nf[k] /= fl; c += no[k] * nf[k]; }
                if (c < 0.0) for (k = 0; k < 3; k++) nf[k] = -nf[k];
                for (k = 0; k < 3; k++) mix[k] = (1.0 - w) * no[k] + w * nf[k];
                ml = sqrt(mix[0] * mix[0] + mix[1] * mix[1] + mix[2] * mix[2]);
                if (ml > 0.5)
                    for (k = 0; k < 3; k++) nrm[i * 3 + (size_t)k] = (float)(mix[k] / ml);
            }
            if (w > 0.01 && (i % W) >= c0 && (i % W) < c1) {
                if (nmoved != NULL) (*nmoved)++;
                if (move2 != NULL) *move2 += mv;
            }
        }
        free(gk); free(gk2); free(mom);
    }
    free(str); free(g); free(gs); free(wgt); free(tmp);
    return 0;
}

/* ------------------------------------------------------------ handedness census */

void RawtexHandedness_free(RawtexHandedness *h)
{
    if (h == NULL) return;
    free(h->pieces);
    h->pieces = NULL;
    h->npieces = 0;
}

/* One neighbour difference along a row or column of a compacted core tile: the
 * painted neighbour's point minus the pixel's, NaN-free, or 0 when the
 * neighbour is unpainted or more than 2 vox away (a 3-D step). */
static int rl_hand_diff(const uint8_t *state, const float *pos, size_t i, size_t j, double d[3])
{
    double l2 = 0.0;
    int k = 0;
    if (state[j] < 1 || state[j] > 3) return 0;
    for (k = 0; k < 3; k++) {
        d[k] = (double)pos[j * 3 + (size_t)k] - (double)pos[i * 3 + (size_t)k];
        l2 += d[k] * d[k];
    }
    return l2 <= 4.0;
}

/* s = (dP/dcol x dP/drow) . n at core pixel i of a tw-wide, H-high compacted
 * tile; 0 when a direction has no usable neighbour.  n need not be unit. */
static double rl_hand_sign(const uint8_t *state, const float *pos, const float *nrm, size_t i, size_t tw, size_t H)
{
    double ex[3] = { 0, 0, 0 }, ey[3] = { 0, 0, 0 }, a[3], b[3], c[3];
    size_t x = i % tw, y = i / tw;
    int ha = 0, hb = 0, k = 0;
    ha = x > 0 && rl_hand_diff(state, pos, i, i - 1, a);
    hb = x + 1 < tw && rl_hand_diff(state, pos, i, i + 1, b);
    if (!ha && !hb) return 0.0;
    for (k = 0; k < 3; k++) ex[k] = ha && hb ? 0.5 * (b[k] - a[k]) : (hb ? b[k] : -a[k]);
    ha = y > 0 && rl_hand_diff(state, pos, i, i - tw, a);
    hb = y + 1 < H && rl_hand_diff(state, pos, i, i + tw, b);
    if (!ha && !hb) return 0.0;
    for (k = 0; k < 3; k++) ey[k] = ha && hb ? 0.5 * (b[k] - a[k]) : (hb ? b[k] : -a[k]);
    c[0] = ex[1] * ey[2] - ex[2] * ey[1];
    c[1] = ex[2] * ey[0] - ex[0] * ey[2];
    c[2] = ex[0] * ey[1] - ex[1] * ey[0];
    return c[0] * (double)nrm[i * 3] + c[1] * (double)nrm[i * 3 + 1] + c[2] * (double)nrm[i * 3 + 2];
}

/* Tiles whose halo-extended column range [t*cols - halo, (t+1)*cols + halo) a face's centre bbox meets. */
static size_t rl_tile_lo(const RlFace *F, size_t halo, int cols)
{
    size_t x = (size_t)F->x0 > halo ? (size_t)F->x0 - halo : 0;
    return x / (size_t)cols;
}

static size_t rl_tile_hi(const RlFace *F, size_t halo, int cols, size_t ntiles)
{
    size_t t = ((size_t)F->x1 + halo) / (size_t)cols;
    return t < ntiles ? t : ntiles - 1;
}

static int rl_write_text(const char *path, const char *text)
{
    FILE *fp = fopen(path, "wb");
    size_t n = strlen(text);
    if (fp == NULL) return -1;
    if (fwrite(text, 1, n, fp) != n) { fclose(fp); return -1; }
    return fclose(fp) == 0 ? 0 : -1;
}

/* Geometry chunk: channels (z, y, x, nz, ny, nx) of core columns [x0, x0+xw)
 * and rows [y0, y0+yh) of a tile (tw columns, npx pixels), NaN-padded to C x C. */
static int rl_write_geom_chunk(const char *path, const float *geo, size_t npx, size_t tw,
                               size_t y0, size_t yh, size_t x0, size_t xw, int C, float *gbuf)
{
    const size_t cc = (size_t)C * (size_t)C;
    size_t ch = 0, yy = 0, j = 0;
    FILE *fp = NULL;
    for (j = 0; j < 6 * cc; j++) gbuf[j] = NAN;
    for (ch = 0; ch < 6; ch++)
        for (yy = 0; yy < yh; yy++)
            memcpy(&gbuf[ch * cc + yy * (size_t)C], &geo[ch * npx + (y0 + yy) * tw + x0], xw * sizeof(float));
    fp = fopen(path, "wb");
    if (fp == NULL) return -1;
    if (fwrite(gbuf, sizeof(float), 6 * cc, fp) != 6 * cc) { fclose(fp); return -1; }
    return fclose(fp) == 0 ? 0 : -1;
}

static int rl_mkdir(const char *path)
{
    if (ves_mkdir(path) == 0) return 0;
    {   /* an existing directory is fine */
        char probe[2300];
        FILE *fp = NULL;
        snprintf(probe, sizeof probe, "%s/.rl_probe", path);
        fp = fopen(probe, "wb");
        if (fp == NULL) return -1;
        fclose(fp);
        remove(probe);
    }
    return 0;
}

int RawtexLayers_write(const char *zarr_dir, const char *mask_path,
                       CubeTable *ct, const float *verts, RawtexUv uv,
                       size_t nv, const int32_t *faces, size_t nf,
                       const float *normals, const uint8_t *face_skip,
                       const RawtexLayersOpts *o, const char *attrs_json,
                       RawtexLayersStats *st)
{
    RawtexPlan plan;
    size_t W = 0, H = 0, ntiles = 0, f = 0, t = 0;
    size_t *tile_start = NULL, *tile_fill = NULL;
    int32_t *tile_faces = NULL;
    uint8_t *paint = NULL, *mask = NULL;
    RlFace *rf = NULL;
    uint32_t *hc_team = NULL, *hc_mirror = NULL;    /* handedness census per column (o->hand) */
    int32_t *hc_r0 = NULL, *hc_r1 = NULL;
    char partial[2300], path[2400], gpartial[2300];
    size_t painted = 0, filled = 0, complete = 0, multi = 0, skip_uv = 0, chunks = 0;
    size_t apron = 0, apron_rej = 0, apron_fits = 0;
    size_t smooth_px = 0, smooth_cap = 0, step_px = 0;
    double smooth_sh2 = 0.0, smooth_turn = 0.0, step_m2 = 0.0;
    size_t halo = 0;
    int failed = 0;
    double t0 = ves_clock_sec();
    FILE *probe = NULL;

    if (st != NULL) memset(st, 0, sizeof *st);
    if (zarr_dir == NULL || mask_path == NULL || ct == NULL || verts == NULL ||
        faces == NULL || normals == NULL || o == NULL || !Rawtex_has_uv(uv) ||
        o->depth < 1 || o->depth > 255 || !(o->step > 0.0) || !isfinite(o->offset) ||
        o->chunk < 1 || o->tile_cols < o->chunk || o->tile_cols % o->chunk != 0 ||
        o->seam_reach < 0 || o->seam_passes < 0 || (o->seam_reach > 0 && !(o->seam_tol > 0.0)) ||
        o->apron_px < 0 || (o->apron_px > 0 && (o->apron_fit_px < 2 || o->apron_cell_px < 1 ||
                                                o->apron_step_px < 1)) ||
        !(o->smooth_sigma_px >= 0.0) || (o->smooth_sigma_px > 0.0 && !(o->smooth_cap_vox > 0.0)) ||
        !(o->step_sigma_px >= 0.0) || (o->step_sigma_px > 0.0 && !(o->step_spread_px > 0.0)) ||
        nv == 0 || nf == 0)
        return -1;
    /* a fill at column c depends on columns within seam_reach per pass: a halo
     * of passes x reach rasterized columns makes the tile width invisible; an
     * apron pixel reads painted pixels up to apron_px + fit radius + one fit
     * cell away */
    halo = o->seam_reach > 0 ? (size_t)o->seam_reach * (size_t)o->seam_passes : 0;
    if (o->smooth_sigma_px > 0.0 || o->step_sigma_px > 0.0) {
        /* fills, then the regularisation kernel, then the step blend (strength, spread, mean), then the
         * apron's reach, chain */
        if (o->smooth_sigma_px > 0.0) halo += (size_t)rl_smooth_radius(o->smooth_sigma_px);
        if (o->step_sigma_px > 0.0)
            halo += 1 + (size_t)rl_smooth_radius(o->step_spread_px) + (size_t)rl_smooth_radius(o->step_sigma_px);
        if (o->apron_px > 0) halo += (size_t)(o->apron_px + o->apron_fit_px + o->apron_cell_px);
    } else if (o->apron_px > 0 && (size_t)(o->apron_px + o->apron_fit_px + o->apron_cell_px) > halo)
        halo = (size_t)(o->apron_px + o->apron_fit_px + o->apron_cell_px);
    if (o->window != NULL) {        /* an absolute window, e.g. columns of a wider sheet's full grid */
        if (Rawtex_plan_field(uv, nv, 1.0, 1.0, RL_MAX_PX, o->window, &plan) != 0 || !plan.ok) {
            fprintf(stderr, "ERROR: window raster %zux%zu exceeds the 1 px/vox caps\n", plan.W, plan.H);
            return -1;
        }
    } else if (Rawtex_plan_field(uv, nv, 1.0, 1.0, RL_MAX_PX, NULL, &plan) != 0 || !plan.ok) {
        fprintf(stderr, "ERROR: layer raster %zux%zu exceeds the 1 px/vox caps; the sheet "
                "cannot be rendered at native pitch\n", plan.W, plan.H);
        return -1;
    }
    W = plan.W; H = plan.H;
    {   /* refuse to replace a published volume */
        const char *dirs[2];
        int d = 0;
        dirs[0] = zarr_dir; dirs[1] = o->geom_dir;
        for (d = 0; d < 2; d++) {
            char zgroup[2300];
            if (dirs[d] == NULL) continue;
            snprintf(zgroup, sizeof zgroup, "%s/.zgroup", dirs[d]);
            probe = fopen(zgroup, "rb");
            if (probe != NULL) {
                fclose(probe);
                fprintf(stderr, "ERROR: refusing to replace existing %s\n", dirs[d]);
                return -1;
            }
        }
    }
    /* a fresh partial directory: stale chunks of a failed run must never
     * survive into a published volume */
    snprintf(partial, sizeof partial, "%s.partial", zarr_dir);
    snprintf(path, sizeof path, "%s/0", partial);
    if (ves_ensure_parent_dir(partial) != 0 || ves_mkdir(partial) != 0 || ves_mkdir(path) != 0) {
        fprintf(stderr, "ERROR: cannot create a fresh %s (remove a stale partial first)\n", path);
        return -1;
    }
    gpartial[0] = 0;
    if (o->geom_dir != NULL) {
        snprintf(gpartial, sizeof gpartial, "%s.partial", o->geom_dir);
        snprintf(path, sizeof path, "%s/0", gpartial);
        if (ves_ensure_parent_dir(gpartial) != 0 || ves_mkdir(gpartial) != 0 || ves_mkdir(path) != 0) {
            fprintf(stderr, "ERROR: cannot create a fresh %s (remove a stale partial first)\n", path);
            return -1;
        }
    }

    /* face prologue once; bucket paintable faces by column tile */
    ntiles = (W + (size_t)o->tile_cols - 1) / (size_t)o->tile_cols;
    rf = (RlFace *)malloc(nf * sizeof *rf);
    paint = (uint8_t *)calloc(nf, 1);
    tile_start = (size_t *)calloc(ntiles + 1, sizeof *tile_start);
    tile_fill = (size_t *)calloc(ntiles, sizeof *tile_fill);
    mask = (uint8_t *)calloc(W * H, 1);
    if (rf == NULL || paint == NULL || tile_start == NULL || tile_fill == NULL || mask == NULL) {
        failed = 1;
        goto done;
    }
    if (o->hand != NULL) {
        size_t c = 0;
        RawtexHandedness_free(o->hand);
        hc_team = (uint32_t *)calloc(W, sizeof *hc_team);
        hc_mirror = (uint32_t *)calloc(W, sizeof *hc_mirror);
        hc_r0 = (int32_t *)malloc(W * sizeof *hc_r0);
        hc_r1 = (int32_t *)malloc(W * sizeof *hc_r1);
        if (hc_team == NULL || hc_mirror == NULL || hc_r0 == NULL || hc_r1 == NULL) { failed = 1; goto done; }
        for (c = 0; c < W; c++) { hc_r0[c] = INT32_MAX; hc_r1[c] = -1; }
    }
    for (f = 0; f < nf; f++) {
        int bad = 0;
        if (face_skip != NULL && face_skip[f]) continue;
        paint[f] = (uint8_t)rl_face(verts, uv, faces, f, plan.umin, plan.vmin, W, H,
                                    o->stretch_ratio, o->stretch_floor, &rf[f], &bad);
        if (bad) skip_uv++;
        if (!paint[f]) continue;
        for (t = rl_tile_lo(&rf[f], halo, o->tile_cols); t <= rl_tile_hi(&rf[f], halo, o->tile_cols, ntiles); t++)
            tile_start[t + 1]++;
    }
    for (t = 0; t < ntiles; t++) tile_start[t + 1] += tile_start[t];
    tile_faces = (int32_t *)malloc((tile_start[ntiles] ? tile_start[ntiles] : 1) * sizeof *tile_faces);
    if (tile_faces == NULL) { failed = 1; goto done; }
    for (f = 0; f < nf; f++) {
        if (!paint[f]) continue;
        for (t = rl_tile_lo(&rf[f], halo, o->tile_cols); t <= rl_tile_hi(&rf[f], halo, o->tile_cols, ntiles); t++)
            tile_faces[tile_start[t] + tile_fill[t]++] = (int32_t)f;
    }

    {
        const int depth = o->depth, C = o->chunk;
        const size_t cols = (size_t)o->tile_cols;
        const size_t cy_n = (H + (size_t)C - 1) / (size_t)C;
        long tt = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) reduction(+:painted,filled,complete,multi,chunks,apron,apron_rej,apron_fits,smooth_px,smooth_cap,smooth_sh2,smooth_turn,step_px,step_m2)
#endif
        for (tt = 0; tt < (long)ntiles; tt++) {
            /* core columns [tx0, tx0+tw) are emitted; the raster (and the seam
             * fill) also covers `halo` columns either side, [hx0, hx0+hw) */
            size_t tx0 = (size_t)tt * cols;
            size_t tw = W - tx0 < cols ? W - tx0 : cols;
            size_t hx0 = tx0 > halo ? tx0 - halo : 0;
            size_t hx1 = tx0 + tw + halo < W ? tx0 + tw + halo : W;
            size_t hw = hx1 - hx0, core = tx0 - hx0;
            size_t hpx = hw * H, npx = tw * H, i = 0, k = 0;
            int pass = 0;
            uint8_t *state = (uint8_t *)calloc(hpx, 1);
            uint8_t *snap = (uint8_t *)malloc(hpx);
            float *pos = (float *)malloc(hpx * 3 * sizeof(float));
            float *nrm = (float *)malloc(hpx * 3 * sizeof(float));
            uint8_t *lay = (uint8_t *)calloc((size_t)depth * npx, 1);
            uint8_t *buf = (uint8_t *)malloc((size_t)depth * (size_t)C * (size_t)C);
            float *geo = o->geom_dir != NULL ? (float *)malloc(6 * npx * sizeof(float)) : NULL;
            float *gbuf = o->geom_dir != NULL ? (float *)malloc(6 * (size_t)C * (size_t)C * sizeof(float)) : NULL;
            if (state == NULL || snap == NULL || pos == NULL || nrm == NULL || lay == NULL || buf == NULL ||
                (o->geom_dir != NULL && (geo == NULL || gbuf == NULL))) {
#ifdef _OPENMP
#pragma omp critical(rl_fail)
#endif
                failed = 1;
                free(state); free(snap); free(pos); free(nrm); free(lay); free(buf); free(geo); free(gbuf);
                continue;
            }
            if (geo != NULL)
                for (i = 0; i < 6 * npx; i++) geo[i] = NAN;
            for (i = tile_start[tt]; i < tile_start[tt + 1]; i++) {
                size_t ff = (size_t)tile_faces[i];
                const RlFace *F = &rf[ff];
                size_t a = (size_t)faces[ff * 3], b = (size_t)faces[ff * 3 + 1], c = (size_t)faces[ff * 3 + 2];
                long xa = F->x0 > (long)hx0 ? F->x0 : (long)hx0;
                long xb = F->x1 < (long)hx1 - 1 ? F->x1 : (long)hx1 - 1;
                long yy = 0, xx = 0;
                for (yy = F->y0; yy <= F->y1; yy++) {
                    for (xx = xa; xx <= xb; xx++) {
                        double pu = (double)xx + 0.5, pv = (double)yy + 0.5;
                        double l1 = ((pu - F->ua) * (F->vc - F->va) - (pv - F->va) * (F->uc - F->ua)) / F->A2;
                        double l2 = ((F->ub - F->ua) * (pv - F->va) - (F->vb - F->va) * (pu - F->ua)) / F->A2;
                        double l0 = 1.0 - l1 - l2;
                        size_t pi = (size_t)yy * hw + ((size_t)xx - hx0);
                        int on_boundary = 0, q = 0;
                        if (l0 < -1e-9 || l1 < -1e-9 || l2 < -1e-9) continue;
                        on_boundary = l0 <= 1e-9 || l1 <= 1e-9 || l2 <= 1e-9;
                        if (on_boundary && state[pi] != 0) continue;
                        if (!on_boundary && state[pi] == 2) {
                            if ((size_t)xx >= tx0 && (size_t)xx < tx0 + tw) multi++;
                            continue;
                        }
                        for (q = 0; q < 3; q++) {
                            pos[pi * 3 + (size_t)q] = (float)(l0 * (double)verts[a * 3 + (size_t)q]
                                + l1 * (double)verts[b * 3 + (size_t)q] + l2 * (double)verts[c * 3 + (size_t)q]);
                            nrm[pi * 3 + (size_t)q] = (float)(l0 * (double)normals[a * 3 + (size_t)q]
                                + l1 * (double)normals[b * 3 + (size_t)q] + l2 * (double)normals[c * 3 + (size_t)q]);
                        }
                        state[pi] = on_boundary ? 1 : 2;
                    }
                }
            }
            for (pass = 0; pass < o->seam_passes && o->seam_reach > 0; pass++)
                if (rl_seam_fill(state, pos, nrm, hw, H, o->seam_reach, o->seam_tol, snap) == 0) break;
            if (o->smooth_sigma_px > 0.0) {
                size_t sp = 0, sc = 0;
                double s2 = 0.0, tu = 0.0;
                if (RawtexLayers_smooth(state, pos, nrm, hw, H, o->smooth_sigma_px, o->smooth_cap_vox,
                                        core, core + tw, &sp, &sc, &s2, &tu) != 0) {
#ifdef _OPENMP
#pragma omp critical(rl_fail)
#endif
                    failed = 1;
                }
                smooth_px += sp; smooth_cap += sc; smooth_sh2 += s2; smooth_turn += tu;
            }
            if (o->step_sigma_px > 0.0) {
                size_t sp = 0;
                double m2 = 0.0;
                if (RawtexLayers_stepblend(state, pos, nrm, hw, H, o->step_sigma_px, o->step_spread_px,
                                           core, core + tw, &sp, &m2) != 0) {
#ifdef _OPENMP
#pragma omp critical(rl_fail)
#endif
                    failed = 1;
                }
                step_px += sp; step_m2 += m2;
            }
            if (o->apron_px > 0) {
                size_t rej = 0, fitc = 0;
                apron += rl_apron(state, pos, nrm, hw, H, hx0, core, core + tw, o, &rej, &fitc);
                apron_rej += rej;
                apron_fits += fitc;
            }
            /* keep only the core columns from here on (in place, row by row) */
            for (k = 0; k < H; k++) {
                memmove(&state[k * tw], &state[k * hw + core], tw);
                memmove(&pos[k * tw * 3], &pos[(k * hw + core) * 3], tw * 3 * sizeof(float));
                memmove(&nrm[k * tw * 3], &nrm[(k * hw + core) * 3], tw * 3 * sizeof(float));
            }
            if (hc_team != NULL)            /* this tile owns core columns [tx0, tx0 + tw) */
                for (i = 0; i < npx; i++) {
                    size_t col = tx0 + i % tw;
                    int32_t row = (int32_t)(i / tw);
                    double s = 0.0, nl = 0.0;
                    if (state[i] < 1 || state[i] > 3) continue;
                    if (row < hc_r0[col]) hc_r0[col] = row;
                    if (row > hc_r1[col]) hc_r1[col] = row;
                    nl = (double)nrm[i * 3] * nrm[i * 3] + (double)nrm[i * 3 + 1] * nrm[i * 3 + 1] +
                         (double)nrm[i * 3 + 2] * nrm[i * 3 + 2];
                    if (nl < 0.25) continue;     /* degenerate normal: no depth axis to judge by */
                    s = rl_hand_sign(state, pos, nrm, i, tw, H);
                    if (s > 0.0) hc_team[col]++;
                    else if (s < 0.0) hc_mirror[col]++;
                }
            for (i = 0; i < npx; i++) {
                const float *p = &pos[i * 3];
                float *n = &nrm[i * 3];
                double nl = 0.0;
                int ok = 1, nok = 1, lk = 0;
                if (!state[i]) continue;
                if (state[i] == 3) filled++;
                else if (state[i] != 4) painted++;      /* 4 = context-only apron, counted by rl_apron */
                nl = sqrt((double)n[0] * n[0] + (double)n[1] * n[1] + (double)n[2] * n[2]);
                if (nl < 0.5) ok = nok = 0;   /* degenerate interpolated normal: no depth axis */
                else { n[0] = (float)(n[0] / nl); n[1] = (float)(n[1] / nl); n[2] = (float)(n[2] / nl); }
                /* `nok` is the depth axis, `ok` the completeness: every layer samples
                 * its own depth even after an earlier layer lacked full RAW support
                 * (the frozen build of 2026-09-29 00:0x let `ok` gate the normal, so
                 * such pixels -- outside the mask, <= 0.06% of a sheet -- repeated the
                 * depth-zero sample in their later layers; fixed 2026-09-29 01:15). */
                for (lk = 0; lk < depth; lk++) {
                    double off = ((double)lk - 0.5 * (double)(depth - 1)) * o->step + o->offset;
                    int full = 0;
                    double s = rl_sample(ct, (double)p[0] + off * (nok ? n[0] : 0.0f),
                                         (double)p[1] + off * (nok ? n[1] : 0.0f),
                                         (double)p[2] + off * (nok ? n[2] : 0.0f), &full);
                    if (s < 0.0) { s = 0.0; full = 0; }
                    if (!full) ok = 0;
                    lay[(size_t)lk * npx + i] = (uint8_t)(s > 255.0 ? 255.0 : floor(s + 0.5));
                }
                if (state[i] == 4) continue;            /* apron: context in the layers only */
                if (ok) {
                    complete++;
                    mask[(i / tw) * W + tx0 + (i % tw)] = 255;
                }
                if (geo != NULL) {
                    geo[i] = p[0]; geo[npx + i] = p[1]; geo[2 * npx + i] = p[2];
                    geo[3 * npx + i] = nok ? n[0] : 0.0f;
                    geo[4 * npx + i] = nok ? n[1] : 0.0f;
                    geo[5 * npx + i] = nok ? n[2] : 0.0f;
                }
            }
            /* chunks: [depth][C][C], written only when a pixel of theirs is painted */
            for (k = 0; k < cy_n && !failed; k++) {
                size_t cx = 0;
                for (cx = 0; cx * (size_t)C < tw; cx++) {
                    size_t y0 = k * (size_t)C, x0 = cx * (size_t)C;
                    size_t yh = H - y0 < (size_t)C ? H - y0 : (size_t)C;
                    size_t xw = tw - x0 < (size_t)C ? tw - x0 : (size_t)C;
                    size_t yy = 0, xx = 0;
                    int any = 0, lk = 0;
                    char cpath[2400];
                    FILE *fp = NULL;
                    for (yy = 0; yy < yh && !any; yy++)
                        for (xx = 0; xx < xw; xx++)
                            if (state[(y0 + yy) * tw + x0 + xx]) { any = 1; break; }
                    if (!any) continue;
                    memset(buf, 0, (size_t)depth * (size_t)C * (size_t)C);
                    for (lk = 0; lk < depth; lk++)
                        for (yy = 0; yy < yh; yy++)
                            memcpy(&buf[((size_t)lk * (size_t)C + yy) * (size_t)C],
                                   &lay[(size_t)lk * npx + (y0 + yy) * tw + x0], xw);
                    snprintf(cpath, sizeof cpath, "%s/0/0.%zu.%zu", partial, k,
                             (tx0 + x0) / (size_t)C);
                    fp = fopen(cpath, "wb");
                    if (fp == NULL || fwrite(buf, 1, (size_t)depth * (size_t)C * (size_t)C, fp)
                                          != (size_t)depth * (size_t)C * (size_t)C) {
                        if (fp != NULL) fclose(fp);
#ifdef _OPENMP
#pragma omp critical(rl_fail)
#endif
                        failed = 1;
                        break;
                    }
                    if (fclose(fp) != 0) {
#ifdef _OPENMP
#pragma omp critical(rl_fail)
#endif
                        failed = 1;
                        break;
                    }
                    if (geo != NULL) {
                        snprintf(cpath, sizeof cpath, "%s/0/0.%zu.%zu", gpartial, k, (tx0 + x0) / (size_t)C);
                        if (rl_write_geom_chunk(cpath, geo, npx, tw, y0, yh, x0, xw, C, gbuf) != 0) {
#ifdef _OPENMP
#pragma omp critical(rl_fail)
#endif
                            failed = 1;
                            break;
                        }
                    }
                    chunks++;
                }
            }
            free(state); free(snap); free(pos); free(nrm); free(lay); free(buf); free(geo); free(gbuf);
        }
    }
    if (failed) { fprintf(stderr, "ERROR: layer render or chunk write failed\n"); goto done; }
    if (hc_team != NULL) {          /* pieces = runs of painted columns */
        size_t c = 0, n = 0, cap = 64;
        RawtexPieceHand *pc = (RawtexPieceHand *)malloc(cap * sizeof *pc);
        if (pc == NULL) { failed = 1; goto done; }
        while (c < W) {
            RawtexPieceHand p;
            if (hc_r1[c] < 0) { c++; continue; }
            memset(&p, 0, sizeof p);
            p.col0 = c; p.row0 = (size_t)hc_r0[c]; p.row1 = (size_t)hc_r1[c] + 1;
            for (; c < W && hc_r1[c] >= 0; c++) {
                if ((size_t)hc_r0[c] < p.row0) p.row0 = (size_t)hc_r0[c];
                if ((size_t)hc_r1[c] + 1 > p.row1) p.row1 = (size_t)hc_r1[c] + 1;
                p.px_team += hc_team[c];
                p.px_mirror += hc_mirror[c];
            }
            p.col1 = c;
            if (n == cap) {
                RawtexPieceHand *g = (RawtexPieceHand *)realloc(pc, 2 * cap * sizeof *pc);
                if (g == NULL) { free(pc); failed = 1; goto done; }
                pc = g; cap *= 2;
            }
            pc[n++] = p;
        }
        o->hand->pieces = pc;
        o->hand->npieces = n;
    }

    snprintf(path, sizeof path, "%s/0/.zarray", partial);
    {
        char text[1024];
        snprintf(text, sizeof text,
                 "{\n  \"zarr_format\": 2,\n  \"shape\": [%d, %zu, %zu],\n  \"chunks\": [%d, %d, %d],\n"
                 "  \"dtype\": \"|u1\",\n  \"compressor\": null,\n  \"fill_value\": 0,\n"
                 "  \"order\": \"C\",\n  \"filters\": null,\n  \"dimension_separator\": \".\"\n}\n",
                 o->depth, H, W, o->depth, o->chunk, o->chunk);
        if (rl_write_text(path, text) != 0) { failed = 1; goto done; }
    }
    snprintf(path, sizeof path, "%s/.zgroup", partial);
    if (rl_write_text(path, "{\n  \"zarr_format\": 2\n}\n") != 0) { failed = 1; goto done; }
    snprintf(path, sizeof path, "%s/.zattrs", partial);
    if (rl_write_text(path, attrs_json != NULL ? attrs_json : "{}\n") != 0) { failed = 1; goto done; }
    if (H > (size_t)INT32_MAX || W > (size_t)INT32_MAX ||
        ves_ensure_parent_dir(mask_path) != 0 ||
        TiffIO_save(mask_path, mask, 1, (int)H, (int)W) != 0) {
        fprintf(stderr, "ERROR: cannot write mask %s\n", mask_path);
        failed = 1;
        goto done;
    }
    if (o->geom_dir != NULL) {
        char text[2048];
        size_t tl = 0;
        int lk = 0;
        snprintf(path, sizeof path, "%s/0/.zarray", gpartial);
        snprintf(text, sizeof text,
                 "{\n  \"zarr_format\": 2,\n  \"shape\": [6, %zu, %zu],\n  \"chunks\": [6, %d, %d],\n"
                 "  \"dtype\": \"<f4\",\n  \"compressor\": null,\n  \"fill_value\": \"NaN\",\n"
                 "  \"order\": \"C\",\n  \"filters\": null,\n  \"dimension_separator\": \".\"\n}\n",
                 H, W, o->chunk, o->chunk);
        if (rl_write_text(path, text) != 0) { failed = 1; goto done; }
        snprintf(path, sizeof path, "%s/.zgroup", gpartial);
        if (rl_write_text(path, "{\n  \"zarr_format\": 2\n}\n") != 0) { failed = 1; goto done; }
        tl += (size_t)snprintf(text + tl, sizeof text - tl,
            "{\n  \"format\": \"scrollfiesta-sheet-geometry-v1\",\n"
            "  \"channels\": [\"z\", \"y\", \"x\", \"nz\", \"ny\", \"nx\"],\n"
            "  \"units\": \"level-0 voxels of the RAW the layers were sampled from\",\n"
            "  \"normal\": \"unit, inward (the layer axis); (0,0,0) where the interpolated normal is degenerate\",\n"
            "  \"layer_rule\": \"layer k samples point + layer_offsets_voxels[k] * normal (trilinear)\",\n"
            "  \"pixel_grid\": \"the layer volume's grid; NaN where no layer column was rendered\",\n"
            "  \"layer_offsets_voxels\": [");
        for (lk = 0; lk < o->depth && tl < sizeof text; lk++)
            tl += (size_t)snprintf(text + tl, sizeof text - tl, "%s%.3f", lk ? ", " : "",
                                   ((double)lk - 0.5 * (double)(o->depth - 1)) * o->step + o->offset);
        if (tl < sizeof text) tl += (size_t)snprintf(text + tl, sizeof text - tl, "]\n}\n");
        snprintf(path, sizeof path, "%s/.zattrs", gpartial);
        if (tl >= sizeof text || rl_write_text(path, text) != 0) { failed = 1; goto done; }
        if (rename(gpartial, o->geom_dir) != 0) {
            fprintf(stderr, "ERROR: cannot publish %s -> %s\n", gpartial, o->geom_dir);
            failed = 1;
            goto done;
        }
    }
    if (rename(partial, zarr_dir) != 0) {
        fprintf(stderr, "ERROR: cannot publish %s -> %s\n", partial, zarr_dir);
        failed = 1;
        goto done;
    }

done:
    if (st != NULL) {
        st->W = W; st->H = H; st->umin = plan.umin; st->vmin = plan.vmin;
        st->painted_px = painted; st->filled_px = filled; st->complete_px = complete; st->multi_px = multi;
        st->apron_px = apron; st->apron_rejected = apron_rej; st->apron_fits = apron_fits;
        st->smooth_px = smooth_px; st->smooth_capped_px = smooth_cap;
        st->smooth_shift_rms = smooth_px ? sqrt(smooth_sh2 / (double)smooth_px) : 0.0;
        st->smooth_turn_mean = smooth_px ? smooth_turn / (double)smooth_px : 0.0;
        st->step_px = step_px;
        st->step_move_rms = step_px ? sqrt(step_m2 / (double)step_px) : 0.0;
        st->skip_uv_faces = skip_uv; st->chunks_written = chunks;
        st->seconds = ves_clock_sec() - t0;
    }
    free(rf); free(paint); free(tile_start); free(tile_fill); free(tile_faces); free(mask);
    free(hc_team); free(hc_mirror); free(hc_r0); free(hc_r1);
    if (failed && o->hand != NULL) RawtexHandedness_free(o->hand);
    return failed ? -1 : 0;
}

/* -------------------------------------------------------------- selftest */

#define RL_CHECK(cond, msg) do { if (!(cond)) { fprintf(stderr, "  FAIL: %s\n", msg); fails++; } } while (0)

/* synthetic RAW: value = 10 + 3z + 2y + x on a 16^3 nested Zarr (chunk 8);
 * trilinear sampling of a linear field is exact */
static int rl_make_zarr(const char *dir, int drop_chunk)
{
    char path[2400];
    int cz = 0, cy = 0, cx = 0;
    uint8_t cube[8 * 8 * 8];
    if (rl_mkdir(dir) != 0) return -1;
    snprintf(path, sizeof path, "%s/0", dir);
    if (rl_mkdir(path) != 0) return -1;
    snprintf(path, sizeof path, "%s/0/.zarray", dir);
    if (rl_write_text(path, "{\"zarr_format\": 2, \"shape\": [16, 16, 16], \"chunks\": [8, 8, 8], "
                            "\"dtype\": \"|u1\", \"compressor\": null, \"fill_value\": 0, "
                            "\"order\": \"C\", \"filters\": null, \"dimension_separator\": \"/\"}") != 0)
        return -1;
    for (cz = 0; cz < 2; cz++)
        for (cy = 0; cy < 2; cy++)
            for (cx = 0; cx < 2; cx++) {
                int z = 0, y = 0, x = 0;
                FILE *fp = NULL;
                if (drop_chunk && cz == 0 && cy == 1 && cx == 0) continue;
                snprintf(path, sizeof path, "%s/0/%d", dir, cz); rl_mkdir(path);
                snprintf(path, sizeof path, "%s/0/%d/%d", dir, cz, cy); rl_mkdir(path);
                for (z = 0; z < 8; z++)
                    for (y = 0; y < 8; y++)
                        for (x = 0; x < 8; x++)
                            cube[(z * 8 + y) * 8 + x] = (uint8_t)(10 + 3 * (cz * 8 + z) + 2 * (cy * 8 + y) + (cx * 8 + x));
                snprintf(path, sizeof path, "%s/0/%d/%d/%d", dir, cz, cy, cx);
                fp = fopen(path, "wb");
                if (fp == NULL) return -1;
                fwrite(cube, 1, sizeof cube, fp);
                fclose(fp);
            }
    return 0;
}

static int rl_read_chunk(const char *zarr, size_t cy, size_t cx, uint8_t *buf, size_t n)
{
    char path[2400];
    FILE *fp = NULL;
    size_t got = 0;
    snprintf(path, sizeof path, "%s/0/0.%zu.%zu", zarr, cy, cx);
    fp = fopen(path, "rb");
    if (fp == NULL) return -1;
    got = fread(buf, 1, n, fp);
    fclose(fp);
    return got == n ? 0 : -1;
}

static void rl_rmtree_zarr(const char *zarr, size_t ncy, size_t ncx)
{
    char path[2400];
    size_t a = 0, b = 0;
    for (a = 0; a < ncy; a++)
        for (b = 0; b < ncx; b++) {
            snprintf(path, sizeof path, "%s/0/0.%zu.%zu", zarr, a, b);
            remove(path);
        }
    snprintf(path, sizeof path, "%s/0/.zarray", zarr); remove(path);
    snprintf(path, sizeof path, "%s/0", zarr); ves_rmdir(path);
    snprintf(path, sizeof path, "%s/.zgroup", zarr); remove(path);
    snprintf(path, sizeof path, "%s/.zattrs", zarr); remove(path);
    ves_rmdir(zarr);
}

static void rl_rm_raw(const char *dir)
{
    char path[2400];
    int cz = 0, cy = 0, cx = 0;
    for (cz = 0; cz < 2; cz++) {
        for (cy = 0; cy < 2; cy++) {
            for (cx = 0; cx < 2; cx++) {
                snprintf(path, sizeof path, "%s/0/%d/%d/%d", dir, cz, cy, cx);
                remove(path);
            }
            snprintf(path, sizeof path, "%s/0/%d/%d", dir, cz, cy); ves_rmdir(path);
        }
        snprintf(path, sizeof path, "%s/0/%d", dir, cz); ves_rmdir(path);
    }
    snprintf(path, sizeof path, "%s/0/.zarray", dir); remove(path);
    snprintf(path, sizeof path, "%s/0", dir); ves_rmdir(path);
    ves_rmdir(dir);
}

int RawtexLayers_selftest(void)
{
    int fails = 0;
    const char *root = "rawtex_layers_selftest";
    char raw[256], raw_gap[256], zarr[256], zarr2[256], zarr3[256], maskp[256], mask2[256], mask3[256];
    char geom[256], zarr6[256], mask6[256];
    /* the plane x = 7.5 carries a 10 x 10 quad in (z,y) at [2,12]^2;
     * uv = (y - 2, z - 2) so pixel (px,py) has centre y = 2.5 + px, z = 2.5 + py */
    float verts[12] = { 2.0f, 2.0f, 7.5f,   2.0f, 12.0f, 7.5f,
                        12.0f, 12.0f, 7.5f, 12.0f, 2.0f, 7.5f };
    double uvd[8] = { 0.0, 0.0,  10.0, 0.0,  10.0, 10.0,  0.0, 10.0 };
    /* winding: cross((0,10,0),(10,10,0)) = (0,0,-100): the winding normal is -x */
    int32_t faces[6] = { 0, 1, 2,  0, 2, 3 };
    float inward[12];
    float *normals = NULL;
    RawtexOrientStats os;
    RawtexLayersStats ls, ls2;
    RawtexLayersOpts opt;
    RawtexUv uv = { NULL, uvd };
    CubeTable ct;
    Arena_T arena = NULL;
    int i = 0;

    snprintf(raw, sizeof raw, "%s/raw.zarr", root);
    snprintf(raw_gap, sizeof raw_gap, "%s/raw_gap.zarr", root);
    snprintf(zarr, sizeof zarr, "%s/layers.zarr", root);
    snprintf(zarr2, sizeof zarr2, "%s/layers_tiled.zarr", root);
    snprintf(zarr3, sizeof zarr3, "%s/layers_gap.zarr", root);
    snprintf(maskp, sizeof maskp, "%s/mask.tif", root);
    snprintf(mask2, sizeof mask2, "%s/mask_tiled.tif", root);
    snprintf(mask3, sizeof mask3, "%s/mask_gap.tif", root);
    snprintf(geom, sizeof geom, "%s/geom.zarr", root);
    snprintf(zarr6, sizeof zarr6, "%s/layers_gap_out.zarr", root);
    snprintf(mask6, sizeof mask6, "%s/mask_gap_out.tif", root);
    rl_mkdir(root);
    {   /* leftovers of an interrupted earlier selftest */
        char stale[300];
        char z3[8][256];
        int j = 0;
        snprintf(z3[0], sizeof z3[0], "%s", zarr); snprintf(z3[1], sizeof z3[1], "%s", zarr2);
        snprintf(z3[2], sizeof z3[2], "%s", zarr3); snprintf(z3[3], sizeof z3[3], "%s/seam16.zarr", root);
        snprintf(z3[4], sizeof z3[4], "%s/seam8.zarr", root); snprintf(z3[5], sizeof z3[5], "%s/seamx.zarr", root);
        snprintf(z3[6], sizeof z3[6], "%s", geom); snprintf(z3[7], sizeof z3[7], "%s", zarr6);
        for (j = 0; j < 8; j++) {
            rl_rmtree_zarr(z3[j], 2, 2);
            snprintf(stale, sizeof stale, "%s.partial", z3[j]);
            rl_rmtree_zarr(stale, 2, 2);
        }
        remove(maskp); remove(mask2); remove(mask3); remove(mask6);
        rl_rm_raw(raw);
        rl_rm_raw(raw_gap);
    }
    if (rl_make_zarr(raw, 0) != 0 || rl_make_zarr(raw_gap, 1) != 0) {
        fprintf(stderr, "[selftest] rawtex_layers: cannot build the synthetic RAW\n");
        return 1;
    }
    /* the axis at x = 0: inward is -x everywhere */
    for (i = 0; i < 4; i++) { inward[i * 3] = 0.0f; inward[i * 3 + 1] = 0.0f; inward[i * 3 + 2] = -1.0f; }

    /* 1. orientation keeps a component that already faces inward ... */
    normals = MeshNormals_compute(verts, 4, faces, 2);
    RL_CHECK(normals != NULL && normals[2] < -0.99f, "winding normal is -x");
    RL_CHECK(RawtexLayers_orient(normals, verts, 4, faces, 2, inward, 0, &os) == 0, "orient runs");
    RL_CHECK(os.components == 1 && os.flipped == 0 && normals[2] < -0.99f, "inward component kept");
    /* ... and flips one facing outward */
    for (i = 0; i < 4; i++) inward[i * 3 + 2] = 1.0f;
    RL_CHECK(RawtexLayers_orient(normals, verts, 4, faces, 2, inward, 2, &os) == 0 &&
             os.flipped == 1 && normals[2] > 0.99f && os.agree_area > 0.999, "outward component flipped");
    /* the decision reads the face winding, so the input must be fresh winding normals */
    free(normals);
    normals = MeshNormals_compute(verts, 4, faces, 2);
    for (i = 0; i < 4; i++) inward[i * 3 + 2] = -1.0f;
    RL_CHECK(normals != NULL && RawtexLayers_orient(normals, verts, 4, faces, 2, inward, 2, &os) == 0 &&
             os.flipped == 0 && normals[2] < -0.99f, "smoothed inward normals stay inward");

    /* 2. render: 5 layers, step 1, offset +0.5 inward: layer k samples x = 7.5 - (k - 2 + 0.5);
     * with the geometry volume on */
    arena = Arena_new();
    memset(&opt, 0, sizeof opt);
    opt.depth = 5; opt.step = 1.0; opt.offset = 0.5; opt.tile_cols = 16; opt.chunk = 8;
    opt.stretch_ratio = 4.0; opt.stretch_floor = 25.0;
    opt.seam_reach = 2; opt.seam_tol = 4.5; opt.seam_passes = 3;
    opt.geom_dir = geom;
    RL_CHECK(cubetable_init(&ct, arena, raw, 8, verts, 4, 5.0) == 0, "cube table");
    cubetable_prewarm_all(&ct);
    RL_CHECK(RawtexLayers_write(zarr, maskp, &ct, verts, uv, 4, faces, 2, normals, NULL, &opt,
                                "{\"test\": 1}\n", &ls) == 0, "render");
    RL_CHECK(ls.W == 10 && ls.H == 10 && ls.painted_px == 100 && ls.complete_px == 100 &&
             ls.chunks_written == 4, "grid, coverage and chunk count");
    opt.geom_dir = NULL;
    RL_CHECK(RawtexLayers_write(zarr, maskp, &ct, verts, uv, 4, faces, 2, normals, NULL, &opt,
                                NULL, &ls2) != 0, "refuses to replace a published volume");
    opt.geom_dir = geom;
    RL_CHECK(RawtexLayers_write(zarr2, mask2, &ct, verts, uv, 4, faces, 2, normals, NULL, &opt,
                                NULL, &ls2) != 0, "refuses to replace a published geometry volume");
    opt.geom_dir = NULL;
    {   /* geometry: point (z,y,x) = (2.5 + py, 2.5 + px, 7.5), normal (0,0,-1); padding NaN */
        float g[6 * 8 * 8];
        char gp[300];
        int px = 0, py = 0, bad = 0, pad_ok = 1;
        size_t cy = 0, cx = 0;
        for (cy = 0; cy < 2; cy++)
            for (cx = 0; cx < 2; cx++) {
                FILE *fp = NULL;
                snprintf(gp, sizeof gp, "%s/0/0.%zu.%zu", geom, cy, cx);
                fp = fopen(gp, "rb");
                if (fp == NULL || fread(g, sizeof(float), 6 * 64, fp) != 6 * 64) {
                    bad++;
                    if (fp != NULL) fclose(fp);
                    continue;
                }
                fclose(fp);
                for (py = 0; py < 8; py++)
                    for (px = 0; px < 8; px++) {
                        int gx = (int)cx * 8 + px, gy = (int)cy * 8 + py, ch = 0;
                        const float want[6] = { 2.5f + (float)gy, 2.5f + (float)gx, 7.5f, 0.0f, 0.0f, -1.0f };
                        for (ch = 0; ch < 6; ch++) {
                            float v = g[(ch * 8 + py) * 8 + px];
                            if (gx >= 10 || gy >= 10) { if (!isnan(v)) pad_ok = 0; }
                            else if (!(fabsf(v - want[ch]) < 1e-5f)) bad++;
                        }
                    }
            }
        RL_CHECK(bad == 0, "geometry = the depth-zero point and unit inward normal of every pixel");
        RL_CHECK(pad_ok, "geometry chunk padding is NaN");
    }
    {
        uint8_t buf[5 * 8 * 8];
        int px = 0, py = 0, k = 0, bad = 0;
        for (py = 0; py < 10; py++)
            for (px = 0; px < 10; px++) {
                if (rl_read_chunk(zarr, (size_t)py / 8, (size_t)px / 8, buf, sizeof buf) != 0) { bad++; continue; }
                for (k = 0; k < 5; k++) {
                    double z = 2.5 + py, y = 2.5 + px, x = 7.5 - ((double)k - 2.0 + 0.5);
                    int want = (int)floor(10.0 + 3.0 * z + 2.0 * y + x + 0.5);
                    if (buf[(k * 8 + py % 8) * 8 + px % 8] != want) bad++;
                }
            }
        RL_CHECK(bad == 0, "layer values = exact trilinear samples, layer index inward");
    }
    {
        Arena_T ma = Arena_new();
        uint8_t *m = NULL;
        int D = 0, Hh = 0, Ww = 0, j = 0, on = 0;
        RL_CHECK(TiffIO_load(ma, maskp, &m, &D, &Hh, &Ww) == 0 && D == 1 && Hh == 10 && Ww == 10,
                 "mask readable at the grid size");
        for (j = 0; m != NULL && j < 100; j++) on += m[j] == 255;
        RL_CHECK(on == 100, "mask complete");
        Arena_dispose(&ma);
    }
    /* 3. tiling does not change a byte: one-chunk tiles vs one tile */
    opt.tile_cols = 8;
    RL_CHECK(RawtexLayers_write(zarr2, mask2, &ct, verts, uv, 4, faces, 2, normals, NULL, &opt,
                                NULL, &ls2) == 0 && ls2.painted_px == 100, "tiled render");
    {
        uint8_t a[5 * 64], b[5 * 64];
        size_t cy = 0, cx = 0;
        int same = 1;
        for (cy = 0; cy < 2; cy++)
            for (cx = 0; cx < 2; cx++)
                if (rl_read_chunk(zarr, cy, cx, a, sizeof a) != 0 ||
                    rl_read_chunk(zarr2, cy, cx, b, sizeof b) != 0 || memcmp(a, b, sizeof a) != 0)
                    same = 0;
        RL_CHECK(same, "tile width is invisible in the output");
    }
    Arena_dispose(&arena);

    /* 4. a missing RAW chunk (z 0-7, y 8-15, x 0-7) never counts as CT: pixels
     * whose samples reach it drop out of the mask, the rest stay */
    arena = Arena_new();
    RL_CHECK(cubetable_init(&ct, arena, raw_gap, 8, verts, 4, 5.0) == 0, "gap cube table");
    cubetable_prewarm_all(&ct);
    opt.tile_cols = 16;
    RL_CHECK(RawtexLayers_write(zarr3, mask3, &ct, verts, uv, 4, faces, 2, normals, NULL, &opt,
                                NULL, &ls2) == 0, "gap render");
    /* layer x samples are 9..5; corners reach the gap (z<=7, y>=8, x<=7) for
     * z corners 2+py <= 7 (py <= 5) and y corners 3+px >= 8 (px >= 5): 6 x 5 */
    RL_CHECK(ls2.painted_px == 100 && ls2.complete_px == 70,
             "missing chunk excluded from the mask only");
    {   /* pixel (5,0) = (z 2.5, y 7.5): layer 2 (x 7) first touches the gap; layers
         * 3 and 4 still sample their own depth (x 6 and 5), renormalized over the
         * y = 7 corners: 10 + 7.5 + 14 + x -> 38 and 37 (depth zero would give 40) */
        uint8_t buf[5 * 8 * 8];
        RL_CHECK(rl_read_chunk(zarr3, 0, 0, buf, sizeof buf) == 0 &&
                 buf[(3 * 8 + 0) * 8 + 5] == 38 && buf[(4 * 8 + 0) * 8 + 5] == 37,
                 "an incomplete pixel's later layers keep their own depth");
    }
    Arena_dispose(&arena);

    /* 6. the geometry volume records the depth axis even where a layer lacks
     * RAW: with the normals flipped to +x, pixels whose z and y corners lie in
     * the missing chunk (py <= 4, px >= 6) drop out of the mask, yet their
     * geometry still holds the unit normal (0,0,+1) and the depth-zero point */
    {
        float pn[12], g[6 * 8 * 8];
        char gp[300];
        int px = 0, py = 0, bad = 0, checked = 0;
        for (i = 0; i < 12; i++) pn[i] = -normals[i];
        arena = Arena_new();
        RL_CHECK(cubetable_init(&ct, arena, raw_gap, 8, verts, 4, 5.0) == 0, "gap cube table (flipped)");
        cubetable_prewarm_all(&ct);
        opt.geom_dir = geom;
        rl_rmtree_zarr(geom, 2, 2);
        RL_CHECK(RawtexLayers_write(zarr6, mask6, &ct, verts, uv, 4, faces, 2, pn, NULL, &opt,
                                    NULL, &ls2) == 0 && ls2.painted_px == 100 && ls2.complete_px < 100,
                 "flipped gap render");
        opt.geom_dir = NULL;
        for (py = 0; py <= 4; py++)
            for (px = 6; px < 8; px++) {
                FILE *fp = NULL;
                snprintf(gp, sizeof gp, "%s/0/0.0.0", geom);
                fp = fopen(gp, "rb");
                if (fp == NULL || fread(g, sizeof(float), 6 * 64, fp) != 6 * 64) {
                    bad++;
                    if (fp != NULL) fclose(fp);
                    continue;
                }
                fclose(fp);
                if (!(fabsf(g[(0 * 8 + py) * 8 + px] - (2.5f + (float)py)) < 1e-5f) ||
                    !(fabsf(g[(1 * 8 + py) * 8 + px] - (2.5f + (float)px)) < 1e-5f) ||
                    !(fabsf(g[(5 * 8 + py) * 8 + px] - 1.0f) < 1e-5f))
                    bad++;
                checked++;
            }
        RL_CHECK(bad == 0 && checked == 10, "geometry keeps the depth axis of masked-out pixels");
        Arena_dispose(&arena);
    }

    /* 5. seam fill: the quad split at u in (7,8) -- the y in (9,10) strip of
     * the plane is missing -- leaves pixel column 7 (centre u = 7.5) unpainted,
     * on the tile boundary of 8-column tiles.  Its anchors (y 8.5 and 10.5,
     * 2 vox apart) fill it with y = 9.5: the plane's own value, at any tile
     * width.  Moving the right part to x = 13 (anchors 5.9 vox apart) keeps
     * the gap open. */
    {
        float sv[24] = { 2.0f, 2.0f, 7.5f,   2.0f, 9.0f, 7.5f,   12.0f, 9.0f, 7.5f,   12.0f, 2.0f, 7.5f,
                         2.0f, 10.0f, 7.5f,  2.0f, 12.0f, 7.5f,  12.0f, 12.0f, 7.5f,  12.0f, 10.0f, 7.5f };
        double suv[16] = { 0.0, 0.0,  7.0, 0.0,  7.0, 10.0,  0.0, 10.0,
                           8.0, 0.0,  10.0, 0.0,  10.0, 10.0,  8.0, 10.0 };
        int32_t sf[12] = { 0, 1, 2,  0, 2, 3,  4, 5, 6,  4, 6, 7 };
        float sin_[24];
        float *sn = NULL;
        RawtexUv suvf = { NULL, suv };
        char za[256], zb[256], zc[256], ma[256], mb[256], mc[256];
        uint8_t a[5 * 64], b[5 * 64];
        size_t cy = 0, cx = 0;
        int same = 1, bad = 0, px = 0, py = 0, k = 0;
        snprintf(za, sizeof za, "%s/seam16.zarr", root); snprintf(ma, sizeof ma, "%s/seam16.tif", root);
        snprintf(zb, sizeof zb, "%s/seam8.zarr", root);  snprintf(mb, sizeof mb, "%s/seam8.tif", root);
        snprintf(zc, sizeof zc, "%s/seamx.zarr", root);  snprintf(mc, sizeof mc, "%s/seamx.tif", root);
        for (i = 0; i < 8; i++) { sin_[i * 3] = 0.0f; sin_[i * 3 + 1] = 0.0f; sin_[i * 3 + 2] = -1.0f; }
        arena = Arena_new();
        RL_CHECK(cubetable_init(&ct, arena, raw, 8, sv, 8, 5.0) == 0, "seam cube table");
        cubetable_prewarm_all(&ct);
        sn = MeshNormals_compute(sv, 8, sf, 4);
        RL_CHECK(sn != NULL && RawtexLayers_orient(sn, sv, 8, sf, 4, sin_, 0, &os) == 0 && os.components == 2,
                 "seam mesh orients as two components");
        opt.tile_cols = 16;
        RL_CHECK(RawtexLayers_write(za, ma, &ct, sv, suvf, 8, sf, 4, sn, NULL, &opt, NULL, &ls) == 0 &&
                 ls.painted_px == 90 && ls.filled_px == 10 && ls.complete_px == 100, "seam column filled");
        opt.tile_cols = 8;
        RL_CHECK(RawtexLayers_write(zb, mb, &ct, sv, suvf, 8, sf, 4, sn, NULL, &opt, NULL, &ls2) == 0 &&
                 ls2.filled_px == 10, "seam fill across a tile boundary");
        for (cy = 0; cy < 2; cy++)
            for (cx = 0; cx < 2; cx++)
                if (rl_read_chunk(za, cy, cx, a, sizeof a) != 0 || rl_read_chunk(zb, cy, cx, b, sizeof b) != 0 ||
                    memcmp(a, b, sizeof a) != 0)
                    same = 0;
        RL_CHECK(same, "seam fill is independent of the tile width");
        for (py = 0; py < 10; py++)
            for (px = 0; px < 10; px++) {
                if (rl_read_chunk(za, (size_t)py / 8, (size_t)px / 8, a, sizeof a) != 0) { bad++; continue; }
                for (k = 0; k < 5; k++) {
                    double z = 2.5 + py, y = 2.5 + px, x = 7.5 - ((double)k - 2.0 + 0.5);
                    if (a[(k * 8 + py % 8) * 8 + px % 8] != (int)floor(10.0 + 3.0 * z + 2.0 * y + x + 0.5)) bad++;
                }
            }
        RL_CHECK(bad == 0, "seam pixels sample the interpolated surface exactly");
        for (i = 4; i < 8; i++) sv[i * 3 + 2] = 13.0f;
        free(sn);
        sn = MeshNormals_compute(sv, 8, sf, 4);
        RL_CHECK(sn != NULL && RawtexLayers_orient(sn, sv, 8, sf, 4, sin_, 0, &os) == 0, "offset seam mesh orients");
        opt.tile_cols = 16;
        RL_CHECK(RawtexLayers_write(zc, mc, &ct, sv, suvf, 8, sf, 4, sn, NULL, &opt, NULL, &ls2) == 0 &&
                 ls2.painted_px == 90 && ls2.filled_px == 0 && ls2.complete_px == 90,
                 "a 3-D discontinuous gap stays open");
        free(sn);
        Arena_dispose(&arena);
        rl_rmtree_zarr(za, 2, 2); rl_rmtree_zarr(zb, 2, 2); rl_rmtree_zarr(zc, 2, 2);
        remove(ma); remove(mb); remove(mc);
    }

    /* 6a. the apron's quadric keeps curvature: a sheet bent on a radius-100
     * cylinder (arc length = pixel x, painted for x < 32), fitted at its border
     * pixel, predicts the unpainted continuation 9 px out to < 0.1 vox and 17 px
     * out to < 0.3 vox (the cubic term; a tangent step would be 0.4 / 1.4 vox off) */
    {
        const long QW = 64, QH = 64;
        uint8_t *qs = (uint8_t *)calloc((size_t)(QW * QH), 1);
        float *qp = (float *)malloc((size_t)(QW * QH) * 3 * sizeof(float));
        double coef[18], err[2] = { 0.0, 0.0 };
        long qx = 0, qy = 0;
        int k = 0, j = 0;
        const long probes[2] = { 40, 48 };
        for (qy = 0; qy < QH; qy++)
            for (qx = 0; qx < QW; qx++) {
                size_t q = (size_t)(qy * QW + qx);
                double th = ((double)qx + 0.5) / 100.0;
                qs[q] = qx < 32 ? 2 : 0;
                qp[q * 3] = (float)qy; qp[q * 3 + 1] = (float)(100.0 * sin(th)); qp[q * 3 + 2] = (float)(100.0 * cos(th));
            }
        RL_CHECK(RawtexLayers_quadric_fit(qs, qp, QW, QH, 31, 32, 16, 1, coef) == 0, "quadric fit on a curved sheet");
        for (k = 0; k < 2; k++) {
            double u = (double)(probes[k] - 31) / 16.0, v = 0.0, phi[6], e = 0.0;
            double th = ((double)probes[k] + 0.5) / 100.0, want[3];
            want[0] = 32.0; want[1] = 100.0 * sin(th); want[2] = 100.0 * cos(th);
            phi[0] = 1.0; phi[1] = u; phi[2] = v; phi[3] = u * u; phi[4] = u * v; phi[5] = v * v;
            for (j = 0; j < 3; j++) {
                double p = 0.0;
                int b = 0;
                for (b = 0; b < 6; b++) p += coef[j * 6 + b] * phi[b];
                e += (p - want[j]) * (p - want[j]);
            }
            err[k] = sqrt(e);
        }
        if (!(err[0] < 0.1 && err[1] < 0.3))
            fprintf(stderr, "  quadric continuation error %.3f vox at 9 px, %.3f at 17 px\n", err[0], err[1]);
        RL_CHECK(err[0] < 0.1 && err[1] < 0.3, "quadric continuation stays on a curved sheet 9 and 17 px out");
        free(qs); free(qp);
    }

    /* 6b. context apron: the plane split at u in (4, 8) leaves a 4 px hole (too
     * wide for the seam fill).  With the apron, the hole's 40 pixels are rendered
     * on the plane's own continuation (layer values = the linear field, within
     * rounding), are NOT in the mask, and do not depend on the tile width;
     * without it they stay zero. */
    {
        float av[24] = { 2.0f, 2.0f, 7.5f,   2.0f, 6.0f, 7.5f,    12.0f, 6.0f, 7.5f,    12.0f, 2.0f, 7.5f,
                         2.0f, 10.0f, 7.5f,  2.0f, 12.0f, 7.5f,  12.0f, 12.0f, 7.5f,  12.0f, 10.0f, 7.5f };
        double auv[16] = { 0.0, 0.0,  4.0, 0.0,  4.0, 10.0,  0.0, 10.0,
                           8.0, 0.0,  10.0, 0.0,  10.0, 10.0,  8.0, 10.0 };
        int32_t af[12] = { 0, 1, 2,  0, 2, 3,  4, 5, 6,  4, 6, 7 };
        float ain[24];
        float *an = NULL;
        RawtexUv auvf = { NULL, auv };
        char zo[256], zp[256], zq[256], mo[256], mp[256], mq[256];
        uint8_t a[5 * 64], b[5 * 64];
        size_t cy = 0, cx = 0;
        int same = 1, bad = 0, off = 0, px = 0, py = 0, k = 0;
        RawtexLayersStats la, lb, lc;
        snprintf(zo, sizeof zo, "%s/apron_off.zarr", root); snprintf(mo, sizeof mo, "%s/apron_off.tif", root);
        snprintf(zp, sizeof zp, "%s/apron16.zarr", root);   snprintf(mp, sizeof mp, "%s/apron16.tif", root);
        snprintf(zq, sizeof zq, "%s/apron8.zarr", root);    snprintf(mq, sizeof mq, "%s/apron8.tif", root);
        for (i = 0; i < 8; i++) { ain[i * 3] = 0.0f; ain[i * 3 + 1] = 0.0f; ain[i * 3 + 2] = -1.0f; }
        arena = Arena_new();
        RL_CHECK(cubetable_init(&ct, arena, raw, 8, av, 8, 5.0) == 0, "apron cube table");
        cubetable_prewarm_all(&ct);
        an = MeshNormals_compute(av, 8, af, 4);
        RL_CHECK(an != NULL && RawtexLayers_orient(an, av, 8, af, 4, ain, 0, &os) == 0, "apron mesh orients");
        opt.tile_cols = 16; opt.apron_px = 0;
        RL_CHECK(RawtexLayers_write(zo, mo, &ct, av, auvf, 8, af, 4, an, NULL, &opt, NULL, &la) == 0 &&
                 la.painted_px == 60 && la.filled_px == 0 && la.apron_px == 0 && la.complete_px == 60,
                 "a 4 px hole stays open without the apron");
        opt.apron_px = 6; opt.apron_fit_px = 6; opt.apron_cell_px = 2; opt.apron_step_px = 1;
        RL_CHECK(RawtexLayers_write(zp, mp, &ct, av, auvf, 8, af, 4, an, NULL, &opt, NULL, &lb) == 0 &&
                 lb.painted_px == 60 && lb.apron_px == 40 && lb.complete_px == 60,
                 "the apron renders the hole as context only (not in the mask)");
        opt.tile_cols = 8;
        RL_CHECK(RawtexLayers_write(zq, mq, &ct, av, auvf, 8, af, 4, an, NULL, &opt, NULL, &lc) == 0 &&
                 lc.apron_px == 40, "apron across a tile boundary");
        for (cy = 0; cy < 2; cy++)
            for (cx = 0; cx < 2; cx++)
                if (rl_read_chunk(zp, cy, cx, a, sizeof a) != 0 || rl_read_chunk(zq, cy, cx, b, sizeof b) != 0 ||
                    memcmp(a, b, sizeof a) != 0)
                    same = 0;
        RL_CHECK(same, "the apron is independent of the tile width");
        for (py = 0; py < 10; py++)
            for (px = 4; px < 8; px++) {
                if (rl_read_chunk(zp, (size_t)py / 8, (size_t)px / 8, a, sizeof a) != 0 ||
                    rl_read_chunk(zo, (size_t)py / 8, (size_t)px / 8, b, sizeof b) != 0) { bad++; continue; }
                for (k = 0; k < 5; k++) {
                    double z = 2.5 + py, y = 2.5 + px, x = 7.5 - ((double)k - 2.0 + 0.5);
                    int want = (int)floor(10.0 + 3.0 * z + 2.0 * y + x + 0.5);
                    int got = a[(k * 8 + py % 8) * 8 + px % 8];
                    if (got < want - 1 || got > want + 1) bad++;
                    if (b[(k * 8 + py % 8) * 8 + px % 8] != 0) off++;
                }
            }
        RL_CHECK(bad == 0 && off == 0, "apron layers sample the plane's continuation; without it they are zero");
        {
            Arena_T ma = Arena_new();
            uint8_t *m = NULL, *m0 = NULL;
            int D = 0, Hh = 0, Ww = 0, D0 = 0, H0 = 0, W0 = 0;
            RL_CHECK(TiffIO_load(ma, mp, &m, &D, &Hh, &Ww) == 0 && TiffIO_load(ma, mo, &m0, &D0, &H0, &W0) == 0 &&
                     Hh == H0 && Ww == W0 && memcmp(m, m0, (size_t)Hh * (size_t)Ww) == 0,
                     "the apron leaves the mask byte-identical");
            Arena_dispose(&ma);
        }
        opt.apron_px = 0; opt.tile_cols = 16;
        free(an);
        Arena_dispose(&arena);
        rl_rmtree_zarr(zo, 2, 2); rl_rmtree_zarr(zp, 2, 2); rl_rmtree_zarr(zq, 2, 2);
        remove(mo); remove(mp); remove(mq);
    }

    /* 7a. regularisation flattens a wobble: the plane x = 7.5 rippled by
     * +-0.5 vox (periods 6 and 7 px), with the ripple's own normals, comes out
     * flat (sigma 4): < 0.02 vox and < 1 deg at >= 12 px from the tile border,
     * < 0.25 vox and < 10 deg everywhere; nothing moves in-plane */
    {
        const size_t SW = 48, SH = 48;
        uint8_t *ss = (uint8_t *)malloc(SW * SH);
        float *sp = (float *)malloc(SW * SH * 3 * sizeof(float));
        float *sn = (float *)malloc(SW * SH * 3 * sizeof(float));
        size_t sx = 0, sy = 0, np = 0, nc = 0;
        double s2 = 0.0, tu = 0.0, in_dx = 0.0, in_da = 0.0, all_dx = 0.0, all_da = 0.0, tang = 0.0;
        const double TWO_PI = 6.28318530717958647692;
        for (sy = 0; sy < SH; sy++)
            for (sx = 0; sx < SW; sx++) {
                size_t q = sy * SW + sx;
                double a = 2.0 * TWO_PI / 12.0 * (double)sx, b = TWO_PI / 7.0 * (double)sy;
                double f = 0.5 * sin(a) * sin(b);
                double fy = 0.5 * (TWO_PI / 6.0) * cos(a) * sin(b), fz = 0.5 * (TWO_PI / 7.0) * sin(a) * cos(b);
                double nl = sqrt(fz * fz + fy * fy + 1.0);
                ss[q] = 2;
                sp[q * 3] = (float)sy; sp[q * 3 + 1] = (float)sx; sp[q * 3 + 2] = (float)(7.5 + f);
                sn[q * 3] = (float)(fz / nl); sn[q * 3 + 1] = (float)(fy / nl); sn[q * 3 + 2] = (float)(-1.0 / nl);
            }
        RL_CHECK(RawtexLayers_smooth(ss, sp, sn, SW, SH, 4.0, 2.0, 0, SW, &np, &nc, &s2, &tu) == 0 &&
                 np == SW * SH && nc == 0, "regularisation runs on a rippled plane");
        for (sy = 0; sy < SH; sy++)
            for (sx = 0; sx < SW; sx++) {
                size_t q = sy * SW + sx;
                double dx = fabs((double)sp[q * 3 + 2] - 7.5);
                double da = acos(fmin(1.0, fabs((double)sn[q * 3 + 2]))) * 180.0 / 3.14159265358979323846;
                int inner = sx >= 12 && sx < SW - 12 && sy >= 12 && sy < SH - 12;
                tang = fmax(tang, fabs((double)sp[q * 3] - (double)sy) + fabs((double)sp[q * 3 + 1] - (double)sx));
                all_dx = fmax(all_dx, dx); all_da = fmax(all_da, da);
                if (inner) { in_dx = fmax(in_dx, dx); in_da = fmax(in_da, da); }
            }
        if (!(in_dx < 0.02 && in_da < 1.0 && all_dx < 0.25 && all_da < 10.0))
            fprintf(stderr, "  ripple after smoothing: inner %.4f vox %.3f deg, all %.4f vox %.3f deg\n",
                    in_dx, in_da, all_dx, all_da);
        RL_CHECK(in_dx < 0.02 && in_da < 1.0 && all_dx < 0.25 && all_da < 10.0,
                 "a +-0.5 vox ripple is regularised flat, normals level");
        RL_CHECK(tang < 0.05, "regularisation moves points along the normal only (ripple normals are near x)");
        free(ss); free(sp); free(sn);
    }

    /* 7b. regularisation keeps real curvature: a radius-100 cylinder stays within
     * 0.12 vox of its radius (sigma 4 shrinks it by ~sigma^2 / 2R = 0.08) and its
     * normals radial within 0.5 deg, away from the tile border */
    {
        const size_t CW = 64, CH = 32;
        uint8_t *cs = (uint8_t *)malloc(CW * CH);
        float *cp = (float *)malloc(CW * CH * 3 * sizeof(float));
        float *cn = (float *)malloc(CW * CH * 3 * sizeof(float));
        size_t cx = 0, cy = 0;
        double worst_r = 0.0, worst_a = 0.0;
        for (cy = 0; cy < CH; cy++)
            for (cx = 0; cx < CW; cx++) {
                size_t q = cy * CW + cx;
                double th = ((double)cx + 0.5) / 100.0;
                cs[q] = 2;
                cp[q * 3] = (float)cy; cp[q * 3 + 1] = (float)(100.0 * sin(th)); cp[q * 3 + 2] = (float)(100.0 * cos(th));
                cn[q * 3] = 0.0f; cn[q * 3 + 1] = (float)(-sin(th)); cn[q * 3 + 2] = (float)(-cos(th));
            }
        RL_CHECK(RawtexLayers_smooth(cs, cp, cn, CW, CH, 4.0, 2.0, 0, CW, NULL, NULL, NULL, NULL) == 0,
                 "regularisation runs on a cylinder");
        for (cy = 12; cy < CH - 12; cy++)
            for (cx = 12; cx < CW - 12; cx++) {
                size_t q = cy * CW + cx;
                double y = cp[q * 3 + 1], x = cp[q * 3 + 2], r = sqrt(y * y + x * x);
                double c = -(cn[q * 3 + 1] * y + cn[q * 3 + 2] * x) / r;
                worst_r = fmax(worst_r, fabs(r - 100.0));
                worst_a = fmax(worst_a, acos(fmin(1.0, c)) * 180.0 / 3.14159265358979323846);
            }
        if (!(worst_r < 0.12 && worst_a < 0.5))
            fprintf(stderr, "  cylinder after smoothing: radius off %.4f vox, normal %.3f deg\n", worst_r, worst_a);
        RL_CHECK(worst_r < 0.12 && worst_a < 0.5, "regularisation keeps a radius-100 turn");
        free(cs); free(cp); free(cn);
    }

    /* 7c. the cap bounds every move: a 5 vox step between two halves of a tile
     * (sigma 2, cap 1) clamps the pixels next to the step, and none moves more */
    {
        const size_t KW = 32, KH = 8;
        uint8_t ks[32 * 8];
        float kp[32 * 8 * 3], kn[32 * 8 * 3], k0[32 * 8 * 3];
        size_t kx = 0, ky = 0, np = 0, nc = 0;
        double worst = 0.0;
        for (ky = 0; ky < KH; ky++)
            for (kx = 0; kx < KW; kx++) {
                size_t q = ky * KW + kx;
                ks[q] = 2;
                kp[q * 3] = (float)ky; kp[q * 3 + 1] = (float)kx; kp[q * 3 + 2] = kx < KW / 2 ? 7.5f : 12.5f;
                kn[q * 3] = 0.0f; kn[q * 3 + 1] = 0.0f; kn[q * 3 + 2] = -1.0f;
            }
        memcpy(k0, kp, sizeof kp);
        RL_CHECK(RawtexLayers_smooth(ks, kp, kn, KW, KH, 2.0, 1.0, 0, KW, &np, &nc, NULL, NULL) == 0 && nc > 0,
                 "a 3-D step clamps");
        for (kx = 0; kx < KW * KH * 3; kx++) worst = fmax(worst, fabs((double)kp[kx] - (double)k0[kx]));
        RL_CHECK(worst <= 1.0 + 1e-5, "no pixel moves more than the cap");
    }

    /* 7d. with regularisation (and the apron) on, the render still does not
     * depend on the tile width, and a flat plane's layers keep the linear field */
    {
        float av[24] = { 2.0f, 2.0f, 7.5f,   2.0f, 6.0f, 7.5f,    12.0f, 6.0f, 7.5f,    12.0f, 2.0f, 7.5f,
                         2.0f, 10.0f, 7.5f,  2.0f, 12.0f, 7.5f,  12.0f, 12.0f, 7.5f,  12.0f, 10.0f, 7.5f };
        double auv[16] = { 0.0, 0.0,  4.0, 0.0,  4.0, 10.0,  0.0, 10.0,
                           8.0, 0.0,  10.0, 0.0,  10.0, 10.0,  8.0, 10.0 };
        int32_t af[12] = { 0, 1, 2,  0, 2, 3,  4, 5, 6,  4, 6, 7 };
        float ain[24];
        float *an = NULL;
        RawtexUv auvf = { NULL, auv };
        char zp[256], zq[256], mp[256], mq[256];
        uint8_t a[5 * 64], b[5 * 64];
        size_t cy = 0, cx = 0;
        int same = 1, bad = 0, px = 0, py = 0, k = 0;
        RawtexLayersStats la, lb;
        snprintf(zp, sizeof zp, "%s/smooth16.zarr", root); snprintf(mp, sizeof mp, "%s/smooth16.tif", root);
        snprintf(zq, sizeof zq, "%s/smooth8.zarr", root);  snprintf(mq, sizeof mq, "%s/smooth8.tif", root);
        rl_rmtree_zarr(zp, 2, 2); rl_rmtree_zarr(zq, 2, 2);
        for (i = 0; i < 8; i++) { ain[i * 3] = 0.0f; ain[i * 3 + 1] = 0.0f; ain[i * 3 + 2] = -1.0f; }
        arena = Arena_new();
        RL_CHECK(cubetable_init(&ct, arena, raw, 8, av, 8, 5.0) == 0, "smooth cube table");
        cubetable_prewarm_all(&ct);
        an = MeshNormals_compute(av, 8, af, 4);
        RL_CHECK(an != NULL && RawtexLayers_orient(an, av, 8, af, 4, ain, 0, &os) == 0, "smooth mesh orients");
        opt.apron_px = 6; opt.apron_fit_px = 6; opt.apron_cell_px = 2; opt.apron_step_px = 1;
        opt.smooth_sigma_px = 2.0; opt.smooth_cap_vox = 1.5;
        opt.tile_cols = 16;
        RL_CHECK(RawtexLayers_write(zp, mp, &ct, av, auvf, 8, af, 4, an, NULL, &opt, NULL, &la) == 0 &&
                 la.painted_px == 60 && la.smooth_px == 60 && la.smooth_capped_px == 0 && la.complete_px == 60,
                 "render with regularisation");
        opt.tile_cols = 8;
        RL_CHECK(RawtexLayers_write(zq, mq, &ct, av, auvf, 8, af, 4, an, NULL, &opt, NULL, &lb) == 0 &&
                 lb.smooth_px == 60, "regularised render across a tile boundary");
        for (cy = 0; cy < 2; cy++)
            for (cx = 0; cx < 2; cx++)
                if (rl_read_chunk(zp, cy, cx, a, sizeof a) != 0 || rl_read_chunk(zq, cy, cx, b, sizeof b) != 0 ||
                    memcmp(a, b, sizeof a) != 0)
                    same = 0;
        RL_CHECK(same, "regularisation is independent of the tile width");
        for (py = 0; py < 10; py++)
            for (px = 0; px < 10; px++) {
                if (rl_read_chunk(zp, (size_t)py / 8, (size_t)px / 8, a, sizeof a) != 0) { bad++; continue; }
                for (k = 0; k < 5; k++) {
                    double z = 2.5 + py, y = 2.5 + px, x = 7.5 - ((double)k - 2.0 + 0.5);
                    int want = (int)floor(10.0 + 3.0 * z + 2.0 * y + x + 0.5);
                    int got = a[(k * 8 + py % 8) * 8 + px % 8];
                    if (got < want - 1 || got > want + 1) bad++;
                }
            }
        RL_CHECK(bad == 0, "a flat plane keeps its layers under regularisation");
        RL_CHECK(la.smooth_shift_rms < 1e-3 && la.smooth_turn_mean < 0.01, "a flat plane does not move");
        opt.apron_px = 0; opt.smooth_sigma_px = 0.0; opt.smooth_cap_vox = 0.0; opt.tile_cols = 16;
        free(an);
        Arena_dispose(&arena);
        rl_rmtree_zarr(zp, 2, 2); rl_rmtree_zarr(zq, 2, 2);
        remove(mp); remove(mq);
    }

    /* 7e. step blending bends a tear: a plane whose right half is misregistered
     * 3 vox along y (adjacent pixels 4 vox apart across column 23|24) comes out
     * with no neighbour step above 1.6 vox; 20+ px from the tear (and never on
     * a continuous plane) nothing moves more than 0.1 vox */
    {
        const size_t TW = 48, TH = 16;
        uint8_t *ts = (uint8_t *)malloc(TW * TH);
        float *tp = (float *)malloc(TW * TH * 3 * sizeof(float));
        float *tn = (float *)malloc(TW * TH * 3 * sizeof(float));
        float *t0 = (float *)malloc(TW * TH * 3 * sizeof(float));
        size_t tx = 0, ty = 0, nm = 0;
        double worst_jump = 0.0, far_move = 0.0, m2 = 0.0;
        int pass = 0;
        for (pass = 0; pass < 2; pass++) {
            for (ty = 0; ty < TH; ty++)
                for (tx = 0; tx < TW; tx++) {
                    size_t q = ty * TW + tx;
                    ts[q] = 2;
                    tp[q * 3] = (float)ty; tp[q * 3 + 1] = (float)tx + (pass == 0 && tx >= 24 ? 3.0f : 0.0f);
                    tp[q * 3 + 2] = 7.5f;
                    tn[q * 3] = 0.0f; tn[q * 3 + 1] = 0.0f; tn[q * 3 + 2] = -1.0f;
                }
            memcpy(t0, tp, TW * TH * 3 * sizeof(float));
            nm = 0; m2 = 0.0;
            RL_CHECK(RawtexLayers_stepblend(ts, tp, tn, TW, TH, 4.0, 8.0, 0, TW, &nm, &m2) == 0, "step blending runs");
            if (pass == 1) {
                double moved = 0.0;
                for (tx = 0; tx < TW * TH * 3; tx++) moved = fmax(moved, fabs((double)tp[tx] - (double)t0[tx]));
                RL_CHECK(nm == 0 && moved == 0.0, "a continuous plane does not move");
                break;
            }
            worst_jump = 0.0; far_move = 0.0;
            for (ty = 0; ty < TH; ty++)
                for (tx = 0; tx < TW; tx++) {
                    size_t q = ty * TW + tx;
                    double mv = 0.0;
                    int k = 0;
                    if (tx + 1 < TW) {
                        double d2 = 0.0;
                        for (k = 0; k < 3; k++) {
                            double e = (double)tp[(q + 1) * 3 + (size_t)k] - (double)tp[q * 3 + (size_t)k];
                            d2 += e * e;
                        }
                        worst_jump = fmax(worst_jump, sqrt(d2));
                    }
                    for (k = 0; k < 3; k++) mv = fmax(mv, fabs((double)tp[q * 3 + (size_t)k] - (double)t0[q * 3 + (size_t)k]));
                    if (tx < 4 || tx >= TW - 4) far_move = fmax(far_move, mv);
                }
            if (!(worst_jump < 1.6 && far_move < 0.1))
                fprintf(stderr, "  step blending: worst neighbour step %.3f vox, far move %.3f vox\n", worst_jump, far_move);
            RL_CHECK(nm > 0 && worst_jump < 1.6, "a 3 vox tear is bent into a continuous surface");
            RL_CHECK(far_move < 0.1, "20+ px from the tear nothing moves");
        }
        free(ts); free(tp); free(tn); free(t0);
    }

    /* 7f. step blending in the render: two quads meeting at u = 5 with a 2 vox
     * tangential misregistration (right quad's y starts at 9, not 7) -- the
     * blend moves pixels and the output does not depend on the tile width */
    {
        float sv[24] = { 2.0f, 2.0f, 7.5f,   2.0f, 7.0f, 7.5f,   12.0f, 7.0f, 7.5f,   12.0f, 2.0f, 7.5f,
                         2.0f, 9.0f, 7.5f,   2.0f, 14.0f, 7.5f,  12.0f, 14.0f, 7.5f,  12.0f, 9.0f, 7.5f };
        double suv[16] = { 0.0, 0.0,  5.0, 0.0,  5.0, 10.0,  0.0, 10.0,
                           5.0, 0.0,  10.0, 0.0,  10.0, 10.0,  5.0, 10.0 };
        int32_t sf[12] = { 0, 1, 2,  0, 2, 3,  4, 5, 6,  4, 6, 7 };
        float sin_[24];
        float *sn = NULL;
        RawtexUv suvf = { NULL, suv };
        char za[256], zb[256], ma[256], mb[256];
        uint8_t a[5 * 64], b[5 * 64];
        size_t cy = 0, cx = 0;
        int same = 1;
        RawtexLayersStats la, lb;
        snprintf(za, sizeof za, "%s/step16.zarr", root); snprintf(ma, sizeof ma, "%s/step16.tif", root);
        snprintf(zb, sizeof zb, "%s/step8.zarr", root);  snprintf(mb, sizeof mb, "%s/step8.tif", root);
        rl_rmtree_zarr(za, 2, 2); rl_rmtree_zarr(zb, 2, 2);
        for (i = 0; i < 8; i++) { sin_[i * 3] = 0.0f; sin_[i * 3 + 1] = 0.0f; sin_[i * 3 + 2] = -1.0f; }
        arena = Arena_new();
        RL_CHECK(cubetable_init(&ct, arena, raw, 8, sv, 8, 5.0) == 0, "step cube table");
        cubetable_prewarm_all(&ct);
        sn = MeshNormals_compute(sv, 8, sf, 4);
        RL_CHECK(sn != NULL && RawtexLayers_orient(sn, sv, 8, sf, 4, sin_, 0, &os) == 0, "step mesh orients");
        opt.step_sigma_px = 2.0; opt.step_spread_px = 3.0; opt.tile_cols = 16;
        RL_CHECK(RawtexLayers_write(za, ma, &ct, sv, suvf, 8, sf, 4, sn, NULL, &opt, NULL, &la) == 0 &&
                 la.painted_px == 100 && la.step_px > 0 && la.step_move_rms > 0.05, "render with step blending");
        opt.tile_cols = 8;
        RL_CHECK(RawtexLayers_write(zb, mb, &ct, sv, suvf, 8, sf, 4, sn, NULL, &opt, NULL, &lb) == 0 &&
                 lb.step_px == la.step_px, "step blending across a tile boundary");
        for (cy = 0; cy < 2; cy++)
            for (cx = 0; cx < 2; cx++)
                if (rl_read_chunk(za, cy, cx, a, sizeof a) != 0 || rl_read_chunk(zb, cy, cx, b, sizeof b) != 0 ||
                    memcmp(a, b, sizeof a) != 0)
                    same = 0;
        RL_CHECK(same, "step blending is independent of the tile width");
        opt.step_sigma_px = 0.0; opt.step_spread_px = 0.0; opt.tile_cols = 16;
        free(sn);
        Arena_dispose(&arena);
        rl_rmtree_zarr(za, 2, 2); rl_rmtree_zarr(zb, 2, 2);
        remove(ma); remove(mb);
    }

    /* 8. handedness census: the plane x = 7.5 with uv = (y - 2, z - 2) and its
     * normal inward (-x) has dP/dcol = +y, dP/drow = +z; in (z, y, x) order
     * (0,1,0) x (1,0,0) = (0,0,-1), and . (0,0,-1) = +1: the team's reading
     * handedness on all 100 px, one piece.  The same plane with the flattening
     * flipped top to bottom (v = 12 - z) is mirror-imaged; so is the original
     * seen with its normals turned outward. */
    {
        double fuv[8] = { 0.0, 10.0,  10.0, 10.0,  10.0, 0.0,  0.0, 0.0 };
        RawtexUv fuvf = { NULL, fuv };
        RawtexHandedness hd;
        char zh[256], mh[256];
        float pn[12];
        int run = 0;
        memset(&hd, 0, sizeof hd);
        snprintf(zh, sizeof zh, "%s/hand.zarr", root); snprintf(mh, sizeof mh, "%s/hand.tif", root);
        for (i = 0; i < 12; i++) pn[i] = -normals[i];
        arena = Arena_new();
        RL_CHECK(cubetable_init(&ct, arena, raw, 8, verts, 4, 5.0) == 0, "handedness cube table");
        cubetable_prewarm_all(&ct);
        opt.hand = &hd; opt.tile_cols = 8;             /* two tiles: the census crosses a tile boundary */
        for (run = 0; run < 3; run++) {
            const float *nn = run == 2 ? pn : normals;
            RawtexUv u = run == 1 ? fuvf : uv;
            rl_rmtree_zarr(zh, 2, 2);
            remove(mh);
            RL_CHECK(RawtexLayers_write(zh, mh, &ct, verts, u, 4, faces, 2, nn, NULL, &opt, NULL, &ls2) == 0 &&
                     hd.npieces == 1, "handedness census runs");
            if (hd.npieces == 1) {
                const RawtexPieceHand *p = &hd.pieces[0];
                RL_CHECK(p->col0 == 0 && p->col1 == 10 && p->row0 == 0 && p->row1 == 10, "one piece, its extent");
                RL_CHECK(run == 0 ? (p->px_team == 100 && p->px_mirror == 0) : (p->px_team == 0 && p->px_mirror == 100),
                         run == 0 ? "the plane shows the recto in the team's handedness"
                                  : (run == 1 ? "a top-to-bottom flipped flattening is mirror-imaged"
                                              : "normals turned outward read as mirror-imaged"));
            }
        }
        RawtexHandedness_free(&hd);
        opt.hand = NULL; opt.tile_cols = 16;
        Arena_dispose(&arena);
        rl_rmtree_zarr(zh, 2, 2);
        remove(mh);
    }

    /* 9. an absolute window renders exactly those columns of the full grid:
     * u in [3, 8) of the 10 x 10 plane gives a 5-wide raster whose layers and
     * mask equal the full render's columns 3..7 */
    {
        RawtexWindow w = { 3.0, 8.0, 0.0, 10.0 };
        char zf[256], zw[256], mf[256], mw[256];
        uint8_t a[5 * 64], b[5 * 64];
        int bad = 0, px = 0, py = 0, k = 0;
        RawtexLayersStats lf, lw;
        snprintf(zf, sizeof zf, "%s/winfull.zarr", root); snprintf(mf, sizeof mf, "%s/winfull.tif", root);
        snprintf(zw, sizeof zw, "%s/win.zarr", root);     snprintf(mw, sizeof mw, "%s/win.tif", root);
        rl_rmtree_zarr(zf, 2, 2); rl_rmtree_zarr(zw, 2, 2);
        arena = Arena_new();
        RL_CHECK(cubetable_init(&ct, arena, raw, 8, verts, 4, 5.0) == 0, "window cube table");
        cubetable_prewarm_all(&ct);
        opt.tile_cols = 8;
        RL_CHECK(RawtexLayers_write(zf, mf, &ct, verts, uv, 4, faces, 2, normals, NULL, &opt, NULL, &lf) == 0,
                 "full render for the window test");
        opt.window = &w;
        RL_CHECK(RawtexLayers_write(zw, mw, &ct, verts, uv, 4, faces, 2, normals, NULL, &opt, NULL, &lw) == 0 &&
                 lw.W == 5 && lw.H == 10 && lw.painted_px == 50 && fabs(lw.umin - 3.0) < 1e-12,
                 "window render: 5 columns, origin u = 3");
        opt.window = NULL; opt.tile_cols = 16;
        for (py = 0; py < 10; py++)
            for (px = 0; px < 5; px++) {
                if (rl_read_chunk(zf, (size_t)py / 8, (size_t)(px + 3) / 8, a, sizeof a) != 0 ||
                    rl_read_chunk(zw, (size_t)py / 8, (size_t)px / 8, b, sizeof b) != 0) { bad++; continue; }
                for (k = 0; k < 5; k++)
                    if (a[(k * 8 + py % 8) * 8 + (px + 3) % 8] != b[(k * 8 + py % 8) * 8 + px % 8]) bad++;
            }
        RL_CHECK(bad == 0, "window layers equal the full render's columns");
        {
            Arena_T ma = Arena_new();
            uint8_t *m0 = NULL, *m1 = NULL;
            int D = 0, H0 = 0, W0 = 0, D1 = 0, H1 = 0, W1 = 0, same = 1, r = 0;
            RL_CHECK(TiffIO_load(ma, mf, &m0, &D, &H0, &W0) == 0 && TiffIO_load(ma, mw, &m1, &D1, &H1, &W1) == 0 &&
                     H1 == H0 && W1 == 5, "window mask readable");
            for (r = 0; m0 != NULL && m1 != NULL && r < H0 && same; r++)
                same = memcmp(m0 + (size_t)r * (size_t)W0 + 3, m1 + (size_t)r * 5, 5) == 0;
            RL_CHECK(same, "window mask equals the full mask's columns");
            Arena_dispose(&ma);
        }
        Arena_dispose(&arena);
        rl_rmtree_zarr(zf, 2, 2); rl_rmtree_zarr(zw, 2, 2);
        remove(mf); remove(mw);
    }

    free(normals);
    rl_rmtree_zarr(zarr, 2, 2);
    rl_rmtree_zarr(zarr2, 2, 2);
    rl_rmtree_zarr(zarr3, 2, 2);
    rl_rmtree_zarr(zarr6, 2, 2);
    rl_rmtree_zarr(geom, 2, 2);
    remove(maskp); remove(mask2); remove(mask3); remove(mask6);
    rl_rm_raw(raw);
    rl_rm_raw(raw_gap);
    ves_rmdir(root);
    fprintf(stderr, "[selftest] rawtex_layers %s (%d failures)\n", fails ? "FAIL" : "PASS", fails);
    return fails ? 1 : 0;
}
