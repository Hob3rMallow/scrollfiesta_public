/* sheet_probe -- per-cube, per-SHEET local-frame extension study.
 *
 * For one cube of a placed grid: split its mesh into sheets (connected
 * components), and for each sheet test "the path" in the sheet's OWN frame:
 *   1. PCA tangent frame (e1,e2 in-plane, n normal) at the sheet centroid;
 *      heightfield h(a,b) on a small lattice (cell vox per cell).
 *   2. Extend h into a --band cell margin: quadratic LSQ predictor (captures
 *      the local curvature, so the extension continues the osculating
 *      cylinder, not the tangent plane) + harmonic residual (Jacobi CG,
 *      sheet cells pinned).
 *   3. Occupancy guard: extension nodes within --guard vox of ANOTHER
 *      sheet's verts are flagged (red) and excluded from the snap.
 *   4. RAW snap: +-n probes for bright papyrus (free-space march blocked by
 *      other sheets), displacement capped at --max-disp.
 * Outputs per sheet: local-frame pre/post RAW bakes (PNG: sheet gray,
 * extension green, guard red), a world-space OBJ of sheet+extension, and a
 * stats line; plus a per-cube vertical montage of all post bakes.
 *
 * This is a STUDY instrument: it shows what parameter extension does when
 * the frame is trustworthy (one sheet, local PCA) instead of inheriting the
 * whole-strip registration. No global u/phi anywhere.
 *
 * Usage:
 *   sheet_probe <placed_dir> <cube_id> --raw <cubes_RAW> [--out DIR]
 *               [--cell F=2.0] [--band N=20] [--min-verts N=200]
 *               [--guard F=3.0] [--reach F=8] [--max-disp F=4]
 *               [--raw-chunk N=128] [--no-snap]
 *   sheet_probe --selftest
 */
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/arena.h"
#include "../common/kdtree.h"
#include "../common/obj_io.h"
#include "../common/pca.h"
#include "../common/raw_sample.h"
#include "../common/tiff_io.h"
#include "../common/union_find.h"
#include "../common/ves_platform.h"
#include "../common/ves_png.h"

#define SP_MAX_GRID (2048u * 2048u)

typedef struct {
    double cell;         /* lattice pitch, vox (2.0) */
    int band;            /* extension margin, cells (20) */
    size_t min_verts;    /* skip smaller sheets (200) */
    double guard_r;      /* foreign clearance, vox (3.0) */
    double reach, step;  /* snap ray (8.0 / 0.25) */
    double occ_thresh, local_r;      /* ray blocking (0.9 / 2.0) */
    double band_frac, min_gain;      /* accept thresholds (0.30 / 25) */
    double max_disp;     /* snap cap, vox (4.0) */
    long raw_chunk;      /* RAW cube edge (128) */
    int do_snap;
} SpOpts;

static void sp_opts_default(SpOpts *o)
{
    memset(o, 0, sizeof(*o));
    o->cell = 2.0;
    o->band = 20;
    o->min_verts = 200;
    o->guard_r = 3.0;
    o->reach = 8.0;
    o->step = 0.25;
    o->occ_thresh = 0.9;
    o->local_r = 2.0;
    o->band_frac = 0.30;
    o->min_gain = 25.0;
    o->max_disp = 4.0;
    o->raw_chunk = 128;
    o->do_snap = 1;
}

static void *sp_xmalloc(size_t n)
{
    void *p = malloc(n > 0 ? n : 1);
    if (p == NULL) { fprintf(stderr, "sheet_probe: OOM %zu\n", n); exit(1); }
    return p;
}

static void *sp_xcalloc(size_t c, size_t s)
{
    void *p = calloc(c > 0 ? c : 1, s);
    if (p == NULL) { fprintf(stderr, "sheet_probe: OOM\n"); exit(1); }
    return p;
}

/* 6x6 Gaussian elimination with partial pivoting; 0 on success */
static int sp_solve6(double A[6][7])
{
    for (int col = 0; col < 6; col++) {
        int piv = col;
        for (int r = col + 1; r < 6; r++)
            if (fabs(A[r][col]) > fabs(A[piv][col])) piv = r;
        if (fabs(A[piv][col]) < 1e-9) return -1;
        if (piv != col)
            for (int j = 0; j < 7; j++) {
                double t = A[col][j];
                A[col][j] = A[piv][j];
                A[piv][j] = t;
            }
        for (int r = 0; r < 6; r++) {
            if (r == col) continue;
            double f = A[r][col] / A[col][col];
            for (int j = col; j < 7; j++) A[r][j] -= f * A[col][j];
        }
    }
    for (int k = 0; k < 6; k++) A[k][6] /= A[k][k];
    return 0;
}

/* Jacobi-CG for (L + diag(w)) x = w*t on a dense wg x hg lattice restricted
 * to alive cells (mirror of param_extend's pe_cg_grid at toy scale) */
static void sp_cg(size_t wg, size_t hg, const uint8_t *alive,
                  const double *w, const double *t, double *x,
                  int iters, double tol)
{
    size_t n = wg * hg;
    double *diag = (double *)sp_xcalloc(n, sizeof(double));
    double *r = (double *)sp_xcalloc(n, sizeof(double));
    double *z = (double *)sp_xcalloc(n, sizeof(double));
    double *p = (double *)sp_xcalloc(n, sizeof(double));
    double *ap = (double *)sp_xcalloc(n, sizeof(double));
    double bn = 0, rz = 0;
    for (size_t i = 0; i < n; i++) {
        x[i] = 0.0;
        if (!alive[i]) continue;
        size_t cx = i % wg, cy = i / wg;
        double deg = 0;
        if (cx > 0 && alive[i - 1]) deg += 1;
        if (cx + 1 < wg && alive[i + 1]) deg += 1;
        if (cy > 0 && alive[i - wg]) deg += 1;
        if (cy + 1 < hg && alive[i + wg]) deg += 1;
        diag[i] = deg + w[i];
        double b = w[i] * t[i];
        r[i] = b;
        bn += b * b;
        z[i] = diag[i] > 1e-12 ? b / diag[i] : 0;
        p[i] = z[i];
        rz += r[i] * z[i];
    }
    if (bn > 0) {
        double btol = tol * tol * bn;
        for (int it = 0; it < iters; it++) {
            double pap = 0;
            for (size_t i = 0; i < n; i++) {
                if (!alive[i]) continue;
                size_t cx = i % wg, cy = i / wg;
                double s = diag[i] * p[i];
                if (cx > 0 && alive[i - 1]) s -= p[i - 1];
                if (cx + 1 < wg && alive[i + 1]) s -= p[i + 1];
                if (cy > 0 && alive[i - wg]) s -= p[i - wg];
                if (cy + 1 < hg && alive[i + wg]) s -= p[i + wg];
                ap[i] = s;
                pap += p[i] * s;
            }
            if (pap <= 0) break;
            double alpha = rz / pap, rr = 0;
            for (size_t i = 0; i < n; i++) {
                if (!alive[i]) continue;
                x[i] += alpha * p[i];
                r[i] -= alpha * ap[i];
                rr += r[i] * r[i];
            }
            if (rr < btol) break;
            double rz2 = 0;
            for (size_t i = 0; i < n; i++) {
                if (!alive[i]) continue;
                z[i] = diag[i] > 1e-12 ? r[i] / diag[i] : 0;
                rz2 += r[i] * z[i];
            }
            double beta = rz2 / rz;
            rz = rz2;
            for (size_t i = 0; i < n; i++)
                if (alive[i]) p[i] = z[i] + beta * p[i];
        }
    }
    free(diag); free(r); free(z); free(p); free(ap);
}

/* +-dir free-space march for bright papyrus (snap_grid ray semantics) */
static int sp_ray(CubeTable *ct, KDTree_T blk, const float *bv,
                  const double p[3], const double d[3], const SpOpts *o,
                  double band, double cv, double *out_t, double *out_s)
{
    double occ2 = o->occ_thresh * o->occ_thresh;
    double loc2 = o->local_r * o->local_r;
    for (double t = o->step; t <= o->reach + 1e-9; t += o->step) {
        float q[3] = { (float)(p[0] + t * d[0]), (float)(p[1] + t * d[1]),
                       (float)(p[2] + t * d[2]) };
        if (blk != NULL) {
            int32_t hit[8];
            size_t nh = KDTree_ball_query(blk, q, (float)occ2, hit, 8), h;
            for (h = 0; h < nh; h++) {
                double ez = (double)bv[hit[h] * 3 + 0] - p[0];
                double ey = (double)bv[hit[h] * 3 + 1] - p[1];
                double ex = (double)bv[hit[h] * 3 + 2] - p[2];
                if (ez * ez + ey * ey + ex * ex > loc2) return 0; /* blocked */
            }
        }
        double s = sample_trilinear(ct, q[0], q[1], q[2]);
        if (s > band && s - cv > o->min_gain) {
            *out_t = t;
            *out_s = s;
            return 1;
        }
    }
    return 0;
}

/* ---- per-sheet probe ------------------------------------------------------ */

typedef struct {
    int id;
    size_t nv, next_cells, n_guard, n_target, n_moved;
    double extent_a, extent_b, h_rms, snap_mean, ext_bright_frac;
    int W, H;              /* bake dims */
    uint8_t *png_pre, *png_post;   /* RGB, malloc'd */
} SheetOut;

static void sp_probe_sheet(const float *verts, size_t nv_all,
                           const int32_t *faces, size_t nf_all,
                           const int32_t *sheet_of, int sid,
                           CubeTable *ct, KDTree_T blk, const float *bv,
                           double win_lo, double win_hi,
                           const SpOpts *o, SheetOut *out)
{
    memset(out, 0, sizeof(*out));
    out->id = sid;

    /* gather sheet verts */
    size_t nv = 0;
    for (size_t i = 0; i < nv_all; i++)
        if (sheet_of[i] == sid) nv++;
    out->nv = nv;
    if (nv < 3) return;
    float *sv = (float *)sp_xmalloc(nv * 3 * sizeof(float));
    {
        size_t w = 0;
        for (size_t i = 0; i < nv_all; i++) {
            if (sheet_of[i] != sid) continue;
            sv[w * 3 + 0] = verts[i * 3 + 0];
            sv[w * 3 + 1] = verts[i * 3 + 1];
            sv[w * 3 + 2] = verts[i * 3 + 2];
            w++;
        }
    }
    (void)faces;
    (void)nf_all;

    /* local frame: PCA normal + in-plane basis, centroid origin */
    float nrm[3], cen[3];
    if (PCA_normal(sv, nv, nrm, cen) != 0) { free(sv); return; }
    float e1f[3], e2f[3];
    PCA_orthonormal_basis(nrm, e1f, e2f);
    double e1[3] = { e1f[0], e1f[1], e1f[2] };
    double e2[3] = { e2f[0], e2f[1], e2f[2] };
    double en[3] = { nrm[0], nrm[1], nrm[2] };
    double c0[3] = { cen[0], cen[1], cen[2] };

    /* project to (a,b,h) */
    double amin = 1e300, amax = -1e300, bmin = 1e300, bmax = -1e300;
    double *abh = (double *)sp_xmalloc(nv * 3 * sizeof(double));
    for (size_t i = 0; i < nv; i++) {
        double d[3] = { sv[i * 3 + 0] - c0[0], sv[i * 3 + 1] - c0[1],
                        sv[i * 3 + 2] - c0[2] };
        double a = d[0] * e1[0] + d[1] * e1[1] + d[2] * e1[2];
        double b = d[0] * e2[0] + d[1] * e2[1] + d[2] * e2[2];
        double h = d[0] * en[0] + d[1] * en[1] + d[2] * en[2];
        abh[i * 3 + 0] = a;
        abh[i * 3 + 1] = b;
        abh[i * 3 + 2] = h;
        if (a < amin) amin = a;
        if (a > amax) amax = a;
        if (b < bmin) bmin = b;
        if (b > bmax) bmax = b;
    }
    out->extent_a = amax - amin;
    out->extent_b = bmax - bmin;

    /* lattice (+band margin) */
    double cell = o->cell > 0.5 ? o->cell : 2.0;
    long W = (long)((amax - amin) / cell) + 1 + 2 * o->band;
    long H = (long)((bmax - bmin) / cell) + 1 + 2 * o->band;
    while (W > 0 && H > 0 && (size_t)W * (size_t)H > SP_MAX_GRID) {
        cell *= 2.0;
        W = (long)((amax - amin) / cell) + 1 + 2 * o->band;
        H = (long)((bmax - bmin) / cell) + 1 + 2 * o->band;
    }
    double a0 = amin - o->band * cell, b0 = bmin - o->band * cell;
    size_t n = (size_t)W * (size_t)H;
    double *hsum = (double *)sp_xcalloc(n, sizeof(double));
    double *hmin = (double *)sp_xmalloc(n * sizeof(double));
    double *hmax = (double *)sp_xmalloc(n * sizeof(double));
    uint32_t *hcnt = (uint32_t *)sp_xcalloc(n, sizeof(uint32_t));
    for (size_t i = 0; i < n; i++) { hmin[i] = 1e300; hmax[i] = -1e300; }
    double hrms = 0.0;
    for (size_t i = 0; i < nv; i++) {
        long cx = (long)((abh[i * 3 + 0] - a0) / cell);
        long cy = (long)((abh[i * 3 + 1] - b0) / cell);
        if (cx < 0 || cy < 0 || cx >= W || cy >= H) continue;
        size_t ci = (size_t)cy * (size_t)W + (size_t)cx;
        double h = abh[i * 3 + 2];
        hsum[ci] += h;
        hcnt[ci]++;
        if (h < hmin[ci]) hmin[ci] = h;
        if (h > hmax[ci]) hmax[ci] = h;
        hrms += h * h;
    }
    out->h_rms = sqrt(hrms / (double)nv);

    /* mask: cells with data and single-valued height (fold spread gate) */
    uint8_t *mask = (uint8_t *)sp_xcalloc(n, 1);
    size_t n_mask = 0;
    for (size_t i = 0; i < n; i++) {
        if (hcnt[i] == 0) continue;
        if (hmax[i] - hmin[i] > 3.0 * cell) continue;   /* folded in frame */
        mask[i] = 1;
        n_mask++;
    }
    if (n_mask < 8) {
        free(sv); free(abh); free(hsum); free(hmin); free(hmax);
        free(hcnt); free(mask);
        return;
    }

    /* chamfer distance to mask; ext = dist in [1, band] */
    uint16_t cap = (uint16_t)(o->band + 1);
    uint16_t *dist = (uint16_t *)sp_xmalloc(n * sizeof(uint16_t));
    for (size_t i = 0; i < n; i++)
        dist[i] = mask[i] ? (uint16_t)0 : cap;
    for (long cy = 0; cy < H; cy++)
        for (long cx = 0; cx < W; cx++) {
            size_t ci = (size_t)cy * (size_t)W + (size_t)cx;
            uint16_t d = dist[ci];
            if (cx > 0 && dist[ci - 1] + 1 < d) d = (uint16_t)(dist[ci - 1] + 1);
            if (cy > 0 && dist[ci - W] + 1 < d) d = (uint16_t)(dist[ci - W] + 1);
            dist[ci] = d;
        }
    for (long cy = H; cy-- > 0; )
        for (long cx = W; cx-- > 0; ) {
            size_t ci = (size_t)cy * (size_t)W + (size_t)cx;
            uint16_t d = dist[ci];
            if (cx + 1 < W && dist[ci + 1] + 1 < d)
                d = (uint16_t)(dist[ci + 1] + 1);
            if (cy + 1 < H && dist[ci + W] + 1 < d)
                d = (uint16_t)(dist[ci + W] + 1);
            dist[ci] = d;
        }

    /* quadratic LSQ predictor over mask cells */
    double G[6] = { 0, 0, 0, 0, 0, 0 };
    {
        double A[6][7];
        memset(A, 0, sizeof(A));
        for (size_t i = 0; i < n; i++) {
            if (!mask[i]) continue;
            double a = a0 + ((double)(i % (size_t)W) + 0.5) * cell;
            double b = b0 + ((double)(i / (size_t)W) + 0.5) * cell;
            double h = hsum[i] / (double)hcnt[i];
            double bas[6] = { 1.0, a, b, a * a, a * b, b * b };
            for (int r = 0; r < 6; r++) {
                for (int cc = 0; cc < 6; cc++) A[r][cc] += bas[r] * bas[cc];
                A[r][6] += bas[r] * h;
            }
        }
        if (sp_solve6(A) == 0)
            for (int k = 0; k < 6; k++) G[k] = A[k][6];
        else
            G[0] = 0.0;   /* fall back to plane h=0 */
    }
#define SP_PRED(a_, b_) (G[0] + G[1]*(a_) + G[2]*(b_) + G[3]*(a_)*(a_) \
                         + G[4]*(a_)*(b_) + G[5]*(b_)*(b_))

    /* harmonic residual over mask+ext */
    uint8_t *alive = (uint8_t *)sp_xcalloc(n, 1);
    double *wgt = (double *)sp_xcalloc(n, sizeof(double));
    double *tgt = (double *)sp_xcalloc(n, sizeof(double));
    double *x = (double *)sp_xmalloc(n * sizeof(double));
    size_t n_ext = 0;
    for (size_t i = 0; i < n; i++) {
        if (mask[i]) {
            double a = a0 + ((double)(i % (size_t)W) + 0.5) * cell;
            double b = b0 + ((double)(i / (size_t)W) + 0.5) * cell;
            alive[i] = 1;
            wgt[i] = 1e4;
            tgt[i] = hsum[i] / (double)hcnt[i] - SP_PRED(a, b);
        } else if (dist[i] <= (uint16_t)o->band) {
            alive[i] = 1;
            n_ext++;
        }
    }
    out->next_cells = n_ext;
    sp_cg((size_t)W, (size_t)H, alive, wgt, tgt, x, 400, 1e-5);

    /* full height per alive cell (mask: data mean; ext: pred+residual) --
     * also feeds the LOCAL surface normal (heightfield gradient), which is
     * what the probes must follow: the frame normal is only right at the
     * sheet centre of a curved sheet */
    double *hval = (double *)sp_xcalloc(n, sizeof(double));
    for (size_t i = 0; i < n; i++) {
        if (!alive[i]) continue;
        double a = a0 + ((double)(i % (size_t)W) + 0.5) * cell;
        double b = b0 + ((double)(i / (size_t)W) + 0.5) * cell;
        hval[i] = mask[i] ? hsum[i] / (double)hcnt[i] : SP_PRED(a, b) + x[i];
    }

    /* extension nodes -> world; guard vs other sheets */
    float *epos = (float *)sp_xmalloc(n * 3 * sizeof(float));
    float *enrm = (float *)sp_xmalloc(n * 3 * sizeof(float));
    uint8_t *eflag = (uint8_t *)sp_xcalloc(n, 1);   /* 1 ext 2 guard */
    double g2 = o->guard_r * o->guard_r;
    for (size_t i = 0; i < n; i++) {
        if (!alive[i] || mask[i]) continue;
        size_t cx = i % (size_t)W, cy = i / (size_t)W;
        double a = a0 + ((double)cx + 0.5) * cell;
        double b = b0 + ((double)cy + 0.5) * cell;
        double h = hval[i];
        float p[3];
        for (int k = 0; k < 3; k++)
            p[k] = (float)(c0[k] + a * e1[k] + b * e2[k] + h * en[k]);
        epos[i * 3 + 0] = p[0];
        epos[i * 3 + 1] = p[1];
        epos[i * 3 + 2] = p[2];
        if (i % 97 == 3)
            fprintf(stderr, "  [dbg] ext i=%zu ab=(%.1f,%.1f) h=%.2f "
                    "p=(%.1f,%.1f,%.1f)\n", i, a, b, h, p[0], p[1], p[2]);
        /* local normal from the height gradient (one-sided at borders) */
        double ha = 0.0, hb = 0.0;
        if (cx > 0 && cx + 1 < (size_t)W && alive[i - 1] && alive[i + 1])
            ha = (hval[i + 1] - hval[i - 1]) / (2.0 * cell);
        else if (cx + 1 < (size_t)W && alive[i + 1])
            ha = (hval[i + 1] - h) / cell;
        else if (cx > 0 && alive[i - 1])
            ha = (h - hval[i - 1]) / cell;
        if (cy > 0 && cy + 1 < (size_t)H && alive[i - (size_t)W]
            && alive[i + (size_t)W])
            hb = (hval[i + (size_t)W] - hval[i - (size_t)W]) / (2.0 * cell);
        else if (cy + 1 < (size_t)H && alive[i + (size_t)W])
            hb = (hval[i + (size_t)W] - h) / cell;
        else if (cy > 0 && alive[i - (size_t)W])
            hb = (h - hval[i - (size_t)W]) / cell;
        double nl[3];
        double nn = 0.0;
        for (int k = 0; k < 3; k++) {
            nl[k] = en[k] - ha * e1[k] - hb * e2[k];
            nn += nl[k] * nl[k];
        }
        nn = sqrt(nn);
        if (nn < 1e-9) { nl[0] = en[0]; nl[1] = en[1]; nl[2] = en[2]; nn = 1; }
        enrm[i * 3 + 0] = (float)(nl[0] / nn);
        enrm[i * 3 + 1] = (float)(nl[1] / nn);
        enrm[i * 3 + 2] = (float)(nl[2] / nn);
        eflag[i] = 1;
        if (blk != NULL) {
            int32_t hit[8];
            if (KDTree_ball_query(blk, p, (float)g2, hit, 8) > 0) {
                eflag[i] = 2;
                out->n_guard++;
            }
        }
    }

    /* snap: +-local-normal probes from unguarded ext nodes */
    double band = win_lo + o->band_frac * (win_hi - win_lo);
    if (o->do_snap && ct != NULL) {
        double msum = 0.0;
        for (size_t i = 0; i < n; i++) {
            if (eflag[i] != 1) continue;
            double p[3] = { epos[i * 3 + 0], epos[i * 3 + 1],
                            epos[i * 3 + 2] };
            double cv = sample_trilinear(ct, (float)p[0], (float)p[1],
                                         (float)p[2]);
            if (cv < 0) continue;
            double bt = 0, t, s;
            int f = 0;
            double dirs[2][3] = {
                { enrm[i * 3 + 0], enrm[i * 3 + 1], enrm[i * 3 + 2] },
                { -enrm[i * 3 + 0], -enrm[i * 3 + 1], -enrm[i * 3 + 2] }
            };
            for (int k = 0; k < 2; k++)
                if (sp_ray(ct, blk, bv, p, dirs[k], o, band, cv, &t, &s)) {
                    if (!f || t < bt) { bt = t; f = k + 1; }
                }
            if (f) {
                out->n_target++;
                double d = bt < o->max_disp ? bt : o->max_disp;
                for (int k = 0; k < 3; k++)
                    epos[i * 3 + k] = (float)(p[k]
                        + d * dirs[f - 1][k]);
                out->n_moved++;
                msum += d;
            }
        }
        out->snap_mean = out->n_moved > 0 ? msum / (double)out->n_moved : 0.0;
        fprintf(stderr, "  [dbg] snap try=%zu nocv=%zu band=%.1f\n",
                dbg_try, dbg_nocv, band);
    }

    /* bakes: pre = surface at prediction; post = after snap */
    out->W = (int)W;
    out->H = (int)H;
    uint8_t *png[2];
    png[0] = (uint8_t *)sp_xmalloc(n * 3);
    png[1] = (uint8_t *)sp_xmalloc(n * 3);
    size_t bright = 0, extpx = 0;
    for (int pass = 0; pass < 2; pass++) {
        for (size_t i = 0; i < n; i++) {
            uint8_t *px = &png[pass][i * 3];
            px[0] = 16; px[1] = 16; px[2] = 48;   /* navy bg */
            double p[3];
            int kind = 0;
            if (mask[i]) {
                double a = a0 + ((double)(i % (size_t)W) + 0.5) * cell;
                double b = b0 + ((double)(i / (size_t)W) + 0.5) * cell;
                double h = hsum[i] / (double)hcnt[i];
                for (int k = 0; k < 3; k++)
                    p[k] = c0[k] + a * e1[k] + b * e2[k] + h * en[k];
                kind = 1;
            } else if (eflag[i]) {
                if (pass == 0 && o->do_snap) {
                    /* pre view = interpolation surface */
                    double a = a0 + ((double)(i % (size_t)W) + 0.5) * cell;
                    double b = b0 + ((double)(i / (size_t)W) + 0.5) * cell;
                    for (int k = 0; k < 3; k++)
                        p[k] = c0[k] + a * e1[k] + b * e2[k]
                             + hval[i] * en[k];
                } else {
                    p[0] = epos[i * 3 + 0];
                    p[1] = epos[i * 3 + 1];
                    p[2] = epos[i * 3 + 2];
                }
                kind = eflag[i] == 2 ? 3 : 2;
            }
            if (!kind) continue;
            double s = ct != NULL
                     ? sample_trilinear(ct, (float)p[0], (float)p[1],
                                        (float)p[2]) : -1.0;
            if (pass == 0 && kind == 1 && i == n / 2)
                fprintf(stderr, "  [dbg] cell %zu p=(%.1f,%.1f,%.1f) s=%.1f "
                        "win=[%.0f,%.0f] ct=%p\n", i, p[0], p[1], p[2], s,
                        win_lo, win_hi, (void *)ct);
            double gy = 0.0;
            if (s >= 0 && win_hi > win_lo) {
                gy = (s - win_lo) / (win_hi - win_lo);
                if (gy < 0) gy = 0;
                if (gy > 1) gy = 1;
            }
            uint8_t g8 = (uint8_t)(gy * 255.0 + 0.5);
            if (kind == 1) {
                px[0] = g8; px[1] = g8; px[2] = g8;
            } else if (kind == 2) {
                px[0] = (uint8_t)(g8 * 3 / 10);
                px[1] = g8 > 36 ? g8 : 36;
                px[2] = (uint8_t)(g8 * 3 / 10);
                if (pass == 1) {
                    extpx++;
                    if (s >= band) bright++;
                }
            } else {
                px[0] = 200; px[1] = 40; px[2] = 40;
            }
        }
    }
    out->ext_bright_frac = extpx > 0 ? (double)bright / (double)extpx : 0.0;
    out->png_pre = png[0];
    out->png_post = png[1];

    free(sv); free(abh); free(hsum); free(hmin); free(hmax); free(hcnt);
    free(mask); free(dist); free(alive); free(wgt); free(tgt); free(x);
    free(hval); free(epos); free(enrm); free(eflag);
#undef SP_PRED
}

/* world OBJ: sheet faces + extension quads (post-snap) -- kept simple by
 * re-running the node math; study artifact, not pipeline data */
/* (extension quads omitted from OBJ v1: the PNG pair is the instrument;
 *  sheets are already in <id>_mesh.obj) */

/* ---- driver ---------------------------------------------------------------- */

static int sp_run_cube(Arena_T arena, const char *placed_dir,
                       const char *cube_id, const char *raw_dir,
                       const char *out_dir, const SpOpts *o)
{
    char path[1024];
    float *verts = NULL;
    int32_t *faces = NULL;
    size_t nv = 0, nf = 0;
    snprintf(path, sizeof(path), "%s/%s_mesh.obj", placed_dir, cube_id);
    if (ObjIO_read(arena, path, &verts, &nv, &faces, &nf) != 0) {
        fprintf(stderr, "sheet_probe: cannot read %s\n", path);
        return -1;
    }
    fprintf(stderr, "[%s] %zu verts %zu faces\n", cube_id, nv, nf);
    if (nv == 0 || nf == 0) return 0;

    /* sheets = connected components over face edges */
    assert(nv < (size_t)INT32_MAX);
    UnionFind uf = UF_new(arena, (int32_t)nv);
    for (size_t f = 0; f < nf; f++) {
        uf_union(&uf, faces[f * 3 + 0], faces[f * 3 + 1]);
        uf_union(&uf, faces[f * 3 + 0], faces[f * 3 + 2]);
    }
    int32_t *sheet_of = (int32_t *)ARENA_ALLOC(arena, nv * sizeof(int32_t));
    int32_t *root_id = (int32_t *)ARENA_ALLOC(arena, nv * sizeof(int32_t));
    size_t *cnt = (size_t *)ARENA_CALLOC(arena, nv, sizeof(size_t));
    for (size_t i = 0; i < nv; i++) root_id[i] = -1;
    int n_sheets = 0;
    for (size_t i = 0; i < nv; i++) {
        int32_t r = uf_find(&uf, (int32_t)i);
        if (root_id[r] < 0) root_id[r] = n_sheets++;
        sheet_of[i] = root_id[r];
        cnt[root_id[r]]++;
    }
    fprintf(stderr, "[%s] %d sheets (>= %zu verts: ", cube_id, n_sheets,
            o->min_verts);
    int n_big = 0;
    for (int s = 0; s < n_sheets; s++)
        if (cnt[s] >= o->min_verts) n_big++;
    fprintf(stderr, "%d)\n", n_big);

    /* RAW table over the cube bbox */
    CubeTable ct;
    int have_ct = 0;
    if (raw_dir != NULL) {
        if (cubetable_init(&ct, arena, raw_dir, o->raw_chunk, verts, nv,
                           o->reach + 4.0) == 0) {
            cubetable_prewarm_all(&ct);
            have_ct = 1;
        } else {
            fprintf(stderr, "[%s] WARN: no RAW table (%s)\n", cube_id,
                    raw_dir);
        }
    }

    /* contrast window: percentiles of samples at a vert subset */
    double win_lo = 0, win_hi = 255;
    if (have_ct) {
        size_t hist[256];
        memset(hist, 0, sizeof(hist));
        size_t m = 0;
        size_t stride = nv > 20000 ? nv / 20000 : 1;
        for (size_t i = 0; i < nv; i += stride) {
            double s = sample_trilinear(&ct, verts[i * 3 + 0],
                                        verts[i * 3 + 1], verts[i * 3 + 2]);
            if (s < 0) continue;
            int b = (int)(s + 0.5);
            if (b < 0) b = 0;
            if (b > 255) b = 255;
            hist[b]++;
            m++;
        }
        if (m > 0) {
            size_t cum = 0, t1 = m / 100, t99 = m - m / 100;
            int lo = 0, hi = 255;
            for (int b = 0; b < 256; b++) {
                cum += hist[b];
                if (cum <= t1) lo = b;
                if (cum < t99) hi = b;
            }
            if (hi > lo) { win_lo = lo; win_hi = hi; }
        }
    }

    char jp[1024];
    snprintf(jp, sizeof(jp), "%s/%s_sheetprobe.json", out_dir, cube_id);
    ves_ensure_parent_dir(jp);
    FILE *jf = fopen(jp, "w");
    if (jf != NULL)
        fprintf(jf, "{ \"cube\": \"%s\", \"sheets\": [\n", cube_id);

    /* montage accumulation */
    size_t mon_w = 0, mon_h = 0;
    SheetOut *outs = (SheetOut *)sp_xcalloc((size_t)n_sheets,
                                            sizeof(SheetOut));
    int n_done = 0;
    for (int s = 0; s < n_sheets; s++) {
        if (cnt[s] < o->min_verts) continue;
        /* blocker tree = all OTHER sheets' verts */
        size_t nb = nv - cnt[s];
        KDTree_T blk = NULL;
        float *bv = NULL;
        if (nb > 0) {
            bv = (float *)ARENA_ALLOC(arena, nb * 3 * sizeof(float));
            size_t w = 0;
            for (size_t i = 0; i < nv; i++) {
                if (sheet_of[i] == s) continue;
                bv[w * 3 + 0] = verts[i * 3 + 0];
                bv[w * 3 + 1] = verts[i * 3 + 1];
                bv[w * 3 + 2] = verts[i * 3 + 2];
                w++;
            }
            blk = KDTree_new(arena, bv, nb);
        }
        SheetOut *so = &outs[n_done];
        sp_probe_sheet(verts, nv, faces, nf, sheet_of, s,
                       have_ct ? &ct : NULL, blk, bv, win_lo, win_hi, o, so);
        if (so->png_post == NULL) continue;
        fprintf(stderr,
                "[%s] sheet %d: nv=%zu ext %.0fx%.0f h_rms=%.1f | ext=%zu "
                "cells guard=%zu | snap tgt=%zu moved=%zu mean=%.2f | "
                "ext bright=%.0f%%\n",
                cube_id, s, so->nv, so->extent_a, so->extent_b, so->h_rms,
                so->next_cells, so->n_guard, so->n_target, so->n_moved,
                so->snap_mean, 100.0 * so->ext_bright_frac);
        if (jf != NULL)
            fprintf(jf, "  %s{ \"id\": %d, \"nv\": %zu, \"ext_cells\": %zu,"
                    " \"guard\": %zu, \"snap_target\": %zu, \"moved\": %zu,"
                    " \"mean_disp\": %.3f, \"h_rms\": %.2f,"
                    " \"ext_bright_frac\": %.4f }",
                    n_done > 0 ? "," : "", s, so->nv, so->next_cells,
                    so->n_guard, so->n_target, so->n_moved, so->snap_mean,
                    so->h_rms, so->ext_bright_frac);
        char pp[1024];
        snprintf(pp, sizeof(pp), "%s/%s_sheet%03d_pre.png", out_dir,
                 cube_id, s);
        VesPng_write_rgb(pp, so->png_pre, so->W, so->H);
        snprintf(pp, sizeof(pp), "%s/%s_sheet%03d_post.png", out_dir,
                 cube_id, s);
        VesPng_write_rgb(pp, so->png_post, so->W, so->H);
        if ((size_t)so->W > mon_w) mon_w = (size_t)so->W;
        mon_h += (size_t)so->H + 6;
        n_done++;
    }
    if (jf != NULL) {
        fprintf(jf, "\n] }\n");
        fclose(jf);
    }

    /* vertical montage of post bakes */
    if (n_done > 0 && mon_w > 0 && mon_h < 30000) {
        uint8_t *mon = (uint8_t *)sp_xcalloc(mon_w * mon_h, 3);
        size_t y = 0;
        for (int k = 0; k < n_done; k++) {
            SheetOut *so = &outs[k];
            for (int r = 0; r < so->H; r++)
                memcpy(&mon[((y + (size_t)r) * mon_w) * 3],
                       &so->png_post[(size_t)r * (size_t)so->W * 3],
                       (size_t)so->W * 3);
            y += (size_t)so->H + 6;
        }
        char mp[1024];
        snprintf(mp, sizeof(mp), "%s/%s_sheets_post.png", out_dir, cube_id);
        VesPng_write_rgb(mp, mon, (int)mon_w, (int)mon_h);
        fprintf(stderr, "[%s] montage %zux%zu -> %s (%d sheets)\n", cube_id,
                mon_w, mon_h, mp, n_done);
        free(mon);
    }
    for (int k = 0; k < n_done; k++) {
        free(outs[k].png_pre);
        free(outs[k].png_post);
    }
    free(outs);
    return 0;
}

/* ---- selftest -------------------------------------------------------------- */

#define SP_ST_DIR "output/_selftest_sheet_probe"

static int sp_check(int cond, const char *what, int *fails)
{
    if (!cond) {
        fprintf(stderr, "[sheet_probe selftest]   FAIL: %s\n", what);
        (*fails)++;
    }
    return cond;
}

static int sp_selftest(void)
{
    int fails = 0;
    Arena_T arena = Arena_new();

    /* RAW: bright shell at r=26 around (32,32) in a 64^3 cube */
    {
        long ch = 64;
        uint8_t *vol = (uint8_t *)sp_xmalloc((size_t)(ch * ch * ch));
        for (long z = 0; z < ch; z++)
            for (long y = 0; y < ch; y++)
                for (long x = 0; x < ch; x++) {
                    double r = sqrt(((double)y - 32.0) * ((double)y - 32.0)
                                    + ((double)x - 32.0) * ((double)x - 32.0));
                    vol[(z * ch + y) * ch + x] =
                        fabs(r - 26.0) <= 0.8 ? (uint8_t)220 : (uint8_t)15;
                }
        char p[512];
        snprintf(p, sizeof(p), SP_ST_DIR "/raw/z00000_y00000_x00000.tif");
        ves_ensure_parent_dir(p);
        if (TiffIO_save(p, vol, (int)ch, (int)ch, (int)ch) != 0)
            sp_check(0, "raw fixture", &fails);
        free(vol);
    }

    /* mesh: an arc sheet ON the r=26 shell (theta 0.15..0.7pi, but 2 vox
     * INSIDE it at r=24 so the snap has work), z 16..48, PLUS a blocker
     * fragment on the arc continuation just past the sheet edge */
    {
        int nth = 41, nz = 17;
        size_t nvfx = (size_t)nth * nz + 16;
        float *v = (float *)sp_xmalloc(nvfx * 3 * sizeof(float));
        int32_t *f = (int32_t *)sp_xmalloc(((size_t)(nth - 1) * (nz - 1) * 6
                                            + 16 * 3) * sizeof(int32_t));
        size_t vi = 0, fi = 0;
        for (int iz = 0; iz < nz; iz++)
            for (int it = 0; it < nth; it++) {
                double th = 0.15 * 3.14159265 + 0.55 * 3.14159265
                          * (double)it / (nth - 1);
                double z = 16.0 + 32.0 * (double)iz / (nz - 1);
                v[vi * 3 + 0] = (float)z;
                v[vi * 3 + 1] = (float)(32.0 + 24.0 * cos(th));
                v[vi * 3 + 2] = (float)(32.0 + 24.0 * sin(th));
                vi++;
            }
        for (int iz = 0; iz + 1 < nz; iz++)
            for (int it = 0; it + 1 < nth; it++) {
                int32_t a = iz * nth + it, b = a + 1, c = a + nth,
                        d = c + 1;
                f[fi * 3 + 0] = a; f[fi * 3 + 1] = b; f[fi * 3 + 2] = d;
                fi++;
                f[fi * 3 + 0] = a; f[fi * 3 + 1] = d; f[fi * 3 + 2] = c;
                fi++;
            }
        /* blocker: fragment on the arc continuation ~4-10 vox past the
         * theta=0.7pi sheet edge, same radius */
        size_t bv0 = vi;
        for (int iz = 0; iz < 4; iz++)
            for (int it = 0; it < 4; it++) {
                double th = (0.72 + 0.03 * it) * 3.14159265;
                double z = 26.0 + 4.0 * iz;
                v[vi * 3 + 0] = (float)z;
                v[vi * 3 + 1] = (float)(32.0 + 24.0 * cos(th));
                v[vi * 3 + 2] = (float)(32.0 + 24.0 * sin(th));
                vi++;
            }
        for (int iz = 0; iz < 3; iz++)
            for (int it = 0; it < 3; it++) {
                int32_t a = (int32_t)(bv0 + (size_t)iz * 4 + it);
                f[fi * 3 + 0] = a; f[fi * 3 + 1] = a + 1;
                f[fi * 3 + 2] = a + 5;
                fi++;
            }
        char p[512];
        snprintf(p, sizeof(p), SP_ST_DIR "/z00000_y00000_x00000_mesh.obj");
        ves_ensure_parent_dir(p);
        sp_check(ObjIO_write(p, v, vi, f, fi) == 0, "mesh fixture", &fails);
        free(v);
        free(f);
    }

    SpOpts o;
    sp_opts_default(&o);
    o.raw_chunk = 64;
    o.band = 8;
    o.min_verts = 100;
    int rc = sp_run_cube(arena, SP_ST_DIR, "z00000_y00000_x00000",
                         SP_ST_DIR "/raw", SP_ST_DIR "/out", &o);
    sp_check(rc == 0, "run rc", &fails);
    /* the extension continues the r=24 arc; the shell at 26 is 2 vox out
     * along the LOCAL normal -> snap must find targets; the blocker on the
     * continuation must trip the guard */
    {
        char p[512];
        snprintf(p, sizeof(p), SP_ST_DIR "/out/"
                 "z00000_y00000_x00000_sheetprobe.json");
        FILE *fp = fopen(p, "rb");
        sp_check(fp != NULL, "json written", &fails);
        if (fp != NULL) {
            char buf[4096];
            size_t got = fread(buf, 1, sizeof(buf) - 1, fp);
            buf[got] = '\0';
            fclose(fp);
            sp_check(strstr(buf, "\"ext_cells\"") != NULL, "json has stats",
                     &fails);
            sp_check(strstr(buf, "\"guard\": 0,") == NULL,
                     "blocker guarded some ext nodes", &fails);
            sp_check(strstr(buf, "\"snap_target\": 0,") == NULL,
                     "snap found bright targets", &fails);
        }
    }

    Arena_dispose(&arena);
    fprintf(stderr, "[sheet_probe selftest] %s (%d failure%s)\n",
            fails == 0 ? "PASSED" : "FAILED", fails, fails == 1 ? "" : "s");
    return fails;
}

static void usage(void)
{
    fprintf(stderr,
        "usage: sheet_probe <placed_dir> <cube_id> --raw <cubes_RAW>\n"
        "       [--out DIR] [--cell F] [--band N] [--min-verts N]\n"
        "       [--guard F] [--reach F] [--max-disp F] [--raw-chunk N]\n"
        "       [--no-snap]\n"
        "       sheet_probe --selftest\n");
}

int main(int argc, char **argv)
{
    if (argc >= 2 && strcmp(argv[1], "--selftest") == 0)
        return sp_selftest() == 0 ? 0 : 1;
    if (argc < 3) {
        usage();
        return 2;
    }
    const char *placed = argv[1];
    const char *cube = argv[2];
    const char *raw = NULL;
    const char *out = "output/tools/sheetprobe";
    SpOpts o;
    sp_opts_default(&o);
    for (int i = 3; i < argc; i++) {
        if (strcmp(argv[i], "--raw") == 0 && i + 1 < argc) raw = argv[++i];
        else if (strcmp(argv[i], "--out") == 0 && i + 1 < argc)
            out = argv[++i];
        else if (strcmp(argv[i], "--cell") == 0 && i + 1 < argc)
            o.cell = atof(argv[++i]);
        else if (strcmp(argv[i], "--band") == 0 && i + 1 < argc)
            o.band = atoi(argv[++i]);
        else if (strcmp(argv[i], "--min-verts") == 0 && i + 1 < argc)
            o.min_verts = (size_t)atoll(argv[++i]);
        else if (strcmp(argv[i], "--guard") == 0 && i + 1 < argc)
            o.guard_r = atof(argv[++i]);
        else if (strcmp(argv[i], "--reach") == 0 && i + 1 < argc)
            o.reach = atof(argv[++i]);
        else if (strcmp(argv[i], "--max-disp") == 0 && i + 1 < argc)
            o.max_disp = atof(argv[++i]);
        else if (strcmp(argv[i], "--raw-chunk") == 0 && i + 1 < argc)
            o.raw_chunk = atol(argv[++i]);
        else if (strcmp(argv[i], "--no-snap") == 0) o.do_snap = 0;
        else fprintf(stderr, "WARN: ignoring %s\n", argv[i]);
    }
    Arena_T arena = Arena_new();
    int rc = sp_run_cube(arena, placed, cube, raw, out, &o);
    Arena_dispose(&arena);
    return rc == 0 ? 0 : 1;
}
