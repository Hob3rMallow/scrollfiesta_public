/* strip_metrics.c -- row-streamed strip QC. See strip_metrics.h. */
#include "strip_metrics.h"
#include "../common/tiff_io.h"
#include "../common/ves_platform.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* covered = painted (ok or multi); 2/3/6 are uncovered-with-reason */
static int sm_covered(uint8_t cls)
{
    return cls == 1 || cls == 4;
}

static int sm_cmp_i32(const void *a, const void *b)
{
    int32_t x = *(const int32_t *)a, y = *(const int32_t *)b;
    return x < y ? -1 : (x > y ? 1 : 0);
}

int StripMetrics_seam_cols(Arena_T arena, const PieceSet *ps,
                           double u0, double du, size_t W,
                           int32_t **out_cols, size_t *out_n)
{
    assert(arena && ps && out_cols && out_n);
    *out_cols = NULL;
    *out_n = 0;
    if (ps->nf == 0 || ps->n_cubes == 0 || du <= 0.0 || W < 8)
        return 0;

    long X0 = (long)floor(u0 / du + 0.5);

    /* verts referenced by kept faces */
    uint8_t *has = (uint8_t *)ARENA_CALLOC(arena, ps->nv + 1, 1);
    for (size_t f = 0; f < ps->nf; f++)
        for (int e = 0; e < 3; e++)
            has[(size_t)ps->faces[f * 3 + e]] = 1;

    /* per (cube, gid) u-interval; gid table grown per cube */
    enum { GID_CAP = 4096 };
    double *gmin = (double *)ARENA_ALLOC(arena, GID_CAP * sizeof(double));
    double *gmax = (double *)ARENA_ALLOC(arena, GID_CAP * sizeof(double));

    size_t cap = 1024, n = 0;
    int32_t *cols = (int32_t *)ARENA_ALLOC(arena, cap * sizeof(int32_t));

    for (size_t c = 0; c < ps->n_cubes; c++) {
        size_t v0 = ps->cube_voff[c], v1 = ps->cube_voff[c + 1];
        int32_t gid_hi = -1;
        for (size_t i = v0; i < v1; i++) {
            if (!has[i]) continue;
            int32_t g = ps->gid != NULL ? ps->gid[i] : 0;
            if (g < 0) g = 0;      /* no sidecar: whole cube = one interval */
            if (g >= GID_CAP) continue;
            if (g > gid_hi) {
                for (int32_t k = gid_hi + 1; k <= g; k++) {
                    gmin[k] = 1e300;
                    gmax[k] = -1e300;
                }
                gid_hi = g;
            }
            double u = (double)ps->uv[i * 2 + 0];
            if (u < gmin[g]) gmin[g] = u;
            if (u > gmax[g]) gmax[g] = u;
        }
        for (int32_t g = 0; g <= gid_hi; g++) {
            if (gmax[g] < gmin[g]) continue;
            long pxs[2];
            pxs[0] = (long)floor(gmin[g] / du + 0.5) - X0;
            pxs[1] = (long)floor(gmax[g] / du + 0.5) - X0;
            for (int e = 0; e < 2; e++) {
                long px = pxs[e];
                if (px < 2 || px + 2 >= (long)W) continue;
                if (n == cap) {
                    size_t nc = cap * 2;
                    int32_t *na = (int32_t *)ARENA_ALLOC(arena,
                                                         nc * sizeof(int32_t));
                    memcpy(na, cols, n * sizeof(int32_t));
                    cols = na;
                    cap = nc;
                }
                cols[n++] = (int32_t)px;
            }
        }
    }
    if (n == 0)
        return 0;

    qsort(cols, n, sizeof(int32_t), sm_cmp_i32);
    size_t w = 0;
    for (size_t i = 0; i < n; i++)
        if (w == 0 || cols[i] != cols[w - 1])
            cols[w++] = cols[i];
    *out_cols = cols;
    *out_n = w;
    return 0;
}

int StripMetrics_compute(Arena_T arena, const char *rawtex_tif,
                         const char *diagclass_tif,
                         const int32_t *seam_cols, size_t ncols,
                         int dark_thresh, double **out_col_dcol,
                         StripMetrics *out)
{
    assert(arena && rawtex_tif && diagclass_tif && out);
    memset(out, 0, sizeof(*out));
    if (out_col_dcol != NULL)
        *out_col_dcol = NULL;

    TiffRowReader ri, rc;
    if (TiffIO_rows_open(arena, rawtex_tif, &ri) != 0)
        return -1;
    if (TiffIO_rows_open(arena, diagclass_tif, &rc) != 0) {
        TiffIO_rows_close(&ri);
        return -1;
    }
    if (ri.W != rc.W || ri.H != rc.H) {
        TiffIO_rows_close(&ri);
        TiffIO_rows_close(&rc);
        return -1;
    }
    int W = ri.W, H = ri.H;

    uint8_t *img = (uint8_t *)ARENA_ALLOC(arena, (size_t)W);
    uint8_t *cls = (uint8_t *)ARENA_ALLOC(arena, (size_t)W);

    double *csum = NULL;
    size_t *ccnt = NULL;
    if (ncols > 0) {
        csum = (double *)ARENA_CALLOC(arena, ncols, sizeof(double));
        ccnt = (size_t *)ARENA_CALLOC(arena, ncols, sizeof(size_t));
    }

    size_t covered = 0, multi = 0, dark = 0, synth = 0;
    double base_sum = 0.0;
    size_t base_cnt = 0;
    int rcod = 0;

    for (int y = 0; y < H && rcod == 0; y++) {
        if (TiffIO_rows_read(&ri, y, 0, W, img) != 0
            || TiffIO_rows_read(&rc, y, 0, W, cls) != 0) {
            rcod = -1;
            break;
        }
        for (int x = 0; x < W; x++) {
            if (cls[x] == 7) { synth++; continue; }
            if (!sm_covered(cls[x])) continue;
            covered++;
            if (cls[x] == 4) multi++;
            if ((int)img[x] < dark_thresh) dark++;
        }
        for (size_t k = 0; k < ncols; k++) {
            int x = (int)seam_cols[k];
            if (x < 0 || x + 1 >= W) continue;
            if (sm_covered(cls[x]) && sm_covered(cls[x + 1])) {
                csum[k] += fabs((double)img[x + 1] - (double)img[x]);
                ccnt[k]++;
            }
        }
        for (int x = 1; x + 1 < W; x += 97) {
            if (sm_covered(cls[x]) && sm_covered(cls[x + 1])) {
                base_sum += fabs((double)img[x + 1] - (double)img[x]);
                base_cnt++;
            }
        }
    }
    TiffIO_rows_close(&ri);
    TiffIO_rows_close(&rc);
    if (rcod != 0)
        return -1;

    out->fill = (double)covered / ((double)W * (double)H);
    out->multi_frac = covered ? (double)multi / (double)covered : 0.0;
    out->dark_frac = covered ? (double)dark / (double)covered : 0.0;
    out->synth_px = synth;
    out->synth_frac = (covered + synth)
                    ? (double)synth / (double)(covered + synth) : 0.0;
    out->fill_total = ((double)covered + (double)synth)
                    / ((double)W * (double)H);
    out->base_dcol_mean = base_cnt ? base_sum / (double)base_cnt : 0.0;

    double ssum = 0.0;
    size_t scols = 0, spairs = 0;
    double *col_dcol = NULL;
    if (ncols > 0 && out_col_dcol != NULL)
        col_dcol = (double *)ARENA_ALLOC(arena, ncols * sizeof(double));
    size_t zoom_floor = (size_t)H / 8;   /* per-col report needs substance:
                                            sliver edges with a few covered
                                            pairs make useless zoom picks */
    for (size_t k = 0; k < ncols; k++) {
        double m = ccnt[k] ? csum[k] / (double)ccnt[k] : -1.0;
        if (col_dcol != NULL)
            col_dcol[k] = ccnt[k] >= zoom_floor ? m : -1.0;
        if (ccnt[k] > 0) {
            ssum += m;
            scols++;
            spairs += ccnt[k];
        }
    }
    out->seam_dcol_mean = scols ? ssum / (double)scols : 0.0;
    out->seam_ratio = (out->base_dcol_mean > 1e-9 && scols)
                    ? out->seam_dcol_mean / out->base_dcol_mean : 0.0;
    out->n_seam_cols = scols;
    out->n_pairs_used = spairs;
    if (out_col_dcol != NULL)
        *out_col_dcol = col_dcol;
    return 0;
}

int StripMetrics_vseams(Arena_T arena, const char *rawtex_tif,
                        const char *diagclass_tif,
                        double v0, long chunk, StripMetrics *inout)
{
    assert(arena && rawtex_tif && diagclass_tif && inout);
    inout->vseam_gap_fill = 0.0;
    inout->vseam_ratio = 0.0;
    inout->n_vseam_rows = 0;
    if (chunk <= 0)
        return 0;

    TiffRowReader ri, rc;
    if (TiffIO_rows_open(arena, rawtex_tif, &ri) != 0)
        return -1;
    if (TiffIO_rows_open(arena, diagclass_tif, &rc) != 0) {
        TiffIO_rows_close(&ri);
        return -1;
    }
    int W = ri.W, H = ri.H;
    uint8_t *ia = (uint8_t *)ARENA_ALLOC(arena, (size_t)W);
    uint8_t *ib = (uint8_t *)ARENA_ALLOC(arena, (size_t)W);
    uint8_t *ca = (uint8_t *)ARENA_ALLOC(arena, (size_t)W);
    uint8_t *cb = (uint8_t *)ARENA_ALLOC(arena, (size_t)W);
    uint8_t *cm = (uint8_t *)ARENA_ALLOC(arena, (size_t)W);

    /* one (rows y-1, y+1 for jump; row y for coverage) measurement */
    double seam_cov = 0.0, ref_cov = 0.0;
    double seam_jump = 0.0, ref_jump = 0.0;
    size_t seam_jn = 0, ref_jn = 0;
    int rc2 = 0;

    long k0 = (long)ceil((v0 + 9.0) / (double)chunk);
    for (long k = k0; rc2 == 0; k++) {
        long y = k * chunk - (long)v0;
        if (y >= H - 9) break;
        inout->n_vseam_rows++;
        /* seam row coverage + jump across it */
        if (TiffIO_rows_read(&rc, (int)y, 0, W, cm) != 0
            || TiffIO_rows_read(&ri, (int)y - 1, 0, W, ia) != 0
            || TiffIO_rows_read(&ri, (int)y + 1, 0, W, ib) != 0
            || TiffIO_rows_read(&rc, (int)y - 1, 0, W, ca) != 0
            || TiffIO_rows_read(&rc, (int)y + 1, 0, W, cb) != 0) {
            rc2 = -1;
            break;
        }
        size_t cov = 0;
        for (int x = 0; x < W; x++) {
            if (sm_covered(cm[x])) cov++;
            if (sm_covered(ca[x]) && sm_covered(cb[x])) {
                seam_jump += fabs((double)ib[x] - (double)ia[x]);
                seam_jn++;
            }
        }
        seam_cov += (double)cov;
        /* reference rows y +- 8 (mean of both) */
        for (int s = -1; s <= 1; s += 2) {
            long yr = y + 8 * s;
            if (TiffIO_rows_read(&rc, (int)yr, 0, W, cm) != 0
                || TiffIO_rows_read(&ri, (int)yr - 1, 0, W, ia) != 0
                || TiffIO_rows_read(&ri, (int)yr + 1, 0, W, ib) != 0
                || TiffIO_rows_read(&rc, (int)yr - 1, 0, W, ca) != 0
                || TiffIO_rows_read(&rc, (int)yr + 1, 0, W, cb) != 0) {
                rc2 = -1;
                break;
            }
            size_t rcov = 0;
            for (int x = 0; x < W; x++) {
                if (sm_covered(cm[x])) rcov++;
                if (sm_covered(ca[x]) && sm_covered(cb[x])) {
                    ref_jump += fabs((double)ib[x] - (double)ia[x]);
                    ref_jn++;
                }
            }
            ref_cov += 0.5 * (double)rcov;
        }
    }
    TiffIO_rows_close(&ri);
    TiffIO_rows_close(&rc);
    if (rc2 != 0)
        return -1;
    inout->vseam_gap_fill = ref_cov > 0.0 ? seam_cov / ref_cov : 0.0;
    double sj = seam_jn ? seam_jump / (double)seam_jn : 0.0;
    double rj = ref_jn ? ref_jump / (double)ref_jn : 0.0;
    inout->vseam_ratio = rj > 1e-9 ? sj / rj : 0.0;
    return 0;
}

/* ============================================================================
 * Self-test.
 * ==========================================================================*/

static int sm_check(int cond, const char *what, int *fails)
{
    if (!cond) {
        fprintf(stderr, "[strip_metrics selftest]   FAIL: %s\n", what);
        (*fails)++;
    }
    return cond;
}

int StripMetrics_selftest(void)
{
    int fails = 0;
    Arena_T arena = Arena_new();
    const char *dir = "output/_selftest_strip_metrics";
    char pi[1024], pc[1024];
    snprintf(pi, sizeof(pi), "%s/img.tif", dir);
    ves_ensure_parent_dir(pi);
    snprintf(pc, sizeof(pc), "%s/cls.tif", dir);

    /* 64x16: left half 100, right half 200 (step at x=32); one dark column;
     * classes: mostly 1, a multi stripe at x=10, uncovered ends */
    enum { TW = 64, TH = 16 };
    uint8_t img[TW * TH], cls[TW * TH];
    for (int y = 0; y < TH; y++) {
        for (int x = 0; x < TW; x++) {
            img[y * TW + x] = (uint8_t)(x < 32 ? 100 : 200);
            cls[y * TW + x] = 1;
        }
        img[y * TW + 20] = 5;            /* dark column (covered) */
        cls[y * TW + 10] = 4;            /* multi column */
        cls[y * TW + 0] = 0;             /* uncovered edge */
        cls[y * TW + TW - 1] = 7;        /* synthetic-fill edge */
    }
    sm_check(TiffIO_save(pi, img, 1, TH, TW) == 0, "write img", &fails);
    sm_check(TiffIO_save(pc, cls, 1, TH, TW) == 0, "write cls", &fails);

    /* seam at the step (x=31 -> pair 31/32 jumps 100), baseline flat */
    int32_t cols[2] = { 31, 45 };
    StripMetrics m;
    double *cd = NULL;
    int rc = StripMetrics_compute(arena, pi, pc, cols, 2, 40, &cd, &m);
    fprintf(stderr, "[strip_metrics selftest] rc=%d fill=%.3f multi=%.3f "
            "dark=%.4f seam=%.1f base=%.2f ratio=%.1f cols=%zu\n",
            rc, m.fill, m.multi_frac, m.dark_frac, m.seam_dcol_mean,
            m.base_dcol_mean, m.seam_ratio, m.n_seam_cols);
    sm_check(rc == 0, "compute rc", &fails);
    sm_check(fabs(m.fill - (double)(TW - 2) / TW) < 1e-6, "fill", &fails);
    sm_check(fabs(m.multi_frac - 1.0 / (TW - 2)) < 1e-6, "multi_frac", &fails);
    sm_check(m.synth_px == TH, "synth_px", &fails);
    sm_check(fabs(m.synth_frac - (double)TH / (double)((TW - 1) * TH)) < 1e-9,
             "synth_frac", &fails);
    sm_check(fabs(m.fill_total - (double)(TW - 1) / TW) < 1e-9, "fill_total",
             &fails);
    sm_check(m.dark_frac > 0.0 && m.dark_frac < 0.05, "dark_frac", &fails);
    sm_check(cd != NULL && fabs(cd[0] - 100.0) < 1e-9 && fabs(cd[1]) < 1e-9,
             "per-col dcol", &fails);
    sm_check(m.n_seam_cols == 2 && m.seam_dcol_mean > 49.0, "seam mean",
             &fails);

    /* v-seams: 64x40 image with an UNCOVERED gap row at y=20 (v0=0,
     * chunk=20 -> seam rows y=20 only, since y=40 is out of range);
     * then a covered variant -> gap_fill jumps to ~1 */
    {
        enum { VW = 64, VH = 40 };
        uint8_t vimg[VW * VH], vcls[VW * VH];
        memset(vimg, 100, sizeof(vimg));
        memset(vcls, 1, sizeof(vcls));
        for (int x = 0; x < VW; x++) vcls[20 * VW + x] = 0;   /* the gap */
        char vp1[1024], vp2[1024];
        snprintf(vp1, sizeof(vp1), "%s/vimg.tif", dir);
        snprintf(vp2, sizeof(vp2), "%s/vcls.tif", dir);
        sm_check(TiffIO_save(vp1, vimg, 1, VH, VW) == 0, "vseam write img",
                 &fails);
        sm_check(TiffIO_save(vp2, vcls, 1, VH, VW) == 0, "vseam write cls",
                 &fails);
        StripMetrics vm;
        memset(&vm, 0, sizeof(vm));
        sm_check(StripMetrics_vseams(arena, vp1, vp2, 0.0, 20, &vm) == 0,
                 "vseam rc", &fails);
        sm_check(vm.n_vseam_rows == 1 && vm.vseam_gap_fill < 0.05,
                 "vseam gap detected", &fails);
        for (int x = 0; x < VW; x++) vcls[20 * VW + x] = 1;   /* closed */
        sm_check(TiffIO_save(vp2, vcls, 1, VH, VW) == 0, "vseam rewrite",
                 &fails);
        sm_check(StripMetrics_vseams(arena, vp1, vp2, 0.0, 20, &vm) == 0
                 && vm.vseam_gap_fill > 0.95 && vm.vseam_ratio <= 1.5,
                 "vseam closed reads ~1", &fails);
    }

    /* seam_cols from a synthetic 2-cube piece set: cube A u=[0,10] gid 0,
     * cube B u=[20,30] gid 0 -> endpoints {0,10,20,30}, canvas-end dropped */
    {
        PieceSet ps;
        memset(&ps, 0, sizeof(ps));
        float verts[8 * 3] = { 0 };
        float uv[8 * 2] = {
            0, 0, 10, 0, 0, 5, 10, 5,     /* cube A */
            20, 0, 30, 0, 20, 5, 30, 5    /* cube B */
        };
        int32_t faces[4 * 3] = { 0, 1, 2, 1, 3, 2, 4, 5, 6, 5, 7, 6 };
        int32_t fc[4] = { 0, 0, 1, 1 };
        int32_t gid[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
        size_t voff[3] = { 0, 4, 8 };
        long org[2][3] = { { 0, 0, 0 }, { 0, 0, 128 } };
        ps.verts = verts; ps.uv = uv; ps.faces = faces; ps.face_cube = fc;
        ps.gid = gid; ps.cube_voff = voff; ps.cube_org = org;
        ps.nv = 8; ps.nf = 4; ps.n_cubes = 2;
        int32_t *sc = NULL;
        size_t nsc = 0;
        sm_check(StripMetrics_seam_cols(arena, &ps, 0.0, 1.0, 40, &sc, &nsc)
                 == 0, "seam_cols rc", &fails);
        /* px 0 dropped (canvas end), expect {10, 20, 30}; 30 kept (40-px W) */
        sm_check(nsc == 3 && sc[0] == 10 && sc[1] == 20 && sc[2] == 30,
                 "seam_cols endpoints", &fails);
    }

    Arena_dispose(&arena);
    fprintf(stderr, "[strip_metrics selftest] %s (%d failure%s)\n",
            fails == 0 ? "PASSED" : "FAILED", fails, fails == 1 ? "" : "s");
    return fails;
}
