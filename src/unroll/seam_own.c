/* seam_own.c -- SeamOwn ownership pass. See seam_own.h.
 *
 * Face-domain GLOBAL pass (no banding: ~10 GB scratch at 52M faces fits the
 * dev box; only the per-region rasters touch RAW and they are grid_cap-
 * bounded). Reuses overlap_quilt's exported subroutines for adjacency,
 * radius, local rasters and the Efros boundary energy; the classification
 * differs from Quilt_run in exactly the ways the whole grid demands:
 *   - pre-weld, cross-cube faces share no verts (adjacency is intra-cube
 *     only); after the v4 step2_join weld they CAN -- mesh-adjacent
 *     cross-cube pairs are the same stitched surface and are exempted in
 *     both the classifier and the multi recount. Consequence for layers-as-
 *     CCs: a welded same-wrap region spanning cubes is ONE component, which
 *     is correct (it IS one sheet); intra-CC folds still hit the multicut;
 *   - same-wrap seam double-paint (|dr| <= gate) is INVISIBLE to the radius
 *     gate -- it gets its own SEAM classification + min-error-boundary cut;
 *   - within a TRUE region, layers are connected components (the multicut is
 *     needed only when one intra-cube CC overlaps itself: a fused fold);
 *   - the energy-only pick is replaced by anchor-radius majority vote first
 *     (energy as tiebreak): a turn-off offender's own fringe can win a pure
 *     energy contest against a small victim chart.
 */
#include "seam_own.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/pca.h"
#include "../common/tiff_io.h"
#include "../common/ves_platform.h"
#include "../flatten/overlap_quilt.h"
#include "../split/multicut_wrap.h"
#include "../whole/cube_register.h"   /* CubeReg_deltaU (header-inline) */

#define SO_TWO_PI 6.28318530717958647692

/* ------------------------------------------------------------------ util */

static void *so_xmalloc(size_t n)
{
    void *p = malloc(n > 0 ? n : 1);
    if (p == NULL) { fprintf(stderr, "seam_own: OOM %zu\n", n); exit(1); }
    return p;
}

static void *so_xcalloc(size_t c, size_t s)
{
    void *p = calloc(c > 0 ? c : 1, s);
    if (p == NULL) { fprintf(stderr, "seam_own: OOM %zux%zu\n", c, s); exit(1); }
    return p;
}

static void *so_xrealloc(void *p, size_t n)
{
    void *q = realloc(p, n > 0 ? n : 1);
    if (q == NULL) { fprintf(stderr, "seam_own: OOM %zu\n", n); exit(1); }
    return q;
}

typedef struct { int32_t *p; } SoUF;
static void so_uf_init(SoUF *uf, size_t n)
{
    uf->p = (int32_t *)so_xmalloc(n * sizeof(int32_t));
    for (size_t i = 0; i < n; i++) uf->p[i] = (int32_t)i;
}
static int32_t so_uf_find(SoUF *uf, int32_t x)
{
    while (uf->p[x] != x) { uf->p[x] = uf->p[uf->p[x]]; x = uf->p[x]; }
    return x;
}
static void so_uf_union(SoUF *uf, int32_t a, int32_t b)
{
    int32_t ra = so_uf_find(uf, a), rb = so_uf_find(uf, b);
    if (ra != rb) uf->p[rb] = ra;
}
static void so_uf_free(SoUF *uf) { free(uf->p); uf->p = NULL; }

static int so_cmp_u64(const void *a, const void *b)
{
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return (x < y) ? -1 : (x > y);
}

void SeamOwnOpts_default(SeamOwnOpts *o)
{
    assert(o);
    memset(o, 0, sizeof(*o));
    o->axis_point[0] = 0.0f;
    o->axis_point[1] = 3405.0f;
    o->axis_point[2] = 2878.0f;
    o->axis_dir[0] = 1.0f;
    o->raw_dir = NULL;
    o->raw_chunk = 128;
    o->ct = NULL;
    o->pitch = 9.5;
    o->radius_gate = 0.0;
    o->cell = 2.0;
    o->seam_gate3d = 8.0;
    o->tex_du = 1.0;
    o->tex_dv = 1.0;
    o->normal_range = 2.0;
    o->normal_samples = 5;
    o->region_cap = 200000;
    o->grid_cap = 3000;
    o->seam_cut = 1;
    o->seam_halfband = 4;
    o->min_overlap_px = 8;
    o->rehome = 1;
    o->rehome_min_faces = 50;
    o->rehome_frac_tol = 0.35;
    o->verbose = 0;
}

/* crude key scan (the placed_index writer emits one key per pattern) */
static const char *so_jfind(const char *s, const char *key)
{
    char pat[128];
    snprintf(pat, sizeof(pat), "\"%s\"", key);
    const char *p = strstr(s, pat);
    if (p == NULL) return NULL;
    p = strchr(p + strlen(pat), ':');
    return p != NULL ? p + 1 : NULL;
}

int SeamOwn_read_index(const char *placed_dir, SeamOwnOpts *o)
{
    assert(placed_dir && o);
    char path[1024];
    snprintf(path, sizeof(path), "%s/placed_index.json", placed_dir);
    FILE *f = fopen(path, "rb");
    if (f == NULL) return -1;
    /* the header keys sit well inside the first 4 KB */
    char buf[4096];
    size_t n = fread(buf, 1, sizeof(buf) - 1, f);
    fclose(f);
    buf[n] = '\0';

    const char *p = NULL;
    double a = 0, b = 0, c = 0;
    if ((p = so_jfind(buf, "axis_point_zyx")) != NULL
        && sscanf(p, " [ %lf , %lf , %lf", &a, &b, &c) == 3) {
        o->axis_point[0] = (float)a;
        o->axis_point[1] = (float)b;
        o->axis_point[2] = (float)c;
    }
    if ((p = so_jfind(buf, "axis_dir_zyx")) != NULL
        && sscanf(p, " [ %lf , %lf , %lf", &a, &b, &c) == 3) {
        o->axis_dir[0] = (float)a;
        o->axis_dir[1] = (float)b;
        o->axis_dir[2] = (float)c;
    }
    if ((p = so_jfind(buf, "pitch")) != NULL && sscanf(p, " %lf", &a) == 1)
        o->pitch = a;
    if ((p = so_jfind(buf, "spiral_a")) != NULL && sscanf(p, " %lf", &a) == 1)
        o->spiral_a = a;
    if ((p = so_jfind(buf, "spiral_b")) != NULL && sscanf(p, " %lf", &a) == 1)
        o->spiral_b = a;
    return 0;
}

/* --------------------------------------------------- cover-only rasterizer */

/* coverage of faces flist[nfl] on the local grid -- Quilt_raster_faces minus
 * the RAW sampling (fill-safety needs cover, not values) */
static void so_cover_faces(const int32_t *flist, size_t nfl,
                           const float *uv, const int32_t *faces,
                           double u0, double v0, double du, double dv,
                           int gw, int gh, uint8_t *cov)
{
    memset(cov, 0, (size_t)gw * (size_t)gh);
    for (size_t fi = 0; fi < nfl; fi++) {
        int32_t f = flist[fi];
        size_t a = (size_t)faces[f * 3 + 0];
        size_t b = (size_t)faces[f * 3 + 1];
        size_t c = (size_t)faces[f * 3 + 2];
        double ua = ((double)uv[a * 2 + 0] - u0) / du;
        double va = ((double)uv[a * 2 + 1] - v0) / dv;
        double ub = ((double)uv[b * 2 + 0] - u0) / du;
        double vb = ((double)uv[b * 2 + 1] - v0) / dv;
        double uc = ((double)uv[c * 2 + 0] - u0) / du;
        double vc = ((double)uv[c * 2 + 1] - v0) / dv;
        double A2 = (ub - ua) * (vc - va) - (vb - va) * (uc - ua);
        if (fabs(A2) < 1e-12) continue;
        double lox = ua < ub ? (ua < uc ? ua : uc) : (ub < uc ? ub : uc);
        double hix = ua > ub ? (ua > uc ? ua : uc) : (ub > uc ? ub : uc);
        double loy = va < vb ? (va < vc ? va : vc) : (vb < vc ? vb : vc);
        double hiy = va > vb ? (va > vc ? va : vc) : (vb > vc ? vb : vc);
        int x0 = (int)ceil(lox - 0.001), x1 = (int)floor(hix + 0.001);
        int y0 = (int)ceil(loy - 0.001), y1 = (int)floor(hiy + 0.001);
        if (x0 < 0) x0 = 0;
        if (y0 < 0) y0 = 0;
        if (x1 >= gw) x1 = gw - 1;
        if (y1 >= gh) y1 = gh - 1;
        for (int yy = y0; yy <= y1; yy++) {
            for (int xx = x0; xx <= x1; xx++) {
                double pu = (double)xx, pv = (double)yy;
                double l1 = ((pu - ua) * (vc - va) - (pv - va) * (uc - ua)) / A2;
                double l2 = ((ub - ua) * (pv - va) - (vb - va) * (pu - ua)) / A2;
                double l0 = 1.0 - l1 - l2;
                if (l0 < -1e-6 || l1 < -1e-6 || l2 < -1e-6) continue;
                cov[(size_t)yy * (size_t)gw + (size_t)xx] = 1;
            }
        }
    }
}

/* ------------------------------------------------------------ pair stores */

typedef struct { int32_t fa, fb; } SoPair;

typedef struct {
    SoPair *p;
    size_t n, cap;
} SoPairList;

static void so_pl_init(SoPairList *l)
{
    l->cap = 256;
    l->n = 0;
    l->p = (SoPair *)so_xmalloc(l->cap * sizeof(SoPair));
}

static void so_pl_push(SoPairList *l, int32_t a, int32_t b)
{
    if (l->n == l->cap) {
        l->cap *= 2;
        l->p = (SoPair *)so_xrealloc(l->p, l->cap * sizeof(SoPair));
    }
    l->p[l->n].fa = a;
    l->p[l->n].fb = b;
    l->n++;
}

/* --------------------------------------------------------------- seam DP */

/* One seam region: rasters of the two dominant cubes on a shared grid, DP
 * min-error boundary through the both-cover band, face drops on the losing
 * side, fill-safety restore. Mutates keep/dec. */
typedef struct {
    size_t dropped, restored;
    int abut_only;
    int skipped;
} SoSeamStat;

static void so_seam_region(const PieceSet *ps, const SeamOwnOpts *o,
                           CubeTable *ct, const float *nrm,
                           const int32_t *rface, size_t nrf,
                           int32_t cube_a, int32_t cube_b,
                           uint8_t *keep, uint8_t *dec, SoSeamStat *st)
{
    memset(st, 0, sizeof(*st));
    double du = o->tex_du > 0 ? o->tex_du : 1.0;
    double dv = o->tex_dv > 0 ? o->tex_dv : 1.0;

    /* A/B face lists + shared uv bbox */
    int32_t *fA = (int32_t *)so_xmalloc(nrf * sizeof(int32_t));
    int32_t *fB = (int32_t *)so_xmalloc(nrf * sizeof(int32_t));
    size_t nA = 0, nB = 0;
    double gu0 = 1e300, gv0 = 1e300, gu1 = -1e300, gv1 = -1e300;
    for (size_t i = 0; i < nrf; i++) {
        int32_t f = rface[i];
        int32_t cb = ps->face_cube[f];
        if (cb == cube_a) fA[nA++] = f;
        else if (cb == cube_b) fB[nB++] = f;
        else continue;
        for (int k = 0; k < 3; k++) {
            size_t vi = (size_t)ps->faces[f * 3 + k];
            double u = (double)ps->uv[vi * 2 + 0];
            double v = (double)ps->uv[vi * 2 + 1];
            if (u < gu0) gu0 = u;
            if (u > gu1) gu1 = u;
            if (v < gv0) gv0 = v;
            if (v > gv1) gv1 = v;
        }
    }
    if (nA == 0 || nB == 0) { free(fA); free(fB); return; }
    int pad = o->seam_halfband > 0 ? o->seam_halfband : 4;
    gu0 -= pad * du; gv0 -= pad * dv;
    gu1 += pad * du; gv1 += pad * dv;
    int gw = (int)((gu1 - gu0) / du) + 1;
    int gh = (int)((gv1 - gv0) / dv) + 1;
    if (gw < 3 || gh < 3 || gw > o->grid_cap || gh > o->grid_cap) {
        st->skipped = 1;
        free(fA); free(fB);
        return;
    }

    size_t np = (size_t)gw * (size_t)gh;
    float *imgA = (float *)so_xmalloc(np * sizeof(float));
    float *imgB = (float *)so_xmalloc(np * sizeof(float));
    uint8_t *covA = (uint8_t *)so_xmalloc(np);
    uint8_t *covB = (uint8_t *)so_xmalloc(np);
    Quilt_raster_faces(fA, nA, ps->verts, ps->uv, nrm, ps->faces, ct,
                       o->normal_range, o->normal_samples,
                       gu0, gv0, du, dv, gw, gh, imgA, covA);
    Quilt_raster_faces(fB, nB, ps->verts, ps->uv, nrm, ps->faces, ct,
                       o->normal_range, o->normal_samples,
                       gu0, gv0, du, dv, gw, gh, imgB, covB);

    /* both-cover band M + its bbox + offset + mean cost */
    size_t nM = 0;
    int mx0 = gw, mx1 = -1, my0 = gh, my1 = -1;
    double sA = 0.0, sB = 0.0;
    for (int y = 0; y < gh; y++) {
        for (int x = 0; x < gw; x++) {
            size_t pi = (size_t)y * (size_t)gw + (size_t)x;
            if (covA[pi] && covB[pi]) {
                nM++;
                sA += imgA[pi];
                sB += imgB[pi];
                if (x < mx0) mx0 = x;
                if (x > mx1) mx1 = x;
                if (y < my0) my0 = y;
                if (y > my1) my1 = y;
            }
        }
    }
    if (nM < (size_t)(o->min_overlap_px > 0 ? o->min_overlap_px : 8)) {
        st->abut_only = 1;
        free(imgA); free(imgB); free(covA); free(covB); free(fA); free(fB);
        return;
    }
    double off = (sA - sB) / (double)nM;
    double csum = 0.0;
    for (int y = my0; y <= my1; y++) {
        for (int x = mx0; x <= mx1; x++) {
            size_t pi = (size_t)y * (size_t)gw + (size_t)x;
            if (covA[pi] && covB[pi])
                csum += fabs((double)imgA[pi] - off - (double)imgB[pi]);
        }
    }
    double c_out = 2.0 * (csum / (double)nM) + 1.0;

    /* orientation: tall band -> vertical path x(y); wide band -> y(x).
     * The horizontal case is the vertical case on a transposed view. */
    int vertical = (my1 - my0) >= (mx1 - mx0);
    int n_along = vertical ? (my1 - my0 + 1) : (mx1 - mx0 + 1);
    int c_lo = vertical ? mx0 : my0;
    int c_hi = vertical ? mx1 : my1;
    int n_cross = c_hi - c_lo + 1;
    double *E = (double *)so_xmalloc((size_t)n_along * (size_t)n_cross
                                     * sizeof(double));
    int32_t *cut = (int32_t *)so_xmalloc((size_t)n_along * sizeof(int32_t));

    for (int t = 0; t < n_along; t++) {
        for (int c = 0; c < n_cross; c++) {
            int x = vertical ? (c_lo + c) : (mx0 + t);
            int y = vertical ? (my0 + t) : (c_lo + c);
            size_t pi = (size_t)y * (size_t)gw + (size_t)x;
            double cost = (covA[pi] && covB[pi])
                        ? fabs((double)imgA[pi] - off - (double)imgB[pi])
                        : c_out;
            double best = 0.0;
            if (t > 0) {
                best = E[(size_t)(t - 1) * (size_t)n_cross + (size_t)c];
                if (c > 0) {
                    double e = E[(size_t)(t - 1) * (size_t)n_cross
                                 + (size_t)(c - 1)];
                    if (e < best) best = e;
                }
                if (c + 1 < n_cross) {
                    double e = E[(size_t)(t - 1) * (size_t)n_cross
                                 + (size_t)(c + 1)];
                    if (e < best) best = e;
                }
            }
            E[(size_t)t * (size_t)n_cross + (size_t)c] = cost + best;
        }
    }
    /* backtrack */
    {
        int t = n_along - 1, bc = 0;
        for (int c = 1; c < n_cross; c++)
            if (E[(size_t)t * (size_t)n_cross + (size_t)c]
                < E[(size_t)t * (size_t)n_cross + (size_t)bc]) bc = c;
        cut[t] = c_lo + bc;
        for (t = n_along - 2; t >= 0; t--) {
            int pc = cut[t + 1] - c_lo, nb = pc;
            for (int d = -1; d <= 1; d++) {
                int c = pc + d;
                if (c < 0 || c >= n_cross) continue;
                if (E[(size_t)t * (size_t)n_cross + (size_t)c]
                    < E[(size_t)t * (size_t)n_cross + (size_t)nb]) nb = c;
            }
            cut[t] = c_lo + nb;
        }
    }
    free(E);

    /* which side is A's? by median cross-coordinate of covered px */
    double mA = 0.0, mB = 0.0;
    size_t cA = 0, cB = 0;
    for (int y = 0; y < gh; y++) {
        for (int x = 0; x < gw; x++) {
            size_t pi = (size_t)y * (size_t)gw + (size_t)x;
            double cc = vertical ? (double)x : (double)y;
            if (covA[pi]) { mA += cc; cA++; }
            if (covB[pi]) { mB += cc; cB++; }
        }
    }
    mA /= (double)(cA > 0 ? cA : 1);
    mB /= (double)(cB > 0 ? cB : 1);
    int a_low = mA <= mB;   /* A owns the low-cross side of the cut */

    /* tentative drops: faces whose bbox meets the M band and whose centroid
     * sits on the OTHER side of the cut */
    for (int pass = 0; pass < 2; pass++) {
        const int32_t *fl = pass == 0 ? fA : fB;
        size_t nl = pass == 0 ? nA : nB;
        int own_low = pass == 0 ? a_low : !a_low;
        for (size_t i = 0; i < nl; i++) {
            int32_t f = fl[i];
            double cu = 0.0, cv = 0.0;
            double bu0 = 1e300, bu1 = -1e300, bv0 = 1e300, bv1 = -1e300;
            for (int k = 0; k < 3; k++) {
                size_t vi = (size_t)ps->faces[f * 3 + k];
                double u = ((double)ps->uv[vi * 2 + 0] - gu0) / du;
                double v = ((double)ps->uv[vi * 2 + 1] - gv0) / dv;
                cu += u; cv += v;
                if (u < bu0) bu0 = u;
                if (u > bu1) bu1 = u;
                if (v < bv0) bv0 = v;
                if (v > bv1) bv1 = v;
            }
            cu /= 3.0; cv /= 3.0;
            if (bu1 < mx0 || bu0 > mx1 || bv1 < my0 || bv0 > my1)
                continue;   /* face never enters the both-cover band */
            int t = vertical ? (int)(cv + 0.5) - my0 : (int)(cu + 0.5) - mx0;
            if (t < 0) t = 0;
            if (t >= n_along) t = n_along - 1;
            double cross = vertical ? cu : cv;
            int on_low = cross < (double)cut[t];
            if (on_low != own_low) {
                keep[f] = 0;
                dec[f] = 5;
                st->dropped++;
            }
        }
    }

    /* fill-safety: the keep set must still cover everything A|B covered */
    if (st->dropped > 0) {
        uint8_t *covK = (uint8_t *)so_xmalloc(np);
        int32_t *fk = (int32_t *)so_xmalloc((nA + nB) * sizeof(int32_t));
        size_t nk = 0;
        for (size_t i = 0; i < nA; i++) if (keep[fA[i]]) fk[nk++] = fA[i];
        for (size_t i = 0; i < nB; i++) if (keep[fB[i]]) fk[nk++] = fB[i];
        so_cover_faces(fk, nk, ps->uv, ps->faces, gu0, gv0, du, dv,
                       gw, gh, covK);
        for (int y = 0; y < gh; y++) {
            for (int x = 0; x < gw; x++) {
                size_t pi = (size_t)y * (size_t)gw + (size_t)x;
                if (!(covA[pi] || covB[pi]) || covK[pi])
                    continue;
                /* restore any dropped face whose bbox contains this px */
                for (int pass = 0; pass < 2; pass++) {
                    const int32_t *fl = pass == 0 ? fA : fB;
                    size_t nl = pass == 0 ? nA : nB;
                    for (size_t i = 0; i < nl; i++) {
                        int32_t f = fl[i];
                        if (keep[f]) continue;
                        double bu0 = 1e300, bu1 = -1e300;
                        double bv0 = 1e300, bv1 = -1e300;
                        for (int k = 0; k < 3; k++) {
                            size_t vi = (size_t)ps->faces[f * 3 + k];
                            double u = ((double)ps->uv[vi * 2 + 0] - gu0) / du;
                            double v = ((double)ps->uv[vi * 2 + 1] - gv0) / dv;
                            if (u < bu0) bu0 = u;
                            if (u > bu1) bu1 = u;
                            if (v < bv0) bv0 = v;
                            if (v > bv1) bv1 = v;
                        }
                        if ((double)x >= bu0 - 0.5 && (double)x <= bu1 + 0.5
                            && (double)y >= bv0 - 0.5
                            && (double)y <= bv1 + 0.5) {
                            keep[f] = 1;
                            dec[f] = 4;
                            st->restored++;
                            st->dropped--;
                        }
                    }
                }
            }
        }
        free(covK);
        free(fk);
    }

    /* survivors of the band get the seam-kept mark for the diag */
    for (size_t i = 0; i < nA; i++)
        if (keep[fA[i]] && dec[fA[i]] == 0) dec[fA[i]] = 4;
    for (size_t i = 0; i < nB; i++)
        if (keep[fB[i]] && dec[fB[i]] == 0) dec[fB[i]] = 4;

    free(imgA); free(imgB); free(covA); free(covB); free(fA); free(fB);
    free(cut);
}

/* ---------------------------------------------------------------- run */

int SeamOwn_run(Arena_T arena, PieceSet *ps,
                const SeamOwnOpts *opts, SeamOwnResult *out)
{
    assert(arena && ps && opts && out);
    memset(out, 0, sizeof(*out));
    double t_start = ves_clock_sec();

    size_t nf = ps->nf, nv = ps->nv;
    out->face_keep = (uint8_t *)ARENA_ALLOC(arena, nf + 1);
    memset(out->face_keep, 1, nf + 1);
    out->face_dec = (uint8_t *)ARENA_CALLOC(arena, nf + 1, 1);
    out->face_region = (int32_t *)ARENA_ALLOC(arena,
                                              (nf + 1) * sizeof(int32_t));
    for (size_t f = 0; f < nf; f++) out->face_region[f] = -1;
    if (nf == 0 || nv == 0)
        return 0;

    const SeamOwnOpts *o = opts;
    double gate = o->radius_gate > 0.0 ? o->radius_gate
                : (o->pitch > 0.0 ? o->pitch * 0.5 : 4.75);
    double gate3d2 = o->seam_gate3d * o->seam_gate3d;

    /* per-face radius + uv/3D centroids */
    float e1f[3], e2f[3];
    float adir[3] = { o->axis_dir[0], o->axis_dir[1], o->axis_dir[2] };
    PCA_orthonormal_basis(adir, e1f, e2f);
    double e1[3] = { e1f[0], e1f[1], e1f[2] };
    double e2[3] = { e2f[0], e2f[1], e2f[2] };

    float *rf = (float *)so_xmalloc(nf * sizeof(float));
    float *fu = (float *)so_xmalloc(nf * sizeof(float));
    float *fv = (float *)so_xmalloc(nf * sizeof(float));
    float *fc3 = (float *)so_xmalloc(nf * 3 * sizeof(float));
    double umin = 1e300, vmin = 1e300, umax = -1e300, vmax = -1e300;
    for (size_t f = 0; f < nf; f++) {
        double r = 0.0, su = 0.0, sv = 0.0, s3[3] = { 0, 0, 0 };
        for (int k = 0; k < 3; k++) {
            size_t vi = (size_t)ps->faces[f * 3 + k];
            r += Quilt_point_radius(&ps->verts[vi * 3], o->axis_point, e1, e2);
            su += ps->uv[vi * 2 + 0];
            sv += ps->uv[vi * 2 + 1];
            for (int d = 0; d < 3; d++) s3[d] += ps->verts[vi * 3 + d];
        }
        rf[f] = (float)(r / 3.0);
        fu[f] = (float)(su / 3.0);
        fv[f] = (float)(sv / 3.0);
        for (int d = 0; d < 3; d++) fc3[f * 3 + d] = (float)(s3[d] / 3.0);
        if (fu[f] < umin) umin = fu[f];
        if (fu[f] > umax) umax = fu[f];
        if (fv[f] < vmin) vmin = fv[f];
        if (fv[f] > vmax) vmax = fv[f];
    }

    /* UV cell CSR over face centroids */
    double cell = o->cell > 0.0 ? o->cell : 2.0;
    size_t gu = (size_t)((umax - umin) / cell) + 1;
    size_t gv = (size_t)((vmax - vmin) / cell) + 1;
    if (gu < 1) gu = 1;
    if (gv < 1) gv = 1;
    size_t ncell = gu * gv;
    int32_t *coff = (int32_t *)so_xcalloc(ncell + 1, sizeof(int32_t));
    int32_t *fcell = (int32_t *)so_xmalloc(nf * sizeof(int32_t));
    for (size_t f = 0; f < nf; f++) {
        size_t cu = (size_t)(((double)fu[f] - umin) / cell);
        size_t cv = (size_t)(((double)fv[f] - vmin) / cell);
        if (cu >= gu) cu = gu - 1;
        if (cv >= gv) cv = gv - 1;
        fcell[f] = (int32_t)(cv * gu + cu);
        coff[fcell[f] + 1]++;
    }
    for (size_t i = 0; i < ncell; i++) coff[i + 1] += coff[i];
    int32_t *ccur = (int32_t *)so_xmalloc(ncell * sizeof(int32_t));
    memcpy(ccur, coff, ncell * sizeof(int32_t));
    int32_t *clist = (int32_t *)so_xmalloc(nf * sizeof(int32_t));
    for (size_t f = 0; f < nf; f++) clist[ccur[fcell[f]]++] = (int32_t)f;
    free(ccur);
    free(fcell);

    /* mesh adjacency + sorted pair keys over ALL faces. Pre-weld, cross-cube
     * faces never share verts (intra-cube only); after step2_join the welded
     * seam faces DO -- their pairs land here and are exempted from the
     * conflict classification (same stitched surface, not double-paint). */
    size_t nadj = 0;
    QuiltAdjEdge *adj = Quilt_build_adjacency(ps->faces, nf, ps->uv, &nadj);
    uint64_t *akeys = (uint64_t *)so_xmalloc((nadj + 1) * sizeof(uint64_t));
    for (size_t i = 0; i < nadj; i++) {
        int32_t a = adj[i].fa, b = adj[i].fb;
        int32_t lo = a < b ? a : b, hi = a < b ? b : a;
        akeys[i] = ((uint64_t)(uint32_t)lo << 32) | (uint32_t)hi;
    }
    qsort(akeys, nadj, sizeof(uint64_t), so_cmp_u64);

    /* ---- classification -------------------------------------------------- */
    SoPairList tp, sp;
    so_pl_init(&tp);
    so_pl_init(&sp);
    uint8_t *is_true = (uint8_t *)so_xcalloc(nf, 1);
    uint8_t *is_seam = (uint8_t *)so_xcalloc(nf, 1);
    size_t multi_before = 0;

    for (size_t c = 0; c < ncell; c++) {
        int32_t a0 = coff[c], a1 = coff[c + 1];
        int cell_multi = 0;
        for (int32_t i = a0; i < a1; i++) {
            for (int32_t j = i + 1; j < a1; j++) {
                int32_t fa = clist[i], fb = clist[j];
                double dr = fabs((double)rf[fa] - (double)rf[fb]);
                if (ps->face_cube[fa] == ps->face_cube[fb]) {
                    if (dr <= gate) continue;   /* one sheet: self-gates */
                    /* adjacency exemption (fold hinge line) */
                    int32_t lo = fa < fb ? fa : fb, hi = fa < fb ? fb : fa;
                    uint64_t key = ((uint64_t)(uint32_t)lo << 32)
                                 | (uint32_t)hi;
                    size_t los = 0, his = nadj;
                    while (los < his) {
                        size_t m = (los + his) / 2;
                        if (akeys[m] < key) los = m + 1;
                        else his = m;
                    }
                    if (los < nadj && akeys[los] == key) continue;
                    so_pl_push(&tp, fa, fb);
                    is_true[fa] = 1;
                    is_true[fb] = 1;
                    cell_multi = 1;
                } else {
                    /* welded meshes: a mesh-adjacent cross-cube pair IS the
                     * same stitched surface -- never a conflict (the early
                     * join stage creates exactly these) */
                    {
                        int32_t lo = fa < fb ? fa : fb;
                        int32_t hi = fa < fb ? fb : fa;
                        uint64_t key = ((uint64_t)(uint32_t)lo << 32)
                                     | (uint32_t)hi;
                        size_t los = 0, his = nadj;
                        while (los < his) {
                            size_t m = (los + his) / 2;
                            if (akeys[m] < key) los = m + 1;
                            else his = m;
                        }
                        if (los < nadj && akeys[los] == key) continue;
                    }
                    if (dr > gate) {
                        so_pl_push(&tp, fa, fb);   /* wrap collision/turn-off */
                        is_true[fa] = 1;
                        is_true[fb] = 1;
                        cell_multi = 1;
                    } else {
                        double d2 = 0.0;
                        for (int d = 0; d < 3; d++) {
                            double dd = (double)fc3[fa * 3 + d]
                                      - (double)fc3[fb * 3 + d];
                            d2 += dd * dd;
                        }
                        if (d2 <= gate3d2) {
                            so_pl_push(&sp, fa, fb);   /* seam double-paint */
                            is_seam[fa] = 1;
                            is_seam[fb] = 1;
                            cell_multi = 1;
                        } else {
                            out->n_mystery_pairs++;
                        }
                    }
                }
            }
        }
        if (cell_multi) multi_before++;
    }
    out->n_true_pairs = tp.n;
    out->n_seam_pairs = sp.n;
    out->multi_cells_before = multi_before;

    /* texture (shared table if given) */
    CubeTable ct_local;
    CubeTable *ct = o->ct;
    float *nrm = NULL;
    int have_ct = 0;
    if (ct != NULL) {
        have_ct = 1;
    } else if (o->raw_dir != NULL
               && cubetable_init(&ct_local, arena, o->raw_dir, o->raw_chunk,
                                 ps->verts, nv, o->normal_range + 2.0) == 0) {
        ct = &ct_local;
        have_ct = 1;
    }
    out->have_texture = have_ct;
    if (have_ct)
        nrm = (float *)ps->normals;   /* PieceSet already carries normals */

    /* ---- TRUE-overlap regions -------------------------------------------- */
    SoUF uf;
    so_uf_init(&uf, nf);
    for (size_t i = 0; i < tp.n; i++)
        so_uf_union(&uf, tp.p[i].fa, tp.p[i].fb);
    for (size_t i = 0; i < nadj; i++)
        if (is_true[adj[i].fa] && is_true[adj[i].fb])
            so_uf_union(&uf, adj[i].fa, adj[i].fb);

    int32_t *region_of = (int32_t *)so_xmalloc(nf * sizeof(int32_t));
    int32_t *root2rid = (int32_t *)so_xmalloc(nf * sizeof(int32_t));
    for (size_t f = 0; f < nf; f++) { region_of[f] = -1; root2rid[f] = -1; }
    int32_t nreg = 0;
    for (size_t f = 0; f < nf; f++) {
        if (!is_true[f]) continue;
        int32_t r = so_uf_find(&uf, (int32_t)f);
        if (root2rid[r] < 0) root2rid[r] = nreg++;
        region_of[f] = root2rid[r];
        out->face_region[f] = region_of[f];
    }
    out->n_regions = (size_t)nreg;

    /* region CSR */
    int32_t *roff = (int32_t *)so_xcalloc((size_t)nreg + 1, sizeof(int32_t));
    for (size_t f = 0; f < nf; f++)
        if (region_of[f] >= 0) roff[region_of[f] + 1]++;
    for (int32_t i = 0; i < nreg; i++) roff[i + 1] += roff[i];
    int32_t *rcur = (int32_t *)so_xmalloc(((size_t)nreg + 1) * sizeof(int32_t));
    memcpy(rcur, roff, ((size_t)nreg + 1) * sizeof(int32_t));
    int32_t *rfaces = (int32_t *)so_xmalloc(
        (roff[nreg] > 0 ? (size_t)roff[nreg] : 1) * sizeof(int32_t));
    for (size_t f = 0; f < nf; f++)
        if (region_of[f] >= 0) rfaces[rcur[region_of[f]]++] = (int32_t)f;
    free(rcur);

    /* true pairs by region (for the intra-CC fold test) */
    int32_t *poff = (int32_t *)so_xcalloc((size_t)nreg + 1, sizeof(int32_t));
    for (size_t i = 0; i < tp.n; i++)
        poff[region_of[tp.p[i].fa] + 1]++;
    for (int32_t i = 0; i < nreg; i++) poff[i + 1] += poff[i];
    int32_t *pcur = (int32_t *)so_xmalloc(((size_t)nreg + 1) * sizeof(int32_t));
    memcpy(pcur, poff, ((size_t)nreg + 1) * sizeof(int32_t));
    int32_t *plist = (int32_t *)so_xmalloc((tp.n > 0 ? tp.n : 1)
                                           * sizeof(int32_t));
    for (size_t i = 0; i < tp.n; i++)
        plist[pcur[region_of[tp.p[i].fa]]++] = (int32_t)i;
    free(pcur);

    /* anchors: non-overlap faces mesh-adjacent to a region */
    int32_t *anoff = (int32_t *)so_xcalloc((size_t)nreg + 1, sizeof(int32_t));
    for (size_t i = 0; i < nadj; i++) {
        int32_t fa = adj[i].fa, fb = adj[i].fb;
        if (region_of[fa] >= 0 && region_of[fb] < 0) anoff[region_of[fa] + 1]++;
        else if (region_of[fb] >= 0 && region_of[fa] < 0)
            anoff[region_of[fb] + 1]++;
    }
    for (int32_t i = 0; i < nreg; i++) anoff[i + 1] += anoff[i];
    int32_t *ancur = (int32_t *)so_xmalloc(((size_t)nreg + 1)
                                           * sizeof(int32_t));
    memcpy(ancur, anoff, ((size_t)nreg + 1) * sizeof(int32_t));
    int32_t *afaces = (int32_t *)so_xmalloc(
        (anoff[nreg] > 0 ? (size_t)anoff[nreg] : 1) * sizeof(int32_t));
    for (size_t i = 0; i < nadj; i++) {
        int32_t fa = adj[i].fa, fb = adj[i].fb;
        if (region_of[fa] >= 0 && region_of[fb] < 0)
            afaces[ancur[region_of[fa]]++] = fb;
        else if (region_of[fb] >= 0 && region_of[fa] < 0)
            afaces[ancur[region_of[fb]]++] = fa;
    }
    free(ancur);

    /* adjacency by region (CSR over edges whose BOTH endpoints share region) */
    int32_t *eoff = (int32_t *)so_xcalloc((size_t)nreg + 1, sizeof(int32_t));
    for (size_t i = 0; i < nadj; i++) {
        int32_t ra = region_of[adj[i].fa];
        if (ra >= 0 && ra == region_of[adj[i].fb]) eoff[ra + 1]++;
    }
    for (int32_t i = 0; i < nreg; i++) eoff[i + 1] += eoff[i];
    int32_t *ecur = (int32_t *)so_xmalloc(((size_t)nreg + 1) * sizeof(int32_t));
    memcpy(ecur, eoff, ((size_t)nreg + 1) * sizeof(int32_t));
    int32_t *elist = (int32_t *)so_xmalloc(
        (eoff[nreg] > 0 ? (size_t)eoff[nreg] : 1) * sizeof(int32_t));
    for (size_t i = 0; i < nadj; i++) {
        int32_t ra = region_of[adj[i].fa];
        if (ra >= 0 && ra == region_of[adj[i].fb])
            elist[ecur[ra]++] = (int32_t)i;
    }
    free(ecur);

    /* per-region resolution */
    int32_t *g2l = (int32_t *)so_xmalloc(nf * sizeof(int32_t));
    for (size_t f = 0; f < nf; f++) g2l[f] = -1;
    /* global layer key per DROPPED face (rehome groups losers by layer) */
    int32_t *face_lkey = (int32_t *)so_xmalloc(nf * sizeof(int32_t));
    for (size_t f = 0; f < nf; f++) face_lkey[f] = -1;
    int32_t layer_base = 0;
    double du = o->tex_du > 0 ? o->tex_du : 1.0;
    double dv = o->tex_dv > 0 ? o->tex_dv : 1.0;
    double e_sum = 0.0;
    size_t e_n = 0;
    size_t scratch_cap = 0;
    int32_t *loclayer = NULL, *lstack = NULL;

    for (int32_t r = 0; r < nreg; r++) {
        int32_t b0 = roff[r], b1 = roff[r + 1];
        size_t nrf = (size_t)(b1 - b0);
        if (nrf > out->max_region_faces) out->max_region_faces = nrf;
        if (nrf == 0) continue;
        if (nrf > scratch_cap) {
            scratch_cap = nrf * 2;
            loclayer = (int32_t *)so_xrealloc(loclayer,
                                              scratch_cap * sizeof(int32_t));
            lstack = (int32_t *)so_xrealloc(lstack,
                                            scratch_cap * sizeof(int32_t));
        }
        for (int32_t i = b0; i < b1; i++) g2l[rfaces[i]] = i - b0;

        /* layers = CCs of intra-region adjacency */
        SoUF lu;
        so_uf_init(&lu, nrf);
        for (int32_t e = eoff[r]; e < eoff[r + 1]; e++) {
            const QuiltAdjEdge *ae = &adj[elist[e]];
            so_uf_union(&lu, g2l[ae->fa], g2l[ae->fb]);
        }
        int nlayers = 0;
        for (size_t i = 0; i < nrf; i++) loclayer[i] = -1;
        for (size_t i = 0; i < nrf; i++) {
            int32_t root = so_uf_find(&lu, (int32_t)i);
            if (loclayer[root] < 0) loclayer[root] = nlayers++;
        }
        for (size_t i = 0; i < nrf; i++)
            lstack[i] = loclayer[so_uf_find(&lu, (int32_t)i)];

        /* fused fold: a TRUE pair inside one CC -> multicut that CC */
        int fold_cc = -1;
        for (int32_t pi = poff[r]; pi < poff[r + 1]; pi++) {
            const SoPair *pp = &tp.p[plist[pi]];
            int la = lstack[g2l[pp->fa]], lb = lstack[g2l[pp->fb]];
            if (la == lb) { fold_cc = la; break; }
        }
        if (fold_cc >= 0 && nrf >= 2) {
            /* local ids of the fold CC */
            size_t ncc = 0;
            int32_t *ccmap = (int32_t *)so_xmalloc(nrf * sizeof(int32_t));
            for (size_t i = 0; i < nrf; i++) {
                ccmap[i] = -1;
                if (lstack[i] == fold_cc) ccmap[i] = (int32_t)ncc++;
            }
            /* subgraph edges */
            size_t na = 0, no = 0;
            for (int32_t e = eoff[r]; e < eoff[r + 1]; e++) {
                const QuiltAdjEdge *ae = &adj[elist[e]];
                if (ccmap[g2l[ae->fa]] >= 0 && ccmap[g2l[ae->fb]] >= 0) na++;
            }
            for (int32_t pi = poff[r]; pi < poff[r + 1]; pi++) {
                const SoPair *pp = &tp.p[plist[pi]];
                if (ccmap[g2l[pp->fa]] >= 0 && ccmap[g2l[pp->fb]] >= 0) no++;
            }
            if (ncc >= 2 && no > 0) {
                int32_t *af = (int32_t *)so_xmalloc((na + 1) * sizeof(int32_t));
                int32_t *at = (int32_t *)so_xmalloc((na + 1) * sizeof(int32_t));
                int32_t *lf2 = (int32_t *)so_xmalloc((na + no + 1)
                                                     * sizeof(int32_t));
                int32_t *lt2 = (int32_t *)so_xmalloc((na + no + 1)
                                                     * sizeof(int32_t));
                double *lw2 = (double *)so_xmalloc((na + no + 1)
                                                   * sizeof(double));
                size_t qa = 0;
                double avg = 0.0;
                for (int32_t e = eoff[r]; e < eoff[r + 1]; e++) {
                    const QuiltAdjEdge *ae = &adj[elist[e]];
                    int32_t la = ccmap[g2l[ae->fa]], lb = ccmap[g2l[ae->fb]];
                    if (la < 0 || lb < 0) continue;
                    af[qa] = la;
                    at[qa] = lb;
                    lw2[qa] = ae->uvlen;
                    avg += ae->uvlen;
                    qa++;
                }
                avg = qa > 0 ? avg / (double)qa : 1.0;
                if (avg < 1e-9) avg = 1.0;
                for (size_t i = 0; i < qa; i++) {
                    lf2[i] = af[i];
                    lt2[i] = at[i];
                    lw2[i] = lw2[i] / avg;
                }
                size_t qo = qa;
                for (int32_t pi = poff[r]; pi < poff[r + 1]; pi++) {
                    const SoPair *pp = &tp.p[plist[pi]];
                    int32_t la = ccmap[g2l[pp->fa]], lb = ccmap[g2l[pp->fb]];
                    if (la < 0 || lb < 0) continue;
                    lf2[qo] = la;
                    lt2[qo] = lb;
                    lw2[qo] = -1.0e6;
                    qo++;
                }
                int32_t *sub = (int32_t *)so_xmalloc(ncc * sizeof(int32_t));
                int32_t nsub = 0;
                if (LiftedMulticut_kernighan_lin((int32_t)ncc, (int32_t)qa,
                                                 af, at, (int32_t)qo,
                                                 lf2, lt2, lw2, 2, sub,
                                                 &nsub) == 0 && nsub >= 2) {
                    for (size_t i = 0; i < nrf; i++)
                        if (ccmap[i] >= 0 && sub[ccmap[i]] > 0)
                            lstack[i] = nlayers + sub[ccmap[i]] - 1;
                    nlayers += nsub - 1;
                    out->n_multicut_fallbacks++;
                }
                free(sub);
                free(af); free(at); free(lf2); free(lt2); free(lw2);
            }
            free(ccmap);
        }
        so_uf_free(&lu);
        out->n_layers_total += (size_t)nlayers;

        if (nlayers < 2) {   /* single layer: cell aliasing, leave alone */
            for (int32_t i = b0; i < b1; i++) g2l[rfaces[i]] = -1;
            continue;
        }

        /* layer sizes + mean radii */
        int32_t *lsz = (int32_t *)so_xcalloc((size_t)nlayers, sizeof(int32_t));
        double *lr = (double *)so_xcalloc((size_t)nlayers, sizeof(double));
        for (size_t i = 0; i < nrf; i++) {
            lsz[lstack[i]]++;
            lr[lstack[i]] += (double)rf[rfaces[b0 + (int32_t)i]];
        }
        for (int L = 0; L < nlayers; L++)
            if (lsz[L] > 0) lr[L] /= (double)lsz[L];
        int largest = 0;
        for (int L = 1; L < nlayers; L++)
            if (lsz[L] > lsz[largest]) largest = L;
        size_t min_layer = (size_t)(0.15 * (double)nrf);
        if (min_layer < 4) min_layer = 4;

        /* tier V: POSITION-MATCHED anchor votes. Radius varies ACROSS a
         * region (axis eccentricity + in-plane slope over a 128-vox cube),
         * so comparing the anchor ring's median against region-wide layer
         * means picks the wrong wrap when the slope exceeds the gate. Each
         * anchor instead votes for the layer whose uv-NEAREST face (3x3
         * uv-cell neighborhood) best continues the anchor's own radius. */
        int keep = -1;
        int32_t n_anchor = anoff[r + 1] - anoff[r];
        size_t votes_cast = 0;
        if (n_anchor >= 4) {
            int32_t *votes = (int32_t *)so_xcalloc((size_t)nlayers,
                                                   sizeof(int32_t));
            double *lbest = (double *)so_xmalloc((size_t)nlayers
                                                 * sizeof(double));
            int32_t stride = n_anchor > 512 ? n_anchor / 512 : 1;
            for (int32_t ai = 0; ai < n_anchor; ai += stride) {
                int32_t af = afaces[anoff[r] + ai];
                for (int L = 0; L < nlayers; L++) lbest[L] = 1e300;
                size_t acu = (size_t)(((double)fu[af] - umin) / cell);
                size_t acv = (size_t)(((double)fv[af] - vmin) / cell);
                if (acu >= gu) acu = gu - 1;
                if (acv >= gv) acv = gv - 1;
                for (int dyc = -1; dyc <= 1; dyc++) {
                    for (int dxc = -1; dxc <= 1; dxc++) {
                        long ccv = (long)acv + dyc, ccu = (long)acu + dxc;
                        if (ccu < 0 || ccu >= (long)gu
                            || ccv < 0 || ccv >= (long)gv) continue;
                        size_t ci = (size_t)ccv * gu + (size_t)ccu;
                        for (int32_t q = coff[ci]; q < coff[ci + 1]; q++) {
                            int32_t g = clist[q];
                            if (region_of[g] != r) continue;
                            int L = lstack[g2l[g]];
                            if ((size_t)lsz[L] < min_layer) continue;
                            double d = fabs((double)rf[g] - (double)rf[af]);
                            if (d < lbest[L]) lbest[L] = d;
                        }
                    }
                }
                int bl = -1;
                for (int L = 0; L < nlayers; L++) {
                    if (lbest[L] > gate) continue;   /* nothing continuing */
                    if (bl < 0 || lbest[L] < lbest[bl]) bl = L;
                }
                if (bl >= 0) { votes[bl]++; votes_cast++; }
            }
            int w0 = 0;
            for (int L = 1; L < nlayers; L++)
                if (votes[L] > votes[w0]) w0 = L;
            int32_t runner = 0;
            for (int L = 0; L < nlayers; L++)
                if (L != w0 && votes[L] > runner) runner = votes[L];
            if (votes_cast >= 4 && votes[w0] >= 2 * (runner > 0 ? runner : 1)
                && votes[w0] > 0) {
                keep = w0;
                out->n_vote_picks++;
            }
            if (o->verbose >= 2) {
                fprintf(stderr, "[seam_own dbg] region %d: nrf=%zu nlayers=%d "
                        "n_anchor=%d votes_cast=%zu keep=%d\n",
                        r, nrf, nlayers, n_anchor, votes_cast, keep);
                for (int L = 0; L < nlayers; L++)
                    fprintf(stderr, "[seam_own dbg]   layer %d: sz=%d r=%.1f "
                            "votes=%d%s\n", L, lsz[L], lr[L], votes[L],
                            (size_t)lsz[L] < min_layer ? " (sliver)" : "");
            }
            free(votes);
            free(lbest);
        }

        /* tier E: boundary energy over substantial layers */
        int used_energy = 0;
        if (keep < 0 && have_ct && nrf <= o->region_cap) {
            double gu0 = 1e300, gv0 = 1e300, gu1 = -1e300, gv1 = -1e300;
            for (int32_t i = b0; i < b1; i++) {
                int32_t f = rfaces[i];
                for (int k = 0; k < 3; k++) {
                    size_t vi = (size_t)ps->faces[f * 3 + k];
                    double u = ps->uv[vi * 2 + 0], v = ps->uv[vi * 2 + 1];
                    if (u < gu0) gu0 = u;
                    if (u > gu1) gu1 = u;
                    if (v < gv0) gv0 = v;
                    if (v > gv1) gv1 = v;
                }
            }
            for (int32_t i = anoff[r]; i < anoff[r + 1]; i++) {
                int32_t f = afaces[i];
                for (int k = 0; k < 3; k++) {
                    size_t vi = (size_t)ps->faces[f * 3 + k];
                    double u = ps->uv[vi * 2 + 0], v = ps->uv[vi * 2 + 1];
                    if (u < gu0) gu0 = u;
                    if (u > gu1) gu1 = u;
                    if (v < gv0) gv0 = v;
                    if (v > gv1) gv1 = v;
                }
            }
            gu0 -= 2 * du; gv0 -= 2 * dv;
            gu1 += 2 * du; gv1 += 2 * dv;
            int gw = (int)((gu1 - gu0) / du) + 1;
            int gh = (int)((gv1 - gv0) / dv) + 1;
            if (gw >= 3 && gh >= 3 && gw <= o->grid_cap && gh <= o->grid_cap
                && n_anchor > 0) {
                size_t np = (size_t)gw * (size_t)gh;
                float *imgA = (float *)so_xmalloc(np * sizeof(float));
                float *imgL = (float *)so_xmalloc(np * sizeof(float));
                uint8_t *covA = (uint8_t *)so_xmalloc(np);
                uint8_t *covL = (uint8_t *)so_xmalloc(np);
                Quilt_raster_faces(&afaces[anoff[r]], (size_t)n_anchor,
                                   ps->verts, ps->uv, nrm, ps->faces, ct,
                                   o->normal_range, o->normal_samples,
                                   gu0, gv0, du, dv, gw, gh, imgA, covA);
                int32_t *lfb = (int32_t *)so_xmalloc(nrf * sizeof(int32_t));
                double bestE = 1e19;
                for (int L = 0; L < nlayers; L++) {
                    if ((size_t)lsz[L] < min_layer) continue;
                    size_t nl = 0;
                    for (size_t i = 0; i < nrf; i++)
                        if (lstack[i] == L)
                            lfb[nl++] = rfaces[b0 + (int32_t)i];
                    Quilt_raster_faces(lfb, nl, ps->verts, ps->uv, nrm,
                                       ps->faces, ct, o->normal_range,
                                       o->normal_samples, gu0, gv0, du, dv,
                                       gw, gh, imgL, covL);
                    double E = Quilt_boundary_energy(imgL, covL, imgA, covA,
                                                     gw, gh);
                    if (E < bestE) { bestE = E; keep = L; }
                }
                free(lfb);
                free(imgA); free(imgL); free(covA); free(covL);
                if (keep >= 0 && bestE < 1e17) {
                    e_sum += bestE;
                    e_n++;
                    out->n_energy_picks++;
                    used_energy = 1;
                } else {
                    keep = -1;
                }
            } else if (gw > o->grid_cap || gh > o->grid_cap) {
                out->n_regions_skipped_grid++;
            }
        }
        if (keep < 0 && nrf > o->region_cap)
            out->n_regions_skipped_cap++;   /* energy skipped, seed = largest */
        if (keep < 0) {
            keep = largest;
            if (!used_energy) out->n_largest_picks++;
        }

        /* CONFLICT PROPAGATION. A region is the TRANSITIVE closure of pair
         * overlaps: in a chain A~B~C, layer C never overlapped A, and a
         * keep-ONE rule would drop it just for sharing A's region. Instead:
         * seed the kept set with the picked layer, then admit every layer
         * (largest first) that has NO direct pair-conflict with an already-
         * kept layer. Drops stay bounded to genuine overlaps of kept
         * material; non-conflicting chain members and same-sheet slivers
         * survive automatically. */
        {
            /* layer-conflict edges from this region's true pairs (dedup'd) */
            int32_t np0 = poff[r], np1 = poff[r + 1];
            size_t nce = 0;
            uint64_t *ce = (uint64_t *)so_xmalloc(
                ((size_t)(np1 - np0) + 1) * sizeof(uint64_t));
            for (int32_t pi = np0; pi < np1; pi++) {
                const SoPair *pp = &tp.p[plist[pi]];
                int la = lstack[g2l[pp->fa]], lb = lstack[g2l[pp->fb]];
                if (la == lb) continue;   /* unsplit fold remainder */
                int lo = la < lb ? la : lb, hi = la < lb ? lb : la;
                ce[nce++] = ((uint64_t)(uint32_t)lo << 32) | (uint32_t)hi;
            }
            qsort(ce, nce, sizeof(uint64_t), so_cmp_u64);
            size_t nu = 0;
            for (size_t i = 0; i < nce; i++)
                if (nu == 0 || ce[i] != ce[nu - 1]) ce[nu++] = ce[i];

            /* layer -> conflict-neighbor CSR (both directions) */
            int32_t *loff2 = (int32_t *)so_xcalloc((size_t)nlayers + 1,
                                                   sizeof(int32_t));
            for (size_t i = 0; i < nu; i++) {
                loff2[(int32_t)(ce[i] >> 32) + 1]++;
                loff2[(int32_t)(uint32_t)ce[i] + 1]++;
            }
            for (int L = 0; L < nlayers; L++) loff2[L + 1] += loff2[L];
            int32_t *lnbr = (int32_t *)so_xmalloc(
                ((size_t)loff2[nlayers] + 1) * sizeof(int32_t));
            {
                int32_t *cur = (int32_t *)so_xmalloc(
                    ((size_t)nlayers + 1) * sizeof(int32_t));
                memcpy(cur, loff2, ((size_t)nlayers + 1) * sizeof(int32_t));
                for (size_t i = 0; i < nu; i++) {
                    int32_t la = (int32_t)(ce[i] >> 32);
                    int32_t lb = (int32_t)(uint32_t)ce[i];
                    lnbr[cur[la]++] = lb;
                    lnbr[cur[lb]++] = la;
                }
                free(cur);
            }
            free(ce);

            /* layers by size desc (stable enough via qsort on packed key) */
            uint64_t *order = (uint64_t *)so_xmalloc((size_t)nlayers
                                                     * sizeof(uint64_t));
            for (int L = 0; L < nlayers; L++)
                order[L] = ((uint64_t)(uint32_t)(INT32_MAX - lsz[L]) << 32)
                         | (uint32_t)L;
            qsort(order, (size_t)nlayers, sizeof(uint64_t), so_cmp_u64);

            uint8_t *lkeep = (uint8_t *)so_xcalloc((size_t)nlayers, 1);
            lkeep[keep] = 1;
            for (int oi = 0; oi < nlayers; oi++) {
                int L = (int)(uint32_t)order[oi];
                if (L == keep) continue;
                int clash = 0;
                for (int32_t q = loff2[L]; q < loff2[L + 1] && !clash; q++)
                    if (lkeep[lnbr[q]]) clash = 1;
                if (!clash) lkeep[L] = 1;
            }
            free(order);
            free(loff2);
            free(lnbr);

            for (size_t i = 0; i < nrf; i++) {
                int32_t gf = rfaces[b0 + (int32_t)i];
                if (lkeep[lstack[i]]) {
                    out->face_dec[gf] = 1;
                    out->n_kept++;
                } else {
                    out->face_dec[gf] = 3;
                    out->face_keep[gf] = 0;
                    out->n_dropped++;
                    face_lkey[gf] = layer_base + lstack[i];
                }
            }
            free(lkeep);
        }
        layer_base += nlayers;
        free(lsz);
        free(lr);
        for (int32_t i = b0; i < b1; i++) g2l[rfaces[i]] = -1;
    }
    free(loclayer);
    free(lstack);
    out->energy_mean = e_n ? e_sum / (double)e_n : 0.0;

    /* ---- SEAM pass (kept faces only) -------------------------------------- */
    if (o->seam_cut && sp.n > 0 && have_ct) {
        SoUF su;
        so_uf_init(&su, nf);
        size_t live_pairs = 0;
        for (size_t i = 0; i < sp.n; i++) {
            int32_t fa = sp.p[i].fa, fb = sp.p[i].fb;
            if (!out->face_keep[fa] || !out->face_keep[fb]) continue;
            so_uf_union(&su, fa, fb);
            live_pairs++;
        }
        for (size_t i = 0; i < nadj; i++) {
            int32_t fa = adj[i].fa, fb = adj[i].fb;
            if (is_seam[fa] && is_seam[fb]
                && out->face_keep[fa] && out->face_keep[fb])
                so_uf_union(&su, fa, fb);
        }
        if (live_pairs > 0) {
            int32_t *sreg = (int32_t *)so_xmalloc(nf * sizeof(int32_t));
            int32_t *sroot = (int32_t *)so_xmalloc(nf * sizeof(int32_t));
            for (size_t f = 0; f < nf; f++) { sreg[f] = -1; sroot[f] = -1; }
            int32_t nsr = 0;
            for (size_t f = 0; f < nf; f++) {
                if (!is_seam[f] || !out->face_keep[f]) continue;
                int32_t rt = so_uf_find(&su, (int32_t)f);
                if (sroot[rt] < 0) sroot[rt] = nsr++;
                sreg[f] = sroot[rt];
            }
            out->n_seam_regions = (size_t)nsr;

            int32_t *soff = (int32_t *)so_xcalloc((size_t)nsr + 1,
                                                  sizeof(int32_t));
            for (size_t f = 0; f < nf; f++)
                if (sreg[f] >= 0) soff[sreg[f] + 1]++;
            for (int32_t i = 0; i < nsr; i++) soff[i + 1] += soff[i];
            int32_t *scur = (int32_t *)so_xmalloc(((size_t)nsr + 1)
                                                  * sizeof(int32_t));
            memcpy(scur, soff, ((size_t)nsr + 1) * sizeof(int32_t));
            int32_t *sfaces = (int32_t *)so_xmalloc(
                (soff[nsr] > 0 ? (size_t)soff[nsr] : 1) * sizeof(int32_t));
            for (size_t f = 0; f < nf; f++)
                if (sreg[f] >= 0) sfaces[scur[sreg[f]]++] = (int32_t)f;
            free(scur);

            for (int32_t r = 0; r < nsr; r++) {
                int32_t b0 = soff[r], b1 = soff[r + 1];
                size_t nrf = (size_t)(b1 - b0);
                if (nrf < 2) continue;
                /* top-2 cubes by face count */
                int32_t ca = -1, cb = -1;
                size_t na = 0, nb = 0;
                for (int32_t i = b0; i < b1; i++) {
                    int32_t cube = ps->face_cube[sfaces[i]];
                    if (cube == ca) na++;
                    else if (cube == cb) nb++;
                    else if (ca < 0) { ca = cube; na = 1; }
                    else if (cb < 0) { cb = cube; nb = 1; }
                    /* other cubes: counted approximately; two dominate seams */
                }
                if (ca < 0 || cb < 0) continue;
                if (nb > na) { int32_t t = ca; ca = cb; cb = t; }
                SoSeamStat sst;
                so_seam_region(ps, o, ct, nrm, &sfaces[b0], nrf, ca, cb,
                               out->face_keep, out->face_dec, &sst);
                out->n_seam_dropped += sst.dropped;
                out->n_seam_restored += sst.restored;
                if (sst.abut_only) out->n_seams_abut_only++;
                if (sst.skipped) out->n_seams_skipped++;
            }
            free(sreg); free(sroot); free(soff); free(sfaces);
        }
        so_uf_free(&su);
    } else if (o->seam_cut && sp.n > 0) {
        out->n_seams_skipped = sp.n;   /* no RAW: value-domain DP impossible */
    }

    /* ---- REHOME pass: find homes for dropped layers ------------------------
     * A losing layer whose radius sits ~k*pitch from the winner at the same
     * (u,v) is the SAME physical sheet mis-registered by k whole turns.
     * Where it is strictly vertex-disjoint from kept faces (relocating a
     * mesh-attached layer would tear it -- those joins belong to the
     * front-half re-registration), move it to its correct turn with the
     * exact spiral map: u += DeltaU(a,b,phi,k), phi += 2pi k, then keep it.
     * Destinations already covered at the same radius are genuine
     * duplicates (stay dropped); destinations holding OTHER-wrap material
     * would mint a fresh conflict (blocked; stay dropped). */
    size_t rehomed_any = 0;
    if (o->rehome && ps->phi != NULL && o->spiral_b != 0.0
        && out->n_dropped > 0) {
        uint8_t *vused = (uint8_t *)so_xcalloc(nv, 1);
        for (size_t f = 0; f < nf; f++) {
            if (!out->face_keep[f]) continue;
            for (int k = 0; k < 3; k++)
                vused[(size_t)ps->faces[f * 3 + k]] = 1;
        }
        /* dropped faces by layer key (CSR over 0..layer_base) */
        int32_t *lcnt = (int32_t *)so_xcalloc((size_t)layer_base + 1,
                                              sizeof(int32_t));
        for (size_t f = 0; f < nf; f++)
            if (face_lkey[f] >= 0) lcnt[face_lkey[f]]++;
        int32_t *loff3 = (int32_t *)so_xcalloc((size_t)layer_base + 2,
                                               sizeof(int32_t));
        for (int32_t k2 = 0; k2 < layer_base; k2++)
            loff3[k2 + 1] = loff3[k2] + lcnt[k2];
        int32_t *lfaces = (int32_t *)so_xmalloc(
            ((size_t)loff3[layer_base] + 1) * sizeof(int32_t));
        {
            int32_t *cur = (int32_t *)so_xmalloc(((size_t)layer_base + 1)
                                                 * sizeof(int32_t));
            memcpy(cur, loff3, (size_t)layer_base * sizeof(int32_t));
            for (size_t f = 0; f < nf; f++)
                if (face_lkey[f] >= 0)
                    lfaces[cur[face_lkey[f]]++] = (int32_t)f;
            free(cur);
        }
        /* layers by size desc (bigger layers claim destinations first) */
        int32_t *lorder = (int32_t *)so_xmalloc(((size_t)layer_base + 1)
                                                * sizeof(int32_t));
        int32_t nlay = 0;
        for (int32_t k2 = 0; k2 < layer_base; k2++)
            if (lcnt[k2] > 0) lorder[nlay++] = k2;
        for (int32_t a2 = 1; a2 < nlay; a2++) {   /* insertion by count desc */
            int32_t key = lorder[a2], b2 = a2 - 1;
            while (b2 >= 0 && lcnt[lorder[b2]] < lcnt[key]) {
                lorder[b2 + 1] = lorder[b2];
                b2--;
            }
            lorder[b2 + 1] = key;
        }
        /* destination-claim map: cells taken by committed rehomes */
        uint8_t *cellblk = (uint8_t *)so_xcalloc(ncell, 1);
        double *drbuf = (double *)so_xmalloc(4096 * sizeof(double));

        for (int32_t li = 0; li < nlay; li++) {
            int32_t key = lorder[li];
            int32_t f0 = loff3[key], f1 = loff3[key + 1];
            int32_t nlf = f1 - f0;
            if (nlf < o->rehome_min_faces) {
                out->n_rehome_small++;
                continue;
            }
            /* tear guard: every vert unused by kept faces */
            int disjoint = 1;
            for (int32_t q = f0; q < f1 && disjoint; q++) {
                int32_t f = lfaces[q];
                for (int k = 0; k < 3; k++)
                    if (vused[(size_t)ps->faces[f * 3 + k]]) {
                        disjoint = 0;
                        break;
                    }
            }
            if (!disjoint) {
                out->n_rehome_adjacent++;
                continue;
            }
            /* position-matched radius offset vs nearest kept face */
            int32_t ndr = 0;
            int32_t stride = nlf > 4096 ? nlf / 4096 : 1;
            for (int32_t q = f0; q < f1 && ndr < 4096; q += stride) {
                int32_t f = lfaces[q];
                size_t acu = (size_t)(((double)fu[f] - umin) / cell);
                size_t acv = (size_t)(((double)fv[f] - vmin) / cell);
                if (acu >= gu) acu = gu - 1;
                if (acv >= gv) acv = gv - 1;
                double bestd = 1e300, bestdr = 0.0;
                for (int dyc = -1; dyc <= 1; dyc++) {
                    for (int dxc = -1; dxc <= 1; dxc++) {
                        long ccv = (long)acv + dyc, ccu = (long)acu + dxc;
                        if (ccu < 0 || ccu >= (long)gu
                            || ccv < 0 || ccv >= (long)gv) continue;
                        size_t ci = (size_t)ccv * gu + (size_t)ccu;
                        for (int32_t w2 = coff[ci]; w2 < coff[ci + 1]; w2++) {
                            int32_t g = clist[w2];
                            if (!out->face_keep[g]) continue;
                            double duu = (double)fu[g] - (double)fu[f];
                            double dvv = (double)fv[g] - (double)fv[f];
                            double d2 = duu * duu + dvv * dvv;
                            if (d2 < bestd) {
                                bestd = d2;
                                bestdr = (double)rf[f] - (double)rf[g];
                            }
                        }
                    }
                }
                if (bestd < 1e300)
                    drbuf[ndr++] = bestdr;
            }
            if (ndr < 4) {
                out->n_rehome_incoherent++;
                continue;
            }
            /* median + mad */
            for (int32_t a2 = 1; a2 < ndr; a2++) {
                double key2 = drbuf[a2];
                int32_t b2 = a2 - 1;
                while (b2 >= 0 && drbuf[b2] > key2) {
                    drbuf[b2 + 1] = drbuf[b2];
                    b2--;
                }
                drbuf[b2 + 1] = key2;
            }
            double med = drbuf[ndr / 2];
            double madr = 0.0;
            {
                double acc = 0.0;
                for (int32_t a2 = 0; a2 < ndr; a2++)
                    acc += fabs(drbuf[a2] - med);
                madr = acc / (double)ndr;   /* mean abs dev: cheap, adequate */
            }
            double kf = med / o->spiral_b;
            int32_t k = (int32_t)lround(kf);
            if (k == 0 || fabs(kf - (double)k) > o->rehome_frac_tol
                || madr > 0.35 * (o->pitch > 0 ? o->pitch : 9.5)) {
                out->n_rehome_incoherent++;
                continue;
            }
            /* destination test: sample faces, look at target cells */
            int32_t n_cov = 0, n_blk = 0, n_test = 0;
            for (int32_t q = f0; q < f1; q += stride) {
                int32_t f = lfaces[q];
                double phif = ((double)ps->phi[(size_t)ps->faces[f * 3 + 0]]
                             + (double)ps->phi[(size_t)ps->faces[f * 3 + 1]]
                             + (double)ps->phi[(size_t)ps->faces[f * 3 + 2]])
                            / 3.0;
                double u2 = (double)fu[f]
                          + CubeReg_deltaU(o->spiral_a, o->spiral_b, phif, k);
                long ccu = (long)((u2 - umin) / cell);
                long ccv = (long)(((double)fv[f] - vmin) / cell);
                n_test++;
                if (ccu < 0 || ccu >= (long)gu || ccv < 0
                    || ccv >= (long)gv)
                    continue;   /* off-canvas: treated as empty */
                int cov = 0, blk = 0;
                for (int dxc = -1; dxc <= 1 && !blk; dxc++) {
                    long cu2 = ccu + dxc;
                    if (cu2 < 0 || cu2 >= (long)gu) continue;
                    size_t ci = (size_t)ccv * gu + (size_t)cu2;
                    if (cellblk[ci]) { blk = 1; break; }
                    for (int32_t w2 = coff[ci]; w2 < coff[ci + 1]; w2++) {
                        int32_t g = clist[w2];
                        if (!out->face_keep[g]) continue;
                        double dr = fabs((double)rf[f] - (double)rf[g]);
                        if (dr <= gate) cov = 1;
                        else { blk = 1; break; }
                    }
                }
                if (blk) n_blk++;
                else if (cov) n_cov++;
            }
            if (n_blk > n_test / 10) {
                out->n_rehome_blocked++;
                continue;
            }
            if (n_cov * 2 >= n_test) {
                out->n_rehome_dup++;
                continue;
            }
            /* COMMIT: exact per-vertex relocation (dedupe verts via mark) */
            for (int32_t q = f0; q < f1; q++) {
                int32_t f = lfaces[q];
                for (int k3 = 0; k3 < 3; k3++) {
                    size_t vi = (size_t)ps->faces[f * 3 + k3];
                    if (vused[vi]) continue;   /* reuse as "moved" mark */
                    double phv = (double)ps->phi[vi];
                    ps->uv[vi * 2 + 0] = (float)((double)ps->uv[vi * 2 + 0]
                        + CubeReg_deltaU(o->spiral_a, o->spiral_b, phv, k));
                    ps->phi[vi] = (float)(phv + SO_TWO_PI * (double)k);
                    vused[vi] = 1;
                }
                out->face_keep[f] = 1;
                out->face_dec[f] = 6;
                out->n_rehomed_faces++;
                out->n_dropped--;
                /* claim the destination cell */
                double phif = ((double)ps->phi[(size_t)ps->faces[f * 3 + 0]]
                             + (double)ps->phi[(size_t)ps->faces[f * 3 + 1]]
                             + (double)ps->phi[(size_t)ps->faces[f * 3 + 2]])
                            / 3.0;
                (void)phif;
                double u2 = ((double)ps->uv[(size_t)ps->faces[f * 3 + 0] * 2]
                           + (double)ps->uv[(size_t)ps->faces[f * 3 + 1] * 2]
                           + (double)ps->uv[(size_t)ps->faces[f * 3 + 2] * 2])
                          / 3.0;
                fu[f] = (float)u2;
                long ccu = (long)((u2 - umin) / cell);
                long ccv = (long)(((double)fv[f] - vmin) / cell);
                if (ccu >= 0 && ccu < (long)gu && ccv >= 0
                    && ccv < (long)gv)
                    cellblk[(size_t)ccv * gu + (size_t)ccu] = 1;
            }
            out->n_rehome_layers++;
            rehomed_any = 1;
        }
        free(vused);
        free(lcnt);
        free(loff3);
        free(lfaces);
        free(lorder);
        free(cellblk);
        free(drbuf);

        /* moved faces changed cells: rebuild the cell CSR so the recount
         * below sees the joined positions */
        if (rehomed_any) {
            free(coff);
            free(clist);
            coff = (int32_t *)so_xcalloc(ncell + 1, sizeof(int32_t));
            int32_t *fcell2 = (int32_t *)so_xmalloc(nf * sizeof(int32_t));
            for (size_t f = 0; f < nf; f++) {
                size_t cu2 = (size_t)(((double)fu[f] - umin) / cell);
                size_t cv2 = (size_t)(((double)fv[f] - vmin) / cell);
                if (cu2 >= gu) cu2 = gu - 1;
                if (cv2 >= gv) cv2 = gv - 1;
                fcell2[f] = (int32_t)(cv2 * gu + cu2);
                coff[fcell2[f] + 1]++;
            }
            for (size_t i = 0; i < ncell; i++) coff[i + 1] += coff[i];
            int32_t *ccur2 = (int32_t *)so_xmalloc(ncell * sizeof(int32_t));
            memcpy(ccur2, coff, ncell * sizeof(int32_t));
            clist = (int32_t *)so_xmalloc(nf * sizeof(int32_t));
            for (size_t f = 0; f < nf; f++)
                clist[ccur2[fcell2[f]]++] = (int32_t)f;
            free(ccur2);
            free(fcell2);
        }
    }

    /* ---- multi-cover cells AFTER (kept faces, same classification) ------- */
    size_t multi_after = 0;
    for (size_t c = 0; c < ncell; c++) {
        int32_t a0 = coff[c], a1 = coff[c + 1];
        int hm = 0;
        for (int32_t i = a0; i < a1 && !hm; i++) {
            int32_t fa = clist[i];
            if (!out->face_keep[fa]) continue;
            for (int32_t j = i + 1; j < a1; j++) {
                int32_t fb = clist[j];
                if (!out->face_keep[fb]) continue;
                double dr = fabs((double)rf[fa] - (double)rf[fb]);
                if (ps->face_cube[fa] == ps->face_cube[fb]) {
                    if (dr <= gate) continue;
                    /* mirror the before-count's adjacency exemption */
                    int32_t lo = fa < fb ? fa : fb, hi = fa < fb ? fb : fa;
                    uint64_t key = ((uint64_t)(uint32_t)lo << 32)
                                 | (uint32_t)hi;
                    size_t los = 0, his = nadj;
                    while (los < his) {
                        size_t m = (los + his) / 2;
                        if (akeys[m] < key) los = m + 1;
                        else his = m;
                    }
                    if (los < nadj && akeys[los] == key) continue;
                    hm = 1;
                    break;
                } else {
                    /* cross-cube adjacency exemption (welded surfaces) */
                    int32_t lo = fa < fb ? fa : fb, hi = fa < fb ? fb : fa;
                    uint64_t key = ((uint64_t)(uint32_t)lo << 32)
                                 | (uint32_t)hi;
                    size_t los = 0, his = nadj;
                    while (los < his) {
                        size_t m = (los + his) / 2;
                        if (akeys[m] < key) los = m + 1;
                        else his = m;
                    }
                    if (los < nadj && akeys[los] == key) continue;
                    if (dr > gate) {
                        hm = 1;
                        break;
                    }
                    double d2 = 0.0;
                    for (int d = 0; d < 3; d++) {
                        double dd = (double)fc3[fa * 3 + d]
                                  - (double)fc3[fb * 3 + d];
                        d2 += dd * dd;
                    }
                    if (d2 <= gate3d2) { hm = 1; break; }
                }
            }
        }
        if (hm) multi_after++;
    }
    out->multi_cells_after = multi_after;

    if (o->verbose) {
        fprintf(stderr,
                "[seam_own] pairs true=%zu seam=%zu mystery=%zu | regions=%zu "
                "layers=%zu max=%zu picks v/e/l=%zu/%zu/%zu mc=%zu "
                "skip cap/grid=%zu/%zu | kept=%zu dropped=%zu | seams=%zu "
                "abut=%zu skip=%zu drop=%zu restore=%zu | cells %zu->%zu "
                "E=%.2f tex=%d gate=%.2f\n",
                out->n_true_pairs, out->n_seam_pairs, out->n_mystery_pairs,
                out->n_regions, out->n_layers_total, out->max_region_faces,
                out->n_vote_picks, out->n_energy_picks, out->n_largest_picks,
                out->n_multicut_fallbacks, out->n_regions_skipped_cap,
                out->n_regions_skipped_grid, out->n_kept, out->n_dropped,
                out->n_seam_regions, out->n_seams_abut_only,
                out->n_seams_skipped, out->n_seam_dropped,
                out->n_seam_restored, out->multi_cells_before,
                out->multi_cells_after, out->energy_mean, have_ct, gate);
    }

    free(rf); free(fu); free(fv); free(fc3);
    free(coff); free(clist);
    free(adj); free(akeys);
    free(tp.p); free(sp.p);
    free(is_true); free(is_seam);
    so_uf_free(&uf);
    free(region_of); free(root2rid);
    free(roff); free(rfaces);
    free(poff); free(plist);
    free(anoff); free(afaces);
    free(eoff); free(elist);
    free(g2l);
    free(face_lkey);
    out->seconds = ves_clock_sec() - t_start;
    return 0;
}

/* ============================================================================
 * Self-test.
 * ==========================================================================*/

static int so_check(int cond, const char *what, int *fails)
{
    if (!cond) {
        fprintf(stderr, "[seam_own selftest]   FAIL: %s\n", what);
        (*fails)++;
    }
    return cond;
}

/* append an n x n grid of quads (2 tris each) lying in the y=Y plane,
 * spanning x=[x0,x0+2n), z=[0,2n), with uv = (u0 + (x-x0), z). Returns the
 * face count added. Vertices/faces/uv/gid/face_cube appended at *pv/*pf. */
static void st_grid(float y, float x0, float u0, int n, int32_t cube,
                    float *verts, float *uv, int32_t *gid,
                    int32_t *faces, int32_t *face_cube,
                    size_t *pv, size_t *pf)
{
    size_t v0 = *pv;
    int side = n + 1;
    for (int iz = 0; iz < side; iz++) {
        for (int ix = 0; ix < side; ix++) {
            size_t vi = v0 + (size_t)iz * (size_t)side + (size_t)ix;
            verts[vi * 3 + 0] = (float)(2 * iz);        /* z */
            verts[vi * 3 + 1] = y;                      /* y */
            verts[vi * 3 + 2] = x0 + (float)(2 * ix);   /* x */
            uv[vi * 2 + 0] = u0 + (float)(2 * ix);      /* u tracks x */
            uv[vi * 2 + 1] = (float)(2 * iz);           /* v tracks z */
            gid[vi] = 0;
        }
    }
    for (int iz = 0; iz < n; iz++) {
        for (int ix = 0; ix < n; ix++) {
            size_t a = v0 + (size_t)iz * (size_t)side + (size_t)ix;
            size_t b = a + 1;
            size_t c = a + (size_t)side;
            size_t d = c + 1;
            faces[*pf * 3 + 0] = (int32_t)a;
            faces[*pf * 3 + 1] = (int32_t)b;
            faces[*pf * 3 + 2] = (int32_t)c;
            face_cube[*pf] = cube;
            (*pf)++;
            faces[*pf * 3 + 0] = (int32_t)b;
            faces[*pf * 3 + 1] = (int32_t)d;
            faces[*pf * 3 + 2] = (int32_t)c;
            face_cube[*pf] = cube;
            (*pf)++;
        }
    }
    *pv += (size_t)side * (size_t)side;
}

/* wire a PieceSet around st_grid-built arrays (2 cubes max) */
static void st_wire(PieceSet *ps, float *verts, float *uv, float *normals,
                    int32_t *gid, int32_t *faces, int32_t *face_cube,
                    size_t nv, size_t nf, size_t *cube_voff, long (*org)[3],
                    size_t n_cubes)
{
    memset(ps, 0, sizeof(*ps));
    ps->verts = verts;
    ps->uv = uv;
    ps->normals = normals;
    ps->gid = gid;
    ps->faces = faces;
    ps->face_cube = face_cube;
    ps->nv = nv;
    ps->nf = nf;
    ps->cube_voff = cube_voff;
    ps->cube_org = org;
    ps->n_cubes = n_cubes;
    for (size_t i = 0; i < nv; i++) {
        normals[i * 3 + 0] = 0.0f;
        normals[i * 3 + 1] = 1.0f;
        normals[i * 3 + 2] = 0.0f;
    }
}

#define ST_RAW_DIR "output/_selftest_seam_own/raw"

static int st_write_raw(void)
{
    static int written = 0;
    if (written) return 0;
    long chunk = 64;
    uint8_t *vol = (uint8_t *)so_xmalloc((size_t)(chunk * chunk * chunk));
    for (long z = 0; z < chunk; z++)
        for (long y = 0; y < chunk; y++)
            for (long x = 0; x < chunk; x++)
                vol[(z * chunk + y) * chunk + x] =
                    (uint8_t)(40 + (x * 3) % 160);
    char path[512];
    snprintf(path, sizeof(path), "%s/z00000_y00000_x00000.tif", ST_RAW_DIR);
    ves_ensure_parent_dir(path);
    int rc = TiffIO_save(path, vol, (int)chunk, (int)chunk, (int)chunk);
    free(vol);
    written = (rc == 0);
    return rc;
}

int SeamOwn_selftest(void)
{
    int fails = 0;
    Arena_T arena = Arena_new();
    st_write_raw();

    enum { GN = 10, VCAP = 2048, FCAP = 2048 };
    float *verts = (float *)so_xmalloc(VCAP * 3 * sizeof(float));
    float *uv = (float *)so_xmalloc(VCAP * 2 * sizeof(float));
    float *normals = (float *)so_xmalloc(VCAP * 3 * sizeof(float));
    int32_t *gid = (int32_t *)so_xmalloc(VCAP * sizeof(int32_t));
    int32_t *faces = (int32_t *)so_xmalloc(FCAP * 3 * sizeof(int32_t));
    int32_t *face_cube = (int32_t *)so_xmalloc(FCAP * sizeof(int32_t));
    size_t cube_voff[3] = { 0, 0, 0 };
    long org[2][3] = { { 0, 0, 0 }, { 0, 0, 0 } };

    SeamOwnOpts o;
    SeamOwnOpts_default(&o);
    o.axis_point[0] = 0.0f;
    o.axis_point[1] = 0.0f;
    o.axis_point[2] = 0.0f;
    o.pitch = 10.0;          /* gate = 5 */
    o.raw_dir = ST_RAW_DIR;
    o.raw_chunk = 64;
    o.rehome = 0;            /* t1-t6 test the drop machinery; t7+ rehome */
    o.verbose = 0;

    /* t1: SEAM double-paint -- two cubes, same 3D slab overlap in u.
     * cube A x=[0,20] u=[0,20]; cube B x=[16,36] u=[16,36]: both paint the
     * SAME 3D points on x=[16,20] (u overlap band 16..20). */
    {
        size_t nv = 0, nf = 0;
        st_grid(20.0f, 0.0f, 0.0f, GN, 0, verts, uv, gid, faces, face_cube,
                &nv, &nf);
        cube_voff[1] = nv;
        st_grid(20.0f, 16.0f, 16.0f, GN, 1, verts, uv, gid, faces, face_cube,
                &nv, &nf);
        cube_voff[2] = nv;
        PieceSet ps;
        st_wire(&ps, verts, uv, normals, gid, faces, face_cube, nv, nf,
                cube_voff, org, 2);
        SeamOwnResult r;
        int rc = SeamOwn_run(arena, &ps, &o, &r);
        so_check(rc == 0, "t1 rc", &fails);
        so_check(r.n_seam_pairs > 0, "t1 seam pairs detected", &fails);
        so_check(r.n_true_pairs == 0, "t1 no true pairs (same radius)",
                 &fails);
        so_check(r.n_seam_regions == 1, "t1 one seam region", &fails);
        so_check(r.n_seam_dropped > 0, "t1 seam DP dropped faces", &fails);
        size_t keptA = 0, keptB = 0;
        for (size_t f = 0; f < nf; f++) {
            if (!r.face_keep[f]) continue;
            if (face_cube[f] == 0) keptA++;
            else keptB++;
        }
        so_check(keptA > 0 && keptB > 0, "t1 both cubes retain faces",
                 &fails);
        so_check(r.multi_cells_after < r.multi_cells_before,
                 "t1 multi cells drop", &fails);
        fprintf(stderr, "[seam_own selftest] t1 seam: pairs=%zu region=%zu "
                "dropped=%zu restored=%zu cells %zu->%zu\n",
                r.n_seam_pairs, r.n_seam_regions, r.n_seam_dropped,
                r.n_seam_restored, r.multi_cells_before, r.multi_cells_after);
    }

    /* t2: TRUE overlap -- cube A slab y=20 spans u=[0,40] (its right half is
     * the clean fringe = anchors), cube B slab y=32 (dr=12 > gate) maps ONTO
     * u=[0,20]: one region, layers split by CC, anchor vote keeps A, B's
     * faces dropped. */
    {
        size_t nv = 0, nf = 0;
        st_grid(20.0f, 0.0f, 0.0f, 2 * GN, 0, verts, uv, gid, faces,
                face_cube, &nv, &nf);
        cube_voff[1] = nv;
        size_t fB0 = nf;
        st_grid(32.0f, 0.0f, 0.0f, GN, 1, verts, uv, gid, faces, face_cube,
                &nv, &nf);
        cube_voff[2] = nv;
        PieceSet ps;
        st_wire(&ps, verts, uv, normals, gid, faces, face_cube, nv, nf,
                cube_voff, org, 2);
        SeamOwnResult r;
        int rc = SeamOwn_run(arena, &ps, &o, &r);
        so_check(rc == 0, "t2 rc", &fails);
        so_check(r.n_true_pairs > 0, "t2 true pairs", &fails);
        so_check(r.n_regions >= 1, "t2 regions", &fails);
        so_check(r.n_dropped > 0, "t2 losing layer dropped", &fails);
        so_check(r.n_vote_picks + r.n_energy_picks >= 1,
                 "t2 vote or energy pick", &fails);
        /* the dropped side must be ONE cube's overlap faces, not both */
        size_t dropA = 0, dropB = 0;
        for (size_t f = 0; f < nf; f++) {
            if (r.face_keep[f]) continue;
            if (f < fB0) dropA++;
            else dropB++;
        }
        so_check((dropA == 0) != (dropB == 0), "t2 one side dropped",
                 &fails);
        fprintf(stderr, "[seam_own selftest] t2 true: pairs=%zu regions=%zu "
                "picks v/e/l=%zu/%zu/%zu dropA=%zu dropB=%zu\n",
                r.n_true_pairs, r.n_regions, r.n_vote_picks,
                r.n_energy_picks, r.n_largest_picks, dropA, dropB);
    }

    /* t3: self-gate -- same radius, disjoint uv: nothing flagged */
    {
        size_t nv = 0, nf = 0;
        st_grid(20.0f, 0.0f, 0.0f, GN, 0, verts, uv, gid, faces, face_cube,
                &nv, &nf);
        cube_voff[1] = nv;
        st_grid(20.0f, 0.0f, 100.0f, GN, 1, verts, uv, gid, faces, face_cube,
                &nv, &nf);
        cube_voff[2] = nv;
        PieceSet ps;
        st_wire(&ps, verts, uv, normals, gid, faces, face_cube, nv, nf,
                cube_voff, org, 2);
        SeamOwnResult r;
        int rc = SeamOwn_run(arena, &ps, &o, &r);
        so_check(rc == 0 && r.n_true_pairs == 0 && r.n_seam_pairs == 0
                 && r.n_regions == 0 && r.n_dropped == 0
                 && r.n_seam_dropped == 0,
                 "t3 disjoint uv leaves everything", &fails);
    }

    /* t4: intra-cube fold -- ONE cube, two slabs at y=20 / y=32 mapped to the
     * same uv, vertex-CONNECTED by two bridge tris (uv'd off to the side):
     * one CC with internal true pairs -> multicut fallback splits it. */
    {
        size_t nv = 0, nf = 0;
        st_grid(20.0f, 0.0f, 0.0f, GN, 0, verts, uv, gid, faces, face_cube,
                &nv, &nf);
        size_t v1 = nv;
        st_grid(32.0f, 0.0f, 0.0f, GN, 0, verts, uv, gid, faces, face_cube,
                &nv, &nf);
        /* bridge: three tris chaining slab1 edge (0,1) to slab2 edge
         * (v1,v1+1) through a mid-radius hinge vert x, with x's uv INSIDE
         * the slab band so the bridge is overlap-flagged too -- the whole
         * fold then lands in ONE region-connected component and only the
         * multicut can split it (a real fused fold). */
        size_t b0v = nv;
        verts[nv * 3 + 0] = 0.0f;
        verts[nv * 3 + 1] = 26.0f;    /* hinge radius between the plies */
        verts[nv * 3 + 2] = 10.0f;
        uv[nv * 2 + 0] = 10.0f;
        uv[nv * 2 + 1] = 10.0f;
        gid[nv] = 0;
        nv++;
        faces[nf * 3 + 0] = 0;
        faces[nf * 3 + 1] = 1;
        faces[nf * 3 + 2] = (int32_t)b0v;
        face_cube[nf] = 0;
        nf++;
        faces[nf * 3 + 0] = 1;
        faces[nf * 3 + 1] = (int32_t)b0v;
        faces[nf * 3 + 2] = (int32_t)v1;
        face_cube[nf] = 0;
        nf++;
        faces[nf * 3 + 0] = (int32_t)b0v;
        faces[nf * 3 + 1] = (int32_t)v1;
        faces[nf * 3 + 2] = (int32_t)(v1 + 1);
        face_cube[nf] = 0;
        nf++;
        cube_voff[1] = nv;
        PieceSet ps;
        st_wire(&ps, verts, uv, normals, gid, faces, face_cube, nv, nf,
                cube_voff, org, 1);
        SeamOwnResult r;
        int rc = SeamOwn_run(arena, &ps, &o, &r);
        so_check(rc == 0, "t4 rc", &fails);
        so_check(r.n_true_pairs > 0, "t4 fold pairs", &fails);
        so_check(r.n_multicut_fallbacks >= 1, "t4 multicut fired", &fails);
        so_check(r.n_dropped > 0, "t4 one ply dropped", &fails);
        fprintf(stderr, "[seam_own selftest] t4 fold: pairs=%zu mc=%zu "
                "layers=%zu dropped=%zu\n", r.n_true_pairs,
                r.n_multicut_fallbacks, r.n_layers_total, r.n_dropped);
    }

    /* t5: no RAW -- t2 geometry, raw_dir NULL: vote still resolves, seam DP
     * would be skipped (none here) */
    {
        size_t nv = 0, nf = 0;
        st_grid(20.0f, 0.0f, 0.0f, 2 * GN, 0, verts, uv, gid, faces,
                face_cube, &nv, &nf);
        cube_voff[1] = nv;
        st_grid(32.0f, 0.0f, 0.0f, GN, 1, verts, uv, gid, faces, face_cube,
                &nv, &nf);
        cube_voff[2] = nv;
        PieceSet ps;
        st_wire(&ps, verts, uv, normals, gid, faces, face_cube, nv, nf,
                cube_voff, org, 2);
        SeamOwnOpts o5 = o;
        o5.raw_dir = NULL;
        SeamOwnResult r;
        int rc = SeamOwn_run(arena, &ps, &o5, &r);
        so_check(rc == 0 && r.have_texture == 0, "t5 no texture", &fails);
        so_check(r.n_dropped > 0, "t5 vote-only still resolves", &fails);
    }

    /* t6: empty input */
    {
        PieceSet ps;
        memset(&ps, 0, sizeof(ps));
        SeamOwnResult r;
        so_check(SeamOwn_run(arena, &ps, &o, &r) == 0 && r.n_regions == 0,
                 "t6 empty input clean", &fails);
    }

    /* t7/t8/t9: REHOME. Same overlap fixture as t2 (cube B one pitch off
     * cube A, fully dropped by the vote); phi channel + spiral pins turn the
     * drop into a relocation. spiral_b = -pitch (sense -1): dr = r_B - r_A =
     * +12 => k = lround(12/-10) = -1 turn. */
    {
        float *phi = (float *)so_xcalloc(VCAP, sizeof(float));
        SeamOwnOpts o7 = o;
        o7.rehome = 1;
        o7.rehome_min_faces = 20;
        o7.spiral_a = 1.0;
        o7.spiral_b = -10.0;

        /* t7: empty destination -> committed, faces kept as dec 6, uv
         * shifted by the exact DeltaU */
        {
            size_t nv = 0, nf = 0;
            st_grid(20.0f, 0.0f, 0.0f, 2 * GN, 0, verts, uv, gid, faces,
                    face_cube, &nv, &nf);
            cube_voff[1] = nv;
            size_t vB0 = nv, fB0 = nf;
            st_grid(32.0f, 0.0f, 0.0f, GN, 1, verts, uv, gid, faces,
                    face_cube, &nv, &nf);
            cube_voff[2] = nv;
            for (size_t i = 0; i < nv; i++)
                phi[i] = i < vB0 ? 2.0f : 10.0f;   /* B: u' lands past umax */
            PieceSet ps;
            st_wire(&ps, verts, uv, normals, gid, faces, face_cube, nv, nf,
                    cube_voff, org, 2);
            ps.phi = phi;
            double u_before = uv[vB0 * 2 + 0];
            SeamOwnResult r;
            int rc = SeamOwn_run(arena, &ps, &o7, &r);
            so_check(rc == 0, "t7 rc", &fails);
            so_check(r.n_rehome_layers == 1 && r.n_rehomed_faces > 0,
                     "t7 layer rehomed", &fails);
            so_check(r.face_dec[fB0] == 6 && r.face_keep[fB0] == 1,
                     "t7 rehomed faces kept as dec 6", &fails);
            double want = u_before
                        + CubeReg_deltaU(1.0, -10.0, 10.0, -1);
            so_check(fabs((double)uv[vB0 * 2 + 0] - want) < 0.05,
                     "t7 exact DeltaU relocation", &fails);
            so_check(fabs((double)phi[vB0] - (10.0 - SO_TWO_PI)) < 1e-3,
                     "t7 phi shifted a turn", &fails);
            fprintf(stderr, "[seam_own selftest] t7 rehome: layers=%zu "
                    "faces=%zu dup=%zu blk=%zu adj=%zu small=%zu inc=%zu\n",
                    r.n_rehome_layers, r.n_rehomed_faces, r.n_rehome_dup,
                    r.n_rehome_blocked, r.n_rehome_adjacent,
                    r.n_rehome_small, r.n_rehome_incoherent);
        }

        /* t8: destination already covered at the SAME radius -> genuine
         * duplicate, stays dropped. Patch C (kept, r=32) sits where B's
         * relocation would land: phi_B tuned so u' falls into C's span. */
        {
            size_t nv = 0, nf = 0;
            st_grid(20.0f, 0.0f, 0.0f, 2 * GN, 0, verts, uv, gid, faces,
                    face_cube, &nv, &nf);
            cube_voff[1] = nv;
            size_t vB0 = nv, fB0 = nf;
            st_grid(32.0f, 0.0f, 0.0f, GN, 1, verts, uv, gid, faces,
                    face_cube, &nv, &nf);
            cube_voff[2] = nv;
            /* patch C: same cube as B, SAME radius ring (same y and x span
             * so rf matches B's), disjoint uv (u 60-80) -- never flagged,
             * always kept; it is the "material already at the destination" */
            size_t vC0 = nv;
            st_grid(32.0f, 0.0f, 60.0f, GN, 1, verts, uv, gid, faces,
                    face_cube, &nv, &nf);
            cube_voff[2] = nv;
            (void)vC0;
            /* B: u' = u + DeltaU(1,-10,phi,-1) must land in C's [60,80]:
             * DeltaU = -2pi + 10*phi - 10pi; u in [0,20] -> want shift ~60:
             * phi = (60 + 2pi + 10pi) / 10 ~= 9.77 */
            for (size_t i = 0; i < nv; i++)
                phi[i] = 2.0f;
            for (size_t i = vB0; i < vC0; i++)
                phi[i] = 9.77f;
            PieceSet ps;
            st_wire(&ps, verts, uv, normals, gid, faces, face_cube, nv, nf,
                    cube_voff, org, 2);
            ps.phi = phi;
            SeamOwnResult r;
            int rc = SeamOwn_run(arena, &ps, &o7, &r);
            so_check(rc == 0, "t8 rc", &fails);
            so_check(r.n_rehome_dup >= 1 && r.n_rehome_layers == 0,
                     "t8 covered destination stays dropped", &fails);
            so_check(r.face_keep[fB0] == 0, "t8 B still dropped", &fails);
        }

        /* t9: layer sharing a vertex with kept faces -> adjacent skip, uv
         * untouched (the tear guard) */
        {
            size_t nv = 0, nf = 0;
            st_grid(20.0f, 0.0f, 0.0f, 2 * GN, 0, verts, uv, gid, faces,
                    face_cube, &nv, &nf);
            cube_voff[1] = nv;
            size_t vB0 = nv, fB0 = nf;
            st_grid(32.0f, 0.0f, 0.0f, GN, 1, verts, uv, gid, faces,
                    face_cube, &nv, &nf);
            cube_voff[2] = nv;
            /* poison one B face with an A vertex reference */
            faces[fB0 * 3 + 0] = 0;
            for (size_t i = 0; i < nv; i++)
                phi[i] = i < vB0 ? 2.0f : 10.0f;
            PieceSet ps;
            st_wire(&ps, verts, uv, normals, gid, faces, face_cube, nv, nf,
                    cube_voff, org, 2);
            ps.phi = phi;
            double uB = uv[(vB0 + 5) * 2 + 0];
            SeamOwnResult r;
            int rc = SeamOwn_run(arena, &ps, &o7, &r);
            so_check(rc == 0, "t9 rc", &fails);
            so_check(r.n_rehome_layers == 0,
                     "t9 attached layer not rehomed", &fails);
            so_check(fabs((double)uv[(vB0 + 5) * 2 + 0] - uB) < 1e-6,
                     "t9 uv untouched", &fails);
        }
        free(phi);
    }

    free(verts); free(uv); free(normals); free(gid);
    free(faces); free(face_cube);
    Arena_dispose(&arena);
    fprintf(stderr, "[seam_own selftest] %s (%d failure%s)\n",
            fails == 0 ? "PASSED" : "FAILED", fails, fails == 1 ? "" : "s");
    return fails;
}
