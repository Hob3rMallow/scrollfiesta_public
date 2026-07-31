/* snap_grid.c -- whole-grid snap. See snap_grid.h.
 *
 * Pass 1 mirrors surface_snap.c + snap_solve.c: direct-at-vertex diagnosis,
 * signed-depth quilting target selection, then the pinned island solve.  It
 * factorizes per cube before welding.  Pass 2 is deliberately different: one
 * global oriented-edge optimization over the complete mesh places the result
 * on the recto and limits 1-ring displacement slope. */
#include "snap_grid.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "../common/csr.h"
#include "../common/kdtree.h"
#include "../common/snap_cg.h"
#include "../common/tiff_io.h"
#include "../common/ves_platform.h"
#include "../flatten/snap_darkcut.h"
#include "../flatten/snap_quilt.h"
#include "../flatten/recto_refine.h"

/* global_mode ceiling: one whole-mesh detect+solve is exact but serial; past
 * this vert count fall back to guarded per-cube (4x21x21 escalation is a
 * separate ticket -- v4 tests 4x5x5 at ~1.5M verts) */
#define SG_GLOBAL_NV_CAP ((size_t)8 * 1024 * 1024)

/* per-vertex / per-region classes (surface_snap.h values, kept local so this
 * module does not drag the flatten header in) */
enum { SG_GOOD = 0, SG_FIXABLE = 1, SG_CRACK = 2 };
enum { SGREG_FIXABLE = 0, SGREG_CRACK = 1, SGREG_ANCHORLESS = 2,
       SGREG_MIXED = 3, SGREG_REJECT = 4 };

static void *sg_xmalloc(size_t n)
{
    void *p = malloc(n > 0 ? n : 1);
    if (p == NULL) { fprintf(stderr, "snap_grid: OOM %zu\n", n); exit(1); }
    return p;
}

static void *sg_xcalloc(size_t c, size_t s)
{
    void *p = calloc(c > 0 ? c : 1, s);
    if (p == NULL) { fprintf(stderr, "snap_grid: OOM %zux%zu\n", c, s); exit(1); }
    return p;
}

void SnapGridOpts_default(SnapGridOpts *o)
{
    assert(o);
    memset(o, 0, sizeof(*o));
    /* detection = SnapOpts_default */
    o->axis_point[1] = 3405.0f;
    o->axis_point[2] = 2878.0f;
    o->axis_dir[0] = 1.0f;
    o->reach = 8.0;
    o->occ_thresh = 0.9;
    o->local_r = 2.0;
    o->gain_bias = 1.0;
    o->tensor_weight = 0.0;
    o->tensor_radius = 2.0;
    o->step = 0.25;
    o->band_frac = 0.30;
    o->min_gain = 25.0;
    o->anchor_frac = 0.25;
    o->region_cap = 4000;
    o->min_region = 8;
    o->pct_lo = 1.0;
    o->pct_hi = 99.0;
    o->dark_thresh = 64;
    o->dark_lambda = 20.0;
    o->dark_sigma = 30.0;
    o->normal_range = 2.0;
    o->normal_samples = 5;
    o->target_mode = 2; /* SNAP_TARGET_QUILT_MRF */
    o->depth_bins = 33;
    o->w_match = 30.0;
    o->w_close = 200.0;
    o->smooth_mu = 2000.0;
    o->smooth_tau = 4.0;
    /* solve = SnapSolveOpts_default */
    o->alpha = 0.85;
    o->w_anchor = 1e4;
    o->cg_max_iter = 96;
    o->cg_tol = 1e-4;
    o->max_disp = 8.0;
    o->recto_iters = 4;
    o->recto_inner_iters = 48;
    o->recto_range = 3.0;
    o->recto_max_total = 3.0;
    o->threads = 32;
    o->verbose = 0;
}

/* mirror of surface_snap.c stretch_window -- keep in sync */
static void sg_stretch_window(const double *val, const uint8_t *has, size_t nv,
                              double pct_lo, double pct_hi,
                              double *out_lo, double *out_hi)
{
    size_t hist[256];
    memset(hist, 0, sizeof(hist));
    size_t n = 0, cum = 0, tlo = 0, thi = 0;
    int b = 0, lo = 0, hi = 255;
    for (size_t i = 0; i < nv; i++) {
        if (!has[i]) continue;
        int bb = (int)(val[i] + 0.5);
        if (bb < 0) bb = 0;
        if (bb > 255) bb = 255;
        hist[bb]++;
        n++;
    }
    *out_lo = 0;
    *out_hi = 255;
    if (n == 0) return;
    tlo = (size_t)(pct_lo / 100.0 * (double)n);
    thi = (size_t)(pct_hi / 100.0 * (double)n);
    cum = 0;
    for (b = 0; b < 256; b++) { cum += hist[b]; if (cum > tlo) { lo = b; break; } }
    cum = 0;
    for (b = 0; b < 256; b++) { cum += hist[b]; if (cum >= thi) { hi = b; break; } }
    if (hi <= lo) { lo = 0; hi = 255; }
    *out_lo = lo;
    *out_hi = hi;
}

/* one cube's detect + solve. Everything indexed relative to [v0,v1) verts and
 * [f0,f1) faces (faces renumbered to local). Mutates ps->verts slice. */
typedef struct SgCubeStats {
    size_t dark, fixable, crack, rejected;
    size_t nreg, reg_fix, reg_crack, reg_anch, reg_mixed, reg_rej;
    size_t max_region;
    size_t moved, reverted;
    double sum_disp, max_disp, sum_quilt_cost, sum_target_dist;
    size_t n_quilt, quilt_fallback;
    int darkcut_fallback;
} SgCubeStats;

static void sg_cube_run(PieceSet *ps, CubeTable *ct, KDTree_T occ,
                        const float *v_orig, const double *cv,
                        const uint8_t *has, size_t v0, size_t v1,
                        size_t f0, size_t f1,
                         double dk, double band, double band_hi,
                         const SnapGridOpts *o, SgCubeStats *st)
{
    memset(st, 0, sizeof(*st));
    size_t nvc = v1 - v0, nfc = f1 - f0;
    if (nvc < 3 || nfc < 1) return;

    double occ2 = (o->occ_thresh > 0 ? o->occ_thresh : 0.9);
    occ2 *= occ2;
    double local2 = (o->local_r > 0 ? o->local_r : 2.0);
    local2 *= local2;
    Arena_T ar = Arena_new();

    /* local faces -- skip any face touching a vert outside [v0,v1): after a
     * cross-cube weld a face can reference another cube's vert; localizing it
     * blindly would corrupt memory (negative index). Global mode (v0=0,
     * v1=nv) keeps every face. */
    int32_t *lf = (int32_t *)sg_xmalloc(nfc * 3 * sizeof(int32_t));
    size_t nlf = 0;
    for (size_t f = 0; f < nfc; f++) {
        size_t a = (size_t)ps->faces[(f0 + f) * 3 + 0];
        size_t b = (size_t)ps->faces[(f0 + f) * 3 + 1];
        size_t cc = (size_t)ps->faces[(f0 + f) * 3 + 2];
        if (a < v0 || a >= v1 || b < v0 || b >= v1 || cc < v0 || cc >= v1)
            continue;
        lf[nlf * 3 + 0] = (int32_t)(a - v0);
        lf[nlf * 3 + 1] = (int32_t)(b - v0);
        lf[nlf * 3 + 2] = (int32_t)(cc - v0);
        nlf++;
    }
    if (nlf < 1) { free(lf); Arena_dispose(&ar); return; }
    nfc = nlf;
    CSR_T adj = CSR_from_faces(ar, lf, nfc, nvc);
    const int32_t *aoff = CSR_offset(adj), *atgt = CSR_target(adj);

    /* dark cut (per-cube GCO is EXACT: the MRF has no cross-cube edges) */
    uint8_t *dark = (uint8_t *)sg_xmalloc(nvc);
    {
        SnapDarkOpts dopt;
        dopt.dk = dk;
        dopt.lambda = o->dark_lambda;
        dopt.sigma = o->dark_sigma;
        dopt.band_hi = band_hi;
        if (o->dark_lambda <= 0.0
            || SnapDark_segment(cv + v0, has + v0, nvc, adj, &dopt, dark) != 0) {
            for (size_t i = 0; i < nvc; i++)
                dark[i] = (has[v0 + i] && cv[v0 + i] < dk) ? (uint8_t)1
                                                           : (uint8_t)0;
            if (o->dark_lambda > 0.0) st->darkcut_fallback = 1;
        }
    }

    /* Pass 1 target selection is shared with surface_snap.c.  The global
     * occupancy snapshot still blocks another wrap, while the signed-depth MRF
     * couples every dark vertex in this exact factorization unit. */
    uint8_t *vclass = (uint8_t *)sg_xcalloc(nvc, 1);
    float *voff = (float *)sg_xcalloc(nvc, sizeof(float));
    float *vdir = (float *)sg_xcalloc(nvc, 3 * sizeof(float));
    float *vgain = (float *)sg_xcalloc(nvc, sizeof(float));
    float *vblock = (float *)sg_xcalloc(nvc, sizeof(float));
    {
        SnapQuiltOpts qo;
        SnapQuiltStats qs;
        memset(&qo, 0, sizeof qo);
        memcpy(qo.axis_point, o->axis_point, sizeof qo.axis_point);
        memcpy(qo.axis_dir, o->axis_dir, sizeof qo.axis_dir);
        qo.reach = o->reach;
        qo.step = o->step;
        qo.occ_thresh = o->occ_thresh;
        qo.local_r = o->local_r;
        qo.band = band;
        qo.min_gain = o->min_gain;
        qo.target_mode = o->target_mode;
        qo.depth_bins = o->depth_bins;
        qo.w_match = o->w_match;
        qo.w_close = o->w_close;
        qo.smooth_mu = o->smooth_mu;
        qo.smooth_tau = o->smooth_tau;
        qo.verbose = 0;
        if (SnapQuilt_select(ar, ct, occ, v_orig, v0, v_orig + v0 * 3,
                             ps->normals != NULL ? ps->normals + v0 * 3 : NULL,
                             nvc, adj, cv + v0, has + v0, dark, &qo,
                             vclass, voff, vdir, vgain, vblock, &qs) != 0) {
            for (size_t i = 0; i < nvc; i++) {
                if (dark[i] && has[v0 + i] && aoff[i] < aoff[i + 1]) {
                    vclass[i] = SG_CRACK;
                    st->dark++;
                    st->crack++;
                }
            }
            st->quilt_fallback++;
        } else {
            st->dark = qs.n_dark;
            st->fixable = qs.n_fixable;
            st->crack = qs.n_crack;
            st->n_quilt = qs.n_fixable;
            st->sum_quilt_cost = qs.mean_quilt_cost * (double)qs.n_fixable;
            st->sum_target_dist = qs.mean_distance * (double)qs.n_fixable;
            st->quilt_fallback = qs.n_gco_fallback;
        }
    }

    /* dark regions (per cube; cannot span cubes) + classification */
    int32_t *vreg = (int32_t *)sg_xmalloc(nvc * sizeof(int32_t));
    for (size_t i = 0; i < nvc; i++) vreg[i] = -1;
    int32_t *ufp = (int32_t *)sg_xmalloc(nvc * sizeof(int32_t));
    for (size_t i = 0; i < nvc; i++) ufp[i] = (int32_t)i;
    for (size_t i = 0; i < nvc; i++) {
        if (vclass[i] == SG_GOOD) continue;
        for (int32_t e = aoff[i]; e < aoff[i + 1]; e++) {
            int32_t j = atgt[e];
            if (vclass[j] == SG_GOOD) continue;
            int32_t ra = (int32_t)i, rb = j;
            while (ufp[ra] != ra) { ufp[ra] = ufp[ufp[ra]]; ra = ufp[ra]; }
            while (ufp[rb] != rb) { ufp[rb] = ufp[ufp[rb]]; rb = ufp[rb]; }
            if (ra != rb) ufp[rb] = ra;
        }
    }
    int32_t nreg = 0;
    int32_t *root2rid = (int32_t *)sg_xmalloc(nvc * sizeof(int32_t));
    for (size_t i = 0; i < nvc; i++) root2rid[i] = -1;
    for (size_t i = 0; i < nvc; i++) {
        if (vclass[i] == SG_GOOD) continue;
        int32_t r = (int32_t)i;
        while (ufp[r] != r) { ufp[r] = ufp[ufp[r]]; r = ufp[r]; }
        if (root2rid[r] < 0) root2rid[r] = nreg++;
        vreg[i] = root2rid[r];
    }
    st->nreg = (size_t)nreg;

    int32_t *rsize = (int32_t *)sg_xcalloc((size_t)nreg, sizeof(int32_t));
    int32_t *rfix = (int32_t *)sg_xcalloc((size_t)nreg, sizeof(int32_t));
    int32_t *rbnd = (int32_t *)sg_xcalloc((size_t)nreg, sizeof(int32_t));
    int32_t *rgood = (int32_t *)sg_xcalloc((size_t)nreg, sizeof(int32_t));
    for (size_t i = 0; i < nvc; i++) {
        int32_t r = vreg[i];
        if (r < 0) continue;
        rsize[r]++;
        if (vclass[i] == SG_FIXABLE) rfix[r]++;
        for (int32_t e = aoff[i]; e < aoff[i + 1]; e++) {
            int32_t j = atgt[e];
            if (vreg[j] != r) {
                rbnd[r]++;
                if (vclass[j] == SG_GOOD) rgood[r]++;
            }
        }
    }
    uint8_t *rclass = (uint8_t *)sg_xmalloc((size_t)(nreg > 0 ? nreg : 1));
    for (int32_t r = 0; r < nreg; r++) {
        double gfrac = rbnd[r] > 0 ? (double)rgood[r] / (double)rbnd[r] : 0.0;
        double ffrac = rsize[r] > 0 ? (double)rfix[r] / (double)rsize[r] : 0.0;
        uint8_t rc;
        if (o->min_region > 1 && rsize[r] < o->min_region) rc = SGREG_REJECT;
        else if (gfrac < o->anchor_frac || rsize[r] > o->region_cap)
            rc = SGREG_ANCHORLESS;
        else if (ffrac >= 0.5) rc = SGREG_FIXABLE;
        else if ((double)(rsize[r] - rfix[r]) / (double)rsize[r] >= 0.5)
            rc = SGREG_CRACK;
        else rc = SGREG_MIXED;
        rclass[r] = rc;
        if ((size_t)rsize[r] > st->max_region)
            st->max_region = (size_t)rsize[r];
        switch (rc) {
        case SGREG_FIXABLE: st->reg_fix++; break;
        case SGREG_CRACK: st->reg_crack++; break;
        case SGREG_ANCHORLESS: st->reg_anch++; break;
        case SGREG_REJECT: st->reg_rej++; break;
        default: st->reg_mixed++; break;
        }
    }
    for (size_t i = 0; i < nvc; i++) {
        int32_t r = vreg[i];
        if (r >= 0 && rclass[r] == SGREG_REJECT && vclass[i] != SG_GOOD) {
            if (vclass[i] == SG_FIXABLE) {
                st->fixable--;
            }
            else st->crack--;
            vclass[i] = SG_GOOD;
            st->dark--;
            st->rejected++;
        }
    }

    /* pinned Laplacian solve over the cube (snap_solve.c semantics): island
     * verts pulled to target, everything else anchored at w_anchor */
    size_t n_isl = 0;
    for (size_t i = 0; i < nvc; i++)
        if (vclass[i] == SG_FIXABLE && vreg[i] >= 0
            && rclass[vreg[i]] == SGREG_FIXABLE) n_isl++;
    if (n_isl > 0) {
        double aw = (o->alpha > 0.0 && o->alpha < 1.0)
                  ? o->alpha / (1.0 - o->alpha) : 1.0;
        double maxd_cap = o->max_disp > 0 ? o->max_disp : 8.0;
        float *w = (float *)sg_xmalloc(nvc * sizeof(float));
        float *target = (float *)sg_xmalloc(nvc * 3 * sizeof(float));
        float *vs = &ps->verts[v0 * 3];
        for (size_t i = 0; i < nvc; i++) {
            int isl = (vclass[i] == SG_FIXABLE && vreg[i] >= 0
                       && rclass[vreg[i]] == SGREG_FIXABLE);
            if (isl) {
                double deg = (double)(aoff[i + 1] - aoff[i]);
                if (deg < 1.0) deg = 1.0;
                w[i] = (float)(aw * deg);
                for (int m = 0; m < 3; m++)
                    target[i * 3 + m] = (float)((double)vs[i * 3 + m]
                                        + (double)voff[i]
                                              * (double)vdir[i * 3 + m]);
            } else {
                w[i] = (float)o->w_anchor;
                target[i * 3 + 0] = vs[i * 3 + 0];
                target[i * 3 + 1] = vs[i * 3 + 1];
                target[i * 3 + 2] = vs[i * 3 + 2];
            }
        }
        SnapCG_solve(ar, adj, w, target, vs, nvc, o->cg_max_iter, o->cg_tol);

        /* cap + occupancy revert against ORIGINAL positions */
        int32_t hit[16];
        for (size_t i = 0; i < nvc; i++) {
            int isl = (vclass[i] == SG_FIXABLE && vreg[i] >= 0
                       && rclass[vreg[i]] == SGREG_FIXABLE);
            if (!isl) continue;
            size_t gi = v0 + i;
            double dz = (double)vs[i * 3 + 0] - v_orig[gi * 3 + 0];
            double dy = (double)vs[i * 3 + 1] - v_orig[gi * 3 + 1];
            double dx = (double)vs[i * 3 + 2] - v_orig[gi * 3 + 2];
            double d = sqrt(dz * dz + dy * dy + dx * dx);
            int revert = (d > maxd_cap);
            if (!revert) {
                float q[3] = { vs[i * 3 + 0], vs[i * 3 + 1], vs[i * 3 + 2] };
                size_t nh = KDTree_ball_query(occ, q, (float)occ2, hit, 16), h;
                if (nh > 16) nh = 16;
                for (h = 0; h < nh; h++) {
                    int32_t bi = hit[h];
                    double ez = (double)v_orig[(size_t)bi * 3 + 0]
                              - v_orig[gi * 3 + 0];
                    double ey = (double)v_orig[(size_t)bi * 3 + 1]
                              - v_orig[gi * 3 + 1];
                    double ex = (double)v_orig[(size_t)bi * 3 + 2]
                              - v_orig[gi * 3 + 2];
                    int own = ((size_t)bi == gi);
                    if ((size_t)bi >= v0 && (size_t)bi < v1) {
                        int32_t li = (int32_t)((size_t)bi - v0);
                        for (int32_t e = aoff[i]; !own && e < aoff[i + 1]; e++)
                            own = (atgt[e] == li);
                    }
                    if (ez * ez + ey * ey + ex * ex > local2 || !own) {
                        revert = 1;
                        break;
                    }
                }
            }
            if (revert) {
                vs[i * 3 + 0] = v_orig[gi * 3 + 0];
                vs[i * 3 + 1] = v_orig[gi * 3 + 1];
                vs[i * 3 + 2] = v_orig[gi * 3 + 2];
                st->reverted++;
            } else if (d > 1e-6) {
                st->moved++;
                st->sum_disp += d;
                if (d > st->max_disp) st->max_disp = d;
            }
        }
        free(w);
        free(target);
    }

    free(lf);
    free(dark);
    free(vclass);
    free(voff);
    free(vdir);
    free(vgain);
    free(vblock);
    free(vreg);
    free(ufp);
    free(root2rid);
    free(rsize);
    free(rfix);
    free(rbnd);
    free(rgood);
    free(rclass);
    Arena_dispose(&ar);
}

int SnapGrid_run(Arena_T arena, PieceSet *ps, CubeTable *ct,
                 const SnapGridOpts *opts, SnapGridStats *out)
{
    assert(arena && ps && opts && out);
    memset(out, 0, sizeof(*out));
    if (ps->nv < 3 || ps->nf < 1 || ct == NULL)
        return ps->nv == 0 ? 0 : -1;
    const SnapGridOpts *o = opts;
    double t0 = ves_clock_sec();

    /* GLOBAL occupancy snapshot: tree + an original-positions copy (the tree
     * copies points internally; rays and the revert guard must judge against
     * pre-move positions even while other cubes mutate theirs) */
    float *v_orig = (float *)sg_xmalloc(ps->nv * 3 * sizeof(float));
    memcpy(v_orig, ps->verts, ps->nv * 3 * sizeof(float));
    KDTree_T occ = KDTree_new(arena, v_orig, ps->nv);
    out->sec_occ = ves_clock_sec() - t0;

    /* GLOBAL direct-at-vertex intensity + window. Normal-max is correct for a
     * forgiving texture bake but wrong for geometric diagnosis: it lets an
     * off-sheet point borrow brightness from the adjacent papyrus layer. */
    t0 = ves_clock_sec();
    double *cv = (double *)sg_xmalloc(ps->nv * sizeof(double));
    uint8_t *has = (uint8_t *)sg_xcalloc(ps->nv, 1);
    {
        /* OpenMP 2.0 (MSVC): signed int loop var; nv < 2^31 checked above
         * implicitly by the CSR nnz bound; assert anyway */
        assert(ps->nv < (size_t)INT32_MAX);
        int n = (int)ps->nv;
        int i = 0;
#ifdef _OPENMP
        omp_set_num_threads(o->threads > 0 ? o->threads : 32);
#pragma omp parallel for schedule(static)
#endif
        for (i = 0; i < n; i++) {
            double s = sample_trilinear(ct,
                                        ps->verts[(size_t)i * 3 + 0],
                                        ps->verts[(size_t)i * 3 + 1],
                                        ps->verts[(size_t)i * 3 + 2]);
            cv[i] = s >= 0 ? s : 0.0;
            has[i] = s >= 0 ? (uint8_t)1 : (uint8_t)0;
        }
    }
    for (size_t i = 0; i < ps->nv; i++)
        if (has[i]) out->n_sampled++;
    double lo = 0, hi = 255;
    sg_stretch_window(cv, has, ps->nv, o->pct_lo, o->pct_hi, &lo, &hi);
    double dk = lo + ((double)o->dark_thresh / 255.0) * (hi - lo);
    double band = lo + o->band_frac * (hi - lo);
    double band_hi = band + 0.5 * (hi - band);
    out->win_lo = lo;
    out->win_hi = hi;
    out->sec_cv = ves_clock_sec() - t0;

    /* per-cube face ranges (faces are contiguous per cube in load order) */
    size_t nc = ps->n_cubes;
    size_t *foff = (size_t *)sg_xmalloc((nc + 1) * sizeof(size_t));
    {
        size_t c = 0;
        foff[0] = 0;
        for (size_t f = 0; f < ps->nf; f++) {
            while (c < (size_t)ps->face_cube[f]) foff[++c] = f;
        }
        while (c < nc) foff[++c] = ps->nf;
    }

    /* detect + solve: one whole-mesh problem (global_mode, exact post-weld)
     * or parallel over cubes (exact pre-weld; face-range guard otherwise) */
    t0 = ves_clock_sec();
    int global = o->global_mode != 0;
    if (global && ps->nv > SG_GLOBAL_NV_CAP) {
        fprintf(stderr,
                "[snap_grid] WARNING: global_mode with nv=%zu > cap %zu; "
                "falling back to per-cube (cross-cube faces skipped)\n",
                ps->nv, (size_t)SG_GLOBAL_NV_CAP);
        global = 0;
    }
    size_t n_units = global ? 1 : nc;
    SgCubeStats *cst = (SgCubeStats *)sg_xcalloc(n_units, sizeof(SgCubeStats));
    if (global) {
        sg_cube_run(ps, ct, occ, v_orig, cv, has, 0, ps->nv, 0, ps->nf,
                    dk, band, band_hi, o, &cst[0]);
    } else {
        int n = (int)nc;
        int c = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
        for (c = 0; c < n; c++) {
            sg_cube_run(ps, ct, occ, v_orig, cv, has,
                        ps->cube_voff[c], ps->cube_voff[c + 1],
                        foff[c], foff[c + 1], dk, band, band_hi,
                        o, &cst[c]);
        }
    }
    out->sec_cubes = ves_clock_sec() - t0;

    size_t n_quilt = 0;
    for (size_t c = 0; c < n_units; c++) {
        const SgCubeStats *s = &cst[c];
        out->n_dark += s->dark;
        out->n_fixable += s->fixable;
        out->n_crack += s->crack;
        out->n_rejected += s->rejected;
        out->nreg += s->nreg;
        out->n_reg_fixable += s->reg_fix;
        out->n_reg_crack += s->reg_crack;
        out->n_reg_anchorless += s->reg_anch;
        out->n_reg_mixed += s->reg_mixed;
        out->n_reg_reject += s->reg_rej;
        if (s->max_region > out->max_region_size)
            out->max_region_size = s->max_region;
        out->n_moved += s->moved;
        out->n_reverted += s->reverted;
        out->mean_disp += s->sum_disp;
        out->mean_quilt_cost += s->sum_quilt_cost;
        out->mean_target_dist += s->sum_target_dist;
        n_quilt += s->n_quilt;
        if (s->max_disp > out->max_disp) out->max_disp = s->max_disp;
        out->n_cubes_run++;
        out->n_cubes_gco_fallback += (size_t)s->darkcut_fallback;
        out->n_quilt_fallback += s->quilt_fallback;
    }
    out->mean_disp = out->n_moved > 0 ? out->mean_disp / (double)out->n_moved
                                       : 0.0;
    out->mean_quilt_cost = n_quilt > 0
                         ? out->mean_quilt_cost / (double)n_quilt : 0.0;
    out->mean_target_dist = n_quilt > 0
                          ? out->mean_target_dist / (double)n_quilt : 0.0;

    /* Pass 2 never factorizes: recto placement is one optimization over every
     * face-connected vertex, including edges introduced by step2_join. */
    if (o->recto_iters > 0) {
        RectoRefineOpts ro;
        RectoRefineStats rs;
        RectoRefineOpts_default(&ro);
        memcpy(ro.axis_point, o->axis_point, sizeof ro.axis_point);
        memcpy(ro.axis_dir, o->axis_dir, sizeof ro.axis_dir);
        ro.window_lo = lo;
        ro.window_hi = hi;
        if (o->recto_range > 0.0) ro.range = o->recto_range;
        if (o->recto_max_total > 0.0) ro.max_total = o->recto_max_total;
        ro.outer_iters = o->recto_iters;
        if (o->recto_inner_iters > 0) ro.inner_iters = o->recto_inner_iters;
        ro.threads = o->threads;
        ro.verbose = o->verbose;
        if (RectoRefine_run(arena, ct, ps->verts, ps->nv, ps->faces, ps->nf,
                            &ro, &rs) != 0) {
            free(v_orig); free(cv); free(has); free(foff); free(cst);
            return -1;
        }
        out->n_recto_supported = rs.n_supported;
        out->n_recto_moved = rs.n_moved;
        out->n_recto_reverted = rs.n_reverted;
        out->n_recto_slope_limited = rs.n_slope_limited;
        out->recto_iterations = rs.iterations;
        out->recto_mean_disp = rs.mean_disp;
        out->recto_max_disp = rs.max_disp;
        out->recto_mean_gradient = rs.mean_recto_gradient;
        out->sec_recto = rs.sec_total;
    }

    if (o->verbose)
        fprintf(stderr,
                "[snap_grid] sampled=%zu dark=%zu fixable=%zu crack=%zu "
                "rejected=%zu | regions=%zu fix/crack/anch/mixed/rej="
                "%zu/%zu/%zu/%zu/%zu max=%zu | moved=%zu reverted=%zu "
                "mean_disp=%.2f max=%.2f | quilt=%.1f target_dist=%.2f "
                "window[%.0f,%.0f] cubes=%zu "
                "dark_gco_fb=%zu quilt_fb=%zu | recto support=%zu moved=%zu "
                "reverted=%zu slope=%zu mean=%.2f max=%.2f iter=%d "
                "(occ %.1fs cv %.1fs pass1 %.1fs recto %.1fs)\n",
                out->n_sampled, out->n_dark, out->n_fixable, out->n_crack,
                out->n_rejected, out->nreg, out->n_reg_fixable,
                out->n_reg_crack, out->n_reg_anchorless, out->n_reg_mixed,
                out->n_reg_reject, out->max_region_size, out->n_moved,
                out->n_reverted, out->mean_disp, out->max_disp,
                out->mean_quilt_cost, out->mean_target_dist, lo, hi,
                out->n_cubes_run, out->n_cubes_gco_fallback,
                out->n_quilt_fallback, out->n_recto_supported,
                out->n_recto_moved, out->n_recto_reverted,
                out->n_recto_slope_limited, out->recto_mean_disp,
                out->recto_max_disp, out->recto_iterations,
                out->sec_occ, out->sec_cv, out->sec_cubes, out->sec_recto);

    free(v_orig);
    free(cv);
    free(has);
    free(foff);
    free(cst);
    return 0;
}

/* ============================================================================
 * Self-test.
 * ==========================================================================*/

static int sg_check(int cond, const char *what, int *fails)
{
    if (!cond) {
        fprintf(stderr, "[snap_grid selftest]   FAIL: %s\n", what);
        (*fails)++;
    }
    return cond;
}

#define SG_RAW_DIR "output/_selftest_snap_grid/raw"

/* fixture volume: a dark POCKET (x in [16,22], y < 32) inside otherwise
 * bright papyrus. A patch at y=28 spanning x=[10,26] then has bright GOOD
 * rims and a dark FIXABLE island whose +y probe reaches bright at y>=32
 * (its -y probe stays inside the dark pocket: no downward target). */
static int sg_write_raw(int bright_blocker)
{
    long chunk = 64;
    uint8_t *vol = (uint8_t *)sg_xmalloc((size_t)(chunk * chunk * chunk));
    for (long z = 0; z < chunk; z++)
        for (long y = 0; y < chunk; y++)
            for (long x = 0; x < chunk; x++)
                vol[(z * chunk + y) * chunk + x] =
                    (y < 32 && x >= 16 && x <= 22
                     && !(bright_blocker && y >= 30 && y <= 31))
                    ? (uint8_t)10 : (uint8_t)200;
    char path[512];
    snprintf(path, sizeof(path), "%s/z00000_y00000_x00000.tif", SG_RAW_DIR);
    ves_ensure_parent_dir(path);
    int rc = TiffIO_save(path, vol, (int)chunk, (int)chunk, (int)chunk);
    free(vol);
    return rc;
}

/* grid patch in the y=Y plane with tiny jitter (perfect lattices are a
 * KD-tree degeneracy; snap_solve.c uses the same trick) */
static void sg_st_grid(float y, int n, double spacing, int32_t cube,
                       float *verts, float *uv,
                       float *normals, int32_t *faces, int32_t *face_cube,
                       size_t *pv, size_t *pf)
{
    size_t v0 = *pv;
    int side = n + 1;
    for (int iz = 0; iz < side; iz++) {
        for (int ix = 0; ix < side; ix++) {
            size_t vi = v0 + (size_t)iz * (size_t)side + (size_t)ix;
            unsigned hh = (unsigned)vi * 2654435761u;
            double jz = ((double)(hh & 255) / 255.0 - 0.5) * 0.03;
            double jy = ((double)((hh >> 8) & 255) / 255.0 - 0.5) * 0.03;
            double jx = ((double)((hh >> 16) & 255) / 255.0 - 0.5) * 0.03;
            verts[vi * 3 + 0] = (float)(10 + spacing * iz + jz);
            verts[vi * 3 + 1] = (float)(y + jy);
            verts[vi * 3 + 2] = (float)(10 + spacing * ix + jx);
            uv[vi * 2 + 0] = (float)(spacing * ix);
            uv[vi * 2 + 1] = (float)(spacing * iz);
            normals[vi * 3 + 0] = 0.0f;
            normals[vi * 3 + 1] = 1.0f;
            normals[vi * 3 + 2] = 0.0f;
        }
    }
    for (int iz = 0; iz < n; iz++) {
        for (int ix = 0; ix < n; ix++) {
            size_t a = v0 + (size_t)iz * (size_t)side + (size_t)ix;
            size_t b = a + 1, c = a + (size_t)side, d = c + 1;
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

int SnapGrid_selftest(void)
{
    int fails = 0;
    Arena_T arena = Arena_new();
    if (sg_write_raw(0) != 0) {
        fprintf(stderr, "snap_grid selftest: cannot write raw fixture\n");
        Arena_dispose(&arena);
        return 1;
    }

    enum { VCAP = 1024, FCAP = 1024 };
    float *verts = (float *)sg_xmalloc(VCAP * 3 * sizeof(float));
    float *uv = (float *)sg_xmalloc(VCAP * 2 * sizeof(float));
    float *normals = (float *)sg_xmalloc(VCAP * 3 * sizeof(float));
    int32_t *gid = (int32_t *)sg_xcalloc(VCAP, sizeof(int32_t));
    int32_t *faces = (int32_t *)sg_xmalloc(FCAP * 3 * sizeof(int32_t));
    int32_t *face_cube = (int32_t *)sg_xmalloc(FCAP * sizeof(int32_t));
    size_t cube_voff[3] = { 0, 0, 0 };
    long org[2][3] = { { 0, 0, 0 }, { 0, 0, 0 } };

    SnapGridOpts o;
    SnapGridOpts_default(&o);
    o.axis_point[0] = 0.0f;
    o.axis_point[1] = 0.0f;
    o.axis_point[2] = 0.0f;
    o.threads = 2;
    o.min_region = 0;      /* tiny fixtures: keep every dark vert */
    o.verbose = 0;

    /* shared PieceSet shell */
    PieceSet ps;

    /* t1: patch at y=28 spanning the dark pocket: bright rim = GOOD, pocket
     * verts = dark FIXABLE island; +y probe reaches bright at y>=32 through
     * free space; solve moves the island up, rim stays. Deliberately make the
     * legacy normal-max reach the slab: direct-at-vertex diagnosis must still
     * see the pocket as dark. */
    {
        o.normal_range = 4.5;
        size_t nv = 0, nf = 0;
        sg_st_grid(28.0f, 8, 2.0, 0, verts, uv, normals, faces, face_cube,
                   &nv, &nf);
        cube_voff[1] = nv;
        memset(&ps, 0, sizeof(ps));
        ps.verts = verts; ps.uv = uv; ps.normals = normals; ps.gid = gid;
        ps.faces = faces; ps.face_cube = face_cube;
        ps.nv = nv; ps.nf = nf; ps.n_cubes = 1;
        ps.cube_voff = cube_voff; ps.cube_org = org;
        CubeTable ct;
        if (cubetable_init(&ct, arena, SG_RAW_DIR, 64, ps.verts, ps.nv,
                           o.normal_range + o.reach + 2.0) != 0) {
            sg_check(0, "t1 cubetable", &fails);
        } else {
            cubetable_prewarm_all(&ct);
            SnapGridStats st;
            int rc = SnapGrid_run(arena, &ps, &ct, &o, &st);
            sg_check(rc == 0, "t1 rc", &fails);
            sg_check(st.n_dark > 0, "t1 dark detected", &fails);
            sg_check(st.n_fixable > 0 && st.n_reg_fixable >= 1,
                     "t1 fixable island", &fails);
            sg_check(st.n_moved > 0 && st.n_reverted == 0,
                     "t1 moved, none reverted", &fails);
            /* island verts (x in the pocket) moved up; rim stayed */
            double ymax_isl = 0.0, ymove_rim = 0.0;
            for (size_t i = 0; i < nv; i++) {
                double x = verts[i * 3 + 2], yv = verts[i * 3 + 1];
                if (x >= 17.0 && x <= 21.0) {
                    if (yv > ymax_isl) ymax_isl = yv;
                } else if (x < 14.0 || x > 24.0) {
                    double dy = fabs(yv - 28.0);
                    if (dy > ymove_rim) ymove_rim = dy;
                }
            }
            sg_check(ymax_isl > 31.0, "t1 island lands close to bright core",
                      &fails);
            sg_check(ymove_rim < 0.5, "t1 rim stays put", &fails);
            fprintf(stderr, "[snap_grid selftest] t1: dark=%zu fix=%zu "
                    "moved=%zu ymax_isl=%.1f rim_dy=%.2f\n", st.n_dark,
                    st.n_fixable, st.n_moved, ymax_isl, ymove_rim);
        }
        o.normal_range = 2.0;
    }

    /* t2: same misplaced patch, but ANOTHER cube's sheet sits at y=30.5
     * between it and the bright side: the +y probe is BLOCKED (occupancy)
     * and -y is dark -> CRACK, nothing moves across the blocker */
    {
        sg_check(sg_write_raw(1) == 0, "t2 bright blocker fixture", &fails);
        size_t nv = 0, nf = 0;
        sg_st_grid(28.0f, 8, 2.0, 0, verts, uv, normals, faces, face_cube,
                   &nv, &nf);
        cube_voff[1] = nv;
        /* Dense blocker, padded beyond the source's outward-radial x extent;
         * 1-vox spacing seals the 0.9-vox occupancy ball. */
        sg_st_grid(30.5f, 20, 1.0, 1, verts, uv, normals, faces, face_cube,
                   &nv, &nf);
        cube_voff[2] = nv;
        memset(&ps, 0, sizeof(ps));
        ps.verts = verts; ps.uv = uv; ps.normals = normals; ps.gid = gid;
        ps.faces = faces; ps.face_cube = face_cube;
        ps.nv = nv; ps.nf = nf; ps.n_cubes = 2;
        ps.cube_voff = cube_voff; ps.cube_org = org;
        CubeTable ct;
        if (cubetable_init(&ct, arena, SG_RAW_DIR, 64, ps.verts, ps.nv,
                           o.normal_range + o.reach + 2.0) != 0) {
            sg_check(0, "t2 cubetable", &fails);
        } else {
            cubetable_prewarm_all(&ct);
            float y0_28 = verts[0 * 3 + 1];
            SnapGridStats st;
            int rc = SnapGrid_run(arena, &ps, &ct, &o, &st);
            /* lower patch verts must NOT have crossed the y=30.5 sheet */
            int crossed = 0;
            for (size_t i = 0; i < cube_voff[1]; i++)
                if (verts[i * 3 + 1] > 30.4f) crossed = 1;
            sg_check(rc == 0, "t2 rc", &fails);
            sg_check(!crossed, "t2 no cross-wrap move through the blocker",
                     &fails);
            (void)y0_28;
            fprintf(stderr, "[snap_grid selftest] t2: dark=%zu fix=%zu "
                    "crack=%zu repair_moved=%zu repair_reverted=%zu "
                    "recto_moved=%zu recto_reverted=%zu recto_max=%.2f\n",
                    st.n_dark, st.n_fixable, st.n_crack, st.n_moved,
                    st.n_reverted, st.n_recto_moved, st.n_recto_reverted,
                    st.recto_max_disp);
        }
        sg_check(sg_write_raw(0) == 0, "restore base RAW fixture", &fails);
    }

    /* t3: empty */
    {
        memset(&ps, 0, sizeof(ps));
        SnapGridStats st;
        sg_check(SnapGrid_run(arena, &ps, NULL, &o, &st) == 0,
                 "t3 empty clean", &fails);
    }

    /* t4: post-weld shape -- ONE t1 patch whose vert range is split across
     * two "cubes" (voff at 40), so mid-grid faces straddle the boundary.
     * (a) per-cube mode must skip the straddlers (face-range guard), not
     * corrupt local indices; (b) global_mode must reproduce t1's fix. */
    {
        size_t nv = 0, nf = 0;
        sg_st_grid(28.0f, 8, 2.0, 0, verts, uv, normals, faces, face_cube,
                   &nv, &nf);
        for (size_t f = 0; f < nf; f++)
            face_cube[f] = f < nf / 2 ? 0 : 1;
        cube_voff[1] = 40;
        cube_voff[2] = nv;
        memset(&ps, 0, sizeof(ps));
        ps.verts = verts; ps.uv = uv; ps.normals = normals; ps.gid = gid;
        ps.faces = faces; ps.face_cube = face_cube;
        ps.nv = nv; ps.nf = nf; ps.n_cubes = 2;
        ps.cube_voff = cube_voff; ps.cube_org = org;
        CubeTable ct;
        if (cubetable_init(&ct, arena, SG_RAW_DIR, 64, ps.verts, ps.nv,
                           o.normal_range + o.reach + 2.0) != 0) {
            sg_check(0, "t4 cubetable", &fails);
        } else {
            cubetable_prewarm_all(&ct);
            SnapGridStats st;
            int rc = SnapGrid_run(arena, &ps, &ct, &o, &st);
            sg_check(rc == 0, "t4a per-cube guarded rc", &fails);
            sg_check(st.n_cubes_run == 2, "t4a two units", &fails);
            int sane = 1;
            for (size_t i = 0; i < nv * 3; i++)
                if (!(verts[i] > -100.0f && verts[i] < 200.0f)) sane = 0;
            sg_check(sane, "t4a verts sane after guarded run", &fails);

            /* fresh geometry, same table coverage; global mode */
            size_t nv2 = 0, nf2 = 0;
            sg_st_grid(28.0f, 8, 2.0, 0, verts, uv, normals, faces,
                       face_cube, &nv2, &nf2);
            for (size_t f = 0; f < nf2; f++)
                face_cube[f] = f < nf2 / 2 ? 0 : 1;
            o.global_mode = 1;
            rc = SnapGrid_run(arena, &ps, &ct, &o, &st);
            o.global_mode = 0;
            sg_check(rc == 0, "t4b global rc", &fails);
            sg_check(st.n_cubes_run == 1, "t4b one unit", &fails);
            sg_check(st.n_dark > 0 && st.n_moved > 0,
                     "t4b global dark+moved", &fails);
            double ymax_isl = 0.0, ymove_rim = 0.0;
            for (size_t i = 0; i < nv2; i++) {
                double x = verts[i * 3 + 2], yv = verts[i * 3 + 1];
                if (x >= 17.0 && x <= 21.0) {
                    if (yv > ymax_isl) ymax_isl = yv;
                } else if (x < 14.0 || x > 24.0) {
                    double dy = fabs(yv - 28.0);
                    if (dy > ymove_rim) ymove_rim = dy;
                }
            }
            sg_check(ymax_isl > 29.5, "t4b island moved (global)", &fails);
            sg_check(ymove_rim < 0.5, "t4b rim stays (global)", &fails);
            fprintf(stderr, "[snap_grid selftest] t4: guarded cubes=2 ok; "
                    "global dark=%zu moved=%zu ymax=%.1f rim=%.2f\n",
                    st.n_dark, st.n_moved, ymax_isl, ymove_rim);
        }
    }

    free(verts); free(uv); free(normals); free(gid);
    free(faces); free(face_cube);
    Arena_dispose(&arena);
    fprintf(stderr, "[snap_grid selftest] %s (%d failure%s)\n",
            fails == 0 ? "PASSED" : "FAILED", fails, fails == 1 ? "" : "s");
    return fails;
}
