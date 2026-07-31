/* ============================================================================
 * ribbon.c -- slice / arc-length / joint parameterization + ribbon fit.
 * See ribbon.h for the pipeline overview. StrokeStrip (Pagurek van Mossel,
 * Liu, Vining, Bessmeltsev, Sheffer 2021) adapted to a scroll sheet:
 *   strokes   = per-slice plane-mesh intersection polylines
 *   strip     = the scroll sheet itself
 *   WTS pairs = cross-slice nearest samples (+ intra-slice continuations)
 *   isolines  = vertical rulings on the papyrus
 *
 * Winding is derived per chain in the SLICE domain: local unwound theta plus
 * an integer wrap index W anchored GLOBALLY by radius vs. the umbilicus
 * (W = round(median(r/pitch - theta/2pi)) per chain). A mesh-graph-integrated
 * angle field (unwrap.c) is NOT used: residual inter-wrap fusion bridges make
 * its BFS shortcut between wraps, collapsing whole turns (verified on
 * grid_4x5x5_v5_full: 27% of UV cells held >=2 wraps). In the slice domain
 * those bridges are visible as radial jumps and are cut geometrically; the
 * radius anchor then places every fragment on its true wrap.
 * ==========================================================================*/
#define _USE_MATH_DEFINES
#include "ribbon.h"

#include "../common/csr.h"
#include "../common/kdtree.h"
#include "../common/pca.h"
#include "../common/union_find.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---- module constants ------------------------------------------------------ */
#define RIB_CONT_GAP_VOX    6.0    /* intra-slice fragment continuation reach */
#define RIB_CONT_DR_MAX     3.0    /* continuation must stay on ITS wrap: a
                                    * genuine fragment continues at ~the same
                                    * radius; a cut fusion-wall remnant jumps
                                    * by ~the pitch and must NOT be re-glued */
#define RIB_CONT_DPHI_MAX   (M_PI / 2.0)
#define RIB_ALIGN_W         0.5    /* alignment weight per pair (x likelihood) */
#define RIB_BRIDGE_FRAC     0.75   /* bridge cut: |dr| > frac*sample_h per step
                                    * (a spiral ramp moves pitch/(2pi r) ~ 0.01
                                    * vox/vox; a fusion wall ~1) */
#define RIB_XFER_R          3.4    /* mesh-vert transfer search radius (vox);
                                    * < 3.5 = half the 7-vox wrap clearance */
#define RIB_VJUMP_MAX       3      /* cross-slice matching stride escalation: a
                                    * sample unmatched at k+1 tries k+2, k+3.
                                    * An EMPTY slice (hole band, or a plane
                                    * coinciding exactly with a vertex row)
                                    * would otherwise sever the pair graph
                                    * vertically and shatter the chart. */
#define RIB_LAMBDA_C        10.0   /* ribbon row tangent smoothness */
#define RIB_LAMBDA_T        10.0   /* ribbon row position-from-tangent weight */
#define RIB_VSMOOTH_MU      0.5    /* ribbon vertical smoothing blend */
#define RIB_CG_TOL          1e-9   /* relative residual */
#define RIB_CG_MAXIT        800
#define RIB_MIN_CHAIN_LEN   0.25   /* drop slivers shorter than this (vox) */
#define RIB_VFILL_MAX       8      /* max rows bridged by vertical grid fill */
#define RIB_UFILL_MAX       8      /* max invalid u-columns bridged WITHIN a
                                    * ribbon row: integrating tangents across
                                    * a long hole extrudes giant free-space
                                    * arcs ("hula hoops") -- long holes stay
                                    * holes and the row splits into segments */
#define RIB_WRAP_GATE       6.0    /* max 3D step (vox) between adjacent ribbon
                                    * grid cells. < the 7-vox inter-wrap
                                    * clearance, so a grid edge can never span
                                    * two wraps: where a row or column would
                                    * jump wraps (winding collapse / a hole
                                    * bridging across wraps), the grid SPLITS
                                    * instead of drawing a spike into space. */
#define RIB_UVFILL_ROUNDS   10     /* neighbor-fill rounds for unmapped verts */
#define RIB_PRIOR_MIN_B     1.0    /* spiral prior usable: |b| above this ... */
#define RIB_PRIOR_MIN_R2    0.15   /* ... and fit r2 above this */
#define RIB_WEDGE_NTHETA    24     /* angular bins for the radial-ordering prior */
#define RIB_PITCH_MIN       3.0    /* auto pitch estimate: ignore radial gaps < */
#define RIB_PITCH_MAX       40.0   /* ... and > these (delamination / multi-hole) */
#define RIB_PITCH_DEFAULT   9.5    /* PHerc0139 fallback pitch (vox/turn) */
#define RIB_MAXW            6      /* max wraps a single radial step may span */

/* ---- small helpers --------------------------------------------------------- */
static int cmp_u64(const void *pa, const void *pb)
{
    uint64_t a = *(const uint64_t *)pa, b = *(const uint64_t *)pb;
    return a < b ? -1 : (a > b ? 1 : 0);
}
static int cmp_dbl(const void *pa, const void *pb)
{
    double a = *(const double *)pa, b = *(const double *)pb;
    return a < b ? -1 : (a > b ? 1 : 0);
}
static int cmp_i32(const void *pa, const void *pb)
{
    int32_t a = *(const int32_t *)pa, b = *(const int32_t *)pb;
    return a < b ? -1 : (a > b ? 1 : 0);
}
static double v3dot(const double a[3], const double b[3])
{ return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }
static double v3norm(double a[3])
{
    double n = sqrt(v3dot(a, a));
    if (n > 1e-12) { a[0] /= n; a[1] /= n; a[2] /= n; }
    return n;
}
static double wrap_pi(double a)
{
    double x = fmod(a + M_PI, 2.0 * M_PI);
    if (x < 0.0) x += 2.0 * M_PI;
    return x - M_PI;
}

void RibbonOpts_default(RibbonOpts *opts)
{
    assert(opts);
    memset(opts, 0, sizeof(*opts));
    opts->axis_dir[0] = 1.0f;   /* Z in (z,y,x) */
    opts->slice_h       = 2.0f;
    opts->sample_h      = 2.0f;
    opts->match_r       = 3.0f;
    opts->match_ang_deg = 20.0f;
    opts->relax_iters   = 5;
    opts->final_iters   = 2;
    opts->grid_u        = 2.0f;
    opts->fit_ribbon    = 1;
    opts->wrap_spacing  = 0.0f;   /* auto-estimate */
}

/* ============================================================================
 * Stage A -- slicing: plane-mesh crossings keyed by mesh edge, chained into
 * per-plane polylines, oriented by ascending theta, resampled, bridge-cut.
 * ==========================================================================*/

typedef struct {
    double p[3];      /* position (z,y,x) */
    double tau[3];    /* unit tangent along the chain (theta-ascending) */
    double c1, c2;    /* in-plane coords relative to the axis */
    double r;         /* in-plane radius */
    double th;        /* local unwound theta along the chain (radians) */
    double phi;       /* th + 2pi*W after winding assignment */
    double s;         /* cumulative arc length along the chain (vox) */
    double u;         /* solved parameter (stage C) */
    int32_t chain;    /* -1 = orphaned by bridge cut */
    int32_t slice;
} Sample;

typedef struct {
    int32_t first;    /* index of first sample */
    int32_t count;    /* number of samples */
    int32_t slice;
    int32_t closed;   /* was a closed loop (cut open) */
    int32_t wind;     /* integer winding index W (stage W) */
    int32_t group;    /* chain-graph group id */
} Chain;

typedef struct {
    Sample *smp;      size_t n_smp;
    Chain  *chn;      size_t n_chn;
    int     nplanes;
    double  tmin, tmax;
    size_t  bridge_cuts;
    int     n_closed;
    int     n_slices_hit;     /* planes with >= 1 chain */
    int     n_multi;          /* planes with > 1 chain */
} SliceSet;

/* Strict-interior plane range for a t-interval [tlo, thi]:
 * planes t_k = tmin + (k + 0.5) h with tlo < t_k < thi.
 * A vertex exactly on a plane counts as ABOVE it (consistent nudge). */
static void plane_range(double tlo, double thi, double tmin, double h,
                        int nplanes, int *k0, int *k1)
{
    double x0 = (tlo - tmin) / h - 0.5;
    double x1 = (thi - tmin) / h - 0.5;
    int a = (int)floor(x0) + 1;
    int b = (int)ceil(x1) - 1;
    if (a < 0) a = 0;
    if (b > nplanes - 1) b = nplanes - 1;
    *k0 = a; *k1 = b;
}

static int slice_mesh(Arena_T arena,
                      const float *verts, size_t nv,
                      const int32_t *faces, size_t nf,
                      const double *t,
                      double tmin, double tmax,
                      const double e1[3], const double e2[3],
                      const float axis_point[3],
                      const RibbonOpts *o, SliceSet *S)
{
    memset(S, 0, sizeof(*S));
    S->tmin = tmin; S->tmax = tmax;
    double h = (double)o->slice_h;
    int nplanes = (int)floor((tmax - tmin) / h);
    if (nplanes < 1) return -1;
    S->nplanes = nplanes;

    /* --- unique undirected edges (key = lo*nv + hi) --- */
    size_t nde = nf * 3;
    uint64_t *keys = (uint64_t *)ARENA_ALLOC(arena, (long)(nde * sizeof(uint64_t)));
    for (size_t f = 0; f < nf; f++) {
        for (int m = 0; m < 3; m++) {
            int32_t a = faces[f*3 + (size_t)m], b = faces[f*3 + (size_t)((m+1)%3)];
            uint64_t lo = (uint64_t)(a < b ? a : b), hi = (uint64_t)(a < b ? b : a);
            keys[f*3 + (size_t)m] = lo * (uint64_t)nv + hi;
        }
    }
    uint64_t *ukeys = (uint64_t *)ARENA_ALLOC(arena, (long)(nde * sizeof(uint64_t)));
    memcpy(ukeys, keys, nde * sizeof(uint64_t));
    qsort(ukeys, nde, sizeof(uint64_t), cmp_u64);
    size_t ne = 0;
    for (size_t i = 0; i < nde; i++)
        if (i == 0 || ukeys[i] != ukeys[i-1]) ukeys[ne++] = ukeys[i];

    /* per-face edge ids (bsearch into the unique list) */
    int32_t *feid = (int32_t *)ARENA_ALLOC(arena, (long)(nde * sizeof(int32_t)));
    for (size_t i = 0; i < nde; i++) {
        uint64_t *hit = (uint64_t *)bsearch(&keys[i], ukeys, ne,
                                            sizeof(uint64_t), cmp_u64);
        assert(hit != NULL);
        feid[i] = (int32_t)(hit - ukeys);
    }

    /* --- per-edge crossing ranges --- */
    int32_t *ek0 = (int32_t *)ARENA_ALLOC(arena, (long)(ne * sizeof(int32_t)));
    int32_t *ecn = (int32_t *)ARENA_ALLOC(arena, (long)(ne * sizeof(int32_t)));
    int32_t *est = (int32_t *)ARENA_ALLOC(arena, (long)((ne + 1) * sizeof(int32_t)));
    size_t total = 0;
    for (size_t e = 0; e < ne; e++) {
        int32_t a = (int32_t)(ukeys[e] / (uint64_t)nv);
        int32_t b = (int32_t)(ukeys[e] % (uint64_t)nv);
        ek0[e] = 0; ecn[e] = 0;
        double ta = t[a], tb = t[b];
        double tlo = ta < tb ? ta : tb, thi = ta < tb ? tb : ta;
        int k0 = 0, k1 = -1;
        plane_range(tlo, thi, tmin, h, nplanes, &k0, &k1);
        if (k1 < k0) continue;
        ek0[e] = k0;
        ecn[e] = k1 - k0 + 1;
        total += (size_t)ecn[e];
    }
    est[0] = 0;
    for (size_t e = 0; e < ne; e++) est[e+1] = est[e] + ecn[e];
    if (total == 0) return -1;
    if (total > (size_t)INT32_MAX / 4) return -1;   /* absurd input */

    /* --- fill crossings --- */
    double  *cpos  = (double *)ARENA_ALLOC(arena, (long)(total * 3 * sizeof(double)));
    double  *cphi  = o->reference_phi != NULL
                   ? (double *)ARENA_ALLOC(arena, (long)(total * sizeof(double)))
                   : NULL;
    int32_t *cpln  = (int32_t *)ARENA_ALLOC(arena, (long)(total * sizeof(int32_t)));
    for (size_t e = 0; e < ne; e++) {
        if (ecn[e] == 0) continue;
        int32_t a = (int32_t)(ukeys[e] / (uint64_t)nv);
        int32_t b = (int32_t)(ukeys[e] % (uint64_t)nv);
        double ta = t[a], tb = t[b];
        for (int32_t k = ek0[e]; k < ek0[e] + ecn[e]; k++) {
            double tk = tmin + ((double)k + 0.5) * h;
            double w = (tk - ta) / (tb - ta);
            if (w < 0.0) w = 0.0;
            if (w > 1.0) w = 1.0;
            size_t ci = (size_t)est[e] + (size_t)(k - ek0[e]);
            for (int d = 0; d < 3; d++)
                cpos[ci*3 + (size_t)d] = (1.0 - w) * (double)verts[(size_t)a*3 + (size_t)d]
                                       +        w  * (double)verts[(size_t)b*3 + (size_t)d];
            if (cphi != NULL)
                cphi[ci] = (1.0 - w) * (double)o->reference_phi[a]
                         +        w  * (double)o->reference_phi[b];
            cpln[ci] = k;
        }
    }

    /* --- segments: link the two crossed edges of each straddling face --- */
    int32_t *nbr = (int32_t *)ARENA_ALLOC(arena, (long)(total * 2 * sizeof(int32_t)));
    for (size_t i = 0; i < total * 2; i++) nbr[i] = -1;
    for (size_t f = 0; f < nf; f++) {
        double f0 = t[faces[f*3]], f1 = t[faces[f*3+1]], f2 = t[faces[f*3+2]];
        double lo = f0, hi = f0;
        if (f1 < lo) lo = f1;
        if (f1 > hi) hi = f1;
        if (f2 < lo) lo = f2;
        if (f2 > hi) hi = f2;
        int k0 = 0, k1 = -1;
        plane_range(lo, hi, tmin, h, nplanes, &k0, &k1);
        for (int32_t k = k0; k <= k1; k++) {
            int32_t slot[3]; int nslot = 0;
            for (int m = 0; m < 3; m++) {
                int32_t e = feid[f*3 + (size_t)m];
                if (ecn[e] > 0 && k >= ek0[e] && k < ek0[e] + ecn[e]) {
                    if (nslot < 3) slot[nslot] = est[e] + (k - ek0[e]);
                    nslot++;
                }
            }
            if (nslot != 2) continue;   /* vertex-on-plane corner */
            int32_t sa = slot[0], sb = slot[1];
            int ok = 1;
            if      (nbr[(size_t)sa*2]   == -1) nbr[(size_t)sa*2]   = sb;
            else if (nbr[(size_t)sa*2+1] == -1) nbr[(size_t)sa*2+1] = sb;
            else ok = 0;
            if (!ok) continue;
            if      (nbr[(size_t)sb*2]   == -1) nbr[(size_t)sb*2]   = sa;
            else if (nbr[(size_t)sb*2+1] == -1) nbr[(size_t)sb*2+1] = sa;
            else { /* undo the sa link */
                if (nbr[(size_t)sa*2+1] == sb) nbr[(size_t)sa*2+1] = -1;
                else if (nbr[(size_t)sa*2] == sb) nbr[(size_t)sa*2] = -1;
            }
        }
    }

    /* crossing theta (raw, wrapped) for orientation + loop cutting */
    double *cth = (double *)ARENA_ALLOC(arena, (long)(total * sizeof(double)));
    for (size_t ci = 0; ci < total; ci++) {
        double d0 = cpos[ci*3+0] - (double)axis_point[0];
        double d1 = cpos[ci*3+1] - (double)axis_point[1];
        double d2 = cpos[ci*3+2] - (double)axis_point[2];
        double a1 = d0*e1[0] + d1*e1[1] + d2*e1[2];
        double a2 = d0*e2[0] + d1*e2[1] + d2*e2[2];
        cth[ci] = atan2(a2, a1);
    }

    /* --- walk chains (open first, then closed loops cut at max |dtheta|) --- */
    uint8_t *vis   = (uint8_t *)ARENA_CALLOC(arena, (long)total, 1);
    int32_t *order = (int32_t *)ARENA_ALLOC(arena, (long)(total * sizeof(int32_t)));
    size_t maxch = total / 2 + 2;
    Chain *chains = (Chain *)ARENA_ALLOC(arena, (long)(maxch * sizeof(Chain)));
    size_t n_chn = 0, opos = 0;

    for (int pass = 0; pass < 2 && n_chn < maxch; pass++) {
        for (size_t s0 = 0; s0 < total && n_chn < maxch; s0++) {
            if (vis[s0]) continue;
            int deg = (nbr[s0*2] != -1) + (nbr[s0*2+1] != -1);
            if (pass == 0 ? (deg != 1) : (deg != 2)) continue;
            size_t start = opos;
            int32_t cur = (int32_t)s0, prev = -1;
            for (;;) {
                order[opos++] = cur;
                vis[cur] = 1;
                int32_t n0 = nbr[(size_t)cur*2], n1 = nbr[(size_t)cur*2+1];
                int32_t nxt = (n0 != prev && n0 != -1) ? n0 : n1;
                if (nxt == -1 || nxt == prev || vis[nxt]) break;
                prev = cur; cur = nxt;
            }
            size_t len = opos - start;
            if (len < 2) { opos = start; continue; }
            if (pass == 1) {
                /* closed loop: rotate so the largest wrapped |dtheta| step
                 * (incl. the wrap-around edge) sits between last and first */
                size_t cutat = 0; double best = -1.0;
                for (size_t i = 0; i < len; i++) {
                    size_t jn = (i + 1) % len;
                    double d = fabs(wrap_pi(cth[order[start+jn]] - cth[order[start+i]]));
                    if (d > best) { best = d; cutat = jn; }
                }
                if (cutat != 0) {
                    Arena_Mark rm = Arena_save(arena);
                    int32_t *tmp = (int32_t *)ARENA_ALLOC(arena, (long)(len * sizeof(int32_t)));
                    for (size_t i = 0; i < len; i++)
                        tmp[i] = order[start + (cutat + i) % len];
                    memcpy(&order[start], tmp, len * sizeof(int32_t));
                    Arena_restore(arena, rm);
                }
                S->n_closed++;
            }
            chains[n_chn].first  = (int32_t)start;
            chains[n_chn].count  = (int32_t)len;
            chains[n_chn].slice  = cpln[order[start]];
            chains[n_chn].closed = (pass == 1);
            n_chn++;
        }
    }

    /* --- chain lengths + net winding (pass 1) --- */
    double *chlen = (double *)ARENA_ALLOC(arena, (long)((n_chn + 1) * sizeof(double)));
    double *chwind = (double *)ARENA_ALLOC(arena, (long)((n_chn + 1) * sizeof(double)));
    double sh = (double)o->sample_h;
    size_t cap = 0;
    for (size_t c = 0; c < n_chn; c++) {
        int32_t len = chains[c].count, first = chains[c].first;
        double tot = 0.0, wnd = 0.0;
        for (int32_t i = 1; i < len; i++) {
            const double *pa = &cpos[(size_t)order[first + i - 1] * 3];
            const double *pb = &cpos[(size_t)order[first + i] * 3];
            double dz = pb[0]-pa[0], dy = pb[1]-pa[1], dx = pb[2]-pa[2];
            tot += sqrt(dz*dz + dy*dy + dx*dx);
            wnd += wrap_pi(cth[order[first + i]] - cth[order[first + i - 1]]);
        }
        chlen[c] = tot;
        chwind[c] = wnd;
        if (tot >= RIB_MIN_CHAIN_LEN)
            cap += (size_t)(tot / sh + 0.5) + 2;
    }
    if (cap < 2) return -1;

    /* --- orient (ascending theta) + resample (pass 2) --- */
    Sample *smp = (Sample *)ARENA_ALLOC(arena, (long)(cap * sizeof(Sample)));
    Chain  *chn0 = (Chain *)ARENA_ALLOC(arena, (long)((n_chn + 1) * sizeof(Chain)));
    size_t n_smp = 0, n_pre_chn = 0;

    for (size_t c = 0; c < n_chn; c++) {
        int32_t len = chains[c].count;
        int32_t first = chains[c].first;
        double totlen = chlen[c];
        if (totlen < RIB_MIN_CHAIN_LEN) continue;
        int rev = chwind[c] < 0.0;
        int nout = (int)(totlen / sh + 0.5) + 1;
        if (nout < 2) nout = 2;

        size_t base = n_smp;
        double step = totlen / (double)(nout - 1);
        double acc = 0.0;
        int32_t seg = 1;
        for (int j = 0; j < nout; j++) {
            double want = (double)j * step;
            if (want > totlen) want = totlen;
            for (;;) {
                int32_t ia = rev ? (first + len - seg)     : (first + seg - 1);
                int32_t ib = rev ? (first + len - seg - 1) : (first + seg);
                const double *pa = &cpos[(size_t)order[ia] * 3];
                const double *pb = &cpos[(size_t)order[ib] * 3];
                double dz = pb[0]-pa[0], dy = pb[1]-pa[1], dx = pb[2]-pa[2];
                double l = sqrt(dz*dz + dy*dy + dx*dx);
                if (acc + l >= want - 1e-12 || seg == len - 1) {
                    double w = l > 1e-12 ? (want - acc) / l : 0.0;
                    if (w < 0.0) w = 0.0;
                    if (w > 1.0) w = 1.0;
                    Sample *sp = &smp[n_smp];
                    memset(sp, 0, sizeof(*sp));
                    for (int d = 0; d < 3; d++)
                        sp->p[d] = (1.0-w)*pa[d] + w*pb[d];
                    if (cphi != NULL)
                        sp->phi = (1.0-w) * cphi[order[ia]]
                                +        w  * cphi[order[ib]];
                    sp->s = want;
                    sp->chain = (int32_t)n_pre_chn;
                    sp->slice = chains[c].slice;
                    n_smp++;
                    break;
                }
                acc += l;
                seg++;
            }
        }
        /* tangents (central differences), in-plane coords, unwound theta */
        for (size_t j = base; j < n_smp; j++) {
            size_t ja = j > base ? j - 1 : j;
            size_t jb = j + 1 < n_smp ? j + 1 : j;
            Sample *sp = &smp[j];
            sp->tau[0] = smp[jb].p[0] - smp[ja].p[0];
            sp->tau[1] = smp[jb].p[1] - smp[ja].p[1];
            sp->tau[2] = smp[jb].p[2] - smp[ja].p[2];
            v3norm(sp->tau);
            double d0 = sp->p[0] - (double)axis_point[0];
            double d1 = sp->p[1] - (double)axis_point[1];
            double d2 = sp->p[2] - (double)axis_point[2];
            sp->c1 = d0*e1[0] + d1*e1[1] + d2*e1[2];
            sp->c2 = d0*e2[0] + d1*e2[1] + d2*e2[2];
            sp->r  = sqrt(sp->c1*sp->c1 + sp->c2*sp->c2);
            double raw = atan2(sp->c2, sp->c1);
            sp->th = (j == base) ? raw
                                 : smp[j-1].th + wrap_pi(raw - (j > base ?
                                       atan2(smp[j-1].c2, smp[j-1].c1) : raw));
        }
        chn0[n_pre_chn].first  = (int32_t)base;
        chn0[n_pre_chn].count  = (int32_t)(n_smp - base);
        chn0[n_pre_chn].slice  = chains[c].slice;
        chn0[n_pre_chn].closed = chains[c].closed;
        chn0[n_pre_chn].wind   = 0;
        chn0[n_pre_chn].group  = -1;
        n_pre_chn++;
    }

    /* --- bridge cut: split chains at radial jumps a spiral ramp cannot make
     * (fusion walls between wraps); samples on sub-min fragments are orphaned */
    double bridge_dr = RIB_BRIDGE_FRAC * sh;
    Chain *chn = (Chain *)ARENA_ALLOC(arena, (long)((n_smp / 2 + 2) * sizeof(Chain)));
    size_t n_out_chn = 0;
    size_t cuts = 0;
    for (size_t c = 0; c < n_pre_chn; c++) {
        int32_t f = chn0[c].first, cnt = chn0[c].count;
        int32_t segstart = 0;
        for (int32_t i = 1; i <= cnt; i++) {
            int cut = (i == cnt) ||
                      (fabs(smp[f+i].r - smp[f+i-1].r) > bridge_dr);
            if (!cut) continue;
            if (i < cnt) cuts++;
            int32_t seglen = i - segstart;
            if (seglen >= 2) {
                chn[n_out_chn].first  = f + segstart;
                chn[n_out_chn].count  = seglen;
                chn[n_out_chn].slice  = chn0[c].slice;
                chn[n_out_chn].closed = chn0[c].closed;
                chn[n_out_chn].wind   = 0;
                chn[n_out_chn].group  = -1;
                for (int32_t q = segstart; q < i; q++)
                    smp[f+q].chain = (int32_t)n_out_chn;
                n_out_chn++;
            } else {
                for (int32_t q = segstart; q < i; q++)
                    smp[f+q].chain = -1;   /* orphan */
            }
            segstart = i;
        }
    }
    S->bridge_cuts = cuts;

    /* re-orient each FINAL chain by its own net winding (a fragment of a
     * fold-back or wall-adjacent chain can run theta-descending even when its
     * parent chain net-ascended; PAVA assumes u ascends along sample order) */
    for (size_t c = 0; c < n_out_chn; c++) {
        int32_t f = chn[c].first, cnt = chn[c].count;
        if (cnt < 2) continue;
        if (smp[f + cnt - 1].th >= smp[f].th) continue;
        for (int32_t i = 0; i < cnt / 2; i++) {
            Sample tmp = smp[f + i];
            smp[f + i] = smp[f + cnt - 1 - i];
            smp[f + cnt - 1 - i] = tmp;
        }
        smp[f].s = 0.0;
        for (int32_t i = 1; i < cnt; i++) {
            double dz = smp[f+i].p[0] - smp[f+i-1].p[0];
            double dy = smp[f+i].p[1] - smp[f+i-1].p[1];
            double dx = smp[f+i].p[2] - smp[f+i-1].p[2];
            smp[f+i].s = smp[f+i-1].s + sqrt(dz*dz + dy*dy + dx*dx);
        }
        for (int32_t i = 0; i < cnt; i++) {
            smp[f+i].tau[0] = -smp[f+i].tau[0];
            smp[f+i].tau[1] = -smp[f+i].tau[1];
            smp[f+i].tau[2] = -smp[f+i].tau[2];
            smp[f+i].chain = (int32_t)c;   /* unchanged, but keep exact */
        }
    }

    S->smp = smp;  S->n_smp = n_smp;
    S->chn = chn;  S->n_chn = n_out_chn;

    /* slice occupancy stats */
    Arena_Mark m2 = Arena_save(arena);
    int32_t *per = (int32_t *)ARENA_CALLOC(arena, (long)nplanes, (long)sizeof(int32_t));
    for (size_t c = 0; c < n_out_chn; c++) per[chn[c].slice]++;
    for (int k = 0; k < nplanes; k++) {
        if (per[k] >= 1) S->n_slices_hit++;
        if (per[k] > 1)  S->n_multi++;
    }
    Arena_restore(arena, m2);
    return (n_smp >= 2 && n_out_chn >= 1) ? 0 : -1;
}

/* ============================================================================
 * Stage B -- alignment pair candidates (winding-index votes come from these).
 * ==========================================================================*/

typedef struct {
    int32_t a, b;     /* sample ids */
    int32_t k;        /* winding relation: W[chain_b] - W[chain_a] */
    double d;         /* tangent-projected rest offset (Eq. 7) */
    double w;         /* base weight */
    double like;      /* reweighted likelihood */
} Pair;

typedef struct {
    Pair  *pairs;  size_t n_pairs;
    size_t n_cross, n_cont;
    size_t matched_a;
    size_t candidates_a;
} PairSet;

/* counting-sort sample ids by slice -> off[] (nplanes+1); orphans excluded */
static int32_t *sort_by_slice(Arena_T arena, const SliceSet *S, int32_t **out_off)
{
    int np = S->nplanes;
    int32_t *off = (int32_t *)ARENA_CALLOC(arena, (long)(np + 1), (long)sizeof(int32_t));
    for (size_t i = 0; i < S->n_smp; i++)
        if (S->smp[i].chain >= 0) off[S->smp[i].slice + 1]++;
    for (int k = 0; k < np; k++) off[k+1] = (int32_t)(off[k+1] + off[k]);
    int32_t *ids = (int32_t *)ARENA_ALLOC(arena, (long)((S->n_smp + 1) * sizeof(int32_t)));
    int32_t *cur = (int32_t *)ARENA_ALLOC(arena, (long)((size_t)np * sizeof(int32_t)));
    memcpy(cur, off, (size_t)np * sizeof(int32_t));
    for (size_t i = 0; i < S->n_smp; i++)
        if (S->smp[i].chain >= 0)
            ids[cur[S->smp[i].slice]++] = (int32_t)i;
    *out_off = off;
    return ids;
}

static int32_t *chains_by_slice(Arena_T arena, const SliceSet *S, int32_t **out_off)
{
    int np = S->nplanes;
    int32_t *off = (int32_t *)ARENA_CALLOC(arena, (long)(np + 1), (long)sizeof(int32_t));
    for (size_t c = 0; c < S->n_chn; c++) off[S->chn[c].slice + 1]++;
    for (int k = 0; k < np; k++) off[k+1] = (int32_t)(off[k+1] + off[k]);
    int32_t *ids = (int32_t *)ARENA_ALLOC(arena, (long)((S->n_chn + 1) * sizeof(int32_t)));
    int32_t *cur = (int32_t *)ARENA_ALLOC(arena, (long)((size_t)np * sizeof(int32_t)));
    memcpy(cur, off, (size_t)np * sizeof(int32_t));
    for (size_t c = 0; c < S->n_chn; c++)
        ids[cur[S->chn[c].slice]++] = (int32_t)c;
    *out_off = off;
    return ids;
}

static int build_pairs(Arena_T arena, const SliceSet *S,
                       const RibbonOpts *o, PairSet *P)
{
    memset(P, 0, sizeof(*P));
    int np = S->nplanes;
    const Sample *smp = S->smp;

    int32_t *off = NULL;
    int32_t *ids = sort_by_slice(arena, S, &off);

    double c1min = 1e300, c1max = -1e300, c2min = 1e300, c2max = -1e300;
    for (size_t i = 0; i < S->n_smp; i++) {
        if (smp[i].chain < 0) continue;
        if (smp[i].c1 < c1min) c1min = smp[i].c1;
        if (smp[i].c1 > c1max) c1max = smp[i].c1;
        if (smp[i].c2 < c2min) c2min = smp[i].c2;
        if (smp[i].c2 > c2max) c2max = smp[i].c2;
    }
    double cell = (double)o->match_r;
    if (cell < 1e-6) cell = 1.0;
    int g1 = (int)((c1max - c1min) / cell) + 1;
    int g2 = (int)((c2max - c2min) / cell) + 1;
    if (g1 < 1) g1 = 1;
    if (g2 < 1) g2 = 1;
    size_t ncell = (size_t)g1 * (size_t)g2;

    double cos_gate = cos((double)o->match_ang_deg * M_PI / 180.0);
    double r2max = (double)o->match_r * (double)o->match_r;

    Pair *pairs = (Pair *)ARENA_ALLOC(arena, (long)((S->n_smp + 16) * sizeof(Pair)));
    size_t n_pairs = 0;

    int32_t *goff = (int32_t *)ARENA_ALLOC(arena, (long)((ncell + 1) * sizeof(int32_t)));
    int32_t *gids = (int32_t *)ARENA_ALLOC(arena, (long)((S->n_smp + 1) * sizeof(int32_t)));
    int32_t *gcur = (int32_t *)ARENA_ALLOC(arena, (long)((ncell + 1) * sizeof(int32_t)));
    uint8_t *matched = (uint8_t *)ARENA_CALLOC(arena, (long)(S->n_smp + 1), 1);

    /* stride escalation: a sample unmatched at k+1 tries k+2, then k+3. An
     * EMPTY slice (hole band, or a plane coinciding exactly with a vertex
     * row) must not sever the pair graph vertically -- one missing slice
     * would otherwise shatter a coherent sheet into stacked fragments. */
    for (int stride = 1; stride <= RIB_VJUMP_MAX; stride++) {
        for (int k = 0; k + stride < np; k++) {
            int32_t a0 = off[k],        a1 = off[k+1];
            int32_t b0 = off[k+stride], b1 = off[k+stride+1];
            if (a1 == a0 || b1 == b0) continue;
            memset(goff, 0, (ncell + 1) * sizeof(int32_t));
            for (int32_t q = b0; q < b1; q++) {
                const Sample *sb = &smp[ids[q]];
                int i1 = (int)((sb->c1 - c1min) / cell), i2 = (int)((sb->c2 - c2min) / cell);
                if (i1 < 0) i1 = 0;
                if (i1 > g1-1) i1 = g1-1;
                if (i2 < 0) i2 = 0;
                if (i2 > g2-1) i2 = g2-1;
                goff[(size_t)i1 * (size_t)g2 + (size_t)i2 + 1]++;
            }
            for (size_t cix = 0; cix < ncell; cix++)
                goff[cix+1] = (int32_t)(goff[cix+1] + goff[cix]);
            memcpy(gcur, goff, ncell * sizeof(int32_t));
            for (int32_t q = b0; q < b1; q++) {
                const Sample *sb = &smp[ids[q]];
                int i1 = (int)((sb->c1 - c1min) / cell), i2 = (int)((sb->c2 - c2min) / cell);
                if (i1 < 0) i1 = 0;
                if (i1 > g1-1) i1 = g1-1;
                if (i2 < 0) i2 = 0;
                if (i2 > g2-1) i2 = g2-1;
                gids[gcur[(size_t)i1 * (size_t)g2 + (size_t)i2]++] = ids[q];
            }
            for (int32_t q = a0; q < a1; q++) {
                int32_t ia = ids[q];
                if (stride == 1) P->candidates_a++;
                if (matched[ia]) continue;
                const Sample *sa = &smp[ia];
                int i1 = (int)((sa->c1 - c1min) / cell), i2 = (int)((sa->c2 - c2min) / cell);
                int32_t bestb = -1; double bestd2 = r2max;
                for (int d1 = -1; d1 <= 1; d1++) for (int d2i = -1; d2i <= 1; d2i++) {
                    int j1 = i1 + d1, j2 = i2 + d2i;
                    if (j1 < 0 || j2 < 0 || j1 >= g1 || j2 >= g2) continue;
                    size_t cix = (size_t)j1 * (size_t)g2 + (size_t)j2;
                    for (int32_t r = goff[cix]; r < goff[cix+1]; r++) {
                        int32_t ib = gids[r];
                        const Sample *sb = &smp[ib];
                        double dd1 = sb->c1 - sa->c1, dd2 = sb->c2 - sa->c2;
                        double d2 = dd1*dd1 + dd2*dd2;
                        if (d2 >= bestd2) continue;
                        if (fabs(wrap_pi(sb->th - sa->th)) > RIB_CONT_DPHI_MAX) continue;
                        if (v3dot(sa->tau, sb->tau) < cos_gate) continue;
                        bestd2 = d2; bestb = ib;
                    }
                }
                if (bestb >= 0) {
                    const Sample *sb = &smp[bestb];
                    double tb[3] = { sa->tau[0] + sb->tau[0],
                                     sa->tau[1] + sb->tau[1],
                                     sa->tau[2] + sb->tau[2] };
                    v3norm(tb);
                    double dp[3] = { sa->p[0]-sb->p[0], sa->p[1]-sb->p[1], sa->p[2]-sb->p[2] };
                    pairs[n_pairs].a = ia;
                    pairs[n_pairs].b = bestb;
                    pairs[n_pairs].k = (int32_t)lround((sa->th - sb->th) / (2.0 * M_PI));
                    pairs[n_pairs].d = v3dot(tb, dp);
                    pairs[n_pairs].w = RIB_ALIGN_W;
                    pairs[n_pairs].like = 1.0;
                    n_pairs++;
                    P->matched_a++;
                    matched[ia] = 1;
                }
            }
        }
    }
    P->n_cross = n_pairs;

    /* --- intra-slice fragment continuation pairs --- */
    size_t cont_cap = S->n_chn + 16;
    Pair *cont = (Pair *)ARENA_ALLOC(arena, (long)(cont_cap * sizeof(Pair)));
    size_t n_cont = 0;
    int32_t *choff = NULL;
    int32_t *chids = chains_by_slice(arena, S, &choff);

    for (int k = 0; k < np; k++) {
        for (int32_t x = choff[k]; x < choff[k+1]; x++) {
            const Chain *ce = &S->chn[chids[x]];
            int32_t iE = ce->first + ce->count - 1;
            const Sample *sE = &smp[iE];
            int32_t bestS = -1; double bestd = RIB_CONT_GAP_VOX;
            for (int32_t y = choff[k]; y < choff[k+1]; y++) {
                if (x == y) continue;
                const Chain *cs = &S->chn[chids[y]];
                int32_t iS = cs->first;
                const Sample *sS = &smp[iS];
                double dth = wrap_pi(sS->th - sE->th);
                if (dth <= 0.0 || dth > RIB_CONT_DPHI_MAX) continue;
                if (fabs(sS->r - sE->r) > RIB_CONT_DR_MAX) continue;
                double dd1 = sS->c1 - sE->c1, dd2 = sS->c2 - sE->c2;
                double d = sqrt(dd1*dd1 + dd2*dd2);
                if (d < bestd) { bestd = d; bestS = iS; }
            }
            if (bestS >= 0 && n_cont < cont_cap) {
                const Sample *sS = &smp[bestS];
                double tb[3] = { sE->tau[0] + sS->tau[0],
                                 sE->tau[1] + sS->tau[1],
                                 sE->tau[2] + sS->tau[2] };
                v3norm(tb);
                double dp[3] = { sE->p[0]-sS->p[0], sE->p[1]-sS->p[1], sE->p[2]-sS->p[2] };
                cont[n_cont].a = iE;
                cont[n_cont].b = bestS;
                /* winding relation from the wrapped local step (the raw th
                 * difference of two chains is offset by their unknown 2pi*W) */
                cont[n_cont].k = (int32_t)lround((sE->th + wrap_pi(sS->th - sE->th)
                                                  - sS->th) / (2.0 * M_PI));
                cont[n_cont].d = v3dot(tb, dp);
                cont[n_cont].w = RIB_ALIGN_W;
                cont[n_cont].like = 1.0;
                n_cont++;
            }
        }
    }
    P->n_cont = n_cont;

    Pair *all = (Pair *)ARENA_ALLOC(arena, (long)((n_pairs + n_cont + 1) * sizeof(Pair)));
    memcpy(all, pairs, n_pairs * sizeof(Pair));
    memcpy(all + n_pairs, cont, n_cont * sizeof(Pair));
    P->pairs = all;
    P->n_pairs = n_pairs + n_cont;
    return 0;
}

/* ============================================================================
 * Stage W -- integer winding index per chain, anchored GLOBALLY by radius vs.
 * the umbilicus (the diffeomorphic-spiral step; no cross-slice vote graph).
 *
 * A scroll is ONE continuous strip and its number of wraps equals the maximal
 * winding number. The previous scheme reconstructed W from cross-slice pair
 * votes and placed vote-disconnected chain groups by a radial-ordering BFS;
 * any group the vote graph could not reach fell to baseW=0 and collapsed onto
 * the anchor, cramming ~30 physical wraps into ~6 u-bands (with empty gaps).
 *
 * Instead anchor every chain directly by radius. On a spiral r = a + b*psi
 * (psi = turn number, b = pitch) the per-sample residual
 *     res_i = r_i/pitch - th_i/2pi
 * is INVARIANT along a chain: for a connected multi-turn chain th (already
 * unwound) carries the turns so res stays ~a/pitch; for a disconnected wrap
 * (a concentric ring / detached fragment) res spans one turn and its median
 * lands mid-wrap. Either way
 *     W[c] = round( median_i( r_i/pitch - th_i/2pi ) )
 * recovers the true relative winding -- adjacent wraps differ by exactly 1 --
 * outlier-robust, O(n), with no groups / votes / BFS. The global a/pitch
 * offset is an irrelevant gauge (fixed later by canonicalization).
 *
 * Delamination trapping falls out: a flap within ~half a pitch of its parent
 * rounds to the SAME W (same UV band); a flap ~one pitch out becomes the next
 * wrap (a genuine ambiguity, deferred to ink). phi* = th + 2pi*W.
 *
 * pitch: caller --wrap-spacing override, else a two-pass median of adjacent
 * radial gaps within (slice, theta-wedge) bins (RIB_PITCH_* clamp/fallback).
 * ==========================================================================*/

/* a referenced sample, keyed for per-wedge radial gap estimation */
typedef struct { int32_t wedge; float r; } RadSample;
static int cmp_radsample(const void *pa, const void *pb)
{
    const RadSample *a = (const RadSample *)pa, *b = (const RadSample *)pb;
    if (a->wedge != b->wedge) return a->wedge < b->wedge ? -1 : 1;
    return a->r < b->r ? -1 : (a->r > b->r ? 1 : 0);
}

/* Radial pitch (vox/turn): two-pass median of adjacent radial gaps within each
 * (slice, theta-wedge) bin. The second pass drops multi-wrap jumps that
 * contaminate a single-pass median when the spacing varies. Returns 0 if it
 * cannot estimate (caller falls back to RIB_PITCH_DEFAULT). */
static double estimate_pitch(Arena_T arena, const SliceSet *S)
{
    Arena_Mark mp = Arena_save(arena);
    double pitch = 0.0;
    size_t nref = 0;
    for (size_t i = 0; i < S->n_smp; i++) if (S->smp[i].chain >= 0) nref++;
    if (nref >= 8) {
        RadSample *rs = (RadSample *)ARENA_ALLOC(arena, (long)((nref + 1) * sizeof(RadSample)));
        size_t q = 0;
        for (size_t i = 0; i < S->n_smp; i++) {
            const Sample *sp = &S->smp[i];
            if (sp->chain < 0) continue;
            double raw = atan2(sp->c2, sp->c1);
            int tb = (int)((raw + M_PI) / (2.0 * M_PI) * RIB_WEDGE_NTHETA);
            if (tb < 0) tb = 0;
            if (tb >= RIB_WEDGE_NTHETA) tb = RIB_WEDGE_NTHETA - 1;
            rs[q].wedge = sp->slice * RIB_WEDGE_NTHETA + tb;
            rs[q].r = (float)sp->r;
            q++;
        }
        qsort(rs, nref, sizeof(RadSample), cmp_radsample);
        double *gaps = (double *)ARENA_ALLOC(arena, (long)((nref + 1) * sizeof(double)));
        size_t ng = 0;
        for (size_t i = 1; i < nref; i++) {
            if (rs[i].wedge != rs[i-1].wedge) continue;
            double dr = (double)rs[i].r - (double)rs[i-1].r;
            if (dr < RIB_PITCH_MIN || dr > RIB_PITCH_MAX) continue;
            gaps[ng++] = dr;
        }
        if (ng >= 8) {
            qsort(gaps, ng, sizeof(double), cmp_dbl);
            double m0 = gaps[ng / 2];
            size_t ng2 = 0;
            for (size_t i = 0; i < ng; i++) if (gaps[i] <= 1.6 * m0) gaps[ng2++] = gaps[i];
            pitch = ng2 >= 8 ? gaps[ng2 / 2] : m0;
        }
    }
    Arena_restore(arena, mp);
    return pitch;
}

/* a chain-pair winding vote (lo,hi ordered), k = W[hi]-W[lo], with tally */
typedef struct { int32_t lo, hi; int32_t k; int32_t votes; } ChainEdge;
static int cmp_chainedge(const void *pa, const void *pb)
{
    const ChainEdge *a = (const ChainEdge *)pa, *b = (const ChainEdge *)pb;
    if (a->lo != b->lo) return a->lo < b->lo ? -1 : 1;
    if (a->hi != b->hi) return a->hi < b->hi ? -1 : 1;
    if (a->k  != b->k)  return a->k  < b->k  ? -1 : 1;
    return 0;
}

static void assign_winding(Arena_T arena, SliceSet *S, const PairSet *P,
                           const RibbonOpts *o, RibbonResult *out)
{
    size_t nc = S->n_chn;
    if (nc == 0) return;
    Arena_Mark mark = Arena_save(arena);

    /* --- radial pitch (vox/turn) --- */
    double pitch = (double)o->wrap_spacing;
    int pitch_source = RIB_PITCH_PINNED;
    if (pitch <= 0.0) {
        pitch = estimate_pitch(arena, S);
        pitch_source = pitch >= RIB_PITCH_MIN
                     ? RIB_PITCH_ESTIMATED : RIB_PITCH_FALLBACK;
    } else if (pitch < RIB_PITCH_MIN) {
        pitch_source = RIB_PITCH_FALLBACK;
    }
    if (pitch < RIB_PITCH_MIN) pitch = RIB_PITCH_DEFAULT;
    double inv2pi = 1.0 / (2.0 * M_PI);

    /* winding SENSE: sign of within-chain cov(th, r), normalized. The radius
     * anchor below uses sense*r/pitch; without the correct sign a negatively-
     * wound scroll (phi decreasing with radius) anchors detached groups a
     * couple wraps off. Robust near zero (rings have r ~constant, corr ~0):
     * default +1, flip to -1 only on a CLEAR negative correlation. */
    double s_thr = 0.0, s_th2 = 0.0, s_r2 = 0.0;
    for (size_t c = 0; c < nc; c++) {
        int32_t f = S->chn[c].first, n = S->chn[c].count;
        if (n < 4) continue;
        double mth = 0.0, mr = 0.0;
        for (int32_t i = 0; i < n; i++) { mth += S->smp[f+i].th; mr += S->smp[f+i].r; }
        mth /= (double)n; mr /= (double)n;
        for (int32_t i = 0; i < n; i++) {
            double dt = S->smp[f+i].th - mth, dr = S->smp[f+i].r - mr;
            s_thr += dt * dr; s_th2 += dt * dt; s_r2 += dr * dr;
        }
    }
    double denom = sqrt(s_th2 * s_r2);
    double sense = (denom > 1e-12 && s_thr / denom < -0.3) ? -1.0 : 1.0;
    /* pinned sense (scroll_whole): a small per-cube patch can have too little
     * radial travel for the covariance to be decisive, and a wrong sense
     * anchors detached groups a couple wraps off -- calibrate once, pin
     * everywhere so every cube of one scroll agrees. */
    if (o->winding_sense != 0) sense = o->winding_sense > 0 ? 1.0 : -1.0;

    /* ============== relative winding via the WTS pair graph ==============
     * Rounding each chain's winding from its own radius INDEPENDENTLY makes
     * neighbouring chains on ONE physical wrap disagree by +-1 wherever the
     * radius sits near a wrap boundary: u then jumps a whole circumference
     * between them and the unroll shatters into streaks. Instead PROPAGATE a
     * smooth relative winding Wrel through the pair graph -- each pair votes
     * W[b]-W[a] = round((th_a - th_b)/2pi) -- so every pair-connected chain
     * shares one consistent winding. Radius then only anchors each connected
     * group's integer OFFSET (next block). Smooth within a group, radius-placed
     * across groups => a continuous sheet with the true wrap count. */
    ChainEdge *ev = (ChainEdge *)ARENA_ALLOC(arena, (long)((P->n_pairs + 1) * sizeof(ChainEdge)));
    size_t nev = 0;
    for (size_t p = 0; p < P->n_pairs; p++) {
        int32_t ca = S->smp[P->pairs[p].a].chain;
        int32_t cb = S->smp[P->pairs[p].b].chain;
        if (ca < 0 || cb < 0 || ca == cb) continue;
        int32_t k = P->pairs[p].k;                 /* W[cb] - W[ca] */
        if (ca < cb) { ev[nev].lo = ca; ev[nev].hi = cb; ev[nev].k = k; }
        else         { ev[nev].lo = cb; ev[nev].hi = ca; ev[nev].k = -k; }
        ev[nev].votes = 1; nev++;
    }
    qsort(ev, nev, sizeof(ChainEdge), cmp_chainedge);
    size_t ne = 0;                                 /* coalesce identical (lo,hi,k) */
    for (size_t i = 0; i < nev; ) {
        size_t j = i;
        while (j < nev && ev[j].lo==ev[i].lo && ev[j].hi==ev[i].hi && ev[j].k==ev[i].k) j++;
        ev[ne] = ev[i]; ev[ne].votes = (int32_t)(j - i); ne++; i = j;
    }
    size_t nmaj = 0, conflicts = 0;                /* majority k per (lo,hi) */
    for (size_t i = 0; i < ne; ) {
        size_t j = i, best = i;
        while (j < ne && ev[j].lo==ev[i].lo && ev[j].hi==ev[i].hi) {
            if (ev[j].votes > ev[best].votes) best = j; j++;
        }
        for (size_t q = i; q < j; q++) if (q != best) conflicts += (size_t)ev[q].votes;
        ev[nmaj++] = ev[best]; i = j;
    }
    int32_t *deg = (int32_t *)ARENA_CALLOC(arena, (long)(nc + 1), (long)sizeof(int32_t));
    for (size_t i = 0; i < nmaj; i++) { deg[ev[i].lo+1]++; deg[ev[i].hi+1]++; }
    for (size_t c = 0; c < nc; c++) deg[c+1] = (int32_t)(deg[c+1] + deg[c]);
    int32_t *adj = (int32_t *)ARENA_ALLOC(arena, (long)(((size_t)deg[nc] + 1) * sizeof(int32_t)));
    int32_t *adk = (int32_t *)ARENA_ALLOC(arena, (long)(((size_t)deg[nc] + 1) * sizeof(int32_t)));
    int32_t *cur = (int32_t *)ARENA_ALLOC(arena, (long)(nc * sizeof(int32_t)));
    memcpy(cur, deg, nc * sizeof(int32_t));
    for (size_t i = 0; i < nmaj; i++) {
        adj[cur[ev[i].lo]] = ev[i].hi; adk[cur[ev[i].lo]++] = ev[i].k;
        adj[cur[ev[i].hi]] = ev[i].lo; adk[cur[ev[i].hi]++] = -ev[i].k;
    }
    int32_t *grp = (int32_t *)ARENA_ALLOC(arena, (long)(nc * sizeof(int32_t)));
    for (size_t c = 0; c < nc; c++) grp[c] = -1;
    int32_t *queue = (int32_t *)ARENA_ALLOC(arena, (long)(nc * sizeof(int32_t)));
    int32_t ngroups = 0;
    for (size_t s0 = 0; s0 < nc; s0++) {
        if (grp[s0] != -1) continue;
        grp[s0] = ngroups; S->chn[s0].wind = 0;
        size_t head = 0, tail = 0; queue[tail++] = (int32_t)s0;
        while (head < tail) {
            int32_t u = queue[head++];
            for (int32_t e = deg[u]; e < deg[u+1]; e++) {
                int32_t v = adj[e];
                int32_t w = S->chn[u].wind + adk[e];
                if (grp[v] == -1) { grp[v] = ngroups; S->chn[v].wind = w; queue[tail++] = v; }
                else if (S->chn[v].wind != w) conflicts++;
            }
        }
        ngroups++;
    }
    for (size_t c = 0; c < nc; c++) S->chn[c].group = grp[c];

    /* --- anchor each group's integer offset by RADIUS ---
     * base[g] = round(mean over the group's samples of
     *                 (sense*r/pitch - th/2pi - Wrel));  W = Wrel + base.
     * The smooth Wrel keeps one wrap contiguous; radius places each group on
     * its true turn, so 30 wraps stay 30 wraps (no unreached->0 collapse). */
    double *baseacc = (double *)ARENA_CALLOC(arena, (long)(ngroups ? ngroups : 1), (long)sizeof(double));
    double *basecnt = (double *)ARENA_CALLOC(arena, (long)(ngroups ? ngroups : 1), (long)sizeof(double));
    for (size_t i = 0; i < S->n_smp; i++) {
        const Sample *sp = &S->smp[i];
        if (sp->chain < 0) continue;
        int32_t g = grp[sp->chain];
        baseacc[g] += sense * sp->r / pitch - sp->th * inv2pi - (double)S->chn[sp->chain].wind;
        basecnt[g] += 1.0;
    }
    int32_t minW = 0x7fffffff, maxW = -0x7fffffff;
    for (size_t c = 0; c < nc; c++) {
        int32_t g = grp[c];
        int32_t base = (basecnt[g] > 0.0) ? (int32_t)lround(baseacc[g] / basecnt[g]) : 0;
        int32_t W = S->chn[c].wind + base;
        S->chn[c].wind = W;
        if (W < minW) minW = W;
        if (W > maxW) maxW = W;
    }
    if (minW > maxW) { minW = 0; maxW = 0; }
    long ambiguous = 0;

    /* final phi per sample. With a graph-traced reference, slice_mesh already
     * interpolated the authoritative absolute field onto every sample. Do not
     * overwrite it with the fragile cross-slice nearest-pair reconstruction:
     * contacting plies can be < match_r apart even though they differ by a
     * whole turn. StrokeStrip consumes the reference only as a turn scaffold;
     * Stage C below still solves its own true arc-length u. */
    if (o->reference_phi == NULL) {
        for (size_t i = 0; i < S->n_smp; i++) {
            Sample *sp = &S->smp[i];
            sp->phi = sp->th + (sp->chain >= 0
                                ? 2.0 * M_PI * (double)S->chn[sp->chain].wind : 0.0);
        }
    } else {
        double pmin = 1e300, pmax = -1e300;
        for (size_t i = 0; i < S->n_smp; i++) {
            const Sample *sp = &S->smp[i];
            if (sp->chain < 0) continue;
            if (sp->phi < pmin) pmin = sp->phi;
            if (sp->phi > pmax) pmax = sp->phi;
        }
        if (pmin <= pmax) {
            minW = (int32_t)floor(pmin / (2.0 * M_PI));
            maxW = (int32_t)ceil (pmax / (2.0 * M_PI));
        }
    }

    /* --- global spiral line r ~= a + b*(phi/2pi). b is PINNED to the robust
     * pitch: an ordinary least-squares slope gets dragged toward ~half the true
     * pitch by the occasional mis-wound fragment (its phi lands a turn off),
     * which then makes u_spiral (the solve init) mis-place every component. The
     * pitch is a robust median of radial gaps, so use it as the slope and fit
     * only the intercept a = median(r - pitch*phi/2pi). sr2 measures how well
     * radius tracks the assigned winding (low = folded core). --- */
    double sa = 0.0, sb = pitch, sr2 = -1.0;
    /* PINNED line (scroll_whole): every cube of one scroll must map phi
     * through the SAME u(phi), so (a,b) come from a one-time seed-cube
     * calibration instead of a per-mesh fit; sr2 stays a per-mesh diagnostic
     * of how well THIS mesh's radius tracks the shared line. */
    int sp_pinned = (o->spiral_b != 0.0);
    if (sp_pinned) { sa = o->spiral_a; sb = o->spiral_b; }
    {
        Arena_Mark ms = Arena_save(arena);
        size_t nn = 0;
        for (size_t i = 0; i < S->n_smp; i++) if (S->smp[i].chain >= 0) nn++;
        if (nn >= 64) {
            /* winding sense: phi may increase inward or outward depending on
             * chirality / basis handedness, so pin |slope| to the robust pitch
             * but take its SIGN from cov(phi, r) (a wrong-sign pin collapses
             * r2 and wrongly disables the spiral placement). */
            double mphi = 0.0, mr = 0.0, cov = 0.0;
            for (size_t i = 0; i < S->n_smp; i++) {
                const Sample *sp = &S->smp[i];
                if (sp->chain < 0) continue;
                mphi += sp->phi; mr += sp->r;
            }
            mphi /= (double)nn; mr /= (double)nn;
            if (!sp_pinned) {
                for (size_t i = 0; i < S->n_smp; i++) {
                    const Sample *sp = &S->smp[i];
                    if (sp->chain < 0) continue;
                    cov += (sp->phi - mphi) * (sp->r - mr);
                }
                sb = (cov >= 0.0 ? 1.0 : -1.0) * pitch;
                double *ic = (double *)ARENA_ALLOC(arena, (nn + 1) * sizeof(double));
                size_t q = 0;
                for (size_t i = 0; i < S->n_smp; i++) {
                    const Sample *sp = &S->smp[i];
                    if (sp->chain < 0) continue;
                    ic[q++] = sp->r - sb * (sp->phi * inv2pi);
                }
                qsort(ic, q, sizeof(double), cmp_dbl);
                sa = ic[q / 2];
            }
            double st = 0.0, sr = 0.0;
            for (size_t i = 0; i < S->n_smp; i++) {
                const Sample *sp = &S->smp[i];
                if (sp->chain < 0) continue;
                double pr = sa + sb * (sp->phi * inv2pi);
                st += (sp->r - mr) * (sp->r - mr);
                sr += (sp->r - pr) * (sp->r - pr);
            }
            sr2 = st > 1e-12 ? 1.0 - sr / st : 1.0;
        }
        Arena_restore(arena, ms);
    }

    /* ambiguity (diagnostic): chains whose mean radius disagrees with the
     * spiral prediction at their phi by more than half a pitch -- a folded
     * core where radius stops tracking winding (data-limited, ink adjudicates);
     * NOT the normal branch-cut wrap, so it does not fire on clean rings. */
    for (size_t c = 0; c < nc; c++) {
        int32_t f = S->chn[c].first, n = S->chn[c].count;
        if (n < 4) continue;
        double mr = 0.0, mphi = 0.0;
        for (int32_t i = 0; i < n; i++) { mr += S->smp[f+i].r; mphi += S->smp[f+i].phi; }
        mr /= (double)n; mphi /= (double)n;
        if (fabs(mr - (sa + sb * (mphi * inv2pi))) > 0.5 * pitch) ambiguous++;
    }

    /* diagnostics (fields repurposed for the radius-anchored scheme) */
    out->w_groups       = (maxW >= minW) ? (int)(maxW - minW + 1) : 0; /* wraps found */
    out->w_prior_groups = (int)ngroups;    /* pair-connected winding groups */
    out->w_unreached    = (int)ambiguous;  /* chains where radius != winding (core) */
    out->w_conflicts    = conflicts;
    out->pitch_used     = pitch;
    out->pitch_source   = pitch_source;
    out->spiral_a       = sa;
    out->spiral_b       = sb;
    out->spiral_r2      = sr2;

    Arena_restore(arena, mark);
}

/* ============================================================================
 * Stage C -- joint parameterization solve (unchanged from the phi-field
 * version: arc-length + alignment quadratic, CG + PAVA + reweighting).
 * ==========================================================================*/

typedef struct {
    const SliceSet *S;
    const PairSet  *P;
    /* one gauge pin PER connected component of the solve graph: without its
     * own pin a disconnected component's absolute u is undetermined (every
     * energy term is a difference) and CG parks it at an arbitrary offset --
     * the source of detached charts landing mid-band. */
    const int32_t *pin_idx; const double *pin_val;
    size_t n_pins; double pin_w;
} SolveCtx;

static void apply_A(const SolveCtx *cx, const double *x, double *y)
{
    const SliceSet *S = cx->S;
    const PairSet  *P = cx->P;
    memset(y, 0, S->n_smp * sizeof(double));
    for (size_t c = 0; c < S->n_chn; c++) {
        int32_t f = S->chn[c].first, n = S->chn[c].count;
        for (int32_t i = 1; i < n; i++) {
            double l = S->smp[f+i].s - S->smp[f+i-1].s;
            if (l < 1e-9) l = 1e-9;
            double w = 1.0 / l;
            double d = x[f+i] - x[f+i-1];
            y[f+i]   += w * d;
            y[f+i-1] -= w * d;
        }
    }
    for (size_t p = 0; p < P->n_pairs; p++) {
        const Pair *pr = &P->pairs[p];
        if (pr->like <= 0.0) continue;
        double w = pr->w * pr->like;
        double d = x[pr->a] - x[pr->b];
        y[pr->a] += w * d;
        y[pr->b] -= w * d;
    }
    for (size_t p = 0; p < cx->n_pins; p++)
        y[cx->pin_idx[p]] += cx->pin_w * x[cx->pin_idx[p]];
}

static void build_b_diag(const SolveCtx *cx, double *b, double *diag)
{
    const SliceSet *S = cx->S;
    const PairSet  *P = cx->P;
    memset(b, 0, S->n_smp * sizeof(double));
    memset(diag, 0, S->n_smp * sizeof(double));
    for (size_t c = 0; c < S->n_chn; c++) {
        int32_t f = S->chn[c].first, n = S->chn[c].count;
        for (int32_t i = 1; i < n; i++) {
            double l = S->smp[f+i].s - S->smp[f+i-1].s;
            if (l < 1e-9) l = 1e-9;
            double w = 1.0 / l;
            b[f+i]   += w * l;
            b[f+i-1] -= w * l;
            diag[f+i]   += w;
            diag[f+i-1] += w;
        }
    }
    for (size_t p = 0; p < P->n_pairs; p++) {
        const Pair *pr = &P->pairs[p];
        if (pr->like <= 0.0) continue;
        double w = pr->w * pr->like;
        b[pr->a] += w * pr->d;
        b[pr->b] -= w * pr->d;
        diag[pr->a] += w;
        diag[pr->b] += w;
    }
    for (size_t p = 0; p < cx->n_pins; p++) {
        b[cx->pin_idx[p]]    += cx->pin_w * cx->pin_val[p];
        diag[cx->pin_idx[p]] += cx->pin_w;
    }
}

static int cg_solve(Arena_T arena, const SolveCtx *cx, double *u)
{
    size_t n = cx->S->n_smp;
    Arena_Mark mark = Arena_save(arena);
    double *b    = (double *)ARENA_ALLOC(arena, (long)(n * sizeof(double)));
    double *diag = (double *)ARENA_ALLOC(arena, (long)(n * sizeof(double)));
    double *r    = (double *)ARENA_ALLOC(arena, (long)(n * sizeof(double)));
    double *z    = (double *)ARENA_ALLOC(arena, (long)(n * sizeof(double)));
    double *pv   = (double *)ARENA_ALLOC(arena, (long)(n * sizeof(double)));
    double *Ap   = (double *)ARENA_ALLOC(arena, (long)(n * sizeof(double)));
    build_b_diag(cx, b, diag);
    for (size_t i = 0; i < n; i++) if (diag[i] < 1e-12) diag[i] = 1e-12;

    apply_A(cx, u, Ap);
    double bnorm = 0.0;
    for (size_t i = 0; i < n; i++) { r[i] = b[i] - Ap[i]; bnorm += b[i]*b[i]; }
    bnorm = sqrt(bnorm);
    if (bnorm < 1e-30) { Arena_restore(arena, mark); return 0; }
    double rz = 0.0;
    for (size_t i = 0; i < n; i++) { z[i] = r[i] / diag[i]; pv[i] = z[i]; rz += r[i]*z[i]; }
    int it = 0;
    for (; it < RIB_CG_MAXIT; it++) {
        double rn = 0.0;
        for (size_t i = 0; i < n; i++) rn += r[i]*r[i];
        if (sqrt(rn) / bnorm < RIB_CG_TOL) break;
        apply_A(cx, pv, Ap);
        double pAp = 0.0;
        for (size_t i = 0; i < n; i++) pAp += pv[i]*Ap[i];
        if (pAp <= 0.0) break;
        double alpha = rz / pAp;
        for (size_t i = 0; i < n; i++) { u[i] += alpha*pv[i]; r[i] -= alpha*Ap[i]; }
        double rz2 = 0.0;
        for (size_t i = 0; i < n; i++) { z[i] = r[i] / diag[i]; rz2 += r[i]*z[i]; }
        double beta = rz2 / rz;
        rz = rz2;
        for (size_t i = 0; i < n; i++) pv[i] = z[i] + beta*pv[i];
    }
    Arena_restore(arena, mark);
    return it;
}

/* Monotone repair: enforce u[i] - u[i-1] >= lb[i] via PAVA on v = u - cumsum(lb).
 * scratch must hold >= 3*n doubles. Returns #samples moved. */
static size_t pava_chain(double *u, const double *lb, int32_t n, double *scratch)
{
    if (n < 2) return 0;
    double  *v  = scratch;
    double  *bm = scratch + n;
    int32_t *bc = (int32_t *)(void *)(scratch + 2 * n);
    double csum = 0.0;
    for (int32_t i = 0; i < n; i++) {
        if (i > 0) csum += lb[i];
        v[i] = u[i] - csum;
    }
    int32_t nb = 0;
    for (int32_t i = 0; i < n; i++) {
        double m = v[i]; int32_t c = 1;
        while (nb > 0 && bm[nb-1] > m + 1e-15) {
            m = (m * (double)c + bm[nb-1] * (double)bc[nb-1]) / (double)(c + bc[nb-1]);
            c += bc[nb-1];
            nb--;
        }
        bm[nb] = m; bc[nb] = c; nb++;
    }
    size_t moved = 0;
    int32_t i = 0;
    csum = 0.0;
    for (int32_t blk = 0; blk < nb; blk++) {
        for (int32_t j = 0; j < bc[blk]; j++, i++) {
            if (i > 0) csum += lb[i];
            double nu = bm[blk] + csum;
            if (fabs(nu - u[i]) > 1e-9) moved++;
            u[i] = nu;
        }
    }
    return moved;
}

/* piecewise phi -> u map built from the MAIN solve component: sorted non-empty
 * phi bins with median u; lookup lerps, ends extrapolate with the edge slope */
typedef struct { double *bphi, *bu; int nb; double end_slope; } PhiUMap;

static double phiu_lookup(const PhiUMap *M, double phi)
{
    if (M->nb == 0) return 0.0;
    if (M->nb == 1) return M->bu[0] + (phi - M->bphi[0]) * M->end_slope;
    if (phi <= M->bphi[0])
        return M->bu[0] - (M->bphi[0] - phi)
               * ((M->bu[1] - M->bu[0]) / (M->bphi[1] - M->bphi[0]));
    if (phi >= M->bphi[M->nb-1])
        return M->bu[M->nb-1] + (phi - M->bphi[M->nb-1])
               * ((M->bu[M->nb-1] - M->bu[M->nb-2])
                  / (M->bphi[M->nb-1] - M->bphi[M->nb-2]));
    int lo = 0, hi = M->nb - 1;
    while (hi - lo > 1) {
        int mid = (lo + hi) / 2;
        if (M->bphi[mid] <= phi) lo = mid; else hi = mid;
    }
    double w = (phi - M->bphi[lo]) / (M->bphi[hi] - M->bphi[lo]);
    return (1.0 - w) * M->bu[lo] + w * M->bu[hi];
}

typedef struct { double phi, u; } PhiU2;
static int cmp_phiu2(const void *pa, const void *pb)
{
    const PhiU2 *a = (const PhiU2 *)pa, *b = (const PhiU2 *)pb;
    return a->phi < b->phi ? -1 : (a->phi > b->phi ? 1 : 0);
}

static void solve_parameterization(Arena_T arena, SliceSet *S, PairSet *P,
                                   double r_ref, const RibbonOpts *o,
                                   RibbonResult *out)
{
    size_t n = S->n_smp;
    size_t nc = S->n_chn;
    Arena_Mark mark = Arena_save(arena);

    /* Init the parameter globally by phi. When the winding is radius-consistent
     * (spiral fit usable) TRUE arc length is the spiral integral
     *   u(phi) = a*phi + b*phi^2/(4pi),   r(phi) = a + b*phi/(2pi),
     * a globally monotone map that places EVERY solve component -- detached
     * fragments, delaminated flaps, outer shells -- on one arc-length chart by
     * its phi, with NO per-component registration (the old binned phi->u map
     * cliff-extrapolated fragments whose phi fell outside the main component).
     * Otherwise fall back to the cylinder approximation phi*r_ref and register
     * detached components afterwards. */
    double sp_a = out->spiral_a, sp_b = out->spiral_b;
    /* Use the global spiral map whenever the fit ran with a valid |pitch|
     * slope (sp_b may be negative -- phi can wind either way). Do NOT gate on a
     * high r2: the pitch+sign+median-intercept line is robust to the outlier
     * fragments that drag an ordinary fit, and the phi*r_ref + registration
     * fallback places detached components WORSE (cliff-extrapolated). r2 == -1
     * marks "fit did not run" (too few samples) -> keep the simple init. */
    /* A PINNED line is trusted unconditionally: its r2 on a folded core cube
     * can be arbitrarily poor, but falling back to phi*r_ref + registration
     * would place that cube in a frame no other cube shares. */
    int spiral_ok = (o->spiral_b != 0.0)
                 || (out->spiral_r2 > -0.99 && fabs(sp_b) > RIB_PRIOR_MIN_B);
    double *u = (double *)ARENA_ALLOC(arena, (long)(n * sizeof(double)));
    for (size_t i = 0; i < n; i++) {
        double phi = S->smp[i].phi;
        u[i] = spiral_ok ? (sp_a * phi + sp_b * phi * phi / (4.0 * M_PI))
                         : phi * r_ref;
    }

    /* Drop CROSS-SLICE pairs whose two ends landed on winding-inconsistent phi
     * (a spurious cross-wrap match). Continuation pairs (the last n_cont, laid
     * out after the n_cross cross-slice pairs) are already radius-gated to the
     * SAME wrap (RIB_CONT_DR_MAX) so they must never be dropped: they stitch
     * the fragments of one wrap (e.g. across a hole/band), and dropping them
     * splits a wrap into detached components that then mis-register. Keeping
     * them lets the solve bridge the gap by the pair's geometric offset even
     * when the two fragments' per-chain phi disagree by a turn. */
    /* A reference field also makes continuation-pair validation possible.
     * Without it, continuations must survive a per-chain whole-turn error so
     * only cross-slice pairs are gated. With it, absolute phi is trustworthy
     * and a continuation that differs by >pi/2 is simply a contacting ply. */
    size_t n_phi_gate = o->reference_phi != NULL ? P->n_pairs : P->n_cross;
    for (size_t p = 0; p < n_phi_gate; p++) {
        Pair *pr = &P->pairs[p];
        if (fabs(S->smp[pr->a].phi - S->smp[pr->b].phi) > RIB_CONT_DPHI_MAX)
            pr->like = -1.0;   /* permanently off */
    }

    /* connected components of the solve graph: chains as nodes, active pairs
     * as edges. Each needs its own gauge pin (see SolveCtx). */
    UnionFind uf = UF_new(arena, (int32_t)(nc ? nc : 1));
    for (size_t p = 0; p < P->n_pairs; p++) {
        const Pair *pr = &P->pairs[p];
        if (pr->like < 0.0) continue;
        int32_t ca = S->smp[pr->a].chain, cb = S->smp[pr->b].chain;
        if (ca >= 0 && cb >= 0 && ca != cb)
            uf_union(&uf, ca, cb);
    }
    /* longest chain per component root -> pin its first sample at its init */
    int32_t *rootbest = (int32_t *)ARENA_ALLOC(arena, (long)((nc ? nc : 1) * sizeof(int32_t)));
    for (size_t c = 0; c < nc; c++) rootbest[c] = -1;
    int32_t bigchain = 0, bestlen = 0;
    for (size_t c = 0; c < nc; c++) {
        int32_t r = uf_find(&uf, (int32_t)c);
        if (rootbest[r] < 0 || S->chn[c].count > S->chn[rootbest[r]].count)
            rootbest[r] = (int32_t)c;
        if (S->chn[c].count > bestlen) { bestlen = S->chn[c].count; bigchain = (int32_t)c; }
    }
    int32_t *pin_idx = (int32_t *)ARENA_ALLOC(arena, (long)((nc ? nc : 1) * sizeof(int32_t)));
    double  *pin_val = (double  *)ARENA_ALLOC(arena, (long)((nc ? nc : 1) * sizeof(double)));
    size_t n_pins = 0;
    for (size_t c = 0; c < nc; c++) {
        if (uf_find(&uf, (int32_t)c) != (int32_t)c) continue;   /* roots only */
        int32_t pc = rootbest[c];
        if (pc < 0) continue;
        pin_idx[n_pins] = S->chn[pc].first;
        pin_val[n_pins] = u[S->chn[pc].first];
        n_pins++;
    }
    out->n_qp_comps = (int)n_pins;
    SolveCtx cx = { S, P, pin_idx, pin_val, n_pins, 1.0 };

    int32_t maxchain = 1;
    for (size_t c = 0; c < S->n_chn; c++)
        if (S->chn[c].count > maxchain) maxchain = S->chn[c].count;
    double *scr = (double *)ARENA_ALLOC(arena, (long)((size_t)maxchain * 3 * sizeof(double) + 64));
    double *lb  = (double *)ARENA_ALLOC(arena, (long)((size_t)maxchain * sizeof(double)));
    double *resbuf = (double *)ARENA_ALLOC(arena, (long)((P->n_pairs + 1) * sizeof(double)));

    int total_rounds = (o->relax_iters > 0 ? o->relax_iters : 1)
                     + (o->final_iters > 0 ? o->final_iters : 0);
    size_t mono_total = 0;
    for (int round = 0; round < total_rounds; round++) {
        cg_solve(arena, &cx, u);

        for (size_t c = 0; c < S->n_chn; c++) {
            int32_t f = S->chn[c].first, cn = S->chn[c].count;
            if (cn < 2) continue;
            lb[0] = 0.0;
            for (int32_t i = 1; i < cn; i++)
                lb[i] = 0.5 * (S->smp[f+i].s - S->smp[f+i-1].s);
            mono_total += pava_chain(&u[f], lb, cn, scr);
        }

        /* IRLS reweight scale from the RESIDUAL distribution (median |res|),
         * NOT the global u-span. With many jointly-solved components that span
         * is huge, so a span-based sigma never down-weights the bad cross-slice
         * matches that locally distort ("streak") an individual chart -- which
         * is why a single small piece parameterized cleanly while the whole
         * sheet streaked. The median tracks the good-pair noise (~vox) no matter
         * how many components share the solve, so outliers are rejected the same
         * whether we solve one piece or the whole scroll. */
        size_t npr = 0;
        for (size_t p = 0; p < P->n_pairs; p++) {
            const Pair *pr = &P->pairs[p];
            if (pr->like < 0.0) continue;
            resbuf[npr++] = fabs(u[pr->a] - u[pr->b] - pr->d);
        }
        double med = 0.0;
        if (npr > 0) { qsort(resbuf, npr, sizeof(double), cmp_dbl); med = resbuf[npr / 2]; }
        double sigma = (round < o->relax_iters ? 3.0 : 1.5) * med;
        if (sigma < 2.0) sigma = 2.0;   /* floor (vox): don't over-reject clean data */
        double inv2s2 = 1.0 / (2.0 * sigma * sigma);
        for (size_t p = 0; p < P->n_pairs; p++) {
            Pair *pr = &P->pairs[p];
            if (pr->like < 0.0) continue;   /* permanently dropped */
            double res = u[pr->a] - u[pr->b] - pr->d;
            pr->like = exp(-res * res * inv2s2);
        }
    }

    /* --- REGISTRATION: place every detached solve component onto the main
     * chart. The pins fix each component's gauge at its cylinder-approx init
     * (phi*r_ref), but the main chart's solved u is TRUE arc length, which
     * drifts from the init by ~15% across the block -- so init-anchored
     * detached charts land mid-band and overlap. Fix: build the main
     * component's monotone phi->u map and rigid-shift every other component
     * by median(f(phi) - u). Preserves each piece's internal arc length;
     * places outer-shell arcs BEYOND the main band (their phi is higher) and
     * genuinely delaminated layers exactly ON the main chart (same phi). */
    out->reg_max_shift = 0.0;
    if (n_pins > 1 && !spiral_ok) {
        Arena_Mark regm = Arena_save(arena);
        int32_t mainroot = uf_find(&uf, bigchain);
        /* per-sample solve-component root (resolved once) */
        int32_t *sroot = (int32_t *)ARENA_ALLOC(arena, (long)((n + 1) * sizeof(int32_t)));
        for (size_t i = 0; i < n; i++) {
            int32_t c = S->smp[i].chain;
            sroot[i] = c >= 0 ? uf_find(&uf, c) : -1;
        }
        /* gather main-component (phi, u) sorted by phi */
        size_t nmain = 0;
        for (size_t i = 0; i < n; i++) if (sroot[i] == mainroot) nmain++;
        PhiU2 *pu = (PhiU2 *)ARENA_ALLOC(arena, (long)((nmain + 1) * sizeof(PhiU2)));
        size_t q = 0;
        for (size_t i = 0; i < n; i++) {
            if (sroot[i] != mainroot) continue;
            pu[q].phi = S->smp[i].phi;
            pu[q].u   = u[i];
            q++;
        }
        qsort(pu, nmain, sizeof(PhiU2), cmp_phiu2);
        /* bin by phi (1/8 turn), median u per bin */
        PhiUMap M; memset(&M, 0, sizeof M);
        M.end_slope = r_ref;
        if (nmain >= 2) {
            double bw = M_PI / 4.0;
            int nb_max = (int)((pu[nmain-1].phi - pu[0].phi) / bw) + 2;
            if (nb_max < 1) nb_max = 1;
            M.bphi = (double *)ARENA_ALLOC(arena, (long)((size_t)nb_max * sizeof(double)));
            M.bu   = (double *)ARENA_ALLOC(arena, (long)((size_t)nb_max * sizeof(double)));
            double *tmpu = (double *)ARENA_ALLOC(arena, (long)((nmain + 1) * sizeof(double)));
            size_t i0 = 0;
            while (i0 < nmain) {
                double b0 = pu[0].phi + floor((pu[i0].phi - pu[0].phi) / bw) * bw;
                size_t i1 = i0;
                double psum = 0.0;
                size_t m2 = 0;
                while (i1 < nmain && pu[i1].phi < b0 + bw) {
                    tmpu[m2++] = pu[i1].u;
                    psum += pu[i1].phi;
                    i1++;
                }
                if (m2 > 0 && M.nb < nb_max) {
                    qsort(tmpu, m2, sizeof(double), cmp_dbl);
                    M.bphi[M.nb] = psum / (double)m2;
                    M.bu[M.nb]   = tmpu[m2 / 2];
                    M.nb++;
                }
                i0 = i1;
            }
        }
        if (M.nb >= 1) {
            /* per detached component: rigid shift to the map */
            double *diffs = (double *)ARENA_ALLOC(arena, (long)((n + 1) * sizeof(double)));
            for (size_t c = 0; c < nc; c++) {
                if (uf_find(&uf, (int32_t)c) != (int32_t)c) continue;
                if ((int32_t)c == mainroot) continue;
                size_t nd = 0;
                double u_lo = 1e300, u_hi = -1e300, ph_lo = 1e300, ph_hi = -1e300;
                for (size_t i = 0; i < n; i++) {
                    if (sroot[i] != (int32_t)c) continue;
                    diffs[nd++] = phiu_lookup(&M, S->smp[i].phi) - u[i];
                    if (u[i] < u_lo) u_lo = u[i];
                    if (u[i] > u_hi) u_hi = u[i];
                    if (S->smp[i].phi < ph_lo) ph_lo = S->smp[i].phi;
                    if (S->smp[i].phi > ph_hi) ph_hi = S->smp[i].phi;
                }
                if (nd == 0) continue;
                qsort(diffs, nd, sizeof(double), cmp_dbl);
                double shift = diffs[nd / 2];
#ifdef RIB_DEBUG_REG
                fprintf(stderr, "[reg] comp root=%d n=%zu phi=[%.2f,%.2f] "
                        "u_presolve=[%.1f,%.1f] f(mid)=%.1f shift=%.1f\n",
                        (int)c, nd, ph_lo, ph_hi, u_lo, u_hi,
                        phiu_lookup(&M, 0.5*(ph_lo+ph_hi)), shift);
#endif
                for (size_t i = 0; i < n; i++)
                    if (sroot[i] == (int32_t)c) u[i] += shift;
                if (fabs(shift) > out->reg_max_shift)
                    out->reg_max_shift = fabs(shift);
            }
        }
        Arena_restore(arena, regm);
    }

    for (size_t i = 0; i < n; i++) S->smp[i].u = u[i];
    out->mono_repairs = mono_total;

    /* --- canonical chart orientation: u increases OUTWARD (with radius).
     * The winding sense of the chart is a gauge (chain-walk dependent); a
     * mirrored chart is internally consistent but puts the outer wraps at
     * LOW u. Canonicalize so detached outer shells always land to the RIGHT
     * and runs are comparable. Mirror = negate u + reverse each chain so
     * per-chain u stays ascending (PAVA's invariant maps onto itself).
     * pin_orient (scroll_whole) SKIPS this: the flip is per-mesh, so cubes of
     * one scroll would mirror independently and never share a frame -- with a
     * pinned spiral the u direction is already deterministic grid-wide. */
    if (!o->pin_orient) {
        double ru = 0.0, rr = 0.0, cuv = 0.0;
        size_t nn = 0;
        for (size_t i = 0; i < n; i++) {
            if (S->smp[i].chain < 0) continue;
            ru += S->smp[i].u; rr += S->smp[i].r; nn++;
        }
        if (nn > 1) {
            ru /= (double)nn; rr /= (double)nn;
            for (size_t i = 0; i < n; i++) {
                if (S->smp[i].chain < 0) continue;
                cuv += (S->smp[i].u - ru) * (S->smp[i].r - rr);
            }
        }
        if (cuv < 0.0) {
            for (size_t i = 0; i < n; i++) S->smp[i].u = -S->smp[i].u;
            for (size_t c = 0; c < nc; c++) {
                int32_t f = S->chn[c].first, cn = S->chn[c].count;
                if (cn < 2) continue;
                for (int32_t i = 0; i < cn / 2; i++) {
                    Sample tmp = S->smp[f + i];
                    S->smp[f + i] = S->smp[f + cn - 1 - i];
                    S->smp[f + cn - 1 - i] = tmp;
                }
                S->smp[f].s = 0.0;
                for (int32_t i = 1; i < cn; i++) {
                    double dz = S->smp[f+i].p[0] - S->smp[f+i-1].p[0];
                    double dy = S->smp[f+i].p[1] - S->smp[f+i-1].p[1];
                    double dx = S->smp[f+i].p[2] - S->smp[f+i-1].p[2];
                    S->smp[f+i].s = S->smp[f+i-1].s + sqrt(dz*dz + dy*dy + dx*dx);
                }
                for (int32_t i = 0; i < cn; i++) {
                    S->smp[f+i].tau[0] = -S->smp[f+i].tau[0];
                    S->smp[f+i].tau[1] = -S->smp[f+i].tau[1];
                    S->smp[f+i].tau[2] = -S->smp[f+i].tau[2];
                }
            }
        }
    }

    double emean = 0.0, emax = 0.0;
    long nedge = 0;
    memset(out->duds_hist, 0, sizeof(out->duds_hist));
    for (size_t c = 0; c < S->n_chn; c++) {
        int32_t f = S->chn[c].first, cn = S->chn[c].count;
        for (int32_t i = 1; i < cn; i++) {
            double l = S->smp[f+i].s - S->smp[f+i-1].s;
            if (l < 1e-9) continue;
            double e = fabs((S->smp[f+i].u - S->smp[f+i-1].u) / l - 1.0);
            emean += e;
            if (e > emax) emax = e;
            nedge++;
            double pc = e * 100.0;
            int bin = pc < 1 ? 0 : pc < 2 ? 1 : pc < 5 ? 2 : pc < 10 ? 3 : pc < 20 ? 4 : 5;
            out->duds_hist[bin]++;
        }
    }
    out->duds_err_mean = nedge ? emean / (double)nedge : 0.0;
    out->duds_err_max  = emax;

    Arena_restore(arena, mark);
}

/* ============================================================================
 * Stage D -- transfer (u, v) to the original mesh vertices via the 3D-nearest
 * slice sample (in-plane search < RIB_XFER_R, inside the no-merger clearance;
 * tangential correction tau . (p_v - p_s)).
 * ==========================================================================*/

static double transfer_uv(Arena_T arena, const float *verts, size_t nv,
                          const int32_t *faces, size_t nf,
                          const double *t,
                          const SliceSet *S, double slice_h,
                          const RibbonOpts *o, RibbonResult *out)
{
    double h = slice_h;
    double tmin = S->tmin;
    int np = S->nplanes;
    const Sample *smp = S->smp;

    /* out->uv arrives PREFILLED by the caller with each vertex's in-plane
     * (c1,c2) coords -- the bucketing key for the per-slice grids below. The
     * final (u,v) overwrites it at the end. */
    out->uv_ok = (uint8_t *)ARENA_CALLOC(arena, (long)nv, 1);

    Arena_Mark mark = Arena_save(arena);

    /* per-vertex bracketing slices + blend weight */
    int32_t *vk0 = (int32_t *)ARENA_ALLOC(arena, (long)(nv * sizeof(int32_t)));
    float   *vw  = (float *)ARENA_ALLOC(arena, (long)(nv * sizeof(float)));
    for (size_t i = 0; i < nv; i++) {
        double fb = (t[i] - tmin) / h - 0.5;
        int k0 = (int)floor(fb);
        double w = fb - (double)k0;
        if (k0 < 0)      { k0 = 0;      w = 0.0; }
        if (k0 > np - 1) { k0 = np - 1; w = 0.0; }
        vk0[i] = k0;
        vw[i] = (float)w;
    }

    /* per-slice sample grids: bucket samples by slice, then a uniform grid on
     * (c1,c2) per slice, queried by vertices bracketing that slice */
    int32_t *soff = NULL;
    int32_t *sids = sort_by_slice(arena, S, &soff);

    double c1min = 1e300, c1max = -1e300, c2min = 1e300, c2max = -1e300;
    for (size_t i = 0; i < S->n_smp; i++) {
        if (smp[i].chain < 0) continue;
        if (smp[i].c1 < c1min) c1min = smp[i].c1;
        if (smp[i].c1 > c1max) c1max = smp[i].c1;
        if (smp[i].c2 < c2min) c2min = smp[i].c2;
        if (smp[i].c2 > c2max) c2max = smp[i].c2;
    }
    double cell = RIB_XFER_R;
    int g1 = (int)((c1max - c1min) / cell) + 1;
    int g2 = (int)((c2max - c2min) / cell) + 1;
    if (g1 < 1) g1 = 1;
    if (g2 < 1) g2 = 1;
    size_t ncell = (size_t)g1 * (size_t)g2;

    double *u0 = (double *)ARENA_ALLOC(arena, (long)(nv * sizeof(double)));
    double *u1 = (double *)ARENA_ALLOC(arena, (long)(nv * sizeof(double)));
    double *ph0 = (double *)ARENA_ALLOC(arena, nv * sizeof(double));
    double *ph1 = (double *)ARENA_ALLOC(arena, nv * sizeof(double));
    int32_t *gr0 = (int32_t *)ARENA_ALLOC(arena, nv * sizeof(int32_t));
    int32_t *gr1 = (int32_t *)ARENA_ALLOC(arena, nv * sizeof(int32_t));
    uint8_t *ok0 = (uint8_t *)ARENA_CALLOC(arena, (long)nv, 1);
    uint8_t *ok1 = (uint8_t *)ARENA_CALLOC(arena, (long)nv, 1);
    uint8_t *mapped = (uint8_t *)ARENA_CALLOC(arena, (long)nv, 1);

    /* vertices bucketed by k0 (counting sort) so each slice's grid is built
     * once and queried by the two vertex populations that bracket it */
    int32_t *voff = (int32_t *)ARENA_CALLOC(arena, (long)(np + 2), (long)sizeof(int32_t));
    for (size_t i = 0; i < nv; i++) voff[vk0[i] + 1]++;
    for (int k = 0; k <= np; k++) voff[k+1] = (int32_t)(voff[k+1] + voff[k]);
    int32_t *vids = (int32_t *)ARENA_ALLOC(arena, (long)(nv * sizeof(int32_t)));
    int32_t *vcur = (int32_t *)ARENA_ALLOC(arena, (long)((size_t)(np + 1) * sizeof(int32_t)));
    memcpy(vcur, voff, (size_t)(np + 1) * sizeof(int32_t));
    for (size_t i = 0; i < nv; i++) vids[vcur[vk0[i]]++] = (int32_t)i;

    int32_t *goff = (int32_t *)ARENA_ALLOC(arena, (long)((ncell + 1) * sizeof(int32_t)));
    int32_t *gids = (int32_t *)ARENA_ALLOC(arena, (long)((S->n_smp + 1) * sizeof(int32_t)));
    int32_t *gcur = (int32_t *)ARENA_ALLOC(arena, (long)((ncell + 1) * sizeof(int32_t)));

    for (int k = 0; k < np; k++) {
        int32_t s0 = soff[k], s1 = soff[k+1];
        if (s1 <= s0) continue;
        /* grid slice k */
        memset(goff, 0, (ncell + 1) * sizeof(int32_t));
        for (int32_t q = s0; q < s1; q++) {
            const Sample *sp = &smp[sids[q]];
            int i1 = (int)((sp->c1 - c1min) / cell), i2 = (int)((sp->c2 - c2min) / cell);
            if (i1 < 0) i1 = 0;
            if (i1 > g1-1) i1 = g1-1;
            if (i2 < 0) i2 = 0;
            if (i2 > g2-1) i2 = g2-1;
            goff[(size_t)i1 * (size_t)g2 + (size_t)i2 + 1]++;
        }
        for (size_t cix = 0; cix < ncell; cix++)
            goff[cix+1] = (int32_t)(goff[cix+1] + goff[cix]);
        memcpy(gcur, goff, ncell * sizeof(int32_t));
        for (int32_t q = s0; q < s1; q++) {
            const Sample *sp = &smp[sids[q]];
            int i1 = (int)((sp->c1 - c1min) / cell), i2 = (int)((sp->c2 - c2min) / cell);
            if (i1 < 0) i1 = 0;
            if (i1 > g1-1) i1 = g1-1;
            if (i2 < 0) i2 = 0;
            if (i2 > g2-1) i2 = g2-1;
            gids[gcur[(size_t)i1 * (size_t)g2 + (size_t)i2]++] = sids[q];
        }
        /* query vertices with k0 == k (their lower slice) and k0 == k-1
         * (their upper slice) */
        for (int side = 0; side < 2; side++) {
            int kv = side == 0 ? k : k - 1;
            if (kv < 0 || kv > np - 1) continue;
            for (int32_t q = voff[kv]; q < voff[kv+1]; q++) {
                int32_t vi = vids[q];
                /* skip second lookup when k1 == k0 (top clamp) */
                if (side == 1 && vk0[vi] + 1 > np - 1) continue;
                double pc1 = (double)out->uv[(size_t)vi*2 + 0];
                double pc2 = (double)out->uv[(size_t)vi*2 + 1];
                int i1 = (int)((pc1 - c1min) / cell), i2 = (int)((pc2 - c2min) / cell);
                int32_t best = -1; double bestd2 = RIB_XFER_R * RIB_XFER_R;
                for (int d1 = -1; d1 <= 1; d1++) for (int d2i = -1; d2i <= 1; d2i++) {
                    int j1 = i1 + d1, j2 = i2 + d2i;
                    if (j1 < 0 || j2 < 0 || j1 >= g1 || j2 >= g2) continue;
                    size_t cix = (size_t)j1 * (size_t)g2 + (size_t)j2;
                    for (int32_t rr = goff[cix]; rr < goff[cix+1]; rr++) {
                        int32_t is = gids[rr];
                        double dd1 = smp[is].c1 - pc1, dd2 = smp[is].c2 - pc2;
                        double d2 = dd1*dd1 + dd2*dd2;
                        if (d2 < bestd2) { bestd2 = d2; best = is; }
                    }
                }
                if (best >= 0) {
                    const Sample *sp = &smp[best];
                    double dp[3] = { (double)verts[(size_t)vi*3+0] - sp->p[0],
                                     (double)verts[(size_t)vi*3+1] - sp->p[1],
                                     (double)verts[(size_t)vi*3+2] - sp->p[2] };
                    double uu = sp->u + v3dot(sp->tau, dp);
                    /* phi/group: nearest-sample values (no tangential
                     * correction -- registration medians round to whole
                     * turns, so the <= sample-spacing error is irrelevant) */
                    int32_t gg = sp->chain >= 0 ? S->chn[sp->chain].group : -1;
                    if (side == 0) {
                        u0[vi] = uu; ph0[vi] = sp->phi; gr0[vi] = gg; ok0[vi] = 1;
                    } else {
                        u1[vi] = uu; ph1[vi] = sp->phi; gr1[vi] = gg; ok1[vi] = 1;
                    }
                }
            }
        }
    }

    double *uvert = (double *)ARENA_ALLOC(arena, (long)(nv * sizeof(double)));
    double *phvert = (double *)ARENA_ALLOC(arena, nv * sizeof(double));
    int32_t *gvert = (int32_t *)ARENA_ALLOC(arena, nv * sizeof(int32_t));
    size_t fallback = 0;
    for (size_t i = 0; i < nv; i++) {
        double w = (double)vw[i];
        if (ok0[i] && ok1[i]) {
            uvert[i]  = (1.0 - w) * u0[i]  + w * u1[i];
            phvert[i] = (1.0 - w) * ph0[i] + w * ph1[i];
            gvert[i]  = w < 0.5 ? gr0[i] : gr1[i];   /* nearer slice's group */
        }
        else if (ok0[i]) { uvert[i] = u0[i]; phvert[i] = ph0[i]; gvert[i] = gr0[i]; }
        else if (ok1[i]) { uvert[i] = u1[i]; phvert[i] = ph1[i]; gvert[i] = gr1[i]; }
        else { uvert[i] = 0.0; phvert[i] = 0.0; gvert[i] = -1; fallback++; continue; }
        out->uv_ok[i] = 1;
        mapped[i] = 1;
    }

    /* neighbor fill for unmapped verts (median of mapped neighbors, robust to
     * a cross-bridge neighbor) */
    if (fallback > 0) {
        CSR_T adj = CSR_from_faces(arena, faces, nf, nv);
        const int32_t *aoff = CSR_offset(adj);
        const int32_t *atgt = CSR_target(adj);
        double nb[64], nbp[64];
        int32_t nbg[64];
        for (int round = 0; round < RIB_UVFILL_ROUNDS && fallback > 0; round++) {
            size_t fixed = 0;
            for (size_t i = 0; i < nv; i++) {
                if (mapped[i]) continue;
                int cnt = 0;
                for (int32_t e = aoff[i]; e < aoff[i+1] && cnt < 64; e++) {
                    int32_t j = atgt[e];
                    if (mapped[j] == 1) {
                        nb[cnt] = uvert[j]; nbp[cnt] = phvert[j];
                        nbg[cnt] = gvert[j]; cnt++;
                    }
                }
                if (cnt > 0) {
                    qsort(nb, (size_t)cnt, sizeof(double), cmp_dbl);
                    qsort(nbp, (size_t)cnt, sizeof(double), cmp_dbl);
                    qsort(nbg, (size_t)cnt, sizeof(int32_t), cmp_i32);
                    uvert[i]  = nb[cnt / 2];
                    phvert[i] = nbp[cnt / 2];
                    gvert[i]  = nbg[cnt / 2];
                    mapped[i] = 2;
                    fixed++;
                }
            }
            for (size_t i = 0; i < nv; i++)
                if (mapped[i] == 2) mapped[i] = 1;
            fallback -= fixed;
            if (fixed == 0) break;
        }
    }
    size_t nfilled = 0, nfb = 0;
    for (size_t i = 0; i < nv; i++) {
        if (!mapped[i]) nfb++;
        else if (!out->uv_ok[i]) nfilled++;
    }
    out->uv_filled = nfilled;
    out->uv_fallback = nfb;

    /* shift so the SAMPLE minimum is 0 (ribbon grid shares this origin).
     * emit_global keeps the absolute frame instead: u in the solve's frame
     * (one pinned spiral => comparable across cubes), v = raw axial t. */
    double U0 = 1e300;
    for (size_t i = 0; i < S->n_smp; i++)
        if (smp[i].chain >= 0 && smp[i].u < U0) U0 = smp[i].u;
    out->u_origin = (U0 < 1e299) ? U0 : 0.0;
    for (size_t i = 0; i < nv; i++) {
        if (o->emit_global) {
            out->uv[i*2 + 0] = (float)uvert[i];
            out->uv[i*2 + 1] = (float)t[i];
        } else {
            out->uv[i*2 + 0] = (float)(uvert[i] - U0);
            out->uv[i*2 + 1] = (float)(t[i] - tmin);
        }
    }
    if (out->phi != NULL)
        for (size_t i = 0; i < nv; i++) out->phi[i] = (float)phvert[i];
    if (out->group != NULL)
        for (size_t i = 0; i < nv; i++) out->group[i] = gvert[i];
    Arena_restore(arena, mark);
    return U0;
}

/* ============================================================================
 * Stage E -- ribbon fit on a regular (u, v) grid (unchanged).
 * ==========================================================================*/

static void thomas(double *a, double *b, double *c, double *d, int n)
{
    for (int i = 1; i < n; i++) {
        double m = a[i] / b[i-1];
        b[i] -= m * c[i-1];
        d[i] -= m * d[i-1];
    }
    d[n-1] /= b[n-1];
    for (int i = n - 2; i >= 0; i--)
        d[i] = (d[i] - c[i] * d[i+1]) / b[i];
}

static void fit_ribbon(Arena_T arena, const SliceSet *S, double U0,
                       const RibbonOpts *o, RibbonResult *out)
{
    double du = (double)o->grid_u;
    if (du < 1e-6) du = 1.0;
    double umax = -1e300;
    for (size_t i = 0; i < S->n_smp; i++) {
        if (S->smp[i].chain < 0) continue;
        double uu = S->smp[i].u - U0;
        if (uu > umax) umax = uu;
    }
    if (umax <= 0.0) return;
    size_t nu = (size_t)(umax / du) + 2;
    size_t nk = (size_t)S->nplanes;
    if (nu < 2 || nk < 1) return;
    if (nu > 500000) return;   /* runaway guard */

    float   *G  = (float *)ARENA_ALLOC(arena, (long)(nk * nu * 3 * sizeof(float)));
    uint8_t *GV = (uint8_t *)ARENA_CALLOC(arena, (long)(nk * nu), 1);
    for (size_t i = 0; i < nk * nu * 3; i++) G[i] = (float)NAN;

    Arena_Mark mark = Arena_save(arena);
    double *P   = (double *)ARENA_ALLOC(arena, (long)(nu * 3 * sizeof(double)));
    double *T   = (double *)ARENA_ALLOC(arena, (long)(nu * 3 * sizeof(double)));
    uint8_t *V  = (uint8_t *)ARENA_ALLOC(arena, (long)nu);
    double *ta_ = (double *)ARENA_ALLOC(arena, (long)(nu * sizeof(double)));
    double *tb_ = (double *)ARENA_ALLOC(arena, (long)(nu * sizeof(double)));
    double *tc_ = (double *)ARENA_ALLOC(arena, (long)(nu * sizeof(double)));
    double *td_ = (double *)ARENA_ALLOC(arena, (long)(nu * sizeof(double)));
    double *X   = (double *)ARENA_ALLOC(arena, (long)(nu * 3 * sizeof(double)));

    int np = S->nplanes;
    int32_t *choff = NULL;
    int32_t *chids = chains_by_slice(arena, S, &choff);

    for (int k = 0; k < np; k++) {
        memset(V, 0, nu);
        int rowvalid = 0;
        for (int32_t x = choff[k]; x < choff[k+1]; x++) {
            const Chain *ch = &S->chn[chids[x]];
            int32_t f = ch->first, cn = ch->count;
            if (cn < 2) continue;
            double ulo = S->smp[f].u - U0, uhi = S->smp[f+cn-1].u - U0;
            int32_t j0 = (int32_t)ceil(ulo / du), j1 = (int32_t)floor(uhi / du);
            if (j0 < 0) j0 = 0;
            if (j1 > (int32_t)nu - 1) j1 = (int32_t)nu - 1;
            int32_t seg = 1;
            for (int32_t j = j0; j <= j1; j++) {
                double uq = (double)j * du + U0;
                while (seg < cn - 1 && S->smp[f+seg].u < uq) seg++;
                double ua = S->smp[f+seg-1].u, ub = S->smp[f+seg].u;
                double w = ub - ua > 1e-12 ? (uq - ua) / (ub - ua) : 0.0;
                if (w < 0.0) w = 0.0;
                if (w > 1.0) w = 1.0;
                for (int d = 0; d < 3; d++) {
                    P[(size_t)j*3 + (size_t)d] = (1.0-w)*S->smp[f+seg-1].p[d] + w*S->smp[f+seg].p[d];
                    T[(size_t)j*3 + (size_t)d] = (1.0-w)*S->smp[f+seg-1].tau[d] + w*S->smp[f+seg].tau[d];
                }
                V[j] = 1;
                rowvalid++;
            }
        }
        if (rowvalid == 0) continue;

        /* Split the row into segments at (a) invalid runs > RIB_UFILL_MAX --
         * integrating tangents across a long hole extrudes giant free-space
         * arcs ("hula hoops") -- and (b) 3D jumps > RIB_WRAP_GATE between
         * consecutive valid columns or across a bridged gap -- a wrap boundary
         * (winding collapse put two wraps at nearby u). Both would otherwise
         * emit a triangle spanning wraps. */
        size_t jseek = 0;
        while (jseek < nu) {
            while (jseek < nu && !V[jseek]) jseek++;
            if (jseek >= nu) break;
            size_t jf = jseek, jl = jseek, scan = jseek + 1;
            for (;;) {
                int wrap_break = 0;
                /* extend over consecutive valid columns; stop at a wrap jump */
                while (scan < nu && V[scan]) {
                    double dz = P[jl*3]   - P[scan*3];
                    double dy = P[jl*3+1] - P[scan*3+1];
                    double dx = P[jl*3+2] - P[scan*3+2];
                    if (sqrt(dz*dz + dy*dy + dx*dx) > RIB_WRAP_GATE) { wrap_break = 1; break; }
                    jl = scan; scan++;
                }
                if (wrap_break) break;
                size_t run0 = scan;
                while (scan < nu && !V[scan]) scan++;
                if (scan >= nu || scan - run0 > RIB_UFILL_MAX) break;
                /* bridge the short gap only if it does not cross a wrap */
                {
                    double dz = P[jl*3]   - P[scan*3];
                    double dy = P[jl*3+1] - P[scan*3+1];
                    double dx = P[jl*3+2] - P[scan*3+2];
                    if (sqrt(dz*dz + dy*dy + dx*dx)
                        > RIB_WRAP_GATE * (double)(scan - jl)) break;
                }
                jl = scan; scan++;
            }
            jseek = jl + 1;
            int n = (int)(jl - jf + 1);
            if (n < 2) {
                for (int d = 0; d < 3; d++)
                    G[(((size_t)k)*nu + jf)*3 + (size_t)d] = (float)P[jf*3 + (size_t)d];
                GV[(size_t)k*nu + jf] = 1;
                continue;
            }

            for (int d = 0; d < 3; d++) {
                for (int j = 0; j < n; j++) {
                    size_t jj = jf + (size_t)j;
                    double wdat = V[jj] ? 1.0 : 0.0;
                    double lap = RIB_LAMBDA_C * ((j > 0 ? 1.0 : 0.0) + (j < n-1 ? 1.0 : 0.0));
                    ta_[j] = j > 0     ? -RIB_LAMBDA_C : 0.0;
                    tc_[j] = j < n - 1 ? -RIB_LAMBDA_C : 0.0;
                    tb_[j] = wdat + lap;
                    td_[j] = wdat * T[jj*3 + (size_t)d];
                }
                thomas(ta_, tb_, tc_, td_, n);
                for (int j = 0; j < n; j++) X[(size_t)j*3 + (size_t)d] = td_[j];
            }
            for (int j = 0; j < n; j++) {
                double *tj = &X[(size_t)j*3];
                if (v3norm(tj) < 1e-9) { tj[0] = 0; tj[1] = 0; tj[2] = 0; }
            }

            for (int d = 0; d < 3; d++) {
                for (int j = 0; j < n; j++) {
                    size_t jj = jf + (size_t)j;
                    double wdat = V[jj] ? 1.0 : 0.0;
                    double diag = wdat + RIB_LAMBDA_T * ((j > 0 ? 1.0 : 0.0) + (j < n-1 ? 1.0 : 0.0));
                    double rhs = wdat * P[jj*3 + (size_t)d];
                    if (j > 0) {
                        double g = du * 0.5 * (X[(size_t)(j-1)*3 + (size_t)d] + X[(size_t)j*3 + (size_t)d]);
                        rhs += RIB_LAMBDA_T * g;
                    }
                    if (j < n - 1) {
                        double g = du * 0.5 * (X[(size_t)j*3 + (size_t)d] + X[(size_t)(j+1)*3 + (size_t)d]);
                        rhs -= RIB_LAMBDA_T * g;
                    }
                    ta_[j] = j > 0     ? -RIB_LAMBDA_T : 0.0;
                    tc_[j] = j < n - 1 ? -RIB_LAMBDA_T : 0.0;
                    tb_[j] = diag;
                    td_[j] = rhs;
                }
                thomas(ta_, tb_, tc_, td_, n);
                for (int j = 0; j < n; j++)
                    G[(((size_t)k)*nu + jf + (size_t)j)*3 + (size_t)d] = (float)td_[j];
            }
            for (int j = 0; j < n; j++)
                GV[(size_t)k*nu + jf + (size_t)j] = V[jf + (size_t)j] ? 1 : 0;
        }
    }

    for (size_t j = 0; j < nu; j++) {
        for (size_t k = 0; k < nk; k++) {
            if (!isnan((double)G[(k*nu + j)*3])) continue;
            long ka = -1, kb = -1;
            for (long q = (long)k - 1; q >= 0 && q >= (long)k - RIB_VFILL_MAX; q--)
                if (!isnan((double)G[((size_t)q*nu + j)*3])) { ka = q; break; }
            for (long q = (long)k + 1; q < (long)nk && q <= (long)k + RIB_VFILL_MAX; q++)
                if (!isnan((double)G[((size_t)q*nu + j)*3])) { kb = q; break; }
            if (ka < 0 || kb < 0) continue;
            /* only bridge if the two anchors are on the same wrap: a per-row
             * step above the gate means the hole spans wraps -- leave it a hole */
            {
                double dz = (double)G[((size_t)ka*nu+j)*3]   - (double)G[((size_t)kb*nu+j)*3];
                double dy = (double)G[((size_t)ka*nu+j)*3+1] - (double)G[((size_t)kb*nu+j)*3+1];
                double dx = (double)G[((size_t)ka*nu+j)*3+2] - (double)G[((size_t)kb*nu+j)*3+2];
                if (sqrt(dz*dz + dy*dy + dx*dx) > RIB_WRAP_GATE * (double)(kb - ka)) continue;
            }
            double w = (double)((long)k - ka) / (double)(kb - ka);
            for (int d = 0; d < 3; d++)
                G[(k*nu + j)*3 + (size_t)d] =
                    (float)((1.0-w) * (double)G[((size_t)ka*nu + j)*3 + (size_t)d]
                          +       w  * (double)G[((size_t)kb*nu + j)*3 + (size_t)d]);
        }
    }
    for (int sweep = 0; sweep < 2; sweep++) {
        for (size_t k = 1; k + 1 < nk; k++) {
            for (size_t j = 0; j < nu; j++) {
                size_t c0 = (k*nu + j)*3, cu = ((k-1)*nu + j)*3, cd = ((k+1)*nu + j)*3;
                if (isnan((double)G[c0]) || isnan((double)G[cu]) || isnan((double)G[cd]))
                    continue;
                for (int d = 0; d < 3; d++) {
                    double m = 0.5 * ((double)G[cu + (size_t)d] + (double)G[cd + (size_t)d]);
                    G[c0 + (size_t)d] = (float)(((double)G[c0 + (size_t)d] + RIB_VSMOOTH_MU * m)
                                                / (1.0 + RIB_VSMOOTH_MU));
                }
            }
        }
    }

    Arena_restore(arena, mark);
    out->grid_pos   = G;
    out->grid_valid = GV;
    out->nu = nu;
    out->nk = nk;
    out->grid_du = (float)du;
    out->grid_dv = o->slice_h;
}

/* ============================================================================
 * Orchestration.
 * ==========================================================================*/

int Ribbon_run(Arena_T arena,
               const float *verts, size_t nv,
               const int32_t *faces, size_t nf,
               const RibbonOpts *opts, RibbonResult *out)
{
    assert(arena && out);
    memset(out, 0, sizeof(*out));
    out->spiral_r2 = -1.0;
    if (nv < 3 || nf < 1 || verts == NULL || faces == NULL)
        return -1;

    RibbonOpts def;
    if (opts == NULL) { RibbonOpts_default(&def); opts = &def; }

    /* axis frame */
    double ad[3] = { (double)opts->axis_dir[0], (double)opts->axis_dir[1],
                     (double)opts->axis_dir[2] };
    double an = sqrt(v3dot(ad, ad));
    if (an < 1e-12) { ad[0] = 1.0; ad[1] = 0.0; ad[2] = 0.0; an = 1.0; }
    ad[0] /= an; ad[1] /= an; ad[2] /= an;
    float axf[3] = { (float)ad[0], (float)ad[1], (float)ad[2] };
    float e1f[3], e2f[3];
    PCA_orthonormal_basis(axf, e1f, e2f);
    double e1[3] = { (double)e1f[0], (double)e1f[1], (double)e1f[2] };
    double e2[3] = { (double)e2f[0], (double)e2f[1], (double)e2f[2] };
    /* canonical handedness: cross(e1,e2) . axis > 0 in stored (z,y,x)
     * components, so theta's sense is deterministic across bases */
    {
        double cr[3] = { e1[1]*e2[2] - e1[2]*e2[1],
                         e1[2]*e2[0] - e1[0]*e2[2],
                         e1[0]*e2[1] - e1[1]*e2[0] };
        if (v3dot(cr, ad) < 0.0) {
            double tswap[3] = { e1[0], e1[1], e1[2] };
            e1[0] = e2[0]; e1[1] = e2[1]; e1[2] = e2[2];
            e2[0] = tswap[0]; e2[1] = tswap[1]; e2[2] = tswap[2];
        }
    }

    /* per-vertex axial t */
    double *t = (double *)ARENA_ALLOC(arena, (long)(nv * sizeof(double)));
    double tmin = 1e300, tmax = -1e300;
    for (size_t i = 0; i < nv; i++) {
        double d0 = (double)verts[i*3+0] - (double)opts->axis_point[0];
        double d1 = (double)verts[i*3+1] - (double)opts->axis_point[1];
        double d2 = (double)verts[i*3+2] - (double)opts->axis_point[2];
        t[i] = d0*ad[0] + d1*ad[1] + d2*ad[2];
        if (t[i] < tmin) tmin = t[i];
        if (t[i] > tmax) tmax = t[i];
    }

    /* A: slice (+ bridge cut) */
    SliceSet S;
    if (slice_mesh(arena, verts, nv, faces, nf, t, tmin, tmax,
                   e1, e2, opts->axis_point, opts, &S) != 0) {
        memset(out, 0, sizeof(*out));
        return -1;
    }
    out->n_slices       = S.n_slices_hit;
    out->n_chains       = (int)S.n_chn;
    out->n_closed       = S.n_closed;
    out->n_multi_slices = S.n_multi;
    out->n_samples      = S.n_smp;
    out->bridge_cuts    = S.bridge_cuts;
    out->v_span         = tmax - tmin;

    /* r_ref = median in-plane sample radius */
    double r_ref = 1.0;
    {
        Arena_Mark m2 = Arena_save(arena);
        double *rad = (double *)ARENA_ALLOC(arena, (long)((S.n_smp + 1) * sizeof(double)));
        size_t q = 0;
        for (size_t i = 0; i < S.n_smp; i++)
            if (S.smp[i].chain >= 0) rad[q++] = S.smp[i].r;
        if (q > 0) {
            qsort(rad, q, sizeof(double), cmp_dbl);
            r_ref = rad[q / 2];
        }
        if (r_ref < 1e-6) r_ref = 1.0;
        Arena_restore(arena, m2);
    }

    /* B: pair candidates */
    PairSet P;
    build_pairs(arena, &S, opts, &P);
    out->n_pairs      = P.n_cross;
    out->n_cont_pairs = P.n_cont;
    out->match_cover  = P.candidates_a ? (double)P.matched_a / (double)P.candidates_a : 0.0;

    /* W: integer winding index per chain -> phi = th + 2pi W */
    assign_winding(arena, &S, &P, opts, out);
    {
        double pmin = 1e300, pmax = -1e300;
        for (size_t i = 0; i < S.n_smp; i++) {
            if (S.smp[i].chain < 0) continue;
            if (S.smp[i].phi < pmin) pmin = S.smp[i].phi;
            if (S.smp[i].phi > pmax) pmax = S.smp[i].phi;
        }
        out->phi_span_turns = (pmax - pmin) / (2.0 * M_PI);
    }

    /* C: solve */
    solve_parameterization(arena, &S, &P, r_ref, opts, out);
    {
        double lo = 1e300, hi = -1e300;
        for (size_t i = 0; i < S.n_smp; i++) {
            if (S.smp[i].chain < 0) continue;
            if (S.smp[i].u < lo) lo = S.smp[i].u;
            if (S.smp[i].u > hi) hi = S.smp[i].u;
        }
        out->u_span = hi - lo;
    }

    /* D: transfer to mesh verts. Pre-fill out->uv with vertex in-plane coords
     * (the transfer's bucketing key); replaced by (u,v) inside. */
    out->uv = (float *)ARENA_ALLOC(arena, (long)(nv * 2 * sizeof(float)));
    if (opts->emit_global) {
        out->phi = (float *)ARENA_CALLOC(arena, nv, sizeof(float));
        out->group = (int32_t *)ARENA_ALLOC(arena, nv * sizeof(int32_t));
        for (size_t i = 0; i < nv; i++) out->group[i] = -1;
    }
    for (size_t i = 0; i < nv; i++) {
        double d0 = (double)verts[i*3+0] - (double)opts->axis_point[0];
        double d1 = (double)verts[i*3+1] - (double)opts->axis_point[1];
        double d2 = (double)verts[i*3+2] - (double)opts->axis_point[2];
        out->uv[i*2+0] = (float)(d0*e1[0] + d1*e1[1] + d2*e1[2]);
        out->uv[i*2+1] = (float)(d0*e2[0] + d1*e2[1] + d2*e2[2]);
    }
    {
        double U0 = transfer_uv(arena, verts, nv, faces, nf, t, &S,
                                (double)opts->slice_h, opts, out);
        /* E: ribbon */
        if (opts->fit_ribbon)
            fit_ribbon(arena, &S, U0, opts, out);
        else
            out->grid_dv = opts->slice_h;
    }
    return 0;
}

/* ============================================================================
 * Coarse -> original UV transfer (simplify-first workflow).
 * ==========================================================================*/

/* closest point on triangle abc to p (Ericson, Real-Time Collision Detection) */
static void closest_on_tri(const double p[3], const double a[3],
                           const double b[3], const double c[3],
                           double out[3], double *out_u, double *out_v)
{
    double ab[3] = { b[0]-a[0], b[1]-a[1], b[2]-a[2] };
    double ac[3] = { c[0]-a[0], c[1]-a[1], c[2]-a[2] };
    double ap[3] = { p[0]-a[0], p[1]-a[1], p[2]-a[2] };
    double d1 = v3dot(ab, ap), d2 = v3dot(ac, ap);
    if (d1 <= 0.0 && d2 <= 0.0) { memcpy(out, a, 3*sizeof(double)); *out_u = 0; *out_v = 0; return; }
    double bp[3] = { p[0]-b[0], p[1]-b[1], p[2]-b[2] };
    double d3 = v3dot(ab, bp), d4 = v3dot(ac, bp);
    if (d3 >= 0.0 && d4 <= d3) { memcpy(out, b, 3*sizeof(double)); *out_u = 1; *out_v = 0; return; }
    double vc = d1*d4 - d3*d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
        double w = d1 / (d1 - d3);
        for (int i = 0; i < 3; i++) out[i] = a[i] + w * ab[i];
        *out_u = w; *out_v = 0; return;
    }
    double cp[3] = { p[0]-c[0], p[1]-c[1], p[2]-c[2] };
    double d5 = v3dot(ab, cp), d6 = v3dot(ac, cp);
    if (d6 >= 0.0 && d5 <= d6) { memcpy(out, c, 3*sizeof(double)); *out_u = 0; *out_v = 1; return; }
    double vb = d5*d2 - d1*d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
        double w = d2 / (d2 - d6);
        for (int i = 0; i < 3; i++) out[i] = a[i] + w * ac[i];
        *out_u = 0; *out_v = w; return;
    }
    double va = d3*d6 - d5*d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
        double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        for (int i = 0; i < 3; i++) out[i] = b[i] + w * (c[i] - b[i]);
        *out_u = 1.0 - w; *out_v = w; return;
    }
    double denom = 1.0 / (va + vb + vc);
    double v_ = vb * denom, w_ = vc * denom;
    for (int i = 0; i < 3; i++) out[i] = a[i] + ab[i]*v_ + ac[i]*w_;
    *out_u = v_; *out_v = w_;
}

int Ribbon_flag_bad_faces(const float *verts, size_t nv,
                          const int32_t *faces, size_t nf,
                          const float *uv, double ratio, double floor_vox,
                          double len_min, uint8_t *out_bad, size_t *out_n)
{
    (void)nv;
    if (out_n) *out_n = 0;
    if (verts == NULL || faces == NULL || uv == NULL || out_bad == NULL) return -1;
    size_t nbad = 0;
    for (size_t f = 0; f < nf; f++) {
        int bad = 0;
        for (int e = 0; e < 3 && !bad; e++) {
            int32_t a = faces[f*3 + e], b = faces[f*3 + (e+1)%3];
            double du = fabs((double)uv[(size_t)a*2] - (double)uv[(size_t)b*2]);
            if (du <= floor_vox) continue;
            double dz = (double)verts[(size_t)a*3]   - (double)verts[(size_t)b*3];
            double dy = (double)verts[(size_t)a*3+1] - (double)verts[(size_t)b*3+1];
            double dx = (double)verts[(size_t)a*3+2] - (double)verts[(size_t)b*3+2];
            double len = sqrt(dz*dz + dy*dy + dx*dx);
            if (len < len_min) continue;      /* short edge: real connection, not a bridge */
            if (du > ratio * len) bad = 1;    /* long edge whose u is stretched: a bad link */
        }
        out_bad[f] = (uint8_t)bad;
        nbad += (size_t)bad;
    }
    if (out_n) *out_n = nbad;
    return 0;
}

int Ribbon_map_uv(Arena_T arena,
                  const float *cverts, size_t cnv,
                  const int32_t *cfaces, size_t cnf,
                  const float *cuv,
                  const float *verts, size_t nv,
                  const int32_t *faces, size_t nf,
                  float **out_uv, uint8_t **out_ok, size_t *out_fallback)
{
    assert(arena && out_uv && out_ok && out_fallback);
    *out_uv = NULL; *out_ok = NULL; *out_fallback = 0;
    if (cnv < 3 || cnf < 1 || nv < 1 || cuv == NULL) return -1;

    float   *uv = (float *)ARENA_ALLOC(arena, (long)(nv * 2 * sizeof(float)));
    uint8_t *ok = (uint8_t *)ARENA_CALLOC(arena, (long)nv, 1);

    Arena_Mark mark = Arena_save(arena);

    /* coarse vertex -> incident faces (CSR) */
    int32_t *voff = (int32_t *)ARENA_CALLOC(arena, (long)(cnv + 1), (long)sizeof(int32_t));
    for (size_t f = 0; f < cnf * 3; f++) voff[cfaces[f] + 1]++;
    for (size_t i = 0; i < cnv; i++) voff[i+1] = (int32_t)(voff[i+1] + voff[i]);
    int32_t *vfac = (int32_t *)ARENA_ALLOC(arena, (long)(cnf * 3 * sizeof(int32_t)));
    int32_t *vcur = (int32_t *)ARENA_ALLOC(arena, (long)(cnv * sizeof(int32_t)));
    memcpy(vcur, voff, cnv * sizeof(int32_t));
    for (size_t f = 0; f < cnf; f++)
        for (int m = 0; m < 3; m++)
            vfac[vcur[cfaces[f*3 + (size_t)m]]++] = (int32_t)f;

    /* v is the AXIAL coordinate = z - tmin (the scroll axis is Z), set on the
     * coarse mesh as (t - tmin) with t = vertex z. It is NOT a free parameter,
     * so the fine vertex's v must come from its OWN z, not the barycentric
     * coarse v (which is the z of the closest point on the DECIMATED surface --
     * that lands on a different slice near folds/steep z, injecting several vox
     * of spurious axial error and anisotropic UV stretch). Recover tmin from
     * the coarse map (cuv.v == cz - tmin is constant across coarse verts; median
     * a sample for robustness) and take v = z - tmin exactly for every fine
     * vertex -- including unmapped ones. The coarse map supplies only u. */
    double tmin = 0.0;
    {
        double samp[1024];
        size_t ns = 0, step = cnv > 1024 ? cnv / 1024 : 1;
        for (size_t i = 0; i < cnv && ns < 1024; i += step)
            samp[ns++] = (double)cverts[i*3+0] - (double)cuv[i*2+1];
        qsort(samp, ns, sizeof(double), cmp_dbl);
        tmin = ns ? samp[ns/2] : 0.0;
    }
    for (size_t i = 0; i < nv; i++)
        uv[i*2+1] = (float)((double)verts[i*3+0] - tmin);   /* v = z - tmin */

    KDTree_T tree = KDTree_new(arena, cverts, cnv);
    int32_t ballbuf[64];

    size_t unmapped = 0;
    double guard2 = RIB_XFER_R * RIB_XFER_R;
    for (size_t i = 0; i < nv; i++) {
        float q[3] = { verts[i*3+0], verts[i*3+1], verts[i*3+2] };
        double p[3] = { (double)q[0], (double)q[1], (double)q[2] };
        /* candidate faces: incident to coarse verts within a ball. The ball
         * radius stays below the 7-vox wrap clearance so candidates can only
         * be on the vertex's own wrap. */
        size_t nball = KDTree_ball_query(tree, q, 25.0f /* 5^2 */, ballbuf, 64);
        double bestd2 = 1e300, bu = 0, bv = 0;
        int32_t bestf = -1;
        double nvd2 = 1e300, nvu = 0.0;   /* nearest coarse vertex u (fallback) */
        for (size_t bq = 0; bq < nball; bq++) {
            int32_t cvid = ballbuf[bq];
            double vz = p[0]-(double)cverts[(size_t)cvid*3];
            double vy = p[1]-(double)cverts[(size_t)cvid*3+1];
            double vx = p[2]-(double)cverts[(size_t)cvid*3+2];
            double vd2 = vz*vz + vy*vy + vx*vx;
            if (vd2 < nvd2) { nvd2 = vd2; nvu = (double)cuv[(size_t)cvid*2+0]; }
            for (int32_t e = voff[cvid]; e < voff[cvid+1]; e++) {
                int32_t f = vfac[e];
                if (f == bestf) continue;
                const int32_t *fc = &cfaces[(size_t)f*3];
                double A[3] = { (double)cverts[(size_t)fc[0]*3], (double)cverts[(size_t)fc[0]*3+1], (double)cverts[(size_t)fc[0]*3+2] };
                double B[3] = { (double)cverts[(size_t)fc[1]*3], (double)cverts[(size_t)fc[1]*3+1], (double)cverts[(size_t)fc[1]*3+2] };
                double C[3] = { (double)cverts[(size_t)fc[2]*3], (double)cverts[(size_t)fc[2]*3+1], (double)cverts[(size_t)fc[2]*3+2] };
                double cp[3], tu, tv;
                double u0 = (double)cuv[(size_t)fc[0]*2+0];
                double u1 = (double)cuv[(size_t)fc[1]*2+0];
                double u2 = (double)cuv[(size_t)fc[2]*2+0];
                double uspan = 0.0, diam = 0.0, e01, e12, e20;
                closest_on_tri(p, A, B, C, cp, &tu, &tv);
                double dz = p[0]-cp[0], dy = p[1]-cp[1], dx = p[2]-cp[2];
                double d2 = dz*dz + dy*dy + dx*dx;
                /* reject faces that STRADDLE a u-discontinuity: a valid coarse
                 * face has u varying like its 3D size (u = arc length), so a
                 * face whose corner u-span far exceeds its 3D diameter bridges
                 * two chains/wraps and would smear when interpolated. Mirrors
                 * the raster smear gate (u-span > max(4*diam, RIB_XFER floor)). */
                uspan = fabs(u0-u1); if (fabs(u1-u2) > uspan) uspan = fabs(u1-u2);
                if (fabs(u2-u0) > uspan) uspan = fabs(u2-u0);
                e01 = (A[0]-B[0])*(A[0]-B[0]) + (A[1]-B[1])*(A[1]-B[1]) + (A[2]-B[2])*(A[2]-B[2]);
                e12 = (B[0]-C[0])*(B[0]-C[0]) + (B[1]-C[1])*(B[1]-C[1]) + (B[2]-C[2])*(B[2]-C[2]);
                e20 = (C[0]-A[0])*(C[0]-A[0]) + (C[1]-A[1])*(C[1]-A[1]) + (C[2]-A[2])*(C[2]-A[2]);
                diam = sqrt(e01 > e12 ? (e01 > e20 ? e01 : e20) : (e12 > e20 ? e12 : e20));
                if (uspan > (4.0*diam > 25.0 ? 4.0*diam : 25.0)) continue;
                if (d2 < bestd2) { bestd2 = d2; bestf = f; bu = tu; bv = tv; }
            }
        }
        if (bestf >= 0 && bestd2 <= guard2) {
            const int32_t *fc = &cfaces[(size_t)bestf*3];
            double w0 = 1.0 - bu - bv, w1 = bu, w2 = bv;
            uv[i*2+0] = (float)(w0 * (double)cuv[(size_t)fc[0]*2+0]
                              + w1 * (double)cuv[(size_t)fc[1]*2+0]
                              + w2 * (double)cuv[(size_t)fc[2]*2+0]);
            /* v already set exactly from z above; keep it */
            ok[i] = 1;
        } else if (nvd2 <= guard2) {
            /* every nearby coarse face straddled a u-jump (or none was in
             * range): fall back to the nearest coarse VERTEX's u, which is a
             * single clean value -- never an interpolation across the jump. */
            uv[i*2+0] = (float)nvu;
            ok[i] = 1;
        } else {
            uv[i*2+0] = 0.0f;
            unmapped++;
        }
    }

    /* neighbor fill (median, robust) on the original mesh -- U ONLY; v stays
     * exact (z - tmin) for unmapped verts too. */
    if (unmapped > 0 && nf > 0) {
        CSR_T adj = CSR_from_faces(arena, faces, nf, nv);
        const int32_t *aoff = CSR_offset(adj);
        const int32_t *atgt = CSR_target(adj);
        double nbu[64];
        for (int round = 0; round < RIB_UVFILL_ROUNDS && unmapped > 0; round++) {
            size_t fixed = 0;
            for (size_t i = 0; i < nv; i++) {
                if (ok[i]) continue;
                int cnt = 0;
                for (int32_t e = aoff[i]; e < aoff[i+1] && cnt < 64; e++) {
                    int32_t j = atgt[e];
                    if (ok[j] == 1) { nbu[cnt] = (double)uv[(size_t)j*2]; cnt++; }
                }
                if (cnt > 0) {
                    qsort(nbu, (size_t)cnt, sizeof(double), cmp_dbl);
                    uv[i*2+0] = (float)nbu[cnt / 2];
                    ok[i] = 2;
                    fixed++;
                }
            }
            for (size_t i = 0; i < nv; i++)
                if (ok[i] == 2) ok[i] = 1;
            unmapped -= fixed;
            if (fixed == 0) break;
        }
    }
    size_t nfb = 0;
    for (size_t i = 0; i < nv; i++) if (!ok[i]) nfb++;

    Arena_restore(arena, mark);
    *out_uv = uv;
    *out_ok = ok;
    *out_fallback = nfb;
    return 0;
}

/* ============================================================================
 * Self-test.
 * ==========================================================================*/

typedef int (*SkipFn)(int i, int j, void *user);

/* Open Archimedean spiral sheet around Z through (cy,cx):
 * r(a) = r0 + g*a, a in [0, 2*pi*turns], nphi x nh grid, verts (z,y,x). */
static void st_build_spiral(Arena_T arena, int nphi, int nh, double r0,
                            double g, double turns, double Hgt,
                            double cy, double cx, SkipFn skip, void *user,
                            float **out_v, size_t *out_nv,
                            int32_t **out_f, size_t *out_nf)
{
    size_t nvv = (size_t)nphi * (size_t)nh;
    float   *v = (float *)ARENA_ALLOC(arena, (long)(nvv * 3 * sizeof(float)));
    int32_t *f = (int32_t *)ARENA_ALLOC(arena, (long)((size_t)(nphi-1) * (size_t)(nh-1) * 6 * sizeof(int32_t)));
    double amax = turns * 2.0 * M_PI;
    for (int j = 0; j < nh; j++) {
        double z = Hgt * (double)j / (double)(nh - 1);
        for (int i = 0; i < nphi; i++) {
            double a = amax * (double)i / (double)(nphi - 1);
            double rr = r0 + g * a;
            size_t idx = (size_t)j * (size_t)nphi + (size_t)i;
            v[idx*3 + 0] = (float)z;
            v[idx*3 + 1] = (float)(cy + rr * sin(a));
            v[idx*3 + 2] = (float)(cx + rr * cos(a));
        }
    }
    size_t fi = 0;
    for (int j = 0; j < nh - 1; j++) {
        for (int i = 0; i < nphi - 1; i++) {
            if (skip && skip(i, j, user)) continue;
            int32_t a = (int32_t)((size_t)j * (size_t)nphi + (size_t)i);
            int32_t b = a + 1;
            int32_t c = (int32_t)((size_t)(j+1) * (size_t)nphi + (size_t)i);
            int32_t d = c + 1;
            f[fi*3+0] = a; f[fi*3+1] = b; f[fi*3+2] = c; fi++;
            f[fi*3+0] = b; f[fi*3+1] = d; f[fi*3+2] = c; fi++;
        }
    }
    *out_v = v; *out_nv = nvv; *out_f = f; *out_nf = fi;
}

static double st_spiral_arclen(double r0, double g, double a)
{
    int n = 4000;
    double sum = 0.0, hstep = a / (double)n;
    for (int i = 0; i < n; i++) {
        double x0 = (double)i * hstep, x1 = x0 + hstep;
        double f0 = sqrt((r0 + g*x0)*(r0 + g*x0) + g*g);
        double f1 = sqrt((r0 + g*x1)*(r0 + g*x1) + g*g);
        sum += 0.5 * (f0 + f1) * hstep;
    }
    return sum;
}

typedef struct { int i0, i1, j0, j1; int bi0, bi1; } StHoles;
static int st_skip_holes(int i, int j, void *user)
{
    const StHoles *H = (const StHoles *)user;
    if (i >= H->i0 && i <= H->i1 && j >= H->j0 && j <= H->j1) return 1;
    if (i >= H->bi0 && i <= H->bi1) return 1;
    return 0;
}

/* skip two single columns (for the fusion-wall case) */
typedef struct { int ia, ib; } StTwoCols;
static int st_skip_twocols(int i, int j, void *user)
{
    const StTwoCols *T = (const StTwoCols *)user;
    (void)j;
    return i == T->ia || i == T->ib;
}

static int st_check(int cond, const char *what, int *fails)
{
    if (!cond) {
        fprintf(stderr, "[ribbon selftest]   FAIL: %s\n", what);
        (*fails)++;
    }
    return cond;
}

int Ribbon_selftest(void)
{
    int fails = 0;
    Arena_T arena = Arena_new();

    /* (1) chain builder on a flat strip: 2 planes -> 2 open chains, u == x. */
    {
        enum { NX = 6 };
        float v[NX*2*3]; int32_t f[(NX-1)*2*3];
        for (int i = 0; i < NX; i++) {
            v[i*3+0] = 0.0f;      v[i*3+1] = 0.0f; v[i*3+2] = (float)i;
            v[(NX+i)*3+0] = 4.0f; v[(NX+i)*3+1] = 0.0f; v[(NX+i)*3+2] = (float)i;
        }
        size_t nfc = 0;
        for (int i = 0; i < NX-1; i++) {
            f[nfc*3+0] = i; f[nfc*3+1] = i+1; f[nfc*3+2] = NX+i; nfc++;
            f[nfc*3+0] = i+1; f[nfc*3+1] = NX+i+1; f[nfc*3+2] = NX+i; nfc++;
        }
        RibbonOpts o; RibbonOpts_default(&o);
        o.axis_point[1] = -60.0f;   /* keep the axis far from the strip */
        RibbonResult R;
        int rc = Ribbon_run(arena, v, NX*2, f, nfc, &o, &R);
        fprintf(stderr, "[ribbon selftest] (1) strip: rc=%d slices=%d chains=%d "
                "samples=%zu duds_max=%.5f u_span=%.3f\n",
                rc, R.n_slices, R.n_chains, R.n_samples, R.duds_err_max, R.u_span);
        st_check(rc == 0, "strip rc", &fails);
        st_check(R.n_slices == 2, "strip 2 slices", &fails);
        st_check(R.n_chains == 2, "strip 2 chains", &fails);
        st_check(R.duds_err_max < 1e-3, "strip u == arc length", &fails);
        st_check(fabs(R.u_span - 5.0) < 0.05, "strip u_span == 5", &fails);
        int okuv = 1;
        for (int i = 1; i < NX; i++) {
            double du_ = fabs((double)R.uv[i*2] - (double)R.uv[(i-1)*2]);
            if (fabs(du_ - 1.0) > 0.08) okuv = 0;
        }
        st_check(okuv, "strip vt |du| == x spacing", &fails);
    }

    /* (2) analytic spiral: u == true arc length within 2%; isolines vertical. */
    {
        double r0 = 15.0, g = 2.0, turns = 3.0, Hgt = 40.0, cy = 100.0, cx = 50.0;
        int nphi = 600, nh = 21;
        float *v; int32_t *f; size_t nvv, nfc;
        st_build_spiral(arena, nphi, nh, r0, g, turns, Hgt, cy, cx, NULL, NULL,
                        &v, &nvv, &f, &nfc);
        RibbonOpts o; RibbonOpts_default(&o);
        o.axis_point[1] = (float)cy; o.axis_point[2] = (float)cx;
        RibbonResult R;
        int rc = Ribbon_run(arena, v, nvv, f, nfc, &o, &R);
        double amax = turns * 2.0 * M_PI;
        double s_total = st_spiral_arclen(r0, g, amax);
        fprintf(stderr, "[ribbon selftest] (2) spiral: rc=%d slices=%d chains=%d "
                "samples=%zu pairs=%zu cover=%.2f groups=%d conf=%zu bridge=%zu "
                "duds_mean=%.4f u_span=%.1f (s_true=%.1f) turns=%.2f grid=%zux%zu\n",
                rc, R.n_slices, R.n_chains, R.n_samples, R.n_pairs, R.match_cover,
                R.w_groups, R.w_conflicts, R.bridge_cuts,
                R.duds_err_mean, R.u_span, s_total,
                R.phi_span_turns, R.nu, R.nk);
        st_check(rc == 0, "spiral rc", &fails);
        st_check(R.match_cover > 0.90, "spiral match cover > 90%", &fails);
        st_check(R.duds_err_mean < 0.02, "spiral duds mean < 2%", &fails);
        st_check(fabs(R.u_span - s_total) / s_total < 0.03, "spiral u_span ~= s_true", &fails);
        st_check(fabs(R.phi_span_turns - turns) < 0.2, "spiral turns ~= 3", &fails);
        if (rc == 0) {
            int j = nh / 2;
            double e_max = 0.0;
            double u_at0 = (double)R.uv[((size_t)j*(size_t)nphi)*2];
            /* winding sense is a gauge: detect the chart's direction first */
            double sdir = (double)R.uv[((size_t)j*(size_t)nphi + (size_t)(nphi-1))*2]
                          > u_at0 ? 1.0 : -1.0;
            for (int i = 0; i < nphi; i += 10) {
                double a = amax * (double)i / (double)(nphi - 1);
                double s = st_spiral_arclen(r0, g, a);
                double uu = sdir * ((double)R.uv[((size_t)j*(size_t)nphi + (size_t)i)*2] - u_at0);
                double e = fabs(uu - s) / s_total;
                if (e > e_max) e_max = e;
            }
            fprintf(stderr, "[ribbon selftest]     mid-row |u - s_true|/s_total max = %.4f\n", e_max);
            st_check(e_max < 0.02, "spiral u tracks s within 2%", &fails);
            double spread_max = 0.0;
            for (int i = 50; i < nphi - 50; i += 100) {
                double lo = 1e300, hi = -1e300;
                for (int jj = 3; jj < nh - 3; jj++) {
                    double uu = (double)R.uv[((size_t)jj*(size_t)nphi + (size_t)i)*2];
                    if (uu < lo) lo = uu;
                    if (uu > hi) hi = uu;
                }
                if (hi - lo > spread_max) spread_max = hi - lo;
            }
            fprintf(stderr, "[ribbon selftest]     isoline u spread across z max = %.3f vox\n", spread_max);
            st_check(spread_max < 1.5, "spiral isolines vertical (<1.5 vox)", &fails);
            st_check(R.nu > 10 && R.nk > 10, "spiral grid nonempty", &fails);
            double worst = 0.0; int checked = 0;
            for (size_t kk = 2; kk < R.nk - 2 && checked < 150; kk += 3) {
                for (size_t jj = 2; jj < R.nu - 2 && checked < 150; jj += 17) {
                    if (!R.grid_valid[kk*R.nu + jj]) continue;
                    const float *gp = &R.grid_pos[(kk*R.nu + jj)*3];
                    if (isnan((double)gp[0])) continue;
                    double best = 1e300;
                    for (size_t vi = 0; vi < nvv; vi++) {
                        double dz = (double)gp[0]-(double)v[vi*3];
                        double dy = (double)gp[1]-(double)v[vi*3+1];
                        double dx = (double)gp[2]-(double)v[vi*3+2];
                        double d2 = dz*dz + dy*dy + dx*dx;
                        if (d2 < best) best = d2;
                    }
                    best = sqrt(best);
                    if (best > worst) worst = best;
                    checked++;
                }
            }
            fprintf(stderr, "[ribbon selftest]     grid-to-mesh max dist over %d cells = %.3f vox\n",
                    checked, worst);
            st_check(checked > 50 && worst < 2.0, "ribbon hugs surface (<2 vox)", &fails);
        }
    }

    /* (2b) MANY-turn spiral GROUND TRUTH: a perfect 12-turn spiral scroll; u
     * must track the analytic arc length across ALL turns and the winding must
     * span every turn -- the regression for "~30 wraps collapsed into ~6". */
    {
        double r0 = 15.0, g = 2.0, turns = 12.0, Hgt = 12.0, cy = 300.0, cx = 300.0;
        int nphi = 1500, nh = 7;
        float *v; int32_t *f; size_t nvv, nfc;
        st_build_spiral(arena, nphi, nh, r0, g, turns, Hgt, cy, cx, NULL, NULL,
                        &v, &nvv, &f, &nfc);
        RibbonOpts o; RibbonOpts_default(&o);
        o.axis_point[1] = (float)cy; o.axis_point[2] = (float)cx;
        RibbonResult R;
        int rc = Ribbon_run(arena, v, nvv, f, nfc, &o, &R);
        double amax = turns * 2.0 * M_PI;
        double s_total = st_spiral_arclen(r0, g, amax);
        fprintf(stderr, "[ribbon selftest] (2b) 12-turn spiral: rc=%d chains=%d "
                "wraps=%d pitch=%.2f turns=%.2f u_span=%.0f (s_true=%.0f) duds=%.4f\n",
                rc, R.n_chains, R.w_groups, R.pitch_used, R.phi_span_turns,
                R.u_span, s_total, R.duds_err_mean);
        st_check(rc == 0, "12-turn rc", &fails);
        st_check(fabs(R.phi_span_turns - turns) < 0.5, "12-turn winding spans all turns", &fails);
        st_check(fabs(R.u_span - s_total) / s_total < 0.05, "12-turn u_span ~= s_true", &fails);
        st_check(R.duds_err_mean < 0.05, "12-turn duds mean < 5%", &fails);
        if (rc == 0) {
            int j = nh / 2;
            double u_at0 = (double)R.uv[((size_t)j*(size_t)nphi)*2];
            double sdir = (double)R.uv[((size_t)j*(size_t)nphi + (size_t)(nphi-1))*2]
                          > u_at0 ? 1.0 : -1.0;
            double e_max = 0.0;
            for (int i = 0; i < nphi; i += 25) {
                double a = amax * (double)i / (double)(nphi - 1);
                double s = st_spiral_arclen(r0, g, a);
                double uu = sdir * ((double)R.uv[((size_t)j*(size_t)nphi + (size_t)i)*2] - u_at0);
                double e = fabs(uu - s) / s_total;
                if (e > e_max) e_max = e;
            }
            fprintf(stderr, "[ribbon selftest]     12-turn |u - s_true|/s_total max = %.4f\n", e_max);
            st_check(e_max < 0.05, "12-turn u tracks analytic arc length within 5%", &fails);
        }
    }

    /* (2e) GLOBAL-FRAME EMISSION (the scroll_whole per-cube contract):
     * emit_global is the SAME parameterization in an absolute frame ((u,v)
     * differ from the classic run only by per-run constants; phi filled and
     * tracking the analytic winding), the frame is invariant to sliding the
     * axis point along the axis (u, phi identical; v shifts by exactly the
     * slide), and pinning sense/pitch/spiral to a seed calibration reproduces
     * the auto result on clean geometry. */
    {
        double r0 = 15.0, g = 2.0, turns = 3.0, Hgt = 40.0, cy = 100.0, cx = 50.0;
        int nphi = 600, nh = 21;
        float *v; int32_t *f; size_t nvv, nfc;
        st_build_spiral(arena, nphi, nh, r0, g, turns, Hgt, cy, cx, NULL, NULL,
                        &v, &nvv, &f, &nfc);
        RibbonOpts oA; RibbonOpts_default(&oA);
        oA.axis_point[1] = (float)cy; oA.axis_point[2] = (float)cx;
        RibbonResult A;
        int rcA = Ribbon_run(arena, v, nvv, f, nfc, &oA, &A);

        RibbonOpts oB = oA; oB.emit_global = 1; oB.pin_orient = 1;
        RibbonResult B;
        int rcB = Ribbon_run(arena, v, nvv, f, nfc, &oB, &B);
        st_check(rcA == 0 && rcB == 0, "glob rc", &fails);
        st_check(B.pitch_source == RIB_PITCH_ESTIMATED,
                 "glob auto pitch was measured", &fails);
        st_check(B.phi != NULL && A.phi == NULL, "glob phi only when emit_global", &fails);
        st_check(B.group != NULL && A.group == NULL, "glob group only when emit_global", &fails);
        if (rcB == 0 && B.group != NULL) {
            int gok = 1;
            for (size_t i = 0; i < nvv; i++) {
                if (B.uv_ok[i] && (B.group[i] < 0 || B.group[i] >= B.w_prior_groups))
                    gok = 0;
            }
            st_check(gok, "glob group ids in [0, w_prior_groups)", &fails);
        }

        /* A supplied graph winding is authoritative, including its absolute
         * gauge. Shift a clean reference by four whole turns and make sure the
         * slice/interpolate/transfer path preserves that shift instead of
         * silently rebuilding winding from the nearest-pair graph. */
        if (rcB == 0 && B.phi != NULL) {
            double ref_shift = 8.0 * M_PI;
            float *ref_phi = (float *)ARENA_ALLOC(
                arena, (long)(nvv * sizeof(float)));
            for (size_t i = 0; i < nvv; i++)
                ref_phi[i] = (float)((double)B.phi[i] + ref_shift);
            RibbonOpts oE = oB;
            oE.reference_phi = ref_phi;
            RibbonResult E;
            int rcE = Ribbon_run(arena, v, nvv, f, nfc, &oE, &E);
            st_check(rcE == 0 && E.phi != NULL, "reference-phi rc", &fails);
            if (rcE == 0 && E.phi != NULL) {
                double de_lo = 1e300, de_hi = -1e300, de_sum = 0.0;
                size_t ne = 0;
                for (size_t i = 0; i < nvv; i++) {
                    if (!B.uv_ok[i] || !E.uv_ok[i]) continue;
                    double de = (double)E.phi[i] - (double)B.phi[i];
                    if (de < de_lo) de_lo = de;
                    if (de > de_hi) de_hi = de;
                    de_sum += de; ne++;
                }
                fprintf(stderr, "[ribbon selftest]     reference-phi shift "
                        "mean=%.4f range=[%.4f,%.4f] (want %.4f)\n",
                        ne ? de_sum / (double)ne : 0.0, de_lo, de_hi, ref_shift);
                st_check(ne > nvv / 2, "reference-phi enough mapped verts", &fails);
                st_check(de_hi - de_lo < 0.1, "reference-phi gauge is constant", &fails);
                st_check(ne > 0 && fabs(de_sum / (double)ne - ref_shift) < 0.1,
                         "reference-phi gauge is preserved", &fails);
            }
        }

        if (rcA == 0 && rcB == 0 && B.phi != NULL) {
            /* sign of B's u relative to A's (A may have canonical-flipped) */
            size_t k0 = 0, k1 = (size_t)(nh / 2) * (size_t)nphi + (size_t)(nphi - 1);
            double sgn = ((double)B.uv[k1*2] - (double)B.uv[k0*2])
                       * ((double)A.uv[k1*2] - (double)A.uv[k0*2]) >= 0.0 ? 1.0 : -1.0;
            double du_lo = 1e300, du_hi = -1e300, dv_lo = 1e300, dv_hi = -1e300;
            double du_sum = 0.0;
            size_t nok = 0;
            for (size_t i = 0; i < nvv; i++) {
                if (!A.uv_ok[i] || !B.uv_ok[i]) continue;
                double du_ = (double)B.uv[i*2]   - sgn * (double)A.uv[i*2];
                double dv_ = (double)B.uv[i*2+1] - (double)A.uv[i*2+1];
                if (du_ < du_lo) du_lo = du_;
                if (du_ > du_hi) du_hi = du_;
                if (dv_ < dv_lo) dv_lo = dv_;
                if (dv_ > dv_hi) dv_hi = dv_;
                du_sum += du_; nok++;
            }
            fprintf(stderr, "[ribbon selftest] (2e) glob-vs-classic: sgn=%+.0f "
                    "du const in [%.4f,%.4f] dv const in [%.4f,%.4f] u_origin=%.3f\n",
                    sgn, du_lo, du_hi, dv_lo, dv_hi, B.u_origin);
            st_check(nok > nvv / 2, "glob enough mapped verts", &fails);
            st_check(du_hi - du_lo < 0.05, "glob u == classic u + const", &fails);
            st_check(dv_hi - dv_lo < 0.05, "glob v == classic v + const", &fails);
            if (sgn > 0.0)
                st_check(fabs(du_sum / (double)nok - B.u_origin) < 0.05,
                         "glob u const == u_origin", &fails);
            /* phi tracks the analytic winding angle up to one global const */
            int j = nh / 2;
            double amax = turns * 2.0 * M_PI;
            double e_phi = 0.0;
            size_t ref = (size_t)j * (size_t)nphi + 50;
            for (int i = 50; i < nphi - 50; i += 50) {
                size_t vi = (size_t)j * (size_t)nphi + (size_t)i;
                if (!B.uv_ok[vi] || !B.uv_ok[ref]) continue;
                double da_true = amax * (double)(i - 50) / (double)(nphi - 1);
                double da_got  = fabs((double)B.phi[vi] - (double)B.phi[ref]);
                double e = fabs(da_got - da_true);
                if (e > e_phi) e_phi = e;
            }
            fprintf(stderr, "[ribbon selftest]     phi-vs-analytic max err = %.3f rad\n", e_phi);
            st_check(e_phi < 0.3, "glob phi tracks analytic winding", &fails);
        }

        /* axis slide: same mesh, axis_point moved +7 along the axis. Slicing is
         * anchored to the mesh's own t-extent, so samples are identical: u and
         * phi must not move; v = t must drop by exactly 7. */
        RibbonOpts oC = oB; oC.axis_point[0] += 7.0f;
        RibbonResult C;
        int rcC = Ribbon_run(arena, v, nvv, f, nfc, &oC, &C);
        st_check(rcC == 0, "glob slide rc", &fails);
        if (rcB == 0 && rcC == 0 && B.phi != NULL && C.phi != NULL) {
            double eu = 0.0, ev = 0.0, ep = 0.0;
            for (size_t i = 0; i < nvv; i++) {
                if (!B.uv_ok[i] || !C.uv_ok[i]) continue;
                double du_ = fabs((double)C.uv[i*2]   - (double)B.uv[i*2]);
                double dv_ = fabs(((double)C.uv[i*2+1] + 7.0) - (double)B.uv[i*2+1]);
                double dp_ = fabs((double)C.phi[i] - (double)B.phi[i]);
                if (du_ > eu) eu = du_;
                if (dv_ > ev) ev = dv_;
                if (dp_ > ep) ep = dp_;
            }
            fprintf(stderr, "[ribbon selftest]     axis-slide max |du|=%.5f |dv-7|=%.5f |dphi|=%.5f\n",
                    eu, ev, ep);
            st_check(eu < 1e-3, "axis slide: u invariant", &fails);
            st_check(ev < 1e-3, "axis slide: v shifts by the slide", &fails);
            st_check(ep < 1e-4, "axis slide: phi invariant", &fails);
        }

        /* seed-calibration pins (the scroll_whole flow: pitch = |spiral_b|,
         * sense = sign(spiral_b), (a,b) = seed fit) reproduce the auto run */
        RibbonOpts oD = oB;
        oD.wrap_spacing  = (float)fabs(B.spiral_b);
        oD.winding_sense = B.spiral_b > 0.0 ? 1 : -1;
        oD.spiral_a = B.spiral_a; oD.spiral_b = B.spiral_b;
        RibbonResult D;
        int rcD = Ribbon_run(arena, v, nvv, f, nfc, &oD, &D);
        st_check(rcD == 0, "glob pinned rc", &fails);
        st_check(D.pitch_source == RIB_PITCH_PINNED,
                 "glob pinned pitch provenance", &fails);
        if (rcB == 0 && rcD == 0 && B.phi != NULL && D.phi != NULL) {
            st_check(fabs(D.pitch_used - fabs(B.spiral_b)) < 1e-9, "pinned pitch echoed", &fails);
            st_check(D.spiral_a == B.spiral_a && D.spiral_b == B.spiral_b,
                     "pinned spiral echoed", &fails);
            double eu = 0.0, ep = 0.0;
            for (size_t i = 0; i < nvv; i++) {
                if (!B.uv_ok[i] || !D.uv_ok[i]) continue;
                double du_ = fabs((double)D.uv[i*2] - (double)B.uv[i*2]);
                double dp_ = fabs((double)D.phi[i] - (double)B.phi[i]);
                if (du_ > eu) eu = du_;
                if (dp_ > ep) ep = dp_;
            }
            fprintf(stderr, "[ribbon selftest]     pinned-vs-auto max |du|=%.4f |dphi|=%.4f\n",
                    eu, ep);
            st_check(eu < 0.5, "pinned calibration reproduces auto u", &fails);
            st_check(ep < 1e-3, "pinned calibration reproduces auto phi", &fails);
        }
    }

    /* (2c) MANY disconnected concentric wraps (16 rings, pitch 10): each is a
     * separate component, so the OLD vote graph left most groups unreached and
     * collapsed them onto W=0. The radius anchor must give all 16 consecutive
     * windings -- wrap count == maximal winding, no collapse, u monotone in r. */
    {
        enum { NU = 96, NH = 5, NW = 16 };
        double R0 = 30.0, DR = 10.0, Hgt = 8.0, cy = 400.0, cx = 400.0;
        size_t nvv = (size_t)NU * NH * NW;
        float *v = (float *)ARENA_ALLOC(arena, (long)(nvv * 3 * sizeof(float)));
        int32_t *f = (int32_t *)ARENA_ALLOC(arena, (long)((size_t)(NU-1)*(NH-1)*NW*2*3 * sizeof(int32_t)));
        size_t nfc = 0;
        for (int w = 0; w < NW; w++) {
            double R_ = R0 + DR * (double)w;
            size_t base = (size_t)w * (size_t)NU * NH;
            for (int j = 0; j < NH; j++)
                for (int i = 0; i < NU; i++) {
                    double a = 2.0 * M_PI * (double)i / (double)NU;   /* open seam */
                    size_t idx = base + (size_t)j * NU + (size_t)i;
                    v[idx*3+0] = (float)(Hgt * (double)j / (double)(NH-1));
                    v[idx*3+1] = (float)(cy + R_ * sin(a));
                    v[idx*3+2] = (float)(cx + R_ * cos(a));
                }
            for (int j = 0; j < NH - 1; j++)
                for (int i = 0; i < NU - 1; i++) {
                    int32_t a0 = (int32_t)(base + (size_t)j*NU + (size_t)i);
                    int32_t b0 = a0 + 1;
                    int32_t c0 = (int32_t)(base + (size_t)(j+1)*NU + (size_t)i);
                    int32_t d0 = c0 + 1;
                    f[nfc*3+0]=a0; f[nfc*3+1]=b0; f[nfc*3+2]=c0; nfc++;
                    f[nfc*3+0]=b0; f[nfc*3+1]=d0; f[nfc*3+2]=c0; nfc++;
                }
        }
        RibbonOpts o; RibbonOpts_default(&o);
        o.axis_point[1] = (float)cy; o.axis_point[2] = (float)cx;
        RibbonResult R;
        int rc = Ribbon_run(arena, v, nvv, f, nfc, &o, &R);
        int mono = (rc == 0);
        double mu[NW];
        if (rc == 0) {
            for (int w = 0; w < NW; w++) {
                size_t base = (size_t)w * (size_t)NU * NH;
                double s = 0.0; long mc = 0;
                for (size_t i = base; i < base + (size_t)NU*NH; i++) { s += (double)R.uv[i*2]; mc++; }
                mu[w] = s / (double)(mc ? mc : 1);
            }
            double d0 = mu[1] - mu[0];
            for (int w = 1; w < NW; w++) {
                double d = mu[w] - mu[w-1];
                if ((d > 0) != (d0 > 0) || fabs(d) < 50.0) mono = 0;
            }
        }
        fprintf(stderr, "[ribbon selftest] (2c) %d disconnected wraps: rc=%d wraps=%d "
                "unreached=%d pitch=%.2f mono=%d\n",
                NW, rc, R.w_groups, R.w_unreached, R.pitch_used, mono);
        st_check(rc == 0, "16wrap rc", &fails);
        st_check(R.w_groups == NW, "16wrap NO collapse (wrap count == 16)", &fails);
        st_check(R.w_unreached == 0, "16wrap all radius-consistent", &fails);
        st_check(mono, "16wrap u monotone in radius, spaced >~1 band", &fails);
    }

    /* (2d) DELAMINATION: a small flap sitting at the PARENT wrap's radius must
     * be trapped in the SAME winding (same u band), NOT promoted to the outer
     * wrap. Three components: parent ring, outer ring (+1 pitch), flap ~parent. */
    {
        enum { NU = 96, NH = 5, FN = 12, FH = 5 };
        double Rp = 100.0, Ro = 110.0, Hgt = 8.0, cy = 500.0, cx = 500.0;
        size_t nring = (size_t)NU * NH, nflap = (size_t)FN * FH;
        size_t nvv = 2*nring + nflap;
        float *v = (float *)ARENA_ALLOC(arena, (long)(nvv * 3 * sizeof(float)));
        int32_t *f = (int32_t *)ARENA_ALLOC(arena,
            (long)(((size_t)(NU-1)*(NH-1)*2*2 + (size_t)(FN-1)*(FH-1)*2) * 3 * sizeof(int32_t)));
        size_t nfc = 0;
        double Rr[2] = { Rp, Ro };
        for (int w = 0; w < 2; w++) {
            size_t base = (size_t)w * nring;
            for (int j = 0; j < NH; j++)
                for (int i = 0; i < NU; i++) {
                    double a = 2.0 * M_PI * (double)i / (double)NU;
                    size_t idx = base + (size_t)j*NU + (size_t)i;
                    v[idx*3+0]=(float)(Hgt*(double)j/(double)(NH-1));
                    v[idx*3+1]=(float)(cy + Rr[w]*sin(a));
                    v[idx*3+2]=(float)(cx + Rr[w]*cos(a));
                }
            for (int j=0;j<NH-1;j++) for (int i=0;i<NU-1;i++){
                int32_t a0=(int32_t)(base+(size_t)j*NU+(size_t)i), b0=a0+1;
                int32_t c0=(int32_t)(base+(size_t)(j+1)*NU+(size_t)i), d0=c0+1;
                f[nfc*3+0]=a0;f[nfc*3+1]=b0;f[nfc*3+2]=c0;nfc++;
                f[nfc*3+0]=b0;f[nfc*3+1]=d0;f[nfc*3+2]=c0;nfc++;
            }
        }
        size_t fbase = 2*nring;
        for (int j=0;j<FH;j++)
            for (int i=0;i<FN;i++){
                double a = 0.6 + 0.5*(double)i/(double)(FN-1);
                size_t idx = fbase + (size_t)j*FN + (size_t)i;
                v[idx*3+0]=(float)(Hgt*(double)j/(double)(FH-1));
                v[idx*3+1]=(float)(cy + (Rp+1.0)*sin(a));
                v[idx*3+2]=(float)(cx + (Rp+1.0)*cos(a));
            }
        for (int j=0;j<FH-1;j++) for (int i=0;i<FN-1;i++){
            int32_t a0=(int32_t)(fbase+(size_t)j*FN+(size_t)i), b0=a0+1;
            int32_t c0=(int32_t)(fbase+(size_t)(j+1)*FN+(size_t)i), d0=c0+1;
            f[nfc*3+0]=a0;f[nfc*3+1]=b0;f[nfc*3+2]=c0;nfc++;
            f[nfc*3+0]=b0;f[nfc*3+1]=d0;f[nfc*3+2]=c0;nfc++;
        }
        RibbonOpts o; RibbonOpts_default(&o);
        o.axis_point[1]=(float)cy; o.axis_point[2]=(float)cx;
        o.wrap_spacing = 10.0f;   /* known pitch */
        RibbonResult R; int rc = Ribbon_run(arena, v, nvv, f, nfc, &o, &R);
        double up_lo=1e300,up_hi=-1e300, uo_lo=1e300,uo_hi=-1e300, uf_lo=1e300,uf_hi=-1e300;
        if (rc==0){
            for (size_t i=0;i<nring;i++){ double u=(double)R.uv[i*2]; if(u<up_lo)up_lo=u; if(u>up_hi)up_hi=u; }
            for (size_t i=nring;i<2*nring;i++){ double u=(double)R.uv[i*2]; if(u<uo_lo)uo_lo=u; if(u>uo_hi)uo_hi=u; }
            for (size_t i=fbase;i<nvv;i++){ double u=(double)R.uv[i*2]; if(u<uf_lo)uf_lo=u; if(u>uf_hi)uf_hi=u; }
        }
        double fmid = 0.5*(uf_lo+uf_hi);
        int in_parent = (rc==0) && (fmid > up_lo - 25.0 && fmid < up_hi + 25.0);
        double dpar = fabs(fmid - 0.5*(up_lo+up_hi));
        double dout = fabs(fmid - 0.5*(uo_lo+uo_hi));
        fprintf(stderr, "[ribbon selftest] (2d) delamination: rc=%d parent_u=[%.0f,%.0f] "
                "outer_u=[%.0f,%.0f] flap_mid=%.0f d_parent=%.0f d_outer=%.0f\n",
                rc, up_lo, up_hi, uo_lo, uo_hi, fmid, dpar, dout);
        st_check(rc==0, "delam rc", &fails);
        st_check(in_parent && dpar < dout, "delamination flap trapped in parent winding", &fails);
    }

    /* (3) punched spiral: fragments stitch, arc length still tracks truth. */
    {
        double r0 = 15.0, g = 2.0, turns = 3.0, Hgt = 40.0, cy = 100.0, cx = 50.0;
        int nphi = 600, nh = 21;
        StHoles H = { 150, 180, 8, 12,  300, 303 };
        float *v; int32_t *f; size_t nvv, nfc;
        st_build_spiral(arena, nphi, nh, r0, g, turns, Hgt, cy, cx,
                        st_skip_holes, &H, &v, &nvv, &f, &nfc);
        RibbonOpts o; RibbonOpts_default(&o);
        o.axis_point[1] = (float)cy; o.axis_point[2] = (float)cx;
        RibbonResult R;
        int rc = Ribbon_run(arena, v, nvv, f, nfc, &o, &R);
        double amax = turns * 2.0 * M_PI;
        fprintf(stderr, "[ribbon selftest] (3) holes: rc=%d chains=%d (slices=%d) "
                "cont_pairs=%zu groups=%d qp_comps=%d duds_mean=%.4f fallback=%zu\n",
                rc, R.n_chains, R.n_slices, R.n_cont_pairs, R.w_groups, R.n_qp_comps,
                R.duds_err_mean, R.uv_fallback);
        st_check(rc == 0, "holes rc", &fails);
        st_check(R.n_chains > R.n_slices, "holes fragmented", &fails);
        st_check(R.n_cont_pairs > 0, "holes continuation pairs found", &fails);
        st_check(R.duds_err_mean < 0.03, "holes duds mean < 3%", &fails);
        st_check(R.uv_fallback < nvv / 20, "holes <5% uv fallback", &fails);
        if (rc == 0) {
            int j = nh / 2, iA = 280, iB = 330;
            double aA = amax * (double)iA / (double)(nphi-1);
            double aB = amax * (double)iB / (double)(nphi-1);
            double sAB = st_spiral_arclen(r0, g, aB) - st_spiral_arclen(r0, g, aA);
            double uA = (double)R.uv[((size_t)j*(size_t)nphi + (size_t)iA)*2];
            double uB = (double)R.uv[((size_t)j*(size_t)nphi + (size_t)iB)*2];
            double err = fabs(fabs(uB - uA) - sAB) / sAB;   /* sign = gauge */
            fprintf(stderr, "[ribbon selftest]     across-band |du|=%.2f (s_true=%.2f, err=%.1f%%)\n",
                    fabs(uB - uA), sAB, err * 100.0);
            st_check(err < 0.10, "holes band bridged within 10%", &fails);
        }
    }

    /* (4) two wraps 7 vox apart, disconnected: radial prior stacks them. */
    {
        enum { NU = 128, NH = 11 };
        double R1 = 20.0, R2 = 27.0, Hgt = 20.0, cy = 60.0, cx = 60.0;
        size_t nvv = (size_t)NU * NH * 2;
        float *v = (float *)ARENA_ALLOC(arena, (long)(nvv * 3 * sizeof(float)));
        int32_t *f = (int32_t *)ARENA_ALLOC(arena, (long)((size_t)(NU-1)*(NH-1)*2*2*3 * sizeof(int32_t)));
        size_t nfc = 0;
        for (int w = 0; w < 2; w++) {
            double R_ = w ? R2 : R1;
            size_t base = (size_t)w * (size_t)NU * NH;
            for (int j = 0; j < NH; j++) {
                for (int i = 0; i < NU; i++) {
                    double a = 2.0 * M_PI * (double)i / (double)NU;
                    size_t idx = base + (size_t)j * NU + (size_t)i;
                    v[idx*3+0] = (float)(Hgt * (double)j / (double)(NH-1));
                    v[idx*3+1] = (float)(cy + R_ * sin(a));
                    v[idx*3+2] = (float)(cx + R_ * cos(a));
                }
            }
            for (int j = 0; j < NH - 1; j++) {
                for (int i = 0; i < NU - 1; i++) {
                    int32_t a = (int32_t)(base + (size_t)j*NU + (size_t)i);
                    int32_t b = a + 1;
                    int32_t c = (int32_t)(base + (size_t)(j+1)*NU + (size_t)i);
                    int32_t d = c + 1;
                    f[nfc*3+0]=a; f[nfc*3+1]=b; f[nfc*3+2]=c; nfc++;
                    f[nfc*3+0]=b; f[nfc*3+1]=d; f[nfc*3+2]=c; nfc++;
                }
            }
        }
        RibbonOpts o; RibbonOpts_default(&o);
        o.axis_point[1] = (float)cy; o.axis_point[2] = (float)cx;
        RibbonResult R;
        int rc = Ribbon_run(arena, v, nvv, f, nfc, &o, &R);
        fprintf(stderr, "[ribbon selftest] (4) two wraps: rc=%d chains=%d groups=%d "
                "prior_groups=%d cont_pairs=%zu duds_mean=%.4f\n",
                rc, R.n_chains, R.w_groups, R.w_prior_groups, R.n_cont_pairs,
                R.duds_err_mean);
        st_check(rc == 0, "wraps rc", &fails);
        st_check(R.w_groups == 2, "wraps two chain groups", &fails);
        st_check(R.n_cont_pairs == 0, "wraps not glued by continuation", &fails);
        if (rc == 0) {
            double arc1 = R1 * 2.0 * M_PI * (double)(NU-1) / (double)NU;
            double arc2 = R2 * 2.0 * M_PI * (double)(NU-1) / (double)NU;
            double lo1 = 1e300, hi1 = -1e300, lo2 = 1e300, hi2 = -1e300;
            for (size_t i = 0; i < (size_t)NU * NH; i++) {
                double uu = (double)R.uv[i*2];
                if (uu < lo1) lo1 = uu;
                if (uu > hi1) hi1 = uu;
            }
            for (size_t i = (size_t)NU * NH; i < nvv; i++) {
                double uu = (double)R.uv[i*2];
                if (uu < lo2) lo2 = uu;
                if (uu > hi2) hi2 = uu;
            }
            fprintf(stderr, "[ribbon selftest]     wrap1 span=%.2f (arc=%.2f) "
                    "wrap2 span=%.2f (arc=%.2f) sep=%.2f\n",
                    hi1 - lo1, arc1, hi2 - lo2, arc2, lo2 - hi1);
            st_check(fabs((hi1-lo1) - arc1) / arc1 < 0.02, "wrap1 arc length", &fails);
            st_check(fabs((hi2-lo2) - arc2) / arc2 < 0.02, "wrap2 arc length", &fails);
            st_check(lo2 > hi1 - 1.0, "wraps sequential in u (not fused)", &fails);
        }
    }

    /* (4b) THREE disconnected concentric wraps (r=25,35,45, pitch 10). The
     * radial-ordering prior must give them consecutive W so u increases with
     * radius: u(wrap0) < u(wrap1) < u(wrap2), spaced by ~one circumference. */
    {
        enum { NU = 128, NH = 9 };
        double Rw[3] = { 25.0, 35.0, 45.0 }, Hgt = 16.0, cy = 70.0, cx = 70.0;
        size_t nvv = (size_t)NU * NH * 3;
        float *v = (float *)ARENA_ALLOC(arena, (long)(nvv * 3 * sizeof(float)));
        int32_t *f = (int32_t *)ARENA_ALLOC(arena, (long)((size_t)(NU-1)*(NH-1)*3*2*3 * sizeof(int32_t)));
        size_t nfc = 0;
        for (int w = 0; w < 3; w++) {
            double R_ = Rw[w];
            size_t base = (size_t)w * (size_t)NU * NH;
            for (int j = 0; j < NH; j++)
                for (int i = 0; i < NU; i++) {
                    double a = 2.0 * M_PI * (double)i / (double)NU;   /* open seam */
                    size_t idx = base + (size_t)j * NU + (size_t)i;
                    v[idx*3+0] = (float)(Hgt * (double)j / (double)(NH-1));
                    v[idx*3+1] = (float)(cy + R_ * sin(a));
                    v[idx*3+2] = (float)(cx + R_ * cos(a));
                }
            for (int j = 0; j < NH - 1; j++)
                for (int i = 0; i < NU - 1; i++) {
                    int32_t a0 = (int32_t)(base + (size_t)j*NU + (size_t)i);
                    int32_t b0 = a0 + 1;
                    int32_t c0 = (int32_t)(base + (size_t)(j+1)*NU + (size_t)i);
                    int32_t d0 = c0 + 1;
                    f[nfc*3+0]=a0; f[nfc*3+1]=b0; f[nfc*3+2]=c0; nfc++;
                    f[nfc*3+0]=b0; f[nfc*3+1]=d0; f[nfc*3+2]=c0; nfc++;
                }
        }
        RibbonOpts o; RibbonOpts_default(&o);
        o.axis_point[1] = (float)cy; o.axis_point[2] = (float)cx;
        RibbonResult R;
        int rc = Ribbon_run(arena, v, nvv, f, nfc, &o, &R);
        /* per-wrap mean u */
        double mu[3] = {0,0,0}; long mc[3] = {0,0,0};
        if (rc == 0)
            for (int w = 0; w < 3; w++) {
                size_t base = (size_t)w * (size_t)NU * NH;
                for (size_t i = base; i < base + (size_t)NU*NH; i++) { mu[w] += (double)R.uv[i*2]; mc[w]++; }
                mu[w] /= (double)(mc[w] ? mc[w] : 1);
            }
        fprintf(stderr, "[ribbon selftest] (4b) three wraps: rc=%d groups=%d unreached=%d "
                "pitch=%.2f  mean u = %.0f, %.0f, %.0f\n",
                rc, R.w_groups, R.w_unreached, R.pitch_used, mu[0], mu[1], mu[2]);
        st_check(rc == 0, "3wrap rc", &fails);
        st_check(R.w_groups == 3, "3wrap 3 groups", &fails);
        st_check(R.w_unreached == 0, "3wrap all placed", &fails);
        /* u strictly ordered by radius (either ascending or descending is fine;
         * the middle must be between the two extremes and gaps ~one circumf) */
        if (rc == 0) {
            double d01 = mu[1]-mu[0], d12 = mu[2]-mu[1];
            st_check((d01 > 0) == (d12 > 0) && fabs(d01) > 50.0 && fabs(d12) > 50.0,
                     "3wrap u monotone in radius, spaced ~1 wrap", &fails);
        }
    }

    /* (4c) DETACHED SHELL ARCS regression (the "party trick" bug): a spiral
     * plus two disconnected outer arc patches (shell regions that failed to
     * join the main body). Expected: each arc registers onto the main chart's
     * phi->u map as its own small rectangular chart to the RIGHT of the
     * spiral's band (their winding is one wrap beyond), ordered by theta,
     * NOT overlapping the spiral chart mid-band. */
    {
        double r0 = 15.0, g = 2.0, turns = 3.0, Hgt = 40.0, cy = 100.0, cx = 50.0;
        int nphi = 600, nh = 21;
        float *sv; int32_t *sf; size_t snv, snf;
        st_build_spiral(arena, nphi, nh, r0, g, turns, Hgt, cy, cx, NULL, NULL,
                        &sv, &snv, &sf, &snf);
        /* arcs continue the spiral one wrap beyond its end (r = r0+g*(th+6pi)):
         * detached shell regions where the main body wasn't meshed that far */
        enum { AN = 40, AH = 9 };
        double a_th0[2] = { 0.3 * M_PI, 1.2 * M_PI };
        double a_th1[2] = { 0.7 * M_PI, 1.6 * M_PI };
        size_t nvv = snv + (size_t)AN * AH * 2;
        size_t nfc_cap = snf + (size_t)(AN-1) * (AH-1) * 2 * 2;
        float *v = (float *)ARENA_ALLOC(arena, (long)(nvv * 3 * sizeof(float)));
        int32_t *f = (int32_t *)ARENA_ALLOC(arena, (long)(nfc_cap * 3 * sizeof(int32_t)));
        memcpy(v, sv, snv * 3 * sizeof(float));
        memcpy(f, sf, snf * 3 * sizeof(int32_t));
        size_t nfc = snf;
        for (int a = 0; a < 2; a++) {
            size_t base = snv + (size_t)a * AN * AH;
            for (int j = 0; j < AH; j++)
                for (int i = 0; i < AN; i++) {
                    double th = a_th0[a] + (a_th1[a] - a_th0[a]) * (double)i / (double)(AN-1);
                    double rr = r0 + g * (th + turns * 2.0 * M_PI);
                    double z = 8.0 + 24.0 * (double)j / (double)(AH-1);
                    size_t idx = base + (size_t)j * AN + (size_t)i;
                    v[idx*3+0] = (float)z;
                    v[idx*3+1] = (float)(cy + rr * sin(th));
                    v[idx*3+2] = (float)(cx + rr * cos(th));
                }
            for (int j = 0; j < AH - 1; j++)
                for (int i = 0; i < AN - 1; i++) {
                    int32_t a0 = (int32_t)(base + (size_t)j*AN + (size_t)i);
                    int32_t b0 = a0 + 1;
                    int32_t c0 = (int32_t)(base + (size_t)(j+1)*AN + (size_t)i);
                    int32_t d0 = c0 + 1;
                    f[nfc*3+0]=a0; f[nfc*3+1]=b0; f[nfc*3+2]=c0; nfc++;
                    f[nfc*3+0]=b0; f[nfc*3+1]=d0; f[nfc*3+2]=c0; nfc++;
                }
        }
        RibbonOpts o; RibbonOpts_default(&o);
        o.axis_point[1] = (float)cy; o.axis_point[2] = (float)cx;
        RibbonResult R;
        int rc = Ribbon_run(arena, v, nvv, f, nfc, &o, &R);
        double sp_max = -1e300, a_lo[2] = {1e300,1e300}, a_hi[2] = {-1e300,-1e300};
        if (rc == 0) {
            for (size_t i = 0; i < snv; i++)
                if ((double)R.uv[i*2] > sp_max) sp_max = (double)R.uv[i*2];
            for (int a = 0; a < 2; a++) {
                size_t base = snv + (size_t)a * AN * AH;
                for (size_t i = base; i < base + (size_t)AN*AH; i++) {
                    double uu = (double)R.uv[i*2];
                    if (uu < a_lo[a]) a_lo[a] = uu;
                    if (uu > a_hi[a]) a_hi[a] = uu;
                }
            }
        }
        fprintf(stderr, "[ribbon selftest] (4c) shell arcs: rc=%d qp_comps=%d "
                "reg_shift=%.1f  spiral u_max=%.0f  arcA=[%.0f,%.0f] arcB=[%.0f,%.0f]\n"
                "[ribbon selftest]     groups=%d unreached=%d pitch=%.2f pairs=%zu "
                "chains=%d turns=%.2f cover=%.3f cont=%zu bridge=%zu smp=%zu\n",
                rc, R.n_qp_comps, R.reg_max_shift, sp_max,
                a_lo[0], a_hi[0], a_lo[1], a_hi[1],
                R.w_groups, R.w_unreached, R.pitch_used, R.n_pairs,
                R.n_chains, R.phi_span_turns, R.match_cover, R.n_cont_pairs,
                R.bridge_cuts, R.n_samples);
        st_check(rc == 0, "arcs rc", &fails);
        st_check(R.n_qp_comps >= 3, "arcs are separate solve components", &fails);
        if (rc == 0) {
            double lenA = (a_th1[0] - a_th0[0])
                          * (r0 + g * (0.5*(a_th0[0]+a_th1[0]) + turns*2.0*M_PI));
            double lenB = (a_th1[1] - a_th0[1])
                          * (r0 + g * (0.5*(a_th0[1]+a_th1[1]) + turns*2.0*M_PI));
            st_check(a_lo[0] > sp_max - 15.0, "arc A chart right of the spiral band", &fails);
            st_check(a_lo[1] > a_hi[0] - 15.0, "arc B right of arc A (theta order)", &fails);
            st_check(fabs((a_hi[0]-a_lo[0]) - lenA) / lenA < 0.15,
                     "arc A span ~= its arc length", &fails);
            st_check(fabs((a_hi[1]-a_lo[1]) - lenB) / lenB < 0.15,
                     "arc B span ~= its arc length", &fails);
        }
    }

    /* (5) FUSION WALL regression: adjacent spiral turns welded by a radial
     * wall -- the defect that collapses a mesh-integrated winding field.
     * Expect: wall crossings bridge-cut, wraps separated in u. */
    {
        double g = 7.0 / (2.0 * M_PI);           /* pitch 7 vox/turn */
        double r0 = 20.0, turns = 3.0, Hgt = 24.0, cy = 80.0, cx = 80.0;
        int nphi = 600, nh = 13;
        int per_turn = (int)((double)(nphi - 1) / turns);   /* cols per 2pi */
        int ia = per_turn / 2;                   /* theta = pi, turn 1 */
        int ib = ia + per_turn;                  /* theta = pi, turn 2 */
        StTwoCols T = { ia, ib };
        float *v; int32_t *f; size_t nvv, nfc0;
        st_build_spiral(arena, nphi, nh, r0, g, turns, Hgt, cy, cx,
                        st_skip_twocols, &T, &v, &nvv, &f, &nfc0);
        /* wall quads: vert line ia <-> vert line ib (both are boundary now) */
        size_t nfc = nfc0;
        int32_t *f2 = (int32_t *)ARENA_ALLOC(arena, (long)((nfc0 + (size_t)(nh-1)*2) * 3 * sizeof(int32_t)));
        memcpy(f2, f, nfc0 * 3 * sizeof(int32_t));
        for (int j = 0; j < nh - 1; j++) {
            int32_t a = (int32_t)((size_t)j*(size_t)nphi + (size_t)ia);
            int32_t b = (int32_t)((size_t)j*(size_t)nphi + (size_t)ib);
            int32_t c = (int32_t)((size_t)(j+1)*(size_t)nphi + (size_t)ia);
            int32_t d = (int32_t)((size_t)(j+1)*(size_t)nphi + (size_t)ib);
            f2[nfc*3+0]=a; f2[nfc*3+1]=b; f2[nfc*3+2]=c; nfc++;
            f2[nfc*3+0]=b; f2[nfc*3+1]=d; f2[nfc*3+2]=c; nfc++;
        }
        RibbonOpts o; RibbonOpts_default(&o);
        o.axis_point[1] = (float)cy; o.axis_point[2] = (float)cx;
        RibbonResult R;
        int rc = Ribbon_run(arena, v, nvv, f2, nfc, &o, &R);
        fprintf(stderr, "[ribbon selftest] (5) fusion wall: rc=%d bridge_cuts=%zu "
                "groups=%d conf=%zu turns=%.2f duds_mean=%.4f\n",
                rc, R.bridge_cuts, R.w_groups, R.w_conflicts,
                R.phi_span_turns, R.duds_err_mean);
        st_check(rc == 0, "wall rc", &fails);
        st_check(R.bridge_cuts > 0, "wall crossings cut", &fails);
        st_check(R.phi_span_turns > 2.5, "wall winding NOT collapsed", &fails);
        if (rc == 0) {
            /* turn 1 vs turn 2 at theta=0: du ~= one turn's arc, not ~0 */
            int j = nh / 2;
            int i1 = 0, i2 = per_turn;
            double u1_ = (double)R.uv[((size_t)j*(size_t)nphi + (size_t)i1)*2];
            double u2_ = (double)R.uv[((size_t)j*(size_t)nphi + (size_t)i2)*2];
            double arc_turn = 2.0 * M_PI * (r0 + g * M_PI);   /* ~ turn-1 arc */
            fprintf(stderr, "[ribbon selftest]     wrap step |du|=%.1f (one-turn arc ~%.1f)\n",
                    fabs(u2_ - u1_), arc_turn);
            st_check(fabs(u2_ - u1_) > 0.5 * arc_turn, "wraps separated in u", &fails);
        }
    }

    /* (6) PAVA unit. */
    {
        double u[4] = { 0.0, 5.0, 3.0, 4.0 };
        double lb[4] = { 0.0, 1.0, 1.0, 1.0 };
        double scr[16];
        size_t moved = pava_chain(u, lb, 4, scr);
        double exp_[4] = { 0.0, 3.0, 4.0, 5.0 };
        int ok = moved > 0;
        for (int i = 0; i < 4; i++) if (fabs(u[i] - exp_[i]) > 1e-9) ok = 0;
        double u2[3] = { 0.0, 2.0, 4.0 };
        double lb2[3] = { 0.0, 1.0, 1.0 };
        size_t moved2 = pava_chain(u2, lb2, 3, scr);
        ok = ok && moved2 == 0;
        fprintf(stderr, "[ribbon selftest] (6) pava: %s\n", ok ? "ok" : "FAIL");
        st_check(ok, "pava unit", &fails);
    }

    /* (8) coarse -> fine UV map: parameterize a coarse spiral, transfer to a
     * fine sampling of the same surface; fine u must still track arc length. */
    {
        double r0 = 15.0, g = 2.0, turns = 3.0, Hgt = 40.0, cy = 100.0, cx = 50.0;
        float *cv, *fv; int32_t *cf, *ff; size_t cnv, cnf, fnv, fnf;
        st_build_spiral(arena, 300, 11, r0, g, turns, Hgt, cy, cx, NULL, NULL,
                        &cv, &cnv, &cf, &cnf);
        st_build_spiral(arena, 600, 21, r0, g, turns, Hgt, cy, cx, NULL, NULL,
                        &fv, &fnv, &ff, &fnf);
        RibbonOpts o; RibbonOpts_default(&o);
        o.axis_point[1] = (float)cy; o.axis_point[2] = (float)cx;
        RibbonResult R;
        int rc = Ribbon_run(arena, cv, cnv, cf, cnf, &o, &R);
        float *fuv = NULL; uint8_t *fok = NULL; size_t nfb = 0;
        int mrc = -1;
        if (rc == 0)
            mrc = Ribbon_map_uv(arena, cv, cnv, cf, cnf, R.uv,
                                fv, fnv, ff, fnf, &fuv, &fok, &nfb);
        fprintf(stderr, "[ribbon selftest] (8) coarse->fine map: rc=%d mrc=%d "
                "fallback=%zu\n", rc, mrc, nfb);
        st_check(rc == 0 && mrc == 0, "map rc", &fails);
        st_check(nfb == 0, "map no fallbacks", &fails);
        if (rc == 0 && mrc == 0) {
            double amax = turns * 2.0 * M_PI;
            double s_total = st_spiral_arclen(r0, g, amax);
            int nphi = 600, nh = 21, j = nh / 2;
            double u_at0 = (double)fuv[((size_t)j*(size_t)nphi)*2];
            double sdir = (double)fuv[((size_t)j*(size_t)nphi + (size_t)(nphi-1))*2]
                          > u_at0 ? 1.0 : -1.0;
            double e_max = 0.0;
            for (int i = 0; i < nphi; i += 10) {
                double a = amax * (double)i / (double)(nphi - 1);
                double s = st_spiral_arclen(r0, g, a);
                double uu = sdir * ((double)fuv[((size_t)j*(size_t)nphi + (size_t)i)*2] - u_at0);
                double e = fabs(uu - s) / s_total;
                if (e > e_max) e_max = e;
            }
            fprintf(stderr, "[ribbon selftest]     fine |u - s_true|/s_total max = %.4f\n", e_max);
            st_check(e_max < 0.03, "mapped fine u tracks s within 3%", &fails);
        }
    }

    /* (9) bad-link flagging: only a PHYSICALLY LONG bridge with a stretched u is
     * cut; a physically SHORT edge with a huge du (winding-collapse artifact) is
     * KEPT (cutting it would sever real geometry -> floaters), and delamination
     * (du ~ 0) is kept. */
    {
        /* v0,v1,v2 normal (~2 vox, du~len); v3 is a real 10-vox inter-wrap
         * bridge from v2 with huge du (bad link); v4 is 2 vox from v1 but du
         * 1000 (collapsed-core short edge -> must be KEPT); v5 delam overlaps v0. */
        float v[]  = { 0,0,0,  0,0,2,  0,2,0,   0,12,0,   0,0,4,   0.3f,0,0 };
        float uv[] = { 0,0,    2,0,    1,0,      900,0,    1000,0,  0.1f,0 };
        int32_t f[] = { 0,1,2,   2,3,1,   1,4,2,   0,5,2 };  /* normal, LONG bridge, SHORT collapse, delam */
        uint8_t bad[4]; size_t nb = 0;
        int rc = Ribbon_flag_bad_faces(v, 6, f, 4, uv, 4.0, 40.0, 5.0, bad, &nb);
        fprintf(stderr, "[ribbon selftest] (9) bad-link: rc=%d flagged=%zu "
                "(normal=%d longbridge=%d shortcollapse=%d delam=%d)\n",
                rc, nb, bad[0], bad[1], bad[2], bad[3]);
        st_check(rc == 0, "flag rc", &fails);
        st_check(nb == 1 && bad[1] == 1, "the long inter-wrap bridge is cut", &fails);
        st_check(bad[0] == 0 && bad[2] == 0 && bad[3] == 0,
                 "short-collapse edge + normal + delam are KEPT", &fails);
    }

    /* (10) neighbour fill is diagnostic, not direct parameterization support.
     * Vertex 1 lies beyond the transfer radius but is mesh-adjacent to two
     * directly mapped vertices, so it receives a finite filled u. It must stay
     * uv_ok=0 or PlacedCube will retain a constant-u smear through it. */
    {
        Arena_Mark mark = Arena_save(arena);
        float v[9] = {1,0,0,  1,10,0,  1,1,0};
        int32_t f[3] = {0,1,2};
        double t[3] = {1,1,1};
        Sample smp[2];
        memset(smp, 0, sizeof(smp));
        smp[0].p[0] = 1; smp[0].c1 = 0; smp[0].u = 0;
        smp[1].p[0] = 1; smp[1].p[1] = 1;
        smp[1].c1 = 1; smp[1].u = 1;
        for (int i = 0; i < 2; i++) {
            smp[i].tau[1] = 1;
            smp[i].chain = 0;
            smp[i].slice = 0;
        }
        Chain chn;
        memset(&chn, 0, sizeof(chn));
        chn.count = 2;
        SliceSet S;
        memset(&S, 0, sizeof(S));
        S.smp = smp; S.n_smp = 2;
        S.chn = &chn; S.n_chn = 1;
        S.nplanes = 1; S.tmin = 0; S.tmax = 2;
        float prefilled[6] = {0,0, 10,0, 1,0};
        RibbonResult R;
        memset(&R, 0, sizeof(R));
        R.uv = prefilled;
        RibbonOpts o;
        RibbonOpts_default(&o);
        o.emit_global = 1;
        (void)transfer_uv(arena, v, 3, f, 1, t, &S, 2.0, &o, &R);
        fprintf(stderr,
            "[ribbon selftest] (10) transfer provenance: ok=%u,%u,%u "
            "filled=%zu fallback=%zu\n",
            (unsigned)R.uv_ok[0], (unsigned)R.uv_ok[1],
            (unsigned)R.uv_ok[2], R.uv_filled, R.uv_fallback);
        st_check(R.uv_ok[0] == 1 && R.uv_ok[1] == 0 && R.uv_ok[2] == 1,
                 "neighbor-filled vertex stays non-direct", &fails);
        st_check(R.uv_filled == 1 && R.uv_fallback == 0,
                 "neighbor-fill provenance counts", &fails);
        Arena_restore(arena, mark);
    }

    /* (7) degenerate inputs return cleanly. */
    {
        RibbonResult R;
        RibbonOpts o; RibbonOpts_default(&o);
        int rc1 = Ribbon_run(arena, NULL, 0, NULL, 0, &o, &R);
        float tv[9] = { 0,0,0, 0.2f,0,1, 0.2f,1,0 };
        int32_t tf[3] = { 0, 1, 2 };
        int rc2 = Ribbon_run(arena, tv, 3, tf, 1, &o, &R);
        fprintf(stderr, "[ribbon selftest] (7) degenerate: rc_empty=%d rc_tiny=%d\n",
                rc1, rc2);
        st_check(rc1 == -1, "empty -> -1", &fails);
        st_check(rc2 == -1, "tiny -> -1", &fails);
    }

    fprintf(stderr, "[ribbon selftest] %s (%d failure%s)\n",
            fails ? "FAILED" : "PASSED", fails, fails == 1 ? "" : "s");
    Arena_dispose(&arena);
    return fails;
}
