/* ============================================================================
 * ribbon.c -- slice / arc-length / joint parameterization + ribbon fit.
 * See ribbon.h for the pipeline overview. StrokeStrip (Pagurek van Mossel,
 * Liu, Vining, Bessmeltsev, Sheffer 2021) adapted to a scroll sheet:
 *   strokes   = per-slice plane-mesh intersection polylines
 *   strip     = the scroll sheet itself
 *   WTS pairs = cross-slice nearest samples (+ intra-slice continuations)
 *   isolines  = vertical rulings on the papyrus
 *
 * Winding is a TOPOLOGICAL SCAFFOLD, not the material U coordinate.  The
 * preferred path supplies unwrap.c's gauge-synchronized lift: local phase is
 * continuous on every disk and integer chart gauges come from same-sheet
 * continuation plus same-ray layer order.  The legacy slice-only fallback may
 * still derive integer turns from radius, but radius never defines U on the
 * authoritative path.  Stage C obtains U from measured slice arclength and
 * cross-slice alignment; winding fixes only its additive chart gauges.
 * ==========================================================================*/
#define _USE_MATH_DEFINES
#include "ribbon.h"

#include "../common/csr.h"
#include "../common/kdtree.h"
#include "../common/mesh_bin.h"   /* selftest reads back world-frame vmesh */
#include "../common/pca.h"
#include "../common/union_find.h"
#include "../common/ves_platform.h"
#include "../remesh/intersection_cleanup.h"
#include "sparse_solve.h"         /* TAUCS Cholesky for the consistency IRLS */

#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

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
/* Stage C uses the StrokeStrip row weights from atlas_strip.c.  RIB_ALIGN_W
 * remains the legacy pair weight used by winding and the reference-U gauge
 * reconciliation path; it is not the material-coordinate discretization. */
#define RIB_STRIP_LENGTH_W  1.0
#define RIB_STRIP_ALIGN_W   1.0
#define RIB_STRIP_LOCAL_W   0.05   /* AtlasStripOptions_default.lambda_local */
#define RIB_STRIP_CONT_W    1.0
#define RIB_STRIP_INITIAL_LIKE 1e-5
#define RIB_STRIP_L1_EPS    0.25
#define RIB_STRIP_LIKE_FLOOR 1e-7  /* robust weight floor, not a geometry epsilon */
#define RIB_STRIP_GEOM_EPS   1e-6
#define RIB_STRIP_CG_TOL     1e-6
/* Scaffold fits exist only to seed claims/winding; stage-3 metric projection
 * re-solves u exactly, so the scaffold's linear solves stop at coordinate
 * noise, not machine noise. */
#define RIB_STRIP_SCAFFOLD_TOL 1e-4
#define RIB_STRIP_ADMM_RHO    0.05
#define RIB_STRIP_ADMM_TOL    1e-2
#define RIB_STRIP_ADMM_MAXIT  60
#define RIB_BRIDGE_FRAC     0.75   /* bridge cut: |dr| > frac*sample_h per step
                                    * (a spiral ramp moves pitch/(2pi r) ~ 0.01
                                    * vox/vox; a fusion wall ~1) */
#define RIB_XFER_R          3.4    /* mesh-vert transfer search radius (vox);
                                    * < 3.5 = half the 7-vox wrap clearance */
#define RIB_XFER_DPHI_MAX   (M_PI / 2.0) /* authoritative lifted-phase gate.
                                    * Unlike wrapped theta, this distinguishes
                                    * contacting plies separated by ~2pi. */
#define RIB_GAUGE_EDGE_MAX  6.0    /* only short, tangent-consistent source
                                    * edges can synchronize metric gauges */
#define RIB_VJUMP_MAX       3      /* cross-slice matching stride escalation: a
                                    * sample unmatched at k+1 tries k+2, k+3.
                                    * An EMPTY slice (hole band, or a plane
                                    * coinciding exactly with a vertex row)
                                    * would otherwise sever the pair graph
                                    * vertically and shatter the chart. */
#define RIB_CG_TOL          1e-5   /* geometric data are sampled at 2 vox; this
                                    * is already materially tighter than the
                                    * legacy 800-step solves (typically 1e-4) */
#define RIB_CG_MAXIT        800
/* Production ribbon systems lose useful conjugacy well before 50 iterations:
 * on both the 4x and 10x likelihood solves the true residual descends for
 * roughly 10--20 steps after replacement, then drifts upward while the
 * recursive residual remains optimistic.  Replacing at 20 costs one extra
 * matrix application per checkpoint and avoids spending most of each old
 * 50-step window moving away from the verified solution. */
#define RIB_CG_RELIABLE_PERIOD 20
#define RIB_MIN_CHAIN_LEN   0.25   /* drop slivers shorter than this (vox) */
#define RIB_UFILL_MAX       8      /* max invalid u-columns bridged WITHIN a
                                    * ribbon row.  The fill is a convex chord
                                    * between real anchors; longer holes stay
                                    * explicitly absent.
                                    * (16 was tried 2026-08-11: it bridges the
                                    * fake-cover-gap gutters the atlas relies
                                    * on and barely moved fill; v-fill is the
                                    * axis that closes the dash dropouts.) */
#define RIB_WRAP_GATE       6.0    /* max 3D step (vox) between adjacent ribbon
                                    * grid cells. < the 7-vox inter-wrap
                                    * clearance, so a grid edge can never span
                                    * two wraps: where a row or column would
                                    * jump wraps (winding collapse / a hole
                                    * bridging across wraps), the grid SPLITS
                                     * instead of drawing a spike into space. */
#define RIB_MERGE_DU_MAX    8.0    /* carried-U agreement at a chain junction.
                                    * Equals chart_graph_layout.py's
                                    * MAX_FEASIBILITY_U_ERROR: the certified
                                    * post-gauge residual bound of a selected
                                    * relation.  A certified continuation's two
                                    * chain ends agree in U within this; a
                                    * wrong-wrap candidate misses it AND the
                                    * 3-D/radial reach gates simultaneously. */
#define RIB_UVFILL_ROUNDS   10     /* neighbor-fill rounds for unmapped verts */
#define RIB_PAIR_CAND_MAX   16     /* current fitted point plus a bounded set
                                    * of real slice claims for block repair */
#define RIB_PRIOR_MIN_B     1.0    /* spiral prior usable: |b| above this ... */
#define RIB_PRIOR_MIN_R2    0.80   /* ... and fit r2 above this.  A weak spiral
                                    * fit is a diagnostic, not a U model. */
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

/* In-place linear-time median selection.  IRLS only needs one order statistic;
 * sorting tens of millions of residuals every round was needless O(n log n)
 * serial work.  Three-way partitioning is important because clean scroll
 * systems contain very many identical (often zero) residuals. */
static double select_order_dbl(double *a, size_t n, size_t kth)
{
    size_t lo = 0, hi = n;
    if (n == 0) return 0.0;
    if (kth >= n) kth = n - 1;
    while (hi - lo > 1) {
        size_t mid = lo + (hi - lo) / 2;
        double x = a[lo], y = a[mid], z = a[hi - 1];
        double pivot = x < y ? (y < z ? y : (x < z ? z : x))
                             : (x < z ? x : (y < z ? z : y));
        size_t lt = lo, i = lo, gt = hi;
        while (i < gt) {
            if (a[i] < pivot) {
                double t = a[lt]; a[lt++] = a[i]; a[i++] = t;
            } else if (a[i] > pivot) {
                double t = a[--gt]; a[gt] = a[i]; a[i] = t;
            } else {
                i++;
            }
        }
        if (kth < lt) hi = lt;
        else if (kth >= gt) lo = gt;
        else return pivot;
    }
    return a[lo];
}

static double select_median_dbl(double *a, size_t n)
{
    return select_order_dbl(a, n, n / 2);
}

static double select_p95_dbl(double *a, size_t n)
{
    /* Nearest-rank 95th percentile, expressed with integers so large audits
     * do not acquire a platform-dependent floating-point index. */
    size_t rank = n > 0 ? (n / 20) * 19
                          + (((n % 20) * 19 + 19) / 20) : 0;
    size_t kth = rank > 0 ? rank - 1 : 0;
    return select_order_dbl(a, n, kth);
}
static int cmp_i32(const void *pa, const void *pb)
{
    int32_t a = *(const int32_t *)pa, b = *(const int32_t *)pb;
    return a < b ? -1 : (a > b ? 1 : 0);
}

/* Keep allocation arithmetic in size_t all the way into Arena_alloc.  Win64
 * uses the LLP64 model, so a legacy cast through `long` truncates any single
 * array larger than 2 GiB even though both size_t and the arena are 64-bit.
 * Large regions routinely exceed that threshold for slice samples/pairs. */
static size_t rib_array_bytes(size_t count, size_t item_size)
{
    if (item_size != 0 && count > SIZE_MAX / item_size) {
        fprintf(stderr, "ribbon: allocation size overflow (%zu * %zu)\n",
                count, item_size);
        RAISE(Arena_Failed);
        return 0;
    }
    return count * item_size;
}

#define RIB_ALLOC_ARRAY(arena, type, count) \
    ((type *)ARENA_ALLOC((arena), \
                        rib_array_bytes((size_t)(count), sizeof(type))))

/* Deterministic work partition used by the large stable radix passes.  A
 * histogram has 65536 size_t counters per worker, so cap small jobs by a
 * reasonably coarse grain instead of blindly allocating for every hardware
 * thread. */
static int rib_work_threads(size_t n, size_t grain)
{
#ifdef _OPENMP
    /* Ribbon_run is also called from the whole-scroll fleet's outer OpenMP
     * loop.  MSVC serializes nested parallel regions by default, so asking an
     * inner region for N workers and then manually partitioning by that
     * requested N leaves all but band zero untouched.  In particular this
     * corrupted the large-mesh radix sort below and made its edge lookup
     * produce a negative index in Release builds.  The fleet already supplies
     * parallelism at this level; keep per-ribbon work serial while nested. */
    if (omp_in_parallel()) return 1;
    int nt = omp_get_max_threads();
    size_t useful = (n + grain - 1) / grain;
    if (useful < 1) useful = 1;
    if ((size_t)nt > useful) nt = (int)useful;
    return nt > 0 ? nt : 1;
#else
    (void)n; (void)grain;
    return 1;
#endif
}

/* Stable 4x16-bit LSD radix sort.  qsort on the 3*nf mesh-edge keys was the
 * largest serial cost on region-scale launches.  Contiguous input bands and
 * per-thread bucket offsets retain exactly the serial radix order (including
 * the order of equal keys), while all scans/scatters use the available cores. */
static int rib_u64_radix(uint64_t *a, size_t n)
{
    enum { NB = 65536 };
    if (n < 2) return 0;
    uint64_t *b = (uint64_t *)malloc(rib_array_bytes(n, sizeof(uint64_t)));
    int nt = rib_work_threads(n, 262144u);
    size_t hist_n = (size_t)nt * (size_t)NB;
    size_t *hist = (size_t *)calloc(hist_n, sizeof(size_t));
    if (b == NULL || hist == NULL) { free(b); free(hist); return -1; }
    uint64_t *src = a, *dst = b;
    for (int pass = 0; pass < 4; pass++) {
        unsigned shift = (unsigned)pass * 16u;
        memset(hist, 0, hist_n * sizeof(size_t));
#ifdef _OPENMP
#pragma omp parallel num_threads(nt)
#endif
        {
            int tid = 0;
#ifdef _OPENMP
            tid = omp_get_thread_num();
#endif
            size_t lo = n * (size_t)tid / (size_t)nt;
            size_t hi = n * (size_t)(tid + 1) / (size_t)nt;
            size_t *h = hist + (size_t)tid * (size_t)NB;
            for (size_t i = lo; i < hi; i++)
                h[(unsigned)((src[i] >> shift) & UINT64_C(0xffff))]++;
        }
        size_t sum = 0;
        for (size_t k = 0; k < (size_t)NB; k++) {
            for (int tid = 0; tid < nt; tid++) {
                size_t *p = &hist[(size_t)tid * (size_t)NB + k];
                size_t c = *p; *p = sum; sum += c;
            }
        }
#ifdef _OPENMP
#pragma omp parallel num_threads(nt)
#endif
        {
            int tid = 0;
#ifdef _OPENMP
            tid = omp_get_thread_num();
#endif
            size_t lo = n * (size_t)tid / (size_t)nt;
            size_t hi = n * (size_t)(tid + 1) / (size_t)nt;
            size_t *h = hist + (size_t)tid * (size_t)NB;
            for (size_t i = lo; i < hi; i++) {
                unsigned k = (unsigned)((src[i] >> shift) & UINT64_C(0xffff));
                dst[h[k]++] = src[i];
            }
        }
        { uint64_t *tmp = src; src = dst; dst = tmp; }
    }
    if (src != a) memcpy(a, src, n * sizeof(uint64_t));
    free(b); free(hist);
    return 0;
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
    opts->metric_iters  = 0;
    opts->metric_weight = 4.0f;
    opts->solve_threads = 0;
    opts->solve_amg     = 1;
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
    int32_t mesh_comp;/* connected source-mesh component that produced chain */
    int32_t source_chart;/* exact chart provenance for direct graph gating */
    int32_t winding_island;/* continuation-only material candidate; a later
                            * branch split may produce several output components */
    int32_t reconstruction_component;/* branch-aware output identity.  Starts
                            * as winding_island, then separates simultaneous
                            * physically distinct claimant paths before atlas
                            * rasterization. */
    double  atlas_shift;/* additive metric-atlas packing shift; subtract to
                         * recover the pre-pack cover-layout solve gauge */
    int32_t solve_comp;/* connected component of the Stage-C metric system */
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
                      const int32_t *vertex_mesh_comp,
                      const int32_t *mesh_comp_island,
                      const int32_t *mesh_comp_chart,
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
    uint64_t *keys = RIB_ALLOC_ARRAY(arena, uint64_t, nde);
    int edge_threads = rib_work_threads(nde, 262144u);
    int ff;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(edge_threads)
#endif
    for (ff = 0; ff < (int)nf; ff++) {
        size_t f = (size_t)ff;
        for (int m = 0; m < 3; m++) {
            int32_t a = faces[f*3 + (size_t)m], b = faces[f*3 + (size_t)((m+1)%3)];
            uint64_t lo = (uint64_t)(a < b ? a : b), hi = (uint64_t)(a < b ? b : a);
            keys[f*3 + (size_t)m] = lo * (uint64_t)nv + hi;
        }
    }
    uint64_t *ukeys = RIB_ALLOC_ARRAY(arena, uint64_t, nde);
    memcpy(ukeys, keys, nde * sizeof(uint64_t));
    if (rib_u64_radix(ukeys, nde) != 0) return -1;
    size_t ne = 0;
    for (size_t i = 0; i < nde; i++)
        if (i == 0 || ukeys[i] != ukeys[i-1]) ukeys[ne++] = ukeys[i];

    /* per-face edge ids (bsearch into the unique list) */
    int32_t *feid = RIB_ALLOC_ARRAY(arena, int32_t, nde);
    int ii;
    int missing_edge_ids = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(edge_threads) reduction(+:missing_edge_ids)
#endif
    for (ii = 0; ii < (int)nde; ii++) {
        size_t i = (size_t)ii;
        uint64_t *hit = (uint64_t *)bsearch(&keys[i], ukeys, ne,
                                            sizeof(uint64_t), cmp_u64);
        if (hit == NULL) {
            feid[i] = -1;
            missing_edge_ids++;
        } else {
            feid[i] = (int32_t)(hit - ukeys);
        }
    }
    if (missing_edge_ids != 0) {
        fprintf(stderr, "ribbon: edge radix/index mismatch (%d of %zu keys)\n",
                missing_edge_ids, nde);
        return -1;
    }

    /* --- per-edge crossing ranges --- */
    int32_t *ek0 = RIB_ALLOC_ARRAY(arena, int32_t, ne);
    int32_t *ecn = RIB_ALLOC_ARRAY(arena, int32_t, ne);
    int32_t *est = RIB_ALLOC_ARRAY(arena, int32_t, ne + 1);
    unsigned long long total64 = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:total64) num_threads(edge_threads)
#endif
    for (ii = 0; ii < (int)ne; ii++) {
        size_t e = (size_t)ii;
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
        total64 += (unsigned long long)ecn[e];
    }
    size_t total = (size_t)total64;
    est[0] = 0;
    for (size_t e = 0; e < ne; e++) est[e+1] = est[e] + ecn[e];
    if (total == 0) return -1;
    if (total > (size_t)INT32_MAX / 4) return -1;   /* absurd input */

    /* --- fill crossings --- */
    double  *cpos  = RIB_ALLOC_ARRAY(arena, double, total * 3);
    double  *cphi  = o->reference_phi != NULL
                   ? RIB_ALLOC_ARRAY(arena, double, total)
                   : NULL;
    double  *cu    = o->reference_u != NULL
                   ? RIB_ALLOC_ARRAY(arena, double, total)
                   : NULL;
    double  *cshift = o->reference_atlas_shift != NULL
                    ? RIB_ALLOC_ARRAY(arena, double, total)
                    : NULL;
    int32_t *cpln  = RIB_ALLOC_ARRAY(arena, int32_t, total);
    int32_t *cmesh = RIB_ALLOC_ARRAY(arena, int32_t, total);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(edge_threads)
#endif
    for (ii = 0; ii < (int)ne; ii++) {
        size_t e = (size_t)ii;
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
            if (cu != NULL)
                cu[ci] = (1.0 - w) * o->reference_u[a]
                       +        w  * o->reference_u[b];
            if (cshift != NULL)
                cshift[ci] = (1.0 - w) * o->reference_atlas_shift[a]
                           +        w  * o->reference_atlas_shift[b];
            cpln[ci] = k;
            cmesh[ci] = vertex_mesh_comp[a];
        }
    }

    /* --- segments: link the two crossed edges of each straddling face --- */
    int32_t *nbr = RIB_ALLOC_ARRAY(arena, int32_t, total * 2);
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
    double *cth = RIB_ALLOC_ARRAY(arena, double, total);
    int crossing_threads = rib_work_threads(total, 65536u);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(crossing_threads)
#endif
    for (ii = 0; ii < (int)total; ii++) {
        size_t ci = (size_t)ii;
        double d0 = cpos[ci*3+0] - (double)axis_point[0];
        double d1 = cpos[ci*3+1] - (double)axis_point[1];
        double d2 = cpos[ci*3+2] - (double)axis_point[2];
        double a1 = d0*e1[0] + d1*e1[1] + d2*e1[2];
        double a2 = d0*e2[0] + d1*e2[1] + d2*e2[2];
        cth[ci] = atan2(a2, a1);
    }

    /* --- walk chains (open first, then closed loops cut at max |dtheta|) --- */
    uint8_t *vis   = (uint8_t *)ARENA_CALLOC(arena, total, 1);
    int32_t *order = RIB_ALLOC_ARRAY(arena, int32_t, total);
    size_t maxch = total / 2 + 2;
    Chain *chains = RIB_ALLOC_ARRAY(arena, Chain, maxch);
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
                    int32_t *tmp = RIB_ALLOC_ARRAY(arena, int32_t, len);
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
            chains[n_chn].mesh_comp = cmesh[order[start]];
            chains[n_chn].source_chart =
                mesh_comp_chart[chains[n_chn].mesh_comp];
            chains[n_chn].winding_island =
                mesh_comp_island[chains[n_chn].mesh_comp];
            chains[n_chn].reconstruction_component =
                chains[n_chn].winding_island;
            chains[n_chn].atlas_shift = cshift != NULL
                                      ? cshift[order[start]] : 0.0;
            n_chn++;
        }
    }

    /* --- chain lengths + net winding (pass 1) --- */
    double *chlen = RIB_ALLOC_ARRAY(arena, double, n_chn + 1);
    double *chwind = RIB_ALLOC_ARRAY(arena, double, n_chn + 1);
    double sh = (double)o->sample_h;
    if (!(sh > 0.0) || !isfinite(sh)) return -1;
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
        if (tot >= RIB_MIN_CHAIN_LEN) {
            double requested = tot / sh + 0.5;
            size_t add = 0;
            if (!isfinite(requested) || requested > (double)(INT32_MAX - 2))
                return -1;
            add = (size_t)requested + 2;
            if (cap > SIZE_MAX - add) return -1;
            cap += add;
        }
    }
    if (cap < 2) return -1;

    /* --- orient (ascending theta) + resample (pass 2) --- */
    Sample *smp = RIB_ALLOC_ARRAY(arena, Sample, cap);
    Chain  *chn0 = RIB_ALLOC_ARRAY(arena, Chain, n_chn + 1);
    size_t n_smp = 0, n_pre_chn = 0;

    for (size_t c = 0; c < n_chn; c++) {
        int32_t len = chains[c].count;
        int32_t first = chains[c].first;
        double totlen = chlen[c];
        if (totlen < RIB_MIN_CHAIN_LEN) continue;
        int rev = cu != NULL
                ? cu[order[first + len - 1]] < cu[order[first]]
                : chwind[c] < 0.0;
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
                    if (cu != NULL)
                        sp->u = (1.0-w) * cu[order[ia]]
                              +        w  * cu[order[ib]];
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
        chn0[n_pre_chn].mesh_comp = chains[c].mesh_comp;
        chn0[n_pre_chn].source_chart = chains[c].source_chart;
        chn0[n_pre_chn].winding_island = chains[c].winding_island;
        chn0[n_pre_chn].reconstruction_component =
            chains[c].reconstruction_component;
        chn0[n_pre_chn].atlas_shift = chains[c].atlas_shift;
        chn0[n_pre_chn].solve_comp = 0;
        n_pre_chn++;
    }

    /* --- bridge cut fallback ------------------------------------------------
     * A radial wall can indicate a fusion between wraps on a circular spiral,
     * but radial slope is NOT a topological invariant: corners and compressed
     * regions of a real scroll can run almost radially for perfectly valid
     * stretches.  In particular, cutting those runs after we have already been
     * given an authoritative graph-winding field shatters one safe chart into
     * thousands of independently solved strips.  Preserve connectivity when a
     * reference field exists; retain the old geometric heuristic only as the
     * no-reference fallback (or when explicitly forced for diagnostics). */
    int do_bridge_cut = o->radial_bridge_cut > 0 ||
                       (o->radial_bridge_cut == 0 &&
                        o->reference_phi == NULL && o->reference_u == NULL);
    double bridge_dr = RIB_BRIDGE_FRAC * sh;
    Chain *chn = RIB_ALLOC_ARRAY(arena, Chain, n_smp / 2 + 2);
    size_t n_out_chn = 0;
    size_t cuts = 0;
    for (size_t c = 0; c < n_pre_chn; c++) {
        int32_t f = chn0[c].first, cnt = chn0[c].count;
        int32_t segstart = 0;
        for (int32_t i = 1; i <= cnt; i++) {
            int cut = (i == cnt) ||
                      (do_bridge_cut &&
                       fabs(smp[f+i].r - smp[f+i-1].r) > bridge_dr);
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
                chn[n_out_chn].mesh_comp = chn0[c].mesh_comp;
                chn[n_out_chn].source_chart = chn0[c].source_chart;
                chn[n_out_chn].winding_island = chn0[c].winding_island;
                chn[n_out_chn].reconstruction_component =
                    chn0[c].reconstruction_component;
                chn[n_out_chn].atlas_shift = chn0[c].atlas_shift;
                chn[n_out_chn].solve_comp = chn0[c].solve_comp;
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
        if (o->reference_u != NULL
                ? smp[f + cnt - 1].u >= smp[f].u
                : smp[f + cnt - 1].th >= smp[f].th)
            continue;
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
    int32_t *per = (int32_t *)ARENA_CALLOC(arena, (size_t)nplanes,
                                           sizeof(int32_t));
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

/* Pair::k is otherwise an integer winding jump.  Metric projection adds
 * face-topology observations only after winding has been assigned, so this
 * impossible jump value can mark evidence that is part of the input surface
 * itself.  Such a relation is still phase/tangent/island gated, but it must
 * not be dismissed merely because the carried U frame contains the block
 * offset that the relation is meant to repair. */
enum { RIB_PAIR_K_TRUSTED_TOPOLOGY = INT32_MIN };

/* A coherent cube/row block offset can easily be a few dozen voxels, so face
 * topology keeps full influence through this residual.  Beyond it the same
 * observation retains Huber influence (it is never Gaussian-killed), but a
 * single leaked/sliver connection cannot drag a whole gauge component by
 * thousands of voxels. */
static const double RIB_METRIC_TOPOLOGY_FULL_RESIDUAL = 64.0;
static const double RIB_METRIC_CERT_ANCHOR_ODDS = 1023.0;

static double rib_metric_topology_likelihood(double residual)
{
    residual = fabs(residual);
    if (!isfinite(residual)) return 1e-6;
    if (residual <= RIB_METRIC_TOPOLOGY_FULL_RESIDUAL) return 1.0;
    double like = RIB_METRIC_TOPOLOGY_FULL_RESIDUAL / residual;
    return like > 1e-6 ? like : 1e-6;
}

typedef struct {
    Pair  *pairs;  size_t n_pairs;
    size_t n_cross, n_cont;
    size_t matched_a;
    size_t candidates_a;
    size_t relation_island_rejects;
    size_t seed_u_rejects;
} PairSet;

/* Robust geometric continuation evidence, quotiented to slice-chain pairs.
 * The raw matcher can emit millions of sample correspondences; this summary
 * keeps only the summed inlier support per chain pair.  It is a gauge-solve
 * diagnostic: emission continuity comes from the merge partition and the
 * serialized direct chart relations, never from this graph. */
typedef struct {
    uint64_t *key;
    double *support;
    size_t capacity;
    size_t n_edges;
} ChainRelationGraph;

static uint64_t chain_relation_key(int32_t a, int32_t b)
{
    uint32_t lo = (uint32_t)(a < b ? a : b);
    uint32_t hi = (uint32_t)(a < b ? b : a);
    return ((uint64_t)lo << 32) | (uint64_t)hi;
}

static size_t chain_relation_slot(uint64_t key, size_t mask)
{
    key ^= key >> 30;
    key *= UINT64_C(0xbf58476d1ce4e5b9);
    key ^= key >> 27;
    key *= UINT64_C(0x94d049bb133111eb);
    key ^= key >> 31;
    return (size_t)key & mask;
}

static void chain_relation_graph_dispose(ChainRelationGraph *graph)
{
    if (graph == NULL) return;
    free(graph->key);
    free(graph->support);
    memset(graph, 0, sizeof(*graph));
}

static int chain_relation_graph_rehash(ChainRelationGraph *graph,
                                       size_t capacity)
{
    uint64_t *key = (uint64_t *)malloc(capacity * sizeof(*key));
    double *support = (double *)calloc(capacity, sizeof(*support));
    if (key == NULL || support == NULL) {
        free(key);
        free(support);
        return -1;
    }
    for (size_t i = 0; i < capacity; i++) key[i] = UINT64_MAX;
    if (graph->key != NULL) {
        size_t mask = capacity - 1;
        for (size_t i = 0; i < graph->capacity; i++) {
            size_t slot;
            if (graph->key[i] == UINT64_MAX) continue;
            slot = chain_relation_slot(graph->key[i], mask);
            while (key[slot] != UINT64_MAX) slot = (slot + 1) & mask;
            key[slot] = graph->key[i];
            support[slot] = graph->support[i];
        }
    }
    free(graph->key);
    free(graph->support);
    graph->key = key;
    graph->support = support;
    graph->capacity = capacity;
    return 0;
}

static int chain_relation_graph_add(ChainRelationGraph *graph,
                                    int32_t a, int32_t b, double support)
{
    uint64_t key;
    size_t mask, slot;
    if (graph == NULL || a < 0 || b < 0 || a == b || !(support > 0.0))
        return 0;
    if (graph->capacity == 0) {
        if (chain_relation_graph_rehash(graph, 1024) != 0) return -1;
    } else if ((graph->n_edges + 1) * 10 >= graph->capacity * 7) {
        if (graph->capacity > SIZE_MAX / 2 ||
            chain_relation_graph_rehash(graph, graph->capacity * 2) != 0)
            return -1;
    }
    key = chain_relation_key(a, b);
    mask = graph->capacity - 1;
    slot = chain_relation_slot(key, mask);
    while (graph->key[slot] != UINT64_MAX && graph->key[slot] != key)
        slot = (slot + 1) & mask;
    if (graph->key[slot] == UINT64_MAX) {
        graph->key[slot] = key;
        graph->n_edges++;
    }
    graph->support[slot] += support;
    return 0;
}

static int direct_chart_relation_allowed(const RibbonOpts *o,
                                         int32_t a, int32_t b)
{
    size_t lo, hi;
    if (a < 0 || b < 0) return 0;
    if (a == b) return 1;
    if (o == NULL || o->reference_relation_a == NULL ||
        o->reference_relation_b == NULL ||
        o->reference_relation_count == 0)
        return o == NULL || o->reference_chart == NULL;
    if (b < a) { int32_t swap = a; a = b; b = swap; }
    lo = 0; hi = o->reference_relation_count;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        int32_t ea = o->reference_relation_a[mid];
        int32_t eb = o->reference_relation_b[mid];
        if (ea < a || (ea == a && eb < b)) lo = mid + 1;
        else hi = mid;
    }
    return lo < o->reference_relation_count &&
           o->reference_relation_a[lo] == a &&
           o->reference_relation_b[lo] == b;
}

/* counting-sort sample ids by slice -> off[] (nplanes+1); orphans excluded */
static int32_t *sort_by_slice(Arena_T arena, const SliceSet *S, int32_t **out_off)
{
    int np = S->nplanes;
    int32_t *off = (int32_t *)ARENA_CALLOC(arena, (size_t)np + 1,
                                           sizeof(int32_t));
    for (size_t i = 0; i < S->n_smp; i++)
        if (S->smp[i].chain >= 0) off[S->smp[i].slice + 1]++;
    for (int k = 0; k < np; k++) off[k+1] = (int32_t)(off[k+1] + off[k]);
    int32_t *ids = RIB_ALLOC_ARRAY(arena, int32_t, S->n_smp + 1);
    int32_t *cur = RIB_ALLOC_ARRAY(arena, int32_t, (size_t)np);
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
    int32_t *off = (int32_t *)ARENA_CALLOC(arena, (size_t)np + 1,
                                           sizeof(int32_t));
    for (size_t c = 0; c < S->n_chn; c++) off[S->chn[c].slice + 1]++;
    for (int k = 0; k < np; k++) off[k+1] = (int32_t)(off[k+1] + off[k]);
    int32_t *ids = RIB_ALLOC_ARRAY(arena, int32_t, S->n_chn + 1);
    int32_t *cur = RIB_ALLOC_ARRAY(arena, int32_t, (size_t)np);
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
    double seed_du_max = fmax((double)o->match_r,
                              3.0 * (double)o->sample_h);

    Pair *pairs = RIB_ALLOC_ARRAY(arena, Pair, S->n_smp + 16);
    size_t n_pairs = 0;

    /* A dense prefix grid was formerly cleared and scanned once per slice.
     * Its O(nplanes * ncell) work dominates tall, wide scroll volumes even
     * though only a small fraction of cells is occupied by either slice.
     * Keep the dense address table (constant-time lookup), but link only the
     * occupied cells and reset just those cells after each slice pair. */
    int32_t *ghead = RIB_ALLOC_ARRAY(arena, int32_t, ncell);
    int32_t *gnext = RIB_ALLOC_ARRAY(arena, int32_t, S->n_smp + 1);
    size_t *gtouched = RIB_ALLOC_ARRAY(arena, size_t, S->n_smp + 1);
    memset(ghead, 0xff, ncell * sizeof(int32_t));
    uint8_t *matched = (uint8_t *)ARENA_CALLOC(arena, S->n_smp + 1, 1);

    /* stride escalation: a sample unmatched at k+1 tries k+2, then k+3. An
     * EMPTY slice (hole band, or a plane coinciding exactly with a vertex
     * row) must not sever the pair graph vertically -- one missing slice
     * would otherwise shatter a coherent sheet into stacked fragments. */
    for (int stride = 1; stride <= RIB_VJUMP_MAX; stride++) {
        for (int k = 0; k + stride < np; k++) {
            int32_t a0 = off[k],        a1 = off[k+1];
            int32_t b0 = off[k+stride], b1 = off[k+stride+1];
            if (a1 == a0 || b1 == b0) continue;
            size_t ntouched = 0;
            for (int32_t q = b0; q < b1; q++) {
                int32_t ib = ids[q];
                const Sample *sb = &smp[ib];
                int i1 = (int)((sb->c1 - c1min) / cell), i2 = (int)((sb->c2 - c2min) / cell);
                if (i1 < 0) i1 = 0;
                if (i1 > g1-1) i1 = g1-1;
                if (i2 < 0) i2 = 0;
                if (i2 > g2-1) i2 = g2-1;
                size_t cix = (size_t)i1 * (size_t)g2 + (size_t)i2;
                if (ghead[cix] < 0) gtouched[ntouched++] = cix;
                gnext[ib] = ghead[cix];
                ghead[cix] = ib;
            }
            for (int32_t q = a0; q < a1; q++) {
                int32_t ia = ids[q];
                if (stride == 1) P->candidates_a++;
                if (matched[ia]) continue;
                const Sample *sa = &smp[ia];
                int i1 = (int)((sa->c1 - c1min) / cell), i2 = (int)((sa->c2 - c2min) / cell);
                int32_t bestb = -1;
                double bestd2 = r2max, best_seed_du = 1e300;
                for (int d1 = -1; d1 <= 1; d1++) for (int d2i = -1; d2i <= 1; d2i++) {
                    int j1 = i1 + d1, j2 = i2 + d2i;
                    if (j1 < 0 || j2 < 0 || j1 >= g1 || j2 >= g2) continue;
                    size_t cix = (size_t)j1 * (size_t)g2 + (size_t)j2;
                    for (int32_t ib = ghead[cix]; ib >= 0; ib = gnext[ib]) {
                        const Sample *sb = &smp[ib];
                        if (sa->chain < 0 || sb->chain < 0 ||
                            S->chn[sa->chain].winding_island !=
                            S->chn[sb->chain].winding_island) {
                            P->relation_island_rejects++;
                            continue;
                        }
                        if (!direct_chart_relation_allowed(
                                o, S->chn[sa->chain].source_chart,
                                S->chn[sb->chain].source_chart)) {
                            P->relation_island_rejects++;
                            continue;
                        }
                        double dd1 = sb->c1 - sa->c1, dd2 = sb->c2 - sa->c2;
                        double d2 = dd1*dd1 + dd2*dd2;
                        if (d2 >= r2max) continue;
                        /* With a registered winding scaffold, compare the
                         * lifted phase itself.  Comparing wrapped theta here
                         * makes phi and phi+2pi indistinguishable and can join
                         * two contacting plies into one StrokeStrip run. */
                        if (o->reference_phi != NULL) {
                            if (!isfinite(sa->phi) || !isfinite(sb->phi) ||
                                fabs(sb->phi - sa->phi) > RIB_CONT_DPHI_MAX)
                                continue;
                        } else if (fabs(wrap_pi(sb->th - sa->th)) >
                                   RIB_CONT_DPHI_MAX) {
                            continue;
                        }
                        /* A parameterized quadribbon already tells us which
                         * material coordinate an intersection carries.  Keep
                         * the geometric nearest-neighbour construction, but do
                         * not let it jump to a contacting ply whose carried U
                         * is far away.  This is correspondence evidence only;
                         * Stage C still solves a new metric U. */
                        double seed_du = 0.0;
                        if (o->reference_u != NULL && o->solve_reference_u) {
                            seed_du = fabs(sb->u - sa->u);
                            if (!isfinite(sa->u) || !isfinite(sb->u) ||
                                seed_du > seed_du_max) {
                                P->seed_u_rejects++;
                                continue;
                            }
                        }
                        if (v3dot(sa->tau, sb->tau) < cos_gate) continue;
                        if (o->reference_u != NULL && o->solve_reference_u) {
                            if (seed_du > best_seed_du + RIB_STRIP_GEOM_EPS)
                                continue;
                            if (fabs(seed_du - best_seed_du) <=
                                    RIB_STRIP_GEOM_EPS && d2 >= bestd2)
                                continue;
                        } else if (d2 >= bestd2) {
                            continue;
                        }
                        best_seed_du = seed_du;
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
            for (size_t q = 0; q < ntouched; q++)
                ghead[gtouched[q]] = -1;
        }
    }
    P->n_cross = n_pairs;

    /* --- intra-slice fragment continuation pairs --- */
    size_t cont_cap = S->n_chn + 16;
    Pair *cont = RIB_ALLOC_ARRAY(arena, Pair, cont_cap);
    size_t n_cont = 0;
    int32_t *choff = NULL;
    int32_t *chids = chains_by_slice(arena, S, &choff);

    int32_t *chn_next = RIB_ALLOC_ARRAY(arena, int32_t, S->n_chn + 1);
    int cont_cell_radius = (int)ceil(RIB_CONT_GAP_VOX / cell);
    if (cont_cell_radius < 1) cont_cell_radius = 1;
    for (int k = 0; k < np; k++) {
        size_t ntouched = 0;
        for (int32_t y = choff[k]; y < choff[k+1]; y++) {
            int32_t cy = chids[y];
            const Sample *sS = &smp[S->chn[cy].first];
            int i1 = (int)((sS->c1 - c1min) / cell);
            int i2 = (int)((sS->c2 - c2min) / cell);
            if (i1 < 0) i1 = 0;
            if (i1 > g1-1) i1 = g1-1;
            if (i2 < 0) i2 = 0;
            if (i2 > g2-1) i2 = g2-1;
            size_t cix = (size_t)i1 * (size_t)g2 + (size_t)i2;
            if (ghead[cix] < 0) gtouched[ntouched++] = cix;
            chn_next[cy] = ghead[cix];
            ghead[cix] = cy;
        }
        for (int32_t x = choff[k]; x < choff[k+1]; x++) {
            int32_t cx = chids[x];
            const Chain *ce = &S->chn[cx];
            int32_t iE = ce->first + ce->count - 1;
            const Sample *sE = &smp[iE];
            int32_t bestS = -1;
            double bestd2 = RIB_CONT_GAP_VOX * RIB_CONT_GAP_VOX;
            int i1 = (int)((sE->c1 - c1min) / cell);
            int i2 = (int)((sE->c2 - c2min) / cell);
            for (int d1 = -cont_cell_radius; d1 <= cont_cell_radius; d1++) {
                int j1 = i1 + d1;
                if (j1 < 0 || j1 >= g1) continue;
                for (int d2i = -cont_cell_radius;
                     d2i <= cont_cell_radius; d2i++) {
                    int j2 = i2 + d2i;
                    if (j2 < 0 || j2 >= g2) continue;
                    size_t cix = (size_t)j1 * (size_t)g2 + (size_t)j2;
                    for (int32_t cy = ghead[cix]; cy >= 0;
                         cy = chn_next[cy]) {
                        if (cy == cx) continue;
                        const Chain *cs = &S->chn[cy];
                        if (ce->winding_island != cs->winding_island) {
                            P->relation_island_rejects++;
                            continue;
                        }
                        if (!direct_chart_relation_allowed(
                                o, ce->source_chart, cs->source_chart)) {
                            P->relation_island_rejects++;
                            continue;
                        }
                        int32_t iS = cs->first;
                        const Sample *sS = &smp[iS];
                        double dth = o->reference_phi != NULL
                                   ? sS->phi - sE->phi
                                   : wrap_pi(sS->th - sE->th);
                        if (!isfinite(dth) || dth <= 0.0 ||
                            dth > RIB_CONT_DPHI_MAX)
                            continue;
                        if (fabs(sS->r - sE->r) > RIB_CONT_DR_MAX) continue;
                        double dd1 = sS->c1 - sE->c1;
                        double dd2 = sS->c2 - sE->c2;
                        double d2 = dd1*dd1 + dd2*dd2;
                        if (d2 < bestd2 ||
                            (d2 == bestd2 && (bestS < 0 || iS < bestS))) {
                            bestd2 = d2;
                            bestS = iS;
                        }
                    }
                }
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
        for (size_t q = 0; q < ntouched; q++)
            ghead[gtouched[q]] = -1;
    }
    P->n_cont = n_cont;

    Pair *all = RIB_ALLOC_ARRAY(arena, Pair, n_pairs + n_cont + 1);
    memcpy(all, pairs, n_pairs * sizeof(Pair));
    memcpy(all + n_pairs, cont, n_cont * sizeof(Pair));
    P->pairs = all;
    P->n_pairs = n_pairs + n_cont;
    return 0;
}

/* Collapse the final robust sample correspondences to chain-level axial
 * continuation support.  This graph is intentionally evidence only: later
 * branch tracking may cut an edge when two supported descendants coexist at
 * the same material U but are physically different sheets. */
static int chain_relation_graph_from_pairs(const SliceSet *S,
                                           const PairSet *P,
                                           ChainRelationGraph *graph)
{
    memset(graph, 0, sizeof(*graph));
    if (S == NULL || P == NULL) return -1;
    for (size_t p = 0; p < P->n_pairs; p++) {
        const Pair *pair = &P->pairs[p];
        int32_t ca, cb;
        if (pair->like < 0.5 || pair->a < 0 || pair->b < 0 ||
            (size_t)pair->a >= S->n_smp ||
            (size_t)pair->b >= S->n_smp)
            continue;
        ca = S->smp[pair->a].chain;
        cb = S->smp[pair->b].chain;
        if (ca < 0 || cb < 0 || ca == cb ||
            (size_t)ca >= S->n_chn || (size_t)cb >= S->n_chn)
            continue;
        if (chain_relation_graph_add(graph, ca, cb, pair->like) != 0) {
            chain_relation_graph_dispose(graph);
            return -1;
        }
    }
    return 0;
}

/* ============================================================================
 * Stage W -- integer winding scaffold.
 *
 * A scroll is ONE continuous strip and its number of wraps equals the maximal
 * winding number. The previous scheme reconstructed W from cross-slice pair
 * votes and placed vote-disconnected chain groups by a radial-ordering BFS;
 * any group the vote graph could not reach fell to baseW=0 and collapsed onto
 * the anchor, cramming ~30 physical wraps into ~6 u-bands (with empty gaps).
 *
 * When no reference_phi is supplied, the legacy fallback anchors every chain
 * directly by radius. On a spiral r = a + b*psi
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
 *
 * With reference_phi, all of the above is diagnostic/fallback bookkeeping:
 * slice_mesh has already interpolated the registered universal-cover lift and
 * the block below preserves it verbatim.  Stage C then measures material U.
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
        RadSample *rs = RIB_ALLOC_ARRAY(arena, RadSample, nref + 1);
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
        double *gaps = RIB_ALLOC_ARRAY(arena, double, nref + 1);
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
    ChainEdge *ev = RIB_ALLOC_ARRAY(arena, ChainEdge, P->n_pairs + 1);
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
    int32_t *deg = (int32_t *)ARENA_CALLOC(arena, nc + 1, sizeof(int32_t));
    for (size_t i = 0; i < nmaj; i++) { deg[ev[i].lo+1]++; deg[ev[i].hi+1]++; }
    for (size_t c = 0; c < nc; c++) deg[c+1] = (int32_t)(deg[c+1] + deg[c]);
    int32_t *adj = RIB_ALLOC_ARRAY(arena, int32_t, (size_t)deg[nc] + 1);
    int32_t *adk = RIB_ALLOC_ARRAY(arena, int32_t, (size_t)deg[nc] + 1);
    int32_t *cur = RIB_ALLOC_ARRAY(arena, int32_t, nc);
    memcpy(cur, deg, nc * sizeof(int32_t));
    for (size_t i = 0; i < nmaj; i++) {
        adj[cur[ev[i].lo]] = ev[i].hi; adk[cur[ev[i].lo]++] = ev[i].k;
        adj[cur[ev[i].hi]] = ev[i].lo; adk[cur[ev[i].hi]++] = -ev[i].k;
    }
    int32_t *grp = RIB_ALLOC_ARRAY(arena, int32_t, nc);
    for (size_t c = 0; c < nc; c++) grp[c] = -1;
    int32_t *queue = RIB_ALLOC_ARRAY(arena, int32_t, nc);
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
    double *baseacc = (double *)ARENA_CALLOC(
        arena, ngroups ? ngroups : 1, sizeof(double));
    double *basecnt = (double *)ARENA_CALLOC(
        arena, ngroups ? ngroups : 1, sizeof(double));
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

    /* winding-scaffold diagnostics */
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
 * Stage C -- StrokeStrip discretization on maximal V-connected runs.
 * ==========================================================================*/

/* A run member is one slice observation of a latent material ruling.  Its
 * longitudinal derivative follows atlas_strip.c exactly: the forward edge is
 * used except at the end of a stroke, where the backward edge is used. */
typedef struct {
    int32_t sample, run;
    int32_t deriv_lo, deriv_hi;
    double deriv_length;
    double dual_width;       /* longitudinal sample dual width */
    double base_weight;      /* normalized transverse (V) dual width */
    double like;             /* robust StrokeStrip membership */
} RibStripMember;

typedef struct {
    size_t first;
    int32_t count;
    double tangent[3];
    double weight;           /* V support in slice intervals */
    double like_sum;         /* sum base_weight * robust membership */
} RibStripRun;

typedef struct {
    int32_t var, run;
    double value;
} RibStripCoeff;

typedef struct {
    RibStripRun *run;
    RibStripMember *member;
    RibStripCoeff *coeff;    /* two signed derivative coeffs per member */
    size_t nrun, nmember, ncoeff;

    int32_t *sample_member;  /* [n_sample], -1 off the run system */
    int32_t *link_pair;      /* selected cross-slice Pair indices */
    int32_t *pruned_pair;    /* valid but non-mutual cross-slice candidates */
    size_t nlink, npruned;

    /* Transpose of coeff[].  Matrix-free PCG computes one row dot followed by
     * one race-free gather per variable; adjacent U rows may share an edge. */
    size_t *var_coeff_off;
    int32_t *var_coeff;

    /* Objective weights are state, not compile-time constants: the final
     * continuation freezes C(u) and ramps the physical stroke metric without
     * changing any correspondence or introducing per-fragment gauges. */
    double length_weight;
    double align_weight;
    double local_weight;
    double cont_weight;
} RibStripSet;

static double rib_strip_sample_dual(const SliceSet *S, int32_t sample)
{
    int32_t chain = S->smp[sample].chain;
    if (chain < 0 || (size_t)chain >= S->n_chn) return 1.0;
    const Chain *st = &S->chn[chain];
    int32_t ordinal = sample - st->first;
    double before = 0.0, after = 0.0;
    if (ordinal > 0)
        before = S->smp[sample].s - S->smp[sample - 1].s;
    if (ordinal + 1 < st->count)
        after = S->smp[sample + 1].s - S->smp[sample].s;
    double width = 0.5 * (before + after);
    if (before == 0.0) width = 0.5 * after;
    if (after == 0.0) width = 0.5 * before;
    return width > RIB_STRIP_GEOM_EPS ? width : 1.0;
}

static double rib_strip_pair_score(const SliceSet *S, const Pair *p,
                                   int score_by_u)
{
    double d1 = S->smp[p->a].c1 - S->smp[p->b].c1;
    double d2 = S->smp[p->a].c2 - S->smp[p->b].c2;
    double geom = d1 * d1 + d2 * d2;
    if (score_by_u) {
        double du = S->smp[p->a].u - S->smp[p->b].u;
        return du * du + 1e-6 * geom;
    }
    return geom;
}

/* Convert the directed nearest-neighbour candidates into non-branching tracks.
 * A link is retained only when it is the best outgoing link of its lower-V
 * sample and the best incoming link of its upper-V sample.  This is the
 * discrete "connected run" contract: no run can contain two observations from
 * one slice, and a gap or ambiguous branch terminates rather than merging two
 * material coordinates.  No sample or geometry is deleted by this operation. */
/* A slice chain can retain coincident samples at a welded/fill degeneracy.
 * Such a sample has no usable local metric derivative and therefore cannot
 * contribute to a StrokeStrip length row.  Prefer the forward edge (the
 * historical convention), fall back to the backward edge, and fail closed
 * only when both are zero. */
static int rib_strip_sample_derivative(const SliceSet *S, int32_t sample,
                                       int32_t *lo, int32_t *hi,
                                       double *length)
{
    if (sample < 0 || (size_t)sample >= S->n_smp) return 0;
    int32_t chain = S->smp[sample].chain;
    if (chain < 0 || (size_t)chain >= S->n_chn) return 0;
    const Chain *st = &S->chn[chain];
    int32_t ordinal = sample - st->first;
    if (ordinal < 0 || ordinal >= st->count) return 0;
    if (ordinal + 1 < st->count) {
        double d = S->smp[sample + 1].s - S->smp[sample].s;
        if (d > RIB_STRIP_GEOM_EPS) {
            if (lo != NULL) *lo = sample;
            if (hi != NULL) *hi = sample + 1;
            if (length != NULL) *length = d;
            return 1;
        }
    }
    if (ordinal > 0) {
        double d = S->smp[sample].s - S->smp[sample - 1].s;
        if (d > RIB_STRIP_GEOM_EPS) {
            if (lo != NULL) *lo = sample - 1;
            if (hi != NULL) *hi = sample;
            if (length != NULL) *length = d;
            return 1;
        }
    }
    return 0;
}

static int rib_strip_build_runs(Arena_T arena, const SliceSet *S,
                                 PairSet *P, int score_by_u,
                                 double initial_like,
                                 RibStripSet *R)
{
    size_t n = S->n_smp, ncross = P->n_cross;
    memset(R, 0, sizeof(*R));
    R->length_weight = RIB_STRIP_LENGTH_W;
    R->align_weight = RIB_STRIP_ALIGN_W;
    R->local_weight = RIB_STRIP_LOCAL_W;
    R->cont_weight = RIB_STRIP_CONT_W;
    int32_t *out_pair = RIB_ALLOC_ARRAY(arena, int32_t, n + 1);
    int32_t *in_pair = RIB_ALLOC_ARRAY(arena, int32_t, n + 1);
    int32_t *next = RIB_ALLOC_ARRAY(arena, int32_t, n + 1);
    int32_t *prev = RIB_ALLOC_ARRAY(arena, int32_t, n + 1);
    double *out_score = RIB_ALLOC_ARRAY(arena, double, n + 1);
    double *in_score = RIB_ALLOC_ARRAY(arena, double, n + 1);
    uint8_t *selected = (uint8_t *)ARENA_CALLOC(
        arena, ncross ? ncross : 1, sizeof(uint8_t));
    for (size_t i = 0; i < n; i++) {
        out_pair[i] = in_pair[i] = next[i] = prev[i] = -1;
        out_score[i] = in_score[i] = 1e300;
    }

    for (size_t pi = 0; pi < ncross; pi++) {
        Pair *pr = &P->pairs[pi];
        if (pr->like < 0.0 || pr->a < 0 || pr->b < 0 ||
            (size_t)pr->a >= n || (size_t)pr->b >= n ||
            S->smp[pr->a].slice >= S->smp[pr->b].slice)
            continue;
        int32_t ca = S->smp[pr->a].chain;
        int32_t cb = S->smp[pr->b].chain;
        if (ca < 0 || cb < 0 || S->chn[ca].count < 2 ||
            S->chn[cb].count < 2 ||
            !rib_strip_sample_derivative(S, pr->a, NULL, NULL, NULL) ||
            !rib_strip_sample_derivative(S, pr->b, NULL, NULL, NULL))
            continue;
        double score = rib_strip_pair_score(S, pr, score_by_u);
        if (score < out_score[pr->a] - RIB_STRIP_GEOM_EPS ||
            (fabs(score - out_score[pr->a]) <= RIB_STRIP_GEOM_EPS &&
             (out_pair[pr->a] < 0 || (int32_t)pi < out_pair[pr->a]))) {
            out_score[pr->a] = score;
            out_pair[pr->a] = (int32_t)pi;
        }
        if (score < in_score[pr->b] - RIB_STRIP_GEOM_EPS ||
            (fabs(score - in_score[pr->b]) <= RIB_STRIP_GEOM_EPS &&
             (in_pair[pr->b] < 0 || (int32_t)pi < in_pair[pr->b]))) {
            in_score[pr->b] = score;
            in_pair[pr->b] = (int32_t)pi;
        }
    }

    size_t nlink = 0;
    for (size_t pi = 0; pi < ncross; pi++) {
        Pair *pr = &P->pairs[pi];
        if (pr->like >= 0.0 && pr->a >= 0 && pr->b >= 0 &&
            (size_t)pr->a < n && (size_t)pr->b < n &&
            out_pair[pr->a] == (int32_t)pi &&
            in_pair[pr->b] == (int32_t)pi) {
            selected[pi] = 1;
            next[pr->a] = pr->b;
            prev[pr->b] = pr->a;
            nlink++;
        }
    }
    R->link_pair = RIB_ALLOC_ARRAY(arena, int32_t, nlink + 1);
    R->pruned_pair = RIB_ALLOC_ARRAY(arena, int32_t, ncross + 1);
    size_t li = 0;
    for (size_t pi = 0; pi < ncross; pi++) {
        if (selected[pi]) {
            R->link_pair[li++] = (int32_t)pi;
        } else if (P->pairs[pi].like >= 0.0) {
            R->pruned_pair[R->npruned] = (int32_t)pi;
            P->pairs[pi].like = -1.0;
            R->npruned++;
        }
    }
    R->nlink = li;

    size_t nrun = 0, nmember = 0;
    for (size_t i = 0; i < n; i++) {
        if (next[i] < 0 || prev[i] >= 0) continue;
        size_t count = 1;
        int32_t at = (int32_t)i;
        while (next[at] >= 0) { at = next[at]; count++; }
        if (count >= 2) { nrun++; nmember += count; }
    }
    R->run = RIB_ALLOC_ARRAY(arena, RibStripRun, nrun + 1);
    R->member = RIB_ALLOC_ARRAY(arena, RibStripMember, nmember + 1);
    R->sample_member = RIB_ALLOC_ARRAY(arena, int32_t, n + 1);
    for (size_t i = 0; i < n; i++) R->sample_member[i] = -1;

    size_t ri = 0, mi = 0;
    for (size_t start = 0; start < n; start++) {
        if (next[start] < 0 || prev[start] >= 0) continue;
        RibStripRun *run = &R->run[ri];
        memset(run, 0, sizeof(*run));
        run->first = mi;
        int32_t at = (int32_t)start;
        while (1) {
            RibStripMember *m = &R->member[mi];
            memset(m, 0, sizeof(*m));
            m->sample = at;
            m->run = (int32_t)ri;
            m->dual_width = rib_strip_sample_dual(S, at);
            m->like = initial_like;
            if (!rib_strip_sample_derivative(
                    S, at, &m->deriv_lo, &m->deriv_hi,
                    &m->deriv_length))
                return -1; /* pair selection above makes this unreachable */
            R->sample_member[at] = (int32_t)mi;
            mi++;
            run->count++;
            if (next[at] < 0) break;
            at = next[at];
        }

        size_t first = run->first;
        int32_t first_slice = S->smp[R->member[first].sample].slice;
        int32_t last_slice =
            S->smp[R->member[first + (size_t)run->count - 1].sample].slice;
        run->weight = (double)(last_slice - first_slice);
        if (!(run->weight > 0.0)) return -1;
        for (int32_t j = 0; j < run->count; j++) {
            RibStripMember *m = &R->member[first + (size_t)j];
            int32_t slice = S->smp[m->sample].slice;
            double vdual;
            if (j == 0) {
                int32_t sn = S->smp[R->member[first + 1].sample].slice;
                vdual = 0.5 * (double)(sn - slice);
            } else if (j + 1 == run->count) {
                int32_t sp = S->smp[R->member[first + (size_t)j - 1].sample].slice;
                vdual = 0.5 * (double)(slice - sp);
            } else {
                int32_t sp = S->smp[R->member[first + (size_t)j - 1].sample].slice;
                int32_t sn = S->smp[R->member[first + (size_t)j + 1].sample].slice;
                vdual = 0.5 * (double)(sn - sp);
            }
            m->base_weight = vdual / run->weight;
            for (int d = 0; d < 3; d++)
                run->tangent[d] += vdual * S->smp[m->sample].tau[d];
        }
        if (v3norm(run->tangent) < RIB_STRIP_GEOM_EPS)
            memcpy(run->tangent, S->smp[R->member[first].sample].tau,
                   sizeof(run->tangent));
        ri++;
    }
    R->nrun = ri;
    R->nmember = mi;

    R->ncoeff = 2 * R->nmember;
    R->coeff = RIB_ALLOC_ARRAY(arena, RibStripCoeff, R->ncoeff + 1);
    R->var_coeff_off = (size_t *)ARENA_CALLOC(
        arena, n + 1, sizeof(size_t));
    for (size_t m = 0; m < R->nmember; m++) {
        R->coeff[2 * m].var = R->member[m].deriv_lo;
        R->coeff[2 * m].run = R->member[m].run;
        R->coeff[2 * m + 1].var = R->member[m].deriv_hi;
        R->coeff[2 * m + 1].run = R->member[m].run;
        R->var_coeff_off[(size_t)R->member[m].deriv_lo + 1]++;
        R->var_coeff_off[(size_t)R->member[m].deriv_hi + 1]++;
    }
    for (size_t i = 0; i < n; i++)
        R->var_coeff_off[i + 1] += R->var_coeff_off[i];
    R->var_coeff = RIB_ALLOC_ARRAY(arena, int32_t, R->ncoeff + 1);
    size_t *cursor = RIB_ALLOC_ARRAY(arena, size_t, n + 1);
    memcpy(cursor, R->var_coeff_off, n * sizeof(size_t));
    for (size_t k = 0; k < R->ncoeff; k++)
        R->var_coeff[cursor[R->coeff[k].var]++] = (int32_t)k;
    return 0;
}

static void rib_strip_update_length_coefficients(const SliceSet *S,
                                                 RibStripSet *R)
{
    for (size_t r = 0; r < R->nrun; r++) {
        RibStripRun *run = &R->run[r];
        double denom = 0.0;
        run->like_sum = 0.0;
        for (int32_t j = 0; j < run->count; j++) {
            RibStripMember *m = &R->member[run->first + (size_t)j];
            /* StrokeStrip Eq. 6 uses the cross-section dual width d_a.
             * Robust WTS likelihood belongs to E_similar (Eq. 13), not to
             * E_length: allowing a bad correspondence to redefine the metric
             * made the length target move during the local/global loop.  The
             * run's normalized transverse quadrature is base_weight. */
            denom += m->base_weight;
            run->like_sum += m->base_weight * m->like;
        }
        if (denom < 1e-30) denom = 1e-30;
        for (int32_t j = 0; j < run->count; j++) {
            size_t mi = run->first + (size_t)j;
            RibStripMember *m = &R->member[mi];
            double projection = v3dot(S->smp[m->sample].tau, run->tangent);
            double a = m->base_weight / denom *
                       projection / m->deriv_length;
            R->coeff[2 * mi].value = -a;
            R->coeff[2 * mi + 1].value = a;
        }
    }
}

static void rib_strip_audit_reference_seed(Arena_T arena, const SliceSet *S,
                                           const PairSet *P,
                                           const RibStripSet *R)
{
    size_t nedge = 0;
    for (size_t c = 0; c < S->n_chn; c++)
        if (S->chn[c].count > 1)
            nedge += (size_t)S->chn[c].count - 1;

    double *edge_err = RIB_ALLOC_ARRAY(arena, double, nedge + 1);
    double *link_du = RIB_ALLOC_ARRAY(arena, double, R->nlink + 1);
    size_t ne = 0, nl = 0, nonpositive = 0, below_half = 0;
    double edge_mean = 0.0, edge_rms = 0.0;
    for (size_t c = 0; c < S->n_chn; c++) {
        int32_t f = S->chn[c].first, count = S->chn[c].count;
        for (int32_t j = 1; j < count; j++) {
            double ds = S->smp[f + j].s - S->smp[f + j - 1].s;
            double du = S->smp[f + j].u - S->smp[f + j - 1].u;
            if (!(ds > RIB_STRIP_GEOM_EPS) || !isfinite(du)) continue;
            double slope = du / ds;
            double err = fabs(slope - 1.0);
            edge_err[ne++] = err;
            edge_mean += err;
            edge_rms += err * err;
            if (slope <= 0.0) nonpositive++;
            if (slope < 0.5) below_half++;
        }
    }
    for (size_t q = 0; q < R->nlink; q++) {
        const Pair *pr = &P->pairs[R->link_pair[q]];
        double du = fabs(S->smp[pr->b].u - S->smp[pr->a].u);
        if (isfinite(du)) link_du[nl++] = du;
    }
    double edge_p95 = ne ? select_p95_dbl(edge_err, ne) : 0.0;
    double link_median = nl ? select_median_dbl(link_du, nl) : 0.0;
    double link_p95 = nl ? select_p95_dbl(link_du, nl) : 0.0;
    fprintf(stderr,
            "  quadribbon U audit: chain edges=%zu |du/ds-1| "
            "mean=%.4g rms=%.4g p95=%.4g; slope<=0=%zu slope<0.5=%zu\n",
            ne, ne ? edge_mean / (double)ne : 0.0,
            ne ? sqrt(edge_rms / (double)ne) : 0.0, edge_p95,
            nonpositive, below_half);
    fprintf(stderr,
            "  quadribbon U audit: selected V links=%zu |delta U| "
            "median=%.4g p95=%.4g max-gate candidates already excluded\n",
            nl, link_median, link_p95);
}

static double rib_strip_update_likelihoods(RibStripSet *R,
                                           const double *u, int l1_phase,
                                           double sigma)
{
    double max_change = 0.0;
    double inv2sigma2 = 1.0 / (2.0 * sigma * sigma);
    for (size_t r = 0; r < R->nrun; r++) {
        const RibStripRun *run = &R->run[r];
        double q = 0.0, den = 0.0;
        for (int32_t j = 0; j < run->count; j++) {
            const RibStripMember *m = &R->member[run->first + (size_t)j];
            double w = m->base_weight * m->like;
            q += w * u[m->sample];
            den += w;
        }
        q = den > 0.0 ? q / den : 0.0;
        for (int32_t j = 0; j < run->count; j++) {
            RibStripMember *m = &R->member[run->first + (size_t)j];
            double residual = u[m->sample] - q;
            double like;
            if (l1_phase) {
                like = RIB_STRIP_INITIAL_LIKE /
                       sqrt(residual * residual +
                            RIB_STRIP_L1_EPS * RIB_STRIP_L1_EPS);
                if (like > 1.0) like = 1.0;
            } else {
                like = exp(-residual * residual * inv2sigma2);
                if (like < RIB_STRIP_LIKE_FLOOR)
                    like = RIB_STRIP_LIKE_FLOOR;
            }
            double change = fabs(like - m->like);
            if (change > max_change) max_change = change;
            m->like = like;
        }
    }
    return max_change;
}

static void rib_strip_measure(const SliceSet *S, const RibStripSet *R,
                              const double *u, int final_mode,
                              double *length_rms, double *align_rms)
{
    double lr2 = 0.0, ar2 = 0.0;
    size_t ln = 0, an = 0;
    for (size_t r = 0; r < R->nrun; r++) {
        const RibStripRun *run = &R->run[r];
        double derivative = 0.0;
        size_t k0 = 2 * run->first;
        size_t k1 = 2 * (run->first + (size_t)run->count);
        for (size_t k = k0; k < k1; k++)
            derivative += R->coeff[k].value * u[R->coeff[k].var];
        double dr = derivative - 1.0;
        lr2 += dr * dr;
        ln++;

        double q = 0.0, den = 0.0;
        for (int32_t j = 0; j < run->count; j++) {
            const RibStripMember *m = &R->member[run->first + (size_t)j];
            double offset = final_mode
                          ? v3dot(run->tangent, S->smp[m->sample].p) : 0.0;
            double w = m->base_weight * m->like;
            q += w * (u[m->sample] - offset);
            den += w;
        }
        q = den > 0.0 ? q / den : 0.0;
        for (int32_t j = 0; j < run->count; j++) {
            const RibStripMember *m = &R->member[run->first + (size_t)j];
            double offset = final_mode
                          ? v3dot(run->tangent, S->smp[m->sample].p) : 0.0;
            double residual = u[m->sample] - offset - q;
            ar2 += residual * residual;
            an++;
        }
    }
    *length_rms = ln ? sqrt(lr2 / (double)ln) : 0.0;
    *align_rms = an ? sqrt(ar2 / (double)an) : 0.0;
}

/* The averaged StrokeStrip length row can hide equal-and-opposite speed errors
 * on individual slices.  Measure the intrinsic per-edge metric explicitly so
 * robust correspondence reweighting cannot be accepted by sacrificing the
 * arc-length parameterization it is supposed to align. */
static void rib_strip_measure_local_metric(const SliceSet *S, const double *u,
                                           double *mean_abs, double *rms)
{
    double abs_sum = 0.0, square_sum = 0.0;
    size_t count = 0;
    for (size_t c = 0; c < S->n_chn; c++) {
        int32_t first = S->chn[c].first, n = S->chn[c].count;
        for (int32_t j = 1; j < n; j++) {
            double length = S->smp[first + j].s - S->smp[first + j - 1].s;
            if (!(length > RIB_STRIP_GEOM_EPS)) continue;
            double error = (u[first + j] - u[first + j - 1]) / length - 1.0;
            abs_sum += fabs(error);
            square_sum += error * error;
            count++;
        }
    }
    *mean_abs = count ? abs_sum / (double)count : 0.0;
    *rms = count ? sqrt(square_sum / (double)count) : 0.0;
}

/* Vertex-row form of the solve graph.  The old edge-scatter matvec was both
 * serial and inherently racy under OpenMP.  Storing each undirected constraint
 * in both endpoint rows makes every output row independent: V bands distribute
 * naturally across cores, with no atomics and no post-hoc seam arithmetic.
 * `pair_ref == -1` denotes an intrinsic consecutive-chain edge; otherwise it
 * indexes PairSet.pairs, whose IRLS likelihood remains live between rounds. */
typedef struct {
    size_t n, n_adj;
    size_t *off;             /* [n+1] */
    int32_t *nbr;            /* [n_adj] */
    int32_t *pair_ref;       /* [n_adj], -1 = chain edge */
    int32_t *active_pair;    /* permanently gated pairs omitted */
    size_t n_active_pair;
} SolveGraph;

/* Geometric multigrid on the sampled ribbon's logical (U,V) grid.  The first
 * stored level uses constant U blocks because the full sample level is smoothed
 * by exact tridiagonal U-line solves.  Every later aggregate is a pair joined by
 * a real U or V topology edge; cuts and rejected branches therefore remain
 * boundaries on every level.  Semi-coarsening doubles the smaller nominal grid
 * spacing first.  Operators and right-hand sides are restricted by P^T A P and
 * P^T b, then full-multigrid nested iteration solves coarse-to-fine. */
#define RIB_AMG_MAX_LEVELS 24
#define RIB_AMG_COARSE_N   256
#define RIB_AMG_SMOOTH     2
#define RIB_GMG_CHEB_RATIO 30.0
#define RIB_AMG_CHAIN_BLOCK 2
#define RIB_AMG_COARSE_SCALE 1.0
#define RIB_AMG_SAMPLE_AT   8000000u
#define RIB_AMG_PAIR_STRIDE 1
#define RIB_GMG_FMG_CYCLES  1
#define RIB_GMG_FINE_CYCLES 1
#define RIB_AMG_MAX_CHOL_VALUES ((size_t)64000000)
#define RIB_AMG_MAX_CHOL_BLOCK  4096
#define RIB_AMG_STAR_ATTACH_RATIO 0.95
#define RIB_AMG_STAR_WIDE_RATIO   0.98
#define RIB_AMG_STAR_FULL_RATIO   0.995
#define RIB_AMG_STAR_CAP          8
#define RIB_AMG_STAR_WIDE_CAP     32

typedef struct { uint64_t key; double w; } RibWEdge;

typedef struct {
    size_t first;
    int32_t count;
    double weight;
} RibAMGHyperRow;

typedef struct {
    int32_t var;
    int32_t row;
    double value;
} RibAMGHyperCoeff;

typedef struct {
    int n;
    size_t ne;
    int32_t *ea, *eb;       /* unique undirected conductance edges */
    double *ew;
    size_t *off;            /* symmetric row graph */
    int32_t *nbr;
    double *aw;
    double *pin, *diag;
    /* Geometry, kept separately from operator strength.  U edges follow slice
     * chains/physical continuations; V edges are the selected non-branching
     * StrokeStrip links. */
    size_t nu_edge, nv_edge;
    int32_t *ua, *ub, *va, *vb;
    double h_u, h_v;
    double lambda_max;
    char from_axis;           /* 'F' on first stored level, then 'U' or 'V' */
    /* Signed rank-one rows retained through Galerkin restriction.  These are
     * the average-before-square StrokeStrip derivative rows; unlike a graph
     * Laplacian they can have positive off-diagonal entries and therefore
     * cannot be faked with conductance edges. */
    RibAMGHyperRow *hrow;
    RibAMGHyperCoeff *hcoeff;
    size_t nhrow, nhcoeff;
    size_t *var_hcoeff_off;
    int32_t *var_hcoeff;
    int32_t *agg;           /* map to next level, NULL on coarsest */
    int n_coarse;
    int smooth_sweeps;      /* pre/post Chebyshev degree for this transition */
    /* The solve graph can contain hundreds of disconnected scroll charts.
     * Its coarsest matrix is therefore block diagonal.  Factor each connected
     * block independently; one dense n*n factor wastes cubic work coupling
     * blocks that have no edge between them. */
    int n_chol_blocks;
    int32_t *chol_block_off, *chol_node;
    size_t *chol_off, chol_values;
    double *chol;           /* packed dense lower factors, one per block */
} RibAMGLevel;

typedef struct {
    int enabled, nlevel, threads, weak_v;
    int block_size, n_sample_agg, n_fine_agg, pair_stride;
    double sample_h, slice_h;
    int32_t *chain_base;    /* [n_chain+1], chain-local block prolongation */
    RibAMGLevel level[RIB_AMG_MAX_LEVELS];
} RibAMG;

typedef struct {
    double *x[RIB_AMG_MAX_LEVELS];
    double *r[RIB_AMG_MAX_LEVELS];
    double *tmp[RIB_AMG_MAX_LEVELS];
    double *smooth_p[RIB_AMG_MAX_LEVELS];
    double *row_dot[RIB_AMG_MAX_LEVELS];
    double *fmg_rhs[RIB_AMG_MAX_LEVELS];
    double *fmg_x[RIB_AMG_MAX_LEVELS];
    double *line_cp;
    double *fine_tmp, *fine_corr;
    int maxchain;
} RibAMGWork;

/* Optional inspection record for the final robust FMG solve.  A stored level
 * has one aggregate correction per fine sample: `node[i]` names that solved
 * level variable and its correction weight is exactly one.  The checkpoint U
 * is the incoming fine reference plus (coarse solution - aggregate mean of the
 * reference).  This full-approximation prolongation preserves the ribbon's
 * known within-strip arclength instead of collapsing every aggregate to zero
 * width.  Keeping the map as well as the resulting U makes the diagnostic
 * auditable instead of emitting only a picture.  Snapshots are malloc-owned
 * because the per-round arena mark is restored immediately after each solve;
 * Ribbon_run writes and releases them after transfer/atlas shifts. */
typedef struct {
    int stored_level;       /* >=0: RibAMG.level index; -1: full sample grid */
    int nodes;
    double h_u, h_v;
    int32_t *node;          /* [n_sample], -1 for samples outside a live chain */
    float *u;               /* [n_sample], exact prolongated level solution */
} RibGMGSnapshot;

typedef struct {
    size_t n_sample;
    int count;
    RibGMGSnapshot snapshot[RIB_AMG_MAX_LEVELS + 1];
} RibGMGSnapshots;

static int rib_dump_solve_round_checkpoint(
    Arena_T arena, SliceSet *S, const PairSet *P, const RibbonOpts *o,
    const RibbonResult *base, int round, int total_rounds,
    const char *phase);

static void rib_gmg_snapshots_dispose(RibGMGSnapshots *C)
{
    if (C == NULL) return;
    for (int l = 0; l < C->count; l++) {
        free(C->snapshot[l].node);
        free(C->snapshot[l].u);
    }
    memset(C, 0, sizeof(*C));
}

static int rib_gmg_capture_stored_level(const RibAMG *H,
                                        const SliceSet *S,
                                        const RibAMGWork *W,
                                        const double *fine_reference,
                                        int lev, RibGMGSnapshots *C)
{
    if (C == NULL) return 0;
    if (lev < 0 || lev >= H->nlevel ||
        C->count >= RIB_AMG_MAX_LEVELS + 1) return -1;
    size_t n = S->n_smp;
    RibGMGSnapshot *D = &C->snapshot[C->count];
    memset(D, 0, sizeof(*D));
    D->node = (int32_t *)malloc(rib_array_bytes(n, sizeof(*D->node)));
    D->u = (float *)malloc(rib_array_bytes(n, sizeof(*D->u)));
    if (D->node == NULL || D->u == NULL) {
        free(D->node); free(D->u);
        memset(D, 0, sizeof(*D));
        return -1;
    }
    for (size_t i = 0; i < n; i++) {
        int32_t chain = S->smp[i].chain;
        if (chain < 0) {
            D->node[i] = -1;
            D->u[i] = (float)(fine_reference != NULL ? fine_reference[i] : 0.0);
            continue;
        }
        int32_t node = H->chain_base[chain] +
                       ((int32_t)i - S->chn[chain].first) / H->block_size;
        for (int l = 0; l < lev; l++) node = H->level[l].agg[node];
        D->node[i] = node;
    }
    for (size_t i = 0; i < n; i++) {
        int32_t node = D->node[i];
        if (node < 0) continue;
        D->u[i] = (float)((fine_reference != NULL ? fine_reference[i] : 0.0) +
                          W->fmg_x[lev][node]);
    }
    D->stored_level = lev;
    D->nodes = H->level[lev].n;
    D->h_u = H->level[lev].h_u;
    D->h_v = H->level[lev].h_v;
    C->n_sample = n;
    C->count++;
    return 0;
}

static int rib_gmg_capture_full_level(const SliceSet *S, const double *u,
                                      double sample_h, double slice_h,
                                      RibGMGSnapshots *C)
{
    if (C == NULL) return 0;
    if (C->count >= RIB_AMG_MAX_LEVELS + 1) return -1;
    size_t n = S->n_smp;
    RibGMGSnapshot *D = &C->snapshot[C->count];
    memset(D, 0, sizeof(*D));
    D->node = (int32_t *)malloc(rib_array_bytes(n, sizeof(*D->node)));
    D->u = (float *)malloc(rib_array_bytes(n, sizeof(*D->u)));
    if (D->node == NULL || D->u == NULL) {
        free(D->node); free(D->u);
        memset(D, 0, sizeof(*D));
        return -1;
    }
    for (size_t i = 0; i < n; i++) {
        D->node[i] = S->smp[i].chain >= 0 ? (int32_t)i : -1;
        D->u[i] = (float)u[i];
    }
    D->stored_level = -1;
    D->nodes = n <= (size_t)INT_MAX ? (int)n : INT_MAX;
    D->h_u = sample_h;
    D->h_v = slice_h;
    C->n_sample = n;
    C->count++;
    return 0;
}

static int rib_wedge_radix(RibWEdge *a, size_t n)
{
    enum { NB = 65536 };
    if (n < 2) return 0;
    RibWEdge *b = (RibWEdge *)malloc(rib_array_bytes(n, sizeof(RibWEdge)));
    int nt = rib_work_threads(n, 262144u);
    size_t hist_n = (size_t)nt * (size_t)NB;
    size_t *hist = (size_t *)calloc(hist_n, sizeof(size_t));
    if (b == NULL || hist == NULL) { free(b); free(hist); return -1; }
    RibWEdge *src = a, *dst = b;
    for (int pass = 0; pass < 4; pass++) {
        unsigned shift = (unsigned)pass * 16u;
        memset(hist, 0, hist_n * sizeof(size_t));
#ifdef _OPENMP
#pragma omp parallel num_threads(nt)
#endif
        {
            int tid = 0;
#ifdef _OPENMP
            tid = omp_get_thread_num();
#endif
            size_t lo = n * (size_t)tid / (size_t)nt;
            size_t hi = n * (size_t)(tid + 1) / (size_t)nt;
            size_t *h = hist + (size_t)tid * (size_t)NB;
            for (size_t i = lo; i < hi; i++)
                h[(unsigned)((src[i].key >> shift) & UINT64_C(0xffff))]++;
        }
        size_t sum = 0;
        for (size_t k = 0; k < (size_t)NB; k++) {
            for (int tid = 0; tid < nt; tid++) {
                size_t *p = &hist[(size_t)tid * (size_t)NB + k];
                size_t c = *p; *p = sum; sum += c;
            }
        }
#ifdef _OPENMP
#pragma omp parallel num_threads(nt)
#endif
        {
            int tid = 0;
#ifdef _OPENMP
            tid = omp_get_thread_num();
#endif
            size_t lo = n * (size_t)tid / (size_t)nt;
            size_t hi = n * (size_t)(tid + 1) / (size_t)nt;
            size_t *h = hist + (size_t)tid * (size_t)NB;
            for (size_t i = lo; i < hi; i++) {
                unsigned k = (unsigned)((src[i].key >> shift) & UINT64_C(0xffff));
                dst[h[k]++] = src[i];
            }
        }
        { RibWEdge *t = src; src = dst; dst = t; }
    }
    if (src != a) memcpy(a, src, n * sizeof(RibWEdge));
    free(b); free(hist);
    return 0;
}

/* Sort and coalesce a temporary edge list into a persistent level. */
static int rib_amg_set_edges(Arena_T arena, RibAMGLevel *L,
                             RibWEdge *tmp, size_t ntmp)
{
    size_t ne = 0;
    if (rib_wedge_radix(tmp, ntmp) != 0) return -1;
    for (size_t i = 0; i < ntmp; ) {
        size_t j = i + 1;
        double w = tmp[i].w;
        while (j < ntmp && tmp[j].key == tmp[i].key) { w += tmp[j].w; j++; }
        if (w > 0.0) { tmp[ne].key = tmp[i].key; tmp[ne].w = w; ne++; }
        i = j;
    }
    L->ne = ne;
    L->ea = RIB_ALLOC_ARRAY(arena, int32_t, ne + 1);
    L->eb = RIB_ALLOC_ARRAY(arena, int32_t, ne + 1);
    L->ew = RIB_ALLOC_ARRAY(arena, double, ne + 1);
    for (size_t e = 0; e < ne; e++) {
        L->ea[e] = (int32_t)(tmp[e].key >> 32);
        L->eb[e] = (int32_t)(tmp[e].key & UINT64_C(0xffffffff));
        L->ew[e] = tmp[e].w;
    }
    return 0;
}

/* Sort/coalesce an unweighted logical-grid edge set.  These edges decide only
 * which geometric neighbors may be aggregated; their operator coefficients
 * never enter that decision. */
static int rib_amg_set_topology(Arena_T arena,
                                int32_t **out_a, int32_t **out_b,
                                size_t *out_n,
                                RibWEdge *tmp, size_t ntmp)
{
    size_t ne = 0;
    if (rib_wedge_radix(tmp, ntmp) != 0) return -1;
    for (size_t i = 0; i < ntmp; ) {
        size_t j = i + 1;
        while (j < ntmp && tmp[j].key == tmp[i].key) j++;
        tmp[ne++].key = tmp[i].key;
        i = j;
    }
    *out_a = RIB_ALLOC_ARRAY(arena, int32_t, ne + 1);
    *out_b = RIB_ALLOC_ARRAY(arena, int32_t, ne + 1);
    *out_n = ne;
    for (size_t e = 0; e < ne; e++) {
        (*out_a)[e] = (int32_t)(tmp[e].key >> 32);
        (*out_b)[e] = (int32_t)(tmp[e].key & UINT64_C(0xffffffff));
    }
    return 0;
}

static int rib_amg_build_rows(Arena_T arena, RibAMGLevel *L)
{
    int n = L->n;
    L->off = (size_t *)ARENA_CALLOC(arena, (size_t)n + 1, sizeof(size_t));
    L->diag = RIB_ALLOC_ARRAY(arena, double, (size_t)n + 1);
    for (int i = 0; i < n; i++) L->diag[i] = L->pin[i];
    for (size_t e = 0; e < L->ne; e++) {
        int32_t a = L->ea[e], b = L->eb[e];
        if (a < 0 || b < 0 || a >= n || b >= n || a == b) return -1;
        L->off[(size_t)a + 1]++;
        L->off[(size_t)b + 1]++;
        L->diag[a] += L->ew[e];
        L->diag[b] += L->ew[e];
    }
    for (int i = 0; i < n; i++) L->off[i + 1] += L->off[i];
    L->nbr = RIB_ALLOC_ARRAY(arena, int32_t, L->off[n] + 1);
    L->aw = RIB_ALLOC_ARRAY(arena, double, L->off[n] + 1);
    size_t *cur = RIB_ALLOC_ARRAY(arena, size_t, (size_t)n + 1);
    memcpy(cur, L->off, (size_t)n * sizeof(size_t));
    for (size_t e = 0; e < L->ne; e++) {
        int32_t a = L->ea[e], b = L->eb[e];
        size_t pa = cur[a]++, pb = cur[b]++;
        L->nbr[pa] = b; L->aw[pa] = L->ew[e];
        L->nbr[pb] = a; L->aw[pb] = L->ew[e];
    }
    return 0;
}

static int32_t rib_amg_sample_node(const RibAMG *H, const SliceSet *S,
                                   int32_t si);

/* Complete a level after hrow/hcoeff have been assigned.  The diagonal is
 * part of the smoother and already contains the graph/pin contribution; add
 * each signed row's exact diagonal and build its transpose for race-free
 * matrix-vector products. */
static int rib_amg_finalize_hyperrows(Arena_T arena, RibAMGLevel *L)
{
    L->var_hcoeff_off = (size_t *)ARENA_CALLOC(
        arena, (size_t)L->n + 1, sizeof(size_t));
    for (size_t k = 0; k < L->nhcoeff; k++) {
        int32_t v = L->hcoeff[k].var;
        if (v < 0 || v >= L->n || !isfinite(L->hcoeff[k].value)) return -1;
        L->var_hcoeff_off[(size_t)v + 1]++;
    }
    for (int i = 0; i < L->n; i++)
        L->var_hcoeff_off[i + 1] += L->var_hcoeff_off[i];
    L->var_hcoeff = RIB_ALLOC_ARRAY(
        arena, int32_t, L->nhcoeff + 1);
    size_t *cursor = RIB_ALLOC_ARRAY(arena, size_t, (size_t)L->n + 1);
    memcpy(cursor, L->var_hcoeff_off, (size_t)L->n * sizeof(size_t));
    for (size_t r = 0; r < L->nhrow; r++) {
        const RibAMGHyperRow *row = &L->hrow[r];
        if (row->count < 1 || !isfinite(row->weight) ||
            fabs(row->weight) < RIB_STRIP_GEOM_EPS ||
            row->first > L->nhcoeff ||
            (size_t)row->count > L->nhcoeff - row->first)
            return -1;
        for (int32_t j = 0; j < row->count; j++) {
            size_t k = row->first + (size_t)j;
            int32_t v = L->hcoeff[k].var;
            L->hcoeff[k].row = (int32_t)r;
            L->var_hcoeff[cursor[v]++] = (int32_t)k;
            L->diag[v] += row->weight *
                          L->hcoeff[k].value * L->hcoeff[k].value;
        }
    }
    /* A fixed 2/3 Jacobi weight is unsafe for a rank-one row spanning many V
     * samples: D^{-1}aa^T has an eigenvalue proportional to row cardinality.
     * Compute a rigorous Gershgorin upper bound for D^{-1}A and use it to scale
     * the Chebyshev smoother instead. */
    double *row_l1 = RIB_ALLOC_ARRAY(arena, double, L->nhrow + 1);
    double *off_abs = (double *)ARENA_CALLOC(
        arena, (size_t)L->n + 1, sizeof(double));
    for (size_t r = 0; r < L->nhrow; r++) {
        const RibAMGHyperRow *row = &L->hrow[r];
        double s = 0.0;
        for (int32_t j = 0; j < row->count; j++)
            s += fabs(L->hcoeff[row->first + (size_t)j].value);
        row_l1[r] = s;
    }
    for (size_t k = 0; k < L->nhcoeff; k++) {
        const RibAMGHyperCoeff *c = &L->hcoeff[k];
        const RibAMGHyperRow *row = &L->hrow[c->row];
        off_abs[c->var] += fabs(row->weight * c->value) *
                           fmax(0.0, row_l1[c->row] - fabs(c->value));
    }
    double rho = 1.0;
    for (int i = 0; i < L->n; i++) {
        if (!(L->diag[i] > 0.0) || !isfinite(L->diag[i])) return -1;
        for (size_t e = L->off[i]; e < L->off[i + 1]; e++)
            off_abs[i] += fabs(L->aw[e]);
        double bound = 1.0 + off_abs[i] / L->diag[i];
        if (!isfinite(bound)) return -1;
        if (bound > rho) rho = bound;
    }
    L->lambda_max = 1.05 * rho;
    return 0;
}

/* Restrict the exact average-derivative rows to the first chain-block space.
 * A derivative wholly inside one constant block cancels, as Galerkin
 * restriction requires; only block-boundary derivatives survive. */
static int rib_amg_build_initial_hyperrows(Arena_T arena,
                                           const SliceSet *S,
                                           const RibStripSet *R,
                                           const RibAMG *H,
                                           RibAMGLevel *L)
{
    if (R == NULL || R->nrun == 0) {
        L->hrow = RIB_ALLOC_ARRAY(arena, RibAMGHyperRow, 1);
        L->hcoeff = RIB_ALLOC_ARRAY(arena, RibAMGHyperCoeff, 1);
        L->nhrow = L->nhcoeff = 0;
        return rib_amg_finalize_hyperrows(arena, L);
    }
    if (R->nrun > (SIZE_MAX - 1) / 2 ||
        R->nmember > (SIZE_MAX - R->ncoeff - 1)) return -1;
    L->hrow = RIB_ALLOC_ARRAY(arena, RibAMGHyperRow, 2 * R->nrun + 1);
    L->hcoeff = RIB_ALLOC_ARRAY(
        arena, RibAMGHyperCoeff, R->ncoeff + R->nmember + 1);
    size_t nr = 0, nk = 0;
    for (size_t r = 0; r < R->nrun; r++) {
        const RibStripRun *run = &R->run[r];
        /* Exact Galerkin image of the average-before-square derivative row. */
        size_t first = nk;
        for (int32_t j = 0; j < run->count; j++) {
            size_t mi = run->first + (size_t)j;
            const RibStripCoeff *lo = &R->coeff[2 * mi];
            const RibStripCoeff *hi = &R->coeff[2 * mi + 1];
            int32_t a = rib_amg_sample_node(H, S, lo->var);
            int32_t b = rib_amg_sample_node(H, S, hi->var);
            if (a < 0 || b < 0) return -1;
            if (a == b) continue;
            if (lo->value != 0.0) {
                L->hcoeff[nk].var = a;
                L->hcoeff[nk].value = lo->value;
                nk++;
            }
            if (hi->value != 0.0) {
                L->hcoeff[nk].var = b;
                L->hcoeff[nk].value = hi->value;
                nk++;
            }
        }
        if (nk - first >= 2) {
            L->hrow[nr].first = first;
            L->hrow[nr].count = (int32_t)(nk - first);
            L->hrow[nr].weight = R->length_weight * run->weight;
            nr++;
        } else {
            nk = first;
        }

        /* The latent isovalue q_R is eliminated exactly rather than inserted
         * as an algebraic graph node:
         *
         *   W [ sum_i m_i u_i^2 - (sum_i m_i u_i)^2 / M ].
         *
         * The positive diagonal W*m_i was added to L->pin before graph rows
         * were built.  This signed rank-one row is the second term.  Keeping
         * it explicitly makes both A_c=P^T A P and a coarse absolute solve
         * exact, including final-pass tangent offsets in the restricted RHS. */
        if (run->like_sum > 0.0) {
            first = nk;
            for (int32_t j = 0; j < run->count; j++) {
                const RibStripMember *member =
                    &R->member[run->first + (size_t)j];
                int32_t a = rib_amg_sample_node(H, S, member->sample);
                double mw = member->base_weight * member->like;
                if (a < 0) return -1;
                if (mw == 0.0) continue;
                L->hcoeff[nk].var = a;
                L->hcoeff[nk].value = mw;
                nk++;
            }
            if (nk > first) {
                L->hrow[nr].first = first;
                L->hrow[nr].count = (int32_t)(nk - first);
                L->hrow[nr].weight =
                    -R->align_weight * run->weight / run->like_sum;
                nr++;
            } else {
                nk = first;
            }
        }
    }
    L->nhrow = nr;
    L->nhcoeff = nk;
    return rib_amg_finalize_hyperrows(arena, L);
}

/* Galerkin restriction of signed rank-one rows.  Coefficients mapped into one
 * aggregate are summed before the row is emitted; cancellation is essential
 * because the derivative rows annihilate constants. */
static int rib_amg_coarsen_hyperrows(Arena_T arena,
                                     const RibAMGLevel *F,
                                     RibAMGLevel *C)
{
    C->hrow = RIB_ALLOC_ARRAY(arena, RibAMGHyperRow, F->nhrow + 1);
    C->hcoeff = RIB_ALLOC_ARRAY(arena, RibAMGHyperCoeff, F->nhcoeff + 1);
    int32_t *seen = RIB_ALLOC_ARRAY(arena, int32_t, (size_t)C->n + 1);
    int32_t *touched = RIB_ALLOC_ARRAY(arena, int32_t, (size_t)C->n + 1);
    double *value = RIB_ALLOC_ARRAY(arena, double, (size_t)C->n + 1);
    for (int i = 0; i < C->n; i++) seen[i] = -1;
    size_t nr = 0, nk = 0;
    for (size_t r = 0; r < F->nhrow; r++) {
        const RibAMGHyperRow *src = &F->hrow[r];
        int32_t tag = (int32_t)r;
        int32_t ntouched = 0;
        for (int32_t j = 0; j < src->count; j++) {
            const RibAMGHyperCoeff *fc =
                &F->hcoeff[src->first + (size_t)j];
            int32_t a = F->agg[fc->var];
            if (seen[a] != tag) {
                seen[a] = tag;
                value[a] = fc->value;
                touched[ntouched++] = a;
            } else {
                value[a] += fc->value;
            }
        }
        double max_abs = 0.0;
        for (int32_t j = 0; j < ntouched; j++) {
            double av = fabs(value[touched[j]]);
            if (av > max_abs) max_abs = av;
        }
        double zero_tol = RIB_STRIP_GEOM_EPS * max_abs;
        size_t first = nk;
        for (int32_t j = 0; j < ntouched; j++) {
            int32_t a = touched[j];
            if (fabs(value[a]) <= zero_tol) continue;
            C->hcoeff[nk].var = a;
            C->hcoeff[nk].value = value[a];
            nk++;
        }
        if (nk > first) {
            C->hrow[nr].first = first;
            C->hrow[nr].count = (int32_t)(nk - first);
            C->hrow[nr].weight = src->weight;
            nr++;
        } else {
            nk = first;
        }
    }
    C->nhrow = nr;
    C->nhcoeff = nk;
    return rib_amg_finalize_hyperrows(arena, C);
}

static int32_t rib_amg_sample_node(const RibAMG *H, const SliceSet *S,
                                   int32_t si)
{
    int32_t c = (si >= 0 && (size_t)si < S->n_smp) ? S->smp[si].chain : -1;
    if (c < 0) return -1;
    return H->chain_base[c] + (si - S->chn[c].first) / H->block_size;
}

static int rib_amg_build_initial_topology(Arena_T arena,
                                          const SliceSet *S,
                                          const PairSet *P,
                                          const SolveGraph *G,
                                          const RibStripSet *R,
                                          const RibAMG *H,
                                          RibAMGLevel *L)
{
    size_t n_cont = 0, n_cross = 0;
    for (size_t q = 0; q < G->n_active_pair; q++) {
        if ((size_t)G->active_pair[q] >= P->n_cross) n_cont++;
        else if (R == NULL) n_cross++;
    }
    size_t u_cap = (size_t)H->n_sample_agg + n_cont + 1;
    size_t v_cap = (R != NULL ? R->nlink : n_cross) + 1;
    RibWEdge *utmp = (RibWEdge *)malloc(
        rib_array_bytes(u_cap, sizeof(RibWEdge)));
    RibWEdge *vtmp = (RibWEdge *)malloc(
        rib_array_bytes(v_cap, sizeof(RibWEdge)));
    if (utmp == NULL || vtmp == NULL) {
        free(utmp); free(vtmp); return -1;
    }
    size_t nu = 0, nv = 0;
    for (size_t c = 0; c < S->n_chn; c++) {
        int32_t begin = H->chain_base[c], end = H->chain_base[c + 1];
        for (int32_t a = begin; a + 1 < end; a++) {
            utmp[nu].key = ((uint64_t)(uint32_t)a << 32) |
                           (uint64_t)(uint32_t)(a + 1);
            utmp[nu++].w = 1.0;
        }
    }
    for (size_t q = 0; q < G->n_active_pair; q++) {
        int32_t pi = G->active_pair[q];
        int is_cont = (size_t)pi >= P->n_cross;
        if (!is_cont && R != NULL) continue;
        const Pair *pr = &P->pairs[pi];
        int32_t a = rib_amg_sample_node(H, S, pr->a);
        int32_t b = rib_amg_sample_node(H, S, pr->b);
        if (a < 0 || b < 0 || a == b) continue;
        uint32_t lo = (uint32_t)(a < b ? a : b);
        uint32_t hi = (uint32_t)(a < b ? b : a);
        RibWEdge *dst = is_cont ? utmp : vtmp;
        size_t *ndst = is_cont ? &nu : &nv;
        dst[*ndst].key = ((uint64_t)lo << 32) | (uint64_t)hi;
        dst[(*ndst)++].w = 1.0;
    }
    if (R != NULL) {
        for (size_t q = 0; q < R->nlink; q++) {
            const Pair *pr = &P->pairs[R->link_pair[q]];
            int32_t a = rib_amg_sample_node(H, S, pr->a);
            int32_t b = rib_amg_sample_node(H, S, pr->b);
            if (a < 0 || b < 0 || a == b) continue;
            uint32_t lo = (uint32_t)(a < b ? a : b);
            uint32_t hi = (uint32_t)(a < b ? b : a);
            vtmp[nv].key = ((uint64_t)lo << 32) | (uint64_t)hi;
            vtmp[nv++].w = 1.0;
        }
    }
    int rc = rib_amg_set_topology(
        arena, &L->ua, &L->ub, &L->nu_edge, utmp, nu);
    if (rc == 0)
        rc = rib_amg_set_topology(
            arena, &L->va, &L->vb, &L->nv_edge, vtmp, nv);
    free(utmp); free(vtmp);
    if (rc != 0) return -1;
    L->h_u = fmax(H->sample_h * (double)H->block_size,
                  RIB_STRIP_GEOM_EPS);
    L->h_v = fmax(H->slice_h, RIB_STRIP_GEOM_EPS);
    L->from_axis = 'F';
    return 0;
}

static int rib_amg_build_initial(Arena_T arena, const SliceSet *S,
                                 const PairSet *P, const SolveGraph *G,
                                 const RibStripSet *R,
                                 const int32_t *pin_idx, size_t n_pins,
                                 double pin_w, const RibAMG *H,
                                 RibAMGLevel *L)
{
    double stage_t0 = ves_clock_sec();
    size_t cap = G->n_active_pair + (size_t)H->n_fine_agg + 1, m = 0;
    RibWEdge *tmp;
    memset(L, 0, sizeof(*L));
    if (H->n_fine_agg <= 0) return -1;
    L->n = H->n_fine_agg;
    L->pin = (double *)ARENA_CALLOC(
        arena, (size_t)L->n + 1, sizeof(double));
    for (size_t p = 0; p < n_pins; p++) {
        int32_t a = rib_amg_sample_node(H, S, pin_idx[p]);
        if (a >= 0) L->pin[a] += pin_w;
    }
    if (R != NULL) {
        /* Positive diagonal of the analytically eliminated alignment term.
         * Its negative rank-one part is emitted by the hyper-row builder. */
        for (size_t r = 0; r < R->nrun; r++) {
            const RibStripRun *run = &R->run[r];
            double w = R->align_weight * run->weight;
            for (int32_t j = 0; j < run->count; j++) {
                const RibStripMember *member =
                    &R->member[run->first + (size_t)j];
                int32_t a = rib_amg_sample_node(H, S, member->sample);
                if (a < 0) return -1;
                L->pin[a] += w * member->base_weight * member->like;
            }
        }
    }
    fprintf(stderr,
            "    GMG build initial: diagonal %.2fs (nodes=%d runs=%zu members=%zu)\n",
            ves_clock_sec() - stage_t0, L->n,
            R != NULL ? R->nrun : 0, R != NULL ? R->nmember : 0);
    stage_t0 = ves_clock_sec();
    tmp = (RibWEdge *)malloc(rib_array_bytes(cap + 1, sizeof(RibWEdge)));
    if (tmp == NULL) return -1;
    /* Galerkin edges where a true-arclength chain crosses an aggregate block. */
    for (size_t c = 0; c < S->n_chn; c++) {
        int32_t f = S->chn[c].first, cn = S->chn[c].count;
        for (int32_t k = H->block_size; k < cn; k += H->block_size) {
            uint32_t a = (uint32_t)(H->chain_base[c] + k / H->block_size - 1);
            uint32_t b = a + 1;
            double l = S->smp[f + k].s - S->smp[f + k - 1].s;
            double length_eps = R != NULL ? RIB_STRIP_GEOM_EPS : 1e-9;
            if (l < length_eps) l = length_eps;
            tmp[m].key = ((uint64_t)a << 32) | (uint64_t)b;
            tmp[m].w = (R != NULL ? R->local_weight : 1.0) / l;
            m++;
        }
    }
    for (size_t q = 0; q < G->n_active_pair; q++) {
        const Pair *pr = &P->pairs[G->active_pair[q]];
        int is_cont = (size_t)G->active_pair[q] >= P->n_cross;
        if (!is_cont && H->pair_stride > 1 &&
            q % (size_t)H->pair_stride != 0) continue;
        int32_t ca = rib_amg_sample_node(H, S, pr->a);
        int32_t cb = rib_amg_sample_node(H, S, pr->b);
        if (ca < 0 || cb < 0 || ca == cb || pr->like <= 0.0) continue;
        uint32_t lo = (uint32_t)(ca < cb ? ca : cb);
        uint32_t hi = (uint32_t)(ca < cb ? cb : ca);
        tmp[m].key = ((uint64_t)lo << 32) | (uint64_t)hi;
        /* Preserve the expected pair/chain conductance ratio when a large
         * fine graph is represented by a deterministic stride sample.  This
         * hierarchy is only a preconditioner; the PCG operator and RHS still
         * contain every constraint exactly.  Continuation edges are sparse
         * and often unique, so they are always retained at their true weight. */
        tmp[m].w = (R != NULL ? R->cont_weight : pr->w) * pr->like *
                   (double)(is_cont ? 1 : H->pair_stride);
        m++;
    }
    if (rib_amg_set_edges(arena, L, tmp, m) != 0) { free(tmp); return -1; }
    free(tmp);
    fprintf(stderr,
            "    GMG build initial: operator edges %.2fs (raw=%zu unique=%zu)\n",
            ves_clock_sec() - stage_t0, m, L->ne);
    stage_t0 = ves_clock_sec();
    if (rib_amg_build_initial_topology(arena, S, P, G, R, H, L) != 0)
        return -1;
    fprintf(stderr,
            "    GMG build initial: topology %.2fs (U=%zu V=%zu)\n",
            ves_clock_sec() - stage_t0, L->nu_edge, L->nv_edge);
    stage_t0 = ves_clock_sec();
    if (rib_amg_build_rows(arena, L) != 0) return -1;
    fprintf(stderr, "    GMG build initial: graph rows %.2fs\n",
            ves_clock_sec() - stage_t0);
    stage_t0 = ves_clock_sec();
    {
        int rc = rib_amg_build_initial_hyperrows(arena, S, R, H, L);
        fprintf(stderr,
                "    GMG build initial: signed rows %.2fs (rows=%zu coeff=%zu)\n",
                ves_clock_sec() - stage_t0, L->nhrow, L->nhcoeff);
        return rc;
    }
}

static int rib_amg_coarsen_topology(Arena_T arena,
                                    const int32_t *fa, const int32_t *fb,
                                    size_t nf, const int32_t *agg,
                                    int32_t **ca, int32_t **cb, size_t *nc)
{
    RibWEdge *tmp = (RibWEdge *)malloc(
        rib_array_bytes(nf + 1, sizeof(RibWEdge)));
    if (tmp == NULL) return -1;
    size_t m = 0;
    for (size_t e = 0; e < nf; e++) {
        int32_t a = agg[fa[e]], b = agg[fb[e]];
        if (a == b) continue;
        uint32_t lo = (uint32_t)(a < b ? a : b);
        uint32_t hi = (uint32_t)(a < b ? b : a);
        tmp[m].key = ((uint64_t)lo << 32) | (uint64_t)hi;
        tmp[m++].w = 1.0;
    }
    int rc = rib_amg_set_topology(arena, ca, cb, nc, tmp, m);
    free(tmp);
    return rc;
}

/* Deterministic unsmoothed aggregation on one logical ribbon axis.  A maximal
 * matching establishes ordinary 2x cells.  Preserve that regular interpolation
 * while it removes at least 5% of the level.  Only after contraction turns the
 * axis graph into stars (pair-only coarse/fine ratio >= 0.8) are unmatched
 * non-isolated nodes attached to the first adjacent matched cell.  Attachment
 * ramps with the stalled pair ratio: cap 8 at 95%, cap 32 at 98%, and uncapped
 * only at 99.5%.  This clears artificial crumbs over several resolvable levels
 * instead of destroying their low-frequency modes in one 20x jump.  Every
 * aggregate remains connected by real edges on the requested axis. */
static int rib_amg_axis_aggregate(const int32_t *ga, const int32_t *gb,
                                  size_t ng, int n, int32_t *agg)
{
    int nc = 0;
    for (int i = 0; i < n; i++) agg[i] = -1;
    for (size_t e = 0; e < ng; e++) {
        int32_t a = ga[e], b = gb[e];
        if (a < 0 || b < 0 || a >= n || b >= n || a == b) return -1;
        if (agg[a] >= 0 || agg[b] >= 0) continue;
        agg[a] = agg[b] = nc++;
    }
    {
        int pair_only_coarse = n - nc;
        double pair_ratio = (double)pair_only_coarse / (double)n;
        if (pair_ratio >= RIB_AMG_STAR_ATTACH_RATIO) {
            int max_agg = pair_ratio >= RIB_AMG_STAR_FULL_RATIO ? INT_MAX
                        : pair_ratio >= RIB_AMG_STAR_WIDE_RATIO
                            ? RIB_AMG_STAR_WIDE_CAP
                            : RIB_AMG_STAR_CAP;
            int32_t *agg_size = (int32_t *)malloc(
                rib_array_bytes((size_t)nc + 1, sizeof(int32_t)));
            if (agg_size == NULL) return -1;
            for (int i = 0; i < nc; i++) agg_size[i] = 2;
            for (size_t e = 0; e < ng; e++) {
                int32_t a = ga[e], b = gb[e];
                if (agg[a] < 0 && agg[b] >= 0 &&
                    agg_size[agg[b]] < max_agg) {
                    agg[a] = agg[b];
                    agg_size[agg[b]]++;
                } else if (agg[b] < 0 && agg[a] >= 0 &&
                           agg_size[agg[a]] < max_agg) {
                    agg[b] = agg[a];
                    agg_size[agg[a]]++;
                }
            }
            free(agg_size);
            fprintf(stderr,
                    "    GMG geometric crumb aggregation: pair_ratio=%.4f cap=%s\n",
                    pair_ratio, max_agg == INT_MAX ? "full"
                                                   : max_agg == RIB_AMG_STAR_WIDE_CAP
                                                       ? "32" : "8");
        }
    }
    for (int i = 0; i < n; i++)
        if (agg[i] < 0) agg[i] = nc++;
    return nc;
}

/* One geometric semi-coarsening step.  Matrix magnitude is deliberately
 * absent: a weak but real ribbon neighbor may aggregate, while a strong
 * spatial contact across a cut can never aggregate. */
static int rib_amg_build_next(Arena_T arena, RibAMGLevel *F,
                               RibAMGLevel *C, char axis)
{
    int n = F->n;
    const int32_t *ga = axis == 'U' ? F->ua : F->va;
    const int32_t *gb = axis == 'U' ? F->ub : F->vb;
    size_t ng = axis == 'U' ? F->nu_edge : F->nv_edge;
    if (ng == 0) return 1;
    int32_t *agg = RIB_ALLOC_ARRAY(arena, int32_t, (size_t)n + 1);
    int nc = rib_amg_axis_aggregate(ga, gb, ng, n, agg);
    if (nc < 0) return -1;
    if (nc >= n) return 1;
    F->agg = agg;
    F->n_coarse = nc;

    memset(C, 0, sizeof(*C));
    C->n = nc;
    C->pin = (double *)ARENA_CALLOC(arena, (size_t)nc + 1, sizeof(double));
    for (int i = 0; i < n; i++) C->pin[agg[i]] += F->pin[i];

    RibWEdge *tmp = (RibWEdge *)malloc(
        rib_array_bytes(F->ne + 1, sizeof(RibWEdge)));
    if (tmp == NULL) return -1;
    size_t m = 0;
    for (size_t e = 0; e < F->ne; e++) {
        int32_t a = agg[F->ea[e]], b = agg[F->eb[e]];
        if (a == b) continue;
        uint32_t lo = (uint32_t)(a < b ? a : b);
        uint32_t hi = (uint32_t)(a < b ? b : a);
        tmp[m].key = ((uint64_t)lo << 32) | (uint64_t)hi;
        tmp[m].w = F->ew[e];
        m++;
    }
    if (rib_amg_set_edges(arena, C, tmp, m) != 0) {
        free(tmp); return -1;
    }
    free(tmp);
    if (rib_amg_coarsen_topology(arena, F->ua, F->ub, F->nu_edge,
                                 agg, &C->ua, &C->ub, &C->nu_edge) != 0)
        return -1;
    if (rib_amg_coarsen_topology(arena, F->va, F->vb, F->nv_edge,
                                 agg, &C->va, &C->vb, &C->nv_edge) != 0)
        return -1;
    C->h_u = F->h_u * (axis == 'U' ? 2.0 : 1.0);
    C->h_v = F->h_v * (axis == 'V' ? 2.0 : 1.0);
    C->from_axis = axis;
    if (rib_amg_build_rows(arena, C) != 0) return -1;
    if (rib_amg_coarsen_hyperrows(arena, F, C) != 0) return -1;
    {
        double ratio = (double)F->n / (double)C->n;
        int sweeps = 2;
        if (ratio > 4.0) sweeps = (int)ceil(2.0 * sqrt(ratio));
        if (sweeps > 16) sweeps = 16;
        F->smooth_sweeps = sweeps;
    }
    return 0;
}

static int rib_amg_factor_coarsest(Arena_T arena, RibAMGLevel *L)
{
    int n = L->n;
    if (n <= 0) return -1;
    UnionFind uf = UF_new(arena, (int32_t)n);
    for (size_t e = 0; e < L->ne; e++)
        uf_union(&uf, L->ea[e], L->eb[e]);
    for (size_t r = 0; r < L->nhrow; r++) {
        const RibAMGHyperRow *row = &L->hrow[r];
        int32_t root = L->hcoeff[row->first].var;
        for (int32_t j = 1; j < row->count; j++)
            uf_union(&uf, root,
                     L->hcoeff[row->first + (size_t)j].var);
    }

    int32_t *root_block = RIB_ALLOC_ARRAY(arena, int32_t, (size_t)n);
    int32_t *block_of = RIB_ALLOC_ARRAY(arena, int32_t, (size_t)n);
    int32_t *count = (int32_t *)ARENA_CALLOC(
        arena, (size_t)n + 1, sizeof(int32_t));
    for (int i = 0; i < n; i++) root_block[i] = -1;
    int nb = 0;
    for (int i = 0; i < n; i++) {
        int32_t root = uf_find(&uf, (int32_t)i);
        if (root_block[root] < 0) root_block[root] = nb++;
        block_of[i] = root_block[root];
        count[block_of[i]]++;
    }
    L->n_chol_blocks = nb;
    L->chol_block_off = RIB_ALLOC_ARRAY(arena, int32_t, (size_t)nb + 1);
    L->chol_node = RIB_ALLOC_ARRAY(arena, int32_t, (size_t)n);
    L->chol_off = RIB_ALLOC_ARRAY(arena, size_t, (size_t)nb + 1);
    L->chol_block_off[0] = 0;
    for (int b = 0; b < nb; b++)
        L->chol_block_off[b + 1] = L->chol_block_off[b] + count[b];
    int32_t *cur = RIB_ALLOC_ARRAY(arena, int32_t, (size_t)nb);
    memcpy(cur, L->chol_block_off, (size_t)nb * sizeof(int32_t));
    for (int i = 0; i < n; i++)
        L->chol_node[cur[block_of[i]]++] = (int32_t)i;

    int32_t *local = RIB_ALLOC_ARRAY(arena, int32_t, (size_t)n);
    L->chol_off[0] = 0;
    size_t max_block = 0;
    for (int b = 0; b < nb; b++) {
        size_t m = (size_t)(L->chol_block_off[b + 1] - L->chol_block_off[b]);
        if (m > max_block) max_block = m;
        if (m != 0 && m > (SIZE_MAX - L->chol_off[b]) / m) return -1;
        L->chol_off[b + 1] = L->chol_off[b] + m * m;
        for (int32_t q = L->chol_block_off[b];
             q < L->chol_block_off[b + 1]; q++)
            local[L->chol_node[q]] = q - L->chol_block_off[b];
    }
    L->chol_values = L->chol_off[nb];
    fprintf(stderr,
            "    GMG coarse preflight: blocks=%d max_block=%zu chol_values=%zu (%.3f GB)\n",
            nb, max_block, L->chol_values,
            (double)L->chol_values * (double)sizeof(double) / 1e9);
    fflush(stderr);
    if (max_block > RIB_AMG_MAX_CHOL_BLOCK ||
        L->chol_values > RIB_AMG_MAX_CHOL_VALUES) {
        fprintf(stderr,
                "    GMG coarse factor refused: limit max_block=%d chol_values=%zu\n",
                RIB_AMG_MAX_CHOL_BLOCK, RIB_AMG_MAX_CHOL_VALUES);
        return -1;
    }
    L->chol = (double *)ARENA_CALLOC(
        arena, L->chol_values ? L->chol_values : 1, sizeof(double));
    if (L->chol == NULL) return -1;

    for (int b = 0; b < nb; b++) {
        int32_t bo = L->chol_block_off[b];
        int m = L->chol_block_off[b + 1] - bo;
        double *A = L->chol + L->chol_off[b];
        for (int i = 0; i < m; i++) {
            int node = L->chol_node[bo + i];
            A[(size_t)i * (size_t)m + (size_t)i] = L->diag[node];
        }
    }
    for (size_t e = 0; e < L->ne; e++) {
        int a = L->ea[e], bnode = L->eb[e], bid = block_of[a];
        int m = L->chol_block_off[bid + 1] - L->chol_block_off[bid];
        int ia = local[a], ib = local[bnode];
        double *A = L->chol + L->chol_off[bid];
        A[(size_t)ia * (size_t)m + (size_t)ib] -= L->ew[e];
        A[(size_t)ib * (size_t)m + (size_t)ia] -= L->ew[e];
    }
    /* L->diag already contains each hyper-row diagonal.  Add only its signed
     * off-diagonal outer-product entries here. */
    for (size_t r = 0; r < L->nhrow; r++) {
        const RibAMGHyperRow *row = &L->hrow[r];
        int32_t first_node = L->hcoeff[row->first].var;
        int bid = block_of[first_node];
        int m = L->chol_block_off[bid + 1] - L->chol_block_off[bid];
        double *A = L->chol + L->chol_off[bid];
        for (int32_t a = 0; a < row->count; a++) {
            const RibAMGHyperCoeff *ca =
                &L->hcoeff[row->first + (size_t)a];
            int ia = local[ca->var];
            for (int32_t b = 0; b < a; b++) {
                const RibAMGHyperCoeff *cb =
                    &L->hcoeff[row->first + (size_t)b];
                int ib = local[cb->var];
                double v = row->weight * ca->value * cb->value;
                A[(size_t)ia * (size_t)m + (size_t)ib] += v;
                A[(size_t)ib * (size_t)m + (size_t)ia] += v;
            }
        }
    }
    for (int bid = 0; bid < nb; bid++) {
        int32_t bo = L->chol_block_off[bid];
        int m = L->chol_block_off[bid + 1] - bo;
        double *A = L->chol + L->chol_off[bid];
        for (int i = 0; i < m; i++) {
            for (int j = 0; j <= i; j++) {
                double s = A[(size_t)i * (size_t)m + (size_t)j];
                for (int k = 0; k < j; k++)
                    s -= A[(size_t)i * (size_t)m + (size_t)k] *
                         A[(size_t)j * (size_t)m + (size_t)k];
                if (i == j) {
                    int node = L->chol_node[bo + i];
                    double floorv = 1e-12 * (1.0 + fabs(L->diag[node]));
                    if (s < floorv) s = floorv;
                    A[(size_t)i * (size_t)m + (size_t)j] = sqrt(s);
                } else {
                    A[(size_t)i * (size_t)m + (size_t)j] =
                        s / A[(size_t)j * (size_t)m + (size_t)j];
                }
            }
            for (int j = i + 1; j < m; j++)
                A[(size_t)i * (size_t)m + (size_t)j] = 0.0;
        }
    }
    return 0;
}

static int rib_amg_build(Arena_T arena, const SliceSet *S, const PairSet *P,
                          const SolveGraph *G, const RibStripSet *R,
                          const int32_t *pin_idx,
                         size_t n_pins, double pin_w,
                         double sample_h, double slice_h, int threads,
                         RibAMG *H)
{
    memset(H, 0, sizeof(*H));
    H->threads = threads;
    H->block_size = RIB_AMG_CHAIN_BLOCK;
    H->sample_h = fmax(sample_h, RIB_STRIP_GEOM_EPS);
    H->slice_h = fmax(slice_h, RIB_STRIP_GEOM_EPS);
    if (R != NULL && R->nrun > 0) {
        double mean_membership = 0.0;
        for (size_t r = 0; r < R->nrun; r++)
            mean_membership += R->run[r].like_sum;
        mean_membership /= (double)R->nrun;
        H->weak_v = mean_membership < 1e-3;
    }
    H->pair_stride = G->n_active_pair + (R != NULL ? R->nlink : 0) >=
                     RIB_AMG_SAMPLE_AT
                   ? RIB_AMG_PAIR_STRIDE : 1;
    H->chain_base = RIB_ALLOC_ARRAY(arena, int32_t, S->n_chn + 1);
    H->chain_base[0] = 0;
    for (size_t c = 0; c < S->n_chn; c++) {
        int32_t nb = (S->chn[c].count + H->block_size - 1) / H->block_size;
        if (nb < 1 || H->chain_base[c] > INT32_MAX - nb) return -1;
        H->chain_base[c + 1] = H->chain_base[c] + nb;
    }
    H->n_sample_agg = H->chain_base[S->n_chn];
    H->n_fine_agg = H->n_sample_agg;
    if (rib_amg_build_initial(arena, S, P, G, R, pin_idx, n_pins, pin_w,
                              H, &H->level[0]) != 0) return -1;
    fprintf(stderr, "    GMG build: initial level complete\n");
    H->nlevel = 1;
    while (H->nlevel < RIB_AMG_MAX_LEVELS &&
           H->level[H->nlevel - 1].n > RIB_AMG_COARSE_N) {
        RibAMGLevel *F = &H->level[H->nlevel - 1];
        RibAMGLevel *C = &H->level[H->nlevel];
        /* True geometric semi-coarsening follows physical grid spacing, not
         * the current IRLS magnitude.  A weak V operator makes slice-to-slice
         * oscillation a LOW-energy mode; refusing ever to coarsen V leaves that
         * mode on an enormous "coarse" grid.  On the 10x ribbon the old U-only
         * rule stopped at 87,139 nodes and allocated 779M dense Cholesky values.
         *
         * Try the smaller-spacing axis first.  A sparse/cut topology can make a
         * nominal axis ineffective, so try the other axis only when the first
         * cannot reduce at all.  Even a small real reduction must be retained:
         * on sparse ribbons it can expose new matches on the other axis at the
         * next level.  The arena mark is essential: a rejected trial must not
         * retain its aggregate arrays. */
        char axis = F->h_u <= F->h_v ? 'U' : 'V';
        int accepted = 0;
        for (int attempt = 0; attempt < 2 && !accepted; attempt++) {
            Arena_Mark step_mark = Arena_save(arena);
            double step_t0 = ves_clock_sec();
            int rc = rib_amg_build_next(arena, F, C, axis);
            if (rc < 0) return -1;
            if (rc == 0 && C->n < F->n) {
                fprintf(stderr,
                        "    GMG build: L%d via %c %d -> %d nodes "
                        "(topoU=%zu topoV=%zu smooth=%d, %.2fs)\n",
                        H->nlevel, axis, F->n, C->n,
                        C->nu_edge, C->nv_edge, F->smooth_sweeps,
                        ves_clock_sec() - step_t0);
                accepted = 1;
                break;
            }
            fprintf(stderr,
                    "    GMG build: rejected L%d via %c %d -> %d nodes (%.2fs)\n",
                    H->nlevel, axis, F->n, rc == 0 ? C->n : F->n,
                    ves_clock_sec() - step_t0);
            F->agg = NULL;
            F->n_coarse = 0;
            Arena_restore(arena, step_mark);
            memset(C, 0, sizeof(*C));
            axis = axis == 'U' ? 'V' : 'U';
        }
        if (!accepted) break;
        H->nlevel++;
    }
    fprintf(stderr, "    GMG build: factoring L%d (%d nodes)\n",
            H->nlevel - 1, H->level[H->nlevel - 1].n);
    double factor_t0 = ves_clock_sec();
    if (rib_amg_factor_coarsest(arena, &H->level[H->nlevel - 1]) != 0)
        return -1;
    fprintf(stderr, "    GMG build: coarse factor complete %.2fs\n",
            ves_clock_sec() - factor_t0);
    H->enabled = 1;
    return 0;
}

static int rib_amg_level_threads(const RibAMG *H, int n)
{
    int nt = H->threads;
    int useful = (n + 8191) / 8192;
    if (useful < 1) useful = 1;
    if (nt > useful) nt = useful;
    return nt;
}

static void rib_amg_matvec(const RibAMG *H, const RibAMGLevel *L,
                           const double *x, double *y, double *row_dot)
{
    int ii, nt = rib_amg_level_threads(H, L->n);
    int rr, nrt = rib_amg_level_threads(H, (int)L->nhrow);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nrt) if(nrt > 1)
#endif
    for (rr = 0; rr < (int)L->nhrow; rr++) {
        const RibAMGHyperRow *row = &L->hrow[rr];
        double s = 0.0;
        for (int32_t j = 0; j < row->count; j++) {
            const RibAMGHyperCoeff *c =
                &L->hcoeff[row->first + (size_t)j];
            s += c->value * x[c->var];
        }
        row_dot[rr] = s;
    }
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nt) if(nt > 1)
#endif
    for (ii = 0; ii < L->n; ii++) {
        double s = L->diag[ii] * x[ii];
        for (size_t e = L->off[ii]; e < L->off[ii + 1]; e++)
            s -= L->aw[e] * x[L->nbr[e]];
        for (size_t q = L->var_hcoeff_off[ii];
             q < L->var_hcoeff_off[ii + 1]; q++) {
            const RibAMGHyperCoeff *c = &L->hcoeff[L->var_hcoeff[q]];
            const RibAMGHyperRow *row = &L->hrow[c->row];
            /* The exact diagonal w*c_i^2 is already in L->diag. */
            s += row->weight * c->value *
                 (row_dot[c->row] - c->value * x[ii]);
        }
        y[ii] = s;
    }
}

static void rib_amg_chol_solve(const RibAMGLevel *L,
                                const double *r, double *x)
{
    for (int bid = 0; bid < L->n_chol_blocks; bid++) {
        int32_t bo = L->chol_block_off[bid];
        int n = L->chol_block_off[bid + 1] - bo;
        const int32_t *node = L->chol_node + bo;
        const double *C = L->chol + L->chol_off[bid];
        for (int i = 0; i < n; i++) {
            double s = r[node[i]];
            for (int k = 0; k < i; k++)
                s -= C[(size_t)i * (size_t)n + (size_t)k] * x[node[k]];
            x[node[i]] = s / C[(size_t)i * (size_t)n + (size_t)i];
        }
        for (int i = n - 1; i >= 0; i--) {
            double s = x[node[i]];
            for (int k = i + 1; k < n; k++)
                s -= C[(size_t)k * (size_t)n + (size_t)i] * x[node[k]];
            x[node[i]] = s / C[(size_t)i * (size_t)n + (size_t)i];
        }
    }
}

static void rib_amg_cheb_smooth(const RibAMG *H, const RibAMGLevel *L,
                                int lev, const double *r, double *x,
                                RibAMGWork *W)
{
    int ii, nt = rib_amg_level_threads(H, L->n);
    double *tmp = W->tmp[lev], *p = W->smooth_p[lev];
    double lambda_hi = L->lambda_max;
    double lambda_lo = lambda_hi / RIB_GMG_CHEB_RATIO;
    double theta = 0.5 * (lambda_hi + lambda_lo);
    double delta = 0.5 * (lambda_hi - lambda_lo);
    double sigma = theta / delta;
    double rho = 1.0 / sigma;

    rib_amg_matvec(H, L, x, tmp, W->row_dot[lev]);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nt) if(nt > 1)
#endif
    for (ii = 0; ii < L->n; ii++) {
        p[ii] = (r[ii] - tmp[ii]) / (theta * L->diag[ii]);
        x[ii] += p[ii];
    }
    int smooth_sweeps = L->smooth_sweeps > 0
                      ? L->smooth_sweeps : RIB_AMG_SMOOTH;
    for (int sweep = 1; sweep < smooth_sweeps; sweep++) {
        double rho_new = 1.0 / (2.0 * sigma - rho);
        double beta = rho_new * rho;
        double alpha = 2.0 * rho_new / delta;
        rib_amg_matvec(H, L, x, tmp, W->row_dot[lev]);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nt) if(nt > 1)
#endif
        for (ii = 0; ii < L->n; ii++) {
            p[ii] = beta * p[ii] +
                    alpha * (r[ii] - tmp[ii]) / L->diag[ii];
            x[ii] += p[ii];
        }
        rho = rho_new;
    }
}

static void rib_amg_vcycle(const RibAMG *H, int lev, RibAMGWork *W)
{
    const RibAMGLevel *L = &H->level[lev];
    double *x = W->x[lev], *r = W->r[lev], *tmp = W->tmp[lev];
    memset(x, 0, (size_t)L->n * sizeof(double));
    if (lev == H->nlevel - 1) {
        rib_amg_chol_solve(L, r, x);
        return;
    }
    int ii, nt = rib_amg_level_threads(H, L->n);
    rib_amg_cheb_smooth(H, L, lev, r, x, W);
    rib_amg_matvec(H, L, x, tmp, W->row_dot[lev]);
    memset(W->r[lev + 1], 0,
           (size_t)H->level[lev + 1].n * sizeof(double));
    /* P^T restriction.  This coarse loop is small relative to the fine line
     * smoother; keeping it serial avoids atomic accumulation and is fully
     * deterministic. */
    for (int i = 0; i < L->n; i++)
        W->r[lev + 1][L->agg[i]] += r[i] - tmp[i];
    rib_amg_vcycle(H, lev + 1, W);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nt) if(nt > 1)
#endif
    for (ii = 0; ii < L->n; ii++) x[ii] += W->x[lev + 1][L->agg[ii]];
    rib_amg_cheb_smooth(H, L, lev, r, x, W);
}

/* Full-multigrid nested iteration on the stored geometric hierarchy.  The
 * hierarchy solves a defect/correction equation.  The caller's fine reference
 * therefore retains the already-good per-slice arc-length ramp while FMG
 * supplies only the smooth cross-slice registration correction represented by
 * the coarse space. */
static int rib_amg_fmg_solve(const RibAMG *H, const SliceSet *S,
                             RibAMGWork *W, const double *rhs0, double *x0,
                             const double *fine_reference,
                             RibGMGSnapshots *capture)
{
    int last = H->nlevel - 1;
    memcpy(W->fmg_rhs[0], rhs0,
           (size_t)H->level[0].n * sizeof(double));
    for (int lev = 0; lev < last; lev++) {
        const RibAMGLevel *F = &H->level[lev];
        const RibAMGLevel *C = &H->level[lev + 1];
        memset(W->fmg_rhs[lev + 1], 0, (size_t)C->n * sizeof(double));
        for (int i = 0; i < F->n; i++)
            W->fmg_rhs[lev + 1][F->agg[i]] += W->fmg_rhs[lev][i];
    }
    memset(W->fmg_x[last], 0,
           (size_t)H->level[last].n * sizeof(double));
    rib_amg_chol_solve(&H->level[last], W->fmg_rhs[last], W->fmg_x[last]);
    if (rib_gmg_capture_stored_level(
            H, S, W, fine_reference, last, capture) != 0) return -1;

    for (int lev = last - 1; lev >= 0; lev--) {
        const RibAMGLevel *F = &H->level[lev];
        for (int i = 0; i < F->n; i++)
            W->fmg_x[lev][i] = W->fmg_x[lev + 1][F->agg[i]];
        for (int cycle = 0; cycle < RIB_GMG_FMG_CYCLES; cycle++) {
            rib_amg_matvec(H, F, W->fmg_x[lev], W->tmp[lev],
                           W->row_dot[lev]);
            for (int i = 0; i < F->n; i++)
                W->r[lev][i] = W->fmg_rhs[lev][i] - W->tmp[lev][i];
            rib_amg_vcycle(H, lev, W);
            for (int i = 0; i < F->n; i++)
                W->fmg_x[lev][i] += W->x[lev][i];
        }
        if (rib_gmg_capture_stored_level(
                H, S, W, fine_reference, lev, capture) != 0) return -1;
    }
    memcpy(x0, W->fmg_x[0],
           (size_t)H->level[0].n * sizeof(double));
    return 0;
}

static void rib_amg_work_init(Arena_T arena, const RibAMG *H,
                              RibAMGWork *W)
{
    memset(W, 0, sizeof(*W));
    for (int l = 0; l < H->nlevel; l++) {
        size_t n = (size_t)H->level[l].n;
        W->x[l] = RIB_ALLOC_ARRAY(arena, double, n + 1);
        W->r[l] = RIB_ALLOC_ARRAY(arena, double, n + 1);
        W->tmp[l] = RIB_ALLOC_ARRAY(arena, double, n + 1);
        W->smooth_p[l] = RIB_ALLOC_ARRAY(arena, double, n + 1);
        W->row_dot[l] = RIB_ALLOC_ARRAY(
            arena, double, H->level[l].nhrow + 1);
        W->fmg_rhs[l] = RIB_ALLOC_ARRAY(arena, double, n + 1);
        W->fmg_x[l] = RIB_ALLOC_ARRAY(arena, double, n + 1);
    }
}

typedef struct {
    const SliceSet  *S;
    const PairSet   *P;
    const SolveGraph *G;
    const RibStripSet *strip; /* NULL for legacy chain-gauge-only solves */
    const RibAMG    *amg;
    /* one gauge pin PER connected component of the solve graph: without its
     * own pin a disconnected component's absolute u is undetermined (every
     * energy term is a difference) and CG parks it at an arbitrary offset --
     * the source of detached charts landing mid-band. */
    const int32_t *pin_idx; const double *pin_val;
    size_t n_pins; double pin_w;
    const double *pin_weight; /* optional [n_pins], else scalar pin_w */
    double *strip_dot;        /* [strip->nrun], matrix-free row scratch */
    double *strip_mean;       /* [strip->nrun], alignment scratch */
    int strip_final;
    double ridge_weight;      /* ADMM u=z split: rho I */
    const double *ridge_target; /* z-y in scaled-dual form */
    int threads;
    double tol;             /* <= 0 selects the legacy size-dependent default */
    /* TAUCS Cholesky of the sparse core (graph edges + length/align diagonal
     * + ridge + pins, with each run's alignment-to-mean term made exact via
     * one auxiliary mean variable).  Used as the PCG preconditioner: the true
     * operator stays matrix-free apply_A, so the per-run rank-1 length terms
     * and any ADMM rho adaptation only cost iterations, never correctness.
     * NULL selects the legacy GMG/Jacobi path. */
    SparseFactor_T chol;
    int chol_attempted;       /* lazy build: rounds with no solves never pay */
    int round_no;             /* current round, for the factor log line */
    double *chol_rhs;         /* [chol_n] scratch: [r ; 0] (malloc'd) */
    double *chol_x;           /* [chol_n] scratch (malloc'd) */
    double *chol_x2;          /* [chol_n] SMW second-solve scratch (malloc'd) */
    size_t chol_n;            /* n_smp + n_aux_runs */
    /* Exact Sherman-Morrison-Woodbury correction for the length rows too
     * long to clique into the factored matrix:
     *   (M + U W U^T)^{-1} = M^{-1} - M^{-1} U C^{-1} U^T M^{-1},
     *   C = W^{-1} + U^T M^{-1} U   (dense k x k, Cholesky-factored).
     * U's columns are the raw merged length rows (x-space, aux rows 0). */
    size_t smw_k;
    int32_t *smw_off;         /* [k+1] CSR into smw_var/smw_val (malloc'd) */
    int32_t *smw_var;
    double *smw_val;
    double *smw_cchol;        /* [k*k] in-place lower Cholesky of C */
    double *smw_s;            /* [k] scratch */
} SolveCtx;

/* In-place dense Cholesky (lower) of the k x k SPD matrix C (row-major).
 * Returns 0 on success, -1 if a pivot fails (not numerically SPD). */
static int rib_dense_chol(double *C, size_t k)
{
    for (size_t j = 0; j < k; j++) {
        double d = C[j * k + j];
        for (size_t p = 0; p < j; p++) d -= C[j * k + p] * C[j * k + p];
        if (!(d > 1e-30)) return -1;
        d = sqrt(d);
        C[j * k + j] = d;
        for (size_t i = j + 1; i < k; i++) {
            double s = C[i * k + j];
            for (size_t p = 0; p < j; p++)
                s -= C[i * k + p] * C[j * k + p];
            C[i * k + j] = s / d;
        }
    }
    return 0;
}

static void rib_dense_chol_solve(const double *C, size_t k, double *s)
{
    for (size_t i = 0; i < k; i++) {
        double v = s[i];
        for (size_t p = 0; p < i; p++) v -= C[i * k + p] * s[p];
        s[i] = v / C[i * k + i];
    }
    for (size_t ii = k; ii > 0; ii--) {
        size_t i = ii - 1;
        double v = s[i];
        for (size_t p = i + 1; p < k; p++) v -= C[p * k + i] * s[p];
        s[i] = v / C[i * k + i];
    }
}

static int solve_thread_count(int requested)
{
#ifdef _OPENMP
    int n = requested > 0 ? requested : omp_get_max_threads();
    if (n < 1) n = 1;
    return n;
#else
    (void)requested;
    return 1;
#endif
}

static int solve_graph_build(Arena_T arena, const SliceSet *S,
                             const PairSet *P, size_t pair_begin,
                             SolveGraph *G)
{
    size_t n = S->n_smp;
    size_t chain_edges = 0, active_pairs = 0;
    memset(G, 0, sizeof(*G));
    G->n = n;
    for (size_t c = 0; c < S->n_chn; c++)
        if (S->chn[c].count > 1) chain_edges += (size_t)S->chn[c].count - 1;
    if (pair_begin > P->n_pairs) return -1;
    for (size_t p = pair_begin; p < P->n_pairs; p++)
        if (P->pairs[p].like >= 0.0) active_pairs++;
    if (P->n_pairs > (size_t)INT32_MAX ||
        chain_edges + active_pairs > SIZE_MAX / 2) return -1;

    G->off = (size_t *)ARENA_CALLOC(arena, n + 1, sizeof(size_t));
    G->active_pair = RIB_ALLOC_ARRAY(arena, int32_t, active_pairs + 1);
    size_t ap = 0;
    for (size_t c = 0; c < S->n_chn; c++) {
        int32_t f = S->chn[c].first, cn = S->chn[c].count;
        for (int32_t k = 1; k < cn; k++) {
            G->off[(size_t)(f + k) + 1]++;
            G->off[(size_t)(f + k - 1) + 1]++;
        }
    }
    for (size_t p = pair_begin; p < P->n_pairs; p++) {
        const Pair *pr = &P->pairs[p];
        if (pr->like < 0.0) continue;
        if (pr->a < 0 || pr->b < 0 || (size_t)pr->a >= n || (size_t)pr->b >= n)
            return -1;
        G->off[(size_t)pr->a + 1]++;
        G->off[(size_t)pr->b + 1]++;
        G->active_pair[ap++] = (int32_t)p;
    }
    for (size_t i = 0; i < n; i++) G->off[i + 1] += G->off[i];
    G->n_adj = G->off[n];
    G->n_active_pair = ap;
    G->nbr = RIB_ALLOC_ARRAY(arena, int32_t, G->n_adj + 1);
    G->pair_ref = RIB_ALLOC_ARRAY(arena, int32_t, G->n_adj + 1);
    size_t *cur = RIB_ALLOC_ARRAY(arena, size_t, n + 1);
    memcpy(cur, G->off, n * sizeof(size_t));
    for (size_t c = 0; c < S->n_chn; c++) {
        int32_t f = S->chn[c].first, cn = S->chn[c].count;
        for (int32_t k = 1; k < cn; k++) {
            int32_t a = f + k, b = a - 1;
            size_t ea = cur[a]++, eb = cur[b]++;
            G->nbr[ea] = b; G->pair_ref[ea] = -1;
            G->nbr[eb] = a; G->pair_ref[eb] = -1;
        }
    }
    for (size_t q = 0; q < ap; q++) {
        int32_t p = G->active_pair[q];
        const Pair *pr = &P->pairs[p];
        size_t ea = cur[pr->a]++, eb = cur[pr->b]++;
        G->nbr[ea] = pr->b; G->pair_ref[ea] = p;
        G->nbr[eb] = pr->a; G->pair_ref[eb] = p;
    }
    G->n_active_pair = ap;
    return 0;
}

static double solve_graph_weight(const SolveCtx *cx, size_t i,
                                 int32_t j, int32_t pref)
{
    if (pref >= 0) {
        const Pair *pr = &cx->P->pairs[pref];
        if (pr->like <= 0.0) return 0.0;
        return cx->strip != NULL ? cx->strip->cont_weight * pr->like
                                 : pr->w * pr->like;
    }
    {
        double l = fabs(cx->S->smp[i].s - cx->S->smp[j].s);
        double length_eps = cx->strip != NULL ? RIB_STRIP_GEOM_EPS : 1e-9;
        if (l < length_eps) l = length_eps;
        return (cx->strip != NULL ? cx->strip->local_weight : 1.0) / l;
    }
}

static void apply_A(const SolveCtx *cx, const double *x, double *y)
{
    const SolveGraph *G = cx->G;
    const RibStripSet *R = cx->strip;
    int nt = cx->threads;
    int ii;
    if (R != NULL) {
        int rr;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nt)
#endif
        for (rr = 0; rr < (int)R->nrun; rr++) {
            const RibStripRun *run = &R->run[rr];
            double dot = 0.0, mean = 0.0;
            size_t k0 = 2 * run->first;
            size_t k1 = 2 * (run->first + (size_t)run->count);
            for (size_t k = k0; k < k1; k++)
                dot += R->coeff[k].value * x[R->coeff[k].var];
            for (int32_t j = 0; j < run->count; j++) {
                const RibStripMember *m =
                    &R->member[run->first + (size_t)j];
                mean += m->base_weight * m->like * x[m->sample];
            }
            cx->strip_dot[rr] = dot;
            cx->strip_mean[rr] = run->like_sum > 0.0
                               ? mean / run->like_sum : 0.0;
        }
    }
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nt)
#endif
    for (ii = 0; ii < (int)G->n; ii++) {
        size_t i = (size_t)ii;
        double yi = 0.0, xi = x[i];
        for (size_t e = G->off[i]; e < G->off[i + 1]; e++) {
            int32_t j = G->nbr[e];
            double w = solve_graph_weight(cx, i, j, G->pair_ref[e]);
            yi += w * (xi - x[j]);
        }
        if (R != NULL) {
            for (size_t q = R->var_coeff_off[i];
                 q < R->var_coeff_off[i + 1]; q++) {
                const RibStripCoeff *c = &R->coeff[R->var_coeff[q]];
                yi += R->length_weight * R->run[c->run].weight *
                      c->value * cx->strip_dot[c->run];
            }
            int32_t mi = R->sample_member[i];
            if (mi >= 0) {
                const RibStripMember *m = &R->member[mi];
                const RibStripRun *run = &R->run[m->run];
                yi += R->align_weight * run->weight *
                      m->base_weight * m->like *
                      (xi - cx->strip_mean[m->run]);
            }
        }
        if (cx->ridge_weight > 0.0)
            yi += cx->ridge_weight * xi;
        y[i] = yi;
    }
    for (size_t p = 0; p < cx->n_pins; p++)
        y[cx->pin_idx[p]] +=
            (cx->pin_weight != NULL ? cx->pin_weight[p] : cx->pin_w) *
            x[cx->pin_idx[p]];
}

static void build_b_diag(const SolveCtx *cx, double *b, double *diag)
{
    const SolveGraph *G = cx->G;
    const RibStripSet *R = cx->strip;
    int nt = cx->threads;
    int ii;
    if (R != NULL) {
        int rr;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nt)
#endif
        for (rr = 0; rr < (int)R->nrun; rr++) {
            const RibStripRun *run = &R->run[rr];
            double mean = 0.0;
            if (cx->strip_final) {
                for (int32_t j = 0; j < run->count; j++) {
                    const RibStripMember *m =
                        &R->member[run->first + (size_t)j];
                    double offset = v3dot(
                        run->tangent, cx->S->smp[m->sample].p);
                    mean += m->base_weight * m->like * offset;
                }
            }
            cx->strip_mean[rr] = run->like_sum > 0.0
                               ? mean / run->like_sum : 0.0;
        }
    }
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nt)
#endif
    for (ii = 0; ii < (int)G->n; ii++) {
        size_t i = (size_t)ii;
        double bi = 0.0, di = 0.0;
        for (size_t e = G->off[i]; e < G->off[i + 1]; e++) {
            int32_t j = G->nbr[e], pref = G->pair_ref[e];
            double w = solve_graph_weight(cx, i, j, pref);
            double d;
            if (pref >= 0) {
                const Pair *pr = &cx->P->pairs[pref];
                d = ((int32_t)i == pr->a) ? pr->d : -pr->d;
            } else {
                d = cx->S->smp[i].s - cx->S->smp[j].s;
            }
            bi += w * d;
            di += w;
        }
        if (R != NULL) {
            for (size_t q = R->var_coeff_off[i];
                 q < R->var_coeff_off[i + 1]; q++) {
                const RibStripCoeff *c = &R->coeff[R->var_coeff[q]];
                double w = R->length_weight * R->run[c->run].weight;
                bi += w * c->value;
                di += w * c->value * c->value;
            }
            int32_t mi = R->sample_member[i];
            if (mi >= 0) {
                const RibStripMember *m = &R->member[mi];
                const RibStripRun *run = &R->run[m->run];
                double mw = m->base_weight * m->like;
                double w = R->align_weight * run->weight;
                double offset = cx->strip_final
                              ? v3dot(run->tangent, cx->S->smp[i].p) : 0.0;
                bi += w * mw * (offset - cx->strip_mean[m->run]);
                if (run->like_sum > 0.0)
                    di += w * mw * (1.0 - mw / run->like_sum);
            }
        }
        if (cx->ridge_weight > 0.0) {
            bi += cx->ridge_weight * cx->ridge_target[i];
            di += cx->ridge_weight;
        }
        b[i] = bi;
        diag[i] = di;
    }
    for (size_t p = 0; p < cx->n_pins; p++) {
        double w = cx->pin_weight != NULL ? cx->pin_weight[p] : cx->pin_w;
        b[cx->pin_idx[p]]    += w * cx->pin_val[p];
        diag[cx->pin_idx[p]] += w;
    }
}

static void rib_line_solve(const SolveCtx *cx, const double *diag,
                           const double *r, double *z, RibAMGWork *W)
{
    size_t n = cx->S->n_smp;
    int ii, nt = cx->threads;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nt)
#endif
    for (ii = 0; ii < (int)n; ii++) z[ii] = r[ii] / diag[ii];
    int cc;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 8) num_threads(nt)
#endif
    for (cc = 0; cc < (int)cx->S->n_chn; cc++) {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        double *cp = W->line_cp + (size_t)tid * (size_t)W->maxchain;
        int32_t f = cx->S->chn[cc].first, cn = cx->S->chn[cc].count;
        if (cn <= 0) continue;
        for (int32_t k = 0; k < cn; k++) {
            int32_t i = f + k;
            double a = 0.0, c = 0.0;
            if (k > 0) {
                double l = cx->S->smp[i].s - cx->S->smp[i - 1].s;
                double length_eps = cx->strip != NULL
                                  ? RIB_STRIP_GEOM_EPS : 1e-9;
                if (l < length_eps) l = length_eps;
                a = -(cx->strip != NULL ? cx->strip->local_weight : 1.0) / l;
            }
            if (k + 1 < cn) {
                double l = cx->S->smp[i + 1].s - cx->S->smp[i].s;
                double length_eps = cx->strip != NULL
                                  ? RIB_STRIP_GEOM_EPS : 1e-9;
                if (l < length_eps) l = length_eps;
                c = -(cx->strip != NULL ? cx->strip->local_weight : 1.0) / l;
            }
            double den = diag[i] - (k > 0 ? a * cp[k - 1] : 0.0);
            if (den < 1e-12) den = 1e-12;
            cp[k] = (k + 1 < cn) ? c / den : 0.0;
            z[i] = (r[i] - (k > 0 ? a * z[i - 1] : 0.0)) / den;
        }
        for (int32_t k = cn - 1; k-- > 0; )
            z[f + k] -= cp[k] * z[f + k + 1];
    }
}

static void rib_apply_preconditioner(const SolveCtx *cx, const double *diag,
                                     const double *r, double *z,
                                     RibAMGWork *W)
{
    size_t n = cx->S->n_smp;
    int ii, cc, nt = cx->threads;
    if (cx->chol != NULL) {
        /* Exact solve of the factored core: copy the residual into the
         * augmented rhs (auxiliary run-mean rows get 0 -- their Schur
         * complement reproduces the alignment variance term exactly), solve,
         * then apply the SMW correction for the un-cliqued long length rows,
         * and keep the x-block of the solution. */
        memcpy(cx->chol_rhs, r, n * sizeof(double));
        if (cx->chol_n > n)
            memset(cx->chol_rhs + n, 0, (cx->chol_n - n) * sizeof(double));
        if (Sparse_factor_solve(cx->chol, cx->chol_rhs, cx->chol_x) == 0) {
            if (cx->smw_k > 0) {
                size_t k = cx->smw_k;
                for (size_t i = 0; i < k; i++) {
                    double dot = 0.0;
                    for (int32_t q = cx->smw_off[i];
                         q < cx->smw_off[i + 1]; q++)
                        dot += cx->smw_val[q] *
                               cx->chol_x[cx->smw_var[q]];
                    cx->smw_s[i] = dot;
                }
                rib_dense_chol_solve(cx->smw_cchol, k, cx->smw_s);
                memset(cx->chol_rhs, 0, cx->chol_n * sizeof(double));
                for (size_t i = 0; i < k; i++) {
                    double si = cx->smw_s[i];
                    if (!(si < 0.0 || si > 0.0)) continue;
                    for (int32_t q = cx->smw_off[i];
                         q < cx->smw_off[i + 1]; q++)
                        cx->chol_rhs[cx->smw_var[q]] +=
                            cx->smw_val[q] * si;
                }
                if (Sparse_factor_solve(cx->chol, cx->chol_rhs,
                                        cx->chol_x2) == 0) {
                    for (size_t i = 0; i < n; i++)
                        z[i] = cx->chol_x[i] - cx->chol_x2[i];
                    return;
                }
                /* second solve failed: uncorrected core is still a valid
                 * (weaker) preconditioner */
            }
            memcpy(z, cx->chol_x, n * sizeof(double));
            return;
        }
        fprintf(stderr,
                "    solve TAUCS: factored solve failed; Jacobi fallback\n");
        /* fall through */
    }
    if (cx->chol != NULL || cx->amg == NULL || !cx->amg->enabled || W == NULL) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nt)
#endif
        for (ii = 0; ii < (int)n; ii++) z[ii] = r[ii] / diag[ii];
        return;
    }

    /* Symmetric multiplicative two-level cycle: exact line pre-smooth,
     * Galerkin coarse correction of its residual, exact line post-smooth. */
    rib_line_solve(cx, diag, r, z, W);
    apply_A(cx, z, W->fine_tmp);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nt)
#endif
    for (ii = 0; ii < (int)n; ii++) W->fine_tmp[ii] = r[ii] - W->fine_tmp[ii];

    /* Algebraic coarse correction, restricted/prolonged by chain blocks. */
    memset(W->r[0], 0,
           (size_t)cx->amg->level[0].n * sizeof(double));
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 8) num_threads(nt)
#endif
    for (cc = 0; cc < (int)cx->S->n_chn; cc++) {
        int32_t f = cx->S->chn[cc].first, cn = cx->S->chn[cc].count;
        int32_t base = cx->amg->chain_base[cc];
        for (int32_t k0 = 0, b = 0; k0 < cn; k0 += cx->amg->block_size, b++) {
            int32_t k1 = k0 + cx->amg->block_size;
            double s = 0.0;
            if (k1 > cn) k1 = cn;
            for (int32_t k = k0; k < k1; k++) s += W->fine_tmp[f + k];
            W->r[0][base + b] = s;
        }
    }
    rib_amg_vcycle(cx->amg, 0, W);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 8) num_threads(nt)
#endif
    for (cc = 0; cc < (int)cx->S->n_chn; cc++) {
        int32_t f = cx->S->chn[cc].first, cn = cx->S->chn[cc].count;
        int32_t base = cx->amg->chain_base[cc];
        for (int32_t k = 0; k < cn; k++)
            z[f + k] += RIB_AMG_COARSE_SCALE *
                        W->x[0][base + k / cx->amg->block_size];
    }
    apply_A(cx, z, W->fine_tmp);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nt)
#endif
    for (ii = 0; ii < (int)n; ii++) W->fine_tmp[ii] = r[ii] - W->fine_tmp[ii];
    rib_line_solve(cx, diag, W->fine_tmp, W->fine_corr, W);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nt)
#endif
    for (ii = 0; ii < (int)n; ii++) z[ii] += W->fine_corr[ii];
}

/* Recompute the true residual instead of extending the floating-point Krylov
 * recurrence indefinitely.  At production scale the latter can lose
 * conjugacy near 1e-6 even though the system is still making geometric
 * progress. */
static double rib_true_residual(const SolveCtx *cx, const double *b,
                                const double *u, double *Ap, double *r,
                                double bnorm)
{
    size_t n = cx->S->n_smp;
    int ii, nt = cx->threads;
    double rn = 0.0;
    apply_A(cx, u, Ap);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:rn) num_threads(nt)
#endif
    for (ii = 0; ii < (int)n; ii++) {
        r[ii] = b[ii] - Ap[ii];
        rn += r[ii] * r[ii];
    }
    return sqrt(rn) / bnorm;
}

/* Full-clique assembly cap: a run's length row expands to ~(distinct vars)^2/2
 * factored entries; the cap ladder lowers the per-run clique size until the
 * total fits this budget.  Rows above the chosen cap take the exact SMW tier
 * (up to RIB_CHOL_SMW_MAX rows; the dense C setup costs one triangular solve
 * per row per round); any excess degrades to per-member springs.
 *
 * MEASURED (4x5x5 soup, 2.04M nodes, 25,663 runs): the run-length tail is so
 * heavy that cliques<=256 + 512 SMW rows still left 57% of length mass on
 * springs while quadrupling factor cost (nnz 8.2M -> 124M, factor 81s ->
 * 340s, trisolve 0.63s -> 3.5s/iter).  Springs-everywhere is the measured
 * operating point; both caps are 0 to select it.  The machinery stays for
 * regimes with short-run-dominated rows. */
#define RIB_CHOL_CLIQUE_CAP 0
#define RIB_CHOL_SMW_MAX    0

static void rib_strip_chol_release(SolveCtx *cx);

/* Collect run rr's length-row coefficients merged by DISTINCT variable
 * (insertion sort by var + in-place merge; chains arrive near-sorted).
 * Returns the merged count. */
static size_t rib_chol_run_merge(const RibStripSet *R, size_t rr,
                                 int32_t *cl_var, double *cl_val)
{
    const RibStripRun *run = &R->run[rr];
    size_t k0 = 2 * run->first;
    size_t k1 = 2 * (run->first + (size_t)run->count);
    size_t nm = 0, out_n = 0;
    for (size_t k = k0; k < k1; k++) {
        const RibStripCoeff *c = &R->coeff[k];
        if (!(c->value < 0.0 || c->value > 0.0)) continue;
        cl_var[nm] = c->var;
        cl_val[nm] = c->value;
        nm++;
    }
    for (size_t a = 1; a < nm; a++) {
        int32_t kv = cl_var[a];
        double kx = cl_val[a];
        size_t b = a;
        while (b > 0 && cl_var[b - 1] > kv) {
            cl_var[b] = cl_var[b - 1];
            cl_val[b] = cl_val[b - 1];
            b--;
        }
        cl_var[b] = kv;
        cl_val[b] = kx;
    }
    for (size_t a = 0; a < nm; a++) {
        if (out_n > 0 && cl_var[out_n - 1] == cl_var[a]) {
            cl_val[out_n - 1] += cl_val[a];
        } else {
            cl_var[out_n] = cl_var[a];
            cl_val[out_n] = cl_val[a];
            out_n++;
        }
    }
    return out_n;
}

/* Assemble and factor the Stage-C operator's sparse form with TAUCS:
 * graph edges (chain + continuation pairs, live IRLS weights), the per-run
 * length row as its full clique over distinct variables (runs are
 * sample-disjoint, so cliques never overlap; springs fallback past the size
 * cap), the per-run alignment term made EXACT through one auxiliary mean
 * variable per run (its Schur complement reproduces the weighted-variance
 * Hessian), the ADMM ridge, and the per-component gauge pins.  With cliques
 * the factorization IS the operator at build-time rho: PCG converges in a
 * few iterations and only mops up ADMM rho adaptation.  Returns 0 on success
 * with cx->chol/chol_rhs/chol_x/chol_n filled; on failure cx->chol stays
 * NULL and the caller keeps the legacy GMG/Jacobi path. */
static int rib_strip_chol_build(SolveCtx *cx, double ridge, int round)
{
    const SolveGraph *G = cx->G;
    const RibStripSet *R = cx->strip;
    size_t n = cx->S->n_smp;
    size_t naux = (R != NULL) ? R->nrun : 0;
    size_t n_aug = n + naux;
    double t0 = ves_clock_sec();
    size_t nt = 0, max_nt = 0, clique_nt = 0, max_run = 0;
    size_t run_cap = 0, smw_support = 0, sk = 0, spring_members = 0;
    int *t_row = NULL, *t_col = NULL;
    int32_t *cl_var = NULL, *m_arr = NULL;
    int32_t *s_off = NULL, *s_var = NULL;
    double *t_val = NULL, *dacc = NULL, *cl_val = NULL, *s_val = NULL;
    double *s_w = NULL;
    int rc = -1;
    cx->chol = NULL;
    cx->chol_rhs = NULL;
    cx->chol_x = NULL;
    cx->chol_x2 = NULL;
    cx->chol_n = 0;
    cx->smw_k = 0;
    cx->smw_off = NULL;
    cx->smw_var = NULL;
    cx->smw_val = NULL;
    cx->smw_cchol = NULL;
    cx->smw_s = NULL;
    if (n_aug == 0 || n_aug > (size_t)INT32_MAX / 2u) return -1;
    max_nt = G->n_adj / 2u + n_aug + 8u;
    if (R != NULL) {
        static const size_t caps[8] = { 1024, 512, 256, 128, 64, 32, 16, 8 };
        max_nt += R->nmember;              /* align aux couplings */
        for (size_t rr = 0; rr < R->nrun; rr++) {
            size_t nm = 2u * (size_t)R->run[rr].count;
            if (nm > max_run) max_run = nm;
        }
        cl_var = (int32_t *)malloc((max_run + 1) * sizeof *cl_var);
        cl_val = (double *)malloc((max_run + 1) * sizeof *cl_val);
        m_arr = (int32_t *)calloc(R->nrun + 1, sizeof *m_arr);
        if (cl_var == NULL || cl_val == NULL || m_arr == NULL) goto done;
        for (size_t rr = 0; rr < R->nrun; rr++) {
            if (!(R->length_weight * R->run[rr].weight > 0.0)) continue;
            m_arr[rr] = (int32_t)rib_chol_run_merge(R, rr, cl_var, cl_val);
        }
        /* cap ladder: the largest per-run clique size whose total fits the
         * memory budget; longer rows go to the SMW tier (exact via dense C),
         * and only past RIB_CHOL_SMW_MAX do rows degrade to springs. */
        for (size_t ci = 0; ci < 8; ci++) {
            run_cap = caps[ci];
            clique_nt = 0;
            for (size_t rr = 0; rr < R->nrun; rr++) {
                size_t m = (size_t)m_arr[rr];
                if (m > 1 && m <= run_cap) clique_nt += m * (m - 1) / 2;
            }
            if (clique_nt <= (size_t)RIB_CHOL_CLIQUE_CAP) break;
        }
        if (clique_nt > (size_t)RIB_CHOL_CLIQUE_CAP) {
            run_cap = 0;
            clique_nt = 0;
        }
        for (size_t rr = 0; rr < R->nrun; rr++) {
            size_t m = (size_t)m_arr[rr];
            if (m == 0 || m <= run_cap) continue;
            if (sk < (size_t)RIB_CHOL_SMW_MAX) {
                sk++;
                smw_support += m;
            } else {
                spring_members += (size_t)R->run[rr].count;
            }
        }
        max_nt += clique_nt + spring_members + 8u;
        if (sk > 0) {
            s_off = (int32_t *)malloc((sk + 1) * sizeof *s_off);
            s_var = (int32_t *)malloc((smw_support + 1) * sizeof *s_var);
            s_val = (double *)malloc((smw_support + 1) * sizeof *s_val);
            s_w = (double *)malloc((sk + 1) * sizeof *s_w);
            if (s_off == NULL || s_var == NULL || s_val == NULL ||
                s_w == NULL)
                goto done;
            s_off[0] = 0;
        }
    }
    if (max_nt > (size_t)INT32_MAX) goto done;
    /* Once-per-round assembly, not a hot path: plain malloc keeps the round
     * arena free for the solver workspace. */
    t_row = (int *)malloc(max_nt * sizeof *t_row);
    t_col = (int *)malloc(max_nt * sizeof *t_col);
    t_val = (double *)malloc(max_nt * sizeof *t_val);
    dacc  = (double *)calloc(n_aug, sizeof *dacc);
    if (t_row == NULL || t_col == NULL || t_val == NULL || dacc == NULL)
        goto done;

    for (size_t i = 0; i < G->n; i++) {
        for (size_t e = G->off[i]; e < G->off[i + 1]; e++) {
            int32_t j = G->nbr[e];
            double w = solve_graph_weight(cx, i, j, G->pair_ref[e]);
            if (!(w > 0.0)) continue;
            dacc[i] += w;
            if ((size_t)j < i) {
                t_row[nt] = (int)i;
                t_col[nt] = (int)j;
                t_val[nt] = -w;
                nt++;
            }
        }
    }
    if (R != NULL) {
        /* Three tiers per length row (runs are sample-disjoint, so cliques
         * never overlap):
         *   m <= run_cap      full clique over the DISTINCT variables -- the
         *                     factored matrix carries the row exactly;
         *   SMW tier          the row is excluded from the factor and applied
         *                     exactly through the dense SMW correction;
         *   overflow          per-member springs (approximate; PCG carries
         *                     the rank-1 remainder). */
        size_t si = 0;
        for (size_t rr = 0; rr < R->nrun; rr++) {
            const RibStripRun *run = &R->run[rr];
            double w = R->length_weight * run->weight;
            size_t m = (size_t)m_arr[rr];
            if (!(w > 0.0) || m == 0) continue;
            if (m <= run_cap) {
                size_t nm = rib_chol_run_merge(R, rr, cl_var, cl_val);
                for (size_t a = 0; a < nm; a++) {
                    double va = cl_val[a];
                    dacc[(size_t)cl_var[a]] += w * va * va;
                    for (size_t b = 0; b < a; b++) {
                        double wv = w * va * cl_val[b];
                        if (!(wv < 0.0 || wv > 0.0)) continue;
                        t_row[nt] = cl_var[a];   /* sorted: var[a] > var[b] */
                        t_col[nt] = cl_var[b];
                        t_val[nt] = wv;
                        nt++;
                    }
                }
            } else if (si < sk) {
                size_t nm = rib_chol_run_merge(R, rr, cl_var, cl_val);
                int32_t base = s_off[si];
                for (size_t a = 0; a < nm; a++) {
                    s_var[base + (int32_t)a] = cl_var[a];
                    s_val[base + (int32_t)a] = cl_val[a];
                }
                s_w[si] = w;
                s_off[si + 1] = base + (int32_t)nm;
                si++;
            } else {
                for (int32_t j = 0; j < run->count; j++) {
                    size_t mm = run->first + (size_t)j;
                    const RibStripCoeff *lo = &R->coeff[2 * mm];
                    const RibStripCoeff *hi = &R->coeff[2 * mm + 1];
                    double wl = w * lo->value * lo->value;
                    double wh = w * hi->value * hi->value;
                    double wo = w * lo->value * hi->value;
                    if (lo->var == hi->var) continue;
                    if (wl > 0.0) dacc[(size_t)lo->var] += wl;
                    if (wh > 0.0) dacc[(size_t)hi->var] += wh;
                    if (wo < 0.0 || wo > 0.0) {
                        t_row[nt] = lo->var > hi->var ? lo->var : hi->var;
                        t_col[nt] = lo->var > hi->var ? hi->var : lo->var;
                        t_val[nt] = wo;
                        nt++;
                    }
                }
            }
        }
        sk = si;   /* rows actually collected */
        for (size_t rr = 0; rr < R->nrun; rr++) {
            const RibStripRun *run = &R->run[rr];
            double w = R->align_weight * run->weight;
            size_t aux = n + rr;
            if (!(w > 0.0)) continue;
            if (run->like_sum > 0.0) {
                for (int32_t j = 0; j < run->count; j++) {
                    const RibStripMember *m =
                        &R->member[run->first + (size_t)j];
                    double wm = w * m->base_weight * m->like;
                    if (!(wm > 0.0)) continue;
                    dacc[(size_t)m->sample] += wm;
                    dacc[aux] += wm;
                    t_row[nt] = (int)aux;
                    t_col[nt] = m->sample;
                    t_val[nt] = -wm;
                    nt++;
                }
            } else {
                /* apply_A treats a dead run's mean as 0: pure diagonal */
                for (int32_t j = 0; j < run->count; j++) {
                    const RibStripMember *m =
                        &R->member[run->first + (size_t)j];
                    double wm = w * m->base_weight * m->like;
                    if (wm > 0.0) dacc[(size_t)m->sample] += wm;
                }
            }
        }
    }
    for (size_t i = 0; i < n; i++) dacc[i] += ridge;
    for (size_t p = 0; p < cx->n_pins; p++)
        dacc[(size_t)cx->pin_idx[p]] +=
            cx->pin_weight != NULL ? cx->pin_weight[p] : cx->pin_w;
    for (size_t i = 0; i < n_aug; i++) {
        double d = dacc[i];
        if (!(d > 1e-12)) d = 1.0; /* decouple dead aux rows, keep SPD */
        t_row[nt] = (int)i;
        t_col[nt] = (int)i;
        t_val[nt] = d;
        nt++;
    }
    assert(nt <= max_nt);
    rc = Sparse_factor_spd((int)n_aug, (int)nt, t_row, t_col, t_val,
                           &cx->chol);
    if (rc == 0 && cx->chol != NULL && sk > 0) {
        /* Dense SMW capacitance C = W^-1 + U^T M^-1 U, Cholesky-factored in
         * place.  One batched multi-rhs solve pass builds U^T M^-1 U; on any
         * failure the correction is dropped and the factored core alone
         * preconditions (still valid, just slower). */
        enum { RIB_SMW_BATCH = 32 };
        double *Bb = (double *)malloc(
            n_aug * (size_t)RIB_SMW_BATCH * sizeof *Bb);
        double *Xb = (double *)malloc(
            n_aug * (size_t)RIB_SMW_BATCH * sizeof *Xb);
        double *C = (double *)calloc(sk * sk, sizeof *C);
        double *ss = (double *)malloc((sk + 1) * sizeof *ss);
        int ok = Bb != NULL && Xb != NULL && C != NULL && ss != NULL;
        for (size_t j0 = 0; ok && j0 < sk; j0 += RIB_SMW_BATCH) {
            size_t jb = sk - j0 < (size_t)RIB_SMW_BATCH ? sk - j0
                                                        : RIB_SMW_BATCH;
            memset(Bb, 0, n_aug * jb * sizeof *Bb);
            for (size_t jj = 0; jj < jb; jj++)
                for (int32_t q = s_off[j0 + jj]; q < s_off[j0 + jj + 1]; q++)
                    Bb[jj * n_aug + (size_t)s_var[q]] = s_val[q];
            if (Sparse_factor_solve_multi(cx->chol, Bb, Xb, (int)jb) != 0) {
                ok = 0;
                break;
            }
            {
                int ci;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 8) num_threads(cx->threads)
#endif
                for (ci = 0; ci < (int)sk; ci++) {
                    for (size_t jj = 0; jj < jb; jj++) {
                        double dot = 0.0;
                        const double *xc = Xb + jj * n_aug;
                        for (int32_t q = s_off[ci]; q < s_off[ci + 1]; q++)
                            dot += s_val[q] * xc[(size_t)s_var[q]];
                        C[(size_t)ci * sk + j0 + jj] = dot;
                    }
                }
            }
        }
        if (ok) {
            for (size_t i = 0; i < sk; i++)
                C[i * sk + i] += 1.0 / s_w[i];
            ok = rib_dense_chol(C, sk) == 0;
        }
        if (ok) {
            cx->smw_k = sk;
            cx->smw_off = s_off;
            cx->smw_var = s_var;
            cx->smw_val = s_val;
            cx->smw_cchol = C;
            cx->smw_s = ss;
            s_off = NULL;   /* ownership moved to cx */
            s_var = NULL;
            s_val = NULL;
        } else {
            fprintf(stderr,
                    "  solve TAUCS: SMW setup failed for %zu long rows; "
                    "factored core only this round\n", sk);
            free(C);
            free(ss);
        }
        free(Bb);
        free(Xb);
    }
done:
    free(t_row);
    free(t_col);
    free(t_val);
    free(dacc);
    free(cl_var);
    free(cl_val);
    free(m_arr);
    free(s_off);
    free(s_var);
    free(s_val);
    free(s_w);
    if (rc != 0 || cx->chol == NULL) {
        cx->chol = NULL;
        return -1;
    }
    cx->chol_rhs = (double *)malloc(n_aug * sizeof *cx->chol_rhs);
    cx->chol_x = (double *)malloc(n_aug * sizeof *cx->chol_x);
    cx->chol_x2 = (double *)malloc(n_aug * sizeof *cx->chol_x2);
    if (cx->chol_rhs == NULL || cx->chol_x == NULL || cx->chol_x2 == NULL) {
        rib_strip_chol_release(cx);
        return -1;
    }
    cx->chol_n = n_aug;
    fprintf(stderr,
            "  solve TAUCS round %d: n=%zu (+%zu aux) nnz=%zu "
            "cliques<=%zu smw_rows=%zu spring_members=%zu factor_sec=%.1f\n",
            round + 1, n, naux, nt, run_cap, cx->smw_k, spring_members,
            ves_clock_sec() - t0);
    return 0;
}

/* Release the round's factorization, scratch, and SMW arrays. */
static void rib_strip_chol_release(SolveCtx *cx)
{
    Sparse_factor_free(&cx->chol);
    free(cx->smw_off);
    free(cx->smw_var);
    free(cx->smw_val);
    free(cx->smw_cchol);
    free(cx->smw_s);
    free(cx->chol_rhs);
    free(cx->chol_x);
    free(cx->chol_x2);
    cx->smw_off = NULL;
    cx->smw_var = NULL;
    cx->smw_val = NULL;
    cx->smw_cchol = NULL;
    cx->smw_s = NULL;
    cx->smw_k = 0;
    cx->chol_rhs = NULL;
    cx->chol_x = NULL;
    cx->chol_x2 = NULL;
    cx->chol_n = 0;
    cx->chol_attempted = 0;
}

static int cg_solve(Arena_T arena, SolveCtx *cx, double *u,
                     double *out_relres, RibGMGSnapshots *capture)
{
    double solve_t0 = ves_clock_sec();
    size_t n = cx->S->n_smp;
    int nt = cx->threads;
    int ii, cc;
    double cg_tol = cx->tol > 0.0 ? cx->tol : RIB_CG_TOL;
    if (cg_tol < RIB_STRIP_CG_TOL) cg_tol = RIB_STRIP_CG_TOL;
    Arena_Mark mark = Arena_save(arena);
    double *b    = RIB_ALLOC_ARRAY(arena, double, n);
    double *diag = RIB_ALLOC_ARRAY(arena, double, n);
    double *r    = RIB_ALLOC_ARRAY(arena, double, n);
    double *z    = RIB_ALLOC_ARRAY(arena, double, n);
    double *pv   = RIB_ALLOC_ARRAY(arena, double, n);
    double *Ap   = RIB_ALLOC_ARRAY(arena, double, n);
    double *best_u = RIB_ALLOC_ARRAY(arena, double, n);
    RibAMGWork amg_work;
    RibAMGWork *amg_wp = NULL;
    if (cx->amg != NULL && cx->amg->enabled) {
        int maxchain = 1;
        rib_amg_work_init(arena, cx->amg, &amg_work);
        for (size_t c = 0; c < cx->S->n_chn; c++)
            if (cx->S->chn[c].count > maxchain) maxchain = cx->S->chn[c].count;
        amg_work.maxchain = maxchain;
        amg_work.line_cp = RIB_ALLOC_ARRAY(
            arena, double, (size_t)maxchain * (size_t)nt + 1);
        amg_work.fine_tmp = RIB_ALLOC_ARRAY(arena, double, n + 1);
        amg_work.fine_corr = RIB_ALLOC_ARRAY(arena, double, n + 1);
        amg_wp = &amg_work;
    }
    build_b_diag(cx, b, diag);
    if (out_relres) *out_relres = 0.0;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nt)
#endif
    for (ii = 0; ii < (int)n; ii++)
        if (diag[ii] < 1e-12) diag[ii] = 1e-12;

    /* Test the actual incoming iterate before nested iteration touches it.
     * IRLS changes the weights only slightly on most rounds, so the previous
     * solution is often already converged.  More importantly, unsmoothed FMG
     * is an approximate correction on an irregular ribbon and is not itself a
     * descent method: it must never be allowed to make a good warm start worse. */
    double bnorm = 0.0, rn0 = 0.0;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:bnorm) num_threads(nt)
#endif
    for (ii = 0; ii < (int)n; ii++) bnorm += b[ii] * b[ii];
    bnorm = sqrt(bnorm);
    if (bnorm < 1e-30) { Arena_restore(arena, mark); return 0; }
    apply_A(cx, u, Ap);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:rn0) num_threads(nt)
#endif
    for (ii = 0; ii < (int)n; ii++) {
        r[ii] = b[ii] - Ap[ii];
        rn0 += r[ii] * r[ii];
    }
    double relres = sqrt(rn0) / bnorm;
    double best_relres = relres;
    memcpy(best_u, u, n * sizeof(*best_u));
    if (relres < cg_tol) {
        if (out_relres) *out_relres = relres;
        Arena_restore(arena, mark);
        return 0;
    }

    /* Lazy TAUCS factor: only a solve that actually iterates pays for the
     * round's factorization (post-IRLS rounds often pass the warm-iterate
     * check above and need none). */
    if (cx->strip != NULL && cx->chol == NULL && !cx->chol_attempted) {
        cx->chol_attempted = 1;
        if (rib_strip_chol_build(cx, RIB_STRIP_ADMM_RHO, cx->round_no) != 0)
            fprintf(stderr,
                    "  solve TAUCS: factor unavailable; Jacobi-PCG fallback\n");
    }

    if (amg_wp != NULL) {
        double fmg_t0 = ves_clock_sec();
        double pre_fmg_relres = relres;
        memcpy(pv, u, n * sizeof(*pv));
        /* Restrict the defect of the current metric parameterization.  The
         * former absolute-FMG path replaced u by a piecewise-constant block
         * field; on the production ribbon PCG hit its iteration cap before it
         * could reconstruct the lost arc-length ramp. */
        memset(amg_work.r[0], 0,
               (size_t)cx->amg->level[0].n * sizeof(double));
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 8) num_threads(nt)
#endif
        for (cc = 0; cc < (int)cx->S->n_chn; cc++) {
            int32_t f = cx->S->chn[cc].first, cn = cx->S->chn[cc].count;
            int32_t base = cx->amg->chain_base[cc];
            for (int32_t k0 = 0, block = 0;
                 k0 < cn; k0 += cx->amg->block_size, block++) {
                int32_t k1 = k0 + cx->amg->block_size;
                double s = 0.0;
                if (k1 > cn) k1 = cn;
                for (int32_t k = k0; k < k1; k++) s += r[f + k];
                amg_work.r[0][base + block] = s;
            }
        }
        if (rib_amg_fmg_solve(cx->amg, cx->S, &amg_work,
                              amg_work.r[0], amg_work.x[0], u,
                              capture) != 0) {
            Arena_restore(arena, mark);
            return -1;
        }
        /* Prolong and add the coarse correction to every ribbon sample. */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 8) num_threads(nt)
#endif
        for (cc = 0; cc < (int)cx->S->n_chn; cc++) {
            int32_t f = cx->S->chn[cc].first, cn = cx->S->chn[cc].count;
            int32_t base = cx->amg->chain_base[cc];
            for (int32_t k = 0; k < cn; k++)
                u[f + k] += amg_work.x[0][base + k / cx->amg->block_size];
        }
        /* The full sample grid is the final FMG level.  A symmetric line/GMG
         * correction restores within-block U detail before Krylov checks the
         * true matrix-free residual. */
        for (int cycle = 0; cycle < RIB_GMG_FINE_CYCLES; cycle++) {
            apply_A(cx, u, Ap);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nt)
#endif
            for (ii = 0; ii < (int)n; ii++) r[ii] = b[ii] - Ap[ii];
            rib_apply_preconditioner(cx, diag, r, z, amg_wp);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nt)
#endif
            for (ii = 0; ii < (int)n; ii++) u[ii] += z[ii];
        }
        double fmg_relres = rib_true_residual(cx, b, u, Ap, r, bnorm);
        if (!(fmg_relres < pre_fmg_relres)) {
            memcpy(u, pv, n * sizeof(*u));
            relres = rib_true_residual(cx, b, u, Ap, r, bnorm);
            if (capture != NULL) rib_gmg_snapshots_dispose(capture);
            fprintf(stderr,
                    "    GMG initial FMG correction REJECTED: relres %.3g -> %.3g "
                    "(%.2fs)\n",
                    pre_fmg_relres, fmg_relres, ves_clock_sec() - fmg_t0);
        } else {
            relres = fmg_relres;
            best_relres = relres;
            memcpy(best_u, u, n * sizeof(*best_u));
            fprintf(stderr,
                    "    GMG initial FMG correction accepted: relres %.3g -> %.3g "
                    "(%.2fs)\n",
                    pre_fmg_relres, relres, ves_clock_sec() - fmg_t0);
        }
    }

    if (relres < cg_tol) {
        if (out_relres) *out_relres = relres;
        Arena_restore(arena, mark);
        return 0;
    }
    double rz = 0.0;
    rib_apply_preconditioner(cx, diag, r, z, amg_wp);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:rz) num_threads(nt)
#endif
    for (ii = 0; ii < (int)n; ii++) {
        pv[ii] = z[ii];
        rz += r[ii] * z[ii];
    }
    int it = 0;
    int residual_reliable = 1;
    double last_reliable_relres = relres;
    for (; it < RIB_CG_MAXIT; it++) {
        double rn = 0.0;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:rn) num_threads(nt)
#endif
        for (ii = 0; ii < (int)n; ii++) rn += r[ii] * r[ii];
        relres = sqrt(rn) / bnorm;
        int restart_direction = 0;
        if (relres < cg_tol && !residual_reliable) {
            relres = rib_true_residual(cx, b, u, Ap, r, bnorm);
            residual_reliable = 1;
            restart_direction = relres >= cg_tol;
        }
        if (it == 0 || (it % 10) == 0) {
            fprintf(stderr,
                    "    Krylov progress: iter=%d relres=%.3g elapsed=%.2fs\n",
                    it, relres, ves_clock_sec() - solve_t0);
        }
        if (relres < cg_tol) break;
        if (restart_direction) {
            rz = 0.0;
            rib_apply_preconditioner(cx, diag, r, z, amg_wp);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:rz) num_threads(nt)
#endif
            for (ii = 0; ii < (int)n; ii++) {
                pv[ii] = z[ii];
                rz += r[ii] * z[ii];
            }
        }
        apply_A(cx, pv, Ap);
        double pAp = 0.0;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:pAp) num_threads(nt)
#endif
        for (ii = 0; ii < (int)n; ii++) pAp += pv[ii] * Ap[ii];
        if (pAp <= 0.0) break;
        double alpha = rz / pAp;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nt)
#endif
        for (ii = 0; ii < (int)n; ii++) {
            u[ii] += alpha * pv[ii];
            r[ii] -= alpha * Ap[ii];
        }
        int reliable_replace = ((it + 1) % RIB_CG_RELIABLE_PERIOD) == 0;
        int restart_after_replace = 0;
        residual_reliable = 0;
        if (reliable_replace) {
            double recursive_rn = 0.0;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:recursive_rn) num_threads(nt)
#endif
            for (ii = 0; ii < (int)n; ii++) recursive_rn += r[ii] * r[ii];
            double recursive_relres = sqrt(recursive_rn) / bnorm;
            /* A reliable residual *check* must not silently replace PCG's
             * recursively updated residual.  Replacing r while retaining p
             * destroys conjugacy (and produced the observed sawtooth robust
             * rounds).  Form the true residual in z, compare the vectors, and
             * copy it into r only when we also restart the direction. */
            double true_relres = rib_true_residual(cx, b, u, Ap, z, bnorm);
            if (!isfinite(true_relres)) {
                fprintf(stderr,
                        "    Krylov BUG: non-finite true residual at iter=%d\n",
                        it + 1);
                if (out_relres) *out_relres = true_relres;
                Arena_restore(arena, mark);
                return -2;
            }
            double residual_delta_n = 0.0;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:residual_delta_n) num_threads(nt)
#endif
            for (ii = 0; ii < (int)n; ii++) {
                double d = z[ii] - r[ii];
                residual_delta_n += d * d;
            }
            double discrepancy = sqrt(residual_delta_n) / bnorm;
            if (true_relres < best_relres) {
                best_relres = true_relres;
                memcpy(best_u, u, n * sizeof(*best_u));
            }
            double scale = fmax(true_relres, RIB_STRIP_CG_TOL);
            restart_after_replace = true_relres < cg_tol ||
                discrepancy > 0.25 * scale;
            if (restart_after_replace) {
                memcpy(r, z, n * sizeof(*r));
                fprintf(stderr,
                        "    Krylov reliable restart: true=%.3g recursive=%.3g "
                        "delta=%.3g previous=%.3g\n",
                        true_relres, recursive_relres, discrepancy,
                        last_reliable_relres);
                residual_reliable = 1;
            }
            relres = true_relres;
            last_reliable_relres = true_relres;
        }
        double rz2 = 0.0;
        rib_apply_preconditioner(cx, diag, r, z, amg_wp);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:rz2) num_threads(nt)
#endif
        for (ii = 0; ii < (int)n; ii++) {
            rz2 += r[ii] * z[ii];
        }
        double beta = restart_after_replace ? 0.0 : rz2 / rz;
        rz = rz2;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nt)
#endif
        for (ii = 0; ii < (int)n; ii++)
            pv[ii] = z[ii] + beta * pv[ii];
    }
    relres = rib_true_residual(cx, b, u, Ap, r, bnorm);
    if (relres < best_relres) {
        best_relres = relres;
        memcpy(best_u, u, n * sizeof(*best_u));
    } else if (relres > best_relres * (1.0 + 1e-12)) {
        fprintf(stderr,
                "    Krylov safeguard: restored best verified iterate "
                "relres %.3g instead of final %.3g\n",
                best_relres, relres);
        memcpy(u, best_u, n * sizeof(*u));
        relres = best_relres;
    }
    if (out_relres) *out_relres = relres;
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
        while (nb > 0 && bm[nb-1] > m + RIB_STRIP_GEOM_EPS) {
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
            if (fabs(nu - u[i]) > RIB_STRIP_GEOM_EPS) moved++;
            u[i] = nu;
        }
    }
    return moved;
}

static size_t rib_project_monotone(const SliceSet *S, double *u,
                                   double *scratch, double *lower_bound,
                                   size_t *thread_moved, int maxchain,
                                   int threads)
{
    memset(thread_moved, 0, (size_t)threads * sizeof(*thread_moved));
    int cc;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 8) num_threads(threads)
#endif
    for (cc = 0; cc < (int)S->n_chn; cc++) {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        int32_t f = S->chn[cc].first, count = S->chn[cc].count;
        if (count < 2) continue;
        double *chain_scratch = scratch +
            (size_t)tid * (size_t)maxchain * 3u;
        double *chain_lb = lower_bound +
            (size_t)tid * (size_t)maxchain;
        chain_lb[0] = 0.0;
        for (int32_t j = 1; j < count; j++)
            chain_lb[j] = 0.5 *
                (S->smp[f + j].s - S->smp[f + j - 1].s);
        thread_moved[tid] +=
            pava_chain(&u[f], chain_lb, count, chain_scratch);
    }
    size_t moved = 0;
    for (int tid = 0; tid < threads; tid++) moved += thread_moved[tid];
    return moved;
}

static const char *rib_strip_round_phase(int round, int relaxed_rounds,
                                         int final_rounds)
{
    if (round >= relaxed_rounds + final_rounds) return "metric";
    if (round >= relaxed_rounds) return "final";
    if (round == 0) return "initial";
    return round <= 3 ? "l1" : "likelihood";
}

/* Publish the intended round sequence before the first long solve.  The live
 * baker can therefore wait for exactly N atomic VMESH publications and turn
 * each one into a PNG while the next solve round is running. */
static int rib_write_solve_round_manifest(const RibbonOpts *o,
                                          int total_rounds,
                                          int relaxed_rounds,
                                          int final_rounds)
{
    char path[2600];
    if (o->solve_level_prefix == NULL) return 0;
    if (snprintf(path, sizeof path, "%s_solve_rounds.json",
                 o->solve_level_prefix) < 0 ||
        strlen(path) + 1 >= sizeof path) return -1;
    if (ves_ensure_parent_dir(path) != 0) return -1;
    FILE *fp = fopen(path, "wb");
    if (fp == NULL) return -1;
    fprintf(fp,
            "{\n  \"schema\": \"ribbon-solve-rounds-v1\",\n"
            "  \"order\": \"solve-order\",\n"
            "  \"count\": %d,\n  \"rounds\": [\n",
            total_rounds);
    for (int round = 0; round < total_rounds; round++) {
        fprintf(fp,
                "    { \"ordinal\": %d, \"phase\": \"%s\", "
                "\"vmesh\": \"%s_solve_round_%02d.vmesh\" }%s\n",
                round, rib_strip_round_phase(
                           round, relaxed_rounds, final_rounds),
                ves_path_basename(o->solve_level_prefix), round,
                round + 1 < total_rounds ? "," : "");
    }
    fprintf(fp, "  ]\n}\n");
    if (fflush(fp) != 0 || fclose(fp) != 0) return -1;
    return 0;
}

/* piecewise phi -> u map built from the MAIN solve component: sorted non-empty
 * phi bins with median u; lookup lerps, ends extrapolate with the edge slope */
typedef struct { double *bphi, *bu; int nb; double end_slope; } PhiUMap;

static double phiu_lookup(const PhiUMap *M, double phi)
{
    if (M->nb == 0) return 0.0;
    if (M->nb == 1) return M->bu[0] + (phi - M->bphi[0]) * M->end_slope;
    if (phi <= M->bphi[0]) {
        double slope = (M->bu[1] - M->bu[0]) /
                       (M->bphi[1] - M->bphi[0]);
        if (!(slope > 1e-9) || !isfinite(slope)) slope = M->end_slope;
        return M->bu[0] - (M->bphi[0] - phi) * slope;
    }
    if (phi >= M->bphi[M->nb-1]) {
        double slope = (M->bu[M->nb-1] - M->bu[M->nb-2]) /
                       (M->bphi[M->nb-1] - M->bphi[M->nb-2]);
        if (!(slope > 1e-9) || !isfinite(slope)) slope = M->end_slope;
        return M->bu[M->nb-1] + (phi - M->bphi[M->nb-1]) * slope;
    }
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

typedef struct {
    int32_t a, b;
    int32_t slice;
    double delta;  /* shift[b] - shift[a] */
} RibGaugeObs;

typedef struct {
    int32_t a, b;
    int32_t slice_lo, slice_hi;
    size_t support;
    double delta, spread80;
} RibGaugeRel;

typedef struct {
    int32_t component, bin;
    double u;
} RibPhaseGaugeSample;

typedef struct {
    int32_t bin;
    double u;
} RibPhaseGaugeBin;

typedef struct {
    int32_t component;
    size_t samples;
    double phi_lo, phi_hi, u_lo, u_hi;
    int32_t slice_lo, slice_hi;
} RibGaugeComponent;

static int cmp_rib_gauge_obs(const void *pa, const void *pb)
{
    const RibGaugeObs *a = (const RibGaugeObs *)pa;
    const RibGaugeObs *b = (const RibGaugeObs *)pb;
    if (a->a != b->a) return a->a < b->a ? -1 : 1;
    if (a->b != b->b) return a->b < b->b ? -1 : 1;
    return a->delta < b->delta ? -1 : (a->delta > b->delta ? 1 : 0);
}

static int cmp_rib_gauge_rel(const void *pa, const void *pb)
{
    const RibGaugeRel *a = (const RibGaugeRel *)pa;
    const RibGaugeRel *b = (const RibGaugeRel *)pb;
    if (a->support != b->support) return a->support > b->support ? -1 : 1;
    if (a->spread80 != b->spread80)
        return a->spread80 < b->spread80 ? -1 : 1;
    if (a->a != b->a) return a->a < b->a ? -1 : 1;
    return a->b < b->b ? -1 : (a->b > b->b ? 1 : 0);
}

static int cmp_rib_phase_gauge_sample(const void *pa, const void *pb)
{
    const RibPhaseGaugeSample *a = (const RibPhaseGaugeSample *)pa;
    const RibPhaseGaugeSample *b = (const RibPhaseGaugeSample *)pb;
    if (a->component != b->component)
        return a->component < b->component ? -1 : 1;
    if (a->bin != b->bin) return a->bin < b->bin ? -1 : 1;
    return a->u < b->u ? -1 : (a->u > b->u ? 1 : 0);
}

static int cmp_rib_gauge_component(const void *pa, const void *pb)
{
    const RibGaugeComponent *a = (const RibGaugeComponent *)pa;
    const RibGaugeComponent *b = (const RibGaugeComponent *)pb;
    if (a->phi_lo != b->phi_lo) return a->phi_lo < b->phi_lo ? -1 : 1;
    if (a->phi_hi != b->phi_hi) return a->phi_hi < b->phi_hi ? -1 : 1;
    return a->component < b->component ? -1 :
           (a->component > b->component ? 1 : 0);
}

/* A mutual-nearest restriction is required to turn adjacent slice samples into
 * non-branching StrokeStrip runs.  The candidates it rejects are nevertheless
 * measured adjacent-slice observations.  Audit whether those observations
 * consistently reconnect additive solve gauges; this function deliberately
 * changes no coordinate. */
static void rib_strip_audit_pruned_gauges(
    Arena_T arena, const SliceSet *S, const PairSet *P,
    const RibStripSet *R, const double *u, int32_t ncomp)
{
    if (ncomp <= 1 || R->npruned == 0) return;
    Arena_Mark mark = Arena_save(arena);
    RibGaugeObs *obs = RIB_ALLOC_ARRAY(arena, RibGaugeObs, R->npruned);
    size_t no = 0;
    for (size_t q = 0; q < R->npruned; q++) {
        const Pair *p = &P->pairs[R->pruned_pair[q]];
        int32_t cha = S->smp[p->a].chain, chb = S->smp[p->b].chain;
        if (cha < 0 || chb < 0) continue;
        int32_t ca = S->chn[cha].solve_comp;
        int32_t cb = S->chn[chb].solve_comp;
        if (ca < 0 || cb < 0 || ca == cb) continue;
        double delta = (u[p->a] - u[p->b]) - p->d;
        if (!isfinite(delta)) continue;
        if (ca < cb) {
            obs[no].a = ca; obs[no].b = cb; obs[no].delta = delta;
        } else {
            obs[no].a = cb; obs[no].b = ca; obs[no].delta = -delta;
        }
        obs[no].slice = S->smp[p->a].slice;
        no++;
    }
    if (no == 0) {
        fprintf(stderr,
                "  pruned-gauge audit: no cross-component observations\n");
        Arena_restore(arena, mark);
        return;
    }
    qsort(obs, no, sizeof(*obs), cmp_rib_gauge_obs);
    RibGaugeRel *rel = RIB_ALLOC_ARRAY(arena, RibGaugeRel, no);
    size_t nr = 0;
    for (size_t i = 0; i < no; ) {
        size_t j = i + 1;
        int32_t slo = obs[i].slice, shi = obs[i].slice;
        while (j < no && obs[j].a == obs[i].a && obs[j].b == obs[i].b) {
            if (obs[j].slice < slo) slo = obs[j].slice;
            if (obs[j].slice > shi) shi = obs[j].slice;
            j++;
        }
        size_t count = j - i;
        rel[nr].a = obs[i].a;
        rel[nr].b = obs[i].b;
        rel[nr].slice_lo = slo;
        rel[nr].slice_hi = shi;
        rel[nr].support = count;
        rel[nr].delta = obs[i + count / 2].delta;
        rel[nr].spread80 = obs[i + (9 * count) / 10].delta -
                           obs[i + count / 10].delta;
        nr++;
        i = j;
    }
    qsort(rel, nr, sizeof(*rel), cmp_rib_gauge_rel);

    UnionFind all = UF_new(arena, ncomp);
    UnionFind stable = UF_new(arena, ncomp);
    size_t nstable = 0;
    for (size_t i = 0; i < nr; i++) {
        uf_union(&all, rel[i].a, rel[i].b);
        if (rel[i].support >= 4 &&
            rel[i].slice_hi - rel[i].slice_lo >= 2 &&
            rel[i].spread80 <= 16.0) {
            uf_union(&stable, rel[i].a, rel[i].b);
            nstable++;
        }
    }
    int all_components = 0, stable_components = 0;
    for (int32_t c = 0; c < ncomp; c++) {
        if (uf_find(&all, c) == c) all_components++;
        if (uf_find(&stable, c) == c) stable_components++;
    }
    fprintf(stderr,
            "  pruned-gauge audit: observations=%zu relations=%zu; "
            "all graph=%d component(s); stable=%zu relation(s) -> %d component(s)\n",
            no, nr, all_components, nstable, stable_components);
    size_t report = nr < 64 ? nr : 64;
    for (size_t i = 0; i < report; i++)
        fprintf(stderr,
                "    gauge %d-%d support=%zu slices=%d..%d "
                "delta=%.3f spread80=%.3f%s\n",
                rel[i].a, rel[i].b, rel[i].support,
                rel[i].slice_lo, rel[i].slice_hi,
                rel[i].delta, rel[i].spread80,
                rel[i].support >= 4 &&
                rel[i].slice_hi - rel[i].slice_lo >= 2 &&
                rel[i].spread80 <= 16.0 ? " stable" : "");
    Arena_restore(arena, mark);
}

/* The lifted phase is a material-identity scaffold, not a metric coordinate.
 * It may nevertheless identify two disconnected solve components observing the
 * same material interval.  Compare their independently solved U curves only at
 * common narrow phase bins and report the component-constant gauge relation.
 * No phase derivative enters the metric and this audit changes no coordinate. */
static void rib_strip_audit_phase_gauges(
    Arena_T arena, const SliceSet *S, const double *u, int32_t ncomp)
{
    const double bin_width = 0.02;
    if (ncomp <= 1) return;
    /* This is a read-only O(C^2) diagnostic, not part of the solve.  Torn
     * quadribbons can legitimately contain hundreds of thousands of tiny
     * metric components; allocating every possible component pair then costs
     * terabytes and used to abort an otherwise completed solve. */
    if (ncomp > 4096) {
        fprintf(stderr,
                "  phase-gauge audit: skipped pairwise diagnostic for %d "
                "components (safe cap 4096)\n",
                ncomp);
        return;
    }
    Arena_Mark mark = Arena_save(arena);
    RibGaugeComponent *component = RIB_ALLOC_ARRAY(
        arena, RibGaugeComponent, (size_t)ncomp);
    for (int32_t c = 0; c < ncomp; c++) {
        component[c].component = c;
        component[c].samples = 0;
        component[c].phi_lo = component[c].u_lo = 1e300;
        component[c].phi_hi = component[c].u_hi = -1e300;
        component[c].slice_lo = INT32_MAX;
        component[c].slice_hi = INT32_MIN;
    }
    RibPhaseGaugeSample *sample = RIB_ALLOC_ARRAY(
        arena, RibPhaseGaugeSample, S->n_smp);
    size_t ns = 0;
    for (size_t i = 0; i < S->n_smp; i++) {
        int32_t ch = S->smp[i].chain;
        if (ch < 0 || !isfinite(S->smp[i].phi) || !isfinite(u[i])) continue;
        double raw_bin = floor(S->smp[i].phi / bin_width);
        if (raw_bin < (double)INT32_MIN || raw_bin > (double)INT32_MAX)
            continue;
        sample[ns].component = S->chn[ch].solve_comp;
        sample[ns].bin = (int32_t)raw_bin;
        sample[ns].u = u[i];
        {
            RibGaugeComponent *g = &component[sample[ns].component];
            double phi = S->smp[i].phi;
            g->samples++;
            if (phi < g->phi_lo) g->phi_lo = phi;
            if (phi > g->phi_hi) g->phi_hi = phi;
            if (u[i] < g->u_lo) g->u_lo = u[i];
            if (u[i] > g->u_hi) g->u_hi = u[i];
            if (S->smp[i].slice < g->slice_lo)
                g->slice_lo = S->smp[i].slice;
            if (S->smp[i].slice > g->slice_hi)
                g->slice_hi = S->smp[i].slice;
        }
        ns++;
    }
    qsort(sample, ns, sizeof(*sample), cmp_rib_phase_gauge_sample);
    size_t *off = (size_t *)ARENA_CALLOC(
        arena, (size_t)ncomp + 1, sizeof(*off));
    RibPhaseGaugeBin *curve = RIB_ALLOC_ARRAY(
        arena, RibPhaseGaugeBin, ns ? ns : 1);
    size_t nc = 0;
    for (size_t i = 0; i < ns; ) {
        size_t j = i + 1;
        while (j < ns &&
               sample[j].component == sample[i].component &&
               sample[j].bin == sample[i].bin)
            j++;
        int32_t comp = sample[i].component;
        if (comp >= 0 && comp < ncomp) {
            curve[nc].bin = sample[i].bin;
            curve[nc].u = sample[i + (j - i) / 2].u;
            off[(size_t)comp + 1]++;
            nc++;
        }
        i = j;
    }
    for (int32_t c = 0; c < ncomp; c++) off[(size_t)c + 1] += off[c];

    qsort(component, (size_t)ncomp, sizeof(*component),
          cmp_rib_gauge_component);
    double sum_metric_width = 0.0;
    fprintf(stderr, "  solve-component intervals in lifted-phase order:\n");
    for (int32_t k = 0; k < ncomp; k++) {
        const RibGaugeComponent *g = &component[k];
        double width = g->u_hi - g->u_lo;
        sum_metric_width += width;
        fprintf(stderr,
                "    solve %d samples=%zu slices=%d..%d "
                "phi=[%.3f,%.3f] U=[%.1f,%.1f] width=%.1f\n",
                g->component, g->samples, g->slice_lo, g->slice_hi,
                g->phi_lo, g->phi_hi, g->u_lo, g->u_hi, width);
    }
    fprintf(stderr,
            "    sum of independently measured component widths: %.1f vox\n",
            sum_metric_width);

    size_t max_bins = 1;
    for (int32_t c = 0; c < ncomp; c++)
        if (off[(size_t)c + 1] - off[c] > max_bins)
            max_bins = off[(size_t)c + 1] - off[c];
    double *delta = RIB_ALLOC_ARRAY(arena, double, max_bins);
    size_t max_rel = (size_t)ncomp * (size_t)(ncomp - 1) / 2;
    RibGaugeRel *rel = RIB_ALLOC_ARRAY(arena, RibGaugeRel, max_rel ? max_rel : 1);
    size_t nr = 0;
    for (int32_t a = 0; a < ncomp; a++) {
        for (int32_t b = a + 1; b < ncomp; b++) {
            size_t ia = off[a], ib = off[b], nd = 0;
            int32_t blo = INT32_MAX, bhi = INT32_MIN;
            while (ia < off[(size_t)a + 1] &&
                   ib < off[(size_t)b + 1]) {
                if (curve[ia].bin < curve[ib].bin) ia++;
                else if (curve[ib].bin < curve[ia].bin) ib++;
                else {
                    delta[nd++] = curve[ia].u - curve[ib].u;
                    if (curve[ia].bin < blo) blo = curve[ia].bin;
                    if (curve[ia].bin > bhi) bhi = curve[ia].bin;
                    ia++; ib++;
                }
            }
            if (nd == 0) continue;
            qsort(delta, nd, sizeof(*delta), cmp_dbl);
            rel[nr].a = a; rel[nr].b = b;
            rel[nr].slice_lo = blo; rel[nr].slice_hi = bhi;
            rel[nr].support = nd;
            rel[nr].delta = delta[nd / 2];
            rel[nr].spread80 = delta[(9 * nd) / 10] - delta[nd / 10];
            nr++;
        }
    }
    qsort(rel, nr, sizeof(*rel), cmp_rib_gauge_rel);
    UnionFind stable = UF_new(arena, ncomp);
    size_t nstable = 0;
    for (size_t i = 0; i < nr; i++) {
        int accepted = rel[i].support >= 8 &&
                       rel[i].slice_hi - rel[i].slice_lo >= 25 &&
                       rel[i].spread80 <= 16.0;
        if (!accepted) continue;
        uf_union(&stable, rel[i].a, rel[i].b);
        nstable++;
    }
    int stable_components = 0;
    for (int32_t c = 0; c < ncomp; c++)
        if (uf_find(&stable, c) == c) stable_components++;
    fprintf(stderr,
            "  phase-gauge audit: samples=%zu bins=%zu relations=%zu; "
            "stable=%zu relation(s) -> %d component(s)\n",
            ns, nc, nr, nstable, stable_components);
    size_t report = nr < 64 ? nr : 64;
    for (size_t i = 0; i < report; i++) {
        int accepted = rel[i].support >= 8 &&
                       rel[i].slice_hi - rel[i].slice_lo >= 25 &&
                       rel[i].spread80 <= 16.0;
        fprintf(stderr,
                "    phase gauge %d-%d support=%zu phi_span=%.3f "
                "delta=%.3f spread80=%.3f%s\n",
                rel[i].a, rel[i].b, rel[i].support,
                bin_width * (double)(rel[i].slice_hi - rel[i].slice_lo),
                rel[i].delta, rel[i].spread80,
                accepted ? " stable" : "");
    }
    Arena_restore(arena, mark);
}

/* Stage C determines every within-component U difference, but disconnected
 * metric components retain arbitrary additive constants.  The authoritative
 * lifted phase supplies their order on the universal cover.  Concatenate the
 * independently measured U intervals in that order.  This is a gauge choice,
 * not a deformation: every component receives one constant U shift; V and all
 * local derivatives are untouched. */
static void rib_strip_stitch_solve_gauges(
    Arena_T arena, SliceSet *S, double *u, int32_t ncomp,
    RibGMGSnapshots *level_capture, RibbonResult *out)
{
    if (ncomp <= 0) return;
    Arena_Mark mark = Arena_save(arena);
    RibGaugeComponent *component = RIB_ALLOC_ARRAY(
        arena, RibGaugeComponent, (size_t)ncomp);
    double *shift = RIB_ALLOC_ARRAY(arena, double, (size_t)ncomp);
    for (int32_t c = 0; c < ncomp; c++) {
        component[c].component = c;
        component[c].samples = 0;
        component[c].phi_lo = component[c].u_lo = 1e300;
        component[c].phi_hi = component[c].u_hi = -1e300;
        component[c].slice_lo = INT32_MAX;
        component[c].slice_hi = INT32_MIN;
        shift[c] = 0.0;
    }
    for (size_t i = 0; i < S->n_smp; i++) {
        int32_t ch = S->smp[i].chain;
        if (ch < 0) continue;
        int32_t c = S->chn[ch].solve_comp;
        if (c < 0 || c >= ncomp || !isfinite(S->smp[i].phi) ||
            !isfinite(u[i]))
            continue;
        RibGaugeComponent *g = &component[c];
        g->samples++;
        if (S->smp[i].phi < g->phi_lo) g->phi_lo = S->smp[i].phi;
        if (S->smp[i].phi > g->phi_hi) g->phi_hi = S->smp[i].phi;
        if (u[i] < g->u_lo) g->u_lo = u[i];
        if (u[i] > g->u_hi) g->u_hi = u[i];
    }
    qsort(component, (size_t)ncomp, sizeof(*component),
          cmp_rib_gauge_component);
    double cursor = 0.0, max_shift = 0.0;
    size_t active = 0;
    for (int32_t k = 0; k < ncomp; k++) {
        RibGaugeComponent *g = &component[k];
        if (g->samples == 0 || !(g->u_hi >= g->u_lo)) continue;
        double ds = cursor - g->u_lo;
        shift[g->component] = ds;
        if (fabs(ds) > max_shift) max_shift = fabs(ds);
        cursor += g->u_hi - g->u_lo;
        active++;
    }
    for (size_t i = 0; i < S->n_smp; i++) {
        int32_t ch = S->smp[i].chain;
        if (ch < 0) continue;
        int32_t c = S->chn[ch].solve_comp;
        if (c < 0 || c >= ncomp) continue;
        u[i] += shift[c];
        if (level_capture != NULL)
            for (int l = 0; l < level_capture->count; l++)
                level_capture->snapshot[l].u[i] += (float)shift[c];
    }
    out->strip_gauge_components = active;
    out->strip_gauge_span = cursor;
    out->strip_gauge_max_shift = max_shift;
    fprintf(stderr,
            "  phase-ordered gauge stitch: %zu component interval(s), "
            "U span=%.1f, max additive shift=%.1f; V unchanged\n",
            active, cursor, max_shift);
    Arena_restore(arena, mark);
}

static void solve_parameterization(Arena_T arena, SliceSet *S, PairSet *P,
                                   double r_ref, const RibbonOpts *o,
                                   RibbonResult *out,
                                   RibGMGSnapshots *level_capture)
{
    size_t n = S->n_smp;
    size_t nc = S->n_chn;
    Arena_Mark mark = Arena_save(arena);

    /* Numerical gauge initialization.  A monotone cover-space integral
     *   u(phi) = a*phi + b*phi^2/(4pi),   r(phi) = a + b*phi/(2pi),
     * places every disconnected metric solve on one injective universal-cover
     * scaffold.  It supplies ONLY one additive pin per solve component: the
     * LS--PAVA rows below still determine every local U difference from measured
     * arclength.  Thus a weak radius fit can be a perfectly valid gauge while
     * remaining a poor material-coordinate model.  The relevant invariant is
     * monotonicity over the observed lifted-phase interval, not radius-fit R2.
     * Without a safe cover gauge, use phi*r_ref as an initializer and register
     * detached components against the measured main chart. */
    double sp_a = out->spiral_a, sp_b = out->spiral_b;
    double phi_lo = 1e300, phi_hi = -1e300;
    for (size_t i = 0; i < n; i++) {
        if (S->smp[i].chain < 0) continue;
        if (S->smp[i].phi < phi_lo) phi_lo = S->smp[i].phi;
        if (S->smp[i].phi > phi_hi) phi_hi = S->smp[i].phi;
    }
    double cover_r0 = sp_a + sp_b * phi_lo / (2.0 * M_PI);
    double cover_r1 = sp_a + sp_b * phi_hi / (2.0 * M_PI);
    int cover_monotone = isfinite(cover_r0) && isfinite(cover_r1) &&
                         cover_r0 > RIB_PRIOR_MIN_B &&
                         cover_r1 > RIB_PRIOR_MIN_B;
    /* A caller-pinned line is an explicit shared-frame contract.  With an
     * authoritative lift, any monotone line is a safe gauge even when R2 is
     * poor.  The legacy slice-only path still demands a genuinely spiral-like
     * fit because its phi itself is less trustworthy. */
    int spiral_ok = (o->spiral_b != 0.0)
                 || (o->reference_phi != NULL &&
                     cover_monotone && fabs(sp_b) > RIB_PRIOR_MIN_B)
                 || (o->reference_phi == NULL &&
                     out->spiral_r2 >= RIB_PRIOR_MIN_R2 &&
                     fabs(sp_b) > RIB_PRIOR_MIN_B);
    int use_reference_seed = o->reference_u != NULL && o->solve_reference_u;
    double *u = RIB_ALLOC_ARRAY(arena, double, n);
    for (size_t i = 0; i < n; i++) {
        double phi = S->smp[i].phi;
        u[i] = use_reference_seed
             ? S->smp[i].u
             : spiral_ok ? (sp_a * phi + sp_b * phi * phi / (4.0 * M_PI))
                         : phi * r_ref;
        if (!isfinite(u[i])) {
            fprintf(stderr,
                    "ribbon: non-finite Stage-C U initializer at sample %zu\n",
                    i);
            RAISE(Arena_Failed);
        }
    }
    if (use_reference_seed)
        fprintf(stderr,
                "  winding-certificate initializer: carried U enters the "
                "global solve unchanged; no per-stroke gauge projection\n");

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

    /* rib_strip_build_runs marks non-selected cross-slice alternatives off in
     * Pair.like.  Keep the post-geometry/post-winding candidate support so the
     * final StrokeStrip step can rebuild C(u) from the relaxed isovalues rather
     * than being trapped in the initial mutual-nearest topology. */
    double *cross_pair_support = RIB_ALLOC_ARRAY(
        arena, double, P->n_cross ? P->n_cross : 1);
    for (size_t p = 0; p < P->n_cross; p++)
        cross_pair_support[p] = P->pairs[p].like;

    RibStripSet strip;
    if (rib_strip_build_runs(arena, S, P, use_reference_seed,
                             RIB_STRIP_INITIAL_LIKE, &strip) != 0) {
        fprintf(stderr, "ribbon: cannot construct V-connected StrokeStrip runs\n");
        RAISE(Arena_Failed);
    }
    rib_strip_update_length_coefficients(S, &strip);
    if (use_reference_seed)
        rib_strip_audit_reference_seed(arena, S, P, &strip);
    out->n_strip_runs = strip.nrun;
    out->n_strip_members = strip.nmember;
    out->n_strip_links = strip.nlink;
    out->n_strip_links_pruned = strip.npruned;

    /* Connected components of the metric system: selected V-run links plus
     * physical same-slice continuations.  Rejected branch alternatives do not
     * share a gauge merely because they were spatially close. */
    UnionFind uf = UF_new(arena, (int32_t)(nc ? nc : 1));
    for (size_t q = 0; q < strip.nlink; q++) {
        const Pair *pr = &P->pairs[strip.link_pair[q]];
        int32_t ca = S->smp[pr->a].chain, cb = S->smp[pr->b].chain;
        if (ca >= 0 && cb >= 0 && ca != cb) uf_union(&uf, ca, cb);
    }
    for (size_t p = P->n_cross; p < P->n_pairs; p++) {
        const Pair *pr = &P->pairs[p];
        if (pr->like < 0.0) continue;
        int32_t ca = S->smp[pr->a].chain, cb = S->smp[pr->b].chain;
        if (ca >= 0 && cb >= 0 && ca != cb)
            uf_union(&uf, ca, cb);
    }
    /* longest chain per component root -> pin its first sample at its init */
    int32_t *rootbest = RIB_ALLOC_ARRAY(arena, int32_t, nc ? nc : 1);
    for (size_t c = 0; c < nc; c++) rootbest[c] = -1;
    int32_t bigchain = 0, bestlen = 0;
    for (size_t c = 0; c < nc; c++) {
        int32_t r = uf_find(&uf, (int32_t)c);
        if (rootbest[r] < 0 || S->chn[c].count > S->chn[rootbest[r]].count)
            rootbest[r] = (int32_t)c;
        if (S->chn[c].count > bestlen) { bestlen = S->chn[c].count; bigchain = (int32_t)c; }
    }
    int32_t *pin_idx = RIB_ALLOC_ARRAY(arena, int32_t, nc ? nc : 1);
    double  *pin_val = RIB_ALLOC_ARRAY(arena, double, nc ? nc : 1);
    size_t n_pins = 0;
    for (size_t c = 0; c < nc; c++) {
        if (uf_find(&uf, (int32_t)c) != (int32_t)c) continue;   /* roots only */
        int32_t pc = rootbest[c];
        if (pc < 0) continue;
        pin_idx[n_pins] = S->chn[pc].first;
        pin_val[n_pins] = u[S->chn[pc].first];
        n_pins++;
    }
    /* Persist a compact solve-component label through Stage D.  The original
     * mesh can bridge sampling holes that disconnect this solve graph; transfer
     * uses those exact chart-continuity constraints to reconcile only the
     * otherwise free additive gauges. */
    int32_t *root_solve = RIB_ALLOC_ARRAY(arena, int32_t, nc ? nc : 1);
    for (size_t c = 0; c < nc; c++) root_solve[c] = -1;
    int32_t nsolve = 0;
    for (size_t c = 0; c < nc; c++) {
        int32_t root = uf_find(&uf, (int32_t)c);
        if (root_solve[root] < 0) root_solve[root] = nsolve++;
        S->chn[c].solve_comp = root_solve[root];
    }
    out->n_qp_comps = (int)nsolve;
    SolveGraph graph;
    double graph_t0 = ves_clock_sec();
    if (solve_graph_build(arena, S, P, P->n_cross, &graph) != 0) {
        fprintf(stderr, "ribbon: solve graph exceeds the supported index range\n");
        RAISE(Arena_Failed);
    }
    double graph_sec = ves_clock_sec() - graph_t0;
    int solve_threads = solve_thread_count(o->solve_threads);
    {
        int useful = (int)((graph.n + 32767u) / 32768u);
        if (useful < 1) useful = 1;
        if (solve_threads > useful) solve_threads = useful;
    }
    /* A rebuilt final C(u) can contain a different number of runs.  n samples
     * is a strict upper bound, and these two arrays are the only O(n) scratch
     * needed to let the topology change without reallocating the solver. */
    double *strip_dot = RIB_ALLOC_ARRAY(arena, double, n + 1);
    double *strip_mean = RIB_ALLOC_ARRAY(arena, double, n + 1);
    SolveCtx cx;
    memset(&cx, 0, sizeof cx);
    cx.S = S;
    cx.P = P;
    cx.G = &graph;
    cx.strip = &strip;
    cx.amg = NULL;
    cx.pin_idx = pin_idx;
    cx.pin_val = pin_val;
    cx.n_pins = n_pins;
    cx.pin_w = 1.0;
    cx.strip_dot = strip_dot;
    cx.strip_mean = strip_mean;
    cx.threads = solve_threads;
    /* The samples themselves are spaced at several voxels.  A 1e-6 relative
     * residual is already below the meaningful coordinate scale; never turn a
     * mesh-processing check into a demand for near-machine-precision CG. */
    cx.tol = RIB_STRIP_CG_TOL;
    fprintf(stderr,
            "  StrokeStrip rows: runs=%zu members=%zu links=%zu "
            "branch_links_pruned=%zu; graph nodes=%zu direct_edges=%zu "
            "continuations=%zu threads=%d build_sec=%.2f\n",
            strip.nrun, strip.nmember, strip.nlink, strip.npruned,
            graph.n, graph.n_adj, graph.n_active_pair, solve_threads, graph_sec);

    int32_t maxchain = 1;
    for (size_t c = 0; c < S->n_chn; c++)
        if (S->chn[c].count > maxchain) maxchain = S->chn[c].count;
    size_t thread_chain_n = (size_t)maxchain * (size_t)solve_threads;
    double *scr = RIB_ALLOC_ARRAY(arena, double, thread_chain_n * 3u + 8u);
    double *lb  = RIB_ALLOC_ARRAY(arena, double, thread_chain_n + 1u);
    size_t *thread_moved = (size_t *)ARENA_CALLOC(
        arena, (size_t)solve_threads, sizeof(size_t));
    int reweight_rounds = o->relax_iters > 0 ? o->relax_iters : 0;
    int relaxed_rounds = 1 + reweight_rounds;
    int final_rounds = o->final_iters > 0 ? o->final_iters : 0;
    int metric_rounds = o->metric_iters > 0 ? o->metric_iters : 0;
    int total_rounds = relaxed_rounds + final_rounds + metric_rounds;
    if (rib_write_solve_round_manifest(
            o, total_rounds, relaxed_rounds, final_rounds) != 0) {
        fprintf(stderr, "ribbon: cannot publish solve-round manifest\n");
        RAISE(Arena_Failed);
    }
    double longest = 0.0;
    for (size_t c = 0; c < S->n_chn; c++) {
        int32_t f = S->chn[c].first, cn = S->chn[c].count;
        if (cn > 1) {
            double span = S->smp[f + cn - 1].s - S->smp[f].s;
            if (span > longest) longest = span;
        }
    }
    double strip_sigma = fmax(longest / 30.0, 1e-3);
    size_t mono_total = 0;
    double *admm_z = RIB_ALLOC_ARRAY(arena, double, n);
    double *admm_y = RIB_ALLOC_ARRAY(arena, double, n);
    double *admm_prev_z = RIB_ALLOC_ARRAY(arena, double, n);
    double *admm_target = RIB_ALLOC_ARRAY(arena, double, n);
    int accepted_final_mode = 0;
    Arena_Mark final_strip_mark = Arena_save(arena);
    int final_rebuilds = 0;
    for (int round = 0; round < total_rounds; round++) {
        int relaxed = round < relaxed_rounds;
        int final_phase = round >= relaxed_rounds &&
                          round < relaxed_rounds + final_rounds;
        int metric_phase = round >= relaxed_rounds + final_rounds;
        int local_round = relaxed ? round : round - relaxed_rounds;
        int metric_round = metric_phase
                         ? round - relaxed_rounds - final_rounds : -1;
        const char *phase = rib_strip_round_phase(
            round, relaxed_rounds, final_rounds);
        double like_change = 0.0;

        if (final_phase || (metric_phase && final_rebuilds == 0)) {
            /* Section 5.3: form new cross-sections from the current parameter
             * isovalues before every final solve.  The dense samples make a
             * direct-sample intersection accurate to sample_h; mutual nearest
             * selection in U is the discrete counterpart of grouping stroke
             * intersections at a common isovalue.  Geometry is the stable tie
             * breaker.  New final cross-sections have unit membership: the
             * relaxed likelihoods classified the initial candidates, whereas
             * these correspondences are induced by the solved field itself. */
            if (final_rebuilds > 0) Arena_restore(arena, final_strip_mark);
            for (size_t i = 0; i < n; i++) S->smp[i].u = u[i];
            for (size_t p = 0; p < P->n_cross; p++)
                P->pairs[p].like = cross_pair_support[p];
            RibStripSet rebuilt;
            if (rib_strip_build_runs(arena, S, P, 1, 1.0, &rebuilt) != 0) {
                fprintf(stderr,
                        "ribbon: cannot rebuild final isovalue cross-sections\n");
                rib_gmg_snapshots_dispose(level_capture);
                RAISE(Arena_Failed);
            }
            strip = rebuilt;
            rib_strip_update_length_coefficients(S, &strip);
            cx.strip = &strip;
            out->n_strip_runs = strip.nrun;
            out->n_strip_members = strip.nmember;
            out->n_strip_links = strip.nlink;
            out->n_strip_links_pruned = strip.npruned;
            final_rebuilds++;
            fprintf(stderr,
                    "  StrokeStrip final C(u) rebuild %d/%d: runs=%zu "
                    "members=%zu links=%zu pruned=%zu\n",
                    final_rebuilds, final_rounds > 0 ? final_rounds : 1,
                    strip.nrun,
                    strip.nmember, strip.nlink, strip.npruned);
        }
        if (relaxed && local_round > 0) {
            like_change = rib_strip_update_likelihoods(
                &strip, u, local_round <= 3, strip_sigma);
            rib_strip_update_length_coefficients(S, &strip);
        }
        if (metric_phase) {
            /* Continuation, not IRLS: C(u), its tangent offsets, and robust
             * memberships are now frozen.  Only the ordinary quadratic weight
             * on measured stroke edges (and certified same-stroke fragment
             * continuations) increases.  Geometric interpolation avoids a
             * sudden ill-conditioned jump from the paper's gentle local term. */
            double target = (double)o->metric_weight;
            if (!isfinite(target) || target < RIB_STRIP_LOCAL_W)
                target = RIB_STRIP_LOCAL_W;
            double alpha = (double)(metric_round + 1) /
                           (double)metric_rounds;
            strip.local_weight = RIB_STRIP_LOCAL_W *
                pow(target / RIB_STRIP_LOCAL_W, alpha);
            strip.cont_weight = RIB_STRIP_CONT_W *
                pow(fmax(1.0, target / RIB_STRIP_CONT_W), alpha);
            fprintf(stderr,
                    "  StrokeStrip metric continuation %d/%d: "
                    "stroke_weight=%.4g continuation_weight=%.4g; "
                    "C(u) frozen\n",
                    metric_round + 1, metric_rounds,
                    strip.local_weight, strip.cont_weight);
        }
        cx.strip_final = !relaxed;
        double round_t0 = ves_clock_sec();
        /* The alignment Schur complement changes with every robust-membership
         * update.  Rebuild its exact signed-row Galerkin hierarchy from the
         * fixed ribbon topology and the live weights; the initial O(1e-5)
         * system and a later O(1) likelihood system are different operators.
         * The round mark releases the hierarchy and solve workspace. */
        Arena_Mark solve_mark = Arena_save(arena);
        RibAMG amg;
        memset(&amg, 0, sizeof amg);
        double amg_t0 = ves_clock_sec();
        /* The TAUCS factor is built LAZILY inside the first cg_solve that
         * actually iterates: a round whose warm iterate already meets the
         * tolerance never pays for a factorization.  The GMG hierarchy is
         * retired from this path (its FMG correction was rejected every
         * round on the 4x5x5 soup); Jacobi-PCG is the fallback when the
         * factor cannot be built. */
        cx.round_no = round;
        cx.chol_attempted = 0;
        double amg_sec = ves_clock_sec() - amg_t0;
        cx.amg = NULL;
        if (round == 0 && amg.enabled) {
            fprintf(stderr,
                    "  solve GMG full grid: nodes=%zu hU=%.3g hV=%.3g; "
                    "first U-block=%d nodes=%d levels=%d policy=%s "
                    "build_sec=%.2f\n",
                    S->n_smp, (double)o->sample_h, (double)o->slice_h,
                    amg.block_size, amg.n_sample_agg, amg.nlevel,
                    amg.weak_v ? "geometric-semi (weak-V operator)"
                               : "geometric-semi",
                    amg_sec);
            for (int l = 0; l < amg.nlevel; l++) {
                const RibAMGLevel *L = &amg.level[l];
                fprintf(stderr,
                        "    GMG L%d via %c: hU=%.3g hV=%.3g nodes=%d "
                        "topoU=%zu topoV=%zu graph=%zu signed_rows=%zu/%zu "
                        "smooth=%d\n",
                        l, L->from_axis, L->h_u, L->h_v, L->n,
                        L->nu_edge, L->nv_edge, L->ne,
                        L->nhrow, L->nhcoeff,
                        L->smooth_sweeps > 0 ? L->smooth_sweeps : 0);
            }
            {
                const RibAMGLevel *coarse = &amg.level[amg.nlevel - 1];
                fprintf(stderr,
                        "    GMG coarsest: blocks=%d chol_values=%zu; "
                        "FMG=%d V-cycle/level + %d fine correction\n",
                        coarse->n_chol_blocks, coarse->chol_values,
                        RIB_GMG_FMG_CYCLES, RIB_GMG_FINE_CYCLES);
            }
        }
        double cg_relres = 0.0;
        double cg_t0 = ves_clock_sec();
        RibGMGSnapshots *round_capture =
            level_capture != NULL && round + 1 == total_rounds
                ? level_capture : NULL;
        /* Solve the actual monotone quadratic, rather than solving an
         * unconstrained system and projecting it once afterward.  The latter
         * is not a QP solve: on the 10x ribbon it moved nearly two million
         * samples and converted smooth parameter lines into broad flat wedges.
         *
         * Scaled ADMM splits the quadratic variable u from one feasible copy z:
         *   u <- argmin E(u) + rho/2 ||u-z+y||^2
         *   z <- projection_monotone(u+y)
         *   y <- y + u-z.
         * The z step is the exact per-chain PAVA projection, while the u step
         * remains SPD and uses the existing matrix-free solver. */
        int cg_iters = 0, admm_iters = 0;
        double primal_rms = 0.0, dual_rms = 0.0;
        size_t mono_round = 0;
        cx.ridge_weight = 0.0;
        cx.ridge_target = NULL;
        /* On a topology-preserving quadribbon, 1e-4 relative residual is
         * already far below the sampled voxel scale.  The 1e-6 research
         * tolerance costs hundreds of Krylov iterations on million-sample
         * systems without a visible or metric benefit.  Feasibility is not
         * relaxed: the published ADMM z remains an exact PAVA projection. */
        int sampled_metric_tol = o->scaffold_solve || o->preserve_input_rows;
        cx.tol = sampled_metric_tol ? RIB_STRIP_SCAFFOLD_TOL
                                    : RIB_STRIP_CG_TOL;
        int base_iters = 0;
        /* A scaffold may skip the unconstrained warm solve: ADMM's first
         * x-update performs the same work with the ridge already in place. */
        if (!o->scaffold_solve) {
            base_iters = cg_solve(arena, &cx, u, &cg_relres, NULL);
            if (base_iters < 0) {
                fprintf(stderr,
                        "ribbon: unconstrained monotone-QP warm solve failed\n");
                rib_strip_chol_release(&cx);
                rib_gmg_snapshots_dispose(level_capture);
                RAISE(Arena_Failed);
            }
            cg_iters += base_iters;
        }
        memcpy(admm_z, u, n * sizeof(*admm_z));
        mono_round = rib_project_monotone(
            S, admm_z, scr, lb, thread_moved, maxchain, solve_threads);
        if (mono_round == 0) {
            memcpy(u, admm_z, n * sizeof(*u));
            fprintf(stderr,
                    "    monotone QP: unconstrained minimizer is feasible "
                    "(cg=%d relres=%.3g)\n",
                    base_iters, cg_relres);
        } else {
        double admm_rho = RIB_STRIP_ADMM_RHO;
        memset(admm_y, 0, n * sizeof(*admm_y));
        memcpy(u, admm_z, n * sizeof(*u));
        cx.ridge_weight = admm_rho;
        cx.ridge_target = admm_target;
        int admm_maxit = o->scaffold_solve ? 3
                       : o->preserve_input_rows ? 8
                       : RIB_STRIP_ADMM_MAXIT;
        for (int admm = 0; admm < admm_maxit; admm++) {
            memcpy(admm_prev_z, admm_z, n * sizeof(*admm_prev_z));
            int si_target;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(solve_threads)
#endif
            for (si_target = 0; si_target < (int)n; si_target++)
                admm_target[si_target] =
                    admm_z[si_target] - admm_y[si_target];
            int inner_iters = cg_solve(arena, &cx, u, &cg_relres, NULL);
            if (inner_iters < 0) {
                fprintf(stderr, "ribbon: monotone ADMM inner solve failed\n");
                rib_strip_chol_release(&cx);
                rib_gmg_snapshots_dispose(level_capture);
                RAISE(Arena_Failed);
            }
            cg_iters += inner_iters;
            int si_project;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(solve_threads)
#endif
            for (si_project = 0; si_project < (int)n; si_project++)
                admm_z[si_project] = u[si_project] + admm_y[si_project];
            mono_round = rib_project_monotone(
                S, admm_z, scr, lb, thread_moved, maxchain, solve_threads);
            double primal2 = 0.0, dual2 = 0.0;
            int si_residual;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:primal2,dual2) num_threads(solve_threads)
#endif
            for (si_residual = 0; si_residual < (int)n; si_residual++) {
                double primal = u[si_residual] - admm_z[si_residual];
                double dual = admm_rho *
                    (admm_z[si_residual] - admm_prev_z[si_residual]);
                primal2 += primal * primal;
                dual2 += dual * dual;
                admm_y[si_residual] += primal;
            }
            primal_rms = sqrt(primal2 / (double)n);
            dual_rms = sqrt(dual2 / (double)n);
            admm_iters = admm + 1;
            if (admm == 0 || ((admm + 1) % 5) == 0 ||
                (primal_rms <= RIB_STRIP_ADMM_TOL &&
                 dual_rms <= RIB_STRIP_ADMM_TOL))
                fprintf(stderr,
                        "    monotone ADMM: iter=%d rho=%.3g primal=%.3g "
                        "dual=%.3g active_projection=%zu inner_cg=%d "
                        "inner_relres=%.3g\n",
                        admm + 1, admm_rho, primal_rms, dual_rms,
                        mono_round, inner_iters, cg_relres);
            if (primal_rms <= RIB_STRIP_ADMM_TOL &&
                dual_rms <= RIB_STRIP_ADMM_TOL)
                break;
            if (primal_rms > 10.0 * dual_rms && admm_rho < 1.0) {
                admm_rho = fmin(1.0, 2.0 * admm_rho);
                for (size_t si = 0; si < n; si++) admm_y[si] *= 0.5;
                cx.ridge_weight = admm_rho;
            } else if (dual_rms > 10.0 * primal_rms &&
                       admm_rho > 1e-4) {
                admm_rho *= 0.5;
                for (size_t si = 0; si < n; si++) admm_y[si] *= 2.0;
                cx.ridge_weight = admm_rho;
            }
        }
        }
        memcpy(u, admm_z, n * sizeof(*u));
        cx.ridge_weight = 0.0;
        cx.ridge_target = NULL;
        cx.tol = sampled_metric_tol ? RIB_STRIP_SCAFFOLD_TOL
                                    : RIB_STRIP_CG_TOL;
        if (admm_iters > 0 && (primal_rms > RIB_STRIP_ADMM_TOL ||
            dual_rms > RIB_STRIP_ADMM_TOL)
           )
            fprintf(stderr,
                    "    monotone ADMM WARNING: iteration cap at primal=%.3g "
                    "dual=%.3g (final z remains exactly feasible)\n",
                    primal_rms, dual_rms);
        double cg_sec = ves_clock_sec() - cg_t0;
        cx.amg = NULL;
        rib_strip_chol_release(&cx);
        Arena_restore(arena, solve_mark);
        mono_total += mono_round;
        if (round_capture != NULL &&
            rib_gmg_capture_full_level(
                S, u, (double)o->sample_h, (double)o->slice_h,
                round_capture) != 0) {
            fprintf(stderr, "ribbon: cannot capture full sample solution\n");
            rib_gmg_snapshots_dispose(level_capture);
            RAISE(Arena_Failed);
        }
        double length_rms = 0.0, align_rms = 0.0;
        double local_mean = 0.0, local_rms = 0.0;
        rib_strip_measure(S, &strip, u, cx.strip_final,
                          &length_rms, &align_rms);
        rib_strip_measure_local_metric(S, u, &local_mean, &local_rms);
        if (!isfinite(length_rms) || !isfinite(align_rms) ||
            !isfinite(local_mean) || !isfinite(local_rms)) {
            fprintf(stderr,
                    "ribbon: non-finite StrokeStrip metric after round %d\n",
                    round + 1);
            rib_gmg_snapshots_dispose(level_capture);
            RAISE(Arena_Failed);
        }
        fprintf(stderr,
                "  StrokeStrip %s round %d/%d: admm_iters=%d "
                "cg_iters=%d relres=%.3g "
                "precond_sec=%.2f cg_sec=%.2f round_sec=%.2f "
                "mono=%zu like_delta=%.3g "
                "length_rms=%.4f local_mean=%.4f local_rms=%.4f "
                "align_rms=%.4f sigma=%.4f\n",
                phase, round + 1, total_rounds, admm_iters,
                cg_iters, cg_relres,
                amg_sec, cg_sec, ves_clock_sec() - round_t0,
                mono_round, like_change,
                length_rms, local_mean, local_rms, align_rms, strip_sigma);

        /* Do not veto a valid StrokeStrip step merely because individual
         * stroke speeds temporarily move away from one.  The paper constrains
         * the cross-section average here.  Per-edge unit speed is recovered by
         * the later metric-continuation phase after C(u) has stabilized. */
        accepted_final_mode = cx.strip_final;
        if (o->solve_level_prefix != NULL) {
            /* A round is useful to a human only when it is visible while the
             * following round is still running.  Copy the just-projected U to
             * the SliceSet and atomically publish the fitted VMESH now.  Robust
             * memberships are solve weights only; they never change geometry
             * existence or reconstruction identity. */
            for (size_t i = 0; i < n; i++) S->smp[i].u = u[i];
            if (rib_dump_solve_round_checkpoint(
                    arena, S, P, o, out, round, total_rounds, phase) != 0) {
                fprintf(stderr,
                        "ribbon: solve-round checkpoint %d/%d failed\n",
                        round + 1, total_rounds);
                rib_gmg_snapshots_dispose(level_capture);
                RAISE(Arena_Failed);
            }
        }
    }

    /* Pair.like remains the fixed geometric support selected when V-runs were
     * constructed.  The alternating robust memberships never leak into
     * downstream topology. */
    rib_strip_measure(S, &strip, u, accepted_final_mode,
                      &out->strip_length_rms, &out->strip_align_rms);
    rib_strip_audit_pruned_gauges(arena, S, P, &strip, u, nsolve);
    rib_strip_audit_phase_gauges(arena, S, u, nsolve);
    if (o->stitch_solve_gauges)
        rib_strip_stitch_solve_gauges(
            arena, S, u, nsolve, level_capture, out);

    /* --- REGISTRATION fallback: place every detached solve component onto the
     * main chart when no monotone cover gauge was available.  Build the main
     * component's monotone phi->u map and robustly rigid-shift every other
     * component by median(f(phi)-u).  This changes one additive constant only,
     * preserving each piece's internal arc length;
     * places outer-shell arcs BEYOND the main band (their phi is higher) and
     * genuinely delaminated layers exactly ON the main chart (same phi).
     * Outside the main chart's observed phi range this is explicitly a smooth
     * atlas extrapolation, not a claim that radius determined physical U. */
    out->reg_max_shift = 0.0;
    if (n_pins > 1 && !spiral_ok && !use_reference_seed) {
        Arena_Mark regm = Arena_save(arena);
        int32_t mainroot = uf_find(&uf, bigchain);
        /* per-sample solve-component root (resolved once) */
        int32_t *sroot = RIB_ALLOC_ARRAY(arena, int32_t, n + 1);
        for (size_t i = 0; i < n; i++) {
            int32_t c = S->smp[i].chain;
            sroot[i] = c >= 0 ? uf_find(&uf, c) : -1;
        }
        /* gather main-component (phi, u) sorted by phi */
        size_t nmain = 0;
        for (size_t i = 0; i < n; i++) if (sroot[i] == mainroot) nmain++;
        PhiU2 *pu = RIB_ALLOC_ARRAY(arena, PhiU2, nmain + 1);
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
            M.bphi = RIB_ALLOC_ARRAY(arena, double, (size_t)nb_max);
            M.bu   = RIB_ALLOC_ARRAY(arena, double, (size_t)nb_max);
            double *tmpu = RIB_ALLOC_ARRAY(arena, double, nmain + 1);
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
            /* Winding order is monotone but noisy axial slices need not give
             * perfectly monotone bin medians.  Isotonic repair prevents the
             * empirical gauge map itself from folding.  It changes only chart
             * offsets; every chart's intrinsic arclength remains untouched. */
            if (M.nb > 1) {
                double *map_lb = RIB_ALLOC_ARRAY(arena, double, (size_t)M.nb);
                double *map_scr = RIB_ALLOC_ARRAY(
                    arena, double, (size_t)M.nb * 3u);
                memset(map_lb, 0, (size_t)M.nb * sizeof(double));
                (void)pava_chain(M.bu, map_lb, (int32_t)M.nb, map_scr);
            }
        }
        if (M.nb >= 1) {
            double *diffs = RIB_ALLOC_ARRAY(arena, double, n + 1);
            /* Per detached component: change only the additive gauge. */
            for (size_t c = 0; c < nc; c++) {
                if (uf_find(&uf, (int32_t)c) != (int32_t)c) continue;
                if ((int32_t)c == mainroot) continue;
                double u_lo = 1e300, u_hi = -1e300, ph_lo = 1e300, ph_hi = -1e300;
                size_t nd = 0;
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
                for (size_t i = 0; i < n; i++) {
                    if (sroot[i] != (int32_t)c) continue;
                    u[i] += shift;
                    if (level_capture != NULL)
                        for (int l = 0; l < level_capture->count; l++)
                            level_capture->snapshot[l].u[i] += (float)shift;
                }
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
            if (level_capture != NULL)
                for (int l = 0; l < level_capture->count; l++)
                    for (size_t i = 0; i < n; i++)
                        level_capture->snapshot[l].u[i] =
                            -level_capture->snapshot[l].u[i];
            for (size_t c = 0; c < nc; c++) {
                int32_t f = S->chn[c].first, cn = S->chn[c].count;
                if (cn < 2) continue;
                for (int32_t i = 0; i < cn / 2; i++) {
                    Sample tmp = S->smp[f + i];
                    S->smp[f + i] = S->smp[f + cn - 1 - i];
                    S->smp[f + cn - 1 - i] = tmp;
                }
                if (level_capture != NULL) {
                    for (int l = 0; l < level_capture->count; l++) {
                        RibGMGSnapshot *D = &level_capture->snapshot[l];
                        for (int32_t i = 0; i < cn / 2; i++) {
                            float tu = D->u[f + i];
                            int32_t tn = D->node[f + i];
                            D->u[f + i] = D->u[f + cn - 1 - i];
                            D->node[f + i] = D->node[f + cn - 1 - i];
                            D->u[f + cn - 1 - i] = tu;
                            D->node[f + cn - 1 - i] = tn;
                        }
                    }
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
 * Stage D -- transfer (u, v) to the original mesh vertices via the nearest
 * same-sheet slice sample.  In-plane distance is only a spatial acceleration
 * gate: contacting plies can be closer than RIB_XFER_R.  When the caller
 * supplies authoritative lifted phase, |sample.phi-reference_phi[vertex]| is
 * therefore also required to be small.  Tangential correction is
 * tau . (p_v - p_s).
 * ==========================================================================*/

typedef struct {
    int32_t a, b;             /* equation shift[b] - shift[a] = delta */
    double delta;
} TransferGaugeObs;

typedef struct {
    int32_t a, b;
    double delta;
    size_t support;
} TransferGaugeRel;

static int cmp_transfer_gauge_obs(const void *pa, const void *pb)
{
    const TransferGaugeObs *a = (const TransferGaugeObs *)pa;
    const TransferGaugeObs *b = (const TransferGaugeObs *)pb;
    if (a->a != b->a) return a->a < b->a ? -1 : 1;
    if (a->b != b->b) return a->b < b->b ? -1 : 1;
    return a->delta < b->delta ? -1 : (a->delta > b->delta ? 1 : 0);
}

static int cmp_transfer_gauge_rel(const void *pa, const void *pb)
{
    const TransferGaugeRel *a = (const TransferGaugeRel *)pa;
    const TransferGaugeRel *b = (const TransferGaugeRel *)pb;
    if (a->support != b->support) return a->support > b->support ? -1 : 1;
    if (a->a != b->a) return a->a < b->a ? -1 : 1;
    return a->b < b->b ? -1 : (a->b > b->b ? 1 : 0);
}

static int transfer_edge_gauge_observation(
    const float *verts, int32_t ia, int32_t ib,
    const float *reference_phi,
    const int32_t *vsolve, const double *vu, const float *vtau,
    TransferGaugeObs *out)
{
    int32_t ca = vsolve[ia], cb = vsolve[ib];
    if (ca < 0 || cb < 0 || ca == cb) return 0;
    if (reference_phi != NULL) {
        double pa = (double)reference_phi[ia];
        double pb = (double)reference_phi[ib];
        if (!isfinite(pa) || !isfinite(pb) ||
            fabs(pa - pb) > RIB_XFER_DPHI_MAX)
            return 0;
    }
    double dp[3] = {
        (double)verts[(size_t)ib*3]   - (double)verts[(size_t)ia*3],
        (double)verts[(size_t)ib*3+1] - (double)verts[(size_t)ia*3+1],
        (double)verts[(size_t)ib*3+2] - (double)verts[(size_t)ia*3+2]
    };
    double l2 = v3dot(dp, dp);
    if (!(l2 > 1e-12) || l2 > RIB_GAUGE_EDGE_MAX * RIB_GAUGE_EDGE_MAX)
        return 0;
    const float *ta = &vtau[(size_t)ia*3];
    const float *tb = &vtau[(size_t)ib*3];
    double tdot = (double)ta[0]*(double)tb[0] +
                  (double)ta[1]*(double)tb[1] +
                  (double)ta[2]*(double)tb[2];
    if (tdot < 0.5) return 0;
    double expected =
        0.5 * ((double)ta[0] + (double)tb[0]) * dp[0] +
        0.5 * ((double)ta[1] + (double)tb[1]) * dp[1] +
        0.5 * ((double)ta[2] + (double)tb[2]) * dp[2];
    /* (vu_b+s_cb) - (vu_a+s_ca) = expected. */
    double delta = expected - (vu[ib] - vu[ia]);
    if (ca < cb) {
        out->a = ca; out->b = cb; out->delta = delta;
    } else {
        out->a = cb; out->b = ca; out->delta = -delta;
    }
    return 1;
}

/* The slice-pair graph can be disconnected across a hole even though the
 * original source chart is continuous.  Each disconnected metric solve then
 * has one free additive U gauge.  Reconcile only those constants using:
 *   (1) two bracketing slice samples evaluated at the SAME source vertex;
 *   (2) short source edges, with their expected du projected on the solved
 *       slice tangents.
 * No radius-to-turn model enters this step.  Pairwise medians reject local
 * transfer noise; a maximum-support forest prevents inconsistent cycles from
 * over-constraining the gauges. */
static void reconcile_transfer_gauges(
    Arena_T arena, SliceSet *S,
    const float *verts, size_t nv, const int32_t *faces, size_t nf,
    const float *reference_phi, const float *vw,
    double *u0, double *u1, const double *ph0, const double *ph1,
    const int32_t *si0, const int32_t *si1,
    const uint8_t *ok0, const uint8_t *ok1,
    RibGMGSnapshots *level_capture, RibbonResult *out)
{
    int32_t ncomp = out->n_qp_comps;
    if (ncomp <= 1) return;
    Arena_Mark mark = Arena_save(arena);

    int32_t *vsolve = RIB_ALLOC_ARRAY(arena, int32_t, nv);
    double *vu = RIB_ALLOC_ARRAY(arena, double, nv);
    float *vtau = RIB_ALLOC_ARRAY(arena, float, nv * 3);
    uint8_t *vok = (uint8_t *)ARENA_CALLOC(arena, nv, 1);
    size_t nobs = 0;
    for (size_t i = 0; i < nv; i++) {
        int choose = -1;
        int32_t c0 = -1, c1 = -1;
        if (ok0[i]) c0 = S->chn[S->smp[si0[i]].chain].solve_comp;
        if (ok1[i]) c1 = S->chn[S->smp[si1[i]].chain].solve_comp;
        if (ok0[i] && ok1[i] && c0 != c1 &&
            fabs(ph0[i] - ph1[i]) <= RIB_XFER_DPHI_MAX)
            nobs++;
        if (ok0[i] && ok1[i]) {
            if (c0 == c1) {
                double w = (double)vw[i];
                vu[i] = (1.0-w)*u0[i] + w*u1[i];
                vsolve[i] = c0;
                for (int d = 0; d < 3; d++)
                    vtau[i*3+(size_t)d] =
                        (float)((1.0-w)*S->smp[si0[i]].tau[d] +
                                      w *S->smp[si1[i]].tau[d]);
                vok[i] = 1;
                {
                    double tn = sqrt((double)vtau[i*3]*(double)vtau[i*3] +
                                     (double)vtau[i*3+1]*(double)vtau[i*3+1] +
                                     (double)vtau[i*3+2]*(double)vtau[i*3+2]);
                    if (tn > 1e-12)
                        for (int d = 0; d < 3; d++)
                            vtau[i*3+(size_t)d] =
                                (float)((double)vtau[i*3+(size_t)d] / tn);
                }
                continue;
            }
            choose = (double)vw[i] < 0.5 ? 0 : 1;
        } else if (ok0[i]) choose = 0;
        else if (ok1[i]) choose = 1;
        if (choose >= 0) {
            int32_t is = choose == 0 ? si0[i] : si1[i];
            vu[i] = choose == 0 ? u0[i] : u1[i];
            vsolve[i] = S->chn[S->smp[is].chain].solve_comp;
            for (int d = 0; d < 3; d++)
                vtau[i*3+(size_t)d] = (float)S->smp[is].tau[d];
            vok[i] = 1;
        } else {
            vu[i] = 0.0;
            vsolve[i] = -1;
            memset(&vtau[i*3], 0, 3*sizeof(float));
        }
    }

    TransferGaugeObs dummy;
    for (size_t f = 0; f < nf; f++) {
        int32_t tri[3] = { faces[f*3], faces[f*3+1], faces[f*3+2] };
        for (int e = 0; e < 3; e++) {
            int32_t a = tri[e], b = tri[(e+1)%3];
            if (!vok[a] || !vok[b]) continue;
            nobs += (size_t)transfer_edge_gauge_observation(
                verts, a, b, reference_phi, vsolve, vu, vtau, &dummy);
        }
    }
    if (nobs == 0) {
        Arena_restore(arena, mark);
        return;
    }

    TransferGaugeObs *obs = RIB_ALLOC_ARRAY(arena, TransferGaugeObs, nobs);
    size_t no = 0;
    for (size_t i = 0; i < nv; i++) {
        if (!ok0[i] || !ok1[i]) continue;
        int32_t c0 = S->chn[S->smp[si0[i]].chain].solve_comp;
        int32_t c1 = S->chn[S->smp[si1[i]].chain].solve_comp;
        if (c0 == c1 || fabs(ph0[i] - ph1[i]) > RIB_XFER_DPHI_MAX)
            continue;
        double delta = u0[i] - u1[i]; /* shift[c1]-shift[c0] */
        if (c0 < c1) {
            obs[no].a = c0; obs[no].b = c1; obs[no].delta = delta;
        } else {
            obs[no].a = c1; obs[no].b = c0; obs[no].delta = -delta;
        }
        no++;
    }
    for (size_t f = 0; f < nf; f++) {
        int32_t tri[3] = { faces[f*3], faces[f*3+1], faces[f*3+2] };
        for (int e = 0; e < 3; e++) {
            int32_t a = tri[e], b = tri[(e+1)%3];
            if (!vok[a] || !vok[b]) continue;
            if (transfer_edge_gauge_observation(
                    verts, a, b, reference_phi, vsolve, vu, vtau, &obs[no]))
                no++;
        }
    }
    assert(no == nobs);
    qsort(obs, nobs, sizeof(*obs), cmp_transfer_gauge_obs);

    TransferGaugeRel *rel = RIB_ALLOC_ARRAY(arena, TransferGaugeRel, nobs);
    size_t nr = 0;
    for (size_t i = 0; i < nobs; ) {
        size_t j = i + 1;
        while (j < nobs && obs[j].a == obs[i].a && obs[j].b == obs[i].b) j++;
        rel[nr].a = obs[i].a;
        rel[nr].b = obs[i].b;
        rel[nr].delta = obs[i + (j-i)/2].delta;
        rel[nr].support = j - i;
        nr++;
        i = j;
    }
    qsort(rel, nr, sizeof(*rel), cmp_transfer_gauge_rel);

    UnionFind guf = UF_new(arena, ncomp);
    TransferGaugeRel *tree = RIB_ALLOC_ARRAY(arena, TransferGaugeRel, (size_t)ncomp);
    size_t nt = 0;
    for (size_t i = 0; i < nr; i++) {
        if (uf_find(&guf, rel[i].a) == uf_find(&guf, rel[i].b)) continue;
        tree[nt++] = rel[i];
        uf_union(&guf, rel[i].a, rel[i].b);
    }
    double *shift = (double *)ARENA_CALLOC(arena, (size_t)ncomp, sizeof(double));
    uint8_t *known = (uint8_t *)ARENA_CALLOC(arena, (size_t)ncomp, 1);
    for (int32_t seed = 0; seed < ncomp; seed++) {
        if (known[seed]) continue;
        known[seed] = 1;
        int progress;
        do {
            progress = 0;
            for (size_t i = 0; i < nt; i++) {
                int32_t a = tree[i].a, b = tree[i].b;
                if (known[a] && !known[b]) {
                    shift[b] = shift[a] + tree[i].delta;
                    known[b] = 1; progress = 1;
                } else if (!known[a] && known[b]) {
                    shift[a] = shift[b] - tree[i].delta;
                    known[a] = 1; progress = 1;
                }
            }
        } while (progress);
    }

    double maxshift = 0.0;
    for (int32_t c = 0; c < ncomp; c++)
        if (fabs(shift[c]) > maxshift) maxshift = fabs(shift[c]);
    for (size_t i = 0; i < S->n_smp; i++) {
        int32_t ch = S->smp[i].chain;
        if (ch >= 0) {
            double ds = shift[S->chn[ch].solve_comp];
            S->smp[i].u += ds;
            if (level_capture != NULL)
                for (int l = 0; l < level_capture->count; l++)
                    level_capture->snapshot[l].u[i] += (float)ds;
        }
    }
    for (size_t i = 0; i < nv; i++) {
        if (ok0[i]) {
            int32_t sc = S->chn[S->smp[si0[i]].chain].solve_comp;
            u0[i] += shift[sc];
        }
        if (ok1[i]) {
            int32_t sc = S->chn[S->smp[si1[i]].chain].solve_comp;
            u1[i] += shift[sc];
        }
    }
    out->uv_gauge_observations = nobs;
    out->uv_gauge_relations = nt;
    out->uv_gauge_max_shift = maxshift;
    Arena_restore(arena, mark);
}

/* The universal-cover phase must keep disconnected relation islands distinct,
 * so winding registration gives them deterministic serial q intervals.  Those
 * intervals do NOT measure material distance: with no relation path, their
 * additive U offset is unobservable.  Once Stage C has measured each island's
 * intrinsic arclength interval, compact the intervals for the final atlas.
 *
 * The gutter is larger than the fitted ribbon's maximum horizontal hole fill,
 * guaranteeing that atlas neighbours cannot become geometric neighbours.  All
 * shifts are additive per island, so every within-island du and every local
 * arclength/monotonicity invariant is preserved exactly. */
static void pack_metric_islands(
    Arena_T arena, SliceSet *S, double grid_u,
    double *u0, double *u1,
    const int32_t *si0, const int32_t *si1,
    const uint8_t *ok0, const uint8_t *ok1, size_t nv,
    RibGMGSnapshots *level_capture, RibbonResult *out)
{
    int32_t max_island = -1;
    for (size_t i = 0; i < S->n_smp; i++) {
        int32_t ch = S->smp[i].chain;
        if (ch < 0) continue;
        int32_t island = S->chn[ch].reconstruction_component;
        if (island > max_island) max_island = island;
    }
    if (max_island < 0) return;

    size_t nisland = (size_t)max_island + 1;
    Arena_Mark mark = Arena_save(arena);
    double *lo = RIB_ALLOC_ARRAY(arena, double, nisland);
    double *hi = RIB_ALLOC_ARRAY(arena, double, nisland);
    double *shift = RIB_ALLOC_ARRAY(arena, double, nisland);
    uint8_t *active = RIB_ALLOC_ARRAY(arena, uint8_t, nisland);
    for (size_t g = 0; g < nisland; g++) {
        lo[g] = 1e300; hi[g] = -1e300; shift[g] = 0.0; active[g] = 0;
    }

    double old_lo = 1e300, old_hi = -1e300;
    for (size_t i = 0; i < S->n_smp; i++) {
        int32_t ch = S->smp[i].chain;
        if (ch < 0) continue;
        int32_t island = S->chn[ch].reconstruction_component;
        if (island < 0 || (size_t)island >= nisland) continue;
        double u = S->smp[i].u;
        if (u < lo[island]) lo[island] = u;
        if (u > hi[island]) hi[island] = u;
        if (u < old_lo) old_lo = u;
        if (u > old_hi) old_hi = u;
        active[island] = 1;
    }

    double du = grid_u > 1e-6 && isfinite(grid_u) ? grid_u : 1.0;
    double gutter = (double)(RIB_UFILL_MAX + 2) * du;
    if (gutter < 16.0) gutter = 16.0;
    size_t nactive = 0;
    double cursor = 0.0, new_lo = 1e300, new_hi = -1e300;
    for (size_t g = 0; g < nisland; g++) {
        if (!active[g]) continue;
        if (nactive == 0) {
            shift[g] = 0.0;
            cursor = hi[g];
        } else {
            shift[g] = cursor + gutter - lo[g];
            cursor = hi[g] + shift[g];
        }
        if (lo[g] + shift[g] < new_lo) new_lo = lo[g] + shift[g];
        if (hi[g] + shift[g] > new_hi) new_hi = hi[g] + shift[g];
        nactive++;
    }

    out->uv_atlas_islands = nactive;
    out->uv_atlas_packed_islands = nactive > 0 ? nactive - 1 : 0;
    out->uv_atlas_gutter = nactive > 1 ? gutter : 0.0;
    if (old_lo < 1e299 && new_lo < 1e299) {
        double old_span = old_hi - old_lo;
        double new_span = new_hi - new_lo;
        out->uv_atlas_pack_saved = old_span > new_span
                                  ? old_span - new_span : 0.0;
    }

    for (size_t c = 0; c < S->n_chn; c++)
        S->chn[c].atlas_shift =
            shift[S->chn[c].reconstruction_component];
    if (nactive > 1) {
        for (size_t i = 0; i < S->n_smp; i++) {
            int32_t ch = S->smp[i].chain;
            if (ch < 0) continue;
            int32_t island = S->chn[ch].reconstruction_component;
            S->smp[i].u += shift[island];
            if (level_capture != NULL)
                for (int l = 0; l < level_capture->count; l++)
                    level_capture->snapshot[l].u[i] += (float)shift[island];
        }
        for (size_t i = 0; i < nv; i++) {
            if (ok0 != NULL && ok0[i] && si0[i] >= 0) {
                int32_t ch = S->smp[si0[i]].chain;
                if (ch >= 0)
                    u0[i] += shift[S->chn[ch].reconstruction_component];
            }
            if (ok1 != NULL && ok1[i] && si1[i] >= 0) {
                int32_t ch = S->smp[si1[i]].chain;
                if (ch >= 0)
                    u1[i] += shift[S->chn[ch].reconstruction_component];
            }
        }
    }
    fprintf(stderr,
            "  metric atlas: islands=%zu packed=%zu gutter=%.1f "
            "cover_span=%.1f atlas_span=%.1f saved=%.1f\n",
            nactive, out->uv_atlas_packed_islands,
            out->uv_atlas_gutter,
            old_lo < 1e299 ? old_hi - old_lo : 0.0,
            new_lo < 1e299 ? new_hi - new_lo : 0.0,
            out->uv_atlas_pack_saved);
    if (nactive <= 32) {
        for (size_t g = 0; g < nisland; g++) {
            if (!active[g]) continue;
            fprintf(stderr,
                    "    island %zu: metric=[%.1f,%.1f] shift=%+.1f "
                    "atlas=[%.1f,%.1f]\n",
                    g, lo[g], hi[g], shift[g],
                    lo[g] + shift[g], hi[g] + shift[g]);
        }
    }
    Arena_restore(arena, mark);
}

static double transfer_uv(Arena_T arena, const float *verts, size_t nv,
                          const int32_t *faces, size_t nf,
                          const int32_t *vertex_mesh_comp,
                          size_t n_mesh_comp,
                           const double *t,
                           SliceSet *S, double slice_h,
                           const RibbonOpts *o,
                           RibGMGSnapshots *level_capture,
                           RibbonResult *out)
{
    double h = slice_h;
    double tmin = S->tmin;
    int np = S->nplanes;
    Sample *smp = S->smp;

    /* out->uv arrives PREFILLED by the caller with each vertex's in-plane
     * (c1,c2) coords -- the bucketing key for the per-slice grids below. The
     * final (u,v) overwrites it at the end. */
    out->uv_ok = (uint8_t *)ARENA_CALLOC(arena, nv, 1);

    Arena_Mark mark = Arena_save(arena);

    size_t *mesh_sample_count = (size_t *)ARENA_CALLOC(
        arena, n_mesh_comp ? n_mesh_comp : 1, sizeof(size_t));
    for (size_t i = 0; i < S->n_smp; i++) {
        int32_t ch = smp[i].chain;
        if (ch < 0) continue;
        int32_t mc = S->chn[ch].mesh_comp;
        if (mc >= 0 && (size_t)mc < n_mesh_comp)
            mesh_sample_count[mc]++;
    }

    /* per-vertex bracketing slices + blend weight */
    int32_t *vk0 = RIB_ALLOC_ARRAY(arena, int32_t, nv);
    float   *vw  = RIB_ALLOC_ARRAY(arena, float, nv);
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

    double *u0 = RIB_ALLOC_ARRAY(arena, double, nv);
    double *u1 = RIB_ALLOC_ARRAY(arena, double, nv);
    double *ph0 = (double *)ARENA_ALLOC(arena, nv * sizeof(double));
    double *ph1 = (double *)ARENA_ALLOC(arena, nv * sizeof(double));
    int32_t *gr0 = (int32_t *)ARENA_ALLOC(arena, nv * sizeof(int32_t));
    int32_t *gr1 = (int32_t *)ARENA_ALLOC(arena, nv * sizeof(int32_t));
    int32_t *si0 = RIB_ALLOC_ARRAY(arena, int32_t, nv);
    int32_t *si1 = RIB_ALLOC_ARRAY(arena, int32_t, nv);
    for (size_t i = 0; i < nv; i++) si0[i] = si1[i] = -1;
    uint8_t *ok0 = (uint8_t *)ARENA_CALLOC(arena, nv, 1);
    uint8_t *ok1 = (uint8_t *)ARENA_CALLOC(arena, nv, 1);
    uint8_t *mapped = (uint8_t *)ARENA_CALLOC(arena, nv, 1);

    /* vertices bucketed by k0 (counting sort) so each slice's grid is built
     * once and queried by the two vertex populations that bracket it */
    int32_t *voff = (int32_t *)ARENA_CALLOC(
        arena, (size_t)np + 2, sizeof(int32_t));
    for (size_t i = 0; i < nv; i++) voff[vk0[i] + 1]++;
    for (int k = 0; k <= np; k++) voff[k+1] = (int32_t)(voff[k+1] + voff[k]);
    int32_t *vids = RIB_ALLOC_ARRAY(arena, int32_t, nv);
    int32_t *vcur = RIB_ALLOC_ARRAY(arena, int32_t, (size_t)np + 1);
    memcpy(vcur, voff, (size_t)(np + 1) * sizeof(int32_t));
    for (size_t i = 0; i < nv; i++) vids[vcur[vk0[i]]++] = (int32_t)i;

    int32_t *goff = RIB_ALLOC_ARRAY(arena, int32_t, ncell + 1);
    int32_t *gids = RIB_ALLOC_ARRAY(arena, int32_t, S->n_smp + 1);
    int32_t *gcur = RIB_ALLOC_ARRAY(arena, int32_t, ncell + 1);
    size_t phase_rejects = 0, phase_blocked = 0;
    size_t component_rejects = 0, component_blocked = 0;

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
                int32_t best = -1;
                double xfer_r2 = RIB_XFER_R * RIB_XFER_R;
                double bestd2 = xfer_r2;
                double bestscore = 1e300;
                int had_geometric = 0;
                int had_component = 0;
                int32_t vertex_comp = vertex_mesh_comp[vi];
                int require_own_comp =
                    vertex_comp >= 0 && (size_t)vertex_comp < n_mesh_comp &&
                    mesh_sample_count[vertex_comp] > 0;
                for (int d1 = -1; d1 <= 1; d1++) for (int d2i = -1; d2i <= 1; d2i++) {
                    int j1 = i1 + d1, j2 = i2 + d2i;
                    if (j1 < 0 || j2 < 0 || j1 >= g1 || j2 >= g2) continue;
                    size_t cix = (size_t)j1 * (size_t)g2 + (size_t)j2;
                    for (int32_t rr = goff[cix]; rr < goff[cix+1]; rr++) {
                        int32_t is = gids[rr];
                        double dd1 = smp[is].c1 - pc1, dd2 = smp[is].c2 - pc2;
                        double d2 = dd1*dd1 + dd2*dd2;
                        if (d2 >= xfer_r2) continue;
                        had_geometric = 1;
                        if (require_own_comp) {
                            int32_t sch = smp[is].chain;
                            int32_t sample_comp = sch >= 0
                                                ? S->chn[sch].mesh_comp : -1;
                            if (sample_comp != vertex_comp) {
                                component_rejects++;
                                continue;
                            }
                        }
                        had_component = 1;
                        double score = d2;
                        if (o->reference_phi != NULL) {
                            double vphi = (double)o->reference_phi[vi];
                            double dphi = fabs(smp[is].phi - vphi);
                            if (!isfinite(vphi) || !isfinite(smp[is].phi) ||
                                dphi > RIB_XFER_DPHI_MAX) {
                                phase_rejects++;
                                continue;
                            }
                            /* Distance alone aliases self-contacting regions.
                             * Lifted phase supplies the missing universal-cover
                             * direction; r*dphi is its local physical scale.
                             * This is only candidate identity, never material U. */
                            double vr = sqrt(pc1*pc1 + pc2*pc2);
                            double lr = 0.5 * (vr + smp[is].r);
                            if (lr < 1.0) lr = 1.0;
                            double phase_arc = lr * dphi;
                            score += phase_arc * phase_arc;
                        }
                        if (score < bestscore ||
                            (score == bestscore && d2 < bestd2)) {
                            bestscore = score;
                            bestd2 = d2;
                            best = is;
                        }
                    }
                }
                if (best < 0 && had_geometric) {
                    if (require_own_comp && !had_component)
                        component_blocked++;
                    else if (had_component && o->reference_phi != NULL)
                        phase_blocked++;
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
                        u0[vi] = uu; ph0[vi] = sp->phi; gr0[vi] = gg;
                        si0[vi] = best; ok0[vi] = 1;
                    } else {
                        u1[vi] = uu; ph1[vi] = sp->phi; gr1[vi] = gg;
                        si1[vi] = best; ok1[vi] = 1;
                    }
                }
            }
        }
    }

    if (!o->stitch_solve_gauges)
        reconcile_transfer_gauges(arena, S, verts, nv, faces, nf,
                                  o->reference_phi, vw,
                                  u0, u1, ph0, ph1, si0, si1, ok0, ok1,
                                  level_capture, out);
    if (!o->stitch_solve_gauges &&
        o->reference_island != NULL && !o->emit_global)
        pack_metric_islands(arena, S, (double)o->grid_u,
                            u0, u1, si0, si1, ok0, ok1, nv,
                            level_capture, out);
    {
        double lo = 1e300, hi = -1e300;
        for (size_t i = 0; i < S->n_smp; i++) {
            if (S->smp[i].chain < 0) continue;
            if (S->smp[i].u < lo) lo = S->smp[i].u;
            if (S->smp[i].u > hi) hi = S->smp[i].u;
        }
        if (lo < 1e299) out->u_span = hi - lo;
    }

    double *uvert = RIB_ALLOC_ARRAY(arena, double, nv);
    double *phvert = (double *)ARENA_ALLOC(arena, nv * sizeof(double));
    int32_t *gvert = (int32_t *)ARENA_ALLOC(arena, nv * sizeof(int32_t));
    size_t fallback = 0;
    for (size_t i = 0; i < nv; i++) {
        double w = (double)vw[i];
        if (ok0[i] && ok1[i]) {
            int32_t ch0 = si0[i] >= 0 ? smp[si0[i]].chain : -1;
            int32_t ch1 = si1[i] >= 0 ? smp[si1[i]].chain : -1;
            int32_t rc0 = ch0 >= 0
                        ? S->chn[ch0].reconstruction_component : -1;
            int32_t rc1 = ch1 >= 0
                        ? S->chn[ch1].reconstruction_component : -1;
            if (rc0 >= 0 && rc1 >= 0 && rc0 != rc1) {
                /* A branch boundary is a real atlas cut.  Interpolating two
                 * independently packed components would place this vertex in
                 * the empty gutter between them and manufacture a third,
                 * unsupported surface.  Keep the nearer measured side. */
                if (w < 0.5) {
                    uvert[i] = u0[i]; phvert[i] = ph0[i]; gvert[i] = gr0[i];
                } else {
                    uvert[i] = u1[i]; phvert[i] = ph1[i]; gvert[i] = gr1[i];
                }
            } else {
                uvert[i]  = (1.0 - w) * u0[i]  + w * u1[i];
                phvert[i] = (1.0 - w) * ph0[i] + w * ph1[i];
                gvert[i]  = w < 0.5 ? gr0[i] : gr1[i];
            }
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
    out->uv_phase_rejects = phase_rejects;
    out->uv_phase_blocked = phase_blocked;
    out->uv_component_rejects = component_rejects;
    out->uv_component_blocked = component_blocked;

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

static uint64_t rib_digest_mix(uint64_t h, uint64_t x)
{
    for (int q = 0; q < 8; q++) {
        h ^= (uint8_t)(x & UINT64_C(0xff));
        h *= UINT64_C(1099511628211);
        x >>= 8;
    }
    return h;
}

static void rib_grid_digest(const char *width_name, const char *stage,
                            const float *G, const uint8_t *present,
                            size_t nk, size_t stride,
                            size_t active_nu)
{
    uint64_t h = UINT64_C(1469598103934665603);
    size_t finite = 0;
    for (size_t k = 0; k < nk; k++) {
        for (size_t j = 0; j < active_nu; j++) {
            const float *p = &G[(k * stride + j) * 3];
            if (!present[k * stride + j]) continue;
            h = rib_digest_mix(h, (uint64_t)(k * active_nu + j));
            for (int d = 0; d < 3; d++) {
                uint32_t bits = 0;
                memcpy(&bits, &p[d], sizeof(bits));
                h = rib_digest_mix(h, (uint64_t)bits);
            }
            finite++;
        }
    }
    fprintf(stderr,
            "  fit-width digest %-7s %-8s: occupied=%zu hash=%016llx\n",
            width_name, stage, finite, (unsigned long long)h);
}

/* ---- hard chain merges + deterministic fusion ------------------------------
 *
 * A fitted-grid cell can be claimed by several slice chains when source
 * charts overlap in their registered metric gauges.  Overlapping cube charts
 * are alternate measurements of the same sheet, never a license to mosaic
 * them cell by cell: the 2026-08-13 claimant row-DP/tabu search left ~139
 * chart runs per raster row with interior quilting energy equal to the
 * V-boundary energy.  Ownership is therefore decided BY CONSTRUCTION, never
 * by search:
 *
 *   1. rib_merge_build joins slice chains into logical chains at certified
 *      continuations only: same material island, a serialized direct chart
 *      relation (or one chart), carried-U agreement within the layout
 *      certificate, and physical seam reach.  Junctions are mutual-best, so
 *      logical chains are simple paths.
 *   2. rib_fuse_select_rows keeps one claimant per cell using a row-invariant
 *      key (global chart support first), so an overlap resolves to the SAME
 *      chart in every row and the ownership frontier sits exactly at the end
 *      of measured support.  Losing claims are classified and counted; a
 *      loser beyond the wrap gate is a reported certificate violation, not a
 *      per-cell decision.
 *   3. rib_fuse_audit_runs re-derives every frozen ownership transition from
 *      the exported chart runs and classifies it (gap / certified /
 *      uncertified); malformed provenance fails the parameterization.
 *
 * A claimant transition can therefore occur only at (i) the end of measured
 * support, (ii) a real source-sheet boundary or hole, or (iii) a directly
 * certified chart-graph continuation.  No raster or texture signal enters
 * any of this. */
typedef struct {
    double p[3];
    double tau[3];
    double phi;       /* authoritative continuous lifted phase (radians) */
    double phi_score; /* |candidate phi - component phi consensus|; claims
                       * (search) ownership mode only, 0 in construction. */
    int32_t mesh_comp;
    int32_t chain;
    int32_t source_chart;
} RibGridClaim;

static double rib_claim_d2(const RibGridClaim *a, const RibGridClaim *b)
{
    double d0 = a->p[0] - b->p[0];
    double d1 = a->p[1] - b->p[1];
    double d2 = a->p[2] - b->p[2];
    return d0*d0 + d1*d1 + d2*d2;
}

typedef struct {
    size_t first;
    size_t width;
    size_t *off;
    RibGridClaim *claim;
} RibClaimRow;

typedef struct {
    size_t first;
    size_t last;
    int32_t chain;
} RibClaimLabelRun;

typedef struct {
    RibClaimLabelRun *run;
    size_t n;
    size_t capacity;
} RibClaimLabelRow;

static int rib_claim_select_labels(const RibGridClaim *claim,
                                   const size_t *off, size_t width,
                                   const RibClaimLabelRow *labels,
                                   size_t absolute_first,
                                   int32_t *selected,
                                   size_t *out_nonfirst);
static int rib_claim_row_build(Arena_T arena, const SliceSet *S,
                               double U0, double du, size_t nu, size_t k,
                               const int32_t *choff, const int32_t *chids,
                               RibClaimRow *row);

static void rib_claim_label_rows_dispose(RibClaimLabelRow *rows, size_t n)
{
    if (rows == NULL) return;
    for (size_t k = 0; k < n; k++) free(rows[k].run);
    free(rows);
}

static int32_t rib_claim_label_at(const RibClaimLabelRow *row, size_t j)
{
    size_t lo = 0, hi = row != NULL ? row->n : 0;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        const RibClaimLabelRun *run = &row->run[mid];
        if (j < run->first) hi = mid;
        else if (j > run->last) lo = mid + 1;
        else return run->chain;
    }
    return -1;
}

static int rib_claim_label_set(RibClaimLabelRow *row,
                               const RibGridClaim *claim,
                               const int32_t *selected,
                               size_t first, size_t width)
{
    size_t need = 0;
    int32_t previous = -1;
    for (size_t j = 0; j < width; j++) {
        int32_t q = selected[j];
        int32_t chain = q >= 0 ? claim[(size_t)q].chain : -1;
        if (chain >= 0 && chain != previous) need++;
        previous = chain;
    }
    if (need > row->capacity) {
        size_t capacity = row->capacity ? row->capacity : 16;
        while (capacity < need) {
            if (capacity > SIZE_MAX / 2) { capacity = need; break; }
            capacity *= 2;
        }
        if (capacity > SIZE_MAX / sizeof(*row->run)) return -1;
        RibClaimLabelRun *grown = (RibClaimLabelRun *)realloc(
            row->run, capacity * sizeof(*row->run));
        if (grown == NULL) return -1;
        row->run = grown;
        row->capacity = capacity;
    }
    row->n = 0;
    previous = -1;
    for (size_t j = 0; j < width; j++) {
        int32_t q = selected[j];
        int32_t chain = q >= 0 ? claim[(size_t)q].chain : -1;
        size_t absolute = first + j;
        if (chain < 0) {
            previous = -1;
            continue;
        }
        if (chain == previous) {
            row->run[row->n - 1].last = absolute;
        } else {
            RibClaimLabelRun *run = &row->run[row->n++];
            run->first = run->last = absolute;
            run->chain = chain;
        }
        previous = chain;
    }
    return 0;
}

static int rib_claim_label_copy(RibClaimLabelRow *dst,
                                const RibClaimLabelRow *src)
{
    if (src->n > dst->capacity) {
        RibClaimLabelRun *grown;
        if (src->n > SIZE_MAX / sizeof(*dst->run)) return -1;
        grown = (RibClaimLabelRun *)realloc(
            dst->run, src->n * sizeof(*dst->run));
        if (grown == NULL) return -1;
        dst->run = grown;
        dst->capacity = src->n;
    }
    if (src->n != 0)
        memcpy(dst->run, src->run, src->n * sizeof(*dst->run));
    dst->n = src->n;
    return 0;
}

/* Freeze the graph decision in a public, compact source-chart map before the
 * internal chain labels are released.  Adjacent chain runs from the same
 * source chart are merged; an unsupported gap always remains a gap.  This map
 * is diagnostic provenance only, but it also makes the key construction
 * invariant auditable after output compaction: measured cells retain exactly
 * the chart selected by the graph. */
static int rib_claim_export_chart_runs(
    Arena_T arena, const SliceSet *S,
    const RibClaimLabelRow *labels, size_t nk,
    size_t **out_offsets, RibbonChartRun **out_runs, size_t *out_count)
{
    size_t *offsets = NULL;
    RibbonChartRun *runs = NULL;
    size_t count = 0;
    *out_offsets = NULL;
    *out_runs = NULL;
    *out_count = 0;
    if (nk > (SIZE_MAX / sizeof(*offsets)) - 1) return -1;
    offsets = RIB_ALLOC_ARRAY(arena, size_t, nk + 1);
    offsets[0] = 0;
    for (size_t k = 0; k < nk; k++) {
        int32_t previous_chart = -1;
        size_t previous_last = 0;
        for (size_t r = 0; r < labels[k].n; r++) {
            const RibClaimLabelRun *src = &labels[k].run[r];
            int32_t chart;
            if (src->chain < 0 || (size_t)src->chain >= S->n_chn) {
                return -1;
            }
            chart = S->chn[src->chain].source_chart;
            if (chart < 0) {
                return -1;
            }
            if (r == 0 || chart != previous_chart ||
                previous_last == SIZE_MAX || src->first != previous_last + 1) {
                if (count == SIZE_MAX) return -1;
                count++;
            }
            previous_chart = chart;
            previous_last = src->last;
        }
        offsets[k + 1] = count;
    }
    if (count > SIZE_MAX / sizeof(*runs)) return -1;
    runs = RIB_ALLOC_ARRAY(arena, RibbonChartRun, count ? count : 1);
    count = 0;
    for (size_t k = 0; k < nk; k++) {
        for (size_t r = 0; r < labels[k].n; r++) {
            const RibClaimLabelRun *src = &labels[k].run[r];
            int32_t chart = S->chn[src->chain].source_chart;
            if (count > offsets[k] &&
                runs[count - 1].source_chart == chart &&
                runs[count - 1].last_col != SIZE_MAX &&
                src->first == runs[count - 1].last_col + 1) {
                runs[count - 1].last_col = src->last;
            } else {
                runs[count].first_col = src->first;
                runs[count].last_col = src->last;
                runs[count].source_chart = chart;
                count++;
            }
        }
        if (count != offsets[k + 1]) {
            return -1;
        }
    }
    *out_offsets = offsets;
    *out_runs = runs;
    *out_count = count;
    return 0;
}

/* Recover the exact graph-selected claimant by its stable chain identity.
 * This is deliberately not a nearest-position lookup: the graph assignment is
 * the decision, while the raster position is only its payload.  Rebuilding a
 * row from the same SliceSet must therefore reproduce every selected identity
 * exactly or fail the parameterization. */
static int rib_claim_select_labels(const RibGridClaim *claim,
                                   const size_t *off, size_t width,
                                   const RibClaimLabelRow *labels,
                                   size_t absolute_first,
                                   int32_t *selected,
                                   size_t *out_nonfirst)
{
    size_t nonfirst = 0;
    for (size_t j = 0; j < width; j++) {
        int32_t wanted = rib_claim_label_at(labels, absolute_first + j);
        selected[j] = -1;
        if (off[j] == off[j+1]) {
            if (wanted >= 0) return -1;
            continue;
        }
        if (wanted < 0) return -1;
        for (size_t q = off[j]; q < off[j+1]; q++) {
            if (claim[q].chain != wanted) continue;
            if (selected[j] >= 0) return -1;
            selected[j] = (int32_t)q;
        }
        if (selected[j] < 0) return -1;
        if ((size_t)selected[j] != off[j]) nonfirst++;
    }
    if (out_nonfirst != NULL) *out_nonfirst = nonfirst;
    return 0;
}

/* Materialize one slice row's claimant graph in a compact U window.  This is
 * deliberately shared by graph assignment and final row fitting so the graph
 * is solved over exactly the candidates which can reach emitted geometry. */
static int rib_claim_row_build(Arena_T arena, const SliceSet *S,
                               double U0, double du, size_t nu, size_t k,
                               const int32_t *choff, const int32_t *chids,
                               RibClaimRow *row)
{
    size_t rjmin = nu, rjmax = 0;
    memset(row, 0, sizeof(*row));
    row->first = nu;
    for (int32_t x = choff[k]; x < choff[k+1]; x++) {
        const Chain *ch = &S->chn[chids[x]];
        int32_t f = ch->first, cn = ch->count;
        if (cn < 2) continue;
        double ulo = S->smp[f].u - U0;
        double uhi = S->smp[f+cn-1].u - U0;
        int32_t j0 = (int32_t)ceil(ulo / du);
        int32_t j1 = (int32_t)floor(uhi / du);
        if (j0 < 0) j0 = 0;
        if (j1 > (int32_t)nu - 1) j1 = (int32_t)nu - 1;
        if (j0 > j1) continue;
        if ((size_t)j0 < rjmin) rjmin = (size_t)j0;
        if ((size_t)j1 > rjmax) rjmax = (size_t)j1;
    }
    if (rjmin > rjmax) return 0;

    size_t rw = rjmax - rjmin + 1;
    size_t *off = (size_t *)ARENA_CALLOC(arena, rw + 1, sizeof(size_t));
    for (int32_t x = choff[k]; x < choff[k+1]; x++) {
        const Chain *ch = &S->chn[chids[x]];
        int32_t f = ch->first, cn = ch->count;
        if (cn < 2) continue;
        double ulo = S->smp[f].u - U0;
        double uhi = S->smp[f+cn-1].u - U0;
        int32_t j0 = (int32_t)ceil(ulo / du);
        int32_t j1 = (int32_t)floor(uhi / du);
        if (j0 < 0) j0 = 0;
        if (j1 > (int32_t)nu - 1) j1 = (int32_t)nu - 1;
        for (int32_t j = j0; j <= j1; j++)
            off[(size_t)j - rjmin + 1]++;
    }
    for (size_t j = 0; j < rw; j++) {
        if (off[j+1] > SIZE_MAX - off[j]) return -1;
        off[j+1] += off[j];
    }
    if (off[rw] == 0) return 0;

    RibGridClaim *claim = RIB_ALLOC_ARRAY(arena, RibGridClaim, off[rw]);
    size_t *cursor = RIB_ALLOC_ARRAY(arena, size_t, rw);
    memcpy(cursor, off, rw * sizeof(size_t));
    for (int32_t x = choff[k]; x < choff[k+1]; x++) {
        const Chain *ch = &S->chn[chids[x]];
        int32_t f = ch->first, cn = ch->count;
        if (cn < 2) continue;
        double ulo = S->smp[f].u - U0;
        double uhi = S->smp[f+cn-1].u - U0;
        int32_t j0 = (int32_t)ceil(ulo / du);
        int32_t j1 = (int32_t)floor(uhi / du);
        if (j0 < 0) j0 = 0;
        if (j1 > (int32_t)nu - 1) j1 = (int32_t)nu - 1;
        int32_t seg = 1;
        for (int32_t j = j0; j <= j1; j++) {
            double uq = (double)j * du + U0;
            while (seg < cn - 1 && S->smp[f+seg].u < uq) seg++;
            double ua = S->smp[f+seg-1].u;
            double ub = S->smp[f+seg].u;
            double w = ub - ua > 1e-12 ? (uq - ua) / (ub - ua) : 0.0;
            RibGridClaim *dst =
                &claim[cursor[(size_t)j - rjmin]++];
            if (w < 0.0) w = 0.0;
            if (w > 1.0) w = 1.0;
            for (int d = 0; d < 3; d++) {
                dst->p[d] = (1.0-w)*S->smp[f+seg-1].p[d]
                          +       w *S->smp[f+seg].p[d];
                dst->tau[d] = (1.0-w)*S->smp[f+seg-1].tau[d]
                            +       w *S->smp[f+seg].tau[d];
            }
            dst->phi = (1.0-w)*S->smp[f+seg-1].phi
                     +       w *S->smp[f+seg].phi;
            dst->phi_score = 0.0;
            dst->mesh_comp = ch->mesh_comp;
            dst->chain = chids[x];
            dst->source_chart = ch->source_chart;
        }
    }
    row->first = rjmin;
    row->width = rw;
    row->off = off;
    row->claim = claim;
    return 0;
}

/* ---- hard chain merges ---------------------------------------------------- */

typedef struct {
    int32_t *logical;       /* [n_chn] logical-chain label, 0..n_logical-1 */
    double  *chart_support; /* [n_charts] total sample count per source chart;
                             * integer-valued doubles, so comparisons between
                             * them are exact */
    size_t   n_charts;
    size_t   n_logical;
    size_t   n_junctions;
    size_t   ambiguous;
    size_t   rejected_island, rejected_relation, rejected_du,
             rejected_gap, rejected_radial, rejected_phase;
} RibMergeSet;

typedef struct {
    double u0, u1;          /* first/last sample u of the chain */
    int32_t chain;
} RibMergeEnd;

static int cmp_rib_merge_end(const void *pa, const void *pb)
{
    const RibMergeEnd *a = (const RibMergeEnd *)pa;
    const RibMergeEnd *b = (const RibMergeEnd *)pb;
    if (a->u0 < b->u0) return -1;
    if (a->u0 > b->u0) return 1;
    return a->chain < b->chain ? -1 : (a->chain > b->chain ? 1 : 0);
}

/* Join slice chains into logical chains at certified continuations.  Within
 * one slice plane, chain c's END may continue into chain d's START iff every
 * gate passes: same material island; a direct serialized chart relation (one
 * chart is intrinsically allowed); carried-U agreement within the layout
 * certificate (RIB_MERGE_DU_MAX); physical seam reach (RIB_CONT_GAP_VOX,
 * below the 7-vox inter-wrap clearance); the same wrap radially
 * (RIB_CONT_DR_MAX); and, when an authoritative winding field exists, lifted
 * phase (RIB_CONT_DPHI_MAX).  Each chain proposes its single best successor
 * under the strict key (|dU|, 3-D gap, chart, start u, chain id); each
 * proposed-to chain keeps its best proposer under the symmetric key.  One
 * successor and one predecessor per chain make every logical chain a simple
 * path.  Near-tied alternatives are counted as `ambiguous` -- the trigger for
 * carrying explicit port intervals in a future layout version, not a license
 * to guess. */
static int rib_merge_build(Arena_T arena, const SliceSet *S,
                           const RibbonOpts *o, RibMergeSet *m)
{
    size_t nc = S->n_chn;
    /* Fast paths without a reference winding field carry phi == u, where the
     * pi/2 gate would be meaningless; the classic solve and the carried
     * registered field both provide a real lifted phase. */
    int phase_gate = o->reference_u == NULL || o->reference_phi != NULL;
    memset(m, 0, sizeof(*m));
    if (nc == 0 || nc > (size_t)INT32_MAX) return nc == 0 ? 0 : -1;

    m->logical = RIB_ALLOC_ARRAY(arena, int32_t, nc);
    m->n_charts = 0;
    for (size_t c = 0; c < nc; c++) {
        int32_t chart = S->chn[c].source_chart;
        if (chart >= 0 && (size_t)chart + 1 > m->n_charts)
            m->n_charts = (size_t)chart + 1;
    }
    m->chart_support = (double *)ARENA_CALLOC(
        arena, m->n_charts ? m->n_charts : 1, sizeof(double));
    for (size_t c = 0; c < nc; c++) {
        int32_t chart = S->chn[c].source_chart;
        if (chart >= 0 && S->chn[c].count > 0)
            m->chart_support[chart] += (double)S->chn[c].count;
    }

    Arena_Mark scratch = Arena_save(arena);
    UnionFind uf = UF_new(arena, (int32_t)nc);
    int32_t *best_succ = RIB_ALLOC_ARRAY(arena, int32_t, nc);
    int32_t *best_pred = RIB_ALLOC_ARRAY(arena, int32_t, nc);
    double *prop_du = RIB_ALLOC_ARRAY(arena, double, nc);
    double *prop_d3 = RIB_ALLOC_ARRAY(arena, double, nc);
    int32_t *choff = NULL;
    int32_t *chids = chains_by_slice(arena, S, &choff);
    RibMergeEnd *ends = RIB_ALLOC_ARRAY(arena, RibMergeEnd, nc + 1);
    for (size_t c = 0; c < nc; c++) {
        best_succ[c] = -1;
        best_pred[c] = -1;
        prop_du[c] = 0.0;
        prop_d3[c] = 0.0;
    }

    for (int k = 0; k < S->nplanes; k++) {
        size_t np = 0;
        for (int32_t x = choff[k]; x < choff[k + 1]; x++) {
            const Chain *ch = &S->chn[chids[x]];
            if (ch->count < 1) continue;
            ends[np].u0 = S->smp[ch->first].u;
            ends[np].u1 = S->smp[ch->first + ch->count - 1].u;
            ends[np].chain = chids[x];
            np++;
        }
        if (np < 2) continue;
        qsort(ends, np, sizeof(*ends), cmp_rib_merge_end);

        for (size_t i = 0; i < np; i++) {
            const Chain *ci = &S->chn[ends[i].chain];
            const Sample *E = &S->smp[ci->first + ci->count - 1];
            double ue = ends[i].u1;
            double best_du = 0.0, best_d3 = 0.0;
            double second_du = 0.0, second_d3 = 0.0;
            int32_t best = -1, second = -1;
            /* lower bound of the diagnostic window in the u0-sorted order */
            size_t lo = 0, hi = np;
            while (lo < hi) {
                size_t mid = lo + (hi - lo) / 2;
                if (ends[mid].u0 < ue - 2.0 * RIB_MERGE_DU_MAX) lo = mid + 1;
                else hi = mid;
            }
            for (size_t j = lo; j < np; j++) {
                const Chain *cj = NULL;
                const Sample *F = NULL;
                double du = 0.0, d3 = 0.0, dz = 0.0, dy = 0.0, dx = 0.0;
                if (ends[j].u0 > ue + 2.0 * RIB_MERGE_DU_MAX) break;
                if (ends[j].chain == ends[i].chain) continue;
                if (!(ends[j].u1 > ue)) continue;   /* must extend beyond */
                cj = &S->chn[ends[j].chain];
                F = &S->smp[cj->first];
                if (ci->winding_island != cj->winding_island) {
                    m->rejected_island++;
                    continue;
                }
                if (!direct_chart_relation_allowed(
                        o, ci->source_chart, cj->source_chart)) {
                    m->rejected_relation++;
                    continue;
                }
                du = ends[j].u0 - ue;
                if (fabs(du) > RIB_MERGE_DU_MAX) {
                    m->rejected_du++;
                    continue;
                }
                dz = E->p[0] - F->p[0];
                dy = E->p[1] - F->p[1];
                dx = E->p[2] - F->p[2];
                d3 = sqrt(dz * dz + dy * dy + dx * dx);
                if (d3 > RIB_CONT_GAP_VOX) {
                    m->rejected_gap++;
                    continue;
                }
                if (fabs(E->r - F->r) > RIB_CONT_DR_MAX) {
                    m->rejected_radial++;
                    continue;
                }
                if (phase_gate && fabs(E->phi - F->phi) > RIB_CONT_DPHI_MAX) {
                    m->rejected_phase++;
                    continue;
                }
                if (best < 0 ||
                    fabs(du) < fabs(best_du) ||
                    (fabs(du) == fabs(best_du) &&
                     (d3 < best_d3 ||
                      (d3 == best_d3 &&
                       (cj->source_chart <
                            S->chn[best].source_chart ||
                        (cj->source_chart ==
                             S->chn[best].source_chart &&
                         ends[j].chain < best)))))) {
                    second = best;
                    second_du = best_du;
                    second_d3 = best_d3;
                    best = ends[j].chain;
                    best_du = du;
                    best_d3 = d3;
                } else if (second < 0 ||
                           fabs(du) < fabs(second_du) ||
                           (fabs(du) == fabs(second_du) && d3 < second_d3)) {
                    second = ends[j].chain;
                    second_du = du;
                    second_d3 = d3;
                }
            }
            if (best < 0) continue;
            if (second >= 0 &&
                fabs(fabs(second_du) - fabs(best_du)) < 0.5 &&
                fabs(second_d3 - best_d3) < 0.5)
                m->ambiguous++;
            best_succ[ends[i].chain] = best;
            prop_du[ends[i].chain] = fabs(best_du);
            prop_d3[ends[i].chain] = best_d3;
        }
    }

    /* Each proposed-to chain keeps its single best proposer; chain-id
     * ascending order plus the strict key makes this deterministic. */
    for (size_t c = 0; c < nc; c++) {
        int32_t d = best_succ[c];
        int32_t cur = 0;
        if (d < 0) continue;
        cur = best_pred[d];
        if (cur < 0 ||
            prop_du[c] < prop_du[cur] ||
            (prop_du[c] == prop_du[cur] &&
             (prop_d3[c] < prop_d3[cur] ||
              (prop_d3[c] == prop_d3[cur] &&
               (S->chn[c].source_chart < S->chn[cur].source_chart ||
                (S->chn[c].source_chart == S->chn[cur].source_chart &&
                 (int32_t)c < cur))))))
            best_pred[d] = (int32_t)c;
    }
    for (size_t d = 0; d < nc; d++) {
        if (best_pred[d] < 0) continue;
        uf_union(&uf, best_pred[d], (int32_t)d);
        m->n_junctions++;
    }

    /* compact logical labels in chain-id order (deterministic) */
    {
        int32_t *root_label = RIB_ALLOC_ARRAY(arena, int32_t, nc);
        int32_t next = 0;
        for (size_t c = 0; c < nc; c++) root_label[c] = -1;
        for (size_t c = 0; c < nc; c++) {
            int32_t root = uf_find(&uf, (int32_t)c);
            if (root_label[root] < 0) root_label[root] = next++;
            m->logical[c] = root_label[root];
        }
        m->n_logical = (size_t)next;
    }
    Arena_restore(arena, scratch);
    return 0;
}

/* ---- branch-aware reconstruction components ------------------------------
 *
 * A material candidate may still contain a genuine branch.  The welded center
 * of PHerc0139 is the concrete counterexample: two physically distinct slice
 * curves belong to mesh component 0 and continuation island 0, and overlap in
 * solved U.  A single raster slot cannot preserve both.
 *
 * Work at the last representation which still contains every observation:
 * slice chains.  Certified within-row endpoint continuations form logical
 * paths.  For each row, two logical paths conflict when they claim the same U
 * cell farther apart than the physical wrap gate.  Paths inherit a lane from
 * robust cross-row continuation support when possible; a supported lane which
 * is simultaneously occupied by a conflicting path is a branch and starts a
 * new lane.  Unsupported, non-conflicting fragments reuse the lowest lane so
 * ordinary holes do not explode into thousands of artifacts.  The resulting
 * (material candidate, lane) pairs are compact reconstruction-component IDs
 * and are packed into disjoint atlas intervals before rasterization.
 *
 * This is deliberately conservative about topology: a rejoin selects the
 * strongest continuing lane but never unions the two histories.  Downstream
 * quilting may reason about their endpoints; the ribbon never guesses a weld. */
typedef struct {
    int32_t group;
    double prior_support;
    size_t samples;
} RibBranchOrder;

static int cmp_rib_branch_order(const void *pa, const void *pb)
{
    const RibBranchOrder *a = (const RibBranchOrder *)pa;
    const RibBranchOrder *b = (const RibBranchOrder *)pb;
    if (a->prior_support > b->prior_support) return -1;
    if (a->prior_support < b->prior_support) return 1;
    if (a->samples > b->samples) return -1;
    if (a->samples < b->samples) return 1;
    return a->group < b->group ? -1 : (a->group > b->group ? 1 : 0);
}

static int rib_branch_lane_allowed(int32_t lane, size_t local_index,
                                   const int32_t *row_groups,
                                   size_t row_group_count,
                                   const int32_t *group_lane,
                                   const uint8_t *conflict_matrix)
{
    for (size_t other = 0; other < row_group_count; other++) {
        int32_t group = row_groups[other];
        if (other == local_index || group_lane[group] != lane) continue;
        if (conflict_matrix[local_index * row_group_count + other]) return 0;
    }
    return 1;
}

static int rib_assign_reconstruction_components(
    Arena_T arena, SliceSet *S, const RibbonOpts *o,
    const ChainRelationGraph *relations, RibbonResult *out)
{
    const double gate2 = RIB_WRAP_GATE * RIB_WRAP_GATE;
    size_t nc = S->n_chn;
    int32_t max_material = -1;
    Arena_Mark mark;
    RibMergeSet merge;
    int32_t *group_slice = NULL, *group_material = NULL;
    int32_t *group_lane = NULL, *chain_lane = NULL;
    size_t *group_samples = NULL, *group_count = NULL, *group_off = NULL;
    int32_t *group_chains = NULL;
    size_t *relation_off = NULL, *relation_cursor = NULL;
    int32_t *relation_other = NULL;
    double *relation_support = NULL;
    int32_t *chain_slice_off = NULL, *chains_by_row = NULL;
    size_t *row_group_off = NULL, *row_group_cursor = NULL;
    int32_t *row_groups = NULL, *group_local = NULL;
    int32_t *lane_count = NULL, *lane_base = NULL;
    double *lane_score = NULL;
    uint32_t *lane_stamp = NULL;
    uint32_t stamp = 0;
    size_t relation_entries = 0, far_conflicts = 0, branch_splits = 0;
    size_t relation_cuts = 0;
    double relation_cut_support = 0.0;
    int rc = -1;

    for (size_t c = 0; c < nc; c++) {
        S->chn[c].reconstruction_component = S->chn[c].winding_island;
        if (S->chn[c].winding_island > max_material)
            max_material = S->chn[c].winding_island;
    }
    if (nc == 0 || max_material < 0 || !o->component_global) return 0;
    if (nc > (size_t)INT32_MAX) return -1;

    mark = Arena_save(arena);
    memset(&merge, 0, sizeof(merge));
    if (rib_merge_build(arena, S, o, &merge) != 0 ||
        merge.n_logical == 0 || merge.n_logical > (size_t)INT32_MAX)
        goto done;

    group_slice = RIB_ALLOC_ARRAY(arena, int32_t, merge.n_logical);
    group_material = RIB_ALLOC_ARRAY(arena, int32_t, merge.n_logical);
    group_lane = RIB_ALLOC_ARRAY(arena, int32_t, merge.n_logical);
    group_samples = (size_t *)ARENA_CALLOC(
        arena, merge.n_logical, sizeof(*group_samples));
    group_count = (size_t *)ARENA_CALLOC(
        arena, merge.n_logical + 1, sizeof(*group_count));
    for (size_t g = 0; g < merge.n_logical; g++) {
        group_slice[g] = -1;
        group_material[g] = -1;
        group_lane[g] = -1;
    }
    for (size_t c = 0; c < nc; c++) {
        int32_t g = merge.logical[c];
        const Chain *chain = &S->chn[c];
        if (g < 0 || (size_t)g >= merge.n_logical) goto done;
        if (group_slice[g] < 0) {
            group_slice[g] = chain->slice;
            group_material[g] = chain->winding_island;
        } else if (group_slice[g] != chain->slice ||
                   group_material[g] != chain->winding_island) {
            fprintf(stderr,
                    "ribbon: logical path %d crosses slice/material identity\n",
                    g);
            goto done;
        }
        group_samples[g] += chain->count > 0 ? (size_t)chain->count : 0;
        group_count[(size_t)g + 1]++;
    }
    for (size_t g = 0; g < merge.n_logical; g++)
        group_count[g + 1] += group_count[g];
    group_off = group_count;
    group_chains = RIB_ALLOC_ARRAY(arena, int32_t, nc);
    {
        size_t *cursor = RIB_ALLOC_ARRAY(arena, size_t, merge.n_logical);
        memcpy(cursor, group_off, merge.n_logical * sizeof(*cursor));
        for (size_t c = 0; c < nc; c++)
            group_chains[cursor[merge.logical[c]]++] = (int32_t)c;
    }

    /* Chain adjacency for robust continuation support. */
    relation_off = (size_t *)ARENA_CALLOC(arena, nc + 1,
                                           sizeof(*relation_off));
    if (relations != NULL && relations->key != NULL) {
        for (size_t slot = 0; slot < relations->capacity; slot++) {
            uint64_t key = relations->key[slot];
            int32_t a, b;
            if (key == UINT64_MAX || !(relations->support[slot] > 0.0))
                continue;
            a = (int32_t)(uint32_t)(key >> 32);
            b = (int32_t)(uint32_t)key;
            if (a < 0 || b < 0 || (size_t)a >= nc || (size_t)b >= nc ||
                a == b || S->chn[a].winding_island !=
                          S->chn[b].winding_island)
                continue;
            relation_off[(size_t)a + 1]++;
            relation_off[(size_t)b + 1]++;
            relation_entries += 2;
        }
    }
    for (size_t c = 0; c < nc; c++)
        relation_off[c + 1] += relation_off[c];
    relation_other = RIB_ALLOC_ARRAY(
        arena, int32_t, relation_entries ? relation_entries : 1);
    relation_support = RIB_ALLOC_ARRAY(
        arena, double, relation_entries ? relation_entries : 1);
    relation_cursor = RIB_ALLOC_ARRAY(arena, size_t, nc);
    memcpy(relation_cursor, relation_off, nc * sizeof(*relation_cursor));
    if (relations != NULL && relations->key != NULL) {
        for (size_t slot = 0; slot < relations->capacity; slot++) {
            uint64_t key = relations->key[slot];
            int32_t a, b;
            double support;
            size_t qa, qb;
            if (key == UINT64_MAX || !(relations->support[slot] > 0.0))
                continue;
            a = (int32_t)(uint32_t)(key >> 32);
            b = (int32_t)(uint32_t)key;
            support = relations->support[slot];
            if (a < 0 || b < 0 || (size_t)a >= nc || (size_t)b >= nc ||
                a == b || S->chn[a].winding_island !=
                          S->chn[b].winding_island)
                continue;
            qa = relation_cursor[a]++;
            qb = relation_cursor[b]++;
            relation_other[qa] = b;
            relation_support[qa] = support;
            relation_other[qb] = a;
            relation_support[qb] = support;
        }
    }

    chain_lane = RIB_ALLOC_ARRAY(arena, int32_t, nc);
    for (size_t c = 0; c < nc; c++) chain_lane[c] = -1;
    chains_by_row = chains_by_slice(arena, S, &chain_slice_off);

    row_group_off = (size_t *)ARENA_CALLOC(
        arena, (size_t)S->nplanes + 1, sizeof(*row_group_off));
    for (size_t g = 0; g < merge.n_logical; g++) {
        int32_t slice = group_slice[g];
        if (slice < 0 || slice >= S->nplanes) goto done;
        row_group_off[(size_t)slice + 1]++;
    }
    for (int k = 0; k < S->nplanes; k++)
        row_group_off[(size_t)k + 1] += row_group_off[k];
    row_groups = RIB_ALLOC_ARRAY(arena, int32_t, merge.n_logical);
    row_group_cursor = RIB_ALLOC_ARRAY(arena, size_t, (size_t)S->nplanes);
    memcpy(row_group_cursor, row_group_off,
           (size_t)S->nplanes * sizeof(*row_group_cursor));
    for (size_t g = 0; g < merge.n_logical; g++)
        row_groups[row_group_cursor[group_slice[g]]++] = (int32_t)g;

    group_local = RIB_ALLOC_ARRAY(arena, int32_t, merge.n_logical);
    for (size_t g = 0; g < merge.n_logical; g++) group_local[g] = -1;
    lane_count = (int32_t *)ARENA_CALLOC(
        arena, (size_t)max_material + 1, sizeof(*lane_count));
    lane_score = (double *)ARENA_CALLOC(arena, nc, sizeof(*lane_score));
    lane_stamp = (uint32_t *)ARENA_CALLOC(arena, nc, sizeof(*lane_stamp));

    {
        double umin = 1e300, umax = -1e300;
        double du = (double)o->grid_u;
        size_t cover_nu;
        if (!(du > 1e-6) || !isfinite(du)) du = 1.0;
        for (size_t i = 0; i < S->n_smp; i++) {
            if (S->smp[i].chain < 0 || !isfinite(S->smp[i].u)) continue;
            if (S->smp[i].u < umin) umin = S->smp[i].u;
            if (S->smp[i].u > umax) umax = S->smp[i].u;
        }
        if (umin >= 1e299 || !(umax >= umin)) goto done;
        {
            double cells = ceil((umax - umin) / du) + 2.0;
            if (!isfinite(cells) || cells > (double)SIZE_MAX) goto done;
            cover_nu = (size_t)cells;
            if (cover_nu < 2) cover_nu = 2;
        }

        for (int k = 0; k < S->nplanes; k++) {
            size_t rb = row_group_off[k], re = row_group_off[k + 1];
            size_t rn = re - rb;
            Arena_Mark row_mark;
            uint8_t *conflict_matrix;
            RibBranchOrder *order;
            RibClaimRow claims;
            if (rn == 0) continue;
            if (rn > SIZE_MAX / rn) goto done;
            row_mark = Arena_save(arena);
            conflict_matrix = (uint8_t *)ARENA_CALLOC(arena, rn * rn, 1);
            order = RIB_ALLOC_ARRAY(arena, RibBranchOrder, rn);
            for (size_t local = 0; local < rn; local++)
                group_local[row_groups[rb + local]] = (int32_t)local;

            if (rib_claim_row_build(arena, S, umin, du, cover_nu,
                                    (size_t)k, chain_slice_off,
                                    chains_by_row, &claims) != 0) {
                Arena_restore(arena, row_mark);
                goto done;
            }
            for (size_t j = 0; j < claims.width; j++) {
                for (size_t qa = claims.off[j]; qa < claims.off[j + 1]; qa++) {
                    int32_t ca = claims.claim[qa].chain;
                    int32_t ga = merge.logical[ca];
                    int32_t la = group_local[ga];
                    for (size_t qb = qa + 1; qb < claims.off[j + 1]; qb++) {
                        int32_t cb = claims.claim[qb].chain;
                        int32_t gb = merge.logical[cb];
                        int32_t lb = group_local[gb];
                        double d2;
                        if (ga == gb || la < 0 || lb < 0 ||
                            group_material[ga] != group_material[gb])
                            continue;
                        d2 = rib_claim_d2(&claims.claim[qa],
                                          &claims.claim[qb]);
                        if (d2 <= gate2) continue;
                        conflict_matrix[(size_t)la * rn + (size_t)lb] = 1;
                        conflict_matrix[(size_t)lb * rn + (size_t)la] = 1;
                    }
                }
            }
            for (size_t a = 0; a < rn; a++)
                for (size_t b = a + 1; b < rn; b++)
                    if (conflict_matrix[a * rn + b]) far_conflicts++;

            for (size_t local = 0; local < rn; local++) {
                int32_t g = row_groups[rb + local];
                double prior = 0.0;
                for (size_t q = group_off[g]; q < group_off[g + 1]; q++) {
                    int32_t chain = group_chains[q];
                    for (size_t e = relation_off[chain];
                         e < relation_off[(size_t)chain + 1]; e++) {
                        int32_t other = relation_other[e];
                        if (S->chn[other].slice >= k || chain_lane[other] < 0)
                            continue;
                        prior += relation_support[e];
                    }
                }
                order[local].group = g;
                order[local].prior_support = prior;
                order[local].samples = group_samples[g];
            }
            qsort(order, rn, sizeof(*order), cmp_rib_branch_order);

            for (size_t oi = 0; oi < rn; oi++) {
                int32_t g = order[oi].group;
                int32_t material = group_material[g];
                size_t local = (size_t)group_local[g];
                int32_t best_lane = -1;
                double best_support = -1.0;
                int had_prior = 0;
                if (++stamp == 0) {
                    memset(lane_stamp, 0, nc * sizeof(*lane_stamp));
                    stamp = 1;
                }
                for (size_t q = group_off[g]; q < group_off[g + 1]; q++) {
                    int32_t chain = group_chains[q];
                    for (size_t e = relation_off[chain];
                         e < relation_off[(size_t)chain + 1]; e++) {
                        int32_t other = relation_other[e];
                        int32_t lane;
                        if (S->chn[other].slice >= k ||
                            S->chn[other].winding_island != material)
                            continue;
                        lane = chain_lane[other];
                        if (lane < 0) continue;
                        had_prior = 1;
                        if (lane_stamp[lane] != stamp) {
                            lane_stamp[lane] = stamp;
                            lane_score[lane] = 0.0;
                        }
                        lane_score[lane] += relation_support[e];
                    }
                }
                for (int32_t lane = 0; lane < lane_count[material]; lane++) {
                    double support = lane_stamp[lane] == stamp
                                   ? lane_score[lane] : 0.0;
                    if (!(support > 0.0) ||
                        !rib_branch_lane_allowed(
                            lane, local, &row_groups[rb], rn,
                            group_lane, conflict_matrix))
                        continue;
                    if (best_lane < 0 || support > best_support ||
                        (support == best_support && lane < best_lane)) {
                        best_lane = lane;
                        best_support = support;
                    }
                }
                if (best_lane < 0 && !had_prior) {
                    for (int32_t lane = 0; lane < lane_count[material]; lane++) {
                        if (rib_branch_lane_allowed(
                                lane, local, &row_groups[rb], rn,
                                group_lane, conflict_matrix)) {
                            best_lane = lane;
                            break;
                        }
                    }
                }
                if (best_lane < 0) {
                    if (lane_count[material] == INT32_MAX) {
                        Arena_restore(arena, row_mark);
                        goto done;
                    }
                    best_lane = lane_count[material]++;
                    if (had_prior) branch_splits++;
                }
                group_lane[g] = best_lane;
                for (size_t q = group_off[g]; q < group_off[g + 1]; q++)
                    chain_lane[group_chains[q]] = best_lane;
            }

            for (size_t local = 0; local < rn; local++)
                group_local[row_groups[rb + local]] = -1;
            Arena_restore(arena, row_mark);
        }
    }

    lane_base = RIB_ALLOC_ARRAY(
        arena, int32_t, (size_t)max_material + 2);
    lane_base[0] = 0;
    for (int32_t material = 0; material <= max_material; material++) {
        int64_t next = (int64_t)lane_base[material] + lane_count[material];
        if (next > INT32_MAX) goto done;
        lane_base[material + 1] = (int32_t)next;
    }
    for (size_t c = 0; c < nc; c++) {
        int32_t material = S->chn[c].winding_island;
        if (material < 0 || chain_lane[c] < 0) goto done;
        S->chn[c].reconstruction_component =
            lane_base[material] + chain_lane[c];
    }

    if (relations != NULL && relations->key != NULL) {
        for (size_t slot = 0; slot < relations->capacity; slot++) {
            uint64_t key = relations->key[slot];
            int32_t a, b;
            if (key == UINT64_MAX || !(relations->support[slot] > 0.0))
                continue;
            a = (int32_t)(uint32_t)(key >> 32);
            b = (int32_t)(uint32_t)key;
            if (a < 0 || b < 0 || (size_t)a >= nc || (size_t)b >= nc)
                continue;
            if (S->chn[a].reconstruction_component !=
                S->chn[b].reconstruction_component) {
                relation_cuts++;
                relation_cut_support += relations->support[slot];
            }
        }
    }

    out->grid_reconstruction_components = (size_t)lane_base[max_material + 1];
    out->grid_branch_conflicts = far_conflicts;
    out->grid_branch_splits = branch_splits;
    out->grid_branch_relation_cuts = relation_cuts;
    out->grid_branch_relation_cut_support = relation_cut_support;
    fprintf(stderr,
            "  reconstruction lanes: material=%d output=%zu "
            "far-conflicts=%zu branch-splits=%zu relation-cuts=%zu "
            "(support=%.1f)\n",
            max_material + 1, out->grid_reconstruction_components,
            far_conflicts, branch_splits, relation_cuts,
            relation_cut_support);
    if (max_material < 32) {
        for (int32_t material = 0; material <= max_material; material++)
            if (lane_count[material] > 1)
                fprintf(stderr, "    material %d: %d branch lane(s)\n",
                        material, lane_count[material]);
    }
    rc = 0;

done:
    Arena_restore(arena, mark);
    return rc;
}

/* ---- deterministic per-cell fusion ----------------------------------------
 * The winner key is ROW-INVARIANT: global chart support first, then chart id,
 * logical chain, chain length, and finally intrinsic claim data.  An overlap
 * therefore resolves to the same chart in every row, which is exactly what
 * removes the per-row mosaic: the ownership frontier can only sit where
 * measured support begins or ends. */
static int rib_fuse_claim_better(const SliceSet *S, const RibMergeSet *m,
                                 const RibGridClaim *a, const RibGridClaim *b)
{
    double sa = 0.0, sb = 0.0;
    int32_t la = 0, lb = 0, ca = 0, cb = 0;
    if (a->source_chart >= 0 && (size_t)a->source_chart < m->n_charts)
        sa = m->chart_support[a->source_chart];
    if (b->source_chart >= 0 && (size_t)b->source_chart < m->n_charts)
        sb = m->chart_support[b->source_chart];
    if (sa != sb) return sa > sb;   /* exact: integer-valued doubles */
    if (a->source_chart != b->source_chart)
        return a->source_chart < b->source_chart;
    la = m->logical[a->chain];
    lb = m->logical[b->chain];
    if (la != lb) return la < lb;
    ca = S->chn[a->chain].count;
    cb = S->chn[b->chain].count;
    if (ca != cb) return ca > cb;
    for (int d = 0; d < 3; d++) {
        if (a->p[d] < b->p[d]) return 1;
        if (a->p[d] > b->p[d]) return 0;
    }
    return a->chain < b->chain;
}

/* Freeze one claimant per supported cell.  Rows are independent pure
 * functions of (frozen merge partition, intrinsic claim data), so the
 * OpenMP schedule cannot change any label.  Losing claims are classified:
 * same logical chain (certified junction overlap), alternate measurement of
 * the same sheet (within the wrap gate), or a reported certificate
 * violation.  Nothing is optimized. */
static int rib_fuse_select_rows(const SliceSet *S, double U0, double du,
                                size_t nu, size_t nk,
                                const int32_t *choff, const int32_t *chids,
                                const RibMergeSet *merge,
                                RibbonResult *out,
                                RibClaimLabelRow **out_labels)
{
    const double gate2 = RIB_WRAP_GATE * RIB_WRAP_GATE;
    RibClaimLabelRow *labels =
        (RibClaimLabelRow *)calloc(nk ? nk : 1, sizeof(*labels));
    size_t fuse_same = 0, fuse_alt = 0, fuse_violations = 0;
    size_t printed = 0;
    double violation_max2 = 0.0;
    int failed = 0;
    *out_labels = NULL;
    if (labels == NULL) return -1;
#ifdef _OPENMP
#pragma omp parallel reduction(+:fuse_same) reduction(+:fuse_alt) \
    reduction(+:fuse_violations) reduction(|:failed)
#endif
    {
        Arena_T tarena = Arena_new();
        int64_t pk = 0;
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 16)
#endif
        for (pk = 0; pk < (int64_t)nk; pk++) {
            size_t k = (size_t)pk;
            Arena_Mark rm = Arena_save(tarena);
            RibClaimRow row;
            int32_t *selected = NULL;
            if (rib_claim_row_build(tarena, S, U0, du, nu, k,
                                    choff, chids, &row) != 0) {
                failed = 1;
                Arena_restore(tarena, rm);
                continue;
            }
            if (row.width == 0) {
                Arena_restore(tarena, rm);
                continue;
            }
            selected = RIB_ALLOC_ARRAY(tarena, int32_t, row.width);
            for (size_t jw = 0; jw < row.width; jw++) {
                size_t begin = row.off[jw], end = row.off[jw + 1];
                size_t win = begin;
                selected[jw] = -1;
                if (begin == end) continue;
                for (size_t q = begin + 1; q < end; q++)
                    if (rib_fuse_claim_better(S, merge, &row.claim[q],
                                              &row.claim[win]))
                        win = q;
                selected[jw] = (int32_t)win;
                for (size_t q = begin; q < end; q++) {
                    double d2 = 0.0;
                    if (q == win) continue;
                    if (merge->logical[row.claim[q].chain] ==
                        merge->logical[row.claim[win].chain]) {
                        fuse_same++;
                        continue;
                    }
                    d2 = rib_claim_d2(&row.claim[q], &row.claim[win]);
                    if (d2 <= gate2) {
                        fuse_alt++;
                        continue;
                    }
                    fuse_violations++;
#ifdef _OPENMP
#pragma omp critical(rib_fuse_report)
#endif
                    {
                        if (d2 > violation_max2) violation_max2 = d2;
                        if (printed < 20) {
                            fprintf(stderr,
                                "    [fit] fusion violation: row=%zu "
                                "col=%zu chart %d (chain %d) kept over "
                                "chart %d (chain %d), step=%.2f vox\n",
                                k, row.first + jw,
                                row.claim[win].source_chart,
                                row.claim[win].chain,
                                row.claim[q].source_chart,
                                row.claim[q].chain, sqrt(d2));
                            printed++;
                        }
                    }
                }
            }
            if (rib_claim_label_set(&labels[k], row.claim, selected,
                                    row.first, row.width) != 0)
                failed = 1;
            Arena_restore(tarena, rm);
        }
        Arena_dispose(&tarena);
    }
    if (failed) {
        rib_claim_label_rows_dispose(labels, nk);
        return -1;
    }
    out->grid_fuse_same_logical += fuse_same;
    out->grid_fuse_alternate += fuse_alt;
    out->grid_fuse_violations += fuse_violations;
    if (sqrt(violation_max2) > out->grid_fuse_violation_max)
        out->grid_fuse_violation_max = sqrt(violation_max2);
    *out_labels = labels;
    return 0;
}

/* Re-derive every frozen ownership transition from the exported chart runs.
 * A contiguous different-chart transition must be a certified chart-graph
 * continuation to count as certified; everything else is reported.  Only a
 * structurally malformed export (overlapping, unsorted, or unmerged runs)
 * fails the parameterization -- physical protection against wrap-spanning
 * geometry already lives in the emission gates. */
static int rib_fuse_audit_runs(const RibbonOpts *o, RibbonResult *out,
                               size_t nk)
{
    size_t gap = 0, certified = 0, uncertified = 0;
    const size_t *offs = out->grid_chart_row_offsets;
    const RibbonChartRun *runs = out->grid_chart_runs;
    if (offs == NULL) return 0;
    for (size_t k = 0; k < nk; k++) {
        for (size_t r = offs[k]; r < offs[k + 1]; r++) {
            if (runs[r].first_col > runs[r].last_col) return -1;
            if (r == offs[k]) continue;
            if (runs[r - 1].last_col >= runs[r].first_col) return -1;
            if (runs[r - 1].last_col + 1 == runs[r].first_col) {
                if (runs[r - 1].source_chart == runs[r].source_chart)
                    return -1;  /* export merges these; malformed */
                if (direct_chart_relation_allowed(
                        o, runs[r - 1].source_chart, runs[r].source_chart))
                    certified++;
                else
                    uncertified++;
            } else {
                gap++;
            }
        }
    }
    out->grid_transition_gap += gap;
    out->grid_transition_certified += certified;
    out->grid_transition_uncertified += uncertified;
    fprintf(stderr,
            "    [fit] ownership transitions: gap=%zu certified=%zu "
            "uncertified=%zu\n", gap, certified, uncertified);
    return 0;
}

static double rib_grid_d2(const float *a, const float *b)
{
    double d0 = (double)a[0] - (double)b[0];
    double d1 = (double)a[1] - (double)b[1];
    double d2 = (double)a[2] - (double)b[2];
    return d0*d0 + d1*d1 + d2*d2;
}

static size_t rib_grid_long_edges(const float *grid, const uint8_t *present,
                                  size_t nk, size_t nu)
{
    const double gate2 = RIB_WRAP_GATE * RIB_WRAP_GATE;
    size_t count = 0;
    for (size_t k = 0; k < nk; k++) {
        for (size_t j = 0; j < nu; j++) {
            const float *p = &grid[(k*nu+j)*3];
            if (!present[k*nu+j]) continue;
            if (j + 1 < nu) {
                const float *q = p + 3;
                if (present[k*nu+j+1] && rib_grid_d2(p,q) > gate2)
                    count++;
            }
            if (k + 1 < nk) {
                const float *q = p + nu*3;
                if (present[(k+1)*nu+j] && rib_grid_d2(p,q) > gate2)
                    count++;
            }
        }
    }
    return count;
}

static void rib_grid_long_edge_directions(const float *grid,
                                          const uint8_t *present,
                                          size_t nk, size_t nu,
                                          size_t *horizontal,
                                          size_t *vertical)
{
    const double gate2 = RIB_WRAP_GATE * RIB_WRAP_GATE;
    size_t h = 0, v = 0;
    for (size_t k = 0; k < nk; k++) {
        for (size_t j = 0; j < nu; j++) {
            const float *p = &grid[(k*nu+j)*3];
            if (!present[k*nu+j]) continue;
            if (j + 1 < nu) {
                const float *q = p + 3;
                if (present[k*nu+j+1] && rib_grid_d2(p,q) > gate2) h++;
            }
            if (k + 1 < nk) {
                const float *q = p + nu*3;
                if (present[(k+1)*nu+j] && rib_grid_d2(p,q) > gate2) v++;
            }
        }
    }
    if (horizontal) *horizontal = h;
    if (vertical) *vertical = v;
}

static double rib_d2f(const float *a, const float *b)
{
    double d0 = (double)a[0] - (double)b[0];
    double d1 = (double)a[1] - (double)b[1];
    double d2 = (double)a[2] - (double)b[2];
    return d0*d0 + d1*d1 + d2*d2;
}

static int rib_grid_triangle_legal(const double a[3], const double b[3],
                                   const double c[3])
{
    const double gate2 = RIB_WRAP_GATE * RIB_WRAP_GATE;
    double ab0 = a[0]-b[0], ab1 = a[1]-b[1], ab2 = a[2]-b[2];
    double ac0 = a[0]-c[0], ac1 = a[1]-c[1], ac2 = a[2]-c[2];
    double bc0 = b[0]-c[0], bc1 = b[1]-c[1], bc2 = b[2]-c[2];
    return ab0*ab0+ab1*ab1+ab2*ab2 <= gate2 &&
           ac0*ac0+ac1*ac1+ac2*ac2 <= gate2 &&
           bc0*bc0+bc1*bc1+bc2*bc2 <= gate2;
}

/* Number of triangles the final adaptive emitter can retain in one cell.
 * `override_slot` substitutes one candidate claim without mutating the grid. */
static int rib_grid_cell_legal(const float *grid, size_t nu,
                               size_t k, size_t j, size_t override_slot,
                               const double override_p[3])
{
    size_t slot[4] = { k*nu+j, k*nu+j+1,
                       (k+1)*nu+j, (k+1)*nu+j+1 };
    double p[4][3];
    int valid[4], nv = 0;
    for (int q = 0; q < 4; q++) {
        const float *g = &grid[slot[q]*3];
        valid[q] = slot[q] == override_slot || isfinite((double)g[0]);
        if (!valid[q]) continue;
        for (int d = 0; d < 3; d++)
            p[q][d] = slot[q] == override_slot
                    ? override_p[d] : (double)g[d];
        nv++;
    }
    if (nv < 3) return 0;
    if (nv == 3) {
        int q[3], n = 0;
        for (int x = 0; x < 4; x++) if (valid[x]) q[n++] = x;
        return rib_grid_triangle_legal(p[q[0]],p[q[1]],p[q[2]]);
    }
    {
        int current = rib_grid_triangle_legal(p[0],p[1],p[2]) +
                      rib_grid_triangle_legal(p[1],p[3],p[2]);
        int alternate = rib_grid_triangle_legal(p[0],p[1],p[3]) +
                        rib_grid_triangle_legal(p[0],p[3],p[2]);
        return current > alternate ? current : alternate;
    }
}

static int rib_grid_local_legal(const float *grid, size_t nk, size_t nu,
                                size_t k, size_t j,
                                const double candidate[3])
{
    int legal = 0;
    size_t slot = k*nu+j;
    if (nk < 2 || nu < 2) return 0;
    size_t k0 = k > 0 ? k-1 : k, k1 = k+1 < nk ? k : k-1;
    size_t j0 = j > 0 ? j-1 : j, j1 = j+1 < nu ? j : j-1;
    for (size_t ck = k0; ck <= k1; ck++)
        for (size_t cj = j0; cj <= j1; cj++)
            legal += rib_grid_cell_legal(grid,nu,ck,cj,slot,candidate);
    return legal;
}

static size_t rib_grid_incident_long(const float *grid, size_t nk, size_t nu,
                                     size_t k, size_t j,
                                     const double candidate[3])
{
    const double gate2 = RIB_WRAP_GATE * RIB_WRAP_GATE;
    size_t bad = 0;
    const long dk[4] = { 0, 0, -1, 1 };
    const long dj[4] = { -1, 1, 0, 0 };
    for (int q = 0; q < 4; q++) {
        long kk = (long)k + dk[q], jj = (long)j + dj[q];
        if (kk < 0 || jj < 0 || (size_t)kk >= nk || (size_t)jj >= nu)
            continue;
        const float *n = &grid[((size_t)kk*nu+(size_t)jj)*3];
        if (!isfinite((double)n[0])) continue;
        double d0=candidate[0]-(double)n[0];
        double d1=candidate[1]-(double)n[1];
        double d2=candidate[2]-(double)n[2];
        if (d0*d0+d1*d1+d2*d2 > gate2) bad++;
    }
    return bad;
}

static int rib_grid_local_legal_current(const float *grid,
                                        size_t nk, size_t nu,
                                        size_t k, size_t j)
{
    int legal = 0;
    if (nk < 2 || nu < 2) return 0;
    size_t k0 = k > 0 ? k-1 : k, k1 = k+1 < nk ? k : k-1;
    size_t j0 = j > 0 ? j-1 : j, j1 = j+1 < nu ? j : j-1;
    for (size_t ck = k0; ck <= k1; ck++)
        for (size_t cj = j0; cj <= j1; cj++)
            legal += rib_grid_cell_legal(
                grid,nu,ck,cj,SIZE_MAX,NULL);
    return legal;
}

typedef struct {
    size_t slot;
    float p[3];
    int supported;
} RibGridLocalRepair;

/* Repair a very specific fitted-grid defect at its geometric source.  If the
 * four axial neighbours of a slot form one tight, already-emittable boundary
 * diamond, but the slot itself is a physical outlier that suppresses at least
 * two local triangles, replace it by their bilinear centre.  The replacement
 * is accepted only when every new edge passes the normal six-voxel inter-wrap
 * gate and the exact local triangle count strictly increases.
 *
 * This is not a generic hole cap: it changes the bad grid parameterization
 * before triangulation.  A larger bay, an incomplete ring, or a ring whose
 * diameter reaches another ply remains untouched. */
static size_t rib_grid_repair_isolated_outliers(
    float *grid, uint8_t *grid_supported, size_t nk, size_t nu,
    size_t *out_supported)
{
    const double gate2 = RIB_WRAP_GATE * RIB_WRAP_GATE;
    const double ring_diameter2 = 6.5 * 6.5;
    RibGridLocalRepair *repair = NULL;
    size_t nrepair = 0, capacity = 0, supported = 0;

    if (out_supported) *out_supported = 0;
    if (!grid || nk < 3 || nu < 3) return 0;

    for (size_t k = 1; k + 1 < nk; k++) {
        for (size_t j = 1; j + 1 < nu; j++) {
            size_t slot = k*nu+j;
            size_t nbr[4] = {
                (k-1)*nu+j, k*nu+j+1, (k+1)*nu+j, k*nu+j-1
            };
            double centre[3] = { 0.0, 0.0, 0.0 };
            double old[3];
            int ring_ok = 1;
            int before, after;
            size_t old_bad, new_bad;

            if (!isfinite((double)grid[slot*3])) continue;
            for (int q = 0; q < 4; q++) {
                const float *p = &grid[nbr[q]*3];
                if (!isfinite((double)p[0])) { ring_ok = 0; break; }
                for (int d = 0; d < 3; d++) centre[d] += 0.25*(double)p[d];
            }
            if (!ring_ok) continue;

            /* All four ring vertices must belong to one compact physical
             * neighbourhood, not merely four pairwise boundary edges. */
            for (int a = 0; a < 4 && ring_ok; a++) {
                const float *pa = &grid[nbr[a]*3];
                for (int b = a + 1; b < 4; b++) {
                    if (rib_grid_d2(pa,&grid[nbr[b]*3]) >
                        ring_diameter2) {
                        ring_ok = 0;
                        break;
                    }
                }
            }
            if (!ring_ok) continue;
            for (int q = 0; q < 4; q++) {
                const float *p = &grid[nbr[q]*3];
                double d0=centre[0]-(double)p[0];
                double d1=centre[1]-(double)p[1];
                double d2=centre[2]-(double)p[2];
                if (d0*d0+d1*d1+d2*d2 > gate2) {
                    ring_ok = 0;
                    break;
                }
            }
            if (!ring_ok) continue;

            for (int d = 0; d < 3; d++)
                old[d] = (double)grid[slot*3+(size_t)d];
            old_bad = rib_grid_incident_long(
                grid,nk,nu,k,j,old);
            if (old_bad < 2) continue;
            new_bad = rib_grid_incident_long(
                grid,nk,nu,k,j,centre);
            if (new_bad != 0) continue;
            before = rib_grid_local_legal_current(grid,nk,nu,k,j);
            after = rib_grid_local_legal(grid,nk,nu,k,j,centre);
            if (after <= before) continue;

            if (nrepair == capacity) {
                size_t next = capacity ? capacity * 2 : 64;
                RibGridLocalRepair *grown;
                if (next < capacity ||
                    next > SIZE_MAX / sizeof(*repair)) {
                    free(repair);
                    return 0;
                }
                grown = (RibGridLocalRepair *)realloc(
                    repair,next*sizeof(*repair));
                if (!grown) {
                    free(repair);
                    return 0;
                }
                repair = grown;
                capacity = next;
            }
            repair[nrepair].slot = slot;
            for (int d = 0; d < 3; d++)
                repair[nrepair].p[d] = (float)centre[d];
            repair[nrepair].supported =
                grid_supported != NULL && grid_supported[slot] != 0;
            nrepair++;
        }
    }

    /* Proposals were all certified against the same unmodified grid.  Their
     * compact-ring precondition makes adjacent proposals impossible: an
     * adjacent outlier would itself violate the other's ring certificate. */
    for (size_t r = 0; r < nrepair; r++) {
        size_t slot = repair[r].slot;
        memcpy(&grid[slot*3],repair[r].p,3*sizeof(float));
        if (grid_supported) grid_supported[slot] = 0;
        if (repair[r].supported) supported++;
    }
    free(repair);
    if (out_supported) *out_supported = supported;
    return nrepair;
}

typedef struct {
    size_t slot[2];
    float p[6];
    size_t supported;
} RibGridPairRepair;

static int rib_grid_region_legal(const float *grid, size_t nu,
                                 size_t k0, size_t k1,
                                 size_t j0, size_t j1)
{
    int legal = 0;
    for (size_t k = k0; k <= k1; k++)
        for (size_t j = j0; j <= j1; j++)
            legal += rib_grid_cell_legal(
                grid,nu,k,j,SIZE_MAX,NULL);
    return legal;
}

/* Construct a simultaneous harmonic replacement for a horizontal or vertical
 * pair.  Six exterior grid vertices form the boundary of the two-slot patch.
 * Solving the two discrete Laplace equations gives:
 *     4*x0-x1=A, -x0+4*x1=B.
 * The exterior ring and every proposed edge are certified below the same
 * inter-wrap thresholds as the single-slot repair. */
static int rib_grid_pair_candidate(float *grid, size_t nk, size_t nu,
                                   size_t k, size_t j, int vertical,
                                   size_t slot[2], float candidate[6],
                                   int *out_before, int *out_after)
{
    const double gate2 = RIB_WRAP_GATE * RIB_WRAP_GATE;
    const double ring_diameter2 = 6.5 * 6.5;
    size_t exterior[6];
    double A[3] = {0.0,0.0,0.0}, B[3] = {0.0,0.0,0.0};
    double C[2][3], old[2][3];
    size_t old_bad = 0;
    size_t k0,k1,j0,j1;

    if (!vertical) {
        if (k == 0 || k + 1 >= nk || j == 0 || j + 2 >= nu) return 0;
        slot[0]=k*nu+j; slot[1]=k*nu+j+1;
        exterior[0]=k*nu+j-1;
        exterior[1]=(k-1)*nu+j;
        exterior[2]=(k+1)*nu+j;
        exterior[3]=k*nu+j+2;
        exterior[4]=(k-1)*nu+j+1;
        exterior[5]=(k+1)*nu+j+1;
        k0=k-1; k1=k; j0=j-1; j1=j+1;
    } else {
        if (k == 0 || k + 2 >= nk || j == 0 || j + 1 >= nu) return 0;
        slot[0]=k*nu+j; slot[1]=(k+1)*nu+j;
        exterior[0]=(k-1)*nu+j;
        exterior[1]=k*nu+j-1;
        exterior[2]=k*nu+j+1;
        exterior[3]=(k+2)*nu+j;
        exterior[4]=(k+1)*nu+j-1;
        exterior[5]=(k+1)*nu+j+1;
        k0=k-1; k1=k+1; j0=j-1; j1=j;
    }
    if (!isfinite((double)grid[slot[0]*3]) ||
        !isfinite((double)grid[slot[1]*3])) return 0;
    for (int q = 0; q < 6; q++)
        if (!isfinite((double)grid[exterior[q]*3])) return 0;
    for (int a = 0; a < 6; a++)
        for (int b = a + 1; b < 6; b++)
            if (rib_grid_d2(&grid[exterior[a]*3],
                            &grid[exterior[b]*3]) > ring_diameter2)
                return 0;

    for (int d = 0; d < 3; d++) {
        A[d]=(double)grid[exterior[0]*3+(size_t)d]+
             (double)grid[exterior[1]*3+(size_t)d]+
             (double)grid[exterior[2]*3+(size_t)d];
        B[d]=(double)grid[exterior[3]*3+(size_t)d]+
             (double)grid[exterior[4]*3+(size_t)d]+
             (double)grid[exterior[5]*3+(size_t)d];
        C[0][d]=(4.0*A[d]+B[d])/15.0;
        C[1][d]=(A[d]+4.0*B[d])/15.0;
        old[0][d]=(double)grid[slot[0]*3+(size_t)d];
        old[1][d]=(double)grid[slot[1]*3+(size_t)d];
        candidate[d]=(float)C[0][d];
        candidate[3+d]=(float)C[1][d];
    }
    {
        const size_t adjacent[2][4] = {
            { exterior[0],exterior[1],exterior[2],slot[1] },
            { slot[0],exterior[3],exterior[4],exterior[5] }
        };
        for (int q = 0; q < 2; q++) {
            for (int n = 0; n < 4; n++) {
                const double *other = NULL;
                double d0,d1,d2;
                if (adjacent[q][n] == slot[0]) other=C[0];
                else if (adjacent[q][n] == slot[1]) other=C[1];
                if (other) {
                    d0=C[q][0]-other[0];
                    d1=C[q][1]-other[1];
                    d2=C[q][2]-other[2];
                } else {
                    const float *p=&grid[adjacent[q][n]*3];
                    d0=C[q][0]-(double)p[0];
                    d1=C[q][1]-(double)p[1];
                    d2=C[q][2]-(double)p[2];
                }
                if (d0*d0+d1*d1+d2*d2 > gate2) return 0;
            }
        }
    }
    old_bad += rib_grid_incident_long(grid,nk,nu,k,j,old[0]);
    old_bad += rib_grid_incident_long(
        grid,nk,nu,k+(vertical?1u:0u),j+(vertical?0u:1u),old[1]);
    if (old_bad < 4) return 0;

    *out_before = rib_grid_region_legal(grid,nu,k0,k1,j0,j1);
    memcpy(&grid[slot[0]*3],candidate,3*sizeof(float));
    memcpy(&grid[slot[1]*3],candidate+3,3*sizeof(float));
    *out_after = rib_grid_region_legal(grid,nu,k0,k1,j0,j1);
    for (int d = 0; d < 3; d++) {
        grid[slot[0]*3+(size_t)d]=(float)old[0][d];
        grid[slot[1]*3+(size_t)d]=(float)old[1][d];
    }
    return *out_after > *out_before;
}

static size_t rib_grid_repair_outlier_pairs(
    float *grid, uint8_t *grid_supported, size_t nk, size_t nu,
    size_t *out_slots, size_t *out_supported)
{
    RibGridPairRepair *repair=NULL;
    uint8_t *reserved=NULL;
    size_t nrepair=0,capacity=0,nslot,supported=0;

    if (out_slots) *out_slots=0;
    if (out_supported) *out_supported=0;
    if (!grid || nk<3 || nu<3 || nk>SIZE_MAX/nu) return 0;
    nslot=nk*nu;
    reserved=(uint8_t *)calloc(nslot,1);
    if (!reserved) return 0;

    for (int vertical=0; vertical<=1; vertical++) {
        size_t kend=vertical ? nk-2 : nk-1;
        size_t jend=vertical ? nu-1 : nu-2;
        for (size_t k=1; k<kend; k++) {
            for (size_t j=1; j<jend; j++) {
                size_t slot[2];
                float candidate[6];
                int before=0,after=0;
                if (!rib_grid_pair_candidate(
                        grid,nk,nu,k,j,vertical,slot,candidate,
                        &before,&after)) continue;
                if (reserved[slot[0]] || reserved[slot[1]]) continue;
                if (nrepair==capacity) {
                    size_t next=capacity ? capacity*2 : 32;
                    RibGridPairRepair *grown;
                    if (next<capacity || next>SIZE_MAX/sizeof(*repair)) {
                        free(repair); free(reserved); return 0;
                    }
                    grown=(RibGridPairRepair *)realloc(
                        repair,next*sizeof(*repair));
                    if (!grown) {
                        free(repair); free(reserved); return 0;
                    }
                    repair=grown; capacity=next;
                }
                repair[nrepair].slot[0]=slot[0];
                repair[nrepair].slot[1]=slot[1];
                memcpy(repair[nrepair].p,candidate,6*sizeof(float));
                repair[nrepair].supported=
                    (grid_supported && grid_supported[slot[0]] ? 1u : 0u)+
                    (grid_supported && grid_supported[slot[1]] ? 1u : 0u);
                reserved[slot[0]]=reserved[slot[1]]=1;
                nrepair++;
            }
        }
    }
    for (size_t r=0; r<nrepair; r++) {
        memcpy(&grid[repair[r].slot[0]*3],repair[r].p,3*sizeof(float));
        memcpy(&grid[repair[r].slot[1]*3],repair[r].p+3,3*sizeof(float));
        if (grid_supported) {
            grid_supported[repair[r].slot[0]]=0;
            grid_supported[repair[r].slot[1]]=0;
        }
        supported+=repair[r].supported;
    }
    free(repair); free(reserved);
    if (out_slots) *out_slots=nrepair*2;
    if (out_supported) *out_supported=supported;
    return nrepair;
}

#define RIB_GRID_LOCAL_RUN_MAX 32

typedef struct {
    size_t n;
    size_t slot[RIB_GRID_LOCAL_RUN_MAX];
    float p[RIB_GRID_LOCAL_RUN_MAX*3];
    size_t supported;
} RibGridRunRepair;

static int rib_grid_edge_slots_safe(const float *grid, size_t a, size_t b)
{
    return rib_grid_d2(&grid[a*3],&grid[b*3]) <=
           RIB_WRAP_GATE*RIB_WRAP_GATE;
}

/* Solve a one-cell-wide horizontal/vertical outlier run with its exterior
 * lattice ring held fixed.  This is the N-variable form of the two-slot
 * harmonic repair above.  It deliberately does not grow into a two-dimensional
 * patch: every exterior ring vertex must already be classified as good. */
static int rib_grid_run_candidate(
    float *grid, const uint8_t *bad, size_t nk, size_t nu,
    size_t k, size_t j, size_t n, int vertical,
    size_t slot[RIB_GRID_LOCAL_RUN_MAX],
    float candidate[RIB_GRID_LOCAL_RUN_MAX*3])
{
    double aa[RIB_GRID_LOCAL_RUN_MAX],bb[RIB_GRID_LOCAL_RUN_MAX];
    double cc[RIB_GRID_LOCAL_RUN_MAX],dd[RIB_GRID_LOCAL_RUN_MAX];
    double old[RIB_GRID_LOCAL_RUN_MAX*3];
    size_t side0[RIB_GRID_LOCAL_RUN_MAX],side1[RIB_GRID_LOCAL_RUN_MAX];
    size_t cap0,cap1,k0,k1,j0,j1;
    size_t old_bad=0;
    int before,after;

    if (n<2 || n>RIB_GRID_LOCAL_RUN_MAX) return 0;
    if (!vertical) {
        if (k==0 || k+1>=nk || j==0 || j+n>=nu) return 0;
        cap0=k*nu+j-1; cap1=k*nu+j+n;
        for (size_t q=0;q<n;q++) {
            slot[q]=k*nu+j+q;
            side0[q]=(k-1)*nu+j+q;
            side1[q]=(k+1)*nu+j+q;
        }
        k0=k-1; k1=k; j0=j-1; j1=j+n-1;
    } else {
        if (k==0 || k+n>=nk || j==0 || j+1>=nu) return 0;
        cap0=(k-1)*nu+j; cap1=(k+n)*nu+j;
        for (size_t q=0;q<n;q++) {
            slot[q]=(k+q)*nu+j;
            side0[q]=(k+q)*nu+j-1;
            side1[q]=(k+q)*nu+j+1;
        }
        k0=k-1; k1=k+n-1; j0=j-1; j1=j;
    }
    if (!isfinite((double)grid[cap0*3]) ||
        !isfinite((double)grid[cap1*3]) ||
        bad[cap0] || bad[cap1]) return 0;
    for (size_t q=0;q<n;q++) {
        if (!isfinite((double)grid[slot[q]*3]) ||
            !isfinite((double)grid[side0[q]*3]) ||
            !isfinite((double)grid[side1[q]*3]) ||
            bad[side0[q]] || bad[side1[q]]) return 0;
        if (q+1<n &&
            (!rib_grid_edge_slots_safe(grid,side0[q],side0[q+1]) ||
             !rib_grid_edge_slots_safe(grid,side1[q],side1[q+1])))
            return 0;
    }
    if (!rib_grid_edge_slots_safe(grid,cap0,side0[0]) ||
        !rib_grid_edge_slots_safe(grid,cap0,side1[0]) ||
        !rib_grid_edge_slots_safe(grid,cap1,side0[n-1]) ||
        !rib_grid_edge_slots_safe(grid,cap1,side1[n-1]))
        return 0;

    for (int axis=0;axis<3;axis++) {
        for (size_t q=0;q<n;q++) {
            aa[q]=q>0 ? -1.0 : 0.0;
            bb[q]=4.0;
            cc[q]=q+1<n ? -1.0 : 0.0;
            dd[q]=(double)grid[side0[q]*3+(size_t)axis]+
                  (double)grid[side1[q]*3+(size_t)axis];
            if (q==0)
                dd[q]+=(double)grid[cap0*3+(size_t)axis];
            if (q+1==n)
                dd[q]+=(double)grid[cap1*3+(size_t)axis];
        }
        thomas(aa,bb,cc,dd,(int)n);
        for (size_t q=0;q<n;q++)
            candidate[q*3+(size_t)axis]=(float)dd[q];
    }

    for (size_t q=0;q<n;q++) {
        size_t neighbour[4]={
            side0[q],side1[q],
            q>0 ? slot[q-1] : cap0,
            q+1<n ? slot[q+1] : cap1
        };
        const float *p=&candidate[q*3];
        for (int z=0;z<4;z++) {
            const float *r;
            if (neighbour[z]==slot[q>0?q-1:0] && q>0)
                r=&candidate[(q-1)*3];
            else if (q+1<n && neighbour[z]==slot[q+1])
                r=&candidate[(q+1)*3];
            else
                r=&grid[neighbour[z]*3];
            if (rib_grid_d2(p,r)>RIB_WRAP_GATE*RIB_WRAP_GATE)
                return 0;
        }
        for (int axis=0;axis<3;axis++)
            old[q*3+(size_t)axis]=(double)grid[slot[q]*3+(size_t)axis];
        {
            double op[3]={old[q*3],old[q*3+1],old[q*3+2]};
            size_t qk=k+(vertical?q:0u), qj=j+(vertical?0u:q);
            old_bad+=rib_grid_incident_long(grid,nk,nu,qk,qj,op);
        }
    }
    if (old_bad<2*n) return 0;
    before=rib_grid_region_legal(grid,nu,k0,k1,j0,j1);
    for (size_t q=0;q<n;q++)
        memcpy(&grid[slot[q]*3],&candidate[q*3],3*sizeof(float));
    after=rib_grid_region_legal(grid,nu,k0,k1,j0,j1);
    for (size_t q=0;q<n;q++)
        for (int axis=0;axis<3;axis++)
            grid[slot[q]*3+(size_t)axis]=(float)old[q*3+(size_t)axis];
    return after>before;
}

static size_t rib_grid_repair_outlier_runs(
    float *grid, uint8_t *grid_supported, size_t nk, size_t nu,
    size_t *out_slots, size_t *out_supported)
{
    RibGridRunRepair *repair=NULL;
    uint8_t *bad=NULL,*reserved=NULL;
    size_t nslot,nrepair=0,capacity=0,total_slots=0,supported=0;

    if (out_slots) *out_slots=0;
    if (out_supported) *out_supported=0;
    if (!grid || nk<3 || nu<3 || nk>SIZE_MAX/nu) return 0;
    nslot=nk*nu;
    bad=(uint8_t *)calloc(nslot,1);
    reserved=(uint8_t *)calloc(nslot,1);
    if (!bad || !reserved) { free(bad); free(reserved); return 0; }
    for (size_t k=0;k<nk;k++) {
        for (size_t j=0;j<nu;j++) {
            size_t s=k*nu+j;
            double p[3];
            if (!isfinite((double)grid[s*3])) continue;
            p[0]=(double)grid[s*3];
            p[1]=(double)grid[s*3+1];
            p[2]=(double)grid[s*3+2];
            if (rib_grid_incident_long(grid,nk,nu,k,j,p)>=2)
                bad[s]=1;
        }
    }

    for (int vertical=0;vertical<=1;vertical++) {
        size_t outer=vertical ? nu : nk;
        size_t inner=vertical ? nk : nu;
        for (size_t a=1;a+1<outer;a++) {
            size_t q=1;
            while (q+1<inner) {
                size_t k=vertical?q:a, j=vertical?a:q;
                size_t s=k*nu+j;
                if (!bad[s]) { q++; continue; }
                size_t begin=q;
                while (q+1<inner) {
                    k=vertical?q:a; j=vertical?a:q;
                    if (!bad[k*nu+j]) break;
                    q++;
                }
                size_t n=q-begin;
                if (n<2 || n>RIB_GRID_LOCAL_RUN_MAX) continue;
                k=vertical?begin:a; j=vertical?a:begin;
                size_t slots[RIB_GRID_LOCAL_RUN_MAX];
                float candidate[RIB_GRID_LOCAL_RUN_MAX*3];
                if (!rib_grid_run_candidate(
                        grid,bad,nk,nu,k,j,n,vertical,slots,candidate))
                    continue;
                int overlap=0;
                for (size_t x=0;x<n;x++)
                    if (reserved[slots[x]]) { overlap=1; break; }
                if (overlap) continue;
                if (nrepair==capacity) {
                    size_t next=capacity?capacity*2:32;
                    RibGridRunRepair *grown;
                    if (next<capacity || next>SIZE_MAX/sizeof(*repair)) {
                        free(repair); free(bad); free(reserved); return 0;
                    }
                    grown=(RibGridRunRepair *)realloc(
                        repair,next*sizeof(*repair));
                    if (!grown) {
                        free(repair); free(bad); free(reserved); return 0;
                    }
                    repair=grown; capacity=next;
                }
                memset(&repair[nrepair],0,sizeof(repair[nrepair]));
                repair[nrepair].n=n;
                memcpy(repair[nrepair].slot,slots,n*sizeof(size_t));
                memcpy(repair[nrepair].p,candidate,n*3*sizeof(float));
                for (size_t x=0;x<n;x++) {
                    reserved[slots[x]]=1;
                    if (grid_supported && grid_supported[slots[x]])
                        repair[nrepair].supported++;
                }
                nrepair++;
            }
        }
    }
    for (size_t r=0;r<nrepair;r++) {
        for (size_t q=0;q<repair[r].n;q++) {
            size_t s=repair[r].slot[q];
            memcpy(&grid[s*3],&repair[r].p[q*3],3*sizeof(float));
            if (grid_supported) grid_supported[s]=0;
        }
        total_slots+=repair[r].n;
        supported+=repair[r].supported;
    }
    free(repair); free(bad); free(reserved);
    if (out_slots) *out_slots=total_slots;
    if (out_supported) *out_supported=supported;
    return nrepair;
}

#define RIB_GRID_LOCAL_PATCH_MAX 256
#define RIB_GRID_LOCAL_PATCH_SURROUND_MIN 64
#define RIB_GRID_LOCAL_PATCH_CLAIM_SNAP 3.0

typedef struct {
    size_t offset;
    size_t n;
} RibGridShortComponent;

typedef struct {
    size_t n;
    size_t slot[RIB_GRID_LOCAL_PATCH_MAX];
    float p[RIB_GRID_LOCAL_PATCH_MAX*3];
    uint8_t source_claim[RIB_GRID_LOCAL_PATCH_MAX];
    size_t supported;
    size_t claim_slots;
} RibGridPatchRepair;

static int rib_grid_patch_local_index(const size_t *slot, size_t n,
                                      size_t query)
{
    for (size_t q=0;q<n;q++)
        if (slot[q]==query) return (int)q;
    return -1;
}

static int rib_grid_patch_nearest_claim(
    const SliceSet *S, double U0, double du,
    const int32_t *choff, const int32_t *chids,
    size_t k, size_t j, const float target[3], float out[3])
{
    double best=RIB_GRID_LOCAL_PATCH_CLAIM_SNAP*
                RIB_GRID_LOCAL_PATCH_CLAIM_SNAP;
    double uq=(double)j*du+U0;
    int found=0;
    if (!S || !choff || !chids || k>=(size_t)S->nplanes) return 0;
    for (int32_t x=choff[k];x<choff[k+1];x++) {
        const Chain *ch=&S->chn[chids[x]];
        int32_t f=ch->first,cn=ch->count;
        int32_t lo,hi,seg;
        double w,p[3],d0,d1,d2,distance;
        if (cn<2 || uq<S->smp[f].u || uq>S->smp[f+cn-1].u) continue;
        lo=1; hi=cn-1;
        while (lo<hi) {
            int32_t mid=lo+(hi-lo)/2;
            if (S->smp[f+mid].u<uq) lo=mid+1;
            else hi=mid;
        }
        seg=lo;
        {
            double ua=S->smp[f+seg-1].u,ub=S->smp[f+seg].u;
            w=ub-ua>1e-12 ? (uq-ua)/(ub-ua) : 0.0;
        }
        if (w<0.0) w=0.0;
        if (w>1.0) w=1.0;
        for (int axis=0;axis<3;axis++)
            p[axis]=(1.0-w)*S->smp[f+seg-1].p[axis]+
                           w *S->smp[f+seg].p[axis];
        d0=p[0]-(double)target[0];
        d1=p[1]-(double)target[1];
        d2=p[2]-(double)target[2];
        distance=d0*d0+d1*d1+d2*d2;
        if (distance>best) continue;
        best=distance;
        for (int axis=0;axis<3;axis++) out[axis]=(float)p[axis];
        found=1;
    }
    /* A chart-registration defect can put the correct physical chain at the
     * wrong solved U, so same-U claims alone are circular evidence.  Project
     * the geometry target onto every measured polyline segment in this slice
     * and retain the nearest point inside the same strict three-voxel snap
     * radius.  The fitted grid supplies the corrected U; the source slice
     * supplies the actual surface position. */
    for (int32_t c=choff[k];c<choff[k+1];c++) {
        const Chain *ch=&S->chn[chids[c]];
        int32_t f=ch->first,cn=ch->count;
        for (int32_t seg=1;seg<cn;seg++) {
            const double *a=S->smp[f+seg-1].p;
            const double *b=S->smp[f+seg].p;
            double v0=b[0]-a[0],v1=b[1]-a[1],v2=b[2]-a[2];
            double w0=(double)target[0]-a[0];
            double w1=(double)target[1]-a[1];
            double w2=(double)target[2]-a[2];
            double vv=v0*v0+v1*v1+v2*v2;
            double t=vv>1e-18 ? (w0*v0+w1*v1+w2*v2)/vv : 0.0;
            double p[3],d0,d1,d2,distance;
            if (t<0.0) t=0.0;
            if (t>1.0) t=1.0;
            p[0]=a[0]+t*v0; p[1]=a[1]+t*v1; p[2]=a[2]+t*v2;
            d0=p[0]-(double)target[0];
            d1=p[1]-(double)target[1];
            d2=p[2]-(double)target[2];
            distance=d0*d0+d1*d1+d2*d2;
            if (distance>best) continue;
            best=distance;
            out[0]=(float)p[0]; out[1]=(float)p[1]; out[2]=(float)p[2];
            found=1;
        }
    }
    return found;
}

static int rib_grid_patch_validate(
    float *grid, const int32_t *label, int32_t cid, size_t nu,
    const size_t *slot, size_t n, const float *candidate,
    size_t kmin, size_t kmax, size_t jmin, size_t jmax,
    int *out_before, int *out_after, size_t *out_touched)
{
    double old[RIB_GRID_LOCAL_PATCH_MAX*3];
    int before=0,after=0;
    size_t touched=0;
    for (size_t ck=kmin-1;ck<=kmax;ck++) {
        for (size_t cj=jmin-1;cj<=jmax;cj++) {
            size_t corner[4]={ck*nu+cj,ck*nu+cj+1,
                              (ck+1)*nu+cj,(ck+1)*nu+cj+1};
            int has_patch=0;
            for (int z=0;z<4;z++)
                if (label[corner[z]]==cid) { has_patch=1; break; }
            if (!has_patch) continue;
            touched++;
            before+=rib_grid_cell_legal(
                grid,nu,ck,cj,SIZE_MAX,NULL);
        }
    }
    if (touched==0) return 0;
    for (size_t q=0;q<n;q++) {
        for (int axis=0;axis<3;axis++)
            old[q*3+(size_t)axis]=(double)grid[slot[q]*3+(size_t)axis];
        memcpy(&grid[slot[q]*3],&candidate[q*3],3*sizeof(float));
    }
    for (size_t ck=kmin-1;ck<=kmax;ck++) {
        for (size_t cj=jmin-1;cj<=jmax;cj++) {
            size_t corner[4]={ck*nu+cj,ck*nu+cj+1,
                              (ck+1)*nu+cj,(ck+1)*nu+cj+1};
            int has_patch=0;
            for (int z=0;z<4;z++)
                if (label[corner[z]]==cid) { has_patch=1; break; }
            if (has_patch)
                after+=rib_grid_cell_legal(
                    grid,nu,ck,cj,SIZE_MAX,NULL);
        }
    }
    for (size_t q=0;q<n;q++)
        for (int axis=0;axis<3;axis++)
            grid[slot[q]*3+(size_t)axis]=
                (float)old[q*3+(size_t)axis];
    if (out_before) *out_before=before;
    if (out_after) *out_after=after;
    if (out_touched) *out_touched=touched;
    return after>before && (size_t)after==2*touched;
}

#define RIB_GRID_LOCAL_PATCH_REFINE_MAX 1.0

static double rib_grid_patch_work_distance(
    const float *grid, const size_t *slot, size_t n,
    const double *work, size_t a, size_t b)
{
    int ia=rib_grid_patch_local_index(slot,n,a);
    int ib=rib_grid_patch_local_index(slot,n,b);
    const double *pa=ia>=0 ? &work[(size_t)ia*3] : NULL;
    const double *pb=ib>=0 ? &work[(size_t)ib*3] : NULL;
    double a0=pa?pa[0]:(double)grid[a*3];
    double a1=pa?pa[1]:(double)grid[a*3+1];
    double a2=pa?pa[2]:(double)grid[a*3+2];
    double b0=pb?pb[0]:(double)grid[b*3];
    double b1=pb?pb[1]:(double)grid[b*3+1];
    double b2=pb?pb[2]:(double)grid[b*3+2];
    double d0=a0-b0,d1=a1-b1,d2=a2-b2;
    return sqrt(d0*d0+d1*d1+d2*d2);
}

static int rib_grid_patch_project_edge(
    const float *grid, const size_t *slot, size_t n,
    double *work, size_t a, size_t b, double limit)
{
    int ia=rib_grid_patch_local_index(slot,n,a);
    int ib=rib_grid_patch_local_index(slot,n,b);
    double pa[3],pb[3],d[3],length2,length,scale;
    for (int axis=0;axis<3;axis++) {
        pa[axis]=ia>=0 ? work[(size_t)ia*3+(size_t)axis]
                       : (double)grid[a*3+(size_t)axis];
        pb[axis]=ib>=0 ? work[(size_t)ib*3+(size_t)axis]
                       : (double)grid[b*3+(size_t)axis];
        d[axis]=pa[axis]-pb[axis];
    }
    length2=d[0]*d[0]+d[1]*d[1]+d[2]*d[2];
    if (length2<=limit*limit) return 1;
    if (ia<0 && ib<0) return 0;
    length=sqrt(length2);
    scale=(length-limit)/length;
    if (ia>=0 && ib>=0) scale*=0.5;
    if (ia>=0)
        for (int axis=0;axis<3;axis++)
            work[(size_t)ia*3+(size_t)axis]-=scale*d[axis];
    if (ib>=0)
        for (int axis=0;axis<3;axis++)
            work[(size_t)ib*3+(size_t)axis]+=scale*d[axis];
    return 1;
}

/* A measured or harmonic target can miss the hard edge gate by a fraction of
 * a voxel even when a nearby coherent patch exists.  Project only the patch
 * variables onto the convex axial-edge balls and one selected diagonal ball
 * per touched cell, while keeping every vertex within one voxel of its target.
 * This is a bounded geometric refinement, not permission to bridge an
 * incompatible exterior: a long fixed-fixed edge fails immediately. */
static int rib_grid_patch_refine_full_cells(
    float *grid, const int32_t *label, int32_t cid, size_t nu,
    const size_t *slot, size_t n, const float *target, float *candidate,
    size_t kmin, size_t kmax, size_t jmin, size_t jmax,
    int *out_before, int *out_after, size_t *out_touched)
{
    const double limit=RIB_WRAP_GATE-0.02;
    double work[RIB_GRID_LOCAL_PATCH_MAX*3];
    for (size_t q=0;q<n;q++)
        for (int axis=0;axis<3;axis++)
            work[q*3+(size_t)axis]=(double)target[q*3+(size_t)axis];

    for (int sweep=0;sweep<300;sweep++) {
        int feasible=1;
        if (sweep<240) {
            for (size_t q=0;q<n;q++)
                for (int axis=0;axis<3;axis++)
                    work[q*3+(size_t)axis]=
                        0.98*work[q*3+(size_t)axis]+
                        0.02*(double)target[q*3+(size_t)axis];
        }
        for (size_t ck=kmin-1;ck<=kmax && feasible;ck++) {
            for (size_t cj=jmin-1;cj<=jmax && feasible;cj++) {
                size_t a=ck*nu+cj,b=a+1,c=a+nu,d=c+1;
                size_t corner[4]={a,b,c,d};
                int has_patch=0;
                for (int z=0;z<4;z++)
                    if (label[corner[z]]==cid) { has_patch=1; break; }
                if (!has_patch) continue;
                feasible=
                    rib_grid_patch_project_edge(
                        grid,slot,n,work,a,b,limit) &&
                    rib_grid_patch_project_edge(
                        grid,slot,n,work,b,d,limit) &&
                    rib_grid_patch_project_edge(
                        grid,slot,n,work,d,c,limit) &&
                    rib_grid_patch_project_edge(
                        grid,slot,n,work,c,a,limit);
                if (!feasible) break;
                {
                    int ad_variable=
                        rib_grid_patch_local_index(slot,n,a)>=0 ||
                        rib_grid_patch_local_index(slot,n,d)>=0;
                    int bc_variable=
                        rib_grid_patch_local_index(slot,n,b)>=0 ||
                        rib_grid_patch_local_index(slot,n,c)>=0;
                    double ad=rib_grid_patch_work_distance(
                        grid,slot,n,work,a,d);
                    double bc=rib_grid_patch_work_distance(
                        grid,slot,n,work,b,c);
                    double ad_score=(!ad_variable && ad>limit)?1e300:ad;
                    double bc_score=(!bc_variable && bc>limit)?1e300:bc;
                    if (ad_score>=1e299 && bc_score>=1e299) {
                        feasible=0;
                        break;
                    }
                    if (ad_score<=bc_score)
                        feasible=rib_grid_patch_project_edge(
                            grid,slot,n,work,a,d,limit);
                    else
                        feasible=rib_grid_patch_project_edge(
                            grid,slot,n,work,b,c,limit);
                }
            }
        }
        if (!feasible) return 0;
        for (size_t q=0;q<n;q++) {
            double delta[3],distance2=0.0;
            for (int axis=0;axis<3;axis++) {
                delta[axis]=work[q*3+(size_t)axis]-
                            (double)target[q*3+(size_t)axis];
                distance2+=delta[axis]*delta[axis];
            }
            if (distance2>RIB_GRID_LOCAL_PATCH_REFINE_MAX*
                         RIB_GRID_LOCAL_PATCH_REFINE_MAX) {
                double scale=RIB_GRID_LOCAL_PATCH_REFINE_MAX/sqrt(distance2);
                for (int axis=0;axis<3;axis++)
                    work[q*3+(size_t)axis]=
                        (double)target[q*3+(size_t)axis]+scale*delta[axis];
            }
        }
    }
    for (size_t q=0;q<n;q++)
        for (int axis=0;axis<3;axis++)
            candidate[q*3+(size_t)axis]=
                (float)work[q*3+(size_t)axis];
    return rib_grid_patch_validate(
        grid,label,cid,nu,slot,n,candidate,kmin,kmax,jmin,jmax,
        out_before,out_after,out_touched);
}

static int rib_grid_patch_refine_neighbourhood(
    float *grid, size_t nk, size_t nu,
    const size_t *core, size_t ncore,
    const size_t *slot, size_t n, const float *target,
    const double *move_limit, float *candidate,
    int *out_before, int *out_after)
{
    const double edge_limit=RIB_WRAP_GATE-0.02;
    double work[RIB_GRID_LOCAL_PATCH_MAX*3];
    float original_pos[RIB_GRID_LOCAL_PATCH_MAX*3];
    size_t kmin=nk,kmax=0,jmin=nu,jmax=0;
    int before=0,after=0;

    if (n==0 || n>RIB_GRID_LOCAL_PATCH_MAX) return 0;
    for (size_t q=0;q<n;q++) {
        size_t k=slot[q]/nu,j=slot[q]-k*nu;
        if (k==0 || k+1>=nk || j==0 || j+1>=nu) return 0;
        if (k<kmin) kmin=k;
        if (k>kmax) kmax=k;
        if (j<jmin) jmin=j;
        if (j>jmax) jmax=j;
        for (int axis=0;axis<3;axis++) {
            if (!isfinite((double)target[q*3+(size_t)axis])) return 0;
            work[q*3+(size_t)axis]=(double)target[q*3+(size_t)axis];
            original_pos[q*3+(size_t)axis]=
                grid[slot[q]*3+(size_t)axis];
        }
    }

    for (int sweep=0;sweep<360;sweep++) {
        int feasible=1;
        if (sweep<300) {
            for (size_t q=0;q<n;q++)
                for (int axis=0;axis<3;axis++)
                    work[q*3+(size_t)axis]=
                        0.985*work[q*3+(size_t)axis]+
                        0.015*(double)target[q*3+(size_t)axis];
        }
        /* Preserve every existing axial connection incident to the movable
         * neighbourhood, including the ring-to-fixed-sheet interface. */
        for (size_t q=0;q<n && feasible;q++) {
            size_t s=slot[q],k=s/nu,j=s-k*nu;
            size_t neighbour[4],nn=0;
            if (j>0) neighbour[nn++]=s-1;
            if (j+1<nu) neighbour[nn++]=s+1;
            if (k>0) neighbour[nn++]=s-nu;
            if (k+1<nk) neighbour[nn++]=s+nu;
            for (size_t z=0;z<nn;z++) {
                size_t t=neighbour[z];
                if (!isfinite((double)grid[t*3])) continue;
                if (!rib_grid_patch_project_edge(
                        grid,slot,n,work,s,t,edge_limit)) {
                    feasible=0;
                    break;
                }
            }
        }
        /* Only cells touching the original bad core are required to become
         * complete.  Outer ring cells are checked for non-regression below. */
        for (size_t ck=kmin-1;ck<=kmax && feasible;ck++) {
            for (size_t cj=jmin-1;cj<=jmax && feasible;cj++) {
                size_t a=ck*nu+cj,b=a+1,c=a+nu,d=c+1;
                size_t corner[4]={a,b,c,d};
                int touches_core=0;
                for (int z=0;z<4;z++)
                    if (rib_grid_patch_local_index(
                            core,ncore,corner[z])>=0) {
                        touches_core=1;
                        break;
                    }
                if (!touches_core) continue;
                for (int z=0;z<4;z++)
                    if (!isfinite((double)grid[corner[z]*3]))
                        feasible=0;
                if (!feasible) break;
                feasible=
                    rib_grid_patch_project_edge(
                        grid,slot,n,work,a,b,edge_limit) &&
                    rib_grid_patch_project_edge(
                        grid,slot,n,work,b,d,edge_limit) &&
                    rib_grid_patch_project_edge(
                        grid,slot,n,work,d,c,edge_limit) &&
                    rib_grid_patch_project_edge(
                        grid,slot,n,work,c,a,edge_limit);
                if (!feasible) break;
                {
                    int ad_variable=
                        rib_grid_patch_local_index(slot,n,a)>=0 ||
                        rib_grid_patch_local_index(slot,n,d)>=0;
                    int bc_variable=
                        rib_grid_patch_local_index(slot,n,b)>=0 ||
                        rib_grid_patch_local_index(slot,n,c)>=0;
                    double ad=rib_grid_patch_work_distance(
                        grid,slot,n,work,a,d);
                    double bc=rib_grid_patch_work_distance(
                        grid,slot,n,work,b,c);
                    double ad_score=(!ad_variable && ad>edge_limit)?1e300:ad;
                    double bc_score=(!bc_variable && bc>edge_limit)?1e300:bc;
                    if (ad_score>=1e299 && bc_score>=1e299) {
                        feasible=0;
                        break;
                    }
                    if (ad_score<=bc_score)
                        feasible=rib_grid_patch_project_edge(
                            grid,slot,n,work,a,d,edge_limit);
                    else
                        feasible=rib_grid_patch_project_edge(
                            grid,slot,n,work,b,c,edge_limit);
                }
            }
        }
        if (!feasible) return 0;
        for (size_t q=0;q<n;q++) {
            double delta[3],distance2=0.0;
            for (int axis=0;axis<3;axis++) {
                delta[axis]=work[q*3+(size_t)axis]-
                            (double)target[q*3+(size_t)axis];
                distance2+=delta[axis]*delta[axis];
            }
            if (distance2>move_limit[q]*move_limit[q]) {
                double scale=move_limit[q]/sqrt(distance2);
                for (int axis=0;axis<3;axis++)
                    work[q*3+(size_t)axis]=
                        (double)target[q*3+(size_t)axis]+scale*delta[axis];
            }
        }
    }

    for (size_t ck=kmin-1;ck<=kmax;ck++) {
        for (size_t cj=jmin-1;cj<=jmax;cj++) {
            size_t corner[4]={ck*nu+cj,ck*nu+cj+1,
                              (ck+1)*nu+cj,(ck+1)*nu+cj+1};
            int touches_variable=0;
            for (int z=0;z<4;z++)
                if (rib_grid_patch_local_index(slot,n,corner[z])>=0) {
                    touches_variable=1;
                    break;
                }
            if (touches_variable)
                before+=rib_grid_cell_legal(
                    grid,nu,ck,cj,SIZE_MAX,NULL);
        }
    }
    for (size_t q=0;q<n;q++) {
        for (int axis=0;axis<3;axis++)
            candidate[q*3+(size_t)axis]=
                (float)work[q*3+(size_t)axis];
        memcpy(&grid[slot[q]*3],&candidate[q*3],3*sizeof(float));
    }
    {
        int valid=1;
        for (size_t ck=kmin-1;ck<=kmax;ck++) {
            for (size_t cj=jmin-1;cj<=jmax;cj++) {
                size_t corner[4]={ck*nu+cj,ck*nu+cj+1,
                                  (ck+1)*nu+cj,(ck+1)*nu+cj+1};
                int touches_variable=0,touches_core=0;
                for (int z=0;z<4;z++) {
                    if (rib_grid_patch_local_index(
                            slot,n,corner[z])>=0) touches_variable=1;
                    if (rib_grid_patch_local_index(
                            core,ncore,corner[z])>=0) touches_core=1;
                }
                if (!touches_variable) continue;
                {
                    int current=rib_grid_cell_legal(
                        grid,nu,ck,cj,SIZE_MAX,NULL);
                    int original;
                    for (size_t q=0;q<n;q++)
                        memcpy(&grid[slot[q]*3],&original_pos[q*3],
                               3*sizeof(float));
                    original=rib_grid_cell_legal(
                        grid,nu,ck,cj,SIZE_MAX,NULL);
                    for (size_t q=0;q<n;q++)
                        memcpy(&grid[slot[q]*3],&candidate[q*3],
                               3*sizeof(float));
                    if ((touches_core && current!=2) ||
                        (!touches_core && current<original))
                        valid=0;
                    after+=current;
                }
            }
        }
        for (size_t q=0;q<n;q++)
            memcpy(&grid[slot[q]*3],&original_pos[q*3],3*sizeof(float));
        if (out_before) *out_before=before;
        if (out_after) *out_after=after;
        return valid && after>before;
    }
}

/* Replace one small physical island in the fitted grid only when the lattice
 * itself supplies an unambiguous certificate.  Short axial edges define
 * physical grid components.  A candidate component must be completely
 * enclosed in finite samples, and every exterior axial neighbour must belong
 * to one much larger short-edge component.  That larger component is the
 * trusted surrounding sheet; this rejects overlaps between independent sheets
 * and bays at a real chart boundary.
 *
 * The candidate positions solve the discrete Dirichlet problem on the island.
 * We then require every grid cell touched by the island to admit both output
 * triangles under the normal six-voxel gate.  Thus an accepted replacement
 * removes a bad fitted-grid selection at its source rather than capping the
 * emitted mesh or changing its vertex topology. */
static int rib_grid_patch_candidate(
    float *grid, const int32_t *label,
    const RibGridShortComponent *component, const int32_t *order,
    const SliceSet *S, double U0, double du,
    const int32_t *choff, const int32_t *chids,
    int32_t cid, size_t nk, size_t nu,
    size_t slot[RIB_GRID_LOCAL_PATCH_MAX],
    float candidate[RIB_GRID_LOCAL_PATCH_MAX*3],
    uint8_t source_claim[RIB_GRID_LOCAL_PATCH_MAX],
    size_t *out_n, size_t *out_claim_slots, int *out_enclosed)
{
    const RibGridShortComponent *C=&component[cid];
    double *A=NULL;
    double rhs[RIB_GRID_LOCAL_PATCH_MAX*3];
    double y[RIB_GRID_LOCAL_PATCH_MAX],x[RIB_GRID_LOCAL_PATCH_MAX];
    float harmonic[RIB_GRID_LOCAL_PATCH_MAX*3];
    float snapped[RIB_GRID_LOCAL_PATCH_MAX*3];
    uint8_t snapped_claim[RIB_GRID_LOCAL_PATCH_MAX];
    int32_t surround=-1;
    size_t kmin=nk,kmax=0,jmin=nu,jmax=0;
    size_t n=C->n;
    int before=0,after=0;
    size_t touched=0;
    size_t claim_slots=0;
    int accepted=0,mode=0;

    if (out_enclosed) *out_enclosed=0;
    if (out_n) *out_n=n;
    if (out_claim_slots) *out_claim_slots=0;
    memset(source_claim,0,n*sizeof(*source_claim));
    if (n==0 || n>RIB_GRID_LOCAL_PATCH_MAX) return 0;
    for (size_t q=0;q<n;q++) {
        size_t s=(size_t)order[C->offset+q];
        size_t k=s/nu,j=s-k*nu;
        size_t neighbour[4];
        if (k==0 || k+1>=nk || j==0 || j+1>=nu) return 0;
        slot[q]=s;
        if (k<kmin) kmin=k;
        if (k>kmax) kmax=k;
        if (j<jmin) jmin=j;
        if (j>jmax) jmax=j;
        neighbour[0]=s-1; neighbour[1]=s+1;
        neighbour[2]=s-nu; neighbour[3]=s+nu;
        for (int z=0;z<4;z++) {
            size_t t=neighbour[z];
            int32_t tc;
            if (!isfinite((double)grid[t*3])) return 0;
            tc=label[t];
            if (tc==cid) continue;
            if (tc<0) return 0;
            if (surround<0) surround=tc;
            else if (surround!=tc) return 0;
        }
    }
    if (surround<0 ||
        component[surround].n<=RIB_GRID_LOCAL_PATCH_SURROUND_MIN ||
        component[surround].n<4*n+8) return 0;
    if (out_enclosed) *out_enclosed=1;

    if (n>SIZE_MAX/n || n*n>SIZE_MAX/sizeof(*A)) goto done;
    A=(double *)calloc(n*n,sizeof(*A));
    if (!A) goto done;
    memset(rhs,0,n*3*sizeof(*rhs));
    for (size_t q=0;q<n;q++) {
        size_t s=slot[q];
        size_t neighbour[4]={s-1,s+1,s-nu,s+nu};
        A[q*n+q]=4.0;
        for (int z=0;z<4;z++) {
            int r=rib_grid_patch_local_index(slot,n,neighbour[z]);
            if (r>=0) {
                A[q*n+(size_t)r]-=1.0;
            } else {
                const float *p=&grid[neighbour[z]*3];
                for (int axis=0;axis<3;axis++)
                    rhs[q*3+(size_t)axis]+=(double)p[axis];
            }
        }
    }

    /* Dense Cholesky is bounded (at most 256 unknowns), deterministic, and avoids
     * making the global sparse parameterization solver part of this local
     * geometric certificate. */
    for (size_t i=0;i<n;i++) {
        for (size_t j=0;j<=i;j++) {
            double sum=A[i*n+j];
            for (size_t z=0;z<j;z++) sum-=A[i*n+z]*A[j*n+z];
            if (i==j) {
                if (!(sum>1e-12) || !isfinite(sum)) goto done;
                A[i*n+j]=sqrt(sum);
            } else {
                A[i*n+j]=sum/A[j*n+j];
            }
        }
    }
    for (int axis=0;axis<3;axis++) {
        for (size_t i=0;i<n;i++) {
            double sum=rhs[i*3+(size_t)axis];
            for (size_t z=0;z<i;z++) sum-=A[i*n+z]*y[z];
            y[i]=sum/A[i*n+i];
        }
        for (size_t ii=n;ii-- > 0;) {
            double sum=y[ii];
            for (size_t z=ii+1;z<n;z++) sum-=A[z*n+ii]*x[z];
            x[ii]=sum/A[ii*n+ii];
        }
        for (size_t q=0;q<n;q++) {
            if (!isfinite(x[q])) goto done;
            harmonic[q*3+(size_t)axis]=(float)x[q];
        }
    }

    memcpy(snapped,harmonic,n*3*sizeof(*snapped));
    memset(snapped_claim,0,n*sizeof(*snapped_claim));
    for (size_t q=0;q<n;q++) {
        size_t k=slot[q]/nu,j=slot[q]-k*nu;
        if (rib_grid_patch_nearest_claim(
                S,U0,du,choff,chids,k,j,&harmonic[q*3],
                &snapped[q*3])) {
            snapped_claim[q]=1;
            claim_slots++;
        }
    }
    /* Large patches are accepted only when measured slice claims support at
     * least three quarters of their interior.  Small patches retain the prior
     * harmonic fallback, but measured claims are always tried first. */
    if (claim_slots>0 &&
        (n<=64 || claim_slots*4>=n*3) &&
        rib_grid_patch_validate(
            grid,label,cid,nu,slot,n,snapped,kmin,kmax,jmin,jmax,
            &before,&after,&touched)) {
        memcpy(candidate,snapped,n*3*sizeof(*candidate));
        memcpy(source_claim,snapped_claim,n*sizeof(*source_claim));
        accepted=1;
        mode=1;
        goto done;
    }
    if (claim_slots>0 && n<=64 &&
        rib_grid_patch_refine_full_cells(
            grid,label,cid,nu,slot,n,snapped,candidate,
            kmin,kmax,jmin,jmax,&before,&after,&touched)) {
        memset(source_claim,0,n*sizeof(*source_claim));
        accepted=1;
        mode=3;
        goto done;
    }
    if (n<=64 && rib_grid_patch_validate(
            grid,label,cid,nu,slot,n,harmonic,kmin,kmax,jmin,jmax,
            &before,&after,&touched)) {
        memcpy(candidate,harmonic,n*3*sizeof(*candidate));
        memset(source_claim,0,n*sizeof(*source_claim));
        accepted=1;
        mode=2;
        goto done;
    }
    if (n<=64 && rib_grid_patch_refine_full_cells(
            grid,label,cid,nu,slot,n,harmonic,candidate,
            kmin,kmax,jmin,jmax,&before,&after,&touched)) {
        memset(source_claim,0,n*sizeof(*source_claim));
        accepted=1;
        mode=4;
    }
    if (!accepted && n<=64) {
        size_t core[RIB_GRID_LOCAL_PATCH_MAX];
        size_t expanded[RIB_GRID_LOCAL_PATCH_MAX];
        float neighbourhood_target[RIB_GRID_LOCAL_PATCH_MAX*3];
        float neighbourhood_candidate[RIB_GRID_LOCAL_PATCH_MAX*3];
        double move_limit[RIB_GRID_LOCAL_PATCH_MAX];
        size_t ne=n;
        int expanded_before=0,expanded_after=0;
        memcpy(core,slot,n*sizeof(*core));
        memcpy(expanded,slot,n*sizeof(*expanded));
        for (size_t q=0;q<n;q++) {
            size_t s=core[q],k=s/nu,j=s-k*nu;
            for (long dk=-1;dk<=1;dk++) {
                for (long dj=-1;dj<=1;dj++) {
                    long kk=(long)k+dk,jj=(long)j+dj;
                    size_t t;
                    if ((dk==0 && dj==0) || kk<0 || jj<0 ||
                        (size_t)kk>=nk || (size_t)jj>=nu) continue;
                    t=(size_t)kk*nu+(size_t)jj;
                    if (label[t]!=surround ||
                        rib_grid_patch_local_index(expanded,ne,t)>=0)
                        continue;
                    if (ne>=RIB_GRID_LOCAL_PATCH_MAX) {
                        ne=0;
                        break;
                    }
                    expanded[ne++]=t;
                }
                if (ne==0) break;
            }
            if (ne==0) break;
        }
        if (ne>n) {
            for (size_t q=0;q<ne;q++) {
                if (q<n) {
                    memcpy(&neighbourhood_target[q*3],&snapped[q*3],
                           3*sizeof(float));
                    move_limit[q]=RIB_GRID_LOCAL_PATCH_REFINE_MAX;
                } else {
                    memcpy(&neighbourhood_target[q*3],
                           &grid[expanded[q]*3],3*sizeof(float));
                    move_limit[q]=0.75;
                }
            }
            if (rib_grid_patch_refine_neighbourhood(
                    grid,nk,nu,core,n,expanded,ne,neighbourhood_target,
                    move_limit,neighbourhood_candidate,
                    &expanded_before,&expanded_after)) {
                memcpy(slot,expanded,ne*sizeof(*slot));
                memcpy(candidate,neighbourhood_candidate,
                       ne*3*sizeof(*candidate));
                memset(source_claim,0,n*sizeof(*source_claim));
                for (size_t q=n;q<ne;q++) source_claim[q]=2;
                if (out_n) *out_n=ne;
                if (out_claim_slots) *out_claim_slots=0;
                before=expanded_before;
                after=expanded_after;
                touched=0;
                accepted=1;
                mode=5;
            }
        }
    }

done:
    fprintf(stderr,
            "  fitted-grid enclosed island: cid=%d slots=%zu "
            "uvbox=[%zu,%zu]x[%zu,%zu] surround=%d/%zu claims=%zu "
            "legal=%d->%d/%zu %s\n",
            cid,n,jmin,jmax,kmin,kmax,surround,
            surround>=0?component[surround].n:0,claim_slots,
            before,after,2*touched,
            !accepted?"rejected":
            mode==1?"accepted-claims":
            mode==2?"accepted-harmonic":
            mode==3?"accepted-refined-claims":
            mode==4?"accepted-refined-harmonic":
                    "accepted-neighbourhood");
    free(A);
    if (accepted && mode==1 && out_claim_slots)
        *out_claim_slots=claim_slots;
    return accepted;
}

static size_t rib_grid_repair_enclosed_outlier_islands(
    float *grid, uint8_t *grid_supported,
    const SliceSet *S, double U0, double du,
    const int32_t *choff, const int32_t *chids,
    size_t nk, size_t nu,
    size_t *out_slots, size_t *out_supported,
    size_t *out_small_components, size_t *out_enclosed_components,
    size_t *out_claim_slots)
{
    RibGridShortComponent *component=NULL;
    RibGridPatchRepair *repair=NULL;
    int32_t *label=NULL,*order=NULL;
    uint8_t *reserved=NULL;
    size_t nslot,ncomponent=0,component_capacity=0,nordered=0;
    size_t nrepair=0,repair_capacity=0,applied_repairs=0;
    size_t total_slots=0,supported=0,claim_slots=0;
    size_t small_components=0,enclosed_components=0;

    if (out_slots) *out_slots=0;
    if (out_supported) *out_supported=0;
    if (out_small_components) *out_small_components=0;
    if (out_enclosed_components) *out_enclosed_components=0;
    if (out_claim_slots) *out_claim_slots=0;
    if (!grid || nk<3 || nu<3 || nk>SIZE_MAX/nu) return 0;
    nslot=nk*nu;
    if (nslot>(size_t)INT32_MAX ||
        nslot>SIZE_MAX/sizeof(*label) ||
        nslot>SIZE_MAX/sizeof(*order)) return 0;
    label=(int32_t *)malloc(nslot*sizeof(*label));
    order=(int32_t *)malloc(nslot*sizeof(*order));
    reserved=(uint8_t *)calloc(nslot,1);
    if (!label || !order || !reserved) goto done;
    memset(label,0xff,nslot*sizeof(*label));

    /* Label components of the grid graph whose axial edges pass the physical
     * emission gate.  `order` doubles as the BFS queue and leaves each
     * component in one contiguous range, avoiding a second large index. */
    for (size_t seed=0;seed<nslot;seed++) {
        size_t begin,head;
        int32_t cid;
        if (!isfinite((double)grid[seed*3]) || label[seed]>=0) continue;
        if (ncomponent==component_capacity) {
            size_t next=component_capacity?component_capacity*2:64;
            RibGridShortComponent *grown;
            if (next<component_capacity ||
                next>SIZE_MAX/sizeof(*component)) goto done;
            grown=(RibGridShortComponent *)realloc(
                component,next*sizeof(*component));
            if (!grown) goto done;
            component=grown; component_capacity=next;
        }
        if (ncomponent>(size_t)INT32_MAX) goto done;
        cid=(int32_t)ncomponent;
        begin=nordered;
        label[seed]=cid;
        order[nordered++]=(int32_t)seed;
        head=begin;
        while (head<nordered) {
            size_t s=(size_t)order[head++];
            size_t k=s/nu,j=s-k*nu;
            size_t neighbour[4],nn=0;
            if (j>0) neighbour[nn++]=s-1;
            if (j+1<nu) neighbour[nn++]=s+1;
            if (k>0) neighbour[nn++]=s-nu;
            if (k+1<nk) neighbour[nn++]=s+nu;
            for (size_t z=0;z<nn;z++) {
                size_t t=neighbour[z];
                if (label[t]>=0 || !isfinite((double)grid[t*3]) ||
                    !rib_grid_edge_slots_safe(grid,s,t)) continue;
                label[t]=cid;
                order[nordered++]=(int32_t)t;
            }
        }
        component[ncomponent].offset=begin;
        component[ncomponent].n=nordered-begin;
        ncomponent++;
    }

    for (size_t c=0;c<ncomponent;c++) {
        size_t slot[RIB_GRID_LOCAL_PATCH_MAX];
        float candidate[RIB_GRID_LOCAL_PATCH_MAX*3];
        uint8_t source_claim[RIB_GRID_LOCAL_PATCH_MAX];
        size_t patch_n=component[c].n;
        size_t patch_claim_slots=0;
        int enclosed=0;
        if (component[c].n==0 ||
            component[c].n>RIB_GRID_LOCAL_PATCH_MAX) continue;
        small_components++;
        if (!rib_grid_patch_candidate(
                grid,label,component,order,S,U0,du,choff,chids,
                (int32_t)c,nk,nu,slot,candidate,source_claim,
                &patch_n,&patch_claim_slots,&enclosed)) {
            if (enclosed) enclosed_components++;
            continue;
        }
        if (enclosed) enclosed_components++;
        {
            int overlap=0;
            for (size_t q=0;q<patch_n;q++)
                if (reserved[slot[q]]) { overlap=1; break; }
            if (overlap) continue;
        }
        if (nrepair==repair_capacity) {
            size_t next=repair_capacity?repair_capacity*2:32;
            RibGridPatchRepair *grown;
            if (next<repair_capacity || next>SIZE_MAX/sizeof(*repair))
                goto done;
            grown=(RibGridPatchRepair *)realloc(
                repair,next*sizeof(*repair));
            if (!grown) goto done;
            repair=grown; repair_capacity=next;
        }
        memset(&repair[nrepair],0,sizeof(repair[nrepair]));
        repair[nrepair].n=patch_n;
        memcpy(repair[nrepair].slot,slot,
               patch_n*sizeof(*slot));
        memcpy(repair[nrepair].p,candidate,
               patch_n*3*sizeof(*candidate));
        memcpy(repair[nrepair].source_claim,source_claim,
               patch_n*sizeof(*source_claim));
        repair[nrepair].claim_slots=patch_claim_slots;
        for (size_t q=0;q<patch_n;q++)
            if (grid_supported && grid_supported[slot[q]])
                repair[nrepair].supported++;
        for (size_t q=0;q<patch_n;q++) reserved[slot[q]]=1;
        nrepair++;
    }

    for (size_t r=0;r<nrepair;r++) {
        for (size_t q=0;q<repair[r].n;q++) {
            size_t s=repair[r].slot[q];
            memcpy(&grid[s*3],&repair[r].p[q*3],3*sizeof(float));
            if (grid_supported && repair[r].source_claim[q]!=2)
                grid_supported[s]=repair[r].source_claim[q];
        }
        total_slots+=repair[r].n;
        supported+=repair[r].supported;
        claim_slots+=repair[r].claim_slots;
        applied_repairs++;
    }

done:
    free(repair); free(component); free(order); free(label); free(reserved);
    if (out_slots) *out_slots=total_slots;
    if (out_supported) *out_supported=supported;
    if (out_small_components) *out_small_components=small_components;
    if (out_enclosed_components) *out_enclosed_components=enclosed_components;
    if (out_claim_slots) *out_claim_slots=claim_slots;
    return applied_repairs;
}

static size_t rib_grid_both_diagonals_long(const float *grid,
                                           const uint8_t *present,
                                           size_t nk, size_t nu)
{
    const double gate2 = RIB_WRAP_GATE * RIB_WRAP_GATE;
    size_t count = 0;
    for (size_t k = 0; k + 1 < nk; k++) {
        for (size_t j = 0; j + 1 < nu; j++) {
            const float *a = &grid[(k*nu+j)*3];
            const float *b = a + 3;
            const float *c = a + nu*3;
            const float *d = c + 3;
            if (!present[k*nu+j] || !present[k*nu+j+1] ||
                !present[(k+1)*nu+j] || !present[(k+1)*nu+j+1])
                continue;
            if (rib_grid_d2(a,d) > gate2 && rib_grid_d2(b,c) > gate2)
                count++;
        }
    }
    return count;
}

/* ==================================================================== *
 *  Claims (search) ownership layer -- restored 2026-08-18 from a05f52a *
 *  (the 2026-08-10 checkpoint that produced pherc0139_4x5x5_0125,      *
 *  deleted by 1559997 "ownership by construction").  Active when       *
 *  component_global && !ownership_construction, which is the DEFAULT.  *
 * ==================================================================== */

typedef struct { double u, phi; } RibUPhi;

static int cmp_rib_uphi(const void *pa, const void *pb)
{
    const RibUPhi *a = (const RibUPhi *)pa;
    const RibUPhi *b = (const RibUPhi *)pb;
    if (a->u < b->u) return -1;
    if (a->u > b->u) return 1;
    return a->phi < b->phi ? -1 : (a->phi > b->phi ? 1 : 0);
}

/* Robust whole-component winding consensus for the fitted grid.
 *
 * A large welded disk can contribute several slice-chain fragments to the
 * same solved-u column.  The historical fitter simply kept the last claim,
 * so chain enumeration order selected different wraps on adjacent rows.  Bin
 * all samples by solved u, take the median authoritative phi in each column,
 * and isotonic-regress those medians into one component-global phi(u) map.
 * The map is used only to choose among competing 3-D observations; it does not
 * move samples, modify topology, or replace the true-arc-length solve. */
static double *phi_consensus_cover(Arena_T arena, RibUPhi *pair, size_t ns,
                                   double du, size_t nu, int winding_sense)
{
    size_t i = 0, j = 0, nb = 0;
    long first = -1, previous = -1;
    double *expect = (double *)ARENA_CALLOC(arena, nu, sizeof(double));
    double *weight = (double *)ARENA_CALLOC(arena, nu, sizeof(double));
    double *scratch = RIB_ALLOC_ARRAY(arena, double, ns);
    double *bval = RIB_ALLOC_ARRAY(arena, double, nu);
    double *bweight = RIB_ALLOC_ARRAY(arena, double, nu);
    size_t *bfirst = RIB_ALLOC_ARRAY(arena, size_t, nu);
    size_t *blast = RIB_ALLOC_ARRAY(arena, size_t, nu);
    qsort(pair, ns, sizeof(RibUPhi), cmp_rib_uphi);

    for (j = 0; j < nu; j++) {
        double lo = ((double)j - 0.5) * du;
        double hi = ((double)j + 0.5) * du;
        size_t begin = 0, end = 0, n = 0;
        while (i < ns && pair[i].u < lo) i++;
        begin = i;
        while (i < ns && (pair[i].u < hi ||
               (j + 1 == nu && pair[i].u <= hi))) i++;
        end = i;
        for (size_t q = begin; q < end; q++) scratch[n++] = pair[q].phi;
        if (n == 0) continue;
        qsort(scratch, n, sizeof(double), cmp_dbl);
        expect[j] = scratch[n / 2];
        weight[j] = (double)n;
        if (first < 0) first = (long)j;
    }
    if (first < 0) return expect;

    for (j = 0; j < (size_t)first; j++) expect[j] = expect[(size_t)first];
    previous = first;
    for (j = (size_t)first + 1; j < nu; j++) {
        if (weight[j] == 0.0) continue;
        if ((long)j > previous + 1) {
            double a = expect[(size_t)previous], b = expect[j];
            for (long q = previous + 1; q < (long)j; q++) {
                double t = (double)(q - previous) /
                           (double)((long)j - previous);
                expect[(size_t)q] = (1.0 - t) * a + t * b;
            }
        }
        previous = (long)j;
    }
    for (j = (size_t)previous + 1; j < nu; j++)
        expect[j] = expect[(size_t)previous];

    {
        double mj = 0.5 * (double)(nu - 1), mp = 0.0, cov = 0.0;
        double direction = winding_sense == -1 ? -1.0 : 1.0;
        for (j = 0; j < nu; j++) mp += expect[j];
        mp /= (double)nu;
        for (j = 0; j < nu; j++)
            cov += ((double)j - mj) * (expect[j] - mp);
        if (winding_sense == 0 && cov < 0.0) direction = -1.0;
        for (j = 0; j < nu; j++) {
            bval[nb] = direction * expect[j];
            bweight[nb] = weight[j] > 0.0 ? weight[j] : 1.0;
            bfirst[nb] = blast[nb] = j;
            nb++;
            while (nb >= 2 && bval[nb - 2] > bval[nb - 1]) {
                double w0 = bweight[nb - 2], w1 = bweight[nb - 1];
                bval[nb - 2] = (w0 * bval[nb - 2] + w1 * bval[nb - 1]) /
                               (w0 + w1);
                bweight[nb - 2] = w0 + w1;
                blast[nb - 2] = blast[nb - 1];
                nb--;
            }
        }
        for (i = 0; i < nb; i++)
            for (j = bfirst[i]; j <= blast[i]; j++)
                expect[j] = direction * bval[i];
    }
    return expect;
}

static double *component_phi_consensus(Arena_T arena, const SliceSet *S,
                                       double U0, double du, size_t nu,
                                       int winding_sense)
{
    size_t ns = 0;
    for (size_t i = 0; i < S->n_smp; i++) {
        int32_t ch = S->smp[i].chain;
        if (ch < 0 || !isfinite(S->smp[i].u) || !isfinite(S->smp[i].phi))
            continue;
        ns++;
    }
    if (ns == 0) return NULL;
    RibUPhi *pair = RIB_ALLOC_ARRAY(arena, RibUPhi, ns);
    ns = 0;
    for (size_t i = 0; i < S->n_smp; i++) {
        int32_t ch = S->smp[i].chain;
        if (ch < 0 || !isfinite(S->smp[i].u) || !isfinite(S->smp[i].phi))
            continue;
        pair[ns].u = S->smp[i].u - U0;
        pair[ns].phi = S->smp[i].phi;
        ns++;
    }
    return phi_consensus_cover(arena, pair, ns, du, nu, winding_sense);
}

static void rib_expect_digest(const char *width_name, const double *expect,
                              size_t active_nu)
{
    uint64_t h = UINT64_C(1469598103934665603);
    size_t count = 0;
    for (size_t j = 0; j < active_nu; j++) {
        uint64_t bits = 0;
        memcpy(&bits, &expect[j], sizeof(bits));
        h = rib_digest_mix(h, (uint64_t)j);
        h = rib_digest_mix(h, bits);
        count++;
    }
    fprintf(stderr,
            "  fit-width digest %-7s expected: values=%zu hash=%016llx\n",
            width_name, count, (unsigned long long)h);
}

static double rib_claim_previous_d2(const RibGridClaim *a, const float *p)
{
    double d0 = a->p[0] - (double)p[0];
    double d1 = a->p[1] - (double)p[1];
    double d2 = a->p[2] - (double)p[2];
    return d0*d0 + d1*d1 + d2*d2;
}

typedef struct {
    size_t breaks;
    size_t horizontal_breaks;
    size_t chain_switches;
    size_t chart_switches;
    double energy;
    int32_t parent;
} RibClaimPath;

static int rib_claim_compare(const void *pa, const void *pb)
{
    const RibGridClaim *a = (const RibGridClaim *)pa;
    const RibGridClaim *b = (const RibGridClaim *)pb;
    for (int d = 0; d < 3; d++) {
        if (a->p[d] < b->p[d]) return -1;
        if (a->p[d] > b->p[d]) return 1;
    }
    if (a->phi_score < b->phi_score) return -1;
    if (a->phi_score > b->phi_score) return 1;
    if (a->phi < b->phi) return -1;
    if (a->phi > b->phi) return 1;
    if (a->mesh_comp < b->mesh_comp) return -1;
    if (a->mesh_comp > b->mesh_comp) return 1;
    if (a->source_chart < b->source_chart) return -1;
    if (a->source_chart > b->source_chart) return 1;
    return a->chain < b->chain ? -1 : (a->chain > b->chain ? 1 : 0);
}

static int rib_claim_path_better(const RibClaimPath *candidate, int32_t tie,
                                 const RibClaimPath *current,
                                 int32_t current_tie)
{
    double eps;
    /* First minimize all physical discontinuities in the U/V quilt.  Among
     * paths with the same total, prefer not to tear the measured row itself;
     * then prefer one source curve over an A->B->A mosaic.  A switch forced by
     * the end of a chain's support remains available because every path through
     * that column pays it. */
    if (candidate->breaks != current->breaks)
        return candidate->breaks < current->breaks;
    if (candidate->horizontal_breaks != current->horizontal_breaks)
        return candidate->horizontal_breaks < current->horizontal_breaks;
    if (candidate->chain_switches != current->chain_switches)
        return candidate->chain_switches < current->chain_switches;
    if (candidate->chart_switches != current->chart_switches)
        return candidate->chart_switches < current->chart_switches;
    eps = 1e-12 * (1.0 + fabs(candidate->energy) +
                   fabs(current->energy));
    if (candidate->energy + eps < current->energy) return 1;
    if (current->energy + eps < candidate->energy) return 0;
    return tie < current_tie;
}

static RibClaimPath rib_claim_path_unary(const RibGridClaim *candidate,
                                         const float *previous,
                                         double gate2)
{
    RibClaimPath path;
    double phase = candidate->phi_score;
    memset(&path, 0, sizeof(path));
    for (int d = 0; d < 3; d++) {
        if (!isfinite(candidate->p[d]) ||
            (previous != NULL && !isfinite((double)previous[d]))) {
            fprintf(stderr,
                    "BUG: non-finite claimant geometry reached ribbon "
                    "ownership\n");
            abort();
        }
    }
    path.energy = phase * phase;
    path.parent = -1;
    if (previous != NULL) {
        double d2 = rib_claim_previous_d2(candidate, previous);
        if (d2 > gate2) path.breaks++;
        path.energy += d2;
    }
    return path;
}

static void rib_claim_path_add_edge(RibClaimPath *path,
                                    const RibGridClaim *left,
                                    const RibGridClaim *right,
                                    double gate2)
{
    double d2 = rib_claim_d2(left, right);
    if (d2 > gate2) {
        path->breaks++;
        path->horizontal_breaks++;
    }
    if (left->chain != right->chain) path->chain_switches++;
    if (left->source_chart != right->source_chart) path->chart_switches++;
    path->energy += d2;
}

/* Exact dynamic-programming solve over every contiguous claimed interval in a
 * row.  The old causal selector committed one column at a time, so a shorter
 * competing curve could win in the middle and produce an A->B->A splice even
 * when A itself covered both endpoints and the complete interval. */
static size_t rib_claim_select_row(Arena_T arena, RibGridClaim *claim,
                                   const size_t *off, size_t nu,
                                   const float *previous_row,
                                   const uint8_t *previous_present,
                                   int component_global,
                                   int32_t *selected)
{
    const double gate2 = RIB_WRAP_GATE * RIB_WRAP_GATE;
    RibClaimPath *state = NULL;
    size_t replacements = 0;

    for (size_t j = 0; j < nu; j++) selected[j] = -1;
    if (!component_global) {
        for (size_t j = 0; j < nu; j++)
            if (off[j + 1] > off[j])
                selected[j] = (int32_t)(off[j + 1] - 1);
        return 0;
    }
    for (size_t j = 0; j < nu; j++)
        if (off[j + 1] > off[j])
            qsort(claim + off[j], off[j + 1] - off[j],
                  sizeof(*claim), rib_claim_compare);
    if (off[nu] == 0) return 0;
    state = RIB_ALLOC_ARRAY(arena, RibClaimPath, off[nu]);

    for (size_t first = 0; first < nu;) {
        size_t last, begin, end;
        while (first < nu && off[first] == off[first + 1]) first++;
        if (first == nu) break;
        last = first;
        begin = off[first];
        end = off[first + 1];
        for (size_t q = begin; q < end; q++)
            state[q] = rib_claim_path_unary(
                &claim[q], previous_row != NULL &&
                             (previous_present == NULL || previous_present[first])
                              ? previous_row + first * 3 : NULL,
                gate2);

        while (last + 1 < nu && off[last + 1] < off[last + 2]) {
            size_t column = last + 1;
            size_t pbegin = off[last], pend = off[last + 1];
            size_t qbegin = off[column], qend = off[column + 1];
            for (size_t q = qbegin; q < qend; q++) {
                RibClaimPath unary = rib_claim_path_unary(
                    &claim[q], previous_row != NULL &&
                                 (previous_present == NULL ||
                                  previous_present[column])
                                  ? previous_row + column * 3 : NULL,
                    gate2);
                RibClaimPath best = {
                    SIZE_MAX, SIZE_MAX, SIZE_MAX, SIZE_MAX, HUGE_VAL, -1
                };
                for (size_t p = pbegin; p < pend; p++) {
                    RibClaimPath next = state[p];
                    rib_claim_path_add_edge(&next, &claim[p], &claim[q],
                                            gate2);
                    next.breaks += unary.breaks;
                    next.energy += unary.energy;
                    next.parent = (int32_t)p;
                    if (rib_claim_path_better(
                            &next, (int32_t)p, &best, best.parent))
                        best = next;
                }
                state[q] = best;
            }
            last = column;
        }

        begin = off[last];
        end = off[last + 1];
        int32_t best = (int32_t)begin;
        for (size_t q = begin + 1; q < end; q++)
            if (rib_claim_path_better(
                    &state[q], (int32_t)q, &state[(size_t)best], best))
                best = (int32_t)q;
        for (size_t j = last;; j--) {
            selected[j] = best;
            best = state[(size_t)best].parent;
            if (j == first) break;
        }
        first = last + 1;
    }
    for (size_t j = 0; j < nu; j++)
        if (selected[j] >= 0 && (size_t)selected[j] != off[j])
            replacements++;
    return replacements;
}

/* Exact global-objective contribution of every cell and grid edge incident to
 * either of two slots.  Callers temporarily substitute both points in the
 * grid; the union (rather than two separately counted neighborhoods) makes a
 * strict improvement monotone in the full emitted-triangle objective. */
static void rib_grid_pair_objective(const float *grid, size_t nk, size_t nu,
                                    size_t slot_a, size_t slot_b,
                                    int *out_legal, size_t *out_bad)
{
    const double gate2 = RIB_WRAP_GATE * RIB_WRAP_GATE;
    size_t ka = slot_a / nu, ja = slot_a - ka*nu;
    size_t kb = slot_b / nu, jb = slot_b - kb*nu;
    size_t kmin = ka < kb ? ka : kb, kmax = ka > kb ? ka : kb;
    size_t jmin = ja < jb ? ja : jb, jmax = ja > jb ? ja : jb;
    size_t ck0 = kmin > 0 ? kmin - 1 : 0;
    size_t cj0 = jmin > 0 ? jmin - 1 : 0;
    size_t ck1 = kmax < nk - 1 ? kmax : nk - 2;
    size_t cj1 = jmax < nu - 1 ? jmax : nu - 2;
    int legal = 0;
    size_t bad = 0;
    const size_t slot[2] = { slot_a, slot_b };
    const long dk[4] = { 0, 0, -1, 1 };
    const long dj[4] = { -1, 1, 0, 0 };

    for (size_t ck = ck0; ck <= ck1; ck++)
        for (size_t cj = cj0; cj <= cj1; cj++)
            legal += rib_grid_cell_legal(
                grid, nu, ck, cj, SIZE_MAX, NULL);

    for (int s = 0; s < 2; s++) {
        size_t k = slot[s] / nu, j = slot[s] - k*nu;
        const float *p = &grid[slot[s]*3];
        for (int q = 0; q < 4; q++) {
            long kk = (long)k + dk[q], jj = (long)j + dj[q];
            size_t neighbor;
            const float *r;
            if (kk < 0 || jj < 0 ||
                (size_t)kk >= nk || (size_t)jj >= nu)
                continue;
            neighbor = (size_t)kk*nu + (size_t)jj;
            if (s == 1 && neighbor == slot_a) continue;
            r = &grid[neighbor*3];
            if (!isfinite((double)r[0])) continue;
            if (rib_grid_d2(p, r) > gate2) bad++;
        }
    }
    *out_legal = legal;
    *out_bad = bad;
}

/* Gather the real slice-chain positions that can claim one fitted-grid slot.
 * The already fitted point is candidate zero, so a two-slot move may retain
 * either endpoint.  Duplicate interpolants are removed after float rounding:
 * the grid stores float, and scoring a higher-precision point which rounds back
 * to the current value creates false progress at a gate threshold. */
static size_t rib_grid_slot_candidates(const SliceSet *S,
                                       double U0, double du,
                                       size_t k, size_t j,
                                       const int32_t *choff,
                                       const int32_t *chids,
                                       const float current[3],
                                       float candidate[RIB_PAIR_CAND_MAX][3])
{
    size_t n = 1;
    double uq = (double)j * du + U0;
    memcpy(candidate[0], current, 3 * sizeof(float));
    for (int32_t x = choff[k]; x < choff[k+1]; x++) {
        const Chain *ch = &S->chn[chids[x]];
        int32_t f = ch->first, cn = ch->count;
        int32_t lo, hi, seg;
        float p[3];
        int duplicate = 0;
        if (cn < 2 || uq < S->smp[f].u ||
            uq > S->smp[f+cn-1].u)
            continue;
        lo = 1; hi = cn - 1;
        while (lo < hi) {
            int32_t mid = lo + (hi - lo) / 2;
            if (S->smp[f+mid].u < uq) lo = mid + 1;
            else hi = mid;
        }
        seg = lo;
        {
            double ua = S->smp[f+seg-1].u;
            double ub = S->smp[f+seg].u;
            double w = ub - ua > 1e-12 ? (uq - ua) / (ub - ua) : 0.0;
            if (w < 0.0) w = 0.0;
            if (w > 1.0) w = 1.0;
            for (int d = 0; d < 3; d++)
                p[d] = (float)((1.0-w)*S->smp[f+seg-1].p[d] +
                                      w *S->smp[f+seg].p[d]);
        }
        for (size_t q = 0; q < n; q++)
            if (p[0] == candidate[q][0] &&
                p[1] == candidate[q][1] &&
                p[2] == candidate[q][2]) {
                duplicate = 1;
                break;
            }
        if (!duplicate && n < RIB_PAIR_CAND_MAX) {
            memcpy(candidate[n], p, 3 * sizeof(float));
            n++;
        }
    }
    return n;
}

/* The row solve is intentionally smooth, but a contested (U,V) cell may have
 * inherited the wrong physical claimant.  Revisit only real slice claims and
 * accept a substitution iff it strictly increases the exact number of locally
 * emittable triangles, or preserves that number while removing a long grid
 * edge.  Each change therefore decreases a global integer defect objective;
 * sweeps cannot oscillate or make the accepted grid worse. */
static size_t rib_grid_repair_claims(const SliceSet *S, double U0, double du,
                                     size_t nu, size_t nk,
                                     const int32_t *choff,
                                     const int32_t *chids, float *grid)
{
    size_t repairs = 0;
    for (int sweep = 0; sweep < 4; sweep++) {
        size_t changed = 0;
        int reverse = sweep & 1;
        for (size_t rk = 0; rk < nk; rk++) {
            size_t k = reverse ? nk-1-rk : rk;
            for (int32_t x = choff[k]; x < choff[k+1]; x++) {
                const Chain *ch = &S->chn[chids[x]];
                int32_t f=ch->first, cn=ch->count;
                if (cn < 2) continue;
                double ulo=S->smp[f].u-U0, uhi=S->smp[f+cn-1].u-U0;
                int32_t j0=(int32_t)ceil(ulo/du), j1=(int32_t)floor(uhi/du);
                if (j0 < 0) j0=0;
                if (j1 > (int32_t)nu-1) j1=(int32_t)nu-1;
                int32_t seg=1;
                for (int32_t jj=j0; jj<=j1; jj++) {
                    size_t j=(size_t)jj, slot=k*nu+j;
                    float *g=&grid[slot*3];
                    if (!isfinite((double)g[0])) continue;
                    double uq=(double)jj*du+U0;
                    while (seg < cn-1 && S->smp[f+seg].u < uq) seg++;
                    double ua=S->smp[f+seg-1].u, ub=S->smp[f+seg].u;
                    double w=ub-ua > 1e-12 ? (uq-ua)/(ub-ua) : 0.0;
                    double candidate[3], current[3] = {
                        (double)g[0], (double)g[1], (double)g[2] };
                    if (w < 0.0) w=0.0;
                    if (w > 1.0) w=1.0;
                    for (int d=0; d<3; d++)
                        candidate[d]=(1.0-w)*S->smp[f+seg-1].p[d] +
                                           w *S->smp[f+seg].p[d];
                    size_t current_bad=rib_grid_incident_long(
                        grid,nk,nu,k,j,current);
                    int current_legal=rib_grid_local_legal(
                        grid,nk,nu,k,j,current);
                    int candidate_legal=rib_grid_local_legal(
                        grid,nk,nu,k,j,candidate);
                    size_t candidate_bad=rib_grid_incident_long(
                        grid,nk,nu,k,j,candidate);
                    if (candidate_legal < current_legal ||
                        (candidate_legal == current_legal &&
                         candidate_bad >= current_bad))
                        continue;
                    for (int d=0; d<3; d++) g[d]=(float)candidate[d];
                    changed++; repairs++;
                }
            }
        }
        if (changed == 0) break;
    }
    return repairs;
}

/* A one-vertex local optimum can still surround a square puncture: two
 * neighboring wrong claimants may need to switch together before either move
 * becomes legal.  Search only defective four-corner cells, enumerate bounded
 * pairs of real slice claims, and accept a simultaneous move iff it increases
 * the exact number of emittable triangles in the union neighborhood (or keeps
 * that number while strictly removing an incident long edge). */
static size_t rib_grid_repair_claim_pairs(const SliceSet *S,
                                          double U0, double du,
                                          size_t nu, size_t nk,
                                          const int32_t *choff,
                                          const int32_t *chids,
                                          float *grid)
{
    static const int corner_dk[4] = { 0, 0, 1, 1 };
    static const int corner_dj[4] = { 0, 1, 0, 1 };
    size_t repairs = 0;
    for (int sweep = 0; sweep < 4; sweep++) {
        size_t changed = 0;
        int reverse = sweep & 1;
        for (size_t rk = 0; rk + 1 < nk; rk++) {
            size_t k = reverse ? nk - 2 - rk : rk;
            for (size_t rj = 0; rj + 1 < nu; rj++) {
                size_t j = reverse ? nu - 2 - rj : rj;
                size_t cell_slot[4] = {
                    k*nu+j, k*nu+j+1,
                    (k+1)*nu+j, (k+1)*nu+j+1
                };
                float candidates[4][RIB_PAIR_CAND_MAX][3];
                size_t ncandidate[4];
                int current_cell;
                int best_delta_legal = 0;
                size_t best_bad_drop = 0;
                int best_a = -1, best_b = -1;
                float best_pa[3] = { 0,0,0 }, best_pb[3] = { 0,0,0 };

                if (!isfinite((double)grid[cell_slot[0]*3]) ||
                    !isfinite((double)grid[cell_slot[1]*3]) ||
                    !isfinite((double)grid[cell_slot[2]*3]) ||
                    !isfinite((double)grid[cell_slot[3]*3]))
                    continue;
                current_cell = rib_grid_cell_legal(
                    grid, nu, k, j, SIZE_MAX, NULL);
                if (current_cell >= 2) continue;
                for (int q = 0; q < 4; q++) {
                    size_t ck = k + (size_t)corner_dk[q];
                    size_t cj = j + (size_t)corner_dj[q];
                    ncandidate[q] = rib_grid_slot_candidates(
                        S,U0,du,ck,cj,choff,chids,
                        &grid[cell_slot[q]*3],candidates[q]);
                }
                for (int a = 0; a < 4; a++) {
                    for (int b = a + 1; b < 4; b++) {
                        size_t sa = cell_slot[a], sb = cell_slot[b];
                        float old_a[3], old_b[3];
                        int old_legal;
                        size_t old_bad;
                        memcpy(old_a,&grid[sa*3],3*sizeof(float));
                        memcpy(old_b,&grid[sb*3],3*sizeof(float));
                        rib_grid_pair_objective(
                            grid,nk,nu,sa,sb,&old_legal,&old_bad);
                        for (size_t ca = 0; ca < ncandidate[a]; ca++) {
                            for (size_t cb = 0; cb < ncandidate[b]; cb++) {
                                int next_legal, delta_legal;
                                size_t next_bad, bad_drop;
                                if (ca == 0 && cb == 0) continue;
                                memcpy(&grid[sa*3],candidates[a][ca],
                                       3*sizeof(float));
                                memcpy(&grid[sb*3],candidates[b][cb],
                                       3*sizeof(float));
                                rib_grid_pair_objective(
                                    grid,nk,nu,sa,sb,
                                    &next_legal,&next_bad);
                                delta_legal = next_legal - old_legal;
                                bad_drop = next_bad < old_bad
                                         ? old_bad - next_bad : 0;
                                if (delta_legal > best_delta_legal ||
                                    (delta_legal == best_delta_legal &&
                                     delta_legal >= 0 &&
                                     bad_drop > best_bad_drop)) {
                                    best_delta_legal = delta_legal;
                                    best_bad_drop = bad_drop;
                                    best_a = a; best_b = b;
                                    memcpy(best_pa,candidates[a][ca],
                                           3*sizeof(float));
                                    memcpy(best_pb,candidates[b][cb],
                                           3*sizeof(float));
                                }
                            }
                        }
                        memcpy(&grid[sa*3],old_a,3*sizeof(float));
                        memcpy(&grid[sb*3],old_b,3*sizeof(float));
                    }
                }
                if (best_a >= 0 &&
                    (best_delta_legal > 0 || best_bad_drop > 0)) {
                    memcpy(&grid[cell_slot[best_a]*3],best_pa,3*sizeof(float));
                    memcpy(&grid[cell_slot[best_b]*3],best_pb,3*sizeof(float));
                    changed++; repairs++;
                }
            }
        }
        fprintf(stderr,
                "  fitted-grid paired claim repair sweep %d: %zu update(s)\n",
                sweep + 1, changed);
        if (changed == 0) break;
    }
    return repairs;
}

static void fit_ribbon(Arena_T arena, const SliceSet *S, double U0,
                       const RibbonOpts *o,
                       const ChainRelationGraph *relations,
                       RibbonResult *out)
{
    double du = (double)o->grid_u;
    if (du < 1e-6) du = 1.0;
    double atlas_umax = -1e300, cover_umax = -1e300;
    for (size_t i = 0; i < S->n_smp; i++) {
        if (S->smp[i].chain < 0) continue;
        int32_t ch = S->smp[i].chain;
        double atlas_u = S->smp[i].u - U0;
        double cover_u = atlas_u - S->chn[ch].atlas_shift;
        if (atlas_u > atlas_umax) atlas_umax = atlas_u;
        if (cover_u > cover_umax) cover_umax = cover_u;
    }
    double umax = o->fit_cover_width ? cover_umax : atlas_umax;
    if (umax <= 0.0) return;
    size_t active_nu = (size_t)(atlas_umax / du) + 2;
    size_t nu = (size_t)(umax / du) + 2;
    size_t nk = (size_t)S->nplanes;
    if (nu < 2 || nk < 1) return;
    {
        /* A CARRIED registered field legitimately spans the whole unrolled
         * scroll (~1.2M vox u on PHerc0139 -> ~600k columns at du=2), so the
         * old fixed 500k-column runaway guard rejected every correct
         * whole-scroll fit (silently, before 2026-08-11).  Guard on the
         * actual allocation instead: the grid is
         * nk*nu*(12B G + 1B GV + 4B lifted phase) plus
         * ~nu*100B of row scratch. */
        double grid_bytes = (double)nk * (double)nu * 17.0
                          + (double)nu * 100.0;
        const double grid_cap = 48.0e9;
        if (grid_bytes > grid_cap) {
            fprintf(stderr,
                "ERROR: fitted grid %zu cols x %zu rows needs %.1f GB "
                "(cap %.0f GB; du=%.3g, u span %.6g, atlas_umax=%.6g, "
                "cover_umax=%.6g). Raise du, tighten island packing/--gap, "
                "or check the u-tail winsorization.\n",
                nu, nk, grid_bytes / 1e9, grid_cap / 1e9, du, umax,
                atlas_umax, cover_umax);
            return;
        }
        if (nu > 500000)
            fprintf(stderr,
                "  fitted grid is whole-scroll scale: %zu cols x %zu rows "
                "(%.1f GB)\n", nu, nk, grid_bytes / 1e9);
    }

    float   *G  = (float *)ARENA_CALLOC(arena, nk * nu * 3, sizeof(float));
    uint8_t *GO = (uint8_t *)ARENA_CALLOC(arena, nk * nu, 1);
    uint8_t *GV = (uint8_t *)ARENA_CALLOC(arena, nk * nu, 1);
    float   *GP = (float *)ARENA_CALLOC(arena, nk * nu, sizeof(float));
    int32_t *GI = RIB_ALLOC_ARRAY(arena, int32_t, nk * nu);
    int32_t *GM = RIB_ALLOC_ARRAY(arena, int32_t, nk * nu);
    for (size_t i = 0; i < nk * nu; i++) {
        GI[i] = -1;
        GM[i] = -1;
    }

    /* The chain-relation graph remains a gauge-solve diagnostic; it is not an
     * emission input.  Continuity enters the fit only through the merge
     * partition and the serialized direct chart relations. */
    (void)relations;

    /* Built before the scratch mark: the merge partition and chart-support
     * table must survive until the post-export transition audit. */
    RibMergeSet merge;
    memset(&merge, 0, sizeof merge);
    double fit_t0 = ves_clock_sec();
    if (o->component_global && o->ownership_construction) {
        if (rib_merge_build(arena, S, o, &merge) != 0) return;
        out->grid_merge_junctions += merge.n_junctions;
        out->grid_merge_ambiguous += merge.ambiguous;
        out->grid_merge_rejected_island += merge.rejected_island;
        out->grid_merge_rejected_relation += merge.rejected_relation;
        out->grid_merge_rejected_du += merge.rejected_du;
        out->grid_merge_rejected_gap += merge.rejected_gap;
        out->grid_merge_rejected_radial += merge.rejected_radial;
        out->grid_merge_rejected_phase += merge.rejected_phase;
        fprintf(stderr,
                "    [fit] chain merges: %zu junction(s), %zu logical "
                "chain(s) from %zu, ambiguous=%zu, rejected "
                "island=%zu relation=%zu du=%zu gap=%zu radial=%zu "
                "phase=%zu (%.2fs)\n",
                merge.n_junctions, merge.n_logical, S->n_chn,
                merge.ambiguous, merge.rejected_island,
                merge.rejected_relation, merge.rejected_du,
                merge.rejected_gap, merge.rejected_radial,
                merge.rejected_phase, ves_clock_sec() - fit_t0);
        fit_t0 = ves_clock_sec();
    }

    Arena_Mark mark = Arena_save(arena);
    int np = S->nplanes;
    int32_t *choff = NULL;
    int32_t *chids = chains_by_slice(arena, S, &choff);
    RibClaimLabelRow *claim_labels = NULL;
    const double *expected_phi = NULL;

    if (o->component_global && o->ownership_construction) {
        /* Freeze ownership by construction before fitting: certified merges
         * define continuity, a row-invariant key resolves overlap, and
         * violations are reported rather than optimized away.  Raster values
         * are payload only: neither the emitted texture nor quilting
         * diagnostics enter this decision. */
        if (rib_fuse_select_rows(
                S, U0, du, nu, nk, choff, chids, &merge, out,
                &claim_labels) != 0) {
            Arena_restore(arena, mark);
            return;
        }
        fprintf(stderr, "    [fit] merge fusion: %.2fs "
                "(same-logical=%zu alternate=%zu violations=%zu)\n",
                ves_clock_sec() - fit_t0, out->grid_fuse_same_logical,
                out->grid_fuse_alternate, out->grid_fuse_violations);
        fit_t0 = ves_clock_sec();
    } else if (o->component_global) {
        /* Claims mode: the component-wide phi consensus is a soft unary for
         * an exact whole-row claimant path. */
        expected_phi = component_phi_consensus(arena, S, U0, du, nu,
                                               o->winding_sense);
        if (o->verify_fit_width && expected_phi != NULL)
            rib_expect_digest(o->fit_cover_width ? "cover" : "compact",
                              expected_phi, active_nu < nu ? active_nu : nu);
        fit_t0 = ves_clock_sec();
    }

    /* ---- v-strips: partition the slice planes into axial strips.  Every
     * claimant is frozen (graph mode) or a pure per-column function (legacy
     * pick-last), so strips own disjoint fitted rows and are semantically
     * invisible; they exist only for thread parallelism. */
    int strip_h = np;
    if (!o->direct_ribbon && o->v_strip_planes > 0) {
        strip_h = o->v_strip_planes;
    } else if (!o->direct_ribbon && o->v_strip_planes < 0) {
        /* Target up to 2 strips per thread but never more than 64 strips
         * total, with a 16-plane floor. */
        int nthreads = 1;
        int n_target = 0;
#ifdef _OPENMP
        nthreads = omp_get_max_threads();
        if (nthreads < 1) nthreads = 1;
#endif
        n_target = 2 * nthreads < 64 ? 2 * nthreads : 64;
        strip_h = np / (n_target > 0 ? n_target : 1);
        if (strip_h < 16) strip_h = 16;
        if (strip_h > 64) strip_h = 64;
    }
    int n_strips = (np + strip_h - 1) / strip_h;
    if (n_strips > 1)
        fprintf(stderr, "  fitted grid v-strips: %d strip(s) of %d plane(s)\n",
                 n_strips, strip_h);
    else if (o->direct_ribbon)
        fprintf(stderr,
                "  fitted grid ownership: one continuous V pass (%d planes)\n",
                np);
    {
    size_t tot_conflicts = 0, tot_conflict_cells = 0;
    size_t tot_discarded = 0, tot_replaced = 0;
    size_t tot_ws = 0, tot_wsm = 0, tot_wss = 0, tot_wsi = 0;
    int claim_handoff_failed = 0;
    int si = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) \
    reduction(+:tot_conflicts) reduction(+:tot_conflict_cells) \
    reduction(+:tot_discarded) reduction(+:tot_replaced) \
    reduction(+:tot_ws) reduction(+:tot_wsm) \
    reduction(+:tot_wss) reduction(+:tot_wsi) \
    reduction(|:claim_handoff_failed)
#endif
    for (si = 0; si < n_strips; si++) {
    Arena_T sarena = Arena_new();
    double *P   = RIB_ALLOC_ARRAY(sarena, double, nu * 3);
    uint8_t *V  = RIB_ALLOC_ARRAY(sarena, uint8_t, nu);
    int32_t *claim_chain = RIB_ALLOC_ARRAY(sarena, int32_t, nu);
    int k0 = si * strip_h;
    int k1 = k0 + strip_h < np ? k0 + strip_h : np;
    /* Empty columns are absent observations, represented only by V.  Keep
     * scratch geometry finite so an accidental read cannot poison a solve. */
    memset(P, 0, nu * 3 * sizeof(*P));

    for (int k = k0; k < k1; k++) {
        Arena_Mark row_mark = Arena_save(sarena);
        float *gpos_row = &G[(size_t)k * nu * 3];
        /* The row's chains occupy a narrow u WINDOW of the whole-scroll grid
         * (a 21x21x21 grid is ~900k columns; one slice plane's chains span a
         * few percent of it).  Every per-row structure and scan below works
         * on [rjmin, rjmax] only -- the old full-width CSR calloc, memsets
         * and column scans were O(nu) per row and dominated the whole fit at
         * scroll scale.  G/GV writes remain at absolute column indices. */
        size_t rjmin = nu, rjmax = 0;
        for (int32_t x = choff[k]; x < choff[k+1]; x++) {
            const Chain *ch = &S->chn[chids[x]];
            int32_t f = ch->first, cn = ch->count;
            if (cn < 2) continue;
            double ulo = S->smp[f].u - U0, uhi = S->smp[f+cn-1].u - U0;
            int32_t j0 = (int32_t)ceil(ulo / du), j1 = (int32_t)floor(uhi / du);
            if (j0 < 0) j0 = 0;
            if (j1 > (int32_t)nu - 1) j1 = (int32_t)nu - 1;
            if (j0 > j1) continue;
            if ((size_t)j0 < rjmin) rjmin = (size_t)j0;
            if ((size_t)j1 > rjmax) rjmax = (size_t)j1;
        }
        if (rjmin > rjmax) {
            Arena_restore(sarena, row_mark);
            continue;
        }
        size_t rw = rjmax - rjmin + 1;
        size_t *off = (size_t *)ARENA_CALLOC(sarena, rw + 1, sizeof(size_t));
        memset(V + rjmin, 0, rw);
        for (size_t j = 0; j < rw; j++) claim_chain[rjmin + j] = -1;

        /* Column CSR (window-shifted), pass 1: count every physical claimant.
         * Gathering first makes conflict accounting explicit while preserving
         * chain order for the established causal continuity selector. */
        for (int32_t x = choff[k]; x < choff[k+1]; x++) {
            const Chain *ch = &S->chn[chids[x]];
            int32_t f = ch->first, cn = ch->count;
            if (cn < 2) continue;
            double ulo = S->smp[f].u - U0, uhi = S->smp[f+cn-1].u - U0;
            int32_t j0 = (int32_t)ceil(ulo / du), j1 = (int32_t)floor(uhi / du);
            if (j0 < 0) j0 = 0;
            if (j1 > (int32_t)nu - 1) j1 = (int32_t)nu - 1;
            for (int32_t j = j0; j <= j1; j++)
                off[(size_t)j - rjmin + 1]++;
        }
        for (size_t j = 0; j < rw; j++) off[j + 1] += off[j];
        size_t nclaim = off[rw];
        if (nclaim == 0) {
            Arena_restore(sarena, row_mark);
            continue;
        }

        RibGridClaim *claim = RIB_ALLOC_ARRAY(sarena, RibGridClaim, nclaim);
        size_t *cursor = RIB_ALLOC_ARRAY(sarena, size_t, rw);
        memcpy(cursor, off, rw * sizeof(size_t));

        /* Pass 2: interpolate each chain at the column isovalue. */
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
                RibGridClaim *dst = &claim[cursor[(size_t)j - rjmin]++];
                if (w < 0.0) w = 0.0;
                if (w > 1.0) w = 1.0;
                for (int d = 0; d < 3; d++) {
                    dst->p[d] = (1.0-w)*S->smp[f+seg-1].p[d]
                              +       w *S->smp[f+seg].p[d];
                    dst->tau[d] = (1.0-w)*S->smp[f+seg-1].tau[d]
                                +       w *S->smp[f+seg].tau[d];
                }
                {
                double candidate_phi = (1.0-w)*S->smp[f+seg-1].phi
                                         +       w *S->smp[f+seg].phi;
                    dst->phi = candidate_phi;
                    dst->phi_score = expected_phi != NULL
                                   ? fabs(candidate_phi -
                                          expected_phi[(size_t)j]) : 0.0;
                }
                dst->mesh_comp = ch->mesh_comp;
                dst->chain = chids[x];
                dst->source_chart = ch->source_chart;
            }
        }

        int32_t *selected = RIB_ALLOC_ARRAY(sarena, int32_t, rw);
        size_t row_replaced = 0;
        if (o->component_global && o->ownership_construction) {
            /* Consume the frozen fusion decision by stable chain identity.
             * Raster coordinates are output data, never a second matching
             * domain. */
            if (rib_claim_select_labels(
                    claim, off, rw, &claim_labels[k], rjmin,
                    selected, &row_replaced) != 0) {
                claim_handoff_failed = 1;
                row_replaced = 0;
                for (size_t jw = 0; jw < rw; jw++) selected[jw] = -1;
            }
        } else if (o->component_global) {
            /* Claims mode: solve each contiguous row interval as one exact
             * claimant path against the fitted row above. */
            const float *previous_row = k > k0
                ? &G[(((size_t)k - 1) * nu + rjmin) * 3] : NULL;
            const uint8_t *previous_present = k > k0
                ? &GO[((size_t)k - 1) * nu + rjmin] : NULL;
            row_replaced = rib_claim_select_row(
                sarena, claim, off, rw, previous_row, previous_present,
                1, selected);
        } else {
            /* Legacy single-mesh path: the historical last-claim-per-column
             * rule (chain enumeration order), unchanged. */
            for (size_t jw = 0; jw < rw; jw++)
                selected[jw] = off[jw + 1] > off[jw]
                             ? (int32_t)(off[jw + 1] - 1) : -1;
        }
        int rowvalid = 0;
        for (size_t jw = 0; jw < rw; jw++) {
            size_t j = rjmin + jw;
            size_t nc = off[jw + 1] - off[jw];
            if (nc > 1) {
                tot_conflicts += nc - 1;
                tot_conflict_cells++;
                if (o->discard_conflicting_claims) {
                    /* A conflict is evidence that this UV location has no
                     * trustworthy XYZ observation.  Do not arbitrate it: the
                     * variational continuation below may bridge the missing
                     * cell, but no claimant is allowed into its data term. */
                    selected[jw] = -1;
                    tot_discarded++;
                }
            }
            if (selected[jw] < 0) continue;
            const RibGridClaim *src = &claim[(size_t)selected[jw]];
            for (int d = 0; d < 3; d++)
                P[j*3 + (size_t)d] = src->p[d];
            V[j] = 1;
            claim_chain[j] = src->chain;
            GP[(size_t)k*nu + j] = (float)src->phi;
            GI[(size_t)k*nu + j] =
                S->chn[src->chain].reconstruction_component;
            GM[(size_t)k*nu + j] = S->chn[src->chain].winding_island;
            if (GI[(size_t)k*nu + j] < 0)
                claim_handoff_failed = 1;
            rowvalid++;
        }
        if (!o->discard_conflicting_claims) tot_replaced += row_replaced;
        if (rowvalid == 0) {
            Arena_restore(sarena, row_mark);
            continue;
        }

        /* A slice may contain several physically separate U runs.  Preserve
         * exact measured claims and bridge only short, geometrically local
         * gaps with a convex chord.  In particular, never integrate tangents
         * across missing data: that was able to leave the convex hull of the
         * input by thousands of voxels.  Any unresolved local gap is handled
         * later by removing its U column globally, never by deleting a patch. */
        {
            size_t jend = rjmax + 1;
            size_t jseek = rjmin;
            const double gate2 = RIB_WRAP_GATE * RIB_WRAP_GATE;
            while (jseek < jend) {
                size_t jf, jl, scan;
                while (jseek < jend && !V[jseek]) jseek++;
                if (jseek >= jend) break;
                jf = jl = jseek;
                scan = jseek + 1;
                for (;;) {
                    int physical_break = 0;
                    while (scan < jend && V[scan]) {
                        double d2 = 0.0;
                        for (int d = 0; d < 3; d++) {
                            double delta = P[jl*3 + (size_t)d]
                                         - P[scan*3 + (size_t)d];
                            d2 += delta * delta;
                        }
                        if (d2 > gate2) {
                            int32_t ca = claim_chain[jl];
                            int32_t cb = claim_chain[scan];
                            tot_ws++;
                            if (ca >= 0 && cb >= 0) {
                                const Chain *a = &S->chn[ca];
                                const Chain *b = &S->chn[cb];
                                if (a->mesh_comp == b->mesh_comp) tot_wsm++;
                                if (a->solve_comp >= 0 &&
                                    a->solve_comp == b->solve_comp) tot_wss++;
                                if (a->reconstruction_component >= 0 &&
                                    a->reconstruction_component ==
                                        b->reconstruction_component) tot_wsi++;
                            }
                            physical_break = 1;
                            break;
                        }
                        jl = scan++;
                    }
                    if (physical_break) break;
                    {
                        size_t gap_first = scan;
                        double d2 = 0.0;
                        double steps;
                        while (scan < jend && !V[scan]) scan++;
                        if (scan >= jend ||
                            scan - gap_first > (size_t)RIB_UFILL_MAX)
                            break;
                        steps = (double)(scan - jl);
                        for (int d = 0; d < 3; d++) {
                            double delta = P[jl*3 + (size_t)d]
                                         - P[scan*3 + (size_t)d];
                            d2 += delta * delta;
                        }
                        if (d2 > gate2 * steps * steps) break;
                    }
                    jl = scan++;
                }
                jseek = jl + 1;

                for (size_t j = jf; j <= jl; j++) {
                    if (!V[j]) continue;
                    for (int d = 0; d < 3; d++)
                        gpos_row[j*3 + (size_t)d] =
                            (float)P[j*3 + (size_t)d];
                    GO[(size_t)k*nu + j] = 1;
                    GV[(size_t)k*nu + j] = 1;
                }
                for (size_t gap = jf; gap <= jl;) {
                    size_t first_gap, right, left;
                    if (V[gap]) { gap++; continue; }
                    first_gap = gap;
                    while (gap <= jl && !V[gap]) gap++;
                    if (first_gap == jf || gap > jl) break;
                    left = first_gap - 1;
                    right = gap;
                    for (size_t q = first_gap; q < right; q++) {
                        double w = (double)(q - left) /
                                   (double)(right - left);
                        for (int d = 0; d < 3; d++)
                            gpos_row[q*3 + (size_t)d] =
                                (float)((1.0 - w) * P[left*3 + (size_t)d] +
                                             w  * P[right*3 + (size_t)d]);
                        GP[(size_t)k*nu + q] =
                            (float)((1.0 - w) * GP[(size_t)k*nu + left] +
                                         w  * GP[(size_t)k*nu + right]);
                        if (GI[(size_t)k*nu + left] >= 0 &&
                            GI[(size_t)k*nu + left] ==
                                GI[(size_t)k*nu + right])
                            GI[(size_t)k*nu + q] =
                                GI[(size_t)k*nu + left];
                        if (GM[(size_t)k*nu + left] >= 0 &&
                            GM[(size_t)k*nu + left] ==
                                GM[(size_t)k*nu + right])
                            GM[(size_t)k*nu + q] =
                                GM[(size_t)k*nu + left];
                        GO[(size_t)k*nu + q] = 1;
                    }
                }
            }
        }
        Arena_restore(sarena, row_mark);
    }
    Arena_dispose(&sarena);
    }
    if (claim_handoff_failed) {
        fprintf(stderr,
                "ERROR: claimant graph handoff changed while rebuilding a "
                "slice; refusing a raster fallback\n");
        rib_claim_label_rows_dispose(claim_labels, nk);
        claim_labels = NULL;
        Arena_restore(arena, mark);
        return;
    }
    out->grid_claim_conflicts += tot_conflicts;
    out->grid_claim_conflict_cells += tot_conflict_cells;
    out->grid_claim_discarded_cells += tot_discarded;
    out->grid_claim_replaced += tot_replaced;
    out->grid_row_wrap_splits += tot_ws;
    out->grid_row_wrap_same_mesh += tot_wsm;
    out->grid_row_wrap_same_solve += tot_wss;
    out->grid_row_wrap_same_island += tot_wsi;
    }
    fprintf(stderr, "    [fit] row loop: %.2fs\n",
            ves_clock_sec() - fit_t0);
    fit_t0 = ves_clock_sec();
    out->grid_long_edges_rowfit = rib_grid_long_edges(G, GO, nk, nu);
    fprintf(stderr, "    [fit] long-edge count A: %.2fs\n",
            ves_clock_sec() - fit_t0);
    fit_t0 = ves_clock_sec();
    if (o->component_global)
        fprintf(stderr,
                "    [fit] claimant repair disabled: %s\n",
                o->discard_conflicting_claims
                    ? "conflicts remain unsupported for minimal continuation"
                    : (o->ownership_construction
                        ? "measured graph ownership is frozen"
                        : "whole-row claimant path is frozen"));
    out->grid_long_edges_repaired = rib_grid_long_edges(G,GO,nk,nu);
    fprintf(stderr, "    [fit] long-edge count B: %.2fs\n",
            ves_clock_sec() - fit_t0);
    fit_t0 = ves_clock_sec();
    if (o->verify_fit_width)
        rib_grid_digest(o->fit_cover_width ? "cover" : "compact", "row-fit",
                        G, GO, nk, nu, active_nu < nu ? active_nu : nu);

    {
        /* v-fill touches exactly one column per iteration: parallel over j is
         * race-free and bit-identical to the serial order. */
        int64_t pj = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1024)
#endif
    for (pj = 0; pj < (int64_t)nu; pj++) {
        size_t j = (size_t)pj;
        long ka = -1;
        for (size_t k = 0; k < nk; k++) {
            if (!GO[k*nu + j]) continue;
            long kb = (long)k;
            if (ka >= 0 && kb > ka + 1) {
                double d2 = 0.0;
                double steps = (double)(kb - ka);
                int32_t ia = GI[(size_t)ka*nu + j];
                int32_t ib = GI[(size_t)kb*nu + j];
                int32_t ma = GM[(size_t)ka*nu + j];
                int32_t mb = GM[(size_t)kb*nu + j];
                for (int d = 0; d < 3; d++) {
                    double delta =
                        (double)G[((size_t)ka*nu + j)*3 + (size_t)d] -
                        (double)G[((size_t)kb*nu + j)*3 + (size_t)d];
                    d2 += delta * delta;
                }
                /* A V continuation is either one convex interpolation between
                 * local anchors or it is not made at all. */
                if (d2 > RIB_WRAP_GATE * RIB_WRAP_GATE * steps * steps) {
                    ka = kb;
                    continue;
                }
                for (long q = ka + 1; q < kb; q++) {
                    double w = (double)(q - ka) / (double)(kb - ka);
                    for (int d = 0; d < 3; d++) {
                        double ga =
                            (double)G[((size_t)ka*nu + j)*3 + (size_t)d];
                        double gb =
                            (double)G[((size_t)kb*nu + j)*3 + (size_t)d];
                        G[((size_t)q*nu + j)*3 + (size_t)d] =
                            (float)((1.0-w) * ga + w * gb);
                    }
                    GP[(size_t)q*nu + j] =
                        (float)((1.0-w) * GP[(size_t)ka*nu + j] +
                                     w  * GP[(size_t)kb*nu + j]);
                    GO[(size_t)q*nu + j] = 1;
                    if (ia >= 0 && ia == ib) GI[(size_t)q*nu + j] = ia;
                    if (ma >= 0 && ma == mb) GM[(size_t)q*nu + j] = ma;
                }
            }
            ka = kb;
        }
    }
    }
    out->grid_long_edges_vfill = rib_grid_long_edges(G, GO, nk, nu);
    fprintf(stderr, "    [fit] v-fill: %.2fs\n", ves_clock_sec() - fit_t0);
    fit_t0 = ves_clock_sec();

    /* This baseline does not infer topology from a local distance threshold.
     * Explicit occupancy is the parameterization domain, every measured point
     * remains present, and the writer emits every occupied four-corner cell.
     * Long edges below are diagnostics only. */
    out->grid_long_edges_local = rib_grid_long_edges(G,GO,nk,nu);
    fit_t0 = ves_clock_sec();
    if (o->verify_fit_width)
        rib_grid_digest(o->fit_cover_width ? "cover" : "compact", "v-fill",
                        G, GO, nk, nu, active_nu < nu ? active_nu : nu);
    /* No post-fit XYZ smoothing: measured cross-sections remain exact. */
    out->grid_long_edges_smooth = rib_grid_long_edges(G, GO, nk, nu);
    out->grid_both_diagonals_long =
        rib_grid_both_diagonals_long(G, GO, nk, nu);
    {
        size_t long_h=0, long_v=0;
        rib_grid_long_edge_directions(G,GO,nk,nu,&long_h,&long_v);
        fprintf(stderr,
                "  fitted-grid diagnostics: long edges row=%zu repaired=%zu "
                "vfill=%zu diagnostic=%zu final=%zu (h=%zu v=%zu) "
                "both-diagonals=%zu claim-updates=%zu "
                "local-repairs=%zu+%zu-pair/%zu-slot+%zu-run/%zu-slot+"
                "%zu-patch/%zu-slot "
                "(patch-small=%zu enclosed=%zu claims=%zu) "
                "(%zu supported) "
                "row-splits=%zu "
                "(same-mesh=%zu same-solve=%zu same-island=%zu)\n",
                out->grid_long_edges_rowfit,out->grid_long_edges_repaired,
                out->grid_long_edges_vfill,out->grid_long_edges_local,
                out->grid_long_edges_smooth,
                long_h,long_v,out->grid_both_diagonals_long,
                out->grid_claim_repairs,out->grid_local_outlier_repairs,
                out->grid_local_pair_repairs,out->grid_local_pair_slots,
                out->grid_local_run_repairs,out->grid_local_run_slots,
                out->grid_local_patch_repairs,out->grid_local_patch_slots,
                out->grid_local_patch_small_components,
                out->grid_local_patch_enclosed_components,
                out->grid_local_patch_claim_slots,
                out->grid_local_supported_replacements,
                out->grid_row_wrap_splits,
                out->grid_row_wrap_same_mesh,out->grid_row_wrap_same_solve,
                out->grid_row_wrap_same_island);
    }
    if (o->verify_fit_width)
        rib_grid_digest(o->fit_cover_width ? "cover" : "compact", "smooth",
                        G, GO, nk, nu, active_nu < nu ? active_nu : nu);

    /* Scratch built after `mark` is dead here.  Restore BEFORE exporting the
     * chart runs: the provenance arrays are handed to the caller and must
     * outlive this call.  (Exporting first and restoring after returned
     * reclaimed arena memory -- a measured use-after-restore.)  The export
     * reads only S (pre-mark) and the heap-allocated labels. */
    Arena_restore(arena, mark);
    if (o->component_global && o->ownership_construction &&
        (rib_claim_export_chart_runs(
             arena, S, claim_labels, nk, &out->grid_chart_row_offsets,
             &out->grid_chart_runs, &out->grid_chart_run_count) != 0 ||
         rib_fuse_audit_runs(o, out, nk) != 0)) {
        fprintf(stderr,
                "ERROR: could not freeze fitted-grid claimant provenance\n");
        rib_claim_label_rows_dispose(claim_labels, nk);
        claim_labels = NULL;
        out->grid_chart_row_offsets = NULL;
        out->grid_chart_runs = NULL;
        out->grid_chart_run_count = 0;
        return;
    }
    rib_claim_label_rows_dispose(claim_labels, nk);
    claim_labels = NULL;
    out->grid_pos   = G;
    out->grid_present = GO;
    out->grid_valid = GV;
    out->grid_phi = GP;
    out->grid_island = GI;
    out->grid_material = GM;
    out->nu = nu;
    out->nk = nk;
    out->grid_du = (float)du;
    out->grid_dv = o->slice_h;
}

/* Turn one completed robust solve round into the same fitted quad-ribbon used
 * by the production endpoint.  Ribbon_write_obj publishes the VMESH through a
 * temporary file + atomic rename; the live RAW baker can safely start the
 * instant the final pathname appears. */
static int rib_dump_solve_round_checkpoint(
    Arena_T arena, SliceSet *S, const PairSet *P, const RibbonOpts *o,
    const RibbonResult *base, int round, int total_rounds,
    const char *phase)
{
    if (o->solve_level_prefix == NULL) return 0;
    Arena_Mark checkpoint_mark = Arena_save(arena);
    double *saved_u = RIB_ALLOC_ARRAY(arena, double, S->n_smp);
    int32_t *saved_component = RIB_ALLOC_ARRAY(arena, int32_t, S->n_chn);
    double *saved_atlas_shift = RIB_ALLOC_ARRAY(arena, double, S->n_chn);
    ChainRelationGraph relations;
    RibbonResult level;
    RibbonWriteStats ws;
    char locator[2600], vmesh[2600];
    int rc = -1;
    memset(&relations, 0, sizeof relations);
    memset(&level, 0, sizeof level);
    memset(&ws, 0, sizeof ws);
    level.source_chart_island = base->source_chart_island;
    level.source_chart_island_count = base->source_chart_island_count;
    for (size_t i = 0; i < S->n_smp; i++) saved_u[i] = S->smp[i].u;
    for (size_t c = 0; c < S->n_chn; c++) {
        saved_component[c] = S->chn[c].reconstruction_component;
        saved_atlas_shift[c] = S->chn[c].atlas_shift;
    }

    if (chain_relation_graph_from_pairs(S, P, &relations) != 0 ||
        rib_assign_reconstruction_components(
            arena, S, o, &relations, &level) != 0)
        goto done;

    /* A solve component has an unobservable additive U gauge.  Direct-ribbon
     * checkpoints used to rasterize those components on top of one another,
     * turning a parameterization diagnostic into a claimant-collision test and
     * deleting most of the visible ribbon.  Pack by additive component shifts
     * before fitting; this preserves every intrinsic derivative exactly. */
    if (o->component_global && !o->emit_global)
        pack_metric_islands(arena, S, (double)o->grid_u,
                            NULL, NULL, NULL, NULL, NULL, NULL, 0,
                            NULL, &level);

    double U0 = 1e300;
    for (size_t i = 0; i < S->n_smp; i++)
        if (S->smp[i].chain >= 0 && S->smp[i].u < U0) U0 = S->smp[i].u;
    if (U0 >= 1e299) goto done;

    RibbonOpts round_opts = *o;
    round_opts.solve_level_prefix = NULL;
    round_opts.solve_level_write_obj = 0;
    fit_ribbon(arena, S, U0, &round_opts, &relations, &level);
    if (level.grid_pos == NULL) goto done;
    if (snprintf(locator, sizeof locator, "%s_solve_round_%02d.obj",
                 o->solve_level_prefix, round) < 0 ||
        snprintf(vmesh, sizeof vmesh, "%s_solve_round_%02d.vmesh",
                 o->solve_level_prefix, round) < 0 ||
        strlen(locator) + 1 >= sizeof locator ||
        strlen(vmesh) + 1 >= sizeof vmesh)
        goto done;
    ves_ensure_parent_dir(locator);
    if (Ribbon_write_obj(locator, &level, o->solve_level_write_obj,
                         0, &ws) != 0)
        goto done;
    fprintf(stderr,
            "  solve-round checkpoint %d/%d (%s): ribbon=%zux%zu "
            "emitted=%zu/%zu vmesh=%s\n",
            round + 1, total_rounds, phase, level.nu, level.nk,
            ws.vertices, ws.faces, vmesh);
    rc = 0;

done:
    for (size_t i = 0; i < S->n_smp; i++) S->smp[i].u = saved_u[i];
    for (size_t c = 0; c < S->n_chn; c++) {
        S->chn[c].reconstruction_component = saved_component[c];
        S->chn[c].atlas_shift = saved_atlas_shift[c];
    }
    chain_relation_graph_dispose(&relations);
    Arena_restore(arena, checkpoint_mark);
    return rc;
}

/* Materialize the final robust solve's FMG ascent without replacing the modern
 * ribbon geometry by source-mesh triangles.  A monotone copy of each level is
 * passed through the ordinary Stage E fitter and direct-grid writer.  The full-
 * approximation prolongation retains the known fine arclength while exposing
 * only the placement correction resolved at that hierarchy level.  The VMESH
 * is therefore the same kind of quad-ribbon surface as the production endpoint
 * and can be baked by sampling RAW at its barycentrically interpolated XYZ
 * positions. */
static int rib_gmg_dump_levels(Arena_T arena, SliceSet *S,
                               const RibbonOpts *o, RibbonResult *base,
                               RibGMGSnapshots *C)
{
    if (C == NULL || C->count == 0 || o->solve_level_prefix == NULL) return 0;
    if (C->n_sample != S->n_smp) return -1;
    size_t n = S->n_smp;
    double *final_u = (double *)malloc(rib_array_bytes(n, sizeof(*final_u)));
    if (final_u == NULL) return -1;
    for (size_t i = 0; i < n; i++) final_u[i] = S->smp[i].u;

    char manifest_path[2600];
    if (snprintf(manifest_path, sizeof manifest_path, "%s_gmg_levels.json",
                 o->solve_level_prefix) < 0 ||
        strlen(manifest_path) + 1 >= sizeof manifest_path) {
        free(final_u); return -1;
    }
    ves_ensure_parent_dir(manifest_path);
    FILE *manifest = fopen(manifest_path, "w");
    if (manifest == NULL) { free(final_u); return -1; }
    fprintf(manifest,
            "{\n  \"schema\": \"ribbon-gmg-levels-v1\",\n"
            "  \"order\": \"coarse-to-fine\",\n"
            "  \"prolongation\": \"fine-reference-plus-coarse-defect-correction\",\n"
            "  \"sample_count\": %zu,\n  \"levels\": [\n", n);

    int rc = 0;
    for (int ordinal = 0; ordinal < C->count; ordinal++) {
        RibGMGSnapshot *D = &C->snapshot[ordinal];
        char locator[2600], vmesh[2600];
        snprintf(locator, sizeof locator, "%s_gmg_level_%02d.obj",
                 o->solve_level_prefix, ordinal);
        if (snprintf(vmesh, sizeof vmesh, "%s_gmg_level_%02d.vmesh",
                     o->solve_level_prefix, ordinal) < 0) {
            rc = -1; break;
        }
        ves_ensure_parent_dir(locator);
        Arena_Mark level_mark = Arena_save(arena);
        for (size_t i = 0; i < n; i++) S->smp[i].u = (double)D->u[i];
        int32_t maxchain = 1;
        for (size_t c = 0; c < S->n_chn; c++)
            if (S->chn[c].count > maxchain) maxchain = S->chn[c].count;
        double *lb = RIB_ALLOC_ARRAY(arena, double, (size_t)maxchain);
        double *scratch = RIB_ALLOC_ARRAY(arena, double, (size_t)maxchain * 3u);
        double *chain_u = RIB_ALLOC_ARRAY(arena, double, (size_t)maxchain);
        for (size_t c = 0; c < S->n_chn; c++) {
            int32_t first = S->chn[c].first, count = S->chn[c].count;
            if (count < 2) continue;
            lb[0] = 0.0;
            chain_u[0] = S->smp[first].u;
            for (int32_t i = 1; i < count; i++) {
                lb[i] = 0.5 * (S->smp[first+i].s - S->smp[first+i-1].s);
                chain_u[i] = S->smp[first+i].u;
            }
            (void)pava_chain(chain_u, lb, count, scratch);
            for (int32_t i = 0; i < count; i++)
                S->smp[first+i].u = chain_u[i];
        }
        double U0 = 1e300;
        for (size_t i = 0; i < n; i++)
            if (S->smp[i].chain >= 0 && S->smp[i].u < U0) U0 = S->smp[i].u;

        RibbonOpts level_opts = *o;
        level_opts.solve_level_prefix = NULL;
        level_opts.solve_level_write_obj = 0;
        RibbonResult level;
        memset(&level, 0, sizeof level);
        level.source_chart_island = base->source_chart_island;
        level.source_chart_island_count = base->source_chart_island_count;
        fit_ribbon(arena, S, U0, &level_opts, NULL, &level);
        if (level.grid_pos == NULL) {
            fprintf(stderr,
                    "  GMG checkpoint %d/%d: stored=%d has no fitted ribbon\n",
                    ordinal + 1, C->count, D->stored_level);
            fprintf(manifest,
                    "    { \"ordinal\": %d, \"stored_level\": %d, "
                    "\"nodes\": %d, \"h_u\": %.17g, \"h_v\": %.17g, "
                    "\"emitted\": false, \"vertices\": 0, \"faces\": 0 }%s\n",
                    ordinal, D->stored_level, D->nodes, D->h_u, D->h_v,
                    ordinal + 1 < C->count ? "," : "");
            Arena_restore(arena, level_mark);
            continue;
        }
        RibbonWriteStats ws;
        memset(&ws, 0, sizeof ws);
        if (Ribbon_write_obj(locator, &level, o->solve_level_write_obj,
                             0, &ws) != 0) {
            fprintf(stderr,
                    "  GMG checkpoint %d/%d: stored=%d has no direct-grid faces\n",
                    ordinal + 1, C->count, D->stored_level);
            fprintf(manifest,
                    "    { \"ordinal\": %d, \"stored_level\": %d, "
                    "\"nodes\": %d, \"h_u\": %.17g, \"h_v\": %.17g, "
                    "\"emitted\": false, \"vertices\": 0, \"faces\": 0 }%s\n",
                    ordinal, D->stored_level, D->nodes, D->h_u, D->h_v,
                    ordinal + 1 < C->count ? "," : "");
            Arena_restore(arena, level_mark);
            continue;
        }
        fprintf(stderr,
                "  GMG checkpoint %d/%d: stored=%d nodes=%d h=(%.3g,%.3g) "
                "ribbon=%zux%zu emitted=%zu/%zu vmesh=%s\n",
                ordinal + 1, C->count, D->stored_level, D->nodes,
                D->h_u, D->h_v, level.nu, level.nk,
                ws.vertices, ws.faces, vmesh);
        fprintf(manifest,
                "    { \"ordinal\": %d, \"stored_level\": %d, "
                "\"nodes\": %d, \"h_u\": %.17g, \"h_v\": %.17g, "
                 "\"locator\": \"%s_gmg_level_%02d.obj\", "
                 "\"vmesh\": \"%s_gmg_level_%02d.vmesh\", "
                 "\"emitted\": true, \"vertices\": %zu, \"faces\": %zu }%s\n",
                ordinal, D->stored_level, D->nodes, D->h_u, D->h_v,
                ves_path_basename(o->solve_level_prefix), ordinal,
                ves_path_basename(o->solve_level_prefix), ordinal,
                ws.vertices, ws.faces,
                ordinal + 1 < C->count ? "," : "");
        Arena_restore(arena, level_mark);
    }
    fprintf(manifest, "  ]\n}\n");
    if (fclose(manifest) != 0) rc = -1;
    for (size_t i = 0; i < n; i++) S->smp[i].u = final_u[i];
    free(final_u);
    return rc;
}

/* The reference-U path already knows every row's metric exactly: after the
 * gauge+s projection, the only remaining freedom is one additive gauge per
 * slice chain.  Solving the original sample graph would therefore spend most
 * of its time rediscovering those fixed within-row differences.  Collapse it
 * to a chain graph instead.  A sample correspondence
 *
 *     u[a] - u[b] = d
 *
 * becomes g[ca] - g[cb] = d - s[a] + s[b].  In carried-registration mode the
 * objective is
 *
 *   sum_p w_p rho(g[a]-g[b]-d_p)
 *       + sum_c n_c (g[c]-g_graph[c])^2,
 *
 * where n_c is the number of graph observations on the slice chain.  The
 * numerical unknown is delta[c] = g[c]-g_graph[c], not absolute g: this removes
 * the O(1e5) universal-cover offset from the CG stopping norm while leaving the
 * objective exactly unchanged.  Thus no connected component has a free
 * additive mode: physical correspondences can reconcile adjacent rows, but
 * cannot replace the globally synchronized graph frame.  Initial robust weights
 * are measured against that frozen frame before the first solve.  No raster or
 * quilting signal enters this objective. */
static int solve_reference_chain_gauges(Arena_T arena, SliceSet *S,
                                         const PairSet *P,
                                         const RibbonOpts *o,
                                         const float *sample_gauge_conf,
                                         ChainRelationGraph *relations,
                                         RibbonResult *out)
{
    size_t nc = S->n_chn;
    if (nc == 0 || nc > (size_t)INT32_MAX) return -1;

    Arena_Mark mark = Arena_save(arena);
    double *g0 = RIB_ALLOC_ARRAY(arena, double, nc);
    double *g  = RIB_ALLOC_ARRAY(arena, double, nc);
    double *chain_conf = NULL;
    if (sample_gauge_conf != NULL) {
        chain_conf = RIB_ALLOC_ARRAY(arena, double, nc);
        for (size_t c = 0; c < nc; c++) {
            int32_t f = S->chn[c].first, n = S->chn[c].count;
            double sum = 0.0;
            for (int32_t i = 0; i < n; i++) {
                double q = (double)sample_gauge_conf[f + i];
                if (!isfinite(q)) q = 0.0;
                if (q < 0.0) q = 0.0;
                if (q > 1.0) q = 1.0;
                sum += q;
            }
            chain_conf[c] = n > 0 ? sum / (double)n : 0.0;
        }
    }
    for (size_t c = 0; c < nc; c++) {
        int32_t f = S->chn[c].first;
        g0[c] = S->smp[f].u - S->chn[c].atlas_shift - S->smp[f].s;
        g[c] = 0.0; /* solve the correction to the graph gauge */
    }

    PairSet Q;
    memset(&Q, 0, sizeof Q);
    Q.pairs = RIB_ALLOC_ARRAY(arena, Pair, P->n_pairs + 1);
    UnionFind uf = UF_new(arena, (int32_t)nc);
    size_t phase_rejected = 0, island_rejected = 0, same_chain = 0;
    size_t trusted_topology = 0;
    for (size_t p = 0; p < P->n_pairs; p++) {
        const Pair *src = &P->pairs[p];
        if (src->a < 0 || src->b < 0 ||
            (size_t)src->a >= S->n_smp || (size_t)src->b >= S->n_smp)
            continue;
        const Sample *sa = &S->smp[src->a];
        const Sample *sb = &S->smp[src->b];
        int32_t ca = sa->chain, cb = sb->chain;
        if (ca < 0 || cb < 0 || (size_t)ca >= nc || (size_t)cb >= nc)
            continue;
        if (ca == cb) { same_chain++; continue; }
        if (S->chn[ca].winding_island != S->chn[cb].winding_island) {
            island_rejected++;
            continue;
        }
        if (o->reference_phi != NULL &&
            fabs(sa->phi - sb->phi) > RIB_CONT_DPHI_MAX) {
            phase_rejected++;
            continue;
        }
        size_t qi = Q.n_pairs++;
        Pair *dst = &Q.pairs[qi];
        dst->a = ca;
        dst->b = cb;
        dst->k = src->k;
        /* Convert the absolute physical gauge equation into a correction
         * equation.  At delta=0, -dst->d is exactly the residual of the frozen
         * graph gauge against this geometric correspondence. */
        dst->d = src->d - sa->s + sb->s - g0[ca] + g0[cb];
        dst->w = src->w > 0.0 ? src->w : RIB_ALIGN_W;
        dst->like = 1.0;
        if (dst->k == RIB_PAIR_K_TRUSTED_TOPOLOGY) trusted_topology++;
        uf_union(&uf, ca, cb);
    }
    if (Q.n_pairs == 0) {
        for (size_t c = 0; c < nc; c++) {
            S->chn[c].solve_comp = (int32_t)c;
            S->chn[c].group = (int32_t)c;
        }
        out->n_qp_comps = (int)nc;
        Arena_restore(arena, mark);
        return 0;
    }

    /* One absolute pin per connected chain-graph component.  Pick its longest
     * physical row, mirroring the full solve's least-fragile gauge choice. */
    int32_t *root_best = RIB_ALLOC_ARRAY(arena, int32_t, nc);
    int32_t *root_label = RIB_ALLOC_ARRAY(arena, int32_t, nc);
    for (size_t c = 0; c < nc; c++) {
        root_best[c] = -1;
        root_label[c] = -1;
    }
    for (size_t c = 0; c < nc; c++) {
        int32_t r = uf_find(&uf, (int32_t)c);
        if (root_best[r] < 0 ||
            S->chn[c].count > S->chn[root_best[r]].count)
            root_best[r] = (int32_t)c;
    }
    int32_t nsolve = 0;
    for (size_t c = 0; c < nc; c++) {
        if (uf_find(&uf, (int32_t)c) != (int32_t)c) continue;
        root_label[c] = nsolve++;
    }
    for (size_t c = 0; c < nc; c++) {
        int32_t root = uf_find(&uf, (int32_t)c);
        S->chn[c].solve_comp = root_label[root];
        S->chn[c].group = root_label[root];
    }
    out->n_qp_comps = (int)nsolve;

    int32_t *pin_idx = RIB_ALLOC_ARRAY(arena, int32_t, nc);
    double *pin_val = RIB_ALLOC_ARRAY(arena, double, nc);
    double *pin_weight = NULL;
    size_t n_pins = 0;
    double anchor_mult_sum = 0.0, anchor_obs_sum = 0.0;
    double anchor_mult_min = INFINITY, anchor_mult_max = 0.0;
    if (o->reference_anchor_gauges) {
        /* One actual graph observation per slice sample.  Using its count as
         * the anchor weight makes both sides of the objective data-counted:
         * a long, well-observed chain cannot be moved by a handful of pairs. */
        pin_weight = RIB_ALLOC_ARRAY(arena, double, nc);
        for (size_t c = 0; c < nc; c++) {
            pin_idx[c] = (int32_t)c;
            pin_val[c] = 0.0;
            double obs = S->chn[c].count > 0
                       ? (double)S->chn[c].count : 1.0;
            double mult = 1.0;
            if (chain_conf != NULL) {
                double q2 = chain_conf[c] * chain_conf[c];
                mult += RIB_METRIC_CERT_ANCHOR_ODDS * q2 * q2;
            }
            pin_weight[c] = obs * mult;
            anchor_mult_sum += obs * mult;
            anchor_obs_sum += obs;
            if (mult < anchor_mult_min) anchor_mult_min = mult;
            if (mult > anchor_mult_max) anchor_mult_max = mult;
        }
        n_pins = nc;
    } else {
        /* Legacy non-registered path: one absolute pin per connected chain
         * component, placed on its longest physical row. */
        for (size_t c = 0; c < nc; c++) {
            if (uf_find(&uf, (int32_t)c) != (int32_t)c) continue;
            int32_t pc = root_best[c];
            if (pc < 0) continue;
            pin_idx[n_pins] = pc;
            pin_val[n_pins] = 0.0;
            n_pins++;
        }
    }

    /* Robust IRLS over the chain-pair graph; each round is an exact TAUCS
     * solve of the likelihood-weighted normal equations (the old matrix-free
     * PCG detour through a fake length-1-chain SliceSet is retired). */
    size_t gq_max_nt = 3u * Q.n_pairs + n_pins;
    if (nc > (size_t)INT32_MAX || gq_max_nt > (size_t)INT32_MAX) {
        Arena_restore(arena, mark);
        return -1;
    }
    int *gq_row = RIB_ALLOC_ARRAY(arena, int, gq_max_nt + 1);
    int *gq_col = RIB_ALLOC_ARRAY(arena, int, gq_max_nt + 1);
    double *gq_val = RIB_ALLOC_ARRAY(arena, double, gq_max_nt + 1);
    double *gq_rhs = RIB_ALLOC_ARRAY(arena, double, nc + 1);
    double *gq_prev = RIB_ALLOC_ARRAY(arena, double, nc + 1);
    size_t stat_n = Q.n_pairs > nc ? Q.n_pairs : nc;
    double *resbuf = RIB_ALLOC_ARRAY(arena, double, stat_n + 1);

    /* Gate the very first linear solve as well as later IRLS rounds.  Starting
     * every pair at likelihood 1 gave one bad edge a full-strength opportunity
     * to pull a disconnected component before robustness was ever measured. */
    for (size_t p = 0; p < Q.n_pairs; p++)
        resbuf[p] = fabs(Q.pairs[p].d);
    double before_med = select_median_dbl(resbuf, Q.n_pairs);
    for (size_t p = 0; p < Q.n_pairs; p++)
        resbuf[p] = fabs(Q.pairs[p].d);
    double before_p95 = select_p95_dbl(resbuf, Q.n_pairs);
    double initial_sigma = 3.0 * before_med;
    if (initial_sigma < 2.0) initial_sigma = 2.0;
    {
        double inv2s2 = 1.0 / (2.0 * initial_sigma * initial_sigma);
    for (size_t p = 0; p < Q.n_pairs; p++) {
        if (Q.pairs[p].k == RIB_PAIR_K_TRUSTED_TOPOLOGY) {
            Q.pairs[p].like = rib_metric_topology_likelihood(
                -Q.pairs[p].d);
            continue;
        }
        double r = -Q.pairs[p].d;
            double like = exp(-r * r * inv2s2);
            Q.pairs[p].like = like > 1e-6 ? like : 1e-6;
        }
    }
    out->chain_gauge_pair_median_before = before_med;
    out->chain_gauge_pair_p95_before = before_p95;
    int total_rounds = (o->relax_iters > 0 ? o->relax_iters : 1)
                     + (o->final_iters > 0 ? o->final_iters : 0);
    if (total_rounds > 7) total_rounds = 7;
    double final_med = 0.0;
    for (int round = 0; round < total_rounds; round++) {
        size_t nt = 0;
        memset(gq_rhs, 0, nc * sizeof(*gq_rhs));
        for (size_t p = 0; p < Q.n_pairs; p++) {
            const Pair *pr = &Q.pairs[p];
            if (pr->like < 0.0) continue;
            double w = pr->w * pr->like;
            if (w <= 0.0) continue;
            int a = pr->a, b = pr->b;
            int lo2 = a < b ? a : b, hi2 = a < b ? b : a;
            gq_row[nt] = a; gq_col[nt] = a; gq_val[nt] = w; nt++;
            gq_row[nt] = b; gq_col[nt] = b; gq_val[nt] = w; nt++;
            gq_row[nt] = hi2; gq_col[nt] = lo2; gq_val[nt] = -w; nt++;
            gq_rhs[a] += w * pr->d;
            gq_rhs[b] -= w * pr->d;
        }
        for (size_t pi = 0; pi < n_pins; pi++) {
            double wp = pin_weight != NULL ? pin_weight[pi] : 1.0;
            gq_row[nt] = pin_idx[pi];
            gq_col[nt] = pin_idx[pi];
            gq_val[nt] = wp;
            gq_rhs[pin_idx[pi]] += wp * pin_val[pi];
            nt++;
        }
        memcpy(gq_prev, g, nc * sizeof(*gq_prev));
        if (Sparse_solve_sym((int)nc, (int)nt, gq_row, gq_col, gq_val,
                             gq_rhs, g, SPARSE_SPD) != 0) {
            fprintf(stderr, "  fast chain-gauge: TAUCS solve FAILED "
                    "(vars=%zu nnz=%zu)\n", nc, nt);
            Arena_restore(arena, mark);
            return -1;
        }
        double gchg = 0.0;
        for (size_t c = 0; c < nc; c++)
            if (fabs(g[c] - gq_prev[c]) > gchg)
                gchg = fabs(g[c] - gq_prev[c]);
        for (size_t p = 0; p < Q.n_pairs; p++)
            resbuf[p] = fabs(g[Q.pairs[p].a] - g[Q.pairs[p].b]
                           - Q.pairs[p].d);
        final_med = select_median_dbl(resbuf, Q.n_pairs);
        double sigma = (round < o->relax_iters ? 3.0 : 1.5) * final_med;
        if (sigma < 2.0) sigma = 2.0;
        double inv2s2 = 1.0 / (2.0 * sigma * sigma);
        for (size_t p = 0; p < Q.n_pairs; p++) {
            if (Q.pairs[p].k == RIB_PAIR_K_TRUSTED_TOPOLOGY) {
                double r = g[Q.pairs[p].a] - g[Q.pairs[p].b]
                         - Q.pairs[p].d;
                Q.pairs[p].like = rib_metric_topology_likelihood(r);
                continue;
            }
            double r = g[Q.pairs[p].a] - g[Q.pairs[p].b] - Q.pairs[p].d;
            double like = exp(-r * r * inv2s2);
            Q.pairs[p].like = like > 1e-6 ? like : 1e-6;
        }
        fprintf(stderr,
                "  fast chain-gauge round %d/%d: vars=%zu pairs=%zu "
                "exact solve, dmax=%.3g median=%.4f sigma=%.4f\n",
                round + 1, total_rounds, nc, Q.n_pairs,
                gchg, final_med, sigma);
        if (round > 0 && gchg < 1e-6) break;
    }

    for (size_t p = 0; p < Q.n_pairs; p++)
        resbuf[p] = fabs(g[Q.pairs[p].a] - g[Q.pairs[p].b]
                       - Q.pairs[p].d);
    double after_med = select_median_dbl(resbuf, Q.n_pairs);
    for (size_t p = 0; p < Q.n_pairs; p++)
        resbuf[p] = fabs(g[Q.pairs[p].a] - g[Q.pairs[p].b]
                       - Q.pairs[p].d);
    double after_p95 = select_p95_dbl(resbuf, Q.n_pairs);
    out->chain_gauge_pair_median_after = after_med;
    out->chain_gauge_pair_p95_after = after_p95;

    if (relations != NULL) {
        size_t retained = 0;
        for (size_t p = 0; p < Q.n_pairs; p++) {
            if (Q.pairs[p].like < 0.5) continue;
            if (chain_relation_graph_add(
                    relations, Q.pairs[p].a, Q.pairs[p].b,
                    Q.pairs[p].like) != 0) {
                Arena_restore(arena, mark);
                return -1;
            }
            retained++;
        }
        fprintf(stderr,
                "  fast chain lineage graph: %zu relation(s), "
                "%zu robust observation(s)\n",
                relations->n_edges, retained);
    }

    double shift2 = 0.0, shift_weight = 0.0, shiftmax = 0.0;
    for (size_t c = 0; c < nc; c++) {
        double shift = g[c];
        double obs = S->chn[c].count > 0 ? (double)S->chn[c].count : 1.0;
        shift2 += obs * shift * shift;
        shift_weight += obs;
        resbuf[c] = fabs(shift);
        if (fabs(shift) > shiftmax) shiftmax = fabs(shift);
        int32_t f = S->chn[c].first, n = S->chn[c].count;
        for (int32_t i = 0; i < n; i++)
            S->smp[f+i].u = g0[c] + g[c] + S->smp[f+i].s
                            + S->chn[c].atlas_shift;
    }
    out->chain_gauge_shift_rms = shift_weight > 0.0
                               ? sqrt(shift2 / shift_weight) : 0.0;
    out->chain_gauge_shift_p95 = select_p95_dbl(resbuf, nc);
    out->chain_gauge_shift_max = shiftmax;
    fprintf(stderr,
            "  fast chain-gauge solve: %s components=%d trusted-topology=%zu "
            "phase-rejected=%zu "
            "same-chain=%zu island-rejected=%zu/%zu "
            "cert-anchor=%.1f[%.1f..%.1f] "
            "pair-med/p95=%.4f/%.4f->%.4f/%.4f "
            "shift-rms/p95/max=%.3f/%.3f/%.3f\n",
            o->reference_anchor_gauges ? "ANCHORED" : "legacy",
            nsolve, trusted_topology, phase_rejected, same_chain, island_rejected,
            P->relation_island_rejects,
            anchor_obs_sum > 0.0 ? anchor_mult_sum / anchor_obs_sum : 1.0,
            isfinite(anchor_mult_min) ? anchor_mult_min : 1.0,
            anchor_mult_max > 0.0 ? anchor_mult_max : 1.0,
            before_med, before_p95, after_med, after_p95,
            out->chain_gauge_shift_rms, out->chain_gauge_shift_p95, shiftmax);
    Arena_restore(arena, mark);
    return 0;
}

/* Metric-projection island placement.  The carried U of a detached geometry
 * component is an atlas-packing artifact (the upstream fit packs disconnected
 * islands side by side), but its carried lifted PHASE is physical: it comes
 * from the contact-resolved winding registration.  Build the dominant
 * component's monotone phi->u map from the projected samples and rigidly
 * place every other component at u = F(phi) with ONE median offset per
 * component.  Direct measurement only -- no solve, no ordering heuristic, no
 * packing; each component's internal metric is untouched.  A component with
 * no finite phase keeps its carried placement and is reported. */
static void rib_project_register_islands(Arena_T arena, SliceSet *S,
                                         double r_ref, double place_pitch,
                                         int solved_u_round,
                                         RibbonResult *out)
{
    Arena_Mark mark = Arena_save(arena);
    if (!(place_pitch > 0.0)) place_pitch = 9.5;
    size_t n = S->n_smp;
    int32_t ncomp = 0;
    for (size_t c = 0; c < S->n_chn; c++)
        if (S->chn[c].mesh_comp >= ncomp) ncomp = S->chn[c].mesh_comp + 1;
    if (ncomp <= 1) {
        fprintf(stderr,
                "  island placement: single geometry component; carried "
                "frame already physical\n");
        Arena_restore(arena, mark);
        return;
    }
    size_t *count = (size_t *)ARENA_CALLOC(arena, (size_t)ncomp,
                                           sizeof(*count));
    size_t nfinite = 0;
    for (size_t i = 0; i < n; i++) {
        int32_t ch = S->smp[i].chain;
        if (ch < 0 || !isfinite(S->smp[i].phi) || !isfinite(S->smp[i].u))
            continue;
        count[S->chn[ch].mesh_comp]++;
        nfinite++;
    }
    int32_t main_comp = 0;
    for (int32_t c = 1; c < ncomp; c++)
        if (count[c] > count[main_comp]) main_comp = c;
    if (count[main_comp] < 2 || nfinite < 2) {
        fprintf(stderr,
                "  island placement: dominant component has no usable "
                "phase; carried frame kept\n");
        Arena_restore(arena, mark);
        return;
    }

    /* The BULK PHASE WINDOW: detached crumbs whose gauge island never
     * connected to the bulk carry lifted phases hundreds of turns away;
     * phi-placing them extrapolates absurd u domains (10x rung: 3.35M vox
     * against ~500k physical while the top gauge island holds 95.5% of
     * samples).  The central-quantile window of ALL finite phases isolates
     * the consistent bulk: islands outside it keep their carried (atlas)
     * frame -- a bounded, honest position -- instead of a fabricated one. */
    double win_lo = 0.0, win_hi = 0.0;
    {
        double *phis = RIB_ALLOC_ARRAY(arena, double, nfinite + 1);
        size_t np = 0;
        for (size_t i = 0; i < n; i++) {
            int32_t ch = S->smp[i].chain;
            if (ch < 0 || !isfinite(S->smp[i].phi) ||
                !isfinite(S->smp[i].u))
                continue;
            phis[np++] = S->smp[i].phi;
        }
        qsort(phis, np, sizeof(double), cmp_dbl);
        {
            /* The window is PHYSICALLY CAPPED: a scroll of outer radius
             * ~2*r_ref at this pitch holds at most ~2*r_ref/pitch turns, so
             * no honest lift spans more.  Take the densest fixed-width
             * stretch of 1.5x that span (two-pointer over the sorted
             * phases): quantile windows cannot survive contaminated tails
             * (10x: 4.5% garbage kept them at 1,643 turns) and even
             * shortest-X% fails when the lift itself DRIFTS (10x: the 95.5%
             * gauge island spans ~944 turns against ~128 physical).  The
             * densest physical-width stretch is the least-drifted coherent
             * mass; everything beyond it keeps the carried frame. */
            double w_cap = 1.5 * (2.0 * M_PI) * (2.0 * r_ref) / place_pitch;
            size_t best_i = 0, best_cnt = 0, j = 0;
            for (size_t i = 0; i < np; i++) {
                if (j < i) j = i;
                while (j + 1 < np && phis[j + 1] - phis[i] <= w_cap) j++;
                if (j - i + 1 > best_cnt) {
                    best_cnt = j - i + 1;
                    best_i = i;
                }
            }
            win_lo = phis[best_i] - 8.0 * M_PI;
            win_hi = phis[best_i] + w_cap + 8.0 * M_PI;
            fprintf(stderr,
                    "  island placement: bulk phase window [%.1f, %.1f] rad "
                    "(densest %.0f-rad stretch holds %zu/%zu = %.1f%% of "
                    "finite-phase samples)\n",
                    win_lo, win_hi, w_cap, best_cnt, np,
                    100.0 * (double)best_cnt / (double)np);
        }
    }

    /* Map source: round 1 uses the dominant component (carried u of other
     * islands is atlas packing, not evidence); round 2 runs on SOLVED u,
     * where every in-window sample is consistency-tied and admissible --
     * this is what gives the map full multi-turn support when the dominant
     * MESH component covers only a sliver of the bulk gauge island. */
    size_t map_cap = solved_u_round ? nfinite : count[main_comp];
    PhiU2 *pu = RIB_ALLOC_ARRAY(arena, PhiU2, map_cap + 1);
    size_t nmain = 0;
    for (size_t i = 0; i < n; i++) {
        int32_t ch = S->smp[i].chain;
        if (ch < 0 || !isfinite(S->smp[i].phi) || !isfinite(S->smp[i].u))
            continue;
        if (!solved_u_round && S->chn[ch].mesh_comp != main_comp)
            continue;
        if (S->smp[i].phi < win_lo || S->smp[i].phi > win_hi)
            continue;
        pu[nmain].phi = S->smp[i].phi;
        pu[nmain].u = S->smp[i].u;
        nmain++;
    }
    if (nmain < 2) {
        fprintf(stderr,
                "  island placement: no in-window map samples; carried "
                "frame kept\n");
        Arena_restore(arena, mark);
        return;
    }
    qsort(pu, nmain, sizeof(PhiU2), cmp_phiu2);
    PhiUMap M;
    memset(&M, 0, sizeof M);
    M.end_slope = r_ref;
    {
        double bw = M_PI / 4.0;
        int nb_max = (int)((pu[nmain-1].phi - pu[0].phi) / bw) + 2;
        if (nb_max < 1) nb_max = 1;
        M.bphi = RIB_ALLOC_ARRAY(arena, double, (size_t)nb_max);
        M.bu = RIB_ALLOC_ARRAY(arena, double, (size_t)nb_max);
        double *tmpu = RIB_ALLOC_ARRAY(arena, double, nmain + 1);
        size_t i0 = 0;
        while (i0 < nmain) {
            double b0 = pu[0].phi +
                        floor((pu[i0].phi - pu[0].phi) / bw) * bw;
            size_t i1 = i0, m2 = 0;
            double psum = 0.0;
            while (i1 < nmain && pu[i1].phi < b0 + bw) {
                tmpu[m2++] = pu[i1].u;
                psum += pu[i1].phi;
                i1++;
            }
            if (m2 > 0 && M.nb < nb_max) {
                qsort(tmpu, m2, sizeof(double), cmp_dbl);
                M.bphi[M.nb] = psum / (double)m2;
                M.bu[M.nb] = tmpu[m2 / 2];
                M.nb++;
            }
            i0 = i1;
        }
        if (M.nb > 1) {
            double *map_lb = RIB_ALLOC_ARRAY(arena, double, (size_t)M.nb);
            double *map_scr = RIB_ALLOC_ARRAY(
                arena, double, (size_t)M.nb * 3u);
            memset(map_lb, 0, (size_t)M.nb * sizeof(double));
            (void)pava_chain(M.bu, map_lb, (int32_t)M.nb, map_scr);
        }
    }
    if (M.nb < 2) {
        fprintf(stderr,
                "  island placement: phi->u map underdetermined (%d bin); "
                "carried frame kept\n", M.nb);
        Arena_restore(arena, mark);
        return;
    }

    double *diffs = RIB_ALLOC_ARRAY(arena, double, n + 1);
    size_t moved = 0, unplaceable = 0, unanchored = 0;
    double max_shift = 0.0;
    /* Bucket samples by mesh component ONCE.  The former per-component full
     * scans were O(ncomp * n): invisible at the 4x5x5's 25 islands, but
     * hours at the 10x rung's tens of thousands (measured ~6 islands/s over
     * 24.7M samples). */
    size_t *comp_off = (size_t *)ARENA_CALLOC(
        arena, (size_t)ncomp + 2, sizeof(size_t));
    int32_t *comp_idx = RIB_ALLOC_ARRAY(arena, int32_t, n + 1);
    for (size_t i = 0; i < n; i++) {
        int32_t ch = S->smp[i].chain;
        int32_t mc = ch >= 0 ? S->chn[ch].mesh_comp : -1;
        if (mc >= 0 && mc < ncomp) comp_off[(size_t)mc + 1]++;
    }
    for (size_t c = 0; c < (size_t)ncomp; c++) comp_off[c + 1] += comp_off[c];
    {
        size_t *cur = RIB_ALLOC_ARRAY(arena, size_t, (size_t)ncomp + 1);
        memcpy(cur, comp_off, (size_t)ncomp * sizeof(size_t));
        for (size_t i = 0; i < n; i++) {
            int32_t ch = S->smp[i].chain;
            int32_t mc = ch >= 0 ? S->chn[ch].mesh_comp : -1;
            if (mc >= 0 && mc < ncomp)
                comp_idx[cur[(size_t)mc]++] = (int32_t)i;
        }
    }
    for (int32_t comp = 0; comp < ncomp; comp++) {
        if (comp == main_comp) continue;
        size_t nd = 0;
        double u_lo = 1e300, u_hi = -1e300;
        double ph_lo = 1e300, ph_hi = -1e300;
        size_t members = 0;
        for (size_t s = comp_off[comp]; s < comp_off[(size_t)comp + 1]; s++) {
            size_t i = (size_t)comp_idx[s];
            members++;
            if (!isfinite(S->smp[i].phi) || !isfinite(S->smp[i].u))
                continue;
            diffs[nd++] = phiu_lookup(&M, S->smp[i].phi) - S->smp[i].u;
            if (S->smp[i].u < u_lo) u_lo = S->smp[i].u;
            if (S->smp[i].u > u_hi) u_hi = S->smp[i].u;
            if (S->smp[i].phi < ph_lo) ph_lo = S->smp[i].phi;
            if (S->smp[i].phi > ph_hi) ph_hi = S->smp[i].phi;
        }
        if (members == 0) continue;
        if (nd == 0) {
            unplaceable++;
            if (unplaceable <= 50 || (unplaceable % 1000) == 0)
                fprintf(stderr,
                        "    island comp %d: %zu sample(s), no finite phase; "
                        "carried placement kept\n", comp, members);
            continue;
        }
        /* Islands whose lifted phase sits outside the bulk window carry an
         * UNANCHORED gauge (a detached crumb's private lift): any phi-based
         * placement of them is fabricated.  They keep the carried atlas
         * frame -- bounded and honestly labeled by provenance. */
        {
            double ph_mid = 0.5 * (ph_lo + ph_hi);
            if (ph_mid < win_lo || ph_mid > win_hi) {
                unanchored++;
                if (unanchored <= 50 || (unanchored % 1000) == 0)
                    fprintf(stderr,
                            "    island comp %d: %zu sample(s) "
                            "phi=[%.2f,%.2f] outside bulk window; carried "
                            "frame kept\n",
                            comp, nd, ph_lo, ph_hi);
                continue;
            }
        }
        qsort(diffs, nd, sizeof(double), cmp_dbl);
        double shift = diffs[nd / 2];
        for (size_t s = comp_off[comp]; s < comp_off[(size_t)comp + 1]; s++)
            S->smp[(size_t)comp_idx[s]].u += shift;
        if (fabs(shift) > max_shift) max_shift = fabs(shift);
        moved++;
        /* tens of thousands of singleton islands at the larger rungs: log
         * the first placements and a 1-in-1000 sample, never all of them */
        if (moved <= 200 || (moved % 1000) == 0)
            fprintf(stderr,
                    "    island comp %d: %zu sample(s) phi=[%.2f,%.2f] "
                    "carried u=[%.1f,%.1f] -> [%.1f,%.1f] (shift %+.1f)\n",
                    comp, nd, ph_lo, ph_hi, u_lo, u_hi,
                    u_lo + shift, u_hi + shift, shift);
    }
    if (max_shift > out->reg_max_shift) out->reg_max_shift = max_shift;
    fprintf(stderr,
            "  island placement: main comp %d, %zu map sample(s) (%s), "
            "phi->u bins=%d; %zu component(s) placed by carried phase, "
            "%zu kept carried (no phase), %zu kept carried (outside bulk "
            "window), max |shift|=%.1f\n",
            main_comp, nmain,
            solved_u_round ? "all in-window, solved u" : "dominant comp",
            M.nb, moved, unplaceable, unanchored, max_shift);
    Arena_restore(arena, mark);
}

/* Metric-projection Stage A: THE GRID IS THE SLICING.  A parameterized
 * quadribbon's param-v rows are its physical slices; cutting geometric
 * planes through a mesh whose components sit on phase-offset row lattices
 * (e.g. the 10x solid ribbon's v4 lattice sliced at 5, or 4x21x21's odd/even
 * interleaved components) splits rows apart and interpolates between them --
 * the shred/seam artifacts.  Build chains directly from each component's
 * rows instead: every finite-U vertex is a sample, carried U/phi are exact
 * vertex values, and the Stage-D transfer becomes distance-zero.  Chains
 * split at measured carried-U gaps (> 4x the row lattice pitch). */
typedef struct { int32_t comp; int64_t q; float u; int32_t vid; } RibRowKey;

static int cmp_rib_row_key(const void *pa, const void *pb)
{
    const RibRowKey *a = (const RibRowKey *)pa;
    const RibRowKey *b = (const RibRowKey *)pb;
    if (a->comp != b->comp) return a->comp < b->comp ? -1 : 1;
    if (a->q != b->q) return a->q < b->q ? -1 : 1;
    if (a->u != b->u) return a->u < b->u ? -1 : 1;
    return a->vid < b->vid ? -1 : (a->vid > b->vid ? 1 : 0);
}

static int cmp_i64(const void *pa, const void *pb)
{
    int64_t a = *(const int64_t *)pa, b = *(const int64_t *)pb;
    return a < b ? -1 : (a > b ? 1 : 0);
}

static int rib_rows_from_quadribbon(Arena_T arena,
                                    const float *verts, size_t nv,
                                    const int32_t *vertex_mesh_comp,
                                    const int32_t *mesh_comp_island,
                                    const int32_t *mesh_comp_chart,
                                    const double *t,
                                    const double e1[3], const double e2[3],
                                    const float axis_point[3],
                                    const RibbonOpts *o, SliceSet *S,
                                    int32_t **out_vertex_sample,
                                    float **out_sample_uref,
                                    float **out_sample_uconf,
                                    double *out_du_lat)
{
    memset(S, 0, sizeof(*S));
    *out_vertex_sample = NULL;
    *out_sample_uref = NULL;
    *out_sample_uconf = NULL;
    *out_du_lat = 1.0;
    if (o->reference_u == NULL || nv == 0) return -1;
    RibRowKey *key = RIB_ALLOC_ARRAY(arena, RibRowKey, nv);
    size_t nk = 0, skipped_nonfinite = 0;
    for (size_t i = 0; i < nv; i++) {
        double uu = o->reference_u[i];
        if (!isfinite(uu) || !isfinite(t[i])) {
            skipped_nonfinite++;
            continue;
        }
        key[nk].comp = vertex_mesh_comp[i];
        key[nk].q = (int64_t)llround(t[i] * 2.0);
        key[nk].u = (float)uu;
        key[nk].vid = (int32_t)i;
        nk++;
    }
    if (nk < 2) return -1;
    qsort(key, nk, sizeof(*key), cmp_rib_row_key);

    /* measured row-lattice pitch -> chain split gate */
    double *dus = RIB_ALLOC_ARRAY(arena, double, nk);
    size_t ndu = 0;
    for (size_t k = 1; k < nk; k++) {
        if (key[k].comp != key[k-1].comp || key[k].q != key[k-1].q) continue;
        double du_ = (double)key[k].u - (double)key[k-1].u;
        if (du_ > 1e-9) dus[ndu++] = du_;
    }
    double du_lat = ndu != 0 ? select_median_dbl(dus, ndu) : 1.0;
    if (du_lat < 0.25) du_lat = 0.25;
    double gap_gate = 4.0 * du_lat + 1.0;

    /* global slice rank per distinct q (component lattices interleave) */
    int64_t *qs = RIB_ALLOC_ARRAY(arena, int64_t, nk);
    for (size_t k = 0; k < nk; k++) qs[k] = key[k].q;
    qsort(qs, nk, sizeof(*qs), cmp_i64);
    size_t nq = 0;
    for (size_t k = 0; k < nk; k++)
        if (k == 0 || qs[k] != qs[nq - 1]) qs[nq++] = qs[k];

    /* count chains, then fill */
    size_t n_chn = 0;
    for (size_t k = 0; k < nk; k++) {
        int newrow = k == 0 || key[k].comp != key[k-1].comp ||
                     key[k].q != key[k-1].q;
        int gap = !newrow &&
                  (double)key[k].u - (double)key[k-1].u > gap_gate;
        if (newrow || gap) n_chn++;
    }
    Sample *smp = (Sample *)ARENA_CALLOC(arena, nk, sizeof(*smp));
    Chain *chn = (Chain *)ARENA_CALLOC(arena, n_chn, sizeof(*chn));
    size_t ci = 0;
    for (size_t k = 0; k < nk; k++) {
        int newrow = k == 0 || key[k].comp != key[k-1].comp ||
                     key[k].q != key[k-1].q;
        int gap = !newrow &&
                  (double)key[k].u - (double)key[k-1].u > gap_gate;
        if (newrow || gap) {
            Chain *c = &chn[ci];
            int64_t *found = (int64_t *)bsearch(
                &key[k].q, qs, nq, sizeof(*qs), cmp_i64);
            c->first = (int32_t)k;
            c->count = 0;
            c->slice = found != NULL ? (int32_t)(found - qs) : 0;
            c->closed = 0;
            c->wind = 0;
            c->group = -1;
            c->mesh_comp = key[k].comp;
            c->source_chart = mesh_comp_chart[key[k].comp];
            c->winding_island = mesh_comp_island[key[k].comp];
            c->reconstruction_component = c->winding_island;
            c->atlas_shift = 0.0;
            c->solve_comp = 0;
            ci++;
        }
        Chain *c = &chn[ci - 1];
        Sample *sp = &smp[k];
        size_t vi = (size_t)key[k].vid;
        sp->p[0] = (double)verts[vi*3+0];
        sp->p[1] = (double)verts[vi*3+1];
        sp->p[2] = (double)verts[vi*3+2];
        double d0 = sp->p[0] - (double)axis_point[0];
        double d1 = sp->p[1] - (double)axis_point[1];
        double d2 = sp->p[2] - (double)axis_point[2];
        sp->c1 = d0*e1[0] + d1*e1[1] + d2*e1[2];
        sp->c2 = d0*e2[0] + d1*e2[1] + d2*e2[2];
        sp->r = sqrt(sp->c1*sp->c1 + sp->c2*sp->c2);
        sp->th = atan2(sp->c2, sp->c1);
        sp->phi = o->reference_phi != NULL &&
                  isfinite((double)o->reference_phi[vi])
                ? (double)o->reference_phi[vi] : sp->th;
        sp->u = o->reference_u[vi];
        sp->chain = (int32_t)(ci - 1);
        sp->slice = c->slice;
        if (c->count == 0) {
            sp->s = 0.0;
        } else {
            const Sample *pp = &smp[k-1];
            double dz = sp->p[0]-pp->p[0], dy = sp->p[1]-pp->p[1];
            double dx = sp->p[2]-pp->p[2];
            sp->s = pp->s + sqrt(dz*dz + dy*dy + dx*dx);
        }
        c->count++;
    }
    /* chain-local tangents from neighbours (one-sided at the ends) */
    for (size_t c = 0; c < n_chn; c++) {
        int32_t f = chn[c].first, cn = chn[c].count;
        for (int32_t i = 0; i < cn; i++) {
            int32_t i0 = i > 0 ? i - 1 : i;
            int32_t i1 = i + 1 < cn ? i + 1 : i;
            double dz = smp[f+i1].p[0] - smp[f+i0].p[0];
            double dy = smp[f+i1].p[1] - smp[f+i0].p[1];
            double dx = smp[f+i1].p[2] - smp[f+i0].p[2];
            double nn = sqrt(dz*dz + dy*dy + dx*dx);
            if (nn < 1e-12) { dz = 0.0; dy = 0.0; dx = 1.0; nn = 1.0; }
            smp[f+i].tau[0] = dz / nn;
            smp[f+i].tau[1] = dy / nn;
            smp[f+i].tau[2] = dx / nn;
        }
    }
    size_t *plane_chains = (size_t *)ARENA_CALLOC(
        arena, nq, sizeof(*plane_chains));
    for (size_t c = 0; c < n_chn; c++)
        plane_chains[chn[c].slice]++;
    int hit = 0, multi = 0;
    for (size_t p = 0; p < nq; p++) {
        if (plane_chains[p] != 0) hit++;
        if (plane_chains[p] > 1) multi++;
    }
    S->smp = smp;
    S->n_smp = nk;
    S->chn = chn;
    S->n_chn = n_chn;
    S->nplanes = (int)nq;
    S->n_slices_hit = hit;
    S->n_multi = multi;
    S->n_closed = 0;
    S->bridge_cuts = 0;
    {
        int32_t *v2s = RIB_ALLOC_ARRAY(arena, int32_t, nv);
        for (size_t i = 0; i < nv; i++) v2s[i] = -1;
        for (size_t k = 0; k < nk; k++) v2s[(size_t)key[k].vid] = (int32_t)k;
        *out_vertex_sample = v2s;
    }
    {
        /* raw carried u per sample, captured BEFORE any monotone remap of
         * u -- the consistency solve classifies grid edges with it */
        float *uref_raw = RIB_ALLOC_ARRAY(arena, float, nk);
        for (size_t k = 0; k < nk; k++) uref_raw[k] = key[k].u;
        *out_sample_uref = uref_raw;
        *out_du_lat = du_lat;
    }
    if (o->reference_u_confidence != NULL) {
        float *uconf = RIB_ALLOC_ARRAY(arena, float, nk);
        for (size_t k = 0; k < nk; k++) {
            double q = (double)o->reference_u_confidence[(size_t)key[k].vid];
            if (!isfinite(q)) q = 0.0;
            if (q < 0.0) q = 0.0;
            if (q > 1.0) q = 1.0;
            uconf[k] = (float)q;
        }
        *out_sample_uconf = uconf;
    }
    fprintf(stderr,
            "  quadribbon row slicing: %zu rows (%zu lattice planes), "
            "%zu chains, %zu samples (= finite-U vertices; %zu skipped), "
            "row pitch %.2f, chain gap gate %.2f\n",
            (size_t)hit, nq, n_chn, nk, skipped_nonfinite, du_lat, gap_gate);
    return 0;
}

/* Metric-projection Stage 2: ONE drift-penalized consistency solve, done
 * as iteratively reweighted least squares over explicitly assembled normal
 * equations, factored exactly by TAUCS (Sparse_solve_sym, supernodal
 * multifrontal Cholesky).  Variables are the sample u values.  Terms, all
 * w*(u_a - u_b - d)^2:
 *   - row metric: consecutive samples of each row chain (d = measured ds,
 *     w = 1/ds; the trusted measurement backbone, never reweighted);
 *   - in-chart vertical drift: cross-row mesh edges classified by the RAW
 *     carried |dU_ref| with per-component pitch (verticals kept, diagonals
 *     excluded, block-offset boundaries rescued when geometrically
 *     vertical), d = 0, w = LAMBDA_V/len;
 *   - cross-chart contact drift: same-row 3D contacts between different
 *     source components, gated so plies are never tied (in-plane
 *     separation, tangent alignment, authoritative-phase agreement),
 *     d = 0, w = LAMBDA_C/len;
 *   - one weight-1 pin per connected component (kills the translation
 *     null space; the component's actual placement is restored after the
 *     solve from the Stage-1 anchoring, see below).
 * IRLS: each round is an exact solve; drift/contact terms are then
 * Huber-reweighted from the solved residuals (a term whose residual stays
 * large is evidence the model tied the wrong samples) and the system is
 * re-solved, at most RIB_PCONS_IRLS_ROUNDS times.  Residual wiggle
 * (drift/contact quantiles) and arc-length error print per round. */
static const double RIB_PCONS_LAMBDA_V = 32.0;
static const double RIB_PCONS_LAMBDA_C = 4.0;
static const double RIB_PCONS_PERP_MAX = 3.5;  /* normal-separation ply gate:
                                * a true ply sits >= ~7 vox away (and the
                                * authoritative phase gate rejects it anyway);
                                * 1.5 rejected REAL hole-rim continuations at
                                * perp 2.8-3.0 with tdot=1.00 (measured 4x),
                                * leaving their charts untied -> visible slips */
static const double RIB_PCONS_HUBER_C = 1.0;   /* vox: vertical drift scale */
/* contacts get a wider Huber: they are SPARSE evidence already filtered by
 * hard geometric/phase gates, and they must pull multi-voxel chart slips
 * home; c=1.0 released real 3-5 px lane seams (measured 4x) while the dense
 * verticals still need the tight scale to reject genuine defects */
static const double RIB_PCONS_HUBER_C_CONTACT = 6.0;
enum { RIB_PCONS_IRLS_ROUNDS = 3 };

typedef struct {
    double value;
    double weight;
} RibMetricGaugeVote;

static int cmp_rib_metric_gauge_vote(const void *pa, const void *pb)
{
    const RibMetricGaugeVote *a = (const RibMetricGaugeVote *)pa;
    const RibMetricGaugeVote *b = (const RibMetricGaugeVote *)pb;
    return a->value < b->value ? -1 : (a->value > b->value ? 1 : 0);
}

static double rib_metric_weighted_median(RibMetricGaugeVote *vote, size_t n)
{
    if (n == 0) return 0.0;
    qsort(vote, n, sizeof(*vote), cmp_rib_metric_gauge_vote);
    double total = 0.0;
    for (size_t i = 0; i < n; i++)
        if (vote[i].weight > 0.0 && isfinite(vote[i].weight))
            total += vote[i].weight;
    if (!(total > 0.0)) return vote[n / 2].value;
    double accum = 0.0;
    for (size_t i = 0; i < n; i++) {
        if (vote[i].weight > 0.0 && isfinite(vote[i].weight))
            accum += vote[i].weight;
        if (accum >= 0.5 * total) return vote[i].value;
    }
    return vote[n - 1].value;
}

/* Project the common certificate map onto the exact metric family of every
 * observed row,
 *
 *                         u_i = g_chain + s_i.
 *
 * Only the additive gauge is estimated.  A posterior-backed, locally metric
 * observation is the best gauge vote; when a short/torn row has none, the
 * selection falls back in explicit stages (posterior, local metric, all)
 * rather than leaving the row unconstrained.  Confidence weights affect only
 * the robust median -- they can never change an intra-row derivative.  This is
 * the certificate-preserving replacement for clipping individual samples
 * after the old full consistency solve, which measured 55% bad derivatives on
 * the 4x5x5 control despite a correct winding order. */
static int rib_project_exact_chain_metric(Arena_T arena, SliceSet *S,
                                           const float *sample_uconf,
                                           RibbonResult *out)
{
    if (S == NULL || S->n_chn == 0) return -1;
    Arena_Mark mark = Arena_save(arena);
    int32_t max_count = 1;
    for (size_t c = 0; c < S->n_chn; c++)
        if (S->chn[c].count > max_count) max_count = S->chn[c].count;
    RibMetricGaugeVote *vote = RIB_ALLOC_ARRAY(
        arena, RibMetricGaugeVote, (size_t)max_count);
    double *displacement = RIB_ALLOC_ARRAY(arena, double, S->n_smp + 1);
    size_t ndisp = 0, confident_votes = 0, fallback_chains = 0;
    size_t exact_edges = 0;

    for (size_t c = 0; c < S->n_chn; c++) {
        int32_t f = S->chn[c].first, count = S->chn[c].count;
        if (count <= 0) continue;
        size_t nv = 0;
        int selected_pass = -1;
        /* pass 0: posterior + locally metric; 1: posterior; 2: locally
         * metric; 3: every finite sample.  Without a posterior, pass 0 is
         * simply the locally-metric tier and pass 1 is skipped. */
        for (int pass = 0; pass < 4 && nv == 0; pass++) {
            if (sample_uconf == NULL && pass == 1) continue;
            for (int32_t i = 0; i < count; i++) {
                size_t q = (size_t)(f + i);
                double u = S->smp[q].u, s = S->smp[q].s;
                if (!isfinite(u) || !isfinite(s)) continue;
                double conf = sample_uconf != NULL
                            ? (double)sample_uconf[q] : 1.0;
                if (!isfinite(conf)) conf = 0.0;
                if (conf < 0.0) conf = 0.0;
                if (conf > 1.0) conf = 1.0;
                int posterior = sample_uconf != NULL && conf >= 0.75;
                int32_t ia = i > 0 ? i - 1 : i;
                int32_t ib = i + 1 < count ? i + 1 : i;
                int metric = 0;
                if (ib > ia) {
                    double ds = S->smp[f + ib].s - S->smp[f + ia].s;
                    double du = S->smp[f + ib].u - S->smp[f + ia].u;
                    metric = ds > 1e-9 && isfinite(du) &&
                             fabs(du / ds - 1.0) <= 0.15;
                }
                int accept = pass == 0 ? (sample_uconf != NULL
                                           ? posterior && metric : metric)
                           : pass == 1 ? posterior
                           : pass == 2 ? metric : 1;
                if (!accept) continue;
                vote[nv].value = u - s;
                /* Posterior odds are deliberately bounded: one near-certain
                 * sample may lead a short row, but cannot overwhelm a rowful
                 * of agreeing observations. */
                vote[nv].weight = sample_uconf != NULL
                                ? 0.05 + conf * conf : 1.0;
                if (metric) vote[nv].weight *= 2.0;
                nv++;
            }
            if (nv != 0) selected_pass = pass;
        }
        if (nv == 0) {
            Arena_restore(arena, mark);
            return -1;
        }
        if (selected_pass > 0) fallback_chains++;
        if (sample_uconf != NULL && selected_pass <= 1)
            confident_votes += nv;
        double gauge = rib_metric_weighted_median(vote, nv);
        for (int32_t i = 0; i < count; i++) {
            size_t q = (size_t)(f + i);
            double next = gauge + S->smp[q].s;
            displacement[ndisp++] = fabs(next - S->smp[q].u);
            S->smp[q].u = next;
            if (i > 0 && S->smp[q].s - S->smp[q - 1].s > 1e-9)
                exact_edges++;
        }
    }

    double p95 = ndisp ? select_p95_dbl(displacement, ndisp) : 0.0;
    double dmax = 0.0;
    for (size_t i = 0; i < ndisp; i++)
        if (displacement[i] > dmax) dmax = displacement[i];
    out->certificate_samples += confident_votes;
    out->certificate_drift_p95 = p95;
    out->certificate_drift_max = dmax;
    fprintf(stderr,
            "  exact chain metric: %zu chains, %zu edges; posterior votes=%zu "
            "fallback chains=%zu; projection |dU| p95/max=%.2f/%.2f vox\n",
            S->n_chn, exact_edges, confident_votes, fallback_chains,
            p95, dmax);
    Arena_restore(arena, mark);
    return 0;
}

/* A carried U coordinate is useful evidence, but it cannot be allowed to veto
 * the very observation that diagnoses a carried block offset.  The input mesh
 * supplies stronger evidence: two vertices joined by a face edge are on the
 * same material surface.  Quotient those edges directly to additive row-gauge
 * relations,
 *
 *   g_a - g_b = dot(tau_bar, p_a-p_b) - s_a + s_b,
 *
 * and robustly combine all observations for each chain pair.  This retains
 * both true verticals and quad diagonals (their tangent projection supplies
 * the correct nonzero rest offset) without consulting U_ref.  The quotient is
 * important: a million-vertex quadribbon normally yields only O(rows) gauge
 * relations, so this stage stays cheap on whole-scroll meshes. */
static const double RIB_METRIC_TOPOLOGY_WEIGHT = 1024.0;
static const double RIB_METRIC_TOPOLOGY_CERT_WEIGHT = 8.0;
static const double RIB_METRIC_TOPOLOGY_HUBER = 1.0;

typedef struct {
    uint64_t key;
    double sum;
    double robust_sum;
    double confidence_sum;
    uint64_t count;
    uint64_t robust_count;
    int32_t sa, sb;
} RibMetricTopologyBin;

typedef struct {
    RibMetricTopologyBin *bin;
    size_t capacity;
    size_t count;
} RibMetricTopologyTable;

static int rib_metric_topology_rehash(RibMetricTopologyTable *tab,
                                      size_t capacity)
{
    RibMetricTopologyBin *next;
    if (capacity < 16 || (capacity & (capacity - 1)) != 0 ||
        capacity > SIZE_MAX / sizeof(*next))
        return -1;
    next = (RibMetricTopologyBin *)malloc(capacity * sizeof(*next));
    if (next == NULL) return -1;
    for (size_t i = 0; i < capacity; i++) next[i].key = UINT64_MAX;
    if (tab->bin != NULL) {
        size_t mask = capacity - 1;
        for (size_t i = 0; i < tab->capacity; i++) {
            if (tab->bin[i].key == UINT64_MAX) continue;
            size_t slot = chain_relation_slot(tab->bin[i].key, mask);
            while (next[slot].key != UINT64_MAX) slot = (slot + 1) & mask;
            next[slot] = tab->bin[i];
        }
    }
    free(tab->bin);
    tab->bin = next;
    tab->capacity = capacity;
    return 0;
}

/* Return one canonical chain-pair observation from a mesh edge. */
static int rib_metric_topology_observation(
    const SliceSet *S, const int32_t *vertex_sample, size_t nv,
    int32_t va, int32_t vb, const RibbonOpts *o,
    int32_t *out_sa, int32_t *out_sb, double *out_target)
{
    if (va < 0 || vb < 0 || (size_t)va >= nv || (size_t)vb >= nv)
        return 0;
    int32_t sa = vertex_sample[va], sb = vertex_sample[vb];
    if (sa < 0 || sb < 0 || sa == sb ||
        (size_t)sa >= S->n_smp || (size_t)sb >= S->n_smp)
        return 0;
    int32_t ca = S->smp[sa].chain, cb = S->smp[sb].chain;
    if (ca < 0 || cb < 0 || ca == cb ||
        (size_t)ca >= S->n_chn || (size_t)cb >= S->n_chn)
        return 0;
    if (S->chn[ca].winding_island != S->chn[cb].winding_island)
        return 0;
    if (o->reference_phi != NULL &&
        (!isfinite(S->smp[sa].phi) || !isfinite(S->smp[sb].phi) ||
         fabs(S->smp[sa].phi - S->smp[sb].phi) > RIB_CONT_DPHI_MAX))
        return 0;

    double cos_gate = cos((double)o->match_ang_deg * M_PI / 180.0);
    double tdot = v3dot(S->smp[sa].tau, S->smp[sb].tau);
    if (!isfinite(tdot) || tdot < cos_gate) return 0;
    if (ca > cb) {
        int32_t tmp = sa; sa = sb; sb = tmp;
        tmp = ca; ca = cb; cb = tmp;
    }
    double tau[3] = { S->smp[sa].tau[0] + S->smp[sb].tau[0],
                      S->smp[sa].tau[1] + S->smp[sb].tau[1],
                      S->smp[sa].tau[2] + S->smp[sb].tau[2] };
    double tn = sqrt(v3dot(tau, tau));
    if (!(tn > 1e-12) || !isfinite(tn)) return 0;
    tau[0] /= tn; tau[1] /= tn; tau[2] /= tn;
    double dp[3] = { S->smp[sa].p[0] - S->smp[sb].p[0],
                     S->smp[sa].p[1] - S->smp[sb].p[1],
                     S->smp[sa].p[2] - S->smp[sb].p[2] };
    double target = v3dot(tau, dp) - S->smp[sa].s + S->smp[sb].s;
    if (!isfinite(target)) return 0;
    *out_sa = sa;
    *out_sb = sb;
    *out_target = target;
    return 1;
}

static int rib_metric_add_topology_gauge_pairs(
    Arena_T arena, const int32_t *faces, size_t nf,
    const int32_t *vertex_sample, size_t nv, const SliceSet *S,
    const float *sample_uconf, const RibbonOpts *o,
    const PairSet *base, PairSet *out)
{
    RibMetricTopologyTable tab;
    memset(&tab, 0, sizeof(tab));
    if (vertex_sample == NULL || faces == NULL || nf == 0) {
        *out = *base;
        return 0;
    }
    if (rib_metric_topology_rehash(&tab, 1024) != 0) return -1;

    size_t observations = 0, same_v = 0, cross_v = 0;
    for (size_t f = 0; f < nf; f++) for (int e = 0; e < 3; e++) {
        int32_t sa, sb;
        double target;
        if (!rib_metric_topology_observation(
                S, vertex_sample, nv, faces[3*f+e], faces[3*f+(e+1)%3],
                o, &sa, &sb, &target))
            continue;
        int32_t ca = S->smp[sa].chain, cb = S->smp[sb].chain;
        uint64_t key = chain_relation_key(ca, cb);
        if ((tab.count + 1) * 10 >= tab.capacity * 7) {
            if (tab.capacity > SIZE_MAX / 2 ||
                rib_metric_topology_rehash(&tab, tab.capacity * 2) != 0) {
                free(tab.bin);
                return -1;
            }
        }
        size_t slot = chain_relation_slot(key, tab.capacity - 1);
        while (tab.bin[slot].key != UINT64_MAX &&
               tab.bin[slot].key != key)
            slot = (slot + 1) & (tab.capacity - 1);
        if (tab.bin[slot].key == UINT64_MAX) {
            memset(&tab.bin[slot], 0, sizeof(tab.bin[slot]));
            tab.bin[slot].key = key;
            tab.bin[slot].sa = sa;
            tab.bin[slot].sb = sb;
            tab.count++;
        }
        tab.bin[slot].sum += target;
        if (sample_uconf != NULL) {
            double qa = (double)sample_uconf[sa];
            double qb = (double)sample_uconf[sb];
            if (!isfinite(qa)) qa = 0.0;
            if (!isfinite(qb)) qb = 0.0;
            if (qa < 0.0) qa = 0.0; if (qa > 1.0) qa = 1.0;
            if (qb < 0.0) qb = 0.0; if (qb > 1.0) qb = 1.0;
            /* A binary relation is only as certificate-backed as its weaker
             * endpoint.  Low-confidence rows deliberately hand authority to
             * topology; two certain rows retain the certificate frame. */
            tab.bin[slot].confidence_sum += qa < qb ? qa : qb;
        }
        tab.bin[slot].count++;
        observations++;
        if (S->smp[sa].slice == S->smp[sb].slice) same_v++;
        else cross_v++;
    }

    /* A second linear mesh pass gives a winsorized mean.  It costs much less
     * than sorting all face edges and prevents one long/sliver triangle from
     * moving an otherwise unanimous chain relation. */
    for (size_t f = 0; f < nf; f++) for (int e = 0; e < 3; e++) {
        int32_t sa, sb;
        double target;
        if (!rib_metric_topology_observation(
                S, vertex_sample, nv, faces[3*f+e], faces[3*f+(e+1)%3],
                o, &sa, &sb, &target))
            continue;
        uint64_t key = chain_relation_key(S->smp[sa].chain,
                                          S->smp[sb].chain);
        size_t slot = chain_relation_slot(key, tab.capacity - 1);
        while (tab.bin[slot].key != key) {
            if (tab.bin[slot].key == UINT64_MAX) break;
            slot = (slot + 1) & (tab.capacity - 1);
        }
        if (tab.bin[slot].key != key || tab.bin[slot].count == 0) continue;
        double mean = tab.bin[slot].sum / (double)tab.bin[slot].count;
        double r = target - mean;
        if (r > RIB_METRIC_TOPOLOGY_HUBER) r = RIB_METRIC_TOPOLOGY_HUBER;
        if (r < -RIB_METRIC_TOPOLOGY_HUBER) r = -RIB_METRIC_TOPOLOGY_HUBER;
        tab.bin[slot].robust_sum += mean + r;
        tab.bin[slot].robust_count++;
    }

    if (base->n_pairs > SIZE_MAX - tab.count) {
        free(tab.bin);
        return -1;
    }
    Pair *pairs = RIB_ALLOC_ARRAY(arena, Pair, base->n_pairs + tab.count + 1);
    if (base->n_pairs != 0)
        memcpy(pairs, base->pairs, base->n_pairs * sizeof(*pairs));
    size_t np = base->n_pairs;
    double support_sum = 0.0, weighted_base_sum = 0.0;
    for (size_t slot = 0; slot < tab.capacity; slot++) {
        RibMetricTopologyBin *b = &tab.bin[slot];
        if (b->key == UINT64_MAX || b->count == 0) continue;
        double target = b->robust_count != 0
                      ? b->robust_sum / (double)b->robust_count
                      : b->sum / (double)b->count;
        Pair *pr = &pairs[np++];
        pr->a = b->sa;
        pr->b = b->sb;
        pr->k = RIB_PAIR_K_TRUSTED_TOPOLOGY;
        pr->d = target + S->smp[b->sa].s - S->smp[b->sb].s;
        /* Interior face edges are normally seen twice, once from each
         * incident triangle.  Half-counting makes support approximately the
         * number of independent vertices while retaining boundary evidence. */
        double support = 0.5 * (double)b->count;
        if (support < 1.0) support = 1.0;
        double confidence = sample_uconf != NULL
                          ? b->confidence_sum / (double)b->count : 0.0;
        double uncertainty = 1.0 - confidence;
        double base_weight = RIB_METRIC_TOPOLOGY_CERT_WEIGHT +
            (RIB_METRIC_TOPOLOGY_WEIGHT - RIB_METRIC_TOPOLOGY_CERT_WEIGHT) *
            uncertainty * uncertainty;
        pr->w = base_weight * support;
        pr->like = 1.0;
        support_sum += support;
        weighted_base_sum += support * base_weight;
    }
    *out = *base;
    out->pairs = pairs;
    out->n_pairs = np;
    fprintf(stderr,
            "  exact-chain topology: %zu face-edge observations -> %zu "
            "chain relations (same-v=%zu cross-v=%zu, effective "
            "support=%.0f, mean weight=%.1f/support, range=%.0f..%.0f)\n",
            observations, tab.count, same_v, cross_v, support_sum,
            support_sum > 0.0 ? weighted_base_sum / support_sum : 0.0,
            RIB_METRIC_TOPOLOGY_CERT_WEIGHT, RIB_METRIC_TOPOLOGY_WEIGHT);
    free(tab.bin);
    return 0;
}

/* ARMING GUARD: the consistency solve is a REFINEMENT stage -- pieces must
 * be placed before it runs (island placement, or an inherently coherent
 * frame).  Contacts entering at packing scale mean the frame is unplaced
 * and the solve would be asked to do placement through millions of stiff
 * long-range ties (measured sec00 failure: +-103k entry drift, Cholesky
 * fill explosion).  The 0143-hybrid's honest wave entered at +-1300 and
 * must pass; an atlas packing at +-100k must refuse loudly. */
static const double RIB_PCONS_ENTRY_DRIFT_MAX = 10000.0;
/* Factor-once PCG: the TAUCS factor sees only the banded grid (rows +
 * verticals + pins, round-constant); contacts and Huber reweights live in
 * the matrix-free operator PCG applies.  The grid ridge keeps the factor
 * SPD: pins are assigned on the FULL solve graph (which includes contacts),
 * so fragments held together only by contacts would otherwise float in the
 * contact-free grid matrix and zero-pivot the Cholesky.  Preconditioner
 * only -- the PCG operator stays exact. */
static const double RIB_PCONS_PCG_TOL = 1e-6;
static const double RIB_PCONS_GRID_RIDGE = 1e-4;
enum { RIB_PCONS_PCG_MAXIT = 600 };
/* Below this size the EXACT per-round factor is cheaper than PCG's matvec
 * budget (measured at the 5%-density operating point: exact 14 s/round at
 * 0.93M vs 200 s of capped PCG); above it the connected multifrontal factor
 * is the pathology and the grid-PCG path wins. */
enum { RIB_PCONS_PCG_MIN_N = 4000000 };

typedef struct { int64_t key; double delta; } RibPairDelta;

static int cmp_rib_pair_delta(const void *pa, const void *pb)
{
    const RibPairDelta *a = (const RibPairDelta *)pa;
    const RibPairDelta *b = (const RibPairDelta *)pb;
    if (a->key != b->key) return a->key < b->key ? -1 : 1;
    return a->delta < b->delta ? -1 : (a->delta > b->delta ? 1 : 0);
}

static void rib_pcons_metrics(const SliceSet *S, const PairSet *Q,
                              size_t nvert_pairs, const double *u,
                              Arena_T arena, const char *tag)
{
    Arena_Mark mark = Arena_save(arena);
    size_t nq = Q->n_pairs;
    double *buf = RIB_ALLOC_ARRAY(arena, double, nq + S->n_smp + 1);
    size_t nd = 0;
    for (size_t p = 0; p < nvert_pairs; p++)
        buf[nd++] = fabs(u[Q->pairs[p].a] - u[Q->pairs[p].b]);
    double v50 = nd ? select_median_dbl(buf, nd) : 0.0;
    for (size_t p = 0; p < nvert_pairs; p++)
        buf[p] = fabs(u[Q->pairs[p].a] - u[Q->pairs[p].b]);
    double v95 = nd ? select_p95_dbl(buf, nd) : 0.0;
    size_t nc = 0;
    for (size_t p = nvert_pairs; p < nq; p++)
        buf[nc++] = u[Q->pairs[p].a] - u[Q->pairs[p].b];
    double c10 = 0.0, c50 = 0.0, c90 = 0.0;
    if (nc != 0) {
        qsort(buf, nc, sizeof(double), cmp_dbl);
        c10 = buf[nc / 10];
        c50 = buf[nc / 2];
        c90 = buf[(9 * nc) / 10];
    }
    size_t ne = 0;
    double esum = 0.0;
    size_t nbig = 0;
    for (size_t c = 0; c < S->n_chn; c++) {
        int32_t f = S->chn[c].first, cn = S->chn[c].count;
        for (int32_t i = 1; i < cn; i++) {
            double l = S->smp[f+i].s - S->smp[f+i-1].s;
            if (l < 1e-9) continue;
            double e = fabs((u[f+i] - u[f+i-1]) / l - 1.0);
            buf[ne++] = e;
            esum += e;
            if (e > 0.05) nbig++;
        }
    }
    double e95 = ne ? select_p95_dbl(buf, ne) : 0.0;
    fprintf(stderr,
            "  consistency %s: wiggle |dU|v p50/p95=%.3f/%.3f  contact "
            "p10/p50/p90=%.2f/%.2f/%.2f (n=%zu)  |du/ds-1| mean/p95=%.4f/"
            "%.4f (>5%%: %.2f%%)\n",
            tag, v50, v95, c10, c50, c90, nc,
            ne ? esum / (double)ne : 0.0, e95,
            ne ? 100.0 * (double)nbig / (double)ne : 0.0);
    Arena_restore(arena, mark);
}

static int rib_project_consistency_solve(Arena_T arena,
                                         const float *verts, size_t nv,
                                         const int32_t *faces, size_t nf,
                                         const int32_t *vertex_sample,
                                         const float *sample_uref,
                                         const float *sample_uconf,
                                         double du_lat,
                                         const RibbonOpts *o,
                                         SliceSet *S, RibbonResult *out)
{
    double t0 = ves_clock_sec();
    size_t n = S->n_smp;
    (void)verts;
    (void)nv;
    if (n == 0) return 0;
    Arena_Mark mark = Arena_save(arena);

    /* --- vertical drift edges: cross-row face edges, deduped, gated by
     * geometry (never by the possibly-defective carried coordinate) --- */
    int64_t *ekey = RIB_ALLOC_ARRAY(arena, int64_t, nf * 3 + 1);
    size_t nek = 0;
    for (size_t f = 0; f < nf; f++) {
        for (int k = 0; k < 3; k++) {
            int32_t a = faces[f*3+k], b = faces[f*3+(k+1)%3];
            int32_t sa = vertex_sample[a], sb = vertex_sample[b];
            if (sa < 0 || sb < 0) continue;
            if (S->smp[sa].slice == S->smp[sb].slice) continue;
            int64_t lo = sa < sb ? sa : sb, hi = sa < sb ? sb : sa;
            ekey[nek++] = (lo << 31) | hi;
        }
    }
    qsort(ekey, nek, sizeof(*ekey), cmp_i64);
    size_t nuniq = 0;
    for (size_t k = 0; k < nek; k++)
        if (k == 0 || ekey[k] != ekey[nuniq - 1]) ekey[nuniq++] = ekey[k];
    double *elen = RIB_ALLOC_ARRAY(arena, double, nuniq + 1);
    for (size_t k = 0; k < nuniq; k++) {
        int32_t sa = (int32_t)(ekey[k] >> 31);
        int32_t sb = (int32_t)(ekey[k] & 0x7fffffff);
        double dz = S->smp[sa].p[0] - S->smp[sb].p[0];
        double dy = S->smp[sa].p[1] - S->smp[sb].p[1];
        double dx = S->smp[sa].p[2] - S->smp[sb].p[2];
        elen[k] = sqrt(dz*dz + dy*dy + dx*dx);
    }
    /* Classify cross-row edges by RAW CARRIED |dU_ref|: a true grid
     * vertical carries exactly 0 (the fit's own claim -> drift target 0);
     * a quad diagonal carries about one column step, and feeding those a
     * target of 0 fights the row metric (measured: it wrecked the spiral
     * selftest); a carried BLOCK-OFFSET boundary carries far more than a
     * column yet is geometrically a normal vertical edge -- rescued (kept)
     * within a tight per-component pitch gate so the solve can pull the
     * block back, and counted loudly.  Pitches are per component (mixed
     * lattices must not cross-gate). */
    int32_t ncomp_v = 0;
    for (size_t c = 0; c < S->n_chn; c++)
        if (S->chn[c].mesh_comp >= ncomp_v) ncomp_v = S->chn[c].mesh_comp + 1;
    if (ncomp_v < 1) ncomp_v = 1;
    size_t *vc_cnt = (size_t *)ARENA_CALLOC(arena, (size_t)ncomp_v,
                                            sizeof(*vc_cnt));
    size_t *vc_fill = (size_t *)ARENA_CALLOC(arena, (size_t)ncomp_v,
                                             sizeof(*vc_fill));
    /* per-component carried column pitch (a detached arc on a 2x lattice
     * must not be classified with the dominant component's pitch) */
    double *comp_du = RIB_ALLOC_ARRAY(arena, double, (size_t)ncomp_v);
    {
        size_t *cd_cnt = (size_t *)ARENA_CALLOC(arena, (size_t)ncomp_v,
                                                sizeof(*cd_cnt));
        size_t *cd_fill = (size_t *)ARENA_CALLOC(arena, (size_t)ncomp_v,
                                                 sizeof(*cd_fill));
        for (size_t c = 0; c < S->n_chn; c++) {
            int32_t cn = S->chn[c].count;
            if (cn > 1)
                cd_cnt[S->chn[c].mesh_comp] += (size_t)cn - 1;
        }
        size_t *cd_off = (size_t *)ARENA_CALLOC(arena, (size_t)ncomp_v + 1,
                                                sizeof(*cd_off));
        for (int32_t c = 0; c < ncomp_v; c++)
            cd_off[c + 1] = cd_off[c] + cd_cnt[c];
        double *cd_val = RIB_ALLOC_ARRAY(arena, double, cd_off[ncomp_v] + 1);
        for (size_t c = 0; c < S->n_chn; c++) {
            int32_t f2 = S->chn[c].first, cn = S->chn[c].count;
            int32_t mc = S->chn[c].mesh_comp;
            for (int32_t i = 1; i < cn; i++)
                cd_val[cd_off[mc] + cd_fill[mc]++] =
                    fabs((double)sample_uref[f2+i] -
                         (double)sample_uref[f2+i-1]);
        }
        for (int32_t c = 0; c < ncomp_v; c++)
            comp_du[c] = cd_fill[c] >= 8
                       ? select_median_dbl(&cd_val[cd_off[c]], cd_fill[c])
                       : du_lat;
    }
    int32_t *ecomp = RIB_ALLOC_ARRAY(arena, int32_t, nuniq + 1);
    int8_t *eclass = RIB_ALLOC_ARRAY(arena, int8_t, nuniq + 1);
    for (size_t k = 0; k < nuniq; k++) {
        int32_t sa = (int32_t)(ekey[k] >> 31);
        int32_t sb = (int32_t)(ekey[k] & 0x7fffffff);
        double duref = fabs((double)sample_uref[sa] -
                            (double)sample_uref[sb]);
        ecomp[k] = S->chn[S->smp[sa].chain].mesh_comp;
        double cdu = comp_du[ecomp[k]];
        eclass[k] = duref <= 0.25 * cdu ? 0
                  : (duref < 1.75 * cdu ? 1 : 2);
        if (eclass[k] == 0) vc_cnt[ecomp[k]]++;
    }
    size_t *vc_off = (size_t *)ARENA_CALLOC(arena, (size_t)ncomp_v + 1,
                                            sizeof(*vc_off));
    for (int32_t c = 0; c < ncomp_v; c++)
        vc_off[c + 1] = vc_off[c] + vc_cnt[c];
    double *vc_len = RIB_ALLOC_ARRAY(arena, double, nuniq + 1);
    for (size_t k = 0; k < nuniq; k++)
        if (eclass[k] == 0)
            vc_len[vc_off[ecomp[k]] + vc_fill[ecomp[k]]++] = elen[k];
    double *vc_gate = RIB_ALLOC_ARRAY(arena, double, (size_t)ncomp_v);
    for (int32_t c = 0; c < ncomp_v; c++)
        vc_gate[c] = vc_fill[c] >= 8
                   ? 1.1 * select_median_dbl(&vc_len[vc_off[c]], vc_fill[c])
                   : 1e300;   /* tiny components: keep their few edges */
    size_t nvert = 0, ndiag = 0, nrescue = 0, nlong = 0;
    for (size_t k = 0; k < nuniq; k++) {
        if (eclass[k] == 1) { ndiag++; eclass[k] = -1; continue; }
        if (eclass[k] == 2) {
            /* an offset-boundary rescue must LOOK like a vertical: pitch
             * length AND perpendicular to the row tangent (a boundary
             * diagonal shares the big dU_ref but carries a tangential
             * component -- feeding it target 0 poisons the solve) */
            int32_t sa = (int32_t)(ekey[k] >> 31);
            int32_t sb = (int32_t)(ekey[k] & 0x7fffffff);
            double ev[3] = { S->smp[sb].p[0] - S->smp[sa].p[0],
                             S->smp[sb].p[1] - S->smp[sa].p[1],
                             S->smp[sb].p[2] - S->smp[sa].p[2] };
            double along = fabs(v3dot(ev, S->smp[sa].tau));
            if (elen[k] <= vc_gate[ecomp[k]] &&
                along <= 0.15 * (elen[k] > 1e-9 ? elen[k] : 1e-9)) {
                nrescue++;
                eclass[k] = 0;
            } else {
                nlong++;
                eclass[k] = -1;
            }
            continue;
        }
        nvert++;
    }

    PairSet Q;
    memset(&Q, 0, sizeof Q);
    Q.pairs = RIB_ALLOC_ARRAY(arena, Pair, nvert + nrescue + 2 * n + 1);
    for (size_t k = 0; k < nuniq; k++) {
        if (eclass[k] != 0) continue;
        Pair *pr = &Q.pairs[Q.n_pairs++];
        pr->a = (int32_t)(ekey[k] >> 31);
        pr->b = (int32_t)(ekey[k] & 0x7fffffff);
        pr->k = 0;
        pr->d = 0.0;
        pr->w = RIB_PCONS_LAMBDA_V / (elen[k] > 0.25 ? elen[k] : 0.25);
        pr->like = 1.0;
    }
    size_t nvert_pairs = Q.n_pairs;

    /* per-sample face-tie coverage toward the previous/next row: a sample
     * missing it sits at a fit crack, where the mesh carries no vertical
     * drift edge and only a contact tie can prevent a visible slip seam */
    uint8_t *vflags = (uint8_t *)ARENA_CALLOC(arena, n + 1, sizeof(*vflags));
    for (size_t k = 0; k < nuniq; k++) {
        int32_t sa = (int32_t)(ekey[k] >> 31);
        int32_t sb = (int32_t)(ekey[k] & 0x7fffffff);
        if (S->smp[sa].slice < S->smp[sb].slice) {
            vflags[sa] |= 1u;          /* has tie toward higher slice */
            vflags[sb] |= 2u;          /* has tie toward lower slice */
        } else if (S->smp[sa].slice > S->smp[sb].slice) {
            vflags[sa] |= 2u;
            vflags[sb] |= 1u;
        }
    }

    /* --- cross-chart contact edges: same-row in-plane hash --- */
    double cs = (double)o->match_r > 1.0 ? (double)o->match_r : 1.0;
    double cos_gate = cos((double)o->match_ang_deg * M_PI / 180.0);
    double ax[3] = { (double)o->axis_dir[0], (double)o->axis_dir[1],
                     (double)o->axis_dir[2] };
    double an = sqrt(v3dot(ax, ax));
    if (an < 1e-12) { ax[0] = 1.0; ax[1] = 0.0; ax[2] = 0.0; an = 1.0; }
    ax[0] /= an; ax[1] /= an; ax[2] /= an;
    RibPairDelta *hk = RIB_ALLOC_ARRAY(arena, RibPairDelta, n + 1);
    size_t nh = 0;
    for (size_t i = 0; i < n; i++) {
        if (S->smp[i].chain < 0) continue;
        int64_t qx = (int64_t)floor(S->smp[i].c1 / cs) + 32768;
        int64_t qy = (int64_t)floor(S->smp[i].c2 / cs) + 32768;
        hk[nh].key = ((int64_t)S->smp[i].slice << 34) | (qx << 17) | qy;
        hk[nh].delta = (double)i;
        nh++;
    }
    qsort(hk, nh, sizeof(*hk), cmp_rib_pair_delta);
    size_t ncontact = 0, gated_out = 0;
    size_t contact_report = 0;
    for (size_t hi_ = 0; hi_ < nh; hi_++) {
        size_t i = (size_t)hk[hi_].delta;
        const Sample *si = &S->smp[i];
        int32_t ci = S->chn[si->chain].mesh_comp;
        int accepted = 0;
        int64_t base = hk[hi_].key;
        /* dslice 0 = classic same-row cross-chart continuation; dslice 1 =
         * the NEXT row bin, which is the only way to tie a horizontal fit
         * crack (no faces there, so no vertical drift edge exists; the
         * same-row hash can never pair v-adjacent samples) */
        for (int dslice = 0; dslice <= 1 && accepted < 2; dslice++) {
        for (int dx = -1; dx <= 1 && accepted < 2; dx++) {
            for (int dy = -1; dy <= 1 && accepted < 2; dy++) {
                int64_t want = base + ((int64_t)dslice << 34) +
                               ((int64_t)dx << 17) + dy;
                size_t lo = 0, hi2 = nh;
                while (lo < hi2) {
                    size_t mid = (lo + hi2) / 2;
                    if (hk[mid].key < want) lo = mid + 1; else hi2 = mid;
                }
                for (size_t e = lo; e < nh && hk[e].key == want &&
                                    accepted < 2; e++) {
                    size_t j = (size_t)hk[e].delta;
                    if (dslice == 0 && j <= i) continue;
                    if (j == i) continue;
                    const Sample *sj = &S->smp[j];
                    if (S->chn[sj->chain].mesh_comp == ci) {
                        /* same component: only a true crack qualifies --
                         * BOTH facing sides must lack any face tie (an
                         * ordinary interior sample is already held by its
                         * vertical drift edges) */
                        if (dslice == 0) continue;
                        if ((vflags[i] & 1u) || (vflags[j] & 2u)) continue;
                    }
                    double dp[3] = { sj->p[0]-si->p[0], sj->p[1]-si->p[1],
                                     sj->p[2]-si->p[2] };
                    double dist = sqrt(v3dot(dp, dp));
                    if (dist > (double)o->match_r) continue;
                    double tdot = fabs(v3dot(si->tau, sj->tau));
                    double along = v3dot(dp, si->tau);
                    double axial = v3dot(dp, ax);
                    double perp2 = v3dot(dp, dp) - along*along - axial*axial;
                    double perp = perp2 > 0.0 ? sqrt(perp2) : 0.0;
                    int ply = perp > RIB_PCONS_PERP_MAX ||
                              tdot < cos_gate ||
                              (o->reference_phi_authoritative &&
                               isfinite(si->phi) && isfinite(sj->phi) &&
                               fabs(si->phi - sj->phi) > RIB_CONT_DPHI_MAX);
                    if (ply) {
                        gated_out++;
                        if (contact_report < 8) {
                            fprintf(stderr,
                                    "    contact gated out at u=%.0f z=%.0f "
                                    "(perp=%.2f tdot=%.2f)\n",
                                    si->u, si->p[0], perp, tdot);
                            contact_report++;
                        }
                        continue;
                    }
                    Pair *pr = &Q.pairs[Q.n_pairs++];
                    pr->a = (int32_t)i;
                    pr->b = (int32_t)j;
                    pr->k = 0;
                    pr->d = 0.0;
                    pr->w = RIB_PCONS_LAMBDA_C / (dist > 1.0 ? dist : 1.0);
                    pr->like = 1.0;
                    ncontact++;
                    accepted++;
                }
            }
        }
        }
    }

    /* --- pins: one per connected component of the solve graph --- */
    UnionFind uf = UF_new(arena, (int32_t)S->n_chn);
    for (size_t p = 0; p < Q.n_pairs; p++)
        uf_union(&uf, S->smp[Q.pairs[p].a].chain,
                 S->smp[Q.pairs[p].b].chain);
    int32_t *root_best = RIB_ALLOC_ARRAY(arena, int32_t, S->n_chn + 1);
    for (size_t c = 0; c < S->n_chn; c++) root_best[c] = -1;
    for (size_t c = 0; c < S->n_chn; c++) {
        int32_t r = uf_find(&uf, (int32_t)c);
        if (root_best[r] < 0 ||
            S->chn[c].count > S->chn[root_best[r]].count)
            root_best[r] = (int32_t)c;
    }
    int32_t *pin_idx = RIB_ALLOC_ARRAY(arena, int32_t, S->n_chn + 1);
    double *pin_val = RIB_ALLOC_ARRAY(arena, double, S->n_chn + 1);
    size_t n_pins = 0;
    for (size_t c = 0; c < S->n_chn; c++) {
        if (uf_find(&uf, (int32_t)c) != (int32_t)c) continue;
        int32_t pc = root_best[c];
        if (pc < 0) continue;
        int32_t ps = S->chn[pc].first;
        pin_idx[n_pins] = ps;
        pin_val[n_pins] = S->smp[ps].u;
        n_pins++;
    }

    double *u = RIB_ALLOC_ARRAY(arena, double, n + 1);
    for (size_t i = 0; i < n; i++) u[i] = S->smp[i].u;
    double *u0 = RIB_ALLOC_ARRAY(arena, double, n + 1);
    memcpy(u0, u, n * sizeof(*u0));

    fprintf(stderr,
            "  consistency solve: %zu row edges (chains), %zu vertical "
            "(+%zu offset-boundary rescued; %zu diagonals and %zu long "
            "excluded), %zu contacts (%zu gated out), %zu pins\n",
            n - S->n_chn, nvert_pairs, nrescue, ndiag, nlong,
            ncontact, gated_out, n_pins);
    rib_pcons_metrics(S, &Q, nvert_pairs, u, arena, "BEFORE");

    /* arming guard: refuse an unplaced frame (see RIB_PCONS_ENTRY_DRIFT_MAX) */
    if (ncontact > 0) {
        Arena_Mark gmark = Arena_save(arena);
        double *cd = RIB_ALLOC_ARRAY(arena, double, ncontact + 1);
        size_t ncd = 0;
        for (size_t p = nvert_pairs; p < Q.n_pairs; p++)
            cd[ncd++] = fabs(u[Q.pairs[p].a] - u[Q.pairs[p].b]);
        qsort(cd, ncd, sizeof(*cd), cmp_dbl);
        {
            double p90 = cd[(ncd * 9) / 10];
            Arena_restore(arena, gmark);
            if (p90 > RIB_PCONS_ENTRY_DRIFT_MAX) {
                fprintf(stderr,
                        "  consistency solve: REFUSED -- contact ties enter "
                        "at p90=%.0f vox (max %.0f): the frame is UNPLACED "
                        "(an atlas packing?).  This solve refines placed "
                        "frames; it must not do placement through stiff "
                        "long-range ties.  Provide the phase sidecar so "
                        "island placement can run (ribbon_phase.f32 beside "
                        "the input; ribbon_sections carves it), or supply a "
                        "coherent carried frame.\n",
                        p90, RIB_PCONS_ENTRY_DRIFT_MAX);
                Arena_restore(arena, mark);
                return -1;
            }
        }
    }

    /* --- factor-once grid preconditioner + contact-carrying PCG (IRLS) ---
     * minimize sum w*(u_b - u_a - d)^2.  The TAUCS factorization sees ONLY
     * the banded grid part (row metric + vertical drift at unit Huber +
     * pins): its topology AND values are round-constant, so ONE factor
     * serves every IRLS round.  Contacts -- and each round's Huber
     * reweights -- live in the matrix-free operator that PCG applies; they
     * cost iterations, never fill.  This keeps the factor banded no matter
     * how far contact ties reach (the former full-matrix factor exploded
     * when ties spanned an atlas packing).  The rhs is round-constant too:
     * every d=0 pair contributes nothing, and chain/pin targets never
     * reweight.  If the grid factor cannot be built, fall back to the exact
     * per-round assembly. */
    size_t ncq = S->n_chn;
    size_t nchain_edges = n - ncq;
    size_t max_nt = 3u * (nchain_edges + Q.n_pairs) + n_pins + n + 8u;
    if (n > (size_t)INT32_MAX || max_nt > (size_t)INT32_MAX) {
        fprintf(stderr, "  consistency solve: system too large (n=%zu "
                "nt=%zu)\n", n, max_nt);
        Arena_restore(arena, mark);
        return -1;
    }
    int *t_row = RIB_ALLOC_ARRAY(arena, int, max_nt + 1);
    int *t_col = RIB_ALLOC_ARRAY(arena, int, max_nt + 1);
    double *t_val = RIB_ALLOC_ARRAY(arena, double, max_nt + 1);
    double *rhs = RIB_ALLOC_ARRAY(arena, double, n + 1);
    double *hub = RIB_ALLOC_ARRAY(arena, double, Q.n_pairs + 1);
    for (size_t p = 0; p < Q.n_pairs; p++) hub[p] = 1.0;

    /* round-constant rhs: chain metric targets + pins */
    memset(rhs, 0, n * sizeof(*rhs));
    for (size_t c = 0; c < ncq; c++) {
        int32_t f2 = S->chn[c].first, cn2 = S->chn[c].count;
        for (int32_t i = 1; i < cn2; i++) {
            int a = f2 + i - 1, b = f2 + i;
            double l = S->smp[b].s - S->smp[a].s;
            if (l < 1e-9) l = 1e-9;
            double w = 1.0 / l;
            rhs[a] -= w * l;
            rhs[b] += w * l;
        }
    }
    for (size_t pi2 = 0; pi2 < n_pins; pi2++)
        rhs[pin_idx[pi2]] += pin_val[pi2];
    double bnorm = 0.0;
    for (size_t i = 0; i < n; i++) bnorm += rhs[i] * rhs[i];
    bnorm = sqrt(bnorm);
    if (bnorm < 1e-30) bnorm = 1.0;

    /* grid-only factor (chains + verticals + pins) -- only where the
     * connected full factor would hurt (see RIB_PCONS_PCG_MIN_N) */
    SparseFactor_T grid_f = NULL;
    if (n > (size_t)RIB_PCONS_PCG_MIN_N) {
        double tf0 = ves_clock_sec();
        size_t nt = 0;
        for (size_t c = 0; c < ncq; c++) {
            int32_t f2 = S->chn[c].first, cn2 = S->chn[c].count;
            for (int32_t i = 1; i < cn2; i++) {
                int a = f2 + i - 1, b = f2 + i;
                double l = S->smp[b].s - S->smp[a].s;
                if (l < 1e-9) l = 1e-9;
                double w = 1.0 / l;
                t_row[nt] = a; t_col[nt] = a; t_val[nt] = w; nt++;
                t_row[nt] = b; t_col[nt] = b; t_val[nt] = w; nt++;
                t_row[nt] = b; t_col[nt] = a; t_val[nt] = -w; nt++;
            }
        }
        for (size_t p = 0; p < nvert_pairs; p++) {
            int a = Q.pairs[p].a, b = Q.pairs[p].b;
            double w = Q.pairs[p].w;
            int lo = a < b ? a : b, hi2 = a < b ? b : a;
            if (w <= 0.0) continue;
            t_row[nt] = a; t_col[nt] = a; t_val[nt] = w; nt++;
            t_row[nt] = b; t_col[nt] = b; t_val[nt] = w; nt++;
            t_row[nt] = hi2; t_col[nt] = lo; t_val[nt] = -w; nt++;
        }
        for (size_t pi2 = 0; pi2 < n_pins; pi2++) {
            t_row[nt] = pin_idx[pi2];
            t_col[nt] = pin_idx[pi2];
            t_val[nt] = 1.0;
            nt++;
        }
        /* PRECONDITIONER-ONLY additions (the operator stays the historical
         * exact system): (a) the contacts' DIAGONAL lumping, so fragments
         * whose only global restraint is a contact are restrained in M at
         * the same magnitude A restrains them (their off-diagonals ride
         * PCG; no fill); (b) a feather ridge for fragments with no
         * restraint at all, so the Cholesky never zero-pivots. */
        for (size_t p = nvert_pairs; p < Q.n_pairs; p++) {
            int a = Q.pairs[p].a, b = Q.pairs[p].b;
            double w = Q.pairs[p].w;
            if (w <= 0.0) continue;
            t_row[nt] = a; t_col[nt] = a; t_val[nt] = w; nt++;
            t_row[nt] = b; t_col[nt] = b; t_val[nt] = w; nt++;
        }
        for (size_t i = 0; i < n; i++) {
            t_row[nt] = (int)i;
            t_col[nt] = (int)i;
            t_val[nt] = RIB_PCONS_GRID_RIDGE;
            nt++;
        }
        if (Sparse_factor_spd((int)n, (int)nt, t_row, t_col, t_val,
                              &grid_f) == 0)
            fprintf(stderr,
                    "  consistency grid factor: n=%zu nnz=%zu "
                    "factor_sec=%.1f (one factor for all IRLS rounds; "
                    "contacts ride PCG)\n",
                    n, nt, ves_clock_sec() - tf0);
        else {
            grid_f = NULL;
            fprintf(stderr,
                    "  consistency grid factor unavailable; exact "
                    "per-round fallback\n");
        }
    }

    double *pcg_r = RIB_ALLOC_ARRAY(arena, double, n + 1);
    double *pcg_z = RIB_ALLOC_ARRAY(arena, double, n + 1);
    double *pcg_p = RIB_ALLOC_ARRAY(arena, double, n + 1);
    double *pcg_ap = RIB_ALLOC_ARRAY(arena, double, n + 1);

    for (int round = 0; round < RIB_PCONS_IRLS_ROUNDS; round++) {
        double tr0 = ves_clock_sec();
        int pcg_iters = 0;
        double relres = 0.0;
        if (grid_f != NULL) {
            /* PCG on the full operator (chains + all pairs*hub + pins),
             * preconditioned by the grid factor, warm-started from u. */
#define RIB_PCONS_APPLY(SRC, DST)                                         \
            do {                                                          \
                memset((DST), 0, n * sizeof(double));                     \
                for (size_t c2 = 0; c2 < ncq; c2++) {                     \
                    int32_t f3 = S->chn[c2].first;                        \
                    int32_t cn3 = S->chn[c2].count;                       \
                    for (int32_t i2 = 1; i2 < cn3; i2++) {                \
                        int a2 = f3 + i2 - 1, b2 = f3 + i2;               \
                        double l2 = S->smp[b2].s - S->smp[a2].s;          \
                        double w2 = 1.0 / (l2 < 1e-9 ? 1e-9 : l2);        \
                        double t2 = w2 * ((SRC)[a2] - (SRC)[b2]);         \
                        (DST)[a2] += t2;                                  \
                        (DST)[b2] -= t2;                                  \
                    }                                                     \
                }                                                         \
                for (size_t p2 = 0; p2 < Q.n_pairs; p2++) {               \
                    int a2 = Q.pairs[p2].a, b2 = Q.pairs[p2].b;           \
                    double w2 = Q.pairs[p2].w * hub[p2];                  \
                    double t2 = 0.0;                                      \
                    if (w2 <= 0.0) continue;                              \
                    t2 = w2 * ((SRC)[a2] - (SRC)[b2]);                    \
                    (DST)[a2] += t2;                                      \
                    (DST)[b2] -= t2;                                      \
                }                                                         \
                for (size_t p2 = 0; p2 < n_pins; p2++)                    \
                    (DST)[pin_idx[p2]] += (SRC)[pin_idx[p2]];             \
            } while (0)
            double rn = 0.0, rz = 0.0;
            RIB_PCONS_APPLY(u, pcg_ap);
            for (size_t i = 0; i < n; i++) {
                pcg_r[i] = rhs[i] - pcg_ap[i];
                rn += pcg_r[i] * pcg_r[i];
            }
            relres = sqrt(rn) / bnorm;
            if (relres > RIB_PCONS_PCG_TOL) {
                if (Sparse_factor_solve(grid_f, pcg_r, pcg_z) != 0)
                    memcpy(pcg_z, pcg_r, n * sizeof(double));
                for (size_t i = 0; i < n; i++) {
                    pcg_p[i] = pcg_z[i];
                    rz += pcg_r[i] * pcg_z[i];
                }
                for (pcg_iters = 0; pcg_iters < RIB_PCONS_PCG_MAXIT;
                     pcg_iters++) {
                    double pap = 0.0, alpha = 0.0, rn2 = 0.0, rz2 = 0.0;
                    RIB_PCONS_APPLY(pcg_p, pcg_ap);
                    for (size_t i = 0; i < n; i++)
                        pap += pcg_p[i] * pcg_ap[i];
                    if (pap <= 0.0) break;
                    alpha = rz / pap;
                    for (size_t i = 0; i < n; i++) {
                        u[i] += alpha * pcg_p[i];
                        pcg_r[i] -= alpha * pcg_ap[i];
                        rn2 += pcg_r[i] * pcg_r[i];
                    }
                    relres = sqrt(rn2) / bnorm;
                    if (relres < RIB_PCONS_PCG_TOL) {
                        pcg_iters++;
                        break;
                    }
                    if (Sparse_factor_solve(grid_f, pcg_r, pcg_z) != 0)
                        memcpy(pcg_z, pcg_r, n * sizeof(double));
                    for (size_t i = 0; i < n; i++)
                        rz2 += pcg_r[i] * pcg_z[i];
                    {
                        double beta = rz2 / rz;
                        rz = rz2;
                        for (size_t i = 0; i < n; i++)
                            pcg_p[i] = pcg_z[i] + beta * pcg_p[i];
                    }
                }
            }
#undef RIB_PCONS_APPLY
        } else {
            /* exact per-round fallback: assemble everything and factor */
            size_t nt = 0;
            for (size_t c = 0; c < ncq; c++) {
                int32_t f2 = S->chn[c].first, cn2 = S->chn[c].count;
                for (int32_t i = 1; i < cn2; i++) {
                    int a = f2 + i - 1, b = f2 + i;
                    double l = S->smp[b].s - S->smp[a].s;
                    if (l < 1e-9) l = 1e-9;
                    double w = 1.0 / l;
                    t_row[nt] = a; t_col[nt] = a; t_val[nt] = w; nt++;
                    t_row[nt] = b; t_col[nt] = b; t_val[nt] = w; nt++;
                    t_row[nt] = b; t_col[nt] = a; t_val[nt] = -w; nt++;
                }
            }
            for (size_t p = 0; p < Q.n_pairs; p++) {
                int a = Q.pairs[p].a, b = Q.pairs[p].b;
                double w = Q.pairs[p].w * hub[p];
                if (w <= 0.0) continue;
                int lo = a < b ? a : b, hi2 = a < b ? b : a;
                t_row[nt] = a; t_col[nt] = a; t_val[nt] = w; nt++;
                t_row[nt] = b; t_col[nt] = b; t_val[nt] = w; nt++;
                t_row[nt] = hi2; t_col[nt] = lo; t_val[nt] = -w; nt++;
            }
            for (size_t pi2 = 0; pi2 < n_pins; pi2++) {
                t_row[nt] = pin_idx[pi2];
                t_col[nt] = pin_idx[pi2];
                t_val[nt] = 1.0;
                nt++;
            }
            if (Sparse_solve_sym((int)n, (int)nt, t_row, t_col, t_val,
                                 rhs, u, SPARSE_SPD) != 0) {
                fprintf(stderr, "  consistency solve: TAUCS factorization "
                        "FAILED (n=%zu nnz=%zu round=%d)\n",
                        n, nt, round + 1);
                Arena_restore(arena, mark);
                return -1;
            }
        }
        /* Huber reweight from the solved residuals: a drift/contact term
         * still far off is evidence the model tied the wrong samples, so
         * its influence decays as c/|r|. */
        size_t ndown = 0;
        double wchg = 0.0;
        for (size_t p = 0; p < Q.n_pairs; p++) {
            double hc = p < nvert_pairs ? RIB_PCONS_HUBER_C
                                        : RIB_PCONS_HUBER_C_CONTACT;
            double r = fabs(u[Q.pairs[p].a] - u[Q.pairs[p].b]);
            double h = r > hc ? hc / r : 1.0;
            if (h < 0.999) ndown++;
            if (fabs(h - hub[p]) > wchg) wchg = fabs(h - hub[p]);
            hub[p] = h;
        }
        char tag[24];
        snprintf(tag, sizeof tag, "round %d", round + 1);
        fprintf(stderr,
                "  consistency IRLS round %d: n=%zu solve=%.2fs "
                "(pcg_iters=%d relres=%.2e), %zu/%zu drift pairs "
                "downweighted (max dw=%.3f)\n",
                round + 1, n, ves_clock_sec() - tr0, pcg_iters, relres,
                ndown, Q.n_pairs, wchg);
        rib_pcons_metrics(S, &Q, nvert_pairs, u, arena, tag);
        if (wchg < 1e-3) break;   /* weights stable: IRLS converged */
    }
    Sparse_factor_free(&grid_f);

    /* Placement restore: every energy term except the pins is a difference,
     * so each component's translation lands exactly where its single pin
     * sample says -- and if that sample happens to sit inside a carried
     * defect band, the exact solve parks the WHOLE component at the
     * defect's offset.  Stage 1's anchoring is the placement authority;
     * restore it per solve-graph component as a rigid translation, the L1
     * (median) estimate of (u - u_init) over the component, which keeps
     * every shape correction the solve made and is independent of pin
     * luck.  Contact-joined charts share a UF root and re-anchor as one
     * unit, so contact-driven placement survives. */
    {
        size_t *rcnt = (size_t *)ARENA_CALLOC(arena, ncq + 1, sizeof(*rcnt));
        for (size_t c = 0; c < ncq; c++) {
            size_t r = (size_t)uf_find(&uf, (int32_t)c);
            rcnt[r] += (size_t)S->chn[c].count;
        }
        size_t *roff = (size_t *)ARENA_CALLOC(arena, ncq + 2, sizeof(*roff));
        for (size_t r = 0; r < ncq; r++) roff[r + 1] = roff[r] + rcnt[r];
        double *rbuf = RIB_ALLOC_ARRAY(arena, double, roff[ncq] + 1);
        size_t *rfill = (size_t *)ARENA_CALLOC(arena, ncq + 1,
                                               sizeof(*rfill));
        for (size_t c = 0; c < ncq; c++) {
            size_t r = (size_t)uf_find(&uf, (int32_t)c);
            int32_t f2 = S->chn[c].first, cn2 = S->chn[c].count;
            for (int32_t k2 = 0; k2 < cn2; k2++)
                rbuf[roff[r] + rfill[r]++] = u[(size_t)(f2 + k2)]
                                           - u0[(size_t)(f2 + k2)];
        }
        double shift_max = 0.0;
        double *rmed = (double *)ARENA_CALLOC(arena, ncq + 1, sizeof(*rmed));
        for (size_t r = 0; r < ncq; r++) {
            if (rcnt[r] == 0) continue;
            double *seg = rbuf + roff[r];
            qsort(seg, rcnt[r], sizeof(*seg), cmp_dbl);
            double med = seg[rcnt[r] / 2];
            if (!isfinite(med)) continue;
            rmed[r] = med;
            if (fabs(med) > fabs(shift_max)) shift_max = med;
        }
        for (size_t c = 0; c < ncq; c++) {
            double med = rmed[(size_t)uf_find(&uf, (int32_t)c)];
            if (med == 0.0) continue;
            int32_t f2 = S->chn[c].first, cn2 = S->chn[c].count;
            for (int32_t k2 = 0; k2 < cn2; k2++)
                u[(size_t)(f2 + k2)] -= med;
        }
        if (fabs(shift_max) > 1e-6)
            fprintf(stderr,
                    "  consistency placement restore: max component "
                    "translation %.2f vox\n", shift_max);
    }

    /* The probabilistic winding certificate is an ORDER constraint, not just
     * another least-squares vote.  The common map u0 is monotone in carried U.
     * On confident samples, keep the metric refinement inside a small fraction
     * of one local turn (2*pi*r).  This still permits seam/metric cleanup but
     * makes an adjacent-ply exchange impossible unless the winding posterior
     * explicitly abstained.  Bounds scale continuously with uncertainty and
     * remain meaningful on cropped scrolls that do not contain the axis. */
    double *cert_lo = NULL, *cert_hi = NULL;
    if (sample_uconf != NULL) {
        enum { RIB_CERT_MIN_CONF_PCT = 75 };
        const double min_conf = (double)RIB_CERT_MIN_CONF_PCT / 100.0;
        cert_lo = RIB_ALLOC_ARRAY(arena, double, n + 1);
        cert_hi = RIB_ALLOC_ARRAY(arena, double, n + 1);
        double *drift = RIB_ALLOC_ARRAY(arena, double, n + 1);
        double *bound = RIB_ALLOC_ARRAY(arena, double, n + 1);
        size_t nc = 0, nclip = 0;
        for (size_t i = 0; i < n; i++) {
            double conf = (double)sample_uconf[i];
            cert_lo[i] = -DBL_MAX;
            cert_hi[i] = DBL_MAX;
            if (!isfinite(conf) || conf < min_conf) continue;
            if (conf > 1.0) conf = 1.0;
            double turn_radius = S->smp[i].r;
            double radius_floor = 4.0 * (double)o->wrap_spacing;
            if (!(radius_floor > 0.0)) radius_floor = 4.0 * 9.5;
            if (turn_radius < radius_floor) turn_radius = radius_floor;
            /* 1% of a turn at confidence 1, rising to 3% at the 0.75
             * arming threshold. */
            double frac = 0.01 + 0.08 * (1.0 - conf);
            double b = frac * (2.0 * M_PI) * turn_radius;
            double bmin = 2.0 * (double)o->sample_h;
            if (bmin < 2.0) bmin = 2.0;
            if (b < bmin) b = bmin;
            cert_lo[i] = u0[i] - b;
            cert_hi[i] = u0[i] + b;
            drift[nc] = fabs(u[i] - u0[i]);
            bound[nc] = b;
            nc++;
            if (u[i] < cert_lo[i]) { u[i] = cert_lo[i]; nclip++; }
            else if (u[i] > cert_hi[i]) { u[i] = cert_hi[i]; nclip++; }
        }
        if (nc != 0) {
            out->certificate_samples += nc;
            out->certificate_clamps += nclip;
            out->certificate_drift_p95 = select_p95_dbl(drift, nc);
            for (size_t i = 0; i < nc; i++)
                if (drift[i] > out->certificate_drift_max)
                    out->certificate_drift_max = drift[i];
            out->certificate_bound_p50 = select_median_dbl(bound, nc);
            fprintf(stderr,
                    "  certificate trust region: %zu confident samples, "
                    "%zu clamped; |drift| p95/max=%.2f/%.2f vox, "
                    "median bound=%.2f vox\n",
                    nc, nclip, out->certificate_drift_p95,
                    out->certificate_drift_max,
                    out->certificate_bound_p50);
        }
    }

    /* Per-row monotonicity plus the certificate boxes.  PAVA supplies the L2
     * monotone projection.  Prefix-lower/suffix-upper envelopes then make the
     * individual boxes monotone-feasible; a final forward pass cannot exceed
     * those nondecreasing upper bounds. */
    size_t mono_bad = 0, mono_chains = 0;
    int32_t max_chain = 2;
    for (size_t c = 0; c < S->n_chn; c++)
        if (S->chn[c].count > max_chain) max_chain = S->chn[c].count;
    double *lb = (double *)ARENA_CALLOC(arena, (size_t)max_chain,
                                        sizeof(*lb));
    double *scr = RIB_ALLOC_ARRAY(arena, double, (size_t)max_chain * 3u);
    for (size_t c = 0; c < S->n_chn; c++) {
        int32_t f = S->chn[c].first, cn = S->chn[c].count;
        int bad = 0;
        for (int32_t i = 1; i < cn; i++)
            if (u[f+i] < u[f+i-1]) { bad = 1; mono_bad++; }
        if (!bad) continue;
        mono_chains++;
        (void)pava_chain(&u[f], lb, cn, scr);
    }
    if (cert_lo != NULL) {
        for (size_t c = 0; c < S->n_chn; c++) {
            int32_t f = S->chn[c].first, cn = S->chn[c].count;
            for (int32_t i = 1; i < cn; i++)
                if (cert_lo[f+i] < cert_lo[f+i-1])
                    cert_lo[f+i] = cert_lo[f+i-1];
            for (int32_t i = cn - 1; i-- > 0; )
                if (cert_hi[f+i] > cert_hi[f+i+1])
                    cert_hi[f+i] = cert_hi[f+i+1];
            for (int32_t i = 0; i < cn; i++) {
                if (cert_lo[f+i] > cert_hi[f+i]) {
                    /* u0 is monotone, so this is only reachable through
                     * numerical noise in an extreme radius envelope. */
                    double mid = 0.5 * (cert_lo[f+i] + cert_hi[f+i]);
                    cert_lo[f+i] = cert_hi[f+i] = mid;
                }
                if (u[f+i] < cert_lo[f+i]) u[f+i] = cert_lo[f+i];
                if (u[f+i] > cert_hi[f+i]) u[f+i] = cert_hi[f+i];
                if (i > 0 && u[f+i] < u[f+i-1]) u[f+i] = u[f+i-1];
            }
        }
    }
    if (mono_bad != 0)
        fprintf(stderr,
                "  consistency monotone repair: %zu edge(s) in %zu chain(s)\n",
                mono_bad, mono_chains);
    out->mono_repairs += mono_bad;

    rib_pcons_metrics(S, &Q, nvert_pairs, u, arena, "AFTER ");
    for (size_t i = 0; i < n; i++) S->smp[i].u = u[i];
    fprintf(stderr, "  consistency solve total: %.2fs\n",
            ves_clock_sec() - t0);
    Arena_restore(arena, mark);
    return 0;
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
    if (nv > (size_t)INT32_MAX)
        return -1;

    RibbonOpts def, resolved;
    if (opts == NULL) { RibbonOpts_default(&def); opts = &def; }
    if (opts->reference_material_island != NULL) {
        if (opts->reference_island != NULL) {
            fprintf(stderr,
                    "ribbon: supply reference_material_island, not both "
                    "material and legacy island fields\n");
            return -1;
        }
        resolved = *opts;
        resolved.reference_island = opts->reference_material_island;
        resolved.reference_island_count =
            opts->reference_material_island_count;
        opts = &resolved;
    }
    if (opts->solve_level_prefix != NULL &&
        (!opts->solve_amg || !opts->fit_ribbon ||
         (opts->reference_u != NULL && !opts->solve_reference_u))) {
        fprintf(stderr,
                "ribbon: solve-level dumps require the FMG Stage-C solve and "
                "fitted ribbon output\n");
        return -1;
    }
    RibGMGSnapshots level_capture;
    memset(&level_capture, 0, sizeof level_capture);
    RibGMGSnapshots *level_capture_ptr =
        opts->solve_level_prefix != NULL ? &level_capture : NULL;
    double stage_t0 = ves_clock_sec();

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

    /* Exact source-chart provenance.  Spatially contacting fragments may be
     * only a voxel apart, so Stage D must know which source component produced
     * each slice sample instead of trying to infer identity from distance. */
    UnionFind mesh_uf = UF_new(arena, (int32_t)nv);
    for (size_t f = 0; f < nf; f++) {
        int32_t a = faces[f*3], b = faces[f*3+1], c = faces[f*3+2];
        if (a < 0 || b < 0 || c < 0 ||
            (size_t)a >= nv || (size_t)b >= nv || (size_t)c >= nv)
            return -1;
        uf_union(&mesh_uf, a, b);
        uf_union(&mesh_uf, b, c);
    }
    int32_t *mesh_root_label = RIB_ALLOC_ARRAY(arena, int32_t, nv);
    int32_t *vertex_mesh_comp = RIB_ALLOC_ARRAY(arena, int32_t, nv);
    for (size_t i = 0; i < nv; i++) mesh_root_label[i] = -1;
    size_t n_mesh_comp = 0;
    for (size_t i = 0; i < nv; i++) {
        int32_t root = uf_find(&mesh_uf, (int32_t)i);
        if (mesh_root_label[root] < 0)
            mesh_root_label[root] = (int32_t)n_mesh_comp++;
        vertex_mesh_comp[i] = mesh_root_label[root];
    }
    fprintf(stderr, "  ribbon stage mesh-components: %.2fs\n",
            ves_clock_sec() - stage_t0);
    int32_t *mesh_comp_island = RIB_ALLOC_ARRAY(
        arena, int32_t, n_mesh_comp ? n_mesh_comp : 1);
    int32_t *mesh_comp_chart = RIB_ALLOC_ARRAY(
        arena, int32_t, n_mesh_comp ? n_mesh_comp : 1);
    size_t reference_island_count = opts->reference_island != NULL
                                  ? opts->reference_island_count : 1;
    if (reference_island_count == 0) {
        fprintf(stderr,
                "ribbon: reference_island_count is required with "
                "reference_island\n");
        return -1;
    }
    for (size_t c = 0; c < n_mesh_comp; c++) {
        mesh_comp_island[c] = -1;
        mesh_comp_chart[c] = -1;
    }
    for (size_t i = 0; i < nv; i++) {
        int32_t mc = vertex_mesh_comp[i];
        int32_t island = opts->reference_island != NULL
                       ? opts->reference_island[i] : 0;
        int32_t source_chart = opts->reference_chart != NULL
                             ? opts->reference_chart[i] : mc;
        if (island < 0 || (size_t)island >= reference_island_count) {
            fprintf(stderr,
                    "ribbon: invalid material-island identity %d at vertex %zu "
                    "(%zu material islands)\n",
                    island, i, reference_island_count);
            return -1;
        }
        if (source_chart < 0 ||
            (opts->reference_chart != NULL &&
             (size_t)source_chart >= opts->reference_chart_count)) {
            fprintf(stderr,
                    "ribbon: invalid source-chart identity %d at vertex %zu\n",
                    source_chart, i);
            return -1;
        }
        if (mesh_comp_island[mc] < 0)
            mesh_comp_island[mc] = island;
        else if (mesh_comp_island[mc] != island) {
            fprintf(stderr,
                    "ribbon: source component %d spans material islands %d/%d\n",
                    mc, mesh_comp_island[mc], island);
            return -1;
        }
        if (mesh_comp_chart[mc] < 0)
            mesh_comp_chart[mc] = source_chart;
        else if (mesh_comp_chart[mc] != source_chart) {
            fprintf(stderr,
                    "ribbon: input component %d spans source charts %d/%d\n",
                    mc, mesh_comp_chart[mc], source_chart);
            return -1;
        }
    }
    if (opts->reference_chart != NULL) {
        out->source_chart_island = RIB_ALLOC_ARRAY(
            arena, int32_t, opts->reference_chart_count);
        out->source_chart_island_count = opts->reference_chart_count;
        for (size_t c = 0; c < out->source_chart_island_count; c++)
            out->source_chart_island[c] = -1;
        for (size_t i = 0; i < nv; i++) {
            int32_t source_chart = opts->reference_chart[i];
            int32_t island = opts->reference_island != NULL
                           ? opts->reference_island[i] : 0;
            int32_t *slot = &out->source_chart_island[source_chart];
            if (*slot < 0) *slot = island;
            else if (*slot != island) {
                fprintf(stderr,
                        "ribbon: source chart %d spans material islands %d/%d\n",
                        source_chart, *slot, island);
                return -1;
            }
        }
    }
    /* per-vertex axial t */
    stage_t0 = ves_clock_sec();
    double *t = RIB_ALLOC_ARRAY(arena, double, nv);
    double tmin = 1e300, tmax = -1e300;
    for (size_t i = 0; i < nv; i++) {
        double d0 = (double)verts[i*3+0] - (double)opts->axis_point[0];
        double d1 = (double)verts[i*3+1] - (double)opts->axis_point[1];
        double d2 = (double)verts[i*3+2] - (double)opts->axis_point[2];
        t[i] = d0*ad[0] + d1*ad[1] + d2*ad[2];
        if (t[i] < tmin) tmin = t[i];
        if (t[i] > tmax) tmax = t[i];
    }

    /* A: slice (+ bridge cut).  A topology-preserving quadribbon solve uses
     * the input's own param-v rows -- geometric planes cannot align with
     * several phase-offset component lattices at once. */
    SliceSet S;
    int32_t *row_vertex_sample = NULL;
    float *row_sample_uref = NULL;
    float *row_sample_uconf = NULL;
    double row_du_lat = 1.0;
    if ((opts->preserve_input_rows || opts->metric_project_only) &&
        opts->reference_u != NULL &&
        opts->solve_reference_u) {
        if (rib_rows_from_quadribbon(arena, verts, nv, vertex_mesh_comp,
                                     mesh_comp_island, mesh_comp_chart, t,
                                     e1, e2, opts->axis_point, opts,
                                     &S, &row_vertex_sample,
                                     &row_sample_uref, &row_sample_uconf,
                                     &row_du_lat) != 0) {
            memset(out, 0, sizeof(*out));
            return -1;
        }
        S.tmin = tmin;
        S.tmax = tmax;
    } else if (slice_mesh(arena, verts, nv, faces, nf, vertex_mesh_comp,
                   mesh_comp_island, mesh_comp_chart,
                   t, tmin, tmax,
                   e1, e2, opts->axis_point, opts, &S) != 0) {
        memset(out, 0, sizeof(*out));
        return -1;
    }
    fprintf(stderr, "  ribbon stage slice: %.2fs\n",
            ves_clock_sec() - stage_t0);
    out->n_slices       = S.n_slices_hit;
    out->n_chains       = (int)S.n_chn;
    out->n_closed       = S.n_closed;
    out->n_multi_slices = S.n_multi;
    out->n_samples      = S.n_smp;
    out->bridge_cuts    = S.bridge_cuts;
    out->v_span         = tmax - tmin;

    /* Metric-projection mode, step 1: the COMMON METRIC MAP.  The carried
     * coordinate is transversely coherent but not physical arclength (3%
     * compressed at 4x, ~50% at 10x), and the compression varies along the
     * sheet.  Estimate the cross-V robust rate du_ref/ds as a function of
     * u_ref per material island, integrate its reciprocal, and re-seed every
     * sample with the same monotone F(u_ref).  This is a measured
     * reparameterization of the carried frame -- no unknowns are solved and
     * every chain sees the identical map, so transverse coherence survives. */
    if (opts->reference_u != NULL && opts->solve_reference_u &&
        opts->metric_project_only) {
        /* ONE map over the whole carried axis, never per island: the carried
         * atlas gives islands disjoint u_ref intervals, so per-bin rates are
         * island-local anyway, while a single cumulative integral preserves
         * the carried frame's relative order and adjacency -- essential when
         * the carried registration itself is the placement authority. */
        Arena_Mark fmark = Arena_save(arena);
        enum { RIB_FMAP_BIN = 32 };
        double lo = 1e300, hi = -1e300;
        size_t nrate = 0;
        for (size_t c = 0; c < S.n_chn; c++) {
            int32_t f = S.chn[c].first, cn = S.chn[c].count;
            for (int32_t i = 0; i < cn; i++) {
                double uu = S.smp[f+i].u;
                if (!isfinite(uu)) continue;
                if (uu < lo) lo = uu;
                if (uu > hi) hi = uu;
                if (i > 0) nrate++;
            }
        }
        if (hi > lo && nrate >= 64) {
            size_t nb = (size_t)((hi - lo) / (double)RIB_FMAP_BIN) + 1;
            if (nb > 65536) nb = 65536;
            if (nb < 2) nb = 2;
            double binw = (hi - lo) / (double)nb;
            size_t *bin_count = (size_t *)ARENA_CALLOC(
                arena, nb, sizeof(*bin_count));
            size_t *bin_fill = (size_t *)ARENA_CALLOC(
                arena, nb, sizeof(*bin_fill));
            for (size_t c = 0; c < S.n_chn; c++) {
                int32_t f = S.chn[c].first, cn = S.chn[c].count;
                for (int32_t i = 1; i < cn; i++) {
                    double du_ = S.smp[f+i].u - S.smp[f+i-1].u;
                    double ds_ = S.smp[f+i].s - S.smp[f+i-1].s;
                    if (!isfinite(du_) || ds_ < 1e-6) continue;
                    double mid = 0.5 * (S.smp[f+i].u + S.smp[f+i-1].u);
                    long bi = (long)((mid - lo) / binw);
                    if (bi < 0) bi = 0;
                    if (bi >= (long)nb) bi = (long)nb - 1;
                    bin_count[bi]++;
                }
            }
            size_t *bin_off = (size_t *)ARENA_CALLOC(
                arena, nb + 1, sizeof(*bin_off));
            for (size_t bi = 0; bi < nb; bi++)
                bin_off[bi + 1] = bin_off[bi] + bin_count[bi];
            double *rates = RIB_ALLOC_ARRAY(arena, double, bin_off[nb] + 1);
            for (size_t c = 0; c < S.n_chn; c++) {
                int32_t f = S.chn[c].first, cn = S.chn[c].count;
                for (int32_t i = 1; i < cn; i++) {
                    double du_ = S.smp[f+i].u - S.smp[f+i-1].u;
                    double ds_ = S.smp[f+i].s - S.smp[f+i-1].s;
                    if (!isfinite(du_) || ds_ < 1e-6) continue;
                    double mid = 0.5 * (S.smp[f+i].u + S.smp[f+i-1].u);
                    long bi = (long)((mid - lo) / binw);
                    if (bi < 0) bi = 0;
                    if (bi >= (long)nb) bi = (long)nb - 1;
                    double rate = du_ / ds_;
                    if (rate < 1e-3) rate = 1e-3;
                    if (rate > 1e3) rate = 1e3;
                    rates[bin_off[bi] + bin_fill[bi]++] = rate;
                }
            }
            /* robust per-bin rate; empty bins inherit the last seen value */
            double *bin_rate = RIB_ALLOC_ARRAY(arena, double, nb);
            double carry = 1.0, rate_sum = 0.0;
            size_t rate_n = 0;
            for (size_t bi = 0; bi < nb; bi++) {
                if (bin_fill[bi] >= 4) {
                    carry = select_median_dbl(&rates[bin_off[bi]],
                                              bin_fill[bi]);
                    rate_sum += carry;
                    rate_n++;
                }
                bin_rate[bi] = carry;
            }
            for (size_t bi = nb; bi-- > 0; ) {   /* leading empties: backfill */
                if (bin_fill[bi] >= 4) carry = bin_rate[bi];
                else bin_rate[bi] = carry;
            }
            /* F(u): anchored at lo, piecewise-linear over bins */
            double *F = RIB_ALLOC_ARRAY(arena, double, nb + 1);
            F[0] = lo;
            for (size_t bi = 0; bi < nb; bi++)
                F[bi + 1] = F[bi] + binw / bin_rate[bi];
            for (size_t c = 0; c < S.n_chn; c++) {
                int32_t f = S.chn[c].first, cn = S.chn[c].count;
                for (int32_t i = 0; i < cn; i++) {
                    double uu = S.smp[f+i].u;
                    if (!isfinite(uu)) continue;
                    double bt = (uu - lo) / binw;
                    long bi = (long)bt;
                    if (bi < 0) { bi = 0; bt = 0.0; }
                    if (bi >= (long)nb) { bi = (long)nb - 1; bt = (double)nb; }
                    S.smp[f+i].u = F[bi] + (bt - (double)bi) *
                                   (F[bi + 1] - F[bi]);
                }
            }
            if (rate_n > 0)
                fprintf(stderr,
                        "  common metric map: u_ref=[%.0f,%.0f] bins=%zu "
                        "mean rate=%.3f -> span %.0f\n",
                        lo, hi, nb, rate_sum / (double)rate_n, F[nb] - F[0]);
        }
        Arena_restore(arena, fmark);
    }

    /* Metric-projection mode, step 2: a chain whose traversal runs against the
     * carried coordinate would be silently mirrored by u = gauge + s.
     * Re-orient it to the carried frame NOW, before any structure (pairs,
     * winding) records sample indices.  A measured re-orientation, not a
     * solve. */
    if (opts->reference_u != NULL && opts->solve_reference_u &&
        opts->metric_project_only) {
        size_t reversed = 0;
        for (size_t c = 0; c < S.n_chn; c++) {
            int32_t f = S.chn[c].first, cn = S.chn[c].count;
            if (cn < 2) continue;
            double mu = 0.0, ms = 0.0, cov = 0.0;
            for (int32_t i = 0; i < cn; i++) {
                mu += S.smp[f+i].u;
                ms += S.smp[f+i].s;
            }
            mu /= (double)cn;
            ms /= (double)cn;
            for (int32_t i = 0; i < cn; i++)
                cov += (S.smp[f+i].u - mu) * (S.smp[f+i].s - ms);
            if (cov >= 0.0) continue;
            for (int32_t i = 0; i < cn / 2; i++) {
                Sample tmp = S.smp[f + i];
                S.smp[f + i] = S.smp[f + cn - 1 - i];
                S.smp[f + cn - 1 - i] = tmp;
            }
            S.smp[f].s = 0.0;
            for (int32_t i = 1; i < cn; i++) {
                double dz = S.smp[f+i].p[0] - S.smp[f+i-1].p[0];
                double dy = S.smp[f+i].p[1] - S.smp[f+i-1].p[1];
                double dx = S.smp[f+i].p[2] - S.smp[f+i-1].p[2];
                S.smp[f+i].s = S.smp[f+i-1].s +
                               sqrt(dz*dz + dy*dy + dx*dx);
            }
            for (int32_t i = 0; i < cn; i++) {
                S.smp[f+i].tau[0] = -S.smp[f+i].tau[0];
                S.smp[f+i].tau[1] = -S.smp[f+i].tau[1];
                S.smp[f+i].tau[2] = -S.smp[f+i].tau[2];
            }
            reversed++;
        }
        if (reversed != 0)
            fprintf(stderr,
                    "  metric projection: %zu chain(s) re-oriented to the "
                    "carried frame\n", reversed);
    }

    /* A caller-provided material coordinate makes Stages B--D redundant.
     * The slice intersections already carry both measured geometry and the
     * authoritative U/phase fields.  Project each slice chain onto exact
     * measured arclength around its robust graph-derived gauge, reconcile those
     * gauges with the anchored chain graph above, then run the same claimant
     * selection and fitted-ribbon Stage E used by the full StrokeStrip solve. */
    if (opts->reference_u != NULL && !opts->solve_reference_u) {
        Arena_Mark fast_mark = Arena_save(arena);
        int32_t max_chain = 2;
        double umin = 1e300, umax = -1e300;
        double pmin = 1e300, pmax = -1e300;
        size_t moved = 0;
        for (size_t c = 0; c < S.n_chn; c++) {
            if (S.chn[c].count > max_chain) max_chain = S.chn[c].count;
            S.chn[c].solve_comp = S.chn[c].mesh_comp;
            S.chn[c].group = S.chn[c].mesh_comp;
        }
        double *work = RIB_ALLOC_ARRAY(arena, double, (size_t)max_chain);
        for (size_t c = 0; c < S.n_chn; c++) {
            int32_t f = S.chn[c].first, n = S.chn[c].count;
            if (n < 2) continue;
            /* reference_u supplies the additive material gauge, while the
             * measured slice chain supplies the metric.  Projecting each row
             * to gauge+s is the closed-form no-pair version of Stage C: it
             * removes chart-seam jumps without a global sparse solve and makes
             * each fitted claim interval equal the physical chain length. */
            for (int32_t i = 0; i < n; i++)
                work[i] = S.smp[f+i].u - S.smp[f+i].s;
            double gauge = select_median_dbl(work, (size_t)n);
            for (int32_t i = 0; i < n; i++) {
                double next = gauge + S.smp[f+i].s;
                if (fabs(next - S.smp[f+i].u) > 1e-9) moved++;
                S.smp[f+i].u = next;
            }
        }
        stage_t0 = ves_clock_sec();
        PairSet fast_pairs;
        ChainRelationGraph fast_relations;
        memset(&fast_relations, 0, sizeof(fast_relations));
        if (build_pairs(arena, &S, opts, &fast_pairs) != 0 ||
            solve_reference_chain_gauges(
                arena, &S, &fast_pairs, opts, NULL,
                &fast_relations, out) != 0) {
            chain_relation_graph_dispose(&fast_relations);
            Arena_restore(arena, fast_mark);
            memset(out, 0, sizeof(*out));
            return -1;
        }
        if (rib_assign_reconstruction_components(
                arena, &S, opts, &fast_relations, out) != 0) {
            chain_relation_graph_dispose(&fast_relations);
            Arena_restore(arena, fast_mark);
            memset(out, 0, sizeof(*out));
            return -1;
        }
        if (opts->component_global && !opts->emit_global)
            pack_metric_islands(arena, &S, (double)opts->grid_u,
                                NULL, NULL, NULL, NULL, NULL, NULL, 0,
                                NULL, out);
        out->n_pairs = fast_pairs.n_cross;
        out->n_cont_pairs = fast_pairs.n_cont;
        out->match_cover = fast_pairs.candidates_a
                         ? (double)fast_pairs.matched_a /
                           (double)fast_pairs.candidates_a
                         : 0.0;
        fprintf(stderr, "  ribbon fast chain alignment: %.2fs\n",
                ves_clock_sec() - stage_t0);
        out->mono_repairs = moved;
        out->w_groups = (int)n_mesh_comp;
        out->w_prior_groups = (int)n_mesh_comp;
        out->pitch_source = opts->wrap_spacing > 0.0f
                          ? RIB_PITCH_PINNED : RIB_PITCH_FALLBACK;
        out->pitch_used = opts->wrap_spacing;
        for (size_t i = 0; i < S.n_smp; i++) {
            if (S.smp[i].chain < 0) continue;
            if (!isfinite(S.smp[i].u)) {
                Arena_restore(arena, fast_mark);
                memset(out, 0, sizeof(*out));
                return -1;
            }
            if (opts->reference_phi == NULL) S.smp[i].phi = S.smp[i].u;
            if (S.smp[i].u < umin) umin = S.smp[i].u;
            if (S.smp[i].u > umax) umax = S.smp[i].u;
            if (S.smp[i].phi < pmin) pmin = S.smp[i].phi;
            if (S.smp[i].phi > pmax) pmax = S.smp[i].phi;
        }
        if (umin >= 1e299 || !(umax > umin)) {
            Arena_restore(arena, fast_mark);
            memset(out, 0, sizeof(*out));
            return -1;
        }
        out->u_origin = umin;
        out->u_span = umax - umin;
        out->phi_span_turns = pmin < 1e299
                            ? (pmax - pmin) / (2.0 * M_PI) : 0.0;
        {
            double error_sum = 0.0, error_max = 0.0;
            long edges = 0;
            memset(out->duds_hist, 0, sizeof(out->duds_hist));
            for (size_t c = 0; c < S.n_chn; c++) {
                int32_t f = S.chn[c].first, n = S.chn[c].count;
                for (int32_t i = 1; i < n; i++) {
                    double ds = S.smp[f+i].s - S.smp[f+i-1].s;
                    if (ds < 1e-9) continue;
                    double e = fabs((S.smp[f+i].u-S.smp[f+i-1].u)/ds-1.0);
                    double pc = e * 100.0;
                    int bin = pc < 1 ? 0 : pc < 2 ? 1 : pc < 5 ? 2
                            : pc < 10 ? 3 : pc < 20 ? 4 : 5;
                    error_sum += e;
                    if (e > error_max) error_max = e;
                    out->duds_hist[bin]++;
                    edges++;
                }
            }
            out->duds_err_mean = edges ? error_sum / (double)edges : 0.0;
            out->duds_err_max = error_max;
        }
        Arena_restore(arena, fast_mark);
        stage_t0 = ves_clock_sec();
        if (opts->fit_ribbon)
            fit_ribbon(arena, &S, umin, opts, &fast_relations, out);
        else
            out->grid_dv = opts->slice_h;
        chain_relation_graph_dispose(&fast_relations);
        fprintf(stderr,
                "  ribbon fast reference-U fit: span=%.1f, metric-projected=%zu, "
                "grid=%zux%zu (%.2fs)\n",
                out->u_span, out->mono_repairs, out->nu, out->nk,
                ves_clock_sec() - stage_t0);
        return opts->fit_ribbon && out->grid_pos == NULL ? -1 : 0;
    }

    /* r_ref = median in-plane sample radius */
    double r_ref = 1.0;
    stage_t0 = ves_clock_sec();
    {
        Arena_Mark m2 = Arena_save(arena);
        double *rad = RIB_ALLOC_ARRAY(arena, double, S.n_smp + 1);
        size_t q = 0;
        for (size_t i = 0; i < S.n_smp; i++)
            if (S.smp[i].chain >= 0) rad[q++] = S.smp[i].r;
        if (q > 0) r_ref = select_median_dbl(rad, q);
        if (r_ref < 1e-6) r_ref = 1.0;
        Arena_restore(arena, m2);
    }
    fprintf(stderr, "  ribbon stage radius-median: %.2fs\n",
            ves_clock_sec() - stage_t0);

    /* B: pair candidates */
    stage_t0 = ves_clock_sec();
    PairSet P;
    build_pairs(arena, &S, opts, &P);
    out->n_pairs      = P.n_cross;
    out->n_cont_pairs = P.n_cont;
    out->match_cover  = P.candidates_a ? (double)P.matched_a / (double)P.candidates_a : 0.0;
    fprintf(stderr, "  ribbon stage pairs: %.2fs\n",
            ves_clock_sec() - stage_t0);
    if (opts->reference_u != NULL && opts->solve_reference_u)
        fprintf(stderr,
                "  quadribbon U correspondence gate: rejected=%zu "
                "candidate contacts with |delta U| > %.3g\n",
                P.seed_u_rejects,
                fmax((double)opts->match_r,
                     3.0 * (double)opts->sample_h));

    /* W: integer winding index per chain -> phi = th + 2pi W */
    stage_t0 = ves_clock_sec();
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
    fprintf(stderr, "  ribbon stage chain-winding: %.2fs\n",
            ves_clock_sec() - stage_t0);

    /* Coarse StrokeStrip level for a parameterized quadribbon.  Each slice is
     * first made exactly metric, u=s+g, then the existing robust chain graph
     * solves only the additive gauges g across V.  This removes the long-wave
     * slice-origin mode before any fine sample solve; asking millions of sample
     * variables to rediscover a few hundred gauges was both ill-conditioned and
     * the source of the diagonal wedge artifacts. */
    if (opts->reference_u != NULL && opts->solve_reference_u &&
        opts->metric_project_only) {
        stage_t0 = ves_clock_sec();
        int32_t max_chain = 2;
        for (size_t c = 0; c < S.n_chn; c++)
            if (S.chn[c].count > max_chain) max_chain = S.chn[c].count;
        double *gauge_sample = RIB_ALLOC_ARRAY(
            arena, double, (size_t)max_chain);
        /* The classic solver path keeps the plain median.  Metric projection
         * uses rib_project_exact_chain_metric below, after the pair/winding
         * graph exists, so its posterior-gated gauges can be reconciled
         * without ever releasing the exact within-row metric. */
        if (!opts->metric_project_only) {
            /* classic solver pre-step: plain per-chain metric gauge */
            for (size_t c = 0; c < S.n_chn; c++) {
                int32_t f = S.chn[c].first, count = S.chn[c].count;
                if (count < 1) continue;
                for (int32_t j = 0; j < count; j++)
                    gauge_sample[j] = S.smp[f + j].u - S.smp[f + j].s;
                double gauge = select_median_dbl(gauge_sample, (size_t)count);
                for (int32_t j = 0; j < count; j++)
                    S.smp[f + j].u = gauge + S.smp[f + j].s;
            }
        }
        if (!opts->metric_project_only) {
            RibbonOpts gauge_opts = *opts;
            gauge_opts.reference_anchor_gauges = 0;
            if (solve_reference_chain_gauges(
                    arena, &S, &P, &gauge_opts, NULL, NULL, out) != 0) {
                fprintf(stderr, "ribbon: coarse chain-gauge solve failed\n");
                rib_gmg_snapshots_dispose(level_capture_ptr);
                return -1;
            }
            fprintf(stderr, "  ribbon stage coarse chain gauges: %.2fs\n",
                    ves_clock_sec() - stage_t0);
        } else {
            /* Metric projection v3: F(u_ref) supplies the common cover frame;
             * rigid phase placement (when authoritative) removes upstream
             * atlas packing; exact per-chain projection supplies geometry;
             * and a chain-only anchored solve reconciles crack/vertical ties.
             * There is no sample-space solve and no post-solve clipping. */
            /* Islands: with an AUTHORITATIVE carried phase (the fit's
             * registration sidecar), carried U is atlas packing and carried
             * phase is the physical placement evidence -- rigid median
             * placement only.  A recomputed fallback phase is unanchored on
             * detached components, so the carried registration stands. */
            if (opts->reference_phi_authoritative)
                rib_project_register_islands(arena, &S, r_ref,
                                             (double)opts->wrap_spacing, 0,
                                             out);
            else
                fprintf(stderr,
                        "  island placement: carried registration kept "
                        "(no authoritative phase sidecar)\n");
            if (rib_project_exact_chain_metric(
                    arena, &S, row_sample_uconf, out) != 0) {
                fprintf(stderr, "ribbon: exact chain metric failed\n");
                rib_gmg_snapshots_dispose(level_capture_ptr);
                return -1;
            }
            {
                PairSet gauge_pairs;
                RibbonOpts gauge_opts = *opts;
                gauge_opts.reference_anchor_gauges = 1;
                if (rib_metric_add_topology_gauge_pairs(
                        arena, faces, nf, row_vertex_sample, nv, &S,
                        row_sample_uconf, &gauge_opts, &P,
                        &gauge_pairs) != 0) {
                    fprintf(stderr,
                            "ribbon: cannot build exact-chain topology "
                            "relations\n");
                    rib_gmg_snapshots_dispose(level_capture_ptr);
                    return -1;
                }
                if (solve_reference_chain_gauges(
                        arena, &S, &gauge_pairs, &gauge_opts,
                        row_sample_uconf, NULL, out) != 0) {
                    fprintf(stderr,
                            "ribbon: anchored chain-gauge solve failed\n");
                    rib_gmg_snapshots_dispose(level_capture_ptr);
                    return -1;
                }
            }
            {
                double u_min = INFINITY;
                for (size_t i = 0; i < S.n_smp; i++)
                    if (isfinite(S.smp[i].u) && S.smp[i].u < u_min)
                        u_min = S.smp[i].u;
                if (isfinite(u_min) && fabs(u_min) > 1e-9) {
                    for (size_t i = 0; i < S.n_smp; i++)
                        S.smp[i].u -= u_min;
                    fprintf(stderr,
                            "  metric projection re-zero: u shifted by "
                            "%+.1f vox\n", -u_min);
                }
            }
            fprintf(stderr,
                    "  ribbon stage metric projection v3: F init + phase "
                    "placement + exact chains + anchored chain gauges "
                    "(%.2fs)\n",
                    ves_clock_sec() - stage_t0);
        }
    }

    /* C: solve.  In metric-projection mode there is nothing to solve: the
     * carried frame is the placement authority and each chain is already
     * exactly metric.  Build the same StrokeStrip rows the solver would use,
     * but only to AUDIT the projected coordinate; no gauge, stitch,
     * registration, or orientation pass may modify it. */
    stage_t0 = ves_clock_sec();
    if (opts->reference_u != NULL && opts->solve_reference_u &&
        opts->metric_project_only) {
        Arena_Mark audit_mark = Arena_save(arena);
        RibStripSet strip;
        for (size_t c = 0; c < S.n_chn; c++) {
            S.chn[c].solve_comp = (int32_t)c;
            S.chn[c].group = (int32_t)c;
        }
        out->n_qp_comps = (int)S.n_chn;
        if (rib_strip_build_runs(arena, &S, &P, 1,
                                 RIB_STRIP_INITIAL_LIKE, &strip) != 0) {
            fprintf(stderr,
                    "ribbon: cannot construct StrokeStrip audit runs\n");
            rib_gmg_snapshots_dispose(level_capture_ptr);
            return -1;
        }
        rib_strip_update_length_coefficients(&S, &strip);
        rib_strip_audit_reference_seed(arena, &S, &P, &strip);
        out->n_strip_runs = strip.nrun;
        out->n_strip_members = strip.nmember;
        out->n_strip_links = strip.nlink;
        out->n_strip_links_pruned = strip.npruned;
        {
            double *uaudit = RIB_ALLOC_ARRAY(arena, double, S.n_smp + 1);
            for (size_t i = 0; i < S.n_smp; i++) uaudit[i] = S.smp[i].u;
            rib_strip_measure(&S, &strip, uaudit, 0,
                              &out->strip_length_rms, &out->strip_align_rms);
        }
        {
            double esum = 0.0, emax = 0.0;
            long nedge = 0;
            memset(out->duds_hist, 0, sizeof(out->duds_hist));
            for (size_t c = 0; c < S.n_chn; c++) {
                int32_t f = S.chn[c].first, cn = S.chn[c].count;
                for (int32_t i = 1; i < cn; i++) {
                    double l = S.smp[f+i].s - S.smp[f+i-1].s;
                    if (l < 1e-9) continue;
                    double e = fabs((S.smp[f+i].u - S.smp[f+i-1].u) / l - 1.0);
                    int bin = e < 0.01 ? 0 : e < 0.02 ? 1 : e < 0.05 ? 2
                            : e < 0.10 ? 3 : e < 0.20 ? 4 : 5;
                    esum += e;
                    if (e > emax) emax = e;
                    out->duds_hist[bin]++;
                    nedge++;
                }
            }
            out->duds_err_mean = nedge ? esum / (double)nedge : 0.0;
            out->duds_err_max = emax;
        }
        fprintf(stderr,
                "  ribbon stage metric-projection audit: runs=%zu links=%zu "
                "align_rms=%.3f length_rms=%.5f (%.2fs)\n",
                out->n_strip_runs, out->n_strip_links, out->strip_align_rms,
                out->strip_length_rms, ves_clock_sec() - stage_t0);
        Arena_restore(arena, audit_mark);
        stage_t0 = ves_clock_sec();
    } else {
        solve_parameterization(arena, &S, &P, r_ref, opts, out,
                               level_capture_ptr);
    }
    if (opts->component_global) {
        ChainRelationGraph branch_relations;
        memset(&branch_relations, 0, sizeof(branch_relations));
        if (chain_relation_graph_from_pairs(&S, &P, &branch_relations) != 0 ||
            rib_assign_reconstruction_components(
                arena, &S, opts, &branch_relations, out) != 0) {
            chain_relation_graph_dispose(&branch_relations);
            rib_gmg_snapshots_dispose(level_capture_ptr);
            memset(out, 0, sizeof(*out));
            return -1;
        }
        chain_relation_graph_dispose(&branch_relations);
    }
    /* A carried reference U is an observed cover coordinate, not an arbitrary
     * solve-component gauge.  Packing reconstruction components here used to
     * translate pieces of that certificate independently after the consistency
     * pass had preserved it, creating new cross-ply order inversions.
     * Unobserved solver gauges still require atlas packing; a carried
     * certificate must retain its single global gauge.
     *
     * metric_project_only and preserve_input_rows are two named instances of
     * that contract, but a quadribbon --scaffold-solve fit of a winding-only
     * VMESH carries the same certificate: winding registration already placed
     * every mesh component on one global U.  Re-packing its several hundred
     * branch-derived reconstruction lanes shredded a 21.8k-vox cover into a
     * 71.2k-vox atlas of confetti.  Gate on the certificate itself. */
    if (opts->component_global && !opts->emit_global) {
        int carried_certificate = opts->reference_u != NULL &&
                                  opts->solve_reference_u;
        if (carried_certificate) {
            fprintf(stderr,
                    "  ribbon certificate contract: retained solved U gauge "
                    "(metric-island packing disabled)\n");
        } else {
            pack_metric_islands(arena, &S, (double)opts->grid_u,
                                NULL, NULL, NULL, NULL, NULL, NULL, 0,
                                level_capture_ptr, out);
        }
    }
    {
        double lo = 1e300, hi = -1e300;
        for (size_t i = 0; i < S.n_smp; i++) {
            if (S.smp[i].chain < 0) continue;
            if (S.smp[i].u < lo) lo = S.smp[i].u;
            if (S.smp[i].u > hi) hi = S.smp[i].u;
        }
        out->u_span = hi - lo;
    }
    fprintf(stderr, "  ribbon stage parameterization: %.2fs\n",
            ves_clock_sec() - stage_t0);

    if (opts->direct_ribbon) {
        /* The quad-ribbon experiment ends in slice space.  Do not map the
         * solved coordinate back onto source-mesh vertices and do not let that
         * legacy transfer reconcile gauges or pack islands afterward. */
        double U0 = 1e300;
        stage_t0 = ves_clock_sec();
        for (size_t i = 0; i < S.n_smp; i++)
            if (S.smp[i].chain >= 0 && S.smp[i].u < U0) U0 = S.smp[i].u;
        if (U0 >= 1e299) {
            rib_gmg_snapshots_dispose(level_capture_ptr);
            return -1;
        }
        out->u_origin = U0;
        if (opts->fit_ribbon) {
            if (rib_gmg_dump_levels(
                    arena, &S, opts, out, level_capture_ptr) != 0) {
                fprintf(stderr, "ERROR: GMG level checkpoint emission failed\n");
                rib_gmg_snapshots_dispose(level_capture_ptr);
                return -1;
            }
            fit_ribbon(arena, &S, U0, opts, NULL, out);
        } else {
            out->grid_dv = opts->slice_h;
        }
        rib_gmg_snapshots_dispose(level_capture_ptr);
        fprintf(stderr, "  ribbon stage direct slice-grid fit: %.2fs\n",
                ves_clock_sec() - stage_t0);
        return opts->fit_ribbon && out->grid_pos == NULL ? -1 : 0;
    }

    /* D: transfer to mesh verts. Pre-fill out->uv with vertex in-plane coords
     * (the transfer's bucketing key); replaced by (u,v) inside. */
    stage_t0 = ves_clock_sec();
    out->uv = RIB_ALLOC_ARRAY(arena, float, nv * 2);
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
        double U0 = 0.0;
        if (row_vertex_sample != NULL) {
            /* Row mode: every finite-U vertex IS its own sample, so the solved
             * coordinate writes back exactly; the KD transfer's uniform
             * slice_h buckets cannot address rank-indexed irregular rows. */
            double umin_ = 1e300;
            for (size_t i = 0; i < S.n_smp; i++)
                if (S.smp[i].chain >= 0 && S.smp[i].u < umin_)
                    umin_ = S.smp[i].u;
            if (umin_ >= 1e299) umin_ = 0.0;
            uint8_t *mapped = (uint8_t *)ARENA_CALLOC(arena, nv, 1);
            size_t direct = 0, filled = 0, fallback = 0;
            for (size_t i = 0; i < nv; i++) {
                int32_t si = row_vertex_sample[i];
                out->uv[i*2+1] = (float)(t[i] - tmin);
                if (si < 0) continue;
                out->uv[i*2+0] = (float)(S.smp[si].u - umin_);
                mapped[i] = 1;
                direct++;
            }
            for (int pass = 0; pass < 8 && direct + filled < nv; pass++) {
                size_t moved = 0;
                for (size_t ff = 0; ff < nf; ff++) {
                    for (int k = 0; k < 3; k++) {
                        int32_t a = faces[ff*3+k];
                        int32_t b = faces[ff*3+(k+1)%3];
                        if (!mapped[a] && mapped[b]) {
                            out->uv[a*2+0] = out->uv[b*2+0];
                            mapped[a] = 1;
                            moved++;
                        }
                    }
                }
                filled += moved;
                if (moved == 0) break;
            }
            for (size_t i = 0; i < nv; i++) {
                if (mapped[i]) continue;
                out->uv[i*2+0] = 0.0f;   /* unreferenced/skipped: finite, counted */
                fallback++;
            }
            out->uv_filled = filled;
            out->uv_fallback = fallback;
            U0 = umin_;
            fprintf(stderr,
                    "  row-mode writeback: %zu exact, %zu neighbour-filled, "
                    "%zu fallback\n", direct, filled, fallback);
        } else {
            U0 = transfer_uv(arena, verts, nv, faces, nf,
                             vertex_mesh_comp, n_mesh_comp, t, &S,
                             (double)opts->slice_h, opts,
                             level_capture_ptr, out);
        }
        /* E: ribbon */
        if (opts->fit_ribbon) {
            if (rib_gmg_dump_levels(
                    arena, &S, opts, out, level_capture_ptr) != 0) {
                fprintf(stderr, "ERROR: GMG level checkpoint emission failed\n");
                rib_gmg_snapshots_dispose(level_capture_ptr);
                return -1;
            }
            fit_ribbon(arena, &S, U0, opts, NULL, out);
            if (opts->verify_fit_width == 1 && out->grid_pos != NULL) {
                RibbonOpts alt_opts = *opts;
                RibbonResult alt;
                memset(&alt, 0, sizeof(alt));
                alt_opts.fit_cover_width = !opts->fit_cover_width;
                alt_opts.verify_fit_width = 2;
                fit_ribbon(arena, &S, U0, &alt_opts, NULL, &alt);
                if (alt.grid_pos != NULL) {
                    size_t shared_nu = out->nu < alt.nu ? out->nu : alt.nu;
                    size_t shared_nk = out->nk < alt.nk ? out->nk : alt.nk;
                    size_t valid_mismatch = 0, value_mismatch = 0, both = 0;
                    double max_delta = 0.0;
                    for (size_t k = 0; k < shared_nk; k++) {
                        for (size_t j = 0; j < shared_nu; j++) {
                            const float *a = &out->grid_pos[(k*out->nu+j)*3];
                            const float *b = &alt.grid_pos[(k*alt.nu+j)*3];
                            int av = out->grid_present != NULL
                                   ? out->grid_present[k*out->nu+j] != 0
                                   : isfinite((double)a[0]);
                            int bv = alt.grid_present != NULL
                                   ? alt.grid_present[k*alt.nu+j] != 0
                                   : isfinite((double)b[0]);
                            if (av != bv) { valid_mismatch++; continue; }
                            if (!av) continue;
                            both++;
                            double delta = 0.0;
                            for (int d = 0; d < 3; d++) {
                                double dd = fabs((double)a[d] - (double)b[d]);
                                if (dd > delta) delta = dd;
                            }
                            if (delta > 1e-6) value_mismatch++;
                            if (delta > max_delta) max_delta = delta;
                        }
                    }
                    fprintf(stderr,
                            "  fit-width compare: compact_nu=%zu cover_nu=%zu "
                            "shared_valid=%zu validity_mismatch=%zu "
                            "value_mismatch=%zu max_delta=%.6g\n",
                            opts->fit_cover_width ? alt.nu : out->nu,
                            opts->fit_cover_width ? out->nu : alt.nu,
                            both, valid_mismatch, value_mismatch, max_delta);
                } else {
                    fprintf(stderr,
                            "  fit-width compare: alternate fit unavailable "
                            "(guard or allocation limit)\n");
                }
            }
        } else
            out->grid_dv = opts->slice_h;
    }
    rib_gmg_snapshots_dispose(level_capture_ptr);
    fprintf(stderr, "  ribbon stage transfer%s: %.2fs\n",
            opts->fit_ribbon ? "+fit" : "", ves_clock_sec() - stage_t0);
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

    float   *uv = RIB_ALLOC_ARRAY(arena, float, nv * 2);
    uint8_t *ok = (uint8_t *)ARENA_CALLOC(arena, nv, 1);

    Arena_Mark mark = Arena_save(arena);

    /* coarse vertex -> incident faces (CSR) */
    int32_t *voff = (int32_t *)ARENA_CALLOC(arena, cnv + 1,
                                            sizeof(int32_t));
    for (size_t f = 0; f < cnf * 3; f++) voff[cfaces[f] + 1]++;
    for (size_t i = 0; i < cnv; i++) voff[i+1] = (int32_t)(voff[i+1] + voff[i]);
    int32_t *vfac = RIB_ALLOC_ARRAY(arena, int32_t, cnf * 3);
    int32_t *vcur = RIB_ALLOC_ARRAY(arena, int32_t, cnv);
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
    float   *v = RIB_ALLOC_ARRAY(arena, float, nvv * 3);
    int32_t *f = RIB_ALLOC_ARRAY(
        arena, int32_t, (size_t)(nphi - 1) * (size_t)(nh - 1) * 6);
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

    /* Pair-only coarsening leaves a contracted star almost unchanged.  The
     * unmatched-attachment pass must consume a bounded one-ring while retaining
     * ordinary 2x aggregation on a path. */
    {
        int32_t star_a[20], star_b[20];
        int32_t severe_a[200], severe_b[200];
        int32_t path_a[7] = {0, 1, 2, 3, 4, 5, 6};
        int32_t path_b[7] = {1, 2, 3, 4, 5, 6, 7};
        int32_t agg[201];
        for (int i = 0; i < 20; i++) {
            star_a[i] = 0;
            star_b[i] = i + 1;
        }
        int star_n = rib_amg_axis_aggregate(star_a, star_b, 20, 21, agg);
        for (int i = 0; i < 200; i++) {
            severe_a[i] = 0;
            severe_b[i] = i + 1;
        }
        int severe_n = rib_amg_axis_aggregate(
            severe_a, severe_b, 200, 201, agg);
        int path_n = rib_amg_axis_aggregate(path_a, path_b, 7, 8, agg);
        int ok = star_n == 14 && severe_n == 1 && path_n == 4;
        st_check(ok, "GMG geometric star/path aggregation", &fails);
        fprintf(stderr,
                "[ribbon selftest] (gmg-aggregate) star=%d severe=%d "
                "path=%d: %s\n",
                star_n, severe_n, path_n, ok ? "ok" : "FAIL");
    }

    /* The ribbon specialization of StrokeStrip has one row per maximal
     * connected run in V.  Removing exactly one adjacent-V link must split
     * only that material ruling: all samples remain members of a row, while
     * unrelated rulings remain connected.  A competing branch is pruned
     * deterministically instead of merging two observations from one slice. */
    {
        enum { NV = 4, NU = 3, NS = NV * NU, MAXPAIR = 10 };
        Sample smp[NS];
        Chain chn[NV];
        SliceSet S;
        RibbonOpts o;
        int all_ok = 1;

        memset(smp, 0, sizeof smp);
        memset(chn, 0, sizeof chn);
        memset(&S, 0, sizeof S);
        memset(&o, 0, sizeof o);
        for (int v = 0; v < NV; v++) {
            chn[v].first = v * NU;
            chn[v].count = NU;
            chn[v].slice = v;
            for (int q = 0; q < NU; q++) {
                int i = v * NU + q;
                smp[i].chain = v;
                smp[i].slice = v;
                smp[i].s = (double)q;
                smp[i].c1 = (double)q;
                smp[i].c2 = 0.0;
                smp[i].p[0] = (double)v;
                smp[i].p[2] = (double)q;
                smp[i].tau[2] = 1.0;
            }
        }
        S.smp = smp;
        S.n_smp = NS;
        S.chn = chn;
        S.n_chn = NV;
        S.nplanes = NV;

        for (int fixture = 0; fixture < 3; fixture++) {
            Pair pair[MAXPAIR];
            PairSet P;
            RibStripSet R;
            int pn = 0;
            memset(pair, 0, sizeof pair);
            memset(&P, 0, sizeof P);
            for (int v = 0; v + 1 < NV; v++) {
                for (int q = 0; q < NU; q++) {
                    if (fixture == 1 && v == 1 && q == 1) continue;
                    pair[pn].a = v * NU + q;
                    pair[pn].b = (v + 1) * NU + q;
                    pair[pn].like = 1.0;
                    pn++;
                }
            }
            if (fixture == 2) {
                pair[pn].a = 0 * NU + 1;
                pair[pn].b = 1 * NU + 0;
                pair[pn].like = 1.0;
                pn++;
            }
            P.pairs = pair;
            P.n_pairs = P.n_cross = (size_t)pn;

            Arena_Mark mark = Arena_save(arena);
            int rc = rib_strip_build_runs(arena, &S, &P, 0,
                                          RIB_STRIP_INITIAL_LIKE, &R);
            int ok = rc == 0;
            size_t expect_runs = fixture == 1 ? 4u : 3u;
            size_t expect_links = fixture == 1 ? 8u : 9u;
            size_t expect_pruned = fixture == 2 ? 1u : 0u;
            if (ok) {
                rib_strip_update_length_coefficients(&S, &R);
                ok = R.nrun == expect_runs && R.nmember == NS &&
                     R.nlink == expect_links &&
                     R.npruned == expect_pruned;
            }
            for (size_t r = 0; ok && r < R.nrun; r++) {
                const RibStripRun *run = &R.run[r];
                unsigned slice_mask = 0;
                double base_sum = 0.0, derivative = 0.0;
                for (int32_t j = 0; j < run->count; j++) {
                    const RibStripMember *m =
                        &R.member[run->first + (size_t)j];
                    unsigned bit = 1u << (unsigned)smp[m->sample].slice;
                    if ((slice_mask & bit) != 0) ok = 0;
                    slice_mask |= bit;
                    base_sum += m->base_weight;
                }
                for (size_t k = 2 * run->first;
                     k < 2 * (run->first + (size_t)run->count); k++)
                    derivative += R.coeff[k].value *
                                  smp[R.coeff[k].var].s;
                if (fabs(base_sum - 1.0) > RIB_STRIP_GEOM_EPS ||
                    fabs(derivative - 1.0) > RIB_STRIP_GEOM_EPS)
                    ok = 0;
                if (fixture == 1) {
                    int32_t first_slice = smp[R.member[run->first].sample].slice;
                    int32_t last_slice = smp[
                        R.member[run->first + (size_t)run->count - 1].sample
                    ].slice;
                    if (first_slice < 2 && last_slice >= 2 &&
                        fabs(smp[R.member[run->first].sample].c1 - 1.0) <=
                            RIB_STRIP_GEOM_EPS)
                        ok = 0;
                }
            }
            for (int i = 0; ok && i < NS; i++)
                if (R.sample_member[i] < 0) ok = 0;
            fprintf(stderr,
                    "[ribbon selftest] (strokestrip-v-runs/%s) "
                    "runs=%zu members=%zu links=%zu pruned=%zu: %s\n",
                    fixture == 0 ? "connected" :
                    fixture == 1 ? "one-gap" : "branch",
                    rc == 0 ? R.nrun : 0u, rc == 0 ? R.nmember : 0u,
                    rc == 0 ? R.nlink : 0u, rc == 0 ? R.npruned : 0u,
                    ok ? "ok" : "FAIL");
            if (!ok) all_ok = 0;
            Arena_restore(arena, mark);
        }
        st_check(all_ok,
                 "StrokeStrip V runs split only at the missing connection",
                 &fails);
    }

    /* The whole-scroll driver unwraps cubes in an outer OpenMP region.  Keep a
     * radix input above the inner parallelism threshold here so a serialized
     * nested team cannot silently process only its first requested band. */
    {
        const size_t n = 300000;
        uint64_t *a = (uint64_t *)malloc(n * sizeof(*a));
        int rc = -1, ok = a != NULL;
        if (a != NULL) {
            for (size_t i = 0; i < n; i++) a[i] = (uint64_t)(n - i);
#ifdef _OPENMP
#pragma omp parallel num_threads(2)
            {
#pragma omp single
                rc = rib_u64_radix(a, n);
            }
#else
            rc = rib_u64_radix(a, n);
#endif
            if (rc != 0) ok = 0;
            for (size_t i = 0; ok && i < n; i++)
                if (a[i] != (uint64_t)(i + 1)) ok = 0;
        }
        fprintf(stderr, "[ribbon selftest] (nested-radix) %s\n",
                ok ? "ok" : "FAIL");
        st_check(ok, "nested OpenMP radix", &fails);
        free(a);
    }

    /* Registered chain gauges are a fully anchored graph problem.  Three
     * slightly inconsistent geometric ladders should improve their pair
     * residual without leaving the graph frame; one gross false edge must be
     * rejected before the first solve.  Reversing pair enumeration must not
     * change the answer. */
    {
        enum { NC = 4, OBS = 20, NS = NC*OBS,
               NGOOD = (NC-1)*OBS, NP = NGOOD+1, NRUN = 3 };
        static const double graph_gauge[NC] = {0.0, 10.8, 20.0, 30.8};
        Sample smp[NRUN][NS];
        Chain chn[NRUN][NC];
        Pair pair[NRUN][NP], canonical[NP];
        SliceSet S[NRUN];
        PairSet P[NRUN];
        RibbonResult R[NRUN];
        RibbonOpts o;
        float phi_present = 0.0f;
        int pn = 0, rc[NRUN] = {-1,-1,-1};

        memset(canonical, 0, sizeof(canonical));
        for (int c = 0; c + 1 < NC; c++) {
            for (int q = 0; q < OBS; q++) {
                Pair *pr = &canonical[pn++];
                pr->a = c*OBS + q;
                pr->b = (c+1)*OBS + q;
                pr->d = -10.0;
                pr->w = RIB_ALIGN_W;
                pr->like = 1.0;
            }
        }
        canonical[pn].a = 0;
        canonical[pn].b = 3*OBS;
        canonical[pn].d = 100.0; /* false: the graph predicts about -30.8 */
        canonical[pn].w = RIB_ALIGN_W;
        canonical[pn].like = 1.0;
        pn++;

        RibbonOpts_default(&o);
        o.reference_phi = &phi_present;
        o.reference_anchor_gauges = 1;
        o.relax_iters = 2;
        o.final_iters = 1;
        o.solve_threads = 1;
        for (int run = 0; run < NRUN; run++) {
            memset(smp[run], 0, sizeof(smp[run]));
            memset(chn[run], 0, sizeof(chn[run]));
            memset(&S[run], 0, sizeof(S[run]));
            memset(&P[run], 0, sizeof(P[run]));
            memset(&R[run], 0, sizeof(R[run]));
            for (int c = 0; c < NC; c++) {
                chn[run][c].first = c*OBS;
                chn[run][c].count = OBS;
                for (int q = 0; q < OBS; q++) {
                    int i = c*OBS + q;
                    smp[run][i].chain = c;
                    smp[run][i].s = (double)q;
                    smp[run][i].u = graph_gauge[c] + (double)q;
                }
            }
            S[run].smp = smp[run];
            S[run].n_smp = NS;
            S[run].chn = chn[run];
            S[run].n_chn = NC;
            memcpy(pair[run], canonical, sizeof(canonical));
            if (run == 1) {
                for (int a = 0, b = NP-1; a < b; a++, b--) {
                    Pair tmp = pair[run][a];
                    pair[run][a] = pair[run][b];
                    pair[run][b] = tmp;
                }
            }
            P[run].pairs = pair[run];
            P[run].n_pairs = run == 2 ? NGOOD : NP;
            rc[run] = solve_reference_chain_gauges(
                arena, &S[run], &P[run], &o, NULL, NULL, &R[run]);
        }

        double order_delta = 0.0, outlier_delta = 0.0;
        int metric_exact = 1;
        for (int c = 0; c < NC; c++) {
            double ga = smp[0][c*OBS].u;
            double gb = smp[1][c*OBS].u;
            double gc = smp[2][c*OBS].u;
            if (fabs(ga-gb) > order_delta) order_delta = fabs(ga-gb);
            if (fabs(ga-gc) > outlier_delta) outlier_delta = fabs(ga-gc);
            for (int q = 1; q < OBS; q++)
                if (fabs((smp[0][c*OBS+q].u -
                          smp[0][c*OBS+q-1].u) - 1.0) > 1e-12)
                    metric_exact = 0;
        }
        int ok = rc[0] == 0 && rc[1] == 0 && rc[2] == 0 &&
                 R[0].chain_gauge_pair_median_after <
                     R[0].chain_gauge_pair_median_before &&
                 R[0].chain_gauge_pair_p95_before < 1.0 &&
                 R[0].chain_gauge_shift_max < 1.0 &&
                 order_delta < 1e-8 && outlier_delta < 1e-3 &&
                 metric_exact;
        fprintf(stderr,
                "[ribbon selftest] (anchored-chain-gauge) "
                "pair-med/p95=%.4f/%.4f->%.4f/%.4f "
                "shift-max=%.4f order-delta=%.3g outlier-delta=%.3g "
                "metric=%s: %s\n",
                R[0].chain_gauge_pair_median_before,
                R[0].chain_gauge_pair_p95_before,
                R[0].chain_gauge_pair_median_after,
                R[0].chain_gauge_pair_p95_after,
                R[0].chain_gauge_shift_max, order_delta, outlier_delta,
                metric_exact ? "exact" : "BROKEN", ok ? "ok" : "FAIL");
        st_check(ok,
                 "anchored chain gauges are bounded, robust and order-invariant",
                 &fails);
    }

    /* Local fitted-grid puncture repair: a single supported centre sample on
     * the wrong ply is replaced from its compact four-neighbour ring.  Moving
     * one ring vertex beyond the inter-wrap certificate must disable repair. */
    {
        float grid[3*3*3];
        uint8_t supported[3*3];
        size_t nsupported = 0;
        for (int k = 0; k < 3; k++) {
            for (int j = 0; j < 3; j++) {
                size_t s = (size_t)k*3+(size_t)j;
                grid[s*3+0] = 0.0f;
                grid[s*3+1] = (float)(2*k);
                grid[s*3+2] = (float)(2*j);
                supported[s] = 1;
            }
        }
        grid[(1*3+1)*3+0] = 20.0f;
        size_t repaired = rib_grid_repair_isolated_outliers(
            grid,supported,3,3,&nsupported);
        int fixed = repaired == 1 && nsupported == 1 &&
                    fabs((double)grid[(1*3+1)*3+0]) < 1e-6 &&
                    fabs((double)grid[(1*3+1)*3+1]-2.0) < 1e-6 &&
                    fabs((double)grid[(1*3+1)*3+2]-2.0) < 1e-6 &&
                    supported[1*3+1] == 0;
        fprintf(stderr,
                "[ribbon selftest] (local-outlier) repaired=%zu "
                "supported=%zu: %s\n",
                repaired,nsupported,fixed ? "ok" : "FAIL");
        st_check(fixed,"isolated fitted-grid outlier repaired",&fails);

        grid[(1*3+1)*3+0] = 20.0f;
        supported[1*3+1] = 1;
        grid[(0*3+1)*3+0] = 8.0f;
        repaired = rib_grid_repair_isolated_outliers(
            grid,supported,3,3,&nsupported);
        st_check(repaired == 0,
                 "wide neighbour ring is not bridged",&fails);
    }

    {
        float grid[3*4*3];
        uint8_t supported[3*4];
        size_t slots=0,nsupported=0;
        for (int k=0;k<3;k++) {
            for (int j=0;j<4;j++) {
                size_t s=(size_t)k*4+(size_t)j;
                grid[s*3]=0.0f;
                grid[s*3+1]=(float)(2*k);
                grid[s*3+2]=(float)(2*j);
                supported[s]=1;
            }
        }
        grid[(1*4+1)*3]=20.0f;
        grid[(1*4+2)*3]=20.0f;
        size_t patches=rib_grid_repair_outlier_pairs(
            grid,supported,3,4,&slots,&nsupported);
        int fixed=patches==1 && slots==2 && nsupported==2 &&
                  fabs((double)grid[(1*4+1)*3])<1e-6 &&
                  fabs((double)grid[(1*4+2)*3])<1e-6 &&
                  supported[1*4+1]==0 && supported[1*4+2]==0;
        fprintf(stderr,
                "[ribbon selftest] (local-pair) patches=%zu slots=%zu "
                "supported=%zu: %s\n",
                patches,slots,nsupported,fixed?"ok":"FAIL");
        st_check(fixed,"adjacent fitted-grid outlier pair repaired",&fails);
    }

    {
        enum { NK=3, NU=7, N=4 };
        float grid[NK*NU*3];
        uint8_t supported[NK*NU];
        size_t slots=0,nsupported=0;
        for (int k=0;k<NK;k++) {
            for (int j=0;j<NU;j++) {
                size_t s=(size_t)k*NU+(size_t)j;
                grid[s*3]=0.0f;
                grid[s*3+1]=(float)(2*k);
                grid[s*3+2]=(float)(2*j);
                supported[s]=1;
            }
        }
        for (int j=1;j<=N;j++) grid[(1*NU+j)*3]=20.0f;
        size_t patches=rib_grid_repair_outlier_runs(
            grid,supported,NK,NU,&slots,&nsupported);
        int fixed=patches==1 && slots==N && nsupported==N;
        for (int j=1;j<=N;j++)
            fixed=fixed && fabs((double)grid[(1*NU+j)*3])<1e-6 &&
                  supported[1*NU+j]==0;
        fprintf(stderr,
                "[ribbon selftest] (local-run) patches=%zu slots=%zu "
                "supported=%zu: %s\n",
                patches,slots,nsupported,fixed?"ok":"FAIL");
        st_check(fixed,"fitted-grid outlier strip repaired",&fails);
    }

    {
        enum { NK=11, NU=11 };
        float grid[NK*NU*3];
        uint8_t supported[NK*NU];
        size_t slots=0,nsupported=0,nsmall=0,nenclosed=0,nclaims=0;
        for (int k=0;k<NK;k++) {
            for (int j=0;j<NU;j++) {
                size_t s=(size_t)k*NU+(size_t)j;
                grid[s*3]=0.0f;
                grid[s*3+1]=(float)(2*k);
                grid[s*3+2]=(float)(2*j);
                supported[s]=1;
            }
        }
        for (int k=4;k<=5;k++)
            for (int j=4;j<=5;j++)
                grid[((size_t)k*NU+(size_t)j)*3]=20.0f;
        size_t patches=rib_grid_repair_enclosed_outlier_islands(
            grid,supported,NULL,0.0,1.0,NULL,NULL,NK,NU,
            &slots,&nsupported,&nsmall,&nenclosed,&nclaims);
        int fixed=patches==1 && slots==4 && nsupported==4 &&
                  nsmall>=1 && nenclosed==1;
        for (int k=4;k<=5;k++) {
            for (int j=4;j<=5;j++) {
                size_t s=(size_t)k*NU+(size_t)j;
                fixed=fixed && fabs((double)grid[s*3])<1e-6 &&
                      supported[s]==0;
            }
        }
        fprintf(stderr,
                "[ribbon selftest] (local-patch) patches=%zu slots=%zu "
                "supported=%zu small=%zu enclosed=%zu: %s\n",
                patches,slots,nsupported,nsmall,nenclosed,
                fixed?"ok":"FAIL");
        st_check(fixed,"enclosed fitted-grid island repaired",&fails);

    }

    {
        enum { NK=3, NU=3 };
        float grid[NK*NU*3],candidate[3];
        float target[3]={5.0f,4.0f,4.0f};
        int32_t label[NK*NU];
        size_t slot[1]={1*NU+1},touched=0;
        int before=0,after=0;
        for (int k=0;k<NK;k++) {
            for (int j=0;j<NU;j++) {
                size_t s=(size_t)k*NU+(size_t)j;
                grid[s*3]=0.0f;
                grid[s*3+1]=(float)(4*k);
                grid[s*3+2]=(float)(4*j);
                label[s]=0;
            }
        }
        grid[slot[0]*3]=20.0f;
        label[slot[0]]=1;
        int refined=rib_grid_patch_refine_full_cells(
            grid,label,1,NU,slot,1,target,candidate,
            1,1,1,1,&before,&after,&touched);
        double dz=(double)candidate[0]-(double)target[0];
        double dy=(double)candidate[1]-(double)target[1];
        double dx=(double)candidate[2]-(double)target[2];
        int bounded=refined && touched==4 && after==8 &&
                    sqrt(dz*dz+dy*dy+dx*dx)<=
                    RIB_GRID_LOCAL_PATCH_REFINE_MAX+1e-4;
        fprintf(stderr,
                "[ribbon selftest] (local-patch-refine) legal=%d->%d/%zu "
                "move=%.4f: %s\n",
                before,after,2*touched,sqrt(dz*dz+dy*dy+dx*dx),
                bounded?"ok":"FAIL");
        st_check(bounded,"bounded full-cell patch refinement",&fails);
    }

    {
        enum { NK=5, NU=5, NV=9 };
        float grid[NK*NU*3],target[NV*3],candidate[NV*3];
        size_t core[1]={2*NU+2};
        size_t slot[NV]={
            2*NU+2,
            1*NU+1,1*NU+2,1*NU+3,
            2*NU+1,       2*NU+3,
            3*NU+1,3*NU+2,3*NU+3
        };
        double limit[NV];
        int before=0,after=0;
        for (int k=0;k<NK;k++) {
            for (int j=0;j<NU;j++) {
                size_t s=(size_t)k*NU+(size_t)j;
                grid[s*3]=0.0f;
                grid[s*3+1]=(float)(4*k);
                grid[s*3+2]=(float)(4*j);
            }
        }
        grid[core[0]*3]=20.0f;
        for (int q=0;q<NV;q++) {
            memcpy(&target[q*3],&grid[slot[q]*3],3*sizeof(float));
            limit[q]=q==0?1.0:0.75;
        }
        target[0]=6.0f;
        int refined=rib_grid_patch_refine_neighbourhood(
            grid,NK,NU,core,1,slot,NV,target,limit,candidate,
            &before,&after);
        int bounded=refined && after>before;
        for (int q=0;q<NV;q++) {
            double d0=(double)candidate[q*3]-(double)target[q*3];
            double d1=(double)candidate[q*3+1]-(double)target[q*3+1];
            double d2=(double)candidate[q*3+2]-(double)target[q*3+2];
            bounded=bounded && sqrt(d0*d0+d1*d1+d2*d2)<=limit[q]+1e-4;
        }
        fprintf(stderr,
                "[ribbon selftest] (local-patch-neighbourhood) "
                "legal=%d->%d: %s\n",
                before,after,bounded?"ok":"FAIL");
        st_check(bounded,
                 "bounded patch neighbourhood preserves outer cells",&fails);
    }

    /* Large-region regression: allocation sizes must never pass through the
     * 32-bit `long` used by Win64.  Exercise the arithmetic without actually
     * reserving the multi-gigabyte block. */
    {
        size_t count = (size_t)INT32_MAX / sizeof(Sample) + 1;
        size_t bytes = rib_array_bytes(count, sizeof(Sample));
        int ok = bytes > (size_t)INT32_MAX &&
                 bytes == count * sizeof(Sample);
        fprintf(stderr, "[ribbon selftest] (alloc64) count=%zu bytes=%zu: %s\n",
                count, bytes, ok ? "ok" : "FAIL");
        st_check(ok, ">2 GiB array arithmetic remains 64-bit", &fails);
    }

    /* (0ab) A connected-component id is only a transitive domain, never a
     * continuation certificate.  Cross-chart matching must consume one direct
     * serialized relation; same-chart continuation remains intrinsic. */
    {
        int32_t edge_a[1] = { 2 }, edge_b[1] = { 7 };
        RibbonOpts graph_opts;
        RibbonOpts_default(&graph_opts);
        graph_opts.reference_chart = edge_a; /* enables strict graph mode */
        graph_opts.reference_relation_a = edge_a;
        graph_opts.reference_relation_b = edge_b;
        graph_opts.reference_relation_count = 1;
        int direct = direct_chart_relation_allowed(&graph_opts, 2, 7);
        int reverse = direct_chart_relation_allowed(&graph_opts, 7, 2);
        int transitive = direct_chart_relation_allowed(&graph_opts, 2, 9);
        int intrinsic = direct_chart_relation_allowed(&graph_opts, 2, 2);
        fprintf(stderr,
                "[ribbon selftest] (0ab) direct chart graph: direct=%d "
                "reverse=%d transitive=%d intrinsic=%d\n",
                direct, reverse, transitive, intrinsic);
        st_check(direct && reverse && !transitive && intrinsic,
                 "direct chart edges gate continuation", &fails);
    }

    /* (m1) Hard chain merges: a certified continuation joins two chains into
     * one logical chain; every gate rejects deterministically; the junction
     * survives chain-enumeration permutation. */
    {
        int32_t edge_a[1] = { 2 }, edge_b[1] = { 7 };
        Sample smp[4];
        Chain chn[2];
        SliceSet S;
        RibbonOpts o;
        RibMergeSet m;
        int accepted = 0, permuted = 0, gated = 1;
        memset(smp, 0, sizeof(smp));
        memset(chn, 0, sizeof(chn));
        memset(&S, 0, sizeof(S));
        RibbonOpts_default(&o);
        o.reference_chart = edge_a;   /* strict graph mode */
        o.reference_relation_a = edge_a;
        o.reference_relation_b = edge_b;
        o.reference_relation_count = 1;

        /* chain 0 (chart 2): u 96..100 ending at (0,0,100); chain 1 (chart 7):
         * u 101.5..140 starting 2 vox away.  Same island, same radius. */
        chn[0].first = 0; chn[0].count = 2; chn[0].source_chart = 2;
        chn[1].first = 2; chn[1].count = 2; chn[1].source_chart = 7;
        smp[0].u = 96.0;  smp[0].p[2] = 96.0;
        smp[1].u = 100.0; smp[1].p[2] = 100.0;
        smp[2].u = 101.5; smp[2].p[2] = 102.0;
        smp[3].u = 140.0; smp[3].p[2] = 140.0;
        for (int i = 0; i < 4; i++) {
            smp[i].r = 50.0;
            smp[i].phi = i < 2 ? 10.0 : 10.3;
            smp[i].chain = i < 2 ? 0 : 1;
        }
        S.smp = smp; S.n_smp = 4;
        S.chn = chn; S.n_chn = 2;
        S.nplanes = 1;

        accepted = rib_merge_build(arena, &S, &o, &m) == 0 &&
                   m.n_junctions == 1 && m.n_logical == 1 &&
                   m.logical[0] == m.logical[1] && m.ambiguous == 0;

        /* permuted enumeration: swap the two chains' identities */
        chn[0].first = 2; chn[0].source_chart = 7;
        chn[1].first = 0; chn[1].source_chart = 2;
        smp[0].chain = 1; smp[1].chain = 1;
        smp[2].chain = 0; smp[3].chain = 0;
        permuted = rib_merge_build(arena, &S, &o, &m) == 0 &&
                   m.n_junctions == 1 && m.logical[0] == m.logical[1];
        chn[0].first = 0; chn[0].source_chart = 2;
        chn[1].first = 2; chn[1].source_chart = 7;
        smp[0].chain = 0; smp[1].chain = 0;
        smp[2].chain = 1; smp[3].chain = 1;

        /* no serialized relation -> rejected_relation */
        o.reference_relation_count = 0;
        gated = gated && rib_merge_build(arena, &S, &o, &m) == 0 &&
                m.n_junctions == 0 && m.rejected_relation == 1;
        o.reference_relation_count = 1;
        /* carried-U disagreement -> rejected_du */
        smp[2].u = 112.0;
        gated = gated && rib_merge_build(arena, &S, &o, &m) == 0 &&
                m.n_junctions == 0 && m.rejected_du == 1;
        smp[2].u = 101.5;
        /* physical gap beyond seam reach -> rejected_gap */
        smp[2].p[1] = 8.0;
        gated = gated && rib_merge_build(arena, &S, &o, &m) == 0 &&
                m.n_junctions == 0 && m.rejected_gap == 1;
        smp[2].p[1] = 0.0;
        /* wrong wrap radially -> rejected_radial */
        smp[2].r = 54.0;
        gated = gated && rib_merge_build(arena, &S, &o, &m) == 0 &&
                m.n_junctions == 0 && m.rejected_radial == 1;
        smp[2].r = 50.0;
        /* lifted-phase flip -> rejected_phase */
        smp[2].phi = 12.0;
        gated = gated && rib_merge_build(arena, &S, &o, &m) == 0 &&
                m.n_junctions == 0 && m.rejected_phase == 1;
        smp[2].phi = 10.3;
        /* island mismatch -> rejected_island */
        chn[1].winding_island = 1;
        gated = gated && rib_merge_build(arena, &S, &o, &m) == 0 &&
                m.n_junctions == 0 && m.rejected_island == 1;
        chn[1].winding_island = 0;

        fprintf(stderr,
                "[ribbon selftest] (m1) chain merge: accepted=%s "
                "permuted=%s gated=%s\n",
                accepted ? "yes" : "NO", permuted ? "yes" : "NO",
                gated ? "yes" : "NO");
        st_check(accepted, "certified junction merges two chains", &fails);
        st_check(permuted, "junctions survive enumeration permutation", &fails);
        st_check(gated, "every junction gate rejects deterministically", &fails);
    }

    /* (m2+m5) Deterministic fusion on a certified overlap: the higher-support
     * chart owns the overlap in every row, the frontier sits at the end of its
     * measured support, the transition audit certifies it, two runs agree
     * byte-for-byte, and the exported provenance survives later arena use
     * (regression: the export used to be reclaimed by Arena_restore). */
    {
        int32_t edge_a[1] = { 2 }, edge_b[1] = { 7 };
        Sample smp[16];
        Chain chn[4];
        SliceSet S;
        RibbonOpts o;
        RibbonResult R1, R2;
        int fused = 0, identical = 1, survives = 1;
        memset(smp, 0, sizeof(smp));
        memset(chn, 0, sizeof(chn));
        memset(&S, 0, sizeof(S));
        memset(&R1, 0, sizeof(R1));
        memset(&R2, 0, sizeof(R2));
        RibbonOpts_default(&o);
        o.grid_u = 1.0f;
        o.slice_h = 2.0f;
        o.component_global = 1;
        o.ownership_construction = 1;
        o.reference_chart = edge_a;
        o.reference_relation_a = edge_a;
        o.reference_relation_b = edge_b;
        o.reference_relation_count = 1;

        /* two planes; per plane: chart 2 with 6 samples over u 0..10 and
         * chart 7 with 2 samples over u 8..18 shifted 0.5 vox in y.  Chart 2
         * carries the larger global support, so it owns the certified
         * junction overlap (columns 8..10) in BOTH rows. */
        for (int row = 0; row < 2; row++) {
            Chain *a = &chn[row * 2], *b = &chn[row * 2 + 1];
            a->first = (int32_t)(row * 8); a->count = 6;
            a->slice = row; a->source_chart = 2;
            b->first = (int32_t)(row * 8 + 6); b->count = 2;
            b->slice = row; b->source_chart = 7;
            for (int q = 0; q < 6; q++) {
                Sample *s = &smp[row * 8 + q];
                s->u = 2.0 * (double)q;
                s->p[0] = 2.0 * (double)row;
                s->p[2] = s->u;
                s->phi = s->u * 0.05;
                s->r = 50.0;
                s->tau[2] = 1.0;
                s->chain = (int32_t)(row * 2);
                s->slice = row;
            }
            for (int q = 0; q < 2; q++) {
                Sample *s = &smp[row * 8 + 6 + q];
                s->u = q ? 18.0 : 8.0;
                s->p[0] = 2.0 * (double)row;
                s->p[1] = 0.5;
                s->p[2] = s->u;
                s->phi = s->u * 0.05;
                s->r = 50.0;
                s->tau[2] = 1.0;
                s->chain = (int32_t)(row * 2 + 1);
                s->slice = row;
            }
        }
        S.smp = smp; S.n_smp = 16;
        S.chn = chn; S.n_chn = 4;
        S.nplanes = 2;
        fit_ribbon(arena, &S, 0.0, &o, NULL, &R1);
        fit_ribbon(arena, &S, 0.0, &o, NULL, &R2);
        fused = R1.grid_pos != NULL && R1.grid_chart_runs != NULL &&
                R1.grid_chart_run_count == 4 &&
                R1.grid_fuse_same_logical > 0 &&
                R1.grid_fuse_violations == 0 &&
                R1.grid_merge_junctions == 2 &&
                R1.grid_transition_certified == 2 &&
                R1.grid_transition_uncertified == 0;
        for (size_t r = 0; fused && r < 4; r++) {
            const RibbonChartRun *run = &R1.grid_chart_runs[r];
            size_t want_first = (r & 1) ? 11u : 0u;
            size_t want_last = (r & 1) ? 18u : 10u;
            int32_t want_chart = (r & 1) ? 7 : 2;
            if (run->first_col != want_first || run->last_col != want_last ||
                run->source_chart != want_chart)
                fused = 0;
        }
        identical = R2.grid_chart_run_count == R1.grid_chart_run_count;
        for (size_t r = 0; identical && r < R1.grid_chart_run_count; r++)
            if (memcmp(&R1.grid_chart_runs[r], &R2.grid_chart_runs[r],
                       sizeof(RibbonChartRun)) != 0)
                identical = 0;
        {   /* m5: the exported provenance must survive later arena use */
            unsigned char *scrub =
                (unsigned char *)ARENA_ALLOC(arena, (size_t)1 << 20);
            memset(scrub, 0xAB, (size_t)1 << 20);
            for (size_t r = 0; r < R2.grid_chart_run_count; r++) {
                const RibbonChartRun *run = &R2.grid_chart_runs[r];
                if (run->source_chart != 2 && run->source_chart != 7)
                    survives = 0;
                if (run->first_col > run->last_col || run->last_col > 18)
                    survives = 0;
            }
            survives = survives && R2.grid_chart_row_offsets != NULL &&
                       R2.grid_chart_row_offsets[R2.nk] ==
                           R2.grid_chart_run_count;
        }
        fprintf(stderr,
                "[ribbon selftest] (m2+m5) fusion: fused=%s identical=%s "
                "export_survives=%s (junctions=%zu same=%zu certified=%zu)\n",
                fused ? "yes" : "NO", identical ? "yes" : "NO",
                survives ? "yes" : "NO", R1.grid_merge_junctions,
                R1.grid_fuse_same_logical, R1.grid_transition_certified);
        st_check(fused,
                 "fusion resolves a certified overlap at end of support",
                 &fails);
        st_check(identical, "fusion is run-to-run identical", &fails);
        st_check(survives,
                 "exported chart runs survive later arena allocation",
                 &fails);
    }

    /* (m4) A one-column support gap is not an ownership-switch license: the
     * winner of a contested column is a row-invariant function of chart
     * support, so the same chart wins on both sides of its own gap, and a
     * higher-support competitor takes exactly its measured columns. */
    {
        int32_t edge_a[1] = { 2 }, edge_b[1] = { 7 };
        Sample smp[10];
        Chain chn[3];
        SliceSet S;
        RibbonOpts o;
        RibbonResult R;
        int low_ok = 0, high_ok = 0;
        memset(smp, 0, sizeof(smp));
        memset(chn, 0, sizeof(chn));
        memset(&S, 0, sizeof(S));
        RibbonOpts_default(&o);
        o.grid_u = 1.0f;
        o.slice_h = 2.0f;
        o.component_global = 1;
        o.ownership_construction = 1;
        o.reference_chart = edge_a;
        o.reference_relation_a = edge_a;
        o.reference_relation_b = edge_b;
        o.reference_relation_count = 1;

        /* chart 2 split across a 1-column hole: chains u 0..4 and 6..10;
         * chart 7 overlaps u 3..7 half a voxel away. */
        chn[0].first = 0; chn[0].count = 2; chn[0].source_chart = 2;
        chn[1].first = 2; chn[1].count = 2; chn[1].source_chart = 2;
        chn[2].first = 4; chn[2].count = 2; chn[2].source_chart = 7;
        smp[0].u = 0.0;  smp[1].u = 4.0;
        smp[2].u = 6.0;  smp[3].u = 10.0;
        smp[4].u = 3.0;  smp[5].u = 7.0;
        for (int i = 0; i < 6; i++) {
            smp[i].p[2] = smp[i].u;
            smp[i].p[1] = i >= 4 ? 0.5 : 0.0;
            smp[i].r = 50.0;
            smp[i].phi = smp[i].u * 0.05;
            smp[i].tau[2] = 1.0;
            smp[i].chain = i / 2;
        }
        S.smp = smp; S.n_smp = 6;
        S.chn = chn; S.n_chn = 3;
        S.nplanes = 1;

        memset(&R, 0, sizeof(R));
        fit_ribbon(arena, &S, 0.0, &o, NULL, &R);
        /* chart 2 support 4 > chart 7 support 2: runs 2|7|2 with the
         * competitor holding ONLY the unmeasured column 5 */
        low_ok = R.grid_chart_run_count == 3 &&
                 R.grid_chart_runs[0].source_chart == 2 &&
                 R.grid_chart_runs[0].last_col == 4 &&
                 R.grid_chart_runs[1].source_chart == 7 &&
                 R.grid_chart_runs[1].first_col == 5 &&
                 R.grid_chart_runs[1].last_col == 5 &&
                 R.grid_chart_runs[2].source_chart == 2 &&
                 R.grid_chart_runs[2].first_col == 6 &&
                 R.grid_fuse_violations == 0;

        /* raise chart 7's support above chart 2's: it takes exactly its
         * measured columns 3..7 as ONE run -- never a mid-support splice */
        chn[2].count = 6;
        {
            static const double bu[6] = {3.0, 3.8, 4.6, 5.4, 6.2, 7.0};
            for (int q = 0; q < 6; q++) {
                Sample *s = &smp[4 + q];
                s->u = bu[q];
                s->p[2] = s->u;
                s->p[1] = 0.5;
                s->p[0] = 0.0;
                s->r = 50.0;
                s->phi = s->u * 0.05;
                s->tau[2] = 1.0;
                s->chain = 2;
                s->slice = 0;
            }
        }
        S.n_smp = 10;
        memset(&R, 0, sizeof(R));
        fit_ribbon(arena, &S, 0.0, &o, NULL, &R);
        high_ok = R.grid_chart_run_count == 3 &&
                  R.grid_chart_runs[0].source_chart == 2 &&
                  R.grid_chart_runs[0].last_col == 2 &&
                  R.grid_chart_runs[1].source_chart == 7 &&
                  R.grid_chart_runs[1].first_col == 3 &&
                  R.grid_chart_runs[1].last_col == 7 &&
                  R.grid_chart_runs[2].source_chart == 2 &&
                  R.grid_chart_runs[2].first_col == 8 &&
                  R.grid_fuse_violations == 0;
        fprintf(stderr,
                "[ribbon selftest] (m4) gap is not a switch license: "
                "low=%s high=%s\n",
                low_ok ? "yes" : "NO", high_ok ? "yes" : "NO");
        st_check(low_ok,
                 "a support gap does not flip ownership across it", &fails);
        st_check(high_ok,
                 "a stronger competitor takes exactly its measured columns",
                 &fails);
    }

    /* (c0) Claims ownership, restored 2026-08-18: a whole welded component may
     * expose several physical wraps at the same solved-u column.  The
     * component-global selector needs a robust, monotone reference rather than
     * chain enumeration order. */
    {
        enum { NBIN = 5, NOBS = 5 };
        Sample smp[NBIN * NOBS];
        SliceSet S;
        memset(smp, 0, sizeof(smp));
        memset(&S, 0, sizeof(S));
        for (int j = 0; j < NBIN; j++) {
            static const double noise[NOBS] = {-50.0, -0.08, 0.0, 0.07, 60.0};
            for (int q = 0; q < NOBS; q++) {
                size_t i = (size_t)j * NOBS + (size_t)q;
                smp[i].u = (double)j;
                smp[i].phi = 1.5 * (double)j + noise[q];
                smp[i].chain = 0;
            }
        }
        S.smp = smp;
        S.n_smp = NBIN * NOBS;
        double *expect = component_phi_consensus(
            arena, &S, 0.0, 1.0, NBIN, 0);
        int robust = expect != NULL, monotone = expect != NULL;
        for (int j = 0; expect != NULL && j < NBIN; j++) {
            if (fabs(expect[j] - 1.5 * (double)j) > 0.1) robust = 0;
            if (j > 0 && expect[j] < expect[j-1]) monotone = 0;
        }
        fprintf(stderr, "[ribbon selftest] (c0) claims phi consensus: "
                "robust=%s monotone=%s\n", robust ? "yes" : "NO",
                monotone ? "yes" : "NO");
        st_check(robust, "claims phi consensus rejects wrap outliers", &fails);
        st_check(monotone, "claims phi consensus is monotone", &fails);
    }

    /* (c0b) Claims ownership, restored 2026-08-18: once the preceding slice
     * identifies a sheet, a competing wrap outside the physical gate must not
     * steal the next row just because its phi residual is smaller.  B has the
     * exact global phi and appears last, reproducing the former overwrite
     * failure. */
    {
        Sample smp[6];
        Chain chn[3];
        SliceSet S;
        RibbonOpts o;
        RibbonResult R;
        memset(smp, 0, sizeof(smp));
        memset(chn, 0, sizeof(chn));
        memset(&S, 0, sizeof(S));
        memset(&R, 0, sizeof(R));
        RibbonOpts_default(&o);
        o.grid_u = 1.0f;
        o.slice_h = 2.0f;
        o.component_global = 1;
        o.ownership_construction = 0;

        for (int c = 0; c < 3; c++) {
            chn[c].first = 2 * c;
            chn[c].count = 2;
            chn[c].slice = c == 0 ? 0 : 1;
            for (int q = 0; q < 2; q++) {
                int i = 2 * c + q;
                smp[i].u = q ? 4.0 : 0.0;
                smp[i].phi = c == 1 ? 5.0 : 0.0;
                smp[i].p[0] = c == 0 ? 0.0 : 2.0;
                smp[i].p[1] = c == 2 ? 10.0 : 0.0;
                smp[i].p[2] = smp[i].u;
                smp[i].tau[2] = 1.0;
                smp[i].chain = c;
                smp[i].slice = chn[c].slice;
            }
        }
        S.smp = smp; S.n_smp = 6;
        S.chn = chn; S.n_chn = 3;
        S.nplanes = 2;
        fit_ribbon(arena, &S, 0.0, &o, NULL, &R);
        int continuous = R.grid_pos != NULL && R.nk == 2 && R.nu >= 5;
        for (size_t j = 0; continuous && j < 5; j++) {
            double y = (double)R.grid_pos[(R.nu + j) * 3 + 1];
            if (!isfinite(y) || fabs(y) > 0.25) continuous = 0;
        }
        fprintf(stderr, "[ribbon selftest] (c0b) claims continuity: "
                "conflicts=%zu stayed_on_sheet=%s\n", R.grid_claim_conflicts,
                continuous ? "yes" : "NO");
        st_check(R.grid_claim_conflicts == 5,
                 "claims continuity exercised all claims", &fails);
        st_check(continuous,
                 "claims path selector does not hop wraps", &fails);
    }

    /* (c0ba) A long source curve and a shorter competing curve can both be
     * physically continuous and use the same source chart.  Cellwise scoring
     * chose the short curve in its middle interval, yielding A->B->A despite A
     * covering both endpoints.  The row path must keep A, while still allowing
     * a transition that is forced by the end of measured support. */
    {
        enum { NCOL = 9, MAX_CLAIM = 14 };
        RibGridClaim claim[MAX_CLAIM];
        RibGridClaim forced[NCOL];
        size_t off[NCOL + 1], forced_off[NCOL + 1];
        float previous[NCOL * 3];
        int32_t selected[NCOL], forced_selected[NCOL];
        size_t nclaim = 0;
        int coherent = 1, forced_handoff = 1;
        memset(claim, 0, sizeof(claim));
        memset(forced, 0, sizeof(forced));
        for (int j = 0; j < NCOL; j++) {
            previous[j*3] = 0.0f;
            previous[j*3 + 1] = 0.0f;
            previous[j*3 + 2] = (float)j;
            off[j] = nclaim;
            claim[nclaim].p[1] = 2.0;
            claim[nclaim].p[2] = (double)j;
            claim[nclaim].phi_score = 1.0;
            claim[nclaim].chain = 101;
            claim[nclaim].source_chart = 7;
            nclaim++;
            if (j >= 2 && j <= 6) {
                claim[nclaim].p[1] = 0.0;
                claim[nclaim].p[2] = (double)j;
                claim[nclaim].chain = 202;
                claim[nclaim].source_chart = 7;
                nclaim++;
            }
            off[j + 1] = nclaim;

            forced_off[j] = (size_t)j;
            forced[j].p[2] = (double)j;
            forced[j].chain = j < 4 ? 303 : 404;
            forced[j].source_chart = 7;
            forced_off[j + 1] = (size_t)j + 1;
        }
        rib_claim_select_row(
            arena, claim, off, NCOL, previous, NULL, 1, selected);
        rib_claim_select_row(
            arena, forced, forced_off, NCOL, previous, NULL, 1,
            forced_selected);
        for (int j = 0; j < NCOL; j++) {
            if (selected[j] < 0 ||
                claim[(size_t)selected[j]].chain != 101)
                coherent = 0;
            if (forced_selected[j] < 0 ||
                forced[(size_t)forced_selected[j]].chain !=
                    (j < 4 ? 303 : 404))
                forced_handoff = 0;
        }
        fprintf(stderr,
                "[ribbon selftest] (c0ba) claimant row path: "
                "contained_overlap=%s forced_handoff=%s\n",
                coherent ? "yes" : "NO",
                forced_handoff ? "yes" : "NO");
        st_check(coherent,
                 "contained claimant cannot create an A-B-A mosaic", &fails);
        st_check(forced_handoff,
                 "claimant path permits support-forced handoff", &fails);
    }

    /* (c0c) Conflict erasure: a doubly claimed middle row contributes no
     * Dirichlet data.  The fitted ribbon may still carry a finite variational
     * continuation from the clean rows above and below, but its support mask
     * must remain zero so downstream fitting can never promote it back to an
     * observation. */
    {
        Sample smp[8];
        Chain chn[4];
        SliceSet S;
        RibbonOpts o;
        RibbonResult R;
        memset(smp, 0, sizeof(smp));
        memset(chn, 0, sizeof(chn));
        memset(&S, 0, sizeof(S));
        memset(&R, 0, sizeof(R));
        RibbonOpts_default(&o);
        o.grid_u = 1.0f;
        o.slice_h = 2.0f;
        o.component_global = 1;
        o.ownership_construction = 0;
        o.discard_conflicting_claims = 1;

        for (int c = 0; c < 4; c++) {
            chn[c].first = 2 * c;
            chn[c].count = 2;
            chn[c].slice = c == 0 ? 0 : (c == 3 ? 2 : 1);
            for (int q = 0; q < 2; q++) {
                int i = 2 * c + q;
                smp[i].u = q ? 4.0 : 0.0;
                smp[i].phi = 0.0;
                smp[i].p[0] = 2.0 * (double)chn[c].slice;
                smp[i].p[1] = c == 2 ? 10.0 : 0.0;
                smp[i].p[2] = smp[i].u;
                smp[i].tau[2] = 1.0;
                smp[i].chain = c;
                smp[i].slice = chn[c].slice;
            }
        }
        S.smp = smp; S.n_smp = 8;
        S.chn = chn; S.n_chn = 4;
        S.nplanes = 3;
        fit_ribbon(arena, &S, 0.0, &o, NULL, &R);
        int generated = R.grid_pos != NULL && R.grid_valid != NULL &&
                        R.nk == 3 && R.nu >= 5;
        for (size_t j = 0; generated && j < 5; j++) {
            size_t s = R.nu + j;
            double y = (double)R.grid_pos[s * 3 + 1];
            if (R.grid_valid[s] || !isfinite(y) || fabs(y) > 0.25)
                generated = 0;
        }
        fprintf(stderr, "[ribbon selftest] (c0c) conflict erasure: "
                "cells=%zu discarded=%zu generated_unsupported=%s\n",
                R.grid_claim_conflict_cells,
                R.grid_claim_discarded_cells,
                generated ? "yes" : "NO");
        st_check(R.grid_claim_conflicts == 5 &&
                 R.grid_claim_conflict_cells == 5 &&
                 R.grid_claim_discarded_cells == 5 &&
                 R.grid_claim_replaced == 0,
                 "conflict erasure discards every multiply claimed cell",
                 &fails);
        st_check(generated,
                 "conflict continuation stays finite but unsupported",
                 &fails);
    }

    /* (c0bb) A material component can contain two genuine reconstruction
     * lanes.  They overlap in solved U on every row, but are farther apart in
     * 3D than the wrap gate and each has its own persistent cross-row
     * continuation.  The branch splitter must keep both histories instead of
     * forcing either one through the single raster claimant slot. */
    {
        enum { NROW = 3, NLANE = 2, NCHAIN = NROW * NLANE };
        Sample smp[NCHAIN * 2];
        Chain chn[NCHAIN];
        SliceSet S;
        RibbonOpts o;
        RibbonResult R;
        ChainRelationGraph graph;
        int split_ok = 1, stable_ok = 1;
        memset(smp, 0, sizeof(smp));
        memset(chn, 0, sizeof(chn));
        memset(&S, 0, sizeof(S));
        memset(&R, 0, sizeof(R));
        memset(&graph, 0, sizeof(graph));
        RibbonOpts_default(&o);
        o.grid_u = 1.0f;
        o.slice_h = 2.0f;
        o.component_global = 1;

        for (int row = 0; row < NROW; row++) {
            for (int lane = 0; lane < NLANE; lane++) {
                int c = row * NLANE + lane;
                chn[c].first = 2 * c;
                chn[c].count = 2;
                chn[c].slice = row;
                chn[c].mesh_comp = 0;
                chn[c].source_chart = 0;
                chn[c].winding_island = 0;
                chn[c].reconstruction_component = 0;
                for (int q = 0; q < 2; q++) {
                    int i = 2 * c + q;
                    smp[i].u = q ? 4.0 : 0.0;
                    smp[i].phi = 0.0;
                    smp[i].p[0] = 2.0 * (double)row;
                    smp[i].p[1] = 10.0 * (double)lane;
                    smp[i].p[2] = smp[i].u;
                    smp[i].r = 50.0 + 10.0 * (double)lane;
                    smp[i].tau[2] = 1.0;
                    smp[i].chain = c;
                    smp[i].slice = row;
                }
            }
        }
        for (int row = 1; row < NROW; row++) {
            chain_relation_graph_add(&graph,
                                     (row - 1) * NLANE,
                                     row * NLANE, 10.0);
            chain_relation_graph_add(&graph,
                                     (row - 1) * NLANE + 1,
                                     row * NLANE + 1, 10.0);
        }
        S.smp = smp; S.n_smp = NCHAIN * 2;
        S.chn = chn; S.n_chn = NCHAIN;
        S.nplanes = NROW;
        if (rib_assign_reconstruction_components(
                arena, &S, &o, &graph, &R) != 0)
            split_ok = stable_ok = 0;
        if (split_ok) {
            split_ok = R.grid_reconstruction_components == 2 &&
                       chn[0].reconstruction_component !=
                           chn[1].reconstruction_component;
            for (int row = 1; row < NROW; row++) {
                if (chn[row * NLANE].reconstruction_component !=
                        chn[0].reconstruction_component ||
                    chn[row * NLANE + 1].reconstruction_component !=
                        chn[1].reconstruction_component)
                    stable_ok = 0;
            }
        }
        fprintf(stderr,
                "[ribbon selftest] (c0bb) branch components: "
                "split=%s stable=%s components=%zu conflicts=%zu\n",
                split_ok ? "yes" : "NO", stable_ok ? "yes" : "NO",
                R.grid_reconstruction_components,
                R.grid_branch_conflicts);
        st_check(split_ok,
                 "overlapping physical lanes become separate components",
                 &fails);
        st_check(stable_ok,
                 "branch component identity persists across rows", &fails);
        chain_relation_graph_dispose(&graph);
    }

    /* (0b) A competing wrap 10 voxels away is a certificate violation, not a
     * per-cell decision: the row-invariant key keeps the row on one sheet and
     * the violation is counted and reported, never optimized away. */
    {
        Sample smp[6];
        Chain chn[3];
        SliceSet S;
        RibbonOpts o;
        RibbonResult R;
        memset(smp, 0, sizeof(smp));
        memset(chn, 0, sizeof(chn));
        memset(&S, 0, sizeof(S));
        memset(&R, 0, sizeof(R));
        RibbonOpts_default(&o);
        o.grid_u = 1.0f;
        o.slice_h = 2.0f;
        o.component_global = 1;
        o.ownership_construction = 1;

        /* row 0: only sheet A; row 1: A first, then a 10-voxel-away wrap B
         * claiming the same solved-u columns (a gauge collision). */
        for (int c = 0; c < 3; c++) {
            chn[c].first = 2 * c;
            chn[c].count = 2;
            chn[c].slice = c == 0 ? 0 : 1;
            for (int q = 0; q < 2; q++) {
                int i = 2 * c + q;
                smp[i].u = q ? 4.0 : 0.0;
                smp[i].phi = c == 1 ? 5.0 : 0.0;
                smp[i].p[0] = c == 0 ? 0.0 : 2.0;
                smp[i].p[1] = c == 2 ? 10.0 : 0.0;
                smp[i].p[2] = smp[i].u;
                smp[i].tau[2] = 1.0;
                smp[i].chain = c;
                smp[i].slice = chn[c].slice;
            }
        }
        S.smp = smp; S.n_smp = 6;
        S.chn = chn; S.n_chn = 3;
        S.nplanes = 2;
        fit_ribbon(arena, &S, 0.0, &o, NULL, &R);
        int continuous = R.grid_pos != NULL && R.nk == 2 && R.nu >= 5;
        for (size_t j = 0; continuous && j < 5; j++) {
            double y = (double)R.grid_pos[(R.nu + j) * 3 + 1];
            if (!isfinite(y) || fabs(y) > 0.25) continuous = 0;
        }
        fprintf(stderr, "[ribbon selftest] (0b) component-global continuity: "
                "conflicts=%zu violations=%zu max=%.1f stayed_on_sheet=%s\n",
                R.grid_claim_conflicts, R.grid_fuse_violations,
                R.grid_fuse_violation_max, continuous ? "yes" : "NO");
        st_check(R.grid_claim_conflicts == 5,
                 "component-global continuity exercised all claims", &fails);
        st_check(continuous, "component-global selector does not hop wraps", &fails);
        st_check(R.grid_fuse_violations == 5 &&
                 R.grid_merge_junctions == 0 &&
                 fabs(R.grid_fuse_violation_max - 10.0) < 0.5,
                 "wrap-gate losers are reported as violations", &fails);
    }

    /* (0c) Unsupported grid width is not an input to the fit.  Exercise both a
     * short bridged U gap (whose P/T slots are deliberately NaN-poisoned) and a
     * detached island whose old cover gauge adds 100 empty columns. */
    {
        Sample smp[12];
        Chain chn[6];
        SliceSet S;
        RibbonOpts compact_opts, cover_opts;
        RibbonResult compact, cover;
        memset(smp, 0, sizeof(smp));
        memset(chn, 0, sizeof(chn));
        memset(&S, 0, sizeof(S));
        memset(&compact, 0, sizeof(compact));
        memset(&cover, 0, sizeof(cover));
        RibbonOpts_default(&compact_opts);
        compact_opts.grid_u = 1.0f;
        compact_opts.slice_h = 2.0f;
        compact_opts.component_global = 1;
        cover_opts = compact_opts;
        cover_opts.fit_cover_width = 1;

        for (int row = 0; row < 2; row++) {
            static const double ulo[3] = {0.0, 6.0, 20.0};
            for (int part = 0; part < 3; part++) {
                int c = row * 3 + part;
                chn[c].first = c * 2;
                chn[c].count = 2;
                chn[c].slice = row;
                chn[c].winding_island = part == 2 ? 1 : 0;
                chn[c].reconstruction_component =
                    chn[c].winding_island;
                chn[c].atlas_shift = part == 2 ? -100.0 : 0.0;
                for (int q = 0; q < 2; q++) {
                    int i = c * 2 + q;
                    smp[i].u = ulo[part] + 2.0 * (double)q;
                    smp[i].phi = smp[i].u;
                    smp[i].p[0] = 2.0 * (double)row;
                    smp[i].p[1] = part == 2 ? 20.0 : 0.0;
                    smp[i].p[2] = smp[i].u;
                    smp[i].tau[2] = 1.0;
                    smp[i].chain = c;
                    smp[i].slice = row;
                }
            }
        }
        S.smp = smp; S.n_smp = 12;
        S.chn = chn; S.n_chn = 6;
        S.nplanes = 2;
        fit_ribbon(arena, &S, 0.0, &compact_opts, NULL, &compact);
        fit_ribbon(arena, &S, 0.0, &cover_opts, NULL, &cover);

        int gap_finite = compact.grid_pos != NULL && compact.nu > 8;
        int invariant = gap_finite && cover.grid_pos != NULL &&
                        cover.nu > compact.nu && compact.nk == cover.nk;
        for (size_t k = 0; gap_finite && k < compact.nk; k++)
            for (size_t j = 3; j <= 5; j++)
                if (!isfinite((double)compact.grid_pos[(k*compact.nu+j)*3]))
                    gap_finite = 0;
        for (size_t k = 0; invariant && k < compact.nk; k++) {
            for (size_t j = 0; j < compact.nu; j++) {
                const float *a = &compact.grid_pos[(k*compact.nu+j)*3];
                const float *b = &cover.grid_pos[(k*cover.nu+j)*3];
                int av = isfinite((double)a[0]), bv = isfinite((double)b[0]);
                if (av != bv || (av && memcmp(a, b, 3*sizeof(float)) != 0)) {
                    invariant = 0;
                    break;
                }
            }
        }
        fprintf(stderr,
                "[ribbon selftest] (0c) fit-width invariance: compact=%zu "
                "cover=%zu gap_finite=%s identical=%s\n",
                compact.nu, cover.nu, gap_finite ? "yes" : "NO",
                invariant ? "yes" : "NO");
        st_check(gap_finite,
                 "short unsupported U gap fits without NaN contamination",
                 &fails);
        st_check(invariant,
                 "unsupported trailing U width cannot change fitted cells",
                 &fails);
    }

    /* (0d) Relation-graph islands need distinct universal-cover identities,
     * but their serial cover gaps are not material arclength.  Compact only the
     * solved metric bounding intervals, preserving every within-island du. */
    {
        Sample smp[6];
        Chain chn[3];
        SliceSet S;
        RibbonResult R;
        static const double base[3] = {0.0, 100000.0, 2500000.0};
        static const double span[3] = {10.0, 20.0, 5.0};
        memset(smp, 0, sizeof(smp));
        memset(chn, 0, sizeof(chn));
        memset(&S, 0, sizeof(S));
        memset(&R, 0, sizeof(R));
        for (int g = 0; g < 3; g++) {
            chn[g].first = 2 * g;
            chn[g].count = 2;
            chn[g].winding_island = g;
            chn[g].reconstruction_component = g;
            for (int q = 0; q < 2; q++) {
                int i = 2 * g + q;
                smp[i].chain = g;
                smp[i].u = base[g] + (q ? span[g] : 0.0);
            }
        }
        S.smp = smp; S.n_smp = 6;
        S.chn = chn; S.n_chn = 3;
        pack_metric_islands(arena, &S, 4.0,
                            NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, &R);
        int metric_exact = 1;
        for (int g = 0; g < 3; g++)
            if (fabs((smp[2*g+1].u - smp[2*g].u) - span[g]) > 1e-9)
                metric_exact = 0;
        double atlas_span = smp[5].u - smp[0].u;
        double expected_span = 10.0 + 40.0 + 20.0 + 40.0 + 5.0;
        int shift_recovers =
            fabs((smp[2].u - chn[1].atlas_shift) - base[1]) < 1e-9 &&
            fabs((smp[5].u - chn[2].atlas_shift) -
                 (base[2] + span[2])) < 1e-9;
        fprintf(stderr,
                "[ribbon selftest] (0d) metric-island atlas: span=%.1f "
                "saved=%.1f exact=%s shift_recovers=%s\n",
                atlas_span, R.uv_atlas_pack_saved,
                metric_exact ? "yes" : "NO",
                shift_recovers ? "yes" : "NO");
        st_check(metric_exact,
                 "metric-island packing preserves every intrinsic du", &fails);
        st_check(fabs(atlas_span - expected_span) < 1e-9 &&
                 R.uv_atlas_islands == 3 &&
                 R.uv_atlas_packed_islands == 2 &&
                 R.uv_atlas_pack_saved > 2000000.0,
                 "fake cover gaps collapse to raster-safe atlas gutters",
                 &fails);
        st_check(shift_recovers,
                 "atlas_shift recovers the pre-pack cover coordinate",
                 &fails);
    }

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
        st_check(R.duds_err_max < 3e-3,
                 "strip u tracks arc length within 0.3%", &fails);
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
            if (R.grid_valid != NULL && R.grid_pos != NULL &&
                R.nu > 4 && R.nk > 4) {
              for (size_t kk = 2; kk + 2 < R.nk && checked < 150; kk += 3) {
                for (size_t jj = 2; jj + 2 < R.nu && checked < 150; jj += 17) {
                    if (!R.grid_valid[kk*R.nu + jj]) continue;
                    const float *gp = &R.grid_pos[(kk*R.nu + jj)*3];
                    if (R.grid_present != NULL &&
                        !R.grid_present[kk*R.nu + jj]) continue;
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
        float *v = RIB_ALLOC_ARRAY(arena, float, nvv * 3);
        int32_t *f = RIB_ALLOC_ARRAY(
            arena, int32_t, (size_t)(NU-1)*(NH-1)*NW*2*3);
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
        float *v = RIB_ALLOC_ARRAY(arena, float, nvv * 3);
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
        float *v = RIB_ALLOC_ARRAY(arena, float, nvv * 3);
        int32_t *f = RIB_ALLOC_ARRAY(
            arena, int32_t, (size_t)(NU-1)*(NH-1)*2*2*3);
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
        float *v = RIB_ALLOC_ARRAY(arena, float, nvv * 3);
        int32_t *f = RIB_ALLOC_ARRAY(
            arena, int32_t, (size_t)(NU-1)*(NH-1)*3*2*3);
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
        float *v = RIB_ALLOC_ARRAY(arena, float, nvv * 3);
        int32_t *f = RIB_ALLOC_ARRAY(arena, int32_t, nfc_cap * 3);
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

    /* (4d) FIXED-TOPOLOGY METRIC PROJECTION: with a carried reference frame,
     * the parameterization is u = s + median(u_ref - s) per slice chain and
     * NOTHING else -- no chain-gauge solve, no Stage-C solve, no interval
     * stitching, no phi->u registration, no orientation flip.  The carried
     * placement of every component must survive exactly, INCLUDING a detached
     * arc whose carried interval lies mid-band inside the spiral's span (the
     * old gauge machinery re-packed or exiled exactly such charts). */
    {
        double r0 = 15.0, g = 2.0, turns = 3.0, Hgt = 40.0, cy = 100.0, cx = 50.0;
        int nphi = 600, nh = 21;
        float *sv; int32_t *sf; size_t snv, snf;
        st_build_spiral(arena, nphi, nh, r0, g, turns, Hgt, cy, cx, NULL, NULL,
                        &sv, &snv, &sf, &snf);
        enum { AN = 40, AH = 9 };
        double a_th0[2] = { 0.3 * M_PI, 1.2 * M_PI };
        double a_th1[2] = { 0.7 * M_PI, 1.6 * M_PI };
        size_t nvv = snv + (size_t)AN * AH * 2;
        size_t nfc_cap = snf + (size_t)(AN-1) * (AH-1) * 2 * 2;
        float *v = RIB_ALLOC_ARRAY(arena, float, nvv * 3);
        int32_t *f = RIB_ALLOC_ARRAY(arena, int32_t, nfc_cap * 3);
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
        /* carried frame: per-row cumulative XYZ arclength; arc A deliberately
         * carried MID-BAND, arc B beyond the spiral span */
        double *uref = RIB_ALLOC_ARRAY(arena, double, nvv);
        float *phir = RIB_ALLOC_ARRAY(arena, float, nvv);
        double amax = turns * 2.0 * M_PI, sp_umax = 0.0;
        for (int j = 0; j < nh; j++) {
            double cum = 0.0;
            for (int i = 0; i < nphi; i++) {
                size_t idx = (size_t)j * (size_t)nphi + (size_t)i;
                if (i > 0) {
                    size_t p = idx - 1;
                    double dz = (double)v[idx*3+0] - (double)v[p*3+0];
                    double dy = (double)v[idx*3+1] - (double)v[p*3+1];
                    double dx = (double)v[idx*3+2] - (double)v[p*3+2];
                    cum += sqrt(dz*dz + dy*dy + dx*dx);
                }
                uref[idx] = cum;
                phir[idx] = (float)(amax * (double)i / (double)(nphi - 1));
                if (cum > sp_umax) sp_umax = cum;
            }
        }
        for (int a = 0; a < 2; a++) {
            size_t base = snv + (size_t)a * AN * AH;
            double off = a == 0 ? 0.35 * sp_umax : 1.02 * sp_umax + 30.0;
            for (int j = 0; j < AH; j++) {
                double cum = 0.0;
                for (int i = 0; i < AN; i++) {
                    size_t idx = base + (size_t)j * AN + (size_t)i;
                    if (i > 0) {
                        size_t p = idx - 1;
                        double dz = (double)v[idx*3+0] - (double)v[p*3+0];
                        double dy = (double)v[idx*3+1] - (double)v[p*3+1];
                        double dx = (double)v[idx*3+2] - (double)v[p*3+2];
                        cum += sqrt(dz*dz + dy*dy + dx*dx);
                    }
                    uref[idx] = off + cum;
                    double th = a_th0[a] + (a_th1[a] - a_th0[a]) * (double)i / (double)(AN-1);
                    phir[idx] = (float)(th + turns * 2.0 * M_PI);
                }
            }
        }
        /* Run three ways: the physical carried frame; a uniformly COMPRESSED
         * frame (rate 0.55, the 10x regime) the common map must decompress;
         * and a frame with a +40 BLOCK OFFSET on spiral rows 8..12 that the
         * consistency solve must pull back through the vertical drift terms
         * (the "should have caught it" regression).  Assertions always
         * compare against the physical frame. */
        double *uref_run = RIB_ALLOC_ARRAY(arena, double, nvv);
        for (int sc = 0; sc < 3; sc++) {
        double scale = sc == 1 ? 0.55 : 1.0;
        for (size_t i = 0; i < nvv; i++) uref_run[i] = scale * uref[i];
        if (sc == 2)
            for (int j = 8; j <= 12; j++)
                for (int i2 = 0; i2 < nphi; i2++)
                    uref_run[(size_t)j * (size_t)nphi + (size_t)i2] += 40.0;
        RibbonOpts o; RibbonOpts_default(&o);
        o.axis_point[1] = (float)cy; o.axis_point[2] = (float)cx;
        o.reference_u = uref_run;
        o.solve_reference_u = 1;
        o.metric_project_only = 1;
        /* Exercise the production branch that reconstructs components.  Its
         * atlas pack must never translate an observed cover certificate. */
        o.component_global = 1;
        o.stitch_solve_gauges = 1;
        o.reference_phi = phir;
        o.reference_phi_authoritative = 1;
        o.fit_ribbon = 0;
        o.direct_ribbon = 0;
        RibbonResult R;
        memset(&R, 0, sizeof R);
        int rc = Ribbon_run(arena, v, nvv, f, nfc, &o, &R);
        size_t finite_uv = 0;
        double med[3] = {0,0,0}, p90[3] = {0,0,0};
        double c_lo[3] = {1e300,1e300,1e300};
        double c_hi[3] = {-1e300,-1e300,-1e300};
        if (rc == 0 && R.uv != NULL) {
            size_t lo[3] = { 0, snv, snv + (size_t)AN*AH };
            size_t hi[3] = { snv, snv + (size_t)AN*AH, nvv };
            double *buf = RIB_ALLOC_ARRAY(arena, double, nvv);
            for (size_t i = 0; i < nvv; i++)
                if (isfinite((double)R.uv[i*2]) &&
                    isfinite((double)R.uv[i*2+1]))
                    finite_uv++;
            for (int cpt = 0; cpt < 3; cpt++) {
                size_t m = 0;
                for (size_t i = lo[cpt]; i < hi[cpt]; i++) {
                    double uu = (double)R.uv[i*2];
                    buf[m++] = uu - uref[i];
                    if (uu < c_lo[cpt]) c_lo[cpt] = uu;
                    if (uu > c_hi[cpt]) c_hi[cpt] = uu;
                }
                med[cpt] = select_median_dbl(buf, m);
                for (size_t i = lo[cpt], q = 0; i < hi[cpt]; i++, q++)
                    buf[q] = fabs(((double)R.uv[i*2] - uref[i]) - med[cpt]);
                qsort(buf, m, sizeof(double), cmp_dbl);
                p90[cpt] = buf[(9 * m) / 10];
            }
        }
        double lenA = (a_th1[0] - a_th0[0])
                      * (r0 + g * (0.5*(a_th0[0]+a_th1[0]) + turns*2.0*M_PI));
        double lenB = (a_th1[1] - a_th0[1])
                      * (r0 + g * (0.5*(a_th0[1]+a_th1[1]) + turns*2.0*M_PI));
        fprintf(stderr,
                "[ribbon selftest] (4d) metric projection: rc=%d finite=%zu/%zu "
                "duds_max=%.2g align_rms=%.3f island_shift_max=%.1f "
                "gauge_comps=%zu\n"
                "[ribbon selftest]     spiral med=%.3f u=[%.1f,%.1f]  "
                "arcA u=[%.1f,%.1f] arcB u=[%.1f,%.1f]  "
                "spread p90 = %.3f/%.3f/%.3f\n",
                rc, finite_uv, nvv, R.duds_err_max, R.strip_align_rms,
                R.reg_max_shift, R.strip_gauge_components,
                med[0], c_lo[0], c_hi[0], c_lo[1], c_hi[1],
                c_lo[2], c_hi[2], p90[0], p90[1], p90[2]);
        st_check(rc == 0, "metric projection rc", &fails);
        st_check(finite_uv == nvv, "metric projection all UV finite", &fails);
        st_check(rc == 0 && R.duds_err_mean < 0.01,
                 "metric projection: chains near-metric (LS consensus)",
                 &fails);
        st_check(rc == 0 && R.strip_gauge_components == 0,
                 "metric projection: no gauge stitch ran", &fails);
        st_check(rc == 0 && fabs(med[0]) < 2.5 && p90[0] < 3.0,
                 "metric projection: main frame keeps the carried placement",
                 &fails);
        st_check(rc == 0 && R.reg_max_shift > 1.0 &&
                 c_lo[1] > c_hi[0] - 15.0 && c_lo[2] > c_hi[1] - 15.0,
                 "metric projection: detached arcs placed by carried phase "
                 "beyond the spiral (mid-band pack corrected)", &fails);
        st_check(rc == 0 &&
                 fabs((c_hi[1]-c_lo[1]) - lenA) / lenA < 0.15 &&
                 fabs((c_hi[2]-c_lo[2]) - lenB) / lenB < 0.15 &&
                 p90[1] < 3.0 && p90[2] < 3.0,
                 "metric projection: placed arcs stay exactly metric and "
                 "coherent", &fails);
        if (sc == 2 && rc == 0 && R.uv != NULL) {
            /* the offset band must be pulled back onto the physical frame */
            double *bb = RIB_ALLOC_ARRAY(arena, double, nvv);
            size_t nb2 = 0;
            for (int j = 8; j <= 12; j++)
                for (int i2 = 0; i2 < nphi; i2++) {
                    size_t idx = (size_t)j * (size_t)nphi + (size_t)i2;
                    bb[nb2++] = (double)R.uv[idx*2] - uref[idx];
                }
            double med_band = select_median_dbl(bb, nb2);
            fprintf(stderr,
                    "[ribbon selftest]     block-offset residual: "
                    "band-med=%.2f vs spiral-med=%.2f\n", med_band, med[0]);
            st_check(fabs(med_band - med[0]) < 1.5,
                     "consistency solve pulls a +40 carried block offset "
                     "back (<1.5 vox)", &fails);
        }
        fprintf(stderr,
                "[ribbon selftest]     (carried variant %d done)\n", sc);
        }
    }

    /* (4e) WARPED-CHART CONTACT regression: strip B abuts strip A in 3D but
     * its carried frame is stretched x1.3, so rigid placement leaves its
     * contact edge ~15 vox off (the measured 4x giant-offset anatomy).  The
     * consistency solve's gated cross-chart contact terms must pull B into
     * agreement, while a stacked ply P (normal-separated duplicate over A)
     * must NOT be tied and keeps its own placement. */
    {
        enum { AC = 61, BC = 50, PC = 41, RR = 21 };
        size_t nvA = (size_t)AC * RR, nvB = (size_t)BC * RR;
        size_t nvP = (size_t)PC * RR;
        size_t nvv = nvA + nvB + nvP;
        float *v = RIB_ALLOC_ARRAY(arena, float, nvv * 3);
        int32_t *f = RIB_ALLOC_ARRAY(
            arena, int32_t,
            ((size_t)(AC-1) + (BC-1) + (PC-1)) * (RR-1) * 6);
        double *uref = RIB_ALLOC_ARRAY(arena, double, nvv);
        size_t nfc = 0;
        /* ply decoy at 7 vox normal separation = the physical minimum
         * inter-wrap clearance (CUT_GAP_DEPTH); anything closer than the
         * perp gate is a hole-rim continuation, not a ply (measured 4x:
         * real rim ties at perp 2.8-3.0 with tdot=1.00) */
        struct { size_t base; int cols; double x0, y, u0, rate; } st[3] = {
            { 0,        AC,  0.0, 0.0,   0.0, 1.0 },
            { 0,        BC, 121.0, 0.0, 121.0, 1.3 },
            { 0,        PC, 30.0, 7.0, 500.0, 1.0 }
        };
        st[1].base = nvA;
        st[2].base = nvA + nvB;
        for (int sidx = 0; sidx < 3; sidx++) {
            for (int j = 0; j < RR; j++)
                for (int i = 0; i < st[sidx].cols; i++) {
                    size_t idx = st[sidx].base + (size_t)j *
                                 (size_t)st[sidx].cols + (size_t)i;
                    double x = st[sidx].x0 + 2.0 * i;
                    v[idx*3+0] = (float)(2.0 * j);
                    v[idx*3+1] = (float)st[sidx].y;
                    v[idx*3+2] = (float)x;
                    uref[idx] = st[sidx].u0 +
                                st[sidx].rate * (x - st[sidx].x0);
                }
            for (int j = 0; j + 1 < RR; j++)
                for (int i = 0; i + 1 < st[sidx].cols; i++) {
                    int32_t a0 = (int32_t)(st[sidx].base +
                        (size_t)j * (size_t)st[sidx].cols + (size_t)i);
                    int32_t b0 = a0 + 1;
                    int32_t c0 = a0 + st[sidx].cols;
                    int32_t d0 = c0 + 1;
                    f[nfc*3+0]=a0; f[nfc*3+1]=b0; f[nfc*3+2]=c0; nfc++;
                    f[nfc*3+0]=b0; f[nfc*3+1]=d0; f[nfc*3+2]=c0; nfc++;
                }
        }
        RibbonOpts o; RibbonOpts_default(&o);
        o.reference_u = uref;
        o.solve_reference_u = 1;
        o.metric_project_only = 1;
        o.stitch_solve_gauges = 1;
        o.fit_ribbon = 0;
        o.direct_ribbon = 0;
        RibbonResult R;
        memset(&R, 0, sizeof R);
        int rc = Ribbon_run(arena, v, nvv, f, nfc, &o, &R);
        double contact_gap = 1e300, span_b = 0.0, ply_sep = 0.0;
        if (rc == 0 && R.uv != NULL) {
            double *buf = RIB_ALLOC_ARRAY(arena, double, nvv);
            size_t m = 0;
            for (int j = 0; j < RR; j++)
                buf[m++] = (double)R.uv[(nvA + (size_t)j * BC) * 2]
                         - (double)R.uv[((size_t)j * AC + AC - 1) * 2];
            contact_gap = select_median_dbl(buf, m);
            double blo = 1e300, bhi = -1e300;
            for (size_t i = nvA; i < nvA + nvB; i++) {
                double uu = (double)R.uv[i*2];
                if (uu < blo) blo = uu;
                if (uu > bhi) bhi = uu;
            }
            span_b = bhi - blo;
            m = 0;
            for (int j = 0; j < RR; j++)
                buf[m++] = (double)R.uv[(nvA + nvB + (size_t)j * PC) * 2]
                         - (double)R.uv[((size_t)j * AC + 15) * 2];
            ply_sep = select_median_dbl(buf, m);
        }
        fprintf(stderr,
                "[ribbon selftest] (4e) warped contact: rc=%d "
                "contact_gap=%.2f (want ~1) span_B=%.1f (want ~98) "
                "ply_sep=%.0f (want > 300)\n",
                rc, contact_gap, span_b, ply_sep);
        st_check(rc == 0, "warped contact rc", &fails);
        st_check(rc == 0 && fabs(contact_gap - 1.0) < 2.5,
                 "contact terms pull a warped chart into agreement", &fails);
        st_check(rc == 0 && fabs(span_b - 98.0) < 5.0,
                 "warped chart stays metric after the pull", &fails);
        st_check(rc == 0 && ply_sep > 300.0,
                 "a normal-separated ply is NOT tied", &fails);
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
        int32_t *f2 = RIB_ALLOC_ARRAY(
            arena, int32_t, (nfc0 + (size_t)(nh-1)*2) * 3);
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
        int32_t vm[3] = {0,0,0};
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
        (void)transfer_uv(arena, v, 3, f, 1, vm, 1,
                          t, &S, 2.0, &o, NULL, &R);
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

    /* (11) transfer must use lifted winding as a sheet label.  The wrong-ply
     * sample is spatially closest but one full turn away; the slightly farther
     * same-ply sample is the only valid source of material U. */
    {
        Arena_Mark mark = Arena_save(arena);
        float v[9] = {1,0,0,  1,1,0,  1,0,1};
        int32_t f[3] = {0,1,2};
        int32_t vm[3] = {0,0,0};
        double t[3] = {1,1,1};
        Sample smp[2];
        memset(smp, 0, sizeof(smp));
        smp[0].c1 = 0.0; smp[0].u = 100.0; smp[0].phi = 2.0 * M_PI;
        smp[1].c1 = 0.5; smp[1].u =   5.0; smp[1].phi = 0.0;
        for (int i = 0; i < 2; i++) {
            smp[i].p[0] = 1.0;
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
        float prefilled[6] = {0,0, 0,0, 0,0};
        float ref_phi[3] = {0,0,0};
        RibbonResult R;
        memset(&R, 0, sizeof(R));
        R.uv = prefilled;
        RibbonOpts o;
        RibbonOpts_default(&o);
        o.emit_global = 1;
        o.reference_phi = ref_phi;
        (void)transfer_uv(arena, v, 3, f, 1, vm, 1,
                          t, &S, 2.0, &o, NULL, &R);
        fprintf(stderr,
            "[ribbon selftest] (11) phase-safe transfer: u=%.3f "
            "rejects=%zu blocked=%zu\n",
            (double)R.uv[0], R.uv_phase_rejects, R.uv_phase_blocked);
        st_check(fabs((double)R.uv[0] - 5.0) < 1e-6,
                 "same-sheet sample wins over closer adjacent ply", &fails);
        st_check(R.uv_phase_rejects >= 3 && R.uv_phase_blocked == 0,
                 "wrong-ply candidates are diagnosed", &fails);
        Arena_restore(arena, mark);
    }

    /* (12) source-chart provenance outranks Euclidean proximity.  Two samples
     * have identical lifted phase; only the farther one was generated by the
     * vertex's own connected input chart. */
    {
        Arena_Mark mark = Arena_save(arena);
        float v[9] = {1,0,0,  1,1,0,  1,0,1};
        int32_t f[3] = {0,1,2};
        int32_t vm[3] = {0,0,0};
        double t[3] = {1,1,1};
        Sample smp[2];
        memset(smp, 0, sizeof(smp));
        smp[0].c1 = 0.0; smp[0].u = 100.0; smp[0].phi = 0.0;
        smp[0].chain = 0; smp[0].slice = 0;
        smp[1].c1 = 0.5; smp[1].u =   5.0; smp[1].phi = 0.0;
        smp[1].chain = 1; smp[1].slice = 0;
        smp[0].p[0] = smp[1].p[0] = 1.0;
        Chain chn[2];
        memset(chn, 0, sizeof(chn));
        chn[0].first = 0; chn[0].count = 1; chn[0].mesh_comp = 1;
        chn[1].first = 1; chn[1].count = 1; chn[1].mesh_comp = 0;
        SliceSet S;
        memset(&S, 0, sizeof(S));
        S.smp = smp; S.n_smp = 2;
        S.chn = chn; S.n_chn = 2;
        S.nplanes = 1; S.tmin = 0; S.tmax = 2;
        float prefilled[6] = {0,0, 0,0, 0,0};
        float ref_phi[3] = {0,0,0};
        RibbonResult R;
        memset(&R, 0, sizeof(R));
        R.uv = prefilled;
        RibbonOpts o;
        RibbonOpts_default(&o);
        o.emit_global = 1;
        o.reference_phi = ref_phi;
        (void)transfer_uv(arena, v, 3, f, 1, vm, 2,
                          t, &S, 2.0, &o, NULL, &R);
        fprintf(stderr,
            "[ribbon selftest] (12) component-safe transfer: u=%.3f "
            "rejects=%zu blocked=%zu\n",
            (double)R.uv[0], R.uv_component_rejects,
            R.uv_component_blocked);
        st_check(fabs((double)R.uv[0] - 5.0) < 1e-6,
                 "own-chart sample wins over closer foreign chart", &fails);
        st_check(R.uv_component_rejects >= 3 &&
                 R.uv_component_blocked == 0,
                 "foreign-chart candidates are diagnosed", &fails);
        Arena_restore(arena, mark);
    }

    /* (13) two disconnected metric solves sampling the same connected source
     * chart have independent additive gauges.  Bracketing observations at the
     * identical source vertices must reconcile them without altering either
     * solve's internal arclength. */
    {
        Arena_Mark mark = Arena_save(arena);
        float v[9] = {2,0,0,  2,1,0,  2,0,1};
        int32_t f[3] = {0,1,2};
        int32_t vm[3] = {0,0,0};
        double t[3] = {2,2,2};
        Sample smp[2];
        memset(smp, 0, sizeof(smp));
        smp[0].p[0] = 1; smp[0].u = 5;   smp[0].phi = 0;
        smp[0].chain = 0; smp[0].slice = 0;
        smp[1].p[0] = 3; smp[1].u = 105; smp[1].phi = 0;
        smp[1].chain = 1; smp[1].slice = 1;
        Chain chn[2];
        memset(chn, 0, sizeof(chn));
        chn[0].first = 0; chn[0].count = 1;
        chn[0].mesh_comp = 0; chn[0].solve_comp = 0;
        chn[1].first = 1; chn[1].count = 1;
        chn[1].mesh_comp = 0; chn[1].solve_comp = 1;
        SliceSet S;
        memset(&S, 0, sizeof(S));
        S.smp = smp; S.n_smp = 2;
        S.chn = chn; S.n_chn = 2;
        S.nplanes = 2; S.tmin = 0; S.tmax = 4;
        float prefilled[6] = {0,0, 0,0, 0,0};
        float ref_phi[3] = {0,0,0};
        RibbonResult R;
        memset(&R, 0, sizeof(R));
        R.uv = prefilled;
        R.n_qp_comps = 2;
        RibbonOpts o;
        RibbonOpts_default(&o);
        o.emit_global = 1;
        o.reference_phi = ref_phi;
        (void)transfer_uv(arena, v, 3, f, 1, vm, 1,
                          t, &S, 2.0, &o, NULL, &R);
        fprintf(stderr,
            "[ribbon selftest] (13) transfer gauge: u=%.3f "
            "obs=%zu rel=%zu max_shift=%.3f\n",
            (double)R.uv[0], R.uv_gauge_observations,
            R.uv_gauge_relations, R.uv_gauge_max_shift);
        st_check(fabs((double)R.uv[0] - 5.0) < 1e-6,
                 "bracketing slices reconcile additive U gauges", &fails);
        st_check(R.uv_gauge_observations >= 3 &&
                 R.uv_gauge_relations == 1 &&
                 fabs(R.uv_gauge_max_shift - 100.0) < 1e-6,
                 "gauge reconciliation diagnostics", &fails);
        Arena_restore(arena, mark);
    }

    /* (14) within one chart and one solve, a self-contact sample may be
     * spatially closest while carrying a different lifted phase.  Universal-
     * cover distance must rank the exact-phase sample first. */
    {
        Arena_Mark mark = Arena_save(arena);
        float v[9] = {1,0,0,  1,1,0,  1,0,1};
        int32_t f[3] = {0,1,2};
        int32_t vm[3] = {0,0,0};
        double t[3] = {1,1,1};
        Sample smp[2];
        memset(smp, 0, sizeof(smp));
        smp[0].c1 = 100; smp[0].r = 100;
        smp[0].u = 100; smp[0].phi = 0.05;
        smp[1].c1 = 101; smp[1].r = 101;
        smp[1].u = 5; smp[1].phi = 0;
        for (int i = 0; i < 2; i++) {
            smp[i].p[0] = 1;
            smp[i].chain = 0;
            smp[i].slice = 0;
        }
        Chain chn;
        memset(&chn, 0, sizeof(chn));
        chn.count = 2; chn.mesh_comp = 0; chn.solve_comp = 0;
        SliceSet S;
        memset(&S, 0, sizeof(S));
        S.smp = smp; S.n_smp = 2;
        S.chn = &chn; S.n_chn = 1;
        S.nplanes = 1; S.tmin = 0; S.tmax = 2;
        float prefilled[6] = {100,0, 100,0, 100,0};
        float ref_phi[3] = {0,0,0};
        RibbonResult R;
        memset(&R, 0, sizeof(R));
        R.uv = prefilled;
        RibbonOpts o;
        RibbonOpts_default(&o);
        o.emit_global = 1;
        o.reference_phi = ref_phi;
        (void)transfer_uv(arena, v, 3, f, 1, vm, 1,
                          t, &S, 2.0, &o, NULL, &R);
        fprintf(stderr,
            "[ribbon selftest] (14) cover-metric transfer: u=%.3f\n",
            (double)R.uv[0]);
        st_check(fabs((double)R.uv[0] - 5.0) < 1e-6,
                 "lifted-phase metric defeats closer self-contact", &fails);
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
