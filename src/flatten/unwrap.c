#define _USE_MATH_DEFINES
#include "unwrap.h"
#include "winding_register.h"

#include "../common/csr.h"
#include "../common/pca.h"
#include "../common/pipeline_constants.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---- small helpers -------------------------------------------------------- */

/* (radius, index) pair so we can seed each component's winding integration at
 * its innermost (smallest-radius) vertex via one ascending sort. */
typedef struct { double r; int32_t idx; } RIdx;

static int cmp_ridx(const void *pa, const void *pb)
{
    const RIdx *a = (const RIdx *)pa;
    const RIdx *b = (const RIdx *)pb;
    if (a->r < b->r) return -1;
    if (a->r > b->r) return 1;
    if (a->idx < b->idx) return -1;
    if (a->idx > b->idx) return 1;
    return 0;
}

/* Wrap an angle into (-pi, pi]. */
static double wrap_to_pi(double a)
{
    double x = fmod(a + M_PI, 2.0 * M_PI);
    if (x < 0.0) x += 2.0 * M_PI;
    return x - M_PI;
}

/* ---- Cross-gauge-island winding sync (overlap-offset solver) -------------
 * The registration's relation forest leaves disconnected GAUGE islands whose
 * lifted phase differs by an unknown integer number of turns; downstream that
 * is exactly why continuation joins radius-reject across islands and why
 * phase placement extrapolates absurd u domains (10x rung: 3.35M vox against
 * a ~500k physical maximum).  Radial overlap is the gauge-free evidence (the
 * directed Villa overlap-offset design, with their snap-gate abstention
 * ported): two samples at the same (z, theta) differ in TRUE lifted turns by
 * their radial order,
 *     k_AB = (w_a - w_b) - alpha * (r_a - r_b) / pitch,   w = sense*Phi/2pi,
 * so k_AB clusters at the integer gauge offset between islands A and B (the
 * radial-order sign alpha is self-calibrated from same-island pairs).  A pair
 * offset is accepted only when >= UNWRAP_SYNC_MIN_OBS observations put >= 70%
 * within 0.25 turn of ONE integer (the snap gate); accepted offsets solve
 * over a spanning tree per connected component with EXACT integer agreement
 * demanded of every non-tree edge (the cycle certificate a scalar Huber
 * cannot give).  Any weaker evidence ABSTAINS and changes nothing. */
enum {
    UNWRAP_SYNC_MIN_OBS = 24,
    UNWRAP_SYNC_KSPAN = 8,           /* histogram half-width around anchor */
    UNWRAP_SYNC_BUCKET_CAP = 3,      /* samples kept per island per bucket */
    UNWRAP_SYNC_PAIR_TAB = 1 << 20   /* open-addressing pair table slots */
};

typedef struct {
    int32_t ga, gb;                  /* ga < gb; empty slot iff gb == 0 */
    int32_t anchor;                  /* rounded k of the first observation */
    int32_t n, n_snap;
    int32_t hist[2 * UNWRAP_SYNC_KSPAN + 1];
} UnwrapSyncPair;

typedef struct {
    int64_t key;                     /* (slab << 8) | theta_bin */
    int32_t island;
    int32_t idx;
} UnwrapSyncKey;

static int cmp_sync_key(const void *pa, const void *pb)
{
    const UnwrapSyncKey *a = (const UnwrapSyncKey *)pa;
    const UnwrapSyncKey *b = (const UnwrapSyncKey *)pb;
    if (a->key != b->key) return a->key < b->key ? -1 : 1;
    if (a->island != b->island) return a->island < b->island ? -1 : 1;
    return a->idx < b->idx ? -1 : (a->idx > b->idx ? 1 : 0);
}

static UnwrapSyncPair *unwrap_sync_pair_slot(UnwrapSyncPair *tab, size_t cap,
                                             int32_t ga, int32_t gb)
{
    size_t h = ((size_t)(uint32_t)ga * 2654435761u ^
                (size_t)(uint32_t)gb * 40503u) % cap;
    for (size_t probe = 0; probe < cap; probe++) {
        UnwrapSyncPair *p = &tab[(h + probe) % cap];
        if (p->gb == 0) {
            p->ga = ga;
            p->gb = gb;
            return p;
        }
        if (p->ga == ga && p->gb == gb) return p;
    }
    return NULL;
}

/* Returns the number of islands whose gauge was corrected (0 = abstained or
 * nothing to do).  Applies corrections to Phi and turn_correction in place. */
static size_t unwrap_sync_gauges(Arena_T arena, size_t nv,
                                 const double *t, const double *r,
                                 double *Phi,
                                 const int32_t *comp, int32_t ncomp,
                                 const int32_t *component_island,
                                 int32_t *turn_correction,
                                 int winding_sense, double pitch)
{
    enum { NTH = 64 };
    const double slab_h = 16.0;
    int32_t nisl = 0;
    size_t synced = 0;
    for (int32_t c = 0; c < ncomp; c++)
        if (component_island[c] >= nisl) nisl = component_island[c] + 1;
    if (nisl <= 1 || !(pitch > 0.0) || nv < 2) return 0;

    Arena_Mark mark = Arena_save(arena);
    UnwrapSyncKey *keys = (UnwrapSyncKey *)ARENA_ALLOC(
        arena, (long)((nv + 1) * sizeof(UnwrapSyncKey)));
    size_t *isl_pop = (size_t *)ARENA_CALLOC(
        arena, (size_t)nisl + 1, sizeof(size_t));
    double t_lo = 1e300;
    size_t nk = 0;
    for (size_t i = 0; i < nv; i++)
        if (t[i] < t_lo) t_lo = t[i];
    for (size_t i = 0; i < nv; i++) {
        int32_t g = comp[i] >= 0 && comp[i] < ncomp
                  ? component_island[comp[i]] : -1;
        double th = 0.0;
        if (g < 0 || g >= nisl) continue;
        if (!isfinite(Phi[i]) || !isfinite(r[i]) || !isfinite(t[i])) continue;
        th = wrap_to_pi(Phi[i]);  /* physical angle from the lift */
        keys[nk].key = ((int64_t)((t[i] - t_lo) / slab_h) << 8) |
                       (int64_t)(((th + M_PI) / (2.0 * M_PI)) * (NTH - 1));
        keys[nk].island = g;
        keys[nk].idx = (int32_t)i;
        isl_pop[g]++;
        nk++;
    }
    {
        /* Coverage tells the campaign whether a global phase sheet is even
         * on the table: the top gauge island's share of samples is the
         * fraction of the scroll the phase authority can honestly place. */
        size_t top1 = 0, top2 = 0, top3 = 0, tot = 0;
        for (int32_t g = 0; g < nisl; g++) {
            size_t p = isl_pop[g];
            tot += p;
            if (p > top1) {
                top3 = top2;
                top2 = top1;
                top1 = p;
            } else if (p > top2) {
                top3 = top2;
                top2 = p;
            } else if (p > top3) {
                top3 = p;
            }
        }
        if (tot > 0)
            fprintf(stderr,
                    "  gauge sync: coverage top3 = %.1f%% / %.1f%% / %.1f%% "
                    "of %zu samples across %d gauge islands\n",
                    100.0 * (double)top1 / (double)tot,
                    100.0 * (double)top2 / (double)tot,
                    100.0 * (double)top3 / (double)tot, tot, nisl);
    }
    if (nk < 2 * (size_t)UNWRAP_SYNC_MIN_OBS) {
        Arena_restore(arena, mark);
        return 0;
    }
    qsort(keys, nk, sizeof *keys, cmp_sync_key);

    /* keep <= BUCKET_CAP samples per (bucket, island); groups are
     * contiguous after the sort, so a run counter against the ORIGINAL
     * previous element suffices */
    size_t kept = 0;
    {
        size_t run = 0;
        for (size_t i = 0; i < nk; i++) {
            if (i > 0 && keys[i].key == keys[i - 1].key &&
                keys[i].island == keys[i - 1].island)
                run++;
            else
                run = 0;
            if (run < (size_t)UNWRAP_SYNC_BUCKET_CAP)
                keys[kept++] = keys[i];
        }
    }

    /* alpha: radial-order sign from same-island pairs */
    double sum_dwdr = 0.0;
    const double inv_two_pi = 1.0 / (2.0 * M_PI);
    for (size_t b0 = 0; b0 < kept;) {
        size_t b1 = b0;
        while (b1 < kept && keys[b1].key == keys[b0].key) b1++;
        for (size_t a = b0; a < b1; a++)
            for (size_t b = a + 1; b < b1; b++) {
                if (keys[a].island != keys[b].island) continue;
                {
                    size_t ia = (size_t)keys[a].idx, ib = (size_t)keys[b].idx;
                    double dw = (double)winding_sense *
                                (Phi[ia] - Phi[ib]) * inv_two_pi;
                    double dr = r[ia] - r[ib];
                    if (fabs(dr) > 0.5 * pitch) sum_dwdr += dw * dr;
                }
            }
        b0 = b1;
    }
    {
        double alpha = sum_dwdr >= 0.0 ? 1.0 : -1.0;
        UnwrapSyncPair *tab = (UnwrapSyncPair *)ARENA_CALLOC(
            arena, (size_t)UNWRAP_SYNC_PAIR_TAB, sizeof(UnwrapSyncPair));
        size_t nobs = 0;
        for (size_t b0 = 0; b0 < kept;) {
            size_t b1 = b0, npairs = 0;
            while (b1 < kept && keys[b1].key == keys[b0].key) b1++;
            for (size_t a = b0; a < b1 && npairs < 512; a++)
                for (size_t b = a + 1; b < b1 && npairs < 512; b++) {
                    int32_t ga = keys[a].island, gb = keys[b].island;
                    size_t ia = (size_t)keys[a].idx;
                    size_t ib = (size_t)keys[b].idx;
                    double dw = 0.0, kf = 0.0;
                    int32_t kr = 0;
                    UnwrapSyncPair *p = NULL;
                    if (ga == gb) continue;
                    if (ga > gb) {
                        int32_t tg = ga;
                        size_t ti = ia;
                        ga = gb; gb = tg;
                        ia = ib; ib = ti;
                    }
                    dw = (double)winding_sense * (Phi[ia] - Phi[ib]) *
                         inv_two_pi;
                    kf = dw - alpha * (r[ia] - r[ib]) / pitch;
                    if (!isfinite(kf) || fabs(kf) > 500.0) continue;
                    kr = (int32_t)lround(kf);
                    p = unwrap_sync_pair_slot(
                        tab, (size_t)UNWRAP_SYNC_PAIR_TAB, ga, gb);
                    if (p == NULL) continue;
                    if (p->n == 0) p->anchor = kr;
                    p->n++;
                    npairs++;
                    nobs++;
                    if (fabs(kf - (double)kr) <= 0.25 &&
                        kr - p->anchor >= -UNWRAP_SYNC_KSPAN &&
                        kr - p->anchor <= UNWRAP_SYNC_KSPAN) {
                        p->hist[kr - p->anchor + UNWRAP_SYNC_KSPAN]++;
                        p->n_snap++;
                    }
                }
            b0 = b1;
        }

        /* accept pairs under the snap gate; solve per component with the
         * exact cycle certificate */
        {
            int32_t *edge_a = (int32_t *)ARENA_ALLOC(
                arena, (long)(((size_t)nisl * 4 + 8) * sizeof(int32_t)));
            int32_t *edge_b = (int32_t *)ARENA_ALLOC(
                arena, (long)(((size_t)nisl * 4 + 8) * sizeof(int32_t)));
            int32_t *edge_k = (int32_t *)ARENA_ALLOC(
                arena, (long)(((size_t)nisl * 4 + 8) * sizeof(int32_t)));
            size_t nedge = 0, printed = 0, naccept = 0;
            for (size_t s = 0; s < (size_t)UNWRAP_SYNC_PAIR_TAB; s++) {
                const UnwrapSyncPair *p = &tab[s];
                int best = 0, bi = -1;
                if (p->gb == 0 || p->n < UNWRAP_SYNC_MIN_OBS) continue;
                for (int h = 0; h < 2 * UNWRAP_SYNC_KSPAN + 1; h++)
                    if (p->hist[h] > best) {
                        best = p->hist[h];
                        bi = h;
                    }
                if (best < UNWRAP_SYNC_MIN_OBS ||
                    (double)best < 0.70 * (double)p->n)
                    continue;
                if (nedge < (size_t)nisl * 4 + 8) {
                    edge_a[nedge] = p->ga;
                    edge_b[nedge] = p->gb;
                    edge_k[nedge] = p->anchor + bi - UNWRAP_SYNC_KSPAN;
                    nedge++;
                }
                naccept++;
                if (printed < 20) {
                    fprintf(stderr,
                            "  gauge sync: island %d ~ island %d offset %d "
                            "(obs=%d snap=%d)\n",
                            p->ga, p->gb,
                            p->anchor + bi - UNWRAP_SYNC_KSPAN, p->n, best);
                    printed++;
                }
            }
            if (nedge > 0) {
                int32_t *off = (int32_t *)ARENA_ALLOC(
                    arena, (long)(((size_t)nisl + 1) * sizeof(int32_t)));
                int32_t *root = (int32_t *)ARENA_ALLOC(
                    arena, (long)(((size_t)nisl + 1) * sizeof(int32_t)));
                uint8_t *have = (uint8_t *)ARENA_CALLOC(
                    arena, (size_t)nisl + 1, 1);
                uint8_t *bad_root = (uint8_t *)ARENA_CALLOC(
                    arena, (size_t)nisl + 1, 1);
                int32_t *queue = (int32_t *)ARENA_ALLOC(
                    arena, (long)(((size_t)nisl + 1) * sizeof(int32_t)));
                /* BFS over accepted edges; the root of each component keeps
                 * offset 0 (sync is relative, so the root choice is
                 * arbitrary) */
                uint8_t *touched = (uint8_t *)ARENA_CALLOC(
                    arena, (size_t)nisl + 1, 1);
                for (size_t e = 0; e < nedge; e++) {
                    touched[edge_a[e]] = 1;
                    touched[edge_b[e]] = 1;
                }
                for (int32_t g = 0; g < nisl; g++) root[g] = -1;
                for (int32_t seed = 0; seed < nisl; seed++) {
                    size_t qh = 0, qt = 0;
                    if (!touched[seed] || root[seed] >= 0) continue;
                    root[seed] = seed;
                    off[seed] = 0;
                    have[seed] = 1;
                    queue[qt++] = seed;
                    while (qh < qt) {
                        int32_t g = queue[qh++];
                        for (size_t e = 0; e < nedge; e++) {
                            int32_t o = -1, k = 0;
                            if (edge_a[e] == g) {
                                o = edge_b[e];
                                /* edge k = off[ga] - off[gb] */
                                k = off[g] - edge_k[e];
                            } else if (edge_b[e] == g) {
                                o = edge_a[e];
                                k = off[g] + edge_k[e];
                            }
                            if (o < 0) continue;
                            if (root[o] < 0) {
                                root[o] = root[g];
                                off[o] = k;
                                have[o] = 1;
                                queue[qt++] = o;
                            } else if (off[o] != k) {
                                bad_root[root[g]] = 1; /* cycle violation */
                            }
                        }
                    }
                }
                /* apply to every island in a clean component */
                {
                    size_t bad_islands = 0;
                    for (int32_t g = 0; g < nisl; g++) {
                        if (!have[g] || root[g] < 0) continue;
                        if (bad_root[root[g]]) {
                            bad_islands++;
                            continue;
                        }
                        if (off[g] != 0) synced++;
                    }
                    if (synced > 0) {
                        for (size_t i = 0; i < nv; i++) {
                            int32_t g = comp[i] >= 0 && comp[i] < ncomp
                                      ? component_island[comp[i]] : -1;
                            if (g < 0 || !have[g] || root[g] < 0 ||
                                bad_root[root[g]] || off[g] == 0)
                                continue;
                            Phi[i] -= (double)winding_sense * 2.0 * M_PI *
                                      (double)off[g];
                        }
                        for (int32_t c = 0; c < ncomp; c++) {
                            int32_t g = component_island[c];
                            if (g < 0 || g >= nisl || !have[g] ||
                                root[g] < 0 || bad_root[root[g]] ||
                                off[g] == 0)
                                continue;
                            turn_correction[c] -= off[g];
                        }
                    }
                    fprintf(stderr,
                            "  gauge sync: %d islands, %zu obs, %zu pair(s) "
                            "accepted -> %zu island(s) re-gauged, %zu "
                            "abstained (cycle), alpha=%+.0f\n",
                            nisl, nobs, naccept, synced, bad_islands, alpha);
                }
            } else {
                fprintf(stderr,
                        "  gauge sync: %d islands, %zu obs, no pair passed "
                        "the snap gate (abstain)\n",
                        nisl, nobs);
            }
        }
    }
    Arena_restore(arena, mark);
    return synced;
}

/* ---- main ----------------------------------------------------------------- */

int Unwrap_run(Arena_T arena,
               const float *verts, size_t nv,
               const int32_t *faces, size_t nf,
               const UnwrapOpts *opts, UnwrapResult *out)
{
    assert(arena);
    assert(out);

    memset(out, 0, sizeof(*out));
    out->spiral_r2 = -1.0;

    if (nv < 3 || nf < 1 || verts == NULL || faces == NULL) {
        return -1;  /* degenerate -- caller falls back to the raw mesh */
    }

    UnwrapAxisMode mode = opts ? opts->axis_mode : UNWRAP_AXIS_AUTO;
    double wrap_spacing = opts ? opts->wrap_spacing : 0.0;

    /* Persistent output (survives the scratch restore below). */
    out->uv = (float *)ARENA_ALLOC(arena, (long)(nv * 2 * sizeof(float)));
    int keep_phi = (opts != NULL && opts->keep_phi != 0);
    if (keep_phi) {
        out->phi = (float *)ARENA_ALLOC(arena, (long)(nv * sizeof(float)));
        out->island = (int32_t *)ARENA_ALLOC(
            arena, (long)(nv * sizeof(int32_t)));
        out->continuation_island = (int32_t *)ARENA_ALLOC(
            arena, (long)(nv * sizeof(int32_t)));
        out->mesh_component = (int32_t *)ARENA_ALLOC(
            arena, (long)(nv * sizeof(int32_t)));
    }

    Arena_Mark mark = Arena_save(arena);

    /* --- 1. Axis + centroid (in (z,y,x) component order). --- */
    float axis[3] = {0, 0, 0};
    float centroid[3] = {0, 0, 0};
    if (mode == UNWRAP_AXIS_EXTERNAL && opts) {
        /* caller-supplied umbilicus line (z,y,x); normalize the direction. */
        double n = sqrt((double)opts->axis_dir[0] * opts->axis_dir[0]
                      +  (double)opts->axis_dir[1] * opts->axis_dir[1]
                      +  (double)opts->axis_dir[2] * opts->axis_dir[2]);
        if (n < 1e-12) n = 1.0;
        axis[0] = (float)(opts->axis_dir[0] / n);
        axis[1] = (float)(opts->axis_dir[1] / n);
        axis[2] = (float)(opts->axis_dir[2] / n);
        centroid[0] = opts->axis_point[0];
        centroid[1] = opts->axis_point[1];
        centroid[2] = opts->axis_point[2];
    } else {
        PCA_principal_axis(verts, nv, axis, centroid);
        if (mode == UNWRAP_AXIS_Z) { axis[0] = 1; axis[1] = 0; axis[2] = 0; }
        else if (mode == UNWRAP_AXIS_Y) { axis[0] = 0; axis[1] = 1; axis[2] = 0; }
        else if (mode == UNWRAP_AXIS_X) { axis[0] = 0; axis[1] = 0; axis[2] = 1; }
    }

    float e1f[3], e2f[3];
    PCA_orthonormal_basis(axis, e1f, e2f);
    double ax[3] = { (double)axis[0], (double)axis[1], (double)axis[2] };
    double e1[3] = { (double)e1f[0], (double)e1f[1], (double)e1f[2] };
    double e2[3] = { (double)e2f[0], (double)e2f[1], (double)e2f[2] };
    double cen[3] = { (double)centroid[0], (double)centroid[1],
                      (double)centroid[2] };

    /* --- 2. Per-vertex axial coord t and in-plane coords (c1, c2). --- */
    double *t  = (double *)ARENA_ALLOC(arena, (long)(nv * sizeof(double)));
    double *c1 = (double *)ARENA_ALLOC(arena, (long)(nv * sizeof(double)));
    double *c2 = (double *)ARENA_ALLOC(arena, (long)(nv * sizeof(double)));
    double tmin = INFINITY, tmax = -INFINITY;
    for (size_t i = 0; i < nv; i++) {
        double d0 = (double)verts[i * 3 + 0] - cen[0];
        double d1 = (double)verts[i * 3 + 1] - cen[1];
        double d2 = (double)verts[i * 3 + 2] - cen[2];
        double ti = d0 * ax[0] + d1 * ax[1] + d2 * ax[2];
        t[i]  = ti;
        c1[i] = d0 * e1[0] + d1 * e1[1] + d2 * e1[2];
        c2[i] = d0 * e2[0] + d1 * e2[1] + d2 * e2[2];
        if (ti < tmin) tmin = ti;
        if (ti > tmax) tmax = ti;
    }
    double span = tmax - tmin;
    if (span < 1e-9) span = 1e-9;

    /* --- 3. Centerline: mean in-plane position per axial slice. --- */
    const int K = FLATTEN_AXIAL_BINS;
    double *cb1 = (double *)ARENA_CALLOC(arena, (long)K, (long)sizeof(double));
    double *cb2 = (double *)ARENA_CALLOC(arena, (long)K, (long)sizeof(double));
    double *cnt = (double *)ARENA_CALLOC(arena, (long)K, (long)sizeof(double));
    char   *fil = (char   *)ARENA_CALLOC(arena, (long)K, (long)sizeof(char));
    for (size_t i = 0; i < nv; i++) {
        int b = (int)((t[i] - tmin) / span * (double)(K - 1));
        if (b < 0) b = 0;
        if (b > K - 1) b = K - 1;
        cb1[b] += c1[i];
        cb2[b] += c2[i];
        cnt[b] += 1.0;
    }
    int f0 = -1;
    for (int b = 0; b < K; b++) {
        if (cnt[b] > 0.0) {
            cb1[b] /= cnt[b];
            cb2[b] /= cnt[b];
            fil[b] = 1;
            if (f0 < 0) f0 = b;
        }
    }
    /* Fill empty bins by carrying the nearest filled center (fwd then back). */
    for (int b = f0 + 1; b < K; b++) {
        if (!fil[b]) { cb1[b] = cb1[b - 1]; cb2[b] = cb2[b - 1]; fil[b] = 1; }
    }
    for (int b = f0 - 1; b >= 0; b--) {
        cb1[b] = cb1[b + 1]; cb2[b] = cb2[b + 1]; fil[b] = 1;
    }
    /* Smooth the centerline (moving average, half-window H). */
    const int H = FLATTEN_CENTERLINE_SMOOTH;
    double *s1 = (double *)ARENA_ALLOC(arena, (long)((size_t)K * sizeof(double)));
    double *s2 = (double *)ARENA_ALLOC(arena, (long)((size_t)K * sizeof(double)));
    for (int b = 0; b < K; b++) {
        double a1 = 0.0, a2 = 0.0;
        int n = 0;
        for (int w = b - H; w <= b + H; w++) {
            if (w < 0 || w >= K) continue;
            a1 += cb1[w]; a2 += cb2[w]; n++;
        }
        s1[b] = a1 / (double)n;
        s2[b] = a2 / (double)n;
    }
    /* External axis: pin the centerline to the supplied umbilicus. c1,c2 are
     * already relative to axis_point, so the per-slice center offset is 0 -- this
     * makes r/theta measure from the TRUE center, not the block's local mean
     * (which sits ~1500 vox off the umbilicus for an off-center sub-block). */
    if (mode == UNWRAP_AXIS_EXTERNAL) {
        for (int b = 0; b < K; b++) { s1[b] = 0.0; s2[b] = 0.0; }
    }

    /* --- 4. Cylindrical depth (r) and angle (theta) per vertex. --- */
    double *r     = (double *)ARENA_ALLOC(arena, (long)(nv * sizeof(double)));
    double *theta = (double *)ARENA_ALLOC(arena, (long)(nv * sizeof(double)));
    for (size_t i = 0; i < nv; i++) {
        double fb = (t[i] - tmin) / span * (double)(K - 1);
        int b0 = (int)floor(fb);
        if (b0 < 0) b0 = 0;
        if (b0 > K - 1) b0 = K - 1;
        int b1 = (b0 + 1 < K) ? b0 + 1 : K - 1;
        double w = fb - (double)b0;
        if (w < 0.0) w = 0.0;
        if (w > 1.0) w = 1.0;
        double ce1 = s1[b0] * (1.0 - w) + s1[b1] * w;
        double ce2 = s2[b0] * (1.0 - w) + s2[b1] * w;
        double p1 = c1[i] - ce1;
        double p2 = c2[i] - ce2;
        r[i]     = sqrt(p1 * p1 + p2 * p2);
        theta[i] = atan2(p2, p1);
    }

    /* --- 5. Graph-integrated winding Phi, per connected component. --- */
    CSR_T adj = CSR_from_faces(arena, faces, nf, nv);
    const int32_t *off = CSR_offset(adj);
    const int32_t *tgt = CSR_target(adj);

    double  *Phi   = (double  *)ARENA_ALLOC(arena, (long)(nv * sizeof(double)));
    int32_t *comp  = (int32_t *)ARENA_ALLOC(arena, (long)(nv * sizeof(int32_t)));
    int32_t *queue = (int32_t *)ARENA_ALLOC(arena, (long)(nv * sizeof(int32_t)));
    for (size_t i = 0; i < nv; i++) comp[i] = -1;

    RIdx *ri = (RIdx *)ARENA_ALLOC(arena, (long)(nv * sizeof(RIdx)));
    for (size_t i = 0; i < nv; i++) { ri[i].r = r[i]; ri[i].idx = (int32_t)i; }
    qsort(ri, nv, sizeof(RIdx), cmp_ridx);

    int32_t ncomp = 0;
    for (size_t k = 0; k < nv; k++) {
        int32_t s = ri[k].idx;
        if (comp[s] != -1) continue;
        comp[s] = ncomp;
        Phi[s]  = theta[s];
        size_t head = 0, tail = 0;
        queue[tail++] = s;
        while (head < tail) {
            int32_t u = queue[head++];
            for (int32_t e = off[u]; e < off[u + 1]; e++) {
                int32_t v = tgt[e];
                if (comp[v] == -1) {
                    comp[v] = ncomp;
                    Phi[v]  = Phi[u] + wrap_to_pi(theta[v] - theta[u]);
                    queue[tail++] = v;
                }
            }
        }
        ncomp++;
    }
    out->n_components = (int)ncomp;

    /* --- 6. Diagnostic local spiral fit on the largest component. ---
     * This estimates pitch sense/scale only.  It is not an absolute winding
     * model for a crushed or eccentric scroll. */
    int32_t *csize = (int32_t *)ARENA_CALLOC(arena, (long)ncomp,
                                             (long)sizeof(int32_t));
    for (size_t i = 0; i < nv; i++) csize[comp[i]]++;
    int32_t big = 0;
    for (int32_t c = 1; c < ncomp; c++) if (csize[c] > csize[big]) big = c;

    double sa = 0.0, sb = 0.0, sr2 = -1.0;
    {
        double n = 0.0, Sx = 0.0, Sy = 0.0, Sxx = 0.0, Sxy = 0.0;
        for (size_t i = 0; i < nv; i++) {
            if (comp[i] != big) continue;
            double x = Phi[i] / (2.0 * M_PI);
            double y = r[i];
            n += 1.0; Sx += x; Sy += y; Sxx += x * x; Sxy += x * y;
        }
        double denom = n * Sxx - Sx * Sx;
        if (n >= (double)FLATTEN_SPIRAL_MIN_PTS && fabs(denom) > 1e-12) {
            sb = (n * Sxy - Sx * Sy) / denom;
            sa = (Sy - sb * Sx) / n;
            double ybar = Sy / n, sstot = 0.0, ssres = 0.0;
            for (size_t i = 0; i < nv; i++) {
                if (comp[i] != big) continue;
                double x = Phi[i] / (2.0 * M_PI);
                double pred = sa + sb * x;
                double dy = r[i] - ybar;
                double dr = r[i] - pred;
                sstot += dy * dy; ssres += dr * dr;
            }
            sr2 = (sstot > 1e-12) ? 1.0 - ssres / sstot : 1.0;
        }
    }
    int winding_sense = sb < 0.0 ? -1 : 1;
    /* A supplied pitch fixes the signed diagnostic slope but does not make the
     * spiral model true.  Measure its actual residual. */
    if (wrap_spacing > 0.0) {
        sb = (double)winding_sense * wrap_spacing;
        double nfit = 0.0, sacc = 0.0;
        for (size_t i = 0; i < nv; i++) {
            if (comp[i] != big) continue;
            sacc += r[i] - sb * (Phi[i] / (2.0 * M_PI));
            nfit += 1.0;
        }
        sa = (nfit > 0.0) ? sacc / nfit : 0.0;
        double ybar = 0.0;
        for (size_t i = 0; i < nv; i++)
            if (comp[i] == big) ybar += r[i];
        if (nfit > 0.0) ybar /= nfit;
        double sstot = 0.0, ssres = 0.0;
        for (size_t i = 0; i < nv; i++) {
            if (comp[i] != big) continue;
            double pred = sa + sb * (Phi[i] / (2.0 * M_PI));
            double dy = r[i] - ybar;
            double dr = r[i] - pred;
            sstot += dy * dy;
            ssres += dr * dr;
        }
        sr2 = sstot > 1e-12 ? 1.0 - ssres / sstot : 1.0;
    }
    out->spiral_a = sa;
    out->spiral_b = sb;
    out->spiral_r2 = sr2;

    /* --- 7. Synchronize the integer gauge of disconnected local lifts. ---
     * q is known inside each disk up to one integer.  Continuation and ray
     * order determine those integers; radius never predicts absolute q. */
    double *qturn = (double *)ARENA_ALLOC(
        arena, (long)(nv * sizeof(double)));
    for (size_t i = 0; i < nv; i++)
        qturn[i] = (double)winding_sense * Phi[i] / (2.0 * M_PI);
    int32_t *turn_correction = NULL, *component_island = NULL;
    int32_t *component_continuation_island = NULL;
    WindingRegisterStats wstats;
    double registration_pitch = wrap_spacing > 0.0
                              ? wrap_spacing : fabs(sb);
    if (WindingRegister_run(
            arena, verts, nv, t, r, theta, qturn, comp, ncomp, csize, big,
            tmin, registration_pitch, winding_sense,
            &turn_correction, &component_island,
            &component_continuation_island, &wstats) != 0) {
        fprintf(stderr, "unwrap: winding gauge registration failed\n");
        Arena_restore(arena, mark);
        return -1;
    }
    for (size_t i = 0; i < nv; i++)
        Phi[i] += (double)winding_sense * 2.0 * M_PI *
                  (double)turn_correction[comp[i]];

    /* --- 7a. Overlap-offset gauge sync across disconnected gauge islands.
     * Runs before the continuation join so joins compare SYNCED lifts. */
    (void)unwrap_sync_gauges(arena, nv, t, r, Phi, comp, ncomp,
                             component_island, turn_correction,
                             winding_sense, registration_pitch);

    /* --- 7b. Join split material islands by lifted-phase continuation. ---
     * The registration's local continuation gate (3-D proximity + helix rate)
     * can never bridge a whole missing arc, so one physical sheet ends up as
     * a CHAIN of material islands tiling consecutive phi intervals with real
     * holes between them (measured 4x: gaps 0.36-1.19 rad at 175-525 vox).
     * The global interval structure IS the continuation evidence: consecutive
     * disjoint phi intervals, overlapping axial support, and a junction
     * radius that continues the same wrap.  Ambiguous gaps (>= a turn's
     * fraction) are left split for an evidence-backed pass. */
    {
        enum { UNWRAP_JOIN_MIN_VERTS = 24 };
        const double join_dphi_max = 1.3;          /* rad; measured 4x chain */
        const double join_dr_tol = 0.5;            /* x pitch at the junction */
        const double join_axial_overlap = 0.3;     /* fraction of min span */
        int32_t nisl = 0;
        for (size_t c = 0; c < (size_t)ncomp; c++)
            if (component_continuation_island[c] >= nisl)
                nisl = component_continuation_island[c] + 1;
        if (nisl > 1) {
            double *ilo = (double *)ARENA_ALLOC(
                arena, (long)((size_t)nisl * sizeof(double)));
            double *ihi = (double *)ARENA_ALLOC(
                arena, (long)((size_t)nisl * sizeof(double)));
            double *tlo = (double *)ARENA_ALLOC(
                arena, (long)((size_t)nisl * sizeof(double)));
            double *thi = (double *)ARENA_ALLOC(
                arena, (long)((size_t)nisl * sizeof(double)));
            double *rhi_sum = (double *)ARENA_ALLOC(
                arena, (long)((size_t)nisl * sizeof(double)));
            double *rlo_sum = (double *)ARENA_ALLOC(
                arena, (long)((size_t)nisl * sizeof(double)));
            size_t *rhi_n = (size_t *)ARENA_ALLOC(
                arena, (long)((size_t)nisl * sizeof(size_t)));
            size_t *rlo_n = (size_t *)ARENA_ALLOC(
                arena, (long)((size_t)nisl * sizeof(size_t)));
            size_t *icount = (size_t *)ARENA_ALLOC(
                arena, (long)((size_t)nisl * sizeof(size_t)));
            int32_t *root = (int32_t *)ARENA_ALLOC(
                arena, (long)((size_t)nisl * sizeof(int32_t)));
            int32_t *order = (int32_t *)ARENA_ALLOC(
                arena, (long)((size_t)nisl * sizeof(int32_t)));
            for (int32_t k = 0; k < nisl; k++) {
                ilo[k] = INFINITY; ihi[k] = -INFINITY;
                tlo[k] = INFINITY; thi[k] = -INFINITY;
                rhi_sum[k] = rlo_sum[k] = 0.0;
                rhi_n[k] = rlo_n[k] = 0;
                icount[k] = 0;
                root[k] = k;
                order[k] = k;
            }
            for (size_t i = 0; i < nv; i++) {
                int32_t k = component_continuation_island[comp[i]];
                double w = (double)winding_sense * Phi[i];
                if (k < 0 || k >= nisl) continue;
                if (w < ilo[k]) ilo[k] = w;
                if (w > ihi[k]) ihi[k] = w;
                if (t[i] < tlo[k]) tlo[k] = t[i];
                if (t[i] > thi[k]) thi[k] = t[i];
                icount[k]++;
            }
            for (size_t i = 0; i < nv; i++) {
                int32_t k = component_continuation_island[comp[i]];
                double w = (double)winding_sense * Phi[i];
                if (k < 0 || k >= nisl) continue;
                if (w > ihi[k] - 0.3) { rhi_sum[k] += r[i]; rhi_n[k]++; }
                if (w < ilo[k] + 0.3) { rlo_sum[k] += r[i]; rlo_n[k]++; }
            }
            /* sort islands by phi_lo (insertion sort: nisl is small) */
            for (int32_t a = 1; a < nisl; a++) {
                int32_t key = order[a];
                int32_t b = a - 1;
                while (b >= 0 && ilo[order[b]] > ilo[key]) {
                    order[b + 1] = order[b];
                    b--;
                }
                order[b + 1] = key;
            }
            {
                size_t joins = 0;
                size_t join_report = 0;
                for (int32_t s = 0; s + 1 < nisl; s++) {
                    int32_t a = order[s], b = order[s + 1];
                    double gap, dr_expect, dr_meas, span_a, span_b, ovl;
                    if (icount[a] < UNWRAP_JOIN_MIN_VERTS ||
                        icount[b] < UNWRAP_JOIN_MIN_VERTS)
                        continue;
                    if (rhi_n[a] == 0 || rlo_n[b] == 0) continue;
                    gap = ilo[b] - ihi[a];
                    span_a = thi[a] - tlo[a];
                    span_b = thi[b] - tlo[b];
                    ovl = (thi[a] < thi[b] ? thi[a] : thi[b]) -
                          (tlo[a] > tlo[b] ? tlo[a] : tlo[b]);
                    dr_expect = registration_pitch *
                                (gap > 0.0 ? gap : 0.0) / (2.0 * M_PI);
                    dr_meas = rlo_sum[b] / (double)rlo_n[b] -
                              rhi_sum[a] / (double)rhi_n[a];
                    if (join_report < 12) {
                        fprintf(stderr,
                                "    join candidate %d(n=%zu phi<=%.2f)+"
                                "%d(n=%zu phi>=%.2f): gap=%.2f ovl=%.0f/"
                                "%.0f dr=%.2f (expect %.2f) -> %s\n",
                                a, icount[a], ihi[a], b, icount[b], ilo[b],
                                gap, ovl,
                                span_a < span_b ? span_a : span_b, dr_meas,
                                dr_expect,
                                !(gap > -0.3) || gap > join_dphi_max
                                    ? "gap-reject"
                                : !(ovl > join_axial_overlap *
                                          (span_a < span_b ? span_a : span_b))
                                    ? "axial-reject"
                                : fabs(dr_meas - dr_expect) >
                                      join_dr_tol * registration_pitch
                                    ? "radius-reject" : "JOIN");
                        join_report++;
                    }
                    /* a small junction overlap (<= 0.3 rad) is a genuine
                     * continuation whose ends double-cover one bin; a
                     * different wrap at the same phi fails the radius gate */
                    if (!(gap > -0.3) || gap > join_dphi_max) continue;
                    if (!(ovl > join_axial_overlap *
                                (span_a < span_b ? span_a : span_b)))
                        continue;
                    if (fabs(dr_meas - dr_expect) >
                        join_dr_tol * registration_pitch)
                        continue;
                    /* union (roots by path compression walk) */
                    {
                        int32_t ra = a, rb = b;
                        while (root[ra] != ra) ra = root[ra];
                        while (root[rb] != rb) rb = root[rb];
                        if (ra != rb) {
                            root[ra > rb ? ra : rb] = ra > rb ? rb : ra;
                            joins++;
                            fprintf(stderr,
                                    "  phase-continuation join: island %d "
                                    "(phi<=%.2f, r~%.1f) + island %d "
                                    "(phi>=%.2f, r~%.1f): gap=%.2f rad, "
                                    "dr=%.2f (expect %.2f)\n",
                                    a, ihi[a],
                                    rhi_sum[a] / (double)rhi_n[a], b, ilo[b],
                                    rlo_sum[b] / (double)rlo_n[b], gap,
                                    dr_meas, dr_expect);
                        }
                    }
                }
                if (joins > 0) {
                    /* compact relabel so downstream label gates stay in
                     * [0, count) */
                    int32_t *newlab = (int32_t *)ARENA_ALLOC(
                        arena, (long)((size_t)nisl * sizeof(int32_t)));
                    int32_t nnew = 0;
                    for (int32_t k = 0; k < nisl; k++) newlab[k] = -1;
                    for (int32_t k = 0; k < nisl; k++) {
                        int32_t rk = k;
                        while (root[rk] != rk) rk = root[rk];
                        if (newlab[rk] < 0) newlab[rk] = nnew++;
                        newlab[k] = newlab[rk];
                    }
                    for (size_t c = 0; c < (size_t)ncomp; c++) {
                        int32_t k = component_continuation_island[c];
                        if (k >= 0 && k < nisl)
                            component_continuation_island[c] = newlab[k];
                    }
                    fprintf(stderr,
                            "  phase-continuation join: %d -> %d material "
                            "islands (%zu join(s), gate dphi<=%.2f rad)\n",
                            nisl, nnew, joins, join_dphi_max);
                } else {
                    fprintf(stderr,
                            "  phase-continuation join: no eligible "
                            "consecutive-interval pairs (%d islands)\n",
                            nisl);
                }
            }
        }
    }
    out->winding_sense = winding_sense;
    out->winding_bins = wstats.bins;
    out->winding_strands = wstats.strands;
    out->continuation_observations = wstats.continuation_observations;
    out->order_observations = wstats.order_observations;
    out->order_observations_suppressed =
        wstats.order_observations_suppressed;
    out->winding_relations = wstats.relations;
    out->winding_eligible_relations = wstats.eligible_relations;
    out->winding_forest_relations = wstats.forest_relations;
    out->winding_order_relations_suppressed =
        wstats.order_relations_suppressed;
    out->winding_continuation_components = wstats.continuation_components;
    out->winding_relation_components = wstats.relation_components;
    out->winding_packed_relation_components =
        wstats.packed_relation_components;
    out->winding_packed_mesh_components = wstats.packed_mesh_components;
    out->winding_relation_conflicts = wstats.relation_conflicts;
    out->winding_observations_dropped = wstats.observations_dropped;
    out->continuation_satisfaction = wstats.continuation_satisfaction;
    out->order_satisfaction = wstats.order_satisfaction;
    out->turn_correction_min = wstats.correction_min;
    out->turn_correction_max = wstats.correction_max;

    /* --- 8. Assemble UV (length-like, each axis shifted to start at 0). ---
     * u is the winding converted to an arc length via a reference radius so the
     * grid is ~isotropic and matches vc_obj2tifxyz's metric (length) UV mode.
     * r_ref = median radius (ri is sorted ascending by r from step 5). */
    double umin = INFINITY, umax = -INFINITY;
    for (size_t i = 0; i < nv; i++) {
        if (Phi[i] < umin) umin = Phi[i];
        if (Phi[i] > umax) umax = Phi[i];
    }
    double r_ref = ri[nv / 2].r;
    if (r_ref < 1e-6) r_ref = 1.0;   /* flat/degenerate -> 1 vox per radian */
    for (size_t i = 0; i < nv; i++) {
        out->uv[i * 2 + 0] = (float)((Phi[i] - umin) * r_ref);
        out->uv[i * 2 + 1] = (float)(t[i] - tmin);
    }
    if (keep_phi) {
        for (size_t i = 0; i < nv; i++) {
            out->phi[i] = (float)Phi[i];
            out->island[i] = component_island[comp[i]];
            out->continuation_island[i] =
                component_continuation_island[comp[i]];
            out->mesh_component[i] = comp[i];
        }
    }

    out->axis[0] = axis[0]; out->axis[1] = axis[1]; out->axis[2] = axis[2];
    out->centroid[0] = centroid[0];
    out->centroid[1] = centroid[1];
    out->centroid[2] = centroid[2];
    out->r_ref  = r_ref;
    out->u_span = (umax - umin) * r_ref;
    out->v_span = tmax - tmin;
    out->turns  = (umax - umin) / (2.0 * M_PI);

    Arena_restore(arena, mark);
    return 0;
}

/* ============================================================================
 * Self-test: synthetic cylinder, spiral, and degenerate inputs.
 * ==========================================================================*/

/* OPEN strip wound ~360 deg around Z (a topological rectangle = one scroll
 * wrap; simply connected, NOT a closed tube). nu angular x nh height samples;
 * faces do NOT wrap the seam. Vertices in (z, y, x). Center offset from origin. */
static void build_cylinder(Arena_T arena, int nu, int nh, double R, double Hgt,
                           float **out_v, size_t *out_nv,
                           int32_t **out_f, size_t *out_nf)
{
    size_t nvv = (size_t)nu * (size_t)nh;
    size_t nff = (size_t)(nu - 1) * (size_t)(nh - 1) * 2;
    float   *v = (float *)ARENA_ALLOC(arena, (long)(nvv * 3 * sizeof(float)));
    int32_t *f = (int32_t *)ARENA_ALLOC(arena, (long)(nff * 3 * sizeof(int32_t)));
    double cy = 100.0, cx = 50.0;
    for (int j = 0; j < nh; j++) {
        double z = Hgt * (double)j / (double)(nh - 1);
        for (int i = 0; i < nu; i++) {
            double ang = 2.0 * M_PI * (double)i / (double)nu;
            size_t idx = (size_t)j * (size_t)nu + (size_t)i;
            v[idx * 3 + 0] = (float)z;
            v[idx * 3 + 1] = (float)(cy + R * sin(ang));
            v[idx * 3 + 2] = (float)(cx + R * cos(ang));
        }
    }
    size_t fi = 0;
    for (int j = 0; j < nh - 1; j++) {
        for (int i = 0; i < nu - 1; i++) {   /* open: no i -> (i+1)%nu wrap */
            int32_t a = (int32_t)((size_t)j * (size_t)nu + (size_t)i);
            int32_t b = (int32_t)((size_t)j * (size_t)nu + (size_t)(i + 1));
            int32_t c = (int32_t)((size_t)(j + 1) * (size_t)nu + (size_t)i);
            int32_t d = (int32_t)((size_t)(j + 1) * (size_t)nu + (size_t)(i + 1));
            f[fi * 3 + 0] = a; f[fi * 3 + 1] = b; f[fi * 3 + 2] = c; fi++;
            f[fi * 3 + 0] = b; f[fi * 3 + 1] = d; f[fi * 3 + 2] = c; fi++;
        }
    }
    *out_v = v; *out_nv = nvv; *out_f = f; *out_nf = fi;
}

/* Open Archimedean spiral sheet wound around Z: nphi x nh grid, NOT wrapped.
 * r(phi) = r0 + kgain*phi, phi in [0, turns*2pi]. Vertices in (z, y, x). */
static void build_spiral(Arena_T arena, int nphi, int nh, double r0,
                         double kgain, double turns, double Hgt,
                         float **out_v, size_t *out_nv,
                         int32_t **out_f, size_t *out_nf)
{
    size_t nvv = (size_t)nphi * (size_t)nh;
    size_t nff = (size_t)(nphi - 1) * (size_t)(nh - 1) * 2;
    float   *v = (float *)ARENA_ALLOC(arena, (long)(nvv * 3 * sizeof(float)));
    int32_t *f = (int32_t *)ARENA_ALLOC(arena, (long)(nff * 3 * sizeof(int32_t)));
    double phimax = turns * 2.0 * M_PI;
    double cy = 100.0, cx = 50.0;
    for (int j = 0; j < nh; j++) {
        double z = Hgt * (double)j / (double)(nh - 1);
        for (int i = 0; i < nphi; i++) {
            double phi = phimax * (double)i / (double)(nphi - 1);
            double rr = r0 + kgain * phi;
            size_t idx = (size_t)j * (size_t)nphi + (size_t)i;
            v[idx * 3 + 0] = (float)z;
            v[idx * 3 + 1] = (float)(cy + rr * sin(phi));
            v[idx * 3 + 2] = (float)(cx + rr * cos(phi));
        }
    }
    size_t fi = 0;
    for (int j = 0; j < nh - 1; j++) {
        for (int i = 0; i < nphi - 1; i++) {
            int32_t a = (int32_t)((size_t)j * (size_t)nphi + (size_t)i);
            int32_t b = (int32_t)((size_t)j * (size_t)nphi + (size_t)(i + 1));
            int32_t c = (int32_t)((size_t)(j + 1) * (size_t)nphi + (size_t)i);
            int32_t d = (int32_t)((size_t)(j + 1) * (size_t)nphi + (size_t)(i + 1));
            f[fi * 3 + 0] = a; f[fi * 3 + 1] = b; f[fi * 3 + 2] = c; fi++;
            f[fi * 3 + 0] = b; f[fi * 3 + 1] = d; f[fi * 3 + 2] = c; fi++;
        }
    }
    *out_v = v; *out_nv = nvv; *out_f = f; *out_nf = fi;
}

int Unwrap_selftest(void)
{
    int fails = WindingRegister_selftest();
    Arena_T arena = Arena_new();

    /* (1) Open ~360-deg strip: auto axis = Z, ~1 turn, axial span ~Hgt. */
    {
        float *v; int32_t *f; size_t nvv, nff;
        build_cylinder(arena, 64, 40, 10.0, 50.0, &v, &nvv, &f, &nff);
        UnwrapResult res;
        UnwrapOpts opt = { UNWRAP_AXIS_AUTO };
        int rc = Unwrap_run(arena, v, nvv, f, nff, &opt, &res);
        int ok = (rc == 0);
        /* axis ~ +Z = (1,0,0) in (z,y,x) */
        ok = ok && (fabs((double)res.axis[0]) > 0.9);
        ok = ok && (res.n_components == 1);
        ok = ok && (res.turns > 0.85 && res.turns < 1.15);
        ok = ok && (res.v_span > 45.0 && res.v_span < 55.0);
        if (!ok) {
            fprintf(stderr, "[unwrap selftest] CYLINDER FAIL "
                    "rc=%d axisZ=%.3f ncomp=%d turns=%.3f vspan=%.2f\n",
                    rc, (double)res.axis[0], res.n_components,
                    res.turns, res.v_span);
            fails++;
        } else {
            fprintf(stderr, "[unwrap selftest] cylinder OK "
                    "(turns=%.3f vspan=%.2f)\n", res.turns, res.v_span);
        }
    }

    /* (1s) Gauge sync: two islands on one spiral, island B lifted 3 turns
     * low.  The overlap-offset solver must recover offset +3 through the
     * snap gate and apply it to Phi and turn_correction. */
    {
        enum { NS = 4000, NA = 2000 };
        double *t = (double *)ARENA_ALLOC(arena, NS * sizeof(double));
        double *r = (double *)ARENA_ALLOC(arena, NS * sizeof(double));
        double *Phi = (double *)ARENA_ALLOC(arena, NS * sizeof(double));
        int32_t *comp = (int32_t *)ARENA_ALLOC(arena, NS * sizeof(int32_t));
        int32_t comp_island[2] = { 0, 1 };
        int32_t tc[2] = { 0, 0 };
        double pitch = 9.5;
        size_t synced = 0;
        int ok = 1;
        for (size_t i = 0; i < NS; i++) {
            int in_b = i >= NA;
            double w_true = in_b
                          ? 3.0 + 6.0 * ((double)(i - NA) / (double)NA)
                          : 6.0 * ((double)i / (double)NA);
            double noise = 0.4 * sin((double)i * 12.9898);
            comp[i] = in_b ? 1 : 0;
            t[i] = (i % 2) ? 5.0 : 21.0;
            r[i] = 50.0 + pitch * w_true + noise;
            Phi[i] = 2.0 * M_PI * (w_true - (in_b ? 3.0 : 0.0));
        }
        synced = unwrap_sync_gauges(arena, NS, t, r, Phi, comp, 2,
                                    comp_island, tc, 1, pitch);
        ok = ok && synced == 1;
        ok = ok && (tc[1] - tc[0] == 3 || tc[0] - tc[1] == -3);
        /* after sync, measured turn difference must match radial order */
        for (size_t i = 0; ok && i < 64; i++) {
            size_t a = i * 7 % NA, b = NA + (i * 13 % NA);
            if (fabs(t[a] - t[b]) > 1.0) continue;
            {
                double dw = (Phi[a] - Phi[b]) / (2.0 * M_PI);
                double dk = dw - (r[a] - r[b]) / pitch;
                if (fabs(dk) > 0.35) ok = 0;
            }
        }
        if (!ok) {
            fprintf(stderr,
                    "[unwrap selftest] GAUGE SYNC FAIL synced=%zu tc=%d/%d\n",
                    synced, tc[0], tc[1]);
            fails++;
        } else {
            fprintf(stderr,
                    "[unwrap selftest] gauge sync OK (offset +3 recovered, "
                    "%zu island re-gauged)\n", synced);
        }
    }

    /* (2) Spiral: forced axis Z, N turns, spiral fit recovers b=kgain*2pi. */
    {
        double r0 = 15.0, kgain = 2.0, turns = 3.0, Hgt = 200.0;
        float *v; int32_t *f; size_t nvv, nff;
        build_spiral(arena, 240, 30, r0, kgain, turns, Hgt, &v, &nvv, &f, &nff);
        UnwrapResult res;
        UnwrapOpts opt = { UNWRAP_AXIS_Z };
        int rc = Unwrap_run(arena, v, nvv, f, nff, &opt, &res);
        int ok = (rc == 0);
        ok = ok && (res.n_components == 1);
        ok = ok && (res.turns > 2.7 && res.turns < 3.3);
        /* Winding direction (sign of Phi/b) is arbitrary, so compare |b|. */
        double bexp = kgain * 2.0 * M_PI;   /* r = r0 + kgain*phi = a + b*(phi/2pi) */
        ok = ok && (fabs(fabs(res.spiral_b) - bexp) < 0.15 * bexp);
        ok = ok && (res.spiral_r2 > 0.95);
        if (!ok) {
            fprintf(stderr, "[unwrap selftest] SPIRAL FAIL "
                    "rc=%d ncomp=%d turns=%.3f a=%.3f b=%.3f (bexp=%.3f) r2=%.4f\n",
                    rc, res.n_components, res.turns, res.spiral_a,
                    res.spiral_b, bexp, res.spiral_r2);
            fails++;
        } else {
            fprintf(stderr, "[unwrap selftest] spiral OK "
                    "(turns=%.3f b=%.3f r2=%.4f)\n",
                    res.turns, res.spiral_b, res.spiral_r2);
        }
    }

    /* (2b) keep_phi: phi is returned iff requested, and u == (phi-min)*r_ref. */
    {
        float *v; int32_t *f; size_t nvv, nff;
        build_spiral(arena, 120, 10, 15.0, 2.0, 2.0, 60.0, &v, &nvv, &f, &nff);
        UnwrapResult res;
        UnwrapOpts opt = { UNWRAP_AXIS_Z };
        int rc = Unwrap_run(arena, v, nvv, f, nff, &opt, &res);
        int ok = (rc == 0) && (res.phi == NULL);   /* default: no phi */
        opt.keep_phi = 1;
        rc = Unwrap_run(arena, v, nvv, f, nff, &opt, &res);
        ok = ok && (rc == 0) && (res.phi != NULL) &&
             (res.island != NULL) && (res.continuation_island != NULL) &&
             (res.mesh_component != NULL);
        if (ok) {
            float pmin = res.phi[0];
            for (size_t i = 1; i < nvv; i++) if (res.phi[i] < pmin) pmin = res.phi[i];
            double err = 0.0;
            for (size_t i = 0; i < nvv; i++) {
                double expect = ((double)res.phi[i] - (double)pmin) * res.r_ref;
                double d = fabs((double)res.uv[i * 2 + 0] - expect);
                if (d > err) err = d;
            }
            ok = (err < 1e-2);   /* float phi round-trip vs double internals */
            if (!ok) fprintf(stderr, "[unwrap selftest] KEEP_PHI FAIL max u err=%.6f\n", err);
        }
        if (!ok) {
            fprintf(stderr, "[unwrap selftest] KEEP_PHI FAIL\n");
            fails++;
        } else {
            fprintf(stderr, "[unwrap selftest] keep_phi OK\n");
        }
    }

    /* (2c) Phase-continuation join: two spiral arc segments split by a real
     * hole.  A 0.5-rad hole joins into one material island (the local
     * continuation gates can never bridge it -- 0.5 rad at r~25 is ~12 vox);
     * a 5-rad hole is ambiguous and stays two. */
    for (int wide = 0; wide < 2; wide++) {
        enum { JNPHI = 60, JNH = 8 };
        double r0 = 20.0, kg = 2.0, hole = wide ? 5.0 : 0.5, hgt2 = 40.0;
        double seg = 2.0 * M_PI;    /* each arc spans one turn */
        size_t nseg = (size_t)JNPHI * (size_t)JNH;
        size_t nvv = nseg * 2, nff = 0;
        float *v = (float *)ARENA_ALLOC(
            arena, (long)(nvv * 3 * sizeof(float)));
        int32_t *f = (int32_t *)ARENA_ALLOC(
            arena, (long)((size_t)(JNPHI - 1) * (JNH - 1) * 4 * 3 *
                          sizeof(int32_t)));
        for (int s = 0; s < 2; s++) {
            double phi0 = (double)s * (seg + hole);
            for (int j = 0; j < JNH; j++)
                for (int i = 0; i < JNPHI; i++) {
                    double phi = phi0 + seg * (double)i / (double)(JNPHI - 1);
                    double rr = r0 + kg * phi;
                    size_t idx = (size_t)s * nseg +
                                 (size_t)j * JNPHI + (size_t)i;
                    v[idx * 3 + 0] =
                        (float)(hgt2 * (double)j / (double)(JNH - 1));
                    v[idx * 3 + 1] = (float)(100.0 + rr * sin(phi));
                    v[idx * 3 + 2] = (float)(50.0 + rr * cos(phi));
                }
            for (int j = 0; j + 1 < JNH; j++)
                for (int i = 0; i + 1 < JNPHI; i++) {
                    int32_t a = (int32_t)((size_t)s * nseg +
                                          (size_t)j * JNPHI + (size_t)i);
                    int32_t b = a + 1, c = a + JNPHI, d = c + 1;
                    f[nff * 3 + 0] = a; f[nff * 3 + 1] = b;
                    f[nff * 3 + 2] = c; nff++;
                    f[nff * 3 + 0] = b; f[nff * 3 + 1] = d;
                    f[nff * 3 + 2] = c; nff++;
                }
        }
        {
            UnwrapResult res;
            UnwrapOpts opt = { UNWRAP_AXIS_Z };
            int rc, ok, nisl = 0;
            opt.keep_phi = 1;
            rc = Unwrap_run(arena, v, nvv, f, nff, &opt, &res);
            ok = (rc == 0) && (res.continuation_island != NULL);
            if (ok) {
                int32_t seen[8];
                for (size_t i = 0; i < nvv; i++) {
                    int32_t k = res.continuation_island[i];
                    int known = 0;
                    for (int q = 0; q < nisl; q++)
                        if (seen[q] == k) known = 1;
                    if (!known && nisl < 8) seen[nisl++] = k;
                }
                ok = wide ? (nisl == 2) : (nisl == 1);
            }
            if (!ok) {
                fprintf(stderr,
                        "[unwrap selftest] PHASE-JOIN FAIL hole=%.1f rad "
                        "rc=%d islands=%d (want %d)\n",
                        hole, rc, nisl, wide ? 2 : 1);
                fails++;
            } else {
                fprintf(stderr,
                        "[unwrap selftest] phase-join OK (hole=%.1f rad -> "
                        "%d island(s))\n", hole, nisl);
            }
        }
    }

    /* (3) Degenerate: empty and single-triangle inputs return cleanly. */
    {
        UnwrapResult res;
        int rc_empty = Unwrap_run(arena, NULL, 0, NULL, 0, NULL, &res);
        int ok = (rc_empty == -1);

        float tri_v[9] = { 0,0,0,  1,0,0,  0,1,0 };
        int32_t tri_f[3] = { 0, 1, 2 };
        int rc_tri = Unwrap_run(arena, tri_v, 3, tri_f, 1, NULL, &res);
        ok = ok && (rc_tri == 0) && (res.uv != NULL);
        if (!ok) {
            fprintf(stderr, "[unwrap selftest] DEGENERATE FAIL "
                    "rc_empty=%d rc_tri=%d\n", rc_empty, rc_tri);
            fails++;
        } else {
            fprintf(stderr, "[unwrap selftest] degenerate OK\n");
        }
    }

    Arena_dispose(&arena);
    return fails;
}
