#define _USE_MATH_DEFINES
#include "unwrap.h"
#include "winding_field.h"
#include "winding_register.h"

#include "../common/csr.h"
#include "../common/pca.h"
#include "../common/pipeline_constants.h"

#include <assert.h>
#include <limits.h>
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

typedef struct {
    int32_t component;
    float value, reliability;
} UnwrapFieldSample;

/* Arc-length u map (step 8).  0.05 rad is ~3 degrees, far finer than a
 * wrap; the bin cap bounds a whole-scroll run; a bin under the sample
 * floor inherits its neighbour rather than inventing a radius. */
#define UNWRAP_ARCLEN_BIN_RAD   0.05
#define UNWRAP_ARCLEN_MAX_BINS  65536
#define UNWRAP_ARCLEN_MIN_BIN   8
#define UNWRAP_ARCLEN_MIN_R     1e-3

/* --- Radial-layer site splitting (R1) --------------------------------------
 * A mesh component whose faces FUSE two wraps (a weld bridge, a thick
 * prediction neck) short-circuits the phase lift: the BFS carries one lifted
 * turn across both wraps, and the register can never separate them because a
 * component is ONE correction site.  Measured 2026-09-01 on the 10x weld
 * certificate: 65.2% of radially adjacent wrap pairs shared a lifted turn
 * (88% at r>640) and the lift spanned 25.6 turns against ~88 physical wraps.
 *
 * The split is GEOMETRIC, on the register's own lattice, so it works even
 * when the lift is already fused: vertices are keyed by (component, axial
 * bin, lifted phase bin); each key run is cut into radial layers at gaps
 * over UNWRAP_SPLIT_LAYER_GAP x pitch (a prediction shell is ~1-2 vox thick
 * and stays whole; the inter-wrap clearance is ~7 vox and always cuts);
 * layers link across neighbouring keys only when their wrap invariant
 * rho = r - pitch*q agrees within UNWRAP_SPLIT_LINK_TOL x pitch.  A part
 * whose inherited lift is internally torn by 2*pi lands in different phase
 * keys with rho a full pitch apart, so it simply splits further -- the
 * construction is self-healing, never wrong-joining.  Classes under
 * UNWRAP_SPLIT_MIN_VERTS merge back into their parent's largest class.
 * The result REFINES the true mesh components (out->mesh_component keeps
 * the unsplit identity). */
#define UNWRAP_SPLIT_AXIAL_H    2.0   /* vox per axial bin (matches WR use) */
#define UNWRAP_SPLIT_PHASE_BINS 256   /* bins per turn of lifted phase */
#define UNWRAP_SPLIT_LAYER_GAP  0.55  /* x pitch: radial cut threshold */
#define UNWRAP_SPLIT_LINK_TOL   0.25  /* x pitch: cross-key union tolerance.
                                       * Must sit BELOW half the layer-gap
                                       * cut: a gradual fusion neck emits a
                                       * merged transition layer at rho
                                       * midway between the two wraps
                                       * (+-gap/2 ~ +-3 vox), and a tolerance
                                       * above that bridges the wraps right
                                       * back together (measured: 13,604
                                       * fused keys yielded only 29 sibling
                                       * pairs at 0.35).  Same-wrap links
                                       * run |drho| 0..1.5 vox. */
#define UNWRAP_SPLIT_MIN_VERTS  64    /* absolute floor; see SITE_BUDGET */
#define UNWRAP_SPLIT_SITE_BUDGET 1500  /* effective min class size is
                                       * nv/SITE_BUDGET: relation evidence is
                                       * only reliable between LARGE sites
                                       * (many independent observations per
                                       * pair).  Measured 2026-09-01: 30,952
                                       * sites at 10x shattered into 7,353
                                       * gauge islands whose pairs carried a
                                       * MEDIAN OF ONE relation -- nothing to
                                       * vote with; the 783-piece pile hit
                                       * the same wall at 8-13%% wrong
                                       * relations. */
#define UNWRAP_SPLIT_STRICT_TOL 1.0   /* vox: tier-1 same-wrap link */
/* Fused-only gate: R1 splits ONLY components whose (axial, phase) keys show
 * two radial layers often enough to be a wrap fusion.  Splitting every piece
 * (2026-09-02, 21x3x3 tube) shattered the gauge -- 14 -> 143 islands, 24 ->
 * 179 lifted turns -- because the flood of sub-sites cannot chain; the
 * 4x5x5 end-to-end run died the same way.  A piece that never stacks on
 * itself stays ONE site. */
#define UNWRAP_SPLIT_FUSED_ONLY     1
#define UNWRAP_SPLIT_FUSED_MIN_KEYS 4     /* two-layer keys needed */
#define UNWRAP_SPLIT_FUSED_MIN_FRAC 0.02  /* ... and share of the piece's keys */
#define UNWRAP_SPLIT_NECK_SPAN  0.60  /* x pitch: a tier-1 class whose loose
                                       * links reach two classes this far
                                       * apart in rho is a fusion NECK: the
                                       * transition band of a gradual merge.
                                       * Its loose links are discarded so it
                                       * can never bridge the wraps it
                                       * touches (min-verts folds it into a
                                       * side at the end). */

static int unwrap_cmp_float(const void *a, const void *b)
{
    float x = *(const float *)a, y = *(const float *)b;
    if (x < y) return -1;
    if (x > y) return 1;
    return 0;
}

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

static int unwrap_compare_field_sample(const void *pa, const void *pb)
{
    const UnwrapFieldSample *a = (const UnwrapFieldSample *)pa;
    const UnwrapFieldSample *b = (const UnwrapFieldSample *)pb;
    if (a->component != b->component)
        return a->component < b->component ? -1 : 1;
    if (a->value != b->value) return a->value < b->value ? -1 : 1;
    return 0;
}

static double unwrap_weighted_field_median(
    const UnwrapFieldSample *sample, size_t first, size_t last)
{
    double total = 0.0;
    for (size_t i = first; i < last; i++) total += sample[i].reliability;
    if (!(total > 0.0)) return NAN;
    double cumulative = 0.0;
    for (size_t i = first; i < last; i++) {
        cumulative += sample[i].reliability;
        if (2.0 * cumulative >= total) return sample[i].value;
    }
    return sample[last - 1].value;
}

/* Build the GWN field over a consistently radial-oriented copy of the mesh,
 * then turn its imperfect values into robust per-component unaries.  For an
 * outward normal, w-minus - w-plus is +1 and -w increases toward larger
 * radius, matching q's registered convention. */
static int unwrap_build_field_unary(
    Arena_T arena,
    const float *verts, size_t nv, const int32_t *faces, size_t nf,
    const double *c1, const double *c2,
    const double basis1[3], const double basis2[3], const double *q,
    const int32_t *component, int32_t ncomponents,
    int32_t anchor_component, double pitch,
    double epsilon_option, double beta_option,
    double **out_center, double **out_sigma, double **out_weight,
    float **out_mean, float **out_jump,
    WindingFieldStats *out_build_stats,
    WindingFieldEvalStats *out_eval_stats,
    size_t *out_samples, size_t *out_supported,
    size_t *out_clean, size_t *out_invalid)
{
    enum { UNWRAP_FIELD_MAX_SAMPLES = 2000000 };
    int32_t *oriented = (int32_t *)ARENA_ALLOC(
        arena, nf * 3 * sizeof *oriented);
    memcpy(oriented, faces, nf * 3 * sizeof *oriented);
    double *vote = (double *)ARENA_CALLOC(
        arena, (size_t)ncomponents, sizeof *vote);
    for (size_t f = 0; f < nf; f++) {
        int32_t ia = faces[f*3], ib = faces[f*3+1], ic = faces[f*3+2];
        const float *a = &verts[(size_t)ia * 3];
        const float *b = &verts[(size_t)ib * 3];
        const float *c = &verts[(size_t)ic * 3];
        double ab[3] = { b[0]-a[0], b[1]-a[1], b[2]-a[2] };
        double ac[3] = { c[0]-a[0], c[1]-a[1], c[2]-a[2] };
        double normal[3] = {
            ab[1]*ac[2] - ab[2]*ac[1],
            ab[2]*ac[0] - ab[0]*ac[2],
            ab[0]*ac[1] - ab[1]*ac[0]
        };
        double rc1 = (c1[ia] + c1[ib] + c1[ic]) / 3.0;
        double rc2 = (c2[ia] + c2[ib] + c2[ic]) / 3.0;
        double radial[3] = {
            rc1*basis1[0] + rc2*basis2[0],
            rc1*basis1[1] + rc2*basis2[1],
            rc1*basis1[2] + rc2*basis2[2]
        };
        double normal_norm = sqrt(normal[0]*normal[0] +
                                  normal[1]*normal[1] +
                                  normal[2]*normal[2]);
        double radial_norm = hypot(rc1, rc2);
        if (normal_norm > 0.0 && radial_norm > 0.0)
            vote[component[ia]] += normal[0]*radial[0] +
                                   normal[1]*radial[1] +
                                   normal[2]*radial[2];
    }
    size_t flipped_faces = 0;
    for (size_t f = 0; f < nf; f++) {
        int32_t c = component[faces[f*3]];
        if (vote[c] < 0.0) {
            int32_t swap = oriented[f*3+1];
            oriented[f*3+1] = oriented[f*3+2];
            oriented[f*3+2] = swap;
            flipped_faces++;
        }
    }

    float *normal = (float *)ARENA_CALLOC(arena, nv * 3, sizeof *normal);
    for (size_t f = 0; f < nf; f++) {
        int32_t ia = oriented[f*3], ib = oriented[f*3+1], ic = oriented[f*3+2];
        const float *a = &verts[(size_t)ia * 3];
        const float *b = &verts[(size_t)ib * 3];
        const float *c = &verts[(size_t)ic * 3];
        double ab[3] = { b[0]-a[0], b[1]-a[1], b[2]-a[2] };
        double ac[3] = { c[0]-a[0], c[1]-a[1], c[2]-a[2] };
        float n[3] = {
            (float)(ab[1]*ac[2] - ab[2]*ac[1]),
            (float)(ab[2]*ac[0] - ab[0]*ac[2]),
            (float)(ab[0]*ac[1] - ab[1]*ac[0])
        };
        for (int k = 0; k < 3; k++) {
            normal[(size_t)ia*3 + (size_t)k] += n[k];
            normal[(size_t)ib*3 + (size_t)k] += n[k];
            normal[(size_t)ic*3 + (size_t)k] += n[k];
        }
    }
    float *mean = (float *)ARENA_ALLOC(arena, nv * sizeof *mean);
    float *jump = (float *)ARENA_ALLOC(arena, nv * sizeof *jump);
    WindingFieldOptions field_options;
    WindingField_default_options(&field_options);
    if (beta_option > 0.0 && isfinite(beta_option))
        field_options.beta = beta_option;
    WindingField_T field = NULL;
    if (WindingField_build(arena, verts, nv, oriented, nf, &field_options,
                           &field, out_build_stats) != 0)
        return -1;
    double epsilon = epsilon_option;
    if (!(epsilon > 0.0) || !isfinite(epsilon)) {
        epsilon = pitch > 0.0 && isfinite(pitch) ? 0.10 * pitch : 1.0;
        if (epsilon < 0.25) epsilon = 0.25;
        if (epsilon > 1.5) epsilon = 1.5;
    }
    if (WindingField_evaluate_sides(
            field, verts, normal, nv, epsilon, WINDING_FIELD_AUTO,
            mean, jump, out_eval_stats) != 0)
        return -1;

    size_t clean = 0, invalid = 0;
    for (size_t i = 0; i < nv; i++) {
        if (!isfinite(mean[i]) || !isfinite(jump[i])) invalid++;
        else if (fabs((double)jump[i] - 1.0) <= 0.5) clean++;
    }
    size_t stride = nv > UNWRAP_FIELD_MAX_SAMPLES
                  ? (nv + UNWRAP_FIELD_MAX_SAMPLES - 1) /
                    UNWRAP_FIELD_MAX_SAMPLES : 1;
    size_t capacity = (nv + stride - 1) / stride + (size_t)ncomponents;
    UnwrapFieldSample *sample = (UnwrapFieldSample *)ARENA_ALLOC(
        arena, capacity * sizeof *sample);
    uint8_t *seen = (uint8_t *)ARENA_CALLOC(
        arena, (size_t)ncomponents, sizeof *seen);
    size_t nsample = 0;
    for (size_t i = 0; i < nv; i++) {
        int32_t c = component[i];
        if (i % stride != 0 && seen[c]) continue;
        if (!isfinite(mean[i]) || !isfinite(jump[i])) continue;
        double z = ((double)jump[i] - 1.0) / 0.35;
        double reliability = exp(-0.5 * z * z);
        if (reliability < 0.05) continue;
        if (nsample >= capacity) return -1;
        sample[nsample].component = c;
        sample[nsample].value = (float)(-(double)mean[i] - q[i]);
        sample[nsample].reliability = (float)reliability;
        nsample++;
        seen[c] = 1;
    }
    qsort(sample, nsample, sizeof *sample, unwrap_compare_field_sample);
    double *center = (double *)ARENA_ALLOC(
        arena, (size_t)ncomponents * sizeof *center);
    double *sigma = (double *)ARENA_ALLOC(
        arena, (size_t)ncomponents * sizeof *sigma);
    double *weight = (double *)ARENA_ALLOC(
        arena, (size_t)ncomponents * sizeof *weight);
    for (int32_t c = 0; c < ncomponents; c++) {
        center[c] = NAN;
        sigma[c] = 1.0;
        weight[c] = 0.0;
    }
    size_t supported = 0;
    for (size_t first = 0; first < nsample;) {
        size_t last = first + 1;
        while (last < nsample &&
               sample[last].component == sample[first].component)
            last++;
        double median = unwrap_weighted_field_median(sample, first, last);
        double effective = 0.0;
        for (size_t i = first; i < last; i++) {
            effective += sample[i].reliability;
            sample[i].value = (float)fabs((double)sample[i].value - median);
        }
        qsort(sample + first, last - first, sizeof *sample,
              unwrap_compare_field_sample);
        double mad = unwrap_weighted_field_median(sample, first, last);
        int32_t c = sample[first].component;
        center[c] = median;
        /* MAD captures sampling noise but not the measured global compression
         * of an incomplete field (about 25 physical wraps -> 8 GWN units on
         * PHerc0139).  Keep a one-turn epistemic floor so a sharp-looking but
         * biased field cannot overturn a well-supported relation. */
        sigma[c] = fmax(1.0, 1.4826 * mad);
        double coverage = effective / (double)(last - first);
        weight[c] = fmin(1.0,
            log1p(effective) * coverage / (1.0 + sigma[c]));
        if (last - first < 3)
            weight[c] *= (double)(last - first) / 3.0;
        if (weight[c] > 0.0) supported++;
        first = last;
    }
    if (isfinite(center[anchor_component])) {
        double anchor = center[anchor_component];
        for (int32_t c = 0; c < ncomponents; c++)
            if (isfinite(center[c])) center[c] -= anchor;
    } else {
        for (int32_t c = 0; c < ncomponents; c++) weight[c] = 0.0;
        supported = 0;
    }
    fprintf(stderr,
        "  winding field: backend=%s nodes=%zu boundary=%zu/%zu "
        "flipped_faces=%zu eps=%.3f clean=%zu/%zu (%.1f%%) "
        "unary_components=%zu/%d samples=%zu\n",
        out_eval_stats->backend_used == WINDING_FIELD_BOUNDARY_EXACT
            ? "boundary-exact" :
        out_eval_stats->backend_used == WINDING_FIELD_DIRECT
            ? "direct" : "fast",
        out_build_stats->bvh_nodes, out_build_stats->exterior_edges,
        out_build_stats->exterior_multiplicity, flipped_faces, epsilon,
        clean, nv, 100.0 * (double)clean / (double)nv,
        supported, ncomponents, nsample);
    *out_center = center;
    *out_sigma = sigma;
    *out_weight = weight;
    *out_mean = mean;
    *out_jump = jump;
    *out_samples = nsample;
    *out_supported = supported;
    *out_clean = clean;
    *out_invalid = invalid;
    return 0;
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
        arena, (size_t)((nv + 1) * sizeof(UnwrapSyncKey)));
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
                arena, (size_t)(((size_t)nisl * 4 + 8) * sizeof(int32_t)));
            int32_t *edge_b = (int32_t *)ARENA_ALLOC(
                arena, (size_t)(((size_t)nisl * 4 + 8) * sizeof(int32_t)));
            int32_t *edge_k = (int32_t *)ARENA_ALLOC(
                arena, (size_t)(((size_t)nisl * 4 + 8) * sizeof(int32_t)));
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
                    arena, (size_t)(((size_t)nisl + 1) * sizeof(int32_t)));
                int32_t *root = (int32_t *)ARENA_ALLOC(
                    arena, (size_t)(((size_t)nisl + 1) * sizeof(int32_t)));
                uint8_t *have = (uint8_t *)ARENA_CALLOC(
                    arena, (size_t)nisl + 1, 1);
                uint8_t *bad_root = (uint8_t *)ARENA_CALLOC(
                    arena, (size_t)nisl + 1, 1);
                int32_t *queue = (int32_t *)ARENA_ALLOC(
                    arena, (size_t)(((size_t)nisl + 1) * sizeof(int32_t)));
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

/* ---- Junction-radius island gauging (goal: weld-equivalent placement) -----
 * Measured 2026-08-31 on BOTH the pre-weld pile fit and the welded fit: the
 * phase-continuation join's candidates radius-reject with dr off by near-
 * integer PITCH multiples (dr 22-83 vox against ~1 expected).  Mechanism:
 * each continuation island's internal lift is fine, but the missing arcs
 * BETWEEN islands collapse out of the global winding, so an island that
 * physically lives k turns further out gets packed phi-adjacent to its
 * predecessor.  Welding never fixed this in the current lane (same rejects
 * on welded input); the claims-unwrapper lineage fixed it with radius-
 * anchored winding.  This pass restores that, from LOCAL evidence only
 * (subset-invariance requirement: chopping out the umbilicus or half the
 * grid must not change the surviving islands' relative gauges):
 *   - interval-adjacent island pairs with axial overlap contribute a
 *     junction observation: the MEDIAN over shared 8-vox axial bins of the
 *     per-bin mean-radius difference (same-z, same-theta-window comparison
 *     -- immune to the scroll's eccentricity, unlike whole-window means);
 *   - k = lround((dr_med - pitch*gap/2pi)/pitch) is accepted only when the
 *     residual snaps within UNWRAP_JGAUGE_SNAP of that integer and
 *     |k| <= UNWRAP_JGAUGE_MAX_TURNS; anything weaker ABSTAINS (k=0);
 *   - accepted k accumulate along the interval-sorted chain; the LARGEST
 *     island anchors (deterministic, locally stable); shifts apply to Phi
 *     and turn_correction, so the downstream join sees corrected lifts,
 *     true k-turn holes open where material is genuinely missing, and
 *     placement lays islands at their physical turns. */
enum {
    UNWRAP_JGAUGE_MIN_VERTS = 24,
    UNWRAP_JGAUGE_MIN_BINS = 4,
    UNWRAP_JGAUGE_MAX_TURNS = 12,
    UNWRAP_JGAUGE_MIN_BIN_N = 3,
    /* MEASURE-ONLY (2026-08-31): applying the chained shifts was measured
     * harmful -- interval-consecutive pairing is not sheet adjacency when
     * lifts are broken (the pile fit's lift compressed 15.1 -> 7.6 turns,
     * the welded fit's overspread to 49 against ~15 physical).  The same
     * junction-radius evidence now travels as MULTI-PITCH ORDER relations
     * through the registration's relation build + forest + loop closure
     * (wr_collect_order), which check cycles instead of chaining blindly.
     * The pass is retained as a per-junction measurement/log. */
    UNWRAP_JGAUGE_APPLY = 0
};
#define UNWRAP_JGAUGE_WINDOW 0.45      /* rad each side of the junction */
#define UNWRAP_JGAUGE_AXIAL_BIN 8.0    /* vox; matches WR_AXIAL_BIN */
#define UNWRAP_JGAUGE_SNAP 0.25        /* x pitch residual gate */

static int unwrap_compare_double(const void *pa, const void *pb)
{
    double a = *(const double *)pa, b = *(const double *)pb;
    return a < b ? -1 : (a > b ? 1 : 0);
}

static void unwrap_gauge_islands_by_junction(
    Arena_T arena, size_t nv, const double *t, const double *r,
    double *Phi, const int32_t *comp, int32_t ncomp,
    const int32_t *continuation_island, int32_t *turn_correction,
    int winding_sense, double pitch,
    size_t *out_shifted, size_t *out_abstained)
{
    Arena_Mark mark = Arena_save(arena);
    int32_t nisl = 0;
    *out_shifted = 0;
    *out_abstained = 0;
    for (int32_t c = 0; c < ncomp; c++)
        if (continuation_island[c] >= nisl)
            nisl = continuation_island[c] + 1;
    if (nisl < 2 || !(pitch > 1e-6)) { Arena_restore(arena, mark); return; }
    {
        double *ilo = (double *)ARENA_ALLOC(
            arena, (size_t)((size_t)nisl * sizeof(double)));
        double *ihi = (double *)ARENA_ALLOC(
            arena, (size_t)((size_t)nisl * sizeof(double)));
        size_t *icount = (size_t *)ARENA_CALLOC(
            arena, (size_t)nisl, sizeof(size_t));
        int32_t *order = (int32_t *)ARENA_ALLOC(
            arena, (size_t)((size_t)nisl * sizeof(int32_t)));
        int32_t *shift = (int32_t *)ARENA_CALLOC(
            arena, (size_t)nisl, sizeof(int32_t));
        double tmin = INFINITY, tmax = -INFINITY;
        size_t nbin = 0;
        for (int32_t k = 0; k < nisl; k++) {
            ilo[k] = INFINITY; ihi[k] = -INFINITY;
            order[k] = k;
        }
        for (size_t i = 0; i < nv; i++) {
            int32_t k = continuation_island[comp[i]];
            double w = (double)winding_sense * Phi[i];
            if (k < 0 || k >= nisl) continue;
            if (w < ilo[k]) ilo[k] = w;
            if (w > ihi[k]) ihi[k] = w;
            if (t[i] < tmin) tmin = t[i];
            if (t[i] > tmax) tmax = t[i];
            icount[k]++;
        }
        if (!(tmax >= tmin)) { Arena_restore(arena, mark); return; }
        nbin = (size_t)((tmax - tmin) / UNWRAP_JGAUGE_AXIAL_BIN) + 2;
        {
            /* per-island junction-window radius accumulators, both ends */
            double *hi_sum = (double *)ARENA_CALLOC(
                arena, (size_t)((size_t)nisl * nbin), sizeof(double));
            double *lo_sum = (double *)ARENA_CALLOC(
                arena, (size_t)((size_t)nisl * nbin), sizeof(double));
            uint32_t *hi_n = (uint32_t *)ARENA_CALLOC(
                arena, (size_t)((size_t)nisl * nbin), sizeof(uint32_t));
            uint32_t *lo_n = (uint32_t *)ARENA_CALLOC(
                arena, (size_t)((size_t)nisl * nbin), sizeof(uint32_t));
            double *delta = (double *)ARENA_ALLOC(
                arena, (size_t)(nbin * sizeof(double)));
            for (size_t i = 0; i < nv; i++) {
                int32_t k = continuation_island[comp[i]];
                double w = (double)winding_sense * Phi[i];
                size_t bin = 0;
                if (k < 0 || k >= nisl) continue;
                bin = (size_t)((t[i] - tmin) / UNWRAP_JGAUGE_AXIAL_BIN);
                if (bin >= nbin) bin = nbin - 1;
                if (w > ihi[k] - UNWRAP_JGAUGE_WINDOW) {
                    hi_sum[(size_t)k * nbin + bin] += r[i];
                    hi_n[(size_t)k * nbin + bin]++;
                }
                if (w < ilo[k] + UNWRAP_JGAUGE_WINDOW) {
                    lo_sum[(size_t)k * nbin + bin] += r[i];
                    lo_n[(size_t)k * nbin + bin]++;
                }
            }
            /* interval-sorted island order (insertion sort; nisl is small) */
            for (int32_t a = 1; a < nisl; a++) {
                int32_t key = order[a];
                int32_t b = a - 1;
                while (b >= 0 && ilo[order[b]] > ilo[key]) {
                    order[b + 1] = order[b];
                    b--;
                }
                order[b + 1] = key;
            }
            /* chain of junction snaps between consecutive ELIGIBLE islands */
            {
                int32_t previous = -1;
                int32_t carried = 0;
                for (int32_t s = 0; s < nisl; s++) {
                    int32_t b = order[s];
                    if (icount[b] < (size_t)UNWRAP_JGAUGE_MIN_VERTS) {
                        shift[b] = carried;   /* inherit the chain state */
                        continue;
                    }
                    if (previous >= 0) {
                        int32_t a = previous;
                        size_t nshared = 0;
                        double gap = ilo[b] - ihi[a];
                        double dr_expect =
                            pitch * (gap > 0.0 ? gap : 0.0) / (2.0 * M_PI);
                        for (size_t bin = 0; bin < nbin; bin++) {
                            uint32_t na = hi_n[(size_t)a * nbin + bin];
                            uint32_t nb2 = lo_n[(size_t)b * nbin + bin];
                            if (na < UNWRAP_JGAUGE_MIN_BIN_N ||
                                nb2 < UNWRAP_JGAUGE_MIN_BIN_N)
                                continue;
                            delta[nshared++] =
                                lo_sum[(size_t)b * nbin + bin] / nb2 -
                                hi_sum[(size_t)a * nbin + bin] / na;
                        }
                        if (nshared >= (size_t)UNWRAP_JGAUGE_MIN_BINS) {
                            double dr_med = 0.0, residual = 0.0;
                            long snap = 0;
                            qsort(delta, nshared, sizeof(double),
                                  unwrap_compare_double);
                            dr_med = delta[nshared / 2];
                            snap = lround((dr_med - dr_expect) / pitch);
                            residual = fabs(dr_med - dr_expect -
                                            (double)snap * pitch);
                            if (labs(snap) <= UNWRAP_JGAUGE_MAX_TURNS &&
                                residual <= UNWRAP_JGAUGE_SNAP * pitch) {
                                carried += (int32_t)snap;
                                if (snap != 0) {
                                    (*out_shifted)++;
                                    fprintf(stderr,
                                            "  island gauge: island %d sits "
                                            "%+ld turn(s) out from island %d "
                                            "(junction dr=%.2f expect %.2f "
                                            "over %zu axial bins, resid "
                                            "%.2f)\n",
                                            b, snap, a, dr_med, dr_expect,
                                            nshared, residual);
                                }
                            } else {
                                (*out_abstained)++;
                                fprintf(stderr,
                                        "  island gauge: abstain %d->%d "
                                        "(dr=%.2f expect %.2f bins=%zu "
                                        "snap=%ld resid=%.2f)\n",
                                        a, b, dr_med, dr_expect, nshared,
                                        snap, residual);
                            }
                        } else {
                            (*out_abstained)++;
                        }
                    }
                    shift[b] = carried;
                    previous = b;
                }
            }
            /* anchor: the largest island keeps its gauge (deterministic and
             * subset-stable -- removing other regions cannot move it) */
            {
                int32_t largest = 0;
                int32_t rebase = 0;
                int applied = 0;
                for (int32_t k = 1; k < nisl; k++)
                    if (icount[k] > icount[largest]) largest = k;
                rebase = shift[largest];
                for (int32_t k = 0; k < nisl; k++) {
                    shift[k] -= rebase;
                    if (shift[k] != 0) applied = 1;
                }
                if (UNWRAP_JGAUGE_APPLY && applied) {
                    for (int32_t c = 0; c < ncomp; c++) {
                        int32_t k = continuation_island[c];
                        if (k < 0 || k >= nisl || shift[k] == 0) continue;
                        turn_correction[c] += shift[k];
                    }
                    for (size_t i = 0; i < nv; i++) {
                        int32_t k = continuation_island[comp[i]];
                        if (k < 0 || k >= nisl || shift[k] == 0) continue;
                        Phi[i] += (double)winding_sense * 2.0 * M_PI *
                                  (double)shift[k];
                    }
                }
                fprintf(stderr,
                        "  island gauge: %d island(s), %zu junction snap(s) "
                        "%s, %zu abstained, anchor=island %d\n",
                        nisl, *out_shifted,
                        UNWRAP_JGAUGE_APPLY ? "applied"
                                            : "measured (apply disarmed)",
                        *out_abstained, largest);
            }
        }
    }
    Arena_restore(arena, mark);
}

/* ---- main ----------------------------------------------------------------- */

/* ---- R1 radial-layer splitter (see the constants block above) ---------- */

typedef struct {
    uint64_t key;
    float    r;
    int32_t  idx;
} UnwrapSplitVert;

static int unwrap_split_cmp(const void *pa, const void *pb)
{
    const UnwrapSplitVert *a = (const UnwrapSplitVert *)pa;
    const UnwrapSplitVert *b = (const UnwrapSplitVert *)pb;
    if (a->key != b->key) return a->key < b->key ? -1 : 1;
    if (a->r != b->r) return a->r < b->r ? -1 : 1;
    return a->idx < b->idx ? -1 : (a->idx > b->idx);
}

typedef struct {
    uint64_t key;      /* the run's lattice key */
    int32_t  first;    /* first index in the sorted vert array */
    int32_t  count;    /* verts in this layer */
    float    rho;      /* median wrap invariant r - pitch*q */
    int32_t  comp;     /* parent component */
} UnwrapSplitNode;

typedef struct {
    int32_t a, b;      /* node ids */
    float   drho;      /* rho difference at the link */
} UnwrapSplitLink;

static int unwrap_split_node_cmp_b(const void *pa, const void *pb)
{
    /* ordering B: (comp, pb, axial) so axial-adjacent runs are adjacent */
    const UnwrapSplitNode *a = (const UnwrapSplitNode *)pa;
    const UnwrapSplitNode *b = (const UnwrapSplitNode *)pb;
    uint64_t ka = ((a->key >> 28) << 28)
                | ((a->key & 0x7fffu) << 13)
                | ((a->key >> 15) & 0x1fffu);
    uint64_t kb = ((b->key >> 28) << 28)
                | ((b->key & 0x7fffu) << 13)
                | ((b->key >> 15) & 0x1fffu);
    if (ka != kb) return ka < kb ? -1 : 1;
    if (a->rho != b->rho) return a->rho < b->rho ? -1 : 1;
    return 0;
}

static int unwrap_split_tri_cmp(const void *pa, const void *pb)
{
    const int32_t *a = (const int32_t *)pa;
    const int32_t *b = (const int32_t *)pb;
    if (a[0] != b[0]) return a[0] < b[0] ? -1 : 1;
    if (a[1] != b[1]) return a[1] < b[1] ? -1 : 1;
    if (a[2] != b[2]) return a[2] < b[2] ? -1 : 1;
    return 0;
}

static int32_t unwrap_split_find(int32_t *par, int32_t x)
{
    while (par[x] != x) {
        par[x] = par[par[x]];
        x = par[x];
    }
    return x;
}

/* Refine comp[] into radially coherent layer classes.  Returns the new
 * component count (>= ncomp), or ncomp unchanged when the input is too small
 * or the pitch is degenerate.  self_conflict_keys counts lattice keys whose
 * run held two or more radial layers -- the direct census of lift fusion. */
static int32_t unwrap_split_radial_layers(
    Arena_T arena, size_t nv, const double *t, double tmin,
    const double *r, const double *qturn, double pitch,
    int32_t *comp, int32_t ncomp, size_t *self_conflict_keys,
    WindingSiblingPair **out_sibling, size_t *out_nsibling)
{
    if (self_conflict_keys != NULL) *self_conflict_keys = 0;
    if (out_sibling != NULL) *out_sibling = NULL;
    if (out_nsibling != NULL) *out_nsibling = 0;
    if (nv < 2 || ncomp <= 0 || !(pitch > 1e-6)) return ncomp;
    Arena_Mark mark = Arena_save(arena);
    UnwrapSplitVert *sv = (UnwrapSplitVert *)ARENA_ALLOC(
        arena, (size_t)(nv * sizeof *sv));
    double cut = UNWRAP_SPLIT_LAYER_GAP * pitch;
    double tol = UNWRAP_SPLIT_LINK_TOL * pitch;
    size_t i = 0, nnode = 0, nrun = 0, conflict = 0;
    int32_t result = ncomp;
    for (i = 0; i < nv; i++) {
        long axial = (long)((t[i] - tmin) / UNWRAP_SPLIT_AXIAL_H);
        double qq = qturn[i];
        long turn = (long)floor(qq + 1e-9);
        long bin = (long)((qq - (double)turn) *
                          (double)UNWRAP_SPLIT_PHASE_BINS);
        long pb = (turn + 64) * UNWRAP_SPLIT_PHASE_BINS + bin;
        if (axial < 0) axial = 0;
        if (axial > 0x1fff) axial = 0x1fff;
        if (pb < 0) pb = 0;
        if (pb > 0x7fff) pb = 0x7fff;
        sv[i].key = ((uint64_t)(uint32_t)comp[i] << 28)
                  | ((uint64_t)axial << 15)
                  | (uint64_t)pb;
        sv[i].r = (float)r[i];
        sv[i].idx = (int32_t)i;
    }
    qsort(sv, nv, sizeof *sv, unwrap_split_cmp);
    for (i = 0; i < nv; i++) {
        int newrun = i == 0 || sv[i].key != sv[i - 1].key;
        int newlayer = newrun ||
            (double)sv[i].r - (double)sv[i - 1].r > cut;
        if (newrun) nrun++;
        if (newlayer) nnode++;
    }
    (void)nrun;
    /* per-component fusion evidence: keys (runs) with >= 2 radial layers */
    uint8_t *fused = (uint8_t *)ARENA_CALLOC(arena, (size_t)ncomp, 1);
    {
        size_t *nkeys = (size_t *)ARENA_CALLOC(arena, (size_t)ncomp,
                                               sizeof(size_t));
        size_t *nconf = (size_t *)ARENA_CALLOC(arena, (size_t)ncomp,
                                               sizeof(size_t));
        size_t layers = 0;
        int32_t rc0 = -1, c2 = 0, nfused = 0;
        for (i = 0; i <= nv; i++) {
            int newrun = i == 0 || i == nv || sv[i].key != sv[i - 1].key;
            if (newrun) {
                if (rc0 >= 0) {
                    nkeys[rc0]++;
                    if (layers >= 2) nconf[rc0]++;
                }
                if (i == nv) break;
                rc0 = comp[sv[i].idx];
                layers = 1;
            } else if ((double)sv[i].r - (double)sv[i - 1].r > cut) {
                layers++;
            }
        }
        for (c2 = 0; c2 < ncomp; c2++) {
            if (!UNWRAP_SPLIT_FUSED_ONLY ||
                (nconf[c2] >= UNWRAP_SPLIT_FUSED_MIN_KEYS &&
                 (double)nconf[c2] >=
                     UNWRAP_SPLIT_FUSED_MIN_FRAC * (double)nkeys[c2]))
                fused[c2] = 1;
            nfused += fused[c2];
        }
        fprintf(stderr, "  radial-site split: fused-only=%d -> %d of %d "
                "components qualify (>= %d two-layer keys and >= %.0f%%)"
                "%c", UNWRAP_SPLIT_FUSED_ONLY, nfused, ncomp,
                UNWRAP_SPLIT_FUSED_MIN_KEYS,
                100.0 * UNWRAP_SPLIT_FUSED_MIN_FRAC, 10);
    }
    {
        UnwrapSplitNode *node = (UnwrapSplitNode *)ARENA_ALLOC(
            arena, (size_t)(nnode * sizeof *node));
        int32_t *vnode = (int32_t *)ARENA_ALLOC(
            arena, (size_t)(nv * sizeof *vnode));
        int32_t *par = NULL, *root_of = NULL;
        size_t at = 0, run_first_node = 0;
        for (i = 0; i < nv; i++) {
            int newrun = i == 0 || sv[i].key != sv[i - 1].key;
            int newlayer = newrun ||
                (double)sv[i].r - (double)sv[i - 1].r > cut;
            if (newlayer) {
                if (newrun) {
                    if (at > run_first_node + 1) conflict++;
                    run_first_node = at;
                }
                node[at].key = sv[i].key;
                node[at].first = (int32_t)i;
                node[at].count = 0;
                node[at].comp = comp[sv[i].idx];
                at++;
            }
            node[at - 1].count++;
            vnode[sv[i].idx] = (int32_t)(at - 1);
        }
        if (at > run_first_node + 1) conflict++;
        assert(at == nnode);
        for (i = 0; i < nnode; i++) {
            int32_t mid = node[i].first + node[i].count / 2;
            int32_t vi = sv[mid].idx;
            node[i].rho = (float)(r[vi] - pitch * qturn[vi]);
        }
        par = (int32_t *)ARENA_ALLOC(arena, (size_t)(nnode * sizeof *par));
        for (i = 0; i < nnode; i++) par[i] = (int32_t)i;
        UnwrapSplitLink *link = (UnwrapSplitLink *)ARENA_ALLOC(
            arena, (size_t)((2 * nnode + 16) * sizeof *link));
        size_t nlink = 0, linkcap = 2 * nnode + 16;
        /* ordering A: nodes already grouped (comp, axial, pb); union layers
         * of pb-adjacent runs whose rho agree */
        {
            size_t a0 = 0, b0 = 0;
            for (a0 = 0; a0 < nnode; ) {
                size_t a1 = a0;
                while (a1 < nnode && node[a1].key == node[a0].key) a1++;
                b0 = a1;
                if (b0 < nnode) {
                    uint64_t ka = node[a0].key, kb = node[b0].key;
                    if ((ka >> 15) == (kb >> 15) &&
                        (kb & 0x7fffu) == (ka & 0x7fffu) + 1u) {
                        size_t b1 = b0, x = a0, y = b0;
                        while (b1 < nnode && node[b1].key == node[b0].key)
                            b1++;
                        while (x < a1 && y < b1) {
                            double dd = (double)node[x].rho -
                                        (double)node[y].rho;
                            if (dd > tol) { y++; continue; }
                            if (dd < -tol) { x++; continue; }
                            if (nlink < linkcap) {
                                link[nlink].a = (int32_t)x;
                                link[nlink].b = (int32_t)y;
                                link[nlink].drho = (float)dd;
                                nlink++;
                            }
                            if ((double)node[x].rho < (double)node[y].rho)
                                x++;
                            else
                                y++;
                        }
                    }
                }
                a0 = a1;
            }
        }
        /* ordering B: axial-adjacent runs */
        {
            UnwrapSplitNode *nb = (UnwrapSplitNode *)ARENA_ALLOC(
                arena, (size_t)(nnode * sizeof *nb));
            int32_t *borig = (int32_t *)ARENA_ALLOC(
                arena, (size_t)(nnode * sizeof *borig));
            size_t a0 = 0;
            memcpy(nb, node, nnode * sizeof *nb);
            for (i = 0; i < nnode; i++) nb[i].first = (int32_t)i;
            qsort(nb, nnode, sizeof *nb, unwrap_split_node_cmp_b);
            for (i = 0; i < nnode; i++) borig[i] = nb[i].first;
            for (a0 = 0; a0 < nnode; ) {
                size_t a1 = a0, b0 = 0, b1 = 0;
                while (a1 < nnode && nb[a1].key == nb[a0].key) a1++;
                b0 = a1;
                if (b0 < nnode) {
                    uint64_t ka = nb[a0].key, kb = nb[b0].key;
                    uint64_t ca = ka >> 28, cb = kb >> 28;
                    uint64_t qa = ka & 0x7fffu, qb = kb & 0x7fffu;
                    uint64_t xa = (ka >> 15) & 0x1fffu;
                    uint64_t xb = (kb >> 15) & 0x1fffu;
                    if (ca == cb && qa == qb && xb == xa + 1u) {
                        size_t x = a0, y = b0;
                        b1 = b0;
                        while (b1 < nnode && nb[b1].key == nb[b0].key) b1++;
                        while (x < a1 && y < b1) {
                            double dd = (double)nb[x].rho -
                                        (double)nb[y].rho;
                            if (dd > tol) { y++; continue; }
                            if (dd < -tol) { x++; continue; }
                            if (nlink < linkcap) {
                                link[nlink].a = borig[x];
                                link[nlink].b = borig[y];
                                link[nlink].drho = (float)dd;
                                nlink++;
                            }
                            if ((double)nb[x].rho < (double)nb[y].rho) x++;
                            else y++;
                        }
                    }
                }
                a0 = a1;
            }
        }
        /* tier 1: strict same-wrap links only */
        for (i = 0; i < nlink; i++) {
            if (fabs((double)link[i].drho) > UNWRAP_SPLIT_STRICT_TOL)
                continue;
            {
                int32_t rx = unwrap_split_find(par, link[i].a);
                int32_t ry = unwrap_split_find(par, link[i].b);
                if (rx != ry) par[ry] = rx;
            }
        }
        /* tier 2: loose links may extend a wrap through rho drift, but a
         * tier-1 class whose loose links span both sides of a fusion is the
         * NECK transition band; its loose links are discarded outright. */
        {
            float *lo2 = (float *)ARENA_ALLOC(
                arena, (size_t)(nnode * sizeof(float)));
            float *hi2 = (float *)ARENA_ALLOC(
                arena, (size_t)(nnode * sizeof(float)));
            double span_gate = UNWRAP_SPLIT_NECK_SPAN * pitch;
            for (i = 0; i < nnode; i++) {
                lo2[i] = 1e30f;
                hi2[i] = -1e30f;
            }
            for (i = 0; i < nlink; i++) {
                double ad = fabs((double)link[i].drho);
                if (ad <= UNWRAP_SPLIT_STRICT_TOL) continue;
                {
                    int32_t ra = unwrap_split_find(par, link[i].a);
                    int32_t rb = unwrap_split_find(par, link[i].b);
                    float pa2 = node[link[i].b].rho;
                    float pb2 = node[link[i].a].rho;
                    if (ra == rb) continue;
                    if (pa2 < lo2[ra]) lo2[ra] = pa2;
                    if (pa2 > hi2[ra]) hi2[ra] = pa2;
                    if (pb2 < lo2[rb]) lo2[rb] = pb2;
                    if (pb2 > hi2[rb]) hi2[rb] = pb2;
                }
            }
            for (i = 0; i < nlink; i++) {
                double ad = fabs((double)link[i].drho);
                if (ad <= UNWRAP_SPLIT_STRICT_TOL) continue;
                {
                    int32_t ra = unwrap_split_find(par, link[i].a);
                    int32_t rb = unwrap_split_find(par, link[i].b);
                    if (ra == rb) continue;
                    if ((double)hi2[ra] - (double)lo2[ra] > span_gate)
                        continue;      /* a-side class is a neck */
                    if ((double)hi2[rb] - (double)lo2[rb] > span_gate)
                        continue;      /* b-side class is a neck */
                    par[unwrap_split_find(par, rb)] =
                        unwrap_split_find(par, ra);
                }
            }
        }
        root_of = (int32_t *)ARENA_ALLOC(
            arena, (size_t)(nnode * sizeof *root_of));
        for (i = 0; i < nnode; i++)
            root_of[i] = unwrap_split_find(par, (int32_t)i);
        {
            size_t *rsize = (size_t *)ARENA_CALLOC(
                arena, nnode, sizeof(size_t));
            int32_t *best_of_comp = (int32_t *)ARENA_ALLOC(
                arena, (size_t)((size_t)ncomp * sizeof(int32_t)));
            int32_t *dense = (int32_t *)ARENA_ALLOC(
                arena, (size_t)(nnode * sizeof(int32_t)));
            int32_t *newcomp = (int32_t *)ARENA_ALLOC(
                arena, (size_t)(nv * sizeof(int32_t)));
            int32_t c = 0, next = 0;
            for (i = 0; i < nnode; i++)
                rsize[root_of[i]] += (size_t)node[i].count;
            for (c = 0; c < ncomp; c++) best_of_comp[c] = -1;
            for (i = 0; i < nnode; i++) {
                int32_t rr = root_of[i];
                int32_t pc = node[i].comp;
                if (best_of_comp[pc] < 0 ||
                    rsize[rr] > rsize[best_of_comp[pc]])
                    best_of_comp[pc] = rr;
            }
            {
                size_t min_verts = nv / UNWRAP_SPLIT_SITE_BUDGET;
                if (min_verts < UNWRAP_SPLIT_MIN_VERTS)
                    min_verts = UNWRAP_SPLIT_MIN_VERTS;
                for (i = 0; i < nnode; i++) {
                    int32_t rr = root_of[i];
                    /* a FUSED piece may split down to the absolute floor:
                     * its second layer is small (a neck's worth), and the
                     * site budget was sized for splitting every piece */
                    size_t lim = fused[node[i].comp]
                                     ? (size_t)UNWRAP_SPLIT_MIN_VERTS
                                     : min_verts;
                    if ((rsize[rr] < lim || !fused[node[i].comp]) &&
                        best_of_comp[node[i].comp] >= 0)
                        root_of[i] = best_of_comp[node[i].comp];
                }
            }
            for (i = 0; i < nnode; i++) dense[i] = -1;
            for (i = 0; i < nnode; i++)
                if (dense[root_of[i]] < 0) dense[root_of[i]] = next++;
            for (i = 0; i < nv; i++)
                newcomp[i] = dense[root_of[vnode[i]]];
            /* sibling votes: consecutive layers inside one key run are
             * geometrically adjacent wraps.  The turn delta is measured per
             * key from the lifted phases; the local wrap spacing is taken
             * from the RUN itself (median consecutive gap) so the vote
             * abstains where a wrap is locally missing -- and the global
             * pitch is only the two-layer fallback, which keeps the guard
             * honest where the true pitch has plateaued wider (the pitch
             * table runs 9.75..18.5 across the 0139 radial range). */
            if (out_sibling != NULL && out_nsibling != NULL) {
                size_t cap = nnode + 1, nt = 0, a0 = 0;
                int32_t *tri = (int32_t *)malloc(cap * 3 * sizeof(int32_t));
                if (tri != NULL) {
                    for (a0 = 0; a0 < nnode; ) {
                        size_t a1 = a0, k2 = 0, ng = 0;
                        double gaps[64];
                        double local = pitch;
                        while (a1 < nnode && node[a1].key == node[a0].key)
                            a1++;
                        for (k2 = a0; k2 + 1 < a1 && ng < 64; k2++) {
                            int32_t m1 = node[k2].first +
                                         node[k2].count / 2;
                            int32_t m2 = node[k2 + 1].first +
                                         node[k2 + 1].count / 2;
                            gaps[ng++] = (double)sv[m2].r -
                                         (double)sv[m1].r;
                        }
                        if (ng >= 2) {
                            size_t g1 = 0, g2 = 0;
                            for (g1 = 1; g1 < ng; g1++) {
                                double tv = gaps[g1];
                                for (g2 = g1;
                                     g2 > 0 && gaps[g2 - 1] > tv; g2--)
                                    gaps[g2] = gaps[g2 - 1];
                                gaps[g2] = tv;
                            }
                            local = gaps[ng / 2];
                        }
                        for (k2 = a0; k2 + 1 < a1; k2++) {
                            int32_t m1 = node[k2].first +
                                         node[k2].count / 2;
                            int32_t m2 = node[k2 + 1].first +
                                         node[k2 + 1].count / 2;
                            int32_t v1 = sv[m1].idx, v2 = sv[m2].idx;
                            double gap = (double)sv[m2].r -
                                         (double)sv[m1].r;
                            int32_t ca = dense[root_of[k2]];
                            int32_t cb = dense[root_of[k2 + 1]];
                            long tgt = 0;
                            if (ca == cb) continue;
                            if (!(gap > 5.5) || gap > 1.45 * local)
                                continue;
                            tgt = lround(qturn[v1] + 1.0 - qturn[v2]);
                            if (tgt < INT32_MIN || tgt > INT32_MAX)
                                continue;
                            if (nt < cap) {
                                tri[nt * 3 + 0] = ca;
                                tri[nt * 3 + 1] = cb;
                                tri[nt * 3 + 2] = (int32_t)tgt;
                                nt++;
                            }
                        }
                        a0 = a1;
                    }
                }
                /* aggregate the votes by (inner, outer): mode target wins */
                if (tri != NULL && nt > 0) {
                    size_t w2 = 0, npair = 0;
                    WindingSiblingPair *sib = NULL;
                    qsort(tri, nt, 3 * sizeof(int32_t),
                          unwrap_split_tri_cmp);
                    for (w2 = 0; w2 < nt; w2++)
                        if (w2 == 0 || tri[w2*3] != tri[w2*3-3] ||
                            tri[w2*3+1] != tri[w2*3-2])
                            npair++;
                    sib = (WindingSiblingPair *)malloc(
                        npair * sizeof *sib);
                    if (sib != NULL) {
                        size_t at2 = 0, f2 = 0;
                        for (f2 = 0; f2 < nt; ) {
                            size_t l2 = f2, best = 0, bc = 0, tc = 0;
                            while (l2 < nt && tri[l2*3] == tri[f2*3] &&
                                   tri[l2*3+1] == tri[f2*3+1])
                                l2++;
                            for (w2 = f2; w2 < l2; ) {
                                size_t m3 = w2;
                                while (m3 < l2 &&
                                       tri[m3*3+2] == tri[w2*3+2])
                                    m3++;
                                if (m3 - w2 > bc) {
                                    bc = m3 - w2;
                                    best = w2;
                                }
                                w2 = m3;
                            }
                            tc = l2 - f2;
                            sib[at2].inner = tri[f2*3];
                            sib[at2].outer = tri[f2*3+1];
                            sib[at2].target = tri[best*3+2];
                            sib[at2].mode_keys = (int32_t)bc;
                            sib[at2].total_keys = (int32_t)tc;
                            at2++;
                            f2 = l2;
                        }
                        *out_sibling = sib;
                        *out_nsibling = at2;
                    }
                    free(tri);
                } else if (tri != NULL) {
                    free(tri);
                }
            }
            for (i = 0; i < nv; i++) comp[i] = newcomp[i];
            if (self_conflict_keys != NULL)
                *self_conflict_keys = conflict;
            result = next;
        }
    }
    Arena_restore(arena, mark);
    return result;
}



/* Frame-independent winding-sense vote.  Bucket vertices by (axial bin of
 * UNWRAP_SENSE_AXIAL_BIN vox, theta bin of 1 degree); within a bucket sort
 * by radius; consecutive SAME-component vertices whose radial gap lies in
 * [0.35, 2.25] x pitch are adjacent wraps of one sheet, and the sign of
 * (Phi_outer - Phi_inner) is the sense.  Needs neither the spiral fit nor
 * the axis position beyond the local radius ordering, so it survives an
 * axis error that scrambles the spiral slope. */
#define UNWRAP_SENSE_AXIAL_BIN 8.0
#define UNWRAP_SENSE_MIN_PAIRS 1000
#define UNWRAP_SENSE_FAIL_AGREE 0.80

typedef struct { uint64_t key; float r; int32_t idx; } UnwrapSenseKey;

static int unwrap_sense_cmp(const void *a, const void *b)
{
    const UnwrapSenseKey *x = (const UnwrapSenseKey *)a;
    const UnwrapSenseKey *y = (const UnwrapSenseKey *)b;
    if (x->key != y->key) return x->key < y->key ? -1 : 1;
    if (x->r != y->r) return x->r < y->r ? -1 : 1;
    return x->idx < y->idx ? -1 : (x->idx > y->idx);
}

static int unwrap_vote_sense(Arena_T arena, size_t nv, const double *t,
                             double tmin, const double *r, const double *theta,
                             const double *Phi, const int32_t *comp,
                             double pitch, double *out_agree, size_t *out_n)
{
    Arena_Mark mark = Arena_save(arena);
    UnwrapSenseKey *k = (UnwrapSenseKey *)ARENA_ALLOC(
        arena, (size_t)(nv * sizeof *k));
    size_t i = 0, n = 0;
    long votes = 0, npairs = 0;
    *out_agree = 0.0;
    *out_n = 0;
    if (k == NULL || !(pitch > 1e-6)) { Arena_restore(arena, mark); return 0; }
    for (i = 0; i < nv; i++) {
        long ab = (long)((t[i] - tmin) / UNWRAP_SENSE_AXIAL_BIN);
        long tb = (long)((theta[i] + M_PI) * (180.0 / M_PI));
        if (ab < 0) ab = 0;
        if (tb < 0) tb = 0;
        if (tb > 359) tb = 359;
        k[n].key = ((uint64_t)ab << 9) | (uint64_t)tb;
        k[n].r = (float)r[i];
        k[n].idx = (int32_t)i;
        n++;
    }
    qsort(k, n, sizeof *k, unwrap_sense_cmp);
    for (i = 1; i < n; i++) {
        double gap = (double)k[i].r - (double)k[i - 1].r;
        double dphi = 0.0;
        if (k[i].key != k[i - 1].key) continue;
        if (comp[k[i].idx] != comp[k[i - 1].idx]) continue;
        if (gap < 0.35 * pitch || gap > 2.25 * pitch) continue;
        dphi = Phi[k[i].idx] - Phi[k[i - 1].idx];
        if (fabs(dphi) < 0.5) continue;      /* same wrap seen twice */
        votes += dphi > 0.0 ? 1 : -1;
        npairs++;
    }
    Arena_restore(arena, mark);
    *out_n = (size_t)npairs;
    if (npairs == 0) return 0;
    *out_agree = fabs((double)votes) / (double)npairs;
    return votes > 0 ? 1 : (votes < 0 ? -1 : 0);
}

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
    out->uv = (float *)ARENA_ALLOC(arena, (size_t)(nv * 2 * sizeof(float)));
    out->winding_index = (float *)ARENA_ALLOC(
        arena, nv * sizeof *out->winding_index);
    out->winding_confidence = (float *)ARENA_ALLOC(
        arena, nv * sizeof *out->winding_confidence);
    int keep_phi = (opts != NULL && opts->keep_phi != 0);
    if (keep_phi) {
        out->phi = (float *)ARENA_ALLOC(arena, (size_t)(nv * sizeof(float)));
        out->island = (int32_t *)ARENA_ALLOC(
            arena, (size_t)(nv * sizeof(int32_t)));
        out->continuation_island = (int32_t *)ARENA_ALLOC(
            arena, (size_t)(nv * sizeof(int32_t)));
        out->mesh_component = (int32_t *)ARENA_ALLOC(
            arena, (size_t)(nv * sizeof(int32_t)));
        out->field_winding = (float *)ARENA_ALLOC(
            arena, nv * sizeof *out->field_winding);
        out->field_jump = (float *)ARENA_ALLOC(
            arena, nv * sizeof *out->field_jump);
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
    double *t  = (double *)ARENA_ALLOC(arena, (size_t)(nv * sizeof(double)));
    double *c1 = (double *)ARENA_ALLOC(arena, (size_t)(nv * sizeof(double)));
    double *c2 = (double *)ARENA_ALLOC(arena, (size_t)(nv * sizeof(double)));
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
    double *cb1 = (double *)ARENA_CALLOC(arena, (size_t)K, sizeof(double));
    double *cb2 = (double *)ARENA_CALLOC(arena, (size_t)K, sizeof(double));
    double *cnt = (double *)ARENA_CALLOC(arena, (size_t)K, sizeof(double));
    char   *fil = (char   *)ARENA_CALLOC(arena, (size_t)K, sizeof(char));
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
    double *s1 = (double *)ARENA_ALLOC(arena, (size_t)((size_t)K * sizeof(double)));
    double *s2 = (double *)ARENA_ALLOC(arena, (size_t)((size_t)K * sizeof(double)));
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
    double *r     = (double *)ARENA_ALLOC(arena, (size_t)(nv * sizeof(double)));
    double *theta = (double *)ARENA_ALLOC(arena, (size_t)(nv * sizeof(double)));
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

    double  *Phi   = (double  *)ARENA_ALLOC(arena, (size_t)(nv * sizeof(double)));
    int32_t *comp  = (int32_t *)ARENA_ALLOC(arena, (size_t)(nv * sizeof(int32_t)));
    int32_t *queue = (int32_t *)ARENA_ALLOC(arena, (size_t)(nv * sizeof(int32_t)));
    for (size_t i = 0; i < nv; i++) comp[i] = -1;

    RIdx *ri = (RIdx *)ARENA_ALLOC(arena, (size_t)(nv * sizeof(RIdx)));
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
    int32_t *csize = (int32_t *)ARENA_CALLOC(arena, (size_t)ncomp,
                                             sizeof(int32_t));
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
    int winding_sense_spiral = sb < 0.0 ? -1 : 1;
    int winding_sense = winding_sense_spiral;
    if (opts != NULL && opts->winding_sense != 0)
        winding_sense = opts->winding_sense > 0 ? 1 : -1;   /* config pin */
    {
        double vote_agree = 0.0;
        size_t vote_n = 0;
        double vote_pitch = wrap_spacing > 0.0 ? wrap_spacing : fabs(sb);
        int vote = unwrap_vote_sense(arena, nv, t, tmin, r, theta, Phi, comp,
                                     vote_pitch, &vote_agree, &vote_n);
        out->winding_sense_vote = vote;
        out->winding_sense_vote_agree = vote_agree;
        out->winding_sense_vote_n = vote_n;
        out->winding_sense_spiral = winding_sense_spiral;
        fprintf(stderr, "  winding sense: %s=%+d spiral=%+d (r2 %.3f) "
                "ray_vote=%+d (agree %.2f, n=%zu)\n",
                (opts != NULL && opts->winding_sense != 0) ? "pinned" : "auto",
                winding_sense, winding_sense_spiral, sr2, vote, vote_agree,
                vote_n);
        if (vote != 0 && vote != winding_sense &&
            vote_n >= UNWRAP_SENSE_MIN_PAIRS &&
            vote_agree >= UNWRAP_SENSE_FAIL_AGREE) {
            if (opts != NULL && opts->winding_sense != 0) {
                fprintf(stderr, "unwrap: pinned winding sense %+d contradicts "
                        "the ray vote %+d (agree %.2f over %zu pairs): "
                        "refusing (fail closed)\n",
                        winding_sense, vote, vote_agree, vote_n);
                return -1;
            }
            fprintf(stderr, "unwrap: WARNING auto winding sense %+d contradicts "
                    "the ray vote %+d (agree %.2f over %zu pairs); using the "
                    "vote\n", winding_sense, vote, vote_agree, vote_n);
            winding_sense = vote;
        }
    }
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
        arena, (size_t)(nv * sizeof(double)));
    for (size_t i = 0; i < nv; i++)
        qturn[i] = (double)winding_sense * Phi[i] / (2.0 * M_PI);
    double *field_center = NULL, *field_sigma = NULL, *field_weight = NULL;
    float *field_mean = NULL, *field_jump = NULL;
    WindingRegisterFieldUnary field_unary;
    const WindingRegisterFieldUnary *field_unary_ptr = NULL;
    WindingFieldStats field_build_stats;
    WindingFieldEvalStats field_eval_stats;
    memset(&field_unary, 0, sizeof field_unary);
    memset(&field_build_stats, 0, sizeof field_build_stats);
    memset(&field_eval_stats, 0, sizeof field_eval_stats);
    int field_mode = opts != NULL ? opts->winding_field_mode : 0;
    double registration_pitch = wrap_spacing > 0.0
                              ? wrap_spacing : fabs(sb);
    /* --- 6b. R1 radial-layer site splitting (see constants block). ---
     * Refines comp[] so the register can correct PER WRAP LAYER; the true
     * mesh-component identity is preserved separately for the result. */
    int32_t *comp_true = (int32_t *)ARENA_ALLOC(
        arena, (size_t)(nv * sizeof(int32_t)));
    WindingSiblingPair *split_sibling = NULL;
    size_t n_split_sibling = 0;
    memcpy(comp_true, comp, nv * sizeof(int32_t));
    {
        size_t split_conflict_keys = 0;
        int32_t ncomp_split = ncomp;
        if (opts != NULL && opts->radial_site_split)
            ncomp_split = unwrap_split_radial_layers(
                arena, nv, t, tmin, r, qturn, registration_pitch,
                comp, ncomp, &split_conflict_keys,
                &split_sibling, &n_split_sibling);
        if (ncomp_split > ncomp) {
            fprintf(stderr,
                    "  unwrap split: %d mesh comps -> %d radial-layer sites "
                    "(%zu fused lattice keys, %zu sibling pairs)",
                    ncomp, ncomp_split, split_conflict_keys,
                    n_split_sibling);
            fputc(10, stderr);
            ncomp = ncomp_split;
            csize = (int32_t *)ARENA_CALLOC(arena, (size_t)ncomp,
                                            sizeof(int32_t));
            for (size_t si = 0; si < nv; si++) csize[comp[si]]++;
            big = 0;
            for (int32_t cs2 = 1; cs2 < ncomp; cs2++)
                if (csize[cs2] > csize[big]) big = cs2;
        }
    }
    /* Exact parent-overlap evidence is aggregated at the actual registration
     * site, never at a relation island.  A relation island may contain many
     * radial layers with different integer transformations (the 21x audit
     * measured only 43--99% agreement per island, but 100% per mesh site).
     * Any disagreement inside one site is therefore a real fused/ambiguous
     * component and is refused instead of being hidden by a median. */
    WindingRegisterBoundary boundary;
    const WindingRegisterBoundary *boundary_ptr = NULL;
    int32_t *boundary_correction = NULL, *boundary_lineage = NULL;
    const int boundary_mode = opts != NULL &&
                              opts->boundary_winding != NULL;
    memset(&boundary, 0, sizeof boundary);
    if (boundary_mode) {
        if (opts->boundary_material == NULL) {
            fprintf(stderr, "unwrap: projective boundary winding requires "
                    "boundary material lineage\n");
            free(split_sibling);
            Arena_restore(arena, mark);
            return -1;
        }
        boundary_correction = (int32_t *)ARENA_ALLOC(
            arena, (size_t)((size_t)ncomp * sizeof *boundary_correction));
        boundary_lineage = (int32_t *)ARENA_ALLOC(
            arena, (size_t)((size_t)ncomp * sizeof *boundary_lineage));
        for (int32_t c = 0; c < ncomp; c++) {
            boundary_correction[c] = INT32_MIN;
            boundary_lineage[c] = -1;
        }
        size_t boundary_vertices = 0;
        for (size_t i = 0; i < nv; i++) {
            double parent_q = (double)opts->boundary_winding[i];
            int32_t parent_lineage = opts->boundary_material[i];
            int32_t c = comp[i];
            if (!isfinite(parent_q)) {
                if (parent_lineage >= 0) {
                    fprintf(stderr, "unwrap: boundary lineage without winding "
                            "at vertex %zu\n", i);
                    free(split_sibling);
                    Arena_restore(arena, mark);
                    return -1;
                }
                continue;
            }
            double delta = parent_q - qturn[i];
            double rounded = nearbyint(delta);
            if (parent_lineage < 0 || !isfinite(delta) ||
                rounded < (double)INT32_MIN || rounded > (double)INT32_MAX ||
                fabs(delta - rounded) > 0.05) {
                fprintf(stderr, "unwrap: non-integral/incomplete parent "
                        "boundary at vertex %zu (dq=%+.9g lineage=%d)\n",
                        i, delta, parent_lineage);
                free(split_sibling);
                Arena_restore(arena, mark);
                return -1;
            }
            int32_t k = (int32_t)rounded;
            if ((boundary_correction[c] != INT32_MIN &&
                 boundary_correction[c] != k) ||
                (boundary_lineage[c] >= 0 &&
                 boundary_lineage[c] != parent_lineage)) {
                fprintf(stderr, "unwrap: parent boundary splits component %d "
                        "at vertex %zu (correction %d/%d lineage %d/%d); "
                        "refusing a median alignment\n", c, i,
                        boundary_correction[c], k,
                        boundary_lineage[c], parent_lineage);
                free(split_sibling);
                Arena_restore(arena, mark);
                return -1;
            }
            boundary_correction[c] = k;
            boundary_lineage[c] = parent_lineage;
            boundary_vertices++;
        }
        if (boundary_vertices == 0) {
            fprintf(stderr, "unwrap: projective boundary contains no finite "
                    "overlap samples\n");
            free(split_sibling);
            Arena_restore(arena, mark);
            return -1;
        }
        int32_t boundary_big = -1;
        size_t boundary_components = 0;
        for (int32_t c = 0; c < ncomp; c++) {
            if (boundary_correction[c] == INT32_MIN) continue;
            boundary_components++;
            if (boundary_big < 0 || csize[c] > csize[boundary_big])
                boundary_big = c;
        }
        if (boundary_big >= 0) big = boundary_big;
        boundary.correction = boundary_correction;
        boundary.lineage = boundary_lineage;
        boundary_ptr = &boundary;
        out->winding_boundary_vertices = boundary_vertices;
        fprintf(stderr, "  projective winding boundary: %zu vertices fix "
                "%zu/%d components; anchor component=%d\n",
                boundary_vertices, boundary_components, ncomp, big);
    }
    if (field_mode > 0 || (field_mode == 0 && ncomp > 1)) {
        size_t field_samples = 0, field_supported = 0;
        size_t field_clean = 0, field_invalid = 0;
        double field_epsilon = opts != NULL
                             ? opts->winding_field_epsilon : 0.0;
        double field_beta = opts != NULL ? opts->winding_field_beta : 0.0;
        if (unwrap_build_field_unary(
                arena, verts, nv, faces, nf, c1, c2, e1, e2, qturn,
                comp, ncomp, big, registration_pitch,
                field_epsilon, field_beta,
                &field_center, &field_sigma, &field_weight,
                &field_mean, &field_jump,
                &field_build_stats, &field_eval_stats,
                &field_samples, &field_supported,
                &field_clean, &field_invalid) == 0) {
            field_unary.center = field_center;
            field_unary.sigma = field_sigma;
            field_unary.weight = field_weight;
            field_unary_ptr = &field_unary;
            out->winding_field_used = 1;
            out->winding_field_backend = (int)field_eval_stats.backend_used;
            out->winding_field_samples = field_samples;
            out->winding_field_supported_components = field_supported;
            out->winding_field_clean_vertices = field_clean;
            out->winding_field_invalid_vertices = field_invalid;
            out->winding_field_clean_fraction =
                (double)field_clean / (double)nv;
        } else {
            fprintf(stderr,
                    "unwrap: winding field unavailable; MRF is relation-only\n");
        }
    }
    int32_t *turn_correction = NULL, *component_island = NULL;
    int32_t *component_continuation_island = NULL;
    float *component_winding_confidence = NULL;
    WindingRegisterStats wstats;
    /* per-vertex unit normals (area-weighted, sign irrelevant) for the
     * register's continuation tangency gate */
    float *vnormal = (float *)ARENA_CALLOC(arena, nv * 3, sizeof(float));
    {
        size_t f2 = 0, v2 = 0;
        for (f2 = 0; f2 < nf; f2++) {
            const int32_t *tri = faces + f2 * 3;
            const float *p0 = verts + (size_t)tri[0] * 3;
            const float *p1 = verts + (size_t)tri[1] * 3;
            const float *p2 = verts + (size_t)tri[2] * 3;
            float e1x = p1[0] - p0[0], e1y = p1[1] - p0[1], e1z = p1[2] - p0[2];
            float e2x = p2[0] - p0[0], e2y = p2[1] - p0[1], e2z = p2[2] - p0[2];
            float n0 = e1y * e2z - e1z * e2y;
            float n1 = e1z * e2x - e1x * e2z;
            float n2 = e1x * e2y - e1y * e2x;
            int k = 0;
            for (k = 0; k < 3; k++) {
                float *dst = vnormal + (size_t)tri[k] * 3;
                dst[0] += n0; dst[1] += n1; dst[2] += n2;
            }
        }
        for (v2 = 0; v2 < nv; v2++) {
            float *n = vnormal + v2 * 3;
            double len = sqrt((double)n[0] * n[0] + (double)n[1] * n[1] +
                              (double)n[2] * n[2]);
            if (len > 1e-12) {
                n[0] = (float)(n[0] / len);
                n[1] = (float)(n[1] / len);
                n[2] = (float)(n[2] / len);
            }
        }
    }
    if (WindingRegister_run_with_field(
            arena, verts, nv, t, r, theta, qturn, comp, ncomp, csize, big,
             tmin, registration_pitch, winding_sense,
             field_unary_ptr,
             boundary_ptr,
             /* The register now carries exact parent locks through its
              * conflict pass; parent presence is not a global disable. */
             opts == NULL || opts->winding_conflict_mode >= 0,
             split_sibling, n_split_sibling,
             opts != NULL ? opts->vertex_cube : NULL,
             vnormal,
             &turn_correction, &component_island,
            &component_continuation_island, &component_winding_confidence,
            &wstats) != 0) {
        fprintf(stderr, "unwrap: winding gauge registration failed\n");
        free(split_sibling);
        Arena_restore(arena, mark);
        return -1;
    }
    free(split_sibling);
    split_sibling = NULL;
    for (size_t i = 0; i < nv; i++)
        Phi[i] += (double)winding_sense * 2.0 * M_PI *
                  (double)turn_correction[comp[i]];

    /* --- 7a. Overlap-offset gauge sync across disconnected gauge islands.
     * Runs before the continuation join so joins compare SYNCED lifts. */
    if (!boundary_mode)
        (void)unwrap_sync_gauges(arena, nv, t, r, Phi, comp, ncomp,
                                 component_island, turn_correction,
                                 winding_sense, registration_pitch);

    /* --- 7a2. Junction-radius island gauging: spread interval-adjacent
     * islands to their radius-implied turns (see the pass's comment). */
    {
        size_t jg_shifted = 0, jg_abstained = 0;
        if (!boundary_mode)
            unwrap_gauge_islands_by_junction(
                arena, nv, t, r, Phi, comp, ncomp,
                component_continuation_island, turn_correction,
                winding_sense, registration_pitch,
                &jg_shifted, &jg_abstained);
        out->island_gauge_shifts = jg_shifted;
        out->island_gauge_abstains = jg_abstained;
    }

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
        if (!boundary_mode && nisl > 1) {
            double *ilo = (double *)ARENA_ALLOC(
                arena, (size_t)((size_t)nisl * sizeof(double)));
            double *ihi = (double *)ARENA_ALLOC(
                arena, (size_t)((size_t)nisl * sizeof(double)));
            double *tlo = (double *)ARENA_ALLOC(
                arena, (size_t)((size_t)nisl * sizeof(double)));
            double *thi = (double *)ARENA_ALLOC(
                arena, (size_t)((size_t)nisl * sizeof(double)));
            double *rhi_sum = (double *)ARENA_ALLOC(
                arena, (size_t)((size_t)nisl * sizeof(double)));
            double *rlo_sum = (double *)ARENA_ALLOC(
                arena, (size_t)((size_t)nisl * sizeof(double)));
            size_t *rhi_n = (size_t *)ARENA_ALLOC(
                arena, (size_t)((size_t)nisl * sizeof(size_t)));
            size_t *rlo_n = (size_t *)ARENA_ALLOC(
                arena, (size_t)((size_t)nisl * sizeof(size_t)));
            size_t *icount = (size_t *)ARENA_ALLOC(
                arena, (size_t)((size_t)nisl * sizeof(size_t)));
            int32_t *root = (int32_t *)ARENA_ALLOC(
                arena, (size_t)((size_t)nisl * sizeof(int32_t)));
            int32_t *order = (int32_t *)ARENA_ALLOC(
                arena, (size_t)((size_t)nisl * sizeof(int32_t)));
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
                        arena, (size_t)((size_t)nisl * sizeof(int32_t)));
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
    out->winding_repair_closers = wstats.repair_closers;
    out->winding_repair_conflicts_pre = wstats.repair_conflicts_pre;
    out->winding_repair_shifts = wstats.repair_shifts;
    out->winding_repair_capped_roots = wstats.repair_capped_roots;
    out->winding_anchor_span_pre_turns = wstats.anchor_span_pre_turns;
    out->winding_anchor_span_turns = wstats.anchor_span_turns;
    out->winding_mrf_rounds = wstats.mrf_rounds;
    out->winding_mrf_label_changes = wstats.mrf_label_changes;
    out->winding_mrf_abstained_sites = wstats.mrf_abstained_sites;
    out->winding_mrf_energy_before = wstats.mrf_energy_before;
    out->winding_mrf_energy_after = wstats.mrf_energy_after;
    out->winding_mrf_mean_confidence = wstats.mrf_mean_confidence;
    out->winding_mrf_field_calibrated_roots =
        wstats.mrf_field_calibrated_roots;
    out->winding_mrf_field_calibrated_sites =
        wstats.mrf_field_calibrated_sites;
    out->winding_mrf_field_calibration_r2 =
        wstats.mrf_field_calibration_r2;
    out->winding_mrf_conflict_rounds = wstats.mrf_conflict_rounds;
    out->winding_mrf_conflict_bins_before =
        wstats.mrf_conflict_bins_before;
    out->winding_mrf_conflict_bins_after =
        wstats.mrf_conflict_bins_after;
    out->winding_mrf_conflict_losing_claims =
        wstats.mrf_conflict_losing_claims;
    out->winding_mrf_conflict_exclusions =
        wstats.mrf_conflict_exclusions;
    out->winding_mrf_conflict_winner_locks =
        wstats.mrf_conflict_winner_locks;
    out->winding_mrf_conflict_label_changes =
        wstats.mrf_conflict_label_changes;
    out->winding_mrf_conflict_converged =
        wstats.mrf_conflict_converged;
    out->winding_boundary_components = wstats.boundary_components;
    out->winding_boundary_relation_cuts = wstats.boundary_relation_cuts;
    out->winding_boundary_lineage_cuts = wstats.boundary_lineage_cuts;
    out->winding_boundary_supported_relation_components =
        wstats.boundary_supported_relation_components;

    /* --- 8. Assemble UV (length-like, each axis shifted to start at 0). ---
     *
     * u is the winding converted to ARC LENGTH.  It used to be
     * (Phi - Phi_min) * r_ref with a single global median radius, which is
     * only correct at that one radius: on the 4x5x5, r_ref = 241.7 against
     * material from r = 19 to r = 457, so the umbilicus was oversampled by
     * ~5x and 2.2% of the emitted triangles carried under a quarter of the
     * sheet's median 3-D area per unit of UV area -- every one of them
     * squashed, none stretched, at a median radius of 45 vox.
     *
     * The fix is a single monotone reparameterization of the SAME lifted
     * phase, so the winding certificate is preserved exactly (every vertex at
     * one Phi still receives one u, and the ordering and turn structure are
     * untouched):
     *
     *     u(Phi) = integral from Phi_min to Phi of rbar(phi) dphi
     *
     * rbar is the empirical median radius at lifted phase phi.  At a fixed
     * LIFTED phase all material lies on one wrap -- that is what lifting
     * means -- so rbar is well defined, strictly positive, and du/dPhi is the
     * local arc length per radian.  A constant rbar reduces to the old
     * behaviour exactly. */
    double umin = INFINITY, umax = -INFINITY;
    for (size_t i = 0; i < nv; i++) {
        if (Phi[i] < umin) umin = Phi[i];
        if (Phi[i] > umax) umax = Phi[i];
    }
    double r_ref = ri[nv / 2].r;
    double arclen_span = 0.0;
    if (r_ref < 1e-6) r_ref = 1.0;   /* flat/degenerate -> 1 vox per radian */
    {
    /* Phase bins for the radius profile.  0.05 rad is ~3 degrees, far finer
     * than a wrap, and the count is capped so a whole-scroll run cannot
     * allocate an unbounded table. */
    size_t nbin = 0;
    double span = umax - umin;
    double bw = 0.0;
    double *cum = NULL, *rbar = NULL;
    size_t *bcount = NULL, *boff = NULL, *bcur = NULL;
    float *brad = NULL;
    if (!(span > 1e-9) || !isfinite(span)) {
        nbin = 0;
    } else {
        double want = span / UNWRAP_ARCLEN_BIN_RAD;
        nbin = want < 1.0 ? 1 : (want > (double)UNWRAP_ARCLEN_MAX_BINS
                                 ? (size_t)UNWRAP_ARCLEN_MAX_BINS
                                 : (size_t)want);
        bw = span / (double)nbin;
    }
    if (nbin > 0) {
        bcount = (size_t *)ARENA_CALLOC(arena, nbin + 1, sizeof *bcount);
        boff = (size_t *)ARENA_CALLOC(arena, nbin + 1, sizeof *boff);
        bcur = (size_t *)ARENA_CALLOC(arena, nbin + 1, sizeof *bcur);
        rbar = (double *)ARENA_ALLOC(arena, (size_t)(nbin * sizeof *rbar));
        cum = (double *)ARENA_ALLOC(arena, (size_t)((nbin + 1) * sizeof *cum));
        brad = (float *)ARENA_ALLOC(arena, (size_t)(nv * sizeof *brad));
    }
    if (nbin > 0 && bcount && boff && bcur && rbar && cum && brad) {
        size_t b = 0;
        /* two-pass bucket fill (count -> prefix sum -> place) */
        for (size_t i = 0; i < nv; i++) {
            double q = (Phi[i] - umin) / bw;
            long bi = (long)q;
            if (bi < 0) bi = 0;
            if ((size_t)bi >= nbin) bi = (long)nbin - 1;
            bcount[(size_t)bi]++;
        }
        for (b = 0; b < nbin; b++) boff[b + 1] = boff[b] + bcount[b];
        for (b = 0; b <= nbin; b++) bcur[b] = boff[b];
        for (size_t i = 0; i < nv; i++) {
            double q = (Phi[i] - umin) / bw;
            long bi = (long)q;
            if (bi < 0) bi = 0;
            if ((size_t)bi >= nbin) bi = (long)nbin - 1;
            brad[bcur[(size_t)bi]++] = (float)r[i];
        }
        for (b = 0; b < nbin; b++) {
            size_t lo = boff[b], hi = boff[b + 1];
            if (hi - lo >= UNWRAP_ARCLEN_MIN_BIN) {
                qsort(brad + lo, hi - lo, sizeof *brad, unwrap_cmp_float);
                rbar[b] = (double)brad[lo + (hi - lo) / 2];
            } else {
                rbar[b] = -1.0;      /* filled from a neighbour below */
            }
        }
        /* Sparse bins inherit the nearest measured profile rather than
         * inventing one; a run with no measured bin at all falls back to
         * r_ref, which is exactly the old behaviour. */
        {
            double last = -1.0;
            for (b = 0; b < nbin; b++) {
                if (rbar[b] > 0.0) last = rbar[b];
                else if (last > 0.0) rbar[b] = last;
            }
            last = -1.0;
            for (b = nbin; b-- > 0; ) {
                if (rbar[b] > 0.0) last = rbar[b];
                else if (last > 0.0) rbar[b] = last;
            }
            for (b = 0; b < nbin; b++)
                if (!(rbar[b] > UNWRAP_ARCLEN_MIN_R)) rbar[b] = r_ref;
        }
        cum[0] = 0.0;
        for (b = 0; b < nbin; b++) cum[b + 1] = cum[b] + rbar[b] * bw;
        for (size_t i = 0; i < nv; i++) {
            double q = (Phi[i] - umin) / bw;
            long bi = (long)q;
            double frac;
            if (bi < 0) bi = 0;
            if ((size_t)bi >= nbin) bi = (long)nbin - 1;
            frac = (Phi[i] - umin) - (double)bi * bw;
            if (frac < 0.0) frac = 0.0;
            if (frac > bw) frac = bw;
            out->uv[i * 2 + 0] = (float)(cum[(size_t)bi] + rbar[(size_t)bi] * frac);
        }
        arclen_span = cum[nbin];
        fprintf(stderr,
                "  unwrap u: arc-length map over %zu phase bin(s) "
                "(rbar %.1f..%.1f vox, r_ref %.1f); span %.1f vox "
                "(constant-r_ref would give %.1f)\n",
                nbin, rbar[0], rbar[nbin - 1], r_ref, arclen_span,
                span * r_ref);
    } else {
        for (size_t i = 0; i < nv; i++)
            out->uv[i * 2 + 0] = (float)((Phi[i] - umin) * r_ref);
        arclen_span = span * r_ref;
    }
    }
    for (size_t i = 0; i < nv; i++) {
        out->uv[i * 2 + 1] = (float)(t[i] - tmin);
        out->winding_index[i] = boundary_mode &&
            isfinite((double)opts->boundary_winding[i])
            ? opts->boundary_winding[i]
            : (float)((double)winding_sense * Phi[i] / (2.0 * M_PI));
        double local_confidence = 1.0;
        if (out->winding_field_used) {
            if (!isfinite(field_jump[i])) local_confidence = 0.0;
            else {
                double z = ((double)field_jump[i] - 1.0) / 0.35;
                local_confidence = exp(-0.5 * z * z);
            }
        }
        double component_confidence = component_winding_confidence != NULL
            ? component_winding_confidence[comp[i]] : 0.0;
        out->winding_confidence[i] = (float)(
            component_confidence * local_confidence);
    }
    if (keep_phi) {
        for (size_t i = 0; i < nv; i++) {
            out->phi[i] = boundary_mode &&
                isfinite((double)opts->boundary_winding[i])
                ? (float)((double)winding_sense * 2.0 * M_PI *
                          (double)opts->boundary_winding[i])
                : (float)Phi[i];
            out->island[i] = component_island[comp[i]];
            out->continuation_island[i] =
                component_continuation_island[comp[i]];
            out->mesh_component[i] = comp_true[i];
            out->field_winding[i] = out->winding_field_used
                                  ? field_mean[i] : NAN;
            out->field_jump[i] = out->winding_field_used
                               ? field_jump[i] : NAN;
        }
    }

    out->axis[0] = axis[0]; out->axis[1] = axis[1]; out->axis[2] = axis[2];
    out->centroid[0] = centroid[0];
    out->centroid[1] = centroid[1];
    out->centroid[2] = centroid[2];
    out->r_ref  = r_ref;
    out->u_span = arclen_span;
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
    float   *v = (float *)ARENA_ALLOC(arena, (size_t)(nvv * 3 * sizeof(float)));
    int32_t *f = (int32_t *)ARENA_ALLOC(arena, (size_t)(nff * 3 * sizeof(int32_t)));
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
    float   *v = (float *)ARENA_ALLOC(arena, (size_t)(nvv * 3 * sizeof(float)));
    int32_t *f = (int32_t *)ARENA_ALLOC(arena, (size_t)(nff * 3 * sizeof(int32_t)));
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

    /* (1g) Junction-radius island gauging: island B's lift packs it
     * phi-adjacent to island A while its radii sit 3 turns further out; the
     * junction snap must shift B by +3 turns.  A noisy half-integer case
     * must abstain. */
    {
        enum { GNB = 4, GPER = 9 };   /* axial bins x samples per bin end */
        size_t nvg = 2u * GNB * GPER * 2u;
        double *gt = (double *)ARENA_ALLOC(arena,
                                           (size_t)(nvg * sizeof(double)));
        double *gr = (double *)ARENA_ALLOC(arena,
                                           (size_t)(nvg * sizeof(double)));
        double *gphi = (double *)ARENA_ALLOC(arena,
                                             (size_t)(nvg * sizeof(double)));
        int32_t *gcomp = (int32_t *)ARENA_ALLOC(arena,
                                                (size_t)(nvg *
                                                       sizeof(int32_t)));
        int32_t gisl[2] = { 0, 1 };
        int32_t gtc[2] = { 0, 0 };
        size_t at = 0;
        size_t shifted = 0, abstained = 0;
        /* island 0: true w in [0,4] turns; island 1: true w in [7,10] but
         * lifted 2.8 turns low so its interval starts 0.2 turns after A */
        for (int isl = 0; isl < 2; isl++) {
            double wlo = isl == 0 ? 0.0 : 7.0 - 2.8;
            double whi = isl == 0 ? 4.0 : 10.0 - 2.8;
            double rlo = isl == 0 ? 50.0 : 50.0 + 9.5 * 7.0;
            double rhi = isl == 0 ? 50.0 + 9.5 * 4.0 : 50.0 + 9.5 * 10.0;
            for (int b = 0; b < GNB; b++) {
                for (int s = 0; s < GPER; s++) {
                    double f = (double)s / (GPER - 1) * 0.05;
                    gt[at] = 4.0 + 8.0 * b;
                    gphi[at] = (wlo + f) * 2.0 * M_PI;
                    gr[at] = rlo + f * 9.5;
                    gcomp[at] = isl;
                    at++;
                    gt[at] = 4.0 + 8.0 * b;
                    gphi[at] = (whi - f) * 2.0 * M_PI;
                    gr[at] = rhi - f * 9.5;
                    gcomp[at] = isl;
                    at++;
                }
            }
        }
        unwrap_gauge_islands_by_junction(
            arena, at, gt, gr, gphi, gcomp, 2, gisl, gtc, 1, 9.5,
            &shifted, &abstained);
        if (shifted != 1 ||
            (UNWRAP_JGAUGE_APPLY ? gtc[1] - gtc[0] != 3
                                 : gtc[1] != gtc[0])) {
            fprintf(stderr,
                    "[unwrap selftest] FAIL: junction gauge expected +3 "
                    "measurement (shifted=%zu tc=%d/%d apply=%d)\n",
                    shifted, gtc[0], gtc[1], UNWRAP_JGAUGE_APPLY);
            fails++;
        } else {
            fprintf(stderr,
                    "[unwrap selftest] junction island gauge OK "
                    "(+3 turns measured%s)\n",
                    UNWRAP_JGAUGE_APPLY ? ", applied" : ", apply disarmed");
        }
        /* abstain: radii midway between integer snaps (offset 2.5 turns) */
        at = 0;
        gtc[0] = gtc[1] = 0;
        for (int isl = 0; isl < 2; isl++) {
            double wlo = isl == 0 ? 0.0 : 7.0 - 2.8;
            double whi = isl == 0 ? 4.0 : 10.0 - 2.8;
            double rlo = isl == 0 ? 50.0 : 50.0 + 9.5 * (7.0 - 2.5 + 2.0);
            double rhi = isl == 0 ? 50.0 + 9.5 * 4.0
                                  : 50.0 + 9.5 * (10.0 - 2.5 + 2.0);
            for (int b = 0; b < GNB; b++) {
                for (int s = 0; s < GPER; s++) {
                    double f = (double)s / (GPER - 1) * 0.05;
                    gt[at] = 4.0 + 8.0 * b;
                    gphi[at] = (wlo + f) * 2.0 * M_PI;
                    gr[at] = rlo + f * 9.5;
                    gcomp[at] = isl;
                    at++;
                    gt[at] = 4.0 + 8.0 * b;
                    gphi[at] = (whi - f) * 2.0 * M_PI;
                    gr[at] = rhi - f * 9.5;
                    gcomp[at] = isl;
                    at++;
                }
            }
        }
        shifted = 0; abstained = 0;
        unwrap_gauge_islands_by_junction(
            arena, at, gt, gr, gphi, gcomp, 2, gisl, gtc, 1, 9.5,
            &shifted, &abstained);
        if (shifted != 0 || abstained == 0 || gtc[1] != gtc[0]) {
            fprintf(stderr,
                    "[unwrap selftest] FAIL: half-integer junction must "
                    "abstain (shifted=%zu abstained=%zu tc=%d/%d)\n",
                    shifted, abstained, gtc[0], gtc[1]);
            fails++;
        } else {
            fprintf(stderr,
                    "[unwrap selftest] junction island gauge abstains on "
                    "half-integer evidence OK\n");
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
            /* u is now an ARC-LENGTH map of phi, not phi*r_ref, so the
             * invariant to check is the one that actually matters: equal phi
             * must give equal u, u must be monotone in phi, and the local
             * slope must track the local radius instead of one global one. */
            double err = 0.0, slope_lo = 1e300, slope_hi = -1e300;
            (void)pmin;
            for (size_t i = 1; i < nvv && ok; i++) {
                size_t probe = i < 64 ? i : 64;
                for (size_t k = i - probe; k < i; k++) {
                    double dphi = (double)res.phi[i] - (double)res.phi[k];
                    double du = (double)res.uv[i*2+0] - (double)res.uv[k*2+0];
                    if (fabs(dphi) < 1e-6) {
                        if (fabs(du) > err) err = fabs(du);
                    } else if (dphi * du < -1e-9) {
                        ok = 0;
                        break;
                    } else if (fabs(dphi) > 0.2) {
                        double sl = du / dphi;
                        if (sl < slope_lo) slope_lo = sl;
                        if (sl > slope_hi) slope_hi = sl;
                    }
                }
            }
            ok = ok && (err < 1e-2);
            if (!ok)
                fprintf(stderr, "[unwrap selftest] KEEP_PHI FAIL "
                        "monotone/equal-phi spread=%.6f\n", err);
            else
                fprintf(stderr, "[unwrap selftest] keep_phi u monotone in phi; "
                        "slope %.1f..%.1f vox/rad (r_ref %.1f)\n",
                        slope_lo, slope_hi, res.r_ref);
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
            arena, (size_t)(nvv * 3 * sizeof(float)));
        int32_t *f = (int32_t *)ARENA_ALLOC(
            arena, (size_t)((size_t)(JNPHI - 1) * (JNH - 1) * 4 * 3 *
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
