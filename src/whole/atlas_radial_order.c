#include "atlas_radial_order.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef ARO_2PI
#define ARO_2PI 6.283185307179586476925286766559
#endif

void AtlasRadialOrderOptions_default(AtlasRadialOrderOptions *opts)
{
    if (opts == NULL) return;
    opts->arc_window = 6.0;
    opts->separation_min_ratio = 0.6;
    opts->turn_tolerance = 0.30;
    opts->core_radius_pitches = 2.0;
    opts->gap_min = 8.0;
    opts->lambda_spacing = 0.10;
    opts->max_turns = 4;
}

/* ---------------------------------------------------------------- geometry */

typedef struct {
    double axis[3];
    double origin[3];
    double e0[3];
    double e1[3];
} AroFrame;

static double aro_dot(const double a[3], const double b[3])
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static void aro_cross(const double a[3], const double b[3], double out[3])
{
    out[0] = a[1] * b[2] - a[2] * b[1];
    out[1] = a[2] * b[0] - a[0] * b[2];
    out[2] = a[0] * b[1] - a[1] * b[0];
}

static int aro_normalize(double v[3])
{
    double n = sqrt(aro_dot(v, v));
    if (!(n > 1e-15)) return -1;
    v[0] /= n;
    v[1] /= n;
    v[2] /= n;
    return 0;
}

/* Right-handed frame with e0,e1 spanning the plane perpendicular to the axis. */
static int aro_frame_build(const ScaffoldCalib *cal, AroFrame *f)
{
    f->axis[0] = (double)cal->axis_dir[0];
    f->axis[1] = (double)cal->axis_dir[1];
    f->axis[2] = (double)cal->axis_dir[2];
    if (aro_normalize(f->axis) != 0) return -1;
    f->origin[0] = (double)cal->axis_point[0];
    f->origin[1] = (double)cal->axis_point[1];
    f->origin[2] = (double)cal->axis_point[2];

    double seed[3] = {0.0, 0.0, 0.0};
    double ax = fabs(f->axis[0]), ay = fabs(f->axis[1]), az = fabs(f->axis[2]);
    if (ax <= ay && ax <= az) seed[0] = 1.0;
    else if (ay <= az) seed[1] = 1.0;
    else seed[2] = 1.0;

    double proj = aro_dot(seed, f->axis);
    for (int i = 0; i < 3; i++) f->e0[i] = seed[i] - proj * f->axis[i];
    if (aro_normalize(f->e0) != 0) return -1;
    aro_cross(f->axis, f->e0, f->e1);
    return aro_normalize(f->e1);
}

static void aro_polar(const AroFrame *f, const double p[3],
                      double *out_r, double *out_theta)
{
    double d[3] = {p[0] - f->origin[0], p[1] - f->origin[1],
                   p[2] - f->origin[2]};
    double a = aro_dot(d, f->e0);
    double b = aro_dot(d, f->e1);
    *out_r = sqrt(a * a + b * b);
    double t = atan2(b, a);
    if (t < 0.0) t += ARO_2PI;
    *out_theta = t;
}

/* ------------------------------------------------------------------ sorting */

typedef struct {
    int32_t plane;
    double  theta;
    double  radius;
    int32_t sample;
} AroSlot;

/* Sorted by plane, then angle: the sweep needs contiguous per-plane runs and
 * a monotone angle inside each so the window can be found by bisection. */
static int aro_slot_cmp(const void *va, const void *vb)
{
    const AroSlot *a = (const AroSlot *)va;
    const AroSlot *b = (const AroSlot *)vb;
    if (a->plane < b->plane) return -1;
    if (a->plane > b->plane) return 1;
    if (a->theta < b->theta) return -1;
    if (a->theta > b->theta) return 1;
    if (a->sample < b->sample) return -1;
    if (a->sample > b->sample) return 1;
    return 0;
}

/* First index in [lo,hi) whose angle is >= target. */
static size_t aro_lower_bound(const AroSlot *slot, size_t lo, size_t hi,
                              double target)
{
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (slot[mid].theta < target) lo = mid + 1;
        else hi = mid;
    }
    return lo;
}

static int aro_double_cmp(const void *va, const void *vb)
{
    double a = *(const double *)va;
    double b = *(const double *)vb;
    if (a < b) return -1;
    if (a > b) return 1;
    return 0;
}

static double aro_wrap_angle(double angle)
{
    while (angle > 0.5 * ARO_2PI) angle -= ARO_2PI;
    while (angle < -0.5 * ARO_2PI) angle += ARO_2PI;
    return angle;
}

static int aro_delam_contains(const uint64_t *keys, size_t nkeys, uint64_t key)
{
    size_t lo = 0, hi = nkeys;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (keys[mid] < key) lo = mid + 1;
        else hi = mid;
    }
    return lo < nkeys && keys[lo] == key;
}

static uint64_t aro_pair_key(int32_t a, int32_t b)
{
    uint32_t lo = (uint32_t)(a < b ? a : b);
    uint32_t hi = (uint32_t)(a < b ? b : a);
    return ((uint64_t)lo << 32) | (uint64_t)hi;
}

/* -------------------------------------------------------------------- build */

int AtlasRadialOrder_build(Arena_T arena,
                           const AtlasStripSample *samples,
                           const int32_t *sample_plane,
                           const int32_t *sample_component,
                           const double *sample_phi,
                           const int32_t *phase_domain,
                           size_t nsamples,
                           const double *base_u,
                           const ScaffoldCalib *cal,
                           const uint64_t *delam_keys,
                           size_t ndelam_keys,
                           const AtlasRadialOrderOptions *input_opts,
                           AtlasRadialOrderSet *out)
{
    if (samples == NULL || sample_plane == NULL || base_u == NULL ||
        sample_phi == NULL || cal == NULL || out == NULL)
        return -1;
    if (nsamples == 0 || nsamples > (size_t)INT32_MAX) return -1;

    AtlasRadialOrderOptions opts;
    if (input_opts != NULL) opts = *input_opts;
    else AtlasRadialOrderOptions_default(&opts);
    if (!(opts.arc_window > 0.0) || opts.max_turns < 1) return -1;
    if (!(opts.gap_min >= 0.0) || !(opts.separation_min_ratio > 0.0))
        return -1;

    double pitch = fabs(cal->pitch);
    if (!(pitch > 1e-6)) return -1;

    AroFrame frame;
    if (aro_frame_build(cal, &frame) != 0) return -1;

    memset(out, 0, sizeof *out);

    /* Everything here is either an output or small and one-shot, so the
     * allocations stay on the caller's arena without a scratch mark: a mark
     * taken before the outputs would release them on restore. */
    AroSlot *slot = (AroSlot *)ARENA_ALLOC(arena, nsamples * sizeof *slot);
    for (size_t i = 0; i < nsamples; i++) {
        double r = 0.0, theta = 0.0;
        aro_polar(&frame, samples[i].p, &r, &theta);
        slot[i].plane = sample_plane[i];
        slot[i].theta = theta;
        slot[i].radius = r;
        slot[i].sample = (int32_t)i;
    }
    qsort(slot, nsamples, sizeof *slot, aro_slot_cmp);

    double core_radius = opts.core_radius_pitches * pitch;
    double separation_min = opts.separation_min_ratio * pitch;
    double max_gap = ((double)opts.max_turns + opts.turn_tolerance) * pitch;

    /* One candidate pair per sample at most, so the outputs can be sized up
     * front and the sweep runs once. */
    AtlasRadialOrderPair *pairs = (AtlasRadialOrderPair *)ARENA_ALLOC(
        arena, nsamples * sizeof *pairs);
    MonotoneQpBound *bounds = (MonotoneQpBound *)ARENA_ALLOC(
        arena, nsamples * sizeof *bounds);
    MonotoneQpRow *rows = (MonotoneQpRow *)ARENA_ALLOC(
        arena, nsamples * sizeof *rows);
    MonotoneQpCoeff *coeff = (MonotoneQpCoeff *)ARENA_ALLOC(
        arena, 2 * nsamples * sizeof *coeff);
    double *gaps = (double *)ARENA_ALLOC(arena, nsamples * sizeof *gaps);

    size_t np = 0;
    size_t rej_core = 0, rej_none = 0, rej_delam = 0, rej_wide = 0;
    size_t rej_same = 0, rej_ambig = 0;
    size_t visits = 0, planes = 0;
    size_t agree_radius = 0, disagree_radius = 0;
    size_t coincident = 0, satisfied = 0, reversed = 0;

    /* Sweep each plane: for every sample, the nearest sample outside it inside
     * the arc window is its next wrap out. */
    size_t pstart = 0;
    while (pstart < nsamples) {
        size_t pend = pstart + 1;
        while (pend < nsamples && slot[pend].plane == slot[pstart].plane) pend++;
        planes++;

        for (size_t i = pstart; i < pend; i++) {
            const AroSlot *a = &slot[i];
            if (a->radius < core_radius) { rej_core++; continue; }
            double half = opts.arc_window / a->radius;
            if (half > ARO_2PI * 0.5) half = ARO_2PI * 0.5;

            int32_t best = -1;
            double best_r = 0.0;
            double want = a->radius + separation_min;
            double limit = a->radius + max_gap;

            /* The window may straddle theta = 0; scan it as up to two runs. */
            double lo_t = a->theta - half, hi_t = a->theta + half;
            size_t runs[2][2];
            int nruns = 0;
            if (lo_t < 0.0) {
                runs[nruns][0] = aro_lower_bound(slot, pstart, pend,
                                                 lo_t + ARO_2PI);
                runs[nruns][1] = pend;
                nruns++;
                runs[nruns][0] = pstart;
                runs[nruns][1] = aro_lower_bound(slot, pstart, pend, hi_t);
                nruns++;
            } else if (hi_t > ARO_2PI) {
                runs[nruns][0] = aro_lower_bound(slot, pstart, pend, lo_t);
                runs[nruns][1] = pend;
                nruns++;
                runs[nruns][0] = pstart;
                runs[nruns][1] = aro_lower_bound(slot, pstart, pend,
                                                 hi_t - ARO_2PI);
                nruns++;
            } else {
                runs[nruns][0] = aro_lower_bound(slot, pstart, pend, lo_t);
                runs[nruns][1] = aro_lower_bound(slot, pstart, pend, hi_t);
                nruns++;
            }

            /*
             * Take the innermost candidate certified as a real wrap out.
             * Within one phase domain, phi differences are valid and retain
             * the important turn-zero delamination test.  Across domains the
             * independently rounded integer gauges do not cancel, so use the
             * local geometric winding differential instead.
             *
             * Testing only the nearest candidate would still be wrong: a
             * delaminated ply may sit between a wrap and its neighbour.
             */
            double best_turns = 0.0;
            int32_t sense = cal->sense != 0 ? cal->sense : 1;
            for (int r = 0; r < nruns; r++) {
                for (size_t j = runs[r][0]; j < runs[r][1]; j++) {
                    visits++;
                    const AroSlot *b = &slot[j];
                    if (b->radius < want || b->radius > limit) continue;
                    if (best >= 0 && b->radius >= best_r) continue;

                    int same_domain =
                        phase_domain == NULL ||
                        phase_domain[b->sample] == phase_domain[a->sample];
                    double turns_real = 0.0;
                    if (same_domain) {
                        double dphi = (sample_phi[b->sample] -
                                       sample_phi[a->sample]) * (double)sense;
                        turns_real = dphi / ARO_2PI;
                    } else {
                        double dtheta =
                            aro_wrap_angle(b->theta - a->theta);
                        double dwind =
                            (double)sense * (b->radius - a->radius) / pitch -
                            dtheta / ARO_2PI;
                        turns_real = fabs(dwind);
                    }
                    double m = floor(turns_real + 0.5);
                    if (fabs(turns_real - m) > opts.turn_tolerance) {
                        rej_ambig++;
                        continue;
                    }
                    if (m < 1.0) { rej_same++; continue; }
                    if (m > (double)opts.max_turns) { rej_wide++; continue; }
                    if (sample_component != NULL && delam_keys != NULL &&
                        ndelam_keys > 0) {
                        int32_t ca = sample_component[a->sample];
                        int32_t cb = sample_component[b->sample];
                        if (ca >= 0 && cb >= 0 && ca != cb &&
                            aro_delam_contains(delam_keys, ndelam_keys,
                                               aro_pair_key(ca, cb))) {
                            rej_delam++;
                            continue;
                        }
                    }
                    best = b->sample;
                    best_r = b->radius;
                    best_turns = m;
                }
            }
            if (best < 0) { rej_none++; continue; }

            double gap = best_r - a->radius;
            double m = best_turns;
            double radius_turns = floor(gap / pitch + 0.5);
            if (radius_turns == m) agree_radius++;
            else disagree_radius++;

            /*
             * One turn costs one circumference at the pair midpoint.  This is
             * the exact Archimedean-spiral difference and, unlike evaluating
             * CubeReg_deltaU at an absolute phi, is immune to per-cube integer
             * gauge errors.
             */
            double radius_mid = 0.5 * (a->radius + best_r);
            double target = (double)sense * m * ARO_2PI * radius_mid;

            double du = base_u[best] - base_u[a->sample];
            double mag = fabs(target);
            if (mag > 0.0) {
                if (fabs(du) < 0.5 * mag) coincident++;
                else if (du * target < 0.0) reversed++;
                else if (fabs(du - target) < 0.5 * mag) satisfied++;
            }

            AtlasRadialOrderPair *pr = &pairs[np];
            pr->inner = a->sample;
            pr->outer = best;
            pr->plane = a->plane;
            pr->turns = (int32_t)m;
            pr->radius_inner = a->radius;
            pr->radius_gap = gap;
            pr->phi_inner = sample_phi[a->sample];
            pr->spacing_target = target;
            pr->weight = fmin(samples[a->sample].confidence,
                              samples[best].confidence);
            if (!(pr->weight > 0.0)) pr->weight = 1.0;
            gaps[np] = gap;
            np++;
        }
        pstart = pend;
    }

    /* Each pair carries its own signed target, so the bound direction is
     * per-pair too: whichever end the spiral says should hold the larger u. */
    size_t npairs = np, nr = 0, nk = 0;
    for (size_t i = 0; i < npairs; i++) {
        const AtlasRadialOrderPair *pr = &pairs[i];
        int forward = pr->spacing_target >= 0.0;
        int32_t lo = forward ? pr->inner : pr->outer;
        int32_t hi = forward ? pr->outer : pr->inner;
        bounds[i].lo = lo;
        bounds[i].hi = hi;
        bounds[i].lower = opts.gap_min;
        bounds[i].stroke = -1;             /* marks a ray bound */
        bounds[i].edge = (int32_t)i;

        MonotoneQpRow *row = &rows[nr++];
        row->first = nk;
        row->count = 2;
        row->target = fabs(pr->spacing_target);
        row->weight = opts.lambda_spacing * pr->weight;
        row->kind = ATLAS_RADIAL_ROW_SPACING;
        row->owner = (int32_t)i;
        coeff[nk].var = hi;
        coeff[nk].value = 1.0;
        nk++;
        coeff[nk].var = lo;
        coeff[nk].value = -1.0;
        nk++;
    }
    size_t nspacing = nr;

    double median_gap = 0.0;
    if (npairs > 0) {
        qsort(gaps, npairs, sizeof *gaps, aro_double_cmp);
        median_gap = gaps[npairs / 2];
    }

    out->pairs = pairs;
    out->npairs = npairs;
    out->bounds = bounds;
    out->nbounds = npairs;
    out->rows = rows;
    out->nrows = nspacing;
    out->coeff = coeff;
    out->ncoeff = nk;

    out->stats.samples_considered = nsamples;
    out->stats.rejected_core = rej_core;
    out->stats.rejected_no_neighbour = rej_none;
    out->stats.rejected_delamination = rej_delam;
    out->stats.rejected_same_wrap = rej_same;
    out->stats.rejected_ambiguous = rej_ambig;
    out->stats.rejected_wide = rej_wide;
    out->stats.admitted = npairs;
    out->stats.planes_used = planes;
    out->stats.window_visits = visits;
    out->stats.turn_agree_radius = agree_radius;
    out->stats.turn_disagree_radius = disagree_radius;
    out->stats.coincident_pairs = coincident;
    out->stats.satisfied_pairs = satisfied;
    out->stats.reversed_pairs = reversed;
    out->stats.median_radius_gap = median_gap;
    return 0;
}

/* ---------------------------------------------------------------- bootstrap */

int AtlasRadialOrder_bootstrap(Arena_T arena,
                               const MonotoneQpBound *bounds, size_t nbounds,
                               const MonotoneQpAnchor *anchors, size_t nanchors,
                               size_t nvar,
                               double tolerance,
                               double *x,
                               uint8_t *disabled,
                               AtlasRadialOrderBootstrapStats *stats)
{
    if (x == NULL || nvar == 0 || nvar > (size_t)INT32_MAX) return -1;
    if (nbounds > 0 && bounds == NULL) return -1;
    if (!(tolerance >= 0.0)) tolerance = 0.0;

    AtlasRadialOrderBootstrapStats st;
    memset(&st, 0, sizeof st);

    Arena_Mark mark = Arena_save(arena);
    uint8_t *anchored = (uint8_t *)ARENA_ALLOC(arena, nvar * sizeof *anchored);
    memset(anchored, 0, nvar * sizeof *anchored);
    for (size_t i = 0; i < nanchors; i++) {
        int32_t v = anchors[i].var;
        if (v >= 0 && (size_t)v < nvar) anchored[v] = 1;
    }

    uint8_t *dead = (uint8_t *)ARENA_ALLOC(
        arena, (nbounds ? nbounds : 1) * sizeof *dead);
    memset(dead, 0, (nbounds ? nbounds : 1) * sizeof *dead);

    /* CSR: outgoing bounds per lo variable. */
    size_t *head = (size_t *)ARENA_ALLOC(arena, (nvar + 1) * sizeof *head);
    memset(head, 0, (nvar + 1) * sizeof *head);
    for (size_t i = 0; i < nbounds; i++) {
        int32_t lo = bounds[i].lo, hi = bounds[i].hi;
        if (lo < 0 || (size_t)lo >= nvar || hi < 0 || (size_t)hi >= nvar) {
            Arena_restore(arena, mark);
            return -1;
        }
        head[(size_t)lo + 1]++;
    }
    for (size_t v = 0; v < nvar; v++) head[v + 1] += head[v];
    int32_t *adj = (int32_t *)ARENA_ALLOC(
        arena, (nbounds ? nbounds : 1) * sizeof *adj);
    size_t *cursor = (size_t *)ARENA_ALLOC(arena, nvar * sizeof *cursor);
    for (size_t v = 0; v < nvar; v++) cursor[v] = head[v];
    for (size_t i = 0; i < nbounds; i++)
        adj[cursor[(size_t)bounds[i].lo]++] = (int32_t)i;

    /* SPFA with a per-variable relaxation counter for cycle detection. */
    int32_t *queue = (int32_t *)ARENA_ALLOC(arena, nvar * sizeof *queue);
    uint8_t *queued = (uint8_t *)ARENA_ALLOC(arena, nvar * sizeof *queued);
    int32_t *raiser = (int32_t *)ARENA_ALLOC(arena, nvar * sizeof *raiser);
    size_t *touched = (size_t *)ARENA_ALLOC(arena, nvar * sizeof *touched);
    double *x0 = (double *)ARENA_ALLOC(arena, nvar * sizeof *x0);
    memset(queued, 0, nvar * sizeof *queued);
    memset(touched, 0, nvar * sizeof *touched);
    for (size_t v = 0; v < nvar; v++) {
        queue[v] = (int32_t)v;
        queued[v] = 1;
        raiser[v] = -1;
        x0[v] = x[v];
    }

    /* Circular buffer holding every variable exactly once: a push happens only
     * when queued[] is clear, so it can never overflow its nvar capacity. */
    size_t cap = nvar, qhead = 0, qtail = 0, qcount = nvar;

    /* The cycle cap must exceed any legitimate chain depth.  Bounds only ever
     * link samples within one slice plane (stroke edges and rays are both
     * in-plane), so the longest chain is a plane's sample count -- a few
     * thousand.  Beyond that, the graph has a positive cycle. */
    size_t touch_cap = nvar < 64 ? 64 : nvar;
    if (touch_cap > 4096) touch_cap = 4096;

    while (qcount > 0) {
        int32_t v = queue[qhead];
        qhead = (qhead + 1) % cap;
        qcount--;
        queued[v] = 0;
        for (size_t e = head[v]; e < head[v + 1]; e++) {
            int32_t bi = adj[e];
            if (dead[bi]) continue;
            int32_t hi = bounds[bi].hi;
            double need = x[v] + bounds[bi].lower;
            if (need <= x[hi] + tolerance) continue;
            if (anchored[hi]) {
                dead[bi] = 1;
                st.anchor_blocked++;
                continue;
            }
            if (++touched[hi] > touch_cap) {
                /* Positive cycle.  Break it on a ray bound, never a stroke. */
                int32_t victim = bounds[bi].stroke < 0 ? bi : raiser[hi];
                if (victim < 0 || bounds[victim].stroke >= 0) victim = bi;
                dead[victim] = 1;
                st.cycle_broken++;
                touched[hi] = 0;
                if (victim == bi) continue;
            }
            x[hi] = need;
            raiser[hi] = bi;
            st.relaxations++;
            if (!queued[hi]) {
                queue[qtail] = hi;
                qtail = (qtail + 1) % cap;
                qcount++;
                queued[hi] = 1;
            }
        }
        st.rounds++;
        if (st.rounds > (size_t)16 * (nvar + nbounds) + 1024) break;
    }

    for (size_t v = 0; v < nvar; v++) {
        double d = fabs(x[v] - x0[v]);
        if (d > 0.0) {
            st.moved_variables++;
            st.total_shift += d;
            if (d > st.max_shift) st.max_shift = d;
        }
    }

    int rc = 0;
    for (size_t i = 0; i < nbounds; i++) {
        if (dead[i]) continue;
        if (x[bounds[i].hi] - x[bounds[i].lo] <
            bounds[i].lower - tolerance - 1e-9) {
            rc = -1;
            break;
        }
    }
    if (disabled != NULL) memcpy(disabled, dead, nbounds * sizeof *dead);
    if (stats != NULL) *stats = st;
    Arena_restore(arena, mark);
    return rc;
}

/* ----------------------------------------------------------------- selftest */

#define ARO_CHECK(cond, msg)                                                  \
    do {                                                                      \
        if (!(cond)) {                                                        \
            fprintf(stderr, "[atlas_radial_order] FAIL %s (%s:%d)\n",         \
                    (msg), __FILE__, __LINE__);                               \
            Arena_dispose(&arena);                                            \
            return -1;                                                        \
        }                                                                     \
    } while (0)

/*
 * Fixture: a 4-wrap Archimedean spiral seen on one slice plane, sampled at a
 * handful of angles.  Wrap w sits at radius r0 + w*pitch, so the ground-truth
 * radial order along every ray is exactly the wrap order.  Case (b) adds a
 * delaminated ply beside wrap 1 (a second sheet at +1.5 vox, far below the
 * cross-wrap threshold) and case (c) removes wrap 2 so the surviving pair must
 * be counted as a two-turn jump.
 */
int AtlasRadialOrder_selftest(void)
{
    Arena_T arena = Arena_new();
    const double pitch = 12.0;
    const double r0 = 60.0;
    const int nwrap = 4;
    const int nang = 8;

    /* Left-handed scroll, as PHerc0139 is: sense = -1, so phi grows more
     * negative outward and u DECREASES outward.  Choosing spiral_a so that
     * CubeReg_deltaU comes out at exactly minus one local circumference keeps
     * the fixture honest about the negative-target path the real data takes. */
    ScaffoldCalib cal;
    memset(&cal, 0, sizeof cal);
    cal.spiral_a = r0;   /* makes one turn cost exactly one circumference */
    cal.spiral_b = -pitch;
    cal.pitch = pitch;
    cal.sense = -1;
    cal.axis_dir[0] = 1.0f;

    AtlasRadialOrderOptions opts;
    AtlasRadialOrderOptions_default(&opts);
    opts.core_radius_pitches = 1.0;   /* r0 = 5 pitches, well outside */
    /* The fixture samples 8 angles per wrap, so the window must be wide enough
     * to reach a neighbour: an eighth of the innermost circumference. */
    opts.arc_window = ARO_2PI * r0 / (double)nang * 0.75;

    size_t cap = (size_t)(nwrap * nang) + 8;
    AtlasStripSample *samples =
        (AtlasStripSample *)ARENA_ALLOC(arena, cap * sizeof *samples);
    int32_t *plane = (int32_t *)ARENA_ALLOC(arena, cap * sizeof *plane);
    int32_t *comp = (int32_t *)ARENA_ALLOC(arena, cap * sizeof *comp);
    int32_t *phase_domain =
        (int32_t *)ARENA_ALLOC(arena, cap * sizeof *phase_domain);
    double *phi = (double *)ARENA_ALLOC(arena, cap * sizeof *phi);
    double *phi_offset = (double *)ARENA_ALLOC(arena, cap * sizeof *phi_offset);
    double *base_u = (double *)ARENA_ALLOC(arena, cap * sizeof *base_u);

#define ARO_PLACE(idx, wrapf, angle, radius, component)                       \
    do {                                                                      \
        double th = ARO_2PI * (double)(angle) / (double)nang;                 \
        double ph = -ARO_2PI * (wrapf);                                       \
        memset(&samples[idx], 0, sizeof samples[idx]);                        \
        samples[idx].p[0] = 0.0;             /* one slice plane */            \
        samples[idx].p[1] = (radius) * cos(th);                               \
        samples[idx].p[2] = (radius) * sin(th);                               \
        samples[idx].confidence = 1.0;                                        \
        plane[idx] = 0;                                                       \
        comp[idx] = (component);                                              \
        phi[idx] = ph;                                                        \
        base_u[idx] = cal.spiral_a * ph +                                     \
                      cal.spiral_b * ph * ph / (2.0 * ARO_2PI);               \
    } while (0)

    size_t k = 0;
    for (int w = 0; w < nwrap; w++)
        for (int a = 0; a < nang; a++) {
            double f = (double)w + (double)a / (double)nang;
            ARO_PLACE(k, f, a, r0 + pitch * f, w);
            k++;
        }

    /* (a) clean spiral: every sample finds the next wrap out, one turn away,
     *     and the target is negative because u decreases outward here. */
    AtlasRadialOrderSet set;
    ARO_CHECK(AtlasRadialOrder_build(arena, samples, plane, comp, phi, NULL, k,
                                     base_u, &cal, NULL, 0, &opts, &set) == 0,
              "build clean spiral");
    ARO_CHECK(set.npairs == (size_t)((nwrap - 1) * nang),
              "one pair per sample that has a wrap outside it");
    ARO_CHECK(set.stats.rejected_no_neighbour == (size_t)nang,
              "the outermost wrap has nothing beyond it");
    ARO_CHECK(set.stats.turn_disagree_radius == 0,
              "radius and phi agree on a clean spiral");
    ARO_CHECK(set.stats.satisfied_pairs == set.npairs,
              "the true field already satisfies its own evidence");
    ARO_CHECK(set.stats.coincident_pairs == 0, "nothing is coincident");
    for (size_t i = 0; i < set.npairs; i++) {
        const AtlasRadialOrderPair *p = &set.pairs[i];
        ARO_CHECK(p->turns == 1, "adjacent wraps are one turn apart");
        ARO_CHECK(p->spacing_target < 0.0, "u decreases outward on this scroll");
        double circ = ARO_2PI * (p->radius_inner + 0.5 * p->radius_gap);
        ARO_CHECK(fabs(fabs(p->spacing_target) - circ) < 1e-6,
                  "one turn costs one local circumference");
        ARO_CHECK(set.bounds[i].lo == p->outer && set.bounds[i].hi == p->inner,
                  "the bound points the way the spiral says");
        ARO_CHECK(set.bounds[i].stroke == -1, "ray bounds are marked");
    }

    /* Independent cube registration may add any whole number of turns to phi.
     * Cross-domain classification and its physical spacing target must be
     * exactly invariant to those arbitrary gauges. */
    for (size_t i = 0; i < k; i++) {
        phase_domain[i] = comp[i];
        phi_offset[i] = phi[i] + ARO_2PI * (double)(3 * comp[i] - 4);
    }
    AtlasRadialOrderSet offset_set;
    ARO_CHECK(AtlasRadialOrder_build(arena, samples, plane, comp, phi_offset,
                                     phase_domain, k, base_u, &cal, NULL, 0,
                                     &opts, &offset_set) == 0,
              "build with independent whole-turn phase gauges");
    ARO_CHECK(offset_set.npairs == set.npairs,
              "phase gauges do not change radial pair admission");
    for (size_t i = 0; i < set.npairs; i++) {
        ARO_CHECK(offset_set.pairs[i].inner == set.pairs[i].inner &&
                  offset_set.pairs[i].outer == set.pairs[i].outer &&
                  offset_set.pairs[i].turns == set.pairs[i].turns,
                  "phase gauges do not change inferred wrap count");
        ARO_CHECK(fabs(offset_set.pairs[i].spacing_target -
                       set.pairs[i].spacing_target) < 1.0e-9,
                  "phase gauges do not change circumference targets");
    }

    /* (b) a delaminated ply shares its parent's phi, so it is the same wrap
     *     and must never be ordered against it -- however far apart the axis
     *     wander makes them look radially. */
    size_t nd = k;
    for (int a = 0; a < nang; a++) {
        double f = 1.0 + (double)a / (double)nang;
        /* Deliberately past the cross-wrap separation, and nearer than the
         * real next wrap: radius alone would call this a wrap, and it also
         * stands between wrap 1 and wrap 2. */
        ARO_PLACE(nd, f, a, r0 + pitch * f + pitch * 0.8, nwrap);
        nd++;
    }
    AtlasRadialOrderSet delam;
    ARO_CHECK(AtlasRadialOrder_build(arena, samples, plane, comp, phi, NULL, nd,
                                     base_u, &cal, NULL, 0, &opts, &delam) == 0,
              "build with delamination");
    ARO_CHECK(delam.stats.rejected_same_wrap >= (size_t)nang,
              "the ply was recognized as the same wrap");
    size_t wrap1_to_wrap2 = 0;
    for (size_t i = 0; i < delam.npairs; i++) {
        int32_t ci = comp[delam.pairs[i].inner];
        int32_t co = comp[delam.pairs[i].outer];
        ARO_CHECK(!(ci == 1 && co == nwrap) && !(ci == nwrap && co == 1),
                  "the ply is never ordered against its own parent wrap");
        if (ci == 1 && co == 2) wrap1_to_wrap2++;
    }
    /* The ply sits between wrap 1 and wrap 2, so this is the real test: the
     * search had to look PAST it rather than give up on the sample. */
    ARO_CHECK(wrap1_to_wrap2 == (size_t)nang,
              "the constraint behind the ply survived");

    /* (c) a missing wrap must be counted as a two-turn jump, not rejected. */
    size_t ng = 0;
    for (int w = 0; w < nwrap; w++) {
        if (w == 2) continue;                    /* wrap 2 is missing */
        for (int a = 0; a < nang; a++) {
            double f = (double)w + (double)a / (double)nang;
            ARO_PLACE(ng, f, a, r0 + pitch * f, w);
            ng++;
        }
    }
    AtlasRadialOrderSet gap;
    ARO_CHECK(AtlasRadialOrder_build(arena, samples, plane, comp, phi, NULL, ng,
                                     base_u, &cal, NULL, 0, &opts, &gap) == 0,
              "build with a missing wrap");
    size_t two_turn = 0;
    for (size_t i = 0; i < gap.npairs; i++)
        if (gap.pairs[i].turns == 2) two_turn++;
    ARO_CHECK(two_turn == (size_t)nang, "the skipped wrap is counted");
#undef ARO_PLACE

    /* (d) bootstrap repairs an infeasible start and reports the shift. */
    {
        double x[4] = {0.0, 5.0, 1.0, 9.0};
        MonotoneQpBound b[3];
        b[0].lo = 0; b[0].hi = 1; b[0].lower = 10.0; b[0].stroke = -1; b[0].edge = 0;
        b[1].lo = 1; b[1].hi = 2; b[1].lower = 10.0; b[1].stroke = -1; b[1].edge = 1;
        b[2].lo = 2; b[2].hi = 3; b[2].lower = 10.0; b[2].stroke = -1; b[2].edge = 2;
        AtlasRadialOrderBootstrapStats bs;
        ARO_CHECK(AtlasRadialOrder_bootstrap(arena, b, 3, NULL, 0, 4, 0.0,
                                             x, NULL, &bs) == 0,
                  "bootstrap chain");
        ARO_CHECK(x[1] >= x[0] + 10.0 - 1e-12 &&
                  x[2] >= x[1] + 10.0 - 1e-12 &&
                  x[3] >= x[2] + 10.0 - 1e-12, "chain is feasible");
        ARO_CHECK(bs.cycle_broken == 0, "an acyclic chain breaks nothing");
        ARO_CHECK(fabs(x[0]) < 1e-12, "the source never moves");
    }

    /* (e) a positive cycle is broken, not spun on. */
    {
        double x[3] = {0.0, 0.0, 0.0};
        MonotoneQpBound b[3];
        b[0].lo = 0; b[0].hi = 1; b[0].lower = 5.0; b[0].stroke = -1; b[0].edge = 0;
        b[1].lo = 1; b[1].hi = 2; b[1].lower = 5.0; b[1].stroke = -1; b[1].edge = 1;
        b[2].lo = 2; b[2].hi = 0; b[2].lower = 5.0; b[2].stroke = -1; b[2].edge = 2;
        AtlasRadialOrderBootstrapStats bs;
        int rc = AtlasRadialOrder_bootstrap(arena, b, 3, NULL, 0, 3, 0.0,
                                            x, NULL, &bs);
        ARO_CHECK(rc == 0, "cycle bootstrap still returns a feasible point");
        ARO_CHECK(bs.cycle_broken >= 1, "the cycle was reported");
    }

    /* (f) an anchor blocks rather than being dragged. */
    {
        double x[2] = {0.0, 0.0};
        MonotoneQpBound b[1];
        b[0].lo = 0; b[0].hi = 1; b[0].lower = 7.0; b[0].stroke = -1; b[0].edge = 0;
        MonotoneQpAnchor a[1];
        a[0].var = 1; a[0].value = 0.0; a[0].component = 0;
        AtlasRadialOrderBootstrapStats bs;
        ARO_CHECK(AtlasRadialOrder_bootstrap(arena, b, 1, a, 1, 2, 0.0,
                                             x, NULL, &bs) == 0,
                  "anchored bootstrap");
        ARO_CHECK(bs.anchor_blocked == 1, "the blocked bound was counted");
        ARO_CHECK(fabs(x[1]) < 1e-12, "the anchor held");
    }

    Arena_dispose(&arena);
    fprintf(stderr, "[atlas_radial_order] selftest OK\n");
    return 0;
}
