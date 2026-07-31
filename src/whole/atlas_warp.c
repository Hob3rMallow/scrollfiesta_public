#include "atlas_warp.h"

#include "../flatten/sparse_solve.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void AtlasWarpOptions_default(AtlasWarpOptions *opts)
{
    if (opts == NULL) return;
    opts->lambda_prior = 0.05;
    opts->lambda_spacing = 0.25;
    opts->irls_rounds = 4;
    opts->irls_scale = 0.35;
    opts->hard_bounds = 0;
    opts->keep_anchors = 0;
}

static int aw_double_cmp(const void *va, const void *vb)
{
    double a = *(const double *)va, b = *(const double *)vb;
    if (a < b) return -1;
    if (a > b) return 1;
    return 0;
}

typedef struct { const uint64_t *key; } AwKeyCtx;
static const uint64_t *aw_sort_key;
static int aw_key_cmp(const void *va, const void *vb)
{
    size_t a = *(const size_t *)va, b = *(const size_t *)vb;
    if (aw_sort_key[a] < aw_sort_key[b]) return -1;
    if (aw_sort_key[a] > aw_sort_key[b]) return 1;
    return 0;
}
static void aw_sort_index(const uint64_t *key, size_t *ord, size_t n)
{
    aw_sort_key = key;
    qsort(ord, n, sizeof *ord, aw_key_cmp);
}

static size_t aw_count_violations(const MonotoneQpBound *bounds, size_t n,
                                  const double *x, double tolerance)
{
    size_t bad = 0;
    for (size_t i = 0; i < n; i++)
        if (x[bounds[i].hi] - x[bounds[i].lo] < bounds[i].lower - tolerance)
            bad++;
    return bad;
}

/* Residual of the radial spacing evidence: |du - one circumference| over the
 * pairs whose wrap count was unambiguous. */
static void aw_spacing_residual(Arena_T arena,
                                const AtlasRadialOrderSet *order,
                                const double *x,
                                double *out_rms, double *out_p95)
{
    *out_rms = 0.0;
    *out_p95 = 0.0;
    if (order->nrows == 0) return;
    Arena_Mark mark = Arena_save(arena);
    double *res = (double *)ARENA_ALLOC(arena, order->nrows * sizeof *res);
    double sum2 = 0.0;
    size_t n = 0;
    for (size_t i = 0; i < order->nrows; i++) {
        const MonotoneQpRow *row = &order->rows[i];
        int32_t hi = order->coeff[row->first].var;
        int32_t lo = order->coeff[row->first + 1].var;
        double d = fabs((x[hi] - x[lo]) - row->target);
        res[n++] = d;
        sum2 += d * d;
    }
    *out_rms = sqrt(sum2 / (double)n);
    qsort(res, n, sizeof *res, aw_double_cmp);
    *out_p95 = res[(size_t)((double)n * 0.95)];
    Arena_restore(arena, mark);
}

int AtlasWarp_solve(Arena_T arena,
                    const AtlasStripProblem *problem,
                    const AtlasStripOptions *strip_opts,
                    const MonotoneQpOptions *qp_opts,
                    const AtlasRadialOrderSet *order,
                    const double *base_sample,
                    const AtlasWarpOptions *input_opts,
                    double *out_sample,
                    AtlasWarpStats *stats)
{
    if (problem == NULL || strip_opts == NULL || order == NULL ||
        base_sample == NULL || out_sample == NULL)
        return -1;

    AtlasWarpOptions opts;
    if (input_opts != NULL) opts = *input_opts;
    else AtlasWarpOptions_default(&opts);
    if (!(opts.lambda_prior > 0.0) && opts.keep_anchors == 0) return -1;

    MonotoneQpOptions qp_defaults;
    MonotoneQpOptions_default(&qp_defaults);
    const MonotoneQpOptions *qp = qp_opts != NULL ? qp_opts : &qp_defaults;

    AtlasWarpStats st;
    memset(&st, 0, sizeof st);

    size_t ns = problem->nsamples;
    size_t nvar = ns + problem->ncross_sections;
    st.nvar = nvar;

    /* The StrokeStrip system that produced HEAD, rebuilt verbatim. */
    AtlasStripOptions final_opts = *strip_opts;
    final_opts.mode = ATLAS_STRIP_FINAL;
    AtlasStripSystem system;
    st.stage_rc = -1;                     /* strip assembly */
    if (AtlasStrip_build(arena, problem, &final_opts, &system) != 0) {
        if (stats != NULL) *stats = st;
        return -1;
    }

    const MonotoneQpProblem *base = &system.qp;
    size_t nrows = base->nrows + ns + order->nrows;
    size_t ncoeff = base->ncoeff + ns + order->ncoeff;
    size_t nbounds = base->nbounds + order->nbounds;

    MonotoneQpRow *rows =
        (MonotoneQpRow *)ARENA_ALLOC(arena, nrows * sizeof *rows);
    MonotoneQpCoeff *coeff =
        (MonotoneQpCoeff *)ARENA_ALLOC(arena, ncoeff * sizeof *coeff);
    MonotoneQpBound *bounds =
        (MonotoneQpBound *)ARENA_ALLOC(arena, nbounds * sizeof *bounds);

    memcpy(rows, base->rows, base->nrows * sizeof *rows);
    memcpy(coeff, base->coeff, base->ncoeff * sizeof *coeff);
    memcpy(bounds, base->bounds, base->nbounds * sizeof *bounds);

    size_t nr = base->nrows, nk = base->ncoeff, nb = base->nbounds;

    /*
     * Proximity to the incoming field, normalized per support component.
     *
     * The obvious form -- one spring of fixed weight per sample -- makes a
     * component's resistance to moving proportional to how finely it happens
     * to be sampled, and because the penalty is quadratic in the displacement
     * the solver then prefers to spread one component's gauge error thinly
     * across all its neighbours rather than fix the component that is wrong.
     * On the three-wrap fixture that recovers 62% of a turn error while
     * dragging the two correct wraps a third of a turn each.
     *
     * Dividing by the component's sample count gives every component the same
     * total prior mass, which is what the prior is actually for: the strip
     * rows already hold each component's internal shape rigid, so the only
     * freedom left is one gauge constant per component, and pinning that is a
     * one-spring job.  It also makes the stage independent of mesh resolution.
     */
    size_t ncomp = 0;
    for (size_t s = 0; s < problem->nstrokes; s++) {
        int32_t c = problem->strokes[s].component;
        if (c >= 0 && (size_t)c + 1 > ncomp) ncomp = (size_t)c + 1;
    }
    size_t *comp_samples =
        (size_t *)ARENA_CALLOC(arena, ncomp ? ncomp : 1, sizeof *comp_samples);
    for (size_t i = 0; i < ns; i++) {
        int32_t s = problem->samples[i].stroke;
        int32_t c = s >= 0 && (size_t)s < problem->nstrokes
                  ? problem->strokes[s].component : -1;
        if (c >= 0 && (size_t)c < ncomp) comp_samples[c]++;
    }
    for (size_t i = 0; i < ns; i++) {
        double conf = problem->samples[i].confidence;
        if (!(conf > 0.0)) conf = 1.0;
        int32_t s = problem->samples[i].stroke;
        int32_t c = s >= 0 && (size_t)s < problem->nstrokes
                  ? problem->strokes[s].component : -1;
        double share = c >= 0 && (size_t)c < ncomp && comp_samples[c] > 0
                     ? (double)comp_samples[c] : 1.0;
        MonotoneQpRow *row = &rows[nr++];
        row->first = nk;
        row->count = 1;
        row->target = base_sample[i];
        row->weight = opts.lambda_prior * conf / share;
        row->kind = ATLAS_WARP_ROW_PRIOR;
        row->owner = (int32_t)i;
        coeff[nk].var = (int32_t)i;
        coeff[nk].value = 1.0;
        nk++;
    }

    /* Radial spacing, re-based onto the appended coefficient array. */
    for (size_t i = 0; i < order->nrows; i++) {
        const MonotoneQpRow *src = &order->rows[i];
        MonotoneQpRow *row = &rows[nr++];
        *row = *src;
        row->first = nk;
        row->weight = src->weight * opts.lambda_spacing;
        for (int32_t c = 0; c < src->count; c++) {
            coeff[nk] = order->coeff[src->first + (size_t)c];
            nk++;
        }
    }

    for (size_t i = 0; i < order->nbounds; i++) bounds[nb++] = order->bounds[i];

    st.stage_rc = -2;                     /* row assembly */
    if (nr != nrows || nk != ncoeff || nb != nbounds) {
        if (stats != NULL) *stats = st;
        return -1;
    }

    MonotoneQpProblem warp;
    memset(&warp, 0, sizeof warp);
    warp.nvar = nvar;
    warp.rows = rows;
    warp.nrows = nrows;
    warp.coeff = coeff;
    warp.ncoeff = ncoeff;
    warp.bounds = bounds;
    warp.nbounds = nbounds;
    if (opts.keep_anchors) {
        warp.anchors = problem->anchors;
        warp.nanchors = problem->nanchors;
    }

    st.nrows = nrows;
    st.nbounds = nbounds;
    st.ray_bounds = order->nbounds;
    st.spacing_rows = order->nrows;

    /* Start from the incoming field and make it feasible. */
    double *x = (double *)ARENA_ALLOC(arena, nvar * sizeof *x);
    memcpy(x, base_sample, ns * sizeof *x);
    AtlasStrip_initialize_intercepts(problem, &final_opts, x);

    st.ray_bounds_violated_before =
        aw_count_violations(order->bounds, order->nbounds, x,
                            qp->feasibility_tolerance);
    aw_spacing_residual(arena, order, x, &st.spacing_rms_before,
                        &st.spacing_p95_before);

    double span_lo = x[0], span_hi = x[0];
    for (size_t i = 1; i < ns; i++) {
        if (x[i] < span_lo) span_lo = x[i];
        if (x[i] > span_hi) span_hi = x[i];
    }
    st.span_before = span_hi - span_lo;

    st.stage_rc = 0;
    if (opts.hard_bounds) {
        uint8_t *disabled =
            (uint8_t *)ARENA_ALLOC(arena, nbounds * sizeof *disabled);
        st.stage_rc = -3;                 /* feasibility bootstrap */
        if (AtlasRadialOrder_bootstrap(arena, bounds, nbounds,
                                       warp.anchors, warp.nanchors, nvar,
                                       qp->feasibility_tolerance, x,
                                       disabled, &st.bootstrap) != 0) {
            if (stats != NULL) *stats = st;
            return -1;
        }
        st.stage_rc = 0;
        /* A bound the bootstrap could not satisfy is contradictory evidence,
         * not a solver limit; drop it so the QP sees a consistent problem and
         * the count is reported rather than hidden inside a failure code. */
        size_t live = 0;
        for (size_t i = 0; i < nbounds; i++) {
            if (disabled[i]) {
                if (bounds[i].stroke < 0) st.ray_bounds_disabled++;
                continue;
            }
            bounds[live++] = bounds[i];
        }
        warp.nbounds = live;
        st.nbounds = live;
    } else {
        /* Bounds off: the strip's own stroke monotonicity goes with them, so
         * the solve is unconstrained least squares and MonotoneQp reduces to
         * a single SPD factorization per round. */
        warp.nbounds = 0;
        st.nbounds = 0;
    }

    AtlasStrip_initialize_intercepts(problem, &final_opts, x);

    /*
     * IRLS over the radial evidence.  A quarter of the pairs disagree with the
     * spiral about which way the separation goes -- turn-off errors upstream,
     * or a neighbour found across a fold -- and plain least squares would
     * split the difference with them.  Reweighting by residual lets the
     * consistent majority win and reports how much was rejected.
     */
    MonotoneQpStats qp_stats;
    memset(&qp_stats, 0, sizeof qp_stats);
    size_t first_spacing = base->nrows + ns;
    int rounds = opts.irls_rounds < 0 ? 0 : opts.irls_rounds;
    double *resid = order->nrows > 0
        ? (double *)ARENA_ALLOC(arena, order->nrows * sizeof *resid) : NULL;

    for (int round = 0; round <= rounds; round++) {
        st.qp_rc = MonotoneQp_solve(arena, &warp, qp, x, NULL, &qp_stats);
        st.irls_rounds_run = round;
        if (getenv("ATLAS_WARP_TRACE") != NULL) {
            double rr = 0.0, pp = 0.0;
            aw_spacing_residual(arena, order, x, &rr, &pp);
            fprintf(stderr, "[atlas_warp] round %d rc=%d it=%d obj %.6g->%.6g"
                            " spacing_rms=%.6g w0=%.6g\n",
                    round, st.qp_rc, qp_stats.iterations,
                    qp_stats.objective_initial, qp_stats.objective_final, rr,
                    order->nrows ? rows[first_spacing].weight : 0.0);
        }
        if (st.qp_rc != 0) break;
        if (round == rounds || order->nrows == 0) break;

        /*
         * The residual is scored as a FRACTION OF THIS PAIR'S OWN TURN, not in
         * voxels.  One circumference spans 1179 to 4782 vox across the 4x5x5
         * crop, so a single voxel-valued scale would tolerate a whole turn of
         * error at the rim while rejecting a tenth of one near the core.
         *
         * The scale is also fixed, not estimated from the spread.  A
         * MAD-derived sigma looks reasonable and is catastrophic here: after a
         * good round every residual is nearly identical, MAD collapses to
         * ~1e-9, and the reweighting reads perfect agreement as zero tolerance
         * and discards all the evidence -- measured, the fixture sprang right
         * back to its broken input with weights at 4e-8.
         */
        double floor_turn = 1.0;
        st.spacing_rows_downweighted = 0;
        double kept = 0.0;
        for (size_t i = 0; i < order->nrows; i++) {
            const MonotoneQpRow *row = &rows[first_spacing + i];
            int32_t hi = coeff[row->first].var;
            int32_t lo = coeff[row->first + 1].var;
            double turn = fabs(row->target);
            if (turn < floor_turn) turn = floor_turn;
            resid[i] = fabs((x[hi] - x[lo]) - row->target);
            double e = resid[i] / (turn * opts.irls_scale);
            double like = exp(-0.5 * e * e);
            if (like < 1e-4) like = 1e-4;
            if (like < 0.5) st.spacing_rows_downweighted++;
            kept += like;
            rows[first_spacing + i].weight =
                opts.lambda_spacing * order->rows[i].weight * like;
        }
        st.irls_sigma = opts.irls_scale;
        st.spacing_weight_retained = kept / (double)order->nrows;
    }
    st.qp_iterations = qp_stats.iterations;
    st.qp_active_final = qp_stats.active_final;
    st.objective_initial = qp_stats.objective_initial;
    st.objective_final = qp_stats.objective_final;
    st.min_slack = qp_stats.min_slack;
    st.stationarity = qp_stats.stationarity_residual;

    if (st.qp_rc == 0) {
        /*
         * Displacement is reported relative to the MEDIAN displacement, not to
         * zero.  A constant added to every u is pure gauge -- it moves no
         * geometry relative to anything -- and a solve that legitimately fixes
         * one component often lands on an equivalent solution offset globally.
         * Measuring raw displacement would read that as a teleport and fire
         * the alarm that is supposed to catch real ones.
         */
        double *shift = (double *)ARENA_ALLOC(arena, ns * sizeof *shift);
        for (size_t i = 0; i < ns; i++) shift[i] = x[i] - base_sample[i];
        qsort(shift, ns, sizeof *shift, aw_double_cmp);
        st.gauge_shift = ns > 0 ? shift[ns / 2] : 0.0;

        double sum2 = 0.0;
        span_lo = x[0];
        span_hi = x[0];
        for (size_t i = 0; i < ns; i++) {
            double d = fabs(x[i] - base_sample[i] - st.gauge_shift);
            shift[i] = d;
            sum2 += d * d;
            if (d > st.shift_max) st.shift_max = d;
            if (x[i] < span_lo) span_lo = x[i];
            if (x[i] > span_hi) span_hi = x[i];
        }
        st.shift_rms = ns > 0 ? sqrt(sum2 / (double)ns) : 0.0;
        st.span_after = span_hi - span_lo;
        qsort(shift, ns, sizeof *shift, aw_double_cmp);
        st.shift_p95 = ns > 0 ? shift[(size_t)((double)ns * 0.95)] : 0.0;
        st.ray_bounds_violated_after =
            aw_count_violations(order->bounds, order->nbounds, x,
                                qp->feasibility_tolerance);
        aw_spacing_residual(arena, order, x, &st.spacing_rms_after,
                            &st.spacing_p95_after);
        memcpy(out_sample, x, ns * sizeof *out_sample);
    }

    if (stats != NULL) *stats = st;
    return st.qp_rc == 0 ? 0 : -1;
}

/* ------------------------------------------------------- per-component gauge */

int AtlasWarp_solve_gauges(Arena_T arena,
                           const AtlasRadialOrderSet *order,
                           const int32_t *sample_component,
                           size_t nsamples,
                           size_t ncomponents,
                           const double *base_sample,
                           const AtlasWarpOptions *input_opts,
                           double *out_shift,
                           AtlasWarpGaugeStats *stats)
{
    if (order == NULL || sample_component == NULL || base_sample == NULL ||
        out_shift == NULL || ncomponents == 0)
        return -1;

    AtlasWarpOptions opts;
    if (input_opts != NULL) opts = *input_opts;
    else AtlasWarpOptions_default(&opts);

    AtlasWarpGaugeStats st;
    memset(&st, 0, sizeof st);
    st.components = ncomponents;

    Arena_Mark mark = Arena_save(arena);
    size_t n = ncomponents;

    /* Gather the cross-component equations once: shift[b] - shift[a] = rhs. */
    size_t cap = order->npairs ? order->npairs : 1;
    int32_t *ea = (int32_t *)ARENA_ALLOC(arena, cap * sizeof *ea);
    int32_t *eb = (int32_t *)ARENA_ALLOC(arena, cap * sizeof *eb);
    double *erhs = (double *)ARENA_ALLOC(arena, cap * sizeof *erhs);
    double *ew = (double *)ARENA_ALLOC(arena, cap * sizeof *ew);
    size_t ne = 0;
    for (size_t i = 0; i < order->npairs; i++) {
        const AtlasRadialOrderPair *p = &order->pairs[i];
        if ((size_t)p->inner >= nsamples || (size_t)p->outer >= nsamples)
            continue;
        int32_t a = sample_component[p->inner];
        int32_t b = sample_component[p->outer];
        if (a < 0 || b < 0 || (size_t)a >= n || (size_t)b >= n) continue;
        if (a == b) { st.intra_skipped++; continue; }
        ea[ne] = a;
        eb[ne] = b;
        /* Observed separation must become the spiral's; the shift makes up
         * the difference. */
        erhs[ne] = p->spacing_target -
                   (base_sample[p->outer] - base_sample[p->inner]);
        ew[ne] = p->weight > 0.0 ? p->weight : 1.0;
        ne++;
    }
    st.equations = ne;

    double *shift = (double *)ARENA_CALLOC(arena, n, sizeof *shift);
    double *diag = (double *)ARENA_ALLOC(arena, n * sizeof *diag);
    double *rhs = (double *)ARENA_ALLOC(arena, n * sizeof *rhs);
    double *like = (double *)ARENA_ALLOC(arena, (ne ? ne : 1) * sizeof *like);
    for (size_t i = 0; i < ne; i++) like[i] = 1.0;

    /*
     * The normal equations of a difference system are a weighted graph
     * Laplacian plus the prior on the diagonal.  Gauss-Seidel on it does NOT
     * converge usefully -- the diagonal is dominated by the tiny prior and the
     * low-frequency modes need far more sweeps than any sane cap (measured:
     * every round hit 2000 sweeps still moving 0.86 per sweep, and the
     * residual ended worse than it started).  The system is only a few hundred
     * square, so it is assembled and factored directly instead.
     */
    int32_t *tri_r = (int32_t *)ARENA_ALLOC(arena, (ne + n) * sizeof *tri_r);
    int32_t *tri_c = (int32_t *)ARENA_ALLOC(arena, (ne + n) * sizeof *tri_c);
    double *tri_v = (double *)ARENA_ALLOC(arena, (ne + n) * sizeof *tri_v);
    int *coo_r = (int *)ARENA_ALLOC(arena, (ne + n) * sizeof *coo_r);
    int *coo_c = (int *)ARENA_ALLOC(arena, (ne + n) * sizeof *coo_c);
    double *coo_v = (double *)ARENA_ALLOC(arena, (ne + n) * sizeof *coo_v);
    uint64_t *key = (uint64_t *)ARENA_ALLOC(arena, (ne ? ne : 1) * sizeof *key);
    size_t *ord = (size_t *)ARENA_ALLOC(arena, (ne ? ne : 1) * sizeof *ord);
    (void)tri_r; (void)tri_c; (void)tri_v;

    int rounds = opts.irls_rounds < 0 ? 0 : opts.irls_rounds;
    for (int round = 0; round <= rounds; round++) {
        for (size_t c = 0; c < n; c++) {
            diag[c] = opts.lambda_prior;
            rhs[c] = 0.0;
        }
        for (size_t i = 0; i < ne; i++) {
            double w = ew[i] * like[i];
            diag[ea[i]] += w;
            diag[eb[i]] += w;
            rhs[ea[i]] -= w * erhs[i];
            rhs[eb[i]] += w * erhs[i];
            uint32_t lo = (uint32_t)(ea[i] < eb[i] ? ea[i] : eb[i]);
            uint32_t hi = (uint32_t)(ea[i] < eb[i] ? eb[i] : ea[i]);
            key[i] = ((uint64_t)hi << 32) | (uint64_t)lo;
            ord[i] = i;
        }
        /* Merge the off-diagonals: many pairs share a component pair. */
        aw_sort_index(key, ord, ne);
        size_t nt = 0;
        for (size_t c = 0; c < n; c++) {
            coo_r[nt] = (int)c;
            coo_c[nt] = (int)c;
            coo_v[nt] = diag[c];
            nt++;
        }
        size_t i = 0;
        while (i < ne) {
            size_t j = i;
            double acc = 0.0;
            while (j < ne && key[ord[j]] == key[ord[i]]) {
                acc += ew[ord[j]] * like[ord[j]];
                j++;
            }
            uint64_t k = key[ord[i]];
            coo_r[nt] = (int)(uint32_t)(k >> 32);
            coo_c[nt] = (int)(uint32_t)(k & 0xffffffffu);
            coo_v[nt] = -acc;
            nt++;
            i = j;
        }
        if (Sparse_solve_sym((int)n, (int)nt, coo_r, coo_c, coo_v,
                             rhs, shift, SPARSE_SPD) != 0) {
            st.solve_failed = 1;
            break;
        }

        double sum2 = 0.0;
        for (size_t e = 0; e < ne; e++) {
            double r = (shift[eb[e]] - shift[ea[e]]) - erhs[e];
            sum2 += r * r;
        }
        double rms = ne ? sqrt(sum2 / (double)ne) : 0.0;
        if (round == 0) st.residual_rms_before = rms;
        st.residual_rms_after = rms;
        st.irls_rounds_run = round;
        if (round == rounds || ne == 0) break;

        /* Same fraction-of-a-turn robustness as the per-sample warp. */
        st.downweighted = 0;
        for (size_t e = 0; e < ne; e++) {
            double r = fabs((shift[eb[e]] - shift[ea[e]]) - erhs[e]);
            double turn = fabs(erhs[e]) > 1.0 ? fabs(erhs[e]) : 1.0;
            double q = r / (turn * opts.irls_scale);
            like[e] = exp(-0.5 * q * q);
            if (like[e] < 1e-4) like[e] = 1e-4;
            if (like[e] < 0.5) st.downweighted++;
        }
    }

    /* Report displacement against the median shift: a global constant is the
     * same atlas, not a move. */
    double *sorted = (double *)ARENA_ALLOC(arena, n * sizeof *sorted);
    memcpy(sorted, shift, n * sizeof *sorted);
    qsort(sorted, n, sizeof *sorted, aw_double_cmp);
    st.gauge_shift = sorted[n / 2];
    double sum2 = 0.0;
    for (size_t c = 0; c < n; c++) {
        double d = fabs(shift[c] - st.gauge_shift);
        sum2 += d * d;
        if (d > st.shift_max) st.shift_max = d;
    }
    st.shift_rms = sqrt(sum2 / (double)n);
    st.final_delta = 0.0;   /* direct solve: no iteration to report */

    /* Components the evidence never links to anything are pinned only by the
     * tiny prior; report them, because their placement is not determined. */
    {
        Arena_Mark m2 = Arena_save(arena);
        int32_t *par = (int32_t *)ARENA_ALLOC(arena, n * sizeof *par);
        for (size_t c = 0; c < n; c++) par[c] = (int32_t)c;
        for (size_t i = 0; i < ne; i++) {
            int32_t ra = ea[i], rb = eb[i];
            while (par[ra] != ra) ra = par[ra];
            while (par[rb] != rb) rb = par[rb];
            if (ra != rb) par[ra] = rb;
        }
        size_t islands = 0;
        for (size_t c = 0; c < n; c++) if (par[c] == (int32_t)c) islands++;
        st.equation_islands = islands;
        Arena_restore(arena, m2);
    }

    memcpy(out_shift, shift, n * sizeof *out_shift);
    if (stats != NULL) *stats = st;
    Arena_restore(arena, mark);
    return 0;
}

/* ----------------------------------------------------------------- selftest */

#define AW_CHECK(cond, msg)                                                   \
    do {                                                                      \
        if (!(cond)) {                                                        \
            fprintf(stderr, "[atlas_warp] FAIL %s (%s:%d)\n",                 \
                    (msg), __FILE__, __LINE__);                               \
            Arena_dispose(&arena);                                            \
            return -1;                                                        \
        }                                                                     \
    } while (0)

/*
 * Fixture: a genuine Archimedean spiral, r(phi) = r0 + pitch*phi/(2*pi), run
 * for three turns across two axial slice planes and then CUT into one stroke
 * per turn.  The cut is what makes this the interesting case: each turn is its
 * own connected component, so the relative gauges are precisely the
 * unobservable that the ordering evidence has to supply.  Within a turn the
 * two planes are tied by cross-sections, as on real data.
 *
 * The spiral has to be real, not three concentric circles.  On the true spiral
 * the exact u-difference between radial neighbours,
 *
 *     u(phi + 2pi) - u(phi) = 2*pi*r0 + pitch*phi + pi*pitch,
 *
 * equals the module's target 2*pi*r_mid identically; on concentric circles it
 * only equals it on average, leaving a +/- pi*pitch wobble that the solver
 * then spreads over the wraps and that no assertion could tell apart from a
 * real defect.
 *
 * Turn 1 then arrives a full circumference short -- the free-floating gauge
 * error the integer-winding stage used to "fix" by teleporting a chart.
 */
int AtlasWarp_selftest(void)
{
    Arena_T arena = Arena_new();
    const double pitch = 12.0, r0 = 60.0, slice = 8.0;
    const int nwrap = 3, nang = 24, nplane = 2;
    const double twopi = 6.283185307179586476925286766559;

    ScaffoldCalib cal;
    memset(&cal, 0, sizeof cal);
    cal.pitch = pitch;
    cal.spiral_a = r0;   /* makes one turn cost exactly one circumference */
    cal.spiral_b = -pitch;
    cal.sense = -1;
    cal.axis_dir[0] = 1.0f;

    size_t nstroke = (size_t)(nwrap * nplane);
    size_t n = nstroke * (size_t)nang;
    AtlasStripSample *samples =
        (AtlasStripSample *)ARENA_ALLOC(arena, n * sizeof *samples);
    AtlasStripStroke *strokes =
        (AtlasStripStroke *)ARENA_ALLOC(arena, nstroke * sizeof *strokes);
    int32_t *plane = (int32_t *)ARENA_ALLOC(arena, n * sizeof *plane);
    int32_t *comp = (int32_t *)ARENA_ALLOC(arena, n * sizeof *comp);
    double *phi = (double *)ARENA_ALLOC(arena, n * sizeof *phi);
    double *truth = (double *)ARENA_ALLOC(arena, n * sizeof *truth);
    double *base = (double *)ARENA_ALLOC(arena, n * sizeof *base);
    double *solved = (double *)ARENA_ALLOC(arena, n * sizeof *solved);

    size_t k = 0;
    for (int p = 0; p < nplane; p++) {
        for (int w = 0; w < nwrap; w++) {
            size_t sid = (size_t)p * (size_t)nwrap + (size_t)w;
            strokes[sid].first = (int32_t)k;
            strokes[sid].count = nang;
            strokes[sid].component = w;
            strokes[sid].resolved = 1;
            /*
             * On a left-handed scroll u decreases outward, so it also
             * decreases with phi along the sheet.  A stroke's arclength must
             * run the way u increases -- which is what AtlasCandidates
             * enforces on real data by flipping any stroke whose arclength
             * anti-correlates with raw u -- so the samples are laid down from
             * the outer end of the turn inward.
             */
            for (int a = 0; a < nang; a++) {
                int step = nang - 1 - a;
                double turn = (double)w + (double)step / (double)nang;
                double phi_signed = -twopi * turn;   /* sense = -1 */
                double theta = twopi * (double)step / (double)nang;
                double r = r0 + pitch * turn;
                memset(&samples[k], 0, sizeof samples[k]);
                samples[k].p[0] = (double)p * slice;
                samples[k].p[1] = r * cos(theta);
                samples[k].p[2] = r * sin(theta);
                samples[k].tangent[0] = 0.0;
                samples[k].tangent[1] = sin(theta);
                samples[k].tangent[2] = -cos(theta);
                /* The pinned spiral map u(phi) = a*phi + b*phi^2/(4pi). */
                truth[k] = cal.spiral_a * phi_signed +
                           cal.spiral_b * phi_signed * phi_signed /
                           (2.0 * twopi);
                samples[k].s = truth[k];   /* monotone along the stroke */
                samples[k].stroke = (int32_t)sid;
                samples[k].ordinal = a;
                samples[k].confidence = 1.0;
                plane[k] = p;
                comp[k] = w;
                phi[k] = phi_signed;
                base[k] = truth[k];
                k++;
            }
        }
    }
    double turn_error = twopi * r0;
    for (int p = 0; p < nplane; p++) {
        size_t sid = (size_t)p * (size_t)nwrap + 1;
        for (int a = 0; a < nang; a++)
            base[(size_t)strokes[sid].first + (size_t)a] -= turn_error;
    }

    /* One cross-section per (wrap, angle) tying the two planes together. */
    size_t ncs = (size_t)(nwrap * nang);
    AtlasStripCrossSection *cs = (AtlasStripCrossSection *)ARENA_ALLOC(
        arena, ncs * sizeof *cs);
    AtlasStripMember *members = (AtlasStripMember *)ARENA_ALLOC(
        arena, ncs * (size_t)nplane * sizeof *members);
    size_t nm = 0, nc = 0;
    for (int w = 0; w < nwrap; w++) {
        for (int a = 0; a < nang; a++) {
            /* Same reversed ordinal the samples were laid down with. */
            int step = nang - 1 - a;
            double theta = twopi * (double)step / (double)nang;
            double r = r0 + pitch * ((double)w + (double)step / (double)nang);
            cs[nc].first = nm;
            cs[nc].count = nplane;
            cs[nc].tangent[0] = 0.0;
            cs[nc].tangent[1] = sin(theta);
            cs[nc].tangent[2] = -cos(theta);
            cs[nc].weight = 1.0;
            cs[nc].id = (int32_t)nc;
            for (int p = 0; p < nplane; p++) {
                size_t sid = (size_t)p * (size_t)nwrap + (size_t)w;
                int32_t first = strokes[sid].first;
                int32_t idx = first + a;
                int32_t lo = a > 0 ? idx - 1 : idx;
                int32_t hi = a > 0 ? idx : idx + 1;
                memset(&members[nm], 0, sizeof members[nm]);
                members[nm].value0 = idx;
                members[nm].value1 = -1;
                members[nm].value_t = 0.0;
                members[nm].deriv_lo = lo;
                members[nm].deriv_hi = hi;
                members[nm].deriv_length =
                    samples[hi].s - samples[lo].s;
                members[nm].p[0] = samples[idx].p[0];
                members[nm].p[1] = samples[idx].p[1];
                members[nm].p[2] = samples[idx].p[2];
                members[nm].tangent[0] = samples[idx].tangent[0];
                members[nm].tangent[1] = samples[idx].tangent[1];
                members[nm].tangent[2] = samples[idx].tangent[2];
                members[nm].dual_width = r * twopi / (double)nang;
                members[nm].membership = 1.0;
                members[nm].base_weight = 1.0;
                members[nm].observation = (int32_t)nm;
                nm++;
            }
            nc++;
        }
    }

    AtlasStripProblem problem;
    memset(&problem, 0, sizeof problem);
    problem.samples = samples;
    problem.nsamples = n;
    problem.strokes = strokes;
    problem.nstrokes = nstroke;
    problem.members = members;
    problem.nmembers = nm;
    problem.cross_sections = cs;
    problem.ncross_sections = nc;

    AtlasStripOptions strip_opts;
    AtlasStripOptions_default(&strip_opts);
    MonotoneQpOptions qp_opts;
    MonotoneQpOptions_default(&qp_opts);

    AtlasRadialOrderOptions order_opts;
    AtlasRadialOrderOptions_default(&order_opts);
    order_opts.core_radius_pitches = 1.0;
    order_opts.arc_window = twopi * r0 / (double)nang * 0.75;

    AtlasRadialOrderSet order;
    AW_CHECK(AtlasRadialOrder_build(arena, samples, plane, comp, phi, NULL,
                                    n, base,
                                    &cal, NULL, 0, &order_opts, &order) == 0,
             "radial order build");
    AW_CHECK(order.npairs == (size_t)((nwrap - 1) * nang * nplane),
             "ray coverage on every plane");

    AtlasWarpOptions warp_opts;
    AtlasWarpOptions_default(&warp_opts);
    AtlasWarpStats st;
    int rc = AtlasWarp_solve(arena, &problem, &strip_opts, &qp_opts, &order,
                             base, &warp_opts, solved, &st);
    if (rc != 0)
        fprintf(stderr, "[atlas_warp] solve rc=%d qp_rc=%d bounds=%zu "
                        "violated_before=%zu disabled=%zu cycle=%zu\n",
                rc, st.qp_rc, st.nbounds, st.ray_bounds_violated_before,
                st.ray_bounds_disabled, st.bootstrap.cycle_broken);
    AW_CHECK(rc == 0, "warp solve");

    /* A whole-turn gauge error clears the 8-vox ordering bound comfortably --
     * only the spacing evidence can see it, which is exactly why the hard
     * bound alone cannot be the mechanism. */
    AW_CHECK(st.spacing_rms_before > 0.5 * turn_error,
             "the displaced wrap really is a turn out of place");
    AW_CHECK(st.ray_bounds_violated_after == 0, "ordering holds");
    AW_CHECK(st.ray_bounds_disabled == 0, "no evidence had to be discarded");

    /*
     * The displaced wrap must come back; the two correct ones must not be
     * dragged along.  That asymmetry is the whole point -- rigid layout moved
     * everything to make room.  Displacement is measured against the solve's
     * own gauge, because adding a constant to every u yields the same atlas.
     */
    double moved_displaced = 0.0, moved_correct = 0.0;
    for (int p = 0; p < nplane; p++) {
        for (int w = 0; w < nwrap; w++) {
            size_t sid = (size_t)p * (size_t)nwrap + (size_t)w;
            for (int a = 0; a < nang; a++) {
                size_t i = (size_t)strokes[sid].first + (size_t)a;
                double d = fabs(solved[i] - base[i] - st.gauge_shift);
                if (w == 1) { if (d > moved_displaced) moved_displaced = d; }
                else if (d > moved_correct) moved_correct = d;
            }
        }
    }
    /* And the recovered field must match the truth, not merely be ordered. */
    double worst = 0.0;
    for (size_t i = 0; i < n; i++) {
        double e = fabs((solved[i] - st.gauge_shift) - truth[i]);
        if (e > worst) worst = e;
    }
    fprintf(stderr, "[atlas_warp] fixture: gauge=%.4g displaced=%.4g "
                    "correct=%.4g turn=%.4g spacing %.4g->%.4g p95 %.4g->%.4g "
                    "worst=%.4g qp(it=%d rc=%d obj %.6g->%.6g) bounds=%zu\n",
            st.gauge_shift, moved_displaced, moved_correct, turn_error,
            st.spacing_rms_before, st.spacing_rms_after,
            st.spacing_p95_before, st.spacing_p95_after, worst,
            st.qp_iterations, st.qp_rc, st.objective_initial,
            st.objective_final, st.nbounds);
    AW_CHECK(moved_displaced > 0.9 * turn_error,
             "the displaced wrap was carried back a full turn");
    AW_CHECK(moved_correct < 0.02 * moved_displaced,
             "the correct wraps stayed put");
    AW_CHECK(worst < 0.10 * turn_error, "the true atlas was recovered");

    /* Rerunning on an already-ordered field must be a near-no-op: this is the
     * "preserve good regions" contract that rigid layout violated. */
    AtlasRadialOrderSet order2;
    AW_CHECK(AtlasRadialOrder_build(arena, samples, plane, comp, phi, NULL,
                                    n, truth,
                                    &cal, NULL, 0, &order_opts, &order2) == 0,
             "radial order on the true field");
    AtlasWarpStats st2;
    double *solved2 = (double *)ARENA_ALLOC(arena, n * sizeof *solved2);
    AW_CHECK(AtlasWarp_solve(arena, &problem, &strip_opts, &qp_opts, &order2,
                             truth, &warp_opts, solved2, &st2) == 0,
             "warp solve on the true field");
    AW_CHECK(st2.spacing_rms_before < 1e-6 * turn_error,
             "the true field already satisfies the spacing evidence");
    AW_CHECK(st2.shift_max < 1e-6 * turn_error,
             "an already-correct field is left alone");

    Arena_dispose(&arena);
    fprintf(stderr, "[atlas_warp] selftest OK\n");
    return 0;
}
