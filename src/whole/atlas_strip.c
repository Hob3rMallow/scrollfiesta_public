#include "atlas_strip.h"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

static double dot3(const double a[3], const double b[3])
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

double AtlasStrip_member_value(const AtlasStripMember *m, const double *u)
{
    if (m == NULL || u == NULL || m->value0 < 0) return 0.0;
    if (m->value1 < 0) return u[m->value0];
    return (1.0 - m->value_t) * u[m->value0] +
           m->value_t * u[m->value1];
}

void AtlasStripOptions_default(AtlasStripOptions *opts)
{
    if (opts == NULL) return;
    opts->mode = ATLAS_STRIP_RELAXED;
    opts->lambda_length = 1.0;
    opts->lambda_align = 1.0;
    opts->lambda_local = 0.05;
    opts->lambda_continuation = 1.0;
    opts->monotone_fraction = 0.5;
    opts->membership_floor = 1e-8;
    opts->lambda_prior = 0.0;
}

void AtlasStripRobustOptions_default(AtlasStripRobustOptions *opts)
{
    if (opts == NULL) return;
    opts->l1_iterations = 3;
    opts->likelihood_iterations = 7;
    opts->initial_likelihood = 1e-5;
    opts->l1_epsilon = 0.25;
    opts->likelihood_sigma = 0.0;
    opts->likelihood_floor = 1e-7;
    opts->convergence_tolerance = 1e-8;
    opts->use_active_set = 0;
}

static int validate_strip_problem(const AtlasStripProblem *p,
                                  const AtlasStripOptions *o)
{
    if (p == NULL || o == NULL || p->samples == NULL || p->nsamples == 0 ||
        p->strokes == NULL || p->nstrokes == 0 ||
        p->members == NULL || p->cross_sections == NULL ||
        p->ncross_sections == 0 ||
        p->nsamples + p->ncross_sections > (size_t)INT_MAX ||
        o->lambda_length < 0.0 || o->lambda_align <= 0.0 ||
        o->lambda_local < 0.0 || o->lambda_continuation < 0.0 ||
        o->monotone_fraction < 0.0 ||
        o->membership_floor < 0.0 ||
        !isfinite(o->lambda_prior) || o->lambda_prior < 0.0) return -1;

    for (size_t s = 0; s < p->nstrokes; s++) {
        const AtlasStripStroke *st = &p->strokes[s];
        if (st->first < 0 || st->count <= 0 ||
            (size_t)st->first >= p->nsamples ||
            (size_t)st->count > p->nsamples - (size_t)st->first)
            return -1;
        for (int32_t i = 0; i < st->count; i++) {
            const AtlasStripSample *sample =
                &p->samples[(size_t)st->first + (size_t)i];
            if (sample->stroke != (int32_t)s || sample->ordinal != i ||
                !isfinite(sample->s) || !isfinite(sample->confidence) ||
                sample->confidence < 0.0)
                return -1;
            for (int k = 0; k < 3; k++)
                if (!isfinite(sample->p[k]) ||
                    !isfinite(sample->tangent[k])) return -1;
            if (i > 0 &&
                sample->s <= p->samples[(size_t)st->first +
                                        (size_t)i - 1].s)
                return -1;
        }
    }
    size_t expected_member = 0;
    for (size_t c = 0; c < p->ncross_sections; c++) {
        const AtlasStripCrossSection *cs = &p->cross_sections[c];
        if (cs->count <= 0 || cs->first != expected_member ||
            cs->first > p->nmembers ||
            (size_t)cs->count > p->nmembers - cs->first ||
            !isfinite(cs->weight) || cs->weight <= 0.0)
            return -1;
        expected_member += (size_t)cs->count;
        double tn = sqrt(dot3(cs->tangent, cs->tangent));
        if (!isfinite(tn) || tn < 0.5 || tn > 1.5) return -1;
        for (int32_t j = 0; j < cs->count; j++) {
            const AtlasStripMember *m =
                &p->members[cs->first + (size_t)j];
            if (m->value0 < 0 || (size_t)m->value0 >= p->nsamples ||
                m->value1 < -1 ||
                (m->value1 >= 0 && (size_t)m->value1 >= p->nsamples) ||
                m->value_t < 0.0 || m->value_t > 1.0 ||
                m->deriv_lo < 0 || m->deriv_hi < 0 ||
                (size_t)m->deriv_lo >= p->nsamples ||
                (size_t)m->deriv_hi >= p->nsamples ||
                m->deriv_lo == m->deriv_hi ||
                !isfinite(m->deriv_length) || m->deriv_length <= 0.0 ||
                !isfinite(m->dual_width) || m->dual_width <= 0.0 ||
                !isfinite(m->membership) || m->membership < 0.0 ||
                m->membership > 1.0 ||
                !isfinite(m->base_weight) || m->base_weight < 0.0)
                return -1;
            for (int k = 0; k < 3; k++)
                if (!isfinite(m->p[k]) || !isfinite(m->tangent[k]))
                    return -1;
        }
    }
    if (expected_member != p->nmembers) return -1;
    if (p->ncontinuations > 0 && p->continuations == NULL) return -1;
    for (size_t i = 0; i < p->ncontinuations; i++) {
        const AtlasStripContinuation *cn = &p->continuations[i];
        if (cn->a < 0 || cn->b < 0 || cn->a == cn->b ||
            (size_t)cn->a >= p->nsamples ||
            (size_t)cn->b >= p->nsamples ||
            !isfinite(cn->target) || !isfinite(cn->weight) ||
            cn->weight < 0.0 || !isfinite(cn->membership) ||
            cn->membership < 0.0 || cn->membership > 1.0)
            return -1;
    }
    for (size_t a = 0; a < p->nanchors; a++)
        if (p->anchors[a].var < 0 ||
            (size_t)p->anchors[a].var >= p->nsamples)
            return -1;
    return 0;
}

int AtlasStrip_build(Arena_T arena,
                     const AtlasStripProblem *p,
                     const AtlasStripOptions *input_opts,
                     AtlasStripSystem *system)
{
    if (arena == NULL || system == NULL) return -1;
    AtlasStripOptions defaults;
    AtlasStripOptions_default(&defaults);
    const AtlasStripOptions *o = input_opts != NULL ? input_opts : &defaults;
    if (validate_strip_problem(p, o) != 0) {
        fprintf(stderr, "[atlas_strip] build: problem validation failed\n");
        return -1;
    }

    size_t nedge = 0;
    for (size_t s = 0; s < p->nstrokes; s++)
        nedge += (size_t)(p->strokes[s].count - 1);

    /*
     * The SOURCE member (j == 0) is always emitted, clamped to a floor
     * weight when the robust loop rejects it.  On an honest field a wholly
     * wrong cross-section drives every member -- source included -- to
     * circumference-scale residuals and sub-floor likelihoods; dropping
     * them all would orphan the latent intercept (singular LS) or, under
     * the old rule, abort the build.  A floor-weight source row couples
     * the source only to its own private intercept, which is a no-op, so
     * a fully rejected cross-section degrades to nothing instead.
     */
    size_t nactive_member = 0, nactive_continuation = 0, nlength = 0;
    size_t ncoeff = nedge * 2;
    for (size_t c = 0; c < p->ncross_sections; c++) {
        const AtlasStripCrossSection *cs = &p->cross_sections[c];
        size_t active_here = 0;
        for (int32_t j = 0; j < cs->count; j++) {
            const AtlasStripMember *m =
                &p->members[cs->first + (size_t)j];
            if (j > 0 &&
                m->membership * m->base_weight <= o->membership_floor)
                continue;
            active_here++;
            nactive_member++;
            ncoeff += (m->value1 >= 0 ? 3u : 2u);
            ncoeff += 2u;
        }
        if (active_here == 0) return -1;
        nlength++;
    }
    for (size_t i = 0; i < p->ncontinuations; i++)
        if (p->continuations[i].membership * p->continuations[i].weight >
            o->membership_floor)
            nactive_continuation++;
    if (nactive_continuation > (SIZE_MAX - ncoeff) / 2) return -1;
    ncoeff += 2 * nactive_continuation;
    size_t nprior = 0;
    if (p->prior_u != NULL && o->lambda_prior > 0.0)
        for (size_t i = 0; i < p->nsamples; i++)
            if (isfinite(p->prior_u[i])) nprior++;
    ncoeff += nprior;
    size_t nrows = nedge + nactive_member + nlength + nactive_continuation +
                   nprior;

    MonotoneQpRow *rows =
        (MonotoneQpRow *)ARENA_ALLOC(arena,
                                     (nrows ? nrows : 1) * sizeof(*rows));
    MonotoneQpCoeff *coeff =
        (MonotoneQpCoeff *)ARENA_ALLOC(arena,
                                       (ncoeff ? ncoeff : 1) * sizeof(*coeff));
    MonotoneQpBound *bounds =
        (MonotoneQpBound *)ARENA_ALLOC(arena,
                                       (nedge ? nedge : 1) * sizeof(*bounds));

    size_t nr = 0, nk = 0, nb = 0;

    /* StrokeStrip Eq. 6: one averaged longitudinal derivative row per C. */
    for (size_t c = 0; c < p->ncross_sections; c++) {
        const AtlasStripCrossSection *cs = &p->cross_sections[c];
        double denom = 0.0;
        for (int32_t j = 0; j < cs->count; j++) {
            const AtlasStripMember *m =
                &p->members[cs->first + (size_t)j];
            double mw = m->membership * m->base_weight;
            if (j == 0) {
                if (mw <= o->membership_floor)
                    mw = 1000.0 * o->membership_floor;
            } else if (mw <= o->membership_floor) continue;
            denom += m->dual_width * mw;
        }
        if (denom <= 0.0) continue;

        MonotoneQpRow *row = &rows[nr++];
        row->first = nk;
        row->target = 1.0;
        row->weight = o->lambda_length * cs->weight;
        row->kind = ATLAS_STRIP_ROW_LENGTH;
        row->owner = (int32_t)c;

        for (int32_t j = 0; j < cs->count; j++) {
            const AtlasStripMember *m =
                &p->members[cs->first + (size_t)j];
            double mw = m->membership * m->base_weight;
            if (j == 0) {
                if (mw <= o->membership_floor)
                    mw = 1000.0 * o->membership_floor;
            } else if (mw <= o->membership_floor) continue;
            double tangent_projection = dot3(m->tangent, cs->tangent);
            double a = (m->dual_width * mw / denom) *
                       tangent_projection / m->deriv_length;
            coeff[nk].var = m->deriv_hi;
            coeff[nk].value = a;
            nk++;
            coeff[nk].var = m->deriv_lo;
            coeff[nk].value = -a;
            nk++;
        }
        row->count = (int32_t)(nk - row->first);
    }

    /*
     * Latent cross-section intercept.  Relaxed mode uses u_a ~= q_C.
     * Final mode uses u_a - tau_C.p_a ~= q_C, equivalent to Eq. 7.
     */
    for (size_t c = 0; c < p->ncross_sections; c++) {
        const AtlasStripCrossSection *cs = &p->cross_sections[c];
        int32_t qvar = (int32_t)(p->nsamples + c);
        for (int32_t j = 0; j < cs->count; j++) {
            const AtlasStripMember *m =
                &p->members[cs->first + (size_t)j];
            double mw = m->membership * m->base_weight;
            if (j == 0) {
                if (mw <= o->membership_floor)
                    mw = 1000.0 * o->membership_floor;
            } else if (mw <= o->membership_floor) continue;
            MonotoneQpRow *row = &rows[nr++];
            row->first = nk;
            row->target = o->mode == ATLAS_STRIP_FINAL
                        ? dot3(cs->tangent, m->p) : 0.0;
            row->weight = o->lambda_align * cs->weight * mw;
            row->kind = ATLAS_STRIP_ROW_ALIGN;
            row->owner = m->observation;

            if (m->value1 >= 0) {
                coeff[nk].var = m->value0;
                coeff[nk].value = 1.0 - m->value_t;
                nk++;
                coeff[nk].var = m->value1;
                coeff[nk].value = m->value_t;
                nk++;
            } else {
                coeff[nk].var = m->value0;
                coeff[nk].value = 1.0;
                nk++;
            }
            coeff[nk].var = qvar;
            coeff[nk].value = -1.0;
            nk++;
            row->count = (int32_t)(nk - row->first);
        }
    }

    /*
     * Weak restriction of the ideal-surface differential energy.  It prevents
     * cross-section averages from hiding equal-and-opposite local speed errors
     * without allowing every noisy observation to organize the global chart.
     */
    for (size_t s = 0; s < p->nstrokes; s++) {
        const AtlasStripStroke *st = &p->strokes[s];
        for (int32_t i = 1; i < st->count; i++) {
            int32_t lo = st->first + i - 1;
            int32_t hi = st->first + i;
            double length = p->samples[hi].s - p->samples[lo].s;
            double conf = fmin(p->samples[lo].confidence,
                               p->samples[hi].confidence);

            MonotoneQpRow *row = &rows[nr++];
            row->first = nk;
            row->count = 2;
            row->target = length;
            row->weight = o->lambda_local * conf / length;
            row->kind = ATLAS_STRIP_ROW_LOCAL;
            row->owner = (int32_t)s;
            coeff[nk].var = hi;
            coeff[nk].value = 1.0;
            nk++;
            coeff[nk].var = lo;
            coeff[nk].value = -1.0;
            nk++;

            bounds[nb].lo = lo;
            bounds[nb].hi = hi;
            bounds[nb].lower = o->monotone_fraction * length;
            bounds[nb].stroke = (int32_t)s;
            bounds[nb].edge = i - 1;
            nb++;
        }
    }

    /* Same-slice fragment continuations retain their physical arc gap in
     * both modes.  These are sparse endpoint rows, not same-isovalue votes. */
    for (size_t i = 0; i < p->ncontinuations; i++) {
        const AtlasStripContinuation *cn = &p->continuations[i];
        double mw = cn->membership * cn->weight;
        if (mw <= o->membership_floor) continue;
        MonotoneQpRow *row = &rows[nr++];
        row->first = nk;
        row->count = 2;
        row->target = cn->target;
        row->weight = o->lambda_continuation * mw;
        row->kind = ATLAS_STRIP_ROW_CONTINUATION;
        row->owner = cn->observation;
        coeff[nk].var = cn->a;
        coeff[nk].value = 1.0;
        nk++;
        coeff[nk].var = cn->b;
        coeff[nk].value = -1.0;
        nk++;
    }

    /*
     * Weak geometric prior: u_i ~= prior_u[i].  Singleton rows; they only
     * add to the normal-equation diagonal, so the cost is negligible.  This
     * is the term that resolves an undetermined fragment gauge toward the
     * spiral instead of toward zero-target stacking, and it removes the
     * unanchored-component nullspace outright.
     */
    if (p->prior_u != NULL && o->lambda_prior > 0.0) {
        for (size_t i = 0; i < p->nsamples; i++) {
            if (!isfinite(p->prior_u[i])) continue;
            MonotoneQpRow *row = &rows[nr++];
            row->first = nk;
            row->count = 1;
            row->target = p->prior_u[i];
            row->weight = o->lambda_prior;
            row->kind = ATLAS_STRIP_ROW_PRIOR;
            row->owner = (int32_t)i;
            coeff[nk].var = (int32_t)i;
            coeff[nk].value = 1.0;
            nk++;
        }
    }

    if (nr != nrows || nk != ncoeff || nb != nedge) {
        fprintf(stderr,
                "[atlas_strip] build: assembly mismatch rows %zu/%zu "
                "coeff %zu/%zu bounds %zu/%zu\n",
                nr, nrows, nk, ncoeff, nb, nedge);
        return -1;
    }

    memset(system, 0, sizeof *system);
    system->nsamples = p->nsamples;
    system->ncross_sections = p->ncross_sections;
    system->qp.nvar = p->nsamples + p->ncross_sections;
    system->qp.rows = rows;
    system->qp.nrows = nrows;
    system->qp.coeff = coeff;
    system->qp.ncoeff = ncoeff;
    system->qp.bounds = bounds;
    system->qp.nbounds = nedge;
    system->qp.anchors = p->anchors;
    system->qp.nanchors = p->nanchors;
    return 0;
}

void AtlasStrip_initialize_intercepts(const AtlasStripProblem *p,
                                      const AtlasStripOptions *input_opts,
                                      double *x)
{
    if (p == NULL || x == NULL) return;
    AtlasStripOptions defaults;
    AtlasStripOptions_default(&defaults);
    const AtlasStripOptions *o = input_opts != NULL ? input_opts : &defaults;

    for (size_t c = 0; c < p->ncross_sections; c++) {
        const AtlasStripCrossSection *cs = &p->cross_sections[c];
        double sum = 0.0, weight = 0.0;
        for (int32_t j = 0; j < cs->count; j++) {
            const AtlasStripMember *m =
                &p->members[cs->first + (size_t)j];
            double w = m->membership * m->base_weight;
            if (w <= o->membership_floor) continue;
            double value = AtlasStrip_member_value(m, x);
            if (o->mode == ATLAS_STRIP_FINAL)
                value -= dot3(cs->tangent, m->p);
            sum += w * value;
            weight += w;
        }
        x[p->nsamples + c] = weight > 0.0 ? sum / weight : 0.0;
    }
}

void AtlasStrip_measure(const AtlasStripProblem *p,
                        const AtlasStripOptions *input_opts,
                        const double *x,
                        AtlasStripMetrics *m)
{
    if (m == NULL) return;
    memset(m, 0, sizeof *m);
    if (p == NULL || x == NULL) return;
    AtlasStripOptions defaults;
    AtlasStripOptions_default(&defaults);
    const AtlasStripOptions *o = input_opts != NULL ? input_opts : &defaults;

    double length_r2 = 0.0, align_r2 = 0.0, speed_r2 = 0.0;
    double continuation_r2 = 0.0;
    size_t nlength = 0, nalign = 0, nspeed = 0, ncontinuation = 0;
    m->min_monotone_ratio = DBL_MAX;

    for (size_t c = 0; c < p->ncross_sections; c++) {
        const AtlasStripCrossSection *cs = &p->cross_sections[c];
        double denom = 0.0, derivative = 0.0;
        for (int32_t j = 0; j < cs->count; j++) {
            const AtlasStripMember *member =
                &p->members[cs->first + (size_t)j];
            double mw = member->membership * member->base_weight;
            if (mw <= o->membership_floor) continue;
            double w = member->dual_width * mw;
            double d = (x[member->deriv_hi] - x[member->deriv_lo]) /
                       member->deriv_length;
            derivative += w * dot3(member->tangent, cs->tangent) * d;
            denom += w;
            m->active_members++;
        }
        if (denom > 0.0) {
            double residual = derivative / denom - 1.0;
            m->energy_length += 0.5 * o->lambda_length * cs->weight *
                                residual * residual;
            length_r2 += residual * residual;
            nlength++;
        }

        int32_t qvar = (int32_t)(p->nsamples + c);
        for (int32_t j = 0; j < cs->count; j++) {
            const AtlasStripMember *member =
                &p->members[cs->first + (size_t)j];
            double mw = member->membership * member->base_weight;
            if (mw <= o->membership_floor) continue;
            double target = o->mode == ATLAS_STRIP_FINAL
                          ? dot3(cs->tangent, member->p) : 0.0;
            double residual = AtlasStrip_member_value(member, x) -
                              x[qvar] - target;
            double weight = o->lambda_align * cs->weight * mw;
            m->energy_align += 0.5 * weight * residual * residual;
            align_r2 += residual * residual;
            nalign++;
            if (fabs(residual) > m->max_align_residual)
                m->max_align_residual = fabs(residual);
        }
    }

    for (size_t s = 0; s < p->nstrokes; s++) {
        const AtlasStripStroke *st = &p->strokes[s];
        for (int32_t i = 1; i < st->count; i++) {
            int32_t lo = st->first + i - 1, hi = st->first + i;
            double length = p->samples[hi].s - p->samples[lo].s;
            double speed = (x[hi] - x[lo]) / length;
            double error = speed - 1.0;
            double conf = fmin(p->samples[lo].confidence,
                               p->samples[hi].confidence);
            m->energy_local += 0.5 * o->lambda_local * conf * length *
                               error * error;
            speed_r2 += error * error;
            nspeed++;
            if (fabs(error) > m->max_local_speed_error)
                m->max_local_speed_error = fabs(error);
            if (speed < m->min_monotone_ratio)
                m->min_monotone_ratio = speed;
        }
    }
    for (size_t i = 0; i < p->ncontinuations; i++) {
        const AtlasStripContinuation *cn = &p->continuations[i];
        double mw = cn->membership * cn->weight;
        if (mw <= o->membership_floor) continue;
        double residual = x[cn->a] - x[cn->b] - cn->target;
        m->energy_continuation += 0.5 * o->lambda_continuation * mw *
                                  residual * residual;
        continuation_r2 += residual * residual;
        if (fabs(residual) > m->max_continuation_residual)
            m->max_continuation_residual = fabs(residual);
        m->active_continuations++;
        ncontinuation++;
    }
    if (p->prior_u != NULL && o->lambda_prior > 0.0) {
        double prior_r2 = 0.0;
        for (size_t i = 0; i < p->nsamples; i++) {
            if (!isfinite(p->prior_u[i])) continue;
            double residual = x[i] - p->prior_u[i];
            m->energy_prior += 0.5 * o->lambda_prior * residual * residual;
            prior_r2 += residual * residual;
            m->active_priors++;
        }
        m->rms_prior_row = m->active_priors
                         ? sqrt(prior_r2 / (double)m->active_priors) : 0.0;
    }
    if (m->min_monotone_ratio == DBL_MAX) m->min_monotone_ratio = 0.0;
    m->rms_length_row = nlength ? sqrt(length_r2 / (double)nlength) : 0.0;
    m->rms_align_row = nalign ? sqrt(align_r2 / (double)nalign) : 0.0;
    m->rms_local_speed_error =
        nspeed ? sqrt(speed_r2 / (double)nspeed) : 0.0;
    m->rms_continuation = ncontinuation
                        ? sqrt(continuation_r2 / (double)ncontinuation) : 0.0;
}

void AtlasStrip_project_monotone(Arena_T arena,
                                 const AtlasStripProblem *p,
                                 const AtlasStripOptions *input_opts,
                                 double *x)
{
    if (arena == NULL || p == NULL || x == NULL) return;
    AtlasStripOptions defaults;
    AtlasStripOptions_default(&defaults);
    const AtlasStripOptions *o = input_opts != NULL ? input_opts : &defaults;
    double beta = o->monotone_fraction;

    Arena_Mark mark = Arena_save(arena);
    uint8_t *anchored = (uint8_t *)ARENA_CALLOC(
        arena, p->nsamples, sizeof(uint8_t));
    for (size_t a = 0; a < p->nanchors; a++)
        if ((size_t)p->anchors[a].var < p->nsamples)
            anchored[p->anchors[a].var] = 1;
    size_t longest = 2;
    for (size_t s = 0; s < p->nstrokes; s++)
        if ((size_t)p->strokes[s].count > longest)
            longest = (size_t)p->strokes[s].count;
    double *value = (double *)ARENA_ALLOC(arena, longest * sizeof(double));
    double *weight = (double *)ARENA_ALLOC(arena, longest * sizeof(double));
    int32_t *count = (int32_t *)ARENA_ALLOC(arena, longest * sizeof(int32_t));

    for (size_t s = 0; s < p->nstrokes; s++) {
        const AtlasStripStroke *st = &p->strokes[s];
        size_t n = (size_t)st->count;
        if (n < 2) continue;
        size_t top = 0;
        for (size_t i = 0; i < n; i++) {
            size_t v = (size_t)st->first + i;
            value[top] = x[v] - beta * p->samples[v].s;
            weight[top] = anchored[v] ? 1e12 : 1.0;
            count[top] = 1;
            top++;
            while (top > 1 && value[top - 1] < value[top - 2]) {
                double w = weight[top - 2] + weight[top - 1];
                value[top - 2] = (weight[top - 2] * value[top - 2] +
                                  weight[top - 1] * value[top - 1]) / w;
                weight[top - 2] = w;
                count[top - 2] += count[top - 1];
                top--;
            }
        }
        size_t i = 0;
        for (size_t b = 0; b < top; b++)
            for (int32_t k = 0; k < count[b]; k++, i++) {
                size_t v = (size_t)st->first + i;
                x[v] = value[b] + beta * p->samples[v].s;
            }
    }
    Arena_restore(arena, mark);
}

static double robust_member_residual(const AtlasStripProblem *p,
                                     const AtlasStripOptions *o,
                                     size_t cross,
                                     const AtlasStripMember *member,
                                     const double *x)
{
    const AtlasStripCrossSection *cs = &p->cross_sections[cross];
    double target = o->mode == ATLAS_STRIP_FINAL
                  ? dot3(cs->tangent, member->p) : 0.0;
    return AtlasStrip_member_value(member, x) -
           x[p->nsamples + cross] - target;
}

static double derived_likelihood_sigma(const AtlasStripProblem *p)
{
    double longest = 0.0;
    for (size_t s = 0; s < p->nstrokes; s++) {
        const AtlasStripStroke *stroke = &p->strokes[s];
        if (stroke->count < 2) continue;
        double span = p->samples[(size_t)stroke->first +
                                 (size_t)stroke->count - 1].s -
                      p->samples[stroke->first].s;
        if (span > longest) longest = span;
    }
    return fmax(longest / 30.0, 1e-3);
}

static void robust_trace_append(AtlasStripRobustTrace *trace,
                                const AtlasStripRobustTraceEntry *entry)
{
    if (trace == NULL || trace->entry == NULL ||
        trace->count >= trace->capacity) return;
    trace->entry[trace->count++] = *entry;
}

int AtlasStrip_solve_robust(Arena_T arena,
                            const AtlasStripProblem *p,
                            const AtlasStripOptions *input_strip_opts,
                            const AtlasStripRobustOptions *input_robust_opts,
                            const MonotoneQpOptions *input_qp_opts,
                            double *x,
                            double *out_membership,
                            AtlasStripRobustTrace *trace,
                            AtlasStripRobustStats *stats)
{
    if (arena == NULL || p == NULL || x == NULL) return -1;

    AtlasStripOptions strip_defaults;
    AtlasStripOptions_default(&strip_defaults);
    const AtlasStripOptions *so = input_strip_opts != NULL
                                ? input_strip_opts : &strip_defaults;
    AtlasStripRobustOptions robust_defaults;
    AtlasStripRobustOptions_default(&robust_defaults);
    const AtlasStripRobustOptions *ro = input_robust_opts != NULL
                                     ? input_robust_opts : &robust_defaults;
    MonotoneQpOptions qp_defaults;
    MonotoneQpOptions_default(&qp_defaults);
    const MonotoneQpOptions *qo = input_qp_opts != NULL
                                ? input_qp_opts : &qp_defaults;

    if (validate_strip_problem(p, so) != 0 ||
        so->mode != ATLAS_STRIP_RELAXED ||
        ro->l1_iterations < 0 || ro->likelihood_iterations < 0 ||
        ro->l1_iterations > 1000 || ro->likelihood_iterations > 1000 ||
        !isfinite(ro->initial_likelihood) ||
        ro->initial_likelihood <= 0.0 || ro->initial_likelihood > 1.0 ||
        !isfinite(ro->l1_epsilon) || ro->l1_epsilon <= 0.0 ||
        !isfinite(ro->likelihood_sigma) ||
        !isfinite(ro->likelihood_floor) ||
        ro->likelihood_floor <= 0.0 || ro->likelihood_floor > 1.0 ||
        !isfinite(ro->convergence_tolerance) ||
        ro->convergence_tolerance < 0.0)
        return -1;

    AtlasStripRobustStats local_stats;
    memset(&local_stats, 0, sizeof local_stats);
    if (trace != NULL) trace->count = 0;

    Arena_Mark outer_mark = Arena_save(arena);
    AtlasStripMember *members = (AtlasStripMember *)ARENA_ALLOC(
        arena, p->nmembers * sizeof(AtlasStripMember));
    double *prior = (double *)ARENA_ALLOC(
        arena, p->nmembers * sizeof(double));
    double *previous_x = (double *)ARENA_ALLOC(
        arena, p->nsamples * sizeof(double));
    memcpy(members, p->members, p->nmembers * sizeof(AtlasStripMember));
    for (size_t i = 0; i < p->nmembers; i++) {
        prior[i] = p->members[i].membership;
        members[i].membership = prior[i] * ro->initial_likelihood;
    }

    AtlasStripProblem work = *p;
    work.members = members;
    AtlasStrip_initialize_intercepts(&work, so, x);

    double sigma = ro->likelihood_sigma > 0.0
                 ? ro->likelihood_sigma : derived_likelihood_sigma(p);
    local_stats.likelihood_sigma = sigma;

    int total_rounds = 1 + ro->l1_iterations + ro->likelihood_iterations;
    int rc = 0;
    for (int round = 0; round < total_rounds; round++) {
        int phase = ATLAS_STRIP_ROBUST_INITIAL;
        if (round > 0 && round <= ro->l1_iterations)
            phase = ATLAS_STRIP_ROBUST_L1;
        else if (round > ro->l1_iterations)
            phase = ATLAS_STRIP_ROBUST_LIKELIHOOD;

        double max_membership_change = 0.0;
        if (round > 0) {
            double inv2sigma2 = 1.0 / (2.0 * sigma * sigma);
            for (size_t c = 0; c < p->ncross_sections; c++) {
                const AtlasStripCrossSection *cs = &p->cross_sections[c];
                for (int32_t j = 0; j < cs->count; j++) {
                    size_t mi = cs->first + (size_t)j;
                    double old = members[mi].membership;
                    double residual =
                        robust_member_residual(&work, so, c, &members[mi], x);
                    double likelihood;
                    if (phase == ATLAS_STRIP_ROBUST_L1) {
                        likelihood = ro->initial_likelihood /
                                     sqrt(residual * residual +
                                          ro->l1_epsilon * ro->l1_epsilon);
                        if (likelihood > 1.0) likelihood = 1.0;
                    } else {
                        likelihood = exp(-residual * residual * inv2sigma2);
                        if (likelihood < ro->likelihood_floor)
                            likelihood = ro->likelihood_floor;
                    }
                    members[mi].membership = prior[mi] * likelihood;
                    double change = fabs(members[mi].membership - old);
                    if (change > max_membership_change)
                        max_membership_change = change;
                }
            }
        }

        memcpy(previous_x, x, p->nsamples * sizeof(double));
        Arena_Mark round_mark = Arena_save(arena);
        AtlasStripSystem system;
        if (AtlasStrip_build(arena, &work, so, &system) != 0) {
            Arena_restore(arena, round_mark);
            rc = -1;
            break;
        }
        MonotoneQpStats qp_stats;
        memset(&qp_stats, 0, sizeof qp_stats);
        int qp_rc = 0;
        if (ro->use_active_set) {
            qp_rc = MonotoneQp_solve(arena, &system.qp, qo, x,
                                     NULL, &qp_stats);
        } else {
            double lr = 0.0;
            qp_rc = MonotoneQp_solve_ls(arena, &system.qp, x, &lr);
            if (qp_rc == 0)
                AtlasStrip_project_monotone(arena, &work, so, x);
            qp_stats.objective_final = MonotoneQp_objective(&system.qp, x);
            qp_stats.max_reduced_linear_residual = lr;
        }

        AtlasStripRobustTraceEntry te;
        memset(&te, 0, sizeof te);
        te.iteration = round;
        te.phase = phase;
        te.qp_rc = qp_rc;
        te.objective = qp_stats.objective_final;
        te.max_membership_change = max_membership_change;
        te.membership_min = DBL_MAX;
        double residual2 = 0.0;
        size_t active_count = 0;
        for (size_t i = 0; i < p->nsamples; i++) {
            double change = fabs(x[i] - previous_x[i]);
            if (change > te.max_u_change) te.max_u_change = change;
        }
        for (size_t c = 0; c < p->ncross_sections; c++) {
            const AtlasStripCrossSection *cs = &p->cross_sections[c];
            for (int32_t j = 0; j < cs->count; j++) {
                size_t mi = cs->first + (size_t)j;
                if (prior[mi] <= 0.0) continue;
                double residual =
                    robust_member_residual(&work, so, c, &members[mi], x);
                residual2 += residual * residual;
                if (fabs(residual) > te.max_residual)
                    te.max_residual = fabs(residual);
                double normalized = members[mi].membership / prior[mi];
                if (normalized < te.membership_min)
                    te.membership_min = normalized;
                if (normalized > te.membership_max)
                    te.membership_max = normalized;
                te.membership_mean += normalized;
                if (normalized < 0.1) te.downweighted_members++;
                active_count++;
            }
        }
        if (active_count > 0) {
            te.residual_rms = sqrt(residual2 / (double)active_count);
            te.membership_mean /= (double)active_count;
        } else {
            te.membership_min = 0.0;
        }
        robust_trace_append(trace, &te);
        if (trace != NULL && trace->iteration_fn != NULL &&
            trace->iteration_fn(trace->iteration_context, &work, x,
                                members, &te) != 0) {
            Arena_restore(arena, round_mark);
            rc = -1;
            break;
        }
        local_stats.total_qp_solves++;
        if (phase == ATLAS_STRIP_ROBUST_INITIAL)
            local_stats.initial_solves++;
        else if (phase == ATLAS_STRIP_ROBUST_L1)
            local_stats.l1_solves++;
        else
            local_stats.likelihood_solves++;
        local_stats.max_u_change = te.max_u_change;
        local_stats.max_membership_change = te.max_membership_change;
        local_stats.downweighted_members = te.downweighted_members;
        local_stats.final_qp = qp_stats;
        Arena_restore(arena, round_mark);

        if (qp_rc != 0) {
            rc = qp_rc;
            break;
        }
        if (phase == ATLAS_STRIP_ROBUST_LIKELIHOOD &&
            te.max_u_change <= ro->convergence_tolerance &&
            te.max_membership_change <= ro->convergence_tolerance)
            break;
    }

    if (out_membership != NULL)
        for (size_t i = 0; i < p->nmembers; i++)
            out_membership[i] = members[i].membership;
    if (stats != NULL) *stats = local_stats;
    Arena_restore(arena, outer_mark);
    return rc;
}
