#include "atlas_ribbon_texture.h"

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/raw_sample.h"

typedef struct {
    size_t sample;
    int32_t column;
} ArtSampleRef;

typedef struct {
    int32_t column;
    double value;
} ArtSignalPoint;

static int art_compare_ref(const void *pa, const void *pb)
{
    const ArtSampleRef *a = (const ArtSampleRef *)pa;
    const ArtSampleRef *b = (const ArtSampleRef *)pb;
    if (a->column != b->column) return a->column < b->column ? -1 : 1;
    return a->sample < b->sample ? -1 : (a->sample > b->sample ? 1 : 0);
}

static int art_compare_double(const void *pa, const void *pb)
{
    double a = *(const double *)pa, b = *(const double *)pb;
    return a < b ? -1 : (a > b ? 1 : 0);
}

void AtlasRibbonTextureOptions_default(AtlasRibbonTextureOptions *o)
{
    assert(o);
    memset(o, 0, sizeof *o);
    o->search_columns = 8;
    o->min_samples = 24;
    o->min_correlation = 0.20;
    o->min_margin = 0.02;
    o->normal_range = 2.0;
    o->normal_steps = 5;
}

static void art_world_point(const AtlasRibbonObservationSet *set,
                            const AtlasRibbonLayerSample *sample,
                            float p[3], float normal[3])
{
    double tangent[3], n[3], nl = 0.0;
    for (int k = 0; k < 3; k++) {
        double xyz = set->axis_point[k] + sample->v * set->axis[k] +
                     sample->p[0] * set->basis0[k] +
                     sample->p[1] * set->basis1[k];
        p[k] = (float)xyz;
        tangent[k] = sample->tangent[0] * set->basis0[k] +
                     sample->tangent[1] * set->basis1[k];
    }
    n[0] = set->axis[1] * tangent[2] - set->axis[2] * tangent[1];
    n[1] = set->axis[2] * tangent[0] - set->axis[0] * tangent[2];
    n[2] = set->axis[0] * tangent[1] - set->axis[1] * tangent[0];
    for (int k = 0; k < 3; k++) nl += n[k] * n[k];
    nl = sqrt(nl);
    for (int k = 0; k < 3; k++) normal[k] =
        nl > 1.0e-9 ? (float)(n[k] / nl) : 0.0f;
}

static int art_correlation(const ArtSignalPoint *a, size_t na,
                           const ArtSignalPoint *b, size_t nb,
                           int lag, size_t minimum,
                           double *out_correlation, uint32_t *out_count)
{
    size_t ia = 0, ib = 0, n = 0;
    double sx = 0.0, sy = 0.0, sxx = 0.0, syy = 0.0, sxy = 0.0;
    while (ia < na && ib < nb) {
        int64_t ca = a[ia].column;
        int64_t cb = (int64_t)b[ib].column + (int64_t)lag;
        if (ca < cb) {
            ia++;
        } else if (ca > cb) {
            ib++;
        } else {
            double x = a[ia].value, y = b[ib].value;
            sx += x; sy += y;
            sxx += x * x; syy += y * y; sxy += x * y;
            n++; ia++; ib++;
        }
    }
    if (n < minimum || n > UINT32_MAX) return -1;
    double nn = (double)n;
    double vx = nn * sxx - sx * sx;
    double vy = nn * syy - sy * sy;
    if (!(vx > 1.0e-8) || !(vy > 1.0e-8)) return -1;
    double corr = (nn * sxy - sx * sy) / sqrt(vx * vy);
    if (!isfinite(corr)) return -1;
    if (corr < -1.0) corr = -1.0;
    if (corr > 1.0) corr = 1.0;
    *out_correlation = corr;
    *out_count = (uint32_t)n;
    return 0;
}

static int art_best_lag(const ArtSignalPoint *a, size_t na,
                        const ArtSignalPoint *b, size_t nb,
                        const AtlasRibbonTextureOptions *o,
                        int *out_lag, uint32_t *out_samples,
                        double *out_best, double *out_zero,
                        double *out_margin)
{
    int search = o->search_columns;
    size_t nscore = (size_t)(2 * search + 1);
    double *score = (double *)malloc(nscore * sizeof(*score));
    uint32_t *count = (uint32_t *)malloc(nscore * sizeof(*count));
    if (score == NULL || count == NULL) {
        free(score); free(count); return -1;
    }
    for (size_t i = 0; i < nscore; i++) {
        score[i] = -DBL_MAX;
        count[i] = 0;
    }
    int best_lag = 0;
    double best = -DBL_MAX;
    uint32_t best_count = 0;
    for (int lag = -search; lag <= search; lag++) {
        size_t at = (size_t)(lag + search);
        if (art_correlation(a, na, b, nb, lag, o->min_samples,
                            &score[at], &count[at]) != 0)
            continue;
        if (score[at] > best + 1.0e-12 ||
            (fabs(score[at] - best) <= 1.0e-12 && abs(lag) < abs(best_lag))) {
            best = score[at];
            best_lag = lag;
            best_count = count[at];
        }
    }
    if (!(best > -DBL_MAX)) {
        free(score); free(count); return -1;
    }
    double second = -DBL_MAX;
    for (int lag = -search; lag <= search; lag++) {
        if (abs(lag - best_lag) <= 1) continue;
        double value = score[(size_t)(lag + search)];
        if (value > second) second = value;
    }
    *out_lag = best_lag;
    *out_samples = best_count;
    *out_best = best;
    *out_zero = score[(size_t)search] > -DBL_MAX
              ? score[(size_t)search] : NAN;
    *out_margin = second > -DBL_MAX ? best - second : 0.0;
    free(score); free(count);
    return 0;
}

int AtlasRibbonTexture_measure(
    Arena_T arena,
    const PieceSet *ps,
    const AtlasSolution *solution,
    const AtlasRibbonObservationSet *set,
    const AtlasRibbonFitResult *fit,
    const char *raw_dir,
    const AtlasRibbonTextureOptions *opts,
    AtlasRibbonTextureResult *out)
{
    if (arena == NULL || ps == NULL || solution == NULL || set == NULL ||
        fit == NULL || raw_dir == NULL || opts == NULL || out == NULL ||
        solution->ncharts == 0 || set->nrows == 0 ||
        set->nrows > SIZE_MAX / solution->ncharts ||
        fit->row_chart_shift == NULL || fit->nrows != set->nrows ||
        opts->search_columns < 0 || opts->search_columns > 128 ||
        opts->min_samples < 3 || !isfinite(opts->min_correlation) ||
        !isfinite(opts->min_margin) || opts->normal_steps < 1)
        return -1;
    memset(out, 0, sizeof *out);
    size_t nfield = set->nrows * solution->ncharts;
    size_t nsample = set->nlayer_samples;
    CubeTable *ct = (CubeTable *)ARENA_ALLOC(arena, sizeof(*ct));
    if (cubetable_init(ct, arena, raw_dir, 128, ps->verts, ps->nv,
                       opts->normal_range + 2.0) != 0)
        return -1;

    double *value = (double *)ARENA_ALLOC(arena, nsample * sizeof(*value));
    size_t *field_count = (size_t *)ARENA_CALLOC(
        arena, nfield, sizeof(*field_count));
    for (size_t i = 0; i < nsample; i++) {
        const AtlasRibbonLayerSample *s = &set->layer_sample[i];
        if (s->row < 0 || (size_t)s->row >= set->nrows ||
            s->chart < 0 || (size_t)s->chart >= solution->ncharts)
            return -1;
        float p[3], normal[3];
        art_world_point(set, s, p, normal);
        value[i] = sample_vertex(ct, p, normal, opts->normal_range,
                                 opts->normal_steps);
        if (value[i] >= 0.0) out->sampled_layer_points++;
        else out->missing_layer_points++;
        field_count[(size_t)s->row * solution->ncharts +
                    (size_t)s->chart]++;
    }
    out->cubes_loaded = ct->n_loaded;
    out->cubes_missing = ct->n_missing;

    size_t *offset = (size_t *)ARENA_ALLOC(
        arena, (nfield + 1) * sizeof(*offset));
    size_t *cursor = (size_t *)ARENA_ALLOC(
        arena, nfield * sizeof(*cursor));
    offset[0] = 0;
    for (size_t field = 0; field < nfield; field++)
        offset[field + 1] = offset[field] + field_count[field];
    if (offset[nfield] != nsample) return -1;
    memcpy(cursor, offset, nfield * sizeof(*cursor));
    ArtSampleRef *ref = (ArtSampleRef *)ARENA_ALLOC(
        arena, nsample * sizeof(*ref));
    for (size_t i = 0; i < nsample; i++) {
        const AtlasRibbonLayerSample *s = &set->layer_sample[i];
        size_t field = (size_t)s->row * solution->ncharts + (size_t)s->chart;
        size_t at = cursor[field]++;
        ref[at].sample = i;
        ref[at].column = s->source_column;
    }

    size_t *signal_start = (size_t *)ARENA_ALLOC(
        arena, nfield * sizeof(*signal_start));
    size_t *signal_count = (size_t *)ARENA_CALLOC(
        arena, nfield, sizeof(*signal_count));
    ArtSignalPoint *signal = (ArtSignalPoint *)ARENA_ALLOC(
        arena, nsample * sizeof(*signal));
    size_t signal_cursor = 0;
    for (size_t field = 0; field < nfield; field++) {
        size_t first = offset[field], last = offset[field + 1];
        signal_start[field] = signal_cursor;
        if (last <= first) continue;
        qsort(&ref[first], last - first, sizeof(*ref), art_compare_ref);
        for (size_t i = first; i < last;) {
            size_t j = i + 1, count = 0;
            double sum = 0.0;
            while (j < last && ref[j].column == ref[i].column) j++;
            for (size_t k = i; k < j; k++) {
                double sample_value = value[ref[k].sample];
                if (sample_value < 0.0) continue;
                sum += sample_value; count++;
            }
            if (count > 0) {
                signal[signal_cursor].column = ref[i].column;
                signal[signal_cursor].value = sum / (double)count;
                signal_cursor++;
                signal_count[field]++;
            }
            i = j;
        }
    }

    size_t pair_capacity = 0;
    for (size_t row = 0; row + 1 < set->nrows; row++)
    for (size_t chart = 0; chart < solution->ncharts; chart++) {
        size_t a = row * solution->ncharts + chart;
        size_t b = a + solution->ncharts;
        if (signal_count[a] >= opts->min_samples &&
            signal_count[b] >= opts->min_samples)
            pair_capacity++;
    }
    out->pair = (AtlasRibbonTexturePair *)ARENA_ALLOC(
        arena, (pair_capacity ? pair_capacity : 1) * sizeof(*out->pair));
    double *accepted_residual = (double *)ARENA_ALLOC(
        arena, (pair_capacity ? pair_capacity : 1) * sizeof(*accepted_residual));
    size_t nresidual = 0;
    double corr_sum = 0.0, margin_sum = 0.0;
    double current2 = 0.0, desired2 = 0.0, residual2 = 0.0;
    for (size_t row = 0; row + 1 < set->nrows; row++)
    for (size_t chart = 0; chart < solution->ncharts; chart++) {
        size_t fa = row * solution->ncharts + chart;
        size_t fb = fa + solution->ncharts;
        if (signal_count[fa] < opts->min_samples ||
            signal_count[fb] < opts->min_samples)
            continue;
        AtlasRibbonTexturePair *pair = &out->pair[out->npairs++];
        memset(pair, 0, sizeof *pair);
        pair->row0 = (int32_t)row;
        pair->row1 = (int32_t)(row + 1);
        pair->chart = (int32_t)chart;
        pair->group = solution->chart[chart].group;
        pair->winding = solution->chart[chart].winding;
        if (art_best_lag(&signal[signal_start[fa]], signal_count[fa],
                         &signal[signal_start[fb]], signal_count[fb], opts,
                         &pair->lag_columns, &pair->samples,
                         &pair->correlation, &pair->zero_correlation,
                         &pair->margin) != 0) {
            pair->correlation = pair->zero_correlation = pair->margin = NAN;
            continue;
        }
        pair->desired_shift_delta =
            (double)pair->lag_columns * set->observation_du;
        pair->current_shift_delta = fit->row_chart_shift[fb] -
                                    fit->row_chart_shift[fa];
        pair->residual = pair->current_shift_delta -
                         pair->desired_shift_delta;
        pair->accepted = pair->correlation >= opts->min_correlation &&
                         pair->margin >= opts->min_margin;
        if (!pair->accepted) continue;
        out->accepted_pairs++;
        corr_sum += pair->correlation;
        margin_sum += pair->margin;
        current2 += pair->current_shift_delta * pair->current_shift_delta;
        desired2 += pair->desired_shift_delta * pair->desired_shift_delta;
        residual2 += pair->residual * pair->residual;
        double magnitude = fabs(pair->residual);
        accepted_residual[nresidual++] = magnitude;
        if (magnitude > out->residual_max) out->residual_max = magnitude;
    }
    if (out->npairs != pair_capacity) return -1;
    if (out->accepted_pairs > 0) {
        double den = (double)out->accepted_pairs;
        out->accepted_correlation_mean = corr_sum / den;
        out->accepted_margin_mean = margin_sum / den;
        out->current_delta_rms = sqrt(current2 / den);
        out->desired_delta_rms = sqrt(desired2 / den);
        out->residual_rms = sqrt(residual2 / den);
        qsort(accepted_residual, nresidual, sizeof(*accepted_residual),
              art_compare_double);
        size_t p95 = (size_t)floor(0.95 * (double)(nresidual - 1));
        out->residual_p95 = accepted_residual[p95];
    }
    return 0;
}

int AtlasRibbonTexture_refresh(
    Arena_T scratch,
    const AtlasSolution *solution,
    const AtlasRibbonObservationSet *set,
    const AtlasRibbonFitResult *fit,
    AtlasRibbonTextureResult *out)
{
    if (scratch == NULL || solution == NULL || set == NULL || fit == NULL ||
        out == NULL || solution->ncharts == 0 || set->nrows == 0 ||
        set->nrows > SIZE_MAX / solution->ncharts ||
        fit->nrows != set->nrows || fit->row_chart_shift == NULL ||
        (out->npairs > 0 && out->pair == NULL))
        return -1;
    size_t ncharts = solution->ncharts;
    double *accepted_residual = (double *)ARENA_ALLOC(
        scratch, (out->npairs ? out->npairs : 1) *
                 sizeof(*accepted_residual));
    size_t nresidual = 0;
    double corr_sum = 0.0, margin_sum = 0.0;
    double current2 = 0.0, desired2 = 0.0, residual2 = 0.0;
    out->accepted_pairs = 0;
    out->accepted_correlation_mean = 0.0;
    out->accepted_margin_mean = 0.0;
    out->current_delta_rms = 0.0;
    out->desired_delta_rms = 0.0;
    out->residual_rms = 0.0;
    out->residual_p95 = 0.0;
    out->residual_max = 0.0;
    for (size_t i = 0; i < out->npairs; i++) {
        AtlasRibbonTexturePair *pair = &out->pair[i];
        if (pair->row0 < 0 || pair->row1 != pair->row0 + 1 ||
            (size_t)pair->row1 >= set->nrows || pair->chart < 0 ||
            (size_t)pair->chart >= ncharts ||
            !isfinite(pair->desired_shift_delta))
            return -1;
        size_t fa = (size_t)pair->row0 * ncharts + (size_t)pair->chart;
        size_t fb = (size_t)pair->row1 * ncharts + (size_t)pair->chart;
        pair->current_shift_delta =
            fit->row_chart_shift[fb] - fit->row_chart_shift[fa];
        pair->residual =
            pair->current_shift_delta - pair->desired_shift_delta;
        if (!pair->accepted) continue;
        if (!isfinite(pair->correlation) || !isfinite(pair->margin) ||
            !isfinite(pair->current_shift_delta) || !isfinite(pair->residual))
            return -1;
        out->accepted_pairs++;
        corr_sum += pair->correlation;
        margin_sum += pair->margin;
        current2 += pair->current_shift_delta * pair->current_shift_delta;
        desired2 += pair->desired_shift_delta * pair->desired_shift_delta;
        residual2 += pair->residual * pair->residual;
        double magnitude = fabs(pair->residual);
        accepted_residual[nresidual++] = magnitude;
        if (magnitude > out->residual_max) out->residual_max = magnitude;
    }
    if (out->accepted_pairs > 0) {
        double den = (double)out->accepted_pairs;
        out->accepted_correlation_mean = corr_sum / den;
        out->accepted_margin_mean = margin_sum / den;
        out->current_delta_rms = sqrt(current2 / den);
        out->desired_delta_rms = sqrt(desired2 / den);
        out->residual_rms = sqrt(residual2 / den);
        qsort(accepted_residual, nresidual, sizeof(*accepted_residual),
              art_compare_double);
        size_t p95 = (size_t)floor(0.95 * (double)(nresidual - 1));
        out->residual_p95 = accepted_residual[p95];
    }
    return 0;
}

int AtlasRibbonTexture_selftest(void)
{
    enum { N = 80, SHIFT = 3 };
    ArtSignalPoint a[N], b[N];
    AtlasRibbonTextureOptions o;
    AtlasRibbonTextureOptions_default(&o);
    o.search_columns = 8;
    o.min_samples = 32;
    for (int i = 0; i < N; i++) {
        a[i].column = i;
        b[i].column = i;
        a[i].value = sin(0.071 * i * i) + 0.37 * cos(0.31 * i);
    }
    for (int i = 0; i < N; i++) {
        int source = i - SHIFT;
        b[i].value = source >= 0
                   ? sin(0.071 * source * source) + 0.37 * cos(0.31 * source)
                   : 0.0;
    }
    int lag = 0;
    uint32_t samples = 0;
    double best = 0.0, zero = 0.0, margin = 0.0;
    int rc = art_best_lag(a, N, b, N, &o, &lag, &samples,
                          &best, &zero, &margin);
    int fail = rc != 0 || lag != -SHIFT || samples < (uint32_t)(N - SHIFT) ||
               best < 0.99 || margin <= 0.0;
    fprintf(stderr,
            "[atlas_ribbon_texture selftest] %s lag=%d best=%.4f margin=%.4f n=%u\n",
            fail ? "FAIL" : "PASS", lag, best, margin, samples);
    return fail;
}
