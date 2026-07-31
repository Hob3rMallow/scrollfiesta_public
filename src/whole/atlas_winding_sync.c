#include "atlas_winding_sync.h"

#include "../common/union_find.h"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef AWS_TWO_PI
#define AWS_TWO_PI 6.283185307179586476925286766559
#endif

typedef struct {
    int32_t axial_bin;
    int32_t phase_bin;
    int32_t island;
    int32_t raw_turn;
    double radius;
    double raw_winding;
} AwsFaceRecord;

typedef struct {
    int32_t axial_bin;
    int32_t phase_bin;
    int32_t island;
    int32_t raw_turn;
    size_t faces;
    double radius;
    double raw_winding;
} AwsStrand;

typedef struct {
    int32_t island0;
    int32_t island1;
    int32_t target;
    double residual;
    uint8_t seam;
} AwsObservation;

typedef struct {
    int32_t other;
    int32_t target;
    int32_t next;
    double weight;
} AwsAdjacency;

static int aws_compare_double(const void *pa, const void *pb)
{
    double a = *(const double *)pa, b = *(const double *)pb;
    return a < b ? -1 : (a > b ? 1 : 0);
}

static double aws_median(double *value, size_t count)
{
    if (count == 0) return NAN;
    qsort(value, count, sizeof(double), aws_compare_double);
    return count & 1u ? value[count / 2]
                      : 0.5 * (value[count / 2 - 1] + value[count / 2]);
}

static int aws_compare_face_record(const void *pa, const void *pb)
{
    const AwsFaceRecord *a = (const AwsFaceRecord *)pa;
    const AwsFaceRecord *b = (const AwsFaceRecord *)pb;
    if (a->axial_bin != b->axial_bin)
        return a->axial_bin < b->axial_bin ? -1 : 1;
    if (a->phase_bin != b->phase_bin)
        return a->phase_bin < b->phase_bin ? -1 : 1;
    if (a->island != b->island) return a->island < b->island ? -1 : 1;
    if (a->raw_turn != b->raw_turn)
        return a->raw_turn < b->raw_turn ? -1 : 1;
    return a->radius < b->radius ? -1 : (a->radius > b->radius ? 1 : 0);
}

static int aws_compare_strand_bin_radius(const void *pa, const void *pb)
{
    const AwsStrand *a = (const AwsStrand *)pa;
    const AwsStrand *b = (const AwsStrand *)pb;
    if (a->axial_bin != b->axial_bin)
        return a->axial_bin < b->axial_bin ? -1 : 1;
    if (a->phase_bin != b->phase_bin)
        return a->phase_bin < b->phase_bin ? -1 : 1;
    if (a->radius != b->radius) return a->radius < b->radius ? -1 : 1;
    if (a->island != b->island) return a->island < b->island ? -1 : 1;
    return a->raw_turn < b->raw_turn ? -1 :
           (a->raw_turn > b->raw_turn ? 1 : 0);
}

static int aws_compare_observation(const void *pa, const void *pb)
{
    const AwsObservation *a = (const AwsObservation *)pa;
    const AwsObservation *b = (const AwsObservation *)pb;
    if (a->island0 != b->island0)
        return a->island0 < b->island0 ? -1 : 1;
    if (a->island1 != b->island1)
        return a->island1 < b->island1 ? -1 : 1;
    if (a->target != b->target) return a->target < b->target ? -1 : 1;
    if (a->seam != b->seam) return a->seam ? -1 : 1;
    return a->residual < b->residual ? -1 :
           (a->residual > b->residual ? 1 : 0);
}

static double aws_face_raw_winding(const AtlasWindingSyncProblem *p,
                                   size_t face)
{
    double value = 0.0;
    for (int k = 0; k < 3; k++) {
        int32_t vertex = p->faces[face * 3 + (size_t)k];
        value += (double)p->sense * (double)p->phi[vertex] / AWS_TWO_PI;
    }
    return value / 3.0;
}

static double aws_face_phase(const AtlasWindingSyncProblem *p, size_t face)
{
    double sine = 0.0, cosine = 0.0;
    for (int k = 0; k < 3; k++) {
        int32_t vertex = p->faces[face * 3 + (size_t)k];
        double w = (double)p->sense * (double)p->phi[vertex];
        sine += sin(w);
        cosine += cos(w);
    }
    double phase = atan2(sine, cosine);
    return phase < 0.0 ? phase + AWS_TWO_PI : phase;
}

static int aws_append_observation(AwsObservation *observation,
                                  size_t capacity, size_t *count,
                                  int32_t island0, int32_t island1,
                                  int32_t target, double residual, int seam)
{
    if (island0 == island1 || *count >= capacity || !isfinite(residual))
        return island0 == island1 ? 0 : -1;
    if (island0 > island1) {
        int32_t t = island0; island0 = island1; island1 = t;
        target = -target;
    }
    observation[*count].island0 = island0;
    observation[*count].island1 = island1;
    observation[*count].target = target;
    observation[*count].residual = residual;
    observation[*count].seam = (uint8_t)(seam != 0);
    (*count)++;
    return 0;
}

static double aws_truncated_square(double value, double limit)
{
    double square = value * value;
    return square < limit ? square : limit;
}

int AtlasWindingSync_solve(
    Arena_T arena, const AtlasWindingSyncProblem *p,
    int32_t **out_face_island, size_t *out_nislands,
    int32_t **out_island_correction,
    double **out_island_prior, double **out_island_prior_mad,
    AtlasWindingSyncRelation **out_relation, size_t *out_nrelation,
    AtlasWindingSyncStats *stats)
{
    if (arena == NULL || p == NULL || p->vertices == NULL ||
        p->nvertices == 0 || p->faces == NULL || p->nfaces == 0 ||
        p->nfaces > (size_t)INT32_MAX || p->phi == NULL ||
        p->face_radius == NULL || p->face_axial == NULL ||
        p->face_chart == NULL || p->ncharts == 0 ||
        p->adjacency_face0 == NULL || p->adjacency_face1 == NULL ||
        p->nintrinsic_adjacency > p->nadjacency ||
        !isfinite(p->spiral_a) || !isfinite(p->spiral_b) ||
        !isfinite(p->pitch) || p->pitch <= 1e-8 ||
        (p->sense != -1 && p->sense != 1) ||
        !isfinite(p->axial_bin_spacing) || p->axial_bin_spacing <= 0.0 ||
        p->phase_bins < 8 || p->phase_bins > 4096 ||
        out_face_island == NULL || out_nislands == NULL ||
        out_island_correction == NULL || out_island_prior == NULL ||
        out_island_prior_mad == NULL || out_relation == NULL ||
        out_nrelation == NULL || stats == NULL)
        return -1;
    memset(stats, 0, sizeof(*stats));

    UnionFind intrinsic = UF_new(arena, (int32_t)p->nfaces);
    for (size_t i = 0; i < p->nintrinsic_adjacency; i++) {
        int32_t a = p->adjacency_face0[i], b = p->adjacency_face1[i];
        if (a < 0 || b < 0 || (size_t)a >= p->nfaces ||
            (size_t)b >= p->nfaces)
            return -1;
        if (p->face_chart[a] == p->face_chart[b]) uf_union(&intrinsic, a, b);
    }
    int32_t *root_island = (int32_t *)ARENA_ALLOC(
        arena, p->nfaces * sizeof(int32_t));
    int32_t *face_island = (int32_t *)ARENA_ALLOC(
        arena, p->nfaces * sizeof(int32_t));
    for (size_t f = 0; f < p->nfaces; f++) root_island[f] = -1;
    size_t nisland = 0;
    for (size_t f = 0; f < p->nfaces; f++) {
        int32_t chart = p->face_chart[f];
        if (chart < 0 || (size_t)chart >= p->ncharts) return -1;
        int32_t root = uf_find(&intrinsic, (int32_t)f);
        if (root_island[root] < 0) root_island[root] = (int32_t)nisland++;
        face_island[f] = root_island[root];
    }
    if (nisland == 0 || nisland > (size_t)INT32_MAX) return -1;

    size_t *island_faces = (size_t *)ARENA_CALLOC(
        arena, nisland, sizeof(size_t));
    double *face_q = (double *)ARENA_ALLOC(
        arena, p->nfaces * sizeof(double));
    double axial_min = DBL_MAX;
    for (size_t f = 0; f < p->nfaces; f++) {
        island_faces[face_island[f]]++;
        face_q[f] = aws_face_raw_winding(p, f);
        if (!isfinite(face_q[f]) || !isfinite(p->face_radius[f]) ||
            !isfinite(p->face_axial[f]))
            return -1;
        if (p->face_axial[f] < axial_min) axial_min = p->face_axial[f];
    }

    size_t *prior_offset = (size_t *)ARENA_CALLOC(
        arena, nisland + 1, sizeof(size_t));
    for (size_t f = 0; f < p->nfaces; f++)
        prior_offset[(size_t)face_island[f] + 1]++;
    for (size_t i = 0; i < nisland; i++)
        prior_offset[i + 1] += prior_offset[i];
    size_t *prior_cursor = (size_t *)ARENA_ALLOC(
        arena, nisland * sizeof(size_t));
    memcpy(prior_cursor, prior_offset, nisland * sizeof(size_t));
    double *prior_value = (double *)ARENA_ALLOC(
        arena, p->nfaces * sizeof(double));
    for (size_t f = 0; f < p->nfaces; f++) {
        int32_t island = face_island[f];
        prior_value[prior_cursor[island]++] =
            (p->face_radius[f] - p->spiral_a) / p->pitch - face_q[f];
    }
    double *prior = (double *)ARENA_ALLOC(arena, nisland * sizeof(double));
    double *prior_mad = (double *)ARENA_ALLOC(
        arena, nisland * sizeof(double));
    double *scratch = (double *)ARENA_ALLOC(
        arena, p->nfaces * sizeof(double));
    for (size_t i = 0; i < nisland; i++) {
        size_t count = prior_offset[i + 1] - prior_offset[i];
        if (count == 0) return -1;
        prior[i] = aws_median(&prior_value[prior_offset[i]], count);
        for (size_t j = 0; j < count; j++)
            scratch[j] = fabs(prior_value[prior_offset[i] + j] - prior[i]);
        prior_mad[i] = aws_median(scratch, count);
    }

    AwsFaceRecord *record = (AwsFaceRecord *)ARENA_ALLOC(
        arena, p->nfaces * sizeof(AwsFaceRecord));
    for (size_t f = 0; f < p->nfaces; f++) {
        double phase = aws_face_phase(p, f);
        int phase_bin = (int)floor(
            phase * (double)p->phase_bins / AWS_TWO_PI);
        if (phase_bin < 0) phase_bin = 0;
        if (phase_bin >= p->phase_bins) phase_bin = p->phase_bins - 1;
        double axial_bin_value =
            floor((p->face_axial[f] - axial_min) / p->axial_bin_spacing);
        if (axial_bin_value < (double)INT32_MIN ||
            axial_bin_value > (double)INT32_MAX)
            return -1;
        record[f].axial_bin = (int32_t)axial_bin_value;
        record[f].phase_bin = phase_bin;
        record[f].island = face_island[f];
        record[f].raw_turn = (int32_t)floor(face_q[f] + 1e-6);
        record[f].radius = p->face_radius[f];
        record[f].raw_winding = face_q[f];
    }
    qsort(record, p->nfaces, sizeof(AwsFaceRecord), aws_compare_face_record);
    AwsStrand *strand = (AwsStrand *)ARENA_ALLOC(
        arena, p->nfaces * sizeof(AwsStrand));
    size_t nstrand = 0;
    for (size_t first = 0; first < p->nfaces;) {
        size_t last = first + 1;
        while (last < p->nfaces &&
               record[last].axial_bin == record[first].axial_bin &&
               record[last].phase_bin == record[first].phase_bin &&
               record[last].island == record[first].island &&
               record[last].raw_turn == record[first].raw_turn)
            last++;
        AwsStrand *s = &strand[nstrand++];
        memset(s, 0, sizeof(*s));
        s->axial_bin = record[first].axial_bin;
        s->phase_bin = record[first].phase_bin;
        s->island = record[first].island;
        s->raw_turn = record[first].raw_turn;
        s->faces = last - first;
        for (size_t i = first; i < last; i++) {
            s->radius += record[i].radius;
            s->raw_winding += record[i].raw_winding;
        }
        s->radius /= (double)s->faces;
        s->raw_winding /= (double)s->faces;
        first = last;
    }
    qsort(strand, nstrand, sizeof(AwsStrand),
          aws_compare_strand_bin_radius);

    if (nstrand > (SIZE_MAX - p->nadjacency) / 3) return -1;
    size_t observation_capacity = 3 * nstrand + p->nadjacency + 1;
    AwsObservation *observation = (AwsObservation *)ARENA_ALLOC(
        arena, observation_capacity * sizeof(AwsObservation));
    size_t nobservation = 0, ncross = 0, nseam = 0, nbins = 0;
    for (size_t first = 0; first < nstrand;) {
        size_t last = first + 1;
        while (last < nstrand &&
               strand[last].axial_bin == strand[first].axial_bin &&
               strand[last].phase_bin == strand[first].phase_bin)
            last++;
        nbins++;
        for (size_t i = first; i < last; i++) {
            size_t stop = i + 4 < last ? i + 4 : last;
            for (size_t j = i + 1; j < stop; j++) {
                const AwsStrand *inner = &strand[i], *outer = &strand[j];
                if (inner->island == outer->island) continue;
                double gap = (outer->radius - inner->radius) / p->pitch;
                if (gap < 0.55 || gap > 4.5) continue;
                long target_long = lround(
                    gap - (outer->raw_winding - inner->raw_winding));
                if (target_long < INT32_MIN || target_long > INT32_MAX)
                    return -1;
                int32_t target = (int32_t)target_long;
                double residual = fabs(
                    gap - ((outer->raw_winding + (double)target) -
                           inner->raw_winding));
                if (residual > 0.35) continue;
                if (aws_append_observation(
                        observation, observation_capacity, &nobservation,
                        inner->island, outer->island, target, residual, 0) != 0)
                    return -1;
                ncross++;
            }
        }
        first = last;
    }
    for (size_t i = p->nintrinsic_adjacency; i < p->nadjacency; i++) {
        int32_t f0 = p->adjacency_face0[i], f1 = p->adjacency_face1[i];
        if (f0 < 0 || f1 < 0 || (size_t)f0 >= p->nfaces ||
            (size_t)f1 >= p->nfaces)
            return -1;
        int32_t a = face_island[f0], b = face_island[f1];
        if (a == b || p->face_chart[f0] != p->face_chart[f1]) continue;
        double correction = face_q[f0] - face_q[f1];
        long target_long = lround(correction);
        if (target_long < INT32_MIN || target_long > INT32_MAX) return -1;
        double residual = fabs(correction - (double)target_long);
        if (residual > 0.25) continue;
        if (aws_append_observation(
                observation, observation_capacity, &nobservation,
                a, b, (int32_t)target_long, residual, 1) != 0)
            return -1;
        nseam++;
    }
    if (nobservation == 0) return -1;
    qsort(observation, nobservation, sizeof(AwsObservation),
          aws_compare_observation);

    AtlasWindingSyncRelation *relation =
        (AtlasWindingSyncRelation *)ARENA_ALLOC(
            arena, nobservation * sizeof(AtlasWindingSyncRelation));
    size_t nrelation = 0;
    for (size_t first = 0; first < nobservation;) {
        size_t last = first + 1;
        while (last < nobservation &&
               observation[last].island0 == observation[first].island0 &&
               observation[last].island1 == observation[first].island1)
            last++;
        AtlasWindingSyncRelation *r = &relation[nrelation++];
        memset(r, 0, sizeof(*r));
        r->island0 = observation[first].island0;
        r->island1 = observation[first].island1;
        r->observations = last - first;
        for (size_t i = first; i < last; i++) {
            if (observation[i].seam) r->seam_observations++;
            else r->cross_section_observations++;
        }
        size_t mode_count = 0;
        int32_t mode = 0;
        for (size_t mfirst = first; mfirst < last;) {
            size_t mlast = mfirst + 1;
            while (mlast < last &&
                   observation[mlast].target == observation[mfirst].target)
                mlast++;
            if (mlast - mfirst > mode_count ||
                (mlast - mfirst == mode_count &&
                 abs(observation[mfirst].target) < abs(mode))) {
                mode_count = mlast - mfirst;
                mode = observation[mfirst].target;
            }
            mfirst = mlast;
        }
        r->target_turn_correction = mode;
        r->mode_agreement = (double)mode_count / (double)r->observations;
        size_t count = 0;
        for (size_t i = first; i < last; i++)
            if (observation[i].target == mode)
                scratch[count++] = observation[i].residual;
        r->residual_median = aws_median(scratch, count);
        int seam_strong = r->seam_observations > 0 &&
                          r->mode_agreement >= 0.80 &&
                          r->residual_median <= 0.20;
        int cross_strong = r->cross_section_observations >= 3 &&
                           r->mode_agreement >= 0.80 &&
                           r->residual_median <= 0.25;
        r->eligible = seam_strong || cross_strong;
        double seam_fraction = (double)r->seam_observations /
                               (double)r->observations;
        r->weight = (1.0 + 3.0 * seam_fraction) *
                    log1p((double)r->observations) * r->mode_agreement /
                    (1.0 + 4.0 * r->residual_median);
        first = last;
    }

    size_t neligible = 0;
    for (size_t i = 0; i < nrelation; i++)
        if (relation[i].eligible) neligible++;
    if (neligible == 0 || neligible > (SIZE_MAX - 1) / 2) return -1;
    int32_t *head = (int32_t *)ARENA_ALLOC(
        arena, nisland * sizeof(int32_t));
    for (size_t i = 0; i < nisland; i++) head[i] = -1;
    AwsAdjacency *adj = (AwsAdjacency *)ARENA_ALLOC(
        arena, (2 * neligible + 1) * sizeof(AwsAdjacency));
    size_t nadj = 0;
    uint8_t *has_relation = (uint8_t *)ARENA_CALLOC(
        arena, nisland, sizeof(uint8_t));
    for (size_t i = 0; i < nrelation; i++) {
        AtlasWindingSyncRelation *r = &relation[i];
        if (!r->eligible) continue;
        int32_t a = r->island0, b = r->island1;
        adj[nadj].other = b;
        adj[nadj].target = r->target_turn_correction;
        adj[nadj].weight = r->weight;
        adj[nadj].next = head[a]; head[a] = (int32_t)nadj++;
        adj[nadj].other = a;
        adj[nadj].target = -r->target_turn_correction;
        adj[nadj].weight = r->weight;
        adj[nadj].next = head[b]; head[b] = (int32_t)nadj++;
        has_relation[a] = has_relation[b] = 1;
    }

    double *prior_weight = (double *)ARENA_ALLOC(
        arena, nisland * sizeof(double));
    double *x = (double *)ARENA_ALLOC(arena, nisland * sizeof(double));
    for (size_t i = 0; i < nisland; i++) {
        double support = sqrt((double)island_faces[i]);
        prior_weight[i] = (1.0 + fmin(8.0, support / 8.0)) /
                          (1.0 + 2.0 * fmin(2.0, prior_mad[i]));
        x[i] = prior[i];
        if (has_relation[i]) stats->islands_with_relations++;
    }
    for (int outer = 0; outer < 8; outer++) {
        stats->robust_outer_iterations = outer + 1;
        for (int sweep = 0; sweep < 100; sweep++) {
            double maximum_change = 0.0;
            for (size_t i = 0; i < nisland; i++) {
                double numerator = prior_weight[i] * prior[i];
                double denominator = prior_weight[i];
                for (int32_t ei = head[i]; ei >= 0; ei = adj[ei].next) {
                    const AwsAdjacency *e = &adj[ei];
                    double predicted = x[e->other] - (double)e->target;
                    double residual = x[i] - predicted;
                    double robust = 1.0 / fmax(1.0, fabs(residual) / 0.75);
                    double weight = e->weight * robust;
                    numerator += weight * predicted;
                    denominator += weight;
                }
                if (denominator <= 0.0) continue;
                double next = numerator / denominator;
                double lower = prior[i] - 6.0, upper = prior[i] + 6.0;
                if (next < lower) next = lower;
                if (next > upper) next = upper;
                if (next > x[i] + 1.0) next = x[i] + 1.0;
                if (next < x[i] - 1.0) next = x[i] - 1.0;
                double change = fabs(next - x[i]);
                if (change > maximum_change) maximum_change = change;
                x[i] = next;
            }
            if (maximum_change < 1e-5) break;
        }
    }

    int32_t *correction = (int32_t *)ARENA_ALLOC(
        arena, nisland * sizeof(int32_t));
    for (size_t i = 0; i < nisland; i++) correction[i] = (int32_t)lround(x[i]);
    for (int sweep = 0; sweep < 20; sweep++) {
        size_t changed = 0;
        for (size_t i = 0; i < nisland; i++) {
            int32_t center = (int32_t)lround(prior[i]);
            int32_t minimum = center - 6, maximum = center + 6;
            int32_t best = correction[i];
            double best_cost = DBL_MAX;
            for (int32_t candidate = minimum; candidate <= maximum;
                 candidate++) {
                double cost = prior_weight[i] *
                    aws_truncated_square((double)candidate - prior[i], 4.0);
                for (int32_t ei = head[i]; ei >= 0; ei = adj[ei].next) {
                    const AwsAdjacency *e = &adj[ei];
                    double residual = (double)correction[e->other] -
                                      (double)candidate -
                                      (double)e->target;
                    cost += e->weight * aws_truncated_square(residual, 4.0);
                }
                if (cost < best_cost - 1e-10 ||
                    (fabs(cost - best_cost) <= 1e-10 &&
                     abs(candidate - center) < abs(best - center))) {
                    best_cost = cost;
                    best = candidate;
                }
            }
            if (best != correction[i]) {
                correction[i] = best;
                changed++;
            }
        }
        stats->discrete_sweeps = sweep + 1;
        if (changed == 0) break;
    }

    UnionFind relation_graph = UF_new(arena, (int32_t)nisland);
    double relation_residual2 = 0.0, prior_residual2 = 0.0;
    stats->minimum_correction = INT32_MAX;
    stats->maximum_correction = INT32_MIN;
    for (size_t i = 0; i < nisland; i++) {
        double residual = (double)correction[i] - prior[i];
        prior_residual2 += residual * residual;
        if (correction[i] != 0) stats->nonzero_corrections++;
        if (correction[i] < stats->minimum_correction)
            stats->minimum_correction = correction[i];
        if (correction[i] > stats->maximum_correction)
            stats->maximum_correction = correction[i];
    }
    for (size_t i = 0; i < nrelation; i++) {
        AtlasWindingSyncRelation *r = &relation[i];
        r->solved_turn_correction =
            correction[r->island1] - correction[r->island0];
        r->final_residual = r->solved_turn_correction -
                            r->target_turn_correction;
        if (!r->eligible) continue;
        uf_union(&relation_graph, r->island0, r->island1);
        relation_residual2 += (double)r->final_residual *
                              (double)r->final_residual;
        if (r->final_residual == 0) stats->satisfied_relations++;
    }
    for (size_t i = 0; i < nisland; i++)
        if (uf_find(&relation_graph, (int32_t)i) == (int32_t)i)
            stats->relation_graph_components++;

    stats->islands = nisland;
    stats->cylindrical_bins = nbins;
    stats->strand_samples = nstrand;
    stats->observations = nobservation;
    stats->cross_section_observations = ncross;
    stats->seam_observations = nseam;
    stats->relations = nrelation;
    stats->eligible_relations = neligible;
    stats->relation_satisfaction_fraction = neligible
        ? (double)stats->satisfied_relations / (double)neligible : 0.0;
    stats->prior_residual_rms =
        sqrt(prior_residual2 / (double)nisland);
    stats->final_relation_residual_rms = neligible
        ? sqrt(relation_residual2 / (double)neligible) : 0.0;

    *out_face_island = face_island;
    *out_nislands = nisland;
    *out_island_correction = correction;
    *out_island_prior = prior;
    *out_island_prior_mad = prior_mad;
    *out_relation = relation;
    *out_nrelation = nrelation;
    return 0;
}
