#include "atlas_boxcut_layout.h"

#include "monotone_qp.h"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    int32_t chart0;
    int32_t chart1;
    double target;
    double separation;
} AbclObservation;

typedef struct {
    size_t relation;
    size_t pairs;
    double score;
} AbclCandidate;

typedef struct {
    int32_t chart;
    double radius;
} AbclChartRank;

typedef struct {
    int32_t inner_chart;
    int32_t outer_chart;
    double lower;
} AbclCollision;

typedef struct {
    int64_t slab;
    int32_t chart;
    size_t radial_rank;
    double umin;
    double umax;
} AbclSlabRecord;

typedef struct {
    int32_t chart0;
    int32_t chart1;
    double weight;
    int seam;
} AbclCoherenceSource;

typedef struct {
    int32_t chart0;
    int32_t chart1;
    double left_upper;
    double right_lower;
} AbclDisjunctionSource;

static int abcl_compare_double(const void *pa, const void *pb)
{
    double a = *(const double *)pa, b = *(const double *)pb;
    return a < b ? -1 : (a > b ? 1 : 0);
}

static double abcl_median(double *value, size_t count)
{
    if (count == 0) return NAN;
    qsort(value, count, sizeof(double), abcl_compare_double);
    return count & 1u ? value[count / 2]
                      : 0.5 * (value[count / 2 - 1] + value[count / 2]);
}

static int abcl_compare_observation(const void *pa, const void *pb)
{
    const AbclObservation *a = (const AbclObservation *)pa;
    const AbclObservation *b = (const AbclObservation *)pb;
    if (a->chart0 != b->chart0) return a->chart0 < b->chart0 ? -1 : 1;
    if (a->chart1 != b->chart1) return a->chart1 < b->chart1 ? -1 : 1;
    return a->target < b->target ? -1 : (a->target > b->target ? 1 : 0);
}

static int abcl_compare_candidate(const void *pa, const void *pb)
{
    const AbclCandidate *a = (const AbclCandidate *)pa;
    const AbclCandidate *b = (const AbclCandidate *)pb;
    if (a->score != b->score) return a->score > b->score ? -1 : 1;
    if (a->pairs != b->pairs) return a->pairs > b->pairs ? -1 : 1;
    return a->relation < b->relation ? -1 :
           (a->relation > b->relation ? 1 : 0);
}

static int abcl_compare_chart_rank(const void *pa, const void *pb)
{
    const AbclChartRank *a = (const AbclChartRank *)pa;
    const AbclChartRank *b = (const AbclChartRank *)pb;
    if (a->radius != b->radius) return a->radius < b->radius ? -1 : 1;
    return a->chart < b->chart ? -1 : (a->chart > b->chart ? 1 : 0);
}

static int abcl_compare_collision(const void *pa, const void *pb)
{
    const AbclCollision *a = (const AbclCollision *)pa;
    const AbclCollision *b = (const AbclCollision *)pb;
    if (a->inner_chart != b->inner_chart)
        return a->inner_chart < b->inner_chart ? -1 : 1;
    if (a->outer_chart != b->outer_chart)
        return a->outer_chart < b->outer_chart ? -1 : 1;
    return a->lower < b->lower ? -1 : (a->lower > b->lower ? 1 : 0);
}

static int abcl_compare_slab_record(const void *pa, const void *pb)
{
    const AbclSlabRecord *a = (const AbclSlabRecord *)pa;
    const AbclSlabRecord *b = (const AbclSlabRecord *)pb;
    if (a->slab != b->slab) return a->slab < b->slab ? -1 : 1;
    if (a->radial_rank != b->radial_rank)
        return a->radial_rank < b->radial_rank ? -1 : 1;
    return a->chart < b->chart ? -1 : (a->chart > b->chart ? 1 : 0);
}

static int abcl_compare_u64(const void *pa, const void *pb)
{
    uint64_t a = *(const uint64_t *)pa;
    uint64_t b = *(const uint64_t *)pb;
    return a < b ? -1 : (a > b ? 1 : 0);
}

static int abcl_compare_coherence_source(const void *pa, const void *pb)
{
    const AbclCoherenceSource *a = (const AbclCoherenceSource *)pa;
    const AbclCoherenceSource *b = (const AbclCoherenceSource *)pb;
    if (a->chart0 != b->chart0) return a->chart0 < b->chart0 ? -1 : 1;
    if (a->chart1 != b->chart1) return a->chart1 < b->chart1 ? -1 : 1;
    if (a->seam != b->seam) return a->seam < b->seam ? -1 : 1;
    return a->weight < b->weight ? -1 : (a->weight > b->weight ? 1 : 0);
}

static int abcl_compare_disjunction_source(const void *pa, const void *pb)
{
    const AbclDisjunctionSource *a = (const AbclDisjunctionSource *)pa;
    const AbclDisjunctionSource *b = (const AbclDisjunctionSource *)pb;
    if (a->chart0 != b->chart0) return a->chart0 < b->chart0 ? -1 : 1;
    if (a->chart1 != b->chart1) return a->chart1 < b->chart1 ? -1 : 1;
    if (a->left_upper != b->left_upper)
        return a->left_upper < b->left_upper ? -1 : 1;
    return a->right_lower < b->right_lower ? -1 :
           (a->right_lower > b->right_lower ? 1 : 0);
}

static uint64_t abcl_face_pair_key(int32_t a, int32_t b)
{
    uint32_t lo = (uint32_t)(a < b ? a : b);
    uint32_t hi = (uint32_t)(a < b ? b : a);
    return ((uint64_t)lo << 32) | (uint64_t)hi;
}

static int abcl_has_u64(const uint64_t *value, size_t count, uint64_t key)
{
    size_t lo = 0, hi = count;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (value[mid] < key) lo = mid + 1;
        else hi = mid;
    }
    return lo < count && value[lo] == key;
}

static double abcl_face_mean_double(const AtlasBoxcutLayoutProblem *p,
                                    size_t face, const double *value)
{
    return (value[p->faces[face * 3]] +
            value[p->faces[face * 3 + 1]] +
            value[p->faces[face * 3 + 2]]) / 3.0;
}

static double abcl_face_mean_float(const AtlasBoxcutLayoutProblem *p,
                                   size_t face, const float *value)
{
    return ((double)value[p->faces[face * 3]] +
            (double)value[p->faces[face * 3 + 1]] +
            (double)value[p->faces[face * 3 + 2]]) / 3.0;
}

static void abcl_face_u_range(const AtlasBoxcutLayoutProblem *p,
                              size_t face, double *out_min, double *out_max)
{
    double minimum = DBL_MAX, maximum = -DBL_MAX;
    for (int k = 0; k < 3; k++) {
        double value = p->base_u[p->faces[face * 3 + (size_t)k]];
        if (value < minimum) minimum = value;
        if (value > maximum) maximum = value;
    }
    *out_min = minimum;
    *out_max = maximum;
}

static int abcl_triangle_u_slice(const AtlasBoxcutLayoutProblem *p,
                                 size_t face, double axial,
                                 double *out_min, double *out_max)
{
    double minimum = DBL_MAX, maximum = -DBL_MAX;
    int intersections = 0;
    double scale = fmax(1.0, fabs(axial));
    for (int edge = 0; edge < 3; edge++) {
        int32_t va = p->faces[face * 3 + (size_t)edge];
        int32_t vb = p->faces[face * 3 + (size_t)((edge + 1) % 3)];
        double aa = p->v[va], ab = p->v[vb];
        double ua = p->base_u[va], ub = p->base_u[vb];
        scale = fmax(scale, fmax(fabs(aa), fabs(ab)));
        double tolerance = 1.0e-12 * scale;
        double dv = ab - aa;
        if (fabs(dv) <= tolerance) {
            if (fabs(axial - aa) <= tolerance) {
                if (ua < minimum) minimum = ua;
                if (ua > maximum) maximum = ua;
                if (ub < minimum) minimum = ub;
                if (ub > maximum) maximum = ub;
                intersections += 2;
            }
            continue;
        }
        double lo = aa < ab ? aa : ab;
        double hi = aa < ab ? ab : aa;
        if (axial < lo - tolerance || axial > hi + tolerance) continue;
        double t = (axial - aa) / dv;
        if (t < 0.0) t = 0.0;
        if (t > 1.0) t = 1.0;
        double u = ua + t * (ub - ua);
        if (u < minimum) minimum = u;
        if (u > maximum) maximum = u;
        intersections++;
    }
    if (intersections == 0 || !isfinite(minimum) || !isfinite(maximum) ||
        minimum > maximum)
        return -1;
    *out_min = minimum;
    *out_max = maximum;
    return 0;
}

static int abcl_face_pair_translation_interval(
    const AtlasBoxcutLayoutProblem *p, size_t face0, size_t face1,
    double margin, double *out_left, double *out_right)
{
    double min_v0 = DBL_MAX, max_v0 = -DBL_MAX;
    double min_v1 = DBL_MAX, max_v1 = -DBL_MAX;
    double candidate[8];
    size_t ncandidate = 0;
    for (int k = 0; k < 3; k++) {
        double value0 = p->v[p->faces[face0 * 3 + (size_t)k]];
        double value1 = p->v[p->faces[face1 * 3 + (size_t)k]];
        if (value0 < min_v0) min_v0 = value0;
        if (value0 > max_v0) max_v0 = value0;
        if (value1 < min_v1) min_v1 = value1;
        if (value1 > max_v1) max_v1 = value1;
    }
    double common_min = fmax(min_v0, min_v1);
    double common_max = fmin(max_v0, max_v1);
    double vscale = fmax(1.0, fmax(fabs(common_min), fabs(common_max)));
    if (common_min > common_max + 1.0e-12 * vscale) return -1;
    candidate[ncandidate++] = common_min;
    candidate[ncandidate++] = common_max;
    for (int k = 0; k < 3; k++) {
        double value0 = p->v[p->faces[face0 * 3 + (size_t)k]];
        double value1 = p->v[p->faces[face1 * 3 + (size_t)k]];
        if (value0 > common_min && value0 < common_max)
            candidate[ncandidate++] = value0;
        if (value1 > common_min && value1 < common_max)
            candidate[ncandidate++] = value1;
    }

    double left = DBL_MAX, right = -DBL_MAX;
    for (size_t i = 0; i < ncandidate; i++) {
        double min0, max0, min1, max1;
        if (abcl_triangle_u_slice(p, face0, candidate[i],
                                  &min0, &max0) != 0 ||
            abcl_triangle_u_slice(p, face1, candidate[i],
                                  &min1, &max1) != 0)
            return -1;
        double local_left = min0 - max1;
        double local_right = max0 - min1;
        if (local_left < left) left = local_left;
        if (local_right > right) right = local_right;
    }
    left -= margin;
    right += margin;
    if (!isfinite(left) || !isfinite(right) || !(left < right)) return -1;
    *out_left = left;
    *out_right = right;
    return 0;
}

static int abcl_collect_collisions(
    Arena_T arena, const AtlasBoxcutLayoutProblem *p,
    const AtlasOverlapAudit *audit, const size_t *chart_rank,
    double average_edge, AbclCollision **out_collision,
    size_t *out_ncollision)
{
    if (arena == NULL || p == NULL || audit == NULL || chart_rank == NULL ||
        out_collision == NULL || out_ncollision == NULL ||
        !audit->broad_phase_complete || audit->indexed_faces != p->nfaces)
        return -1;
    AbclCollision *collision = (AbclCollision *)ARENA_ALLOC(
        arena, (audit->npairs ? audit->npairs : 1) *
                   sizeof(AbclCollision));
    size_t nobservation = 0;
    for (size_t i = 0; i < audit->npairs; i++) {
        const AtlasOverlapPair *pair = &audit->pairs[i];
        if (pair->face0 < 0 || pair->face1 < 0 ||
            (size_t)pair->face0 >= p->nfaces ||
            (size_t)pair->face1 >= p->nfaces)
            return -1;
        if (abcl_has_u64(
                p->allowed_overlap_key, p->allowed_overlap_keys,
                abcl_face_pair_key(pair->face0, pair->face1)))
            continue;
        size_t face0 = (size_t)pair->face0;
        size_t face1 = (size_t)pair->face1;
        int32_t chart0 = p->face_label[face0];
        int32_t chart1 = p->face_label[face1];
        if (chart0 < 0 || chart1 < 0 ||
            (size_t)chart0 >= p->ncharts ||
            (size_t)chart1 >= p->ncharts)
            return -1;
        if (chart0 == chart1) continue;
        int chart0_is_inner = chart_rank[chart0] < chart_rank[chart1];
        int32_t inner_chart = chart0_is_inner ? chart0 : chart1;
        int32_t outer_chart = chart0_is_inner ? chart1 : chart0;
        size_t inner_face = chart0_is_inner ? face0 : face1;
        size_t outer_face = chart0_is_inner ? face1 : face0;
        double inner_min, inner_max, outer_min, outer_max;
        abcl_face_u_range(p, inner_face, &inner_min, &inner_max);
        abcl_face_u_range(p, outer_face, &outer_min, &outer_max);
        double lower = inner_max - outer_min + 0.05 * average_edge;
        /* A negative lower bound is still active information.  For example,
         * a relation currently pulled to -100 may collide even though -5 is
         * sufficient separation in the preserved initial ordering.  Dropping
         * that row made the lazy exact loop declare no updates while exact
         * overlaps remained. */
        if (!isfinite(lower)) continue;
        collision[nobservation].inner_chart = inner_chart;
        collision[nobservation].outer_chart = outer_chart;
        collision[nobservation].lower = lower;
        nobservation++;
        (void)inner_min;
        (void)outer_max;
    }
    qsort(collision, nobservation, sizeof(AbclCollision),
          abcl_compare_collision);
    size_t ncollision = 0;
    for (size_t first = 0; first < nobservation;) {
        size_t last = first + 1;
        double lower = collision[first].lower;
        while (last < nobservation &&
               collision[last].inner_chart == collision[first].inner_chart &&
               collision[last].outer_chart == collision[first].outer_chart) {
            if (collision[last].lower > lower) lower = collision[last].lower;
            last++;
        }
        collision[ncollision].inner_chart = collision[first].inner_chart;
        collision[ncollision].outer_chart = collision[first].outer_chart;
        collision[ncollision].lower = lower;
        ncollision++;
        first = last;
    }
    *out_collision = collision;
    *out_ncollision = ncollision;
    return 0;
}

static int abcl_merge_collisions(
    Arena_T arena, const AbclCollision *old_collision, size_t nold,
    const AbclCollision *new_collision, size_t nnew,
    AbclCollision **out_collision, size_t *out_count,
    size_t *out_updates)
{
    if (arena == NULL || out_collision == NULL || out_count == NULL ||
        out_updates == NULL || nold > SIZE_MAX - nnew)
        return -1;
    AbclCollision *merged = (AbclCollision *)ARENA_ALLOC(
        arena, (nold + nnew ? nold + nnew : 1) * sizeof(AbclCollision));
    size_t i = 0, j = 0, count = 0, updates = 0;
    while (i < nold || j < nnew) {
        int take_old = 0, take_new = 0;
        if (j >= nnew) take_old = 1;
        else if (i >= nold) take_new = 1;
        if (!take_old && !take_new) {
            if (old_collision[i].inner_chart != new_collision[j].inner_chart) {
                take_old = old_collision[i].inner_chart <
                           new_collision[j].inner_chart;
                take_new = !take_old;
            } else if (old_collision[i].outer_chart !=
                     new_collision[j].outer_chart) {
                take_old = old_collision[i].outer_chart <
                           new_collision[j].outer_chart;
                take_new = !take_old;
            }
        }
        if (take_old) {
            merged[count++] = old_collision[i++];
        } else if (take_new) {
            merged[count++] = new_collision[j++];
            updates++;
        } else {
            merged[count] = old_collision[i];
            if (new_collision[j].lower > merged[count].lower + 1.0e-9) {
                merged[count].lower = new_collision[j].lower;
                updates++;
            }
            count++;
            i++;
            j++;
        }
    }
    *out_collision = merged;
    *out_count = count;
    *out_updates = updates;
    return 0;
}

static int abcl_slab_index(double value, double origin, double slab_height,
                           int64_t *out)
{
    long double scaled = ((long double)value - (long double)origin) /
                         (long double)slab_height;
    if (!isfinite(scaled) || scaled < 0.0L ||
        scaled > (long double)INT64_MAX)
        return -1;
    *out = (int64_t)floorl(scaled);
    return 0;
}

static int abcl_collect_slab_collisions(
    Arena_T arena, const AtlasBoxcutLayoutProblem *p,
    const size_t *chart_rank,
    double average_edge, AbclCollision **out_collision,
    size_t *out_ncollision, size_t *out_records,
    size_t *out_slabs, double *out_slab_height)
{
    if (arena == NULL || p == NULL || chart_rank == NULL ||
        out_collision == NULL ||
        out_ncollision == NULL || out_records == NULL ||
        out_slabs == NULL || out_slab_height == NULL ||
        p->nfaces > SIZE_MAX / 32)
        return -1;
    double vmin = DBL_MAX, vmax = -DBL_MAX;
    for (size_t f = 0; f < p->nfaces; f++) {
        for (int k = 0; k < 3; k++) {
            int32_t vertex = p->faces[f * 3 + (size_t)k];
            double value = p->v[vertex];
            if (!isfinite(value)) return -1;
            if (value < vmin) vmin = value;
            if (value > vmax) vmax = value;
        }
    }
    double slab_height = fmax(8.0, 2.0 * average_edge);
    int64_t *slab0 = (int64_t *)ARENA_ALLOC(
        arena, p->nfaces * sizeof(int64_t));
    int64_t *slab1 = (int64_t *)ARENA_ALLOC(
        arena, p->nfaces * sizeof(int64_t));
    size_t record_count = 0;
    int ready = 0;
    size_t record_target = 32 * p->nfaces;
    for (int attempt = 0; attempt < 64; attempt++) {
        long double estimate = 0.0L;
        int overflow = 0;
        for (size_t f = 0; f < p->nfaces; f++) {
            double face_vmin = DBL_MAX, face_vmax = -DBL_MAX;
            for (int k = 0; k < 3; k++) {
                int32_t vertex = p->faces[f * 3 + (size_t)k];
                double value = p->v[vertex];
                if (value < face_vmin) face_vmin = value;
                if (value > face_vmax) face_vmax = value;
            }
            if (abcl_slab_index(face_vmin, vmin, slab_height,
                                &slab0[f]) != 0 ||
                abcl_slab_index(face_vmax, vmin, slab_height,
                                &slab1[f]) != 0) {
                overflow = 1;
                break;
            }
            estimate += (long double)slab1[f] -
                        (long double)slab0[f] + 1.0L;
        }
        if (!overflow && estimate <= (long double)record_target) {
            record_count = (size_t)estimate;
            ready = 1;
            break;
        }
        long double growth = 2.0L;
        if (overflow) {
            long double span = (long double)vmax - (long double)vmin;
            long double required = span / ((long double)INT64_MAX / 4.0L);
            if (required > (long double)slab_height)
                growth = 2.0L * required / (long double)slab_height;
        } else if (estimate > 0.0L) {
            growth = estimate / (long double)record_target;
            if (growth < 2.0L) growth = 2.0L;
        }
        long double next = (long double)slab_height * growth;
        if (!isfinite(next)) return -1;
        slab_height = (double)next;
        if (!isfinite(slab_height) || !(slab_height > 0.0)) return -1;
    }
    if (!ready || record_count > SIZE_MAX / sizeof(AbclSlabRecord))
        return -1;

    AbclSlabRecord *record = (AbclSlabRecord *)ARENA_ALLOC(
        arena, (record_count ? record_count : 1) *
                   sizeof(AbclSlabRecord));
    size_t cursor = 0;
    for (size_t f = 0; f < p->nfaces; f++) {
        int32_t chart = p->face_label[f];
        double face_umin, face_umax;
        abcl_face_u_range(p, f, &face_umin, &face_umax);
        int64_t slab = slab0[f];
        for (;;) {
            if (cursor >= record_count) return -1;
            record[cursor].slab = slab;
            record[cursor].chart = chart;
            record[cursor].radial_rank = chart_rank[chart];
            record[cursor].umin = face_umin;
            record[cursor].umax = face_umax;
            cursor++;
            if (slab == slab1[f]) break;
            slab++;
        }
    }
    if (cursor != record_count) return -1;
    qsort(record, record_count, sizeof(AbclSlabRecord),
          abcl_compare_slab_record);

    size_t nenvelope = 0;
    for (size_t first = 0; first < record_count;) {
        size_t last = first + 1;
        double umin = record[first].umin, umax = record[first].umax;
        while (last < record_count &&
               record[last].slab == record[first].slab &&
               record[last].chart == record[first].chart) {
            if (record[last].umin < umin) umin = record[last].umin;
            if (record[last].umax > umax) umax = record[last].umax;
            last++;
        }
        record[nenvelope] = record[first];
        record[nenvelope].umin = umin;
        record[nenvelope].umax = umax;
        nenvelope++;
        first = last;
    }

    AbclCollision *collision = (AbclCollision *)ARENA_ALLOC(
        arena, (nenvelope ? nenvelope : 1) * sizeof(AbclCollision));
    size_t ncollision_observation = 0, nslab = 0;
    for (size_t first = 0; first < nenvelope;) {
        size_t last = first + 1;
        while (last < nenvelope && record[last].slab == record[first].slab)
            last++;
        nslab++;
        for (size_t i = first; i + 1 < last; i++) {
            size_t j = i + 1;
            double lower = record[i].umax - record[j].umin +
                           0.05 * average_edge;
            /* Keep already-satisfied (negative) bounds too.  They are what
             * prevent another active constraint from pulling this ordered
             * neighbour pair through one another later in the solve. */
            if (!isfinite(lower)) continue;
            collision[ncollision_observation].inner_chart = record[i].chart;
            collision[ncollision_observation].outer_chart = record[j].chart;
            collision[ncollision_observation].lower = lower;
            ncollision_observation++;
        }
        first = last;
    }
    qsort(collision, ncollision_observation, sizeof(AbclCollision),
          abcl_compare_collision);
    size_t ncollision = 0;
    for (size_t first = 0; first < ncollision_observation;) {
        size_t last = first + 1;
        double lower = collision[first].lower;
        while (last < ncollision_observation &&
               collision[last].inner_chart == collision[first].inner_chart &&
               collision[last].outer_chart == collision[first].outer_chart) {
            if (collision[last].lower > lower) lower = collision[last].lower;
            last++;
        }
        collision[ncollision] = collision[first];
        collision[ncollision].lower = lower;
        ncollision++;
        first = last;
    }
    *out_collision = collision;
    *out_ncollision = ncollision;
    *out_records = nenvelope;
    *out_slabs = nslab;
    *out_slab_height = slab_height;
    return 0;
}

static int abcl_collect_coherence(
    Arena_T arena, const AtlasBoxcutLayoutProblem *p,
    AtlasBoxcutLayoutCoherence **out_coherence, size_t *out_count,
    AtlasBoxcutLayoutStats *stats)
{
    if (arena == NULL || p == NULL || out_coherence == NULL ||
        out_count == NULL || stats == NULL ||
        p->adjacency_edges < p->intrinsic_adjacency_edges ||
        (p->adjacency_edges > 0 &&
         (p->adjacency_face0 == NULL || p->adjacency_face1 == NULL ||
          p->adjacency_weight == NULL)))
        return -1;

    uint64_t *overlap_key = (uint64_t *)ARENA_ALLOC(
        arena, (p->overlap->npairs ? p->overlap->npairs : 1) *
                   sizeof(uint64_t));
    size_t noverlap_key = 0;
    for (size_t i = 0; i < p->overlap->npairs; i++) {
        int32_t f0 = p->overlap->pairs[i].face0;
        int32_t f1 = p->overlap->pairs[i].face1;
        if (f0 < 0 || f1 < 0 || (size_t)f0 >= p->nfaces ||
            (size_t)f1 >= p->nfaces)
            return -1;
        uint64_t key = abcl_face_pair_key(f0, f1);
        if (abcl_has_u64(p->allowed_overlap_key,
                         p->allowed_overlap_keys, key))
            continue;
        overlap_key[noverlap_key++] = key;
    }
    qsort(overlap_key, noverlap_key, sizeof(uint64_t),
          abcl_compare_u64);

    AbclCoherenceSource *source = (AbclCoherenceSource *)ARENA_ALLOC(
        arena, (p->adjacency_edges ? p->adjacency_edges : 1) *
                   sizeof(AbclCoherenceSource));
    size_t nsource = 0;
    for (size_t i = 0; i < p->adjacency_edges; i++) {
        int32_t f0 = p->adjacency_face0[i];
        int32_t f1 = p->adjacency_face1[i];
        if (f0 < 0 || f1 < 0 || (size_t)f0 >= p->nfaces ||
            (size_t)f1 >= p->nfaces ||
            !isfinite(p->adjacency_weight[i]) ||
            !(p->adjacency_weight[i] > 0.0))
            return -1;
        int32_t c0 = p->face_label[f0];
        int32_t c1 = p->face_label[f1];
        if (c0 < 0 || c1 < 0 || (size_t)c0 >= p->ncharts ||
            (size_t)c1 >= p->ncharts)
            return -1;
        if (c0 == c1) continue;

        /* An adjacent pair that has positive-area overlap is a measured local
         * fold, not evidence that its two sides should remain registered.
         * This exclusion uses only the parameterization audit; the XYZ/BPA
         * weld remains completely held out from the optimizer. */
        if (abcl_has_u64(overlap_key, noverlap_key,
                         abcl_face_pair_key(f0, f1))) {
            stats->coherence_direct_overlap_edges_excluded++;
            continue;
        }
        AbclCoherenceSource *s = &source[nsource++];
        s->chart0 = c0 < c1 ? c0 : c1;
        s->chart1 = c0 < c1 ? c1 : c0;
        s->weight = p->adjacency_weight[i];
        s->seam = i >= p->intrinsic_adjacency_edges;
    }
    qsort(source, nsource, sizeof(AbclCoherenceSource),
          abcl_compare_coherence_source);

    AtlasBoxcutLayoutCoherence *coherence =
        (AtlasBoxcutLayoutCoherence *)ARENA_ALLOC(
            arena, (nsource ? nsource : 1) *
                       sizeof(AtlasBoxcutLayoutCoherence));
    size_t ncoherence = 0;
    for (size_t first = 0; first < nsource;) {
        size_t last = first + 1;
        while (last < nsource &&
               source[last].chart0 == source[first].chart0 &&
               source[last].chart1 == source[first].chart1)
            last++;
        AtlasBoxcutLayoutCoherence *r = &coherence[ncoherence++];
        memset(r, 0, sizeof(*r));
        r->chart0 = source[first].chart0;
        r->chart1 = source[first].chart1;
        for (size_t i = first; i < last; i++) {
            r->source_weight += source[i].weight;
            if (source[i].seam) r->seam_edges++;
            else r->intrinsic_edges++;
        }
        r->qp_weight = p->coherence_weight_scale * r->source_weight;
        if (!isfinite(r->qp_weight) || !(r->qp_weight > 0.0)) return -1;
        stats->coherence_intrinsic_edges += r->intrinsic_edges;
        stats->coherence_seam_edges += r->seam_edges;
        first = last;
    }
    stats->coherence_adjacency_edges = nsource;
    stats->coherence_relations = ncoherence;
    *out_coherence = coherence;
    *out_count = ncoherence;
    return 0;
}

int AtlasBoxcutLayout_collect_disjunctions(
    Arena_T arena, const AtlasBoxcutLayoutProblem *p,
    const AtlasOverlapAudit *overlap,
    AtlasBoxcutLayoutDisjunction **out_disjunction, size_t *out_count)
{
    if (arena == NULL || p == NULL || out_disjunction == NULL ||
        out_count == NULL || overlap == NULL ||
        !isfinite(p->average_uv_edge_length) ||
        !(p->average_uv_edge_length > 0.0) ||
        !overlap->broad_phase_complete || overlap->indexed_faces != p->nfaces)
        return -1;
    AbclDisjunctionSource *source =
        (AbclDisjunctionSource *)ARENA_ALLOC(
            arena, (overlap->npairs ? overlap->npairs : 1) *
                       sizeof(AbclDisjunctionSource));
    size_t nsource = 0;
    double margin = 0.05 * p->average_uv_edge_length;
    for (size_t i = 0; i < overlap->npairs; i++) {
        int32_t f0 = overlap->pairs[i].face0;
        int32_t f1 = overlap->pairs[i].face1;
        if (f0 < 0 || f1 < 0 || (size_t)f0 >= p->nfaces ||
            (size_t)f1 >= p->nfaces)
            return -1;
        if (abcl_has_u64(
                p->allowed_overlap_key, p->allowed_overlap_keys,
                abcl_face_pair_key(f0, f1)))
            continue;
        int32_t c0 = p->face_label[f0];
        int32_t c1 = p->face_label[f1];
        if (c0 < 0 || c1 < 0 || (size_t)c0 >= p->ncharts ||
            (size_t)c1 >= p->ncharts)
            return -1;
        if (c0 == c1) continue;
        size_t face0 = (size_t)f0, face1 = (size_t)f1;
        if (c1 < c0) {
            int32_t swap_chart = c0;
            c0 = c1;
            c1 = swap_chart;
            size_t swap_face = face0;
            face0 = face1;
            face1 = swap_face;
        }
        AbclDisjunctionSource *s = &source[nsource++];
        s->chart0 = c0;
        s->chart1 = c1;
        /* With delta = shift[chart1] - shift[chart0], horizontal
         * bounding-box separation is guaranteed on either side by:
         *   delta <= left_upper  OR  delta >= right_lower.
         * These are deliberately exported as a disjunction; choosing one
         * side prematurely is exactly the fixed-order failure under study. */
        if (abcl_face_pair_translation_interval(
                p, face0, face1, margin,
                &s->left_upper, &s->right_lower) != 0)
            return -1;
    }
    qsort(source, nsource, sizeof(AbclDisjunctionSource),
          abcl_compare_disjunction_source);
    AtlasBoxcutLayoutDisjunction *disjunction =
        (AtlasBoxcutLayoutDisjunction *)ARENA_ALLOC(
            arena, (nsource ? nsource : 1) *
                       sizeof(AtlasBoxcutLayoutDisjunction));
    size_t count = 0;
    for (size_t group_first = 0; group_first < nsource;) {
        size_t group_last = group_first + 1;
        while (group_last < nsource &&
               source[group_last].chart0 == source[group_first].chart0 &&
               source[group_last].chart1 == source[group_first].chart1)
            group_last++;

        /* Preserve the connected components of the forbidden-delta union.
         * Collapsing every interval for a chart pair to one convex hull forces
         * a global left/right order and destroys valid interleavings.  Each
         * disjoint union component is instead its own overlap block and gets
         * an independent binary choice in the MILP prototype. */
        size_t first = group_first;
        while (first < group_last) {
            double left_upper = source[first].left_upper;
            double right_lower = source[first].right_lower;
            size_t last = first + 1;
            while (last < group_last) {
                double scale = fmax(1.0, fmax(fabs(right_lower),
                                             fabs(source[last].left_upper)));
                if (source[last].left_upper >
                    right_lower + 1.0e-12 * scale)
                    break;
                if (source[last].right_lower > right_lower)
                    right_lower = source[last].right_lower;
                last++;
            }
            AtlasBoxcutLayoutDisjunction *d = &disjunction[count++];
            d->chart0 = source[first].chart0;
            d->chart1 = source[first].chart1;
            d->overlap_pairs = last - first;
            d->chart1_left_upper = left_upper;
            d->chart1_right_lower = right_lower;
            first = last;
        }
        group_first = group_last;
    }
    *out_disjunction = disjunction;
    *out_count = count;
    return 0;
}

int AtlasBoxcutLayout_solve(
    Arena_T arena, const AtlasBoxcutLayoutProblem *p,
    double **out_chart_shift,
    AtlasBoxcutLayoutRelation **out_relation, size_t *out_nrelation,
    AtlasBoxcutLayoutCoherence **out_coherence, size_t *out_ncoherence,
    AtlasBoxcutLayoutDisjunction **out_disjunction,
    size_t *out_ndisjunction,
    AtlasBoxcutLayoutStats *stats)
{
    if (arena == NULL || p == NULL || p->faces == NULL || p->nfaces == 0 ||
        p->nfaces > (size_t)INT32_MAX || p->nvertices == 0 ||
        p->base_u == NULL || p->v == NULL || p->registered_u == NULL ||
        p->phi == NULL ||
        p->face_radius == NULL || p->face_label == NULL ||
        p->ncharts == 0 || p->ncharts > (size_t)INT32_MAX ||
        p->overlap == NULL || !isfinite(p->radius_gate) ||
        p->radius_gate < 0.0 ||
        !isfinite(p->average_uv_edge_length) ||
        p->average_uv_edge_length <= 0.0 ||
        !isfinite(p->coherence_weight_scale) ||
        p->coherence_weight_scale <= 0.0 || out_chart_shift == NULL ||
        out_relation == NULL || out_nrelation == NULL ||
        out_coherence == NULL || out_ncoherence == NULL ||
        out_disjunction == NULL || out_ndisjunction == NULL || stats == NULL)
        return -1;
    memset(stats, 0, sizeof(*stats));
    *out_chart_shift = NULL;
    *out_relation = NULL;
    *out_nrelation = 0;
    *out_coherence = NULL;
    *out_ncoherence = 0;
    *out_disjunction = NULL;
    *out_ndisjunction = 0;

    size_t *chart_faces = (size_t *)ARENA_CALLOC(
        arena, p->ncharts, sizeof(size_t));
    double *chart_radius = (double *)ARENA_CALLOC(
        arena, p->ncharts, sizeof(double));
    double *chart_base_u = (double *)ARENA_CALLOC(
        arena, p->ncharts, sizeof(double));
    for (size_t f = 0; f < p->nfaces; f++) {
        int32_t chart = p->face_label[f];
        if (chart < 0 || (size_t)chart >= p->ncharts) return -1;
        if (!isfinite(p->face_radius[f])) return -1;
        chart_faces[chart]++;
        chart_radius[chart] += p->face_radius[f];
        chart_base_u[chart] += abcl_face_mean_double(p, f, p->base_u);
        for (int k = 0; k < 3; k++) {
            int32_t v = p->faces[f * 3 + (size_t)k];
            if (v < 0 || (size_t)v >= p->nvertices) return -1;
        }
    }
    for (size_t chart = 0; chart < p->ncharts; chart++) {
        if (chart_faces[chart] == 0) return -1;
        chart_radius[chart] /= (double)chart_faces[chart];
        chart_base_u[chart] /= (double)chart_faces[chart];
    }

    AbclObservation *observation = (AbclObservation *)ARENA_ALLOC(
        arena, (p->overlap->npairs ? p->overlap->npairs : 1) *
                   sizeof(AbclObservation));
    size_t nobservation = 0;
    const double two_pi = 6.283185307179586476925286766559;
    const double average_edge = p->average_uv_edge_length;
    for (size_t i = 0; i < p->overlap->npairs; i++) {
        const AtlasOverlapPair *pair = &p->overlap->pairs[i];
        if (pair->face0 < 0 || pair->face1 < 0 ||
            (size_t)pair->face0 >= p->nfaces ||
            (size_t)pair->face1 >= p->nfaces)
            return -1;
        if (abcl_has_u64(
                p->allowed_overlap_key, p->allowed_overlap_keys,
                abcl_face_pair_key(pair->face0, pair->face1)))
            continue;
        size_t f0 = (size_t)pair->face0, f1 = (size_t)pair->face1;
        double r0 = p->face_radius[f0], r1 = p->face_radius[f1];
        if (!isfinite(r0) || !isfinite(r1)) return -1;
        if (fabs(r1 - r0) <= p->radius_gate) continue;
        int32_t l0 = p->face_label[f0], l1 = p->face_label[f1];
        if (l0 == l1) continue;

        int first_is_inner = r0 < r1;
        size_t inner_face = first_is_inner ? f0 : f1;
        size_t outer_face = first_is_inner ? f1 : f0;
        int32_t inner_chart = first_is_inner ? l0 : l1;
        int32_t outer_chart = first_is_inner ? l1 : l0;
        double base_du = abcl_face_mean_double(p, outer_face, p->base_u) -
                         abcl_face_mean_double(p, inner_face, p->base_u);
        double desired_du = fabs(
            abcl_face_mean_float(p, outer_face, p->registered_u) -
            abcl_face_mean_float(p, inner_face, p->registered_u));
        double predicted_du = fabs((double)pair->turn) * two_pi *
                              0.5 * (r0 + r1);
        if (desired_du < 4.0 * average_edge && predicted_du > desired_du)
            desired_du = predicted_du;
        double correction = desired_du - base_du;
        double clearance = 2.0 * average_edge - base_du;
        if (correction < clearance) correction = clearance;
        if (!isfinite(correction) || correction <= 0.0) continue;

        double inner_min, inner_max, outer_min, outer_max;
        abcl_face_u_range(p, inner_face, &inner_min, &inner_max);
        abcl_face_u_range(p, outer_face, &outer_min, &outer_max);
        double separation = inner_max - outer_min + 0.05 * average_edge;
        if (!isfinite(separation) || separation <= 0.0) continue;
        (void)inner_min;
        (void)outer_max;

        int32_t chart0 = l0 < l1 ? l0 : l1;
        int32_t chart1 = l0 < l1 ? l1 : l0;
        AbclObservation *o = &observation[nobservation++];
        o->chart0 = chart0;
        o->chart1 = chart1;
        o->target = outer_chart == chart1 ? correction : -correction;
        o->separation = outer_chart == chart1 ? separation : -separation;
        (void)inner_chart;
    }
    if (nobservation == 0) return -1;
    qsort(observation, nobservation, sizeof(AbclObservation),
          abcl_compare_observation);

    AtlasBoxcutLayoutRelation *relation =
        (AtlasBoxcutLayoutRelation *)ARENA_ALLOC(
            arena, nobservation * sizeof(AtlasBoxcutLayoutRelation));
    double *scratch = (double *)ARENA_ALLOC(
        arena, nobservation * sizeof(double));
    size_t nrelation = 0;
    for (size_t first = 0; first < nobservation;) {
        size_t last = first + 1;
        while (last < nobservation &&
               observation[last].chart0 == observation[first].chart0 &&
               observation[last].chart1 == observation[first].chart1)
            last++;
        AtlasBoxcutLayoutRelation *r = &relation[nrelation++];
        memset(r, 0, sizeof(*r));
        r->chart0 = observation[first].chart0;
        r->chart1 = observation[first].chart1;
        r->pairs = last - first;
        size_t count = 0;
        for (size_t i = first; i < last; i++) {
            scratch[count++] = observation[i].target;
            if (observation[i].target > 0.0) r->positive_pairs++;
            else r->negative_pairs++;
        }
        r->target_median = abcl_median(scratch, count);
        count = 0;
        for (size_t i = first; i < last; i++)
            scratch[count++] = fabs(observation[i].target -
                                      r->target_median);
        r->target_mad = abcl_median(scratch, count);
        size_t dominant = r->positive_pairs > r->negative_pairs
            ? r->positive_pairs : r->negative_pairs;
        r->sign_agreement = (double)dominant / (double)r->pairs;
        int dominant_positive = r->positive_pairs >= r->negative_pairs;
        r->separation_lower = 0.0;
        for (size_t i = first; i < last; i++) {
            if ((observation[i].target > 0.0) != dominant_positive) continue;
            double lower = fabs(observation[i].separation);
            if (lower > r->separation_lower) r->separation_lower = lower;
        }
        double mad_limit = fmax(8.0 * average_edge,
                                0.35 * fabs(r->target_median));
        r->eligible = r->sign_agreement >= 0.90 &&
                      fabs(r->target_median) >= average_edge &&
                      r->target_mad <= mad_limit &&
                      r->separation_lower > 0.0;
        if (r->sign_agreement < 0.90)
            stats->mixed_direction_relations++;
        first = last;
    }

    AbclChartRank *radial_order = (AbclChartRank *)ARENA_ALLOC(
        arena, p->ncharts * sizeof(AbclChartRank));
    size_t *radial_rank = (size_t *)ARENA_ALLOC(
        arena, p->ncharts * sizeof(size_t));
    AbclChartRank *placement_order = (AbclChartRank *)ARENA_ALLOC(
        arena, p->ncharts * sizeof(AbclChartRank));
    size_t *placement_rank = (size_t *)ARENA_ALLOC(
        arena, p->ncharts * sizeof(size_t));
    for (size_t chart = 0; chart < p->ncharts; chart++) {
        radial_order[chart].chart = (int32_t)chart;
        radial_order[chart].radius = chart_radius[chart];
        placement_order[chart].chart = (int32_t)chart;
        placement_order[chart].radius = chart_base_u[chart];
    }
    qsort(radial_order, p->ncharts, sizeof(AbclChartRank),
          abcl_compare_chart_rank);
    qsort(placement_order, p->ncharts, sizeof(AbclChartRank),
          abcl_compare_chart_rank);
    for (size_t i = 0; i < p->ncharts; i++) {
        radial_rank[radial_order[i].chart] = i;
        placement_rank[placement_order[i].chart] = i;
    }

    AbclCandidate *candidate = (AbclCandidate *)ARENA_ALLOC(
        arena, (nrelation ? nrelation : 1) * sizeof(AbclCandidate));
    size_t ncandidate = 0;
    for (size_t i = 0; i < nrelation; i++) {
        AtlasBoxcutLayoutRelation *r = &relation[i];
        if (!r->eligible) continue;
        int target_positive = r->target_median > 0.0;
        int chart1_is_outer = radial_rank[r->chart1] > radial_rank[r->chart0];
        if (target_positive != chart1_is_outer) {
            r->eligible = 0;
            stats->collision_relations_rejected++;
            continue;
        }
        double scale = fmax(average_edge, fabs(r->target_median));
        candidate[ncandidate].relation = i;
        candidate[ncandidate].pairs = r->pairs;
        candidate[ncandidate].score =
            (double)r->pairs * r->sign_agreement /
            (1.0 + r->target_mad / scale);
        r->qp_weight = candidate[ncandidate].score;
        r->selected = 1;
        ncandidate++;
    }
    qsort(candidate, ncandidate, sizeof(AbclCandidate),
          abcl_compare_candidate);

    AtlasBoxcutLayoutCoherence *coherence = NULL;
    size_t ncoherence = 0;
    if (abcl_collect_coherence(
            arena, p, &coherence, &ncoherence, stats) != 0)
        return -1;
    AtlasBoxcutLayoutDisjunction *disjunction = NULL;
    size_t ndisjunction = 0;
    if (AtlasBoxcutLayout_collect_disjunctions(
            arena, p, p->overlap, &disjunction, &ndisjunction) != 0)
        return -1;
    stats->collision_disjunctions = ndisjunction;

    AbclCollision *collision = NULL;
    size_t ncollision = 0;
    if (abcl_collect_collisions(
            arena, p, p->overlap, placement_rank, average_edge,
            &collision, &ncollision) != 0)
        return -1;
    AbclCollision *slab_collision = NULL;
    size_t nslab_collision = 0;
    /* The slab envelope is a useful sufficient-condition diagnostic, but it
     * is not a valid default placement constraint: serializing every chart
     * present in an axial slab tears physically adjacent pieces and creates a
     * scroll-width cascade.  The complete exact audit below supplies the
     * necessary collision rows lazily instead. */
    AbclCollision *seed_collision = NULL;
    size_t nseed_collision = 0, seed_updates = 0;
    if (abcl_merge_collisions(
            arena, collision, ncollision, slab_collision, nslab_collision,
            &seed_collision, &nseed_collision, &seed_updates) != 0)
        return -1;
    collision = seed_collision;
    ncollision = nseed_collision;
    stats->slab_bounds = nslab_collision;
    stats->collision_bounds_added = seed_updates;

    if (ncandidate > SIZE_MAX / 2 || ncoherence > SIZE_MAX / 2 ||
        p->ncharts > SIZE_MAX - ncandidate ||
        p->ncharts + ncandidate > SIZE_MAX - ncoherence ||
        p->ncharts > SIZE_MAX - 2 * ncandidate ||
        p->ncharts + 2 * ncandidate > SIZE_MAX - 2 * ncoherence)
        return -1;
    size_t nrows = p->ncharts + ncandidate + ncoherence;
    size_t ncoeff = p->ncharts + 2 * ncandidate + 2 * ncoherence;
    MonotoneQpRow *row = (MonotoneQpRow *)ARENA_ALLOC(
        arena, nrows * sizeof(MonotoneQpRow));
    MonotoneQpCoeff *coeff = (MonotoneQpCoeff *)ARENA_ALLOC(
        arena, ncoeff * sizeof(MonotoneQpCoeff));
    double average_chart_faces = (double)p->nfaces / (double)p->ncharts;
    size_t coeff_cursor = 0;
    for (size_t chart = 0; chart < p->ncharts; chart++) {
        row[chart].first = coeff_cursor;
        row[chart].count = 1;
        row[chart].target = 0.0;
        row[chart].weight = 0.25 * (double)chart_faces[chart] /
                            average_chart_faces;
        row[chart].kind = 0;
        row[chart].owner = (int32_t)chart;
        coeff[coeff_cursor].var = (int32_t)chart;
        coeff[coeff_cursor].value = 1.0;
        coeff_cursor++;
    }
    for (size_t i = 0; i < ncandidate; i++) {
        size_t ri = candidate[i].relation;
        AtlasBoxcutLayoutRelation *r = &relation[ri];
        MonotoneQpRow *rr = &row[p->ncharts + i];
        rr->first = coeff_cursor;
        rr->count = 2;
        rr->target = r->target_median;
        rr->weight = fmax(1.0e-6, r->qp_weight);
        rr->kind = 1;
        rr->owner = ri <= (size_t)INT32_MAX ? (int32_t)ri : INT32_MAX;
        coeff[coeff_cursor].var = r->chart0;
        coeff[coeff_cursor].value = -1.0;
        coeff_cursor++;
        coeff[coeff_cursor].var = r->chart1;
        coeff[coeff_cursor].value = 1.0;
        coeff_cursor++;
    }
    for (size_t i = 0; i < ncoherence; i++) {
        AtlasBoxcutLayoutCoherence *r = &coherence[i];
        MonotoneQpRow *rr = &row[p->ncharts + ncandidate + i];
        rr->first = coeff_cursor;
        rr->count = 2;
        rr->target = 0.0;
        rr->weight = r->qp_weight;
        rr->kind = 2;
        rr->owner = i <= (size_t)INT32_MAX ? (int32_t)i : INT32_MAX;
        coeff[coeff_cursor].var = r->chart0;
        coeff[coeff_cursor].value = -1.0;
        coeff_cursor++;
        coeff[coeff_cursor].var = r->chart1;
        coeff[coeff_cursor].value = 1.0;
        coeff_cursor++;
    }
    if (coeff_cursor != ncoeff ||
        p->nfaces > (size_t)INT32_MAX / 3) return -1;

    size_t ncorner = 3 * p->nfaces;
    int32_t *audit_faces = (int32_t *)ARENA_ALLOC(
        arena, ncorner * sizeof(int32_t));
    double *audit_u = (double *)ARENA_ALLOC(
        arena, ncorner * sizeof(double));
    double *audit_v = (double *)ARENA_ALLOC(
        arena, ncorner * sizeof(double));
    float *audit_registered = (float *)ARENA_ALLOC(
        arena, ncorner * sizeof(float));
    float *audit_phi = (float *)ARENA_ALLOC(
        arena, ncorner * sizeof(float));
    int32_t *audit_component = (int32_t *)ARENA_ALLOC(
        arena, ncorner * sizeof(int32_t));
    for (size_t f = 0; f < p->nfaces; f++) {
        int32_t chart = p->face_label[f];
        for (int k = 0; k < 3; k++) {
            size_t corner = f * 3 + (size_t)k;
            int32_t source = p->faces[corner];
            audit_faces[corner] = (int32_t)corner;
            audit_v[corner] = p->v[source];
            audit_registered[corner] = p->registered_u[source];
            audit_phi[corner] = p->phi[source];
            audit_component[corner] = chart;
        }
    }

    stats->initial_exact_overlap_pairs = p->overlap->exact_face_pairs;
    for (size_t i = 0; i < p->overlap->npairs; i++) {
        const AtlasOverlapPair *pair = &p->overlap->pairs[i];
        if (p->face_label[pair->face0] != p->face_label[pair->face1])
            stats->initial_cross_overlap_pairs++;
    }

    double *shift = (double *)ARENA_ALLOC(
        arena, p->ncharts * sizeof(double));
    int solved = 0;
    const int maximum_outer_iterations = 128;
    for (int outer_iteration = 0;
         outer_iteration < maximum_outer_iterations; outer_iteration++) {
        Arena_T work = Arena_new();
        MonotoneQpBound *bound = (MonotoneQpBound *)ARENA_ALLOC(
            work, (ncollision ? ncollision : 1) * sizeof(MonotoneQpBound));
        for (size_t i = 0; i < ncollision; i++) {
            bound[i].lo = collision[i].inner_chart;
            bound[i].hi = collision[i].outer_chart;
            bound[i].lower = collision[i].lower;
            bound[i].stroke = -1;
            bound[i].edge = i <= (size_t)INT32_MAX
                          ? (int32_t)i : INT32_MAX;
        }

        memset(shift, 0, p->ncharts * sizeof(double));
        double initial_margin = 1.0e-4 * fmax(1.0, average_edge);
        for (size_t rank_index = 0; rank_index < p->ncharts; rank_index++) {
            int32_t outer = placement_order[rank_index].chart;
            for (size_t e = 0; e < ncollision; e++) {
                if (bound[e].hi != outer) continue;
                double minimum = shift[bound[e].lo] + bound[e].lower +
                                 initial_margin;
                if (shift[outer] < minimum) shift[outer] = minimum;
            }
        }
        for (size_t e = 0; e < ncollision; e++) {
            if (shift[bound[e].hi] - shift[bound[e].lo] <
                bound[e].lower - 1.0e-8) {
                fprintf(stderr,
                        "[atlas_boxcut_layout] infeasible initializer: "
                        "outer=%d bounds=%zu\n",
                        outer_iteration, ncollision);
                Arena_dispose(&work);
                return -1;
            }
        }

        MonotoneQpProblem qp;
        memset(&qp, 0, sizeof(qp));
        qp.nvar = p->ncharts;
        qp.rows = row;
        qp.nrows = nrows;
        qp.coeff = coeff;
        qp.ncoeff = ncoeff;
        qp.bounds = bound;
        qp.nbounds = ncollision;
        MonotoneQpOptions qp_options;
        MonotoneQpOptions_default(&qp_options);
        size_t iteration_base = p->ncharts + ncollision;
        size_t iteration_limit = iteration_base > (SIZE_MAX - 64) / 4
            ? SIZE_MAX : 4 * iteration_base + 64;
        qp_options.max_active_iterations = iteration_limit > (size_t)INT_MAX
            ? INT_MAX : (int)iteration_limit;
        MonotoneQpStats qp_stats;
        int qp_rc = MonotoneQp_solve(
            work, &qp, &qp_options, shift, NULL, &qp_stats);
        if (qp_rc != 0) {
            fprintf(stderr,
                    "[atlas_boxcut_layout] QP failed: outer=%d rc=%d "
                    "bounds=%zu limit=%d spd=%d multiplier=%d active=%d\n",
                    outer_iteration, qp_rc, ncollision,
                    qp_options.max_active_iterations, qp_stats.spd_solves,
                    qp_stats.multiplier_solves, qp_stats.active_final);
            Arena_dispose(&work);
            return -1;
        }
        solved = 1;
        stats->collision_outer_iterations = (size_t)outer_iteration + 1;
        stats->collision_bounds = ncollision;
        if (qp_stats.iterations > INT_MAX - stats->qp_iterations)
            stats->qp_iterations = INT_MAX;
        else
            stats->qp_iterations += qp_stats.iterations;
        stats->qp_active_bounds = qp_stats.active_final;
        stats->qp_objective = qp_stats.objective_final;
        stats->qp_stationarity = qp_stats.stationarity_residual;
        stats->qp_max_violation = qp_stats.max_violation;

        for (size_t f = 0; f < p->nfaces; f++) {
            int32_t chart = p->face_label[f];
            for (int k = 0; k < 3; k++) {
                size_t corner = f * 3 + (size_t)k;
                int32_t source = p->faces[corner];
                audit_u[corner] = p->base_u[source] + shift[chart];
            }
        }
        AtlasOverlapAudit current_overlap;
        if (AtlasOverlapAudit_build(
                work, audit_faces, p->nfaces, ncorner, audit_u, audit_v,
                audit_registered, audit_phi, audit_component, p->ncharts,
                &current_overlap) != 0 ||
            !current_overlap.broad_phase_complete ||
            current_overlap.indexed_faces != p->nfaces) {
            Arena_dispose(&work);
            return -1;
        }
        stats->final_exact_overlap_pairs = current_overlap.exact_face_pairs;
        stats->final_cross_overlap_pairs =
            current_overlap.cross_component_pairs;
        if (current_overlap.cross_component_pairs == 0) {
            Arena_dispose(&work);
            break;
        }

        AbclCollision *new_collision = NULL;
        size_t nnew_collision = 0;
        if (abcl_collect_collisions(
                work, p, &current_overlap, placement_rank, average_edge,
                &new_collision, &nnew_collision) != 0) {
            Arena_dispose(&work);
            return -1;
        }
        AbclCollision *merged = NULL;
        size_t nmerged = 0, updates = 0;
        if (abcl_merge_collisions(
                arena, collision, ncollision, new_collision, nnew_collision,
                &merged, &nmerged, &updates) != 0) {
            Arena_dispose(&work);
            return -1;
        }
        Arena_dispose(&work);
        if (updates == 0) break;
        collision = merged;
        ncollision = nmerged;
        stats->collision_bounds_added += updates;
    }
    if (!solved) return -1;

    double relation_residual2 = 0.0;
    for (size_t i = 0; i < nrelation; i++) {
        AtlasBoxcutLayoutRelation *r = &relation[i];
        r->solved_delta = shift[r->chart1] - shift[r->chart0];
        r->residual = r->solved_delta - r->target_median;
        relation_residual2 += r->residual * r->residual;
        if (fabs(r->residual) > stats->relation_residual_max)
            stats->relation_residual_max = fabs(r->residual);
    }
    stats->relation_residual_rms = nrelation
        ? sqrt(relation_residual2 / (double)nrelation) : 0.0;

    double coherence_residual2 = 0.0;
    double coherence_weighted_residual2 = 0.0;
    double coherence_weight_sum = 0.0;
    double *coherence_scratch = (double *)ARENA_ALLOC(
        arena, (ncoherence ? ncoherence : 1) * sizeof(double));
    for (size_t i = 0; i < ncoherence; i++) {
        AtlasBoxcutLayoutCoherence *r = &coherence[i];
        r->solved_delta = shift[r->chart1] - shift[r->chart0];
        r->residual = r->solved_delta;
        double absolute = fabs(r->residual);
        coherence_scratch[i] = absolute;
        coherence_residual2 += r->residual * r->residual;
        coherence_weighted_residual2 +=
            r->source_weight * r->residual * r->residual;
        coherence_weight_sum += r->source_weight;
        if (absolute <= average_edge) {
            stats->coherence_relations_preserved++;
            stats->coherence_edges_preserved +=
                r->intrinsic_edges + r->seam_edges;
        }
        if (absolute > stats->coherence_residual_max)
            stats->coherence_residual_max = absolute;
    }
    stats->coherence_residual_rms = ncoherence
        ? sqrt(coherence_residual2 / (double)ncoherence) : 0.0;
    stats->coherence_weighted_residual_rms = coherence_weight_sum > 0.0
        ? sqrt(coherence_weighted_residual2 / coherence_weight_sum) : 0.0;
    if (ncoherence > 0) {
        qsort(coherence_scratch, ncoherence, sizeof(double),
              abcl_compare_double);
        stats->coherence_residual_median = ncoherence & 1u
            ? coherence_scratch[ncoherence / 2]
            : 0.5 * (coherence_scratch[ncoherence / 2 - 1] +
                     coherence_scratch[ncoherence / 2]);
        size_t p95 = ncoherence - ncoherence / 20;
        stats->coherence_residual_p95 = coherence_scratch[p95 - 1];
    }

    for (size_t i = 0; i < nobservation; i++) {
        const AbclObservation *o = &observation[i];
        double solved_delta = shift[o->chart1] - shift[o->chart0];
        if (solved_delta * o->target > 0.0)
            stats->radial_order_pairs_satisfied++;
        if ((o->target > 0.0 &&
             solved_delta >= o->target - average_edge) ||
            (o->target < 0.0 &&
             solved_delta <= o->target + average_edge))
            stats->target_pairs_satisfied++;
    }
    for (size_t i = 0; i < p->ncharts; i++) {
        if (fabs(shift[i]) > 1e-10) {
            stats->moved_charts++;
            stats->moved_faces += chart_faces[i];
        }
        if (fabs(shift[i]) > stats->maximum_abs_shift)
            stats->maximum_abs_shift = fabs(shift[i]);
    }
    stats->observations = nobservation;
    stats->relations = nrelation;
    stats->eligible_relations = ncandidate;
    stats->selected_relations = ncandidate;
    stats->radial_order_pair_fraction =
        (double)stats->radial_order_pairs_satisfied / (double)nobservation;
    stats->target_pair_fraction =
        (double)stats->target_pairs_satisfied / (double)nobservation;

    *out_chart_shift = shift;
    *out_relation = relation;
    *out_nrelation = nrelation;
    *out_coherence = coherence;
    *out_ncoherence = ncoherence;
    *out_disjunction = disjunction;
    *out_ndisjunction = ndisjunction;
    return 0;
}

int AtlasBoxcutLayout_selftest(void)
{
    const int32_t faces[9] = {
        0, 1, 2,
        3, 4, 5,
        6, 7, 8
    };
    const double u[9] = {
        0.0, 2.0, 0.0,
        10.0, 12.0, 10.0,
        0.0, 2.0, 0.0
    };
    const double v[9] = {
        0.0, 0.0, 2.0,
        0.0, 0.0, 2.0,
        3.0, 3.0, 5.0
    };
    AtlasBoxcutLayoutProblem problem;
    memset(&problem, 0, sizeof(problem));
    problem.faces = faces;
    problem.nfaces = 3;
    problem.nvertices = 9;
    problem.base_u = u;
    problem.v = v;
    double left = 0.0, right = 0.0;
    if (abcl_face_pair_translation_interval(
            &problem, 0, 1, 0.0, &left, &right) != 0 ||
        fabs(left + 12.0) > 1.0e-12 || fabs(right + 8.0) > 1.0e-12)
        return 1;
    if (abcl_face_pair_translation_interval(
            &problem, 0, 1, 0.25, &left, &right) != 0 ||
        fabs(left + 12.25) > 1.0e-12 || fabs(right + 7.75) > 1.0e-12)
        return 1;
    if (abcl_face_pair_translation_interval(
            &problem, 0, 2, 0.0, &left, &right) == 0)
        return 1;
    return 0;
}
